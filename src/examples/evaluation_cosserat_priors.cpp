// Evaluation driver: Cosserat prior vs. baseline estimator on an S-shape under
// an unknown tip-force disturbance.
//
// Pipeline:
//   1. Load YAML config (topology + hyperparams + options).
//   2. Cosserat solve A (prior, f_ext = 0) -> ContinuumRodPriors fed to the estimator.
//   3. Cosserat solve B (real,  f_ext != 0) -> ground-truth disk frames.
//   4. Build the EstimatorPriors package (deterministic, one-shot).
//   5. For each of N trials (different noise seed each), synthesize a noisy
//      tip-position measurement and run the estimator three ways:
//        Method 1 (Baseline)              : noisy tip only, Straight init, no priors
//        Method 2 (Strain as measurement) : noisy tip + Cosserat strain measurements
//                                           (ControlInputMode::None)
//        Method 3 (Force as input)        : noisy tip + acceleration inputs
//                                           -K^-1 [f; l] + junction strain jump
//                                           (ControlInputMode::ForceAsInput)
//   6. Report per-method RMSE / MaxErr / Time stats (mean ± std over N trials).
//   7. Optionally dump the LAST trial's per-node positions to CSV for plotting.
//
// Requires USE_LOCAL_TDCR=ON. See doc/evaluation.md for the experiment
// design, results and open issues.

#include "config_loader.h"
#include "continuum_robot_state_estimator.h"
#include "continuum_rod_priors.h"
#include "cosserat_priors_adapter.h"
#include "cosseratrodmodel.h"

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using Estimator       = ContinuumRobotStateEstimator;
using SensorMeas      = Estimator::SensorMeasurement;
using ControlInputT   = Estimator::ControlInput;
using SystemStateT    = Estimator::SystemState;

namespace {

// Cosserat body (z-forward) -> estimator body (x-forward).
// Same cyclic permutation the adapter uses internally; redefined here because
// the adapter keeps it in an anonymous namespace.
const Eigen::Matrix3d& rotConv()
{
    static const Eigen::Matrix3d R = (Eigen::Matrix3d() <<
        0, 0, 1,
        1, 0, 0,
        0, 1, 0).finished();
    return R;
}

Eigen::Matrix4d nearestDiskFrame(double s_k,
                                 const Eigen::VectorXd& s,
                                 const Eigen::MatrixXd& diskFrames)
{
    const Eigen::Index N = s.size();
    Eigen::Index best = 0;
    double best_dist = std::abs(s(0) - s_k);
    for (Eigen::Index j = 1; j < N; ++j) {
        const double d = std::abs(s(j) - s_k);
        if (d < best_dist) { best_dist = d; best = j; }
    }
    return diskFrames.block(0, 4 * best, 4, 4);
}

Eigen::Matrix4d cosseratToInertial(const Eigen::Matrix4d& T_cos,
                                   const Eigen::Matrix4d& Ti0)
{
    const Eigen::Matrix3d& R_conv = rotConv();
    const Eigen::Matrix3d  R_convT = R_conv.transpose();

    Eigen::Matrix4d T_body = Eigen::Matrix4d::Identity();
    T_body.block(0, 0, 3, 3) = R_conv * T_cos.block(0, 0, 3, 3) * R_convT;
    T_body.block(0, 3, 3, 1) = R_conv * T_cos.block(0, 3, 3, 1);
    return Ti0 * T_body;
}

std::vector<Eigen::Vector3d> cosseratPositionsAtNodes(
    const Eigen::MatrixXd& diskFrames,
    const Eigen::VectorXd& s_disk,
    const Eigen::Matrix4d& Ti0,
    unsigned int K,
    double L)
{
    std::vector<Eigen::Vector3d> p(K);
    for (unsigned int k = 0; k < K; ++k) {
        const double s_k = (static_cast<double>(k) / (K - 1)) * L;
        const Eigen::Matrix4d T_cos = nearestDiskFrame(s_k, s_disk, diskFrames);
        const Eigen::Matrix4d T_in  = cosseratToInertial(T_cos, Ti0);
        p[k] = T_in.block(0, 3, 3, 1);
    }
    return p;
}

std::pair<double, double> positionErrors(
    const SystemStateT& state,
    const std::vector<Eigen::Vector3d>& p_gt,
    unsigned int robot_idx)
{
    const auto& nodes = state.robots.at(robot_idx).estimation_nodes;
    const std::size_t K = std::min<std::size_t>(nodes.size(), p_gt.size());
    double sse = 0.0;
    double max_err = 0.0;
    for (std::size_t k = 0; k < K; ++k) {
        const Eigen::Vector3d p_est = nodes[k].pose.block(0, 3, 3, 1);
        const double e = (p_est - p_gt[k]).norm();
        sse += e * e;
        if (e > max_err) max_err = e;
    }
    const double rmse = (K > 0) ? std::sqrt(sse / static_cast<double>(K)) : 0.0;
    return { rmse, max_err };
}

// Per-node position error [m] between an estimate and the ground truth.
std::vector<double> nodeErrors(
    const SystemStateT& state,
    const std::vector<Eigen::Vector3d>& p_gt,
    unsigned int robot_idx)
{
    const auto& nodes = state.robots.at(robot_idx).estimation_nodes;
    const std::size_t K = std::min<std::size_t>(nodes.size(), p_gt.size());
    std::vector<double> e(K);
    for (std::size_t k = 0; k < K; ++k)
        e[k] = (Eigen::Vector3d(nodes[k].pose.block(0, 3, 3, 1)) - p_gt[k]).norm();
    return e;
}

std::string fmtSci(double v, int width = 12, int prec = 4)
{
    std::ostringstream os;
    os << std::setw(width) << std::scientific << std::setprecision(prec) << v;
    return os.str();
}

std::string fmtFix(double v, int width = 10, int prec = 2)
{
    std::ostringstream os;
    os << std::setw(width) << std::fixed << std::setprecision(prec) << v;
    return os.str();
}

void printBanner(const std::string& title)
{
    std::cout << "============================================================\n  " << title
              << "\n============================================================\n";
}

// "(   161.2,   -0.0,  -64.9) mm" — fixed-point mm tuple for tip positions.
std::string fmtVec3Mm(const Eigen::Vector3d& v, int width = 7, int prec = 1)
{
    std::ostringstream os;
    os << std::fixed << std::setprecision(prec) << "("
       << std::setw(width) << v(0) * 1000.0 << ", "
       << std::setw(width) << v(1) * 1000.0 << ", "
       << std::setw(width) << v(2) * 1000.0 << ") mm";
    return os.str();
}

// Format a (mean, std) pair as e.g. "  14.92 ±   1.40" right-aligned to a
// target visual width. We do the padding here instead of via std::setw at the
// call site because std::setw counts bytes, while the ± character (U+00B1) is
// 2 bytes in UTF-8 but renders as 1 visual column — so std::setw would
// silently misalign by 1 char in every row that contains ±.
//
// The core "{6w} ± {6w}" format always produces exactly 15 visual columns
// (6 + space + ± + space + 6); we pad with leading spaces to the target.
std::string fmtMmStd(double mean_m, double std_m, int visual_width = 16, int prec = 2)
{
    std::ostringstream os;
    os << std::fixed << std::setprecision(prec)
       << std::setw(6) << mean_m * 1000.0 << " ± "
       << std::setw(6) << std_m  * 1000.0;
    std::string s = os.str();
    const int pad = std::max(0, visual_width - 15);
    return std::string(pad, ' ') + s;
}
std::string fmtMsStd(double mean_ms, double std_ms, int visual_width = 16, int prec = 2)
{
    std::ostringstream os;
    os << std::fixed << std::setprecision(prec)
       << std::setw(6) << mean_ms << " ± "
       << std::setw(6) << std_ms;
    std::string s = os.str();
    const int pad = std::max(0, visual_width - 15);
    return std::string(pad, ' ') + s;
}

// Backbone RMSE of the (no-force) Cosserat prior shape against the (with-force)
// ground-truth shape. This is the "what if we used the prior directly, with no
// estimator at all?" baseline that all three methods are competing against.
double priorTruthRmse(const std::vector<Eigen::Vector3d>& p_prior,
                      const std::vector<Eigen::Vector3d>& p_gt)
{
    if (p_prior.empty() || p_prior.size() != p_gt.size()) return 0.0;
    double sse = 0.0;
    for (std::size_t k = 0; k < p_prior.size(); ++k) {
        const double e = (p_prior[k] - p_gt[k]).norm();
        sse += e * e;
    }
    return std::sqrt(sse / static_cast<double>(p_prior.size()));
}

struct Args
{
    std::string  config_path;
    std::string  csv_path = "";              // if non-empty, dump last trial here
    std::string  trials_csv_path = "";       // if non-empty, append one row per trial and method
    unsigned int n_trials = 1;               // multi-seed averaging
    double       sigma   = 0.001;            // tip-position noise std [m]
    double       fx      = 0.0;              // tip force [N], base frame
    double       fy      = 0.0;
    double       fz      = 0.3;              // disturbs the tip ~43 mm in z
    unsigned int seed = 42;
    // Tendon tensions tuned for a balanced S-shape in the XZ plane:
    //   seg-1 pulls tendon 0          (curves rod upward in +z, first half of S)
    //   seg-2 pulls tendons 1 and 2   (curves rod downward in -z, second half)
    // Visually classic S — segment-1 bend and segment-2 curl are similar size.
    double q0 = 6.0, q1 = 0.0, q2 = 0.0;
    double q3 = 0.0, q4 = 4.0, q5 = 4.0;
};

bool parseArgs(int argc, char* argv[], Args& a)
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0]
                  << " <config.yaml>"
                  << " [--n-trials N] [--csv path] [--trials-csv path]"
                  << " [--sigma S] [--fx X] [--fy Y] [--fz Z] [--seed N]"
                  << " [--q F1,F2,F3,F4,F5,F6]\n";
        return false;
    }
    a.config_path = argv[1];
    for (int i = 2; i < argc; ++i) {
        const std::string s = argv[i];
        auto nextDbl = [&](const char* flag, double& dst) {
            if (i + 1 >= argc) {
                std::cerr << flag << " requires a value\n";
                std::exit(1);
            }
            dst = std::atof(argv[++i]);
        };
        auto nextStr = [&](const char* flag, std::string& dst) {
            if (i + 1 >= argc) {
                std::cerr << flag << " requires a value\n";
                std::exit(1);
            }
            dst = argv[++i];
        };
        if      (s == "--sigma")    nextDbl("--sigma",  a.sigma);
        else if (s == "--fx")       nextDbl("--fx",     a.fx);
        else if (s == "--fy")       nextDbl("--fy",     a.fy);
        else if (s == "--fz")       nextDbl("--fz",     a.fz);
        else if (s == "--csv")          nextStr("--csv",         a.csv_path);
        else if (s == "--trials-csv")   nextStr("--trials-csv",  a.trials_csv_path);
        else if (s == "--n-trials") {
            if (i + 1 >= argc) { std::cerr << "--n-trials requires a value\n"; return false; }
            a.n_trials = static_cast<unsigned int>(std::atoi(argv[++i]));
            if (a.n_trials == 0) { std::cerr << "--n-trials must be >= 1\n"; return false; }
        }
        else if (s == "--seed") {
            if (i + 1 >= argc) { std::cerr << "--seed requires a value\n"; return false; }
            a.seed = static_cast<unsigned int>(std::atoi(argv[++i]));
        }
        else if (s == "--q") {
            if (i + 1 >= argc) { std::cerr << "--q requires 6 comma-separated values\n"; return false; }
            std::string vals = argv[++i];
            for (char& c : vals) if (c == ',') c = ' ';
            std::istringstream is(vals);
            if (!(is >> a.q0 >> a.q1 >> a.q2 >> a.q3 >> a.q4 >> a.q5)) {
                std::cerr << "--q must be 6 comma-separated numbers\n"; return false;
            }
        }
        else {
            std::cerr << "Unknown flag: " << s << "\n"; return false;
        }
    }
    return true;
}

// ----- Per-trial result --------------------------------------------------------

struct TrialMetrics
{
    bool        ok_old, ok_strain, ok_force;
    double      rmse_old, rmse_strain, rmse_force;
    double      max_old,  max_strain,  max_force;
    double      time_old, time_strain, time_force;
    std::size_t iters_old, iters_strain, iters_force;
    unsigned int seed;
    std::vector<double> err_old, err_strain, err_force;   // per-node position error [m]
};

struct TrialState
{
    Eigen::Matrix4d T_tip_noisy;
    SystemStateT    state_old, state_strain, state_force;
};

auto clock_now() { return std::chrono::steady_clock::now(); }
double ms_between(std::chrono::steady_clock::time_point a,
                  std::chrono::steady_clock::time_point b)
{
    return std::chrono::duration<double, std::milli>(b - a).count();
}

TrialMetrics runTrial(
    const Estimator::RobotTopology& topology,
    const Estimator::Hyperparameters& params,
    const Estimator::Options& base_opts_robust,
    const EstimatorPriors& ep_strain,
    const EstimatorPriors& ep_force,
    const Eigen::MatrixXd& diskFrames_real,
    const Eigen::VectorXd& s_disk_real,
    const Eigen::Matrix4d& Ti0,
    const std::vector<Eigen::Vector3d>& p_gt,
    unsigned int robot_idx,
    unsigned int K,
    double L,
    double sigma,
    unsigned int seed,
    TrialState* opt_state_out)
{
    std::mt19937 rng(seed);
    std::normal_distribution<double> noise(0.0, sigma);

    // Noisy tip measurement (position only, mask = [1,1,1,0,0,0])
    const Eigen::Matrix4d T_cos_tip = nearestDiskFrame(L, s_disk_real, diskFrames_real);
    Eigen::Matrix4d T_tip = cosseratToInertial(T_cos_tip, Ti0);
    T_tip(0, 3) += noise(rng);
    T_tip(1, 3) += noise(rng);
    T_tip(2, 3) += noise(rng);

    SensorMeas tip_meas;
    tip_meas.type      = SensorMeas::Pose;
    tip_meas.value     = T_tip;
    tip_meas.mask      << 1, 1, 1, 0, 0, 0;
    tip_meas.idx_robot = robot_idx;
    tip_meas.idx_node  = K - 1;

    // Method 1: Baseline -- noisy tip only, Straight init, no priors
    auto opts_old = base_opts_robust;
    opts_old.init_guess_type = Estimator::Options::Straight;
    Estimator est_old(topology, params, opts_old);
    SystemStateT state_old;
    std::vector<double> cost_old;
    const auto t1_0 = clock_now();
    const bool ok_old = est_old.computeStateEstimate(
        state_old, cost_old, { tip_meas }, {}, false);
    const double t_old = ms_between(t1_0, clock_now());

    // Method 2: Strain as measurement -- noisy tip + strain measurements, Cosserat init
    auto opts_strain = base_opts_robust;
    opts_strain.init_guess_type = Estimator::Options::Custom;
    opts_strain.custom_guess    = ep_strain.initial_guess;
    std::vector<SensorMeas> meas_strain = ep_strain.measurements;
    meas_strain.push_back(tip_meas);
    Estimator est_strain(topology, params, opts_strain);
    SystemStateT state_strain;
    std::vector<double> cost_strain;
    const auto t2_0 = clock_now();
    const bool ok_strain = est_strain.computeStateEstimate(
        state_strain, cost_strain, meas_strain, {}, false);
    const double t_strain = ms_between(t2_0, clock_now());

    // Method 3: Force as input -- noisy tip + load-driven strain derivative
    // -K^-1 [f; l] plus the junction strain jump on the acceleration channel,
    // Cosserat init.
    auto opts_force = base_opts_robust;
    opts_force.init_guess_type = Estimator::Options::Custom;
    opts_force.custom_guess    = ep_force.initial_guess;
    Estimator est_force(topology, params, opts_force);
    SystemStateT state_force;
    std::vector<double> cost_force;
    const auto t3_0 = clock_now();
    const bool ok_force = est_force.computeStateEstimate(
        state_force, cost_force, { tip_meas }, ep_force.control_inputs, false);
    const double t_force = ms_between(t3_0, clock_now());

    const auto [rmse_old,    max_old]    = positionErrors(state_old,    p_gt, robot_idx);
    const auto [rmse_strain, max_strain] = positionErrors(state_strain, p_gt, robot_idx);
    const auto [rmse_force,  max_force]  = positionErrors(state_force,  p_gt, robot_idx);

    if (opt_state_out) {
        opt_state_out->T_tip_noisy  = T_tip;
        opt_state_out->state_old    = state_old;
        opt_state_out->state_strain = state_strain;
        opt_state_out->state_force  = state_force;
    }

    TrialMetrics m;
    m.ok_old     = ok_old;     m.ok_strain     = ok_strain;     m.ok_force     = ok_force;
    m.rmse_old   = rmse_old;   m.rmse_strain   = rmse_strain;   m.rmse_force   = rmse_force;
    m.max_old    = max_old;    m.max_strain    = max_strain;    m.max_force    = max_force;
    m.time_old   = t_old;      m.time_strain   = t_strain;      m.time_force   = t_force;
    m.iters_old  = cost_old.size();
    m.iters_strain = cost_strain.size();
    m.iters_force  = cost_force.size();
    m.seed       = seed;
    m.err_old    = nodeErrors(state_old,    p_gt, robot_idx);
    m.err_strain = nodeErrors(state_strain, p_gt, robot_idx);
    m.err_force  = nodeErrors(state_force,  p_gt, robot_idx);
    return m;
}

// ----- Aggregation -------------------------------------------------------------

struct Stat { double mean, std; };

Stat meanStd(const std::vector<double>& v)
{
    if (v.empty()) return { 0.0, 0.0 };
    double s = 0.0;
    for (double x : v) s += x;
    const double m = s / static_cast<double>(v.size());
    double sq = 0.0;
    for (double x : v) sq += (x - m) * (x - m);
    const double sd = (v.size() > 1)
                          ? std::sqrt(sq / static_cast<double>(v.size() - 1))
                          : 0.0;
    return { m, sd };
}

// ----- CSV ---------------------------------------------------------------------

void writeCSV(const std::string& path,
              const std::vector<Eigen::Vector3d>& p_prior,
              const std::vector<Eigen::Vector3d>& p_gt,
              const TrialState& last,
              unsigned int robot_idx,
              double L)
{
    const std::size_t K = p_gt.size();
    std::ofstream f(path);
    if (!f) {
        std::cerr << "Failed to open CSV path: " << path << "\n";
        return;
    }
    f << "node,s,"
      << "prior_x,prior_y,prior_z,"
      << "truth_x,truth_y,truth_z,"
      << "m1_old_x,m1_old_y,m1_old_z,"
      << "m2_strain_x,m2_strain_y,m2_strain_z,"
      << "m3_force_x,m3_force_y,m3_force_z\n";

    const auto& nodes_old    = last.state_old.robots.at(robot_idx).estimation_nodes;
    const auto& nodes_strain = last.state_strain.robots.at(robot_idx).estimation_nodes;
    const auto& nodes_force  = last.state_force.robots.at(robot_idx).estimation_nodes;

    f << std::setprecision(8) << std::fixed;
    for (std::size_t k = 0; k < K; ++k) {
        const double s_k = (K > 1) ? (static_cast<double>(k) / (K - 1)) * L : 0.0;
        const Eigen::Vector3d p_p = p_prior[k];
        const Eigen::Vector3d p_t = p_gt[k];
        const Eigen::Vector3d p_1 = nodes_old   [k].pose.block(0, 3, 3, 1);
        const Eigen::Vector3d p_2 = nodes_strain[k].pose.block(0, 3, 3, 1);
        const Eigen::Vector3d p_3 = nodes_force [k].pose.block(0, 3, 3, 1);
        f << k << "," << s_k << ","
          << p_p(0) << "," << p_p(1) << "," << p_p(2) << ","
          << p_t(0) << "," << p_t(1) << "," << p_t(2) << ","
          << p_1(0) << "," << p_1(1) << "," << p_1(2) << ","
          << p_2(0) << "," << p_2(1) << "," << p_2(2) << ","
          << p_3(0) << "," << p_3(1) << "," << p_3(2) << "\n";
    }
    // Append the noisy tip measurement as a special row (node = -1)
    f << "-1,"     // node
      << L  << "," // arclength = tip
      << "0,0,0," // prior placeholder
      << "0,0,0," // truth placeholder
      << last.T_tip_noisy(0, 3) << "," << last.T_tip_noisy(1, 3) << "," << last.T_tip_noisy(2, 3) << ","
      << "0,0,0,"
      << "0,0,0\n";

    std::cout << "Wrote CSV: " << path << " (" << K << " nodes + 1 tip-measurement row)\n";
}

} // anonymous namespace


int main(int argc, char* argv[])
{
    Args args;
    if (!parseArgs(argc, argv, args)) return 1;

    // --- 1. Load YAML config -------------------------------------------------
    ConfigLoader config(args.config_path);
    const auto topology     = config.getTopology();
    const auto params       = config.getHyperparameters();
    const auto base_options = config.getOptions();

    const unsigned int robot_idx = 0;
    const unsigned int K         = topology.K.at(robot_idx);
    const double       L         = topology.L.at(robot_idx);
    const Eigen::Matrix4d Ti0    = topology.Ti0.at(robot_idx);

    printBanner("Cosserat-prior evaluation  -  S-shape + tip-force disturbance");
    std::cout << "\n"
              << "Scenario\n"
              << "  Config               : " << args.config_path << "\n"
              << "  Tendon tensions      : q = ["
              << args.q0 << ", " << args.q1 << ", " << args.q2 << " | "
              << args.q3 << ", " << args.q4 << ", " << args.q5 << "] N\n"
              << "  Tip-force disturbance: F = ("
              << std::fixed << std::setprecision(2)
              << args.fx << ", " << args.fy << ", " << args.fz << ") N\n"
              << "  Sensor noise std     : sigma = "
              << std::fixed << std::setprecision(2) << args.sigma * 1000.0 << " mm\n"
              << "  Noise seeds          : " << args.n_trials << " trials"
              << "  (seed_start = " << args.seed << ")\n"
              << "  Estimation nodes     : K = " << K << " along L = "
              << std::fixed << std::setprecision(1) << L * 1000.0 << " mm rod\n\n";

    Eigen::Matrix<double, 6, 1> q;
    q << args.q0, args.q1, args.q2, args.q3, args.q4, args.q5;
    const Eigen::Vector3d F_disturb(args.fx, args.fy, args.fz);
    const Eigen::Vector3d L_disturb = Eigen::Vector3d::Zero();
    const double L1 = 0.1, L2 = 0.1;

    // --- 2. Cosserat solves --------------------------------------------------
    CosseratRodModel model_prior;
    Eigen::MatrixXd diskFrames_prior;
    if (!model_prior.forwardKinematics(diskFrames_prior, q,
                                       Eigen::Vector3d::Zero(),
                                       Eigen::Vector3d::Zero())
        || !model_prior.hasAuxOutputs()) {
        std::cerr << "Cosserat prior solve FAILED (residual = "
                  << model_prior.getFinalResudial() << ")\n";
        return 1;
    }
    const ContinuumRodPriors priors = priorsFromCosseratModel(model_prior, L1, L2, diskFrames_prior);
    const Eigen::VectorXd s_disk_prior = model_prior.getArclengthSamples();

    CosseratRodModel model_real;
    // Warm-start the real solve from the prior's converged shooting state.
    // The default initial guess (straight rod) diverges for tip forces above
    // ~0.02 N when the rod is already heavily bent by tendon tensions; using
    // the prior's converged (u0, v0) as the starting point lets the shooting
    // method converge cleanly for much larger disturbances.
    model_real.setDefaultInitValues(model_prior.getFinalInitValues());
    Eigen::MatrixXd diskFrames_real;
    if (!model_real.forwardKinematics(diskFrames_real, q, F_disturb, L_disturb)
        || !model_real.hasAuxOutputs()) {
        std::cerr << "Cosserat real solve FAILED (residual = "
                  << model_real.getFinalResudial() << ")\n";
        return 1;
    }
    const Eigen::VectorXd s_disk_real = model_real.getArclengthSamples();

    const auto p_prior = cosseratPositionsAtNodes(diskFrames_prior, s_disk_prior, Ti0, K, L);
    const auto p_gt    = cosseratPositionsAtNodes(diskFrames_real,  s_disk_real,  Ti0, K, L);
    const double prior_truth_rmse = priorTruthRmse(p_prior, p_gt);

    // --- 3. Adapter (once each, deterministic) ------------------------------
    // M2 uses the model strain as measurements (mode None: no control inputs);
    // M3 uses only the ForceAsInput control inputs. Initial guesses are
    // identical in the two outputs (velocity input is zero in both).
    const auto t_adapter0 = clock_now();
    const EstimatorPriors ep_strain = cosseratPriorsToEstimator(
        priors, topology, robot_idx, diskFrames_prior,
        ControlInputMode::None);
    const EstimatorPriors ep_force  = cosseratPriorsToEstimator(
        priors, topology, robot_idx, diskFrames_prior,
        ControlInputMode::ForceAsInput);
    const double adapter_ms = ms_between(t_adapter0, clock_now());

    auto base_opts_robust = base_options;
    base_opts_robust.solver = Estimator::Options::NewtonLineSearch;

    // --- 4. Trial loop ------------------------------------------------------
    std::vector<TrialMetrics> trials;
    trials.reserve(args.n_trials);
    TrialState last_state;

    for (unsigned int t = 0; t < args.n_trials; ++t) {
        const unsigned int seed = args.seed + t;
        const bool capture_state = (t == args.n_trials - 1);
        trials.push_back(runTrial(
            topology, params, base_opts_robust, ep_strain, ep_force,
            diskFrames_real, s_disk_real, Ti0, p_gt,
            robot_idx, K, L, args.sigma, seed,
            capture_state ? &last_state : nullptr));
    }

    // --- 5. Aggregate stats -------------------------------------------------
    auto collect = [&](double TrialMetrics::*member) {
        std::vector<double> v;
        v.reserve(trials.size());
        for (const auto& m : trials) v.push_back(m.*member);
        return meanStd(v);
    };
    auto convergedCount = [&](bool TrialMetrics::*member) {
        std::size_t n = 0;
        for (const auto& m : trials) if (m.*member) ++n;
        return n;
    };

    const auto s_rmse_old    = collect(&TrialMetrics::rmse_old);
    const auto s_rmse_strain = collect(&TrialMetrics::rmse_strain);
    const auto s_rmse_force  = collect(&TrialMetrics::rmse_force);
    const auto s_max_old     = collect(&TrialMetrics::max_old);
    const auto s_max_strain  = collect(&TrialMetrics::max_strain);
    const auto s_max_force   = collect(&TrialMetrics::max_force);
    const auto s_time_old    = collect(&TrialMetrics::time_old);
    const auto s_time_strain = collect(&TrialMetrics::time_strain);
    const auto s_time_force  = collect(&TrialMetrics::time_force);

    const std::size_t conv_old    = convergedCount(&TrialMetrics::ok_old);
    const std::size_t conv_strain = convergedCount(&TrialMetrics::ok_strain);
    const std::size_t conv_force  = convergedCount(&TrialMetrics::ok_force);

    // --- 6. Reporting -------------------------------------------------------
    //
    // Layout (each block separated by a blank line):
    //   1. Cosserat ground truth      — tip positions + prior-truth gap.
    //   2. Headline                   — one-paragraph plain-English verdict.
    //   3. Per-method results table   — mm-based, fixed-point.
    //   4. Comparisons                — relative numbers (vs. prior, vs. method 1).
    //
    // The headline is computed AFTER the trials run; it picks winners on each
    // axis (accuracy / speed) among methods that produced a usable shape
    // (RMSE < L/2).

    const double L_mm = L * 1000.0;
    const double half_L_mm = L_mm * 0.5;
    const Stat rmse_means[3] = { s_rmse_old, s_rmse_strain, s_rmse_force };
    const Stat time_means[3] = { s_time_old, s_time_strain, s_time_force };
    const std::size_t conv_counts[3] = { conv_old, conv_strain, conv_force };
    const char* method_names[3] = { "Baseline (tip only)",
                                    "Strain as measurement",
                                    "Force as input" };

    int accuracy_winner = -1, speed_winner = -1;
    double best_rmse_m = std::numeric_limits<double>::infinity();
    double best_time_ms = std::numeric_limits<double>::infinity();
    for (int i = 0; i < 3; ++i) {
        const bool valid = (rmse_means[i].mean * 1000.0 < half_L_mm);
        if (!valid) continue;
        if (rmse_means[i].mean < best_rmse_m) { best_rmse_m  = rmse_means[i].mean; accuracy_winner = i; }
        if (time_means[i].mean < best_time_ms) { best_time_ms = time_means[i].mean; speed_winner    = i; }
    }
    auto isFailed = [&](int i) {
        return rmse_means[i].mean * 1000.0 >= half_L_mm || conv_counts[i] == 0;
    };

    // --- Block 1: Cosserat ground truth ------------------------------------
    std::cout << "Cosserat ground truth\n"
              << "  Prior tip position  : " << fmtVec3Mm(p_prior.back()) << "\n"
              << "  Truth tip position  : " << fmtVec3Mm(p_gt.back())   << "\n"
              << "  Prior<->truth gap   : "
              << std::fixed << std::setprecision(2) << prior_truth_rmse * 1000.0
              << " mm RMSE   (baseline without estimator)\n\n";

    // --- Block 2: Headline -------------------------------------------------
    std::cout << "[*] Headline\n";
    if (accuracy_winner < 0) {
        std::cout << "  No method produced a usable shape (all RMSE > L/2 = "
                  << std::fixed << std::setprecision(0) << half_L_mm << " mm).\n";
    } else if (accuracy_winner == speed_winner) {
        const int w = accuracy_winner;
        std::cout << "  " << method_names[w]
                  << " wins on BOTH accuracy and speed: "
                  << std::fixed << std::setprecision(2)
                  << rmse_means[w].mean * 1000.0 << " mm RMSE in "
                  << time_means[w].mean << " ms.\n";
        // Show the speedup against the slowest other method.
        int slow_i = (w == 0) ? 1 : 0;
        for (int i = 0; i < 3; ++i)
            if (i != w && time_means[i].mean > time_means[slow_i].mean) slow_i = i;
        if (rmse_means[slow_i].mean * 1000.0 < half_L_mm) {
            std::cout << "  That is "
                      << std::fixed << std::setprecision(1)
                      << rmse_means[slow_i].mean / rmse_means[w].mean
                      << "x more accurate and "
                      << std::setprecision(0) << time_means[slow_i].mean / time_means[w].mean
                      << "x faster than " << method_names[slow_i] << ".\n";
        }
    } else {
        std::cout << "  Mixed result (we are past the crossover where the prior stops helping on accuracy):\n"
                  << "    " << method_names[accuracy_winner] << " wins on accuracy: "
                  << std::fixed << std::setprecision(2)
                  << rmse_means[accuracy_winner].mean * 1000.0 << " mm RMSE.\n"
                  << "    " << method_names[speed_winner]    << " wins on speed:    "
                  << std::fixed << std::setprecision(2)
                  << time_means[speed_winner].mean << " ms ("
                  << std::setprecision(0)
                  << time_means[accuracy_winner].mean / time_means[speed_winner].mean
                  << "x faster than " << method_names[accuracy_winner] << ").\n";
    }
    for (int i = 0; i < 3; ++i) {
        if (isFailed(i)) {
            std::cout << "  " << method_names[i] << " FAILS: "
                      << std::fixed << std::setprecision(1)
                      << rmse_means[i].mean * 1000.0 << " mm RMSE";
            if (conv_counts[i] == 0)
                std::cout << ", " << conv_counts[i] << "/" << args.n_trials << " converged";
            std::cout << ".\n";
        }
    }
    std::cout << "\n";

    // --- Block 3: per-method results table ---------------------------------
    //
    // Column layout (every width below is in *visual* columns, not bytes):
    //
    //   "  [mark] [name padded W_NAME]  [RMSE W_VAL]  [MaxErr W_VAL]  [Time W_VAL]  [Conv W_CONV]"
    //      ^^^^                          ^^^^                ^^^^                ^^^^         ^^^^^
    //      W_MARK                         16                 16                 16            8
    //
    // The fmt*Std helpers do their own visual-width padding (they handle the
    // ± byte/visual mismatch internally), so we do NOT wrap them in std::setw
    // at the call site — only header strings (pure ASCII) use std::setw.
    constexpr int W_MARK = 3;
    constexpr int W_NAME = 24;
    constexpr int W_VAL  = 16;     // applies to RMSE, MaxErr, Time columns
    constexpr int W_CONV = 8;
    const std::string SEP = "  ";  // 2-space column separator

    const int line_width =
        2 + W_MARK + 1 + W_NAME
          + static_cast<int>(SEP.size()) + W_VAL
          + static_cast<int>(SEP.size()) + W_VAL
          + static_cast<int>(SEP.size()) + W_VAL
          + static_cast<int>(SEP.size()) + W_CONV;
    const std::string hline = std::string(line_width, '-');

    std::cout << "Per-method results  (mean " << "\xC2\xB1" << " 1 std over "
              << args.n_trials << " trials)\n"
              << hline << "\n"
              << "  " << std::setw(W_MARK) << "" << " "
              << std::left  << std::setw(W_NAME) << "Method"
              << SEP << std::right << std::setw(W_VAL)  << "RMSE [mm]"
              << SEP <<               std::setw(W_VAL)  << "MaxErr [mm]"
              << SEP <<               std::setw(W_VAL)  << "Time [ms]"
              << SEP <<               std::setw(W_CONV) << "Conv"
              << "\n" << hline << "\n";

    auto printRow = [&](int idx, const char* short_name, Stat rmse, Stat maxe, Stat tms, std::size_t conv) {
        // Markers: [*] = winner on at least one axis; [X] = failed/invalid.
        const char* marker = "   ";
        if (isFailed(idx))                                       marker = "[X]";
        else if (idx == accuracy_winner || idx == speed_winner)  marker = "[*]";

        const std::string conv_str = std::to_string(conv) + "/" + std::to_string(args.n_trials);

        std::cout << "  " << std::left  << std::setw(W_MARK) << marker << " "
                  <<        std::setw(W_NAME) << short_name
                  << SEP << fmtMmStd(rmse.mean, rmse.std, W_VAL)
                  << SEP << fmtMmStd(maxe.mean, maxe.std, W_VAL)
                  << SEP << fmtMsStd(tms.mean,  tms.std,  W_VAL)
                  << SEP << std::right << std::setw(W_CONV) << conv_str
                  << "\n";
    };
    printRow(0, "1. Baseline (tip only)",   s_rmse_old,    s_max_old,    s_time_old,    conv_old);
    printRow(1, "2. Strain as measurement", s_rmse_strain, s_max_strain, s_time_strain, conv_strain);
    printRow(2, "3. Force as input",        s_rmse_force,  s_max_force,  s_time_force,  conv_force);
    std::cout << hline << "\n"
              << "  [*] = best on at least one axis     [X] = failed (RMSE > L/2 or 0 converged)\n"
              << "  Adapter prep (one-shot, not included in per-trial timings): "
              << std::fixed << std::setprecision(2) << adapter_ms << " ms\n\n";

    // --- Block 4: comparisons ----------------------------------------------
    std::cout << "Comparisons\n";
    const double gap_mm = prior_truth_rmse * 1000.0;
    std::cout << "  vs. prior" << "\xE2\x86\x94" << "truth gap ("
              << std::fixed << std::setprecision(2)
              << gap_mm << " mm - the no-estimator baseline):\n";
    for (int i = 0; i < 3; ++i) {
        const double rmse_mm = rmse_means[i].mean * 1000.0;
        std::ostringstream tag;
        if (gap_mm > 0 && rmse_mm < gap_mm) {
            const double pct = (gap_mm - rmse_mm) / gap_mm * 100.0;
            tag << std::fixed << std::setprecision(0) << pct << "% better than prior";
        } else if (gap_mm > 0) {
            const double ratio = rmse_mm / gap_mm;
            tag << std::fixed << std::setprecision(1) << ratio << "x WORSE than prior";
        }
        std::cout << "    " << std::left << std::setw(28) << method_names[i] << std::right
                  << " : " << std::fixed << std::setprecision(2)
                  << std::setw(7) << rmse_mm << " mm   (" << tag.str() << ")\n";
    }
    std::cout << "\n"
              << "  Strain as measurement vs Baseline:\n"
              << "    Accuracy: " << std::fixed << std::setprecision(2)
              << s_rmse_strain.mean * 1000.0 << " vs " << s_rmse_old.mean * 1000.0 << " mm";
    if (s_rmse_old.mean > 0) {
        const double pct = (s_rmse_strain.mean - s_rmse_old.mean) / s_rmse_old.mean * 100.0;
        std::cout << "  ("
                  << std::fixed << std::setprecision(0) << std::abs(pct) << "% "
                  << (pct < 0 ? "better" : "worse") << ")";
    }
    std::cout << "\n    Speed:    "
              << std::fixed << std::setprecision(2) << s_time_strain.mean << " vs "
              << s_time_old.mean << " ms   ("
              << std::fixed << std::setprecision(1)
              << s_time_old.mean / s_time_strain.mean << "x faster)\n\n";

    // --- 6a. Optional per-trial results: one row per trial and method, plus
    // one "0_model" row (trial = -1) for the no-force model on its own.
    // Errors in metres; e0..e{K-1} are the per-node position errors.
    if (!args.trials_csv_path.empty()) {
        const bool exists = std::ifstream(args.trials_csv_path).good();
        std::ofstream f(args.trials_csv_path, std::ios::app);
        if (!f) {
            std::cerr << "Failed to open trials-csv path: " << args.trials_csv_path << "\n";
        } else {
            if (!exists) {
                f << "fx,fy,fz,sigma,L,trial,seed,method,rmse,max_err,time_ms,iterations,converged";
                for (unsigned int k = 0; k < K; ++k) f << ",e" << k;
                f << "\n";
            }
            f << std::setprecision(9);
            auto writeTrial = [&](long trial, long seed, const char* method, double rmse,
                                  double max_err, double time_ms, std::size_t iters,
                                  bool converged, const std::vector<double>& err) {
                f << args.fx << "," << args.fy << "," << args.fz << "," << args.sigma << ","
                  << L << "," << trial << "," << seed << "," << method << ","
                  << rmse << "," << max_err << "," << time_ms << "," << iters << ","
                  << (converged ? 1 : 0);
                for (double e : err) f << "," << e;
                f << "\n";
            };
            std::vector<double> err_prior(K);
            for (unsigned int k = 0; k < K; ++k) err_prior[k] = (p_prior[k] - p_gt[k]).norm();
            writeTrial(-1, -1, "0_model", prior_truth_rmse,
                       *std::max_element(err_prior.begin(), err_prior.end()), 0.0, 0, true, err_prior);
            for (std::size_t t = 0; t < trials.size(); ++t) {
                const auto& m = trials[t];
                writeTrial(t, m.seed, "1_old",    m.rmse_old,    m.max_old,    m.time_old,    m.iters_old,    m.ok_old,    m.err_old);
                writeTrial(t, m.seed, "2_strain", m.rmse_strain, m.max_strain, m.time_strain, m.iters_strain, m.ok_strain, m.err_strain);
                writeTrial(t, m.seed, "3_force",  m.rmse_force,  m.max_force,  m.time_force,  m.iters_force,  m.ok_force,  m.err_force);
            }
        }
    }

    // --- 6b. Optional per-trial CSV dump (last trial) ----------------------
    if (!args.csv_path.empty()) {
        writeCSV(args.csv_path, p_prior, p_gt, last_state, robot_idx, L);
    }

    printBanner("Evaluation OK");

    // Exit non-zero if any method failed to converge on the LAST trial — useful
    // for CI checks but harmless for normal use.
    const bool any_fail = !(trials.back().ok_old &&
                            trials.back().ok_strain &&
                            trials.back().ok_force);
    return any_fail ? 2 : 0;
}
