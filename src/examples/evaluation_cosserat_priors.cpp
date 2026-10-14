// Evaluation driver: four ways of giving the Cosserat model to the estimator,
// compared on an S-shaped rod pushed by an unknown tip force.
//
// Pipeline:
//   1. Load YAML config (topology + hyperparams + options).
//   2. Cosserat solve without the force -> model prediction, fed to the estimator.
//   3. Cosserat solve with the force    -> ground truth (hidden from the estimator).
//   4. Turn the model prediction into estimator inputs, once per method.
//   5. For each of N trials (different noise seed each), synthesize a noisy
//      tip-position measurement and run the estimator four ways:
//        Baseline              : noisy tip only, straight start, nothing from the model
//        Strain as measurement : noisy tip + model strain as measurements
//                                (ControlInputMode::None)
//        Force as input        : noisy tip + acceleration inputs -K^-1 [f; l]
//                                + junction strain jump (ControlInputMode::ForceAsInput)
//        Strain as input       : noisy tip + velocity inputs = model strain
//                                (ControlInputMode::StrainAsInput)
//   6. Report per-method RMSE / MaxErr / Time stats (mean ± std over N trials).
//      Time = computeStateEstimate() only; --timing-csv times the stages before it
//      (Cosserat solve, conversion, estimator set-up), see timePipelineStages().
//   7. Optionally write per-run results and the last trial's shapes to CSV.
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
#include <array>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using Estimator    = ContinuumRobotStateEstimator;
using SensorMeas   = Estimator::SensorMeasurement;
using ControlInput = Estimator::ControlInput;
using SystemState  = Estimator::SystemState;
using Positions    = std::vector<Eigen::Vector3d>;   // one position per estimation node

namespace {

// ----- The compared methods ----------------------------------------------------

enum Method { Baseline, StrainMeas, ForceIn, StrainIn, N_METHODS };

const char* kMethodName[N_METHODS] = { "Baseline (tip only)",
                                       "Strain as measurement",
                                       "Force as input",
                                       "Strain as input" };
const char* kTableLabel[N_METHODS] = { "1. Baseline (tip only)",
                                       "2. Strain as measurement",
                                       "3. Force as input",
                                       "4. Strain as input" };
// Method keys in the CSV files (scripts/evaluate.py reads them).
const char* kKey[N_METHODS]        = { "baseline", "strain_meas", "force_in", "strain_in" };

// What the estimator gets from the model for one method: solver options
// (starting shape), model measurements and control inputs. Built once.
struct MethodSetup
{
    Estimator::Options        options;
    std::vector<SensorMeas>   model_measurements;
    std::vector<ControlInput> control_inputs;
};

// ----- Frames and positions ----------------------------------------------------

// Cosserat body axes (z-forward) -> estimator body axes (x-forward).
// Same cyclic permutation the adapter uses internally; redefined here because
// the adapter keeps it in an anonymous namespace.
const Eigen::Matrix3d& cosseratToEstimatorAxes()
{
    static const Eigen::Matrix3d R = (Eigen::Matrix3d() <<
        0, 0, 1,
        1, 0, 0,
        0, 1, 0).finished();
    return R;
}

Eigen::Matrix4d nearestDiskFrame(double s,
                                 const Eigen::VectorXd& s_disks,
                                 const Eigen::MatrixXd& diskFrames)
{
    const Eigen::Index N = s_disks.size();
    Eigen::Index best = 0;
    double best_dist = std::abs(s_disks(0) - s);
    for (Eigen::Index j = 1; j < N; ++j) {
        const double d = std::abs(s_disks(j) - s);
        if (d < best_dist) { best_dist = d; best = j; }
    }
    return diskFrames.block(0, 4 * best, 4, 4);
}

// Cosserat disk frame -> world frame of the estimator (base pose Ti0).
Eigen::Matrix4d cosseratFrameToWorld(const Eigen::Matrix4d& T_cos,
                                     const Eigen::Matrix4d& Ti0)
{
    const Eigen::Matrix3d& R_conv = cosseratToEstimatorAxes();
    const Eigen::Matrix3d  R_convT = R_conv.transpose();

    Eigen::Matrix4d T_body = Eigen::Matrix4d::Identity();
    T_body.block(0, 0, 3, 3) = R_conv * T_cos.block(0, 0, 3, 3) * R_convT;
    T_body.block(0, 3, 3, 1) = R_conv * T_cos.block(0, 3, 3, 1);
    return Ti0 * T_body;
}

// Cosserat shape sampled at the K estimation nodes, in the world frame.
Positions nodePositionsFromCosserat(const Eigen::MatrixXd& diskFrames,
                                    const Eigen::VectorXd& s_disks,
                                    const Eigen::Matrix4d& Ti0,
                                    unsigned int K,
                                    double L)
{
    Positions p(K);
    for (unsigned int k = 0; k < K; ++k) {
        const double s_k = (static_cast<double>(k) / (K - 1)) * L;
        const Eigen::Matrix4d T_cos = nearestDiskFrame(s_k, s_disks, diskFrames);
        p[k] = cosseratFrameToWorld(T_cos, Ti0).block(0, 3, 3, 1);
    }
    return p;
}

Positions nodePositionsOf(const SystemState& state, unsigned int robot_idx)
{
    const auto& nodes = state.robots.at(robot_idx).estimation_nodes;
    Positions p(nodes.size());
    for (std::size_t k = 0; k < nodes.size(); ++k) p[k] = nodes[k].pose.block(0, 3, 3, 1);
    return p;
}

// ----- Errors ------------------------------------------------------------------

// Distance [m] between two shapes at each node.
std::vector<double> errorPerNode(const Positions& p, const Positions& p_truth)
{
    const std::size_t K = std::min(p.size(), p_truth.size());
    std::vector<double> e(K);
    for (std::size_t k = 0; k < K; ++k) e[k] = (p[k] - p_truth[k]).norm();
    return e;
}

double rmseOf(const std::vector<double>& e)
{
    if (e.empty()) return 0.0;
    double sse = 0.0;
    for (double x : e) sse += x * x;
    return std::sqrt(sse / static_cast<double>(e.size()));
}

double maxOf(const std::vector<double>& e)
{
    return e.empty() ? 0.0 : *std::max_element(e.begin(), e.end());
}

// ----- Formatting --------------------------------------------------------------

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
std::string fmtMeanStd(double mean, double std, int visual_width = 16, int prec = 2)
{
    std::ostringstream os;
    os << std::fixed << std::setprecision(prec)
       << std::setw(6) << mean << " ± "
       << std::setw(6) << std;
    const int pad = std::max(0, visual_width - 15);
    return std::string(pad, ' ') + os.str();
}

// ----- Command line ------------------------------------------------------------

struct Args
{
    std::string  config_path;
    std::string  csv_path = "";              // if non-empty, write the last trial's shapes here
    std::string  trials_csv_path = "";       // if non-empty, append one row per trial and method
    std::string  timing_csv_path = "";       // if non-empty, time every pipeline stage (see timePipelineStages)
    unsigned int timing_reps = 300;          // repetitions per stage for --timing-csv
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
                  << " [--timing-csv path] [--timing-reps N]"
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
        else if (s == "--timing-csv")   nextStr("--timing-csv",  a.timing_csv_path);
        else if (s == "--timing-reps") {
            if (i + 1 >= argc) { std::cerr << "--timing-reps requires a value\n"; return false; }
            char* end = nullptr;
            errno = 0;
            const long reps = std::strtol(argv[++i], &end, 10);
            if (end == argv[i] || *end != '\0' || errno == ERANGE || reps < 1
                || static_cast<unsigned long>(reps) > std::numeric_limits<unsigned int>::max()) {
                std::cerr << "--timing-reps must be an integer >= 1\n"; return false;
            }
            a.timing_reps = static_cast<unsigned int>(reps);
        }
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

// ----- One estimate, one trial -------------------------------------------------

using Clock = std::chrono::steady_clock;
double msSince(Clock::time_point t0)
{
    return std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
}

// Result of one estimator run.
struct MethodResult
{
    bool                converged = false;
    double              rmse = 0.0;          // [m]
    double              max_err = 0.0;       // [m]
    double              time_ms = 0.0;       // computeStateEstimate() only
    std::size_t         iterations = 0;
    std::vector<double> node_err;            // [m], one per node
    Positions           shape;               // estimated node positions
};

struct TrialResult
{
    unsigned int                          seed;
    Eigen::Matrix4d                       T_tip_noisy;
    std::array<MethodResult, N_METHODS>   method;
};

MethodResult runEstimator(const Estimator::RobotTopology& topology,
                          const Estimator::Hyperparameters& params,
                          const MethodSetup& setup,
                          const SensorMeas& tip_meas,
                          const Positions& p_truth,
                          unsigned int robot_idx)
{
    std::vector<SensorMeas> measurements = setup.model_measurements;
    measurements.push_back(tip_meas);

    Estimator estimator(topology, params, setup.options);
    SystemState state;
    std::vector<double> cost;
    const auto t0 = Clock::now();
    const bool converged = estimator.computeStateEstimate(
        state, cost, measurements, setup.control_inputs, false);

    MethodResult r;
    r.time_ms    = msSince(t0);
    r.converged  = converged;
    r.iterations = cost.size();
    r.shape      = nodePositionsOf(state, robot_idx);
    r.node_err   = errorPerNode(r.shape, p_truth);
    r.rmse       = rmseOf(r.node_err);
    r.max_err    = maxOf(r.node_err);
    return r;
}

// Noisy tip position measurement (position only, mask = [1,1,1,0,0,0]).
SensorMeas noisyTipMeasurement(const Eigen::Matrix4d& T_tip_true,
                               double sigma, unsigned int seed,
                               unsigned int robot_idx, unsigned int K)
{
    std::mt19937 rng(seed);
    std::normal_distribution<double> noise(0.0, sigma);
    Eigen::Matrix4d T_tip = T_tip_true;
    T_tip(0, 3) += noise(rng);
    T_tip(1, 3) += noise(rng);
    T_tip(2, 3) += noise(rng);

    SensorMeas m;
    m.type      = SensorMeas::Pose;
    m.value     = T_tip;
    m.mask      << 1, 1, 1, 0, 0, 0;
    m.idx_robot = robot_idx;
    m.idx_node  = K - 1;
    return m;
}

// ----- Statistics --------------------------------------------------------------

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

struct MethodStats
{
    Stat        rmse, max_err, time_ms;   // rmse/max_err in [m]
    std::size_t converged;
};

MethodStats statsOf(const std::vector<TrialResult>& trials, int method)
{
    std::vector<double> rmse, max_err, time_ms;
    std::size_t converged = 0;
    for (const auto& t : trials) {
        const MethodResult& r = t.method[method];
        rmse.push_back(r.rmse);
        max_err.push_back(r.max_err);
        time_ms.push_back(r.time_ms);
        if (r.converged) ++converged;
    }
    return { meanStd(rmse), meanStd(max_err), meanStd(time_ms), converged };
}

// ----- CSV output --------------------------------------------------------------

// Shapes of one trial, one row per node: the model, the truth, the noisy tip
// measurement (same on every row) and each method's estimate. Positions in m.
void writeShapeCsv(const std::string& path,
                   const Positions& p_model,
                   const Positions& p_truth,
                   const TrialResult& trial,
                   double L)
{
    const std::size_t K = p_truth.size();
    std::ofstream f(path);
    if (!f) {
        std::cerr << "Failed to open CSV path: " << path << "\n";
        return;
    }
    auto header = [&](const char* name) { f << "," << name << "_x," << name << "_y," << name << "_z"; };
    f << "node,s";
    header("model");
    header("truth");
    header("tip_meas");
    for (int m = 0; m < N_METHODS; ++m) header(kKey[m]);
    f << "\n";

    auto writeXyz = [&](const Eigen::Vector3d& p) { f << "," << p(0) << "," << p(1) << "," << p(2); };
    const Eigen::Vector3d tip_meas = trial.T_tip_noisy.block(0, 3, 3, 1);

    f << std::setprecision(8) << std::fixed;
    for (std::size_t k = 0; k < K; ++k) {
        const double s_k = (K > 1) ? (static_cast<double>(k) / (K - 1)) * L : 0.0;
        f << k << "," << s_k;
        writeXyz(p_model[k]);
        writeXyz(p_truth[k]);
        writeXyz(tip_meas);
        for (int m = 0; m < N_METHODS; ++m) writeXyz(trial.method[m].shape[k]);
        f << "\n";
    }
    std::cout << "Wrote CSV: " << path << " (" << K << " nodes)\n";
}

// One row per trial and method, plus one "model" row (trial = -1) for the
// model on its own. Appends; writes the header only when the file is new.
// Errors in metres; e0..e{K-1} are the per-node position errors.
void appendTrialsCsv(const std::string& path, const Args& args, double L,
                     const std::vector<double>& model_err,
                     const std::vector<TrialResult>& trials)
{
    const std::size_t K = model_err.size();
    const bool exists = std::ifstream(path).good();
    std::ofstream f(path, std::ios::app);
    if (!f) {
        std::cerr << "Failed to open trials-csv path: " << path << "\n";
        return;
    }
    if (!exists) {
        f << "fx,fy,fz,sigma,L,trial,seed,method,rmse,max_err,time_ms,iterations,converged";
        for (std::size_t k = 0; k < K; ++k) f << ",e" << k;
        f << "\n";
    }
    f << std::setprecision(9);
    auto writeRow = [&](long trial, long seed, const char* method, double rmse, double max_err,
                        double time_ms, std::size_t iters, bool converged,
                        const std::vector<double>& err) {
        f << args.fx << "," << args.fy << "," << args.fz << "," << args.sigma << ","
          << L << "," << trial << "," << seed << "," << method << ","
          << rmse << "," << max_err << "," << time_ms << "," << iters << ","
          << (converged ? 1 : 0);
        for (double e : err) f << "," << e;
        f << "\n";
    };
    writeRow(-1, -1, "model", rmseOf(model_err), maxOf(model_err), 0.0, 0, true, model_err);
    for (std::size_t t = 0; t < trials.size(); ++t)
        for (int m = 0; m < N_METHODS; ++m) {
            const MethodResult& r = trials[t].method[m];
            writeRow(t, trials[t].seed, kKey[m], r.rmse, r.max_err, r.time_ms,
                     r.iterations, r.converged, r.node_err);
        }
}

// Wall-clock time of every stage that runs BEFORE computeStateEstimate(), each
// repeated `reps` times; one row per repetition (stage,rep,ms,ok). The per-trial
// time_ms covers computeStateEstimate() only, so a model-based method's time per
// estimate is model_solve + dto + its adapter + estimator_construct + time_ms.
//   model_solve       Cosserat solve without the tip force, as the experiment does it:
//                     new object, default (straight) initial guess. Includes the tendon
//                     loads (computeInternalLoads). Used by every method except the Baseline.
//   model_solve_warm  information only: one object, initial guess = previous solution,
//                     tensions alternating between 1.00 and 1.01 x nominal (a 1 % change
//                     per frame, as for a robot that moves smoothly).
//   truth_solve       Cosserat solve WITH the tip force (simulated reality), warm-started
//                     from the model solution as in main(). It is not part of any method;
//                     reported only so it is not mistaken for one.
//   dto                 priorsFromCosseratModel (strains + loads -> ContinuumRodPriors).
//   adapter_<method>    cosseratPriorsToEstimator for that method.
//   estimator_construct_<method>  the ContinuumRobotStateEstimator constructor.
void timePipelineStages(const std::string& path, unsigned int reps,
                        const Eigen::Matrix<double, 6, 1>& q,
                        const Eigen::Vector3d& tip_force, const Eigen::Vector3d& tip_moment,
                        double L1, double L2,
                        const Estimator::RobotTopology& topology,
                        const Estimator::Hyperparameters& params,
                        const std::array<MethodSetup, N_METHODS>& setup,
                        unsigned int robot_idx)
{
    std::ofstream f(path);
    if (!f) {
        std::cerr << "Failed to open timing-csv path: " << path << "\n";
        return;
    }
    f << "stage,rep,ms,ok\n" << std::setprecision(9);
    auto row = [&](const std::string& stage, unsigned int r, double ms, bool ok) {
        f << stage << "," << r << "," << ms << "," << (ok ? 1 : 0) << "\n";
    };
    const Eigen::Vector3d zero = Eigen::Vector3d::Zero();

    CosseratRodModel model;                 // converged model solve, for dto/adapter
    Eigen::MatrixXd  diskFrames;
    if (!model.forwardKinematics(diskFrames, q, zero, zero) || !model.hasAuxOutputs()) {
        std::cerr << "Timing: Cosserat prior solve FAILED\n";
        return;
    }

    for (unsigned int r = 0; r < reps; ++r) {
        const auto t0 = Clock::now();
        CosseratRodModel m;
        Eigen::MatrixXd df;
        const bool ok = m.forwardKinematics(df, q, zero, zero) && m.hasAuxOutputs();
        row("model_solve", r, msSince(t0), ok);
    }
    {
        CosseratRodModel m;
        Eigen::MatrixXd df;
        m.forwardKinematics(df, q, zero, zero);
        for (unsigned int r = 0; r < reps; ++r) {
            const Eigen::Matrix<double, 6, 1> q_next = q * ((r % 2) ? 1.00 : 1.01);
            m.setDefaultInitValues(m.getFinalInitValues());
            const auto t0 = Clock::now();
            const bool ok = m.forwardKinematics(df, q_next, zero, zero) && m.hasAuxOutputs();
            row("model_solve_warm", r, msSince(t0), ok);
        }
    }
    for (unsigned int r = 0; r < reps; ++r) {
        CosseratRodModel m;
        m.setDefaultInitValues(model.getFinalInitValues());   // as in main()
        Eigen::MatrixXd df;
        const auto t0 = Clock::now();
        const bool ok = m.forwardKinematics(df, q, tip_force, tip_moment);
        row("truth_solve", r, msSince(t0), ok);
    }

    const ControlInputMode mode[N_METHODS] = { ControlInputMode::None, ControlInputMode::None,
                                               ControlInputMode::ForceAsInput,
                                               ControlInputMode::StrainAsInput };
    for (unsigned int r = 0; r < reps; ++r) {
        auto t0 = Clock::now();
        const ContinuumRodPriors priors = priorsFromCosseratModel(model, L1, L2, diskFrames);
        row("dto", r, msSince(t0), true);
        for (int m = StrainMeas; m < N_METHODS; ++m) {
            t0 = Clock::now();
            const EstimatorPriors ep = cosseratPriorsToEstimator(priors, topology, robot_idx,
                                                                 diskFrames, mode[m]);
            row(std::string("adapter_") + kKey[m], r, msSince(t0),
                !ep.measurements.empty() || !ep.control_inputs.empty());
        }
        for (int m = 0; m < N_METHODS; ++m) {
            t0 = Clock::now();
            Estimator estimator(topology, params, setup[m].options);
            row(std::string("estimator_construct_") + kKey[m], r, msSince(t0), true);
        }
    }
    std::cout << "Wrote timing CSV: " << path << " (" << reps << " repetitions per stage)\n";
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
    const Eigen::Vector3d tip_force(args.fx, args.fy, args.fz);
    const Eigen::Vector3d tip_moment = Eigen::Vector3d::Zero();
    const double L1 = 0.1, L2 = 0.1;

    // --- 2. Cosserat solves: model (no force) and ground truth (with force) --
    CosseratRodModel rod_model;
    Eigen::MatrixXd diskFrames_model;
    if (!rod_model.forwardKinematics(diskFrames_model, q,
                                     Eigen::Vector3d::Zero(),
                                     Eigen::Vector3d::Zero())
        || !rod_model.hasAuxOutputs()) {
        std::cerr << "Cosserat prior solve FAILED (residual = "
                  << rod_model.getFinalResudial() << ")\n";
        return 1;
    }
    const ContinuumRodPriors priors = priorsFromCosseratModel(rod_model, L1, L2, diskFrames_model);
    const Eigen::VectorXd s_disks_model = rod_model.getArclengthSamples();

    CosseratRodModel rod_truth;
    // Warm-start the real solve from the model's converged shooting state.
    // With the original Cosserat solver (before 2026-10-14) the default initial
    // guess (straight rod) was unreliable for tip forces above ~0.02 N when the
    // rod is already heavily bent by tendon tensions; starting from the model's
    // converged (u0, v0) is also faster.
    rod_truth.setDefaultInitValues(rod_model.getFinalInitValues());
    Eigen::MatrixXd diskFrames_truth;
    if (!rod_truth.forwardKinematics(diskFrames_truth, q, tip_force, tip_moment)
        || !rod_truth.hasAuxOutputs()) {
        std::cerr << "Cosserat real solve FAILED (residual = "
                  << rod_truth.getFinalResudial() << ")\n";
        return 1;
    }
    const Eigen::VectorXd s_disks_truth = rod_truth.getArclengthSamples();

    const Positions p_model = nodePositionsFromCosserat(diskFrames_model, s_disks_model, Ti0, K, L);
    const Positions p_truth = nodePositionsFromCosserat(diskFrames_truth, s_disks_truth, Ti0, K, L);
    const std::vector<double> model_err = errorPerNode(p_model, p_truth);
    const double model_rmse = rmseOf(model_err);
    const Eigen::Matrix4d T_tip_true =
        cosseratFrameToWorld(nearestDiskFrame(L, s_disks_truth, diskFrames_truth), Ti0);

    // --- 3. What each method gets from the model (once, deterministic) ------
    // Strain as measurement: model strain as measurements, no control inputs.
    // Force / Strain as input: control inputs only; the starting shape is the
    // model's (with zero strain bias for Strain as input).
    const auto t_adapter0 = Clock::now();
    const EstimatorPriors from_model[N_METHODS] = {
        {},
        cosseratPriorsToEstimator(priors, topology, robot_idx, diskFrames_model, ControlInputMode::None),
        cosseratPriorsToEstimator(priors, topology, robot_idx, diskFrames_model, ControlInputMode::ForceAsInput),
        cosseratPriorsToEstimator(priors, topology, robot_idx, diskFrames_model, ControlInputMode::StrainAsInput),
    };
    const double adapter_ms = msSince(t_adapter0);

    auto solver_options = base_options;
    solver_options.solver = Estimator::Options::NewtonLineSearch;

    std::array<MethodSetup, N_METHODS> setup;
    for (int m = 0; m < N_METHODS; ++m) {
        setup[m].options = solver_options;
        if (m == Baseline) {
            setup[m].options.init_guess_type = Estimator::Options::Straight;
            continue;
        }
        setup[m].options.init_guess_type = Estimator::Options::Custom;
        setup[m].options.custom_guess    = from_model[m].initial_guess;
        if (m == StrainMeas) setup[m].model_measurements = from_model[m].measurements;
        else                 setup[m].control_inputs     = from_model[m].control_inputs;
    }

    if (!args.timing_csv_path.empty())
        timePipelineStages(args.timing_csv_path, args.timing_reps, q, tip_force, tip_moment,
                           L1, L2, topology, params, setup, robot_idx);

    // --- 4. Trials: same noisy tip measurement for every method -------------
    std::vector<TrialResult> trials;
    trials.reserve(args.n_trials);
    for (unsigned int t = 0; t < args.n_trials; ++t) {
        TrialResult trial;
        trial.seed = args.seed + t;
        const SensorMeas tip_meas = noisyTipMeasurement(T_tip_true, args.sigma, trial.seed, robot_idx, K);
        trial.T_tip_noisy = tip_meas.value;
        for (int m = 0; m < N_METHODS; ++m)
            trial.method[m] = runEstimator(topology, params, setup[m], tip_meas, p_truth, robot_idx);
        trials.push_back(std::move(trial));
    }

    // --- 5. Aggregate stats -------------------------------------------------
    std::array<MethodStats, N_METHODS> stats;
    for (int m = 0; m < N_METHODS; ++m) stats[m] = statsOf(trials, m);

    // --- 6. Reporting -------------------------------------------------------
    //
    // Layout (each block separated by a blank line):
    //   1. Cosserat ground truth      — tip positions + model-truth gap.
    //   2. Headline                   — one-paragraph plain-English verdict.
    //   3. Per-method results table   — mm-based, fixed-point.
    //   4. Comparisons                — relative numbers (vs. model, vs. Baseline).
    //
    // The headline is computed AFTER the trials run; it picks winners on each
    // axis (accuracy / speed) among methods that produced a usable shape
    // (RMSE < L/2).

    const double half_L_mm = L * 1000.0 * 0.5;
    auto rmseMm = [&](int m) { return stats[m].rmse.mean * 1000.0; };
    auto timeMs = [&](int m) { return stats[m].time_ms.mean; };
    auto isFailed = [&](int m) { return rmseMm(m) >= half_L_mm || stats[m].converged == 0; };

    int accuracy_winner = -1, speed_winner = -1;
    for (int m = 0; m < N_METHODS; ++m) {
        if (rmseMm(m) >= half_L_mm) continue;
        if (accuracy_winner < 0 || rmseMm(m) < rmseMm(accuracy_winner)) accuracy_winner = m;
        if (speed_winner    < 0 || timeMs(m) < timeMs(speed_winner))    speed_winner    = m;
    }

    // --- Block 1: Cosserat ground truth ------------------------------------
    std::cout << "Cosserat ground truth\n"
              << "  Prior tip position  : " << fmtVec3Mm(p_model.back()) << "\n"
              << "  Truth tip position  : " << fmtVec3Mm(p_truth.back()) << "\n"
              << "  Prior<->truth gap   : "
              << std::fixed << std::setprecision(2) << model_rmse * 1000.0
              << " mm RMSE   (baseline without estimator)\n\n";

    // --- Block 2: Headline -------------------------------------------------
    std::cout << "[*] Headline\n";
    if (accuracy_winner < 0) {
        std::cout << "  No method produced a usable shape (all RMSE > L/2 = "
                  << std::fixed << std::setprecision(0) << half_L_mm << " mm).\n";
    } else if (accuracy_winner == speed_winner) {
        const int w = accuracy_winner;
        std::cout << "  " << kMethodName[w]
                  << " wins on BOTH accuracy and speed: "
                  << std::fixed << std::setprecision(2)
                  << rmseMm(w) << " mm RMSE in " << timeMs(w) << " ms.\n";
        // Show the speedup against the slowest other method.
        int slow = (w == 0) ? 1 : 0;
        for (int m = 0; m < N_METHODS; ++m)
            if (m != w && timeMs(m) > timeMs(slow)) slow = m;
        if (rmseMm(slow) < half_L_mm) {
            std::cout << "  That is "
                      << std::fixed << std::setprecision(1) << rmseMm(slow) / rmseMm(w)
                      << "x more accurate and "
                      << std::setprecision(0) << timeMs(slow) / timeMs(w)
                      << "x faster than " << kMethodName[slow] << ".\n";
        }
    } else {
        std::cout << "  Mixed result (we are past the crossover where the prior stops helping on accuracy):\n"
                  << "    " << kMethodName[accuracy_winner] << " wins on accuracy: "
                  << std::fixed << std::setprecision(2) << rmseMm(accuracy_winner) << " mm RMSE.\n"
                  << "    " << kMethodName[speed_winner]    << " wins on speed:    "
                  << std::fixed << std::setprecision(2) << timeMs(speed_winner) << " ms ("
                  << std::setprecision(0) << timeMs(accuracy_winner) / timeMs(speed_winner)
                  << "x faster than " << kMethodName[accuracy_winner] << ").\n";
    }
    for (int m = 0; m < N_METHODS; ++m) {
        if (!isFailed(m)) continue;
        std::cout << "  " << kMethodName[m] << " FAILS: "
                  << std::fixed << std::setprecision(1) << rmseMm(m) << " mm RMSE";
        if (stats[m].converged == 0)
            std::cout << ", " << stats[m].converged << "/" << args.n_trials << " converged";
        std::cout << ".\n";
    }
    std::cout << "\n";

    // --- Block 3: per-method results table ---------------------------------
    //
    // Column layout (every width below is in *visual* columns, not bytes):
    //
    //   "  [mark] [name padded W_NAME]  [RMSE W_VAL]  [MaxErr W_VAL]  [Time W_VAL]  [Conv W_CONV]"
    //
    // fmtMeanStd does its own visual-width padding (it handles the ± byte/visual
    // mismatch internally), so we do NOT wrap it in std::setw at the call site —
    // only header strings (pure ASCII) use std::setw.
    constexpr int W_MARK = 3;
    constexpr int W_NAME = 24;
    constexpr int W_VAL  = 16;     // applies to RMSE, MaxErr, Time columns
    constexpr int W_CONV = 8;
    const std::string SEP = "  ";  // 2-space column separator

    const int line_width = 2 + W_MARK + 1 + W_NAME
                         + 3 * (static_cast<int>(SEP.size()) + W_VAL)
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

    for (int m = 0; m < N_METHODS; ++m) {
        // Markers: [*] = winner on at least one axis; [X] = failed/invalid.
        const char* marker = "   ";
        if (isFailed(m))                                     marker = "[X]";
        else if (m == accuracy_winner || m == speed_winner)  marker = "[*]";
        const MethodStats& s = stats[m];
        const std::string conv = std::to_string(s.converged) + "/" + std::to_string(args.n_trials);
        std::cout << "  " << std::left  << std::setw(W_MARK) << marker << " "
                  <<        std::setw(W_NAME) << kTableLabel[m]
                  << SEP << fmtMeanStd(s.rmse.mean * 1000.0,    s.rmse.std * 1000.0,    W_VAL)
                  << SEP << fmtMeanStd(s.max_err.mean * 1000.0, s.max_err.std * 1000.0, W_VAL)
                  << SEP << fmtMeanStd(s.time_ms.mean,          s.time_ms.std,          W_VAL)
                  << SEP << std::right << std::setw(W_CONV) << conv
                  << "\n";
    }
    std::cout << hline << "\n"
              << "  [*] = best on at least one axis     [X] = failed (RMSE > L/2 or 0 converged)\n"
              << "  Adapter prep (one-shot, not included in per-trial timings): "
              << std::fixed << std::setprecision(2) << adapter_ms << " ms\n\n";

    // --- Block 4: comparisons ----------------------------------------------
    std::cout << "Comparisons\n";
    const double gap_mm = model_rmse * 1000.0;
    std::cout << "  vs. prior" << "\xE2\x86\x94" << "truth gap ("
              << std::fixed << std::setprecision(2)
              << gap_mm << " mm - the no-estimator baseline):\n";
    for (int m = 0; m < N_METHODS; ++m) {
        std::ostringstream tag;
        if (gap_mm > 0 && rmseMm(m) < gap_mm) {
            tag << std::fixed << std::setprecision(0) << (gap_mm - rmseMm(m)) / gap_mm * 100.0
                << "% better than prior";
        } else if (gap_mm > 0) {
            tag << std::fixed << std::setprecision(1) << rmseMm(m) / gap_mm << "x WORSE than prior";
        }
        std::cout << "    " << std::left << std::setw(28) << kMethodName[m] << std::right
                  << " : " << std::fixed << std::setprecision(2)
                  << std::setw(7) << rmseMm(m) << " mm   (" << tag.str() << ")\n";
    }
    std::cout << "\n"
              << "  Strain as measurement vs Baseline:\n"
              << "    Accuracy: " << std::fixed << std::setprecision(2)
              << rmseMm(StrainMeas) << " vs " << rmseMm(Baseline) << " mm";
    if (rmseMm(Baseline) > 0) {
        const double pct = (rmseMm(StrainMeas) - rmseMm(Baseline)) / rmseMm(Baseline) * 100.0;
        std::cout << "  ("
                  << std::fixed << std::setprecision(0) << std::abs(pct) << "% "
                  << (pct < 0 ? "better" : "worse") << ")";
    }
    std::cout << "\n    Speed:    "
              << std::fixed << std::setprecision(2) << timeMs(StrainMeas) << " vs "
              << timeMs(Baseline) << " ms   ("
              << std::fixed << std::setprecision(1)
              << timeMs(Baseline) / timeMs(StrainMeas) << "x faster)\n\n";

    // --- 7. Optional CSV output ---------------------------------------------
    if (!args.trials_csv_path.empty())
        appendTrialsCsv(args.trials_csv_path, args, L, model_err, trials);
    if (!args.csv_path.empty())
        writeShapeCsv(args.csv_path, p_model, p_truth, trials.back(), L);

    printBanner("Evaluation OK");

    // Exit non-zero if any method failed to converge on the LAST trial — useful
    // for CI checks but harmless for normal use.
    for (const auto& r : trials.back().method)
        if (!r.converged) return 2;
    return 0;
}
