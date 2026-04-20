// test_cosserat_adapter.cpp
// Validation tests for the Cosserat -> estimator integration (Phase 5).
//
// Requires USE_LOCAL_TDCR=ON (links tdcr_modeling for tests C / D).
//
// Tests:
//   A. Frame permutation correctness (Cosserat body -> estimator body)
//   B. Arclength resampling accuracy on a linear strain field
//   C. End-to-end convergence with realistic Cosserat priors
//   D. Prior-informed estimate vs no-prior baseline
//
// Run from the project root (the adapter binary is placed in examples/):
//   ./examples/test_cosserat_adapter [path_to_config.yaml]
// Default config path: config/6_cosserat_priors.yaml (relative to the repo root).

#include "config_loader.h"
#include "continuum_robot_state_estimator.h"
#include "continuum_rod_priors.h"
#include "cosserat_priors_adapter.h"
#include "cosseratrodmodel.h"

#include <Eigen/Core>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

int g_pass = 0;
int g_fail = 0;

#define CHECK(cond, msg)                                                     \
    do {                                                                     \
        if (!(cond)) {                                                       \
            std::cerr << "  FAIL: " << msg << " (line " << __LINE__ << ")\n";\
            return false;                                                    \
        }                                                                    \
    } while (0)

void runTest(const std::string& name, const std::function<bool()>& fn)
{
    std::cout << "[RUN]  " << name << "\n";
    const bool ok = fn();
    std::cout << (ok ? "[OK]   " : "[FAIL] ") << name << "\n\n";
    if (ok) ++g_pass; else ++g_fail;
}

bool vecNearEq(const Eigen::MatrixXd& a, const Eigen::MatrixXd& b, double tol)
{
    if (a.rows() != b.rows() || a.cols() != b.cols()) return false;
    return (a - b).cwiseAbs().maxCoeff() <= tol;
}

// Build a minimal RobotTopology for tests.
ContinuumRobotStateEstimator::RobotTopology makeTopology(unsigned int K, double L,
                                                         bool lock_first_pose = false,
                                                         bool lock_first_strain = false)
{
    ContinuumRobotStateEstimator::RobotTopology t;
    t.N                   = 1;
    t.K                   = {K};
    t.M                   = {1};
    t.L                   = {L};
    t.Ti0                 = {Eigen::Matrix4d::Identity()};
    t.common_end_effector = false;
    t.lock_first_pose     = {lock_first_pose};
    t.lock_last_pose      = {false};
    t.lock_first_strain   = {lock_first_strain};
    t.lock_last_strain    = {false};
    t.fbg_theta_offset    = {0};
    t.fbg_core_distance   = {0};
    return t;
}

// Build a (4, 4*N) matrix of identity disk frames.
Eigen::MatrixXd makeIdentityDiskFrames(int N)
{
    Eigen::MatrixXd frames = Eigen::MatrixXd::Zero(4, 4 * N);
    for (int i = 0; i < N; ++i)
        frames.block(0, 4 * i, 4, 4) = Eigen::Matrix4d::Identity();
    return frames;
}

// Build a ContinuumRodPriors with constant v, u across two samples at s=[0, L].
ContinuumRodPriors makeConstantPriors(const Eigen::Vector3d& v,
                                      const Eigen::Vector3d& u,
                                      double L = 0.1)
{
    ContinuumRodPriors p;
    p.s = Eigen::Vector2d(0.0, L);
    p.v.resize(2, 3); p.v.row(0) = v.transpose(); p.v.row(1) = v.transpose();
    p.u.resize(2, 3); p.u.row(0) = u.transpose(); p.u.row(1) = u.transpose();
    return p;
}

// ============================================================================
// Test A — Frame permutation correctness
// ============================================================================
bool testA_framePermutation()
{
    const auto topo   = makeTopology(2, 0.1);
    const auto frames = makeIdentityDiskFrames(2);
    const double tol  = 1e-12;

    struct Case {
        const char*                 label;
        Eigen::Vector3d             v_cos;
        Eigen::Vector3d             u_cos;
        Eigen::Matrix<double, 6, 1> expected;
    };

    const std::vector<Case> cases = {
        { "straight (v=[0,0,1], u=0)",
          {0,0,1}, {0,0,0}, (Eigen::Matrix<double,6,1>() << 1,0,0, 0,0,0).finished() },
        { "Cosserat bend about x (u=[1,0,0]) -> estimator om_2=1",
          {0,0,1}, {1,0,0}, (Eigen::Matrix<double,6,1>() << 1,0,0, 0,1,0).finished() },
        { "Cosserat bend about y (u=[0,1,0]) -> estimator om_3=1",
          {0,0,1}, {0,1,0}, (Eigen::Matrix<double,6,1>() << 1,0,0, 0,0,1).finished() },
        { "Cosserat torsion about z (u=[0,0,1]) -> estimator om_1=1",
          {0,0,1}, {0,0,1}, (Eigen::Matrix<double,6,1>() << 1,0,0, 1,0,0).finished() },
        { "Cosserat shear x (v=[1,0,0]) -> estimator nu_2=1",
          {1,0,0}, {0,0,0}, (Eigen::Matrix<double,6,1>() << 0,1,0, 0,0,0).finished() },
        { "Cosserat shear y (v=[0,1,0]) -> estimator nu_3=1",
          {0,1,0}, {0,0,0}, (Eigen::Matrix<double,6,1>() << 0,0,1, 0,0,0).finished() },
    };

    for (const auto& c : cases) {
        const auto priors = makeConstantPriors(c.v_cos, c.u_cos);
        const auto ep     = cosseratPriorsToEstimator(priors, topo, 0, frames);
        CHECK(!ep.measurements.empty(), c.label);
        const Eigen::Matrix<double, 6, 1> got = ep.measurements[0].value;
        if (!vecNearEq(got, c.expected, tol)) {
            std::cerr << "  FAIL: " << c.label
                      << "\n         expected: [" << c.expected.transpose() << "]"
                      << "\n         got:      [" << got.transpose()      << "]\n";
            return false;
        }
    }
    return true;
}

// ============================================================================
// Test B — Arclength resampling on a linear strain field
// ============================================================================
bool testB_resamplingAccuracy()
{
    // Cosserat samples on a 5-point grid with linear v(s) = [0, 0, 1 + s]
    const int N_cos = 5;
    ContinuumRodPriors priors;
    priors.s.resize(N_cos);
    priors.v.resize(N_cos, 3);
    priors.u = Eigen::MatrixXd::Zero(N_cos, 3);
    for (int i = 0; i < N_cos; ++i) {
        const double s = 0.05 * i;        // 0, 0.05, 0.10, 0.15, 0.20
        priors.s(i)   = s;
        priors.v.row(i) << 0.0, 0.0, 1.0 + s;
    }

    // Resample to K=9 nodes (finer uniform grid, interior points fall between Cosserat samples)
    const unsigned int K = 9;
    const double       L = 0.20;
    const auto topo      = makeTopology(K, L);
    const auto frames    = makeIdentityDiskFrames(N_cos);

    const auto ep = cosseratPriorsToEstimator(priors, topo, 0, frames);

    CHECK(ep.measurements.size() == K, "measurements count");

    // After frame permutation: estimator nu_1 = v_cos(2) = 1 + s_k
    const double tol = 1e-12;
    for (unsigned int k = 0; k < K; ++k) {
        const double s_k         = (static_cast<double>(k) / (K - 1)) * L;
        const double nu1_expected = 1.0 + s_k;
        const double nu1_got      = ep.measurements[k].value(0);
        if (std::abs(nu1_got - nu1_expected) > tol) {
            std::cerr << "  FAIL: k=" << k << "  s=" << s_k
                      << "  expected nu_1 = " << nu1_expected
                      << "  got " << nu1_got << "\n";
            return false;
        }
    }
    return true;
}

// ============================================================================
// Helpers for tests C / D
// ============================================================================
struct PipelineResult {
    ContinuumRobotStateEstimator::SystemState state;
    std::vector<double>                       cost;
    bool                                      converged;
};

// Runs the estimator with the given setup. Returns the result.
PipelineResult runEstimator(
    const ContinuumRobotStateEstimator::RobotTopology&                    topology,
    const ContinuumRobotStateEstimator::Hyperparameters&                  params,
    ContinuumRobotStateEstimator::Options                                 options,
    const std::vector<ContinuumRobotStateEstimator::SensorMeasurement>&   measurements,
    const std::vector<ContinuumRobotStateEstimator::ControlInput>&        control_inputs)
{
    ContinuumRobotStateEstimator estimator(topology, params, options);
    PipelineResult r;
    r.converged = estimator.computeStateEstimate(
        r.state, r.cost, measurements, control_inputs, /*verbose=*/false);
    return r;
}

// Mean absolute strain difference between estimator and Cosserat priors, per-node.
double meanAbsStrainDiff(
    const ContinuumRobotStateEstimator::SystemState&                     state,
    const std::vector<ContinuumRobotStateEstimator::SensorMeasurement>&  cosserat_priors,
    unsigned int                                                          robot_idx)
{
    const auto& nodes = state.robots.at(robot_idx).estimation_nodes;
    const size_t K    = std::min(nodes.size(), cosserat_priors.size());
    if (K == 0) return 0.0;

    double sum = 0.0;
    for (size_t k = 0; k < K; ++k) {
        const Eigen::Matrix<double, 6, 1>& s_cos = cosserat_priors[k].value;
        const Eigen::Matrix<double, 6, 1>& s_est = nodes[k].strain;
        sum += (s_cos - s_est).cwiseAbs().sum() / 6.0;
    }
    return sum / static_cast<double>(K);
}

// Run Cosserat FK with the default-driver tensions and return the populated DTO (Data Transfer Object).
bool runCosseratFk(CosseratRodModel& model,
                   Eigen::MatrixXd& diskFrames,
                   ContinuumRodPriors& priors)
{
    Eigen::Matrix<double, 6, 1> q;
    q << 0.5, 0.2, 0.0, 0.3, 0.0, 0.1;
    const Eigen::Vector3d f_ext = Eigen::Vector3d::Zero();
    const Eigen::Vector3d l_ext = Eigen::Vector3d::Zero();

    if (!model.forwardKinematics(diskFrames, q, f_ext, l_ext)) return false;
    if (!model.hasAuxOutputs())                                return false;
    priors = priorsFromCosseratModel(model, /*L1=*/0.1, /*L2=*/0.1);
    return true;
}

// ============================================================================
// Test C — End-to-end convergence
// Uses the YAML config so hyperparameters match the driver exactly.
// ============================================================================
bool testC_endToEndConvergence(const std::string& config_path)
{
    ConfigLoader    config(config_path);
    const auto      topology    = config.getTopology();
    const auto      params      = config.getHyperparameters();
    auto            options     = config.getOptions();
    const unsigned  robot_idx   = 0;

    CosseratRodModel   model;
    Eigen::MatrixXd    diskFrames;
    ContinuumRodPriors priors;
    CHECK(runCosseratFk(model, diskFrames, priors), "Cosserat FK converges");

    const auto ep = cosseratPriorsToEstimator(priors, topology, robot_idx, diskFrames);
    options.custom_guess = ep.initial_guess;

    const auto r = runEstimator(topology, params, options, ep.measurements, ep.control_inputs);
    CHECK(r.converged, "estimator converges");
    CHECK(r.cost.size() >= 1, "cost trajectory populated");

    // Per-node max diff between estimator strain and the Cosserat prior we injected.
    double max_diff = 0.0;
    const auto& nodes = r.state.robots.at(robot_idx).estimation_nodes;
    CHECK(nodes.size() == ep.measurements.size(),
          "estimation_nodes count matches measurement count");
    for (size_t k = 0; k < nodes.size(); ++k) {
        const double d = (ep.measurements[k].value - nodes[k].strain).cwiseAbs().maxCoeff();
        max_diff       = std::max(max_diff, d);
    }
    std::cout << "  max |diff|_inf = " << max_diff << "\n";
    CHECK(max_diff < 1.0e-2, "max strain diff < 1% (expected ~5e-3 from junction smoothing)");
    return true;
}

// ============================================================================
// Test D — Prior vs no-prior comparison
// Run 1: Straight guess, no measurements, no control inputs.
// Run 2: Cosserat priors (measurements + control inputs + custom initial guess).
// Run 2 must be closer to the Cosserat prediction than Run 1.
// ============================================================================
bool testD_priorVsNoPrior(const std::string& config_path)
{
    ConfigLoader   config(config_path);
    auto           topology  = config.getTopology();
    const auto     params    = config.getHyperparameters();
    const unsigned robot_idx = 0;

    // Lock both the base pose AND the base strain so the no-prior run
    // (empty measurements) has a well-posed system. Both runs use the same
    // topology, so the only variable under test is whether Cosserat priors
    // are injected as measurements + control inputs + custom initial guess.
    topology.lock_first_strain.at(robot_idx) = true;

    CosseratRodModel   model;
    Eigen::MatrixXd    diskFrames;
    ContinuumRodPriors priors;
    CHECK(runCosseratFk(model, diskFrames, priors), "Cosserat FK converges");
    const auto ep = cosseratPriorsToEstimator(priors, topology, robot_idx, diskFrames);

    // Run 1 — no prior. Straight initial guess, empty measurements + control inputs.
    ContinuumRobotStateEstimator::Options opts_none = config.getOptions();
    opts_none.init_guess_type = ContinuumRobotStateEstimator::Options::Straight;
    const auto r_none = runEstimator(topology, params, opts_none, {}, {});
    CHECK(r_none.converged, "no-prior run converges (to the Straight initial guess)");

    // Run 2 — Cosserat priors.
    ContinuumRobotStateEstimator::Options opts_prior = config.getOptions();
    opts_prior.init_guess_type = ContinuumRobotStateEstimator::Options::Custom;
    opts_prior.custom_guess    = ep.initial_guess;
    const auto r_prior = runEstimator(topology, params, opts_prior, ep.measurements, ep.control_inputs);
    CHECK(r_prior.converged, "prior-informed run converges");

    const double diff_none  = meanAbsStrainDiff(r_none.state,  ep.measurements, robot_idx);
    const double diff_prior = meanAbsStrainDiff(r_prior.state, ep.measurements, robot_idx);

    std::cout << "  mean |diff|  no-prior = " << diff_none  << "\n";
    std::cout << "  mean |diff|  prior    = " << diff_prior << "\n";
    CHECK(diff_prior < diff_none, "prior-informed estimate must be closer to Cosserat");
    // Sanity: the margin should be substantial (Cosserat bending dominates Run 1's diff).
    CHECK(diff_prior * 10.0 < diff_none,
          "prior-informed diff should be at least 10x smaller than no-prior diff");
    return true;
}

}  // anonymous namespace


int main(int argc, char* argv[])
{
    const std::string config_path =
        (argc >= 2) ? argv[1] : "config/6_cosserat_priors.yaml";

    std::cout << "=== Cosserat adapter + integration tests ===\n";
    std::cout << "Config: " << config_path << "\n\n";

    runTest("A. Frame permutation",              testA_framePermutation);
    runTest("B. Resampling (linear strain)",     testB_resamplingAccuracy);
    runTest("C. End-to-end convergence",
            [&config_path] { return testC_endToEndConvergence(config_path); });
    runTest("D. Prior vs no-prior",
            [&config_path] { return testD_priorVsNoPrior(config_path); });

    std::cout << "=== Results: " << g_pass << " passed, " << g_fail << " failed ===\n";
    return g_fail == 0 ? 0 : 1;
}
