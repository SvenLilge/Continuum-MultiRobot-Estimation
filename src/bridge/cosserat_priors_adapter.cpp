#include "cosserat_priors_adapter.h"
#include "cosseratrodmodel.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

ContinuumRodPriors priorsFromCosseratModel(const CosseratRodModel& model,
                                           double L1, double L2)
{
    ContinuumRodPriors p;

    p.s          = model.getArclengthSamples();
    p.v          = model.getStrainV();
    p.u          = model.getStrainU();
    p.v_dot      = model.getStrainVDot();
    p.u_dot      = model.getStrainUDot();
    p.n_internal = model.getInternalForce();
    p.m_internal = model.getInternalMoment();
    p.f_dist     = model.getDistributedForce();
    p.l_dist     = model.getDistributedMoment();

    Eigen::Vector3d F_junction, L_junction, F_tip, L_tip;
    model.getDiscreteLoads(F_junction, L_junction, F_tip, L_tip);
    p.s_discrete = { L1, L1 + L2 };
    p.F_discrete = { F_junction, F_tip };
    p.L_discrete = { L_junction, L_tip };

    p.validate();
    return p;
}


// ============================================================================
// cosseratPriorsToEstimator — frame-permute + resample a ContinuumRodPriors
// into estimator-ready measurements, control inputs, and initial guess.
// Math / conventions: see doc/bridge_adapter_guide.md §4.2, §4.3, §4.4.
// ============================================================================

namespace {

// Cosserat body (z-forward) -> estimator body (x-forward).
// Cyclic permutation; the matrix form is equivalent to the elementwise
// reorder [v0, v1, v2] -> [v2, v0, v1] for any 3-vector.
const Eigen::Matrix3d& rotConv()
{
    static const Eigen::Matrix3d R = (Eigen::Matrix3d() <<
        0, 0, 1,
        1, 0, 0,
        0, 1, 0).finished();
    return R;
}

// Permute a Cosserat body-frame strain (v, u) into the estimator's 6x1 strain
// [nu_1, nu_2, nu_3, omega_1, omega_2, omega_3].
Eigen::Matrix<double,6,1> permuteStrain(
    const Eigen::Vector3d& v_cos, const Eigen::Vector3d& u_cos)
{
    Eigen::Matrix<double,6,1> s;
    s << v_cos(2), v_cos(0), v_cos(1),
         u_cos(2), u_cos(0), u_cos(1);
    return s;
}

// Linearly interpolate an (N x 3) matrix at arclength s_k along the s vector.
// Clamps at endpoints. s must be sorted non-decreasingly.
Eigen::Vector3d interpolateRowAt(
    double s_k, const Eigen::VectorXd& s, const Eigen::MatrixXd& M)
{
    const Eigen::Index N = s.size();
    if (N == 0 || M.rows() != N) return Eigen::Vector3d::Zero();
    if (s_k <= s(0))       return M.row(0).transpose();
    if (s_k >= s(N - 1))   return M.row(N - 1).transpose();

    Eigen::Index j = 0;
    while (j + 1 < N && s(j + 1) < s_k) ++j;
    const double denom = s(j + 1) - s(j);
    if (denom <= 0.0) return M.row(j).transpose();
    const double alpha = (s_k - s(j)) / denom;
    return ((1.0 - alpha) * M.row(j) + alpha * M.row(j + 1)).transpose();
}

// Nearest-neighbor lookup: the 4x4 disk frame whose arclength sample is
// closest to s_k. Used for the initial-guess pose (simpler than SLERP; the
// estimator refines it anyway).
Eigen::Matrix4d nearestDiskFrame(
    double s_k, const Eigen::VectorXd& s, const Eigen::MatrixXd& diskFrames)
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

} // anonymous namespace


EstimatorPriors cosseratPriorsToEstimator(
    const ContinuumRodPriors&                            priors,
    const ContinuumRobotStateEstimator::RobotTopology&   topology,
    unsigned int                                         robot_idx,
    const Eigen::MatrixXd&                               diskFrames)
{
    priors.validate();

    if (robot_idx >= topology.N)
        throw std::invalid_argument(
            "cosseratPriorsToEstimator: robot_idx >= topology.N");

    if (priors.s.size() == 0)
        throw std::invalid_argument(
            "cosseratPriorsToEstimator: priors.s is empty");

    const Eigen::Index N_disk = priors.s.size();
    if (diskFrames.rows() != 4 || diskFrames.cols() != 4 * N_disk)
        throw std::invalid_argument(
            "cosseratPriorsToEstimator: diskFrames shape must be 4 x (4 * priors.s.size())");

    const unsigned int K = topology.K.at(robot_idx);
    const double       L = topology.L.at(robot_idx);
    if (K < 2)
        throw std::invalid_argument(
            "cosseratPriorsToEstimator: topology.K[robot_idx] must be >= 2");

    const Eigen::Matrix3d& R_conv = rotConv();
    const Eigen::Matrix3d  R_convT = R_conv.transpose();

    EstimatorPriors ep;
    ep.initial_guess.robots.resize(topology.N);
    auto& target = ep.initial_guess.robots[robot_idx];
    target.estimation_nodes.resize(K);

    // Per-node resampled values, reused across measurements + control inputs
    // + initial-guess construction.
    std::vector<Eigen::Vector3d>           v_at(K), u_at(K), vdot_at(K), udot_at(K);
    std::vector<Eigen::Matrix<double,6,1>> strain_at(K);
    std::vector<double>                    s_at(K);

    for (unsigned int k = 0; k < K; ++k) {
        const double s_k = (static_cast<double>(k) / (K - 1)) * L;
        s_at[k]      = s_k;
        v_at[k]      = interpolateRowAt(s_k, priors.s, priors.v);
        u_at[k]      = interpolateRowAt(s_k, priors.s, priors.u);
        vdot_at[k]   = (priors.v_dot.size() == 0)
                           ? Eigen::Vector3d::Zero()
                           : interpolateRowAt(s_k, priors.s, priors.v_dot);
        udot_at[k]   = (priors.u_dot.size() == 0)
                           ? Eigen::Vector3d::Zero()
                           : interpolateRowAt(s_k, priors.s, priors.u_dot);
        strain_at[k] = permuteStrain(v_at[k], u_at[k]);
    }

    // Strain measurements — one per estimator node.
    // The user-facing convention is that m.value is the target strain directly:
    // internally the estimator does strain_des = -m.value (l. 869) and then
    // convertStateMeanBodyInertial() negates the strain again on the way out
    // (l. 1815), so the two negations cancel and the returned state's strain
    // equals m.value at convergence. See doc/bridge_adapter_guide.md §4.3.
    for (unsigned int k = 0; k < K; ++k) {
        ContinuumRobotStateEstimator::SensorMeasurement m;
        m.type      = ContinuumRobotStateEstimator::SensorMeasurement::Strain;
        m.value     = strain_at[k];
        m.mask      << 1, 1, 1, 1, 1, 1;
        m.idx_robot = robot_idx;
        m.idx_node  = k;
        ep.measurements.push_back(m);
    }

    // Control inputs — one per segment between consecutive nodes.
    // 12x1 layout: [velocity_se3(6); acceleration_se3(6)].
    // getTransitionFunction() negates the values internally, so no sign flip.
    for (unsigned int k = 0; k + 1 < K; ++k) {
        const Eigen::Matrix<double,6,1> strain_dot_k =
            permuteStrain(vdot_at[k], udot_at[k]);

        ContinuumRobotStateEstimator::ControlInput ci;
        ci.type        = ContinuumRobotStateEstimator::ControlInput::Constant;
        ci.idx_robot   = static_cast<int>(robot_idx);
        ci.idx_segment = static_cast<int>(k);
        ci.values.resize(1);
        ci.values[0].resize(12, 1);
        ci.values[0].block(0, 0, 6, 1) = strain_at[k];
        ci.values[0].block(6, 0, 6, 1) = strain_dot_k;
        ep.control_inputs.push_back(ci);
    }

    // Initial guess — pose from the nearest Cosserat disk frame (rotated to
    // the estimator's body convention), strain from the resampled permuted
    // strain vector.
    const Eigen::Matrix4d& Ti0 = topology.Ti0.at(robot_idx);
    for (unsigned int k = 0; k < K; ++k) {
        const Eigen::Matrix4d T_cos = nearestDiskFrame(s_at[k], priors.s, diskFrames);
        const Eigen::Matrix3d R_cos = T_cos.block(0, 0, 3, 3);
        const Eigen::Vector3d p_cos = T_cos.block(0, 3, 3, 1);

        Eigen::Matrix4d T_est_body = Eigen::Matrix4d::Identity();
        T_est_body.block(0, 0, 3, 3) = R_conv * R_cos * R_convT;
        T_est_body.block(0, 3, 3, 1) = R_conv * p_cos;

        auto& node     = target.estimation_nodes[k];
        node.arclength = s_at[k];
        node.pose      = Ti0 * T_est_body;
        node.strain    = strain_at[k];
    }

    return ep;
}
