#include "cosserat_priors_adapter.h"
#include "cosseratrodmodel.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

ContinuumRodPriors priorsFromCosseratModel(const CosseratRodModel& model,
                                           double L1, double L2,
                                           const Eigen::MatrixXd& diskFrames)
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

    // Load-driven strain derivative: -K^-1 * [f_body; l_body]
    // (from K eps' = -N(eps) - [R^T f; R^T l]). This is the acceleration
    // input in the physical (API) sign convention; the estimator's internal
    // convention, where strain has the opposite sign, uses +K^-1 f_in. Needs
    // disk-frame rotations R(s_i) to bring world-frame loads into the body
    // frame. If diskFrames is empty we leave epsilon_in_accel empty so
    // ForceAsInput consumers can detect the missing data.
    const Eigen::Index N = p.s.size();
    const bool have_frames =
        diskFrames.rows() == 4 && diskFrames.cols() == 4 * N;
    if (p.f_dist.size() != 0 && have_frames) {
        const Eigen::Matrix3d Kse_inv = model.getKse().inverse();
        const Eigen::Matrix3d Kbt_inv = model.getKbt().inverse();
        p.epsilon_in_accel.resize(N, 6);
        const bool have_l = (p.l_dist.size() != 0);
        for (Eigen::Index i = 0; i < N; ++i) {
            const Eigen::Matrix3d R_i = diskFrames.block(0, 4 * i, 3, 3);
            const Eigen::Vector3d f_w = p.f_dist.row(i).transpose();
            const Eigen::Vector3d f_b = R_i.transpose() * f_w;
            const Eigen::Vector3d l_b = have_l
                ? Eigen::Vector3d(R_i.transpose() * p.l_dist.row(i).transpose())
                : Eigen::Vector3d::Zero();
            p.epsilon_in_accel.row(i).head<3>() = -(Kse_inv * f_b).transpose();
            p.epsilon_in_accel.row(i).tail<3>() = -(Kbt_inv * l_b).transpose();
        }

        // Strain jump across each tendon termination:
        // eps(s+) - eps(s-) = -K^-1 [R^T F; R^T L], R at the load's sample.
        for (std::size_t j = 0; j < p.s_discrete.size(); ++j) {
            Eigen::Index i_near = 0;
            for (Eigen::Index i = 1; i < N; ++i)
                if (std::abs(p.s(i) - p.s_discrete[j]) < std::abs(p.s(i_near) - p.s_discrete[j]))
                    i_near = i;
            const Eigen::Matrix3d R_j = diskFrames.block(0, 4 * i_near, 3, 3);
            Eigen::Matrix<double,6,1> jump;
            jump << -(Kse_inv * R_j.transpose() * p.F_discrete[j]),
                    -(Kbt_inv * R_j.transpose() * p.L_discrete[j]);
            p.epsilon_jump_discrete.push_back(jump);
        }
    }

    p.validate();
    return p;
}


// ============================================================================
// cosseratPriorsToEstimator — frame-permute + resample a ContinuumRodPriors
// into estimator-ready measurements, control inputs, and initial guess.
// Math / conventions: see doc/cosserat_integration.md §3.2.
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
    const Eigen::MatrixXd&                               diskFrames,
    ControlInputMode                                     mode)
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
    std::vector<Eigen::Vector3d>           v_at(K), u_at(K);
    std::vector<Eigen::Matrix<double,6,1>> strain_at(K);
    std::vector<double>                    s_at(K);

    for (unsigned int k = 0; k < K; ++k) {
        const double s_k = (static_cast<double>(k) / (K - 1)) * L;
        s_at[k]      = s_k;
        v_at[k]      = interpolateRowAt(s_k, priors.s, priors.v);
        u_at[k]      = interpolateRowAt(s_k, priors.s, priors.u);
        strain_at[k] = permuteStrain(v_at[k], u_at[k]);
    }

    ep.model_strain = strain_at;

    // State strain implied by the mode. With a velocity input the estimator's
    // state strain is only the bias (total - velocity input); StrainAsInput puts
    // the whole model strain on the input, so the bias target is zero.
    // See the comment on EstimatorPriors in the header.
    const bool strain_on_input = (mode == ControlInputMode::StrainAsInput);
    std::vector<Eigen::Matrix<double,6,1>> state_strain_at = strain_at;
    if (strain_on_input)
        for (auto& e : state_strain_at) e.setZero();

    // Strain measurements — one per estimator node.
    // The user-facing convention is that m.value is the target strain directly:
    // internally the estimator does strain_des = -m.value (l. 869) and then
    // convertStateMeanBodyInertial() negates the strain again on the way out
    // (l. 1815), so the two negations cancel and the returned state's strain
    // equals m.value at convergence. See doc/cosserat_integration.md §3.2.
    for (unsigned int k = 0; k < K; ++k) {
        ContinuumRobotStateEstimator::SensorMeasurement m;
        m.type      = ContinuumRobotStateEstimator::SensorMeasurement::Strain;
        m.value     = state_strain_at[k];
        m.mask      << 1, 1, 1, 1, 1, 1;
        m.idx_robot = robot_idx;
        m.idx_node  = k;
        ep.measurements.push_back(m);
    }

    // Control inputs — one per segment between consecutive nodes.
    // 12x1 layout: [velocity_se3(6); acceleration_se3(6)].
    // getTransitionFunction() negates the values internally, so no sign flip.
    //
    // Modes:
    //   None          : no control inputs.
    //   StrainAsInput : velocity = model strain at the segment midpoint,
    //                   acceleration = 0. The acceleration input is
    //                   NOT the derivative of the velocity input, so putting
    //                   strain_dot there would count the change twice.
    //   ForceAsInput  : velocity = 0, acceleration = K^-1 * f_in
    //                   (= priors.epsilon_in_accel, resampled & permuted).

    // Pre-resample epsilon_in_accel (Cosserat body frame) at each estimator node
    // for ForceAsInput mode; permute to estimator body frame at the same time.
    std::vector<Eigen::Matrix<double,6,1>> eps_accel_at(
        K, Eigen::Matrix<double,6,1>::Zero());
    if (mode == ControlInputMode::ForceAsInput && priors.epsilon_in_accel.size() != 0) {
        const Eigen::MatrixXd eps_f = priors.epsilon_in_accel.leftCols(3);
        const Eigen::MatrixXd eps_l = priors.epsilon_in_accel.rightCols(3);
        for (unsigned int k = 0; k < K; ++k) {
            const Eigen::Vector3d top = interpolateRowAt(s_at[k], priors.s, eps_f);
            const Eigen::Vector3d bot = interpolateRowAt(s_at[k], priors.s, eps_l);
            eps_accel_at[k] = permuteStrain(top, bot);
        }

        // Concentrated loads inside the rod make the strain jump. A segment
        // input is constant, so spread each jump over the segment that starts
        // at (or contains) the load: input += jump / segment length. A load at
        // the tip has no segment after it; it only sets the tip boundary
        // condition, which an input on the strain derivative cannot express.
        for (std::size_t j = 0; j < priors.epsilon_jump_discrete.size(); ++j) {
            const double s_j = priors.s_discrete[j];
            for (unsigned int k = 0; k + 1 < K; ++k) {
                if (s_j >= s_at[k] - 1e-9 && s_j < s_at[k + 1] - 1e-9) {
                    const auto& d = priors.epsilon_jump_discrete[j];
                    eps_accel_at[k] += permuteStrain(d.head<3>(), d.tail<3>())
                                       / (s_at[k + 1] - s_at[k]);
                    break;
                }
            }
        }
    }

    for (unsigned int k = 0; mode != ControlInputMode::None && k + 1 < K; ++k) {
        ContinuumRobotStateEstimator::ControlInput ci;
        ci.type        = ContinuumRobotStateEstimator::ControlInput::Constant;
        ci.idx_robot   = static_cast<int>(robot_idx);
        ci.idx_segment = static_cast<int>(k);
        ci.values.resize(1);
        ci.values[0].resize(12, 1);
        if (mode == ControlInputMode::ForceAsInput) {
            ci.values[0].block(0, 0, 6, 1).setZero();
            ci.values[0].block(6, 0, 6, 1) = eps_accel_at[k];
        } else { // StrainAsInput
            // Constant over the segment, so sample at its midpoint. This also
            // picks the correct side of a strain jump at a node (e.g. the
            // segment junction, where the node sample holds the pre-jump value).
            const double s_mid = 0.5 * (s_at[k] + s_at[k + 1]);
            ci.values[0].block(0, 0, 6, 1) = permuteStrain(
                interpolateRowAt(s_mid, priors.s, priors.v),
                interpolateRowAt(s_mid, priors.s, priors.u));
            ci.values[0].block(6, 0, 6, 1).setZero();
        }
        ep.control_inputs.push_back(ci);
    }

    // Initial guess — pose from the nearest Cosserat disk frame (rotated to
    // the estimator's body convention), strain = the state strain implied by
    // the mode (model strain, or zero bias for StrainAsInput).
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
        node.strain    = state_strain_at[k];
    }

    return ep;
}
