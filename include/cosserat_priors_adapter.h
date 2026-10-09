#ifndef COSSERAT_PRIORS_ADAPTER_H
#define COSSERAT_PRIORS_ADAPTER_H

#include "continuum_rod_priors.h"
#include "continuum_robot_state_estimator.h"

// Forward-declaration — keeps this header free of tdcr-modeling types so the
// estimator's public API stays model-agnostic. The .cpp is the only place
// that includes cosseratrodmodel.h, and it is compiled only when
// USE_LOCAL_TDCR=ON.
class CosseratRodModel;

// Extract auxiliary outputs from a CosseratRodModel (post-FK-convergence) and
// pack them into a ContinuumRodPriors DTO (Data Transfer Object) for consumption by the estimator.
//
// Throws std::runtime_error (propagated from the model's getters) if the last
// forwardKinematics() call did not converge or was never issued.
// Throws std::invalid_argument (from validate()) on shape inconsistency.
//
// L1, L2 are the physical lengths of the two CosseratRodModel segments —
// needed to populate the s_discrete arclengths for the junction and tip
// concentrated loads.
//
// diskFrames is the 4 x (4*N) disk-frame matrix from the same converged
// forwardKinematics() call. Each 4x4 block's 3x3 upper-left is the rod's
// body-frame orientation R(s_i) used to rotate the world-frame distributed
// loads into body frame for the K^-1 * f_in computation (priors.epsilon_in_accel).
// If the producing model does not expose disk frames, pass an empty
// Eigen::MatrixXd() — epsilon_in_accel will then be left empty and downstream
// ForceAsInput consumers will fall back to zero acceleration input.
ContinuumRodPriors priorsFromCosseratModel(const CosseratRodModel& model,
                                           double L1, double L2,
                                           const Eigen::MatrixXd& diskFrames);


// Bundle of estimator-ready inputs produced by cosseratPriorsToEstimator().
//
// IMPORTANT — what "strain" means once control inputs are used: the strain
// stored in the estimator's state is only the *bias*, i.e. the part NOT
// already carried by the velocity input. The strain that actually shapes the
// rod is (state strain) + (velocity input). The estimator never adds the
// input back, so `measurements`, the initial guess and the returned state all
// refer to the bias. `model_strain` keeps the rod model's total strain at
// each node for comparisons.
struct EstimatorPriors
{
    std::vector<ContinuumRobotStateEstimator::SensorMeasurement> measurements;
    std::vector<ContinuumRobotStateEstimator::ControlInput>       control_inputs;
    ContinuumRobotStateEstimator::SystemState                     initial_guess;
    std::vector<Eigen::Matrix<double,6,1>>                        model_strain;
};


// How the rod model's prediction is handed to the estimator:
//
//   None           : model strain as strain measurements (value = model strain),
//                    no control inputs. Initial-guess strain = model strain.
//                    (Evaluation "Method 2".)
//   StrainAsInput  : velocity input = model strain at the segment midpoint,
//                    acceleration input = 0. The state strain
//                    is then the bias, so measurements and initial-guess strain
//                    are set to 0 ("the model strain is fully explained").
//                    Do NOT also feed model-strain measurements: the rod would
//                    be built from twice the strain.
//   ForceAsInput   : velocity input = 0, acceleration input = load-driven
//                    strain derivative -K^-1 [f; l] (priors.epsilon_in_accel)
//                    plus each interior strain jump -K^-1 [F; L]
//                    (priors.epsilon_jump_discrete) spread over the segment
//                    that starts at the load. This
//                    only tells the estimator how the strain CHANGES along the
//                    rod; its level must come from measurements.
//                    Velocity input is zero, so bias = total strain:
//                    measurements / initial-guess strain = model strain.
//                    (Evaluation "Method 3", which uses only the inputs.)
enum class ControlInputMode {
    None,
    StrainAsInput,
    ForceAsInput
};

// Convert a ContinuumRodPriors (stored in the producing model's body frame)
// into estimator-ready measurements, control inputs, and an initial guess.
//
// Handles:
//   - Arclength resampling from the rod-model grid to the estimator's K nodes.
//   - Body-frame convention rotation (Cosserat z-forward -> estimator x-forward)
//     via the cyclic permutation R_conv = [[0,0,1],[1,0,0],[0,1,0]].
//   - Sign conventions required by assembleMeasurementTerms() and
//     getTransitionFunction() — see doc/cosserat_integration.md §3.2.
//
// Parameters:
//   priors     — DTO (Data Transfer Object) produced by a rod-model adapter (e.g. priorsFromCosseratModel).
//   topology   — estimator RobotTopology (K, L, Ti0, ...).
//   robot_idx  — which robot in the topology this prior drives.
//   diskFrames — the 4 x (4*N) disk-frame matrix from the rod model; 4x4 blocks
//                correspond 1-to-1 with priors.s samples.
//
// Multi-robot note: only robots[robot_idx].estimation_nodes is populated in
// the returned initial_guess. If topology.N > 1 and you assign the result to
// options.custom_guess directly, the other robots will have empty
// estimation_nodes — the estimator's convertStateMeanBodyInertial() will
// access them out-of-bounds. For multi-robot setups call this once per robot
// and merge the robots[] vectors yourself before assigning to custom_guess.
//
// `mode` selects how the prediction is packed (see ControlInputMode above).
// Default None: strain measurements only, no control inputs.
//
// Throws std::invalid_argument on shape mismatches or out-of-range indices.
EstimatorPriors cosseratPriorsToEstimator(
    const ContinuumRodPriors&                            priors,
    const ContinuumRobotStateEstimator::RobotTopology&   topology,
    unsigned int                                         robot_idx,
    const Eigen::MatrixXd&                               diskFrames,
    ControlInputMode                                     mode = ControlInputMode::None);

#endif // COSSERAT_PRIORS_ADAPTER_H
