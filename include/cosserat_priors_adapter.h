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
ContinuumRodPriors priorsFromCosseratModel(const CosseratRodModel& model,
                                           double L1, double L2);


// Bundle of estimator-ready inputs produced by cosseratPriorsToEstimator().
struct EstimatorPriors
{
    std::vector<ContinuumRobotStateEstimator::SensorMeasurement> measurements;
    std::vector<ContinuumRobotStateEstimator::ControlInput>       control_inputs;
    ContinuumRobotStateEstimator::SystemState                     initial_guess;
};

// Convert a ContinuumRodPriors (stored in the producing model's body frame)
// into estimator-ready measurements, control inputs, and an initial guess.
//
// Handles:
//   - Arclength resampling from the rod-model grid to the estimator's K nodes.
//   - Body-frame convention rotation (Cosserat z-forward -> estimator x-forward)
//     via the cyclic permutation R_conv = [[0,0,1],[1,0,0],[0,1,0]].
//   - Sign conventions required by assembleMeasurementTerms() and
//     getTransitionFunction() — see doc/bridge_adapter_guide.md §4.3 and §4.4.
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
// Throws std::invalid_argument on shape mismatches or out-of-range indices.
EstimatorPriors cosseratPriorsToEstimator(
    const ContinuumRodPriors&                            priors,
    const ContinuumRobotStateEstimator::RobotTopology&   topology,
    unsigned int                                         robot_idx,
    const Eigen::MatrixXd&                               diskFrames);

#endif // COSSERAT_PRIORS_ADAPTER_H
