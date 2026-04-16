#ifndef COSSERAT_PRIORS_ADAPTER_H
#define COSSERAT_PRIORS_ADAPTER_H

#include "continuum_rod_priors.h"

// Forward-declaration — keeps this header free of tdcr-modeling types so the
// estimator's public API stays model-agnostic. The .cpp is the only place
// that includes cosseratrodmodel.h, and it is compiled only when
// USE_LOCAL_TDCR=ON.
class CosseratRodModel;

// Extract auxiliary outputs from a CosseratRodModel (post-FK-convergence) and
// pack them into a ContinuumRodPriors DTO for consumption by the estimator.
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

#endif // COSSERAT_PRIORS_ADAPTER_H
