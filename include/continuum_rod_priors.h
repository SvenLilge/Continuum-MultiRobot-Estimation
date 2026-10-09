#ifndef CONTINUUM_ROD_PRIORS_H
#define CONTINUUM_ROD_PRIORS_H

#include <Eigen/Core>
#include <vector>

// Model-agnostic prior / snapshot fed into ContinuumRobotStateEstimator.
// Pure Eigen + STL; no dependency on any specific rod-model library.
// Any model (Cosserat, PCC, sub-segment Cosserat, analytical) can produce one
// — the consumer (estimator) never needs to know which model generated it.
//
// Dimensions are co-aligned: all N-sample matrices share the same row index
// into `s`. Discrete loads are separate, keyed by their own s_discrete entries.
//
// Optional fields (the `_dot`, `_internal`, `_dist`, and discrete arrays) may
// be left default-constructed (empty) by a model that cannot produce them.
// `validate()` treats those as acceptable; if populated, sizes must match.
struct ContinuumRodPriors
{
    // Arclength discretization (N samples).
    Eigen::VectorXd s;                         // (N,)   monotonic non-decreasing

    // Kinematic state (body frame). Required.
    Eigen::MatrixXd v;                         // (N,3)  linear strain
    Eigen::MatrixXd u;                         // (N,3)  angular strain / curvature

    // Derivatives along s (body frame). Optional.
    Eigen::MatrixXd v_dot;                     // (N,3)
    Eigen::MatrixXd u_dot;                     // (N,3)

    // Internal force / moment resultants (world frame). Optional.
    Eigen::MatrixXd n_internal;                // (N,3)
    Eigen::MatrixXd m_internal;                // (N,3)

    // Distributed applied load densities (world frame). Optional.
    Eigen::MatrixXd f_dist;                    // (N,3)
    Eigen::MatrixXd l_dist;                    // (N,3)

    // Strain derivative caused by the distributed loads:
    // -K^-1 * [f_body; l_body], K = diag(K_se, K_bt). This is the
    // acceleration-level control input, written in the physical sign
    // convention the estimator API expects.
    // Producer body frame (e.g. Cosserat z-forward); the adapter to the
    // estimator applies the body-convention permutation. Optional; empty or
    // (N,6). Populated by a producer that exposes both K and distributed
    // wrenches (e.g. priorsFromCosseratModel).
    Eigen::MatrixXd epsilon_in_accel;          // (N,6)

    // Concentrated loads at specific arclengths (world frame). Optional.
    // For tendon-driven rods, typical entries are at segment junctions and
    // at the tip where tendons terminate. All three arrays must have the
    // same length M.
    std::vector<double>          s_discrete;
    std::vector<Eigen::Vector3d> F_discrete;
    std::vector<Eigen::Vector3d> L_discrete;

    // Strain jump across each concentrated load:
    // eps(s+) - eps(s-) = -K^-1 * [R^T F; R^T L].
    // Producer body frame. Optional; empty or one entry per s_discrete.
    std::vector<Eigen::Matrix<double,6,1>> epsilon_jump_discrete;

    // Minimal shape-consistency check. Throws std::invalid_argument on mismatch.
    // Called at the adapter boundary so garbage never reaches the estimator.
    void validate() const;
};

#endif // CONTINUUM_ROD_PRIORS_H
