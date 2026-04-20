# Cosserat-Informed State Estimation — Integration Report

**Author:** Matteo Guidi
**Date:** 2026-04-19
**Scope:** Feeding physically-consistent outputs from a Cosserat rod model
(sibling repo [`tdcr-modeling`](../../tdcr-modeling)) into the SE(3)
continuum-robot state estimator in this repo, and quantifying the resulting
improvement.

---

## 1. Executive summary

The estimator in [`continuum_robot_state_estimator.{h,cpp}`](../include/continuum_robot_state_estimator.h)
solves a MAP problem over a Gaussian-process prior on SE(3) backbone states,
using user-supplied measurements and control inputs. Before this work, the
prior was physics-agnostic (straight-rod Gaussian) and the measurements had to
come from sensors (FBG, pose markers, ...). For sensor-sparse or
sensor-absent configurations, the estimator had no way to incorporate the
well-known tendon-force / elasticity physics of a continuum robot.

This integration lets the estimator consume the **full auxiliary output of a
Cosserat rod model** (strains, strain rates, disk frames) as (i) strain
measurements, (ii) GP-prior control inputs, and (iii) an initial guess. With
nothing else changed in the estimator, the mean per-node strain error vs the
Cosserat reference drops from **0.224 rad/m (no priors)** to
**2.34 × 10⁻⁴ rad/m (with priors)** — a **≈ 957× reduction**. The full
6-DOF per-node strain field is reproduced sign-for-sign.

The integration is **additive**: no lines in `continuum_robot_state_estimator.{h,cpp}`
were changed (one `#include <cassert>` aside, required to make Release builds
compile). Standalone builds still run without GSL, without `tdcr-modeling`, and
without any behaviour change. The Cosserat path is opt-in via
`-DUSE_LOCAL_TDCR=ON`.

---

## 2. What the relevant commits contain

Two commits land the core of this work — one on each side of the bridge.

### 2.1 `tdcr-modeling` @ [`86f3b90`](../../tdcr-modeling) — auxiliary outputs from the Cosserat model

> *Add auxiliary output getters and validation tests for CosseratRodModel*

Before this commit, `CosseratRodModel::forwardKinematics()` returned only
disk frames. The physics state computed inside the shooting loop (body-frame
strains `v, u`, their s-derivatives `v̇, u̇`, internal resultants `n, m`,
distributed applied loads `f_dist, l_dist`, and discrete tendon loads at the
segment junction / tip) was discarded on return.

The commit adds:

| New getter | Shape | Meaning (Rucker 2011 §II) |
|---|---|---|
| `getArclengthSamples()` | `(N,)` | Sample arclengths `s_j` (irregular, 21 pts for the default two-segment rod) |
| `getStrainV()`, `getStrainU()` | `(N, 3)` | Linear & angular strains in the body frame |
| `getStrainVDot()`, `getStrainUDot()` | `(N, 3)` | Arclength derivatives (needed by the GP prior) |
| `getInternalForce()`, `getInternalMoment()` | `(N, 3)` | Static equilibrium resultants (world frame) |
| `getDistributedForce()`, `getDistributedMoment()` | `(N, 3)` | Applied load densities per unit arclength |
| `getDiscreteLoads(F_junction, L_junction, F_tip, L_tip)` | 4 × 3-vectors | Concentrated tendon terminations |
| `hasAuxOutputs()` | `bool` | Flip to `true` only on a converged FK solve |

Populating these members required rerunning the shooting loop's local
computations *once more* on the converged solution. The commit also adds a
validation test suite (`test_cosseratrodmodel_aux.cpp`, 322 lines) covering
zero-tension, single-tendon, equilibrium, and continuation-mode edge cases.

Net effect: the Cosserat model went from "black-box IK/FK" to **"physics-aware
data producer"** — exactly what the estimator needs to regularize its MAP
problem.

### 2.2 `Continuum-MultiRobot-Estimation` @ [`068fd4e`](../CMakeLists.txt) — model-agnostic pipeline scaffolding

> *feat: Implement Cosserat rod model integration with continuum state estimator*

This commit lays the integration substrate on the estimator side:

| New / changed | Role |
|---|---|
| [`include/continuum_rod_priors.h`](../include/continuum_rod_priors.h) | **DTO (Data Transfer Object)** carrying `s, v, u, v̇, u̇, n, m, f_dist, l_dist, s_discrete, F_discrete, L_discrete`. Pure Eigen + STL, header-only. Any rod model (Cosserat, PCC, analytical, ...) can produce one. |
| [`src/continuum_rod_priors.cpp`](../src/continuum_rod_priors.cpp) | `validate()` — shape consistency check called at every adapter boundary. |
| [`include/cosserat_priors_adapter.h`](../include/cosserat_priors_adapter.h) | Forward-declares `CosseratRodModel`; exposes `priorsFromCosseratModel()`. Public header, stays tdcr-free. |
| [`src/bridge/cosserat_priors_adapter.cpp`](../src/bridge/cosserat_priors_adapter.cpp) | Includes `cosseratrodmodel.h`; implements the Cosserat → DTO (Data Transfer Object) copy. Lives in `src/bridge/` **deliberately outside** `file(GLOB src/*.cpp)` so it is never compiled when `USE_LOCAL_TDCR=OFF`. |
| [`src/examples/cosserat_estimator_driver.cpp`](../src/examples/cosserat_estimator_driver.cpp) | End-to-end smoke test — runs Cosserat FK and prints the DTO (Data Transfer Object). Proves the cross-repo boundary works. |
| [`CMakeLists.txt`](../CMakeLists.txt) | Adds `option(USE_LOCAL_TDCR)` (default OFF), conditional `add_subdirectory(../tdcr-modeling/c++)`, conditional new target. No effect on standalone builds. |

At this point the estimator could **see** Cosserat priors but did not yet
**use** them — the driver printed the DTO (Data Transfer Object) and exited.

### 2.3 Current branch work — connecting the DTO (Data Transfer Object) to `computeStateEstimate()`

The uncommitted changes on [`add-control-inputs`](../CMakeLists.txt) complete
the consumption side:

| File | Change |
|---|---|
| [`include/cosserat_priors_adapter.h`](../include/cosserat_priors_adapter.h) | New `EstimatorPriors` struct (measurements + control inputs + initial guess) and new `cosseratPriorsToEstimator()` function. |
| [`src/bridge/cosserat_priors_adapter.cpp`](../src/bridge/cosserat_priors_adapter.cpp) | `cosseratPriorsToEstimator()` body: body-frame permutation, arclength resampling to the estimator's K-node grid, initial-guess pose construction. |
| [`src/examples/cosserat_estimator_driver.cpp`](../src/examples/cosserat_estimator_driver.cpp) | Rewritten: actually runs `computeStateEstimate()`, prints a per-node strain comparison table, honors `--no-control-inputs` and `--visualize`. |
| [`config/6_cosserat_priors.yaml`](../config/6_cosserat_priors.yaml) | New scenario: single robot, `K=11` nodes over `L=0.20 m`, strain noise `R_v = R_u = 0.01`, `initial_guess: Custom`, empty measurement / control-input lists (filled at runtime). |
| [`src/tests/test_cosserat_adapter.cpp`](../src/tests/test_cosserat_adapter.cpp) | Four validation tests: frame-permutation unit tests, resampling accuracy, end-to-end convergence, and a prior-vs-no-prior A/B. |
| [`CMakeLists.txt`](../CMakeLists.txt) | Gates `test_cosserat_adapter` and `cosserat_estimator_driver` inside `if(USE_LOCAL_TDCR)`. |

**The estimator itself was not modified.** Integration is purely additive
through the existing public API (`SensorMeasurement`, `ControlInput`,
`Options::custom_guess`).

---

## 3. Why this is a meaningful improvement

### 3.1 What the estimator does (Lilge et al. 2024)

`ContinuumRobotStateEstimator::computeStateEstimate()` performs MAP state
estimation over a factor graph with three kinds of factors:

1. **GP prior factor** — a Gauss-Markov prior on the backbone state evolution
   along `s` parameterized by `Q_c` and a per-segment expected velocity
   (control input).
2. **Measurement factors** — sensor constraints (pose, strain, FBG
   wavelengths, ...).
3. **Coupling factors** — multi-robot kinematic constraints.

Without Cosserat priors the user has only two levers: hand-tune `Q_c`, or
install more sensors. Both are unsatisfying — `Q_c` isn't physics, and
sensors are expensive / invasive.

### 3.2 What Cosserat priors add

Given tendon tensions and segment geometry, the Cosserat model solves the
full quasi-static equilibrium of the rod and produces the physically-correct
strain field `v(s), u(s)`. Feeding these into the estimator:

| Estimator input | Source | Physical meaning |
|---|---|---|
| `SensorMeasurement::Strain` at each of K nodes | Permuted `v(s_k), u(s_k)` from the Cosserat model | Pseudo-measurements "the rod *should* have this strain here". Pins each node's SE(3) body velocity. |
| `ControlInput::Constant` on each of K−1 segments | `[v_k, u_k, v̇_k, u̇_k]` (12×1) | Sets the GP prior's expected pose-tangent per unit arclength + its expected rate of change. Makes the prior's mean trajectory match the physics. |
| `Options::custom_guess` | Cosserat disk frames (permuted to estimator body convention), plus strain at each node | Seeds Newton's method near the physics solution, avoiding local minima / slow start from a straight rod. |

The net effect is that the estimator's MAP problem now has a **physics-consistent
mean** and **dense pseudo-observations** along the rod, not just a
straight-rod mean.

### 3.3 Two subtle correctness obligations

Body-frame convention differs between the two codebases. The adapter corrects
this once, before data enters the estimator:

| | Cosserat | Estimator |
|---|---|---|
| Backbone axis | local **z** | local **x** |
| Straight rod strain | `v* = [0, 0, 1]` | `ν = [1, 0, 0]` |

The mapping is the cyclic permutation `R_conv = [[0,0,1],[1,0,0],[0,1,0]]`
applied to all 3-vectors (strains, disk-frame positions, and both sides of
rotation matrices). Unit tests A and B in
[`test_cosserat_adapter.cpp`](../src/tests/test_cosserat_adapter.cpp) verify
this on six canonical cases (bend-about-x, bend-about-y, torsion, axial,
two shears) to `1e-12`.

The estimator also has **two internal sign flips on strain** that cancel
(`strain_des = -m.value` at [l. 869](../src/continuum_robot_state_estimator.cpp#L869)
and a second negation in [`convertStateMeanBodyInertial()`](../src/continuum_robot_state_estimator.cpp#L1815)),
and **one sign flip on control inputs** that does not cancel
([l. 2416](../src/continuum_robot_state_estimator.cpp#L2416)). The adapter
respects these — `m.value` is set to the user-facing target strain directly,
and `ControlInput.values` are stored without any manual negation. A first
draft of the adapter got the measurement sign wrong; the strain field came
out exactly negated, which is how we caught the audit.

---

## 4. Experimental protocol

All runs use [`config/6_cosserat_priors.yaml`](../config/6_cosserat_priors.yaml):

- **Geometry**: single robot, two 0.10 m segments (`L = 0.20 m`), `K = 11`
  estimator nodes uniformly spaced at `Δs = 0.02 m`. Cosserat sample
  arclengths align with estimator node arclengths so resampling never
  interpolates across the junction discontinuity at `s = 0.10`.
- **Boundary conditions**: base pose locked (world identity). Base strain
  locked only in the A/B comparison so that the "no-prior" run is well-posed
  with no measurements.
- **Loading**: tendon tensions `q = [0.5, 0.2, 0.0, 0.3, 0.0, 0.1] N`
  (non-trivial 3D bending in both segments), no external point / distributed
  loads.
- **Noise model**: `R_v = R_u = 0.01` (strain), `R_p = 0.002, R_o = 0.05`
  (unused here, no pose measurements). Low strain noise = estimator trusts
  the Cosserat pseudo-measurements.
- **Optimizer**: Newton, `max_iterations = 100`, `tol = 1e-3`,
  `kirchhoff_rods = true`.

**Three runs are compared**, all using the same topology and hyperparameters:

| Run | Measurements | Control inputs | Initial guess |
|---|---|---|---|
| **A — No priors** | `[]` | `[]` | `Straight` (built-in) |
| **B — Priors, no control inputs** (`--no-control-inputs`) | Cosserat strains (11 nodes) | `[]` | Cosserat disk frames |
| **C — Full priors** | Cosserat strains (11 nodes) | Cosserat `[v, u, v̇, u̇]` (10 segments) | Cosserat disk frames |

Run A is the *no-physics baseline*. Run B isolates the contribution of the
pseudo-measurements + initial guess. Run C adds the GP-prior control inputs
(the "full Cosserat MAP").

Ground truth is the Cosserat strain field itself. The metric is the mean
per-node L¹ strain error (6 components averaged, 11 nodes averaged):

```
mean_abs_strain_diff = (1 / (6 K)) Σ_k Σ_c |ν_est[k, c] − ν_cos[k, c]|
```

---

## 5. Results

### 5.1 Headline — prior vs no-prior (Test D from `test_cosserat_adapter`)

| Run | Mean per-node strain error | Ratio vs no-prior |
|---|---|---|
| A — no priors (straight-rod MAP, no measurements) | **2.24 × 10⁻¹ rad/m** | 1.00× (baseline) |
| C — full Cosserat priors | **2.34 × 10⁻⁴ rad/m** | **957× smaller** |

Raw log (from `./examples/test_cosserat_adapter config/6_cosserat_priors.yaml`):

```
[RUN]  D. Prior vs no-prior
  mean |diff|  no-prior = 0.223755
  mean |diff|  prior    = 0.000233825
[OK]   D. Prior vs no-prior
```

The test's built-in sanity gate (`diff_prior * 10 < diff_none`) passes with
**96× headroom**: the margin is not driven by the test's tolerance but by
a real physics gap.

### 5.2 Per-component strain agreement (full driver run, Run C)

From `./examples/cosserat_estimator_driver config/6_cosserat_priors.yaml`:

```
Per-node strain comparison (Cosserat permuted to estimator body frame)

  Cosserat prior
    k   s[m]        nu_1      nu_2      nu_3      om_1      om_2      om_3
    0  0.000    1.0000    0.0000   -0.0000   -0.0000   -1.5958    0.2940
    1  0.020    1.0000    0.0000   -0.0000   -0.0000   -1.5958    0.2940
    ...
    5  0.100    1.0000    0.0000   -0.0000   -0.0000   -1.5958    0.2940
    6  0.120    1.0000    0.0000   -0.0000   -0.0000   -0.5093   -0.1764
    ...
   10  0.200    1.0000    0.0000   -0.0000    0.0000   -0.5093   -0.1764

  Estimator result
    k   s[m]        nu_1      nu_2      nu_3      om_1      om_2      om_3
    0  0.000    1.0000    0.0000   -0.0000   -0.0000   -1.5958    0.2940
    ...
    5  0.100    1.0000    0.0000   -0.0000    0.0000   -1.5904    0.2917
    6  0.120    1.0000    0.0000   -0.0000    0.0000   -0.5147   -0.1741
    ...
   10  0.200    1.0000    0.0000   -0.0000    0.0000   -0.5093   -0.1764

  Max |diff| per component across all nodes:
    nu_1=  0.00e+00  nu_2=  0.00e+00  nu_3=  0.00e+00
    om_1=  2.81e-08  om_2=  5.36e-03  om_3=  2.32e-03
  Max |diff|_inf overall:    5.358e-03
```

Observations:

- **Linear strain (ν₁, ν₂, ν₃) agreement is exact to machine precision** — the
  Kirchhoff-rod constraint (`kirchhoff_rods: true`) pins ν to `[1, 0, 0]` and
  both sides agree.
- **Torsion (ω₁) matches to 3 × 10⁻⁸** — essentially noise.
- **Bending (ω₂, ω₃) matches to 6 × 10⁻³ rad/m at worst**, and the worst two
  nodes are the two nodes straddling the segment junction at `s = 0.10`
  (k = 5 and k = 6). At those nodes, the Cosserat strain field is *discontinuous*
  (tendons terminating at the junction change the curvature step-wise), while
  the GP prior in the estimator prefers smoothness. The smoothing bleeds a
  few percent into the two adjacent nodes — everywhere else, agreement is at
  ~10⁻⁸ rad/m.

### 5.3 Convergence characteristics

| Run | Initial cost | Final cost | Newton iterations |
|---|---|---|---|
| B — priors, no control inputs | 1.402 × 10² | 1.402 × 10² | 1 |
| C — priors + control inputs | 3.579 × 10⁴ | 1.090 × 10² | 2 |

Run B's initial-cost ≈ final-cost is telling: the Cosserat disk-frame initial
guess is *already* the estimator's minimum under this measurement set, so
Newton has nothing to refine beyond one step's worth of strain matching.

Run C's initial cost is larger because the GP-prior cost now penalizes
deviation from the Cosserat-specified tangent field, which the straight-guess
component of the initial state violates for two iterations — still trivial
convergence.

In Run A (no priors, no measurements), the cost is flat and the estimator
stays at the straight initial guess, producing the 0.224 rad/m error above.

### 5.4 Qualitative physical meaning of the error reduction

- The Cosserat-predicted bending curvature in segment 1 is
  `‖ω‖ = √(1.5958² + 0.2940²) ≈ 1.62 rad/m`. Integrated over 0.10 m that is
  **9.3°** of bend in segment 1 alone — the rod clearly bends. The no-prior
  run predicts a **straight rod** and therefore misses this deformation
  entirely; the ≈ 0.22 rad/m mean error quantifies exactly that miss.
- The prior-informed estimator recovers the full bending field to **sub-0.01
  rad/m average error**, which corresponds to a sub-degree deviation over
  a 10 cm segment — well below the physical repeatability of any real tendon
  actuation.

### 5.5 Unit-test suite status

Full output of `./examples/test_cosserat_adapter config/6_cosserat_priors.yaml`:

```
=== Cosserat adapter + integration tests ===

[RUN]  A. Frame permutation
[OK]   A. Frame permutation

[RUN]  B. Resampling (linear strain)
[OK]   B. Resampling (linear strain)

[RUN]  C. End-to-end convergence
  max |diff|_inf = 0.00535762
[OK]   C. End-to-end convergence

[RUN]  D. Prior vs no-prior
  mean |diff|  no-prior = 0.223755
  mean |diff|  prior    = 0.000233825
[OK]   D. Prior vs no-prior

=== Results: 4 passed, 0 failed ===
```

| Test | What it validates |
|---|---|
| **A** | The six canonical body-frame permutations (straight, bend-x, bend-y, torsion-z, shear-x, shear-y) map to the correct estimator strain to `1e-12`. |
| **B** | A linear strain field `v(s) = [0, 0, 1+s]` at 5 Cosserat samples is resampled onto a denser 9-node estimator grid to `1e-12` — confirms interpolation math is right. |
| **C** | End-to-end Cosserat FK + adapter + estimator reproduce the Cosserat strain field to `max |diff|_inf = 5.4e-3` on the real two-segment 6-tendon scenario. |
| **D** | Prior-informed estimate is **957× closer** to the Cosserat reference than the uninformed one (baseline: 0.224 rad/m, with priors: 0.000234 rad/m). |

---

## 6. Design properties preserved

| Invariant | Status |
|---|---|
| Standalone build (`USE_LOCAL_TDCR=OFF`) produces `continuum_robot_viewer`, `test_config_loader`, `test_estimation` with **identical behavior** to pre-integration. | Verified. |
| Public estimator headers do not `#include` any `tdcr-modeling` type. | Verified — adapter header forward-declares `CosseratRodModel`; only the `.cpp` in `src/bridge/` includes `cosseratrodmodel.h`. |
| `continuum_robot_state_estimator.{h,cpp}` logic is untouched. | Verified — only a missing `#include <cassert>` was added, which was a pre-existing Release-build transitive-include issue, fixed with explicit authorization. |
| `src/bridge/cosserat_priors_adapter.cpp` is never compiled in standalone builds. | Verified — it lives outside `file(GLOB src/*.cpp)` and is only added to the combined-build targets inside the `if (USE_LOCAL_TDCR)` block. |
| No git submodules, no copies of tdcr sources, no `find_package(tdcr_modeling)` (no upstream Config.cmake exists). | Verified — integration is via `add_subdirectory(../tdcr-modeling/c++)` only when opted in. |
| DTO (Data Transfer Object — `ContinuumRodPriors`) is model-agnostic — a PCC or analytical adapter could produce one with no estimator-side changes. | Verified — the DTO (Data Transfer Object) is pure Eigen + STL, no Cosserat-specific fields. |

---

## 7. How to reproduce

From the repo root:

```bash
# 1. Combined build (opt-in). Ninja is preferred (see CLAUDE.md §2).
rm -rf build && mkdir build && cd build
cmake -G Ninja -DCMAKE_BUILD_TYPE=Release -DUSE_LOCAL_TDCR=ON ..
cmake --build .
cd ..

# 2. Headline A/B/C experiment and unit tests:
./examples/test_cosserat_adapter config/6_cosserat_priors.yaml
# Expect: 4/4 passed, including "mean |diff| no-prior = 0.223755,
#         mean |diff| prior = 0.000233825"

# 3. Full per-node strain comparison table (Run C):
./examples/cosserat_estimator_driver config/6_cosserat_priors.yaml

# 4. Run B (priors but no GP-prior control inputs):
./examples/cosserat_estimator_driver config/6_cosserat_priors.yaml --no-control-inputs

# 5. Visualization (requires a display):
./examples/cosserat_estimator_driver config/6_cosserat_priors.yaml --visualize
```

Standalone regression check (no GSL, no tdcr):

```bash
rm -rf build && mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build . -j
./examples/test_config_loader && ./examples/test_estimation
# Expect: all pass, no GSL probe in the CMake output.
```

---

## 8. Reference paths

- Integration plan (phased, with verified numbers): [`doc/cosserat_integration_plan.md`](cosserat_integration_plan.md)
- Pipeline-phase summary: [`doc/changes_and_new_files.md`](changes_and_new_files.md)
- Cosserat model extension plan (sibling repo): [`../../tdcr-modeling/doc/cosseratrodmodel_extension_plan.md`](../../tdcr-modeling/doc/cosseratrodmodel_extension_plan.md)
- Adapter header: [`include/cosserat_priors_adapter.h`](../include/cosserat_priors_adapter.h)
- Adapter body (only built with `USE_LOCAL_TDCR=ON`): [`src/bridge/cosserat_priors_adapter.cpp`](../src/bridge/cosserat_priors_adapter.cpp)
- Combined driver: [`src/examples/cosserat_estimator_driver.cpp`](../src/examples/cosserat_estimator_driver.cpp)
- Validation tests: [`src/tests/test_cosserat_adapter.cpp`](../src/tests/test_cosserat_adapter.cpp)
- Scenario config: [`config/6_cosserat_priors.yaml`](../config/6_cosserat_priors.yaml)
