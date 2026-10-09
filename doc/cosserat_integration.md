# Feeding a Physics Model into the State Estimator

This document explains how the Cosserat rod model from the sibling repository [`tdcr-modeling`](../../tdcr-modeling) is connected to the state estimator in this repository, and how the connection was checked. It assumes no prior knowledge of state estimation.

Reading order:

1. **This document**: background, how the connection works, sanity checks.
2. [`evaluation.md`](evaluation.md): the experiment comparing three ways of using the model, with results and open issues.
3. The physics of the rod model itself: [`tdcr-modeling/doc/cosserat_rod_model.md`](../../tdcr-modeling/doc/cosserat_rod_model.md).

---

## 1. The big picture

A continuum robot bends smoothly along its whole length. Knowing its exact shape matters for control, but measuring it directly is hard: a typical setup has only one position sensor at the tip.

Two sources of information are available:

- **A sensor** gives real but sparse and noisy data, for example the tip position ±1 mm.
- **A physics model** (the Cosserat rod model) predicts the full shape from the tendon tensions. It is dense but imperfect: it does not know about unexpected forces, such as someone pushing on the robot.

**State estimation** combines the two. It looks for the shape that best agrees with both the sensor data and the model's expectations, weighting each source by how much it is trusted. This work lets the estimator use the Cosserat model's predictions. Before this, the estimator only had a generic assumption that "the rod bends smoothly."

### Key terms

| Term | Meaning |
|---|---|
| Arclength $s$ | Distance along the rod from the base (m). The rod is $L = 0.2$ m long. |
| Node | One of $K$ evenly spaced points along the rod where the estimator stores the shape ($K = 11$ here, every 0.02 m). |
| Pose $T$ | Position and orientation of the cross-section, as a 4×4 matrix. |
| Strain $\varepsilon = [\nu;\ \omega]$ | How the pose changes per unit length, in the rod's own (body) frame. $\nu$ = linear strain (stretch/shear), $\omega$ = angular strain (bending curvature and twist, rad/m). |
| Prior | What the estimator believes about the shape *before* seeing any measurement. |
| Gaussian process (GP) prior | The prior used here: strain is treated as a random curve that varies smoothly along $s$. |
| Measurement | A sensor reading, here a tip position or a strain value, with an assumed noise level. |
| Control input | Known information that shifts the prior, for example "the strain should follow this curve." |
| Noise standard deviation | How uncertain a measurement is. Smaller means more trusted. The estimator weights each error by $1/\sigma^2$. |
| MAP estimate | *Maximum a posteriori*: the single most likely shape given the prior and the measurements. |
| DTO | *Data transfer object*: a plain container of numbers passed between the two codebases. |

## 2. How the estimator works (just enough to follow the rest)

The estimator (Lilge, Barfoot, Burgner-Kahrs, *State Estimation for Continuum Multi-Robot Systems on SE(3)*, T-RO 2024) stores a pose $T_k$ and a strain $\varepsilon_k$ at every node $k$. It finds them by minimizing a weighted sum of squared errors with Newton's method:

$$
J = \tfrac12 \sum_{k} e_{\text{prior},k}^\top\, Q_k^{-1}\, e_{\text{prior},k} \;+\; \tfrac12 \sum_{j} e_{\text{meas},j}^\top\, R_j^{-1}\, e_{\text{meas},j}.
$$

- **Prior terms** $e_{\text{prior},k}$ measure how much the shape between node $k$ and node $k+1$ deviates from what the GP prior expects. $Q_k$ grows with the process-noise setting $Q_c$: a larger $Q_c$ means a more flexible prior.
- **Measurement terms** $e_{\text{meas},j}$ measure how far the estimate is from each measurement. $R_j$ is that measurement's noise covariance.

**The default prior** (no control inputs) is called *white noise on acceleration* (WNOA). It says the strain's derivative is pure random noise. On average the strain therefore stays constant: an arc of constant curvature costs nothing, and every change in curvature costs something. The prior has no idea *which* curvature the rod has; only measurements can tell it.

**Control inputs** (Lilge and Barfoot, *Incorporating Control Inputs in Continuous-Time GP State Estimation*, 2025, eq. 34) give the prior a known trend. Each input is 12 numbers per segment between two nodes:

- the first 6 (the *velocity channel*) are a known strain $\varepsilon_{in}$;
- the last 6 (the *acceleration channel*) are a known strain derivative $a_{in}$.

The prior then becomes: total strain $\varepsilon = \varepsilon_{in} + \varepsilon_b$, where the unknown part $\varepsilon_b$ satisfies $\varepsilon_b' = a_{in} + \text{noise}$. With both channels at zero this reduces to WNOA.

**Conventions.** The estimator's body $x$-axis points along the rod, so a straight rod has $\nu = [1, 0, 0]$ and $\omega = 0$. This config sets `kirchhoff_rods: true`, which fixes $\nu = [1,0,0]$ (no stretch or shear) and leaves only $\omega$ free.

## 3. The bridge: from Cosserat outputs to estimator inputs

```
CosseratRodModel::forwardKinematics(q, f_ext, l_ext)        [tdcr-modeling]
        │  getters (strains, loads, disk frames)
        ▼
priorsFromCosseratModel()  ──►  ContinuumRodPriors (DTO)    [this repo, src/bridge/]
        │
        ▼
cosseratPriorsToEstimator() ──► EstimatorPriors {measurements, control_inputs, initial_guess}
        │
        ▼
ContinuumRobotStateEstimator::computeStateEstimate(...)     [unchanged estimator]
```

Both bridge functions are in [`src/bridge/cosserat_priors_adapter.cpp`](../src/bridge/cosserat_priors_adapter.cpp) and declared in [`include/cosserat_priors_adapter.h`](../include/cosserat_priors_adapter.h). **The estimator itself was not modified.** Everything enters through its existing public interface. The one exception is a missing `#include <cassert>`, added so Release builds compile.

### 3.1 Step 1: copy the model outputs into a plain container

`priorsFromCosseratModel(model, L1, L2, diskFrames)` copies the model's outputs (see the rod-model doc, Section 7) into `ContinuumRodPriors` ([`include/continuum_rod_priors.h`](../include/continuum_rod_priors.h)). That container depends only on Eigen, so any rod model could fill it.

| Field | Content |
|---|---|
| `s` | 21 arclength samples |
| `v`, `u`, `v_dot`, `u_dot` | Strains and their derivatives (Cosserat body frame) |
| `n_internal`, `m_internal`, `f_dist`, `l_dist` | Internal force/moment and tendon loads per unit length (world frame) |
| `s_discrete`, `F_discrete`, `L_discrete` | Point loads where tendons end (junction and tip) |
| `epsilon_in_accel` | $-[K_{se}^{-1} R^\top f_t;\ K_{bt}^{-1} R^\top l_t]$ at each sample: the strain derivative caused by the distributed tendon loads (R&W eq. 5). Used by Method 3 in the evaluation. |
| `epsilon_jump_discrete` | $-[K_{se}^{-1} R^\top F;\ K_{bt}^{-1} R^\top L]$ per point load: the strain jump where tendons end (R&W eq. 20). Used by Method 3. |

`validate()` checks all array sizes and throws `std::invalid_argument` on a mismatch. The model's getters throw if the Cosserat solve did not converge, so bad data never reaches the estimator.

### 3.2 Step 2: convert to what the estimator understands

`cosseratPriorsToEstimator(priors, topology, robot_idx, diskFrames, mode)` does three things.

**(a) Rotate the axis convention.** The Cosserat code points the rod along body $z$; the estimator points it along body $x$. Relabelling the axes with the cyclic permutation

$$
P = \begin{bmatrix} 0 & 0 & 1 \\ 1 & 0 & 0 \\ 0 & 1 & 0 \end{bmatrix}
\qquad (P\,e_z = e_x,\ \ P\,e_x = e_y,\ \ P\,e_y = e_z,\ \ \det P = +1)
$$

converts everything:

- strains: $\nu = P\,v$ and $\omega = P\,u$;
- disk positions: $p_{est} = P\,p$;
- disk orientations: $R_{est} = P\,R\,P^\top$. Both the world axes and the body axes are relabelled, so $P$ appears on both sides.

The result is then placed in the world with the robot's base pose: $T = T_{i0}\, T_{est}$.

**(b) Resample to the estimator's nodes.** The model gives 21 samples (every 0.01 m). The estimator wants $K$ evenly spaced nodes at $s_k = kL/(K-1)$. Strains are linearly interpolated. Poses use the nearest disk, because rotation matrices cannot be averaged linearly; Newton refines them anyway. With $K = 11$ every node lands exactly on a sample, so nothing is ever interpolated across the strain jump at the junction ($s = 0.10$).

**(c) Build three kinds of input.**

| Output | How many | Content |
|---|---|---|
| Strain measurements | $K$ (one per node) | Target for the strain stored in the estimator's state, all 6 components active. Their noise is set by `R_v` (linear part) and `R_u` (angular part) in the YAML. |
| Control inputs | $K-1$ (one per segment) | Depends on the mode (table below). |
| Initial guess | 1 | Cosserat poses, plus the state strain implied by the mode, as the starting point for Newton. |
| `model_strain` | $K$ | The Cosserat strain at each node, kept for comparisons |

**The strain in the estimator's state.** Once control inputs are used, the estimator's state strain is only the *bias*: the part not already carried by the velocity input. The strain that shapes the rod is (state strain) + (velocity input). This follows from the 2025 control-input paper, eqs. 1, 3 and 13; the estimator never adds the input back. Measurements and the initial guess must therefore target the bias.

The `mode` argument picks one consistent combination:

| `ControlInputMode` | Velocity input | Acceleration input | Strain measurements and initial strain | Used by |
|---|---|---|---|---|
| `None` (default) | none | none | Cosserat strain | Run B; evaluation Method 2 |
| `StrainAsInput` | Cosserat strain at the segment midpoint | 0 | 0 (bias) | Run C |
| `ForceAsInput` | 0 | `epsilon_in_accel`, plus each interior strain jump spread over the segment that starts at it | Cosserat strain (bias = total, because the velocity input is 0) | Evaluation Method 3 |

The acceleration input is 0 in `StrainAsInput` because Lilge & Barfoot (2025) define it as independent of the velocity input, not as its derivative. The velocity input is sampled at the segment midpoint because it is constant over the segment. That also picks the correct side of the strain jump at the junction node.

**Signs.** Internally the estimator stores the inverse pose and therefore the negative strain, and converts on the way in and on the way out. As a user you always pass the physical body-frame strain and strain derivative; the adapter never flips signs by hand.

**Multi-robot caveat.** The initial guess fills only robot `robot_idx`. With several robots, call the adapter once per robot and merge the results before passing the guess to the estimator; otherwise the estimator reads empty data for the other robots.

## 4. Software design

Two independent repositories must work together without tangling. The rules, and why each exists:

| Rule | Reason |
|---|---|
| The container (`ContinuumRodPriors`) depends only on Eigen. | The estimator never depends on a specific rod model; a different model only needs a new adapter. |
| The adapter header only *forward-declares* `CosseratRodModel`; only the `.cpp` includes `cosseratrodmodel.h`. | No public estimator header pulls in `tdcr-modeling`. |
| The adapter `.cpp` lives in `src/bridge/`. | The main build globs `src/*.cpp`. Putting the adapter there would break builds that don't have `tdcr-modeling`. |
| The CMake option `USE_LOCAL_TDCR` is `OFF` by default. | The estimator builds exactly as before, with no GSL needed. `ON` runs `add_subdirectory(../tdcr-modeling/c++)` and adds the bridge targets. |
| No copied sources, no git submodule, no `find_package`. | Copies drift apart; submodules add friction; `tdcr-modeling` has no CMake package config. |

Files added on this side:

| File | Role |
|---|---|
| [`include/continuum_rod_priors.h`](../include/continuum_rod_priors.h), [`src/continuum_rod_priors.cpp`](../src/continuum_rod_priors.cpp) | The container and its `validate()` |
| [`include/cosserat_priors_adapter.h`](../include/cosserat_priors_adapter.h), [`src/bridge/cosserat_priors_adapter.cpp`](../src/bridge/cosserat_priors_adapter.cpp) | The two bridge functions |
| [`src/examples/cosserat_estimator_driver.cpp`](../src/examples/cosserat_estimator_driver.cpp) | Demo: runs model + bridge + estimator, prints per-node strain and position comparisons; returns an error if the shape is off by ≥ 2 mm (`--no-control-inputs`, `--visualize`) |
| [`src/tests/test_cosserat_adapter.cpp`](../src/tests/test_cosserat_adapter.cpp) | 6 tests (Section 5) |
| [`config/6_cosserat_priors.yaml`](../config/6_cosserat_priors.yaml) | Scenario for the demo and the tests |
| [`src/examples/evaluation_cosserat_priors.cpp`](../src/examples/evaluation_cosserat_priors.cpp), [`config/7_evaluation_s_shape.yaml`](../config/7_evaluation_s_shape.yaml), [`scripts/evaluate.py`](../scripts/evaluate.py) | The evaluation: runs the sweeps and makes the figures (see [`evaluation.md`](evaluation.md)) |

## 5. Sanity check: is the plumbing correct?

**Scenario** ([`config/6_cosserat_priors.yaml`](../config/6_cosserat_priors.yaml)): one robot, $L = 0.2$ m, $K = 11$, tendon tensions $q = [0.5, 0.2, 0, 0.3, 0, 0.1]$ N, no external load, strain noise `R_v = R_u = 0.01`, base pose fixed. The rod bends by about 1.6 rad/m in segment 1, roughly 9° over that segment.

**Question:** if the estimator is given the Cosserat strain everywhere, does it return the Cosserat **shape**? If frames, signs, resampling or the use of control inputs were wrong, it would not. Comparing node positions is essential. Strain alone is not enough, because with control inputs the state strain is only the bias.

| Run | Mode | Starting cost → final | Strain check | Max position error |
|---|---|---|---|---|
| A | nothing: no measurements, no inputs, straight start | – | mean strain error **0.224 rad/m** (stays straight) | – |
| B | `None`: strain measurements | 140 → 34.7 (2 iterations) | max 5.4 × 10⁻³ rad/m vs Cosserat | **1.1 mm** (at the tip) |
| C | `StrainAsInput`: velocity input | 2.5 × 10⁻¹⁸ (1 iteration) | bias ≈ 10⁻²⁷ (target 0) | **≈ 0** (2 × 10⁻⁹ mm) |

Mean strain error is the average of $|\varepsilon_{est} - \varepsilon_{Cosserat}|$ over all 6 components and 11 nodes. Position errors compare estimated node positions with the Cosserat disk positions.

What this shows:

- **The bridge is correct in both modes.** Run C reproduces the Cosserat shape exactly: each segment's input is exactly the Cosserat strain, which is constant within each Cosserat segment here.
- **Run B is slightly off at the junction.** There the true strain jumps, but run B's strain is a state variable that the GP prior wants to be smooth, so it blurs the jump at nodes 5–6. That gives the 1.1 mm at the tip.
- **What A vs B does not show:** run A has no information at all, so it can only return a straight rod. The 0.224 → 0.00023 rad/m drop (≈ 957×) confirms the information arrives; it says nothing about accuracy against an independent ground truth. That question is answered in [`evaluation.md`](evaluation.md).

**Unit tests** (`./examples/test_cosserat_adapter config/6_cosserat_priors.yaml`, 6/6 pass):

| Test | Checks |
|---|---|
| A. Axis permutation | Six basic cases (straight, bend about x, bend about y, twist, two shears) map correctly, to 1e-12 |
| B. Resampling | A linear strain field on 5 samples is interpolated exactly onto 9 nodes, to 1e-12 |
| C. End-to-end, mode `None` | Strain within 1e-2 of Cosserat (actual 5.4 × 10⁻³ rad/m) **and** positions within 2 mm (actual 1.1 mm) |
| D. With vs without priors | Run B is at least 10× closer to the Cosserat strain than run A (actual ≈ 957×) |
| E. Strain as velocity input | Bias strain stays ≈ 0 **and** positions within 2 mm (actual 2 × 10⁻⁹ mm). Guards against counting the model strain twice (once as input, once as measurements). |
| F. Force as acceleration input | Given only the true base strain, the shape built from the acceleration inputs stays within 2 mm (actual 1.1 mm, from spreading the junction jump over one segment). A sign error on the inputs gives 10.7 mm and fails. |

## 6. Build and run

```bash
# Combined build (Ninja is more reliable than Make for this target on macOS)
cmake -G Ninja -DCMAKE_BUILD_TYPE=Release -DUSE_LOCAL_TDCR=ON -S . -B build
cmake --build build

./examples/test_cosserat_adapter     config/6_cosserat_priors.yaml   # tests A–F
./examples/cosserat_estimator_driver config/6_cosserat_priors.yaml   # run C: strain + position tables
./examples/cosserat_estimator_driver config/6_cosserat_priors.yaml --no-control-inputs   # run B
```

Standalone build (no `tdcr-modeling`, no GSL), which must keep working:

```bash
cmake -DCMAKE_BUILD_TYPE=Release -S . -B build-standalone
cmake --build build-standalone -j
./examples/test_config_loader . && ./examples/test_estimation .   # argument = repo root
```

If CMake fails with *"Could NOT find Boost: Found unsuitable version"* or a program crashes with *"dyld: Library not loaded"*, Homebrew has upgraded VTK's dependencies without rebuilding VTK. Run `brew upgrade vtk`, then reconfigure with `cmake --fresh ...`.
