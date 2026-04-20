# The Bridge Adapter — How Cosserat Outputs Become Estimator Inputs

**Audience:** anyone who wants to understand what the integration code actually does at runtime — what data flows where, what gets transformed, and how the estimator's MAP problem ends up with Cosserat physics baked in. Read this alongside [`changes_and_new_files.md`](changes_and_new_files.md) (the *files*) and [`cosserat_integration_results.md`](cosserat_integration_results.md) (the *proof*).

---

## 1. What the bridge is, in one sentence

The bridge is a pair of pure-function adapters that convert the auxiliary outputs of a Cosserat rod model into the three public-API objects that `ContinuumRobotStateEstimator::computeStateEstimate()` consumes — **without modifying the estimator** and **without leaking Cosserat types into any public estimator header**.

```
 CosseratRodModel
       │
       │ (aux outputs via getters)
       ▼
 ┌──────────────────────────────────┐
 │ priorsFromCosseratModel()        │   file: src/bridge/cosserat_priors_adapter.cpp
 └──────────────────────────────────┘
       │
       │ ContinuumRodPriors DTO (Data Transfer Object)
       ▼
 ┌──────────────────────────────────┐
 │ cosseratPriorsToEstimator()      │   file: src/bridge/cosserat_priors_adapter.cpp
 └──────────────────────────────────┘
       │
       │ EstimatorPriors { measurements, control_inputs, initial_guess }
       ▼
 ContinuumRobotStateEstimator::computeStateEstimate(...)
```

Both functions live in [`src/bridge/cosserat_priors_adapter.cpp`](../src/bridge/cosserat_priors_adapter.cpp), declared in [`include/cosserat_priors_adapter.h`](../include/cosserat_priors_adapter.h). The directory name `bridge/` is deliberate: it is excluded from `file(GLOB src/*.cpp)`, so these files are compiled **only** when `USE_LOCAL_TDCR=ON`.

---

## 2. Data flow, end-to-end

Using the default driver inputs as a concrete example:

| Stage | Object | Shape / cardinality | Source |
|---|---|---|---|
| 1 | Tendon tensions `q` | 6 × 1 | driver hard-codes `[0.5, 0.2, 0, 0.3, 0, 0.1]` |
| 2 | Cosserat FK | `diskFrames (4, 4·21)` + internal aux state | `CosseratRodModel::forwardKinematics(q, …)` |
| 3 | `ContinuumRodPriors` (DTO) | `s (21,)`, `v, u (21,3)`, `v̇, u̇ (21,3)`, `n, m, f_dist, l_dist (21,3)`, 2 discrete loads | `priorsFromCosseratModel(model, L1, L2)` |
| 4 | `EstimatorPriors` | 11 strain measurements, 10 control inputs, 1 initial-guess SystemState | `cosseratPriorsToEstimator(priors, topology, 0, diskFrames)` |
| 5 | MAP state estimate | `SystemState` (K=11 nodes) | `estimator.computeStateEstimate(...)` |

The bridge handles stages **3 and 4**. Stage 5 is the existing estimator, unchanged.

---

## 3. Bridge function #1 — `priorsFromCosseratModel`

**Signature** ([`include/cosserat_priors_adapter.h:23`](../include/cosserat_priors_adapter.h#L23))

```cpp
ContinuumRodPriors priorsFromCosseratModel(const CosseratRodModel& model,
                                           double L1, double L2);
```

**Role:** copy the auxiliary output getters of a converged Cosserat FK solve into a model-agnostic Eigen struct. Nothing physical happens here — it's a data marshalling step.

### 3.1 Inputs

| Parameter | Meaning | Preconditions |
|---|---|---|
| `model` | A `CosseratRodModel` on which `forwardKinematics()` has **already converged** and whose `hasAuxOutputs()` is `true`. | Throws `std::runtime_error` (from the model-side getters) if not converged. |
| `L1`, `L2` | Physical lengths of the two rod segments (meters). Needed to populate the `s_discrete` arclengths for the junction and tip concentrated loads. | `L1, L2 > 0`. |

### 3.2 What it does, step by step

Lines [`src/bridge/cosserat_priors_adapter.cpp:8-31`](../src/bridge/cosserat_priors_adapter.cpp#L8-L31):

```cpp
ContinuumRodPriors p;

p.s          = model.getArclengthSamples();    // (N,)   N = 21 for the default 2-seg rod
p.v          = model.getStrainV();             // (N,3)  body-frame linear strain
p.u          = model.getStrainU();             // (N,3)  body-frame angular strain
p.v_dot      = model.getStrainVDot();          // (N,3)  ∂v/∂s
p.u_dot      = model.getStrainUDot();          // (N,3)  ∂u/∂s
p.n_internal = model.getInternalForce();       // (N,3)  world-frame resultant force
p.m_internal = model.getInternalMoment();      // (N,3)  world-frame resultant moment
p.f_dist     = model.getDistributedForce();    // (N,3)  distributed load density
p.l_dist     = model.getDistributedMoment();   // (N,3)  distributed moment density

Eigen::Vector3d F_junction, L_junction, F_tip, L_tip;
model.getDiscreteLoads(F_junction, L_junction, F_tip, L_tip);
p.s_discrete = { L1, L1 + L2 };
p.F_discrete = { F_junction, F_tip };
p.L_discrete = { L_junction, L_tip };

p.validate();        // throws std::invalid_argument on shape mismatch
return p;
```

Each getter returns `const Eigen::MatrixXd&`; the assignment `=` does a deep copy, so the DTO owns its data and the model can be destroyed / reused afterwards.

### 3.3 What does **not** happen here

- **No frame conversion.** The DTO still holds Cosserat z-forward body-frame data.
- **No resampling.** The DTO still lives on the Cosserat 21-sample arclength grid.
- **No sign flips.** Values are passed through verbatim.

These conversions happen in the second bridge function (§4) so the DTO stays model-agnostic. A PCC or analytical adapter producing the same DTO will feed the same downstream code untouched.

### 3.4 `ContinuumRodPriors::validate()`

Defined in [`src/continuum_rod_priors.cpp`](../src/continuum_rod_priors.cpp). Called at every adapter boundary. Checks:

- `s` is size `N ≥ 2` and monotonic non-decreasing.
- Every populated `(N,3)` matrix has `N` rows and 3 columns.
- `s_discrete`, `F_discrete`, `L_discrete` have the same length, and every `s_discrete` is inside `[s[0], s[N-1]]`.

Any violation throws `std::invalid_argument` with a message naming the offending field. Optional fields (derivatives, internal resultants, distributed loads) may be left default-constructed (empty); if empty, no shape check runs on them.

---

## 4. Bridge function #2 — `cosseratPriorsToEstimator`

This is where the real work happens: frame-permutation, arclength resampling, and packing into the three estimator-API structs.

**Signature** ([`include/cosserat_priors_adapter.h:60`](../include/cosserat_priors_adapter.h#L60))

```cpp
EstimatorPriors cosseratPriorsToEstimator(
    const ContinuumRodPriors&                          priors,
    const ContinuumRobotStateEstimator::RobotTopology& topology,
    unsigned int                                       robot_idx,
    const Eigen::MatrixXd&                             diskFrames);
```

**Return type** ([`include/cosserat_priors_adapter.h:28-33`](../include/cosserat_priors_adapter.h#L28-L33)):

```cpp
struct EstimatorPriors {
    std::vector<ContinuumRobotStateEstimator::SensorMeasurement> measurements;
    std::vector<ContinuumRobotStateEstimator::ControlInput>      control_inputs;
    ContinuumRobotStateEstimator::SystemState                    initial_guess;
};
```

### 4.1 Inputs — what gets passed in and why

| Parameter | Meaning |
|---|---|
| `priors` | The DTO from §3. Still in Cosserat body-frame convention, still on the 21-sample grid. |
| `topology` | The estimator's `RobotTopology` — needs `N`, `K[robot_idx]`, `L[robot_idx]`, `Ti0[robot_idx]`. Typically loaded from the YAML config. |
| `robot_idx` | Which robot in a multi-robot topology this prior drives. Must be `< topology.N`. |
| `diskFrames` | The `4 × 4N` matrix of disk-frame homogeneous transforms **returned by `CosseratRodModel::forwardKinematics()`**. Each 4×4 block corresponds 1-to-1 with the row of the same index in `priors.s` / `priors.v` / `priors.u`. Used to seed the initial-guess poses. |

Validation at entry ([`src/bridge/cosserat_priors_adapter.cpp:108-127`](../src/bridge/cosserat_priors_adapter.cpp#L108-L127)):

- `priors.validate()` re-runs.
- `robot_idx < topology.N`.
- `priors.s.size() > 0`.
- `diskFrames` shape is `4 × (4·N)` with `N = priors.s.size()`.
- `topology.K[robot_idx] ≥ 2`.

Every violation throws `std::invalid_argument` before any computation starts.

### 4.2 The key transformations

The rest of this section walks through the algorithm. The file-local helpers it uses — `rotConv()`, `permuteStrain()`, `interpolateRowAt()`, `nearestDiskFrame()` — are defined in the anonymous namespace of [`src/bridge/cosserat_priors_adapter.cpp:40-98`](../src/bridge/cosserat_priors_adapter.cpp#L40-L98).

#### 4.2.1 Body-frame permutation

The Cosserat code uses a z-forward body convention (straight rod has strain `v* = [0,0,1]`). The estimator uses x-forward (straight rod has strain `ν = [1,0,0]`). The mapping is a right-handed cyclic permutation:

```
R_conv = [ 0  0  1 ]         z_cosserat → x_estimator
         [ 1  0  0 ]         x_cosserat → y_estimator
         [ 0  1  0 ]         y_cosserat → z_estimator
```

Three applications:

1. **Strain vectors** — use `permuteStrain(v, u)` which elementwise maps `[v0, v1, v2, u0, u1, u2] → [v2, v0, v1, u2, u0, u1]`. This is the same permutation as `R_conv @ v`.
2. **Disk-frame positions** — `p_est = R_conv * p_cosserat`.
3. **Disk-frame rotations** — `R_est = R_conv * R_cosserat * R_conv^T` (similarity transform, because the rotation is expressed *in* the body frame which itself is being rotated).

Unit-test coverage (6 canonical cases to `1e-12`): [`src/tests/test_cosserat_adapter.cpp`](../src/tests/test_cosserat_adapter.cpp) Test A.

#### 4.2.2 Arclength resampling

The Cosserat grid is irregular and has `N = 21` samples (11 in segment 1 + 10 in segment 2, with a duplicate at the junction). The estimator has `K` uniform nodes at `s_k = (k/(K-1)) · L`. For every `s_k`, we linearly interpolate each of `v, u, v̇, u̇` from the two bracketing Cosserat samples — `interpolateRowAt()` at [`src/bridge/cosserat_priors_adapter.cpp:67-81`](../src/bridge/cosserat_priors_adapter.cpp#L67-L81).

Edge cases:
- `s_k ≤ s[0]` → return `row(0)` (clamp).
- `s_k ≥ s[N-1]` → return `row(N-1)` (clamp).
- `s_k` exactly on a Cosserat sample → return that row.
- `v̇` or `u̇` absent (empty matrix) → return `Zero()`.

**Why K = 11 in the default config**: `Δs = 0.02 m`, so `s_k ∈ {0, 0.02, …, 0.20}`. These coincide exactly with the Cosserat sample arclengths, so interpolation across the discontinuity at `s = 0.10` (the segment junction) never happens — avoiding smoothing of the physical curvature step.

For disk frames, we do **not** linearly interpolate (rotation matrices are not a vector space). Instead, `nearestDiskFrame()` picks the 4×4 block whose sample arclength is closest to `s_k` — [`src/bridge/cosserat_priors_adapter.cpp:86-97`](../src/bridge/cosserat_priors_adapter.cpp#L86-L97). This is crude but fine because the estimator refines poses through Newton anyway; the nearest-neighbor seed just has to be in the attractor of the true optimum.

### 4.3 Output #1 — strain measurements

One per estimator node — `K` measurements total.

```cpp
for k in 0..K-1:
    m.type      = SensorMeasurement::Strain;
    m.value     = permuteStrain(v_at[k], u_at[k]);   // 6x1, user-facing target strain
    m.mask      << 1, 1, 1, 1, 1, 1;                 // all 6 components active
    m.idx_robot = robot_idx;
    m.idx_node  = k;
```

**Sign convention — the critical bit.** Internally the estimator does *two* sign flips that cancel:

1. [`continuum_robot_state_estimator.cpp` l. 869](../src/continuum_robot_state_estimator.cpp#L869) — inside `assembleMeasurementTerms()` on the T_bi convention: `strain_des = -measurement.value; e = strain_des - strain_cur;`
2. [`continuum_robot_state_estimator.cpp` l. 1815](../src/continuum_robot_state_estimator.cpp#L1815) — inside `convertStateMeanBodyInertial()` on the way out: `state.robots[n].estimation_nodes[k].strain = -state.robots[n].estimation_nodes[k].strain;`

Since they cancel, the user-facing strain the estimator *reports* equals `+m.value`. So **the adapter stores `m.value` as the target strain directly, no sign flip** ([`src/bridge/cosserat_priors_adapter.cpp:163-171`](../src/bridge/cosserat_priors_adapter.cpp#L163-L171)). An earlier draft of the adapter negated it and produced a sign-inverted estimator output — that is how the audit turned up the bug.

### 4.4 Output #2 — control inputs

One per inter-node segment — `K−1` control inputs total.

```cpp
for k in 0..K-2:
    ci.type        = ControlInput::Constant;
    ci.idx_robot   = robot_idx;
    ci.idx_segment = k;
    ci.values      = { 12x1 };
    ci.values[0].topRows(6)    = permuteStrain(v_at[k], u_at[k]);      // velocity
    ci.values[0].bottomRows(6) = permuteStrain(v̇_at[k], u̇_at[k]);     // acceleration
```

The 12×1 vector packs a pose-tangent in se(3) (6D velocity) followed by its arclength derivative (6D acceleration). The estimator's GP prior uses `velocity` as the *expected* body-frame pose-tangent per unit arclength; setting it to the Cosserat-predicted strain makes the prior's expected evolution match the physics.

**Sign convention.** The estimator has **one** internal sign flip here ([`continuum_robot_state_estimator.cpp` l. 2416](../src/continuum_robot_state_estimator.cpp#L2416)):
```cpp
v_in = -1 * input.values[0].topRows(6);
a_in = -1 * input.values[0].bottomRows(6);
```
It does *not* cancel with anything downstream. The adapter therefore stores the values **without** manual negation ([`src/bridge/cosserat_priors_adapter.cpp:176-189`](../src/bridge/cosserat_priors_adapter.cpp#L176-L189)) — user-facing convention is "state-like, no sign flip", same as measurements.

### 4.5 Output #3 — initial guess

A full `SystemState` with `topology.N` robots, of which only `robot_idx` has its `estimation_nodes` populated. The other robots' `estimation_nodes` vectors are default-constructed (empty) — see §4.6 for the multi-robot caveat.

```cpp
for k in 0..K-1:
    T_cos           = nearestDiskFrame(s_k, priors.s, diskFrames);        // 4x4, Cosserat body frame
    R_est           = R_conv * T_cos.R * R_conv.T;                        // similarity transform
    p_est           = R_conv * T_cos.p;
    T_est_body      = [R_est, p_est; 0 0 0 1];

    node.arclength  = s_k;
    node.pose       = topology.Ti0[robot_idx] * T_est_body;               // inertial frame T_ib
    node.strain     = permuteStrain(v_at[k], u_at[k]);
```

The estimator expects `pose` in the inertial frame `T_ib` ([`continuum_robot_state_estimator.h:87`](../include/continuum_robot_state_estimator.h#L87)); left-multiplying by `Ti0` handles that conversion. `strain` fills the `Node::strain` field in the estimator's body-frame-consumption convention (same as measurements — no sign flip).

Lines: [`src/bridge/cosserat_priors_adapter.cpp:194-208`](../src/bridge/cosserat_priors_adapter.cpp#L194-L208).

### 4.6 Multi-robot caveat

If `topology.N > 1`, the returned `initial_guess.robots` has `N` entries but **only `robot_idx` is populated**. Assigning this directly to `options.custom_guess` will crash later in the estimator — `convertStateMeanBodyInertial()` indexes into every robot's `estimation_nodes`.

For multi-robot setups, call the adapter once per Cosserat-driven robot and merge the `robots[]` vectors yourself before assigning to `options.custom_guess`. Documented in the header at [`include/cosserat_priors_adapter.h:52-57`](../include/cosserat_priors_adapter.h#L52-L57).

---

## 5. How each output is consumed by `computeStateEstimate()`

The adapter produces three things, which all enter the estimator through its existing public API. Nothing about the estimator was modified for this integration.

### 5.1 Strain measurements → measurement factors in the MAP problem

`computeStateEstimate(state, cost, measurements, control_inputs, verbose)` is called with the adapter's `measurements` list. Internally:

1. Measurements are validated (`validateMeasurements`).
2. `assembleMeasurementTerms()` ([`continuum_robot_state_estimator.h:266`](../include/continuum_robot_state_estimator.h#L266)) iterates over each measurement. For `Strain` type, at l. 869 it converts `m.value` into the T_bi convention (`strain_des = -m.value`), forms the error `e = strain_des - strain_cur`, and contributes a cost term `e.T · R_strain⁻¹ · e` to the MAP objective — plus the corresponding linearized block in the sparse Hessian.
3. `R_strain` is the 6×6 covariance loaded from the YAML hyperparameters (`R_v, R_u`, plus the `R_strain_scale` multiplier). Smaller `R_v, R_u` = the estimator trusts the Cosserat pseudo-measurements more.
4. On the way out, `convertStateMeanBodyInertial()` negates strain again at l. 1815, so what comes back in `state.robots[i].estimation_nodes[k].strain` matches the user-facing `m.value` convention.

**Effect on optimization**: 11 dense strain measurements overdetermine the 11-node strain field, so the optimizer essentially reproduces the Cosserat prediction at every node. The only residual is ~5e-3 rad/m at nodes straddling the junction (smoothness vs discontinuity trade-off).

### 5.2 Control inputs → GP-prior transition function

`assemblePriorTerms()` ([`continuum_robot_state_estimator.h:264`](../include/continuum_robot_state_estimator.h#L264)) iterates over each inter-node segment and computes the GP prior factor:

```
e_prior[k] = Φ(Δs, u_k) · x_k − x_{k+1}       // with control input u_k
cost       += e_prior.T · Q⁻¹ · e_prior
```

where `Φ` is the transition function returned by `getTransitionFunction(ControlInput, delta_s)` ([`continuum_robot_state_estimator.h:256`](../include/continuum_robot_state_estimator.h#L256), impl at l. 2416). When `input.type == None`, this is the standard WNOA (white-noise-on-acceleration) transition — a straight-line mean. When `input.type == Constant`, the transition function bakes the 6D velocity (and optionally 6D acceleration) into the prior's *mean* — so the GP's expected pose-tangent over that segment becomes the Cosserat-predicted strain, not zero.

The negation on l. 2416 (`v_in = -input.values[0].topRows(6)`) is the one that motivates the adapter to store values without manual sign flip.

**Effect on optimization**: with dense strain measurements, control inputs add redundancy not information — the final strain field is the same whether `--no-control-inputs` is on or off. Their payoff is in sensor-sparse regimes: with only a few measurements, control inputs propagate the Cosserat-informed mean between them instead of the straight-rod default.

### 5.3 Custom initial guess → starting point for Newton

Inside `computeStateEstimate`, `constructInitialGuess(Options::Custom)` simply returns `m_options.custom_guess` ([`continuum_robot_state_estimator.h:225-227`](../include/continuum_robot_state_estimator.h#L225-L227)). The rest of the optimization takes it as the starting `SystemState` and iterates from there.

Because the adapter's `initial_guess` is already at (or very near) the Cosserat minimum:
- With control inputs on: Newton converges in 2 iterations.
- With control inputs off: Newton converges in 1 iteration.

Compare to `initial_guess: Straight` with no measurements — Newton has nothing to push against, so it returns the straight rod unchanged with a mean per-node strain error of 0.224 rad/m.

### 5.4 Summary table

| Adapter output | Enters estimator via | Consumed in | Contributes to |
|---|---|---|---|
| `measurements` (11 × `SensorMeasurement::Strain`) | `computeStateEstimate(..., measurements, ...)` | `assembleMeasurementTerms()` (l. 869) | Measurement factors in the MAP cost |
| `control_inputs` (10 × `ControlInput::Constant`) | `computeStateEstimate(..., inputs)` | `getTransitionFunction()` (l. 2416), `assemblePriorTerms()` | GP prior mean on each segment |
| `initial_guess` (`SystemState`) | `options.custom_guess`, then `Options::Custom` | `constructInitialGuess(Custom)` | Newton start point |

---

## 6. Sign and frame conventions — cheat sheet

| Concern | Cosserat convention | Estimator external convention (what adapter writes) | Estimator internal convention | Resolution |
|---|---|---|---|---|
| Body axis | z-forward (`v* = [0,0,1]`) | x-forward (`ν = [1,0,0]`) | x-forward | Cyclic permutation `R_conv` (§4.2.1) |
| Strain measurement value | body-frame strain | body-frame strain (**no sign flip**) | `-m.value` internally, negated again on output | Both negations cancel — adapter stores target strain directly |
| Control input value | — | body-frame `[v; v̇]` (**no sign flip**) | `-input.values[0]` internally | Adapter stores without negation; estimator's single internal flip is absorbed |
| Pose | body-frame `T_cos` | inertial-frame `T_ib = Ti0 · R_conv·T_cos·R_conv^T` | `T_ib` / `T_bi` internally | Handled in initial-guess construction (§4.5) |

**Rule of thumb:** the adapter writes values in the estimator's **external user-facing convention** — the same convention the existing YAML configs use. Any internal sign games inside the estimator are its own business and cancel out by the time data crosses the API boundary in either direction.

---

## 7. What the estimator does **not** know

- It never sees `CosseratRodModel` — the type isn't even declared in any estimator-side public header.
- It never sees the Cosserat 21-sample grid — the adapter resamples to the estimator's K uniform nodes before data enters the API.
- It never sees the z-forward body convention — the adapter permutes before the DTO's values reach `SensorMeasurement` / `ControlInput` / `custom_guess`.
- It never sees distributed loads, internal resultants, or discrete tendon loads. Those are *in* the DTO (populated by `priorsFromCosseratModel`) but the current adapter does not project them onto any estimator input. They are there for a future extension that could turn them into additional cost terms or coupling constraints without touching the DTO schema.

---

## 8. Concrete walkthrough — the default driver, end to end

Inputs fed to the adapter by [`src/examples/cosserat_estimator_driver.cpp`](../src/examples/cosserat_estimator_driver.cpp):

```
q = [0.5, 0.2, 0.0, 0.3, 0.0, 0.1] N       # tendon tensions, two 3-tendon sets
L1 = L2 = 0.10 m                            # segment lengths
topology: N=1, K=11, L=0.20 m, Ti0=Identity # single robot, 11 uniform nodes at Δs=0.02
robot_idx = 0
```

What the bridge produces:

```
priors:
  s        : 21 samples, [0, 0.01, ..., 0.09, 0.10, 0.10, 0.11, ..., 0.20]
  v, u     : (21,3) Cosserat body-frame strains
  v̇, u̇   : (21,3) arclength derivatives
  n, m     : (21,3) internal resultants  (unused downstream)
  f_dist   : (21,3) distributed force    (unused downstream)
  l_dist   : (21,3) distributed moment   (unused downstream)
  s_discrete : [0.10, 0.20]              (unused downstream)
  F_discrete : [F_junction, F_tip]       (unused downstream)
  L_discrete : [L_junction, L_tip]       (unused downstream)

ep (EstimatorPriors):
  measurements   : 11 × Strain, value = permuted Cosserat strain at s_k ∈ {0, 0.02, ..., 0.20}
  control_inputs : 10 × Constant, values = [v, v̇] permuted, one per segment (k=0..9)
  initial_guess  : SystemState { robots[0].estimation_nodes[0..10] = { pose, strain } }
```

At node k = 5 (s = 0.10, exactly on the junction):
- Cosserat raw:        `v = [0, 0, 1], u = [-1.5958, 0.2940, 0]` (in z-forward body frame)
- After permutation:   `strain = [1, 0, 0, 0, -1.5958, 0.2940]` (estimator ν₁, ν₂, ν₃, ω₁, ω₂, ω₃)
- `m.value` stored directly as above; estimator's two internal negations cancel.

After `computeStateEstimate`, the estimator reports back `state.robots[0].estimation_nodes[5].strain = [1, 0, 0, 0, -1.5904, 0.2917]` — within 5e-3 rad/m of the Cosserat prediction. The small residual on ω₂, ω₃ is the junction-smoothing artifact discussed in §7 of [`cosserat_integration_results.md`](cosserat_integration_results.md).

---

## 9. Error handling

All failure modes are loud and early:

| Failure | Where it throws | Exception type |
|---|---|---|
| Cosserat FK did not converge | `priorsFromCosseratModel` (getters internally) | `std::runtime_error` |
| DTO shape mismatch | `priors.validate()` (called in both bridge functions) | `std::invalid_argument` |
| `robot_idx ≥ topology.N` | `cosseratPriorsToEstimator` entry check | `std::invalid_argument` |
| `diskFrames` shape wrong | `cosseratPriorsToEstimator` entry check | `std::invalid_argument` |
| `K[robot_idx] < 2` | `cosseratPriorsToEstimator` entry check | `std::invalid_argument` |
| Estimator rejects measurements / options | `computeStateEstimate` internal validators | `std::runtime_error` |

The driver does not catch any of these — a malformed pipeline aborts the process immediately, which is what you want during development. For production use, wrap the two bridge calls in a try/catch and log the message.

---

## 10. Further reading

- [`../../tdcr-modeling/doc/integration_design_for_supervisor.md`](../../tdcr-modeling/doc/integration_design_for_supervisor.md) — the cross-repo design proposal (rationale, constraints, alternatives considered).
- [`cosserat_integration_results.md`](cosserat_integration_results.md) — measured improvement (957× error reduction) and experimental protocol.
- [`changes_and_new_files.md`](changes_and_new_files.md) — file-by-file summary of the pipeline phase.
- Sibling repo: [`../../tdcr-modeling/doc/paper_equation_map.md`](../../tdcr-modeling) — Rucker 2011 sign audit of the Cosserat side.
