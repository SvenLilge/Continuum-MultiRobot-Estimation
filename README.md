# Continuum-MultiRobot-Estimation

This repository is part of the following publication:

> **State Estimation for Continuum Multi-Robot Systems on SE(3)**
> Sven Lilge, Timothy D. Barfoot, Jessica Burgner-Kahrs
> *IEEE Transactions on Robotics 2024*

A detailed documentation for the code [can be found here](https://github.com/SvenLilge/Continuum-MultiRobot-Estimation/wiki/Documentation).

---

## Dependencies

| Library | Purpose |
|---|---|
| [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) | Matrix arithmetic |
| [VTK](https://vtk.org/) | 3D rendering |
| [yaml-cpp](https://github.com/jbeder/yaml-cpp) | Configuration file parsing |

**macOS** (Homebrew):

```bash
brew install vtk eigen yaml-cpp
```

**Ubuntu/Debian**:

```bash
sudo apt install libeigen3-dev libvtk9-dev libyaml-cpp-dev
```

---

## Building

```bash
mkdir build && cd build
cmake ..
cmake --build .
```

The code compiles in Debug mode by default (useful asserts and error messages). For faster execution:

```bash
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build .
```

Executables are placed in `examples/`.

---

## Running the Examples

A single viewer executable takes any YAML configuration file as argument:

```bash
./examples/continuum_robot_viewer config/1_continuum_robot.yaml
```

### Available configurations

| Config | Description |
|---|---|
| `1_continuum_robot.yaml` | Single robot with strain measurements |
| `2_parallel_continuum_robot.yaml` | Two parallel robots with common end-effector |
| `3_continuous_stewart_gough.yaml` | Six-legged continuous Stewart-Gough platform |
| `4_collaborative_continuum_robots.yaml` | Three collaborating robots with coupling constraints |
| `5_fbg_measurements.yaml` | Two robots with FBG sensors (data from CSV) |

Config `6_cosserat_priors.yaml` is **not** runnable through the viewer — its measurements, control inputs, and initial guess are populated at runtime by `cosserat_estimator_driver`. See [Running the combined driver](#running-the-combined-driver) below.

### Interactive controls

The viewer opens a VTK window where you can adjust strain and force inputs in real time:

| Key | Action |
|---|---|
| `Up` / `Down` | Select component |
| `Left` / `Right` | Adjust value |
| `Tab` | Switch between strain and force columns |
| `[` / `]` | Halve / double step size |
| `0` | Zero selected component |
| `r` | Reset all to zero |

---

## Optional: Combined Build with Cosserat Rod Model

The estimator can be driven by **physics-based priors** from the `CosseratRodModel` in the [tdcr-modeling](https://github.com/SvenLilge/tdcr-modeling) repo — strains, strain derivatives, and disk-frame poses predicted from tendon tensions. This lets the estimator infer the rod's full deformation field without needing physical sensors on every node. The integration is **fully optional** and opt-in; the default standalone build is unchanged (no GSL, no tdcr checkout required).

### What the combined pipeline does

Given tendon tensions, the Cosserat model solves the quasi-static rod equilibrium and produces strains `(v, u)`, their s-derivatives `(v̇, u̇)`, and disk frames along the backbone. A thin **adapter** (in `src/bridge/`) converts these into three estimator inputs: (1) per-node strain *pseudo-measurements*, (2) per-segment GP-prior *control inputs*, and (3) a *custom initial guess* built from the disk frames. The estimator then runs its usual MAP optimization — now seeded with physics. 

### Directory layout

Both repos must be **sibling directories** (no submodule, no special setup):

```
your-workspace/
  tdcr-modeling/                     # clone of tdcr-modeling
  Continuum-MultiRobot-Estimation/   # clone of this repo
```

### Additional dependency

The combined build requires [GSL](https://www.gnu.org/software/gsl/) (used by the Cosserat model's ODE solver — not needed by the standalone build):

```bash
# macOS
brew install gsl

# Ubuntu/Debian
sudo apt install libgsl-dev
```

### Building the combined pipeline

Ninja is strongly preferred — the Unix Makefiles generator has been observed to occasionally fail creating intermediate output directories for sources under `src/examples/` and `src/bridge/`.

```bash
cd Continuum-MultiRobot-Estimation
mkdir build && cd build
cmake -G Ninja -DCMAKE_BUILD_TYPE=Release -DUSE_LOCAL_TDCR=ON ..
cmake --build .
```

If the two repos are **not** in sibling directories, point CMake at the tdcr path explicitly:

```bash
cmake -G Ninja -DCMAKE_BUILD_TYPE=Release -DUSE_LOCAL_TDCR=ON \
      -DTDCR_ROOT=/path/to/tdcr-modeling/c++ ..
```

The combined build produces two **additional** executables on top of the standalone targets:

| Executable | Role |
|---|---|
| `cosserat_estimator_driver` | End-to-end demo: Cosserat FK → adapter → state estimator → per-node strain comparison table. Optional VTK visualization. |
| `test_cosserat_adapter` | Four validation tests (frame permutation, arclength resampling, end-to-end convergence, prior-vs-no-prior A/B). |

### Running the combined driver

Run from the repo root (not from `build/`) — executables are placed in `examples/`:

```bash
./examples/cosserat_estimator_driver config/6_cosserat_priors.yaml
```

You will see:
1. Cosserat FK convergence (residual ≈ 1e-19).
2. Adapter summary (how many measurements / control inputs / initial-guess nodes it produced).
3. Estimator Newton iterations and final cost.
4. A per-node strain-comparison table: **Cosserat prior vs Estimator result** for all 6 strain components (ν₁–ν₃, ω₁–ω₃).
5. `max |diff|` per component and overall — expect ~5e-3 rad/m on the worst node (junction) and ~1e-8 everywhere else.

**Flags** (both optional):

| Flag | Effect |
|---|---|
| `--visualize` | After the estimator run, open a VTK window with the backbone, coordinate frames, and covariance ellipsoids. Mouse drag rotates, scroll zooms, close to exit. |
| `--no-control-inputs` | Skip the adapter's GP-prior control inputs. Useful for A/B comparing "measurements + initial guess only" against the full-prior run. |

Example: estimator run with visualization enabled:

```bash
./examples/cosserat_estimator_driver config/6_cosserat_priors.yaml --visualize
```

### Running the validation tests

```bash
./examples/test_cosserat_adapter                             # defaults to config/6_cosserat_priors.yaml
./examples/test_cosserat_adapter config/6_cosserat_priors.yaml
```

Expected output: `=== Results: 4 passed, 0 failed ===`, including the prior-vs-no-prior line:

```
mean |diff|  no-prior = 0.223755
mean |diff|  prior    = 0.000233825
```

### How it works (architecture)

The `ContinuumRodPriors` DTO (Data Transfer Object — pure Eigen + STL, header [`include/continuum_rod_priors.h`](include/continuum_rod_priors.h)) is the **data boundary** between the two repos — no estimator-side header `#include`s anything from tdcr-modeling. The Cosserat-specific adapter lives in `src/bridge/` so it is compiled **only** when `USE_LOCAL_TDCR=ON`. Any future rod model (PCC, analytical, sub-segment Cosserat, …) can produce a `ContinuumRodPriors` with no estimator-side changes.

When `USE_LOCAL_TDCR=OFF` (the default), none of the above applies — no GSL, no tdcr checkout, no extra targets, and the standalone build is byte-identical to the pre-integration behavior.

---

## Running the Tests

```bash
./examples/test_config_loader        # ConfigLoader unit tests
./examples/test_estimation           # Estimation integration tests
```

Both accept an optional argument to override the project root path (default: `..`):

```bash
./examples/test_config_loader /path/to/project
```

The combined build (`USE_LOCAL_TDCR=ON`) adds a third test executable:

```bash
./examples/test_cosserat_adapter     # Cosserat adapter + integration tests (4 checks)
```

---

## Configuration File Structure

All parameters live in `config/`. Edit a YAML file and re-run the viewer -- no recompilation needed.

### `topology` -- robot geometry and structure

- `N` -- number of robots
- `K` -- estimation nodes per robot
- `M` -- interpolation nodes between estimation nodes
- `L` -- robot lengths (meters)
- `Ti0` -- base frame per robot (`identity`, `translation`, `matrix`, `csv_file`)
- `robot_coupling` -- coupling constraints between robots / end-effector
- `lock_first_pose`, `lock_last_pose`, `lock_first_strain`, `lock_last_strain`

### `hyperparameters` -- covariance tuning

- `noise_std` -- measurement noise standard deviations (`R_p`, `R_o`, `R_v`, `R_u`, `R_fbg`)
- `R_pose_scale`, `R_strain_scale`, etc. -- scale factors for covariance matrices
- `Qc_diagonal` -- process noise diagonal (controls prior stiffness)

### `options` -- solver settings

- `solver`: `Newton` or `NewtonLineSearch`
- `initial_guess`: `Straight`, `Last`, or `Custom`
- `max_iterations`, `convergence_threshold`, `kirchhoff_rods`

### `measurements` -- sensor inputs

- `type`: `Strain`, `Pose`, or `FBGStrain`
- Single node, node range (`idx_node_range`), or CSV (`source: csv_file`)

### `control_inputs` -- external inputs

- `type`: `None`, `Constant`, or `PiecewiseLinear`
- 12-element vectors: 6 strain components + 6 force components

### `visualization` -- rendering settings

- `window_width`, `window_height`, `render_frames`, `render_covariance`, `covariance_n_std`, `verbose`

---

## References

If you found the provided continuum robot state estimation implementation helpful or used parts of it yourself, please refer to it using the following BibTeX entries to cite our work:

[1] State Estimation for Continuum Multi-Robot Systems on SE(3)

	@article{Lilge2024,
		  author={Lilge, Sven and Barfoot, Timothy D. and Burgner-Kahrs, Jessica},
		  journal={IEEE Transactions on Robotics},
  	 	  title={State Estimation for Continuum Multi-Robot Systems on SE(3)},
		  year={2024},
	   	  volume={},
  		  number={},
  		  pages={1-20}
	}

[2] Continuum Robot State Estimation using Gaussian Process Regression on SE(3)

	@article{Lilge2022,
		title={Continuum Robot State Estimation using Gaussian Process Regression on SE (3)},
		author={Lilge, Sven and Barfoot, Timothy D and Burgner-Kahrs, Jessica},
		journal={The International Journal of Robotics Research},
		volume={41},
		number={13-14},
		pages={1099--1120},
		year={2022},
		publisher={SAGE Publications Sage UK: London, England}
	}
