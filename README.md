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

## Running the Tests

```bash
./examples/test_config_loader        # ConfigLoader unit tests
./examples/test_estimation           # Estimation integration tests
```

Both accept an optional argument to override the project root path (default: `..`):

```bash
./examples/test_config_loader /path/to/project
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
