# Continuum-MultiRobot-Estimation

This repository is part of the following publication:

State Estimation for Continuum Multi-Robot Systems on SE(3)\
Sven Lilge, Timothy D. Barfoot, Jessica Burgner-Kahrs\
IEEE Transactions on Robotics 2024

**A detailed documentation for the code [can be found here](https://github.com/SvenLilge/Continuum-MultiRobot-Estimation/wiki/Documentation)!**

### Dependencies (C++)

The C++ implementation requires the following libraries:

- [Eigen Library](http://eigen.tuxfamily.org/index.php?title=Main_Page) — matrix arithmetic
- [Visualization Toolkit (VTK)](https://vtk.org/) — 3D rendering
- [yaml-cpp](https://github.com/jbeder/yaml-cpp) — runtime configuration file parsing

On macOS with Homebrew:

	brew install vtk eigen yaml-cpp

On Ubuntu/Debian:

	sudo apt install libeigen3-dev libvtk9-dev libyaml-cpp-dev

### Building the Code (CMake)

In the root directory of the repository:

	mkdir build
	cd build
	cmake ..
	cmake --build .

The code compiles in Debug mode by default, which enables useful asserts (input validation, descriptive error messages). For faster execution, compile in Release mode:

	mkdir build
	cd build
	cmake -DCMAKE_BUILD_TYPE=Release ..
	cmake --build .

The compiled executables are placed in the `examples/` folder.

### Running the Examples

Each example reads its parameters from a YAML configuration file in the `config/` folder. Run from the `examples/` directory:

	cd examples

	./1_continuum_robot
	./2_parallel_continuum_robot
	./3_continuous_stewart_gough
	./4_collaborative_continuum_robots
	./5_fbg_measurements

By default each executable loads its corresponding file from `../config/`. You can pass a different config file as an argument:

	./1_continuum_robot path/to/my_config.yaml

### Modifying Simulation Parameters (No Recompilation)

All simulation parameters live in the `config/` folder. To change a parameter:

1. Open the corresponding YAML file, e.g. `config/1_continuum_robot.yaml`
2. Edit the desired parameter
3. Re-run the executable — no recompilation needed

Example: change the position noise in example 1:

	# config/1_continuum_robot.yaml
	hyperparameters:
	  noise_std:
	    R_p: 0.005   # was 0.002 — increase position noise

Then simply run:

	./examples/1_continuum_robot

### Configuration File Structure

Each YAML config file has four main sections:

**`topology`** — robot geometry and structure
- `N`: number of robots
- `K`: estimation nodes per robot
- `M`: interpolation nodes between estimation nodes (1 = no interpolation)
- `L`: robot lengths in meters
- `Ti0`: base frame of each robot (supports `identity`, `translation`, `matrix`, `csv_file` types)
- `robot_coupling`: coupling constraints between robots or to a common end-effector
- boundary condition flags: `lock_first_pose`, `lock_last_pose`, `lock_first_strain`, `lock_last_strain`

**`hyperparameters`** — probabilistic tuning (covariance matrices)
- `noise_std`: measurement noise standard deviations (`R_p`, `R_o`, `R_v`, `R_u`, `R_fbg`)
- `R_pose_scale`, `R_strain_scale`, etc.: scale factors applied to each covariance matrix
- `Qc_diagonal`: process noise diagonal (controls stiffness of the prior)

**`options`** — solver settings
- `solver`: `Newton` or `NewtonLineSearch`
- `initial_guess`: `Straight`, `Last`, or `Custom`
- `max_iterations`, `convergence_threshold`, `kirchhoff_rods`

**`measurements`** — sensor inputs
- `type`: `Strain`, `Pose`, or `FBGStrain`
- single node, node range (`idx_node_range`), or loaded from CSV (`source: csv_file`)

**`visualization`** — rendering settings
- `window_width`, `window_height`, `render_frames`, `render_covariance`, `covariance_n_std`, `verbose`

### Examples Overview

| Executable | Config file | Description |
|---|---|---|
| `1_continuum_robot` | `config/1_continuum_robot.yaml` | Single robot with strain measurements |
| `2_parallel_continuum_robot` | `config/2_parallel_continuum_robot.yaml` | Two parallel robots with a common end-effector and pose measurement |
| `3_continuous_stewart_gough` | `config/3_continuous_stewart_gough.yaml` | Six-legged continuous Stewart-Gough platform (no measurements) |
| `4_collaborative_continuum_robots` | `config/4_collaborative_continuum_robots.yaml` | Three collaborating robots with coupling constraints |
| `5_fbg_measurements` | `config/5_fbg_measurements.yaml` | Two robots with Fiber Bragg Grating sensors; data loaded from CSV |

### References

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
