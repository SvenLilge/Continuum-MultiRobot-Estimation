# Continuum-MultiRobot-Estimation

This repository is part of the following publication:

State Estimation for Continuum Multi-Robot Systems on SE(3)\
Sven Lilge, Timothy D. Barfoot, Jessica Burgner-Kahrs\
IEEE Transactions on Robotics 2024 (To appear)

Detailed instructions on how to use this code will be added shortly!

### Dependencies (C++)

The C++ implementation makes use of both the Eigen library for matrix arithmetic and the Visualization Toolkit for rendering, which need to be installed in order to compile the provided C++ code.

- [Eigen Library](http://eigen.tuxfamily.org/index.php?title=Main_Page)
- [Visualization Toolkit](https://vtk.org/)
 
### Installation Instructions using CMake (C++)

In the root directory of the C++ code run the following commands:

	mkdir build
	cd build
	cmake ..
	make

The code should compile in Release mode by default (for performance), otherwise this can be forced with:
	
	mkdir build
	cd build
	cmake -DCMAKE_BUILD_TYPE=Release ..
	make

Alternatively, the code can be compiled in Debug mode, enabling useful asserts for debugging:
	
	mkdir build
	cd build
	cmake -DCMAKE_BUILD_TYPE=Debug ..
	make


Afterwards, you can find and execute compiled examples in the 'examples' folder.

### References

If you found the provided continuum robot state estimation implementation helpful or used parts of it yourself, please refer to it using the following BibTeX entries to cite our work:

[1] State Estimation for Continuum Multi-Robot Systems on SE(3)
	 
	@article{Lilge2024,
		title={State Estimation for Continuum Multi-Robot Systems on SE (3)},
		author={Lilge, Sven and Barfoot, Timothy D and Burgner-Kahrs, Jessica},
		journal={arXiv preprint arXiv:2401.13540},
		year={2024}
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
