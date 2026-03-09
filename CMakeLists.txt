cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(continuum-multirobot-estimation)

#Important CMAKE stuff
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

#VTK
find_package(VTK REQUIRED)

#Eigen
find_package (Eigen3 REQUIRED NO_MODULE)

#yaml-cpp
find_package(yaml-cpp REQUIRED)

#include directories
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR})

#include source code
file(GLOB H_FILES ${PROJECT_SOURCE_DIR}/include/*.h)
file(GLOB CPP_FILES ${PROJECT_SOURCE_DIR}/src/*.cpp)

# Define the executable
add_executable(1_continuum_robot src/examples/1_continuum_robot.cpp ${H_FILES} ${CPP_FILES})
add_executable(2_parallel_continuum_robot src/examples/2_parallel_continuum_robot.cpp ${H_FILES} ${CPP_FILES})
add_executable(3_continuous_stewart_gough src/examples/3_continuous_stewart_gough.cpp ${H_FILES} ${CPP_FILES})
add_executable(4_collaborative_continuum_robots src/examples/4_collaborative_continuum_robots.cpp ${H_FILES} ${CPP_FILES})
add_executable(5_fbg_measurements src/examples/5_fbg_measurements.cpp ${H_FILES} ${CPP_FILES})

# Link the necessary libraries
target_link_libraries(1_continuum_robot ${VTK_LIBRARIES} Eigen3::Eigen yaml-cpp::yaml-cpp)
target_link_libraries(2_parallel_continuum_robot ${VTK_LIBRARIES} Eigen3::Eigen yaml-cpp::yaml-cpp)
target_link_libraries(3_continuous_stewart_gough ${VTK_LIBRARIES} Eigen3::Eigen yaml-cpp::yaml-cpp)
target_link_libraries(4_collaborative_continuum_robots ${VTK_LIBRARIES} Eigen3::Eigen yaml-cpp::yaml-cpp)
target_link_libraries(5_fbg_measurements ${VTK_LIBRARIES} Eigen3::Eigen yaml-cpp::yaml-cpp)

# Set the output directory for the executable
set_target_properties(1_continuum_robot PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../examples)
set_target_properties(2_parallel_continuum_robot PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../examples)
set_target_properties(3_continuous_stewart_gough PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../examples)
set_target_properties(4_collaborative_continuum_robots PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../examples)
set_target_properties(5_fbg_measurements PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../examples)

# ── Automated tests (no VTK / no GUI) ────────────────────────────────────────

# T1: ConfigLoader unit tests — only needs yaml-cpp + Eigen
add_executable(test_config_loader
    src/tests/test_config_loader.cpp
    ${H_FILES}
    src/config_loader.cpp
    src/utilities.cpp
)
target_link_libraries(test_config_loader Eigen3::Eigen yaml-cpp::yaml-cpp)
set_target_properties(test_config_loader PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../examples)

# T2: Estimation integration tests — needs estimator but NOT VTK
add_executable(test_estimation
    src/tests/test_estimation.cpp
    ${H_FILES}
    src/continuum_robot_state_estimator.cpp
    src/config_loader.cpp
    src/utilities.cpp
)
target_link_libraries(test_estimation Eigen3::Eigen yaml-cpp::yaml-cpp)
set_target_properties(test_estimation PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../examples)
