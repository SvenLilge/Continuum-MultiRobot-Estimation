cmake_minimum_required(VERSION 3.2)
set(CMAKE_CXX_STANDARD 17)

project(continuum-multirobot-estimation)

# Force Release mode unless specified otherwise
if (NOT EXISTS ${CMAKE_BINARY_DIR}/CMakeCache.txt)
  if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
  endif()
endif()

#Important CMAKE stuff
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

#VTK
find_package(VTK REQUIRED)

#Eigen
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

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
target_link_libraries(1_continuum_robot ${VTK_LIBRARIES} Eigen3::Eigen)
target_link_libraries(2_parallel_continuum_robot ${VTK_LIBRARIES} Eigen3::Eigen)
target_link_libraries(3_continuous_stewart_gough ${VTK_LIBRARIES} Eigen3::Eigen)
target_link_libraries(4_collaborative_continuum_robots ${VTK_LIBRARIES} Eigen3::Eigen)
target_link_libraries(5_fbg_measurements ${VTK_LIBRARIES} Eigen3::Eigen)

# Set the output directory for the executable
set_target_properties(1_continuum_robot PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../examples)
set_target_properties(2_parallel_continuum_robot PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../examples)
set_target_properties(3_continuous_stewart_gough PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../examples)
set_target_properties(4_collaborative_continuum_robots PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../examples)
set_target_properties(5_fbg_measurements PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../examples)
