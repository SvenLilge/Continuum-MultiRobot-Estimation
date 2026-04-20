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

# --- Optional: link against a sibling tdcr-modeling checkout -----------------
# Default OFF so the estimator builds standalone (no GSL required). Enable with
# -DUSE_LOCAL_TDCR=ON to build the combined Cosserat-driven pipeline against
# a sibling tdcr-modeling repo checkout.
option(USE_LOCAL_TDCR "Link against sibling tdcr-modeling checkout" OFF)

if (USE_LOCAL_TDCR)
    set(TDCR_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../tdcr-modeling/c++" CACHE PATH
        "Path to the c++ subdir of a tdcr-modeling checkout")
    if (NOT EXISTS "${TDCR_ROOT}/CMakeLists.txt")
        message(FATAL_ERROR
            "USE_LOCAL_TDCR=ON but TDCR_ROOT='${TDCR_ROOT}' is not a tdcr-modeling c++ dir. "
            "Pass -DTDCR_ROOT=<path-to-tdcr-modeling/c++> to override.")
    endif()
    add_subdirectory("${TDCR_ROOT}" "${CMAKE_BINARY_DIR}/tdcr_modeling_build")
endif()

#include directories
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR})

#include source code
file(GLOB H_FILES ${PROJECT_SOURCE_DIR}/include/*.h)
file(GLOB CPP_FILES ${PROJECT_SOURCE_DIR}/src/*.cpp)

# Unified viewer executable (pass any YAML config as argument)
add_executable(continuum_robot_viewer src/examples/continuum_robot_viewer.cpp ${H_FILES} ${CPP_FILES})
target_link_libraries(continuum_robot_viewer ${VTK_LIBRARIES} Eigen3::Eigen yaml-cpp::yaml-cpp)
set_target_properties(continuum_robot_viewer PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../examples)

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

# --- Optional combined driver — USE_LOCAL_TDCR=ON only ----------------------
# Cosserat FK → aux outputs → ContinuumRodPriors DTO (Data Transfer Object).
# End-to-end build smoke
# test of the cross-repo adapter boundary. The adapter .cpp lives in
# src/bridge/ deliberately — outside file(GLOB src/*.cpp) — so it is NEVER
# compiled into the standalone targets above.
if (USE_LOCAL_TDCR)
    # Adapter + end-to-end validation tests (frame permutation, resampling,
    # estimator convergence, prior vs no-prior).  No VTK.
    add_executable(test_cosserat_adapter
        src/tests/test_cosserat_adapter.cpp
        ${H_FILES}
        src/continuum_robot_state_estimator.cpp
        src/continuum_rod_priors.cpp
        src/bridge/cosserat_priors_adapter.cpp
        src/config_loader.cpp
        src/utilities.cpp
    )
    target_link_libraries(test_cosserat_adapter
        PRIVATE Eigen3::Eigen yaml-cpp::yaml-cpp tdcr_modeling)
    set_target_properties(test_cosserat_adapter PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../examples)

    add_executable(cosserat_estimator_driver
        src/examples/cosserat_estimator_driver.cpp
        ${H_FILES}
        src/continuum_robot_state_estimator.cpp
        src/continuum_rod_priors.cpp
        src/bridge/cosserat_priors_adapter.cpp
        src/config_loader.cpp
        src/utilities.cpp
        src/visualizervtk.cpp
    )
    target_link_libraries(cosserat_estimator_driver
        PRIVATE ${VTK_LIBRARIES} Eigen3::Eigen yaml-cpp::yaml-cpp tdcr_modeling)
    set_target_properties(cosserat_estimator_driver PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../examples)
endif()
