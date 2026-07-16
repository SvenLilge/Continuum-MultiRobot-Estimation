// test_config_loader.cpp
// Automated unit tests for ConfigLoader (no VTK, no GUI).
//
// Compile alongside: config_loader.cpp, utilities.cpp
// Link: yaml-cpp, Eigen3
//
// Run from the project root:
//   ./examples/test_config_loader

#include "config_loader.h"
#include "utilities.h"

#include <Eigen/Dense>
#include <Eigen/LU>

#include <cmath>
#include <fstream>
#include <iostream>

// helpers

static int g_pass = 0;
static int g_fail = 0;

#define ASSERT_TRUE(cond, msg)                                  \
  do                                                            \
  {                                                             \
    if (!(cond))                                                \
    {                                                           \
      std::cerr << "  FAIL [" << msg << "]: condition false\n"; \
      ++g_fail;                                                 \
    }                                                           \
    else                                                        \
    {                                                           \
      std::cout << "  PASS [" << msg << "]\n";                  \
      ++g_pass;                                                 \
    }                                                           \
  } while (0)

#define ASSERT_THROWS(expr, msg)                                                \
  do                                                                            \
  {                                                                             \
    bool caught = false;                                                        \
    try                                                                         \
    {                                                                           \
      expr;                                                                     \
    }                                                                           \
    catch (const std::exception &)                                              \
    {                                                                           \
      caught = true;                                                            \
    }                                                                           \
    if (!caught)                                                                \
    {                                                                           \
      std::cerr << "  FAIL [" << msg << "]: expected exception, none thrown\n"; \
      ++g_fail;                                                                 \
    }                                                                           \
    else                                                                        \
    {                                                                           \
      std::cout << "  PASS [" << msg << "]\n";                                  \
      ++g_pass;                                                                 \
    }                                                                           \
  } while (0)

#define ASSERT_THROWS_MSG(expr, substr, msg)                                                              \
  do                                                                                                      \
  {                                                                                                       \
    bool caught = false;                                                                                  \
    std::string caught_what;                                                                              \
    try                                                                                                   \
    {                                                                                                     \
      expr;                                                                                               \
    }                                                                                                     \
    catch (const std::exception &e)                                                                       \
    {                                                                                                     \
      caught = true;                                                                                      \
      caught_what = e.what();                                                                             \
    }                                                                                                     \
    if (!caught)                                                                                          \
    {                                                                                                     \
      std::cerr << "  FAIL [" << msg << "]: expected exception, none thrown\n";                           \
      ++g_fail;                                                                                           \
    }                                                                                                     \
    else if (caught_what.find(substr) == std::string::npos)                                               \
    {                                                                                                     \
      std::cerr << "  FAIL [" << msg << "]: exception message '" << caught_what << "' does not contain '" \
                << substr << "'\n";                                                                       \
      ++g_fail;                                                                                           \
    }                                                                                                     \
    else                                                                                                  \
    {                                                                                                     \
      std::cout << "  PASS [" << msg << "]\n";                                                            \
      ++g_pass;                                                                                           \
    }                                                                                                     \
  } while (0)

static std::string writeTmpYaml(const std::string &content)
{
  static int counter = 0;
  std::string path = "/tmp/test_cfg_" + std::to_string(counter++) + ".yaml";
  std::ofstream f(path);
  f << content;
  return path;
}

// Test01–Test06 (T01-T06): valid config files

static void test_valid_configs(const std::string &root)
{
  std::cout << "\n[T01-T06] Valid config loading\n";

  struct Expected
  {
    std::string path;
    unsigned int N;
    std::vector<unsigned int> K;
    int meas_count; // -1 = don't check
  };

  std::vector<Expected> configs = {
      {root + "/config/1_continuum_robot.yaml", 1, {21}, 21},
      {root + "/config/2_parallel_continuum_robot.yaml", 2, {21, 21}, 1},
      {root + "/config/3_continuous_stewart_gough.yaml", 6, {15, 15, 15, 15, 15, 15}, 0},
      {root + "/config/4_collaborative_continuum_robots.yaml", 3, {31, 16, 16}, -1},
      {root + "/config/5_fbg_measurements.yaml", 2, {25, 21}, 46},
  };

  for (auto &cfg : configs)
  {
    try
    {
      ConfigLoader loader(cfg.path);
      auto topo = loader.getTopology();
      auto hp = loader.getHyperparameters();
      auto opts = loader.getOptions();
      auto meas = loader.getMeasurements();

      // T01: N correct
      ASSERT_TRUE(topo.N == cfg.N, "T01 N=" + std::to_string(cfg.N) + " in " + cfg.path);

      // T02: array sizes match N
      ASSERT_TRUE(topo.K.size() == cfg.N && topo.M.size() == cfg.N && topo.L.size() == cfg.N &&
                      topo.Ti0.size() == cfg.N && topo.lock_first_pose.size() == cfg.N &&
                      topo.lock_last_pose.size() == cfg.N,
                  "T02 per-robot array sizes for " + cfg.path);

      // T02b: K values match
      bool k_ok = (topo.K == cfg.K);
      ASSERT_TRUE(k_ok, "T02b K values for " + cfg.path);

      // T03: Ti0 are valid SE(3)
      bool ti0_ok = true;
      for (auto &T : topo.Ti0)
      {
        try
        {
          validate_transformation_matrix(T);
        }
        catch (...)
        {
          ti0_ok = false;
          break;
        }
      }
      ASSERT_TRUE(ti0_ok, "T03 Ti0 valid SE(3) for " + cfg.path);

      // T04: hyperparameter diagonals all positive
      bool hp_ok = true;
      for (int i = 0; i < 6; i++)
      {
        if (hp.Qc(i, i) <= 0)
        {
          hp_ok = false;
          break;
        }
      }
      ASSERT_TRUE(hp_ok, "T04 Qc positive diagonal for " + cfg.path);

      // T05/T06: measurement count
      if (cfg.meas_count >= 0)
      {
        ASSERT_TRUE((int)meas.size() == cfg.meas_count, "T05/T06 meas count=" + std::to_string(cfg.meas_count) +
                                                            " (got " + std::to_string(meas.size()) + ") for " +
                                                            cfg.path);
      }
    }
    catch (const std::exception &e)
    {
      std::cerr << "  FAIL [load " << cfg.path << "]: " << e.what() << "\n";
      ++g_fail;
    }
  }
}

// T07: translation transform

static void test_translation_transform()
{
  std::cout << "\n[T07] parseTransformMatrix translation\n";
  std::string yaml = R"(
topology:
  N: 1
  K: [5]
  M: [1]
  L: [0.1]
  lock_first_pose:   [true]
  lock_last_pose:    [false]
  lock_first_strain: [false]
  lock_last_strain:  [false]
  fbg_theta_offset:  [0]
  fbg_core_distance: [0]
  common_end_effector: false
  Ti0:
    - type: translation
      xyz: [0.1, 0.2, 0.3]
  robot_coupling: []
hyperparameters:
  noise_std:
    R_p: 0.001
    R_o: 0.01
    R_fbg: 1e-5
  R_pose_scale: 1.0
  R_fbg_strain_scale: 1.0
  R_coupling_scale: 1.0
  R_coupling_diagonal: [1,1,1,1,1,1]
  Qc_scale: 1.0
  Qc_diagonal: [1,1,1,1,1,1]
options:
  initial_guess: Straight
  solver: Newton
)";
  std::string path = writeTmpYaml(yaml);
  try
  {
    ConfigLoader loader(path);
    auto topo = loader.getTopology();
    auto &T = topo.Ti0[0];
    ASSERT_TRUE(std::abs(T(0, 3) - 0.1) < 1e-9 && std::abs(T(1, 3) - 0.2) < 1e-9 && std::abs(T(2, 3) - 0.3) < 1e-9,
                "T07 translation xyz values correct");
    // Rotation part should be identity
    Eigen::Matrix3d R = T.block<3, 3>(0, 0);
    ASSERT_TRUE((R - Eigen::Matrix3d::Identity()).norm() < 1e-9, "T07 translation rotation is identity");
  }
  catch (const std::exception &e)
  {
    std::cerr << "  FAIL [T07]: " << e.what() << "\n";
    ++g_fail;
  }
}

// T08: csv_file + orthonormalize

static void test_csv_transform(const std::string &root)
{
  std::cout << "\n[T08] parseTransformMatrix csv_file + orthonormalize\n";
  // Example 5 loads robot 1's base from CSV with orthonormalize:true
  try
  {
    ConfigLoader loader(root + "/config/5_fbg_measurements.yaml");
    auto topo = loader.getTopology();
    // Ti0[1] was loaded from CSV with orthonormalization
    Eigen::Matrix3d R = topo.Ti0[1].block<3, 3>(0, 0);
    double det = R.determinant();
    double ortho_err = (R * R.transpose() - Eigen::Matrix3d::Identity()).norm();
    ASSERT_TRUE(std::abs(det - 1.0) < 1e-6, "T08 csv R det=1 after orthonorm");
    ASSERT_TRUE(ortho_err < 1e-6, "T08 csv R^T R = I after orthonorm");
  }
  catch (const std::exception &e)
  {
    std::cerr << "  FAIL [T08]: " << e.what() << "\n";
    ++g_fail;
  }
}

// T09: matrix transform

static void test_matrix_transform()
{
  std::cout << "\n[T09] parseTransformMatrix matrix type\n";
  std::string yaml = R"(
topology:
  N: 1
  K: [5]
  M: [1]
  L: [0.1]
  lock_first_pose:   [true]
  lock_last_pose:    [false]
  lock_first_strain: [false]
  lock_last_strain:  [false]
  fbg_theta_offset:  [0]
  fbg_core_distance: [0]
  common_end_effector: false
  Ti0:
    - type: matrix
      data:
        - [1, 0, 0, 0.5]
        - [0, 1, 0, 0.0]
        - [0, 0, 1, 0.0]
        - [0, 0, 0, 1.0]
  robot_coupling: []
hyperparameters:
  noise_std:
    R_p: 0.001
    R_o: 0.01
    R_fbg: 1e-5
  R_pose_scale: 1.0
  R_fbg_strain_scale: 1.0
  R_coupling_scale: 1.0
  R_coupling_diagonal: [1,1,1,1,1,1]
  Qc_scale: 1.0
  Qc_diagonal: [1,1,1,1,1,1]
options:
  initial_guess: Straight
  solver: Newton
)";
  std::string path = writeTmpYaml(yaml);
  try
  {
    ConfigLoader loader(path);
    auto topo = loader.getTopology();
    ASSERT_TRUE(std::abs(topo.Ti0[0](0, 3) - 0.5) < 1e-9, "T09 matrix type translation x=0.5");
  }
  catch (const std::exception &e)
  {
    std::cerr << "  FAIL [T09]: " << e.what() << "\n";
    ++g_fail;
  }
}

// T10: missing file

static void test_missing_file()
{
  std::cout << "\n[T10] Missing config file\n";
  ASSERT_THROWS_MSG(ConfigLoader("/tmp/does_not_exist_xyz.yaml"), "Cannot open",
                    "T10 missing file throws 'Cannot open'");
}

// T11: missing topology section

static void test_missing_topology()
{
  std::cout << "\n[T11] Missing 'topology' section\n";
  std::string yaml = R"(
hyperparameters:
  Qc_diagonal: [1,1,1,1,1,1]
options:
  solver: Newton
)";
  ASSERT_THROWS_MSG(ConfigLoader(writeTmpYaml(yaml)), "topology", "T11 missing topology section");
}

// T12: N mismatch 

static void test_n_mismatch()
{
  std::cout << "\n[T12] N mismatch with per-robot arrays\n";
  std::string yaml = R"(
topology:
  N: 2
  K: [5]
  M: [1, 1]
  L: [0.1, 0.1]
  lock_first_pose:   [true, true]
  lock_last_pose:    [false, false]
  lock_first_strain: [false, false]
  lock_last_strain:  [false, false]
  fbg_theta_offset:  [0, 0]
  fbg_core_distance: [0, 0]
  common_end_effector: false
  Ti0:
    - type: identity
    - type: identity
  robot_coupling: []
hyperparameters:
  noise_std:
    R_p: 0.001
    R_o: 0.01
    R_fbg: 1e-5
  R_pose_scale: 1.0
  R_fbg_strain_scale: 1.0
  R_coupling_scale: 1.0
  R_coupling_diagonal: [1,1,1,1,1,1]
  Qc_scale: 1.0
  Qc_diagonal: [1,1,1,1,1,1]
options:
  solver: Newton
)";
  // K has only 1 element but N=2 → should throw
  ASSERT_THROWS(ConfigLoader(writeTmpYaml(yaml)), "T12 N mismatch K array throws");
}

// T13: invalid initial_guess string

static void test_invalid_initial_guess()
{
  std::cout << "\n[T13] Invalid initial_guess value\n";
  std::string yaml = R"(
topology:
  N: 1
  K: [5]
  M: [1]
  L: [0.1]
  lock_first_pose:   [true]
  lock_last_pose:    [false]
  lock_first_strain: [false]
  lock_last_strain:  [false]
  fbg_theta_offset:  [0]
  fbg_core_distance: [0]
  common_end_effector: false
  Ti0:
    - type: identity
  robot_coupling: []
hyperparameters:
  noise_std:
    R_p: 0.001
    R_o: 0.01
    R_fbg: 1e-5
  R_pose_scale: 1.0
  R_fbg_strain_scale: 1.0
  R_coupling_scale: 1.0
  R_coupling_diagonal: [1,1,1,1,1,1]
  Qc_scale: 1.0
  Qc_diagonal: [1,1,1,1,1,1]
options:
  initial_guess: BadValue
  solver: Newton
)";
  ASSERT_THROWS_MSG(ConfigLoader(writeTmpYaml(yaml)), "unknown value", "T13 invalid initial_guess");
}

// T14: translation missing xyz

static void test_translation_missing_xyz()
{
  std::cout << "\n[T14] Translation type missing xyz field\n";
  std::string yaml = R"(
topology:
  N: 1
  K: [5]
  M: [1]
  L: [0.1]
  lock_first_pose:   [true]
  lock_last_pose:    [false]
  lock_first_strain: [false]
  lock_last_strain:  [false]
  fbg_theta_offset:  [0]
  fbg_core_distance: [0]
  common_end_effector: false
  Ti0:
    - type: translation
  robot_coupling: []
hyperparameters:
  noise_std:
    R_p: 0.001
    R_o: 0.01
    R_fbg: 1e-5
  R_pose_scale: 1.0
  R_fbg_strain_scale: 1.0
  R_coupling_scale: 1.0
  R_coupling_diagonal: [1,1,1,1,1,1]
  Qc_scale: 1.0
  Qc_diagonal: [1,1,1,1,1,1]
options:
  solver: Newton
)";
  ASSERT_THROWS_MSG(ConfigLoader(writeTmpYaml(yaml)), "xyz", "T14 translation missing xyz gives descriptive error");
}

// T15: mask with wrong number of elements

static void test_mask_size_validation()
{
  std::cout << "\n[T15] Mask with wrong number of elements\n";
  std::string yaml = R"(
topology:
  N: 1
  K: [5]
  M: [1]
  L: [0.1]
  lock_first_pose:   [true]
  lock_last_pose:    [false]
  lock_first_strain: [false]
  lock_last_strain:  [false]
  fbg_theta_offset:  [0]
  fbg_core_distance: [0]
  common_end_effector: false
  Ti0:
    - type: identity
  robot_coupling: []
hyperparameters:
  noise_std:
    R_p: 0.001
    R_o: 0.01
    R_fbg: 1e-5
  R_pose_scale: 1.0
  R_fbg_strain_scale: 1.0
  R_coupling_scale: 1.0
  R_coupling_diagonal: [1,1,1,1,1,1]
  Qc_scale: 1.0
  Qc_diagonal: [1,1,1,1,1,1]
options:
  solver: Newton
measurements:
  - type: Strain
    idx_robot: 0
    idx_node: 0
    mask: [1, 1, 1]
    value: [0, 0, 0, 0, 1, 0]
)";
  ASSERT_THROWS_MSG(ConfigLoader(writeTmpYaml(yaml)), "6 elements",
                    "T15 mask with 3 elements gives descriptive error");
}

// T_extra: missing R_pose_scale gives descriptive error

static void test_missing_scale_field()
{
  std::cout << "\n[T_extra] Missing R_pose_scale gives descriptive error\n";
  std::string yaml = R"(
topology:
  N: 1
  K: [5]
  M: [1]
  L: [0.1]
  lock_first_pose:   [true]
  lock_last_pose:    [false]
  lock_first_strain: [false]
  lock_last_strain:  [false]
  fbg_theta_offset:  [0]
  fbg_core_distance: [0]
  common_end_effector: false
  Ti0:
    - type: identity
  robot_coupling: []
hyperparameters:
  noise_std:
    R_p: 0.001
    R_o: 0.01
    R_fbg: 1e-5
  # R_pose_scale intentionally omitted
  R_fbg_strain_scale: 1.0
  R_coupling_scale: 1.0
  R_coupling_diagonal: [1,1,1,1,1,1]
  Qc_scale: 1.0
  Qc_diagonal: [1,1,1,1,1,1]
options:
  solver: Newton
)";
  ASSERT_THROWS_MSG(ConfigLoader(writeTmpYaml(yaml)), "R_pose_scale",
                    "T_extra missing R_pose_scale gives descriptive error");
}

// VisualizationSettings defaults

static void test_vis_defaults()
{
  std::cout << "\n[T_vis] VisualizationSettings defaults when section absent\n";
  std::string yaml = R"(
topology:
  N: 1
  K: [5]
  M: [1]
  L: [0.1]
  lock_first_pose:   [true]
  lock_last_pose:    [false]
  lock_first_strain: [false]
  lock_last_strain:  [false]
  fbg_theta_offset:  [0]
  fbg_core_distance: [0]
  common_end_effector: false
  Ti0:
    - type: identity
  robot_coupling: []
hyperparameters:
  noise_std:
    R_p: 0.001
    R_o: 0.01
    R_fbg: 1e-5
  R_pose_scale: 1.0
  R_fbg_strain_scale: 1.0
  R_coupling_scale: 1.0
  R_coupling_diagonal: [1,1,1,1,1,1]
  Qc_scale: 1.0
  Qc_diagonal: [1,1,1,1,1,1]
options:
  solver: Newton
)";
  try
  {
    ConfigLoader loader(writeTmpYaml(yaml));
    auto vis = loader.getVisualizationSettings();
    ASSERT_TRUE(vis.window_width == 1280, "T_vis default window_width=1280");
    ASSERT_TRUE(vis.window_height == 720, "T_vis default window_height=720");
    ASSERT_TRUE(vis.render_frames == true, "T_vis default render_frames=true");
    ASSERT_TRUE(vis.covariance_n_std == 3, "T_vis default covariance_n_std=3");
  }
  catch (const std::exception &e)
  {
    std::cerr << "  FAIL [T_vis]: " << e.what() << "\n";
    ++g_fail;
  }
}

// main

int main(int argc, char *argv[])
{
  // Accept optional project root directory as first argument
  std::string root = (argc > 1) ? argv[1] : "..";

  std::cout << "=== ConfigLoader Unit Tests ===\n";
  std::cout << "Project root: " << root << "\n";

  test_valid_configs(root);
  test_translation_transform();
  test_csv_transform(root);
  test_matrix_transform();
  test_missing_file();
  test_missing_topology();
  test_n_mismatch();
  test_invalid_initial_guess();
  test_translation_missing_xyz();
  test_mask_size_validation();
  test_missing_scale_field();
  test_vis_defaults();

  std::cout << "\n=== Results: " << g_pass << " passed, " << g_fail << " failed ===\n";
  return (g_fail == 0) ? 0 : 1;
}
