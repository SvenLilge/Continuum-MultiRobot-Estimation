// test_estimation.cpp
// Headless integration tests: runs computeStateEstimate() for all 5 examples.
// No VTK / no GUI required.
//
// Compile alongside: config_loader.cpp, continuum_robot_state_estimator.cpp, utilities.cpp
// Link: yaml-cpp, Eigen3
//
// Run from the project root:
//   ./examples/test_estimation [project_root]

#include "config_loader.h"
#include "continuum_robot_state_estimator.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

//  helpers

static int g_pass = 0;
static int g_fail = 0;

#define ASSERT_TRUE(cond, msg)                       \
    do                                               \
    {                                                \
        if (!(cond))                                 \
        {                                            \
            std::cerr << "  FAIL [" << msg << "]\n"; \
            ++g_fail;                                \
        }                                            \
        else                                         \
        {                                            \
            std::cout << "  PASS [" << msg << "]\n"; \
            ++g_pass;                                \
        }                                            \
    } while (0)

// Returns true if all finite (no NaN / Inf) in the matrix.
static bool allFinite(const Eigen::MatrixXd &m)
{
    return m.allFinite();
}

// Check all poses and strains in a robot state for NaN/Inf.
static bool stateIsFinite(const ContinuumRobotStateEstimator::SystemState &state)
{
    for (auto &robot : state.robots)
    {
        for (auto &node : robot.estimation_nodes)
        {
            if (!allFinite(node.pose))
                return false;
            if (!allFinite(node.strain))
                return false;
            if (!allFinite(node.pose_std))
                return false;
            if (!allFinite(node.strain_std))
                return false;
        }
        for (auto &node : robot.interpolation_nodes)
        {
            if (!allFinite(node.pose))
                return false;
            if (!allFinite(node.strain))
                return false;
        }
    }
    return true;
}

// Run one estimation case and apply standard checks.
static void runCase(const std::string &label, const std::string &config_path, bool expect_convergence = true)
{
    std::cout << "\n[" << label << "] " << config_path << "\n";
    try
    {
        ConfigLoader cfg(config_path);
        ContinuumRobotStateEstimator estimator(cfg.getTopology(), cfg.getHyperparameters(), cfg.getOptions());

        ContinuumRobotStateEstimator::SystemState state;
        std::vector<double> cost;
        bool converged = estimator.computeStateEstimate(state, cost, cfg.getMeasurements(), /*verbose=*/false);

        // Convergence
        if (expect_convergence)
        {
            ASSERT_TRUE(converged, label + " converged");
        }
        else
        {
            std::cout << "  INFO [" << label << "] converged=" << converged << "\n";
        }

        // Cost decreased
        ASSERT_TRUE(!cost.empty(), label + " cost vector non-empty");
        if (cost.size() >= 2)
        {
            ASSERT_TRUE(cost.back() <= cost.front(), label + " cost decreased over iterations");
        }

        // No NaN / Inf
        ASSERT_TRUE(stateIsFinite(state), label + " all state values finite");

        // At least one robot in the state
        ASSERT_TRUE(!state.robots.empty(), label + " state has robots");

        // Estimation nodes count matches K
        auto topo = cfg.getTopology();
        bool nodes_ok = true;
        for (unsigned int r = 0; r < topo.N; r++)
        {
            if (state.robots[r].estimation_nodes.size() != topo.K[r])
                nodes_ok = false;
        }
        ASSERT_TRUE(nodes_ok, label + " estimation node counts match K");

        // Uncertainties positive
        bool std_ok = true;
        for (auto &robot : state.robots)
        {
            for (auto &node : robot.estimation_nodes)
            {
                for (int i = 0; i < 6; i++)
                {
                    if (node.pose_std(i) < 0 || node.strain_std(i) < 0)
                        std_ok = false;
                }
            }
        }
        ASSERT_TRUE(std_ok, label + " uncertainties non-negative");
    }
    catch (const std::exception &e)
    {
        std::cerr << "  FAIL [" << label << "]: exception: " << e.what() << "\n";
        ++g_fail;
    }
}

// Extra checks for specific examples

// Example 1: single robot — tip should move from base (non-zero curvature)
static void check_example1_shape(const std::string &config_path)
{
    std::cout << "\n[Ex1-shape] Tip displacement > 0\n";
    try
    {
        ConfigLoader cfg(config_path);
        ContinuumRobotStateEstimator estimator(cfg.getTopology(), cfg.getHyperparameters(), cfg.getOptions());
        ContinuumRobotStateEstimator::SystemState state;
        std::vector<double> cost;
        estimator.computeStateEstimate(state, cost, cfg.getMeasurements(), false);

        auto &tip = state.robots[0].estimation_nodes.back();
        double tip_norm = tip.pose.block<3, 1>(0, 3).norm();
        ASSERT_TRUE(tip_norm > 0.01, "Ex1 tip position norm > 0.01 m (robot is not straight)");
    }
    catch (const std::exception &e)
    {
        std::cerr << "  FAIL [Ex1-shape]: " << e.what() << "\n";
        ++g_fail;
    }
}

// Example 2: two parallel robots — tips separated in y by ~0.1 m
static void check_example2_separation(const std::string &config_path)
{
    std::cout << "\n[Ex2-sep] Two robot tips y-separation\n";
    try
    {
        ConfigLoader cfg(config_path);
        ContinuumRobotStateEstimator estimator(cfg.getTopology(), cfg.getHyperparameters(), cfg.getOptions());
        ContinuumRobotStateEstimator::SystemState state;
        std::vector<double> cost;
        estimator.computeStateEstimate(state, cost, cfg.getMeasurements(), false);

        Eigen::Vector3d t0 = state.robots[0].estimation_nodes.back().pose.block<3, 1>(0, 3);
        Eigen::Vector3d t1 = state.robots[1].estimation_nodes.back().pose.block<3, 1>(0, 3);
        double dy = std::abs(t0(1) - t1(1));
        ASSERT_TRUE(dy > 0.05, "Ex2 robots separated in y by > 0.05 m (got " + std::to_string(dy) + " m)");
    }
    catch (const std::exception &e)
    {
        std::cerr << "  FAIL [Ex2-sep]: " << e.what() << "\n";
        ++g_fail;
    }
}

// Example 5: FBG — both robots have curvature (non-trivial strain)
static void check_example5_fbg(const std::string &config_path)
{
    std::cout << "\n[Ex5-fbg] FBG-driven robots have curvature\n";
    try
    {
        ConfigLoader cfg(config_path);
        ContinuumRobotStateEstimator estimator(cfg.getTopology(), cfg.getHyperparameters(), cfg.getOptions());
        ContinuumRobotStateEstimator::SystemState state;
        std::vector<double> cost;
        estimator.computeStateEstimate(state, cost, cfg.getMeasurements(), false);

        for (unsigned int r = 0; r < 2; r++)
        {
            double max_omega = 0;
            for (auto &node : state.robots[r].estimation_nodes)
                max_omega = std::max(max_omega, node.strain.tail<3>().norm());
            ASSERT_TRUE(max_omega > 1e-4, "Ex5 robot " + std::to_string(r) + " has non-zero curvature");
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << "  FAIL [Ex5-fbg]: " << e.what() << "\n";
        ++g_fail;
    }
}

// queryAdditionalStates — interpolation at mid-length
static void test_query_additional_states(const std::string &config_path)
{
    std::cout << "\n[T_query] queryAdditionalStates\n";
    try
    {
        ConfigLoader cfg(config_path);
        auto topo = cfg.getTopology();
        ContinuumRobotStateEstimator estimator(topo, cfg.getHyperparameters(), cfg.getOptions());
        ContinuumRobotStateEstimator::SystemState state;
        std::vector<double> cost;
        estimator.computeStateEstimate(state, cost, cfg.getMeasurements(), false);

        double mid = topo.L[0] / 2.0;
        estimator.queryAdditionalStates(state, {{0, mid}});

        ASSERT_TRUE(!state.robots[0].queried_nodes.empty(), "T_query queried_nodes non-empty after call");
        ASSERT_TRUE(allFinite(state.robots[0].queried_nodes[0].pose), "T_query queried node pose is finite");
        ASSERT_TRUE(allFinite(state.robots[0].queried_nodes[0].strain), "T_query queried node strain is finite");
    }
    catch (const std::exception &e)
    {
        std::cerr << "  FAIL [T_query]: " << e.what() << "\n";
        ++g_fail;
    }
}

// max_iterations: 1 — should NOT converge
static void test_no_convergence(const std::string &config_path)
{
    std::cout << "\n[T_noconv] max_iterations=1 returns false\n";
    try
    {
        ConfigLoader cfg(config_path);
        auto opts = cfg.getOptions();
        opts.max_optimization_iterations = 1;

        ContinuumRobotStateEstimator estimator(cfg.getTopology(), cfg.getHyperparameters(), opts);
        ContinuumRobotStateEstimator::SystemState state;
        std::vector<double> cost;
        bool converged = estimator.computeStateEstimate(state, cost, cfg.getMeasurements(), false);

        ASSERT_TRUE(!converged, "T_noconv max_iterations=1 returns false");
    }
    catch (const std::exception &e)
    {
        std::cerr << "  FAIL [T_noconv]: " << e.what() << "\n";
        ++g_fail;
    }
}

// kirchhoff_rods: false — 6D strain, should still run
static void test_kirchhoff_off(const std::string &config_path)
{
    std::cout << "\n[T_kirchhoff] kirchhoff_rods=false still converges\n";
    try
    {
        ConfigLoader cfg(config_path);
        auto opts = cfg.getOptions();
        opts.kirchhoff_rods = false;

        ContinuumRobotStateEstimator estimator(cfg.getTopology(), cfg.getHyperparameters(), opts);
        ContinuumRobotStateEstimator::SystemState state;
        std::vector<double> cost;
        bool converged = estimator.computeStateEstimate(state, cost, cfg.getMeasurements(), false);

        ASSERT_TRUE(converged, "T_kirchhoff converged with kirchhoff_rods=false");
        ASSERT_TRUE(stateIsFinite(state), "T_kirchhoff state finite");
    }
    catch (const std::exception &e)
    {
        std::cerr << "  FAIL [T_kirchhoff]: " << e.what() << "\n";
        ++g_fail;
    }
}

// warm start: 'Last' initial guess after one solve — second solve ≤ first
static void test_warm_start(const std::string &config_path)
{
    std::cout << "\n[T_warm] 'Last' initial guess (warm start)\n";
    try
    {
        ConfigLoader cfg(config_path);
        auto opts = cfg.getOptions();
        opts.init_guess_type = ContinuumRobotStateEstimator::Options::Last;

        ContinuumRobotStateEstimator estimator(cfg.getTopology(), cfg.getHyperparameters(), opts);
        auto meas = cfg.getMeasurements();

        ContinuumRobotStateEstimator::SystemState state;
        std::vector<double> cost1, cost2;

        estimator.computeStateEstimate(state, cost1, meas, false);
        estimator.computeStateEstimate(state, cost2, meas, false);

        // With warm start the second solve should need ≤ iterations than the first
        ASSERT_TRUE(cost2.size() <= cost1.size(), "T_warm second solve needs <= iterations than first");
    }
    catch (const std::exception &e)
    {
        std::cerr << "  FAIL [T_warm]: " << e.what() << "\n";
        ++g_fail;
    }
}

// main

int main(int argc, char *argv[])
{
    std::string root = (argc > 1) ? argv[1] : "..";

    std::cout << "=== Estimation Integration Tests ===\n";
    std::cout << "Project root: " << root << "\n";

    std::string cfg1 = root + "/config/1_continuum_robot.yaml";
    std::string cfg2 = root + "/config/2_parallel_continuum_robot.yaml";
    std::string cfg3 = root + "/config/3_continuous_stewart_gough.yaml";
    std::string cfg4 = root + "/config/4_collaborative_continuum_robots.yaml";
    std::string cfg5 = root + "/config/5_fbg_measurements.yaml";

    // Standard convergence + sanity checks for all 5 examples
    runCase("Ex1", cfg1);
    runCase("Ex2", cfg2);
    runCase("Ex3", cfg3);
    runCase("Ex4", cfg4);
    runCase("Ex5", cfg5);

    // Example-specific shape checks
    check_example1_shape(cfg1);
    check_example2_separation(cfg2);
    check_example5_fbg(cfg5);

    // Cross-cutting feature tests (run on example 1 for speed)
    test_query_additional_states(cfg1);
    test_no_convergence(cfg1);
    test_kirchhoff_off(cfg1);
    test_warm_start(cfg1);

    std::cout << "\n=== Results: " << g_pass << " passed, " << g_fail << " failed ===\n";
    return (g_fail == 0) ? 0 : 1;
}
