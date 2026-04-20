// Combined driver: CosseratRodModel FK -> ContinuumRodPriors -> estimator.
// Requires USE_LOCAL_TDCR=ON. See doc/bridge_adapter_guide.md for the runtime
// walkthrough and doc/cosserat_integration_results.md for measured A/B results.
//
// Pipeline:
//   1. Load YAML config.
//   2. Run Cosserat FK.
//   3. Extract ContinuumRodPriors DTO (Data Transfer Object).
//   4. Convert to estimator-ready measurements + control inputs + initial guess.
//   5. Run computeStateEstimate().
//   6. Print per-node strain comparison: Cosserat prediction vs estimator result.
//   7. If --visualize is passed, render the estimator's SystemState in VTK.
//
// Flags:
//   --visualize            open a VTK window showing the backbone + frames + cov.
//   --no-control-inputs    skip the adapter's ControlInput injection (for
//                          A/B comparison of prior vs measurement influence).

#include "config_loader.h"
#include "continuum_robot_state_estimator.h"
#include "continuum_rod_priors.h"
#include "cosserat_priors_adapter.h"
#include "cosseratrodmodel.h"
#include "visualizervtk.h"

#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include <vtkAutoInit.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkRenderingFreeType);
VTK_MODULE_INIT(vtkInteractionStyle);

namespace {

std::string fmt(double v, int width = 10, int prec = 4)
{
    std::ostringstream os;
    os << std::setw(width) << std::fixed << std::setprecision(prec) << v;
    return os.str();
}

std::string fmtSci(double v, int width = 12, int prec = 3)
{
    std::ostringstream os;
    os << std::setw(width) << std::scientific << std::setprecision(prec) << v;
    return os.str();
}

void printBanner(const std::string& title)
{
    std::cout << "========================================\n  " << title
              << "\n========================================\n";
}

} // anonymous namespace


int main(int argc, char* argv[])
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0]
                  << " <config.yaml> [--visualize] [--no-control-inputs]\n";
        return 1;
    }
    const std::string config_path = argv[1];
    bool use_control_inputs = true;
    bool visualize          = false;
    for (int i = 2; i < argc; ++i) {
        const std::string a = argv[i];
        if      (a == "--no-control-inputs") use_control_inputs = false;
        else if (a == "--visualize")         visualize          = true;
        else {
            std::cerr << "Unknown flag: " << a << "\n";
            return 1;
        }
    }

    // --- 1. Load YAML config -------------------------------------------------
    ConfigLoader config(config_path);
    auto topology       = config.getTopology();
    auto params         = config.getHyperparameters();
    auto options        = config.getOptions();
    auto measurements   = config.getMeasurements();
    auto control_inputs = config.getControlInputs();
    const auto vis_settings = config.getVisualizationSettings();

    printBanner("Cosserat + Estimator Combined Driver");
    std::cout << "\nConfig: " << config_path << "\n";
    std::cout << "  Robots: N = " << topology.N
              << ",   K[0] = " << topology.K.at(0)
              << ",   L[0] = " << topology.L.at(0) << " m\n\n";

    // --- 2. Run Cosserat FK --------------------------------------------------
    CosseratRodModel model;
    Eigen::MatrixXd diskFrames;
    Eigen::Matrix<double, 6, 1> q;
    q << 0.5, 0.2, 0.0, 0.3, 0.0, 0.1;          // default tendon tensions [N]
    const Eigen::Vector3d f_ext = Eigen::Vector3d::Zero();
    const Eigen::Vector3d l_ext = Eigen::Vector3d::Zero();
    const double L1 = 0.1, L2 = 0.1;            // default Cosserat segment lengths

    std::cout << "Cosserat FK\n";
    std::cout << "  Tendons [N]:       ["
              << q(0) << ", " << q(1) << ", " << q(2) << " | "
              << q(3) << ", " << q(4) << ", " << q(5) << "]\n";
    std::cout << "  Segment lengths:   L1 = " << L1 << ",  L2 = " << L2 << "\n";

    const bool fk_ok = model.forwardKinematics(diskFrames, q, f_ext, l_ext);
    if (!fk_ok || !model.hasAuxOutputs()) {
        std::cerr << "  FAILED (residual = " << model.getFinalResudial() << ")\n";
        return 1;
    }
    std::cout << "  Converged:         residual = " << fmtSci(model.getFinalResudial()) << "\n\n";

    // --- 3. Extract DTO (Data Transfer Object) ------------------------------
    const ContinuumRodPriors priors = priorsFromCosseratModel(model, L1, L2);

    // --- 4. Convert to estimator inputs -------------------------------------
    const unsigned int robot_idx = 0;
    const EstimatorPriors ep = cosseratPriorsToEstimator(priors, topology, robot_idx, diskFrames);

    std::cout << "Adapter (Cosserat -> Estimator)\n";
    std::cout << "  Strain measurements:  " << ep.measurements.size() << "\n";
    std::cout << "  Control inputs:       " << ep.control_inputs.size() << "\n";
    std::cout << "  Initial-guess nodes:  "
              << ep.initial_guess.robots.at(robot_idx).estimation_nodes.size() << "\n\n";

    // Merge adapter outputs with whatever the config already had (both empty in
    // config/6_cosserat_priors.yaml, but user may add FBG / pose measurements later).
    measurements.insert(measurements.end(),
                        ep.measurements.begin(), ep.measurements.end());
    if (use_control_inputs) {
        control_inputs.insert(control_inputs.end(),
                              ep.control_inputs.begin(), ep.control_inputs.end());
    } else {
        std::cout << "  (--no-control-inputs: skipping adapter control inputs)\n";
    }

    // Seed the estimator's optimizer with the Cosserat-predicted poses.
    options.custom_guess = ep.initial_guess;

    // --- 5. Run estimator ---------------------------------------------------
    ContinuumRobotStateEstimator estimator(topology, params, options);
    ContinuumRobotStateEstimator::SystemState state;
    std::vector<double> cost;

    const bool est_ok = estimator.computeStateEstimate(
        state, cost, measurements, control_inputs, vis_settings.verbose);

    std::cout << "Estimator\n";
    std::cout << "  Iterations:        " << cost.size() << "\n";
    if (!cost.empty()) {
        std::cout << "  Initial cost:      " << fmtSci(cost.front()) << "\n";
        std::cout << "  Final   cost:      " << fmtSci(cost.back())  << "\n";
    }
    std::cout << "  Converged:         " << (est_ok ? "yes" : "NO") << "\n\n";

    // --- 6. Per-node strain comparison (all 6 components) -------------------
    const auto& nodes = state.robots.at(robot_idx).estimation_nodes;
    const int K = static_cast<int>(nodes.size());

    const char* labels[6] = { "nu_1", "nu_2", "nu_3", "om_1", "om_2", "om_3" };

    std::cout << "Per-node strain comparison (Cosserat permuted to estimator body frame)\n\n";

    // Cosserat side table
    std::cout << "  Cosserat prior\n";
    std::cout << "    k   s[m]  ";
    for (int c = 0; c < 6; ++c) std::cout << std::setw(10) << labels[c];
    std::cout << "\n";
    for (int k = 0; k < K; ++k) {
        const Eigen::Matrix<double,6,1> s_cos = ep.measurements.at(k).value;
        std::cout << "   " << std::setw(2) << k << "  "
                  << fmt(nodes.at(k).arclength, 5, 3);
        for (int c = 0; c < 6; ++c) std::cout << " " << fmt(s_cos(c), 9, 4);
        std::cout << "\n";
    }

    // Estimator side table
    std::cout << "\n  Estimator result\n";
    std::cout << "    k   s[m]  ";
    for (int c = 0; c < 6; ++c) std::cout << std::setw(10) << labels[c];
    std::cout << "\n";
    for (int k = 0; k < K; ++k) {
        const Eigen::Matrix<double,6,1> s_est = nodes.at(k).strain;
        std::cout << "   " << std::setw(2) << k << "  "
                  << fmt(nodes.at(k).arclength, 5, 3);
        for (int c = 0; c < 6; ++c) std::cout << " " << fmt(s_est(c), 9, 4);
        std::cout << "\n";
    }

    // Per-component max-abs diff across all nodes
    Eigen::Matrix<double,6,1> max_diff_per_component = Eigen::Matrix<double,6,1>::Zero();
    double max_diff = 0.0;
    for (int k = 0; k < K; ++k) {
        const Eigen::Matrix<double,6,1> s_cos = ep.measurements.at(k).value;
        const Eigen::Matrix<double,6,1> s_est = nodes.at(k).strain;
        const Eigen::Matrix<double,6,1> d = (s_cos - s_est).cwiseAbs();
        max_diff_per_component = max_diff_per_component.cwiseMax(d);
        max_diff = std::max(max_diff, d.maxCoeff());
    }

    std::cout << "\n  Max |diff| per component across all nodes:\n    ";
    for (int c = 0; c < 6; ++c)
        std::cout << labels[c] << "=" << fmtSci(max_diff_per_component(c), 10, 2) << "  ";
    std::cout << "\n  Max |diff|_inf overall: " << fmtSci(max_diff) << "\n\n";

    printBanner("Pipeline OK");

    // --- 7. Optional visualization -----------------------------------------
    if (visualize) {
        std::cout << "\nOpening VTK viewer. Close the window to exit.\n";
        Visualizer vis(topology, vis_settings);
        vis.update(state,
                   vis_settings.render_frames,
                   vis_settings.render_covariance,
                   vis_settings.covariance_n_std);

        auto style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
        auto iren  = vtkSmartPointer<vtkRenderWindowInteractor>::New();
        iren->SetRenderWindow(vis.getRenderWindow());
        iren->UpdateSize(vis_settings.window_width, vis_settings.window_height);
        iren->SetInteractorStyle(style);
        iren->Initialize();
        iren->Start();
    }

    return est_ok ? 0 : 2;
}
