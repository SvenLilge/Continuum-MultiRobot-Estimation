// Combined driver: CosseratRodModel FK -> ContinuumRodPriors -> estimator.
// Requires USE_LOCAL_TDCR=ON. See doc/cosserat_integration.md for how the bridge
// works and for the measured A/B results.
//
// Pipeline:
//   1. Load YAML config.
//   2. Run Cosserat FK.
//   3. Extract ContinuumRodPriors DTO (Data Transfer Object).
//   4. Convert to estimator-ready measurements + control inputs + initial guess.
//   5. Run computeStateEstimate().
//   6. Print per-node strain and position comparison: Cosserat vs estimator.
//      The position check is the decisive one: a strain-only check cannot see
//      a wrong shape when control inputs are used (state strain = bias only).
//   7. If --visualize is passed, render the estimator's SystemState in VTK.
//
// Flags:
//   --visualize            open a VTK window showing the backbone + frames + cov.
//   --no-control-inputs    feed the model strain as strain measurements
//                          (ControlInputMode::None) instead of as a velocity
//                          input (ControlInputMode::StrainAsInput, default).

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
                  << " <config.yaml> [--visualize] [--no-control-inputs]"
                  << " [--q F1,F2,F3,F4,F5,F6]\n";
        return 1;
    }
    const std::string config_path = argv[1];
    bool use_control_inputs = true;
    bool visualize          = false;
    // Default tendon tensions; overridable via --q so the viewer can show the
    // same S-shape the evaluation driver uses (q = 6,0,0,0,8,8 N).
    Eigen::Matrix<double, 6, 1> q;
    q << 0.5, 0.2, 0.0, 0.3, 0.0, 0.1;
    for (int i = 2; i < argc; ++i) {
        const std::string a = argv[i];
        if      (a == "--no-control-inputs") use_control_inputs = false;
        else if (a == "--visualize")         visualize          = true;
        else if (a == "--q") {
            if (i + 1 >= argc) {
                std::cerr << "--q requires 6 comma-separated values\n";
                return 1;
            }
            std::string vals = argv[++i];
            for (char& c : vals) if (c == ',') c = ' ';
            std::istringstream is(vals);
            if (!(is >> q(0) >> q(1) >> q(2) >> q(3) >> q(4) >> q(5))) {
                std::cerr << "--q must be 6 comma-separated numbers\n";
                return 1;
            }
        }
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
    // q was parsed from CLI above (default = [0.5, 0.2, 0, 0.3, 0, 0.1] N).
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
    const ContinuumRodPriors priors = priorsFromCosseratModel(model, L1, L2, diskFrames);

    // --- 4. Convert to estimator inputs -------------------------------------
    const unsigned int robot_idx = 0;
    const ControlInputMode mode =
        use_control_inputs ? ControlInputMode::StrainAsInput : ControlInputMode::None;
    const EstimatorPriors ep = cosseratPriorsToEstimator(priors, topology, robot_idx, diskFrames, mode);

    std::cout << "Adapter (Cosserat -> Estimator)\n";
    std::cout << "  Strain measurements:  " << ep.measurements.size() << "\n";
    std::cout << "  Control inputs:       " << ep.control_inputs.size() << "\n";
    std::cout << "  Initial-guess nodes:  "
              << ep.initial_guess.robots.at(robot_idx).estimation_nodes.size() << "\n\n";

    // Merge adapter outputs with whatever the config already had (both empty in
    // config/6_cosserat_priors.yaml, but user may add FBG / pose measurements later).
    measurements.insert(measurements.end(),
                        ep.measurements.begin(), ep.measurements.end());
    control_inputs.insert(control_inputs.end(),
                          ep.control_inputs.begin(), ep.control_inputs.end());
    std::cout << (use_control_inputs
        ? "  Mode: model strain as velocity input; state strain = bias, target 0\n"
        : "  Mode: model strain as strain measurements (--no-control-inputs)\n");

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
    std::cout << "  Cosserat model strain\n";
    std::cout << "    k   s[m]  ";
    for (int c = 0; c < 6; ++c) std::cout << std::setw(10) << labels[c];
    std::cout << "\n";
    for (int k = 0; k < K; ++k) {
        const Eigen::Matrix<double,6,1> s_cos = ep.model_strain.at(k);
        std::cout << "   " << std::setw(2) << k << "  "
                  << fmt(nodes.at(k).arclength, 5, 3);
        for (int c = 0; c < 6; ++c) std::cout << " " << fmt(s_cos(c), 9, 4);
        std::cout << "\n";
    }

    // Estimator side table
    std::cout << (use_control_inputs
        ? "\n  Estimator state strain (bias; total = bias + velocity input, target bias = 0)\n"
        : "\n  Estimator state strain\n");
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

    // Per-component max-abs diff between state strain and its target across all nodes
    Eigen::Matrix<double,6,1> max_diff_per_component = Eigen::Matrix<double,6,1>::Zero();
    double max_diff = 0.0;
    for (int k = 0; k < K; ++k) {
        const Eigen::Matrix<double,6,1> s_cos = ep.measurements.at(k).value;  // target
        const Eigen::Matrix<double,6,1> s_est = nodes.at(k).strain;
        const Eigen::Matrix<double,6,1> d = (s_cos - s_est).cwiseAbs();
        max_diff_per_component = max_diff_per_component.cwiseMax(d);
        max_diff = std::max(max_diff, d.maxCoeff());
    }

    std::cout << "\n  Max |state strain - target| per component across all nodes:\n    ";
    for (int c = 0; c < 6; ++c)
        std::cout << labels[c] << "=" << fmtSci(max_diff_per_component(c), 10, 2) << "  ";
    std::cout << "\n  Max |diff|_inf overall: " << fmtSci(max_diff) << "\n\n";

    // Position comparison: estimated node positions vs the Cosserat shape
    // (the initial guess holds the Cosserat poses at the estimator nodes).
    const auto& ref_nodes = ep.initial_guess.robots.at(robot_idx).estimation_nodes;
    std::cout << "Per-node position comparison [mm]\n";
    std::cout << "    k   s[m]          Cosserat p                 estimate p           |err|\n";
    double max_pos_err_mm = 0.0;
    for (int k = 0; k < K; ++k) {
        const Eigen::Vector3d p_cos = ref_nodes.at(k).pose.block<3,1>(0,3) * 1e3;
        const Eigen::Vector3d p_est = nodes.at(k).pose.block<3,1>(0,3) * 1e3;
        const double err = (p_cos - p_est).norm();
        max_pos_err_mm = std::max(max_pos_err_mm, err);
        std::cout << "   " << std::setw(2) << k << "  " << fmt(nodes.at(k).arclength, 5, 3) << "  ";
        for (int c = 0; c < 3; ++c) std::cout << fmt(p_cos(c), 8, 2);
        std::cout << "  ";
        for (int c = 0; c < 3; ++c) std::cout << fmt(p_est(c), 8, 2);
        std::cout << fmt(err, 9, 3) << "\n";
    }
    std::cout << "\n  Max position error: " << fmt(max_pos_err_mm, 0, 3) << " mm\n\n";

    // 2 mm allows for the junction smoothing (~1.1 mm at the tip in mode None).
    const bool shape_ok = max_pos_err_mm < 2.0;
    printBanner(shape_ok ? "Pipeline OK" : "Pipeline FAILED: estimated shape does not match Cosserat");

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

    if (!shape_ok) return 3;
    return est_ok ? 0 : 2;
}
