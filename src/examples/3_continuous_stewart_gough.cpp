#include "config_loader.h"
#include "visualizervtk.h"
#include "continuum_robot_state_estimator.h"

#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>

// VTK Factory initialisation (for VTK version above 6)
#include <vtkAutoInit.h>
VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkRenderingFreeType);
VTK_MODULE_INIT(vtkInteractionStyle);

int main(int argc, char *argv[])
{
    // Load configuration from YAML file
    std::string config_path = (argc > 1) ? argv[1] : "../config/3_continuous_stewart_gough.yaml";
    ConfigLoader config(config_path);

    auto topology = config.getTopology();
    auto params = config.getHyperparameters();
    auto options = config.getOptions();
    auto measurements = config.getMeasurements();
    auto vis_settings = config.getVisualizationSettings();

    // Create estimator and compute state estimate
    ContinuumRobotStateEstimator state_estimator(topology, params, options);

    ContinuumRobotStateEstimator::SystemState state;
    std::vector<double> cost;
    state_estimator.computeStateEstimate(state, cost, measurements, vis_settings.verbose);

    state_estimator.printStateMean(state);

    // Visualization
    Visualizer vis(topology);
    vis.update(state, vis_settings.render_frames, vis_settings.render_covariance, vis_settings.covariance_n_std);

    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(vis.getRenderWindow());
    renderWindowInteractor->UpdateSize(vis_settings.window_width, vis_settings.window_height);
    vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
    renderWindowInteractor->SetInteractorStyle(style);
    renderWindowInteractor->Initialize();
    renderWindowInteractor->Start();

    return 0;
}
