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

// Interactive context shared with keyboard handler
struct InteractiveContext {
    ContinuumRobotStateEstimator* estimator;
    std::vector<ContinuumRobotStateEstimator::SensorMeasurement> measurements;
    std::vector<ContinuumRobotStateEstimator::ControlInput> inputs;
    Visualizer* vis;
    ConfigLoader::VisualizationSettings vis_settings;
    double scale = 1.0;        // multiplier applied to all control input values
    double scale_step = 0.1;   // increment per key press
};

// Custom interactor style that intercepts arrow keys to interactively
// scale the control inputs and observe their effect on the estimated shape.
//
// Arrow Up/Down:   Increase/decrease the control input scale factor.
//                  At scale=1 the original inputs are used; at scale=0 the
//                  estimator runs with no control inputs (pure WNOA prior).
//                  This is useful to visualize how much the control inputs
//                  influence the estimated robot shape compared to measurements alone.
// Arrow Left/Right: Halve/double the step size for finer or coarser adjustment.
// 'r':             Reset scale to 1.0 (original control inputs).
// '0':             Set scale to 0.0 (no control inputs).
class KeyInteractorStyle : public vtkInteractorStyleTrackballCamera
{
public:
    static KeyInteractorStyle* New() { return new KeyInteractorStyle; }
    vtkTypeMacro(KeyInteractorStyle, vtkInteractorStyleTrackballCamera);

    InteractiveContext* ctx = nullptr;

    void OnKeyPress() override
    {
        std::string key = this->Interactor->GetKeySym();

        bool changed = false;
        if (key == "Up") {
            ctx->scale += ctx->scale_step;
            changed = true;
        } else if (key == "Down") {
            ctx->scale -= ctx->scale_step;
            changed = true;
        } else if (key == "Right") {
            ctx->scale_step *= 2.0;
            std::cout << "[Step size: " << ctx->scale_step << "]" << std::endl;
        } else if (key == "Left") {
            ctx->scale_step *= 0.5;
            std::cout << "[Step size: " << ctx->scale_step << "]" << std::endl;
        } else if (key == "r") {
            ctx->scale = 1.0;
            changed = true;
        } else if (key == "0") {
            ctx->scale = 0.0;
            changed = true;
        }

        if (changed) {
            std::cout << "Scale: " << ctx->scale << std::endl;

            // Scale all control input values
            auto scaled_inputs = ctx->inputs;
            for (auto& ci : scaled_inputs)
                for (auto& v : ci.values)
                    v *= ctx->scale;

            // Re-run estimator with scaled inputs and update the visualization
            ContinuumRobotStateEstimator::SystemState state;
            std::vector<double> cost;
            ctx->estimator->computeStateEstimate(state, cost, ctx->measurements, scaled_inputs, false);

            ctx->vis->update(state, ctx->vis_settings.render_frames,
                             ctx->vis_settings.render_covariance, ctx->vis_settings.covariance_n_std);
            this->Interactor->GetRenderWindow()->Render();
        }

        // Forward to parent for camera controls
        vtkInteractorStyleTrackballCamera::OnKeyPress();
    }
};

int main(int argc, char *argv[])
{
    // Load configuration from YAML file
    std::string config_path = (argc > 1) ? argv[1] : "../config/1_continuum_robot.yaml";
    ConfigLoader config(config_path);

    auto topology = config.getTopology();
    auto params = config.getHyperparameters();
    auto options = config.getOptions();
    auto measurements = config.getMeasurements();
    auto inputs = config.getControlInputs();
    auto vis_settings = config.getVisualizationSettings();

    // Create estimator and compute initial state estimate
    ContinuumRobotStateEstimator state_estimator(topology, params, options);

    ContinuumRobotStateEstimator::SystemState state;
    std::vector<double> cost;
    state_estimator.computeStateEstimate(state, cost, measurements, inputs, vis_settings.verbose);

    state_estimator.printStateMean(state);

    // Visualization
    Visualizer vis(topology);
    vis.update(state, vis_settings.render_frames, vis_settings.render_covariance, vis_settings.covariance_n_std);

    // Setup interactive context
    InteractiveContext ctx;
    ctx.estimator = &state_estimator;
    ctx.measurements = measurements;
    ctx.inputs = inputs;
    ctx.vis = &vis;
    ctx.vis_settings = vis_settings;

    std::cout << "\n--- Interactive Control Input Viewer ---" << std::endl;
    std::cout << "  Up/Down    : increase/decrease input scale" << std::endl;
    std::cout << "  Left/Right : halve/double step size" << std::endl;
    std::cout << "  r          : reset scale to 1.0" << std::endl;
    std::cout << "  0          : set scale to 0.0 (no inputs)" << std::endl;
    std::cout << "  Mouse      : rotate/pan/zoom as usual" << std::endl;
    std::cout << "------------------------------------------\n" << std::endl;

    vtkSmartPointer<KeyInteractorStyle> style = vtkSmartPointer<KeyInteractorStyle>::New();
    style->ctx = &ctx;

    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(vis.getRenderWindow());
    renderWindowInteractor->UpdateSize(vis_settings.window_width, vis_settings.window_height);
    renderWindowInteractor->SetInteractorStyle(style);
    renderWindowInteractor->Initialize();
    renderWindowInteractor->Start();

    return 0;
}
