#include "config_loader.h"
#include "visualizervtk.h"
#include "continuum_robot_state_estimator.h"

#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRendererCollection.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>

#include <iomanip>
#include <sstream>

// VTK Factory initialisation (for VTK version above 6)
#include <vtkAutoInit.h>
VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkRenderingFreeType);
VTK_MODULE_INIT(vtkInteractionStyle);

// All labels exactly 20 chars for alignment
static const char* STRAIN_LABELS[6] = {
    "v_x [elongation]    ", "v_y [shear y]       ", "v_z [shear z]       ",
    "w_x [torsion]       ", "w_y [bending y]     ", "w_z [bending z]     "
};
static const char* FORCE_LABELS[6] = {
    "f_x [axial force]   ", "f_y [shear force y] ", "f_z [shear force z] ",
    "m_x [torque]        ", "m_y [moment y]      ", "m_z [moment z]      "
};

struct InteractiveContext {
    ContinuumRobotStateEstimator* estimator;
    std::vector<ContinuumRobotStateEstimator::SensorMeasurement> measurements;
    Visualizer* vis;
    ConfigLoader::VisualizationSettings vis_settings;

    int num_segments;
    Eigen::Matrix<double,12,1> input_values = Eigen::Matrix<double,12,1>::Zero();

    int selected_component = 4;       // start on w_y
    double step = 1.0;

    vtkSmartPointer<vtkTextActor> hud_main;       // gray text (non-selected rows)
    vtkSmartPointer<vtkTextActor> hud_highlight;  // red text (selected row only)
};

static std::string fmtVal(double v)
{
    std::ostringstream ss;
    ss << std::showpos << std::fixed << std::setprecision(1) << std::setw(7) << v;
    return ss.str();
}

// Build a single row string for both columns
static std::string buildRow(const InteractiveContext& ctx, int i, bool show_marker)
{
    std::ostringstream ss;
    bool sel_s = show_marker && (ctx.selected_component == i);
    bool sel_f = show_marker && (ctx.selected_component == i + 6);

    ss << (sel_s ? ">> " : "   ")
       << STRAIN_LABELS[i] << fmtVal(ctx.input_values(i))
       << (sel_s ? " <<" : "   ")
       << "   "
       << (sel_f ? ">> " : "   ")
       << FORCE_LABELS[i] << fmtVal(ctx.input_values(i + 6))
       << (sel_f ? " <<" : "   ");
    return ss.str();
}

// Number of lines above the data rows (header + separator)
static const int HEADER_LINES = 2;
// Number of lines below data rows (separator + step + blank + 6 keybinds)
static const int FOOTER_LINES = 9;
// Total HUD lines = 2 + 6 + 9 = 17
static const int HUD_TOTAL_LINES = 17;

static std::string buildMainHUD(const InteractiveContext& ctx)
{
    std::ostringstream ss;

    ss << "   Strain (eps_in)                      Forces (f_in)\n";
    ss << "   ---------------------------------------------------------------\n";

    for (int i = 0; i < 6; ++i)
        ss << buildRow(ctx, i, false) << "\n";

    ss << "   ---------------------------------------------------------------\n";
    ss << "   Left/Right changes value by: " << ctx.step << "\n\n";
    ss << "   Up/Down      select component\n";
    ss << "   Left/Right   adjust value\n";
    ss << "   Tab          switch column\n";
    ss << "   [ / ]        halve/double step size\n";
    ss << "   0            zero selected component\n";
    ss << "   r            reset all to zero\n";

    return ss.str();
}

static std::string buildHighlightHUD(const InteractiveContext& ctx)
{
    std::ostringstream ss;

    int sel_row = (ctx.selected_component >= 6)
                  ? ctx.selected_component - 6
                  : ctx.selected_component;

    for (int line = 0; line < HUD_TOTAL_LINES; ++line) {
        int data_line = line - 2;  // data rows start after 2 header lines
        if (data_line >= 0 && data_line < 6 && data_line == sel_row) {
            ss << buildRow(ctx, sel_row, true) << "\n";
        } else {
            ss << "\n";
        }
    }

    return ss.str();
}

static void updateDisplay(InteractiveContext& ctx)
{
    ctx.hud_main->SetInput(buildMainHUD(ctx).c_str());
    ctx.hud_highlight->SetInput(buildHighlightHUD(ctx).c_str());
}

static void printStartupMessage(const InteractiveContext& ctx)
{
    std::cout << "\n\033[1;36m  Live Control Input Viewer\033[0m\n";
    std::cout << "  " << ctx.num_segments << " segments\n\n";
    std::cout << "  All controls are in the VTK window.\n";
    std::cout << "  Keep the VTK window focused to use keyboard.\n\n";
    std::cout << "\033[0;90m";
    std::cout << "  Up/Down      select component\n";
    std::cout << "  Left/Right   adjust value\n";
    std::cout << "  Tab          switch column\n";
    std::cout << "  [ / ]        halve/double step\n";
    std::cout << "  0            zero component\n";
    std::cout << "  r            reset all\n";
    std::cout << "\033[0m" << std::endl;
}

static void rerunEstimator(InteractiveContext& ctx)
{
    std::vector<ContinuumRobotStateEstimator::ControlInput> inputs;
    for (int seg = 0; seg < ctx.num_segments; ++seg) {
        ContinuumRobotStateEstimator::ControlInput ci;
        ci.type = ContinuumRobotStateEstimator::ControlInput::Constant;
        ci.idx_robot = 0;
        ci.idx_segment = seg;
        ci.values.push_back(ctx.input_values);
        inputs.push_back(ci);
    }

    ContinuumRobotStateEstimator::SystemState state;
    std::vector<double> cost;
    ctx.estimator->computeStateEstimate(state, cost, ctx.measurements, inputs, false);

    ctx.vis->update(state, ctx.vis_settings.render_frames,
                    ctx.vis_settings.render_covariance, ctx.vis_settings.covariance_n_std);
}

class LiveInputInteractorStyle : public vtkInteractorStyleTrackballCamera
{
public:
    static LiveInputInteractorStyle* New() { return new LiveInputInteractorStyle; }
    vtkTypeMacro(LiveInputInteractorStyle, vtkInteractorStyleTrackballCamera);

    InteractiveContext* ctx = nullptr;

    void OnKeyPress() override
    {
        std::string key = this->Interactor->GetKeySym();
        bool changed = false;

        if (key == "Up") {
            int col = (ctx->selected_component >= 6) ? 6 : 0;
            int row = ctx->selected_component - col;
            row = (row - 1 + 6) % 6;
            ctx->selected_component = col + row;
        } else if (key == "Down") {
            int col = (ctx->selected_component >= 6) ? 6 : 0;
            int row = ctx->selected_component - col;
            row = (row + 1) % 6;
            ctx->selected_component = col + row;
        } else if (key == "Right") {
            ctx->input_values(ctx->selected_component) += ctx->step;
            changed = true;
        } else if (key == "Left") {
            ctx->input_values(ctx->selected_component) -= ctx->step;
            changed = true;
        } else if (key == "Tab") {
            ctx->selected_component = (ctx->selected_component + 6) % 12;
        } else if (key == "bracketleft") {
            ctx->step *= 0.5;
        } else if (key == "bracketright") {
            ctx->step *= 2.0;
        } else if (key == "0") {
            ctx->input_values(ctx->selected_component) = 0.0;
            changed = true;
        } else if (key == "r") {
            ctx->input_values.setZero();
            changed = true;
        } else {
            vtkInteractorStyleTrackballCamera::OnKeyPress();
            return;
        }

        updateDisplay(*ctx);

        if (changed)
            rerunEstimator(*ctx);

        this->Interactor->GetRenderWindow()->Render();
        vtkInteractorStyleTrackballCamera::OnKeyPress();
    }
};

int main(int argc, char *argv[])
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <config.yaml>\n";
        return 1;
    }

    std::string config_path = argv[1];
    ConfigLoader config(config_path);

    auto topology = config.getTopology();
    auto params = config.getHyperparameters();
    auto options = config.getOptions();
    auto measurements = config.getMeasurements();
    auto control_inputs = config.getControlInputs();
    auto vis_settings = config.getVisualizationSettings();

    ContinuumRobotStateEstimator state_estimator(topology, params, options);

    ContinuumRobotStateEstimator::SystemState state;
    std::vector<double> cost;
    state_estimator.computeStateEstimate(state, cost, measurements, control_inputs, vis_settings.verbose);
    state_estimator.printStateMean(state);

    Visualizer vis(topology, vis_settings);
    vis.update(state, vis_settings.render_frames, vis_settings.render_covariance, vis_settings.covariance_n_std);

    // Get renderer
    vtkRendererCollection* renderers = vis.getRenderWindow()->GetRenderers();
    renderers->InitTraversal();
    vtkRenderer* renderer = renderers->GetNextItem();

    // Main HUD (dark gray -- non-selected rows + controls)
    vtkSmartPointer<vtkTextActor> hud_main = vtkSmartPointer<vtkTextActor>::New();
    hud_main->GetTextProperty()->SetFontFamilyToCourier();
    hud_main->GetTextProperty()->SetFontSize(24);
    hud_main->GetTextProperty()->SetBold(true);
    hud_main->GetTextProperty()->SetColor(0.2, 0.2, 0.2);
    hud_main->GetTextProperty()->SetOpacity(0.85);
    hud_main->SetPosition(20, 20);
    renderer->AddViewProp(hud_main);

    // Highlight HUD (red -- selected row only, overlaid at same position)
    vtkSmartPointer<vtkTextActor> hud_highlight = vtkSmartPointer<vtkTextActor>::New();
    hud_highlight->GetTextProperty()->SetFontFamilyToCourier();
    hud_highlight->GetTextProperty()->SetFontSize(24);
    hud_highlight->GetTextProperty()->SetBold(true);
    hud_highlight->GetTextProperty()->SetColor(0.85, 0.15, 0.15);
    hud_highlight->GetTextProperty()->SetOpacity(1.0);
    hud_highlight->SetPosition(20, 20);
    renderer->AddViewProp(hud_highlight);

    // Setup context
    InteractiveContext ctx;
    ctx.estimator = &state_estimator;
    ctx.measurements = measurements;
    ctx.vis = &vis;
    ctx.vis_settings = vis_settings;
    ctx.num_segments = topology.K.at(0) - 1;
    ctx.hud_main = hud_main;
    ctx.hud_highlight = hud_highlight;

    // Initialize control input values from YAML (use first Constant input if available)
    for (const auto& ci : control_inputs) {
        if (ci.type == ContinuumRobotStateEstimator::ControlInput::Constant && !ci.values.empty()) {
            ctx.input_values = ci.values[0];
            break;
        }
    }

    printStartupMessage(ctx);
    updateDisplay(ctx);

    vtkSmartPointer<LiveInputInteractorStyle> style = vtkSmartPointer<LiveInputInteractorStyle>::New();
    style->ctx = &ctx;

    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(vis.getRenderWindow());
    renderWindowInteractor->UpdateSize(vis_settings.window_width, vis_settings.window_height);
    renderWindowInteractor->SetInteractorStyle(style);
    renderWindowInteractor->Initialize();
    renderWindowInteractor->Start();

    return 0;
}
