#pragma once

#include <vector>

#include <vtkSmartPointer.h>
#include <vtkActor.h>
#include <vtkAxesActor.h>
#include <vtkPoints.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>

#include "continuum_robot_state_estimator.h"
#include "config_loader.h"

// VTK-based visualization utility for the estimator state.
//
// Typical usage:
// 1) Construct once with the same topology used by the estimator.
// 2) Call update(...) whenever a new SystemState is available.
// 3) Attach getRenderWindow() to a vtkRenderWindowInteractor.
class Visualizer
{

public:
    // Build and initialize the rendering scene (actors, camera, buffers).
    // The topology defines how many robots/couplings/frames must be created.
    // Camera position/orientation is read from vis_settings.
    Visualizer(ContinuumRobotStateEstimator::RobotTopology topology,
               const ConfigLoader::VisualizationSettings& vis_settings = {});

    // Smart pointers handle memory, so no manual cleanup is needed here.
    ~Visualizer();

    // Update all rendered objects from the current estimated system state.
    // render_frames: show/hide local coordinate frames at estimation nodes.
    // render_covariance: show/hide position covariance ellipsoids.
    // n_std: ellipsoid scale in number of standard deviations.
    void update(ContinuumRobotStateEstimator::SystemState state, bool render_frames = true, bool render_covariance = false, int n_std = 3);

    // Return the VTK render window so external code can create an interactor.
    vtkSmartPointer<vtkRenderWindow> getRenderWindow();

private:
    // Allocate actors/mappers/points once and add them to the renderer.
    // Runtime updates only move/reshape these existing objects.
    void InitScene();

    // Copy of estimator topology used to size and index all visual objects.
    ContinuumRobotStateEstimator::RobotTopology m_topology;

    // Visualization settings (camera position, etc.)
    ConfigLoader::VisualizationSettings m_vis_settings;

    // Top-level VTK rendering objects.
    vtkSmartPointer<vtkRenderWindow> mp_renWin;
    vtkSmartPointer<vtkRenderer>     mp_ren;

    // One axis actor per estimation node (plus optional end-effector frame).
    std::vector<vtkSmartPointer<vtkAxesActor>> mp_axes;

    // Robot backbone geometry:
    // - points store sampled centerline positions
    // - actors render each backbone as a tube/line object
    std::vector<vtkSmartPointer<vtkPoints>> mp_backbone_points;
    std::vector<vtkSmartPointer<vtkActor>> mp_backbone_actors;

    // Coupling link geometry between robots/end-effector.
    std::vector<vtkSmartPointer<vtkPoints>> mp_coupling_points;
    std::vector<vtkSmartPointer<vtkActor>> mp_coupling_actors;

    // Joint marker actors (cube for rigid, sphere for spherical/other).
    std::vector<vtkSmartPointer<vtkActor>> mp_joint_actors;

    // Position uncertainty actors (scaled/rotated spheres -> ellipsoids).
    std::vector<vtkSmartPointer<vtkActor>> mp_ellipsoid_actors;
};
