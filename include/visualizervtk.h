#pragma once

#include <vector>

#include <vtkSmartPointer.h>
#include <vtkActor.h>
#include <vtkAxesActor.h>
#include <vtkPoints.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>

#include "continuum_robot_state_estimator.h"

class Visualizer
{

public:

    Visualizer(ContinuumRobotStateEstimator::RobotTopology topology);
    ~Visualizer();

    void update(ContinuumRobotStateEstimator::SystemState state, bool render_frames = true, bool render_covariance = false, int n_std = 3);
	vtkSmartPointer<vtkRenderWindow> getRenderWindow();


private:

    void InitScene();
    ContinuumRobotStateEstimator::RobotTopology m_topology;
    vtkSmartPointer<vtkRenderWindow> mp_renWin;
    vtkSmartPointer<vtkRenderer>     mp_ren;




    //Stuff to visualize
    std::vector<vtkSmartPointer<vtkAxesActor>> mp_axes;

    std::vector<vtkSmartPointer<vtkPoints>> mp_backbone_points;
    std::vector<vtkSmartPointer<vtkActor>> mp_backbone_actors;

    std::vector<vtkSmartPointer<vtkPoints>> mp_coupling_points;
    std::vector<vtkSmartPointer<vtkActor>> mp_coupling_actors;

    std::vector<vtkSmartPointer<vtkActor>> mp_joint_actors;

    std::vector<vtkSmartPointer<vtkActor>> mp_ellipsoid_actors;


};


