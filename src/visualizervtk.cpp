#include "visualizervtk.h"

#include <Eigen/Eigenvalues>

#include <vtkCamera.h>
#include <vtkMatrix4x4.h>
#include <vtkSphereSource.h>
#include <vtkCubeSource.h>
#include <vtkProperty.h>
#include <vtkMapper.h>
#include <vtkPolyDataMapper.h>
#include <vtkTubeFilter.h>
/**
 * Constructor.
 * @brief Visualizer::Visualizer
 */
Visualizer::Visualizer(ContinuumRobotStateEstimator::RobotTopology topology)
{

    mp_renWin = vtkSmartPointer<vtkRenderWindow>::New();
    mp_ren = vtkSmartPointer<vtkRenderer>::New();

    mp_renWin->AddRenderer(mp_ren);

    m_topology = topology;

    InitScene();
}

Visualizer::~Visualizer()
{

}

vtkSmartPointer<vtkRenderWindow> Visualizer::getRenderWindow() {

	return mp_renWin;
}

void Visualizer::update(ContinuumRobotStateEstimator::SystemState state, bool render_frames, bool render_covariance, int n_std)
{
    //Backbones
    for(unsigned int n = 0; n < m_topology.N; n++)
    {
        //Run through all interpolation nodes

        for(unsigned int m = 0; m < state.robots[n].interpolation_nodes.size(); m++)
        {
            Eigen::Vector3d pos = state.robots[n].interpolation_nodes[m].pose.block(0,3,3,1);

            mp_backbone_points[n]->SetPoint(m,pos(0),pos(1),pos(2)); //Set the point of the line
        }
        mp_backbone_points[n]->Modified();
    }

    //Coupling Links and Joints
    for(unsigned int c = 0; c < m_topology.robot_coupling.size(); c++)
    {
        //LINKS
        ContinuumRobotStateEstimator::RobotTopology::Coupling coupling = m_topology.robot_coupling[c];

        Eigen::Matrix4d pose_node_A = state.robots[coupling.idxA].estimation_nodes[coupling.coupling_node_robot_A].pose;



        Eigen::Matrix4d pose_node_B;
        if(coupling.idxB == m_topology.N)
        {
            pose_node_B = state.end_effector.pose;
        }
        else
        {
            pose_node_B = state.robots[coupling.idxB].estimation_nodes[coupling.coupling_node_robot_B].pose;
        }

        Eigen::Matrix4d pose_coupling = pose_node_A*coupling.T_bA_c;


        //Set point 1 (node of robot A)
        Eigen::Vector3d pos1 = pose_node_A.block(0,3,3,1);

        mp_coupling_points[c]->SetPoint(0,pos1(0),pos1(1),pos1(2));

        //Set point 2 (x-offset of robot A)
        Eigen::Vector3d pos2 = pose_coupling.block(0,3,3,1);

        mp_coupling_points[c]->SetPoint(1,pos2(0),pos2(1),pos2(2));

        //Set point 3 (node of robot B)
        Eigen::Vector3d pos3 = pose_node_B.block(0,3,3,1);

        mp_coupling_points[c]->SetPoint(2,pos3(0),pos3(1),pos3(2));

        mp_coupling_points[c]->Modified();

        //JOINTS
        vtkSmartPointer<vtkMatrix4x4> joint_frame_vtk = vtkSmartPointer<vtkMatrix4x4>::New();
        for(int i = 0; i < 4; i++)
        {
            for(int j = 0; j < 4; j++)
            {
                joint_frame_vtk->SetElement(i,j,pose_coupling(i,j));
            }
        }
        mp_joint_actors[c]->SetUserMatrix(joint_frame_vtk);


    }

    //Axes (set visibility to true or false depending on input)
    if(render_frames)
    {
        //Show all axes
        for(unsigned int i = 0; i < mp_axes.size(); i++)
        {
            mp_axes[i]->SetVisibility(true);
        }

        //Update Robot frames
        int k_offset = 0;
        for(unsigned int n = 0; n < m_topology.N; n++)
        {
            for(unsigned int k = 0; k < m_topology.K[n]; k++)
            {
                Eigen::Matrix4d axes_pose = state.robots[n].estimation_nodes[k].pose;

                vtkSmartPointer<vtkMatrix4x4> axes_frame_vtk = vtkSmartPointer<vtkMatrix4x4>::New();
                for(int i = 0; i < 4; i++)
                {
                    for(int j = 0; j < 4; j++)
                    {
                        axes_frame_vtk->SetElement(i,j,axes_pose(i,j));
                    }
                }
                mp_axes[k+k_offset]->SetUserMatrix(axes_frame_vtk);
            }
            k_offset = k_offset + m_topology.K[n];
        }

        //Update EE frame
        if(m_topology.common_end_effector)
        {
            Eigen::Matrix4d ee_pose = state.end_effector.pose;

            vtkSmartPointer<vtkMatrix4x4> axes_frame_vtk = vtkSmartPointer<vtkMatrix4x4>::New();
            for(int i = 0; i < 4; i++)
            {
                for(int j = 0; j < 4; j++)
                {
                    axes_frame_vtk->SetElement(i,j,ee_pose(i,j));
                }
            }
            mp_axes.back()->SetUserMatrix(axes_frame_vtk);
        }

    }
    else
    {
        //Hide all axes
        for(unsigned int i = 0; i < mp_axes.size(); i++)
        {
            mp_axes[i]->SetVisibility(false);
        }
    }


    if(render_covariance)
    {
        Eigen::EigenSolver<Eigen::MatrixXd> solver;
        int m_offset = 0;
        for(unsigned int n = 0; n < m_topology.N; n++)
        {
            for(unsigned int m = 0; m < state.robots[n].interpolation_nodes.size(); m++)
            {

                Eigen::Vector3d pos = state.robots[n].interpolation_nodes[m].pose.block(0,3,3,1);
                Eigen::Matrix3d cov = state.robots[n].interpolation_nodes[m].position_covariance;
                solver.compute(cov);
                Eigen::MatrixXd eigen_vectors = solver.eigenvectors().real();
                Eigen::VectorXd eigen_values = solver.eigenvalues().real();
                eigen_values = n_std*eigen_values.cwiseSqrt();

                Eigen::Matrix3d R = eigen_vectors;
                Eigen::Vector3d s = eigen_values;


                //Make sure that R resembles a rotation matrix
                if((R.col(0).cross(R.col(1))).dot(R.col(2)) < 0)
                {
                    R << eigen_vectors.col(1), eigen_vectors.col(0), eigen_vectors.col(2);
                    s << eigen_values(1),
                            eigen_values(0),
                            eigen_values(2);
                }

                //Create sphere pose
                Eigen::Matrix4d sphere_pose = Eigen::Matrix4d::Identity();

                sphere_pose.block(0,0,3,3) = R;
                sphere_pose.block(0,3,3,1) = pos;

                //Check if one of the singular values is below a treshold (important to not scale the sphere's dimensions to zero)
                for(int i = 0; i < s.size(); i++)
                {
                    if(s(i) < 1e-8)
                        s(i) = 1e-8;
                }


                mp_ellipsoid_actors[m+m_offset]->SetVisibility(true);

                vtkSmartPointer<vtkMatrix4x4> sphere_frame_vtk = vtkSmartPointer<vtkMatrix4x4>::New();
                for(int i = 0; i < 4; i++)
                {
                    for(int j = 0; j < 4; j++)
                    {
                        sphere_frame_vtk->SetElement(i,j,sphere_pose(i,j));
                    }
                }
                mp_ellipsoid_actors[m+m_offset]->SetUserMatrix(sphere_frame_vtk);
                mp_ellipsoid_actors[m+m_offset]->SetScale(s(0),s(1),s(2));

            }
            m_offset = m_offset + state.robots[n].interpolation_nodes.size();
        }

        if(m_topology.common_end_effector)
        {

            Eigen::Vector3d pos = state.end_effector.pose.block(0,3,3,1);
            Eigen::Matrix3d cov = state.end_effector.position_covariance;
            solver.compute(cov);
            Eigen::MatrixXd eigen_vectors = solver.eigenvectors().real();
            Eigen::VectorXd eigen_values = solver.eigenvalues().real();
            eigen_values = n_std*eigen_values.cwiseSqrt();

            Eigen::Matrix3d R = eigen_vectors;
            Eigen::Vector3d s = eigen_values;


            //Make sure that R resembles a rotation matrix
            if((R.col(0).cross(R.col(1))).dot(R.col(2)) < 0)
            {
                R << eigen_vectors.col(1), eigen_vectors.col(0), eigen_vectors.col(2);
                s << eigen_values(1),
                        eigen_values(0),
                        eigen_values(2);
            }

            //Create sphere pose
            Eigen::Matrix4d sphere_pose = Eigen::Matrix4d::Identity();

            sphere_pose.block(0,0,3,3) = R;
            sphere_pose.block(0,3,3,1) = pos;

            //Check if one of the singular values is below a treshold (important to not scale the sphere's dimensions to zero)
            for(int i = 0; i < s.size(); i++)
            {
                if(s(i) < 1e-8)
                    s(i) = 1e-8;
            }


            mp_ellipsoid_actors.back()->SetVisibility(true);

            vtkSmartPointer<vtkMatrix4x4> sphere_frame_vtk = vtkSmartPointer<vtkMatrix4x4>::New();
            for(int i = 0; i < 4; i++)
            {
                for(int j = 0; j < 4; j++)
                {
                    sphere_frame_vtk->SetElement(i,j,sphere_pose(i,j));
                }
            }
            mp_ellipsoid_actors.back()->SetUserMatrix(sphere_frame_vtk);
            mp_ellipsoid_actors.back()->SetScale(s(0),s(1),s(2));

        }
    }
    else
    {
        //Hide all ellipsoids
        for(unsigned int i = 0; i < mp_ellipsoid_actors.size(); i++)
        {
            mp_ellipsoid_actors[i]->SetVisibility(false);
        }
    }

    //Update scene
    mp_renWin->Render();

}

void Visualizer::InitScene()
{

    //Background
    mp_ren->SetBackground(1.,1.,1.);

    //Backbones
    for(unsigned int n = 0; n < m_topology.N; n++)
    {
        //Store the points for each robot
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

        for(unsigned int m = 0; m < (m_topology.K[n]-1)*m_topology.M[n] + 1; m++)
        {
            points->InsertPoint(m,0,0,0); //Set the point of the line
        }
        mp_backbone_points.push_back(points);

        vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
        polydata->SetPoints(points);
        polydata->Allocate();

        for(unsigned int m = 0; m < (m_topology.K[n]-1)*m_topology.M[n]; m++)
        {
            vtkIdType connectivity[2];
            connectivity[0] = m;
            connectivity[1] = m+1;
            polydata->InsertNextCell(VTK_LINE,2,connectivity);
        }

        vtkSmartPointer<vtkTubeFilter> tubeFilter = vtkSmartPointer<vtkTubeFilter>::New();

        tubeFilter->SetInputData(polydata);
        tubeFilter->SetRadius(0.0005);
        tubeFilter->SetNumberOfSides(50);

        vtkSmartPointer<vtkPolyDataMapper> backboneMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        backboneMapper->SetInputConnection(tubeFilter->GetOutputPort());
        mp_backbone_actors.push_back(vtkSmartPointer<vtkActor>::New());
        mp_backbone_actors.back()->SetMapper(backboneMapper);
        mp_backbone_actors.back()->GetProperty()->SetColor(0,0,0);
        mp_backbone_actors.back()->GetProperty()->SetAmbient(0.3);
        mp_backbone_actors.back()->GetProperty()->SetDiffuse(0.5);
        mp_backbone_actors.back()->GetProperty()->SetSpecular(0.1);
        mp_ren->AddActor(mp_backbone_actors.back());
    }

    //Coupling Links and Joints

    for(unsigned int c = 0; c < m_topology.robot_coupling.size(); c++)
    {

        //LINKS


        //Store the points for each robot
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

        //Three points for each coupling position
        for(unsigned int m = 0; m < 3; m++)
        {
            points->InsertPoint(m,0,0,0); //Set the point of the line
        }
        mp_coupling_points.push_back(points);

        vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
        polydata->SetPoints(points);
        polydata->Allocate();

        for(unsigned int m = 0; m < 2; m++)
        {
            vtkIdType connectivity[2];
            connectivity[0] = m;
            connectivity[1] = m+1;
            polydata->InsertNextCell(VTK_LINE,2,connectivity);
        }

        vtkSmartPointer<vtkTubeFilter> tubeFilter = vtkSmartPointer<vtkTubeFilter>::New();

        tubeFilter->SetInputData(polydata);
        tubeFilter->SetRadius(0.0005);
        tubeFilter->SetNumberOfSides(50);

        vtkSmartPointer<vtkPolyDataMapper> backboneMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        backboneMapper->SetInputConnection(tubeFilter->GetOutputPort());
        mp_coupling_actors.push_back(vtkSmartPointer<vtkActor>::New());
        mp_coupling_actors.back()->SetMapper(backboneMapper);
        mp_coupling_actors.back()->GetProperty()->SetColor(0.75,0.75,0.75);
        mp_coupling_actors.back()->GetProperty()->SetAmbient(0.3);
        mp_coupling_actors.back()->GetProperty()->SetDiffuse(0.5);
        mp_coupling_actors.back()->GetProperty()->SetSpecular(0.1);
        mp_ren->AddActor(mp_coupling_actors.back());

        //JOINTS

        ContinuumRobotStateEstimator::RobotTopology::Coupling coupling = m_topology.robot_coupling[c];

        Eigen::MatrixXi mask = coupling.mask;

        vtkSmartPointer<vtkPolyDataMapper> jointMapper = vtkSmartPointer<vtkPolyDataMapper>::New();

        if(mask(0,0) == 1 && mask(1,0) == 1 && mask(2,0) == 1 && mask(3,0) == 1 && mask(4,0) == 1 && mask(5,0) == 1) //Rigid joint
        {
            vtkSmartPointer<vtkCubeSource> cubeSource = vtkSmartPointer<vtkCubeSource>::New();
            cubeSource->SetXLength(0.003);
            cubeSource->SetYLength(0.003);
            cubeSource->SetZLength(0.003);
            cubeSource->Update();

            //Connect to Mapper
            jointMapper->SetInputConnection(cubeSource->GetOutputPort());

        }
        else//Spherical Joint and anything else
        {
            vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
            sphereSource->SetRadius(0.0015);
            sphereSource->SetPhiResolution(100);
            sphereSource->SetThetaResolution(100);
            sphereSource->Update();

            //Connect to Mapper
            jointMapper->SetInputConnection(sphereSource->GetOutputPort());
        }
        mp_joint_actors.push_back(vtkSmartPointer<vtkActor>::New());
        mp_joint_actors.back()->SetMapper(jointMapper);
        mp_joint_actors.back()->GetProperty()->SetColor(0.75,0.75,0.75);
        mp_joint_actors.back()->GetProperty()->SetAmbient(0.3);
        mp_joint_actors.back()->GetProperty()->SetDiffuse(0.5);
        mp_joint_actors.back()->GetProperty()->SetSpecular(0.1);
        mp_ren->AddActor(mp_joint_actors.back());

    }


    //Axes


    //Robot frames
    for(unsigned int n = 0; n < m_topology.N; n++)
    {
        for(unsigned int k = 0; k < m_topology.K[n]; k++)
        {
            vtkSmartPointer<vtkAxesActor> robot_axes = vtkSmartPointer<vtkAxesActor>::New();
            robot_axes->SetXAxisLabelText("");
            robot_axes->SetYAxisLabelText("");
            robot_axes->SetZAxisLabelText("");
            robot_axes->SetShaftTypeToCylinder();
            robot_axes->SetCylinderRadius(0.025);
            robot_axes->SetTotalLength(0.01,0.01,0.01);
            robot_axes->SetVisibility(true);
            mp_axes.push_back(robot_axes);
            mp_ren->AddActor(robot_axes);
        }
    }

    //End-effector frame (if there)
    if(m_topology.common_end_effector)
    {
        vtkSmartPointer<vtkAxesActor> ee_frame = vtkSmartPointer<vtkAxesActor>::New();
        ee_frame->SetXAxisLabelText("");
        ee_frame->SetYAxisLabelText("");
        ee_frame->SetZAxisLabelText("");
        ee_frame->SetShaftTypeToCylinder();
        ee_frame->SetCylinderRadius(0.025);
        ee_frame->SetTotalLength(0.01,0.01,0.01);
        ee_frame->SetVisibility(true);
        mp_axes.push_back(ee_frame);
        mp_ren->AddActor(ee_frame);
    }



    //Covariance ellipsoids
    for(unsigned int n = 0; n < m_topology.N; n++)
    {
        // First ellipsoid at root of robot
        vtkSmartPointer<vtkSphereSource> ellipsoid_source = vtkSmartPointer<vtkSphereSource>::New();
        ellipsoid_source->SetCenter(0.0, 0.0, 0.0);
        ellipsoid_source->SetRadius(1);
        ellipsoid_source->SetPhiResolution(100);
        ellipsoid_source->SetThetaResolution(100);

        vtkSmartPointer<vtkPolyDataMapper> ellipsoid_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        ellipsoid_mapper->SetInputConnection(ellipsoid_source->GetOutputPort());

        vtkSmartPointer<vtkActor> ellipsoid_actor_init = vtkSmartPointer<vtkActor>::New();
        ellipsoid_actor_init->SetMapper(ellipsoid_mapper);
        ellipsoid_actor_init->GetProperty()->SetColor(0,0,1);
        ellipsoid_actor_init->GetProperty()->SetOpacity(0.1);
        ellipsoid_actor_init->GetProperty()->SetAmbient(0.3);
        ellipsoid_actor_init->GetProperty()->SetDiffuse(0.5);
        ellipsoid_actor_init->GetProperty()->SetSpecular(0.1);
        ellipsoid_actor_init->SetVisibility(false);

        mp_ellipsoid_actors.push_back(ellipsoid_actor_init);
        mp_ren->AddActor(ellipsoid_actor_init);


        //Remaining ellipsoids
        for(unsigned int k = 0; k < m_topology.K[n]; k++)
        {
            for(unsigned int m = 0; m < m_topology.M[n]; m++)
            {

                vtkSmartPointer<vtkActor> ellipsoid_actor = vtkSmartPointer<vtkActor>::New();
                ellipsoid_actor->SetMapper(ellipsoid_mapper);
                ellipsoid_actor->GetProperty()->SetColor(0,0,1);
                ellipsoid_actor->GetProperty()->SetOpacity(0.1);
                ellipsoid_actor->GetProperty()->SetAmbient(0.3);
                ellipsoid_actor->GetProperty()->SetDiffuse(0.5);
                ellipsoid_actor->GetProperty()->SetSpecular(0.1);
                ellipsoid_actor->SetVisibility(false);

                mp_ellipsoid_actors.push_back(ellipsoid_actor);
                mp_ren->AddActor(ellipsoid_actor);
            }
        }
    }

    if(m_topology.common_end_effector)
    {
        vtkSmartPointer<vtkSphereSource> ellipsoid_source = vtkSmartPointer<vtkSphereSource>::New();
        ellipsoid_source->SetCenter(0.0, 0.0, 0.0);
        ellipsoid_source->SetRadius(1);
        ellipsoid_source->SetPhiResolution(100);
        ellipsoid_source->SetThetaResolution(100);

        vtkSmartPointer<vtkPolyDataMapper> ellipsoid_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        ellipsoid_mapper->SetInputConnection(ellipsoid_source->GetOutputPort());

        vtkSmartPointer<vtkActor> ellipsoid_actor_init = vtkSmartPointer<vtkActor>::New();
        ellipsoid_actor_init->SetMapper(ellipsoid_mapper);
        ellipsoid_actor_init->GetProperty()->SetColor(0,0,1);
        ellipsoid_actor_init->GetProperty()->SetOpacity(0.1);
        ellipsoid_actor_init->GetProperty()->SetAmbient(0.3);
        ellipsoid_actor_init->GetProperty()->SetDiffuse(0.5);
        ellipsoid_actor_init->GetProperty()->SetSpecular(0.1);
        ellipsoid_actor_init->SetVisibility(false);

        mp_ellipsoid_actors.push_back(ellipsoid_actor_init);
        mp_ren->AddActor(ellipsoid_actor_init);
    }




    //Camera
    mp_ren->GetActiveCamera()->SetPosition(0.2,0.3,0.5);
    mp_ren->GetActiveCamera()->SetFocalPoint(0.1,0,0);
    mp_ren->GetActiveCamera()->SetViewUp(1,0,0);


    //Update scene
    mp_renWin->Render();











}



