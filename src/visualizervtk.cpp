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
#include <vtkAssembly.h>
#include <vtkLine.h>
#include <vtkLineSource.h>
#include <vtkTransform.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkDelaunay3D.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkGeometryFilter.h>
#include <vtkCleanPolyData.h>
#include <vtkTriangleFilter.h>
#include <vtkSignedDistance.h>
#include <vtkPCANormalEstimation.h>
#include <vtkExtractSurface.h>
#include <vtkPolyDataNormals.h>
#include <vtkPolyDataMapper.h>
#include <vtkPointData.h>
#include <vtkPolygon.h>



#include <utilities.h>
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
    double transparency = 1; // Set transparency level (0.0 = fully transparent, 1.0 = fully opaque)

    vtkSmartPointer<vtkAssembly> assembly = vtkSmartPointer<vtkAssembly>::New();

    for (int n = 0; n < state.robots.size(); ++n) {
        for (int i = 0; i < state.robots.at(n).estimation_nodes.size(); ++i) {
            Eigen::Matrix4d transform = state.robots.at(n).estimation_nodes.at(i).pose.inverse();
            // Invert transform exploiting the fact that its a rigid body transform
            Eigen::Matrix3d rotation = transform.block<3, 3>(0, 0);
            Eigen::Matrix3d rotation_transpose = rotation.transpose();
            transform.block<3, 3>(0, 0) = rotation_transpose;
            transform.block<3, 1>(0, 3) = -rotation_transpose * transform.block<3, 1>(0, 3);

            vtkSmartPointer<vtkMatrix4x4> vtk_transform = vtkSmartPointer<vtkMatrix4x4>::New();
            for (int r = 0; r < 4; ++r) {
                for (int c = 0; c < 4; ++c) {
                    vtk_transform->SetElement(r, c, transform(r, c));
                }
            }

            if(i == state.robots.at(n).estimation_nodes.size() - 1)
            {
                vtkSmartPointer<vtkSphereSource> half_sphere = vtkSmartPointer<vtkSphereSource>::New();
                half_sphere->SetRadius(0.005);
                half_sphere->SetPhiResolution(50);
                half_sphere->SetThetaResolution(50);
                half_sphere->SetStartTheta(0.0);
                half_sphere->SetEndTheta(180.0); // Only generate half sphere
                half_sphere->SetCenter(0.0, 0.0, 0.0);

                vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
                mapper->SetInputConnection(half_sphere->GetOutputPort());

                vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
                actor->SetMapper(mapper);
                vtkSmartPointer<vtkTransform> rotation_transform = vtkSmartPointer<vtkTransform>::New();
                rotation_transform->RotateZ(-90); // Rotate half sphere to align with x-axis
                vtkSmartPointer<vtkTransform> final_transform = vtkSmartPointer<vtkTransform>::New();
                final_transform->Concatenate(vtk_transform);
                final_transform->Concatenate(rotation_transform);
                actor->SetUserTransform(final_transform);

                vtkSmartPointer<vtkProperty> property = vtkSmartPointer<vtkProperty>::New();
                property->SetColor(0.5, 0.5, 0.5); // Grey color for plastic
                property->SetColor(0.3,0.3,0.3); // Grey color for plastic
                property->SetSpecular(0.1); // Lower specular to make it look less shiny
                property->SetSpecularPower(5); // Lower specular power for a more matte finish
                property->SetInterpolationToPhong(); // Use Phong interpolation for a more realistic look
                property->SetRepresentationToSurface(); // Ensure the representation is set to surface
                property->SetOpacity(transparency); // Apply transparency
                actor->SetProperty(property);

                assembly->AddPart(actor);
            }

            // Draw connecting line to the next state
            if (i < state.robots.at(n).estimation_nodes.size() - 1) {
                Eigen::Vector3d start_point = transform.block<3, 1>(0, 3);
                Eigen::Matrix4d next_transform = state.robots.at(n).estimation_nodes.at(i + 1).pose.inverse();
                Eigen::Matrix3d next_rotation = next_transform.block<3, 3>(0, 0);
                Eigen::Matrix3d next_rotation_transpose = next_rotation.transpose();
                next_transform.block<3, 3>(0, 0) = next_rotation_transpose;
                next_transform.block<3, 1>(0, 3) = -next_rotation_transpose * next_transform.block<3, 1>(0, 3);
                Eigen::Vector3d end_point = next_transform.block<3, 1>(0, 3);
                // get start and end strain from the robot_state (last 6 columns of robot_state row)
                Eigen::VectorXd strain = state.robots.at(n).estimation_nodes.at(i).strain;
                Eigen::VectorXd next_strain = state.robots.at(n).estimation_nodes.at(i + 1).strain;

                Eigen::Vector3d previous_point = start_point;
                for (int j = 1; j <= 20; ++j) {
                    double a = j / 20.0;
                    double s_interval = 0.03;

                    // Operate in the lie algebra
                    Eigen::MatrixXd xi_k1 = tran_to_vec(next_transform.inverse() * transform);

                    Eigen::MatrixXd Jinv = vec_to_jac_inverse(xi_k1);
                    Eigen::Matrix<double, 6, 1> xi_k_dot = -1 * strain;
                    Eigen::Matrix<double, 6, 1> xi_k1_dot = -1 * Jinv * next_strain;

                    Eigen::MatrixXd xi_tau = (a * a * a - 2 * a * a + a) * s_interval * xi_k_dot + (3 * a * a - 2 * a * a * a) * xi_k1 + (a * a * a - a * a) * s_interval * xi_k1_dot;
                    Eigen::MatrixXd xi_dot_tau = ((3 * a * a - 4 * a + 1) * s_interval * xi_k_dot + (6 * a - 6 * a * a) * xi_k1 + (3 * a * a - 2 * a) * s_interval * xi_k1_dot) / s_interval;

                    Eigen::Matrix4d T_diff_step = vec_to_tran(xi_tau);

                    Eigen::Matrix4d T_cur = transform * T_diff_step.inverse();

                    Eigen::Vector3d interpolated_point = T_cur.block<3, 1>(0, 3);

                    vtkSmartPointer<vtkLineSource> line_source = vtkSmartPointer<vtkLineSource>::New();
                    line_source->SetPoint1(previous_point(0), previous_point(1), previous_point(2));
                    line_source->SetPoint2(interpolated_point(0), interpolated_point(1), interpolated_point(2));

                    vtkSmartPointer<vtkPolyDataMapper> line_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
                    line_mapper->SetInputConnection(line_source->GetOutputPort());

                    vtkSmartPointer<vtkTubeFilter> tube_filter = vtkSmartPointer<vtkTubeFilter>::New();
                    tube_filter->SetInputConnection(line_source->GetOutputPort());
                    tube_filter->SetRadius(0.00075); // Set the radius of the tube
                    tube_filter->SetRadius(0.005); // Set the radius of the tube
                    tube_filter->SetNumberOfSides(50); // Set the number of sides for the tube to make it smooth
                    tube_filter->Update();

                    vtkSmartPointer<vtkPolyDataMapper> tube_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
                    tube_mapper->SetInputConnection(tube_filter->GetOutputPort());

                    vtkSmartPointer<vtkActor> tube_actor = vtkSmartPointer<vtkActor>::New();
                    tube_actor->SetMapper(tube_mapper);
                    tube_actor->GetProperty()->SetColor(0.8, 0.8, 0.8); // Black color for the tube
                    if(j <= 3)
                    {
                        tube_actor->GetProperty()->SetColor(0.3, 0.3, 0.3); // Black color for the tube
                    }
                    tube_actor->GetProperty()->SetSpecular(1.0); // Set specular to make it metallic
                    tube_actor->GetProperty()->SetSpecularPower(50); // Set specular power for shininess
                    tube_actor->GetProperty()->SetOpacity(transparency); // Apply transparency

                    assembly->AddPart(tube_actor);

                    previous_point = interpolated_point;
                }

                // for (int k = 0; k < 4; ++k) {
                //     double angle = k * M_PI / 2.0;
                //     Eigen::Vector3d offset(0.0, cos(angle) * 0.007, sin(angle) * 0.007);
                //     Eigen::Vector3d start_offset_point = start_point + rotation_transpose * offset;
                //     Eigen::Vector3d end_offset_point = end_point + next_rotation_transpose * offset;

                //     vtkSmartPointer<vtkLineSource> line_source = vtkSmartPointer<vtkLineSource>::New();
                //     line_source->SetPoint1(start_offset_point(0), start_offset_point(1), start_offset_point(2));
                //     line_source->SetPoint2(end_offset_point(0), end_offset_point(1), end_offset_point(2));

                //     vtkSmartPointer<vtkPolyDataMapper> line_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
                //     line_mapper->SetInputConnection(line_source->GetOutputPort());

                //     vtkSmartPointer<vtkTubeFilter> tube_filter = vtkSmartPointer<vtkTubeFilter>::New();
                //     tube_filter->SetInputConnection(line_source->GetOutputPort());
                //     tube_filter->SetRadius(0.00015); // Set the radius of the tube
                //     tube_filter->SetNumberOfSides(50); // Set the number of sides for the tube to make it smooth
                //     tube_filter->Update();

                //     vtkSmartPointer<vtkPolyDataMapper> tube_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
                //     tube_mapper->SetInputConnection(tube_filter->GetOutputPort());

                //     vtkSmartPointer<vtkActor> tube_actor = vtkSmartPointer<vtkActor>::New();
                //     tube_actor->SetMapper(tube_mapper);
                //     tube_actor->GetProperty()->SetSpecular(0.0); // Set specular to zero to make it non-shiny
                //     tube_actor->GetProperty()->SetDiffuse(0.8); // Increase diffuse to make it look more like a rope
                //     tube_actor->GetProperty()->SetAmbient(0.2); // Set ambient to give it a softer look
                //     tube_actor->GetProperty()->SetColor(0.5, 0.2, 0.1); // Darker brown color with more red for the rope
                //     tube_actor->GetProperty()->SetOpacity(transparency); // Apply transparency

                //     assembly->AddPart(tube_actor);
                // }
            }
        }
    }
    
    mp_ren->AddActor(assembly);

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
            for(unsigned int m = 0; m < state.robots[n].estimation_nodes.size(); m = m+2)
            {

                Eigen::Vector3d pos = state.robots[n].estimation_nodes[m].pose.block(0,3,3,1);
                Eigen::Matrix3d cov = state.robots[n].estimation_nodes[m].position_covariance;
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
            m_offset = m_offset + state.robots[n].estimation_nodes.size();
        }
        
        vtkSmartPointer<vtkPoints> all_points = vtkSmartPointer<vtkPoints>::New();
        for(unsigned int n = 0; n < m_topology.N; n++)
        {
            vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
            vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
            vtkSmartPointer<vtkCellArray> polys = vtkSmartPointer<vtkCellArray>::New();
            std::vector<vtkIdType> previous_ring_ids;
            for(unsigned int m = 0; m < state.robots[n].interpolation_nodes.size(); m = m+1)
            {
                Eigen::Vector3d pos = state.robots[n].interpolation_nodes[m].pose.block(0,3,3,1);
                Eigen::Matrix3d rot = state.robots[n].interpolation_nodes[m].pose.block(0,0,3,3);
                Eigen::Matrix3d cov = state.robots[n].interpolation_nodes[m].position_covariance;
                solver.compute(cov);
                if (solver.info() != Eigen::Success) {
                    continue; // Skip this covariance if the solver fails
                }
                Eigen::MatrixXd eigen_vectors = solver.eigenvectors().real();
                Eigen::VectorXd eigen_values = solver.eigenvalues().real();
                if ((eigen_values.array() < 0).any()) {
                    continue; // Skip this covariance if any eigenvalue is negative
                }
                eigen_values = n_std * eigen_values.cwiseSqrt();

                Eigen::Matrix3d R = eigen_vectors;
                Eigen::Vector3d s = eigen_values;

                // Make sure that R resembles a rotation matrix
                if ((R.col(0).cross(R.col(1))).dot(R.col(2)) < 0)
                {
                    R << eigen_vectors.col(1), eigen_vectors.col(0), eigen_vectors.col(2);
                    s << eigen_values(1), eigen_values(0), eigen_values(2);
                }

                // Identify the column of R that has the biggest angle to the first column of rot
                int max_col = 0;
                double max_dot = std::abs(R.col(0).dot(rot.col(0)));
                double max_dot_sign = R.col(0).dot(rot.col(0));
                for (int col = 1; col < 3; ++col) {
                    double dot = std::abs(R.col(col).dot(rot.col(0)));
                    if (dot > max_dot) {
                        max_dot = dot;
                        max_dot_sign = R.col(col).dot(rot.col(0));
                        max_col = col;
                    }
                }

                int num_points_per_circle = 100; // Number of points per circle

                if(m == state.robots[n].interpolation_nodes.size() - 1)
                {
                    // Sample points on the surface of the ellipsoid
                    int num_latitude = 20;
                    int num_longitude = num_points_per_circle;

                    int step = 1;
                    int goal = num_latitude;

                    //if(max_dot_sign < 0)
                    //{
                    //    step = -1;
                    //    goal = 0;
                    //}

                    for (int i = num_latitude/2.0; i <= goal; i = i + step)
                    {
                        double phi = vtkMath::Pi() * i / num_latitude;
                        std::vector<vtkIdType> ring_ids;
                        for (int j = 0; j < num_longitude; ++j)
                        {
                            double theta = 2.0 * vtkMath::Pi() * j / num_longitude;

                            Eigen::Vector3d point;
                            if(max_col == 2)
                            {
                                double x = s(0) * sin(phi) * cos(theta);
                                double y = s(1) * sin(phi) * sin(theta);
                                double z = s(2) * cos(phi);
                                point = R * Eigen::Vector3d(x, y, z) + pos;
                            }
                            else if(max_col == 1)
                            {
                                double x = s(0) * sin(phi) * cos(theta);
                                double y = s(1) * cos(phi);
                                double z = s(2) * sin(phi) * sin(theta);
                                point = R * Eigen::Vector3d(x, y, z) + pos;
                            }
                            else if(max_col == 0)
                            {
                                double x = s(0) * cos(phi);
                                double y = s(1) * sin(phi) * cos(theta);
                                double z = s(2) * sin(phi) * sin(theta);
                                point = R * Eigen::Vector3d(x, y, z) + pos;
                            }

                            vtkIdType id = points->InsertNextPoint(point(0), point(1), point(2));
                            ring_ids.push_back(id);
                        }


                        for (int j = 0; j < num_longitude; ++j)
                        {
                            if (!previous_ring_ids.empty())
                            {
                                vtkSmartPointer<vtkPolygon> polygon = vtkSmartPointer<vtkPolygon>::New();
                                polygon->GetPointIds()->SetNumberOfIds(4);

                                polygon->GetPointIds()->SetId(0, previous_ring_ids[j]);
                                polygon->GetPointIds()->SetId(1, previous_ring_ids[(j + 1) % num_longitude]);
                                polygon->GetPointIds()->SetId(2, ring_ids[(j + 1) % num_longitude]);
                                polygon->GetPointIds()->SetId(3, ring_ids[j]);
                                polys->InsertNextCell(polygon);
                            }
                        }

                        previous_ring_ids = ring_ids;
                    }
                }
                else
                {
                    // Sample points along the three main circles/ellipses composing the ellipsoid
                    std::vector<vtkIdType> ring_ids;

                    int sign = 1;
                    if (max_dot_sign < 0) {
                        sign = -1;
                        }

                    // Circle in the XY plane
                    if(max_col == 2)
                    {
                        

                        for (int i = 0; i < num_points_per_circle; ++i)
                        {
                            double theta = 2.0 * vtkMath::Pi() * i / num_points_per_circle;
                            double x = sign*s(0) * cos(theta);
                            double y = sign*s(1) * sin(theta);
                            double z = 0.0;

                            Eigen::Vector3d point = R * Eigen::Vector3d(x, y, z) + pos;
                            vtkIdType id = points->InsertNextPoint(point(0), point(1), point(2));
                            ring_ids.push_back(id);
                        }
                    }

                    // Circle in the XZ plane
                    if(max_col == 1)
                    {
                        for (int i = 0; i < num_points_per_circle; ++i)
                        {
                            double theta = 2.0 * vtkMath::Pi() * i / num_points_per_circle;
                            double x = sign*s(0) * sin(theta);
                            double y = 0.0;
                            double z = sign*s(2) * cos(theta);

                            Eigen::Vector3d point = R * Eigen::Vector3d(x, y, z) + pos;
                            vtkIdType id = points->InsertNextPoint(point(0), point(1), point(2));
                            ring_ids.push_back(id);
                        }
                    }

                    // Circle in the YZ plane
                    if(max_col == 0)
                    {
                        for (int i = 0; i < num_points_per_circle; ++i)
                        {
                            double theta = 2.0 * vtkMath::Pi() * i / num_points_per_circle;
                            double x = 0.0;
                            double y = sign*s(1) * cos(theta);
                            double z = sign*s(2) * sin(theta);

                            Eigen::Vector3d point = R * Eigen::Vector3d(x, y, z) + pos;
                            vtkIdType id = points->InsertNextPoint(point(0), point(1), point(2));
                            ring_ids.push_back(id);
                        }
                    }

                    for (int i = 0; i < num_points_per_circle; ++i)
                    {
                        if (!previous_ring_ids.empty())
                        {
                            vtkSmartPointer<vtkPolygon> polygon = vtkSmartPointer<vtkPolygon>::New();
                            polygon->GetPointIds()->SetNumberOfIds(4);
                            polygon->GetPointIds()->SetId(0, previous_ring_ids[i]);
                            polygon->GetPointIds()->SetId(1, previous_ring_ids[(i + 1) % num_points_per_circle]);
                            polygon->GetPointIds()->SetId(2, ring_ids[(i + 1) % num_points_per_circle]);
                            polygon->GetPointIds()->SetId(3, ring_ids[i]);
                            polys->InsertNextCell(polygon);
                        }
                    }

                    previous_ring_ids = ring_ids;
                }
            }

            polyData->SetPoints(points);
            polyData->SetPolys(polys);

            vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
            mapper->SetInputData(polyData);

            vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
            actor->SetMapper(mapper);
            actor->GetProperty()->SetColor(0, 0, 1);
            actor->GetProperty()->SetOpacity(0.1); // Increase opacity for better 3D effect
            actor->GetProperty()->SetAmbient(0.2); // Adjust ambient lighting
            actor->GetProperty()->SetDiffuse(0.7); // Increase diffuse lighting for better shading
            actor->GetProperty()->SetSpecular(0.5); // Increase specular for shininess
            actor->GetProperty()->SetSpecularPower(20); // Increase specular power for sharper highlights

            mp_ren->AddActor(actor);
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
    //Background
    mp_ren->SetBackground(1.,1.,1.);

    //Camera
    mp_ren->GetActiveCamera()->SetPosition(0.2,0.3,0.5);
    mp_ren->GetActiveCamera()->SetFocalPoint(0.1,0,0);
    mp_ren->GetActiveCamera()->SetViewUp(1,0,0);

        //plot one coordinate frame at the origin
    vtkSmartPointer<vtkAxesActor> axes = vtkSmartPointer<vtkAxesActor>::New();
    axes->SetTotalLength(0.025, 0.025, 0.025);
    axes->SetShaftType(0);
    axes->SetAxisLabels(0);
    mp_ren->AddActor(axes);


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
        ellipsoid_actor_init->GetProperty()->SetColor(0, 0, 1);
        ellipsoid_actor_init->GetProperty()->SetOpacity(0.1); // Increase opacity for better 3D effect
        ellipsoid_actor_init->GetProperty()->SetAmbient(0.2); // Adjust ambient lighting
        ellipsoid_actor_init->GetProperty()->SetDiffuse(0.7); // Increase diffuse lighting for better shading
        ellipsoid_actor_init->GetProperty()->SetSpecular(0.5); // Increase specular for shininess
        ellipsoid_actor_init->GetProperty()->SetSpecularPower(20); // Increase specular power for sharper highlights
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
                ellipsoid_actor->GetProperty()->SetColor(0, 0, 1);
                ellipsoid_actor->GetProperty()->SetOpacity(0.1); // Increase opacity for better 3D effect
                ellipsoid_actor->GetProperty()->SetAmbient(0.2); // Adjust ambient lighting
                ellipsoid_actor->GetProperty()->SetDiffuse(0.7); // Increase diffuse lighting for better shading
                ellipsoid_actor->GetProperty()->SetSpecular(0.5); // Increase specular for shininess
                ellipsoid_actor->GetProperty()->SetSpecularPower(20); // Increase specular power for sharper highlights
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



