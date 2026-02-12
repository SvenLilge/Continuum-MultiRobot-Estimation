#include "visualizervtk.h"
#include "continuum_robot_state_estimator.h"
#include "utilities.h"

#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>


// VTK Factory initialisation (for VTK version above 6)
#include <vtkAutoInit.h>
VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkRenderingFreeType);
VTK_MODULE_INIT(vtkInteractionStyle);

int main(int argc, char *argv[])
{

    // Load data

    std::string file_name = "../data/dataTongjia.csv";
    Eigen::MatrixXd data = load_csv<Eigen::MatrixXd>(file_name, true);

    // Data format
    // timestamp, disk idx 0, x, y, z, rotX, rotY, rotZ, disk idx 1, x, y, z, rotX, rotY, rotZ, etc

    int sample = 3000;

    // Load data from csv file (disk indices 0 to 6)
    std::vector<Eigen::Matrix4d> T_disks;
    for(unsigned int i = 0; i < 7; i++)
    {
        // get the row indicated by sample from the csv file and get the transformation matrix for disk i
        Eigen::Matrix4d T_disk = Eigen::Matrix4d::Identity();
        // position
        T_disk.block(0,3,3,1) << data(sample, 1 + i*7 + 1), data(sample, 1 + i*7 + 2), data(sample, 1 + i*7 + 3);
        // rotation from euler angles (eul2rotm(angles,“XYZ”))
        double rotX = data(sample, 1 + i*7 + 4);
        double rotY = data(sample, 1 + i*7 + 5);
        double rotZ = data(sample, 1 + i*7 + 6);

        Eigen::Matrix3d Rx, Ry, Rz;
        double cx = std::cos(rotX), sx = std::sin(rotX);
        double cy = std::cos(rotY), sy = std::sin(rotY);
        double cz = std::cos(rotZ), sz = std::sin(rotZ);

        // Rotation about X
        Rx << 1,  0,   0,
            0, cx, -sx,
            0, sx,  cx;

        // Rotation about Y
        Ry <<  cy, 0, sy,
              0, 1,  0,
            -sy, 0, cy;

        // Rotation about Z
        Rz << cz, -sz, 0,
            sz,  cz, 0,
             0,   0, 1;

        Eigen::Matrix3d R = Rx * Ry * Rz;
        T_disk.block(0,0,3,3) = R;

        T_disks.push_back(T_disk);
    }

    // Now for all transformations, I want:
    // the z axis to be the x axis
    // the x axis to be the y axis
    // the y axis to be the z axis
    // Update just the rotation part of the transformation matrix accordingly
    for(unsigned int i = 0; i < T_disks.size(); i++)
    {
        Eigen::Matrix3d R = T_disks.at(i).block(0,0,3,3);
        Eigen::Matrix3d R_new;
        R_new.block(0,0,3,1) = R.block(0,2,3,1); // new x axis is old z axis
        R_new.block(0,1,3,1) = R.block(0,0,3,1); // new y axis is old x axis
        R_new.block(0,2,3,1) = R.block(0,1,3,1); // new z axis is old y axis
        T_disks.at(i).block(0,0,3,3) = R_new;
    }


    //Define robot topology
    ContinuumRobotStateEstimator::RobotTopology topology;

    // Number of robots
    topology.N =2;
    // Number of total estimation nodes per robot (including the node at the root of each robot)
    topology.K = std::vector<unsigned int>{11,11};
    // Number of interpolated states between estimation nodes per robot
    // M=1 results in no interpolation and the interpolation nodes will be equal to the estimation nodes
    // M=2 results in one additional interpolated node between each estimation node etc
    topology.M = std::vector<unsigned int>{3,3};
    // Lengths of robots
    topology.L = std::vector<double>{0.54,0.54};
    //Define if we lock the pose of the robots' ends
    topology.lock_first_pose = std::vector<bool>{true,true};
    topology.lock_last_pose = std::vector<bool>{false,false};
    topology.lock_first_strain = std::vector<bool>{false,false};
    topology.lock_last_strain = std::vector<bool>{false,false};

    topology.fbg_core_distance = std::vector<double>{0,0};
    topology.fbg_theta_offset = std::vector<double>{0,0};


    topology.Ti0.clear();

    topology.Ti0.push_back(T_disks.at(3));
    topology.Ti0.push_back(T_disks.at(0));



    //Define coupling

    topology.common_end_effector = false;

    topology.robot_coupling.clear();

    ContinuumRobotStateEstimator::RobotTopology::Coupling coupling;
    coupling.idxA = 1;
    coupling.idxB = 0;
    coupling.coupling_node_robot_A = (topology.K.at(1)-1); //not needed for EE
    coupling.coupling_node_robot_B = (topology.K.at(0)-1)/2; //not needed for EE
    coupling.T_bA_c = Eigen::Matrix4d::Identity();
    coupling.T_bB_c = Eigen::Matrix4d::Identity();
    coupling.T_bA_c.block(0,0,3,3) << 0,1,0,
                                      -1,0,0,
                                      0,0,1;

    Eigen::Matrix<int,6,1> mask_coupling;
    mask_coupling << 1,1,1,1,1,1; // first three are position, last three orientation
    coupling.mask = mask_coupling;

    topology.robot_coupling.push_back(coupling);


    //Noise on measurements
    double R_p = 2*1e-3;
    double R_o = 1*0.05;
    double R_v = 1*0.05;
    double R_u = 1*0.05;

    double R_fbg = 1*1e-5;


    //Define hyperparameters
    ContinuumRobotStateEstimator::Hyperparameters params;

    Eigen::Matrix<double,6,1> R_pose;
    R_pose << R_p*R_p, R_p*R_p, R_p*R_p, R_o*R_o, R_o*R_o, R_o*R_o;

    Eigen::Matrix<double,6,1> R_strain;
    R_strain << R_v*R_v, R_v*R_v, R_v*R_v, R_u*R_u, R_u*R_u, R_u*R_u;

    Eigen::Matrix<double,4,1> R_fbg_strain;
    R_fbg_strain << R_fbg*R_fbg, R_fbg*R_fbg, R_fbg*R_fbg, R_fbg*R_fbg;


    Eigen::Matrix<double,6,1> R_coupling;
    R_coupling << 1, 1, 1, 1, 1, 1;


    Eigen::Matrix<double,6,1> Qc;
    Qc << 1e-1, 1e-1, 1e-1, 1e2, 1e2, 1e2;



    params.R_pose = 1*R_pose.asDiagonal();
    params.R_strain = 10*R_strain.asDiagonal();

    params.R_fbg_strain = 20e0*R_fbg_strain.asDiagonal();

    params.R_coupling = 1e-8*R_coupling.asDiagonal();

    params.Qc = 1e0*Qc.asDiagonal();





    //Define solver options
    ContinuumRobotStateEstimator::Options options;

    options.init_guess_type = ContinuumRobotStateEstimator::Options::InitialGuessType::Straight;
    options.solver = ContinuumRobotStateEstimator::Options::Solver::NewtonLineSearch;
    options.max_optimization_iterations = 200;
    options.kirchhoff_rods = true;
    options.convergence_threshold = 5e-1;


    ContinuumRobotStateEstimator state_estimator(topology, params, options);


    //Define measurements
    std::vector<ContinuumRobotStateEstimator::SensorMeasurement> measurements;

    //Pose
    ContinuumRobotStateEstimator::SensorMeasurement meas;
    Eigen::Matrix4d pose = T_disks.at(6);

    meas.type = ContinuumRobotStateEstimator::SensorMeasurement::Pose;
    meas.idx_robot = 0;
    meas.idx_node = (topology.K.at(0)-1); //not needed for EE
    meas.mask = Eigen::Matrix<int,6,1>(1,1,1,1,1,1);
    meas.value = pose;
    measurements.push_back(meas);


    //Run state estimator
    ContinuumRobotStateEstimator::SystemState state;
    std::vector<double> cost;

    //Compute the state estimate at the estimation node and interpolate intermediate nodes as specified by M
    state_estimator.computeStateEstimate(state,cost,measurements,true);

    //state_estimator.printStateMean(state);

    Visualizer vis(topology);

    //Update the visualizer with the state
    vis.update(state,true,true,3);


    //Create Window Interactor
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(vis.getRenderWindow());
	
	//Set up and start main loop
	renderWindowInteractor->UpdateSize(1280,720);
	vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
    renderWindowInteractor->SetInteractorStyle(style);
	renderWindowInteractor->Initialize();
	renderWindowInteractor->Start();

    return 1;
}




