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

    //Define robot topology
    ContinuumRobotStateEstimator::RobotTopology topology;

    // Number of robots
    topology.N =3;
    // Number of total estimation nodes per robot (including the node at the root of each robot)
    topology.K = std::vector<unsigned int>{31,16,16};
    // Number of interpolated states between estimation nodes per robot
    // M=1 results in no interpolation and the interpolation nodes will be equal to the estimation nodes
    // M=2 results in one additional interpolated node between each estimation node etc
    topology.M = std::vector<unsigned int>{3,3,3};
    // Lengths of robots
    topology.L = std::vector<double>{0.30,0.15,0.15};
    //Define if we lock the pose of the robots' ends
    topology.lock_first_pose = std::vector<bool>{true,true,true};
    topology.lock_last_pose = std::vector<bool>{false,false,false};
    topology.lock_first_strain = std::vector<bool>{false,false,false};
    topology.lock_last_strain = std::vector<bool>{false,false,false};

    topology.fbg_core_distance = std::vector<double>{0,0,0};
    topology.fbg_theta_offset = std::vector<double>{0,0,0};


    topology.Ti0.clear();

    Eigen::Matrix4d T1 = Eigen::Matrix4d::Identity();

    topology.Ti0.push_back(T1);


    Eigen::Matrix4d T2 = Eigen::Matrix4d::Identity();
    T2.block(0,0,3,3) << 0,1,0,
                         -1,0,0,
                         0,0,1;

    T2.block(0,3,3,1) << 0.1,0.1,0;

    topology.Ti0.push_back(T2);



    Eigen::Matrix4d T3 = Eigen::Matrix4d::Identity();
    T3.block(0,0,3,3) << 0,-1,0,
                         1,0,0,
                         0,0,1;

    T3.block(0,3,3,1) << 0.2,-0.1,0;

    topology.Ti0.push_back(T3);



    //Define coupling

    topology.common_end_effector = false;

    topology.robot_coupling.clear();

    ContinuumRobotStateEstimator::RobotTopology::Coupling coupling;
    coupling.idxA = 1;
    coupling.idxB = 0;
    coupling.coupling_node_robot_A = (topology.K.at(1)-1); //not needed for EE
    coupling.coupling_node_robot_B = (topology.K.at(0)-4)/2; //not needed for EE
    coupling.T_bA_c = Eigen::Matrix4d::Identity();
    coupling.T_bB_c = Eigen::Matrix4d::Identity();
    coupling.T_bB_c.block(0,0,3,3) << 0,1,0,
                                      -1,0,0,
                                      0,0,1;

    Eigen::Matrix<int,6,1> mask_coupling;
    mask_coupling << 1,1,1,1,1,1; // first three are position, last three orientation
    coupling.mask = mask_coupling;

    topology.robot_coupling.push_back(coupling);



    coupling.idxA = 2;
    coupling.idxB = 0;
    coupling.coupling_node_robot_A = (topology.K.at(2)-1); //not needed for EE
    coupling.coupling_node_robot_B = (topology.K.at(0)-4); //not needed for EE
    coupling.T_bA_c = Eigen::Matrix4d::Identity();
    coupling.T_bB_c = Eigen::Matrix4d::Identity();
    coupling.T_bB_c.block(0,0,3,3) << 0,-1,0,
                                      1,0,0,
                                      0,0,1;

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

    params.R_coupling = 0.5*1e-6*R_coupling.asDiagonal();

    params.Qc = 1e0*Qc.asDiagonal();





    //Define solver options
    ContinuumRobotStateEstimator::Options options;

    options.init_guess_type = ContinuumRobotStateEstimator::Options::InitialGuessType::Straight;
    options.solver = ContinuumRobotStateEstimator::Options::Solver::NewtonLineSearch;
    options.max_optimization_iterations = 200;
    options.kirchhoff_rods = true;
    options.convergence_threshold = 5e1;


    ContinuumRobotStateEstimator state_estimator(topology, params, options);


    //Define measurements
    std::vector<ContinuumRobotStateEstimator::SensorMeasurement> measurements;

    //Pose
    ContinuumRobotStateEstimator::SensorMeasurement mes;
    Eigen::Matrix4d pose = Eigen::Matrix4d::Identity();

    pose.block(0,3,3,1) << 0.2, 0.05, 0;

    mes.type = ContinuumRobotStateEstimator::SensorMeasurement::Pose;
    mes.idx_robot = 0;
    mes.idx_node = 30;
    mes.mask = Eigen::Matrix<int,6,1>(1,1,1,0,0,0);
    mes.value = pose;
    measurements.push_back(mes);


    //Run state estimator
    ContinuumRobotStateEstimator::SystemState state;
    std::vector<double> cost;

    //Compute the state estimate at the estimation node and interpolate intermediate nodes as specified by M
    state_estimator.computeStateEstimate(state,cost,measurements,true);

    state_estimator.printStateMean(state);

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




