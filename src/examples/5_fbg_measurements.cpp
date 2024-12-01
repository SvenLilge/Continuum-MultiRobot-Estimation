#include "visualizervtk.h"
#include "continuum_robot_state_estimator.h"
#include "utilities.h"

#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>

#include <Eigen/Dense>

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
    topology.N = 2;
    // Number of total estimation nodes per robot (including the node at the root of each robot)
    topology.K = std::vector<unsigned int>{25,21};
    // Number of interpolated states between estimation nodes per robot
    // M=1 results in no interpolation and the interpolation nodes will be equal to the estimation nodes
    // M=2 results in one additional interpolated node between each estimation node etc
    topology.M = std::vector<unsigned int>{3,3};
    // Lengths of robots
    topology.L = std::vector<double>{0.24,0.2};
    //Define if we lock the pose of the robots' ends
    topology.lock_first_pose = std::vector<bool>{true,true};
    topology.lock_last_pose = std::vector<bool>{false,false};
    topology.lock_first_strain = std::vector<bool>{false,false};
    topology.lock_last_strain = std::vector<bool>{false,false};

    topology.fbg_core_distance = std::vector<double>{37.534162e-6,37.534162e-6};
    topology.fbg_theta_offset = std::vector<double>{-0.2516,-0.7194};

    Eigen::Matrix4d T1 = Eigen::Matrix4d::Identity();

    std::string file_name_t2 = "../data/base_2.csv";

    Eigen::Matrix4d T2 = load_csv<Eigen::Matrix4d>(file_name_t2);

    Eigen::Vector3d Rx = T2.block(0,0,3,1);
    Rx = Rx/Rx.norm();
    Eigen::Vector3d Ry = T2.block(0,1,3,1);
    Ry = Ry/Ry.norm();
    Eigen::Vector3d Rz = Rx.cross(Ry);
    Rz = Rz/Rz.norm();
    Ry = Rz.cross(Rx);
    Ry = Ry/Ry.norm();

    T2.block(0,0,3,3) << Rx,Ry,Rz;



    topology.Ti0.clear();
    topology.Ti0.push_back(T1);
    topology.Ti0.push_back(T2);

    //Define coupling

    topology.common_end_effector = true;

    topology.robot_coupling.clear();

    ContinuumRobotStateEstimator::RobotTopology::Coupling coupling;
    coupling.idxA = 0;
    coupling.idxB = 2;
    coupling.coupling_node_robot_A = (topology.K.at(0)-1);
    coupling.coupling_node_robot_B = 0; //not needed for EE
    coupling.T_bA_c = Eigen::Matrix4d::Identity();
    coupling.T_bB_c <<  Eigen::Matrix4d::Identity();
    coupling.T_bB_c(2,3) = 0.05;

    Eigen::Matrix<int,6,1> mask_coupling;
    mask_coupling << 1,1,1,1,1,1; // first three are position, last three orientation
    coupling.mask = mask_coupling;

    topology.robot_coupling.push_back(coupling);


    coupling.idxA = 1;
    coupling.idxB = 2;
    coupling.coupling_node_robot_A = (topology.K.at(1)-1);
    coupling.coupling_node_robot_B = 0; //not needed for EE
    coupling.T_bA_c = Eigen::Matrix4d::Identity();
    coupling.T_bB_c <<  Eigen::Matrix4d::Identity();
    coupling.T_bB_c(2,3) = -0.05;

    mask_coupling << 1,1,1,1,1,1; // first three are position, last three orientation
    coupling.mask = mask_coupling;

    topology.robot_coupling.push_back(coupling);



    //Noise on measurements
    double R_p = 2*1e-3;
    double R_o = 1*0.05;
    double R_v = 10*0.05;
    double R_u = 10*0.05;

    double R_fbg = 10*1e-5;


    //Define hyperparameters
    ContinuumRobotStateEstimator::Hyperparameters params;

    Eigen::Matrix<double,6,1> R_pose;
    R_pose << R_p*R_p, R_p*R_p, R_p*R_p, R_o*R_o, R_o*R_o, R_o*R_o;

    Eigen::Matrix<double,6,1> R_strain;
    R_strain << R_v*R_v, R_v*R_v, R_v*R_v, R_u*R_u, R_u*R_u, R_u*R_u;

    Eigen::Matrix<double,4,1> R_fbg_strain;
    R_fbg_strain << R_fbg*R_fbg, R_fbg*R_fbg, R_fbg*R_fbg, R_fbg*R_fbg;


    Eigen::Matrix<double,6,1> R_coupling;
    R_coupling << 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6;


    Eigen::Matrix<double,6,1> Qc;
    Qc << 1e-2, 1e-2, 1e-2, 1e3, 1e3, 1e3;

    params.R_pose = 0.1*R_pose.asDiagonal();

    params.R_fbg_strain = 2e0*R_fbg_strain.asDiagonal();

    params.R_coupling = 0.1*R_coupling.asDiagonal();

    params.Qc = 0.1e0*Qc.asDiagonal();





    //Define solver options
    ContinuumRobotStateEstimator::Options options;

    options.init_guess_type = ContinuumRobotStateEstimator::Options::InitialGuessType::Straight;
    options.solver = ContinuumRobotStateEstimator::Options::Solver::NewtonLineSearch;
    options.max_optimization_iterations = 50;
    options.kirchhoff_rods = true;
    options.convergence_threshold = 5e1;

    ContinuumRobotStateEstimator state_estimator(topology, params, options);


    //Define measurements
    std::vector<ContinuumRobotStateEstimator::SensorMeasurement> measurements;


    //Load in data from experiments

    std::string file_name_fbgs_1 = "../data/fbg_strains_1.csv";
    std::string file_name_fbgs_2 = "../data/fbg_strains_2.csv";

    Eigen::MatrixXd data_fbgs_1 = load_csv<Eigen::MatrixXd>(file_name_fbgs_1);
    Eigen::MatrixXd data_fbgs_2 = load_csv<Eigen::MatrixXd>(file_name_fbgs_2);

    //Strains

    for(unsigned int n = 0; n < topology.N; n++)
    {
        for(unsigned int k = 0; k < topology.K.at(n); k++)
        {
            ContinuumRobotStateEstimator::SensorMeasurement meas;

            meas.type = ContinuumRobotStateEstimator::SensorMeasurement::FBGStrain;
            meas.idx_robot = n;
            meas.idx_node = k;
            meas.mask = Eigen::Matrix<int,6,1>(1,1,1,1,1,1);

            Eigen::Matrix<double,4,1> fbg_strain;
            if(n == 0)
            {
                fbg_strain = data_fbgs_1.col(k);
            }
            else
            {
                fbg_strain = data_fbgs_2.col(k);
            }

            meas.value = fbg_strain;
            measurements.push_back(meas);

        }
    }

    //Run state estimator
    ContinuumRobotStateEstimator::SystemState state;
    std::vector<double> cost;

    //Compute the state estimate at the estimation node and interpolate intermediate nodes as specified by M
    state_estimator.computeStateEstimate(state,cost,measurements,true);

    // state_estimator.printStateMean(state);

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




