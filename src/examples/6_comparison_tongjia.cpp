#include "visualizervtk.h"
#include "continuum_robot_state_estimator.h"
#include "utilities.h"

#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkCallbackCommand.h>

#include <atomic>
#include <chrono>
#include <mutex>
#include <thread>


// VTK Factory initialisation (for VTK version above 6)
#include <vtkAutoInit.h>
VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkRenderingFreeType);
VTK_MODULE_INIT(vtkInteractionStyle);

int main(int argc, char *argv[])
{
    // Load data
    std::string file_name = "../data/RSS2026/base 1/multiCR oscillation 1/dataVicon.csv";
    //std::string file_name = "../data/dataTongjia.csv";
    Eigen::MatrixXd data = load_csv<Eigen::MatrixXd>(file_name, true);

    // Data format:
    // timestamp, disk idx 0, x, y, z, rotX, rotY, rotZ, disk idx 1, x, y, z, rotX, rotY, rotZ, etc

    // Sample range to iterate over
    int sample_start = 0;
    int sample_end = static_cast<int>(data.rows()) - 1;
    sample_end = sample_start; 
    if(sample_end < 0) return 1;
    int sample_step = 10;

    // Prepare topology based on the first sample (sample_start)
    int sample = sample_start;

    // Build T_disks for the initial sample (same code as before)
    std::vector<Eigen::Matrix4d> T_disks_init;
    for(unsigned int i = 0; i < 7; i++)
    {
        Eigen::Matrix4d T_disk = Eigen::Matrix4d::Identity();
        T_disk.block(0,3,3,1) << data(sample, 1 + i*7 + 1), data(sample, 1 + i*7 + 2), data(sample, 1 + i*7 + 3);
        double rotX = data(sample, 1 + i*7 + 4);
        double rotY = data(sample, 1 + i*7 + 5);
        double rotZ = data(sample, 1 + i*7 + 6);

        Eigen::Matrix3d Rx, Ry, Rz;
        double cx = std::cos(rotX), sx = std::sin(rotX);
        double cy = std::cos(rotY), sy = std::sin(rotY);
        double cz = std::cos(rotZ), sz = std::sin(rotZ);

        Rx << 1,  0,   0,
              0, cx, -sx,
              0, sx,  cx;

        Ry <<  cy, 0, sy,
               0, 1,  0,
             -sy, 0, cy;

        Rz << cz, -sz, 0,
              sz,  cz, 0,
               0,   0, 1;

        Eigen::Matrix3d R = Rx * Ry * Rz;
        T_disk.block(0,0,3,3) = R;
        T_disks_init.push_back(T_disk);
    }
    // Axis permutation
    for(unsigned int i = 0; i < T_disks_init.size(); i++)
    {
        Eigen::Matrix3d R = T_disks_init.at(i).block(0,0,3,3);
        Eigen::Matrix3d R_new;
        R_new.block(0,0,3,1) = R.block(0,2,3,1); // new x axis is old z axis
        R_new.block(0,1,3,1) = R.block(0,0,3,1); // new y axis is old x axis
        R_new.block(0,2,3,1) = R.block(0,1,3,1); // new z axis is old y axis
        T_disks_init.at(i).block(0,0,3,3) = R_new;
    }

    // Define robot topology (use initial transforms)
    ContinuumRobotStateEstimator::RobotTopology topology;

    topology.N =2;
    topology.K = std::vector<unsigned int>{13,13};
    topology.M = std::vector<unsigned int>{2,2};
    topology.L = std::vector<double>{0.54,0.54};
    topology.lock_first_pose = std::vector<bool>{true,true};
    topology.lock_last_pose = std::vector<bool>{false,false};
    topology.lock_first_strain = std::vector<bool>{false,false};
    topology.lock_last_strain = std::vector<bool>{false,false};
    topology.fbg_core_distance = std::vector<double>{0,0};
    topology.fbg_theta_offset = std::vector<double>{0,0};

    topology.Ti0.clear();
    topology.Ti0.push_back(T_disks_init.at(3));
    topology.Ti0.push_back(T_disks_init.at(0));

    topology.common_end_effector = false;
    topology.robot_coupling.clear();

    ContinuumRobotStateEstimator::RobotTopology::Coupling coupling;
    coupling.idxA = 1;
    coupling.idxB = 0;
    coupling.coupling_node_robot_A = (topology.K.at(1)-1);
    coupling.coupling_node_robot_B = (topology.K.at(0)-1)/2;
    coupling.T_bA_c = Eigen::Matrix4d::Identity();
    coupling.T_bB_c = Eigen::Matrix4d::Identity();
    coupling.T_bA_c.block(0,0,3,3) << 0,1,0,
                                      -1,0,0,
                                       0,0,1;
    Eigen::Matrix<int,6,1> mask_coupling;
    mask_coupling << 1,1,1,1,1,1;
    coupling.mask = mask_coupling;
    topology.robot_coupling.push_back(coupling);

    // Noise on measurements
    double R_p = 2*1e-3;
    double R_o = 1*0.05;
    double R_v = 1*0.05;
    double R_u = 1*0.05;
    double R_fbg = 1*1e-5;

    // Hyperparameters
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
    Qc << 1e-1, 1e-1, 1e-1, 1e1, 1e1, 1e1;

    params.R_pose = 1e-1*R_pose.asDiagonal();
    params.R_strain = 10*R_strain.asDiagonal();
    params.R_fbg_strain = 20e0*R_fbg_strain.asDiagonal();
    params.R_coupling = 1e-8*R_coupling.asDiagonal();
    params.Qc = 1e0*Qc.asDiagonal();

    // Solver options
    ContinuumRobotStateEstimator::Options options;
    options.init_guess_type = ContinuumRobotStateEstimator::Options::InitialGuessType::Last;
    options.solver = ContinuumRobotStateEstimator::Options::Solver::Newton;
    options.max_optimization_iterations = 20;
    options.kirchhoff_rods = true;
    options.convergence_threshold = 5e-1;

    // Create estimator
    ContinuumRobotStateEstimator state_estimator(topology, params, options);

    // Visualizer
    Visualizer vis(topology);

    // Create window interactor and attach renderer
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(vis.getRenderWindow());
    renderWindowInteractor->UpdateSize(1280,720);
    vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
    renderWindowInteractor->SetInteractorStyle(style);

    // Prepare measurement template (only pose measurement used in example)
    ContinuumRobotStateEstimator::SensorMeasurement measTemplate;
    measTemplate.type = ContinuumRobotStateEstimator::SensorMeasurement::Pose;
    measTemplate.idx_robot = 0;
    measTemplate.idx_node = (topology.K.at(0)-1);
    measTemplate.mask = Eigen::Matrix<int,6,1>(1,1,1,1,1,1);

    // Shared state for communication between worker thread and UI thread
    struct SharedState {
        std::mutex mtx;
        ContinuumRobotStateEstimator::SystemState latest_state;
        std::vector<double> latest_cost;
        std::atomic<bool> new_state{false};
        std::atomic<bool> finished{false};
        Visualizer* vis{nullptr};
    } shared;
    shared.vis = &vis;

    std::vector<Eigen::MatrixXd> accuracy; // For each sample, store error between estimate and measurement for each disk


    std::vector<Eigen::MatrixXd> estimates; // For each sample, store the estimated T_disks


    std::vector<double> computation_times; // For each sample, store the computation time of the state estimation


    // Worker thread: computes estimates off the UI thread and publishes latest result
    std::thread worker([&](){
        for(int s = sample_start; s <= sample_end; s += sample_step)
        {
            // Build T_disks for sample s
            std::vector<Eigen::Matrix4d> T_disks;
            for(unsigned int i = 0; i < 7; i++)
            {
                Eigen::Matrix4d T_disk = Eigen::Matrix4d::Identity();
                T_disk.block(0,3,3,1) << data(s, 1 + i*7 + 1), data(s, 1 + i*7 + 2), data(s, 1 + i*7 + 3);
                double rotX = data(s, 1 + i*7 + 4);
                double rotY = data(s, 1 + i*7 + 5);
                double rotZ = data(s, 1 + i*7 + 6);

                Eigen::Matrix3d Rx, Ry, Rz;
                double cx = std::cos(rotX), sx = std::sin(rotX);
                double cy = std::cos(rotY), sy = std::sin(rotY);
                double cz = std::cos(rotZ), sz = std::sin(rotZ);

                Rx << 1,  0,   0,
                      0, cx, -sx,
                      0, sx,  cx;

                Ry <<  cy, 0, sy,
                       0, 1,  0,
                     -sy, 0, cy;

                Rz << cz, -sz, 0,
                      sz,  cz, 0,
                       0,   0, 1;

                Eigen::Matrix3d R = Rx * Ry * Rz;
                T_disk.block(0,0,3,3) = R;
                T_disks.push_back(T_disk);
            }
            // Axis permutation
            for(unsigned int i = 0; i < T_disks.size(); i++)
            {
                Eigen::Matrix3d R = T_disks.at(i).block(0,0,3,3);
                Eigen::Matrix3d R_new;
                R_new.block(0,0,3,1) = R.block(0,2,3,1); // new x axis is old z axis
                R_new.block(0,1,3,1) = R.block(0,0,3,1); // new y axis is old x axis
                R_new.block(0,2,3,1) = R.block(0,1,3,1); // new z axis is old y axis
                T_disks.at(i).block(0,0,3,3) = R_new;
            }

            // Prepare measurement for this sample
            ContinuumRobotStateEstimator::SensorMeasurement meas = measTemplate;
            meas.value = T_disks.at(6);

            std::vector<ContinuumRobotStateEstimator::SensorMeasurement> measurements;
            measurements.push_back(meas);

            // Compute state estimate (heavy work off UI thread)
            ContinuumRobotStateEstimator::SystemState state;
            std::vector<double> cost;

            auto start = std::chrono::high_resolution_clock::now();
            state_estimator.computeStateEstimate(state, cost, measurements, false);
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            computation_times.push_back(duration.count() / 1000.0); // Store in milliseconds

            // Store estimate of each disk for analysis
            Eigen::MatrixXd T_disk_estimated;
            T_disk_estimated.resize(4*7,4); // 7 disks, each with a 4x4 transform
            // First disk is first node of second robot
            T_disk_estimated.block(0,0,4,4) = state.robots.at(1).estimation_nodes.at(0).pose;
            // Second disk is fifth node of second robot
            T_disk_estimated.block(4,0,4,4) = state.robots.at(1).estimation_nodes.at(4).pose;
            // Third disk is ninth node of second robot
            T_disk_estimated.block(8,0,4,4) = state.robots.at(1).estimation_nodes.at(8).pose;
            // Fourth disk is first node of first robot
            T_disk_estimated.block(12,0,4,4) = state.robots.at(0).estimation_nodes.at(0).pose;
            // Fifth disk is fifth node of first robot
            T_disk_estimated.block(16,0,4,4) = state.robots.at(0).estimation_nodes.at(4).pose;
            // Sixth disk is ninth node of first robot
            T_disk_estimated.block(20,0,4,4) = state.robots.at(0).estimation_nodes.at(8).pose;
            // Seventh disk is thirteenth node of first robot
            T_disk_estimated.block(24,0,4,4) = state.robots.at(0).estimation_nodes.at(12).pose;
            estimates.push_back(T_disk_estimated);

            // Compare estimate to measurement and store accuracy
            Eigen::MatrixXd acc;
            acc.resize(7,2); // 7 disks, each with 2 error components (pos and rot)
            for(unsigned int d = 0; d < 7; d++)
            {
                Eigen::Matrix4d T_est = T_disk_estimated.block(d*4,0,4,4);
                Eigen::Matrix4d T_meas = T_disks.at(d);

                Eigen::Matrix4d T_err = T_est * invert_transformation(T_meas);
                Eigen::Matrix<double,6,1> err_vec = tran_to_vec(T_err);

                acc(d,0) = err_vec.block(0,0,3,1).norm(); // position error
                acc(d,1) = err_vec.block(3,0,3,1).norm(); // orientation error
            }

            accuracy.push_back(acc);

            // Publish result (wait if UI hasn't consumed previous)
            while(shared.new_state.load()) std::this_thread::sleep_for(std::chrono::milliseconds(5));
            {
                std::lock_guard<std::mutex> lk(shared.mtx);
                shared.latest_state = state;
                shared.latest_cost = cost;
                shared.new_state.store(true);
            }

            // Optional small delay to avoid saturating CPU and allow UI time to render
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
        }
        shared.finished.store(true);

        // After finishing all samples, print accuracy results (average all samples)
        Eigen::MatrixXd acc_sum;
        acc_sum.resize(7,2);
        acc_sum.setZero();
        for(int i = 0; i < accuracy.size(); i++)        {
            acc_sum = acc_sum + accuracy.at(i);
        }
        Eigen::MatrixXd acc_avg = acc_sum / static_cast<double>(accuracy.size());
        std::cout << "Average accuracy over " << accuracy.size() << " samples:" << std::endl;
        for(unsigned int d = 0; d < 7; d++)        {
            std::cout << "Disk " << d << ": Position error = " << acc_avg(d,0) << " m, Orientation error = " << acc_avg(d,1) << " rad" << std::endl;
        }

        // Print average computation time
        double comp_time_sum = 0.0;
        for(int i = 0; i < computation_times.size(); i++)        {
            comp_time_sum += computation_times.at(i);
        }
        double comp_time_avg = comp_time_sum / static_cast<double>(computation_times.size());
        std::cout << "Average computation time: " << comp_time_avg << " ms" << std::endl;

    });

    // Timer callback: if worker published a new state, update visualizer and render
    vtkSmartPointer<vtkCallbackCommand> timerCallback = vtkSmartPointer<vtkCallbackCommand>::New();
    timerCallback->SetClientData(&shared);
    timerCallback->SetCallback([](vtkObject* caller, unsigned long, void* clientData, void*){
        SharedState* s = static_cast<SharedState*>(clientData);
        if(s->new_state.load())
        {
            ContinuumRobotStateEstimator::SystemState state_copy;
            {
                std::lock_guard<std::mutex> lk(s->mtx);
                state_copy = s->latest_state;
                s->new_state.store(false);
            }
            // update visualizer and render
            if(s->vis)
            {
                s->vis->update(state_copy, true, true, 3);
                s->vis->getRenderWindow()->Render();
            }
        }
    });

    renderWindowInteractor->AddObserver(vtkCommand::TimerEvent, timerCallback.GetPointer());
    int timerId = renderWindowInteractor->CreateRepeatingTimer(100); // 100 ms tick

    // Start interaction (UI remains responsive). Worker thread computes in background and timer callback applies updates.
    renderWindowInteractor->Initialize();
    renderWindowInteractor->Start();

    // Interactor stopped (window closed) - join worker
    if(worker.joinable()) worker.join();

    // Destroy timer
    if(timerId > 0) renderWindowInteractor->DestroyTimer(timerId);

    return 1;
}




