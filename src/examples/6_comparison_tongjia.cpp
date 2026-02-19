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
#include <fstream>
#include <iomanip>


// VTK Factory initialisation (for VTK version above 6)
#include <vtkAutoInit.h>
VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkRenderingFreeType);
VTK_MODULE_INIT(vtkInteractionStyle);

// Structure to hold sample information for processing
struct SampleInfo {
    int row1, row2;  // surrounding rows for interpolation
    double alpha;     // interpolation factor for refresh rate mode
};

// Helper function: Build T_disks from position and Euler angles (rotX, rotY, rotZ)
Eigen::Matrix4d buildDiskTransformFromEulerData(double x, double y, double z, 
                                                       double rotX, double rotY, double rotZ)
{
    std::vector<Eigen::Matrix4d> T_disks;
    
    // Build Euler angle rotation matrices
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
    
    
    Eigen::Matrix4d T_disk = Eigen::Matrix4d::Identity();
    T_disk.block(0,3,3,1) << x, y, z;
    T_disk.block(0,0,3,3) = R;
    
    // Axis permutation
    Eigen::Matrix3d R_rot = T_disk.block(0,0,3,3);
    Eigen::Matrix3d R_new;
    R_new.block(0,0,3,1) = R_rot.block(0,2,3,1); // new x axis is old z axis
    R_new.block(0,1,3,1) = R_rot.block(0,0,3,1); // new y axis is old x axis
    R_new.block(0,2,3,1) = R_rot.block(0,1,3,1); // new z axis is old y axis
    T_disk.block(0,0,3,3) = R_new;
    
    return T_disk;
}

// Helper function: Extract position and Euler angles from a data row for a specific disk
void extractDiskData(const Eigen::MatrixXd& data, int row, int disk_idx,
                     double& x, double& y, double& z, double& rotX, double& rotY, double& rotZ)
{
    x = data(row, 1 + disk_idx*7 + 1);
    y = data(row, 1 + disk_idx*7 + 2);
    z = data(row, 1 + disk_idx*7 + 3);
    rotX = data(row, 1 + disk_idx*7 + 4);
    rotY = data(row, 1 + disk_idx*7 + 5);
    rotZ = data(row, 1 + disk_idx*7 + 6);
}

// Helper function: Interpolate disk data between two rows
void interpolateDiskData(const Eigen::MatrixXd& data, int row1, int row2, double alpha,
                        int disk_idx, double& x, double& y, double& z, 
                        double& rotX, double& rotY, double& rotZ)
{
    double x1, y1, z1, rotX1, rotY1, rotZ1;
    double x2, y2, z2, rotX2, rotY2, rotZ2;
    
    extractDiskData(data, row1, disk_idx, x1, y1, z1, rotX1, rotY1, rotZ1);
    extractDiskData(data, row2, disk_idx, x2, y2, z2, rotX2, rotY2, rotZ2);
    
    // Linear interpolation for position and Euler angles
    x = (1.0 - alpha) * x1 + alpha * x2;
    y = (1.0 - alpha) * y1 + alpha * y2;
    z = (1.0 - alpha) * z1 + alpha * z2;
    rotX = (1.0 - alpha) * rotX1 + alpha * rotX2;
    rotY = (1.0 - alpha) * rotY1 + alpha * rotY2;
    rotZ = (1.0 - alpha) * rotZ1 + alpha * rotZ2;
}

// Helper function: Convert rotation matrix back to Euler angles (XYZ convention)
// This reverses the construction: R = Rx * Ry * Rz
void rotationMatrixToEulerXYZ(const Eigen::Matrix3d& R, double& rotX, double& rotY, double& rotZ)
{
    // Extract Euler angles from rotation matrix using XYZ convention
    // R = Rx(rotX) * Ry(rotY) * Rz(rotZ)
    
    rotY = std::asin(-R(2,0));
    
    if(std::abs(std::cos(rotY)) > 1e-6)
    {
        rotX = std::atan2(R(2,1), R(2,2));
        rotZ = std::atan2(R(1,0), R(0,0));
    }
    else
    {
        // Gimbal lock case
        rotX = 0.0;
        rotZ = std::atan2(-R(0,1), R(1,1));
    }
}

// Helper function: Revert axis permutation and extract Euler angles from T_disk
void extractPoseFromTransform(const Eigen::Matrix4d& T_disk, double& x, double& y, double& z,
                              double& rotX, double& rotY, double& rotZ)
{
    // Extract position
    x = T_disk(0,3);
    y = T_disk(1,3);
    z = T_disk(2,3);
    
    // Revert axis permutation: the transform has permuted axes, we need to reverse it
    Eigen::Matrix3d R_permuted = T_disk.block(0,0,3,3);
    Eigen::Matrix3d R_original;
    R_original.block(0,0,3,1) = R_permuted.block(0,1,3,1); // old x axis was new y axis
    R_original.block(0,1,3,1) = R_permuted.block(0,2,3,1); // old y axis was new z axis
    R_original.block(0,2,3,1) = R_permuted.block(0,0,3,1); // old z axis was new x axis
    
    // Convert rotation matrix to Euler angles
    rotationMatrixToEulerXYZ(R_original, rotX, rotY, rotZ);
}

// Helper function: Export results to CSV file
void exportResultsToCSV(const std::string& output_file,
                        const std::vector<SampleInfo>& samples_to_process,
                        const Eigen::MatrixXd& data,
                        const std::vector<std::vector<Eigen::Matrix4d>>& estimated_T_disks,
                        const std::vector<double>& computation_times)
{
    std::ofstream outfile(output_file);
    
    // Write header
    outfile << "timestamp";
    for(unsigned int d = 0; d < 7; d++)
    {
        outfile << ",disk_" << d << "_x,disk_" << d << "_y,disk_" << d << "_z,"
                << "disk_" << d << "_rotX,disk_" << d << "_rotY,disk_" << d << "_rotZ";
    }
    outfile << ",computation_time_ms" << std::endl;
    
    // Write data rows
    for(size_t sample_idx = 0; sample_idx < samples_to_process.size(); ++sample_idx)
    {
        const auto& sample_info = samples_to_process[sample_idx];
        
        // Interpolate timestamp
        double t_interp;
        if(sample_info.alpha > 1e-6)
        {
            double t1 = data(sample_info.row1, 0);
            double t2 = data(sample_info.row2, 0);
            t_interp = (1.0 - sample_info.alpha) * t1 + sample_info.alpha * t2;
        }
        else
        {
            t_interp = data(sample_info.row1, 0);
        }
        
        outfile << std::fixed << std::setprecision(6) << t_interp;
        
        // Write pose for each disk
        for(unsigned int d = 0; d < 7 && d < estimated_T_disks[sample_idx].size(); ++d)
        {
            double x, y, z, rotX, rotY, rotZ;
            extractPoseFromTransform(estimated_T_disks[sample_idx][d], x, y, z, rotX, rotY, rotZ);
            
            outfile << "," << std::fixed << std::setprecision(6) << x
                    << "," << std::fixed << std::setprecision(6) << y
                    << "," << std::fixed << std::setprecision(6) << z
                    << "," << std::fixed << std::setprecision(6) << rotX
                    << "," << std::fixed << std::setprecision(6) << rotY
                    << "," << std::fixed << std::setprecision(6) << rotZ;
        }
        
        // Write computation time
        outfile << "," << std::fixed << std::setprecision(3) << computation_times[sample_idx];
        outfile << std::endl;
    }
    
    outfile.close();
    std::cout << "Results exported to " << output_file << std::endl;
}

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
    if(sample_end < 0) return 1;
    
    // Sampling mode: either sample_step or refresh_rate
    bool use_refresh_rate = true;
    double refresh_rate = 40.0; // Hz, default
    int sample_step = 10;
    
    // Parse command line arguments
    for(int i = 1; i < argc; i++)
    {
        std::string arg = argv[i];
        if(arg == "--refresh-rate" && i + 1 < argc)
        {
            use_refresh_rate = true;
            refresh_rate = std::stod(argv[i + 1]);
            i++;
            std::cout << "Using refresh rate mode: " << refresh_rate << " Hz" << std::endl;
        }
        else if(arg == "--sample-step" && i + 1 < argc)
        {
            use_refresh_rate = false;
            sample_step = std::stoi(argv[i + 1]);
            i++;
            std::cout << "Using sample step mode with step size: " << sample_step << std::endl;
        }
    }
    
    // Generate sample indices and interpolation parameters based on mode
    std::vector<SampleInfo> samples_to_process;
    
    if(use_refresh_rate)
    {
        // Generate time-based samples with interpolation
        double t_start = data(sample_start, 0);
        double t_end = data(sample_end, 0);
        double dt = 1.0 / refresh_rate;
        
        for(double t = t_start; t <= t_end; t += dt)
        {
            // Find the two surrounding samples in the data
            int idx_lower = sample_start;
            int idx_upper = sample_end;
            
            for(int s = sample_start; s <= sample_end; s++)
            {
                if(data(s, 0) <= t)
                    idx_lower = s;
                if(data(s, 0) >= t && idx_upper == sample_end)
                {
                    idx_upper = s;
                    break;
                }
            }
            
            double alpha = 0.0;
            if(idx_lower != idx_upper)
            {
                double t1 = data(idx_lower, 0);
                double t2 = data(idx_upper, 0);
                if(t2 > t1)
                    alpha = (t - t1) / (t2 - t1);
            }
            
            samples_to_process.push_back({idx_lower, idx_upper, alpha});
        }
        
        std::cout << "Generated " << samples_to_process.size() << " samples at " << refresh_rate << " Hz" << std::endl;
    }
    else
    {
        // Use sample step
        for(int s = sample_start; s <= sample_end; s += sample_step)
        {
            samples_to_process.push_back({s, s, 0.0});
        }
        std::cout << "Using sample step of " << sample_step << ", processing " << samples_to_process.size() << " samples" << std::endl;
    }

    // Prepare topology based on the first sample (sample_start)
    double x, y, z, rotX, rotY, rotZ;
    std::vector<Eigen::Matrix4d> T_disks_init;
    for(int disk_idx = 0; disk_idx < 7; disk_idx++) // Assuming 2 disks for this example
    {
        extractDiskData(data, sample_start, disk_idx, x, y, z, rotX, rotY, rotZ);
        Eigen::Matrix4d T_disk = buildDiskTransformFromEulerData(x, y, z, rotX, rotY, rotZ);
        T_disks_init.push_back(T_disk);
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
        std::vector<Eigen::Matrix4d> latest_T_disks; // Ground truth T_disks
        std::atomic<bool> new_state{false};
        std::atomic<bool> finished{false};
        Visualizer* vis{nullptr};
    } shared;
    shared.vis = &vis;

    std::vector<Eigen::MatrixXd> accuracy; // For each sample, store error between estimate and measurement for each disk

    std::vector<Eigen::MatrixXd> estimates; // For each sample, store the estimated T_disks

    std::vector<double> computation_times; // For each sample, store the computation time of the state estimation
    
    std::vector<std::vector<Eigen::Matrix4d>> estimated_T_disks; // Store estimated T_disks for each sample


    // Worker thread: computes estimates off the UI thread and publishes latest result
    std::thread worker([&](){
        for(const auto& sample_info : samples_to_process)
        {
            // Build T_disks for this sample (with interpolation if needed)
            std::vector<Eigen::Matrix4d> T_disks;
            for(unsigned int disk_idx = 0; disk_idx < 7; disk_idx++)
            {
                double x, y, z, rotX, rotY, rotZ;
                
                if(sample_info.alpha > 1e-6)  // Need interpolation
                {
                    interpolateDiskData(data, sample_info.row1, sample_info.row2, sample_info.alpha,
                                       disk_idx, x, y, z, rotX, rotY, rotZ);
                }
                else  // No interpolation needed, use row1 directly
                {
                    extractDiskData(data, sample_info.row1, disk_idx, x, y, z, rotX, rotY, rotZ);
                }
                
                // Build T_disk with the extracted/interpolated data
                Eigen::Matrix4d T_disk = buildDiskTransformFromEulerData(x, y, z, rotX, rotY, rotZ);
                T_disks.push_back(T_disk);
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

            // Store estimated T_disks as a vector for export
            std::vector<Eigen::Matrix4d> estimated_disks_vec;
            for(unsigned int d = 0; d < 7; d++)
            {
                estimated_disks_vec.push_back(T_disk_estimated.block(d*4,0,4,4));
            }
            estimated_T_disks.push_back(estimated_disks_vec);

            // Compare estimate to measurement and store accuracy
            Eigen::MatrixXd acc;
            acc.resize(7,2); // 7 disks, each with 2 error components (pos and rot)
            for(unsigned int d = 0; d < 7; d++)
            {
                Eigen::Matrix4d T_est = T_disk_estimated.block(d*4,0,4,4);
                Eigen::Matrix4d T_meas = T_disks.at(d);

                Eigen::Matrix4d T_err = invert_transformation(T_est) * T_meas;
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
                shared.latest_T_disks = T_disks;
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

        // Export results to CSV file
        // define filename from input file (output file = input file_GP_estimates)
        std::string output_filename = file_name.substr(0, file_name.find_last_of('/')) + "/GP_estimates.csv";
        exportResultsToCSV(output_filename, samples_to_process, data, estimated_T_disks, computation_times);

    });

    // Timer callback: if worker published a new state, update visualizer and render
    vtkSmartPointer<vtkCallbackCommand> timerCallback = vtkSmartPointer<vtkCallbackCommand>::New();
    timerCallback->SetClientData(&shared);
    timerCallback->SetCallback([](vtkObject* caller, unsigned long, void* clientData, void*){
        SharedState* s = static_cast<SharedState*>(clientData);
        if(s->new_state.load())
        {
            ContinuumRobotStateEstimator::SystemState state_copy;
            std::vector<Eigen::Matrix4d> T_disks_copy;
            {
                std::lock_guard<std::mutex> lk(s->mtx);
                state_copy = s->latest_state;
                T_disks_copy = s->latest_T_disks;
                s->new_state.store(false);
            }
            // update visualizer and render
            if(s->vis)
            {
                s->vis->update(state_copy, true, true, 3, &T_disks_copy);
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




