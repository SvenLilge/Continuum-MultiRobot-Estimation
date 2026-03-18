#ifndef CONTINUUM_ROBOT_STATE_ESTIMATOR_H
#define CONTINUUM_ROBOT_STATE_ESTIMATOR_H

#include <utility>
#include <vector>
#include <Eigen/Core>
#include <Eigen/SparseCholesky>

// Main estimator for one or more continuum robots on SE(3).
//
// Workflow:
// 1) Define RobotTopology, Hyperparameters, and Options.
// 2) Construct estimator and provide sensor measurements.
// 3) Call computeStateEstimate(...) to get means + uncertainties.
// 4) Optionally query additional arclength states via interpolation.
class ContinuumRobotStateEstimator
{
public:
    // Geometric/system configuration that remains fixed during estimation.
    struct RobotTopology
    {
        struct Coupling
        {

            unsigned int idxA; // ID of continuum robot (A) being coupled (0 to N-1)
            unsigned int idxB; // ID of continuum robot (B) to which robot (A) is coupled to (0 to N, where N is a common end-effector, e.g. a platform)
            unsigned int coupling_node_robot_A; // Estimation node of robot (A) considered for the coupling
            unsigned int coupling_node_robot_B; // Estimation node of robot (B) to which the node of robot (A) is coupled
            Eigen::Matrix4d T_bA_c; // Transformation between the body frame of robot (A) and coupling frame c
            Eigen::Matrix4d T_bB_c; // Transformation between the body frame of robot (B) (or common end-effector) and coupling frame c
            Eigen::Matrix<int,6,1> mask; // Mask for valid coupling components, first three entries are position/translational constraint, last three entries are orientation constraint (e.g. 1,1,1,0,0,0 is only position constraint)
        };

        unsigned int N; // Number of robots
        std::vector<unsigned int> K; // Number of estimation node states along each robot (must be dimension N)
        std::vector<unsigned int> M; // Number of interpolated notes between estimation nodes (must be dimension N)
        std::vector<double> L; // Total arclength of each robot (must be dimension N)
        std::vector<Eigen::Matrix4d> Ti0; // Base frame of each robot T_i0 (i is static/interial frame, while 0 is body frame b of the first node)

        bool common_end_effector; // Indicates whether the topology features a common end-effector


        // These values can be set to enforce boundary conditions
        // If set to true, the first and/or last poses/strains of each robot will be locked to the initial guess
        std::vector<bool> lock_first_pose;
        std::vector<bool> lock_last_pose;
        std::vector<bool> lock_first_strain;
        std::vector<bool> lock_last_strain;

        // If the robots are equipped with FBG sensors, the following parameters can be set
        // They indicate the distance of the central core of the fiber to the remaining cores
        // and the offset angle of the fiber w.r.t. the robot's backbone
        std::vector<double> fbg_theta_offset;
        std::vector<double> fbg_core_distance;

        std::vector<Coupling> robot_coupling; // Information on robot coupling

    };

    // Probabilistic tuning parameters (measurement and process covariances).
    // Ideally one should be able to set those individually for each cost term
    // For now, the same covariance is used for all measurements/cost term of one type
    struct Hyperparameters 
    {
        // Covanriance matrices for pose and strain measurements
        Eigen::Matrix<double,6,6> R_pose;
        Eigen::Matrix<double,6,6> R_strain;
        Eigen::Matrix<double,4,4> R_fbg_strain;

        // Covariance matrix for coupling constraints
        Eigen::Matrix<double,6,6> R_coupling;

        // Covariance matrix for the prior process noise
        Eigen::Matrix<double,6,6> Qc;
    };

    // Full estimated system state (mean + uncertainty) for all robots.
    struct SystemState
    {
        struct RobotState
        {
            struct Node
            {
                double arclength;

                // State mean
                Eigen::Matrix4d pose; //Expressed in inertial frame T_ik
                Eigen::Matrix<double,6,1> strain;


                // Standard deviation
                Eigen::Matrix<double,6,1> pose_std;
                Eigen::Matrix<double,6,1> strain_std;

                //Covariance matrices
                Eigen::Matrix3d position_covariance;
                Eigen::Matrix3d orientation_covariance;
                Eigen::Matrix3d nu_covariance;
                Eigen::Matrix3d omega_covariance;

            };

            // Nodes for estimation and interpolation
            std::vector<Node> estimation_nodes;
            std::vector<Node> interpolation_nodes;
            std::vector<Node> queried_nodes;
        };

        struct EndEffectorState
        {
            Eigen::Matrix4d pose; //Expressed in inertial frame T_iee

            Eigen::Matrix<double,6,1> pose_std;

            Eigen::Matrix3d position_covariance;
            Eigen::Matrix3d orientation_covariance;
        };

        EndEffectorState end_effector;
        std::vector<RobotState> robots;
    };

    // Numerical settings that control optimizer behavior.
    struct Options
    {
         // Choice of the initial guess of the optimization problem
         // (straight robots, last known system state - if applicable, custom guess - e.g. we can use this to get our guess from strain measurements)
        enum InitialGuessType {Straight, Last, PriorMean, Custom};
        enum Solver {Newton, NewtonLineSearch}; // Solver for the optimization problem

        InitialGuessType init_guess_type;
        Solver solver;
        SystemState custom_guess; // Custom initial guess for solving the system (Expressed with T_ib frames)
        unsigned int max_optimization_iterations;
        double convergence_threshold;
        bool kirchhoff_rods;
    };

    // One sensor input term used by the optimization.
    struct SensorMeasurement //Expressed with T_ib frames
    {
        enum Type {Pose, Strain, FBGStrain};

        Type type; // Pose or strain measurement
        Eigen::MatrixXd value; // 4x4 Transformation matrix or 6x1 strain vector
        // Mask for valid measurement components, first three entries are position/translational strain,
        // last three entries are orientation/rotational strain (e.g. 1,1,1,0,0,0 is only position/translational strain)
        // Ignroed for FBG strain measurements
        Eigen::Matrix<int,6,1> mask; 
        unsigned int idx_robot; // ID of continuum robot (or end-effector) the measurement belongs to (0 to N, where N is a common end-effector, e.g. a platform)
        unsigned int idx_node; // ID of estimatation node the measurement belongs to (ignored if end-effector)
    };

    // Control input for a specific robot segment, used in the GP prior.
    // When type is None, the standard WNOA prior is recovered.
    struct ControlInput
    {
        enum Type {None, Constant, PiecewiseLinear};

        Type type = None;
        std::vector<Eigen::MatrixXd> values = {}; // 12x1 vectors (velocity + acceleration in se(3))

        int idx_robot = -1;   // Robot index (0 to N-1)
        int idx_segment = -1; // Segment index (0 to K-2)
    };

    // Empty constructor; configuration must be provided before estimation.
    ContinuumRobotStateEstimator();
    // Fully configured constructor for immediate use.
    ContinuumRobotStateEstimator(RobotTopology topology, Hyperparameters parameters, Options options);

    // Set/get static robot geometry and coupling configuration.
    void setRobotTopology(RobotTopology topology);
    RobotTopology getRobotTopology();

    // Set/get covariance and noise tuning parameters.
    void setHyperparameters(Hyperparameters parameters);
    Hyperparameters getHyperparameters();

    // Set/get solver options.
    void setOptions(Options options);
    Options getOptions();

    //Get system state expressed with T_ib frames
    SystemState getSystemState();

    //Prints the mean of a state handed to the function
    void printStateMean(SystemState state);

    //Prints the data of a particular robot node handed to the function
    void printNodeInfo(SystemState::RobotState::Node node);

    // Computes the state estimate from a batch of measurements.
    // Returns true if optimization converged; false otherwise.
    // cost stores the optimization objective value history per iteration.
    //Set verbose to true for additional terminal outputs (useful for debugging etc)
    bool computeStateEstimate(SystemState &state, std::vector<double> &cost, std::vector<SensorMeasurement> measurements, std::vector<ControlInput> inputs = {}, bool verbose_mode = false);

    //Returns the state estimate of the last estimation computation with additional queried nodes
    //Careful: Will use the last known state (i.e. the last state computed and returned from computeStateEstimate)
    void queryAdditionalStates(SystemState &state, std::vector<std::pair<unsigned int,double>> arclengths, std::vector<ControlInput> inputs = {});


private:
    RobotTopology m_robot_topology;
    Hyperparameters m_hyperparameters;
    SystemState m_state; //Expressed with T_ib frames
    Options m_options;

    Eigen::SparseMatrix<double> m_P; // Projection matrix

    Eigen::SparseMatrix<double> m_covariance; // Covariance matrix


    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> m_solver;

    //Validates the parameters within each structure and ensures they are set correctly
    void validateRobotTopology(RobotTopology topology);
    void validateHyperparameters(Hyperparameters parameters);
    void validateOptions(Options options);
    void validateMeasurements(std::vector<SensorMeasurement> measurements);

    //Returns a system state to be used as an initial guess in the optimization based on the chosen type
    //Expressed with T_ib frames
    SystemState constructInitialGuess(Options::InitialGuessType type, std::vector<ControlInput> inputs = {});

    // Cached precomputed transition function values per (robot, segment) pair.
    // Populated by getTransitionFunction and reused by getTransitionFunctionIntegral/getQInv.
    struct PrecomputedValues
    {
        struct TransitionFunction {
            Eigen::Matrix<double,12,12> phi_total;
            std::vector<Eigen::Matrix<double,12,12>> phi_segments;
            std::vector<Eigen::Matrix<double,12,12>> phi_segments_transpose;
        };
        struct TransitionFunctionIntegral {
            Eigen::Matrix<double,12,1> phi_integral_total;
            std::vector<Eigen::Matrix<double,12,1>> phi_integral_segments;
        };
        struct CovarianceMatrix {
            Eigen::Matrix<double,12,12> inverse_covariance_total;
            std::vector<Eigen::Matrix<double,12,12>> covariance_segments;
        };
        TransitionFunction trans;
        TransitionFunctionIntegral trans_int;
        CovarianceMatrix Q;
        int idx_robot;
        int idx_segment;
    };

    std::vector<PrecomputedValues> m_precomputed_values;
    Eigen::Matrix<double,12,12> m_LQLT; // = L * Qc * L^T where L = [0; I]

    // Transition function and covariance computation, parameterized by control input.
    // When ControlInput::None, these return the standard WNOA expressions.
    Eigen::MatrixXd getTransitionFunction(ControlInput input, double delta_s);
    Eigen::MatrixXd getTransitionFunctionIntegral(ControlInput input, double delta_s);
    Eigen::MatrixXd getQInv(ControlInput input, double delta_s);
    Eigen::MatrixXd getTransitionFunctionPartial(ControlInput input, double t_k, double t_k1, double delta_s);
    Eigen::MatrixXd getTransitionFunctionIntegralPartial(ControlInput input, double tau, double delta_s);
    Eigen::MatrixXd getQPartial(ControlInput input, double tau, double delta_s);

    //Returns matrix A, vector b and the cost for the prior, coupling and measurement terms based on robot topology, current state and measurements
    void assemblePriorTerms(std::vector<Eigen::Triplet<double>> &A_tripletList, std::vector<Eigen::Triplet<double>> &b_tripletList, double &cost, SystemState state, std::vector<ControlInput> inputs);
    void assembleCouplingTerms(std::vector<Eigen::Triplet<double>> &A_tripletList, std::vector<Eigen::Triplet<double>> &b_tripletList, double &cost, SystemState state);
    void assembleMeasurementTerms(std::vector<Eigen::Triplet<double>> &A_tripletList, std::vector<Eigen::Triplet<double>> &b_tripletList, double &cost, SystemState state, std::vector<SensorMeasurement> measurements, std::vector<ControlInput> inputs);

    //Only computes the prior, coupling and measurement costs based on robot topology, current state and measurements
    //Used for linesearch
    double getPriorCost(SystemState state, std::vector<ControlInput> inputs);
    double getCouplingCost(SystemState state);
    double getMeasurementCost(SystemState state, std::vector<SensorMeasurement> measurements, std::vector<ControlInput> inputs);

    // FBG Sensor model and derivative
    Eigen::Matrix<double,4,1> computeFBGSensorModel(Eigen::Matrix<double,6,1> curvature_strains, double theta_offset, double core_distance);
    Eigen::Matrix<double,4,6> computeFBGSensorModelDerivative(Eigen::Matrix<double,6,1> curvature_strains, double theta_offset, double core_distance);

    //Constructs the projection matrix of the system based on robot toplogy
    Eigen::SparseMatrix<double> constructProjectionMatrix();

    //Solves the linear system in each iteration Ax=b, while considering the projection matrix M and the chosen solving method
    Eigen::MatrixXd solveLinearSystem(Eigen::SparseMatrix<double> A, Eigen::SparseMatrix<double> b, Eigen::SparseMatrix<double> P, bool initialize);

    //Updates the State Variables based on dx
    void updateStateVariables(SystemState &state, Eigen::MatrixXd dx);
    //Updates the State Uncertainties based on system's covariance matrix
    void updateStateUncertainties(SystemState &state, Eigen::SparseMatrix<double> covariance);

    //Interpolate between the estimation nodes
    void interpolateStates(SystemState &state, Eigen::SparseMatrix<double> covariance, std::vector<ControlInput> inputs);

    //Converts a state (only the mean values) expressed in T_ib frames to state expressed in T_bi frames and vice versa
    void convertStateMeanBodyInertial(SystemState &state);

    // Diagnostic helper for optimization debugging.
    // Prints where the system matrix A has non-zero structure (X) vs zeros (-).
    // This is useful to verify that factors are connected as expected:
    // - prior terms should form banded local connections along each robot,
    // - measurement terms should add local diagonal blocks,
    // - coupling terms should add off-diagonal links between coupled states.
    // Note: this is not for visualization of robot shape; it is for checking
    // correctness and conditioning of the linear system assembly.
    void printSparsity(Eigen::MatrixXd A);
};

#endif // CONTINUUM_ROBOT_STATE_ESTIMATOR_H
