#ifndef CONTINUUM_ROBOT_STATE_ESTIMATOR_H
#define CONTINUUM_ROBOT_STATE_ESTIMATOR_H

#include <utility>
#include <vector>
#include <Eigen/Core>
#include <Eigen/SparseCholesky>

class ContinuumRobotStateEstimator
{
public:
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
        std::vector<bool> lock_first_pose;
        std::vector<bool> lock_last_pose;
        std::vector<bool> lock_first_strain;
        std::vector<bool> lock_last_strain;
        std::vector<double> fbg_theta_offset; //TODO: Is this the best spot to store this?
        std::vector<double> fbg_core_distance; //TODO: Is this the best spot to store this?
        bool common_end_effector;
        std::vector<Coupling> robot_coupling; // Information on robot coupling
    };

    struct Hyperparameters
    {
        //TODO: Move these to the individual measurement and coupling structures, such that each can have an independent uncertainty/covariance
        Eigen::Matrix<double,6,6> R_pose;
        Eigen::Matrix<double,6,6> R_strain;
        Eigen::Matrix<double,4,4> R_fbg_strain;

        Eigen::Matrix<double,6,6> R_coupling;

        Eigen::Matrix<double,6,6> Qc;
    };

    struct SystemState
    {
        struct RobotState
        {
            struct Node
            {
                Eigen::Matrix4d pose; //Expressed in inertial frame T_ik
                Eigen::Matrix<double,6,1> strain;

                Eigen::Matrix<double,6,1> pose_std;
                Eigen::Matrix<double,6,1> strain_std;

                double arclength;

                //Save whole covariance
                Eigen::Matrix3d position_covariance;
                Eigen::Matrix3d orientation_covariance;
                Eigen::Matrix3d nu_covariance;
                Eigen::Matrix3d omega_covariance;

            };

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

    struct Options
    {
        enum InitialGuessType {Straight, Last, Custom}; // Choice of the initial guess of the optimization problem (straight robots, last known system state - if applicable, custom guess - e.g. we can use this to get our guess from strain measurements)
        enum Solver {Newton, NewtonLineSearch}; // Solver for the optimization problem

        InitialGuessType init_guess_type;
        Solver solver;
        SystemState custom_guess; // Custom initial guess for solving the system (Expressed with T_ib frames)
        unsigned int max_optimization_iterations;
        double convergence_threshold;
        bool kirchhoff_rods;
    };

    struct SensorMeasurement //Expressed with T_ib frames
    {
        enum Type {Pose, Strain, FBGStrain};

        Type type; // Pose or strain measurement
        Eigen::MatrixXd value; // 4x4 Transformation matrix or 6x1 strain vector
        Eigen::Matrix<int,6,1> mask; // Mask for valid measurement components, first three entries are position/translational strain, last three entries are orientation/rotational strain (e.g. 1,1,1,0,0,0 is only position/translational strain)

        unsigned int idx_robot; // ID of continuum robot (or end-effector) the measurement belongs to (0 to N, where N is a common end-effector, e.g. a platform)
        unsigned int idx_node; // ID of estimatation node the measurement belongs to (ignored if end-effector)
    };

    ContinuumRobotStateEstimator();
    ContinuumRobotStateEstimator(RobotTopology topology, Hyperparameters parameters, Options options);

    void setRobotTopology(RobotTopology topology);
    RobotTopology getRobotTopology();

    void setHyperparameters(Hyperparameters parameters);
    Hyperparameters getHyperparameters();

    void setOptions(Options options);
    Options getOptions();

    //Get system state expressed with T_ib frames
    SystemState getSystemState();

    //Prints the mean of a state handed to the function
    void printStateMean(SystemState state);

    //Prints the data of a particular robot node handed to the function
    void printNodeInfo(SystemState::RobotState::Node node);

    //Computes the and returns state estimate given a set of sensor measurements
    bool computeStateEstimate(SystemState &state, std::vector<double> &cost, std::vector<SensorMeasurement> measurements, bool verbose_mode = false);

    //Returns the state estimate of the last estimation computation with additional queried nodes
    //Careful: Will use the last known state (i.e. the last state computed and returned from computeStateEstimate)
    void queryAdditionalStates(SystemState &state, std::vector<std::pair<unsigned int,double>> arclengths);

    //TODO: Should probably move to private
    static Eigen::Matrix<double,4,1> computeFBGSensorModel(Eigen::Matrix<double,6,1> curvature_strains, double theta_offset, double core_distance);
    static Eigen::Matrix<double,4,6> computeFBGSensorModelDerivative(Eigen::Matrix<double,6,1> curvature_strains, double theta_offset, double core_distance);

private:
    RobotTopology m_robot_topology;
    Hyperparameters m_hyperparameters;
    SystemState m_state; //Expressed with T_ib frames
    Options m_options;

    Eigen::SparseMatrix<double> m_P; // Projection matrix

    Eigen::SparseMatrix<double> m_covariance; // Covariance matrix


    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> m_solver;

    //Validates the parameters within each structure and ensures they are set correctly
    //TODO: Need to check newly added topology and options parameters
    //TODO: Need to check FBG strain measurements
    void validateRobotTopology(RobotTopology topology);
    void validateHyperparameters(Hyperparameters parameters);
    void validateOptions(Options options);
    void validateMeasurements(std::vector<SensorMeasurement> measurements);

    //Returns a system state to be used as an initial guess in the optimization based on the chosen type
    //Expressed with T_ib frames
    SystemState constructInitialGuess(Options::InitialGuessType type);

    //Returns matrix A, vector b and the cost for the prior, coupling and measurement terms based on robot topology, current state and measurements
    void assemblePriorTerms(std::vector<Eigen::Triplet<double>> &A_tripletList, std::vector<Eigen::Triplet<double>> &b_tripletList, double &cost, SystemState state);
    void assembleCouplingTerms(std::vector<Eigen::Triplet<double>> &A_tripletList, std::vector<Eigen::Triplet<double>> &b_tripletList, double &cost, SystemState state);
    void assembleMeasurementTerms(std::vector<Eigen::Triplet<double>> &A_tripletList, std::vector<Eigen::Triplet<double>> &b_tripletList, double &cost, SystemState state, std::vector<SensorMeasurement> measurements);

    //Only computes the prior, coupling and measurement costs based on robot topology, current state and measurements
    //Used for linesearch
    double getPriorCost(SystemState state);
    double getCouplingCost(SystemState state);
    double getMeasurementCost(SystemState state, std::vector<SensorMeasurement> measurements);

    //Constructs the projection matrix of the system based on robot toplogy
    Eigen::SparseMatrix<double> constructProjectionMatrix();

    //Solves the linear system in each iteration Ax=b, while considering the projection matrix M and the chosen solving method
    Eigen::MatrixXd solveLinearSystem(Eigen::SparseMatrix<double> A, Eigen::SparseMatrix<double> b, Eigen::SparseMatrix<double> P, bool initialize);

    //Updates the State Variables based on dx
    void updateStateVariables(SystemState &state, Eigen::MatrixXd dx);
    //Updates the State Uncertainties based on system's covariance matrix
    void updateStateUncertainties(SystemState &state, Eigen::SparseMatrix<double> covariance);

    //Interpolate between the estimation nodes
    void interpolateStates(SystemState &state, Eigen::SparseMatrix<double> covariance);

    //Converts a state (only the mean values) expressed in T_ib frames to state expressed in T_bi frames and vice versa
    void convertStateMeanBodyInertial(SystemState &state);

    //Prints the sparsity pattern of A
    void printSparsity(Eigen::MatrixXd A);
};

#endif // CONTINUUM_ROBOT_STATE_ESTIMATOR_H
