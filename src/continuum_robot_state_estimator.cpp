#include "continuum_robot_state_estimator.h"
#include "utilities.h"

#include <numeric>
#include <algorithm>
#include <iostream>
#include <Eigen/LU>


ContinuumRobotStateEstimator::ContinuumRobotStateEstimator()
{

}

ContinuumRobotStateEstimator::ContinuumRobotStateEstimator(RobotTopology topology, Hyperparameters parameters, Options options)
{
    //Asserts for inputs to make sure they are valid
    validateRobotTopology(topology);
    m_robot_topology = topology;
    validateHyperparameters(parameters);
    m_hyperparameters = parameters;
    validateOptions(options);
    m_options = options;

    //Constrcut projection matrix
    m_P = constructProjectionMatrix();

    // Initialize System State as straight rods based on topology
    m_state = constructInitialGuess(Options::InitialGuessType::Straight);

}

void ContinuumRobotStateEstimator::setRobotTopology(RobotTopology topology)
{
    // Asserts for inputs to make sure they are valid
    validateRobotTopology(topology);

    m_robot_topology = topology;

    //Constrcut projection matrix
    m_P = constructProjectionMatrix();

    // Initialize System State as straight rods based on topology
    m_state = constructInitialGuess(Options::InitialGuessType::Straight);

}

ContinuumRobotStateEstimator::RobotTopology ContinuumRobotStateEstimator::getRobotTopology()
{
    return m_robot_topology;
}

void ContinuumRobotStateEstimator::setHyperparameters(Hyperparameters parameters)
{
    //Asserts for inputs to make sure they are valid
    validateHyperparameters(parameters);

    m_hyperparameters = parameters;
}

ContinuumRobotStateEstimator::Hyperparameters ContinuumRobotStateEstimator::getHyperparameters()
{
    return m_hyperparameters;
}

void ContinuumRobotStateEstimator::setOptions(Options options)
{
    //Asserts for inputs to make sure they are valid
    validateOptions(options);

    m_options = options;
}

ContinuumRobotStateEstimator::Options ContinuumRobotStateEstimator::getOptions()
{
    return m_options;
}

ContinuumRobotStateEstimator::SystemState ContinuumRobotStateEstimator::getSystemState()
{
    return m_state;
}

void ContinuumRobotStateEstimator::validateRobotTopology(RobotTopology topology)
{
    //Go through every variable in structure and make sure their values are valid
    assert((topology.N > 0) && "Topology needs to feature at least one robot");

    assert((topology.K.size() == topology.N) && "Vector size must be equal to number of robots");
    assert((topology.M.size() == topology.N) && "Vector size must be equal to number of robots");
    assert((topology.L.size() == topology.N) && "Vector size must be equal to number of robots");
    assert((topology.Ti0.size() == topology.N) && "Vector size must be equal to number of robots");
    assert((topology.lock_first_pose.size() == topology.N) && "Vector size must be equal to number of robots");
    assert((topology.lock_first_strain.size() == topology.N) && "Vector size must be equal to number of robots");
    assert((topology.lock_last_pose.size() == topology.N) && "Vector size must be equal to number of robots");
    assert((topology.lock_first_strain.size() == topology.N) && "Vector size must be equal to number of robots");

    for(unsigned int i = 0; i < topology.N; i++)
    {
        validate_transformation_matrix(topology.Ti0[i]);
    }

    bool coupling_to_ee_exists = false;
    for(unsigned int i = 0; i < topology.robot_coupling.size(); i++)
    {
        assert((topology.robot_coupling[i].idxA <= topology.N-1) && "Index needs to be between 0 and N-1");
        assert((topology.robot_coupling[i].idxB <= topology.N) && "Index needs to be between 0 and N");
        if(topology.common_end_effector == false)
        {
            assert((topology.robot_coupling[i].idxB != topology.N) && "Can't be coupled to end-effector if end-effector is set to false");
        }
        else
        {
            if(topology.robot_coupling[i].idxB == topology.N)
            {
                coupling_to_ee_exists = true;
            }
        }
        validate_transformation_matrix(topology.robot_coupling[i].T_bB_c);

        for(unsigned j = 0; j < 6; j++)
        {
            assert((topology.robot_coupling[i].mask(j,0) == 1 || topology.robot_coupling[i].mask(j,0) == 0) && "Mask can only hold entries with 1 (enabled) or 0 (disabled)");
        }

    }

    if(topology.common_end_effector == true)
    {
        assert(coupling_to_ee_exists && "At least one robot has to be coupled to common end-effector");
    }

}

void ContinuumRobotStateEstimator::validateHyperparameters(Hyperparameters parameters)
{
    //Go through every variable in structure and make sure their values are valid
    validate_diagonal_matrix(parameters.R_pose);
    validate_diagonal_matrix(parameters.R_strain);
    validate_diagonal_matrix(parameters.R_fbg_strain);
    validate_diagonal_matrix(parameters.R_coupling);
    validate_diagonal_matrix(parameters.Qc);

}

void ContinuumRobotStateEstimator::validateOptions(Options options)
{
    //Go through every variable in structure and make sure their values are valid
    if(options.init_guess_type == Options::InitialGuessType::Custom)
    {
        if(m_robot_topology.common_end_effector == true)
        {
            validate_transformation_matrix(options.custom_guess.end_effector.pose);
        }

        assert((options.custom_guess.robots.size() == m_robot_topology.N) && "Number of robots in initial guess must match number of robots in topology");

        for(unsigned int i = 0; i < options.custom_guess.robots.size(); i++)
        {
            assert((options.custom_guess.robots[i].estimation_nodes.size() == m_robot_topology.K[i]) && "Number of estimation nodes must match number of estimation nodes in topology");

            for(unsigned int j = 0; j < options.custom_guess.robots[i].estimation_nodes.size(); j++)
            {
                validate_transformation_matrix(options.custom_guess.robots[i].estimation_nodes[j].pose);

                if(j == 0)
                {
                    assert((options.custom_guess.robots[i].estimation_nodes[j].pose - m_robot_topology.Ti0[i]).isZero() && "Pose of first node of each robot in initial guess must be identical to its base frame");
                }
            }
        }

    }

    assert(options.max_optimization_iterations > 0 && "Maximum number of iterations must be greater than 0");

}

void ContinuumRobotStateEstimator::validateMeasurements(std::vector<SensorMeasurement> measurements)
{
    //Go through every variable in structure and make sure their values are valid
    for(unsigned int i = 0; i < measurements.size(); i++)
    {
        assert((measurements[i].idx_robot <= m_robot_topology.N) && "Index needs to be between 0 and N");
        if(measurements[i].idx_robot == m_robot_topology.N)
        {
            assert(m_robot_topology.common_end_effector && "Measurement for end-effector only possible, if end-effector exists in topology");
            assert((measurements[i].type == SensorMeasurement::Type::Pose) && "End-effector can only have a pose measurement");
            validate_transformation_matrix(measurements[i].value);
        }
        else
        {
            assert((measurements[i].idx_node < m_robot_topology.K[measurements[i].idx_robot]) && "Node index must be within number of estimation nodes set in topology");
            if(measurements[i].type == SensorMeasurement::Type::Pose)
            {
                validate_transformation_matrix(measurements[i].value);
            }
            else if(measurements[i].type == SensorMeasurement::Type::Strain)
            {
                assert((measurements[i].value.rows() == 6 && measurements[i].value.cols() == 1) && "Vector of strain measurements has to be of dimension 6x1");
            }
            else if(measurements[i].type == SensorMeasurement::Type::FBGStrain)
            {
                assert((measurements[i].value.rows() == 4 && measurements[i].value.cols() == 1) && "Vector of strain measurements has to be of dimension 4x1");
                assert((topology.fbg_theta_offset.size() == topology.N) && "Need to set FBG theta offset for each robot");
                assert((topology.fbg_core_distance.size() == topology.N) && "Need to set FBG core distance for each robot");
            }
        }

        for(unsigned j = 0; j < 6; j++)
        {
            assert((measurements[i].mask(j,0) == 1 || measurements[i].mask(j,0) == 0) && "Mask can only hold entries with 1 (enabled) or 0 (disabled)");
        }

    }

}

ContinuumRobotStateEstimator::SystemState ContinuumRobotStateEstimator::constructInitialGuess(Options::InitialGuessType type)
{

    SystemState initial_guess;
    initial_guess.robots = {};

    if(type == Options::InitialGuessType::Straight)
    {
        if(m_robot_topology.common_end_effector)
        {
            initial_guess.end_effector.pose = Eigen::Matrix4d::Identity();
        }

        for(unsigned int n = 0; n < m_robot_topology.N; n++)
        {
            SystemState::RobotState robot_state;
            robot_state.estimation_nodes = {};
            robot_state.interpolation_nodes = {};

            for(unsigned int k = 0; k < m_robot_topology.K[n]; k++)
            {
                SystemState::RobotState::Node node;

                //Compute arc-length at node k
                double s = (double)k/((double)(m_robot_topology.K[n]-1))*m_robot_topology.L[n];
                node.arclength = s;

                //Define pose of node k in base frame of robot
                Eigen::Matrix4d T_0k = Eigen::Matrix4d::Identity();
                T_0k(0,3) = s; //Set x-component to arclength of node

                //Transform pose to be expressed in inertial frame
                Eigen::Matrix4d T_ik = m_robot_topology.Ti0[n]*T_0k;
                node.pose = T_ik;

                //Define strain
                Eigen::Matrix<double,6,1> strain;
                strain << 1, 0, 0, 0, 0, 0;
                node.strain = strain;

                robot_state.estimation_nodes.push_back(node);

            }

            initial_guess.robots.push_back(robot_state);

        }

    }
    else if(type == Options::InitialGuessType::Custom)
    {
        initial_guess = m_options.custom_guess;
    }
    else if(type == Options::InitialGuessType::Last)
    {
        initial_guess = m_state;
    }

    return initial_guess;


}

void ContinuumRobotStateEstimator::assemblePriorTerms(std::vector<Eigen::Triplet<double> > &A_tripletList, std::vector<Eigen::Triplet<double> > &b_tripletList, double &cost, ContinuumRobotStateEstimator::SystemState state)
{
    //Save total number of nodes
    int total_nodes = std::accumulate(m_robot_topology.K.begin(),m_robot_topology.K.end(),0);

    //Set cost to zero
    cost = 0;


    // Loop over all robots
    int k_offset = 0;
    for(unsigned int n = 0; n < m_robot_topology.N; n++)
    {
        //Loop over all prior terms between two consecutive estimation nodes along the length of each robot
        for(unsigned int k = 0; k < m_robot_topology.K[n] - 1; k++)
        {
            //Get pose of current and next node
            Eigen::Matrix4d T_k = state.robots[n].estimation_nodes[k].pose;
            Eigen::Matrix4d T_k1 = state.robots[n].estimation_nodes[k+1].pose;

            //Get arclength of current and next node
            double s_k = state.robots[n].estimation_nodes[k].arclength;
            double s_k1 = state.robots[n].estimation_nodes[k+1].arclength;

            //Get strain of current and next node
            Eigen::MatrixXd varpi_k = state.robots[n].estimation_nodes[k].strain;
            Eigen::MatrixXd varpi_k1 = state.robots[n].estimation_nodes[k+1].strain;

            //Compute quantities used below
            Eigen::MatrixXd xi = tran_to_vec(T_k1*invert_transformation(T_k));
            Eigen::MatrixXd J_inv = vec_to_jac_inverse(xi);
            Eigen::MatrixXd T_ad = tran_adjoint(T_k1*invert_transformation(T_k));
            double delta_s = s_k1 - s_k;


            Eigen::MatrixXd Qc_inv = invert_diagonal(m_hyperparameters.Qc);

            //Compute error vector at operating point
            Eigen::Matrix<double, 12,1> e;
            e << xi - delta_s*varpi_k,
                    J_inv*varpi_k1 - varpi_k;

            //Jacobian of the error
            Eigen::Matrix<double,12,24> F;
            F << -J_inv*T_ad, -delta_s*Eigen::Matrix<double,6,6>::Identity(), J_inv, Eigen::Matrix<double,6,6>::Zero(),
                    -0.5*curlyhat(varpi_k1)*J_inv*T_ad, -Eigen::Matrix<double,6,6>::Identity(), 0.5*curlyhat(varpi_k1)*J_inv, J_inv;

            //Covariance for prior terms
            Eigen::Matrix<double,12,12> Q_inv;
            Q_inv << 12/(delta_s*delta_s*delta_s)*Qc_inv, -6/(delta_s*delta_s)*Qc_inv,
                    -6/(delta_s*delta_s)*Qc_inv, 4/(delta_s)*Qc_inv;

            //Assemble terms into correct spot in matrix
            int k_start = 12*k + 12*k_offset;

            Eigen::MatrixXd A_block = F.transpose()*Q_inv*F;
            Eigen::MatrixXd b_block = -1*F.transpose()*Q_inv*e;

            //Save the coefficients and entries in A and b
            int row_block = 0;
            int col_block = 0;
            for(int idx_row = k_start; idx_row < k_start + 24; idx_row++)
            {
                for(int idx_col = k_start; idx_col < k_start + 24; idx_col++)
                {
                    A_tripletList.push_back(Eigen::Triplet<double>(idx_row,idx_col,A_block(row_block,col_block)));

                    col_block++;
                }

                b_tripletList.push_back(Eigen::Triplet<double>(idx_row,0,b_block(row_block,0)));

                col_block = 0;
                row_block++;
            }

            //Save the cost
            cost = cost + 0.5*e.transpose()*Q_inv*e;


        }
        k_offset = k_offset + m_robot_topology.K[n];

    }



}

void ContinuumRobotStateEstimator::assembleCouplingTerms(std::vector<Eigen::Triplet<double> > &A_tripletList, std::vector<Eigen::Triplet<double> > &b_tripletList, double &cost, ContinuumRobotStateEstimator::SystemState state)
{
    // Save total number of nodes
    int total_nodes = std::accumulate(m_robot_topology.K.begin(),m_robot_topology.K.end(),0);

    //Set cost to zero
    cost = 0;

    //Run through all of the defined coupling terms in the robot topology
    for(unsigned int c = 0; c < m_robot_topology.robot_coupling.size(); c++)
    {
        RobotTopology::Coupling coupling = m_robot_topology.robot_coupling[c];

        //Get indices of robots to couple (or end-effector)
        unsigned int n1 = coupling.idxA;
        unsigned int n2 = coupling.idxB;

        //Coupling node for first robot is the last node (tip of robot)
        unsigned int k1 = coupling.coupling_node_robot_A;
        unsigned int k2 = coupling.coupling_node_robot_B;

        //Get the pose of the first robot's node
        Eigen::Matrix4d T_1s = state.robots[n1].estimation_nodes[k1].pose;

        //Get the pose of the second robot's node (or the end-effector platform)
        Eigen::Matrix4d T_2s;
        if(n2 == m_robot_topology.N) //In this case, we are coupling to the end-effector
        {
            T_2s = state.end_effector.pose;
        }
        else
        {
            T_2s = state.robots[n2].estimation_nodes[k2].pose;
        }

        //Get the transformations from each of the nodes to the coupling locations
        Eigen::Matrix4d T_1c = coupling.T_bA_c; // Transformation from the body (node) frame to the coupling frame

        Eigen::Matrix4d T_2c = coupling.T_bB_c; // Transformation from the body (node) frame to the coupling frame

        // This is the current pose error in the coupling (should be identity matrix when coupling constraint is fullfilled)
        // T_c1*T_1s*T_s2*T_2c = T_cc
        Eigen::Matrix4d T_cur = invert_transformation(T_1c)*T_1s*invert_transformation(T_2s)*T_2c;

        //Compute the error from the above
        Eigen::Matrix<double,6,1> e = tran_to_vec(T_cur);

        //Compute the Jacobian of the error (6x24, because we consider to full states)
        Eigen::Matrix<double,6,24> F;
        F << vec_to_jac(e).inverse()*tran_adjoint(invert_transformation(T_1c)), Eigen::Matrix<double,6,6>::Zero(), -vec_to_jac(e).inverse()*tran_adjoint(invert_transformation(T_1c)*T_1s*invert_transformation(T_2s)), Eigen::Matrix<double,6,6>::Zero();

        //Offset for the first robot
        int k_offset_1 = 0;
        for(unsigned int n = 0; n < n1; n++)
        {
            k_offset_1 += m_robot_topology.K[n];
        }

        //Offset for second robot (if not end-effector)
        int k_offset_2 = 0;
        for(unsigned int n = 0; n < n2; n++)
        {
            k_offset_2 += m_robot_topology.K[n];
        }



        //CONSIDER MASKING
        for(unsigned int i = 0; i < 6; i++)
        {
            if(coupling.mask(i,0) == 0)
            {
                e(i,0) = 0;
                F.row(i) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
            }
        }


        Eigen::MatrixXd R_c_inv = invert_diagonal(m_hyperparameters.R_coupling);


        Eigen::MatrixXd A_block = F.transpose()*R_c_inv*F;
        Eigen::MatrixXd b_block = -1*F.transpose()*R_c_inv*e;

        cost = cost + 0.5*e.transpose()*R_c_inv*e;


        if(n2 == m_robot_topology.N) //end-effector case
        {
            int k_start_1 = 12*k1 + 12*k_offset_1;
            int k_start_2 = 12*total_nodes;

            // Top left block belongs to node 1 (on diagonal of A)
            int row_block = 0;
            int col_block = 0;
            for(int idx_row = k_start_1; idx_row < k_start_1 + 12; idx_row++)
            {
                for(int idx_col = k_start_1; idx_col < k_start_1 + 12; idx_col++)
                {
                    A_tripletList.push_back(Eigen::Triplet<double>(idx_row,idx_col,A_block(row_block,col_block)));

                    col_block++;
                }

                //Also update the first half of b
                b_tripletList.push_back(Eigen::Triplet<double>(idx_row,0,b_block(row_block,0)));

                col_block = 0;
                row_block++;
            }

            //Bottom right block (minus strain) belongs to end-effector (on diagonal of A)
            row_block = 12;
            col_block = 12;
            for(int idx_row = k_start_2; idx_row < k_start_2 + 6; idx_row++)
            {
                for(int idx_col = k_start_2; idx_col < k_start_2 + 6; idx_col++)
                {
                    A_tripletList.push_back(Eigen::Triplet<double>(idx_row,idx_col,A_block(row_block,col_block)));

                    col_block++;
                }

                //Also update the second part of b
                b_tripletList.push_back(Eigen::Triplet<double>(idx_row,0,b_block(row_block,0)));

                col_block = 12;
                row_block++;
            }

            //Top right block belongs to rows of node 1 and columns of ee
            row_block = 0;
            col_block = 12;
            for(int idx_row = k_start_1; idx_row < k_start_1 + 12; idx_row++)
            {
                for(int idx_col = k_start_2; idx_col < k_start_2 + 6; idx_col++)
                {
                    A_tripletList.push_back(Eigen::Triplet<double>(idx_row,idx_col,A_block(row_block,col_block)));

                    col_block++;
                }

                col_block = 12;
                row_block++;
            }

            //Bottom left block belongs to columns of node 1 and rows of ee
            row_block = 12;
            col_block = 0;
            for(int idx_row = k_start_2; idx_row < k_start_2 + 6; idx_row++)
            {
                for(int idx_col = k_start_1; idx_col < k_start_1 + 12; idx_col++)
                {
                    A_tripletList.push_back(Eigen::Triplet<double>(idx_row,idx_col,A_block(row_block,col_block)));

                    col_block++;
                }

                col_block = 0;
                row_block++;
            }


        }
        else
        {
            int k_start_1 = 12*k1 + 12*k_offset_1;
            int k_start_2 = 12*k2 + 12*k_offset_2;

            // Top left block belongs to node 1 (on diagonal of A)
            int row_block = 0;
            int col_block = 0;
            for(int idx_row = k_start_1; idx_row < k_start_1 + 12; idx_row++)
            {
                for(int idx_col = k_start_1; idx_col < k_start_1 + 12; idx_col++)
                {
                    A_tripletList.push_back(Eigen::Triplet<double>(idx_row,idx_col,A_block(row_block,col_block)));

                    col_block++;
                }

                //Also update the first half of b
                b_tripletList.push_back(Eigen::Triplet<double>(idx_row,0,b_block(row_block,0)));

                col_block = 0;
                row_block++;
            }

            //Bottom right block belongs to node 2 (on diagonal of A)
            row_block = 12;
            col_block = 12;
            for(int idx_row = k_start_2; idx_row < k_start_2 + 12; idx_row++)
            {
                for(int idx_col = k_start_2; idx_col < k_start_2 + 12; idx_col++)
                {
                    A_tripletList.push_back(Eigen::Triplet<double>(idx_row,idx_col,A_block(row_block,col_block)));

                    col_block++;
                }

                //Also update the second part of b
                b_tripletList.push_back(Eigen::Triplet<double>(idx_row,0,b_block(row_block,0)));

                col_block = 12;
                row_block++;
            }

            //Top right block belongs to rows of node 1 and columns of node 2
            row_block = 0;
            col_block = 12;
            for(int idx_row = k_start_1; idx_row < k_start_1 + 12; idx_row++)
            {
                for(int idx_col = k_start_2; idx_col < k_start_2 + 12; idx_col++)
                {
                    A_tripletList.push_back(Eigen::Triplet<double>(idx_row,idx_col,A_block(row_block,col_block)));

                    col_block++;
                }

                col_block = 12;
                row_block++;
            }

            //Bottom left block belongs to columns of node 1 and rows of node 2
            row_block = 12;
            col_block = 0;
            for(int idx_row = k_start_2; idx_row < k_start_2 + 12; idx_row++)
            {
                for(int idx_col = k_start_1; idx_col < k_start_1 + 12; idx_col++)
                {
                    A_tripletList.push_back(Eigen::Triplet<double>(idx_row,idx_col,A_block(row_block,col_block)));

                    col_block++;
                }

                col_block = 0;
                row_block++;
            }
        }


    }

}

void ContinuumRobotStateEstimator::assembleMeasurementTerms(std::vector<Eigen::Triplet<double> > &A_tripletList, std::vector<Eigen::Triplet<double> > &b_tripletList, double &cost, ContinuumRobotStateEstimator::SystemState state, std::vector<ContinuumRobotStateEstimator::SensorMeasurement> measurements)
{
    // Resize matrices according to robot topology
    int total_nodes = std::accumulate(m_robot_topology.K.begin(),m_robot_topology.K.end(),0);

    //Set cost to zero
    cost = 0;


    //Run through all of the measurements
    for(unsigned int m = 0; m < measurements.size(); m++)
    {
        SensorMeasurement measurement = measurements[m];

        //Measurement is for end-effector
        if(measurement.idx_robot == m_robot_topology.N)
        {
            //Measurement can only be pose (we checked this before when we validate the measurements)
            Eigen::Matrix4d T_des = invert_transformation(measurement.value); // Invert T_des so that it is expressed as T_bi
            Eigen::Matrix4d T_ee = state.end_effector.pose; //Get end-effector pose of state

            Eigen::Matrix<double,6,1> e = tran_to_vec(T_ee*invert_transformation(T_des));
            Eigen::Matrix<double,6,6> F = vec_to_jac(e).inverse();

            //CONSIDER MASKING
            for(unsigned int i = 0; i < 6; i++)
            {
                if(measurement.mask(i,0) == 0)
                {
                    e(i,0) = 0;
                    F.row(i) << 0, 0, 0, 0, 0, 0;
                }
            }

            Eigen::Matrix<double,6,6> R_pose_inv = invert_diagonal(m_hyperparameters.R_pose);

            //Assemble terms into correct spot in matrix
            int k_start = 12*total_nodes;

            Eigen::MatrixXd A_block = F.transpose()*R_pose_inv*F;
            Eigen::MatrixXd b_block = -1*F.transpose()*R_pose_inv*e;

            //Save the coefficients and entries in A and b
            int row_block = 0;
            int col_block = 0;
            for(int idx_row = k_start; idx_row < k_start + 6; idx_row++)
            {
                for(int idx_col = k_start; idx_col < k_start + 6; idx_col++)
                {
                    A_tripletList.push_back(Eigen::Triplet<double>(idx_row,idx_col,A_block(row_block,col_block)));

                    col_block++;
                }

                b_tripletList.push_back(Eigen::Triplet<double>(idx_row,0,b_block(row_block,0)));

                col_block = 0;
                row_block++;
            }

            cost = cost + 0.5*e.transpose()*R_pose_inv*e;

        }
        else //Measurement is for one of the continuum robots
        {
            Eigen::MatrixXd e;
            Eigen::MatrixXd F;
            Eigen::MatrixXd R_inv;

            if(measurement.type == SensorMeasurement::Type::Strain)
            {
                e.resize(6,1);
                F.resize(6,12);
                R_inv.resize(6,6);

                Eigen::Matrix<double,6,1> strain_des = -measurement.value;
                Eigen::Matrix<double,6,1> strain_cur = state.robots[measurement.idx_robot].estimation_nodes[measurement.idx_node].strain;

                e = strain_des - strain_cur;
                F << Eigen::Matrix<double,6,6>::Zero(), -Eigen::Matrix<double,6,6>::Identity();

                R_inv = invert_diagonal(m_hyperparameters.R_strain);

                //CONSIDER MASKING
                for(unsigned int i = 0; i < 6; i++)
                {
                    if(measurement.mask(i,0) == 0)
                    {
                        e(i,0) = 0;
                        F.row(i) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                    }
                }
            }
            else if(measurement.type == SensorMeasurement::Type::Pose)
            {
                e.resize(6,1);
                F.resize(6,12);
                R_inv.resize(6,6);

                Eigen::Matrix4d T_des = invert_transformation(measurement.value); // Invert T_des so that it is expressed as T_bi
                Eigen::Matrix4d T_cur = state.robots[measurement.idx_robot].estimation_nodes[measurement.idx_node].pose; //Get node pose of state

                e = tran_to_vec(T_cur*invert_transformation(T_des));
                F << vec_to_jac(e).inverse(), Eigen::Matrix<double,6,6>::Zero();

                R_inv = invert_diagonal(m_hyperparameters.R_pose);

                //CONSIDER MASKING
                for(unsigned int i = 0; i < 6; i++)
                {
                    if(measurement.mask(i,0) == 0)
                    {
                        e(i,0) = 0;
                        F.row(i) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                    }
                }
            }
            else if(measurement.type == SensorMeasurement::Type::FBGStrain)
            {
                e.resize(4,1);
                F.resize(4,12);
                R_inv.resize(4,4);

                Eigen::Matrix<double,4,1> fbg_strain_des = measurement.value;

                Eigen::Matrix<double,6,1> curvature_strain_cur = -1*state.robots[measurement.idx_robot].estimation_nodes[measurement.idx_node].strain;
                Eigen::Matrix<double,4,1> fbg_strain_cur = computeFBGSensorModel(curvature_strain_cur,m_robot_topology.fbg_theta_offset[measurement.idx_robot],m_robot_topology.fbg_core_distance[measurement.idx_robot]);

                e = fbg_strain_des - fbg_strain_cur;
                F << Eigen::Matrix<double,4,6>::Zero(), computeFBGSensorModelDerivative(curvature_strain_cur,m_robot_topology.fbg_theta_offset[measurement.idx_robot],m_robot_topology.fbg_core_distance[measurement.idx_robot]);

                R_inv = invert_diagonal(m_hyperparameters.R_fbg_strain);

                //CONSIDER MASKING
                for(unsigned int i = 0; i < 4; i++)
                {
                    if(measurement.mask(i,0) == 0)
                    {
                        e(i,0) = 0;
                        F.row(i) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                    }
                }

            }

            //Define indexing depending on robot and node ID
            int k_offset = 0;
            for(unsigned int n = 0; n < measurement.idx_robot; n++)
            {
                k_offset = k_offset + m_robot_topology.K[n];
            }
            int k_start = 12*measurement.idx_node + 12*k_offset;

            Eigen::MatrixXd A_block = F.transpose()*R_inv*F;
            Eigen::MatrixXd b_block = -1*F.transpose()*R_inv*e;

            //Save the coefficients and entries in A and b
            int row_block = 0;
            int col_block = 0;
            for(int idx_row = k_start; idx_row < k_start + 12; idx_row++)
            {
                for(int idx_col = k_start; idx_col < k_start + 12; idx_col++)
                {
                    A_tripletList.push_back(Eigen::Triplet<double>(idx_row,idx_col,A_block(row_block,col_block)));

                    col_block++;
                }

                b_tripletList.push_back(Eigen::Triplet<double>(idx_row,0,b_block(row_block,0)));

                col_block = 0;
                row_block++;
            }

            //Update cost
            cost = cost + 0.5*(e.transpose()*R_inv*e)(0,0);

        }


    }


}

double ContinuumRobotStateEstimator::getPriorCost(ContinuumRobotStateEstimator::SystemState state)
{
    //Set cost to zero
    double cost = 0;


    // Loop over all robots
    for(unsigned int n = 0; n < m_robot_topology.N; n++)
    {
        //Loop over all prior terms between two consecutive estimation nodes along the length of each robot
        for(unsigned int k = 0; k < m_robot_topology.K[n] - 1; k++)
        {
            //Get pose of current and next node
            Eigen::Matrix4d T_k = state.robots[n].estimation_nodes[k].pose;
            Eigen::Matrix4d T_k1 = state.robots[n].estimation_nodes[k+1].pose;

            //Get arclength of current and next node
            double s_k = state.robots[n].estimation_nodes[k].arclength;
            double s_k1 = state.robots[n].estimation_nodes[k+1].arclength;

            //Get strain of current and next node
            Eigen::MatrixXd varpi_k = state.robots[n].estimation_nodes[k].strain;
            Eigen::MatrixXd varpi_k1 = state.robots[n].estimation_nodes[k+1].strain;

            //Compute quantities used below
            Eigen::MatrixXd xi = tran_to_vec(T_k1*invert_transformation(T_k));
            Eigen::MatrixXd J_inv = vec_to_jac_inverse(xi);
            double delta_s = s_k1 - s_k;


            Eigen::MatrixXd Qc_inv = invert_diagonal(m_hyperparameters.Qc);

            //Compute error vector at operating point
            Eigen::Matrix<double, 12,1> e;
            e << xi - delta_s*varpi_k,
                    J_inv*varpi_k1 - varpi_k;


            //Covariance for prior terms
            Eigen::Matrix<double,12,12> Q_inv;
            Q_inv << 12/(delta_s*delta_s*delta_s)*Qc_inv, -6/(delta_s*delta_s)*Qc_inv,
                    -6/(delta_s*delta_s)*Qc_inv, 4/(delta_s)*Qc_inv;

            cost = cost + 0.5*e.transpose()*Q_inv*e;


        }

    }

    return cost;

}

double ContinuumRobotStateEstimator::getCouplingCost(ContinuumRobotStateEstimator::SystemState state)
{

    //Set cost to zero
    double cost = 0;

    //Run through all of the defined coupling terms in the robot topology
    for(unsigned int c = 0; c < m_robot_topology.robot_coupling.size(); c++)
    {
        RobotTopology::Coupling coupling = m_robot_topology.robot_coupling[c];

        //Get indices of robots to couple (or end-effector)
        unsigned int n1 = coupling.idxA;
        unsigned int n2 = coupling.idxB;

        //Coupling node for first robot is the last node (tip of robot)
        unsigned int k1 = coupling.coupling_node_robot_A;
        unsigned int k2 = coupling.coupling_node_robot_B;

        //Get the pose of the first robot's node
        Eigen::Matrix4d T_1s = state.robots[n1].estimation_nodes[k1].pose;

        //Get the pose of the second robot's node (or the end-effector platform)
        Eigen::Matrix4d T_2s;
        if(n2 == m_robot_topology.N) //In this case, we are coupling to the end-effector
        {
            T_2s = state.end_effector.pose;
        }
        else
        {
            T_2s = state.robots[n2].estimation_nodes[k2].pose;
        }

        //Get the transformations from each of the nodes to the coupling locations
        Eigen::Matrix4d T_1c = coupling.T_bA_c; // Transformation from the body (node) frame to the coupling frame

        Eigen::Matrix4d T_2c = coupling.T_bB_c; // Transformation from the body (node) frame to the coupling frame

        // This is the current pose error in the coupling (should be identity matrix when coupling constraint is fullfilled)
        // T_c1*T_1s*T_s2*T_2c = T_cc
        Eigen::Matrix4d T_cur = invert_transformation(T_1c)*T_1s*invert_transformation(T_2s)*T_2c;

        //Compute the error from the above
        Eigen::Matrix<double,6,1> e = tran_to_vec(T_cur);

        //CONSIDER MASKING
        for(unsigned int i = 0; i < 6; i++)
        {
            if(coupling.mask(i,0) == 0)
            {
                e(i,0) = 0;
            }
        }


        Eigen::MatrixXd R_c_inv = invert_diagonal(m_hyperparameters.R_coupling);



        cost = cost + 0.5*e.transpose()*R_c_inv*e;



    }

    return cost;

}

double ContinuumRobotStateEstimator::getMeasurementCost(ContinuumRobotStateEstimator::SystemState state, std::vector<ContinuumRobotStateEstimator::SensorMeasurement> measurements)
{

    //Set cost to zero
    double cost = 0;


    //Run through all of the measurements
    for(unsigned int m = 0; m < measurements.size(); m++)
    {
        SensorMeasurement measurement = measurements[m];

        //Measurement is for end-effector
        if(measurement.idx_robot == m_robot_topology.N)
        {
            //Measurement can only be pose (we checked this before when we validate the measurements)
            Eigen::Matrix4d T_des = invert_transformation(measurement.value); // Invert T_des so that it is expressed as T_bi
            Eigen::Matrix4d T_ee = state.end_effector.pose; //Get end-effector pose of state

            Eigen::Matrix<double,6,1> e = tran_to_vec(T_ee*invert_transformation(T_des));

            //CONSIDER MASKING
            for(unsigned int i = 0; i < 6; i++)
            {
                if(measurement.mask(i,0) == 0)
                {
                    e(i,0) = 0;
                }
            }

            Eigen::Matrix<double,6,6> R_pose_inv = invert_diagonal(m_hyperparameters.R_pose);

            cost = cost + 0.5*e.transpose()*R_pose_inv*e;

        }
        else //Measurement is for one of the continuum robots
        {
            Eigen::MatrixXd e;
            Eigen::MatrixXd R_inv;

            if(measurement.type == SensorMeasurement::Type::Strain)
            {
                e.resize(6,1);
                R_inv.resize(6,6);

                Eigen::Matrix<double,6,1> strain_des = -measurement.value;
                Eigen::Matrix<double,6,1> strain_cur = state.robots[measurement.idx_robot].estimation_nodes[measurement.idx_node].strain;

                e = strain_des - strain_cur;

                R_inv = invert_diagonal(m_hyperparameters.R_strain);

                //CONSIDER MASKING
                for(unsigned int i = 0; i < 6; i++)
                {
                    if(measurement.mask(i,0) == 0)
                    {
                        e(i,0) = 0;
                    }
                }
            }
            else if(measurement.type == SensorMeasurement::Type::Pose)
            {
                e.resize(6,1);
                R_inv.resize(6,6);

                Eigen::Matrix4d T_des = invert_transformation(measurement.value); // Invert T_des so that it is expressed as T_bi
                Eigen::Matrix4d T_cur = state.robots[measurement.idx_robot].estimation_nodes[measurement.idx_node].pose; //Get node pose of state

                e = tran_to_vec(T_cur*invert_transformation(T_des));

                R_inv = invert_diagonal(m_hyperparameters.R_pose);

                //CONSIDER MASKING
                for(unsigned int i = 0; i < 6; i++)
                {
                    if(measurement.mask(i,0) == 0)
                    {
                        e(i,0) = 0;
                    }
                }
            }
            else if(measurement.type == SensorMeasurement::Type::FBGStrain)
            {
                e.resize(4,1);
                R_inv.resize(4,4);

                Eigen::Matrix<double,4,1> fbg_strain_des = measurement.value;

                Eigen::Matrix<double,6,1> curvature_strain_cur = -1*state.robots[measurement.idx_robot].estimation_nodes[measurement.idx_node].strain;
                Eigen::Matrix<double,4,1> fbg_strain_cur = computeFBGSensorModel(curvature_strain_cur,m_robot_topology.fbg_theta_offset[measurement.idx_robot],m_robot_topology.fbg_core_distance[measurement.idx_robot]);

                e = fbg_strain_des - fbg_strain_cur;

                R_inv = invert_diagonal(m_hyperparameters.R_fbg_strain);

                //CONSIDER MASKING
                for(unsigned int i = 0; i < 4; i++)
                {
                    if(measurement.mask(i,0) == 0)
                    {
                        e(i,0) = 0;
                    }
                }
            }


            //Update cost
            cost = cost + 0.5*(e.transpose()*R_inv*e)(0,0);

        }
    }

    return cost;
}

Eigen::SparseMatrix<double> ContinuumRobotStateEstimator::constructProjectionMatrix()
{
    Eigen::MatrixXd P;
    Eigen::MatrixXd Id;
    int num_cols;
    int num_rows;
    bool kirchhoff = m_options.kirchhoff_rods; //Option to also constraint the last strain of each robot (should be disable for the most general assumptions)

    int k_total = 0;
    for(unsigned int n = 0; n < m_robot_topology.N; n++)
    {
        k_total = k_total + m_robot_topology.K[n];
    }



    //Define the size of the projection matrix
    num_cols = 12*k_total;
    //Total states
    num_rows = 12*k_total;

    //For each robot, check if first and/or last pose are locked (thus, won't be updated)
    for(unsigned int i = 0; i < m_robot_topology.N; i++)
    {
        if(m_robot_topology.lock_first_pose[i])
        {
            num_rows = num_rows - 6; //Substract initial pose of robot i (6 states)
        }

        if(m_robot_topology.lock_last_pose[i])
        {
            num_rows = num_rows - 6; //Substract last pose of robot i (6 states)
        }

        if(m_robot_topology.lock_first_strain[i])
        {
            if(kirchhoff)
                num_rows = num_rows - 3; //Substract initial strain of robot i (3 additional states for kirchoff)
            else
                num_rows = num_rows - 6; //Substract initial strain of robot i (6 additional states for general)
        }

        if(m_robot_topology.lock_last_strain[i])
        {
            if(kirchhoff)
                num_rows = num_rows - 3; //Substract initial strain of robot i (3 additional states for kirchoff)
            else
                num_rows = num_rows - 6; //Substract initial strain of robot i (6 additional states for general)
        }


    }

    //If we consider Kirchhoff rods, the translational strain variables of each node are locked won't be updated
    if(kirchhoff)
    {
        num_rows = num_rows - 3*k_total;
    }

    //If we have a common end-effector, include 6 additional states for its pose
    if(m_robot_topology.common_end_effector)
    {
        num_cols = num_cols + 6;
        num_rows = num_rows + 6;
    }



    P.resize(num_rows,num_cols);
    P.setZero();

    Id.resize(num_cols,num_cols);
    Id.setIdentity();

    //Fill in the projection matrix

    //Start with end-effector state, if included
    if(m_robot_topology.common_end_effector)
    {
        //Set the bottom right corner to the identity to include the full pose of end-effector
        P.bottomRightCorner(6,6) = Eigen::Matrix<double,6,6>::Identity();
    }

    int k_offset = 0;
    for(unsigned int n = 0; n < m_robot_topology.N; n++) //Run through all robots
    {
        //For each robot add the corresponding Identity block to P


        int Id_idx = 12*k_offset;

        int P_idx  = 12*k_offset;

        //Check if we are locking translational strains
        if(kirchhoff)
        {
            P_idx = P_idx - 3*k_offset;
        }

        //Check if the first and last poses of previous robots in the topology are locked
        for(unsigned int n_t = 0; n_t < n; n_t++)
        {
            if(m_robot_topology.lock_first_pose[n_t])
            {
                P_idx = P_idx - 6;
            }
            if(m_robot_topology.lock_last_pose[n_t])
            {
                P_idx = P_idx - 6;
            }

            if(m_robot_topology.lock_first_strain[n_t])
            {
                if(kirchhoff)
                    P_idx = P_idx - 3; //Substract initial strain of robot i (3 additional states for kirchoff)
                else
                    P_idx = P_idx - 6; //Substract initial strain of robot i (6 additional states for general)
            }

            if(m_robot_topology.lock_last_strain[n_t])
            {
                if(kirchhoff)
                    P_idx = P_idx - 3; //Substract initial strain of robot i (3 additional states for kirchoff)
                else
                    P_idx = P_idx - 6; //Substract initial strain of robot i (6 additional states for general)
            }
        }

        //Fill in the projection matrix accordingly
        for(unsigned int k = 0; k < m_robot_topology.K.at(n); k++)
        {
            unsigned int p_offset = 0;
            unsigned int id_offset = 0;



            bool lock_strain = ((k == 0 && m_robot_topology.lock_first_strain[n]) || ((k == m_robot_topology.K[n] -1) && m_robot_topology.lock_last_strain[n]));
            bool lock_pose = ((k == 0 && m_robot_topology.lock_first_pose[n]) || ((k == m_robot_topology.K[n] -1) && m_robot_topology.lock_last_pose[n]));

            if(!lock_pose)
            {
                //Update pose
                P.block(P_idx,0,6,num_cols) = Id.block(Id_idx,0,6,num_cols);
                p_offset += 6;
            }

            if(lock_strain)
            {
                //Update nothing
            }
            else if(kirchhoff)
            {
                //Update only rot strain
                P.block(P_idx+p_offset,0,3,num_cols) = Id.block(Id_idx+9,0,3,num_cols);
                p_offset += 3;
            }
            else
            {
                //Update full strain
                P.block(P_idx+p_offset,0,6,num_cols) = Id.block(Id_idx + 6,0,6,num_cols);
                p_offset += 6;
            }
            id_offset = 12;

            P_idx += p_offset;
            Id_idx += id_offset;
        }



        k_offset = k_offset + m_robot_topology.K[n];
    }

    return P.sparseView();

}

Eigen::MatrixXd ContinuumRobotStateEstimator::solveLinearSystem(Eigen::SparseMatrix<double> A, Eigen::SparseMatrix<double> b, Eigen::SparseMatrix<double> P, bool initialize)
{
    //Apply projection matrix
    A = P*A*P.transpose();
    b = P*b;

    if(initialize)
    {
        m_solver.analyzePattern(A);
        m_solver.factorize(A);

        assert(m_solver.vectorD().minCoeff() > 1e-15 && "System matrix is not invertible, state estimation problem not well defined (i.e. not enough constraints or measurements)");
    }
    else
    {
        m_solver.factorize(A);
    }

    Eigen::MatrixXd dx = m_solver.solve(b);



    //Undo projection matrix
    dx = P.transpose()*dx;

    return dx;


}

Eigen::Matrix<double,4,1> ContinuumRobotStateEstimator::computeFBGSensorModel(Eigen::Matrix<double, 6, 1> curvature_strains, double theta_offset, double core_distance)
{
    Eigen::Matrix<double,4,1> y;

    double nu_1 = curvature_strains(0);
    double omega_1 = curvature_strains(3);
    double omega_2 = curvature_strains(4);
    double omega_3 = curvature_strains(5);

    //Compute strain of innermost core
    y(0) = nu_1 - 1;

    //Compute remaining strain variables
    for(unsigned int i = 0; i < 3; i++)
    {
        double theta_i = -1.0*(double)i*2.0*M_PI/3.0 + theta_offset;
        double tmp1 = nu_1 - core_distance*(omega_3*std::cos(theta_i) - omega_2*std::sin(theta_i));
        double tmp2 = core_distance*omega_1;
        double epsilon = std::sqrt(tmp1*tmp1+tmp2*tmp2) - 1;
        y(i+1) = epsilon;
    }

    return y;
}

Eigen::Matrix<double,4,6> ContinuumRobotStateEstimator::computeFBGSensorModelDerivative(Eigen::Matrix<double, 6, 1> curvature_strains, double theta_offset, double core_distance)
{
    Eigen::Matrix<double,4,6> G;


    double nu_1 = curvature_strains(0);
    double omega_1 = curvature_strains(3);
    double omega_2 = curvature_strains(4);
    double omega_3 = curvature_strains(5);

    //Compute derivative for innermost core
    G.row(0) << 1, 0, 0, 0, 0, 0;


    //Compute remaining derivatives
    for(unsigned int i = 0; i < 3; i++)
    {
        double theta_i = -1.0*(double)i*2.0*M_PI/3.0 + theta_offset;

        double tmp1 = nu_1 - core_distance*(omega_3*std::cos(theta_i) - omega_2*std::sin(theta_i));
        double tmp2 = core_distance*omega_1;
        double divisor = std::sqrt(tmp1*tmp1+tmp2*tmp2);

        double d_nu_1 = tmp1/divisor;
        double d_omega_1 = (omega_1*core_distance*core_distance)/divisor;
        double d_omega_2 = core_distance*std::sin(theta_i)*tmp1/divisor;
        double d_omega_3 = -1*core_distance*std::cos(theta_i)*tmp1/divisor;

        G.row(i+1) << d_nu_1, 0, 0, d_omega_1, d_omega_2, d_omega_3;
    }

    return G;
}


void ContinuumRobotStateEstimator::updateStateVariables(ContinuumRobotStateEstimator::SystemState &state, Eigen::MatrixXd dx)
{
    int k_offset = 0;
    for(unsigned int n = 0; n < m_robot_topology.N; n++)
    {
        for(unsigned int k = 0; k < m_robot_topology.K[n]; k++)
        {
            int k_start = 12*k + 12*k_offset;
            Eigen::MatrixXd dx_k = dx.block(k_start,0,12,1);

            state.robots[n].estimation_nodes[k].pose = vec_to_tran(dx_k.block(0,0,6,1))*state.robots[n].estimation_nodes[k].pose;
            state.robots[n].estimation_nodes[k].strain = state.robots[n].estimation_nodes[k].strain + dx_k.block(6,0,6,1);

        }
        k_offset = k_offset + m_robot_topology.K[n];
    }

    if(m_robot_topology.common_end_effector)
    {
        Eigen::MatrixXd dx_ee = dx.bottomRows(6);

        state.end_effector.pose = vec_to_tran(dx_ee)*state.end_effector.pose;
    }
}

void ContinuumRobotStateEstimator::updateStateUncertainties(ContinuumRobotStateEstimator::SystemState &state, Eigen::SparseMatrix<double> covariance)
{
    int k_offset = 0;
    for(unsigned int n = 0; n < m_robot_topology.N; n++)
    {
        for(unsigned int k = 0; k < m_robot_topology.K[n]; k++)
        {
            int k_start = 12*k + 12*k_offset;

            Eigen::Matrix3d rot = state.robots[n].estimation_nodes[k].pose.block(0,0,3,3);

            //We are already transforming in the world frame here (but do it later for the state itself)
            Eigen::Matrix3d position_covariance = rot.transpose()*covariance.block(k_start,k_start,3,3)*rot;
            Eigen::Matrix3d orientation_covariance = rot.transpose()*covariance.block(k_start+3,k_start+3,3,3)*rot;


            Eigen::Matrix3d nu_covariance = covariance.block(k_start+6,k_start+6,3,3);
            Eigen::Matrix3d omega_covariance = covariance.block(k_start+9,k_start+9,3,3);

            // We are saving the position covariance so that we can plot ellipsoids from it
            state.robots[n].estimation_nodes[k].position_covariance = position_covariance;
            state.robots[n].estimation_nodes[k].orientation_covariance = orientation_covariance;
            state.robots[n].estimation_nodes[k].nu_covariance = nu_covariance;
            state.robots[n].estimation_nodes[k].omega_covariance = omega_covariance;

            //Save standard deviations
            Eigen::Matrix<double,6,1> pose_std;
            pose_std << std::sqrt(position_covariance(0,0)),
                    std::sqrt(position_covariance(1,1)),
                    std::sqrt(position_covariance(2,2)),
                    std::sqrt(orientation_covariance(0,0)),
                    std::sqrt(orientation_covariance(1,1)),
                    std::sqrt(orientation_covariance(2,2));

            Eigen::Matrix<double,6,1> strain_std;
            strain_std << std::sqrt(nu_covariance(0,0)),
                    std::sqrt(nu_covariance(1,1)),
                    std::sqrt(nu_covariance(2,2)),
                    std::sqrt(omega_covariance(0,0)),
                    std::sqrt(omega_covariance(1,1)),
                    std::sqrt(omega_covariance(2,2));

            state.robots[n].estimation_nodes[k].pose_std = pose_std;
            state.robots[n].estimation_nodes[k].strain_std = strain_std;

        }
        k_offset = k_offset + m_robot_topology.K[n];
    }

    if(m_robot_topology.common_end_effector)
    {
        Eigen::Matrix3d rot = state.end_effector.pose.block(0,0,3,3);

        //We are already transforming in the world frame here (but do it later for the state itself)

        Eigen::Matrix3d position_covariance = rot.transpose()*covariance.bottomRightCorner(6,6).block(0,0,3,3)*rot;
        Eigen::Matrix3d orientation_covariance = rot.transpose()*covariance.bottomRightCorner(3,3)*rot;

        state.end_effector.position_covariance = position_covariance;
        state.end_effector.orientation_covariance = orientation_covariance;


        Eigen::Matrix<double,6,1> pose_std;
        pose_std << std::sqrt(position_covariance(0,0)),
                std::sqrt(position_covariance(1,1)),
                std::sqrt(position_covariance(2,2)),
                std::sqrt(orientation_covariance(0,0)),
                std::sqrt(orientation_covariance(1,1)),
                std::sqrt(orientation_covariance(2,2));


        state.end_effector.pose_std = pose_std;

    }
}

void ContinuumRobotStateEstimator::interpolateStates(ContinuumRobotStateEstimator::SystemState &state, Eigen::SparseMatrix<double> covariance)
{
    int k_offset = 0;
    for(unsigned int n = 0; n < m_robot_topology.N; n++)
    {
        state.robots[n].interpolation_nodes.clear();


        //Define first node for each robot (copy from estimation node)
        SystemState::RobotState::Node node = state.robots[n].estimation_nodes[0];

        //Push back the node into interpolation node vector
        state.robots[n].interpolation_nodes.push_back(node);

        for(unsigned int k = 0; k < m_robot_topology.K[n] - 1; k++)
        {
            //Define variables that only depend on K

            //Get pose of current and next node
            Eigen::Matrix4d T_k = state.robots[n].estimation_nodes[k].pose;
            Eigen::Matrix4d T_k1 = state.robots[n].estimation_nodes[k+1].pose;

            //Get arclength of current and next node
            double s_k = state.robots[n].estimation_nodes[k].arclength;
            double s_k1 = state.robots[n].estimation_nodes[k+1].arclength;

            //Get strain of current and next node
            Eigen::MatrixXd varpi_k = state.robots[n].estimation_nodes[k].strain;
            Eigen::MatrixXd varpi_k1 = state.robots[n].estimation_nodes[k+1].strain;

            //Compute quantities used below
            Eigen::MatrixXd xi = tran_to_vec(T_k1*invert_transformation(T_k));
            Eigen::MatrixXd J_inv = vec_to_jac_inverse(xi);

            double delta_s = s_k1 - s_k;


            Eigen::MatrixXd Qc = m_hyperparameters.Qc;
            Eigen::MatrixXd Qc_inv = invert_diagonal(m_hyperparameters.Qc);

            //Covariance for prior terms
            Eigen::Matrix<double,12,12> Q_inv;
            Q_inv << 12/(delta_s*delta_s*delta_s)*Qc_inv, -6/(delta_s*delta_s)*Qc_inv,
                    -6/(delta_s*delta_s)*Qc_inv, 4/(delta_s)*Qc_inv;


            Eigen::Matrix<double,12,12> phi;
            phi << Eigen::Matrix<double,6,6>::Identity(), delta_s*Eigen::Matrix<double,6,6>::Identity(),
                    Eigen::Matrix<double,6,6>::Zero(), Eigen::Matrix<double,6,6>::Identity();


            for(unsigned int m = 1; m <= m_robot_topology.M[n]; m++)
            {
                //Create a new node
                SystemState::RobotState::Node node;

                //Interpolation coefficients
                double s_q = s_k +  (double)m/((double)m_robot_topology.M[n])*delta_s;
                node.arclength = s_q; //arclength

                double delta_s_q = s_q - s_k;

                Eigen::MatrixXd Q;
                Q.resize(2*Qc.rows(),2*Qc.cols());

                Q << delta_s_q*delta_s_q*delta_s_q/3*Qc, delta_s_q*delta_s_q/2*Qc,
                        delta_s_q*delta_s_q/2*Qc, delta_s_q*Qc;

                Eigen::Matrix<double,12,12> psi_temp;

                psi_temp << Eigen::Matrix<double,6,6>::Identity(), (s_k1 - s_q)*Eigen::Matrix<double,6,6>::Identity(),
                        Eigen::Matrix<double,6,6>::Zero(), Eigen::Matrix<double,6,6>::Identity();

                Eigen::MatrixXd psi = Q*psi_temp.transpose()*Q_inv;


                Eigen::Matrix<double,12,12> lambda;

                lambda << Eigen::Matrix<double,6,6>::Identity(), delta_s_q*Eigen::Matrix<double,6,6>::Identity(),
                        Eigen::Matrix<double,6,6>::Zero(), Eigen::Matrix<double,6,6>::Identity();
                lambda = lambda - psi*phi;

                //Interpolate the mean
                Eigen::Matrix<double,12,1> gamma_1;
                gamma_1 << Eigen::Matrix<double,6,1>::Zero(),
                        varpi_k;

                Eigen::Matrix<double,12,1> gamma_2;
                gamma_2 << xi,
                        J_inv*varpi_k1;

                Eigen::MatrixXd gamma = lambda*gamma_1 + psi*gamma_2;

                Eigen::Matrix4d T_q = vec_to_tran(gamma.topRows(6)) * T_k;
                node.pose = T_q; //pose

                Eigen::MatrixXd J_tau = vec_to_jac(gamma.topRows(6));

                Eigen::MatrixXd varpi_q = J_tau*gamma.bottomRows(6);
                node.strain = varpi_q; // strain

                //Interpolate the covariance
                Eigen::Matrix<double,12,12> Ga_1;
                Ga_1 << Eigen::Matrix<double,6,6>::Identity(), Eigen::Matrix<double,6,6>::Zero(),
                        0.5*curlyhat(varpi_k),  Eigen::Matrix<double,6,6>::Identity();


                Eigen::Matrix<double,12,12> Ga_2;
                Ga_2 << J_inv, Eigen::Matrix<double,6,6>::Zero(),
                        0.5*curlyhat(varpi_k1)*J_inv, J_inv;


                Eigen::Matrix<double,12,12> Xi_1;
                Xi_1 << Eigen::Matrix<double,6,6>::Identity(), Eigen::Matrix<double,6,6>::Zero(),
                        Eigen::Matrix<double,6,6>::Zero(), Eigen::Matrix<double,6,6>::Zero();

                Eigen::Matrix<double,12,12> Xi_2;
                Xi_2 << tran_adjoint(vec_to_tran(xi)), Eigen::Matrix<double,6,6>::Zero(),
                        Eigen::Matrix<double,6,6>::Zero(), Eigen::Matrix<double,6,6>::Zero();


                int k_start = 12*k + 12*k_offset;

                Eigen::MatrixXd P11 = covariance.block(k_start,k_start,12,12);
                Eigen::MatrixXd P12 = covariance.block(k_start,k_start+12,12,12);
                Eigen::MatrixXd P22 = covariance.block(k_start+12,k_start+12,12,12);

                Eigen::MatrixXd P11k = Ga_1*(P11 - Xi_1*P11*Xi_1.transpose())*Ga_1.transpose();
                Eigen::MatrixXd P12k = Ga_1*(P12 - Xi_1*P11*Xi_2.transpose())*Ga_2.transpose();
                Eigen::MatrixXd P22k = Ga_2*(P22 - Xi_2*P11*Xi_2.transpose())*Ga_2.transpose();

                Eigen::Matrix<double,24,24> Pk;
                Pk << P11k, P12k,
                        P12k.transpose(), P22k;


                Eigen::MatrixXd Q_tau;
                Eigen::Matrix<double,12,12> Q_tau_temp;
                Q_tau_temp << Eigen::Matrix<double,6,6>::Identity(), (s_k1-s_q)*Eigen::Matrix<double,6,6>::Identity(),
                        Eigen::Matrix<double,6,6>::Zero(), Eigen::Matrix<double,6,6>::Identity();

                Q_tau = Q - psi*Q_tau_temp*Q;

                Eigen::Matrix<double,12,24> lambda_psi;
                lambda_psi << lambda, psi;

                Eigen::MatrixXd Pq_local = lambda_psi*Pk*lambda_psi.transpose() + Q_tau;

                Eigen::Matrix<double,12,12> Ga_inv;
                Ga_inv << J_tau, Eigen::Matrix<double,6,6>::Zero(),
                        -0.5*J_tau*curlyhat(varpi_q), J_tau;

                Eigen::Matrix<double,12,12> Xi;
                Xi << tran_adjoint(vec_to_tran(gamma.topRows(6))), Eigen::Matrix<double,6,6>::Zero(),
                        Eigen::Matrix<double,6,6>::Zero(),Eigen::Matrix<double,6,6>::Zero();

                Eigen::MatrixXd covariance_q = Ga_inv*Pq_local*Ga_inv.transpose() + Xi*covariance.block(k_start,k_start,12,12)*Xi.transpose();

                //Save out the relevant variables

                Eigen::Matrix3d rot = node.pose.block(0,0,3,3);

                //We are already transforming in the world frame here (but do it later for the state itself)
                Eigen::Matrix3d position_covariance = rot.transpose()*covariance_q.block(0,0,3,3)*rot;
                Eigen::Matrix3d orientation_covariance = rot.transpose()*covariance_q.block(3,3,3,3)*rot;
                Eigen::Matrix3d nu_covariance = covariance_q.block(6,6,3,3);
                Eigen::Matrix3d omega_covariance = covariance_q.block(9,9,3,3);

                // We are saving the position covariance so that we can plot ellipsoids from it
                node.position_covariance = position_covariance;
                node.orientation_covariance = orientation_covariance;
                node.nu_covariance = nu_covariance;
                node.omega_covariance = omega_covariance;

                //Save standard deviations
                Eigen::Matrix<double,6,1> pose_std;
                pose_std << std::sqrt(position_covariance(0,0)),
                        std::sqrt(position_covariance(1,1)),
                        std::sqrt(position_covariance(2,2)),
                        std::sqrt(orientation_covariance(0,0)),
                        std::sqrt(orientation_covariance(1,1)),
                        std::sqrt(orientation_covariance(2,2));

                Eigen::Matrix<double,6,1> strain_std;
                strain_std << std::sqrt(nu_covariance(0,0)),
                        std::sqrt(nu_covariance(1,1)),
                        std::sqrt(nu_covariance(2,2)),
                        std::sqrt(omega_covariance(0,0)),
                        std::sqrt(omega_covariance(1,1)),
                        std::sqrt(omega_covariance(2,2));

                node.pose_std = pose_std;
                node.strain_std = strain_std;

                //Push back the node
                state.robots[n].interpolation_nodes.push_back(node);


            }



        }
        k_offset = k_offset + m_robot_topology.K[n];
    }

}

//Converts the state from body to inertial frame and vice versa
//Careful: Only done for the means, not the covariances
void ContinuumRobotStateEstimator::convertStateMeanBodyInertial(ContinuumRobotStateEstimator::SystemState &state)
{
    //Invert all the frames from T_ib to T_bi or vice versa
    //Change the sign of the strain
    if(m_robot_topology.common_end_effector)
    {
        state.end_effector.pose = invert_transformation(state.end_effector.pose);
    }

    //Run through all robots
    for(unsigned int n = 0; n < m_robot_topology.N; n++)
    {
        //Run through all estimation nodes
        for(unsigned int k = 0; k < m_robot_topology.K[n]; k++)
        {
            state.robots[n].estimation_nodes[k].pose = invert_transformation(state.robots[n].estimation_nodes[k].pose);
            state.robots[n].estimation_nodes[k].strain = -state.robots[n].estimation_nodes[k].strain;
        }

        //Run through all interpolation nodes
        for(unsigned int m = 0; m < state.robots[n].interpolation_nodes.size(); m++)
        {

            state.robots[n].interpolation_nodes[m].pose = invert_transformation(state.robots[n].interpolation_nodes[m].pose);
            state.robots[n].interpolation_nodes[m].strain = -state.robots[n].interpolation_nodes[m].strain;
        }

        //Run through all queried nodes
        for(unsigned int q = 0; q < state.robots[n].queried_nodes.size(); q++)
        {

            state.robots[n].queried_nodes[q].pose = invert_transformation(state.robots[n].queried_nodes[q].pose);
            state.robots[n].queried_nodes[q].strain = -state.robots[n].queried_nodes[q].strain;
        }
    }

}

void ContinuumRobotStateEstimator::printSparsity(Eigen::MatrixXd A)
{
    std::cout  << std::endl << "Sparsity pattern of A (each entry is a 6x6 block):" << std::endl;
    for(unsigned int i = 0; i < A.rows()/6; i++)
    {
        for(unsigned int j = 0; j < A.cols()/6; j++)
        {
            Eigen::MatrixXd block = A.block(6*i,6*j,6,6);

            if(block.isZero())
            {
                std::cout << "- ";

            }
            else
            {
                std::cout << "X ";
            }

        }
        std::cout << std::endl;
    }
    std::cout  << std::endl;
}

void ContinuumRobotStateEstimator::printStateMean(ContinuumRobotStateEstimator::SystemState state)
{


    //Run through all robots
    for(unsigned int n = 0; n < m_robot_topology.N; n++)
    {
        std::cout << "Robot " << n+1 << ":" << std::endl << std::endl;

        //Run through all estimation nodes
        for(unsigned int k = 0; k < m_robot_topology.K[n]; k++)
        {

            std::cout << "Node " << k+1 << ":" << std::endl;

            std::cout << "Arclength:" << std::endl;

            std::cout <<state.robots[n].estimation_nodes[k].arclength << std::endl;

            std::cout << "Pose:" << std::endl;

            std::cout <<state.robots[n].estimation_nodes[k].pose << std::endl;

            std::cout << "Strain:" << std::endl;

            std::cout <<state.robots[n].estimation_nodes[k].strain.transpose() << std::endl << std::endl;
        }
    }

    if(m_robot_topology.common_end_effector)
    {
        std::cout << "End-Effector:" << std::endl << std::endl;

        std::cout << state.end_effector.pose << std::endl << std::endl;
    }

}


void ContinuumRobotStateEstimator::printNodeInfo(ContinuumRobotStateEstimator::SystemState::RobotState::Node node)
{
    std::cout << "Arclegnth:" << std::endl;

    std::cout <<node.arclength << std::endl;

    std::cout << "Pose:" << std::endl;

    std::cout <<node.pose << std::endl;

    std::cout << "Strain:" << std::endl;

    std::cout <<node.strain.transpose() << std::endl;

    std::cout << "Pose StD:" << std::endl;

    std::cout <<node.pose_std.transpose() << std::endl;

    std::cout << "Strain StD:" << std::endl;

    std::cout <<node.strain_std.transpose() << std::endl;

    std::cout << "Position Covariance:" << std::endl;

    std::cout <<node.position_covariance << std::endl;

    std::cout << "Orientation Covariance:" << std::endl;

    std::cout <<node.orientation_covariance << std::endl;

    std::cout << "Nu Covariance:" << std::endl;

    std::cout <<node.nu_covariance << std::endl;

    std::cout << "Omega Covariance:" << std::endl;

    std::cout <<node.omega_covariance << std::endl << std::endl;

}

bool ContinuumRobotStateEstimator::computeStateEstimate(SystemState &state, std::vector<double> &cost, std::vector<SensorMeasurement> measurements, bool verbose_mode)
{
    // Asserts for inputs to make sure they are valid
    validateMeasurements(measurements);

    // Setup the initial guess depending on chosen option
    SystemState current_state = constructInitialGuess(m_options.init_guess_type);

    // Convert the state to be expressed with T_bi frames
    // This is necessary as all of our derivations throughout the paper are based on this convention
    convertStateMeanBodyInertial(current_state);

    //Clear cost
    cost.clear();



    // Create and resize matrices according to robot topology
    int total_nodes = std::accumulate(m_robot_topology.K.begin(),m_robot_topology.K.end(),0);

    int end_effector_state = 0;
    if(m_robot_topology.common_end_effector)
    {
        end_effector_state = 6;
    }

    Eigen::SparseMatrix<double> A(12*total_nodes+end_effector_state,12*total_nodes+end_effector_state);

    Eigen::SparseMatrix<double> b(12*total_nodes+end_effector_state,1);


    //Create the triplet lists
    std::vector<Eigen::Triplet<double>> A_tripletList;
    std::vector<Eigen::Triplet<double>> b_tripletList;

    // Solve the system iteratively
    bool max_iter = false;
    for(unsigned int iter = 0; iter < m_options.max_optimization_iterations; iter++)
    {

        //Clear triplet lists
        A_tripletList.clear();
        b_tripletList.clear();


        // Assemble prior terms
        double cost_p;
        assemblePriorTerms(A_tripletList,b_tripletList,cost_p,current_state);

        // Assemble coupling terms
        double cost_c;

        assembleCouplingTerms(A_tripletList,b_tripletList,cost_c,current_state);

        // Assemble measurement terms
        double cost_m;

        assembleMeasurementTerms(A_tripletList,b_tripletList,cost_m,current_state,measurements);



        //Add terms together
        A.setFromTriplets(A_tripletList.begin(), A_tripletList.end());
        b.setFromTriplets(b_tripletList.begin(), b_tripletList.end());

        //Save cost for current iteration
        cost.push_back(cost_p + cost_c + cost_m);

        //Print information about current iteration
        if(verbose_mode)
        {
            if(iter == 0)
            {
                printSparsity(A);
            }
            std::cout << "Iteration: " << iter << std::endl;
            std::cout << "Cost: " << cost.back() << std::endl;
            std::cout << "Cost Prior: " << cost_p << std::endl;
            std::cout << "Cost Coupling: " << cost_c << std::endl;
            std::cout << "Cost Measurements: " << cost_m << std::endl;
            std::cout << std::endl;
        }

        // Solve for dx
        bool initialize = false;
        if(iter == 0)
        {
            initialize = true;
        }

        Eigen::MatrixXd dx = solveLinearSystem(A,b,m_P, initialize);

        Eigen::MatrixXd p = dx;
        double m;
        m = (-dx.transpose()*p)(0,0);
        double alpha = 1;


        // Update the state variables
        if(m_options.solver == Options::Solver::NewtonLineSearch) //Do linesearch if selected
        {
            SystemState state_check;
            double new_cost;

            Eigen::MatrixXd step;
            double c = 0.5;
            double tau = 0.5;
            int ls_iter = 0;
            int max_ls_iter = 10;

            do// As long as the cost of the next step are larger than the current cost
            {
                //Reset state_check to current state
                state_check = current_state;

                //Define step
                step = alpha*p;

                //Apply current step to state_check
                updateStateVariables(state_check,step);

                //Compute the cost for the new state
                new_cost = getPriorCost(state_check) + getMeasurementCost(state_check,measurements) + getCouplingCost(state_check);
                
                alpha = tau*alpha; //Half the step to make for next iteration

                ls_iter++;



            } while((new_cost > (cost.back() + m*c*alpha) || std::isnan(new_cost)) && ls_iter < max_ls_iter);

            //Set the current state to the state that minimzed the cost
            current_state = state_check;

        }
        else //Otherwise just apply the whole step
        {
            updateStateVariables(current_state,dx);
        }



        //Check if the cost during the last couple iterations are still changing
        bool detect_convergence = false;

        if(m > -1*m_options.convergence_threshold)
        {
            detect_convergence = true;
        }

        if(detect_convergence)
        {
            if(verbose_mode)
            {
                std::cout << "Detected convergence of cost function." << std::endl << std::endl;
            }
            break;
        }

        if(iter == m_options.max_optimization_iterations - 1)
        {
            if(verbose_mode)
            {
                std::cout << "Maximum number of iterations reached." << std::endl << std::endl;
            }
            max_iter = true;
        }
    }


    // Extract covariances and uncertainties and update state
    Eigen::SparseMatrix<double> tmp_covariance = m_P*A*m_P.transpose();
    Eigen::SparseMatrix<double> Identity(tmp_covariance.rows(), tmp_covariance.cols());
    Identity.setIdentity();

    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double> > solver;
    Eigen::SparseMatrix<double> tmp_covariance_inv = solver.compute(tmp_covariance).solve(Identity);

    m_covariance = m_P.transpose()*tmp_covariance_inv*m_P; //Save the covariance of the estimate at the final iteration
    updateStateUncertainties(current_state, m_covariance);

    //Fill in the interpolation nodes between the estimation nodes
    interpolateStates(current_state, m_covariance);


    //Convert the state mean back to inertial frame (more intuitive to work with)
    convertStateMeanBodyInertial(current_state);

    // Set member variable m_state and return it
    m_state = current_state;
    state = current_state;


    if(max_iter)
    {
        return false;
    }
    return true;
}

//  Careful: uses the last saved system state
//  Need to solve for state estimate before calling this function
void ContinuumRobotStateEstimator::queryAdditionalStates(ContinuumRobotStateEstimator::SystemState &state, std::vector<std::pair<unsigned int, double>> arclengths)
{
    state = m_state;



    //Convert state to body frames (our equations are written that way)
    convertStateMeanBodyInertial(state);

    //Sorts pairs in increasing order of their first value (or second value if first value is the same)
    std::sort(arclengths.begin(), arclengths.end());

    for(unsigned int i = 0; i < arclengths.size(); i++)
    {
        unsigned int robot_idx = arclengths[i].first;
        double arclength = arclengths[i].second;

        //Make sure the entries in the vector are correct
        assert((robot_idx < m_robot_topology.N) && "Index needs to be between 0 and N-1");
        assert((arclength >= 0) && (arclength <= m_robot_topology.L[robot_idx]) && "Arclength needs to be between 0 and total length L");

        //Create a new node
        SystemState::RobotState::Node node;
        node.arclength = arclength;

        //Get the estimation node before and after the interpolated node
        SystemState::RobotState::Node node_k;
        SystemState::RobotState::Node node_k1;
        int k_idx = 0;
        for(unsigned int k = 0; k < m_robot_topology.K[robot_idx]; k++)
        {
            node_k = state.robots[robot_idx].estimation_nodes[k];
            node_k1 = state.robots[robot_idx].estimation_nodes[k+1];
            k_idx = k;
            // Check if the requested arclength agrees with the arclength of the estimation node
            // Add a small epsilon disturbance to be sure we don't get numerical issues with comparing two doubles
            // If the requested arclength is L, the loop will end with node_k1 being the last node
            // Requested arclength can't be bigger than L, we took care about that in the assert above
            if(node.arclength + 1e-10 >= node_k.arclength && node.arclength - 1e-10 <= node_k1.arclength)
            {
                break;
            }

        }

        //Interpolate between the found nodes

        //Get pose of current and next node
        Eigen::Matrix4d T_k = node_k.pose;
        Eigen::Matrix4d T_k1 = node_k1.pose;

        //Get arclength of current and next node
        double s_k = node_k.arclength;
        double s_k1 = node_k1.arclength;

        //Get strain of current and next node
        Eigen::MatrixXd varpi_k = node_k.strain;
        Eigen::MatrixXd varpi_k1 = node_k1.strain;

        //Compute quantities used below
        Eigen::MatrixXd xi = tran_to_vec(T_k1*invert_transformation(T_k));
        Eigen::MatrixXd J_inv = vec_to_jac_inverse(xi);

        double delta_s = s_k1 - s_k;

        Eigen::MatrixXd Qc = m_hyperparameters.Qc;
        Eigen::MatrixXd Qc_inv = invert_diagonal(m_hyperparameters.Qc);

        //Covariance for prior terms
        Eigen::Matrix<double,12,12> Q_inv;
        Q_inv << 12/(delta_s*delta_s*delta_s)*Qc_inv, -6/(delta_s*delta_s)*Qc_inv,
                -6/(delta_s*delta_s)*Qc_inv, 4/(delta_s)*Qc_inv;


        Eigen::Matrix<double,12,12> phi;
        phi << Eigen::Matrix<double,6,6>::Identity(), delta_s*Eigen::Matrix<double,6,6>::Identity(),
                Eigen::Matrix<double,6,6>::Zero(), Eigen::Matrix<double,6,6>::Identity();

        //Interpolation coefficients
        double s_q = arclength;

        double delta_s_q = s_q - s_k;

        Eigen::MatrixXd Q;
        Q.resize(2*Qc.rows(),2*Qc.cols());

        Q << delta_s_q*delta_s_q*delta_s_q/3*Qc, delta_s_q*delta_s_q/2*Qc,
                delta_s_q*delta_s_q/2*Qc, delta_s_q*Qc;

        Eigen::Matrix<double,12,12> psi_temp;

        psi_temp << Eigen::Matrix<double,6,6>::Identity(), (s_k1 - s_q)*Eigen::Matrix<double,6,6>::Identity(),
                Eigen::Matrix<double,6,6>::Zero(), Eigen::Matrix<double,6,6>::Identity();

        Eigen::MatrixXd psi = Q*psi_temp.transpose()*Q_inv;


        Eigen::Matrix<double,12,12> lambda;

        lambda << Eigen::Matrix<double,6,6>::Identity(), delta_s_q*Eigen::Matrix<double,6,6>::Identity(),
                Eigen::Matrix<double,6,6>::Zero(), Eigen::Matrix<double,6,6>::Identity();
        lambda = lambda - psi*phi;

        //Interpolate the mean
        Eigen::Matrix<double,12,1> gamma_1;
        gamma_1 << Eigen::Matrix<double,6,1>::Zero(),
                varpi_k;

        Eigen::Matrix<double,12,1> gamma_2;
        gamma_2 << xi,
                J_inv*varpi_k1;

        Eigen::MatrixXd gamma = lambda*gamma_1 + psi*gamma_2;

        Eigen::Matrix4d T_q = vec_to_tran(gamma.topRows(6)) * T_k;
        node.pose = T_q; //pose

        Eigen::MatrixXd J_tau = vec_to_jac(gamma.topRows(6));

        Eigen::MatrixXd varpi_q = J_tau*gamma.bottomRows(6);
        node.strain = varpi_q; // strain

        //Interpolate the covariance
        Eigen::Matrix<double,12,12> Ga_1;
        Ga_1 << Eigen::Matrix<double,6,6>::Identity(), Eigen::Matrix<double,6,6>::Zero(),
                0.5*curlyhat(varpi_k),  Eigen::Matrix<double,6,6>::Identity();


        Eigen::Matrix<double,12,12> Ga_2;
        Ga_2 << J_inv, Eigen::Matrix<double,6,6>::Zero(),
                0.5*curlyhat(varpi_k1)*J_inv, J_inv;


        Eigen::Matrix<double,12,12> Xi_1;
        Xi_1 << Eigen::Matrix<double,6,6>::Identity(), Eigen::Matrix<double,6,6>::Zero(),
                Eigen::Matrix<double,6,6>::Zero(), Eigen::Matrix<double,6,6>::Zero();

        Eigen::Matrix<double,12,12> Xi_2;
        Xi_2 << tran_adjoint(vec_to_tran(xi)), Eigen::Matrix<double,6,6>::Zero(),
                Eigen::Matrix<double,6,6>::Zero(), Eigen::Matrix<double,6,6>::Zero();


        int k_offset = 0;
        for(unsigned int n = 0; n < robot_idx; n++)
        {
            k_offset += m_robot_topology.K[n];
        }

        int k_start = 12*k_idx + 12*k_offset;

        Eigen::MatrixXd P11 = m_covariance.block(k_start,k_start,12,12);
        Eigen::MatrixXd P12 = m_covariance.block(k_start,k_start+12,12,12);
        Eigen::MatrixXd P22 = m_covariance.block(k_start+12,k_start+12,12,12);

        Eigen::MatrixXd P11k = Ga_1*(P11 - Xi_1*P11*Xi_1.transpose())*Ga_1.transpose();
        Eigen::MatrixXd P12k = Ga_1*(P12 - Xi_1*P11*Xi_2.transpose())*Ga_2.transpose();
        Eigen::MatrixXd P22k = Ga_2*(P22 - Xi_2*P11*Xi_2.transpose())*Ga_2.transpose();

        Eigen::Matrix<double,24,24> Pk;
        Pk << P11k, P12k,
                P12k.transpose(), P22k;


        Eigen::MatrixXd Q_tau;
        Eigen::Matrix<double,12,12> Q_tau_temp;
        Q_tau_temp << Eigen::Matrix<double,6,6>::Identity(), (s_k1-s_q)*Eigen::Matrix<double,6,6>::Identity(),
                Eigen::Matrix<double,6,6>::Zero(), Eigen::Matrix<double,6,6>::Identity();

        Q_tau = Q - psi*Q_tau_temp*Q;

        Eigen::Matrix<double,12,24> lambda_psi;
        lambda_psi << lambda, psi;

        Eigen::MatrixXd Pq_local = lambda_psi*Pk*lambda_psi.transpose() + Q_tau;

        Eigen::Matrix<double,12,12> Ga_inv;
        Ga_inv << J_tau, Eigen::Matrix<double,6,6>::Zero(),
                -0.5*J_tau*curlyhat(varpi_q), J_tau;

        Eigen::Matrix<double,12,12> Xi;
        Xi << tran_adjoint(vec_to_tran(gamma.topRows(6))), Eigen::Matrix<double,6,6>::Zero(),
                Eigen::Matrix<double,6,6>::Zero(),Eigen::Matrix<double,6,6>::Zero();

        Eigen::MatrixXd covariance_q = Ga_inv*Pq_local*Ga_inv.transpose() + Xi*m_covariance.block(k_start,k_start,12,12)*Xi.transpose();

        //Save out the relevant variables

        Eigen::Matrix3d rot = node.pose.block(0,0,3,3);

        //We are already transforming in the world frame here (but do it later for the state itself)
        Eigen::Matrix3d position_covariance = rot.transpose()*covariance_q.block(0,0,3,3)*rot;
        Eigen::Matrix3d orientation_covariance = rot.transpose()*covariance_q.block(3,3,3,3)*rot;
        Eigen::Matrix3d nu_covariance = covariance_q.block(6,6,3,3);
        Eigen::Matrix3d omega_covariance = covariance_q.block(9,9,3,3);

        // We are saving the position covariance so that we can plot ellipsoids from it
        node.position_covariance = position_covariance;
        node.orientation_covariance = orientation_covariance;
        node.nu_covariance = nu_covariance;
        node.omega_covariance = omega_covariance;

        //Save standard deviations
        Eigen::Matrix<double,6,1> pose_std;
        pose_std << std::sqrt(position_covariance(0,0)),
                std::sqrt(position_covariance(1,1)),
                std::sqrt(position_covariance(2,2)),
                std::sqrt(orientation_covariance(0,0)),
                std::sqrt(orientation_covariance(1,1)),
                std::sqrt(orientation_covariance(2,2));

        Eigen::Matrix<double,6,1> strain_std;
        strain_std << std::sqrt(nu_covariance(0,0)),
                std::sqrt(nu_covariance(1,1)),
                std::sqrt(nu_covariance(2,2)),
                std::sqrt(omega_covariance(0,0)),
                std::sqrt(omega_covariance(1,1)),
                std::sqrt(omega_covariance(2,2));

        node.pose_std = pose_std;
        node.strain_std = strain_std;

        //Push back the node to the corresponding robot (will be ordered by arclength, since we sorted in the beginning)
        state.robots[robot_idx].queried_nodes.push_back(node);

    }


    //Convert state to inertial frames (more intuitive to work with)
    convertStateMeanBodyInertial(state);

    //Set saved system state to the state
    m_state = state;

}


