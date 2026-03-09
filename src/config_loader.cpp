// This file contains the actual code behind every method declared in
// config_loader.h. Other files never include this directly; they only see
// the header. The compiler links the two together automatically.
//
// Sections in this file mirror the top-level sections of the YAML config:
//   [topology]         → parseTopology()
//   [hyperparameters]  → parseHyperparameters()
//   [options]          → parseOptions()
//   [measurements]     → parseMeasurements()      (optional)
//   [visualization]    → parseVisualization()     (optional)

#include "config_loader.h"
#include "utilities.h"         // load_csv<>() helper

#include <yaml-cpp/yaml.h>     // YAML parser
#include <Eigen/Geometry>      // Eigen matrix types
#include <filesystem>          // std::filesystem::path  (C++17)
#include <stdexcept>           // std::runtime_error

// The parse methods receive YAML::Node as a void* pointer so that yaml-cpp
// does not appear in the public header (keeping compile dependencies minimal).
// This macro casts it back to a const YAML::Node& for convenient use.
#define AS_NODE(ptr) (*static_cast<const YAML::Node*>(ptr))

// These all do the same job for different types:
//   1. Look up a named field in a YAML node.
//   2. Throw a descriptive error if it does not exist.
//   3. Convert the YAML sequence to a std::vector.
//   4. Optionally verify the element count matches `expected`.
static std::vector<double> readDoubleArray(const YAML::Node& node, const std::string& field, size_t expected = 0)
{
    if (!node[field])
        throw std::runtime_error("Missing required field: " + field);
    auto vec = node[field].as<std::vector<double>>();
    if (expected > 0 && vec.size() != expected)
        throw std::runtime_error(field + ": expected " + std::to_string(expected) +
                                 " elements, got " + std::to_string(vec.size()));
    return vec;
}

static std::vector<unsigned int> readUintArray(const YAML::Node& node, const std::string& field, size_t expected = 0)
{
    if (!node[field])
        throw std::runtime_error("Missing required field: " + field);
    auto vec = node[field].as<std::vector<unsigned int>>();
    if (expected > 0 && vec.size() != expected)
        throw std::runtime_error(field + ": expected " + std::to_string(expected) +
                                 " elements, got " + std::to_string(vec.size()));
    return vec;
}

static std::vector<bool> readBoolArray(const YAML::Node& node, const std::string& field, size_t expected = 0)
{
    if (!node[field])
        throw std::runtime_error("Missing required field: " + field);
    auto vec = node[field].as<std::vector<bool>>();
    if (expected > 0 && vec.size() != expected)
        throw std::runtime_error(field + ": expected " + std::to_string(expected) +
                                 " elements, got " + std::to_string(vec.size()));
    return vec;
}

static std::vector<int> readIntArray(const YAML::Node& node, const std::string& field, size_t expected = 0)
{
    if (!node[field])
        throw std::runtime_error("Missing required field: " + field);
    auto vec = node[field].as<std::vector<int>>();
    if (expected > 0 && vec.size() != expected)
        throw std::runtime_error(field + ": expected " + std::to_string(expected) +
                                 " elements, got " + std::to_string(vec.size()));
    return vec;
}

// Constructor
ConfigLoader::ConfigLoader(const std::string& config_path)
{
    // Remember the directory containing the config file so that relative paths
    // written in the YAML (e.g. "data/fbg.csv") can be resolved to absolute
    // paths later in parseMeasurements / parseTransformMatrix.
    std::filesystem::path p(config_path);
    m_config_dir = p.parent_path().string();
    if (m_config_dir.empty()) m_config_dir = ".";

    // Load and parse the entire YAML file into memory.
    // This is the only point where disk I/O happens.
    YAML::Node root;
    try {
        root = YAML::LoadFile(config_path);
    } catch (const YAML::BadFile&) {
        throw std::runtime_error("Cannot open config file: " + config_path);
    } catch (const YAML::ParserException& e) {
        throw std::runtime_error("YAML parse error in " + config_path + ": " + e.what());
    }

    // Validate that the three mandatory top-level sections exist.
    if (!root["topology"])
        throw std::runtime_error("Config file missing 'topology' section");
    if (!root["hyperparameters"])
        throw std::runtime_error("Config file missing 'hyperparameters' section");
    if (!root["options"])
        throw std::runtime_error("Config file missing 'options' section");

    // Extract each section as a sub-node and parse it.
    // We pass a pointer to avoid copying large YAML trees.
    YAML::Node topology_node = root["topology"];
    YAML::Node hp_node       = root["hyperparameters"];
    YAML::Node opt_node      = root["options"];

    parseTopology(&topology_node);
    parseHyperparameters(&hp_node);
    parseOptions(&opt_node);

    // Optional sections — parsed only if present in the file.
    if (root["measurements"]) {
        YAML::Node meas_node = root["measurements"];
        parseMeasurements(&meas_node);
    }

    if (root["visualization"]) {
        YAML::Node vis_node = root["visualization"];
        parseVisualization(&vis_node);
    }
}

// Accessors
// Simple getters — return copies of the already-parsed member variables.
// Calling these is cheap; no YAML parsing happens here.

ContinuumRobotStateEstimator::RobotTopology ConfigLoader::getTopology() const { return m_topology; }
ContinuumRobotStateEstimator::Hyperparameters ConfigLoader::getHyperparameters() const { return m_hyperparameters; }
ContinuumRobotStateEstimator::Options ConfigLoader::getOptions() const { return m_options; }
std::vector<ContinuumRobotStateEstimator::SensorMeasurement> ConfigLoader::getMeasurements() const { return m_measurements; }
ConfigLoader::VisualizationSettings ConfigLoader::getVisualizationSettings() const { return m_vis_settings; }

// Path resolution
// Converts a path written in the YAML file (possibly relative, e.g. "data/fbg.csv")
// into an absolute path by prepending the config file's own directory.
// Absolute paths (starting with '/') are returned unchanged.
std::string ConfigLoader::resolvePath(const std::string& relative_path) const
{
    std::filesystem::path p(relative_path);
    if (p.is_absolute()) return relative_path;
    return (std::filesystem::path(m_config_dir) / p).string();
}

// Transform matrix parsing

// Reads a 4x4 homogeneous transformation matrix from a YAML node.
// A homogeneous matrix encodes both rotation (top-left 3x3) and translation
// (right column) in a single 4x4 matrix. It is used here to describe the
// base frame of each robot (Ti0) and inter-robot coupling transforms.
//
// Supported YAML formats:
//
//   type: identity          → returns the 4x4 identity matrix
//
//   type: translation       → pure translation, no rotation
//   xyz: [x, y, z]
//
//   type: matrix            → explicit 4x4 matrix
//   data:
//     - [r00, r01, r02, tx]
//     - [r10, r11, r12, ty]
//     - [r20, r21, r22, tz]
//     - [0,   0,   0,   1 ]
//
//   type: csv_file          → load 4x4 matrix from a CSV file
//   path: relative/path.csv
//   orthonormalize: true    → (optional) fix floating-point drift in rotation
Eigen::Matrix4d ConfigLoader::parseTransformMatrix(const void* ptr)
{
    const YAML::Node& node = AS_NODE(ptr);

    if (!node["type"])
        throw std::runtime_error("Transform matrix entry missing 'type' field");

    std::string type = node["type"].as<std::string>();

    if (type == "identity") {
        return Eigen::Matrix4d::Identity();
    }
    else if (type == "translation") {
        if (!node["xyz"])
            throw std::runtime_error("Transform 'translation' missing required 'xyz' field");
        auto xyz = node["xyz"].as<std::vector<double>>();
        if (xyz.size() != 3)
            throw std::runtime_error("Transform 'translation' xyz must have 3 elements");
        Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
        T(0,3) = xyz[0];
        T(1,3) = xyz[1];
        T(2,3) = xyz[2];
        return T;
    }
    else if (type == "matrix") {
        if (!node["data"])
            throw std::runtime_error("Transform 'matrix' requires 'data' field");
        Eigen::Matrix4d T;
        auto rows = node["data"];
        if (rows.size() != 4)
            throw std::runtime_error("Transform matrix data must have 4 rows");
        for (int i = 0; i < 4; i++) {
            auto row = rows[i].as<std::vector<double>>();
            if (row.size() != 4)
                throw std::runtime_error("Transform matrix row must have 4 elements");
            for (int j = 0; j < 4; j++)
                T(i,j) = row[j];
        }
        return T;
    }
    else if (type == "csv_file") {
        if (!node["path"])
            throw std::runtime_error("Transform 'csv_file' requires 'path' field");
        std::string path = resolvePath(node["path"].as<std::string>());
        Eigen::Matrix4d T = load_csv<Eigen::Matrix4d>(path);

        // Optional: re-orthonormalize the rotation part of the matrix.
        // Useful when the matrix was measured/exported with small numerical errors
        // that break the orthonormality constraint (R^T R = I).
        if (node["orthonormalize"] && node["orthonormalize"].as<bool>()) {
            Eigen::Vector3d Rx = T.block(0,0,3,1).normalized();
            Eigen::Vector3d Ry = T.block(0,1,3,1).normalized();
            Eigen::Vector3d Rz = Rx.cross(Ry).normalized();
            Ry = Rz.cross(Rx).normalized();
            T.block(0,0,3,3) << Rx, Ry, Rz;
        }
        return T;
    }
    else {
        throw std::runtime_error("Unknown transform type: " + type);
    }
}

// Topology

// Parses the [topology] section which defines the physical structure of the
// multi-robot system:
//
//   N                    : number of robots
//   K                    : discretization nodes per robot (array of N values)
//   M                    : segments per robot (array of N values)
//   L                    : arc length of each robot [m]
//   lock_first/last_pose : fix boundary pose (e.g. clamped base)
//   lock_first/last_strain: fix boundary strain (e.g. zero strain at tip)
//   fbg_core_distance    : distance between FBG fiber cores [m]
//   fbg_theta_offset     : angular offset of the FBG core array [rad]
//   common_end_effector  : true if all robots share a common tip frame
//   Ti0                  : list of N base-frame transforms (one per robot)
//   robot_coupling       : optional kinematic constraints between robots
void ConfigLoader::parseTopology(const void* ptr)
{
    const YAML::Node& t = AS_NODE(ptr);

    if (!t["N"])
        throw std::runtime_error("topology: missing 'N' (number of robots)");
    m_topology.N = t["N"].as<unsigned int>();
    unsigned int N = m_topology.N;

    // Read per-robot arrays — each must have exactly N elements.
    m_topology.K = readUintArray(t, "K", N);
    m_topology.M = readUintArray(t, "M", N);
    m_topology.L = readDoubleArray(t, "L", N);

    m_topology.lock_first_pose   = readBoolArray(t, "lock_first_pose",   N);
    m_topology.lock_last_pose    = readBoolArray(t, "lock_last_pose",    N);
    m_topology.lock_first_strain = readBoolArray(t, "lock_first_strain", N);
    m_topology.lock_last_strain  = readBoolArray(t, "lock_last_strain",  N);

    m_topology.fbg_core_distance = readDoubleArray(t, "fbg_core_distance", N);
    m_topology.fbg_theta_offset  = readDoubleArray(t, "fbg_theta_offset",  N);

    // Default to false if the field is absent (single-robot configs omit it).
    m_topology.common_end_effector = t["common_end_effector"] ? t["common_end_effector"].as<bool>() : false;

    // Parse base frames Ti0 — one 4x4 transform per robot describing where
    // each robot's base is mounted in the world frame.
    m_topology.Ti0.clear();
    if (!t["Ti0"])
        throw std::runtime_error("topology: missing 'Ti0' (base frames)");
    auto ti0_seq = t["Ti0"];
    if (ti0_seq.size() != N)
        throw std::runtime_error("topology.Ti0: expected " + std::to_string(N) +
                                 " entries, got " + std::to_string(ti0_seq.size()));
    for (size_t i = 0; i < ti0_seq.size(); i++) {
        YAML::Node entry = ti0_seq[i];
        m_topology.Ti0.push_back(parseTransformMatrix(&entry));
    }

    // Parse optional robot_coupling entries.
    // Each entry constrains two robots to share a common point (e.g. shared
    // end-effector or a physical joint). The mask selects which of the 6 DOFs
    // (3 translation + 3 rotation) the constraint applies to.
    m_topology.robot_coupling.clear();
    if (t["robot_coupling"] && t["robot_coupling"].IsSequence()) {
        for (size_t i = 0; i < t["robot_coupling"].size(); i++) {
            YAML::Node c = t["robot_coupling"][i];
            ContinuumRobotStateEstimator::RobotTopology::Coupling coupling;

            coupling.idxA = c["idxA"].as<unsigned int>();   // Index of robot A
            coupling.idxB = c["idxB"].as<unsigned int>();   // Index of robot B
            coupling.coupling_node_robot_A = c["coupling_node_robot_A"].as<unsigned int>();
            coupling.coupling_node_robot_B = c["coupling_node_robot_B"].as<unsigned int>();

            YAML::Node tbac = c["T_bA_c"];
            coupling.T_bA_c = parseTransformMatrix(&tbac);  // Transform: body-A to coupling point

            YAML::Node tbbc = c["T_bB_c"];
            coupling.T_bB_c = parseTransformMatrix(&tbbc);  // Transform: body-B to coupling point

            // 6-element mask: 1 = constrain this DOF, 0 = ignore it
            auto mask_vec = c["mask"].as<std::vector<int>>();
            if (mask_vec.size() != 6)
                throw std::runtime_error("Coupling mask must have 6 elements");
            for (int j = 0; j < 6; j++)
                coupling.mask(j) = mask_vec[j];

            m_topology.robot_coupling.push_back(coupling);
        }
    }
}

// Hyperparameters

// Parses the [hyperparameters] section which defines the noise covariance
// matrices used by the Gaussian-process state estimator.
//
// These matrices control how much the estimator trusts each sensor vs the
// physical model:
//   - Small noise  → estimator trusts that sensor more
//   - Large noise  → estimator relies more on the physical model
//
// Two input formats are supported so that the YAML can be written concisely:
//
// FORMAT 1  —  noise_std (scalar standard deviations + scale factors)
//   noise_std:
//     R_p:   0.001   # position noise std [m]
//     R_o:   0.01    # orientation noise std [rad]
//     R_fbg: 0.005   # FBG strain noise std
//   R_pose_scale: 1.0
//   ...
//   The code builds diagonal covariance matrices as σ² × scale × I.
//
// FORMAT 2  —  direct diagonal (explicit covariance diagonal)
//   R_pose_diagonal:       [1e-6, 1e-6, 1e-6, 1e-4, 1e-4, 1e-4]
//   R_fbg_strain_diagonal: [1e-5, 1e-5, 1e-5, 1e-5]
//   ...
void ConfigLoader::parseHyperparameters(const void* ptr)
{
    const YAML::Node& hp = AS_NODE(ptr);

    // Zero-initialize all matrices so any un-set matrix stays at zero.
    m_hyperparameters.R_pose.setZero();
    m_hyperparameters.R_strain.setZero();
    m_hyperparameters.R_fbg_strain.setZero();
    m_hyperparameters.R_coupling.setZero();
    m_hyperparameters.Qc.setZero();

    if (hp["noise_std"]) {
        // ---- FORMAT 1: build covariances from scalar standard deviations ----
        auto ns = hp["noise_std"];
        double R_p   = ns["R_p"].as<double>();    // position std [m]
        double R_o   = ns["R_o"].as<double>();    // orientation std [rad]
        double R_fbg = ns["R_fbg"].as<double>();  // FBG strain std

        // R_pose: 6x6 diagonal with [σ_p², σ_p², σ_p², σ_o², σ_o², σ_o²]
        Eigen::Matrix<double,6,1> R_pose_diag;
        R_pose_diag << R_p*R_p, R_p*R_p, R_p*R_p, R_o*R_o, R_o*R_o, R_o*R_o;
        double R_pose_scale = hp["R_pose_scale"].as<double>();
        m_hyperparameters.R_pose = R_pose_scale * R_pose_diag.asDiagonal();

        // R_strain: optional (example 5 uses only FBG, not direct strain)
        if (ns["R_v"] && ns["R_u"] && hp["R_strain_scale"]) {
            double R_v = ns["R_v"].as<double>();  // linear  strain std
            double R_u = ns["R_u"].as<double>();  // angular strain std
            Eigen::Matrix<double,6,1> R_strain_diag;
            R_strain_diag << R_v*R_v, R_v*R_v, R_v*R_v, R_u*R_u, R_u*R_u, R_u*R_u;
            double R_strain_scale = hp["R_strain_scale"].as<double>();
            m_hyperparameters.R_strain = R_strain_scale * R_strain_diag.asDiagonal();
        }

        // R_fbg_strain: 4x4 diagonal (FBG has 4 strain components)
        Eigen::Matrix<double,4,1> R_fbg_diag;
        R_fbg_diag << R_fbg*R_fbg, R_fbg*R_fbg, R_fbg*R_fbg, R_fbg*R_fbg;
        double R_fbg_scale = hp["R_fbg_strain_scale"].as<double>();
        m_hyperparameters.R_fbg_strain = R_fbg_scale * R_fbg_diag.asDiagonal();

        // R_coupling: 6x6 diagonal for inter-robot coupling constraint noise
        auto R_coupling_vec = readDoubleArray(hp, "R_coupling_diagonal", 6);
        double R_coupling_scale = hp["R_coupling_scale"].as<double>();
        Eigen::Matrix<double,6,1> R_coupling_diag;
        for (int i = 0; i < 6; i++) R_coupling_diag(i) = R_coupling_vec[i];
        m_hyperparameters.R_coupling = R_coupling_scale * R_coupling_diag.asDiagonal();

        // Qc: 6x6 process noise matrix (controls stiffness of the GP prior)
        // Larger Qc → more flexible robot model
        auto Qc_vec = readDoubleArray(hp, "Qc_diagonal", 6);
        double Qc_scale = hp["Qc_scale"].as<double>();
        Eigen::Matrix<double,6,1> Qc_diag;
        for (int i = 0; i < 6; i++) Qc_diag(i) = Qc_vec[i];
        m_hyperparameters.Qc = Qc_scale * Qc_diag.asDiagonal();
    }
    else {
        // ---- FORMAT 2: direct diagonal vectors ----
        // Helper lambda: converts a 6-element vector to a 6x6 diagonal matrix.
        auto toMat6 = [](const std::vector<double>& v) {
            Eigen::Matrix<double,6,1> d;
            for (int i = 0; i < 6; i++) d(i) = v[i];
            return Eigen::Matrix<double,6,6>(d.asDiagonal());
        };

        m_hyperparameters.R_pose = toMat6(readDoubleArray(hp, "R_pose_diagonal", 6));

        if (hp["R_strain_diagonal"])
            m_hyperparameters.R_strain = toMat6(readDoubleArray(hp, "R_strain_diagonal", 6));

        auto fbg_vec = readDoubleArray(hp, "R_fbg_strain_diagonal", 4);
        Eigen::Matrix<double,4,1> fbg_d;
        for (int i = 0; i < 4; i++) fbg_d(i) = fbg_vec[i];
        m_hyperparameters.R_fbg_strain = Eigen::Matrix<double,4,4>(fbg_d.asDiagonal());

        m_hyperparameters.R_coupling = toMat6(readDoubleArray(hp, "R_coupling_diagonal", 6));
        m_hyperparameters.Qc         = toMat6(readDoubleArray(hp, "Qc_diagonal", 6));
    }
}

// Options

// Parses the [options] section which controls solver behaviour:
//
//   initial_guess        : how to initialise the state before optimising
//                          "Straight"  → assume all robots are straight
//                          "Last"      → reuse the previous solution (warm start)
//                          "Custom"    → caller provides the initial state
//   solver               : "Newton" or "NewtonLineSearch"
//                          NewtonLineSearch is more robust but slightly slower
//   max_iterations       : stop after this many Newton steps (default 200)
//   convergence_threshold: stop when the update norm falls below this (default 0.5)
//   kirchhoff_rods       : if true, shear and extension strains are ignored
//                          (Kirchhoff rod model — bending and torsion only)
void ConfigLoader::parseOptions(const void* ptr)
{
    const YAML::Node& opt = AS_NODE(ptr);

    // Initial guess type (default: Straight)
    std::string ig = opt["initial_guess"] ? opt["initial_guess"].as<std::string>() : "Straight";
    if (ig == "Straight")
        m_options.init_guess_type = ContinuumRobotStateEstimator::Options::Straight;
    else if (ig == "Last")
        m_options.init_guess_type = ContinuumRobotStateEstimator::Options::Last;
    else if (ig == "Custom")
        m_options.init_guess_type = ContinuumRobotStateEstimator::Options::Custom;
    else
        throw std::runtime_error("options.initial_guess: unknown value '" + ig + "' (expected Straight, Last, or Custom)");

    // Solver type (default: NewtonLineSearch)
    std::string sol = opt["solver"] ? opt["solver"].as<std::string>() : "NewtonLineSearch";
    if (sol == "Newton")
        m_options.solver = ContinuumRobotStateEstimator::Options::Newton;
    else if (sol == "NewtonLineSearch")
        m_options.solver = ContinuumRobotStateEstimator::Options::NewtonLineSearch;
    else
        throw std::runtime_error("options.solver: unknown value '" + sol + "' (expected Newton or NewtonLineSearch)");

    m_options.max_optimization_iterations = opt["max_iterations"]        ? opt["max_iterations"].as<unsigned int>()  : 200;
    m_options.convergence_threshold        = opt["convergence_threshold"] ? opt["convergence_threshold"].as<double>() : 0.5;
    m_options.kirchhoff_rods               = opt["kirchhoff_rods"]        ? opt["kirchhoff_rods"].as<bool>()          : true;
}

// Measurements

// Parses the optional [measurements] section which injects sensor observations
// directly from the YAML file (useful for offline / replay scenarios).
//
// Each entry in the list specifies:
//   type       : "Pose", "Strain", or "FBGStrain"
//   idx_robot  : which robot this measurement belongs to (0-indexed)
//   idx_node   : the discretization node where the sensor is located
//                (alternatively: idx_node_range: [start, end] for a batch)
//   mask       : 6-element binary vector selecting which DOFs are observed
//   value      : the measured value
//   source     : optional "csv_file" to load many measurements from a CSV
//
// Mask example:  [1,1,1,0,0,0]  → observe only position (x,y,z), not rotation
void ConfigLoader::parseMeasurements(const void* ptr)
{
    const YAML::Node& meas_list = AS_NODE(ptr);
    if (!meas_list.IsSequence()) return;

    for (size_t i = 0; i < meas_list.size(); i++) {
        YAML::Node mj = meas_list[i];

        std::string type_str = mj["type"].as<std::string>();

        // ---- Special case: FBG measurements loaded from a CSV file ----
        // Each column in the CSV becomes one FBGStrain measurement at a
        // successive node index (column 0 → node 0, column 1 → node 1, ...).
        if (mj["source"] && mj["source"].as<std::string>() == "csv_file") {
            std::string path = resolvePath(mj["path"].as<std::string>());
            unsigned int idx_robot = mj["idx_robot"].as<unsigned int>();

            Eigen::Matrix<int,6,1> mask;
            auto mask_vec = mj["mask"].as<std::vector<int>>();
            for (int j = 0; j < 6; j++) mask(j) = mask_vec[j];

            Eigen::MatrixXd data = load_csv<Eigen::MatrixXd>(path);
            for (int k = 0; k < data.cols(); k++) {
                ContinuumRobotStateEstimator::SensorMeasurement m;
                m.type      = ContinuumRobotStateEstimator::SensorMeasurement::FBGStrain;
                m.idx_robot = idx_robot;
                m.idx_node  = k;
                m.mask      = mask;
                m.value     = data.col(k);
                m_measurements.push_back(m);
            }
            continue;  // move to the next measurement entry
        }

        // ---- Batch (node range) vs single node ----
        // idx_node_range: [2, 5]  → creates measurements at nodes 2, 3, 4, 5
        // idx_node: 3             → creates a single measurement at node 3
        bool is_range = mj["idx_node_range"].IsDefined();

        unsigned int start_node = 0, end_node = 0;
        if (is_range) {
            auto range = mj["idx_node_range"].as<std::vector<unsigned int>>();
            if (range.size() != 2)
                throw std::runtime_error("measurements: idx_node_range must have 2 elements [start, end]");
            start_node = range[0];
            end_node   = range[1];
        } else {
            start_node = mj["idx_node"].as<unsigned int>();
            end_node   = start_node;
        }

        unsigned int idx_robot = mj["idx_robot"].as<unsigned int>();
        Eigen::Matrix<int,6,1> mask;
        auto mask_vec = mj["mask"].as<std::vector<int>>();
        for (int j = 0; j < 6; j++) mask(j) = mask_vec[j];

        // Create one SensorMeasurement struct per node in the range.
        for (unsigned int n = start_node; n <= end_node; n++) {
            ContinuumRobotStateEstimator::SensorMeasurement m;
            m.idx_robot = idx_robot;
            m.idx_node  = n;
            m.mask      = mask;

            if (type_str == "Strain") {
                // Direct strain measurement: 6-element vector [v1,v2,v3,u1,u2,u3]
                m.type = ContinuumRobotStateEstimator::SensorMeasurement::Strain;
                auto val = mj["value"].as<std::vector<double>>();
                if (val.size() != 6)
                    throw std::runtime_error("Strain measurement value must have 6 elements");
                Eigen::Matrix<double,6,1> strain;
                for (int j = 0; j < 6; j++) strain(j) = val[j];
                m.value = strain;
            }
            else if (type_str == "Pose") {
                // Pose measurement: stored as a 4x4 homogeneous transform
                m.type = ContinuumRobotStateEstimator::SensorMeasurement::Pose;
                YAML::Node val_node = mj["value"];
                m.value = parseTransformMatrix(&val_node);
            }
            else if (type_str == "FBGStrain") {
                // FBG strain: 4-element vector from a fiber Bragg grating sensor
                m.type = ContinuumRobotStateEstimator::SensorMeasurement::FBGStrain;
                auto val = mj["value"].as<std::vector<double>>();
                if (val.size() != 4)
                    throw std::runtime_error("FBGStrain measurement value must have 4 elements");
                Eigen::Matrix<double,4,1> fbg;
                for (int j = 0; j < 4; j++) fbg(j) = val[j];
                m.value = fbg;
            }
            else {
                throw std::runtime_error("Unknown measurement type: " + type_str);
            }

            m_measurements.push_back(m);
        }
    }
}

// Visualization

// Parses the optional [visualization] section.
// All fields are optional; missing fields keep their default values
// (set in the VisualizationSettings struct definition in the header).
void ConfigLoader::parseVisualization(const void* ptr)
{
    const YAML::Node& vis = AS_NODE(ptr);

    if (vis["window_width"])      m_vis_settings.window_width      = vis["window_width"].as<int>();
    if (vis["window_height"])     m_vis_settings.window_height     = vis["window_height"].as<int>();
    if (vis["render_frames"])     m_vis_settings.render_frames     = vis["render_frames"].as<bool>();
    if (vis["render_covariance"]) m_vis_settings.render_covariance = vis["render_covariance"].as<bool>();
    if (vis["covariance_n_std"])  m_vis_settings.covariance_n_std  = vis["covariance_n_std"].as<int>();
    if (vis["verbose"])           m_vis_settings.verbose           = vis["verbose"].as<bool>();
}
