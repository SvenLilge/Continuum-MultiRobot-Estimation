// This file declares the ConfigLoader class. It tells other files WHAT the
// class can do, without revealing HOW it does it (that lives in config_loader.cpp).

// Usage:
// ConfigLoader cfg("path/to/config.yaml");
// auto topology = cfg.getTopology();
// auto hp       = cfg.getHyperparameters();

#ifndef CONFIG_LOADER_H   // "Include guard": prevents this file from being
#define CONFIG_LOADER_H   // processed more than once by the compiler.

#include "continuum_robot_state_estimator.h"
#include <string>
#include <vector>

// Reads a YAML configuration file once (in the constructor) and stores all
// parsed parameters internally. Callers then use the get*() methods to
// retrieve the data. The YAML file is NOT re-read on every get() call.
class ConfigLoader
{
public:

    // Optional rendering parameters read from the [visualization] section of
    // the YAML file. All fields have sensible defaults so the section can be
    // omitted entirely.
    struct VisualizationSettings {
        int  window_width      = 1280;  // Render window width  (pixels)
        int  window_height     = 720;   // Render window height (pixels)
        bool render_frames     = true;  // Draw coordinate frames on each node
        bool render_covariance = true;  // Draw uncertainty ellipsoids
        int  covariance_n_std  = 3;     // Ellipsoid radius in # of std deviations
        bool verbose           = true;  // Print extra info to the terminal
    };

    // Opens and fully parses the YAML file at config_path. Throws
    // std::runtime_error if the file cannot be opened, is malformed, or is
    // missing required fields.
    //
    // After construction all get*() methods are safe to call at any time.
    explicit ConfigLoader(const std::string& config_path);

    // Physical robot structure (number of robots, lengths, base frames, ...)
    ContinuumRobotStateEstimator::RobotTopology getTopology()const;

    // Noise covariance matrices (R_pose, R_strain, Qc, ...)
    ContinuumRobotStateEstimator::Hyperparameters getHyperparameters() const;

    // Solver settings (Newton / NewtonLineSearch, max iterations, ...)
    ContinuumRobotStateEstimator::Options getOptions() const;

    // Optional list of sensor measurements defined in the YAML file
    std::vector<ContinuumRobotStateEstimator::SensorMeasurement> getMeasurements() const;

    // Optional list of control inputs defined in the YAML file
    std::vector<ContinuumRobotStateEstimator::ControlInput> getControlInputs() const;

    // Optional visualization / rendering settings
    VisualizationSettings getVisualizationSettings() const;

private:

    // Internal storage  (m_ prefix = "member variable")
    // These hold the parsed data. They are filled once by the constructor and
    // are read-only after that. External code cannot access them directly.

    std::string m_config_dir;   // Directory of the config file (used to resolve
                                // relative paths to CSV / mesh files)

    ContinuumRobotStateEstimator::RobotTopology   m_topology;
    ContinuumRobotStateEstimator::Hyperparameters m_hyperparameters;
    ContinuumRobotStateEstimator::Options         m_options;
    std::vector<ContinuumRobotStateEstimator::SensorMeasurement> m_measurements;
    std::vector<ContinuumRobotStateEstimator::ControlInput> m_control_inputs;
    VisualizationSettings m_vis_settings;

    // Internal parse methods
    // Each method handles one top-level section of the YAML file.
    // They receive a raw pointer to a YAML::Node cast to void* so that
    // yaml-cpp does not need to be exposed in this header (keeping compile
    // times low and the public API clean).
    void parseTopology        (const void* node);
    void parseHyperparameters (const void* node);
    void parseOptions         (const void* node);
    void parseMeasurements    (const void* node);
    void parseControlInputs   (const void* node);
    void parseVisualization   (const void* node);

    // Parses a 4x4 homogeneous transform from a YAML node.
    // Supported formats: "identity", "translation", "matrix", "csv_file".
    Eigen::Matrix4d parseTransformMatrix(const void* node);

    // Converts a path that may be relative (to the config file's directory)
    // into an absolute path. Absolute paths are returned unchanged.
    std::string resolvePath(const std::string& relative_path) const;
};

#endif  // CONFIG_LOADER_H
