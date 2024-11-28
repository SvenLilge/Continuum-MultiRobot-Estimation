#ifndef UTILITIES_H
#define UTILITIES_H

#include <Eigen/Core>
#include <fstream>

//Generate the nth Bernoulli number (hard coded for performance)
double bernoulli_number(unsigned int n);

// Builds the 6x6 curly hat matrix from the 6x1 inputs
//
// input:
//   vec: 6x1 vector xi
//
// output:
//   vec_curly_hat: the 6x6 curly hat matrix
//
// From: Timothy D Barfoot and Paul T Furgale,
//       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
//		 DOI: 10.1109/TRO.2014.2298059
Eigen::MatrixXd curlyhat(Eigen::MatrixXd vec);

// Builds the 3x3 skew symmetric matrix from the 3x1 input or 4x4 from 6x1 input
//
// input:
//   vec: 3x1 or 6x1 vector xi
//
// output:
//   vec_hat: 3x3 or 4x4 matrix
//
// From: Timothy D Barfoot and Paul T Furgale,
//       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
//		 DOI: 10.1109/TRO.2014.2298059
Eigen::MatrixXd hat(Eigen::MatrixXd vec);

// Computes the matrix log of the rotation matrix C
//
// input:
//   C: a 3x3 rotation matrix
//
// output:
//   phi: a 3x1 vector (axis * angle) computed from C
//
// From: Timothy D Barfoot and Paul T Furgale,
//       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
//		 DOI: 10.1109/TRO.2014.2298059
Eigen::MatrixXd rot_to_vec(Eigen::MatrixXd C);

// Throws an error if the rotation matrix is not orthonormal
//
// input:
//   C: a 3x3 matrix
//
// From: Timothy D Barfoot and Paul T Furgale,
//       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
//		 DOI: 10.1109/TRO.2014.2298059
void validate_rotation_matrix(Eigen::Matrix3d C);

// Computes the matrix log of the transformation matrix T
//
// input:
//   T: a 4x4 transformation matrix
//
// output:
//   p: a 6x1 vector in tangent coordinates computed from T
//
// From: Timothy D Barfoot and Paul T Furgale,
//       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
//		 DOI: 10.1109/TRO.2014.2298059
Eigen::MatrixXd tran_to_vec(Eigen::MatrixXd T);

// Computes the 6x6 adjoint mapping of the 4x4 transformation matrix
//
// input:
//   T: a 4x4 transformation matrix
//
// output:
//   Ad_T: a 6x6 adjoint transformation
//
// From: Timothy D Barfoot and Paul T Furgale,
//       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
//		 DOI: 10.1109/TRO.2014.2298059
Eigen::MatrixXd tran_adjoint(Eigen::MatrixXd T);

// Throws an error if the transformation matrix is not in SE(3)
//
// input:
//   T: a 4x4 matrix
//
// From: Timothy D Barfoot and Paul T Furgale,
//       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
//		 DOI: 10.1109/TRO.2014.2298059
void validate_transformation_matrix(Eigen::Matrix4d T);

// Construction of the 3x3 J matrix or 6x6 J matrix
//
// input:
//   vec: 3x1 or 6x1 vector xi
//
// output:
//   J: 3x3 or 6x6 jacobian matrix
//
// From: Timothy D Barfoot and Paul T Furgale,
//       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
//		 DOI: 10.1109/TRO.2014.2298059
Eigen::MatrixXd vec_to_jac(Eigen::MatrixXd vec);

// Construction of the 3x3 J^-1 matrix or 6x6 J^-1 matrix
//
// input:
//   vec: 3x1 or 6x1 vector xi
//
// output:
//   inv_J: 3x3 or 6x6 inverse jacobian matrix
//
// From: Timothy D Barfoot and Paul T Furgale,
//       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
//		 DOI: 10.1109/TRO.2014.2298059
Eigen::MatrixXd vec_to_jac_inverse(Eigen::MatrixXd vec);

// Construction of the 3x3 J^-1 matrix or 6x6 J^-1 matrix using series representation
//
// input:
//   vec: 3x1 or 6x1 vector xi
//   N:   number of terms to include in the series
//
// output:
//   inv_J: 3x3 or 6x6 inverse jacobian matrix
//
// From: Timothy D Barfoot and Paul T Furgale,
//       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
//		 DOI: 10.1109/TRO.2014.2298059
Eigen::MatrixXd vec_to_jac_inverse_series(Eigen::MatrixXd vec, unsigned int N);

// Construction of the 3x3 J matrix or 6x6 J matrix using series representation
//
// input:
//   vec: 3x1 or 6x1 vector xi
//   N:   number of terms to include in the series
//
// output:
//   inv_J: 3x3 or 6x6 jacobian matrix
//
// From: Timothy D Barfoot and Paul T Furgale,
//       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
//		 DOI: 10.1109/TRO.2014.2298059
Eigen::MatrixXd vec_to_jac_series(Eigen::MatrixXd vec, unsigned int N);

// Construction of the 3x3 Q matrix
//
// input:
//   vec: 6x1 vector xi
//
// output:
//   Q: 3x3 Q matrix
//
// From: Timothy D Barfoot and Paul T Furgale,
//       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
//		 DOI: 10.1109/TRO.2014.2298059
Eigen::MatrixXd vec_to_Q(Eigen::MatrixXd vec);

// Build a rotation matrix using the exponential map
//
// input:
//   phi: 3x1 vector
//
// output:
//   C: 3x3 rotation matrix
//
// From: Timothy D Barfoot and Paul T Furgale,
//       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
//		 DOI: 10.1109/TRO.2014.2298059
Eigen::MatrixXd vec_to_rot(Eigen::MatrixXd phi);

// Build a rotation matrix using the exponential map using series representation
//
// input:
//   phi: 3x1 vector
//   N:   number of terms to include in the series
//
// output:
//   C: 3x3 rotation matrix
//
// From: Timothy D Barfoot and Paul T Furgale,
//       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
//		 DOI: 10.1109/TRO.2014.2298059
Eigen::MatrixXd vec_to_rot_series(Eigen::MatrixXd phi, unsigned int N);

// Build a transformation matrix using the exponential map
//
// input:
//   vec: 6x1 vector
//
// output:
//   T: 4x4 transformation matrix
//
// From: Timothy D Barfoot and Paul T Furgale,
//       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
//		 DOI: 10.1109/TRO.2014.2298059
Eigen::MatrixXd vec_to_tran(Eigen::MatrixXd vec);

// Build a transformation matrix using the exponential map with series representation
//
// input:
//   vec: 6x1 vector
//   N:   number of terms to include in the series
//
// output:
//   T: 4x4 transformation matrix
//
// From: Timothy D Barfoot and Paul T Furgale,
//       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
//		 DOI: 10.1109/TRO.2014.2298059
Eigen::MatrixXd vec_to_tran_series(Eigen::MatrixXd vec, unsigned int N);

void validate_diagonal_matrix(Eigen::MatrixXd matrix);

Eigen::Matrix4d invert_transformation(Eigen::Matrix4d T);

Eigen::MatrixXd invert_diagonal(Eigen::MatrixXd matrix);

// Load CSV files to Eigen Matrix
// From: https://stackoverflow.com/questions/34247057/how-to-read-csv-file-and-assign-to-eigen-matrix
template<typename M> M load_csv(const std::string & path) {
    std::ifstream indata;
    indata.open(path);
    std::string line;
    std::vector<double> values;
    uint rows = 0;
    while (std::getline(indata, line)) {
        std::stringstream lineStream(line);
        std::string cell;
        while (std::getline(lineStream, cell, ',')) {
            values.push_back(std::stod(cell));
        }
        ++rows;
    }
    return Eigen::Map<const Eigen::Matrix<typename M::Scalar, M::RowsAtCompileTime, M::ColsAtCompileTime, Eigen::RowMajor>>(values.data(), rows, values.size()/rows);
}

#endif // UTILITIES_H
