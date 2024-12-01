#include "utilities.h"

#include <assert.h>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/MatrixFunctions>
#include <cmath>
#include <iostream>

double bernoulli_number(unsigned int n)
{
    switch(n)
    {
    case 0:
        return 1;
    case 1:
        return -0.5;
    case 2:
        return 0.166666666;
    case 3:
        return 0;
    case 4:
        return -0.033333333;
    case 5:
        return 0;
    case 6:
        return +0.023809523;
    case 7:
        return 0;
    case 8:
        return -0.033333333;
    case 9:
        return 0;
    case 10:
        return +0.075757575;
    case 11:
        return 0;
    case 12:
        return -0.253113553;
    case 13:
        return 0;
    case 14:
        return 1.166666666;
    case 15:
        return 0;
    }

    std::cout << "Warning: Only the first 15 Bernoulli numbers are supported!" << std::endl;
    return 0;

}

Eigen::MatrixXd curlyhat(Eigen::MatrixXd vec)
{
    assert((vec.cols() == 1 && vec.rows() == 6) && "Input needs to be 6x1 column vector");

    Eigen::Vector3d phi = vec.block(3,0,3,1);
    Eigen::Vector3d v = vec.block(0,0,3,1);

    Eigen::Matrix3d phi_hat = hat(phi);
    Eigen::Matrix3d v_hat = hat(v);

    Eigen::Matrix<double,6,6> vec_curly_hat;

    vec_curly_hat << phi_hat, v_hat,
            Eigen::Matrix3d::Zero(), phi_hat;

    return vec_curly_hat;

}

Eigen::MatrixXd hat(Eigen::MatrixXd vec)
{
    assert((vec.cols() == 1 && (vec.rows() == 3 || vec.rows() == 6)) && "Input needs to be 6x1 or 3x1 column vector");

    if(vec.rows() == 3)
    {
        Eigen::Matrix3d vec_hat;

        vec_hat << 0, -vec(2), vec(1),
                vec(2), 0, -vec(0),
                -vec(1), vec(0), 0;

        return vec_hat;

    }
    else if (vec.rows() == 6)
    {

        Eigen::Matrix3d vec_hat;

        vec_hat << 0, -vec(5), vec(4),
                vec(5), 0, -vec(3),
                -vec(4), vec(3), 0;

        Eigen::Matrix4d vec_hat_full;

        vec_hat_full << vec_hat, vec.block(0,0,3,1),
                0, 0, 0, 0;

        return vec_hat_full;

    }

    // This line will never be reached
    return Eigen::Matrix4d::Zero();
}

Eigen::MatrixXd rot_to_vec(Eigen::MatrixXd C)
{
    validate_rotation_matrix(C);

    // From ASRL lgmath
    const double phi_ba = acos(std::clamp(0.5 * (C.trace() - 1.0), -1.0, 1.0));
      const double sinphi_ba = sin(phi_ba);

      if (fabs(sinphi_ba) > 1e-5) {
        // General case, angle is NOT near 0, pi, or 2*pi
        Eigen::Vector3d axis;
        axis << C(2, 1) - C(1, 2), C(0, 2) - C(2, 0),
            C(1, 0) - C(0, 1);
        return (0.5 * phi_ba / sinphi_ba) * axis;

      } else if (fabs(phi_ba) > 1e-9) {
        // Angle is near pi or 2*pi
        // ** Note with this method we do not know the sign of 'phi', however since
        // we know phi is
        //    close to pi or 2*pi, the sign is unimportant..

        // Find the eigenvalues and eigenvectors
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigenSolver(C);

        // Try each eigenvalue
        for (int i = 0; i < 3; i++) {
          // Check if eigen value is near +1.0
          if (fabs(eigenSolver.eigenvalues()[i] - 1.0) < 1e-3) {
            // Get corresponding angle-axis
            Eigen::Vector3d aaxis_ba = phi_ba * eigenSolver.eigenvectors().col(i);
            return aaxis_ba;
          }
        }

        // Runtime error
        throw std::runtime_error(
            "so3 logarithmic map failed to find an axis-angle, "
            "angle was near pi, or 2*pi, but no eigenvalues were near 1");

      } else {
        // Angle is near zero
        return Eigen::Vector3d::Zero();
      }
}

void validate_rotation_matrix(Eigen::Matrix3d C)
{
    assert((C.cols() == 3 && C.rows() == 3) && "Input needs to be 3x3 matrix");

    Eigen::Matrix3d CtC = C.transpose()*C;
    Eigen::Matrix3d E = CtC - Eigen::Matrix3d::Identity();
    E = E.array().abs();

    double error = E.maxCoeff();

    assert((error < 1e-6) && "C must be a valid rotation matrix");

}

Eigen::MatrixXd tran_to_vec(Eigen::MatrixXd T)
{
    validate_transformation_matrix(T);

    Eigen::Matrix3d C = T.block(0,0,3,3);
    Eigen::Matrix<double,3,1> r = T.block(0,3,3,1);


    Eigen::MatrixXd phi = rot_to_vec(C);
    Eigen::MatrixXd inv_J = vec_to_jac_inverse(phi);


    Eigen::MatrixXd rho = inv_J*r;


    Eigen::Matrix<double,6,1> vec;

    vec << rho,
            phi;

    return vec;
}

Eigen::MatrixXd tran_adjoint(Eigen::MatrixXd T)
{
    validate_transformation_matrix(T);

    Eigen::Matrix3d C = T.block(0,0,3,3);
    Eigen::Matrix<double,3,1> r = T.block(0,3,3,1);

    Eigen::Matrix<double,6,6> Ad_T;
    Ad_T << C, hat(r)*C,
            Eigen::Matrix3d::Zero(), C;

    return Ad_T;
}

void validate_transformation_matrix(Eigen::Matrix4d T)
{
    assert((T.cols() == 4 && T.rows() == 4) && "Input needs to be 4x4 matrix");

    Eigen::Matrix3d C = T.block(0,0,3,3);
    validate_rotation_matrix(C);

    Eigen::Matrix<double,1,4> bottom_T;
    bottom_T << 0, 0, 0, 1;

    Eigen::MatrixXd diff = bottom_T - T.bottomRows(1);
    diff = diff.array().abs();

    double error = diff.maxCoeff();

    assert((error < 1e-10) && "T must be a valid transformation matrix");
}

Eigen::MatrixXd vec_to_jac(Eigen::MatrixXd vec)
{
    assert((vec.cols() == 1 && (vec.rows() == 3 || vec.rows() == 6)) && "Input needs to be 6x1 or 3x1 column vector");

    double tolerance = 1e-12;


    Eigen::MatrixXd J;

    if(vec.rows() == 3)
    {
        Eigen::MatrixXd phi = vec;

        double phi_norm = phi.norm();

        if(phi_norm < tolerance)
        {
            //Fall back on series representation
            J = vec_to_jac_series(phi,10);
        }
        else
        {
            Eigen::MatrixXd axis = phi/phi_norm;
            double cph = (1-std::cos(phi_norm))/phi_norm;
            double sph = std::sin(phi_norm)/phi_norm;
            J = sph*Eigen::Matrix3d::Identity() + (1-sph)*axis*axis.transpose()+cph*hat(axis);
        }

    }
    else if (vec.rows() == 6)
    {
        Eigen::MatrixXd rho = vec.block(0,0,3,1);
        Eigen::MatrixXd phi = vec.block(3,0,3,1);

        J.resize(6,6);
        Eigen::MatrixXd J_small;
        Eigen::MatrixXd Q;

        double phi_norm = phi.norm();

        if(phi_norm < tolerance)
        {
            //Fall back on series representation
            Q = 0.5*hat(rho) + (hat(phi)*hat(rho) + hat(rho)*hat(phi))/6;
            J_small = vec_to_jac_series(phi,10);
        }
        else
        {
            J_small = vec_to_jac(phi);
            Q = vec_to_Q(vec);
        }

        J << J_small, Q,
                Eigen::Matrix3d::Zero(), J_small;



    }

    return J;
}

Eigen::MatrixXd vec_to_jac_inverse(Eigen::MatrixXd vec)
{

    assert((vec.cols() == 1 && (vec.rows() == 3 || vec.rows() == 6)) && "Input needs to be 6x1 or 3x1 column vector");

    double tolerance = 1e-12;

    Eigen::MatrixXd inv_J;

    if(vec.rows() == 3)
    {
        Eigen::MatrixXd phi = vec;

        double phi_norm = phi.norm();

        if(phi_norm < tolerance)
        {
            //Fall back on series representation
            inv_J = vec_to_jac_inverse_series(phi,10);
        }
        else
        {
            Eigen::MatrixXd axis = phi/phi_norm;
            double phi_norm_half = 0.5*phi_norm;
            double cot_phi_norm_half = std::cos(phi_norm_half)/std::sin(phi_norm_half);

            inv_J = phi_norm_half*cot_phi_norm_half*Eigen::Matrix3d::Identity() + (1-phi_norm_half*cot_phi_norm_half)*axis*axis.transpose() - phi_norm_half*hat(axis);
        }

    }
    else if (vec.rows() == 6)
    {
        Eigen::MatrixXd rho = vec.block(0,0,3,1);
        Eigen::MatrixXd phi = vec.block(3,0,3,1);

        inv_J.resize(6,6);

        double phi_norm = phi.norm();

        if(phi_norm < tolerance)
        {
            //Fall back on series representation
            inv_J = vec_to_jac_inverse_series(vec,10);
        }
        else
        {
            Eigen::MatrixXd inv_J_small = vec_to_jac_inverse(phi);
            Eigen::MatrixXd Q = vec_to_Q(vec);

            inv_J << inv_J_small, -inv_J_small*Q*inv_J_small,
                    Eigen::Matrix3d::Zero(), inv_J_small;
        }

    }

    return inv_J;

}


Eigen::MatrixXd vec_to_jac_inverse_series(Eigen::MatrixXd vec, unsigned int N)
{
    assert((vec.cols() == 1 && (vec.rows() == 3 || vec.rows() == 6)) && "Input needs to be 6x1 or 3x1 column vector");

    Eigen::MatrixXd inv_J;

    if(vec.rows() == 3)
    {
        inv_J = Eigen::Matrix3d::Identity();
        Eigen::Matrix3d pxn = Eigen::Matrix3d::Identity();
        Eigen::MatrixXd px = hat(vec);

        for(unsigned int n = 1; n <= N; n++)
        {
            pxn = pxn*px/n;
            inv_J = inv_J + bernoulli_number(n)*pxn;
        }


    }
    else if (vec.rows() == 6)
    {
        inv_J = Eigen::Matrix<double,6,6>::Identity();
        Eigen::Matrix<double,6,6> pxn = Eigen::Matrix<double,6,6>::Identity();
        Eigen::MatrixXd px = curlyhat(vec);

        for(unsigned int n = 1; n <= N; n++)
        {
            pxn = pxn*px/n;
            inv_J = inv_J + bernoulli_number(n)*pxn;
        }

    }

    return inv_J;
}

Eigen::MatrixXd vec_to_jac_series(Eigen::MatrixXd vec, unsigned int N)
{
    assert((vec.cols() == 1 && (vec.rows() == 3 || vec.rows() == 6)) && "Input needs to be 6x1 or 3x1 column vector");

    Eigen::MatrixXd J;

    if(vec.rows() == 3)
    {
        J = Eigen::Matrix3d::Identity();
        Eigen::Matrix3d pxn = Eigen::Matrix3d::Identity();
        Eigen::MatrixXd px = hat(vec);

        for(unsigned int n = 1; n <= N; n++)
        {
            pxn = pxn*px/(n+1);
            J = J + pxn;
        }

    }
    else if (vec.rows() == 6)
    {
        J = Eigen::Matrix<double,6,6>::Identity();
        Eigen::Matrix<double,6,6> pxn = Eigen::Matrix<double,6,6>::Identity();
        Eigen::MatrixXd px = curlyhat(vec);

        for(unsigned int n = 1; n <= N; n++)
        {
            pxn = pxn*px/(n+1);
            J = J + pxn;
        }

    }

    return J;
}

Eigen::MatrixXd vec_to_Q(Eigen::MatrixXd vec)
{
    assert((vec.cols() == 1 && vec.rows() == 6) && "Input needs to be 6x1 column vector");

    Eigen::MatrixXd rho = vec.block(0,0,3,1);
    Eigen::MatrixXd phi = vec.block(3,0,3,1);

    double ph = phi.norm();
    double ph2 = ph*ph;
    double ph3 = ph2*ph;
    double ph4 = ph3*ph;
    double ph5 = ph4*ph;

    double cph = std::cos(ph);
    double sph = std::sin(ph);

    Eigen::Matrix3d rx = hat(rho);
    Eigen::Matrix3d px = hat(phi);

    Eigen::MatrixXd t1 = 0.5*rx;
    Eigen::MatrixXd t2 = ((ph-sph)/ph3)*(px*rx + rx*px + px*rx*px);
    double m3 = (1 - 0.5*ph2 - cph)/ph4;
    Eigen::MatrixXd t3 = -m3*(px*px*rx + rx*px*px - 3*px*rx*px);
    double m4 = 0.5*(m3 - 3*(ph-sph-ph3/6)/ph5);
    Eigen::MatrixXd t4 = -m4*(px*rx*px*px + px*px*rx*px);

    Eigen::MatrixXd Q = t1 + t2 + t3 + t4;

    return Q;
}

Eigen::MatrixXd vec_to_rot(Eigen::MatrixXd phi)
{
    assert((phi.rows() == 3 && phi.cols() == 1) && "Input needs to be 3x1 column vector");

    double tolerance = 1e-12;

    double angle = phi.norm();

    Eigen::MatrixXd C;

    if(angle < tolerance)
    {
        // Fall back to series representation
        C = vec_to_rot_series(phi,10);
    }
    else
    {
        Eigen::MatrixXd axis = phi/angle;

        double cp = std::cos(angle);
        double sp = std::sin(angle);

        C = cp*Eigen::Matrix3d::Identity() + (1-cp)*axis*axis.transpose() + sp*hat(axis);

    }

    return C;

}

Eigen::MatrixXd vec_to_rot_series(Eigen::MatrixXd phi, unsigned int N)
{
    assert((phi.rows() == 3 && phi.cols() == 1) && "Input needs to be 3x1 column vector");
    Eigen::Matrix3d C = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d xM = Eigen::Matrix3d::Identity();
    Eigen::MatrixXd cmPhi = hat(phi);

    for(unsigned int n = 1; n <=N; n++)
    {
        xM = xM*(cmPhi/n);
        C = C + xM;
    }

    C = C *(C.transpose()*C).sqrt().inverse();

    validate_rotation_matrix(C);

    return C;

}

Eigen::MatrixXd vec_to_tran(Eigen::MatrixXd vec)
{
    assert((vec.rows() == 6 && vec.cols() == 1) && "Input needs to be 6x1 column vector");

    Eigen::MatrixXd rho = vec.block(0,0,3,1);
    Eigen::MatrixXd phi = vec.block(3,0,3,1);

    Eigen::MatrixXd C = vec_to_rot(phi);
    Eigen::MatrixXd J = vec_to_jac(phi);

    Eigen::Matrix4d T = Eigen::Matrix4d::Identity();

    T.block(0,0,3,3) = C;
    T.block(0,3,3,1) = J*rho;

    return T;

}

Eigen::MatrixXd vec_to_tran_series(Eigen::MatrixXd vec, unsigned int N)
{
    assert((vec.rows() == 6 && vec.cols() == 1) && "Input needs to be 6x1 column vector");

    Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
    Eigen::Matrix4d xM = Eigen::Matrix4d::Identity();
    Eigen::MatrixXd bpr = hat(vec);

    for(unsigned int n = 1; n <= N; n++)
    {
        xM = xM*(bpr/n);
        T = T + xM;
    }

    validate_transformation_matrix(T);

    return T;

}

void validate_diagonal_matrix(Eigen::MatrixXd matrix)
{
    assert((matrix.rows() == matrix.cols()) && "Input needs to be square matrix");

    bool is_diagonal = true;
    for(unsigned int i = 0; i < matrix.rows(); i++)
    {
        for(unsigned int j = 0; j < matrix.cols(); j++)
        {
            if(i != j && std::abs(matrix(i,j)) > 1e-10)
            {
                is_diagonal = false;
            }

        }
    }

    assert(is_diagonal && "Matrix can only have non-zero entries on diagonal");
}

Eigen::Matrix4d invert_transformation(Eigen::Matrix4d T)
{
    validate_transformation_matrix(T);

    Eigen::Matrix4d inv_T = Eigen::Matrix4d::Identity();

    inv_T.block(0,0,3,3) = T.block(0,0,3,3).transpose();
    inv_T.block(0,3,3,1) = -T.block(0,0,3,3).transpose()*T.block(0,3,3,1);

    return inv_T;
}

Eigen::MatrixXd invert_diagonal(Eigen::MatrixXd matrix)
{
    validate_diagonal_matrix(matrix);

    Eigen::MatrixXd matrix_inv = matrix;

    for(unsigned int i = 0; i < matrix_inv.rows(); i++)
    {
        matrix_inv(i,i) = 1/matrix_inv(i,i);
    }

    return matrix_inv;
}


