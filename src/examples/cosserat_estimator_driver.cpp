// Combined driver: CosseratRodModel FK -> ContinuumRodPriors DTO.
// Requires USE_LOCAL_TDCR=ON.

#include "cosseratrodmodel.h"
#include "continuum_robot_state_estimator.h"
#include "continuum_rod_priors.h"
#include "cosserat_priors_adapter.h"

#include <Eigen/Core>
#include <iomanip>
#include <iostream>

namespace {

// Fixed-width number formatting for aligned columns.
std::string fmt(double v, int width = 10, int prec = 4)
{
    std::ostringstream os;
    os << std::setw(width) << std::fixed << std::setprecision(prec) << v;
    return os.str();
}

std::string fmtSci(double v, int width = 12, int prec = 3)
{
    std::ostringstream os;
    os << std::setw(width) << std::scientific << std::setprecision(prec) << v;
    return os.str();
}

void printVec3(const char* label, const Eigen::Vector3d& v)
{
    std::cout << "    " << label << " = ["
              << fmt(v(0)) << ", " << fmt(v(1)) << ", " << fmt(v(2)) << " ]\n";
}

} // anonymous namespace

int main()
{
    // --- Configuration -------------------------------------------------------
    CosseratRodModel model;
    Eigen::MatrixXd diskFrames;
    Eigen::Matrix<double, 6, 1> q;
    q << 0.5, 0.2, 0.0, 0.3, 0.0, 0.1;           // tendon tensions [N]
    Eigen::Vector3d f_ext = Eigen::Vector3d::Zero(); // external tip force [N]
    Eigen::Vector3d l_ext = Eigen::Vector3d::Zero(); // external tip moment [Nm]
    const double L1 = 0.1, L2 = 0.1;                // segment lengths [m]

    // --- Forward kinematics --------------------------------------------------
    std::cout << "========================================\n";
    std::cout << "  Cosserat + Estimator Combined Driver\n";
    std::cout << "========================================\n\n";

    std::cout << "Input\n";
    std::cout << "  Tendon tensions [N]:  ["
              << q(0) << ", " << q(1) << ", " << q(2) << " | "
              << q(3) << ", " << q(4) << ", " << q(5) << "]\n";
    std::cout << "  Segment lengths [m]:  L1 = " << L1 << ",  L2 = " << L2 << "\n";
    std::cout << "  External tip force:   " << f_ext.transpose() << " N\n";
    std::cout << "  External tip moment:  " << l_ext.transpose() << " Nm\n\n";

    const bool ok = model.forwardKinematics(diskFrames, q, f_ext, l_ext);
    if (!ok || !model.hasAuxOutputs()) {
        std::cerr << "FK did not converge (residual = "
                  << model.getFinalResudial() << ")\n";
        return 1;
    }

    std::cout << "Forward kinematics converged  (residual = "
              << fmtSci(model.getFinalResudial()) << ")\n\n";

    // --- Extract priors ------------------------------------------------------
    const ContinuumRodPriors priors = priorsFromCosseratModel(model, L1, L2);

    // --- DTO shape -----------------------------------------------------------
    const int N = static_cast<int>(priors.s.size());
    std::cout << "ContinuumRodPriors\n";
    std::cout << "  Arclength samples:    N = " << N
              << "   s in [" << priors.s(0) << ", " << priors.s(N - 1) << "] m\n";
    std::cout << "  Strains v, u:         " << N << " x 3  (body frame)\n";
    std::cout << "  Distributed f, l:     " << N << " x 3  (world frame)\n";
    std::cout << "  Internal    n, m:     " << N << " x 3  (world frame)\n";
    std::cout << "  Discrete loads:       " << priors.s_discrete.size() << " entries\n\n";

    // --- Per-arclength ranges ------------------------------------------------
    std::cout << "Per-arclength ranges\n";
    std::cout << "  -----------------------------------------------\n";
    std::cout << "  Quantity              Max norm       Unit\n";
    std::cout << "  -----------------------------------------------\n";
    std::cout << "  ||v - e3||  (shear)  " << fmtSci((priors.v.rowwise() - Eigen::RowVector3d(0, 0, 1)).rowwise().norm().maxCoeff()) << "   [-]\n";
    std::cout << "  ||u||    (curvature) " << fmtSci(priors.u.rowwise().norm().maxCoeff()) << "   [rad/m]\n";
    std::cout << "  ||f_dist|| (tendon)  " << fmtSci(priors.f_dist.rowwise().norm().maxCoeff()) << "   [N/m]\n";
    std::cout << "  ||l_dist|| (tendon)  " << fmtSci(priors.l_dist.rowwise().norm().maxCoeff()) << "   [Nm/m]\n";
    std::cout << "  ||n||   (int force)  " << fmtSci(priors.n_internal.rowwise().norm().maxCoeff()) << "   [N]\n";
    std::cout << "  ||m||  (int moment)  " << fmtSci(priors.m_internal.rowwise().norm().maxCoeff()) << "   [Nm]\n";
    std::cout << "  -----------------------------------------------\n\n";

    // --- Discrete loads ------------------------------------------------------
    std::cout << "Discrete tendon-termination loads (world frame)\n\n";

    std::cout << "  Junction  s = " << priors.s_discrete[0] << " m  (seg-1 tendons terminate)\n";
    printVec3("F [N] ", priors.F_discrete[0]);
    printVec3("L [Nm]", priors.L_discrete[0]);
    std::cout << "    ||F|| = " << fmt(priors.F_discrete[0].norm(), 8, 4)
              << " N,   ||L|| = " << fmtSci(priors.L_discrete[0].norm()) << " Nm\n\n";

    std::cout << "  Tip       s = " << priors.s_discrete[1] << " m  (seg-2 tendons terminate)\n";
    printVec3("F [N] ", priors.F_discrete[1]);
    printVec3("L [Nm]", priors.L_discrete[1]);
    std::cout << "    ||F|| = " << fmt(priors.F_discrete[1].norm(), 8, 4)
              << " N,   ||L|| = " << fmtSci(priors.L_discrete[1].norm()) << " Nm\n\n";

    // --- Done ----------------------------------------------------------------
    std::cout << "========================================\n";
    std::cout << "  Pipeline OK\n";
    std::cout << "========================================\n";
    return 0;
}
