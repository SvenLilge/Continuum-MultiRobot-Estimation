#include "continuum_rod_priors.h"

#include <sstream>
#include <stdexcept>

namespace {

// Optional field: either empty (not populated) OR exactly (N x 3).
void require_Nx3_or_empty(const Eigen::MatrixXd& m, Eigen::Index N, const char* name)
{
    if (m.size() == 0) return;
    if (m.rows() != N || m.cols() != 3) {
        std::ostringstream oss;
        oss << "ContinuumRodPriors::validate: field \"" << name
            << "\" has shape (" << m.rows() << "x" << m.cols()
            << "), expected (" << N << "x3).";
        throw std::invalid_argument(oss.str());
    }
}

} // anonymous namespace

void ContinuumRodPriors::validate() const
{
    const Eigen::Index N = s.size();
    if (N < 1) {
        throw std::invalid_argument("ContinuumRodPriors::validate: s is empty");
    }

    // Required fields.
    if (v.rows() != N || v.cols() != 3) {
        std::ostringstream oss;
        oss << "ContinuumRodPriors::validate: v has shape (" << v.rows()
            << "x" << v.cols() << "), expected (" << N << "x3).";
        throw std::invalid_argument(oss.str());
    }
    if (u.rows() != N || u.cols() != 3) {
        std::ostringstream oss;
        oss << "ContinuumRodPriors::validate: u has shape (" << u.rows()
            << "x" << u.cols() << "), expected (" << N << "x3).";
        throw std::invalid_argument(oss.str());
    }

    // Optional fields: either empty OR (N x 3).
    require_Nx3_or_empty(v_dot,      N, "v_dot");
    require_Nx3_or_empty(u_dot,      N, "u_dot");
    require_Nx3_or_empty(n_internal, N, "n_internal");
    require_Nx3_or_empty(m_internal, N, "m_internal");
    require_Nx3_or_empty(f_dist,     N, "f_dist");
    require_Nx3_or_empty(l_dist,     N, "l_dist");

    // s must be monotonically non-decreasing.
    for (Eigen::Index i = 1; i < N; ++i) {
        if (s(i) < s(i - 1)) {
            std::ostringstream oss;
            oss << "ContinuumRodPriors::validate: s is not monotonic at index "
                << i << " (s[" << (i - 1) << "]=" << s(i - 1)
                << " > s[" << i << "]=" << s(i) << ").";
            throw std::invalid_argument(oss.str());
        }
    }

    // Discrete-load arrays must all have the same length.
    const std::size_t M = s_discrete.size();
    if (F_discrete.size() != M || L_discrete.size() != M) {
        std::ostringstream oss;
        oss << "ContinuumRodPriors::validate: discrete arrays have mismatched "
            << "sizes (s=" << M << ", F=" << F_discrete.size()
            << ", L=" << L_discrete.size() << ").";
        throw std::invalid_argument(oss.str());
    }
}
