#include "cosserat_priors_adapter.h"
#include "cosseratrodmodel.h"

ContinuumRodPriors priorsFromCosseratModel(const CosseratRodModel& model,
                                           double L1, double L2)
{
    ContinuumRodPriors p;

    p.s          = model.getArclengthSamples();
    p.v          = model.getStrainV();
    p.u          = model.getStrainU();
    p.v_dot      = model.getStrainVDot();
    p.u_dot      = model.getStrainUDot();
    p.n_internal = model.getInternalForce();
    p.m_internal = model.getInternalMoment();
    p.f_dist     = model.getDistributedForce();
    p.l_dist     = model.getDistributedMoment();

    Eigen::Vector3d F_junction, L_junction, F_tip, L_tip;
    model.getDiscreteLoads(F_junction, L_junction, F_tip, L_tip);
    p.s_discrete = { L1, L1 + L2 };
    p.F_discrete = { F_junction, F_tip };
    p.L_discrete = { L_junction, L_tip };

    p.validate();
    return p;
}
