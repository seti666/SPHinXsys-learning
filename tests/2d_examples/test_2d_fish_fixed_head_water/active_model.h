#ifndef ACTIVE_MODEL_H
#define ACTIVE_MODEL_H

#include "complex_solid.h"
#include "elastic_dynamics.h"

namespace SPH
{
class ActiveModelSolid : public SaintVenantKirchhoffSolid
{
    Matd *active_strain_;

  public:
    explicit ActiveModelSolid(Real rho0, Real youngs_modulus, Real poisson_ratio);
    virtual ~ActiveModelSolid(){};

    virtual void initializeLocalParameters(BaseParticles *base_particles) override;
    virtual Matd StressPK1(Matd &deformation, size_t particle_index_i) override;
};
} // namespace SPH
#endif // ACTIVE_MODEL_H
