#ifndef NEPTUNE_FLUID_BOUNDARY_CONDITION_HPP
#define NEPTUNE_FLUID_BOUNDARY_CONDITION_HPP

#include "neptune/common/utils.hpp"

namespace neptune {

/**
 * A fluid boundary condition
 * This is just a robin boundary condition:
 *  alpha * f + beta * df/dn = gamma
 */
class FluidBoundaryCondition {
   public:
    FluidBoundaryCondition(real_t alpha, real_t beta, real_t gamma)
        : alpha(alpha), beta(beta), gamma(gamma)
    {}

    real_t alpha = 0.0;  // The absorption coefficient
    real_t beta = 0.0;   // The reflection coefficient
    real_t gamma = 0.0;
};

}  // namespace neptune

#endif