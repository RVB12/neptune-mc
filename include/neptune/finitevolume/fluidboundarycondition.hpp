#ifndef NEPTUNE_FLUID_BOUNDARY_CONDITION_HPP
#define NEPTUNE_FLUID_BOUNDARY_CONDITION_HPP

#include "neptune/common/utils.hpp"

namespace neptune {

/**
 * A fluid boundary condition
 */
class FluidBoundaryCondition {
   public:
    FluidBoundaryCondition() = default;
    ~FluidBoundaryCondition() = default;

    virtual char get_type() const = 0;
};

/**
 * Purely absorbing fluid boundary condition
 */
class AbsorbingFluidBoundaryCondition : public FluidBoundaryCondition {
   public:
    AbsorbingFluidBoundaryCondition() : FluidBoundaryCondition() {}

    char get_type() const override { return 'a'; }
};

/**
 * No flux fluid boundary condition.
 * This is the same as a reflective boundary condition
 */
class NoFluxFluidBoundaryCondition : public FluidBoundaryCondition {
   public:
    NoFluxFluidBoundaryCondition() : FluidBoundaryCondition() {}
    char get_type() const override { return 'n'; }
};

/**
 *
 */
class RobinFluidBoundaryCondition : public FluidBoundaryCondition {
   public:
    RobinFluidBoundaryCondition(real_t alpha, real_t beta, real_t gamma)
        : FluidBoundaryCondition(), alpha(alpha), beta(beta), gamma(gamma)
    {}

    real_t alpha;
    real_t beta;
    real_t gamma;

    char get_type() const override { return 'r'; }
};

}  // namespace neptune

#endif