#ifndef NEPTUNE_BOUNDARYCONDITION_HPP
#define NEPTUNE_BOUNDARYCONDITION_HPP

#include <cassert>

#include "neptune/common/utils.hpp"
#include "neptune/geometry/vec3.hpp"
#include "neptune/physics/particle.hpp"

namespace neptune {

/**
 * Abstract class for a boundary condition
 */
class BoundaryCondition {
   public:
    BoundaryCondition() = default;
    ~BoundaryCondition() = default;

    /**
     * Apply the boundary condition to a particle
     * 
     * @param p      The particle to apply the boundary condition to
     *               It is assumed that the particle is on the boundary
     * @param normal The normal vector of the boundary
     */
    virtual void apply(Particle& p, Vec3 normal) = 0;
};

/**
 * Reflective boundary condition
 */
class SpecularBC : public BoundaryCondition {
   public:
    SpecularBC() = default;
    ~SpecularBC() = default;

    void apply(Particle& p, Vec3 normal) override
    {
        // Reflect the velocity of the particle
        p.geometry.velocity -= 2 * (p.geometry.velocity.dot(normal)) * normal;
    }
};

/**
 * Absorbing boundary condition
 */
class AbsorbingBC : public BoundaryCondition {
   public:
    AbsorbingBC() = default;
    ~AbsorbingBC() = default;

    void apply(Particle& p, Vec3 normal) override
    {
        p.weight = 0;  // Set the weight to zero to absorb the particle
    }
};

/**
 * Partially reflecting boundary condition
 */
class PartiallySpecularBC : public BoundaryCondition {
   public:
    PartiallySpecularBC(real_t absorption_coefficient)
        : absorption_coefficient(absorption_coefficient)
    {}
    ~PartiallySpecularBC() = default;

    void apply(Particle& p, Vec3 normal) override
    {
        // Reflect the velocity of the particle
        p.geometry.velocity -= 2 * (p.geometry.velocity.dot(normal)) * normal;
        // Absorb a fraction of the particle's weight
        p.weight *= (1 - absorption_coefficient);
    }

   private:
    real_t absorption_coefficient;
};

}  // namespace neptune

#endif