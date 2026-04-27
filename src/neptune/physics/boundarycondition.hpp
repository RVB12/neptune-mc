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

    real_t absorption_coefficient;
};

/**
 * Mixed reflection. Probability alpha of having diffuse reflection,
 * probability (1-alpha) of having specular reflection.
 *
 * Also absorbs a fraction of the particle's weight given by
 * absorption_coefficient.
 */
class MixedDiffuseReflectiveBC : public BoundaryCondition {
   public:
    MixedDiffuseReflectiveBC(real_t alpha, real_t absorption_coefficient)
        : alpha(alpha), absorption_coefficient(absorption_coefficient)
    {}
    ~MixedDiffuseReflectiveBC() = default;

    void apply(Particle& p, Vec3 normal) override
    {
        std::uniform_real_distribution<real_t> dist_u(0, 1);

        real_t u = dist_u(p.rng);

        if (u < alpha) {
            // Diffuse reflection
            Vec3 tangent1;
            if (std::abs(normal.z) < 0.9) {
                tangent1 = Vec3(0, 0, 1).cross(normal);
            }
            else {
                tangent1 = Vec3(1, 0, 0).cross(normal);
            }
            tangent1 /= tangent1.norm();
            Vec3 tangent2 = normal.cross(tangent1);

            real_t phi = 2 * M_PI * dist_u(p.rng);
            real_t cos_theta = std::sqrt(dist_u(p.rng));
            real_t sin_theta = std::sqrt(1 - cos_theta * cos_theta);

            Vec3 direction = -cos_theta * normal +
                             sin_theta * std::cos(phi) * tangent1 +
                             sin_theta * std::sin(phi) * tangent2;
            direction /= direction.norm();

            p.geometry.velocity = direction * p.geometry.velocity.norm();
        }
        else {
            // Specular reflection
            p.geometry.velocity -=
                2 * (p.geometry.velocity.dot(normal)) * normal;
        }

        // Absorb a fraction of the particle's weight
        p.weight *= (1 - absorption_coefficient);
    }

   private:
    real_t alpha;
    real_t absorption_coefficient;
};

}  // namespace neptune

#endif