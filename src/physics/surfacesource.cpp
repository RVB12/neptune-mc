#include "neptune/physics/surfacesource.hpp"

#include <iostream>

#include "neptune/physics/physicalconstants.hpp"

namespace neptune {

Particle SurfaceSource::generate_particle(index_t i_particle, const Mesh& mesh, const Background& background)
{
    if (i_particle >= num_particles()) {
        throw std::runtime_error(
            "SurfaceSource::generate_particle: "
            "Tried to generate more particles than "
            "allowed by num_particles");
    }
    Particle p(i_particle + particle_id_offset_);

    p.geometry.cell_id = cell_idx;
    p.geometry.boundary_id = surface_id;

    p.geometry.position =
        mesh.generate_random_surface_position(cell_idx, surface_id, p.rng);

    // Calculate normal of the surface and two tangent vectors
    Vec3 normal =
        mesh.get_surface_normal(cell_idx, surface_id, p.geometry.position);
    Vec3 tangent1;
    if (std::abs(normal.z) < 0.9) {
        tangent1 = Vec3(0, 0, 1).cross(normal);
    }
    else {
        tangent1 = Vec3(1, 0, 0).cross(normal);
    }
    tangent1 /= tangent1.norm();
    Vec3 tangent2 = normal.cross(tangent1);
    // The cross of two normalized vectors is already normalized
    // so this next line is not needed
    // tangent2 /= tangent2.norm();

    // Now generate a direction from the cosine distribution
    std::uniform_real_distribution<real_t> dist_u(0, 1);

    real_t phi = 2 * M_PI * dist_u(p.rng);
    real_t cos_theta = std::sqrt(dist_u(p.rng));
    real_t sin_theta = std::sqrt(1 - cos_theta * cos_theta);

    Vec3 direction = -cos_theta * normal +
                     sin_theta * std::cos(phi) * tangent1 +
                     sin_theta * std::sin(phi) * tangent2;
    direction /= direction.norm();
    
    p.geometry.velocity = direction * v_norm;

    std::exponential_distribution<real_t> dist_t(1.0);
    p.u = dist_t(p.rng);

    return p;
}

}  // namespace neptune