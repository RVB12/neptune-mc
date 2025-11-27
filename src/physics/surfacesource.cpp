#include "neptune/physics/surfacesource.hpp"
#include "neptune/physics/physicalconstants.hpp"

#include <iostream>

namespace neptune {

Particle SurfaceSource::generate_particle(index_t i_particle, const Mesh& mesh)
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

    p.geometry.position = mesh.generate_random_surface_position(cell_idx, surface_id, p.rng);

    std::normal_distribution<real_t> dist_vx(vel_mean.x, vel_stddev);
    std::normal_distribution<real_t> dist_vy(vel_mean.y, vel_stddev);
    std::normal_distribution<real_t> dist_vz(vel_mean.z, vel_stddev);

    p.geometry.velocity = Vec3(dist_vx(p.rng), dist_vy(p.rng), dist_vz(p.rng));

    // Make sure that the velocity is pointing away from the surface
    Vec3 normal = mesh.get_surface_normal(cell_idx, surface_id, p.geometry.position);
    if (p.geometry.velocity.dot(normal) > 0) {
        // Reflect the velocity vector
        p.geometry.velocity = p.geometry.velocity - normal * (2 * p.geometry.velocity.dot(normal));
    }

    return p;
}

}  // namespace neptune