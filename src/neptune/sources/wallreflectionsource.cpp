#include "neptune/sources/wallreflectionsource.hpp"

#include <iostream>

#include "neptune/common/utils.hpp"
#include "neptune/physics/physicalconstants.hpp"

namespace neptune {

WallReflectionSource::WallReflectionSource(index_t num_particles_,
                                           index_t cell_idx_,
                                           index_t surface_id_,
                                           real_t sheath_potential_,
                                           real_t franck_condon_temp_,
                                           real_t Rt_,
                                           real_t Rr_,
                                           Species background_species,
                                           Species particle_species,
                                           const Mesh& mesh,
                                           const Background& background)
    : Source(num_particles_),
      cell_idx(cell_idx_),
      surface_id(surface_id_),
      sheath_potential(sheath_potential_),
      franck_condon_temp(franck_condon_temp_),
      Rt(Rt_),
      Rr(Rr_),
      background_species(background_species),
      particle_species(particle_species)
{
    // Do some calculations to determine the intensity
    Vec3 normal = mesh.get_surface_normal(cell_idx, surface_id, Vec3(1, 0, 0));
    Vec3 boundary_drift_3d = background.velocity(background_species, cell_idx);
    real_t boundary_drift = dot(boundary_drift_3d, normal);

    real_t ion_temp = background.temperature(background_species, cell_idx);
    real_t boundary_bgk_stddev =
        std::sqrt(ion_temp * constants::eV_to_Joules / constants::m_p *
                  SPECIES_MASS[static_cast<int>(background_species)]);

    real_t ion_density = background.density(background_species, cell_idx);

    real_t flux = boundary_bgk_stddev / std::sqrt(2 * M_PI) *
                      std::exp(-0.5 * (boundary_drift / boundary_bgk_stddev) *
                               (boundary_drift / boundary_bgk_stddev)) +
                  boundary_drift * 0.5 *
                      (1 + std::erf(boundary_drift /
                                    (boundary_bgk_stddev * std::sqrt(2))));
    flux *= ion_density;

    real_t surface_area = mesh.get_surface_area(cell_idx, surface_id);

    intensity_ = flux * Rt * surface_area;
}

Particle WallReflectionSource::generate_particle(uint64_t particle_id,
                                                 rng_key_t key,
                                                 const Mesh& mesh,
                                                 const Background& background)
{
    Particle p(particle_id, key, particle_species);

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

    Vec3 boundary_drift_3d = background.velocity(background_species, cell_idx);
    real_t boundary_drift = dot(boundary_drift_3d, normal);

    real_t ion_temp = background.temperature(background_species, cell_idx);
    real_t boundary_bgk_stddev =
        std::sqrt(ion_temp * constants::eV_to_Joules /
                  (constants::m_p *
                   SPECIES_MASS[static_cast<int>(background_species)]));

    p.intensity = intensity() / num_particles();

    // Draw a random number to determine if the particle is backscattered or not
    std::uniform_real_distribution<real_t> dist_u(0, 1);
    real_t r = dist_u(p.rng);
    if (r < Rr) {
        // The particle is backscattered, draw a velocity from the maxwellian of
        // the ions but only consider particles that are moving towards the wall
        // (this is a sample from a truncated normal distribution)
        real_t x = truncated_normal_right_tail(
            -boundary_drift / boundary_bgk_stddev, p.rng);
        real_t vel_to_wall = x * boundary_bgk_stddev + boundary_drift;

        real_t drift_tangent_1 = dot(boundary_drift_3d, tangent1);
        real_t drift_tangent_2 = dot(boundary_drift_3d, tangent2);

        std::normal_distribution<real_t> dist_v1(drift_tangent_1,
                                                 boundary_bgk_stddev);
        std::normal_distribution<real_t> dist_v2(drift_tangent_2,
                                                 boundary_bgk_stddev);

        real_t vel_tangent_1 = dist_v1(p.rng);
        real_t vel_tangent_2 = dist_v2(p.rng);

        real_t alpha = 0.5;
        vel_to_wall =
            std::sqrt(alpha * (vel_to_wall * vel_to_wall +
                               2 * sheath_potential * constants::eV_to_Joules /
                                   constants::m_p));

        p.geometry.velocity = vel_to_wall * normal + vel_tangent_1 * tangent1 +
                              vel_tangent_2 * tangent2;

        // Now reflect the velocity specularly
        p.geometry.velocity -= 2 * (p.geometry.velocity.dot(normal)) * normal;
    }
    else {
        // The particle undergoes Franck-Condon dissociation, draw a new
        // velocity from a maxwellian distribution with the franck-condon
        // temperature
        real_t bgk_stddev =
            std::sqrt(franck_condon_temp * constants::eV_to_Joules /
                      (constants::m_p *
                       SPECIES_MASS[static_cast<int>(particle_species)]));
        std::normal_distribution<real_t> dist_vx(0, bgk_stddev);
        std::normal_distribution<real_t> dist_vy(0, bgk_stddev);
        std::normal_distribution<real_t> dist_vz(0, bgk_stddev);
        p.geometry.velocity =
            Vec3(dist_vx(p.rng), dist_vy(p.rng), dist_vz(p.rng));

        // Make sure the velocity is pointing into the domain
        if (p.geometry.velocity.dot(normal) > 0) {
            p.geometry.velocity -=
                2 * (p.geometry.velocity.dot(normal)) * normal;
        }
    }

    std::exponential_distribution<real_t> dist_t(1.0);
    p.u = dist_t(p.rng);

    return p;
}

}  // namespace neptune