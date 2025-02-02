#include <omp.h>

#include "eiron2/geometry/voxelmesh.hpp"
#include "eiron2/physics/linesource.hpp"
#include "eiron2/simulator/simulator_kdmc.hpp"

#include "eiron2/physics/physicalconstants.hpp"
namespace eiron2 {

void simulate_kdmc(Mesh& mesh,
                   PlasmaBackground& background,
                   vector<std::unique_ptr<Source>>& sources,
                   vector<real_t>& density,
                   vector<BoundaryCondition>& boundary_conditions,
                   SimulationParams params,
                   ExecutionStatistics& stats,
                   real_t dt)
{
    std::uniform_real_distribution<real_t> r_unif(0, 1);
    std::normal_distribution<real_t> r_normal(0, 1);

    // initialize the density array to zero
    for (index_t i = 0; i < mesh.n_cells(); i++) {
        density[i] = 0;
    }

    // History based simulation (because its easy)
    for (int i_source = 0; i_source < sources.size(); i_source++) {
#pragma omp parallel for
        for (int i_particle = 0;
             i_particle < sources[i_source]->num_particles(); i_particle++) {
            real_t intensity_multiplier = sources[i_source]->intensity() /
                                          sources[i_source]->num_particles();
            index_t num_crossings = 0;
            index_t num_collisions = 0;

            // Generate particle
            Particle particle =
                sources[i_source]->generate_particle(i_particle);
            particle.geometry.cell_id =
                mesh.get_cell_index(particle.geometry.position);

            if (particle.geometry.cell_id == -1) {
                throw std::runtime_error("Particle generated outside of mesh");
            }

            for (index_t i_step = 0; i_step < params.max_steps; i_step++) {
                real_t weight_old = particle.weight;
                // real_t t_old = particle.t;
                index_t cell_id_old = particle.geometry.cell_id;

                Vec2 drift_term = Vec2(background.plasma_vel_x[cell_id_old], background.plasma_vel_y[cell_id_old]) * dt;
                real_t D0 = background.plasma_temp[cell_id_old] / background.plasma_dens[cell_id_old] / constants::K_cx;
                Vec2 diffusion_term = Vec2(r_normal(particle.rng), r_normal(particle.rng)) * std::sqrt(2.0 * dt * D0);

                Vec2 new_position = particle.geometry.position + drift_term + diffusion_term;

                // If the particle travels for time dt at the assumed velocity, it ends up at new_position
                Vec2 assumed_velocity = (new_position - particle.geometry.position) / dt;

                particle.geometry.velocity = assumed_velocity;

                real_t time_left = dt;
                while (true) {
                    // Calculate the next cell crossing
                    CellCrossing crossing = mesh.get_next_crossing(particle.geometry);

                    if (crossing.t < time_left) {
                        // We will hit the boundary before the time is up
                        num_crossings++;
                        time_left -= crossing.t;
                        particle.geometry.position = crossing.position;

                        // Tally
                        real_t v = mesh.get_volume(particle.geometry.cell_id);
#pragma omp atomic
                        density[particle.geometry.cell_id] +=
                            crossing.t * weight_old / v * intensity_multiplier;

                        if (crossing.next_cell_id != -1) {
                            // Regular grid crossing, move the particle to the
                            // next cell
                            particle.geometry.cell_id = crossing.next_cell_id;
                        }
                        else {
                            // Boundary crossing
                            BoundaryCondition bc =
                                boundary_conditions[crossing.boundary_id];

                            if (bc == BoundaryCondition::REFLECTIVE) {
                                // Reflect the particle
                                Vec2 normal = mesh.get_surface_normal(
                                    particle.geometry.cell_id,
                                    crossing.surface_id);
                                particle.geometry.velocity =
                                    particle.geometry.velocity -
                                    2 * dot(particle.geometry.velocity, normal) *
                                        normal;
                            }
                            else if (bc == BoundaryCondition::ABSORBING) {
                                // Absorb the particle
                                particle.weight = 0;
                                break;
                            }
                            else {
                                throw std::runtime_error(
                                    "Using a boundary condition that is not "
                                    "(yet) implemented");
                            }
                        }
                    }
                    else {
                        // We will not hit the boundary before the time is up
                        // This is the final iteration
                        particle.geometry.position +=
                            time_left * particle.geometry.velocity;

                        // Tally
                        real_t v = mesh.get_volume(particle.geometry.cell_id);
#pragma omp atomic
                        density[particle.geometry.cell_id] +=
                            time_left * weight_old / v * intensity_multiplier;

                        break;
                    }
                }

                // Check for charge exchange
                // if(kinetic_flight_time < dt){
                //     num_collisions++;
                //     // Generate new velocity from plasma background (BGK
                //     // operator)
                //     real_t stddev = std::sqrt(background.plasma_temp[cell_id_old]);
                //     std::normal_distribution<real_t> dist_vx(
                //         background.plasma_vel_x[cell_id_old], stddev);
                //     std::normal_distribution<real_t> dist_vy(
                //         background.plasma_vel_y[cell_id_old], stddev);
                //     particle.geometry.velocity = Vec2(dist_vx(particle.rng), dist_vy(particle.rng));
                // }

                // Ionisation
                // particle.weight = particle.weight * std::exp(-dt * constants::K_i * background.plasma_dens[cell_id_old]);

                if (particle.weight < 1.0e-4) {
                    break;
                }
            }

#pragma omp atomic
            stats.num_crossings += num_crossings;
#pragma omp atomic
            stats.num_collisions += num_collisions;
        }
    }
}

}  // namespace eiron2