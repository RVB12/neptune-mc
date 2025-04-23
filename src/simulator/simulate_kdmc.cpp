#include <omp.h>

#include "neptune/geometry/voxelmesh.hpp"
#include "neptune/physics/linesource.hpp"
#include "neptune/physics/physicalconstants.hpp"
#include "neptune/physics/update_particle.hpp"
#include "neptune/simulator/simulator.hpp"
namespace neptune {

/**
 * @brief Simulate a particle using a diffusive fluid approximation
 */
void kdmc_fluid_step(const Mesh& mesh,
                     const PlasmaBackground& background,
                     const vector<std::unique_ptr<BoundaryCondition>>& boundary_conditions,
                     Particle& p,
                     real_t dt)
{
    std::uniform_real_distribution<real_t> r_unif(0, 1);
    std::normal_distribution<real_t> r_normal(0, 1);

    index_t cell_id = p.geometry.cell_id;  // Current cell id of the particle


    real_t K_cx = constants::calculate_K_cx(background.plasma_temp[cell_id]);
    real_t R_cx = K_cx * background.plasma_dens[cell_id];
    real_t R_cx_theta = R_cx * dt;

    Vec3 u =
        Vec3(background.plasma_vel_x[cell_id], background.plasma_vel_y[cell_id],
             background.plasma_vel_z[cell_id]);

    Vec3 drift_term = u * dt;

    real_t D = background.plasma_temp[cell_id] * constants::eV_to_Joules /
               constants::m_p / (R_cx * R_cx_theta);
    D *= (std::exp(-R_cx_theta) - 1 + R_cx_theta);
    if (D < 0) {
        D = 0;
    }
    D = std::sqrt(2 * D * dt);

    Vec3 diffusion_term;
    diffusion_term.x = D * r_normal(p.rng);
    diffusion_term.y = D * r_normal(p.rng);
    diffusion_term.z = D * r_normal(p.rng);

    Vec3 new_pos = p.geometry.position + drift_term + diffusion_term;

    index_t new_cell_id = mesh.get_cell_index(new_pos);

    // Now handle boundary conditions
    // This only works for VoxelMesh for now
    if (new_cell_id == -1) {
        return;
        VoxelMesh vmesh = dynamic_cast<const VoxelMesh&>(mesh);

        real_t A;
        // TODO: fix boundary conditions
        // if (boundary_conditions[0] == BoundaryCondition::ABSORBING) {
        //     A = -1.;
        // }
        // else {
        //     A = 1.;
        // }

        // Reflect the particle back inside the mesh
        // This may have to be done multiple times
        while (true) {
            if (new_pos.x < 0) {
                new_pos.x = -new_pos.x;
                p.weight *= A;
            }
            else if (new_pos.x > vmesh.Lx) {
                new_pos.x = 2 * vmesh.Lx - new_pos.x;
                p.weight *= A;
            }
            if (new_pos.y < 0) {
                new_pos.y = -new_pos.y;
                p.weight *= A;
            }
            else if (new_pos.y > vmesh.Ly) {
                new_pos.y = 2 * vmesh.Ly - new_pos.y;
                p.weight *= A;
            }
            if (new_pos.z < 0) {
                new_pos.z = -new_pos.z;
                p.weight *= A;
            }
            else if (new_pos.z > vmesh.Lz) {
                new_pos.z = 2 * vmesh.Lz - new_pos.z;
                p.weight *= A;
            }

            new_cell_id = mesh.get_cell_index(new_pos);
            if (new_cell_id != -1) {
                break;
            }
        }
    }

    // Update the particle
    p.geometry.position = new_pos;
    p.geometry.cell_id = new_cell_id;
    p.geometry.boundary_id = -1;
    p.t += dt;
    p.weight *= std::exp(-dt * constants::calculate_K_i(
                                  background.plasma_temp[p.geometry.cell_id]) *
                         background.plasma_dens[p.geometry.cell_id]);
}

void simulate_kdmc(Mesh& mesh,
                    PlasmaBackground& background,
                    vector<std::unique_ptr<Source>>& sources,
                    Estimator& estimator,
                    vector<std::unique_ptr<BoundaryCondition>>& boundary_conditions,
                    SimulationParams params,
                    ExecutionStatistics& stats,
                    real_t dt)
{
    std::uniform_real_distribution<real_t> r_unif(0, 1);
    std::normal_distribution<real_t> r_normal(0, 1);

    // History based simulation (because its easy)
#pragma omp parallel
    {
        index_t num_crossings_local = 0;
        index_t num_collisions_local = 0;

        Estimator estimator_local(mesh.n_cells());
        for (int i_source = 0; i_source < sources.size(); i_source++) {
#pragma omp for
            for (int i_particle = 0;
                 i_particle < sources[i_source]->num_particles();
                 i_particle++) {
                // Generate particle
                Particle particle =
                    sources[i_source]->generate_particle(i_particle);
                particle.geometry.cell_id =
                    mesh.get_cell_index(particle.geometry.position);

                particle.intensity = sources[i_source]->intensity() /
                                     sources[i_source]->num_particles();

                if (particle.geometry.cell_id == -1) {
                    throw std::runtime_error(
                        "Particle generated outside of mesh");
                }

                for (index_t i_step = 0; i_step < params.max_steps; i_step++) {
                    Particle particle_old = particle;

                    // Simulate particle kinetically until a charge exchange
                    // collision occurs

                    real_t t_max = particle.t + dt;
                    if (params.t_final > 0) {
                        t_max = std::min(t_max, params.t_final);
                    }

                    while (true) {
                        particle_old = particle;
                        EventType event = update_particle(mesh, background,
                                                          boundary_conditions,
                                                          particle, t_max);

                        estimator_local.add_track_length(particle, particle_old,
                                                         mesh);
                        if (particle.weight < 1.0e-4) {
                            break;
                        }
                        if (event == EventType::FINAL_TIME) {
                            break;
                        }
                        else if (event == EventType::CHARGE_EXCHANGE) {
                            num_collisions_local++;
                            break;
                        }
                        else if (event == EventType::CELL_CROSSING) {
                            num_crossings_local++;
                        }
                    }

                    // Now simulate the particle fluidly for the rest of the
                    // time step

                    if (particle.t < t_max) {
                        kdmc_fluid_step(mesh, background, boundary_conditions,
                                        particle, t_max - particle.t);

                        // TODO: add fluid estimator here
                    }

                    if (std::abs(particle.weight) < 1.0e-4) {
                        break;
                    }

                    // Tally
                    if (particle.t >= params.t_final) {
                        estimator_local.add_final(particle, mesh);
                        break;
                    }
                }
            }
        }
#pragma omp critical
        {
            // Add the local results to the global results
            estimator.combine(estimator_local);
        }
#pragma omp atomic
        stats.num_crossings += num_crossings_local;
#pragma omp atomic
        stats.num_collisions += num_collisions_local;
    }
}

}  // namespace neptune