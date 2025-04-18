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
                     const vector<BoundaryCondition>& boundary_conditions,
                     Particle& p,
                     real_t dt)
{
    std::uniform_real_distribution<real_t> r_unif(0, 1);
    std::normal_distribution<real_t> r_normal(0, 1);

    //
    // During a diffusive step, the particle behaves according to the SDE:
    // dX = V dt
    // dV = gamma(X) * (mu(X) - V) dt + sqrt(2 * gamma(x)) * sigma(X) dW
    //
    // here, gamma is the charge exchange rate, mu is the velocity of the plasma
    // background and sigma standard deviation of the velocity distribution
    //
    // If we assume that gamma, mu and sigma are constant during the time step,
    // we can use the analytical solution to the SDE:
    //
    // Mean of X(t) = X(0) + mu * t + (V(0) - mu) * (1 - e^(-gamma * t)) / gamma
    // Variance of X(t) = 2 * sigma^2 / gamma * (t - 2/gamma*e^(-gamma * t) +
    // 1/(2*gamma)*(1-e^(-2*gamma*t))) Mean of V(t) = V(0) * e^(-gamma * t) + mu
    // * (1 - e^(-gamma * t)) Variance of V(t) = sigma^2 * (1 - e^(-2*gamma *
    // t)) Covariance of X(t) and V(t) = sigma^2 / gamma * (1 - 2*e^(-gamma * t)
    // + e^(-2*gamma * t))
    //
    // Of course, the true gamma, mu and sigma are not constant.
    // We could use smaller timesteps, higher order methods or ad hoc
    // corrections but this is still in the works.
    //

    real_t gamma =
        constants::calculate_K_cx(background.plasma_temp[p.geometry.cell_id]) *
        background.plasma_dens[p.geometry.cell_id];
    Vec3 mu = Vec3(background.plasma_vel_x[p.geometry.cell_id],
                   background.plasma_vel_y[p.geometry.cell_id],
                   background.plasma_vel_z[p.geometry.cell_id]);
    real_t sigma = constants::calculate_BGK_stddev(
        background.plasma_temp[p.geometry.cell_id]);

    Vec3 mean_pos =
        p.geometry.position + mu * dt +
        (p.geometry.velocity - mu) * (1 - std::exp(-gamma * dt)) / gamma;
    Vec3 mean_vel = p.geometry.velocity * std::exp(-gamma * dt) +
                    mu * (1 - std::exp(-gamma * dt));
    real_t var_pos = 2 * sigma * sigma / gamma *
                     (dt - 2 / gamma * std::exp(-gamma * dt) +
                      1 / (2 * gamma) * (1 - std::exp(-2 * gamma * dt)));
    if(var_pos < 0) {
        var_pos = 0;
    }
    real_t var_vel = sigma * sigma * (1 - std::exp(-2 * gamma * dt));
    real_t cov_pos_vel =
        sigma * sigma / gamma *
        (1 - 2 * std::exp(-gamma * dt) + std::exp(-2 * gamma * dt));

    Vec3 n1 = Vec3(r_normal(p.rng), r_normal(p.rng), r_normal(p.rng));
    Vec3 n2 = Vec3(r_normal(p.rng), r_normal(p.rng), r_normal(p.rng));

    // Calculate the new position and velocity
    Vec3 new_vel = mean_vel + std::sqrt(var_vel) * n1;
    real_t corrected_var = var_pos - cov_pos_vel * cov_pos_vel / var_vel;
    if (corrected_var < 0) {
        corrected_var = 0;
    }
    Vec3 new_pos =
        mean_pos + cov_pos_vel / std::sqrt(var_vel) * n1 +
        std::sqrt(corrected_var) * n2;

    index_t new_cell_id = mesh.get_cell_index(new_pos);

    // Now handle boundary conditions
    // This only works for VoxelMesh for now
    if (new_cell_id == -1) {
        return;
        // VoxelMesh vmesh = dynamic_cast<const VoxelMesh&>(mesh);

        // real_t A;
        // if (boundary_conditions[0] == BoundaryCondition::ABSORBING) {
        //     A = -1.;
        // }
        // else {
        //     A = 1.;
        // }

        // // Reflect the particle back inside the mesh
        // // This may have to be done multiple times
        // while (true) {
        //     if (new_pos.x < 0) {
        //         new_pos.x = -new_pos.x;
        //         new_vel.x = -new_vel.x;
        //         p.weight *= A;
        //     }
        //     else if (new_pos.x > vmesh.Lx) {
        //         new_pos.x = 2 * vmesh.Lx - new_pos.x;
        //         new_vel.x = -new_vel.x;
        //         p.weight *= A;
        //     }
        //     if (new_pos.y < 0) {
        //         new_pos.y = -new_pos.y;
        //         new_vel.y = -new_vel.y;
        //         p.weight *= A;
        //     }
        //     else if (new_pos.y > vmesh.Ly) {
        //         new_pos.y = 2 * vmesh.Ly - new_pos.y;
        //         new_vel.y = -new_vel.y;
        //         p.weight *= A;
        //     }
        //     if (new_pos.z < 0) {
        //         new_pos.z = -new_pos.z;
        //         new_vel.z = -new_vel.z;
        //         p.weight *= A;
        //     }
        //     else if (new_pos.z > vmesh.Lz) {
        //         new_pos.z = 2 * vmesh.Lz - new_pos.z;
        //         new_vel.z = -new_vel.z;
        //         p.weight *= A;
        //     }

        //     new_cell_id = mesh.get_cell_index(new_pos);
        //     if (new_cell_id != -1) {
        //         break;
        //     }
        // }
    }

    // Update the particle
    p.geometry.position = new_pos;
    p.geometry.cell_id = new_cell_id;
    p.geometry.velocity = new_vel;
    p.geometry.boundary_id = -1;
    p.t += dt;
    p.weight *= std::exp(-dt * constants::calculate_K_i(
                                  background.plasma_temp[p.geometry.cell_id]) *
                         background.plasma_dens[p.geometry.cell_id]);
}

void simulate_kdmc2(Mesh& mesh,
                    PlasmaBackground& background,
                    vector<std::unique_ptr<Source>>& sources,
                    Estimator& estimator,
                    vector<BoundaryCondition>& boundary_conditions,
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