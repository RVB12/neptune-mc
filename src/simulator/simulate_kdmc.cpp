#include <omp.h>

#include "neptune/geometry/voxelmesh.hpp"
#include "neptune/physics/linesource.hpp"
#include "neptune/physics/physicalconstants.hpp"
#include "neptune/physics/update_particle.hpp"
#include "neptune/simulator/simulator.hpp"
namespace neptune {

void update_particle_kdmc(const Mesh& mesh,
                          const PlasmaBackground& background,
                          const vector<BoundaryCondition>& boundary_conditions,
                          Particle& p,
                          index_t& num_crossings,
                          index_t& num_collisions,
                          real_t t_final,
                          real_t dt)
{
    std::uniform_real_distribution<real_t> r_unif(0, 1);
    std::normal_distribution<real_t> r_normal(0, 1);

    index_t cell_id = p.geometry.cell_id;  // Current cell id of the particle

    // Limit the time step to the final time
    real_t dt2 = std::min(dt, t_final - p.t);

    // Charge exchange rate coefficient at the starting cell
    real_t K_cx0 = constants::calculate_K_cx(background.plasma_temp[cell_id]);

    // Charge exchange rate at the starting cell
    real_t R_cx0 = K_cx0 * background.plasma_dens[cell_id];

    // Kinetic flight time: part of the timestep using normal kinetic model
    real_t kinetic_flight_time = -std::log(r_unif(p.rng)) / (R_cx0);

    if (kinetic_flight_time > dt2) {
        kinetic_flight_time = dt2;
    }

    // Diffusive flight time: part of the timestep using diffusive model
    // For now, this is a simple and relatively inefficient way to sample
    // Should be improved in the future
    // real_t diffusive_flight_time = 0;
    // while (true) {
    //     real_t time2 = -std::log(r_unif(p.rng)) / (R_cx0);
    //     if (kinetic_flight_time + diffusive_flight_time + time2 > dt2) {
    //         break;
    //     }
    //     diffusive_flight_time += time2;
    // }

    // real_t kinetic_flight_time2 =
    //     dt2 - kinetic_flight_time - diffusive_flight_time;

    real_t kinetic_flight_time2 = -std::log(r_unif(p.rng)) / (R_cx0);
    if (kinetic_flight_time2 > dt2 - kinetic_flight_time) {
        kinetic_flight_time2 = dt2 - kinetic_flight_time;
    }
    // real_t kinetic_flight_time2 = 0;
    real_t diffusive_flight_time =
        dt2 - kinetic_flight_time - kinetic_flight_time2;

    Vec3 pos2 = p.geometry.position + kinetic_flight_time * p.geometry.velocity;

    // Check if the new position is outside the mesh
    index_t cell_id2 = mesh.get_cell_index(pos2);
    bool do_mc = false;
    if (cell_id2 == -1) {
        do_mc = true;
        // TODO, this is just a kinetic step, we could properly handle
        // the boundary here
    }

    // std::cout << "dt2: " << dt2 << std::endl;
    // std::cout << "kinetic_flight_time: " << kinetic_flight_time << std::endl;
    // std::cout << "diffusive_flight_time: " << diffusive_flight_time <<
    // std::endl; std::cout << "kinetic_flight_time2: " << kinetic_flight_time2
    // << std::endl; std::cout << std::endl;

    // real_t diffusive_flight_time = dt2 - kinetic_flight_time;
    // real_t kinetic_flight_time2 = 0;

    //
    // Calculate the drift and diffusion terms.
    // This is equation 3.1 -- 3.6 in the KDMC paper.
    // We don't do scaling with epsilon here, it should be safe
    // to just remove epsilon from the equations.
    //

    // Charge exchange rate coefficient at the intermediate cell
    real_t K_cx = constants::calculate_K_cx(background.plasma_temp[cell_id]);
    real_t R_cx = K_cx * background.plasma_dens[cell_id];
    real_t R_cx_theta = R_cx * diffusive_flight_time;

    Vec3 u =
        Vec3(background.plasma_vel_x[cell_id], background.plasma_vel_y[cell_id],
             background.plasma_vel_z[cell_id]);

    Vec3 drift_term = u * diffusive_flight_time;

    // TODO: add spatial derivative term

    real_t D = background.plasma_temp[cell_id] * constants::eV_to_Joules /
               constants::m_p / (R_cx * R_cx_theta);
    D *= (std::exp(-R_cx_theta) - 1 + R_cx_theta);

    // real_t D_diffusive = background.plasma_temp[cell_id] *
    //                      constants::eV_to_Joules / constants::m_p / R_cx;

    // real_t D_kinetic = background.plasma_temp[cell_id] *
    //                    constants::eV_to_Joules / constants::m_p *
    //                    std::sqrt(diffusive_flight_time) * 0.5;

    // real_t D = std::exp(-R_cx_theta*0.2) * D_kinetic +
    //            (1 - std::exp(-R_cx_theta*0.2)) * D_diffusive;
    // real_t D = D_kinetic;
    if (D < 0) {
        D = 0;
    }
    D = std::sqrt(2 * D * diffusive_flight_time);

    // real_t D = background.plasma_temp[cell_id] *
    //            constants::eV_to_Joules / constants::m_p / (R_cx * R_cx);
    // D *= (2 * std::exp(-R_cx_theta) + R_cx_theta + R_cx_theta *
    // std::exp(-R_cx_theta) - 2);

    Vec3 diffusion_term;
    diffusion_term.x = D * r_normal(p.rng);
    diffusion_term.y = D * r_normal(p.rng);
    diffusion_term.z = D * r_normal(p.rng);

    Vec3 new_position =
        p.geometry.position + kinetic_flight_time * p.geometry.velocity;

    if (diffusive_flight_time > 0) {
        new_position += drift_term + diffusion_term;
    }

    Vec3 new_velocity;
    if (kinetic_flight_time2 > 0) {
        // Draw new velocity (todo, update the cell id)
        real_t stddev =
            constants::calculate_BGK_stddev(background.plasma_temp[cell_id]);
        std::normal_distribution<real_t> dist_vx(
            background.plasma_vel_x[cell_id], stddev);
        std::normal_distribution<real_t> dist_vy(
            background.plasma_vel_y[cell_id], stddev);
        std::normal_distribution<real_t> dist_vz(
            background.plasma_vel_z[cell_id], stddev);

        new_velocity = Vec3(dist_vx(p.rng), dist_vy(p.rng), dist_vz(p.rng));

        new_position += new_velocity * kinetic_flight_time2;
    }

    // Don't just check if the new position is outside the mesh
    // but check if any position that could realistically be drawn
    // is outside the mesh
    // we check the 99.9% confidence interval of the position
    VoxelMesh vmesh = dynamic_cast<const VoxelMesh&>(mesh);

    Vec3 temp_pos = p.geometry.position + drift_term +
                    kinetic_flight_time * p.geometry.velocity;
    // std::cout << "D" << D << std::endl;
    Vec3 temp_diff(D, D, D);
    Vec3 temp_min = temp_pos - 3.291 * temp_diff;
    Vec3 temp_max = temp_pos + 3.291 * temp_diff;
    if (temp_min.x < 0 or temp_max.x > vmesh.Lx or temp_min.y < 0 or
        temp_max.y > vmesh.Ly or temp_min.z < 0 or temp_max.z > vmesh.Lz) {
        do_mc = true;
    }

    // // now do the same including the final flight
    // temp_pos = p.geometry.position + drift_term + kinetic_flight_time *
    // p.geometry.velocity + kinetic_flight_time2 * new_velocity; temp_min =
    // temp_pos - 2.57 * temp_diff; temp_max = temp_pos + 2.57 * temp_diff;
    // if(temp_min.x < 0 or temp_max.x > vmesh.Lx or
    //    temp_min.y < 0 or temp_max.y > vmesh.Ly or
    //    temp_min.z < 0 or temp_max.z > vmesh.Lz)
    // {
    //     do_mc = true;
    // }

    // Check if the new position is outside the mesh
    index_t new_cell_id = mesh.get_cell_index(new_position);
    if (new_cell_id == -1) {
        do_mc = true;
    }

    if (do_mc) {
        // Particle is outside the mesh
        // Good handling of the boundary is still in the works
        // For now, we just switch to normal MC

        real_t t_next = p.t + dt2;
        while (p.t < t_next) {
            EventType event = update_particle(mesh, background, boundary_conditions, p, t_next);
            if(event == EventType::FINAL_TIME) {
                break;
            }
            else if (event == EventType::CELL_CROSSING) {
                num_crossings++;
            }
            else if (event == EventType::CHARGE_EXCHANGE) {
                num_collisions++;
            }
            if (p.weight < 1.0e-4) {
                break;
            }
        }
    }
    else {
        p.geometry.position = new_position;
        p.geometry.cell_id = new_cell_id;
        p.t += dt2;
    }

    if (kinetic_flight_time < dt2) {
        num_collisions++;
        p.geometry.velocity = new_velocity;
    }
}

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
#pragma omp parallel
        {
            index_t num_crossings_local = 0;
            index_t num_collisions_local = 0;

            // This should really be changed to a reduction clause instead of
            // separate local variables
            vector<real_t> density_local(mesh.n_cells(), 0);
#pragma omp for
            for (int i_particle = 0;
                 i_particle < sources[i_source]->num_particles();
                 i_particle++) {
                real_t intensity_multiplier =
                    sources[i_source]->intensity() /
                    sources[i_source]->num_particles();

                // Generate particle
                Particle particle =
                    sources[i_source]->generate_particle(i_particle);
                particle.geometry.cell_id =
                    mesh.get_cell_index(particle.geometry.position);

                if (particle.geometry.cell_id == -1) {
                    throw std::runtime_error(
                        "Particle generated outside of mesh");
                }

                for (index_t i_step = 0; i_step < params.max_steps; i_step++) {
                    // Simulate the particle up to the next event
                    update_particle_kdmc(mesh, background, boundary_conditions,
                                         particle, num_crossings_local,
                                         num_collisions_local, params.t_final,
                                         dt);

                    // Tally
                    if (particle.t >= params.t_final) {
                        real_t v = mesh.get_volume(particle.geometry.cell_id);
                        density_local[particle.geometry.cell_id] +=
                            intensity_multiplier / v;
                        break;
                    }

                    if (particle.weight < 1.0e-4) {
                        break;
                    }
                }

#ifdef _OPENMP
                auto thread_id = omp_get_thread_num();
                if (thread_id == 0) {
                    if (i_particle % 1000 == 0) {
                        std::cout.flush();
                        std::cout
                            << "Progress: "
                            << i_particle /
                                   ((float)sources[i_source]->num_particles()) *
                                   100
                            << "% \r";
                    }
                }
#endif
            }

            for (index_t i = 0; i < mesh.n_cells(); i++) {
#pragma omp atomic
                density[i] += density_local[i];
            }
#pragma omp atomic
            stats.num_crossings += num_crossings_local;
#pragma omp atomic
            stats.num_collisions += num_collisions_local;
        }
    }
}

// void simulate_kdmc(Mesh& mesh,
//                    PlasmaBackground& background,
//                    vector<std::unique_ptr<Source>>& sources,
//                    vector<real_t>& density,
//                    vector<BoundaryCondition>& boundary_conditions,
//                    SimulationParams params,
//                    ExecutionStatistics& stats,
//                     real_t t_final,
//                    real_t dt)
// {
//     std::uniform_real_distribution<real_t> r_unif(0, 1);
//     std::normal_distribution<real_t> r_normal(0, 1);

//     // initialize the density array to zero
//     for (index_t i = 0; i < mesh.n_cells(); i++) {
//         density[i] = 0;
//     }

//     // History based simulation (because its easy)
//     for (int i_source = 0; i_source < sources.size(); i_source++) {
// #pragma omp parallel for
//         for (int i_particle = 0;
//              i_particle < sources[i_source]->num_particles(); i_particle++) {
//             real_t intensity_multiplier = sources[i_source]->intensity() /
//                                           sources[i_source]->num_particles();
//             index_t num_crossings = 0;
//             index_t num_collisions = 0;

//             // Generate particle
//             Particle particle =
//                 sources[i_source]->generate_particle(i_particle);
//             particle.geometry.cell_id =
//                 mesh.get_cell_index(particle.geometry.position);

//             if (particle.geometry.cell_id == -1) {
//                 throw std::runtime_error("Particle generated outside of
//                 mesh");
//             }

//             for (index_t i_step = 0; i_step < params.max_steps; i_step++) {
//                 real_t weight_old = particle.weight;
//                 // real_t t_old = particle.t;
//                 index_t cell_id_old = particle.geometry.cell_id;

//                 // Charge exchange rate coefficient
//                 real_t K_cx = constants::calculate_K_cx(
//                     background.plasma_temp[cell_id_old]);

//                 // Charge exchange rate
//                 real_t R_cx = K_cx * background.plasma_dens[cell_id_old];

//                 // Calculate time until next charge exchange
//                 real_t kinetic_flight_time =
//                     -std::log(r_unif(particle.rng)) / (R_cx);

//                 bool charge_exchange_occurred = false;
//                 Vec2 new_velocity = particle.geometry.velocity;
//                 if (kinetic_flight_time <= dt) {
//                     charge_exchange_occurred = true;
//                     real_t stddev = constants::calculate_BGK_stddev(
//                         background.plasma_temp[cell_id_old]);
//                     std::normal_distribution<real_t> dist_vx(
//                         background.plasma_vel_x[cell_id_old], stddev);
//                     std::normal_distribution<real_t> dist_vy(
//                         background.plasma_vel_y[cell_id_old], stddev);
//                     new_velocity =
//                         Vec2(dist_vx(particle.rng), dist_vy(particle.rng));
//                 } else {
//                     kinetic_flight_time = dt;
//                     new_velocity = particle.geometry.velocity;
//                 }

//                 // Fill up remaining time with diffusion
//                 real_t diffusive_flight_time =
//                     std::max(0.0, dt - kinetic_flight_time);

//                 // Correspondence with formulas (3.1--3.2) in KDMC paper:
//                 // T / sigma = plasma_background_variance
//                 // D = Variance(Delta X| V) / 2 / diffusive_flight_time
//                 // charge_exchange_rate = sigma / epsilon^2
//                 Vec2 u = Vec2(background.plasma_vel_x[cell_id_old],
//                               background.plasma_vel_y[cell_id_old]);

//                 Vec2 drift_term = u * diffusive_flight_time;
//                 // drift_term += (new_velocity - u) / R_cx *
//                 //               (1 - std::exp(-R_cx *
//                 diffusive_flight_time));
//                 // TODO: add spatial derivative term

//                 // Note, this still needs to be updated to take the velocity
//                 // correlation into account
//                 real_t D0 = background.plasma_temp[cell_id_old] *
//                             constants::eV_to_Joules / constants::m_p / R_cx;

//                 // D0 *= (2 * std::exp(-R_cx * diffusive_flight_time) +
//                 //        R_cx * diffusive_flight_time +
//                 //        R_cx * diffusive_flight_time *
//                 //            std::exp(-R_cx * diffusive_flight_time) -
//                 //        2);

//                 // D0 += (new_velocity - u).norm2() / (R_cx * R_cx) *
//                 //       (1 -
//                 //        2 * R_cx * diffusive_flight_time *
//                 //            std::exp(-R_cx * diffusive_flight_time) -
//                 //        std::exp(-2 * R_cx * diffusive_flight_time));

//                 // assert(D0 >= 0);
//                 if(D0 < 0){
//                     D0 = 0;
//                 }
//                 Vec2 diffusion_term =
//                     Vec2(r_normal(particle.rng), r_normal(particle.rng)) *
//                     std::sqrt(2.0 * diffusive_flight_time * D0);

//                 Vec2 new_position =
//                     particle.geometry.position +
//                     kinetic_flight_time * particle.geometry.velocity +
//                     drift_term + diffusion_term;

//                 Vec2 old_velocity = particle.geometry.velocity;

//                 // We have calculated, now need to trace through the mesh
//                 // for tallying. We do this by pretending the particles moves
//                 // to the new position at a constant velocity
//                 Vec2 assumed_velocity =
//                     (new_position - particle.geometry.position) / dt;
//                 particle.geometry.velocity = assumed_velocity;
//                 real_t time_left = dt;
//                 while (true) {
//                     // Calculate the next cell crossing
//                     CellCrossing crossing =
//                         mesh.get_next_crossing(particle.geometry);

//                     if (crossing.t < time_left) {
//                         // We will hit the boundary before the time is up
//                         num_crossings++;
//                         time_left -= crossing.t;
//                         particle.geometry.position = crossing.position;

//                         // Tally
//                         real_t v =
//                         mesh.get_volume(particle.geometry.cell_id);
// #pragma omp atomic
//                         density[particle.geometry.cell_id] +=
//                             crossing.t * weight_old / v *
//                             intensity_multiplier;

//                         if (crossing.next_cell_id != -1) {
//                             // Regular grid crossing, move the particle to
//                             the
//                             // next cell
//                             particle.geometry.cell_id =
//                             crossing.next_cell_id;
//                             particle.geometry.boundary_id =
//                                 crossing.next_surface_id;
//                         }
//                         else {
//                             // Boundary crossing
//                             BoundaryCondition bc =
//                                 boundary_conditions[crossing
//                                                         .boundary_condition_id];
//                             particle.geometry.boundary_id =
//                             crossing.surface_id;

//                             if (bc == BoundaryCondition::REFLECTIVE) {
//                                 // Reflect the particle
//                                 Vec2 normal = mesh.get_surface_normal(
//                                     particle.geometry.cell_id,
//                                     crossing.surface_id);
//                                 particle.geometry.velocity =
//                                     particle.geometry.velocity -
//                                     2 *
//                                         dot(particle.geometry.velocity,
//                                             normal) *
//                                         normal;
//                             }
//                             else if (bc == BoundaryCondition::ABSORBING) {
//                                 // Absorb the particle
//                                 particle.weight = 0;
//                                 break;
//                             }
//                             else {
//                                 throw std::runtime_error(
//                                     "Using a boundary condition that is not "
//                                     "(yet) implemented");
//                             }
//                         }
//                     }
//                     else {
//                         // We will not hit the boundary before the time is up
//                         // This is the final iteration
//                         particle.geometry.position +=
//                             time_left * particle.geometry.velocity;
//                         particle.geometry.boundary_id = -1;

//                         // Tally
//                         real_t v =
//                         mesh.get_volume(particle.geometry.cell_id);
// #pragma omp atomic
//                         density[particle.geometry.cell_id] +=
//                             time_left * weight_old / v *
//                             intensity_multiplier;

//                         break;
//                     }
//                 }

//                 particle.geometry.velocity = old_velocity;

//                 // Check for charge exchange
//                 if (charge_exchange_occurred) {
//                     num_collisions++;
//                     // Generate new velocity from plasma background (BGK
//                     // operator)
//                     particle.geometry.velocity = new_velocity;
//                 }

//                 // Ionisation
//                 particle.weight = particle.weight *
//                            std::exp(-dt *
//                                     constants::calculate_K_i(
//                                         background.plasma_temp[cell_id_old])
//                                         *
//                                     background.plasma_dens[cell_id_old]);

//                 particle.t += dt;
//                 if (particle.weight < 1.0e-4) {
//                     // std::cout << "final time: " << particle.t <<
//                     std::endl; break;
//                 }
//             }

// #pragma omp atomic
//             stats.num_crossings += num_crossings;
// #pragma omp atomic
//             stats.num_collisions += num_collisions;
//         }
//     }
// }

}  // namespace neptune