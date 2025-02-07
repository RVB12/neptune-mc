#include <omp.h>

#include "eigrene/geometry/voxelmesh.hpp"
#include "eigrene/physics/linesource.hpp"
#include "eigrene/physics/physicalconstants.hpp"
#include "eigrene/physics/update_particle.hpp"
#include "eigrene/simulator/simulator.hpp"
namespace eigrene {

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

    // Charge exchange rate coefficient
    real_t K_cx = constants::calculate_K_cx(background.plasma_temp[cell_id]);

    // Charge exchange rate
    real_t R_cx = K_cx * background.plasma_dens[cell_id];

    // Kinetic flight time: part of the timestep using normal kinetic model
    real_t kinetic_flight_time = -std::log(r_unif(p.rng)) / (R_cx);

    if (kinetic_flight_time > dt2) {
        kinetic_flight_time = dt2;
    }

    real_t stddev =
        constants::calculate_BGK_stddev(background.plasma_temp[cell_id]);
    std::normal_distribution<real_t> dist_vx(background.plasma_vel_x[cell_id],
                                             stddev);
    std::normal_distribution<real_t> dist_vy(background.plasma_vel_y[cell_id],
                                             stddev);
    Vec2 new_velocity(dist_vx(p.rng), dist_vy(p.rng));

    // Diffusive flight time: part of the timestep using diffusive model
    real_t diffusive_flight_time = std::max(0.0, dt2 - kinetic_flight_time);

    // Correspondence with formulas (3.1--3.2) in KDMC paper:
    // T / sigma = plasma_background_variance
    // D = Variance(Delta X| V) / 2 / diffusive_flight_time
    // charge_exchange_rate = sigma / epsilon^2
    Vec2 u = Vec2(background.plasma_vel_x[cell_id],
                  background.plasma_vel_y[cell_id]);

    Vec2 drift_term = u * diffusive_flight_time;
    drift_term += (new_velocity - u) / R_cx * (1 - std::exp(-R_cx * dt));

    // Note, this still needs to be updated to take the velocity
    // correlation into account
    real_t D0 = background.plasma_temp[cell_id] * constants::eV_to_Joules /
                constants::m_p / R_cx;

    D0 *= (2 * std::exp(-R_cx * diffusive_flight_time) +
           R_cx * diffusive_flight_time +
           R_cx * diffusive_flight_time *
               std::exp(-R_cx * diffusive_flight_time) -
           2);

    D0 += (new_velocity - u).norm2() / (R_cx * R_cx) *
          (1 -
           2 * R_cx * diffusive_flight_time *
               std::exp(-R_cx * diffusive_flight_time) -
           std::exp(-2 * R_cx * diffusive_flight_time));

    // assert(D0 >= 0);
    if (D0 < 0) {
        D0 = 0;
    }

    Vec2 diffusion_term = Vec2(r_normal(p.rng), r_normal(p.rng)) *
                          std::sqrt(2.0 * diffusive_flight_time * D0);

    Vec2 new_position = p.geometry.position + drift_term + diffusion_term +
                        kinetic_flight_time * p.geometry.velocity;

    // Check if the new position is outside the mesh
    index_t new_cell_id = mesh.get_cell_index(new_position);
    if (new_cell_id == -1) {
        // Particle is outside the mesh
        // Good handling of the boundary is still in the works
        // For now, we just switch to normal MC

        real_t t_next = p.t + dt2;
        while (p.t < t_next) {
            update_particle_fixed_time(mesh, background, boundary_conditions, p,
                                       num_crossings, num_collisions, t_next);
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

    if(kinetic_flight_time < dt2){
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
                   real_t t_final,
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
                                         num_collisions_local, t_final, dt);

                    // Tally
                    if (particle.t >= t_final) {
                        real_t v = mesh.get_volume(particle.geometry.cell_id);
                        density_local[particle.geometry.cell_id] +=
                            intensity_multiplier / v;
                        break;
                    }

                    if (particle.weight < 1.0e-4) {
                        break;
                    }
                }
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

}  // namespace eigrene