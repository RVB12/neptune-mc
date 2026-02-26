#include <omp.h>

#include "neptune/common/utils.hpp"
#include "neptune/geometry/trianglemesh.hpp"
#include "neptune/geometry/voxelmesh.hpp"
#include "neptune/physics/physicalconstants.hpp"
#include "neptune/physics/update_particle.hpp"
#include "neptune/simulator/simulator.hpp"
namespace neptune {

EventType update_particle_arcmc(
    const Mesh& mesh,
    const vector<std::unique_ptr<BGKScatteringProcess>>& bgk_processes,
    const vector<std::unique_ptr<IonisationProcess>>& ion_processes,
    const vector<std::unique_ptr<BoundaryCondition>>& boundary_conditions,
    Particle& p,
    real_t t_final,
    index_t M,
    real_t& V2_eff,
    real_t& V2_eff_x)
{
    std::uniform_real_distribution<real_t> r_unif(0, 1);
    std::normal_distribution<real_t> r_normal(0, 1);

    index_t cell_id = p.geometry.cell_id;  // Current cell id of the particle

    // Calculate total collision rate (charge exchange + neutral-neutral)

    real_t collision_rate = 0;
    for (int i = 0; i < bgk_processes.size(); i++) {
        collision_rate += bgk_processes[i]->rate[cell_id];
    }

    // random number to determine time till collision
    real_t eps = r_unif(p.rng);
    // time till collision
    real_t tc = p.u / collision_rate;

    // Calculate the next cell crossing
    CellCrossing crossing = mesh.get_next_crossing(p.geometry);

    // Check which one is earlier, charge exchange, cell crossing or final time
    if (t_final > 0) {
        real_t time_left = t_final - p.t;
        if (time_left < tc && time_left < crossing.t) {
            // We will not hit the boundary before the time is up
            // This is the final iteration
            p.geometry.position += time_left * p.geometry.velocity;
            p.geometry.boundary_id = -1;
            p.t = t_final;
            for (int i = 0; i < ion_processes.size(); i++) {
                p.weight *=
                    std::exp(-time_left * ion_processes[i]->rate[cell_id]);
            }
            return EventType::FINAL_TIME;
        }
    }
    real_t dt;
    EventType event_type;
    if (tc < crossing.t) {
        event_type = EventType::CHARGE_EXCHANGE;
        // Move the particle until the collision
        p.geometry.position += tc * p.geometry.velocity;

        // A collision has occurred, normally, we would now determine which
        // collision has occurred but ARC-MC currently only supports one
        // collision process.

        // Draw a new "velocity" according to the ARCMC method

        double S, W;
        if (M <= 4) {
            std::exponential_distribution<double> exp_dist(1.);
            S = 0.0;
            W = 0.0;
            for (uint64_t j = 0; j < M; ++j) {
                double tau = exp_dist(p.rng);
                S += tau;
                W += tau * tau;
            }
            W /= S * S;
        }
        else {
            std::gamma_distribution<double> gamma_dist_S((double)M, 1.0);
            double shape_w = (M - 1.) * (M + 2.) * (M + 3.) / (4. * M * M);
            double rate_w = (M + 1.) * (M + 2.) * (M + 3.) / (4. * M);
            double scale_w = 1.0 / rate_w;
            std::gamma_distribution<double> gamma_dist_W(shape_w, scale_w);

            S = gamma_dist_S(p.rng);
            W = 1.0 / M + gamma_dist_W(p.rng);
        }

        std::normal_distribution<double> normal_dist(0., 1.);
        double xi_x = normal_dist(p.rng);
        double xi_y = normal_dist(p.rng);
        double xi_z = normal_dist(p.rng);

        Vec3 xi(xi_x, xi_y, xi_z);

        // Now we will pretend that the particle moves in a straight line
        // with an assumed velocity
        Vec3 drift(bgk_processes[0]->bgk_mean_x[cell_id],
                   bgk_processes[0]->bgk_mean_y[cell_id],
                   bgk_processes[0]->bgk_mean_z[cell_id]);
        real_t bgk_std = bgk_processes[0]->bgk_stddev[cell_id];

        p.geometry.velocity = drift + std::sqrt(W) * bgk_std * xi;

        // V2_eff and V2_eff_x are used for the energy estimator,
        // since we can't just use v_eff^2 for the energy estimate, we need to calculate
        // which is a standin for v_eff^2 that accounts for the fact that the actual velocity varies more

        // Very strange behavior here, the following two lines should be better,
        // since it corrects the mean while taking the actual velocity into account,
        // but it seems to lead to worse estimates.
        // Just setting V2_eff to its expected value seems to work better, maybe it has to do with some covariance or something?
        // V2_eff = p.geometry.velocity.norm2() + 3 * bgk_std * bgk_std * (1.0 - W);
        // V2_eff_x = p.geometry.velocity.x * p.geometry.velocity.x + bgk_std * bgk_std * (1.0 - W);
        
        V2_eff = drift.norm2() + 3 * bgk_std * bgk_std;
        V2_eff_x = drift.x * drift.x + bgk_std * bgk_std;

        p.u = S;

        // Reset boundary id (we are no longer at a boundary)
        p.geometry.boundary_id = -1;

        dt = tc;
    }
    else {
        event_type = EventType::CELL_CROSSING;
        dt = crossing.t;
        p.u -=
            collision_rate * dt;  // Update the scaled time till next collision
        if (crossing.next_cell_id != -1) {
            // Regular grid crossing, move the particle to the next
            // cell
            p.geometry.position = crossing.position;
            p.geometry.cell_id = crossing.next_cell_id;
            p.geometry.boundary_id = crossing.next_surface_id;
        }
        else {
            // Boundary crossing
            p.geometry.position = crossing.position;
            p.geometry.boundary_id = crossing.surface_id;

            Vec3 normal = mesh.get_surface_normal(cell_id, crossing.surface_id,
                                                  crossing.position);

            boundary_conditions[crossing.boundary_condition_id]->apply(p,
                                                                       normal);
        }
    }

    // Ionisation
    for (int i = 0; i < ion_processes.size(); i++) {
        p.weight *= std::exp(-dt * ion_processes[i]->rate[cell_id]);
    }
    p.t += dt;

    return event_type;
}

void simulate_armc(
    Mesh& mesh,
    const vector<std::unique_ptr<BGKScatteringProcess>>& bgk_processes,
    const vector<std::unique_ptr<IonisationProcess>>& ion_processes,
    vector<std::unique_ptr<Source>>& sources,
    ArcmcEstimator& estimator,
    vector<std::unique_ptr<BoundaryCondition>>& boundary_conditions,
    SimulationParams params,
    ExecutionStatistics& stats,
    index_t M)
{
    std::uniform_real_distribution<real_t> r_unif(0, 1);
    std::normal_distribution<real_t> r_normal(0, 1);

    // Calculate total number of particles
    for (int i_source = 0; i_source < sources.size(); i_source++) {
        stats.total_particles += sources[i_source]->num_particles();
    }

    index_t num_crossings = 0;
    index_t num_collisions = 0;

    // History based simulation (because its easy)
#pragma omp parallel reduction(+ : num_crossings, num_collisions) \
    reduction(ArcmcEstimatorReduction : estimator)
    for (int i_source = 0; i_source < sources.size(); i_source++) {
#pragma omp for
        for (int i_particle = 0;
             i_particle < sources[i_source]->num_particles(); i_particle++) {
            // Generate particle
            Particle particle =
                sources[i_source]->generate_particle(i_particle, mesh);
            real_t V2_eff = particle.geometry.velocity.norm();
            real_t V2_eff_x = particle.geometry.velocity.x * particle.geometry.velocity.x;

            particle.intensity = sources[i_source]->intensity() /
                                 sources[i_source]->num_particles();

            for (index_t i_step = 0; i_step < params.max_steps; i_step++) {
                Particle particle_old = particle;
                // Simulate the particle up to the next event
                EventType event = update_particle_arcmc(
                    mesh, bgk_processes, ion_processes, boundary_conditions,
                    particle, params.t_final, M, V2_eff, V2_eff_x);

                estimator.add_track_length(particle, particle_old, mesh, V2_eff, V2_eff_x);

                if (particle.weight < 1.0e-4) {
                    break;
                }

                if (event == EventType::CELL_CROSSING) {
                    num_crossings++;
                }
                else if (event == EventType::CHARGE_EXCHANGE) {
                    // TODO: this assume we always do M collisions in ARCMC, but
                    // in the first iteration, we usually do just 1 collision.
                    num_collisions += M;
                }
                else if (event == EventType::FINAL_TIME) {
                    estimator.add_final(particle, mesh, V2_eff, V2_eff_x);
                    particle.weight = 0;
                    break;
                }
            }
        }
    }

    stats.num_crossings = num_crossings;
    stats.num_collisions = num_collisions;
}

}  // namespace neptune