#include "neptune/simulator/simulator.hpp"

#include <omp.h>

#include "neptune/estimation/estimators.hpp"
#include "neptune/geometry/voxelmesh.hpp"
#include "neptune/physics/collision_processes.hpp"
#include "neptune/physics/update_particle.hpp"
namespace neptune {

void simulate(
    Mesh& mesh,
    const vector<std::unique_ptr<BGKScatteringProcess>>& bgk_processes,
    const vector<std::unique_ptr<IonisationProcess>>& ion_processes,
    vector<std::unique_ptr<Source>>& sources,
    Estimator& estimator,
    vector<std::unique_ptr<BoundaryCondition>>& boundary_conditions,
    SimulationParams params,
    ExecutionStatistics& stats)
{
    std::uniform_real_distribution<real_t> r_unif(0, 1);
    std::normal_distribution<real_t> r_normal(0, 1);

    // Calculate total number of particles
    for (int i_source = 0; i_source < sources.size(); i_source++) {
        stats.total_particles += sources[i_source]->num_particles();
    }

    index_t num_crossings = 0;
    index_t num_collisions = 0;

    Background background;

    // History based simulation (because its easy)
#pragma omp parallel reduction(+ : num_crossings, num_collisions) \
    reduction(EstimatorReduction : estimator)
    for (int i_source = 0; i_source < sources.size(); i_source++) {
        #pragma omp for
        for (int i_particle = 0;
             i_particle < sources[i_source]->num_particles(); i_particle++) {
            // Generate particle
            Particle particle =
                sources[i_source]->generate_particle(i_particle, mesh, background);

            for (index_t i_step = 0; i_step < params.max_steps; i_step++) {
                Particle particle_old = particle;
                // Simulate the particle up to the next event
                EventType event = update_particle(
                    mesh, bgk_processes, ion_processes, boundary_conditions,
                    particle, params.t_final);

                real_t R_i = 0.0;
                if (ion_processes.size() > 0) {
                    R_i = ion_processes[0]->rate[particle_old.geometry.cell_id];
                }
                estimator.add_track_length(particle, particle_old, mesh, R_i);

                if (particle.weight < 1.0e-4) {
                    break;
                }

                if (event == EventType::CELL_CROSSING) {
                    num_crossings++;
                }
                else if (event == EventType::CHARGE_EXCHANGE) {
                    num_collisions++;
                }
                else if (event == EventType::FINAL_TIME) {
                    estimator.add_final(particle, mesh);
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