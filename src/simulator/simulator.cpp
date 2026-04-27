#include "neptune/simulator/simulator.hpp"

#include <omp.h>

#include "neptune/collisions/abstract_collisions.hpp"
#include "neptune/estimation/estimator.hpp"
#include "neptune/geometry/voxelmesh.hpp"
#include "neptune/physics/update_particle.hpp"
namespace neptune {

void simulate(
    const Mesh& mesh,
    const Background& background,
    const vector<std::unique_ptr<ChargeExchangeCollision>>& cx_collisions,
    const vector<std::unique_ptr<IonisationCollision>>& ion_collisions,
    vector<std::unique_ptr<Source>>& sources,
    EstimatorGroup& estimator_group,
    vector<std::unique_ptr<BoundaryCondition>>& boundary_conditions,
    SimulationParams params,
    ExecutionStatistics& stats)
{
    std::uniform_real_distribution<real_t> r_unif(0, 1);
    std::normal_distribution<real_t> r_normal(0, 1);

    // Calculate total number of particles
    index_t total_particles = 0;
    for (int i_source = 0; i_source < sources.size(); i_source++) {
        total_particles += sources[i_source]->num_particles();
    }

    index_t num_crossings = 0;
    index_t num_collisions = 0;

    stats.total_particles = total_particles * params.n_batches;
    vector<size_t> cumulative_particles_per_source(sources.size() + 1, 0);
    for (int i_source = 0; i_source < sources.size(); i_source++) {
        cumulative_particles_per_source[i_source + 1] =
            cumulative_particles_per_source[i_source] +
            sources[i_source]->num_particles();
    }

    rng_key_t base_key{static_cast<uint32_t>(params.seed >> 32), static_cast<uint32_t>(params.seed)};

    // History based simulation (because its easy)
#pragma omp parallel reduction(+ : num_crossings, num_collisions) \
    reduction(EstimatorGroupReduction : estimator_group)
    for(index_t i_batch = 0; i_batch < params.n_batches; i_batch++) {
#pragma omp for schedule(guided, 10)
        for (index_t global_particle_index = 0;
             global_particle_index < total_particles; global_particle_index++) {
            // Find source and local particles index
            index_t i_source = 0;
            // binary search to find the source index
            index_t left = 0;
            index_t right = sources.size();
            while (left < right) {
                index_t mid = left + (right - left) / 2;
                if (cumulative_particles_per_source[mid] <=
                    global_particle_index) {
                    left = mid + 1;
                }
                else {
                    right = mid;
                }
            }
            i_source = left - 1;

            Particle particle = sources[i_source]->generate_particle(
                global_particle_index * (i_batch + 1), base_key, mesh, background);

            for (index_t i_step = 0; i_step < params.max_steps; i_step++) {
                Particle particle_old = particle;
                // Simulate the particle up to the next event
                EventType event = update_particle(
                    mesh, background, cx_collisions, ion_collisions,
                    boundary_conditions, particle, params.t_final);

                estimator_group.tally_track_length(particle, particle_old, mesh,
                                                   background, ion_collisions);

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
                    estimator_group.tally_final(particle, mesh);
                    particle.weight = 0;
                    break;
                }
            }
        }
        estimator_group.finalize_batch();
    }
    estimator_group.finalize();


    stats.num_crossings = num_crossings;
    stats.num_collisions = num_collisions;
}

}  // namespace neptune