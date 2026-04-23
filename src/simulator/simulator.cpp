#include "neptune/simulator/simulator.hpp"

#include <omp.h>

#include "neptune/estimation/estimator.hpp"
#include "neptune/geometry/voxelmesh.hpp"
#include "neptune/collisions/abstract_collisions.hpp"
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
    for (int i_source = 0; i_source < sources.size(); i_source++) {
        stats.total_particles += sources[i_source]->num_particles();
    }

    index_t num_crossings = 0;
    index_t num_collisions = 0;

    // History based simulation (because its easy)
#pragma omp parallel reduction(+ : num_crossings, num_collisions) \
    reduction(EstimatorGroupReduction : estimator_group)
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
                    mesh, background, cx_collisions, ion_collisions, boundary_conditions,
                    particle, params.t_final);

                estimator_group.tally_track_length(particle, particle_old, mesh, background, ion_collisions);

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
    }

    stats.num_crossings = num_crossings;
    stats.num_collisions = num_collisions;
}

}  // namespace neptune