#include "neptune/simulator/simulator.hpp"

#include <omp.h>

#include "neptune/estimation/estimators.hpp"
#include "neptune/geometry/voxelmesh.hpp"
#include "neptune/physics/linesource.hpp"
#include "neptune/physics/update_particle.hpp"
namespace neptune {

void simulate(Mesh& mesh,
              PlasmaBackground& background,
              vector<std::unique_ptr<Source>>& sources,
              Estimator& estimator,
              vector<BoundaryCondition>& boundary_conditions,
              SimulationParams params,
              ExecutionStatistics& stats)
{
    std::uniform_real_distribution<real_t> r_unif(0, 1);
    std::normal_distribution<real_t> r_normal(0, 1);

    // History based simulation (because its easy)
    for (int i_source = 0; i_source < sources.size(); i_source++) {
#pragma omp parallel
        {
            index_t num_crossings_local = 0;
            index_t num_collisions_local = 0;

            Estimator estimator_local(mesh.n_cells());

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
                    // Simulate the particle up to the next event
                    EventType event = update_particle(mesh, background, boundary_conditions,
                                    particle, params.t_final);

                    estimator_local.add_track_length(particle, particle_old,
                                                     mesh);

                    if (particle.weight < 1.0e-4) {
                        break;
                    }

                    if(event == EventType::CELL_CROSSING) {
                        num_crossings_local++;
                    }
                    else if (event == EventType::CHARGE_EXCHANGE) {
                        num_collisions_local++;
                    }
                    else if (event == EventType::FINAL_TIME) {
                        estimator_local.add_final(particle, mesh);
                        particle.weight = 0;
                        break;
                    }
                }
            }

// Combine the results of the local estimators
#pragma omp critical
            {
                estimator.combine(estimator_local);
            }

#pragma omp atomic
            stats.num_crossings += num_crossings_local;
#pragma omp atomic
            stats.num_collisions += num_collisions_local;
        }
    }
}

}  // namespace neptune