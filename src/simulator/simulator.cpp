#include "eiron2/simulator/simulator.hpp"

#include <omp.h>

#include "eiron2/geometry/voxelmesh.hpp"
#include "eiron2/physics/linesource.hpp"
#include "eiron2/physics/update_particle.hpp"
namespace eiron2 {

void simulate(Mesh& mesh,
              PlasmaBackground& background,
              vector<std::unique_ptr<Source>>& sources,
              vector<real_t>& density,
              vector<BoundaryCondition>& boundary_conditions,
              SimulationParams params,
              ExecutionStatistics& stats)
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

            // This should really be changed to a reduction clause instead of separate local variables
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
                    real_t weight_old = particle.weight;
                    real_t t_old = particle.t;
                    index_t cell_id_old = particle.geometry.cell_id;
                    // Simulate the particle up to the next event
                    update_particle(mesh, background, boundary_conditions,
                                    particle, num_crossings_local,
                                    num_collisions_local);

                    // Tally
                    real_t v = mesh.get_volume(cell_id_old);
                    density_local[cell_id_old] += (particle.t - t_old) *
                                                  weight_old / v *
                                                  intensity_multiplier;

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

}  // namespace eiron2