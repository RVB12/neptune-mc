#include "eiron2/simulator/simulator.hpp"
#include "eiron2/physics/update_particle.hpp"

#include "eiron2/geometry/voxelmesh.hpp"
#include "eiron2/physics/linesource.hpp"
namespace eiron2 {

void simulate(Mesh& mesh,
              PlasmaBackground& background,
              vector<std::unique_ptr<Source>>& sources,
              vector<real_t>& density,
              SimulationParams params)
{

    std::cout << "Background:" << std::endl;
    for (index_t i = 0; i < background.plasma_dens.size(); i++) {
        std::cout << background.plasma_dens[i] << std::endl;
    }
    for (index_t i = 0; i < background.plasma_dens.size(); i++) {
        std::cout << background.plasma_temp[i] << std::endl;
    }
    for (index_t i = 0; i < background.plasma_dens.size(); i++) {
        std::cout << background.plasma_vel_x[i] << std::endl;
    }
    for (index_t i = 0; i < background.plasma_dens.size(); i++) {
        std::cout << background.plasma_vel_y[i] << std::endl;
    }

    VoxelMesh& mesh_voxel = dynamic_cast<VoxelMesh&>(mesh);
    std::cout << "Mesh:" << std::endl;
    std::cout << mesh_voxel.Lx << std::endl;
    std::cout << mesh_voxel.Ly << std::endl;
    std::cout << mesh_voxel.nx << std::endl;
    std::cout << mesh_voxel.ny << std::endl;


    for (index_t i_source = 0; i_source < sources.size(); i_source++) {
        std::cout << "Source:" << std::endl;
        LineSource& line = dynamic_cast<LineSource&>(*sources[i_source]);
        std::cout << line.pos_a.x << std::endl;
        std::cout << line.pos_a.y << std::endl;
        std::cout << line.pos_b.x << std::endl;
        std::cout << line.pos_b.y << std::endl;
        std::cout << line.angle_min << std::endl;
        std::cout << line.angle_max << std::endl;
    }

    std::uniform_real_distribution<real_t> r_unif(0, 1);
    std::normal_distribution<real_t> r_normal(0, 1);

    // calculate the total number of particles (important for normalization)
    index_t n_particles_total = 0;
    for (index_t i_source = 0; i_source < sources.size(); i_source++) {
        n_particles_total += sources[i_source]->max_particles();
    }

    // initialize the density array to zero
    for (index_t i = 0; i < mesh.n_cells(); i++) {
        density[i] = 0;
    }

    // History based simulation (because its easy)
    for (int i_source = 0; i_source < sources.size(); i_source++) {
        while (sources[i_source]->particles_left() > 0) {
            // Generate particle
            Particle particle = sources[i_source]->generate_particle();
            particle.geometry.cell_id =
                mesh.get_cell_index(particle.geometry.position);

            for (index_t i_step = 0; i_step < params.max_steps; i_step++) {
                real_t weight_old = particle.weight;
                real_t t_old = particle.t;
                index_t cell_id_old = particle.geometry.cell_id;
                // Simulate the particle up to the next event
                update_particle(mesh, background, particle);

                // Tally
                real_t sfn = 0.5 * background.plasma_dens[0] * 0.1 * 1.0e4 /
                             n_particles_total;
                real_t v = mesh.get_volume(cell_id_old);
                density[cell_id_old] += (particle.t - t_old) * weight_old / v * sfn;

                if (particle.weight < 1.0e-8) {
                    break;
                }
            }
        }
    }
}

}  // namespace eiron2