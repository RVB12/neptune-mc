#include "eiron2/simulator/simulator.hpp"
namespace eiron2 {

void simulate(Mesh& mesh,
              PlasmaBackground& background,
              vector<Source*>& sources,
              vector<real_t>& density,
              SimulationParams params)
{
    real_t T_i = 40;
    real_t T_e = T_i;
    real_t T_n = T_i;
    real_t K_cx = 3.2e-15 * std::sqrt(T_i / 0.026);
    real_t K_i = 2.0e-13 / (6.0 + (T_e / 13.6)) * std::sqrt(T_e / 13.6) *
                 std::exp(-13.6 / T_e);

        std::uniform_real_distribution<real_t> r_unif(0, 1);
        std::normal_distribution<real_t> r_normal(0, 1);

    index_t n_particles_total = 0;
    for (index_t i_source = 0; i_source < sources.size(); i_source++) {
        n_particles_total += sources[i_source]->max_particles();
    }

    // History based simulation (because its easy)
    for (int i_source = 0; i_source < sources.size(); i_source++) {
        while (sources[i_source]->particles_left() > 0) {
            Particle particle = sources[i_source]->generate_particle();
            particle.geometry.cell_id = mesh.get_cell_index(particle.geometry.position);

            for (index_t i_step = 0; i_step < params.max_steps; i_step++) {
                index_t cell_id =
                    particle.geometry
                        .cell_id;  // Current cell id of the particle
                real_t weight_old = particle.weight;

                // random number to determine time till collision
                real_t eps = r_unif(particle.rng);
                // time till collision
                real_t tc =
                    -std::log(eps) / (K_cx * background.plasma_dens[cell_id]);

                // Calculate the next cell crossing
                CellCrossing crossing =
                    mesh.get_next_crossing(particle.geometry);

                // Check which one is earlier, charge exchange or cell crossing
                real_t dt;
                if (tc < crossing.t) {
                    // Move the particle until the collision
                    particle.geometry.position +=
                        tc * particle.geometry.velocity;
                    // Generate new velocity from plasma background (BGK
                    // operator)
                    real_t stddev = std::sqrt(background.plasma_temp[cell_id]);
                    std::normal_distribution<real_t> dist_vx(
                        background.plasma_vel_x[cell_id], stddev);
                    std::normal_distribution<real_t> dist_vy(
                        background.plasma_vel_y[cell_id], stddev);
                    particle.geometry.velocity =
                        Vec2(dist_vx(particle.rng), dist_vy(particle.rng));

                    dt = tc;
                }
                else {
                    dt = crossing.t;
                    if (crossing.next_cell_id != -1) {
                        // Regular grid crossing, move the particle to the next
                        // cell
                        particle.geometry.position = crossing.position;
                        particle.geometry.cell_id = crossing.next_cell_id;
                    }
                    else {
                        // Boundary crossing
                        // For now, we just reflect particles at boundaries
                        // 0,1,2,3,4 and absorb particles at boundary 3 (x ==
                        // Xl)

                        // Todo: Implement proper boundary conditions
                        // maybe implement surface norm in the mesh

                        particle.geometry.position = crossing.position;
                        if (crossing.boundary_id == 0) {
                            particle.geometry.velocity.x =
                                -particle.geometry.velocity.x;
                        }
                        else if (crossing.boundary_id == 1 ||
                                 crossing.boundary_id == 3) {
                            particle.geometry.velocity.y =
                                -particle.geometry.velocity.y;
                        }
                        else if (crossing.boundary_id == 2) {
                            particle.weight = 0;
                        }
                        else {
                            assert(false);
                        }
                    }
                }

                // Ionisation
                particle.weight =
                    particle.weight *
                    std::exp(-dt * K_i * background.plasma_dens[cell_id]);

                particle.t += dt;

                real_t sfn = 0.5 * background.plasma_dens[0] * 0.1 * 1.0e4 /
                            n_particles_total;
                real_t v = mesh.get_volume(cell_id);
                density[cell_id] += dt * weight_old / v * sfn;
                // density[cell_id] += dt * (weight_old + particle.weight) / 2 /
                // v * sfn;

                if (particle.weight < 0.000000001) {
                    break;
                }
            }
        }
    }

}

}  // namespace eiron2