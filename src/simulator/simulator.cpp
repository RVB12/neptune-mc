#include "eiron2/simulator/simulator.hpp"
namespace eiron2 {

void simulate(Mesh& mesh,
              PlasmaBackground& background,
              vectorview<double> density,
              SimulationParams params)
{
    double T_i = 40;
    double T_e = T_i;
    double T_n = T_i;
    double K_cx = 3.2e-15 * std::sqrt(T_i / 0.026);
    double K_i = 2.0e-13 / (6.0 + (T_e / 13.6)) *
                           std::sqrt(T_e / 13.6) * std::exp(-13.6 / T_e);

    // History based simulation (because its easy)
    #pragma omp parallel for num_threads(4)
    for (int i_particle = 0; i_particle < params.n_particles; i_particle++) {
        std::uniform_real_distribution<double> r_unif(0, 1);
        std::normal_distribution<double> r_normal(0, 1);

        // By basing the generator seed on the particle id, we can get
        // consistent results even if we parallelize.
        std::mt19937 gen(i_particle);

        // Generate a new particle
        Particle particle;

        particle.weight = 1.0;

        // Currently, these are just some ad hoc initial conditions
        // I should look into simulating a more realistic line source
        // particle.geometry.position.x = 0;
        // particle.geometry.position.y = r_unif(gen);
        particle.geometry.position.x = 0.5;
        particle.geometry.position.y = 0.5;

        particle.geometry.cell_id =
            mesh.get_cell_index(particle.geometry.position);

        particle.geometry.velocity.x = 0.5 + 0.25 * r_unif(gen);
        particle.geometry.velocity.y = 0.25 * r_normal(gen);

        // Full simulation loop
        for (int i_step = 0; i_step < params.max_steps; i_step++) {
            int cell_id =
                particle.geometry.cell_id;  // Current cell id of the particle
            double weight_old = particle.weight;

            // random number to determine time till collision
            double eps = r_unif(gen);
            // time till collision
            double tc =
                -std::log(eps) / (K_cx * background.plasma_dens[cell_id]);

            // Calculate the next cell crossing
            CellCrossing crossing = mesh.get_next_crossing(particle.geometry);

            // Check which one is earlier, charge exchange or cell crossing
            double dt;
            if (tc < crossing.t) {
                // Move the particle until the collision
                particle.geometry.position += tc * particle.geometry.velocity;
                // Generate new velocity from plasma background (BGK operator)
                double stddev = std::sqrt(background.plasma_temp[cell_id]);
                std::normal_distribution<double> dist_vx(
                    background.plasma_vel_x[cell_id], stddev);
                std::normal_distribution<double> dist_vy(
                    background.plasma_vel_y[cell_id], stddev);
                particle.geometry.velocity =
                    Vec2(dist_vx(gen), dist_vy(gen));

                dt = tc;
            }
            else {
                dt = crossing.t;
                if (crossing.next_cell_id != -1) {
                    // Regular grid crossing, move the particle to the next cell
                    particle.geometry.position = crossing.position;
                    particle.geometry.cell_id = crossing.next_cell_id;
                }
                else {
                    // Boundary crossing
                    // For now, we just reflect particles at boundaries
                    // 0,1,2,3,4 and absorb particles at boundary 3 (x == Xl)

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

            double sfn = 0.5 * background.plasma_dens[0] * 0.1 * 1.0e4 / params.n_particles;
            double v = mesh.get_volume(cell_id);
            // Note, this interpolation is not in the literature, make sure
            // it doesn't cause issues.
            #pragma omp atomic
            density[cell_id] += dt * weight_old / v * sfn;
            // density[cell_id] += dt * (weight_old + particle.weight) / 2 / v * sfn;

            if (particle.weight < 0.000000001) {
                break;
            }
        }
    }
}

}  // namespace eiron2