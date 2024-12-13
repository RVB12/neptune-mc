#ifndef EIRON2_SIMULATOR
#define EIRON2_SIMULATOR

#include <cstdint>
#include <random>
#include <vector>

#include "eiron2/geometry/cellcrossing.hpp"
#include "eiron2/geometry/mesh.hpp"
#include "eiron2/physics/particle.hpp"

namespace eiron2 {

struct PlasmaBackground {
    std::vector<double>& plasma_dens;
    std::vector<double>& plasma_temp;
    std::vector<double>& plasma_vel_x;
    std::vector<double>& plasma_vel_y;
    std::vector<double>& plasma_vel_z;
};

struct SimulationParams {
    int64_t seed = 0;
    int n_particles = -1;
    int max_steps = 1000;
};

/**
 * @brief The main function of Eiron2,
 * Simulate a number of particles on a given grid.
 *
 */
void simulate(Mesh& mesh,
              PlasmaBackground& background,
              std::vector<double>& density,
              SimulationParams params)
{
    double Kcx = 0.1;
    double Ki = 1;

    std::uniform_real_distribution<double> r_unif(0, 1);
    std::normal_distribution<double> r_normal(0, 1);

    // History based simulation (because its easy)
    for (int i_particle; i_particle < params.n_particles; i_particle++) {
        // By basing the generator seed on the particle id, we can get
        // consistent results even if we parallelize.
        std::mt19937 gen{i_particle};

        // Generate a new particle
        Particle particle;

        // Currently, these are just some ad hoc initial conditions
        // I should look into simulating a more realistic line source
        particle.geometry.position.x = 0;
        particle.geometry.position.z = r_unif(gen);
        particle.geometry.position.y = r_unif(gen);

        particle.geometry.cell_id = mesh.get_cell_index(particle.geometry.position);

        particle.geometry.velocity.x = 0.5 + 0.25 * r_unif(gen);
        particle.geometry.velocity.y = 0.25 * r_normal(gen);
        particle.geometry.velocity.z = 0.25 * r_normal(gen);

        // Full simulation loop
        for (int i_step; i_step < params.max_steps; i_step++) {
            int cell_id =
                particle.geometry.cell_id;  // Current cell id of the particle

            // random number to determine time till collision
            double eps = r_unif(gen);
            // time till collision
            double tc =
                -std::log(eps) / (Kcx * background.plasma_dens[cell_id]);

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
                std::normal_distribution<double> dist_vz(
                    background.plasma_vel_z[cell_id], stddev);
                particle.geometry.velocity =
                    Vec3(dist_vx(gen), dist_vy(gen), dist_vz(gen));

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
                    particle.geometry.position = crossing.position;
                    if (crossing.boundary_id == 0) {
                        particle.geometry.velocity.x =
                            -particle.geometry.velocity.x;
                    }
                    else if (crossing.boundary_id == 1 ||
                             crossing.boundary_id == 4) {
                        particle.geometry.velocity.y =
                            -particle.geometry.velocity.y;
                    }
                    else if (crossing.boundary_id == 2 ||
                             crossing.boundary_id == 5) {
                        particle.geometry.velocity.z =
                            -particle.geometry.velocity.z;
                    }
                    else if (crossing.boundary_id == 3) {
                        particle.weight = 0;
                    }
                    else {
                        assert(false);
                    }
                }
            }

            // Ionisation
            particle.weight = particle.weight * std::exp(-dt*Ki*background.plasma_dens[cell_id]);

            particle.t += dt;

            density[cell_id] += dt * particle.weight;

            if (particle.weight < 0.0001) break;
        }
    }
}

}  // namespace eiron2

#endif