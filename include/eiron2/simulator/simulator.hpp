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
    int max_steps = 100000;
};

/**
 * @brief The main function of Eiron2,
 * Simulate a number of particles on a given grid.
 *
 */
void simulate(Mesh& mesh,
              PlasmaBackground& background,
              std::vector<double>& density,
              SimulationParams params);

}  // namespace eiron2

#endif