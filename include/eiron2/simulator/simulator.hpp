#ifndef EIRON2_SIMULATOR
#define EIRON2_SIMULATOR

#include <cstdint>
#include <random>

#include "eiron2/common/vector.hpp"
#include "eiron2/geometry/cellcrossing.hpp"
#include "eiron2/geometry/mesh.hpp"
#include "eiron2/physics/particle.hpp"

namespace eiron2 {

class PlasmaBackground {
   public:
    PlasmaBackground(vectorview<double> plasma_dens_,
                     vectorview<double> plasma_temp_,
                     vectorview<double> plasma_vel_x_,
                     vectorview<double> plasma_vel_y_)
        : plasma_dens(plasma_dens_),
          plasma_temp(plasma_temp_),
          plasma_vel_x(plasma_vel_x_),
          plasma_vel_y(plasma_vel_y_)
    {}

    vectorview<double> plasma_dens;
    vectorview<double> plasma_temp;
    vectorview<double> plasma_vel_x;
    vectorview<double> plasma_vel_y;
};

struct SimulationParams {
    int64_t seed = 0;
    int n_particles = 1000;
    int max_steps = 100000;
};

/**
 * @brief The main function of Eiron2,
 * Simulate a number of particles on a given grid.
 *
 */
void simulate(Mesh& mesh,
              PlasmaBackground& background,
              vectorview<double> density,
              SimulationParams params);

}  // namespace eiron2

#endif