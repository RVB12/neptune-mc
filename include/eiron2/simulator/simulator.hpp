#ifndef EIRON2_SIMULATOR
#define EIRON2_SIMULATOR

#include <cstdint>
#include <random>

#include "eiron2/common/vector.hpp"
#include "eiron2/geometry/cellcrossing.hpp"
#include "eiron2/geometry/mesh.hpp"
#include "eiron2/physics/particle.hpp"
#include "eiron2/common/utils.hpp"

namespace eiron2 {

class PlasmaBackground {
   public:
    PlasmaBackground(vectorview<real_t> plasma_dens_,
                     vectorview<real_t> plasma_temp_,
                     vectorview<real_t> plasma_vel_x_,
                     vectorview<real_t> plasma_vel_y_)
        : plasma_dens(plasma_dens_),
          plasma_temp(plasma_temp_),
          plasma_vel_x(plasma_vel_x_),
          plasma_vel_y(plasma_vel_y_)
    {}

    vectorview<real_t> plasma_dens;
    vectorview<real_t> plasma_temp;
    vectorview<real_t> plasma_vel_x;
    vectorview<real_t> plasma_vel_y;
};

struct SimulationParams {
    int64_t seed = 0;
    index_t n_particles = 1000;
    index_t max_steps = 100000;
};

/**
 * @brief The main function of Eiron2,
 * Simulate a number of particles on a given grid.
 *
 */
void simulate(Mesh& mesh,
              PlasmaBackground& background,
              vectorview<real_t> density,
              SimulationParams params);

}  // namespace eiron2

#endif