#ifndef EIRON2_SIMULATOR
#define EIRON2_SIMULATOR

#include <cstdint>
#include <random>
#include <vector>

#include "eiron2/common/vector.hpp"
#include "eiron2/geometry/cellcrossing.hpp"
#include "eiron2/geometry/mesh.hpp"
#include "eiron2/physics/particle.hpp"
#include "eiron2/physics/source.hpp"
#include "eiron2/common/utils.hpp"
#include "eiron2/physics/boundarycondition.hpp"

namespace eiron2 {

class PlasmaBackground {
   public:
    PlasmaBackground() = default;

    PlasmaBackground(index_t n_cells)
        : plasma_dens(n_cells),
          plasma_temp(n_cells),
          plasma_vel_x(n_cells),
          plasma_vel_y(n_cells)
    {}

    vector<real_t> plasma_dens;
    vector<real_t> plasma_temp;
    vector<real_t> plasma_vel_x;
    vector<real_t> plasma_vel_y;
};

struct SimulationParams {
    int64_t seed = 0;
    index_t max_steps = 100000;
};

/**
 * @brief The main function of Eiron2,
 * Simulate a number of particles on a given grid.
 *
 */
void simulate(Mesh& mesh,
              PlasmaBackground& background,
              vector<std::unique_ptr<Source>>& sources,
              vector<real_t>& density,
              vector<BoundaryCondition>& boundary_conditions,
              SimulationParams params);

}  // namespace eiron2

#endif