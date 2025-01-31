#ifndef EIRON2_SIMULATOR_KDMC
#define EIRON2_SIMULATOR_KDMC

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
#include "eiron2/simulator/execution_statistics.hpp"
#include "eiron2/simulator/simulator.hpp"

namespace eiron2 {

/**
 * @brief The main function of Eiron2,
 * Simulate a number of particles on a given grid.
 *
 */
void simulate_kdmc(Mesh& mesh,
              PlasmaBackground& background,
              vector<std::unique_ptr<Source>>& sources,
              vector<real_t>& density,
              vector<BoundaryCondition>& boundary_conditions,
              SimulationParams params,
              ExecutionStatistics& stats,
              real_t dt);

}  // namespace eiron2

#endif