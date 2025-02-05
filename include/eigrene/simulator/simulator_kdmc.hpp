#ifndef EIGRENE_SIMULATOR_KDMC
#define EIGRENE_SIMULATOR_KDMC

#include <cstdint>
#include <random>
#include <vector>

#include "eigrene/common/vector.hpp"
#include "eigrene/geometry/cellcrossing.hpp"
#include "eigrene/geometry/mesh.hpp"
#include "eigrene/physics/particle.hpp"
#include "eigrene/physics/source.hpp"
#include "eigrene/common/utils.hpp"
#include "eigrene/physics/boundarycondition.hpp"
#include "eigrene/simulator/execution_statistics.hpp"
#include "eigrene/simulator/simulator.hpp"

namespace eigrene {

/**
 * @brief The main function of Eigrene,
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

}  // namespace eigrene

#endif