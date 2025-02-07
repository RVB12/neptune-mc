#ifndef EIGRENE_SIMULATOR
#define EIGRENE_SIMULATOR

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
#include "eigrene/physics/plasmabackground.hpp"

namespace eigrene {


struct SimulationParams {
    int64_t seed = 0;
    index_t max_steps = 100000;
    real_t t_final = 1.0;
};

/**
 * @brief The main function of Eigrene,
 * Simulate a number of particles on a given grid.
 *
 */
void simulate(Mesh& mesh,
              PlasmaBackground& background,
              vector<std::unique_ptr<Source>>& sources,
              vector<real_t>& density,
              vector<BoundaryCondition>& boundary_conditions,
              SimulationParams params,
              ExecutionStatistics& stats);

void simulate_fixed_time(Mesh& mesh,
              PlasmaBackground& background,
              vector<std::unique_ptr<Source>>& sources,
              vector<real_t>& density,
              vector<BoundaryCondition>& boundary_conditions,
              SimulationParams params,
              ExecutionStatistics& stats,
              real_t t_final);

void simulate_kdmc(Mesh& mesh,
              PlasmaBackground& background,
              vector<std::unique_ptr<Source>>& sources,
              vector<real_t>& density,
              vector<BoundaryCondition>& boundary_conditions,
              SimulationParams params,
              ExecutionStatistics& stats,
              real_t t_final,
              real_t dt);


}  // namespace eigrene

#endif