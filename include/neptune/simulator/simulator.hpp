#ifndef NEPTUNE_SIMULATOR
#define NEPTUNE_SIMULATOR

#include <cstdint>
#include <random>
#include <vector>

#include "neptune/common/utils.hpp"
#include "neptune/common/vector.hpp"
#include "neptune/geometry/cellcrossing.hpp"
#include "neptune/geometry/mesh.hpp"
#include "neptune/physics/boundarycondition.hpp"
#include "neptune/physics/particle.hpp"
#include "neptune/physics/plasmabackground.hpp"
#include "neptune/physics/source.hpp"
#include "neptune/simulator/execution_statistics.hpp"
#include "neptune/estimation/estimators.hpp"

namespace neptune {

struct SimulationParams {
    // This value is added to the particle id to determine an rng seed
    int64_t seed = 0;
    // The maximum number of simulation steps a particle can take
    index_t max_steps = 100000;
    // The maximum time a particle can take
    // if it is negative, no
    real_t t_final = -1.0;
};

/**
 * @brief The main function of Neptune,
 * Simulate a number of particles on a given grid.
 *
 */
void simulate(Mesh& mesh,
              PlasmaBackground& background,
              vector<std::unique_ptr<Source>>& sources,
              Estimator& estimator,
              vector<BoundaryCondition>& boundary_conditions,
              SimulationParams params,
              ExecutionStatistics& stats);

void simulate_kdmc(Mesh& mesh,
                   PlasmaBackground& background,
                   vector<std::unique_ptr<Source>>& sources,
                   Estimator& estimator,
                   vector<BoundaryCondition>& boundary_conditions,
                   SimulationParams params,
                   ExecutionStatistics& stats,
                   real_t dt);

void simulate_kdmc2(Mesh& mesh,
                    PlasmaBackground& background,
                    vector<std::unique_ptr<Source>>& sources,
                    Estimator& estimator,
                    vector<BoundaryCondition>& boundary_conditions,
                    SimulationParams params,
                    ExecutionStatistics& stats,
                    real_t dt);

}  // namespace neptune

#endif