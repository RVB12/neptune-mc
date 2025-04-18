#ifndef EIGRENE_SIMULATOR
#define EIGRENE_SIMULATOR

#include <cstdint>
#include <random>
#include <vector>

#include "eigrene/common/utils.hpp"
#include "eigrene/common/vector.hpp"
#include "eigrene/geometry/cellcrossing.hpp"
#include "eigrene/geometry/mesh.hpp"
#include "eigrene/physics/boundarycondition.hpp"
#include "eigrene/physics/particle.hpp"
#include "eigrene/physics/plasmabackground.hpp"
#include "eigrene/physics/source.hpp"
#include "eigrene/simulator/execution_statistics.hpp"
#include "eigrene/estimation/estimators.hpp"

namespace eigrene {

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
 * @brief The main function of Eigrene,
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

}  // namespace eigrene

#endif