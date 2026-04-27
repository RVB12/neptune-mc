#ifndef NEPTUNE_SIMULATOR
#define NEPTUNE_SIMULATOR

#include <cstdint>
#include <random>
#include <vector>

#include "neptune/collisions/abstract_collisions.hpp"
#include "neptune/common/utils.hpp"
#include "neptune/common/vector.hpp"
#include "neptune/estimation/estimators.hpp"
#include "neptune/geometry/cellcrossing.hpp"
#include "neptune/geometry/mesh.hpp"
#include "neptune/physics/boundarycondition.hpp"
#include "neptune/physics/collision_processes.hpp"
#include "neptune/physics/particle.hpp"
#include "neptune/physics/plasmabackground.hpp"
#include "neptune/physics/source.hpp"
#include "neptune/simulator/execution_statistics.hpp"
#include "neptune/physics/background.hpp"
#include "neptune/estimation/estimator.hpp"
#include "neptune/simulator/simulationparams.hpp"

namespace neptune {

/**
 * @brief The main function of Neptune,
 * Simulate a number of particles on a given grid.
 *
 */
void simulate(
    const Mesh& mesh,
    const Background& background,
    const vector<std::unique_ptr<ChargeExchangeCollision>>& cx_collisions,
    const vector<std::unique_ptr<IonisationCollision>>& ion_collisions,
    vector<std::unique_ptr<Source>>& sources,
    EstimatorGroup& estimator_group,
    vector<std::unique_ptr<BoundaryCondition>>& boundary_conditions,
    SimulationParams params,
    ExecutionStatistics& stats);

void simulate_kdmc(
    Mesh& mesh,
    const vector<std::unique_ptr<BGKScatteringProcess>>& bgk_processes,
    const vector<std::unique_ptr<IonisationProcess>>& ion_processes,
    vector<std::unique_ptr<Source>>& sources,
    Estimator_old& estimator,
    FluidEstimator& fluid_estimator,
    vector<std::unique_ptr<BoundaryCondition>>& boundary_conditions,
    SimulationParams params,
    ExecutionStatistics& stats,
    real_t dt,
    bool use_new_version = true);

void simulate_armc(
    Mesh& mesh,
    const vector<std::unique_ptr<BGKScatteringProcess>>& bgk_processes,
    const vector<std::unique_ptr<IonisationProcess>>& ion_processes,
    vector<std::unique_ptr<Source>>& sources,
    ArcmcEstimator& estimator,
    vector<std::unique_ptr<BoundaryCondition>>& boundary_conditions,
    SimulationParams params,
    ExecutionStatistics& stats,
    index_t M);

}  // namespace neptune

#endif