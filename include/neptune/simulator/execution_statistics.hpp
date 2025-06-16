#ifndef NEPTUNE_SIMULATOR_EXECUTION_STATISTICS_HPP
#define NEPTUNE_SIMULATOR_EXECUTION_STATISTICS_HPP

#include "neptune/common/utils.hpp"

namespace neptune {

struct ExecutionStatistics {
    // The number of cell crossings
    index_t num_crossings = 0;
    // The number of charge exchange collisions
    index_t num_collisions = 0;
    // The number of particles that have been simulated
    index_t total_particles = 0;
    // The time spent building the spatial index and neighbor list (in seconds)
    real_t grid_building_time = 0;
    // The time spent tracing the particles (in seconds)
    real_t trace_time = 0;
    // The total runtime of the simulation (in seconds)
    real_t total_time = 0;
    // The total number of kd iterations
    index_t num_kd_steps = 0;
    // The number of times the kd simulation switched to kinetic
    index_t num_kd_to_kinetic_switches = 0;

};

}  // namespace neptune

#endif