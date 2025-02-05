#ifndef EIGRENE_SIMULATOR_EXECUTION_STATISTICS_HPP
#define EIGRENE_SIMULATOR_EXECUTION_STATISTICS_HPP

#include "eigrene/common/utils.hpp"

namespace eigrene {

struct ExecutionStatistics {
    // The number of cell crossings
    index_t num_crossings = 0;
    // The number of charge exchange collisions
    index_t num_collisions = 0;
    // The time spent building the spatial index and neighbor list (in seconds)
    real_t grid_building_time = 0;
    // The time spent tracing the particles (in seconds)
    real_t trace_time = 0;
    // The total runtime of the simulation (in seconds)
    real_t total_time = 0;
};

}  // namespace eigrene

#endif