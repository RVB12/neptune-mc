#ifndef NEPTUNE_SIMULATIONPARAMS_HPP
#define NEPTUNE_SIMULATIONPARAMS_HPP

#include "neptune/common/utils.hpp"

namespace neptune {

enum class AlgorithmType {
    MC,  // Standard Monte Carlo
    KDMC,  // Kinetic-Diffusion Monte Carlo with new boundary conditions
    KDMC_OLD,  // Kinetic-Diffusion Monte Carlo with old boundary conditions
    ARCMC  // Artifially Reduced Collisionality Monte Carlo
};

struct SimulationParams {
    // This value is used to seed the rng
    uint64_t seed = 1302;
    // The maximum number of simulation steps a particle can take
    index_t max_steps = 100000;
    // The maximum time a particle can take
    // if it is negative, no
    real_t t_final = -1.0;
    // The number of batches to run (for variance estimation)
    index_t n_batches = 1;
    // The algorithm to use
    AlgorithmType algorithm = AlgorithmType::MC;
    // Time step to use for KDMC
    real_t kdmc_dt;
    // Parameter for ARCMC
    index_t armc_M;
};
}  // namespace neptune

#endif