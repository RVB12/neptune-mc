#ifndef NEPTUNE_UPDATE_PARTICLE_HPP
#define NEPTUNE_UPDATE_PARTICLE_HPP

#include <cassert>

#include "neptune/common/utils.hpp"
#include "neptune/geometry/vec2.hpp"
#include "neptune/geometry/mesh.hpp"
#include "neptune/physics/particle.hpp"
#include "neptune/simulator/simulator.hpp"

namespace neptune {

/**
 * @brief Simulate a particle until an event occurs.
 * An event can be a collision with a cell boundary or a charge exchange.
 *
 * @param p
 */
void update_particle(const Mesh& mesh,
                     const PlasmaBackground& background,
                     const vector<BoundaryCondition>& boundary_conditions,
                     Particle& p,
                     index_t& num_crossings,
                     index_t& num_collisions,
                     real_t t_final = -1);

}  // namespace neptune

#endif  // NEPTUNE_UPDATE_PARTICLE_HPP