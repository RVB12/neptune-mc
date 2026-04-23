#ifndef NEPTUNE_UPDATE_PARTICLE_HPP
#define NEPTUNE_UPDATE_PARTICLE_HPP

#include <cassert>

#include "neptune/common/utils.hpp"
#include "neptune/geometry/mesh.hpp"
#include "neptune/geometry/vec2.hpp"
#include "neptune/physics/event.hpp"
#include "neptune/physics/particle.hpp"
#include "neptune/simulator/simulator.hpp"

namespace neptune {

/**
 * @brief Simulate a particle until an event occurs.
 * An event can be a collision with a cell boundary, a charge exchange collision
 * or that the particle has reached the final time.
 *
 * @param p
 */
EventType update_particle(
    const Mesh& mesh,
    const Background& background,
    const vector<std::unique_ptr<ChargeExchangeCollision>>& cx_collisions,
    const vector<std::unique_ptr<IonisationCollision>>& ion_collisions,
    const vector<std::unique_ptr<BoundaryCondition>>& boundary_conditions,
    Particle& p,
    real_t t_final = -1);

}  // namespace neptune

#endif  // NEPTUNE_UPDATE_PARTICLE_HPP