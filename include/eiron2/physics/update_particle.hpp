#ifndef EIRON2_UPDATE_PARTICLE_HPP
#define EIRON2_UPDATE_PARTICLE_HPP

#include <cassert>

#include "eiron2/common/utils.hpp"
#include "eiron2/geometry/vec2.hpp"
#include "eiron2/geometry/mesh.hpp"
#include "eiron2/physics/particle.hpp"
#include "eiron2/simulator/simulator.hpp"

namespace eiron2 {

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
                     index_t& num_collisions);

}  // namespace eiron2

#endif  // EIRON2_UPDATE_PARTICLE_HPP