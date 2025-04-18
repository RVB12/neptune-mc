#ifndef EIGRENE_UPDATE_PARTICLE_HPP
#define EIGRENE_UPDATE_PARTICLE_HPP

#include <cassert>

#include "eigrene/common/utils.hpp"
#include "eigrene/geometry/vec2.hpp"
#include "eigrene/geometry/mesh.hpp"
#include "eigrene/physics/particle.hpp"
#include "eigrene/simulator/simulator.hpp"

namespace eigrene {

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

}  // namespace eigrene

#endif  // EIGRENE_UPDATE_PARTICLE_HPP