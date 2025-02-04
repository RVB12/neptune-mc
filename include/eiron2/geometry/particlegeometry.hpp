#ifndef EIRON2_PARTICLEGEOMETRY_HPP
#define EIRON2_PARTICLEGEOMETRY_HPP

#include <cassert>

#include "eiron2/geometry/vec2.hpp"
#include "eiron2/common/utils.hpp"

namespace eiron2 {

/**
 * @brief Class containing geometric information about a particle,
 * this consists of the particle's position, its velocity and
 * the id of the cell the particle is in.
 */
class ParticleGeometry {
   public:
   ParticleGeometry() = default;
    // Contructors
    ParticleGeometry(Vec2 position_, Vec2 velocity_, index_t cell_id_)
        : position(position_), velocity(velocity_), cell_id(cell_id_)
    {}

    // Member variables
    Vec2 position;     // Current position
    Vec2 velocity;     // Current velocity
    index_t cell_id = -1;  // id of the cell the particle is in
    index_t boundary_id = -1; // id of the boundary the particle is on
                              // -1 if the particle is not on a boundary
};

}  // namespace eiron2

#endif