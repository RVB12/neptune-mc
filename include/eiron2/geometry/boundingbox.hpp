#ifndef EIRON2_GEOMETRY_BOUNDINGBOX_HPP
#define EIRON2_GEOMETRY_BOUNDINGBOX_HPP

#include <cassert>

#include "eiron2/common/utils.hpp"
#include "eiron2/common/vec2.hpp"
#include "eiron2/common/vec3.hpp"

namespace eiron2 {

/**
 * @brief A simple axis-aligned bounding box
 */
class BoundingBox
{

public:
    BoundingBox() : min(Vec2(0, 0)), max(Vec2(0, 0)) {}

    BoundingBox(Vec2 min_, Vec2 max_) : min(min_), max(max_) {}

    inline BoundingBox merge(BoundingBox other)
    {
        return BoundingBox(
            Vec2(std::min(min.x, other.min.x), std::min(min.y, other.min.y)),
            Vec2(std::max(max.x, other.max.x), std::max(max.y, other.max.y)));
    }

    Vec2 min;
    Vec2 max;
};

}  // namespace eiron2

#endif