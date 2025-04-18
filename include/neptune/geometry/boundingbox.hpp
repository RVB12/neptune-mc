#ifndef NEPTUNE_GEOMETRY_BOUNDINGBOX_HPP
#define NEPTUNE_GEOMETRY_BOUNDINGBOX_HPP

#include <cassert>

#include "neptune/common/utils.hpp"
#include "neptune/geometry/vec2.hpp"
#include "neptune/geometry/vec3.hpp"

namespace neptune {

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

}  // namespace neptune

#endif