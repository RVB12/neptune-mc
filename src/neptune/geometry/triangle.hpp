#ifndef NEPTUNE_GEOMETRY_TRIANGLE_HPP
#define NEPTUNE_GEOMETRY_TRIANGLE_HPP

#include "neptune/common/utils.hpp"
#include "neptune/geometry/boundingbox.hpp"
#include "neptune/geometry/vec2.hpp"

namespace neptune {

/**
 * @brief A triangle, defined by three indices into the vertex array
 *        The indices should be listed in counter clockwise order.
 */
class TriangleIndices {
   public:
    index_t ind1;
    index_t ind2;
    index_t ind3;
};

/**
 * @brief A triangle, defined by three vertices
 *        The vertices should be listed in counter clockwise order.
 */
class Triangle {
   public:
    Triangle(const Vec2& p1, const Vec2& p2, const Vec2& p3)
        : p1(p1), p2(p2), p3(p3)
    {}

    bool contains(const Vec2& pos) const;

    BoundingBox bounding_box() const;

    const Vec2& p1;
    const Vec2& p2;
    const Vec2& p3;
};

}  // namespace neptune

#endif