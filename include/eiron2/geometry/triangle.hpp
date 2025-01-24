#ifndef EIRON2_GEOMETRY_TRIANGLE_HPP
#define EIRON2_GEOMETRY_TRIANGLE_HPP

#include "eiron2/common/utils.hpp"
#include "eiron2/geometry/vec2.hpp"
#include "eiron2/geometry/boundingbox.hpp"

namespace eiron2 {

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
    Triangle(const Vec2& p1, const Vec2& p2, const Vec2& p3) : p1(p1), p2(p2), p3(p3) {}

    bool contains(const Vec2& pos) const;

    BoundingBox bounding_box() const;

    const Vec2& p1;
    const Vec2& p2;
    const Vec2& p3;
};

}  // namespace eiron2

#endif