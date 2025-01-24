#include "eiron2/geometry/triangle.hpp"

#include <algorithm>

namespace eiron2 {

real_t triangle_sign_helper(const Vec2& p1, const Vec2& p2, const Vec2& p3)
{
    return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y);
}

// Source:
// https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle
bool Triangle::contains(const Vec2& pos) const
{
    // TODO: update this so that it has a tolerance

    // Half-plane checks
    real_t d1 = triangle_sign_helper(pos, p1, p2);
    real_t d2 = triangle_sign_helper(pos, p2, p3);
    real_t d3 = triangle_sign_helper(pos, p3, p1);

    // Check if the point is on the same side of all triangle edges
    bool has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
    bool has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);

    return !(has_neg && has_pos);
}

BoundingBox Triangle::bounding_box() const
{
    real_t min_x = std::min(std::min(p1.x, p2.x), p3.x);
    real_t min_y = std::min(std::min(p1.y, p2.y), p3.y);

    real_t max_x = std::max(std::max(p1.x, p2.x), p3.x);
    real_t max_y = std::max(std::max(p1.y, p2.y), p3.y);

    return BoundingBox(Vec2(min_x, min_y), Vec2(max_x, max_y));
}

}  // namespace eiron2
