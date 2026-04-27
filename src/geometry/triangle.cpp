#include "neptune/geometry/triangle.hpp"

#include <algorithm>

namespace neptune {

/**
 * Determine the side of a line that a point is on.
 *
 *
 * @param p1 Beginning of the line
 * @param p2 End of the line
 * @param p3 The point to check
 *
 * @return real_t If positive, the point is on the left side of the line, if
 * negative, the point is on the right side of the line, if zero, the point is
 * on the line (but not necessarily between p1 and p2)
 *
 * @note How it works:
 * 1. Determine a vector orthogonal to the line (p1, p2), this vector is: (
 * -(p2.y-p1.y) , p2.x-p1.x )
 * 2. Calculate the dot product of point - p1 and the orthogonal vector
 */
real_t triangle_sign_helper(const Vec2& p1, const Vec2& p2, const Vec2& point)
{
    // Vec2 orth(  -(p2.y - p1.y), (p2.x - p1.x) );
    // Vec2 diff = point - p1;
    // return orth.dot(diff);
    return (point.x - p1.x) * (p1.y - p2.y) + (point.y - p1.y) * (p2.x - p1.x);
}

// Source:
// http://totologic.blogspot.com/2014/01/accurate-point-in-triangle-test.html
bool Triangle::contains(const Vec2& pos) const
{
    // center of the triangle
    Vec2 center = (p1 + p2 + p3) / 3.0;

    // Make the triangle slightly larger to account for floating point errors
    real_t eps = 1.0e-10;
    Vec2 p1_tol = p1 + eps * (p1 - center);
    Vec2 p2_tol = p2 + eps * (p2 - center);
    Vec2 p3_tol = p3 + eps * (p3 - center);

    // Check which side of the lines the point is on
    real_t d1 = triangle_sign_helper(p1_tol, p2_tol, pos);
    real_t d2 = triangle_sign_helper(p2_tol, p3_tol, pos);
    real_t d3 = triangle_sign_helper(p3_tol, p1_tol, pos);

    // If all three sides are positive (i.e. it is always on the left side)
    // the point is inside the triangle, otherwise it is outside
    return (d1 >= 0) && (d2 >= 0) && (d3 >= 0);
}

BoundingBox Triangle::bounding_box() const
{
    real_t min_x = std::min(std::min(p1.x, p2.x), p3.x);
    real_t min_y = std::min(std::min(p1.y, p2.y), p3.y);

    real_t max_x = std::max(std::max(p1.x, p2.x), p3.x);
    real_t max_y = std::max(std::max(p1.y, p2.y), p3.y);

    return BoundingBox(Vec2(min_x, min_y), Vec2(max_x, max_y));
}

}  // namespace neptune
