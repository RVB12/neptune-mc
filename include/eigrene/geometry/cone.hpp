#ifndef EIGRENE_GEOMETRY_CONE_HPP
#define EIGRENE_GEOMETRY_CONE_HPP

#include "eigrene/geometry/vec2.hpp"
#include "eigrene/geometry/vec3.hpp"

namespace eigrene {

/**
 * @brief Find the intersection of a ray with position pos and velocity vel with
 * the cone defined by the line segment a-b rotated around the z-axis
 *
 * @param pos The position of the ray
 * @param vel The direction of the ray
 * @param r_a The radius of the cone at point a
 * @param z_a The z-coordinate of point a
 * @param r_b The radius of the cone at point b
 * @param z_b The z-coordinate of point b
 *
 * @return The time untill the intersection at point pos + t * vel
 */
real_t ray_cone_intersection(
    Vec3 pos, Vec3 vel, real_t r_a, real_t z_a, real_t r_b, real_t z_b);
}  // namespace eigrene

#endif