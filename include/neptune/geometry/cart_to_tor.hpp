#ifndef CART_TO_TOR_HPP
#define CART_TO_TOR_HPP

#include "neptune/geometry/vec3.hpp"

namespace neptune {

/**
 * @brief Convert a vector toroidal coordinates to cartesian coordinates
 *
 * This is mostly used for directions, where e.g. the mean velocity
 * is represented in toroidal coordinates, but the particle geometry is
 * represented in cartesian coordinates.
 *
 * @param position The position of the particle, used to determine the angle for
 * the transformation
 * @param direction The direction in toroidal coordinates
 *                  x: the radial coordinate
 *                  y: the toroidal coordinate
 *                  z: the vertical coordinate
 *
 * @return Vec3 The direction in cartesian coordinates
 */
inline Vec3 transform_tor_dir_to_cart(Vec3 position, Vec3 direction)
{
    // Get cos and sine of the angle of the position
    real_t r = std::sqrt(position.x * position.x + position.y * position.y);
    real_t c = position.x / r;
    real_t s = position.y / r;

    // Rotate the direction vector
    real_t x_new = direction.x * c + direction.y * s;
    real_t y_new = -direction.x * s + direction.y * c;
    real_t z_new = direction.z;

    return Vec3(x_new, y_new, z_new);
}

/**
 * @brief Convert a vector in cartesian coordinates to toroidal coordinates
 *
 * This is mostly used for directions, where e.g. the mean velocity
 * is represented in toroidal coordinates, but the particle geometry is
 * represented in cartesian coordinates.
 *
 * @param position The position of the particle, used to determine the angle for
 * the transformation
 * @param direction The direction in toroidal coordinates
 *                  x: the radial coordinate
 *                  y: the toroidal coordinate
 *                  z: the vertical coordinate
 *
 * @return Vec3 The direction in cartesian coordinates
 */
inline Vec3 transform_cart_to_tor(Vec3 position, Vec3 direction)
{
    // Get cos and sine of the angle of the position
    real_t r = std::sqrt(position.x * position.x + position.y * position.y);
    real_t c = position.x / r;
    real_t s = position.y / r;

    // Inverse rotation of the direction vector
    real_t x_new = direction.x * c - direction.y * s;
    real_t y_new = direction.x * s + direction.y * c;
    real_t z_new = direction.z;

    return Vec3(x_new, y_new, z_new);
}

}  // namespace neptune

#endif