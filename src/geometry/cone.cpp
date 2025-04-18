#include "neptune/geometry/cone.hpp"

namespace neptune {

real_t ray_cone_intersection(
    Vec3 pos, Vec3 vel, real_t r_a, real_t z_a, real_t r_b, real_t z_b)
{
    // Precompute slope terms
    real_t k_r = r_b - r_a;
    real_t k_z = z_b - z_a;

    // Compute the coefficients of the quadratic equation
    real_t a =
        k_z * k_z * (vel.x * vel.x + vel.y * vel.y) - k_r * k_r * vel.z * vel.z;
    real_t b = 2 * k_z * k_z * (pos.x * vel.x + pos.y * vel.y) -
               2 * k_r * vel.z * (k_z * r_a + k_r * (pos.z - z_a));
    real_t c = k_z * k_z * (pos.x * pos.x + pos.y * pos.y) -
               std::pow(k_z * r_a + k_r * (pos.z - z_a), 2);

    // Solve the quadratic equation
    real_t discriminant = b * b - 4 * a * c;
    real_t t1;
    real_t t2;
    if (discriminant < 0.) {
        // No intersection, return negative number because it will be ignored.
        return -1;
    }
    else if (std::abs(a) < 1.0e-10) {
        // Only one solution
        // Use a different formula to avoid rounding errors
        t1 = -c / b;
        t2 = -c / b;
    }
    else if (std::abs(c) < 1.0e-10) {
        // The initial position is on the cone (within tolerance)
        // One solution is at the initial position
        // Use a different formula to avoid rounding errors
        t1 = -b / a;
        t2 = 0;
    }
    else {
        // Normal case, just calculate both solutions
        if (b >= 0.) {
            t1 = (-b - std::sqrt(discriminant)) / (2 * a);
            t2 = c / (a * t1);
        }
        else {
            t1 = 2 * c / (-b + std::sqrt(discriminant));
            t2 = c / (a * t1);
        }
    }

    // Check if the intersections are within the bounds of the cone
    real_t z_new = pos.z + t1 * vel.z;
    if (z_new < z_a || z_new > z_b) {
        t1 = -1;
    }
    z_new = pos.z + t2 * vel.z;
    if (z_new < z_a || z_new > z_b) {
        t2 = -1;
    }

    // Return the smallest positive solution
    if (t1 < 0 && t2 < 0) {
        return -1;
    }
    else if (t1 < 0) {
        return t2;
    }
    else if (t2 < 0) {
        return t1;
    }
    else {
        return std::min(t1, t2);
    }
}

}  // namespace neptune