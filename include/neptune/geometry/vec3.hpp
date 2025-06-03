#ifndef NEPTUNE_VEC3_HPP
#define NEPTUNE_VEC3_HPP

#include <cassert>
#include <ostream>

#include "neptune/common/utils.hpp"

namespace neptune {

/**
 * @brief Simple class that represents a 3D vector
 */
class Vec3 {
   public:
    // Contructors
    Vec3() = default;
    Vec3(real_t x_, real_t y_, real_t z_) : x(x_), y(y_), z(z_){};

    // Member variables
    real_t x = 0.;
    real_t y = 0.;
    real_t z = 0.;

    // Operators
    real_t& operator[](index_t);
    const real_t& operator[](index_t) const;

    Vec3& operator+=(Vec3);
    Vec3& operator-=(Vec3);
    Vec3& operator*=(Vec3);
    Vec3& operator/=(Vec3);

    Vec3& operator+=(real_t);
    Vec3& operator-=(real_t);
    Vec3& operator*=(real_t);
    Vec3& operator/=(real_t);

    Vec3 operator-();

    // Calculate the dot product of two vectors
    real_t dot(Vec3) const;

    // Calculate the distance between two vectors
    real_t dist(Vec3) const;

    // Calculate the norm of the vector
    real_t norm() const;

    // Calculate the squared norm of the vector
    real_t norm2() const;

    // Calculate the cross product of two vectors
    Vec3 cross(Vec3) const;
};

// Define some wrappers around the operators
inline Vec3 operator+(Vec3 a, Vec3 b) { return a += b; }
inline Vec3 operator-(Vec3 a, Vec3 b) { return a -= b; }
inline Vec3 operator*(Vec3 a, Vec3 b) { return a *= b; }
inline Vec3 operator/(Vec3 a, Vec3 b) { return a /= b; }

inline Vec3 operator+(Vec3 a, real_t s) { return a += s; }
inline Vec3 operator-(Vec3 a, real_t s) { return a -= s; }
inline Vec3 operator*(Vec3 a, real_t s) { return a *= s; }
inline Vec3 operator/(Vec3 a, real_t s) { return a /= s; }

inline Vec3 operator+(real_t s, Vec3 a) { return a += s; }
inline Vec3 operator-(real_t s, Vec3 a) { return a -= s; }
inline Vec3 operator*(real_t s, Vec3 a) { return a *= s; }
inline Vec3 operator/(real_t s, Vec3 a) { return a /= s; }

inline real_t dot(Vec3 a, Vec3 b) { return a.dot(b); }

inline real_t dist(Vec3 a, Vec3 b) { return a.dist(b); }

}  // namespace neptune

#endif