#ifndef EIRON2_Vec2_HPP
#define EIRON2_Vec2_HPP

#include <cassert>
#include <ostream>

#include "eiron2/common/utils.hpp"

namespace eiron2 {

/**
 * @brief Simple class that represents a 2D vector
 */
class Vec2 {
   public:
    // Contructors
    Vec2() = default;
    Vec2(real_t x_, real_t y_) : x(x_), y(y_){};

    // Member variables
    real_t x = 0.;
    real_t y = 0.;

    // Operators
    real_t& operator[](index_t);
    const real_t& operator[](index_t) const;

    Vec2& operator+=(Vec2);
    Vec2& operator-=(Vec2);
    Vec2& operator*=(Vec2);
    Vec2& operator/=(Vec2);

    Vec2& operator+=(real_t);
    Vec2& operator-=(real_t);
    Vec2& operator*=(real_t);
    Vec2& operator/=(real_t);

    Vec2 operator-();

    // Calculate the dot product of two vectors
    real_t dot(Vec2) const;

    // Calculate the distance between two vectors
    real_t dist(Vec2) const;

    // Calculate the norm of the vector
    real_t norm() const;

    // Calculate the squared norm of the vector
    real_t norm2() const;
};

std::ostream& operator<<(std::ostream& os, Vec2 const& v);

// Define some wrappers around the operators
inline Vec2 operator+(Vec2 a, Vec2 b) { return a += b; }
inline Vec2 operator-(Vec2 a, Vec2 b) { return a -= b; }
inline Vec2 operator*(Vec2 a, Vec2 b) { return a *= b; }
inline Vec2 operator/(Vec2 a, Vec2 b) { return a /= b; }

inline Vec2 operator+(Vec2 a, real_t s) { return a += s; }
inline Vec2 operator-(Vec2 a, real_t s) { return a -= s; }
inline Vec2 operator*(Vec2 a, real_t s) { return a *= s; }
inline Vec2 operator/(Vec2 a, real_t s) { return a /= s; }

inline Vec2 operator+(real_t s, Vec2 a) { return a += s; }
inline Vec2 operator-(real_t s, Vec2 a) { return a -= s; }
inline Vec2 operator*(real_t s, Vec2 a) { return a *= s; }
inline Vec2 operator/(real_t s, Vec2 a) { return a /= s; }

inline real_t dot(Vec2 a, Vec2 b) { return a.dot(b); }

inline real_t dist(Vec2 a, Vec2 b) { return a.dist(b); }

}  // namespace eiron2

#endif