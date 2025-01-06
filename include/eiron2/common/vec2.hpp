#ifndef EIRON2_Vec2_HPP
#define EIRON2_Vec2_HPP

#include <cassert>

namespace eiron2 {

/**
 * @brief Simple class that represents a 3D vector
 */
class Vec2 {
   public:
    // Contructors
    Vec2() = default;
    Vec2(double x_, double y_) : x(x_), y(y_){};

    // Member variables
    double x = 0.;
    double y = 0.;

    // Operators
    double& operator[](int);
    const double& operator[](int) const;

    Vec2& operator+=(Vec2);
    Vec2& operator-=(Vec2);
    Vec2& operator*=(Vec2);
    Vec2& operator/=(Vec2);

    Vec2& operator+=(double);
    Vec2& operator-=(double);
    Vec2& operator*=(double);
    Vec2& operator/=(double);

    Vec2 operator-();

    // Calculate the dot product of two vectors
    double dot(Vec2) const;

    // Calculate the distance between two vectors
    double dist(Vec2) const;

    // Calculate the norm of the vector
    double norm() const;

    // Calculate the squared norm of the vector
    double norm2() const;
};

// Define some wrappers around the operators
inline Vec2 operator+(Vec2 a, Vec2 b) { return a += b; }
inline Vec2 operator-(Vec2 a, Vec2 b) { return a -= b; }
inline Vec2 operator*(Vec2 a, Vec2 b) { return a *= b; }
inline Vec2 operator/(Vec2 a, Vec2 b) { return a /= b; }

inline Vec2 operator+(Vec2 a, double s) { return a += s; }
inline Vec2 operator-(Vec2 a, double s) { return a -= s; }
inline Vec2 operator*(Vec2 a, double s) { return a *= s; }
inline Vec2 operator/(Vec2 a, double s) { return a /= s; }

inline Vec2 operator+(double s, Vec2 a) { return a += s; }
inline Vec2 operator-(double s, Vec2 a) { return a -= s; }
inline Vec2 operator*(double s, Vec2 a) { return a *= s; }
inline Vec2 operator/(double s, Vec2 a) { return a /= s; }

inline double dot(Vec2 a, Vec2 b) { return a.dot(b); }

inline double dist(Vec2 a, Vec2 b) { return a.dist(b); }

}  // namespace eiron2

#endif