#ifndef EIRON2_VEC3_HPP
#define EIRON2_VEC3_HPP

#include <cassert>

namespace eiron2 {

/**
 * @brief Simple class that represents a 3D vector
 */
class Vec3 {
   public:
    // Contructors
    Vec3() = default;
    Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_){};

    // Member variables
    double x = 0.;
    double y = 0.;
    double z = 0.;

    // Operators
    double& operator[](int);
    const double& operator[](int) const;

    Vec3& operator+=(Vec3);
    Vec3& operator-=(Vec3);
    Vec3& operator*=(Vec3);
    Vec3& operator/=(Vec3);

    Vec3& operator+=(double);
    Vec3& operator-=(double);
    Vec3& operator*=(double);
    Vec3& operator/=(double);

    Vec3 operator-();

    // Calculate the dot product of two vectors
    double dot(Vec3) const;

    // Calculate the distance between two vectors
    double dist(Vec3) const;

    // Calculate the norm of the vector
    double norm() const;

    // Calculate the squared norm of the vector
    double norm2() const;
};

// Define some wrappers around the operators
inline Vec3 operator+(Vec3 a, Vec3 b) { return a += b; }
inline Vec3 operator-(Vec3 a, Vec3 b) { return a -= b; }
inline Vec3 operator*(Vec3 a, Vec3 b) { return a *= b; }
inline Vec3 operator/(Vec3 a, Vec3 b) { return a /= b; }

inline Vec3 operator+(Vec3 a, double s) { return a += s; }
inline Vec3 operator-(Vec3 a, double s) { return a -= s; }
inline Vec3 operator*(Vec3 a, double s) { return a *= s; }
inline Vec3 operator/(Vec3 a, double s) { return a /= s; }

inline Vec3 operator+(double s, Vec3 a) { return a += s; }
inline Vec3 operator-(double s, Vec3 a) { return a -= s; }
inline Vec3 operator*(double s, Vec3 a) { return a *= s; }
inline Vec3 operator/(double s, Vec3 a) { return a /= s; }

inline double dot(Vec3 a, Vec3 b) { return a.dot(b); }

inline double dist(Vec3 a, Vec3 b) { return a.dist(b); }

}  // namespace eiron2

#endif