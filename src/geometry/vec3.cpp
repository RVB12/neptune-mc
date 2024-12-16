

#include <cmath>
#include <cassert>

#include "eiron2/common/utils.hpp"
#include "eiron2/common/vec3.hpp"

namespace eiron2 {

double& Vec3::operator[](int i)
{
    assert(i >= 0);
    assert(i < 3);
    switch (i) {
        case 0:
            return x;
        case 1:
            return y;
        case 2:
            return z;
    }
}

const double& Vec3::operator[](int i) const
{
    assert(i >= 0);
    assert(i < 3);
    switch (i) {
        case 0:
            return x;
        case 1:
            return y;
        case 2:
            return z;
    }
}

Vec3& Vec3::operator+=(Vec3 other)
{
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
}

Vec3& Vec3::operator-=(Vec3 other)
{
    x -= other.x;
    y -= other.y;
    z -= other.z;
    return *this;
}

Vec3& Vec3::operator*=(Vec3 other)
{
    x *= other.x;
    y *= other.y;
    z *= other.z;
    return *this;
}

Vec3& Vec3::operator/=(Vec3 other)
{
    x /= other.x;
    y /= other.y;
    z /= other.z;
    return *this;
}

Vec3& Vec3::operator+=(double s)
{
    x += s;
    y += s;
    z += s;
    return *this;
}

Vec3& Vec3::operator-=(double s)
{
    x -= s;
    y -= s;
    z -= s;
    return *this;
}

Vec3& Vec3::operator*=(double s)
{
    x *= s;
    y *= s;
    z *= s;
    return *this;
}

Vec3& Vec3::operator/=(double s)
{
    x /= s;
    y /= s;
    z /= s;
    return *this;
}


Vec3 Vec3::operator-()
{
    return Vec3(-x,-y,-z);
}

double Vec3::dot(Vec3 other) const { return x * other.x + y * other.y + z * other.z; }

double Vec3::dist(Vec3 other) const {
    double d2 = square(x - other.x) + square(y - other.y) + square(z - other.z);
    return std::sqrt(d2);    
}

double Vec3::norm2() const {
    return square(x) + square(y) + square(z);
}

double Vec3::norm() const {
    return std::sqrt(norm2());    
}

}  // namespace eiron2