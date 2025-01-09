

#include <cmath>
#include <cassert>

#include "eiron2/common/utils.hpp"
#include "eiron2/common/vec2.hpp"

namespace eiron2 {

double& Vec2::operator[](int i)
{
    assert(i >= 0);
    assert(i < 2);
    switch (i) {
        case 0:
            return x;
        case 1:
            return y;
    }
}

const double& Vec2::operator[](int i) const
{
    assert(i >= 0);
    assert(i < 2);
    switch (i) {
        case 0:
            return x;
        case 1:
            return y;
    }
}

Vec2& Vec2::operator+=(Vec2 other)
{
    x += other.x;
    y += other.y;
    return *this;
}

Vec2& Vec2::operator-=(Vec2 other)
{
    x -= other.x;
    y -= other.y;
    return *this;
}

Vec2& Vec2::operator*=(Vec2 other)
{
    x *= other.x;
    y *= other.y;
    return *this;
}

Vec2& Vec2::operator/=(Vec2 other)
{
    x /= other.x;
    y /= other.y;
    return *this;
}

Vec2& Vec2::operator+=(double s)
{
    x += s;
    y += s;
    return *this;
}

Vec2& Vec2::operator-=(double s)
{
    x -= s;
    y -= s;
    return *this;
}

Vec2& Vec2::operator*=(double s)
{
    x *= s;
    y *= s;
    return *this;
}

Vec2& Vec2::operator/=(double s)
{
    x /= s;
    y /= s;
    return *this;
}


Vec2 Vec2::operator-()
{
    return Vec2(-x,-y);
}

double Vec2::dot(Vec2 other) const { return x * other.x + y * other.y; }

double Vec2::dist(Vec2 other) const {
    double d2 = square(x - other.x) + square(y - other.y);
    return std::sqrt(d2);    
}

double Vec2::norm2() const {
    return square(x) + square(y);
}

double Vec2::norm() const {
    return std::sqrt(norm2());    
}

std::ostream& operator<<(std::ostream& os, Vec2 const& v)
{
    os << "[2]" "](" << v.x << ", " << v.y << ")";
    return os;
}

}  // namespace eiron2