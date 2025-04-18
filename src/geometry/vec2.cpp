

#include <cmath>
#include <cassert>

#include "neptune/common/utils.hpp"
#include "neptune/geometry/vec2.hpp"

namespace neptune {

real_t& Vec2::operator[](index_t i)
{
    assert(i >= 0);
    assert(i < 2);
    switch (i) {
        case 0:
            return x;
        case 1:
            return y;
        default:
            throw std::runtime_error("Index out of bounds");
    }
}

const real_t& Vec2::operator[](index_t i) const
{
    assert(i >= 0);
    assert(i < 2);
    switch (i) {
        case 0:
            return x;
        case 1:
            return y;
        default:
            throw std::runtime_error("Index out of bounds");
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

Vec2& Vec2::operator+=(real_t s)
{
    x += s;
    y += s;
    return *this;
}

Vec2& Vec2::operator-=(real_t s)
{
    x -= s;
    y -= s;
    return *this;
}

Vec2& Vec2::operator*=(real_t s)
{
    x *= s;
    y *= s;
    return *this;
}

Vec2& Vec2::operator/=(real_t s)
{
    x /= s;
    y /= s;
    return *this;
}


Vec2 Vec2::operator-()
{
    return Vec2(-x,-y);
}

real_t Vec2::dot(Vec2 other) const { return x * other.x + y * other.y; }

real_t Vec2::dist(Vec2 other) const {
    real_t d2 = square(x - other.x) + square(y - other.y);
    return std::sqrt(d2);    
}

real_t Vec2::norm2() const {
    return square(x) + square(y);
}

real_t Vec2::norm() const {
    return std::sqrt(norm2());    
}

std::ostream& operator<<(std::ostream& os, Vec2 const& v)
{
    os << "[2]" "](" << v.x << ", " << v.y << ")";
    return os;
}

}  // namespace neptune