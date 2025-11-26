#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include "neptune/geometry/rotatedtrianglemesh.hpp"
#include "neptune/geometry/vec2.hpp"

using namespace neptune;

void check_cone_intersection(const Vec2& p1,
                             const Vec2& p2,
                             const Vec3& position,
                             const Vec3& velocity,
                             bool is_on_surface,
                             bool expect_intersection = true)
{
    real_t t = find_cone_intersection(p1, p2, position, velocity, false);

    if(t == -1) {
        REQUIRE(!expect_intersection);
        return;
    }

    REQUIRE(t > 0);

    Vec3 intersection_point = position + velocity * t;

    // The intersection point should lie on the cone surface
    real_t r = std::sqrt(intersection_point.x * intersection_point.x +
                         intersection_point.y * intersection_point.y);
    real_t z = intersection_point.z;
    real_t r_cone = std::abs((p2.x - p1.x) / (p2.y - p1.y) * (z - p1.y) + p1.x);
    REQUIRE(std::abs(r - r_cone) < 1e-14);
}

TEST_CASE("Cone intersection:Simple cases", "[RotatedTriangleMesh]")
{
    Vec2 p1(1.0, 0.0);
    Vec2 p2(0.0, 1.0);
    SECTION("Starting at zero")
    {
        Vec3 position(0.0, 0.0, 0.0);
        SECTION("Along x-axis")
        {
            Vec3 velocity(1.0, 0.0, 0.0);
            check_cone_intersection(p1, p2, position, velocity, false);
        }

        SECTION("Along y-axis")
        {
            Vec3 velocity(0.0, 1.0, 0.0);
            check_cone_intersection(p1, p2, position, velocity, false);
        }

        SECTION("Along z-axis")
        {
            Vec3 velocity(0.0, 0.0, 1.0);
            check_cone_intersection(p1, p2, position, velocity, false);
        }
    }

    SECTION("On surface")
    {
        Vec3 position(-1.0, 0.0, 0.0);
        SECTION("Along x-axis")
        {
            Vec3 velocity(1.0, 0.0, 0.0);
            check_cone_intersection(p1, p2, position, velocity, true, false);
        }

        SECTION("Along y-axis")
        {
            Vec3 velocity(0.0, 1.0, 0.0);
            check_cone_intersection(p1, p2, position, velocity, true, false);
        }

        SECTION("Along z-axis")
        {
            Vec3 velocity(0.0, 0.0, 1.0);
            check_cone_intersection(p1, p2, position, velocity, true, false);
        }
    }

    SECTION("Almost on surface")
    {
        Vec3 position(-0.9999, 0.0, 0.0);
        SECTION("Along x-axis")
        {
            Vec3 velocity(1.0, 0.0, 0.0);
            check_cone_intersection(p1, p2, position, velocity, false, false);
        }

        SECTION("Along y-axis")
        {
            Vec3 velocity(0.0, 1.0, 0.0);
            check_cone_intersection(p1, p2, position, velocity, false, false);
        }

        SECTION("Along z-axis")
        {
            Vec3 velocity(0.0, 0.0, 1.0);
            check_cone_intersection(p1, p2, position, velocity, false, false);
        }
    }
}


TEST_CASE("Cone test case that failed", "[RotatedTriangleMesh]")
{
    Vec2 p1(1.267430427433943, -0.016825296355507825);
    Vec2 p2(1.2840816326530613, 0);
    Vec3 position(1.2641767797715269, 0.0069106210230192432, -0.0014997451417587322);
    Vec3 velocity(-7062.1497659581173, 11068.823255313708, 6191.2515955083018);

    check_cone_intersection(p1, p2, position, velocity, false);
}