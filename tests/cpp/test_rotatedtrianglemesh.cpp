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
    real_t t =
        find_cone_intersection(p1, p2, position, velocity, is_on_surface);

    if (t == -1) {
        REQUIRE(!expect_intersection);
        return;
    }
    REQUIRE(expect_intersection);

    REQUIRE(t > 0);

    Vec3 intersection_point = position + velocity * t;

    real_t r = std::sqrt(intersection_point.x * intersection_point.x +
                         intersection_point.y * intersection_point.y);
    real_t z = intersection_point.z;
    // Check if cone is a plane perpendicular to the z axis
    if (std::abs(p2.y - p1.y) < 1.0e-10) {
        REQUIRE(std::abs(z - p1.y) < 1e-14);
        return;
    }
    else {
        real_t r_cone =
            std::abs((p2.x - p1.x) / (p2.y - p1.y) * (z - p1.y) + p1.x);
        REQUIRE(std::abs(r - r_cone) < 1e-14);
    }
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
            check_cone_intersection(p1, p2, position, Vec3(1.0, 0.0, 0.0),
                                    false);
        }

        SECTION("Along y-axis")
        {
            check_cone_intersection(p1, p2, position, Vec3(0.0, 1.0, 0.0),
                                    false);
        }

        SECTION("Along z-axis")
        {
            check_cone_intersection(p1, p2, position, Vec3(0.0, 0.0, 1.0),
                                    false);
        }

        SECTION("Diagonal toward the cone surface")
        {
            check_cone_intersection(p1, p2, position, Vec3(1.0, 0.0, 1.0),
                                    false);
        }

        SECTION("Diagonal away from the cone surface")
        {
            check_cone_intersection(p1, p2, position, Vec3(-1.0, 0.0, 1.0),
                                    false);
        }
    }

    SECTION("Starting outside the cone")
    {
        SECTION("Moving inward along x-axis")
        {
            check_cone_intersection(p1, p2, Vec3(-2.0, 0.0, 0.0),
                                    Vec3(1.0, 0.0, 0.0), false);
        }

        SECTION("Moving inward along y-axis")
        {
            check_cone_intersection(p1, p2, Vec3(0.0, 2.0, 0.0),
                                    Vec3(0.0, -1.0, 0.0), false);
        }

        SECTION("Moving inward while increasing z")
        {
            check_cone_intersection(p1, p2, Vec3(2.0, 0.0, 0.0),
                                    Vec3(-1.0, 0.0, 0.25), false);
        }

        SECTION("Moving away from the cone")
        {
            check_cone_intersection(p1, p2, Vec3(2.0, 0.0, 0.0),
                                    Vec3(1.0, 0.0, 0.0), false, false);
        }

        SECTION("Parallel to the cone surface")
        {
            check_cone_intersection(p1, p2, Vec3(2.0, 0.0, 0.0),
                                    Vec3(-1.0, 0.0, 1.0), false, false);
        }
    }

    SECTION("On surface")
    {
        Vec3 position(-1.0, 0.0, 0.0);

        SECTION("Along x-axis")
        {
            check_cone_intersection(p1, p2, position, Vec3(1.0, 0.0, 0.0), true,
                                    true);
            check_cone_intersection(p1, p2, position, Vec3(-1.0, 0.0, 0.0),
                                    true, false);
        }

        SECTION("Along y-axis")
        {
            check_cone_intersection(p1, p2, position, Vec3(0.0, 1.0, 0.0), true,
                                    false);
        }

        SECTION("Along z-axis")
        {
            check_cone_intersection(p1, p2, position, Vec3(0.0, 0.0, 1.0), true,
                                    false);
        }
    }

    SECTION("Almost on surface")
    {
        Vec3 position(-0.9999, 0.0, 0.0);

        SECTION("Along x-axis")
        {
            check_cone_intersection(p1, p2, position, Vec3(1.0, 0.0, 0.0),
                                    false, true);
            check_cone_intersection(p1, p2, position, Vec3(-1.0, 0.0, 0.0),
                                    false, true);
        }

        SECTION("Along y-axis")
        {
            check_cone_intersection(p1, p2, position, Vec3(0.0, 1.0, 0.0),
                                    false, true);
            check_cone_intersection(p1, p2, position, Vec3(0.0, -1.0, 0.0),
                                    false, true);
        }

        SECTION("Along z-axis")
        {
            check_cone_intersection(p1, p2, position, Vec3(0.0, 0.0, 1.0),
                                    false, true);
            check_cone_intersection(p1, p2, position, Vec3(0.0, 0.0, -1.0),
                                    false, false);
        }
    }
}

TEST_CASE("Cone intersection:Different cone profiles", "[RotatedTriangleMesh]")
{
    SECTION("Translated and widened cone")
    {
        Vec2 p1(0.5, -1.0);
        Vec2 p2(1.5, 1.0);

        check_cone_intersection(p1, p2, Vec3(2.0, 0.0, 0.0),
                                Vec3(-1.0, 0.0, 0.0), false);
    }

    SECTION("Diagonal intersection")
    {
        Vec2 p1(0.5, -1.0);
        Vec2 p2(1.5, 1.0);

        check_cone_intersection(p1, p2, Vec3(2.0, 0.0, -0.5),
                                Vec3(-1.0, 0.0, 0.5), false);
    }

    SECTION("No intersection outside the profile")
    {
        Vec2 p1(0.5, -1.0);
        Vec2 p2(1.5, 1.0);

        check_cone_intersection(p1, p2, Vec3(2.0, 0.0, 0.0),
                                Vec3(0.0, 0.0, 1.0), false, false);
    }

    SECTION("Horizontal cone")
    {
        Vec2 p1(0.5, 1.0);
        Vec2 p2(1.5, 1.0);

        check_cone_intersection(p1, p2, Vec3(1.0, 0.0, 0.0),
                                Vec3(0.0, 0.0, 1.0), false);
    }
}
