#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include "eigrene/geometry/vec3.hpp"

using namespace eigrene;

TEST_CASE("Default vec3 constructor works", "[vector]")
{
    Vec3 a;
    REQUIRE(a.x == 0);
    REQUIRE(a.y == 0);
    REQUIRE(a.z == 0);
}

TEST_CASE("vec3 constructor works", "[vector]")
{
    Vec3 a(3.4,8.2,9.7);
    REQUIRE(a.x == 3.4);
    REQUIRE(a.y == 8.2);
    REQUIRE(a.z == 9.7);
}


TEST_CASE("vec3 operators work", "[vector]")
{
    Vec3 a(1,2,3);
    Vec3 b(4,8,9);

    SECTION("addition"){
        auto c = a + b;
        REQUIRE(c.x == 5);
        REQUIRE(c.y == 10);
        REQUIRE(c.z == 12);

        a += b;
        REQUIRE(a.x == 5);
        REQUIRE(a.y == 10);
        REQUIRE(a.z == 12);
    }

    SECTION("subtraction"){
        auto c = a - b;
        REQUIRE(c.x == -3);
        REQUIRE(c.y == -6);
        REQUIRE(c.z == -6);

        a -= b;
        REQUIRE(a.x == -3);
        REQUIRE(a.y == -6);
        REQUIRE(a.z == -6);
    }

    SECTION("multiplication"){
        auto c = a * b;
        REQUIRE(c.x == 4);
        REQUIRE(c.y == 16);
        REQUIRE(c.z == 27);

        a *= b;
        REQUIRE(a.x == 4);
        REQUIRE(a.y == 16);
        REQUIRE(a.z == 27);
    }

    SECTION("division"){
        auto c = a / b;
        REQUIRE(c.x == 1./4.);
        REQUIRE(c.y == 2./8.);
        REQUIRE(c.z == 3./9.);

        a /= b;
        REQUIRE(a.x == 1./4.);
        REQUIRE(a.y == 2./8.);
        REQUIRE(a.z == 3./9.);
    }
}
