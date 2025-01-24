#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include "eiron2/geometry/vec2.hpp"

using namespace eiron2;

TEST_CASE("Default vec2 constructor works", "[vector]")
{
    Vec2 a;
    REQUIRE(a.x == 0);
    REQUIRE(a.y == 0);
}

TEST_CASE("vec2 constructor works", "[vector]")
{
    Vec2 a(3.4,8.2);
    REQUIRE(a.x == 3.4);
    REQUIRE(a.y == 8.2);
}


TEST_CASE("vec2 operators work", "[vector]")
{
    Vec2 a(1,2);
    Vec2 b(4,8);

    SECTION("addition"){
        auto c = a + b;
        REQUIRE(c.x == 5);
        REQUIRE(c.y == 10);

        a += b;
        REQUIRE(a.x == 5);
        REQUIRE(a.y == 10);
    }

    SECTION("subtraction"){
        auto c = a - b;
        REQUIRE(c.x == -3);
        REQUIRE(c.y == -6);

        a -= b;
        REQUIRE(a.x == -3);
        REQUIRE(a.y == -6);
    }

    SECTION("multiplication"){
        auto c = a * b;
        REQUIRE(c.x == 4);
        REQUIRE(c.y == 16);

        a *= b;
        REQUIRE(a.x == 4);
        REQUIRE(a.y == 16);
    }

    SECTION("division"){
        auto c = a / b;
        REQUIRE(c.x == 1./4.);
        REQUIRE(c.y == 2./8.);

        a /= b;
        REQUIRE(a.x == 1./4.);
        REQUIRE(a.y == 2./8.);
    }
}
