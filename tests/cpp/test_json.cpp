#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include "neptune/common/json_conversions.hpp"
#include "neptune/geometry/vec2.hpp"
#include "neptune/geometry/vec3.hpp"

using namespace neptune;

TEST_CASE("Vec2 json conversion works", "[vector][json]")
{
    Vec2 a(1.5, 2.5);

    nlohmann::json j = a;

    Vec2 b = j;

    REQUIRE(a.x == b.x);
    REQUIRE(a.y == b.y);

    std::string s_expected = R"([1.5,2.5])";
    REQUIRE(j.dump() == s_expected);
}

TEST_CASE("Vec3 json conversion works", "[vector][json]")
{
    Vec3 a(2.1, 3.7, 9.1);

    nlohmann::json j = a;

    Vec3 b = j;

    REQUIRE(a.x == b.x);
    REQUIRE(a.y == b.y);
    REQUIRE(a.z == b.z);

    std::string s_expected = R"([2.1,3.7,9.1])";
    REQUIRE(j.dump() == s_expected);
}

TEST_CASE("Vector json conversion works", "[vector][json]")
{
    vector<double> a = {1.5, 2.5, 3.5};

    nlohmann::json j = a;

    vector<double> b = j;

    REQUIRE(a.size() == b.size());
    for (index_t i = 0; i < a.size(); i++) {
        REQUIRE(a[i] == b[i]);
    }

    std::string s_expected = R"([1.5,2.5,3.5])";
    REQUIRE(j.dump() == s_expected);
}

// TEST_CASE("LineSource json conversion works", "[vector][json]")
// {
//     LineSource a(10, 100, Vec3(1.5, 2.5, 15.4), Vec3(3.5, 4.5, 1.5), 0.1, 0.2);

//     nlohmann::json j = a;

//     LineSource b = j;

//     REQUIRE(a.pos_a.x == b.pos_a.x);
//     REQUIRE(a.pos_a.y == b.pos_a.y);
//     REQUIRE(a.pos_b.x == b.pos_b.x);
//     REQUIRE(a.pos_b.y == b.pos_b.y);
//     REQUIRE(a.angle_min == b.angle_min);
//     REQUIRE(a.angle_max == b.angle_max);
//     REQUIRE(a.num_particles() == b.num_particles());

//     std::string s_expected =
//         R"({"angle_max":0.2,"angle_min":0.1,"n_particles":10,"pos_a":[1.5,2.5],"pos_b":[3.5,4.5],"type":"line"})";
//     REQUIRE(j.dump() == s_expected);
// }

// TEST_CASE("Source json conversion works", "[vector][json]")
// {
//     LineSource a(10, 100, Vec3(1.5, 2.5, 15.4), Vec3(3.5, 4.5, 1.5), 0.1, 0.2);

//     nlohmann::json j = a;

//     std::unique_ptr<Source> b = j;

//     LineSource* c = dynamic_cast<LineSource*>(b.get());
//     REQUIRE(c != nullptr);
//     REQUIRE(a.pos_a.x == c->pos_a.x);
//     REQUIRE(a.pos_a.y == c->pos_a.y);
//     REQUIRE(a.pos_b.x == c->pos_b.x);
//     REQUIRE(a.pos_b.y == c->pos_b.y);
//     REQUIRE(a.angle_min == c->angle_min);
//     REQUIRE(a.angle_max == c->angle_max);
//     REQUIRE(a.num_particles() == c->num_particles());

//     std::string s_expected =
//         R"({"angle_max":0.2,"angle_min":0.1,"n_particles":10,"pos_a":[1.5,2.5],"pos_b":[3.5,4.5],"type":"line"})";
//     REQUIRE(j.dump() == s_expected);
// }

// TEST_CASE("VoxelMesh json conversion works", "[vector][json]")
// {
//     VoxelMesh a(1.5, 2.5, 3.5, 10, 20, 30);

//     nlohmann::json j = a;

//     VoxelMesh b = j;

//     REQUIRE(a.Lx == b.Lx);
//     REQUIRE(a.Ly == b.Ly);
//     REQUIRE(a.Lz == b.Lz);
//     REQUIRE(a.nx == b.nx);
//     REQUIRE(a.ny == b.ny);
//     REQUIRE(a.nz == b.nz);
//     REQUIRE(a.dx == b.dx);
//     REQUIRE(a.dy == b.dy);
//     REQUIRE(a.dz == b.dz);

//     std::string s_expected =
//         R"({"Lx":1.5,"Ly":2.5,"nx":10,"ny":20,"type":"voxel"})";
//     REQUIRE(j.dump() == s_expected);
// }