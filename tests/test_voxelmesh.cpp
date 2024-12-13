#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include "eiron2/geometry/vec3.hpp"
#include "eiron2/geometry/voxelmesh.hpp"

using namespace eiron2;

TEST_CASE("voxelmesh constructor works", "[voxelmesh][mesh]")
{
    VoxelMesh mesh(1., 1., 1., 2, 2, 2);
}

TEST_CASE("voxelmesh index conversion works for square mesh",
          "[voxelmesh][mesh]")
{
    VoxelMesh mesh(1., 1., 1., 3, 3, 3);

    for (int i = 0; i < mesh.n_cells(); i++) {
        Index3D ind = mesh.linear_to_triplet(i);

        int i2 = mesh.triplet_to_linear(ind);

        REQUIRE(i == i2);
    }
}

TEST_CASE("voxelmesh index conversion works for rectangular mesh",
          "[voxelmesh][mesh]")
{
    VoxelMesh mesh(1., 1., 1., 3, 6, 5);

    for (int i = 0; i < mesh.n_cells(); i++) {
        Index3D ind = mesh.linear_to_triplet(i);

        int i2 = mesh.triplet_to_linear(ind);

        REQUIRE(i == i2);
    }
}

TEST_CASE("check the cell crossings with some manually verified examples",
          "[voxelmesh][mesh]")
{
    VoxelMesh mesh(1., 1.6, 2.1, 10, 15, 13);

    Vec3 pos(0.2, 0.3, 0.4);
    Vec3 vel(0.1, 0.4, -0.8);

    Particle p(pos, vel, mesh.find_cell(pos));

    auto crossing = mesh.get_next_crossing(p);
}
