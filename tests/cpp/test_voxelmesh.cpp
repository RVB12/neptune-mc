#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include "neptune/geometry/vec2.hpp"
#include "neptune/geometry/voxelmesh.hpp"

using namespace neptune;

TEST_CASE("voxelmesh constructor works", "[voxelmesh][mesh]")
{
    VoxelMesh mesh(1., 1., 1., 2, 2, 2);
}

TEST_CASE("voxelmesh index conversion works for square mesh",
          "[voxelmesh][mesh]")
{
    VoxelMesh mesh(1., 1., 1., 3, 3, 3);

    for (index_t i = 0; i < mesh.n_cells(); i++) {
        Index3D ind = mesh.linear_to_tuple(i);

        index_t i2 = mesh.tuple_to_linear(ind);

        CHECK(i == i2);
    }
}

TEST_CASE("voxelmesh index conversion works for rectangular mesh",
          "[voxelmesh][mesh]")
{
    VoxelMesh mesh(1., 1., 1., 3, 6, 9);

    for (index_t i = 0; i < mesh.n_cells(); i++) {
        Index3D ind = mesh.linear_to_tuple(i);

        index_t i2 = mesh.tuple_to_linear(ind);

        CHECK(i == i2);
    }
}

TEST_CASE("check the cell crossings with some manually verified examples",
          "[voxelmesh][mesh]")
{
    VoxelMesh mesh(1., 1., 1., 10, 5, 5);
    Vec3 cell_size = mesh.cell_size();

    Index3D ind(3, 2, 0);
    Vec3 pos_in_cell(0.05, 0.04, 0.);
    Vec3 pos_edge(ind.i * cell_size.x, ind.j * cell_size.y, ind.k * cell_size.z);

    Vec3 pos = pos_edge + pos_in_cell;

    SECTION("vel 1")
    {
        Vec3 vel(0.1, 0.1, 0.0);
        ParticleGeometry p(pos, vel, mesh.tuple_to_linear(ind));
        auto crossing = mesh.get_next_crossing(p);

        real_t t_expected = 0.5;
        CHECK(std::abs(crossing.t - t_expected) <= 0.0001);

        Index3D expected_ind = ind;
        expected_ind.i += 1;
        CHECK(crossing.next_cell_id == mesh.tuple_to_linear(expected_ind));

        Vec3 expected_pos = pos + t_expected * vel;
        CHECK((expected_pos - crossing.position).norm() <= 0.0001);

        CHECK(mesh.is_in_cell(crossing.position, crossing.next_cell_id));
    }

    SECTION("vel 2")
    {
        Vec3 vel(-0.2, -0.8, 0.0);
        ParticleGeometry p(pos, vel, mesh.tuple_to_linear(ind));
        auto crossing = mesh.get_next_crossing(p);

        real_t t_expected = 0.05;
        CHECK(std::abs(crossing.t - t_expected) <= 0.0001);

        Index3D expected_ind = ind;
        expected_ind.j -= 1;
        CHECK(crossing.next_cell_id == mesh.tuple_to_linear(expected_ind));

        Vec3 expected_pos = pos + t_expected * vel;
        CHECK((expected_pos - crossing.position).norm() <= 0.0001);

        CHECK(mesh.is_in_cell(crossing.position, crossing.next_cell_id));
    }

    // TODO: check the boundary crossings
}
