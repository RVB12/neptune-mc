#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include "eiron2/geometry/vec2.hpp"
#include "eiron2/geometry/voxelmesh.hpp"

using namespace eiron2;

TEST_CASE("voxelmesh constructor works", "[voxelmesh][mesh]")
{
    VoxelMesh mesh(1., 1., 2, 2);
}

TEST_CASE("voxelmesh index conversion works for square mesh",
          "[voxelmesh][mesh]")
{
    VoxelMesh mesh(1., 1., 3, 3);

    for (index_t i = 0; i < mesh.n_cells(); i++) {
        Index2D ind = mesh.linear_to_tuple(i);

        index_t i2 = mesh.tuple_to_linear(ind);

        CHECK(i == i2);
    }
}

TEST_CASE("voxelmesh index conversion works for rectangular mesh",
          "[voxelmesh][mesh]")
{
    VoxelMesh mesh(1., 1., 3, 6);

    for (index_t i = 0; i < mesh.n_cells(); i++) {
        Index2D ind = mesh.linear_to_tuple(i);

        index_t i2 = mesh.tuple_to_linear(ind);

        CHECK(i == i2);
    }
}

TEST_CASE("check the cell crossings with some manually verified examples",
          "[voxelmesh][mesh]")
{
    VoxelMesh mesh(1., 1., 10, 5);
    Vec2 cell_size = mesh.cell_size();

    Index2D ind(3, 2);
    Vec2 pos_in_cell(0.05, 0.04);
    Vec2 pos_edge(ind.i * cell_size.x, ind.j * cell_size.y);

    Vec2 pos = pos_edge + pos_in_cell;

    SECTION("vel 1")
    {
        Vec2 vel(0.1, 0.1);
        ParticleGeometry p(pos, vel, mesh.tuple_to_linear(ind));
        auto crossing = mesh.get_next_crossing(p);

        real_t t_expected = 0.5;
        CHECK(std::abs(crossing.t - t_expected) <= 0.0001);

        Index2D expected_ind = ind;
        expected_ind.i += 1;
        CHECK(crossing.next_cell_id == mesh.tuple_to_linear(expected_ind));

        Vec2 expected_pos = pos + t_expected * vel;
        CHECK((expected_pos - crossing.position).norm() <= 0.0001);

        CHECK(mesh.is_in_cell(crossing.position, crossing.next_cell_id));
    }

    SECTION("vel 2")
    {
        Vec2 vel(-0.2, -0.8);
        ParticleGeometry p(pos, vel, mesh.tuple_to_linear(ind));
        auto crossing = mesh.get_next_crossing(p);

        real_t t_expected = 0.05;
        CHECK(std::abs(crossing.t - t_expected) <= 0.0001);

        Index2D expected_ind = ind;
        expected_ind.j -= 1;
        CHECK(crossing.next_cell_id == mesh.tuple_to_linear(expected_ind));

        Vec2 expected_pos = pos + t_expected * vel;
        CHECK((expected_pos - crossing.position).norm() <= 0.0001);

        CHECK(mesh.is_in_cell(crossing.position, crossing.next_cell_id));
    }

    // TODO: check the boundary crossings
}
