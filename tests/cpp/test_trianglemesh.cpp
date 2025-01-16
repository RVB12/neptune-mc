#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include "eiron2/common/vec2.hpp"
#include "eiron2/geometry/trianglemesh.hpp"

using namespace eiron2;
TEST_CASE("TriangleMesh::n_cells", "[TriangleMesh]") {
    vector<Vec2> vertices = {Vec2(0, 0), Vec2(1, 0), Vec2(0, 1)};
    vector<Triangle> triangles = {Triangle(0, 1, 2)};
    TriangleMesh mesh(vertices, triangles);

    REQUIRE(mesh.n_cells() == 1);
}

TEST_CASE("TriangleMesh::get_triangle", "[TriangleMesh]") {
    vector<Vec2> vertices = {Vec2(0, 0), Vec2(1, 0), Vec2(0, 1)};
    vector<Triangle> triangles = {Triangle(0, 1, 2)};
    TriangleMesh mesh(vertices, triangles);

    Triangle t = mesh.get_triangle(0);
    REQUIRE(t.ind1 == 0);
    REQUIRE(t.ind2 == 1);
    REQUIRE(t.ind3 == 2);
}

TEST_CASE("TriangleMesh::get_cell_index", "[TriangleMesh]") {
    vector<Vec2> vertices = {Vec2(0, 0), Vec2(1, 0), Vec2(0, 1), Vec2(1, 1)};
    vector<Triangle> triangles = {Triangle(0, 1, 2), Triangle(1, 3, 2)};
    TriangleMesh mesh(vertices, triangles);

    REQUIRE(mesh.get_cell_index(Vec2(0.1, 0.1)) == 0);
    REQUIRE(mesh.get_cell_index(Vec2(0.5, 0.5)) == 0);
    REQUIRE(mesh.get_cell_index(Vec2(0.9, 0.9)) == 1);
    REQUIRE(mesh.get_cell_index(Vec2(0.5, 0.1)) == 0);
    REQUIRE(mesh.get_cell_index(Vec2(0.1, 0.5)) == 0);
    REQUIRE(mesh.get_cell_index(Vec2(0.9, 0.12)) == 1);
    REQUIRE(mesh.get_cell_index(Vec2(0.1, 0.9)) == 1);
    REQUIRE(mesh.get_cell_index(Vec2(0.5, 0.9)) == 1);
    REQUIRE(mesh.get_cell_index(Vec2(0.9, 0.5)) == 1);
}

TEST_CASE("TriangleMesh::is_in_cell", "[TriangleMesh]") {
    vector<Vec2> vertices = {Vec2(0, 0), Vec2(1, 0), Vec2(0, 1), Vec2(1, 1)};
    vector<Triangle> triangles = {Triangle(0, 1, 2), Triangle(1, 3, 2)};
    TriangleMesh mesh(vertices, triangles);

    REQUIRE(mesh.is_in_cell(Vec2(0.1, 0.1), 0));
    REQUIRE(mesh.is_in_cell(Vec2(0.5, 0.5), 0));
    REQUIRE(mesh.is_in_cell(Vec2(0.9, 0.9), 1));
    REQUIRE(mesh.is_in_cell(Vec2(0.5, 0.1), 0));
    REQUIRE(mesh.is_in_cell(Vec2(0.1, 0.5), 0));
    REQUIRE(mesh.is_in_cell(Vec2(0.9, 0.1), 1));
    REQUIRE(mesh.is_in_cell(Vec2(0.1, 0.9), 1));
    REQUIRE(mesh.is_in_cell(Vec2(0.5, 0.9), 1));
    REQUIRE(mesh.is_in_cell(Vec2(0.9, 0.5), 1));
}

TEST_CASE("TriangleMesh::get_volume", "[TriangleMesh]") {
    vector<Vec2> vertices = {Vec2(0, 0), Vec2(2, 0), Vec2(1, 2), Vec2(3, 1)};
    vector<Triangle> triangles = {Triangle(0, 1, 2), Triangle(1, 3, 2)};
    TriangleMesh mesh(vertices, triangles);

    REQUIRE(mesh.get_volume(0) == 2.0);
    REQUIRE(mesh.get_volume(1) == 1.5);
}

TEST_CASE("TriangleMesh::get_next_crossing", "[TriangleMesh]") {
    vector<Vec2> vertices = {Vec2(0, 0), Vec2(1, 0), Vec2(0, 1), Vec2(1, 1)};
    vector<Triangle> triangles = {Triangle(0, 1, 2), Triangle(1, 3, 2)};
    // vector<index_t> neighbors = {-1, 1, -1, -1, -1, 0};
    // TriangleMesh mesh(vertices, triangles, neighbors);
    TriangleMesh mesh(vertices, triangles);

    ParticleGeometry p(Vec2(0.1, 0.1), Vec2(1, 1), 0);
    CellCrossing c = mesh.get_next_crossing(p);

    Vec2 expected_pos = Vec2(0.5, 0.5);
    REQUIRE(c.position.x == expected_pos.x);
    REQUIRE(c.position.y == expected_pos.y);
    REQUIRE(c.t == 0.4);
    REQUIRE(c.next_cell_id == 1);
    REQUIRE(c.boundary_id == -1);
}