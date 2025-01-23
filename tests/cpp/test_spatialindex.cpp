#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers.hpp>

#include "eiron2/common/vec2.hpp"
#include "eiron2/geometry/2Dspatialindex.hpp"

using namespace eiron2;

TEST_CASE("Spatial index constructor works (only one value)",
          "[geometry][2Dspatialindex]")
{
    index_t nx = 10;
    index_t ny = 10;

    vector<BoundingBox> bboxes = {
        BoundingBox(Vec2(0.25, 0.26), Vec2(0.7, 0.8))};

    auto bbox_func = [&](index_t i) { return bboxes[i]; };

    SpatialIndex2D index(nx, ny, bboxes.size(), bbox_func);

    REQUIRE(index.nx == nx);
    REQUIRE(index.ny == ny);

    // There is only one value, so the bounding box should be the same as the
    // value
    REQUIRE(bboxes[0].min.x == index.bbox.min.x);
    REQUIRE(bboxes[0].min.y == index.bbox.min.y);

    // There is only one value, so each cell should contain that value
    REQUIRE(index.cell_idx.size() == nx * ny + 1);
    for(index_t i = 0; i < nx * ny + 1; i++) {
        REQUIRE(index.cell_idx[i] == i);
    }

    REQUIRE(index.data.size() == nx * ny);
    for(index_t i = 0; i < nx * ny; i++) {
        REQUIRE(index.data[i] == 0);
    }

    // Query the value
    auto values1 = index.query(Vec2(0.5, 0.5));
    REQUIRE(values1.size() == 1);
    REQUIRE(values1[0] == 0);

    // Query the value with a bounding box
    auto values2 = index.query(BoundingBox(Vec2(0.5, 0.5), Vec2(0.7, 0.7)));
    REQUIRE(values2.size() == 1);
    REQUIRE(values2[0] == 0);

}


TEST_CASE("Spatial index constructor works (two values)",
          "[geometry][2Dspatialindex]")
{
    index_t nx = 2;
    index_t ny = 2;

    vector<BoundingBox> bboxes = {
        BoundingBox(Vec2(0, 0), Vec2(0.75, 0.25)),
        BoundingBox(Vec2(0.6, 0), Vec2(1, 1))};

    auto bbox_func = [&](index_t i) { return bboxes[i]; };

    SpatialIndex2D index(nx, ny, bboxes.size(), bbox_func);

    REQUIRE(index.nx == nx);
    REQUIRE(index.ny == ny);

    // There are two values, so the bounding box should be the union of the two
    // values
    REQUIRE(index.bbox.min.x == 0);
    REQUIRE(index.bbox.min.y == 0);
    REQUIRE(index.bbox.max.x == 1);
    REQUIRE(index.bbox.max.y == 1);

    // Query the first cell
    auto values1 = index.query(Vec2(0.25, 0.25));
    REQUIRE(values1.size() == 1);
    REQUIRE(values1[0] == 0);

    // Query the second cell
    auto values2 = index.query(Vec2(0.25, 0.75));
    REQUIRE(values2.size() == 0);

    // Query the third cell
    auto values3 = index.query(Vec2(0.75, 0.25));
    REQUIRE(values3.size() == 2);
    REQUIRE(values3[0] == 0);
    REQUIRE(values3[1] == 1);

    // Query the fourth cell
    auto values4 = index.query(Vec2(0.75, 0.75));
    REQUIRE(values4.size() == 1);
    REQUIRE(values4[0] == 1);

}