#ifndef EIRON2_SPATIALINDEX_HPP
#define EIRON2_SPATIALINDEX_HPP

#include <cassert>
#include <functional>
#include <set>

#include "eiron2/common/utils.hpp"
#include "eiron2/common/vec2.hpp"
#include "eiron2/common/vector.hpp"

namespace eiron2 {

struct BoundingBox {
    BoundingBox() : min(Vec2(0, 0)), max(Vec2(0, 0)) {}

    BoundingBox(Vec2 min_, Vec2 max_) : min(min_), max(max_) {}

    inline BoundingBox merge(BoundingBox other)
    {
        return BoundingBox(
            Vec2(std::min(min.x, other.min.x), std::min(min.y, other.min.y)),
            Vec2(std::max(max.x, other.max.x), std::max(max.y, other.max.y)));
    }

    Vec2 min;
    Vec2 max;
};

/**
 * @brief A very simple grid-based spatial index
 *
 * The class stores a list of values sorted according to their position in a 2D
 * grid. The grid is defined by a minimum and maximum corner and the number of
 * cells in the x and y direction. You can easily query which values are in a
 * cell because a seperate list also keeps track of the start and end index of
 * each cell.
 *
 */
class SpatialIndex2D {
   public:
    /**
     * @brief Construct a new Spatial Index from a list of values
     *
     * @param nx_ The number of cells in the x direction
     * @param ny_  The number of cells in the y direction
     * @param n_values The number of values to be indexed
     * @param bbox_func A function that returns the bounding box of a value
     *
     */
    SpatialIndex2D(index_t nx_,
                   index_t ny_,
                   index_t n_values,
                   std::function<BoundingBox(index_t)> bbox_func)
        : nx(nx_), ny(ny_)
    {
        // Make three passes over the data.
        // First, determine the bounding box of the data
        // Second, determine how many values are in each cell
        // Third, allocate data array and fill it with the values

        // First pass, bounding box
        bbox = bbox_func(0);
        for (index_t i = 1; i < n_values; ++i) {
            bbox = bbox.merge(bbox_func(i));
        }

        cell_size = (bbox.max - bbox.min) / Vec2(nx, ny);

        // Second pass, count values in each cell
        vector<index_t> cell_count(nx * ny, 0);

        for (index_t i = 0; i < n_values; ++i) {
            BoundingBox box = bbox_func(i);
            index_t x0 = std::max<index_t>(
                0, std::min<index_t>(nx - 1,
                                     (box.min.x - bbox.min.x) / cell_size.x));
            index_t x1 = std::max<index_t>(
                0, std::min<index_t>(nx - 1,
                                     (box.max.x - bbox.min.x) / cell_size.x));
            index_t y0 = std::max<index_t>(
                0, std::min<index_t>(ny - 1,
                                     (box.min.y - bbox.min.y) / cell_size.y));
            index_t y1 = std::max<index_t>(
                0, std::min<index_t>(ny - 1,
                                     (box.max.y - bbox.min.y) / cell_size.y));

            for (index_t y = y0; y <= y1; ++y) {
                for (index_t x = x0; x <= x1; ++x) {
                    cell_count[x + y * nx]++;
                }
            }
        }

        // Construct the cell_idx array
        cell_idx.resize(nx * ny + 1);
        cell_idx[0] = 0;
        for (index_t i = 0; i < nx * ny; ++i) {
            cell_idx[i + 1] = cell_idx[i] + cell_count[i];
        }

        // The final value in cell_idx is also the full size of the data array
        // so we can allocate it now
        data.resize(cell_idx[nx * ny]);

        // Reset the cell_count array (we'll use it to keep track of where to
        // put the next value)
        for (index_t i = 0; i < nx * ny; ++i) {
            cell_count[i] = 0;
        }

        // Third pass, fill the data array
        for (index_t i = 0; i < n_values; ++i) {
            BoundingBox box = bbox_func(i);
            index_t x0 = std::max<index_t>(
                0, std::min<index_t>(nx - 1,
                                     (box.min.x - bbox.min.x) / cell_size.x));
            index_t x1 = std::max<index_t>(
                0, std::min<index_t>(nx - 1,
                                     (box.max.x - bbox.min.x) / cell_size.x));
            index_t y0 = std::max<index_t>(
                0, std::min<index_t>(ny - 1,
                                     (box.min.y - bbox.min.y) / cell_size.y));
            index_t y1 = std::max<index_t>(
                0, std::min<index_t>(ny - 1,
                                     (box.max.y - bbox.min.y) / cell_size.y));

            for (index_t y = y0; y <= y1; ++y) {
                for (index_t x = x0; x <= x1; ++x) {
                    index_t idx = cell_idx[x + y * nx] + cell_count[x + y * nx];
                    data[idx] = i;
                    cell_count[x + y * nx]++;
                }
            }
        }
    }

    // Copy constructor
    SpatialIndex2D(const SpatialIndex2D& other)
        : bbox(other.bbox),
          nx(other.nx),
          ny(other.ny),
          cell_size(other.cell_size),
          cell_idx(other.cell_idx),
          data(other.data)
    {
    }

    // Move constructor
    SpatialIndex2D(SpatialIndex2D&& other)
        : bbox(other.bbox),
          nx(other.nx),
          ny(other.ny),
          cell_size(other.cell_size),
          cell_idx(std::move(other.cell_idx)),
          data(std::move(other.data))
    {
    }

    /**
     * @brief Query the values in the same cell as a given position
     *
     * @param pos The position to query
     * 
     * @return vectorview<index_t> The values in the cell
     */
    inline vectorview<index_t> query(Vec2 position) const
    {
        assert(position.x >= bbox.min.x);
        assert(position.y >= bbox.min.y);
        assert(position.x <= bbox.max.x);
        assert(position.y <= bbox.max.y);

        index_t ix = std::max<index_t>(0, std::min<index_t>(nx - 1, (position.x - bbox.min.x) / cell_size.x));
        index_t iy = std::max<index_t>(0, std::min<index_t>(ny - 1, (position.y - bbox.min.y) / cell_size.y));

        index_t start = cell_idx[ix + iy * nx];
        index_t end = cell_idx[ix + iy * nx + 1];

        return data.subvector(start, end);
    }

    /**
     * @brief Query the values that are in the same cells as a given bounding box
     * 
     * @param box The bounding box to query
     */
    inline vector<index_t> query(BoundingBox query_box) const
    {
        assert(query_box.min.x >= bbox.min.x);
        assert(query_box.min.y >= bbox.min.y);
        assert(query_box.max.x <= bbox.max.x);
        assert(query_box.max.y <= bbox.max.y);

        index_t x0 = std::max<index_t>(0, std::min<index_t>(nx - 1, (query_box.min.x - bbox.min.x) / cell_size.x));
        index_t x1 = std::max<index_t>(0, std::min<index_t>(nx - 1, (query_box.max.x - bbox.min.x) / cell_size.x));
        index_t y0 = std::max<index_t>(0, std::min<index_t>(ny - 1, (query_box.min.y - bbox.min.y) / cell_size.y));
        index_t y1 = std::max<index_t>(0, std::min<index_t>(ny - 1, (query_box.max.y - bbox.min.y) / cell_size.y));

        std::set<index_t> result;
        for (index_t y = y0; y <= y1; ++y) {
            for (index_t x = x0; x <= x1; ++x) {
                index_t start = cell_idx[x + y * nx];
                index_t end = cell_idx[x + y * nx + 1];
                for (index_t i = start; i < end; ++i) {
                    result.insert(data[i]);
                }
            }
        }

        vector<index_t> result_vector(result.size());
        std::copy(result.begin(), result.end(), result_vector.begin());

        return result_vector;
    }

    // The bounding box of the spatial index
    BoundingBox bbox;
    // The number of cells in the x direction
    index_t nx;
    // The number of cells in the y direction
    index_t ny;
    // The size of one cell in the x and y direction
    Vec2 cell_size;
    // Column-major list of the start index of each cell
    // in the data array. If cell_idx[i + j * nx] = k
    // and cell_idx[i + 1 + j * nx] = l, then the values
    // in the cell (i, j) are stored in data[k:l].
    vector<index_t> cell_idx;
    // The list of the actual values indexed by this spatial index
    // We only store integers, these are meant to be indices into
    // a different (external) array.
    vector<index_t> data;
};

}  // namespace eiron2

#endif