#ifndef NEPTUNE_RotatedTriangleMesh_HPP
#define NEPTUNE_RotatedTriangleMesh_HPP

#include "neptune/common/utils.hpp"
#include "neptune/common/vector.hpp"
#include "neptune/geometry/cellcrossing.hpp"
#include "neptune/geometry/mesh.hpp"
#include "neptune/geometry/particlegeometry.hpp"
#include "neptune/geometry/spatialindex.hpp"
#include "neptune/geometry/triangle.hpp"
#include "neptune/geometry/trianglemesh.hpp"
#include "neptune/geometry/vec2.hpp"

namespace neptune {

/**
 * @brief An unstructured mesh consisting of triangles that are rotated around
 * the z-axis to form 3D cells.
 *
 * The is an implicit assumption that the x-coordinate of the vertices (the
 * radial coordinate) is non-negative.
 *
 */
class RotatedTriangleMesh : public Mesh {
   public:
    RotatedTriangleMesh(vector<Vec2>& vertices_,
                        vector<TriangleIndices>& triangles_)
        : Mesh(),
          vertices(vertices_),
          triangles(triangles_),
          spatialindex(vertices, triangles)
    {
        make_neighborlist(vertices, triangles, spatialindex, neighbors);
    }

    RotatedTriangleMesh(vector<Vec2>& vertices_,
                        vector<TriangleIndices>& triangles_,
                        vector<index_t>& neighbors_)
        : Mesh(),
          vertices(vertices_),
          triangles(triangles_),
          spatialindex(vertices, triangles),
          neighbors(neighbors_)
    {
    }

    RotatedTriangleMesh(vector<Vec2>& vertices_,
                        vector<TriangleIndices>& triangles_,
                        SpatialIndex& spatialindex,
                        vector<index_t>& neighbors_)
        : Mesh(),
          vertices(vertices_),
          triangles(triangles_),
          spatialindex(spatialindex),
          neighbors(neighbors_)
    {}

    RotatedTriangleMesh(vector<Vec2>&& vertices_,
                        vector<TriangleIndices>&& triangles_,
                        SpatialIndex&& spatialindex,
                        vector<index_t>&& neighbors_)
        : Mesh(),
          vertices(std::move(vertices_)),
          triangles(std::move(triangles_)),
          spatialindex(std::move(spatialindex)),
          neighbors(std::move(neighbors_))
    {}

    RotatedTriangleMesh(const RotatedTriangleMesh& other)
        : Mesh(),
          vertices(other.vertices),
          triangles(other.triangles),
          spatialindex(other.spatialindex),
          neighbors(other.neighbors)
    {}

    RotatedTriangleMesh(RotatedTriangleMesh&& other)
        : Mesh(),
          vertices(std::move(other.vertices)),
          triangles(std::move(other.triangles)),
          spatialindex(std::move(other.spatialindex)),
          neighbors(std::move(other.neighbors))
    {}

    index_t n_cells() const;

    CellCrossing get_next_crossing(ParticleGeometry& p) const;

    /**
     * @brief Return the id of the cell that a given position is in.
     */
    index_t get_cell_index(Vec3 position) const;

    /**
     * @brief Return the triangle indices corresponding to the given index
     */
    TriangleIndices get_triangle_indices(index_t cell_index) const;

    /**
     * @brief Return the triangle corresponding to the given index
     */
    Triangle get_triangle(index_t cell_index) const;

    /**
     * Check if a given position is inside of the cell with the given index
     */
    bool is_in_cell(Vec3 p, index_t cell_index) const;

    /**
     * @brief Return the volume of a given cell
     */
    real_t get_volume(index_t cell_index) const;

    /**
     * @brief Return the normal of a given surface
     */
    Vec3 get_surface_normal(index_t cell_index,
                            index_t surface_id,
                            Vec3 position) const;

    /**
     * @brief Return the index of the neighbour cell that shares a given surface
     */
    index_t get_neighbour(index_t cell_index, index_t surface_id) const;

    /**
     * @brief Return the number of surfaces of a given cell
     */
    index_t get_n_surfaces(index_t cell_index) const;

    /**
     * @brief Return the area of a given surface
     */
    real_t get_surface_area(index_t cell_index, index_t surface_id) const;

    /**
     * @brief Return the center of a given cell
     */
    Vec3 get_center(index_t cell_index) const;

    /**
     * @brief Return the center of a given surface of a cell
     */
    Vec3 get_surface_center(index_t cell_index, index_t surface_id) const;

    /**
     * @brief Return the first boundary crossing (i.e. an actual boundary, not a
     * cell boundary) That a particle will cross.
     */
    std::pair<index_t, index_t> get_boundary_crossing(Vec3 position,
                                                      Vec3 direction) const;

    /**
     * @brief Generate a random position within a given cell
     */
    Vec3 generate_random_position(index_t cell_index, rng_t& rng) const;

    /**
     * @brief Generate a random position on a given surface of a cell
     */
    Vec3 generate_random_surface_position(index_t cell_index,
                                          index_t surface_id,
                                          rng_t& rng) const;

    vector<index_t> get_neighbor_list() const { return neighbors; }

    /**
     * @brief Return whether the mesh uses toroidal coordinates or carthesian
     * coordinates
     */
    bool is_toroidal() const { return true; }

   private:
    // Member variables

    // List of all the triangles in the mesh
    vector<Vec2> vertices;

    // List of all the triangles in the mesh
    vector<TriangleIndices> triangles;

    // Spatial index for fast spatial queries
    SpatialIndex spatialindex;

    // Map from triangle edge to the triangle that shares the edge
    // used for fast traversal
    vector<index_t> neighbors;
};

// Return the time until the intersection of particle with position pos and
// velocity vel with a cone defined by rotating points p1 and p2 around the
// z-axis. Note, this function does not verify whether the returned time is
// positive or whether the intersection lies between the two points. This should
// be done later.
real_t find_cone_intersection(const Vec2& p1,
                              const Vec2& p2,
                              const Vec3& pos,
                              const Vec3& vel,
                              bool is_on_surface);

}  // namespace neptune

#endif