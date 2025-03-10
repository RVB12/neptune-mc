#ifndef EIGRENE_TRIANGLEMESH_HPP
#define EIGRENE_TRIANGLEMESH_HPP

#include "eigrene/common/utils.hpp"
#include "eigrene/geometry/vec2.hpp"
#include "eigrene/common/vector.hpp"
#include "eigrene/geometry/cellcrossing.hpp"
#include "eigrene/geometry/mesh.hpp"
#include "eigrene/geometry/particlegeometry.hpp"
#include "eigrene/geometry/spatialindex.hpp"
#include "eigrene/geometry/triangle.hpp"

namespace eigrene {

/**
 * @brief Construct a neighbor list for a list of triangles
 *
 * If two triangles share an edge, the neighbor list will contain the index of
 * the neighboring triangle. If no neighbor is found, the corresponding entry
 * will be -1. If an edge has entry -1, it is a boundary edge (or the mesh is
 * ill defined)
 *
 * After this function, you should still fill change the boundary edges to the
 * correct boundary id.
 *
 * This version is pretty slow and dus not use any spatial index.
 *
 * @param triangles
 * @param neighbors
 */
void make_neighborlist(const vector<TriangleIndices>& triangles,
                       vector<index_t>& neighbors);
/**
 * @brief Construct a neighbor list for a list of triangles
 *
 * If two triangles share an edge, the neighbor list will contain the index of
 * the neighboring triangle. If no neighbor is found, the corresponding entry
 * will be -1. If an edge has entry -1, it is a boundary edge (or the mesh is
 * ill defined)
 *
 * After this function, you should still fill change the boundary edges to the
 * correct boundary id.
 *
 * This version uses a spatial index for efficiency.
 */
void make_neighborlist(const vector<Vec2>& vertices,
                       const vector<TriangleIndices>& triangles,
                       const SpatialIndex& spatialindex,
                       vector<index_t>& neighbors);

/**
 * @brief An unstructured mesh consisting of triangles
 *
 * The mesh is stored as a list of vertices and a list of triangles,
 * where each each triangle is defined by 3 indices, indexing into the vertex
 * list.
 *
 * To allow fast traversal, a neighbor list is also maintained.
 *
 */
class TriangleMesh : public Mesh {
   public:
    TriangleMesh(vector<Vec2>& vertices_, vector<TriangleIndices>& triangles_)
        : Mesh(),
          vertices(vertices_),
          triangles(triangles_),
          spatialindex(vertices, triangles)
    {
        make_neighborlist(vertices, triangles, spatialindex, neighbors);
    }

    TriangleMesh(vector<Vec2>& vertices_,
                 vector<TriangleIndices>& triangles_,
                 SpatialIndex& spatialindex,
                 vector<index_t>& neighbors_)
        : Mesh(),
          vertices(vertices_),
          triangles(triangles_),
          spatialindex(spatialindex),
          neighbors(neighbors_)
    {}

    TriangleMesh(vector<Vec2>&& vertices_,
                 vector<TriangleIndices>&& triangles_,
                 SpatialIndex&& spatialindex,
                 vector<index_t>&& neighbors_)
        : Mesh(),
          vertices(std::move(vertices_)),
          triangles(std::move(triangles_)),
          spatialindex(std::move(spatialindex)),
          neighbors(std::move(neighbors_))
    {}

    TriangleMesh(const TriangleMesh& other)
        : Mesh(),
          vertices(other.vertices),
          triangles(other.triangles),
          spatialindex(other.spatialindex),
          neighbors(other.neighbors)
    {}

    TriangleMesh(TriangleMesh&& other)
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
    Vec3 get_surface_normal(index_t cell_index, index_t surface_id) const;

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

}  // namespace eigrene

#endif