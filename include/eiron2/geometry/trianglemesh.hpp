#ifndef EIRON2_TRIANGLEMESH_HPP
#define EIRON2_TRIANGLEMESH_HPP

#include "eiron2/geometry/cellcrossing.hpp"
#include "eiron2/geometry/mesh.hpp"
#include "eiron2/geometry/particlegeometry.hpp"
#include "eiron2/common/vec2.hpp"
#include "eiron2/common/vector.hpp"
#include "eiron2/common/utils.hpp"

namespace eiron2 {

/**
 * @brief A triangle, defined by three indices into the vertex array
 *        The indices should be listed in counter clockwise order.
 */
class Triangle{
    public:
    index_t ind1;
    index_t ind2;
    index_t ind3;
};

/**
 * @brief An unstructured mesh consisting of triangles
 * 
 * The mesh is stored as a list of vertices and a list of triangles,
 * where each each triangle is defined by 3 indices, indexing into the vertex list.
 * 
 * To allow fast traversal, a neighbor list is also maintained.
 *
 */
class TriangleMesh : public Mesh {
   public:

    TriangleMesh(vector<Vec2>& vertices_, vector<Triangle>& triangles_)
        : Mesh(),
         vertices(vertices_),
         triangles(triangles_),
         neighbors(3*triangles.size())
    {
        make_neighborlist();
    }

    TriangleMesh(vector<Vec2>& vertices_, vector<Triangle>& triangles_, vector<index_t>& neighbors_)
        : Mesh(),
         vertices(vertices_),
         triangles(triangles_),
         neighbors(neighbors_)
    {
        make_neighborlist();
    }

    TriangleMesh(vectorview<real_t> vertices_, vectorview<index_t> triangles_)
        : Mesh(),
         vertices(vertices_.size() / 2),
         triangles(triangles_.size() / 3),
         neighbors(3*triangles_.size())
    {
        for (index_t i = 0; i < vertices.size(); i++) {
            vertices[i] = Vec2(vertices_[2*i], vertices_[2*i+1]);
        }
        for(index_t i = 0; i < triangles.size(); i++) {
            triangles[i] = Triangle(triangles_[3*i], triangles_[3*i+1], triangles_[3*i+2]);
        }
        make_neighborlist();
    }

    index_t n_cells() const;

    CellCrossing get_next_crossing(ParticleGeometry& p) const;

    /**
     * @brief Return the id of the cell that a given position is in.
     */
    index_t get_cell_index(Vec2 position) const;

    /**
     * @brief Return the triangle corresponding to the given index
     */
    Triangle get_triangle(index_t cell_index) const;

    /**
     * Check if a given position is inside of the cell with the given index
     */
    bool is_in_cell(Vec2 p, index_t cell_index) const;

    /**
     * @brief Return the volume of a given cell
     */
    real_t get_volume(index_t cell_index) const;

    /**
     * @brief Return the normal of a given surface
     */
    Vec2 get_surface_normal(index_t cell_index, index_t surface_id) const;

   private:
    // Go through all the triangles and construct the neighborlist
    void make_neighborlist();

    // Member variables
    vector<Vec2> vertices; // List of all the triangles in the mesh
    vector<Triangle> triangles; // List of all the triangles in the mesh
    vector<index_t> neighbors; // Map of the shared edges of the triangles
};

}  // namespace eiron2

#endif