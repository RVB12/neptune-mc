#ifndef NEPTUNE_NONUNIFORMMESH_HPP
#define NEPTUNE_NONUNIFORMMESH_HPP

#include "neptune/common/utils.hpp"
#include "neptune/common/vector.hpp"
#include "neptune/geometry/vec2.hpp"
#include "neptune/geometry/cellcrossing.hpp"
#include "neptune/geometry/mesh.hpp"
#include "neptune/geometry/particlegeometry.hpp"
#include "neptune/geometry/index3D.hpp"

namespace neptune {

/**
 * @brief A structured rectangular mesh consisting of rectangular voxels
 *
 * The mesh is of size (Lx x Ly x Lz) and is split into (nx x ny x nz) cells.
 * The cells are ordered sequentially in the x, then the y and then the z
 * direction.
 *
 * Unlike the VoxelMesh, the NonUniformMesh allows for varying cell sizes.
 *
 */
class NonUniformMesh : public Mesh {
   public:
    NonUniformMesh() = default;
    NonUniformMesh(vector<real_t>& cb_x, vector<real_t>& cb_y, vector<real_t>& cb_z)
        : Mesh(), cb_x(cb_x), cb_y(cb_y), cb_z(cb_z)
    {}

    NonUniformMesh(vector<real_t>&& cb_x, vector<real_t>&& cb_y, vector<real_t>&& cb_z)
        : Mesh(), cb_x(std::move(cb_x)), cb_y(std::move(cb_y)), cb_z(std::move(cb_z))
    {}

    NonUniformMesh(NonUniformMesh&& other) noexcept
        : Mesh(std::move(other)), cb_x(std::move(other.cb_x)), cb_y(std::move(other.cb_y)), cb_z(std::move(other.cb_z))
    {}

    NonUniformMesh(const NonUniformMesh& other) = default;

    NonUniformMesh& operator=(NonUniformMesh&& other) noexcept {
        if (this != &other) {
            Mesh::operator=(std::move(other));
            cb_x = std::move(other.cb_x);
            cb_y = std::move(other.cb_y);
            cb_z = std::move(other.cb_z);
        }
        return *this;
    }

    index_t n_cells() const;

    CellCrossing get_next_crossing(ParticleGeometry& p) const;

    /**
     * @brief Convert a linear index to triplet index
     */
    Index3D linear_to_tuple(index_t index) const;

    /**
     * @brief Convert a triplet index to linear index
     */
    index_t tuple_to_linear(Index3D ij) const;

    /**
     * @brief Return the id of the cell that a given position is in.
     */
    index_t get_cell_index(Vec3 position) const;

    /**
     * @brief Return the size of the cell in each direction
     */
    Vec3 cell_size() const;

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

    /**
     * @brief Return the number of surfaces of a given cell
     */
    index_t get_n_surfaces(index_t cell_index) const;

    /**
     * @brief Return the index of the neighbour cell that shares a given surface
     */
    index_t get_neighbour(index_t cell_index, index_t surface_id) const;

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
     * @brief Return the first boundary crossing (i.e. an actual boundary, not a cell boundary)
     * That a particle will cross.
     */
    std::pair<index_t, index_t> get_boundary_crossing(Vec3 position, Vec3 direction) const;

    /**
     * @brief Generate a random position within a given cell
     */
    Vec3 generate_random_position(index_t cell_index, rng_t& rng) const;

    /**
     * @brief Generate a random position on a given surface of a cell
     */
    Vec3 generate_random_surface_position(index_t cell_index, index_t surface_id, rng_t& rng) const;

    /**
     * @brief Transform a direction vector to take the mesh and the local position of a particle into account
     * 
     * (This is useful for meshes that represent curved geometries, such as rotated triangle meshes)
     * 
     * @param position The current position of the particle
     * @param direction The direction to be transformed
     * 
     * @return Vec3 The transformed direction
     */
    Vec3 transform_direction(Vec3 position, Vec3 direction) const;

    // Member variables
    // Note, these vectors must be sorted!
    vector<real_t> cb_x;  // x-coordinates of the cell boundaries
    vector<real_t> cb_y;  // y-coordinates of the cell boundaries
    vector<real_t> cb_z;  // z-coordinates of the cell boundaries
};

}  // namespace neptune

#endif