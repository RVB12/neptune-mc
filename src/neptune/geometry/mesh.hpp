#ifndef NEPTUNE_MESH_HPP
#define NEPTUNE_MESH_HPP

#include "neptune/common/utils.hpp"
#include "neptune/common/vector.hpp"
#include "neptune/geometry/cellcrossing.hpp"
#include "neptune/geometry/particlegeometry.hpp"
#include "neptune/geometry/vec2.hpp"
#include "neptune/rng/rng.hpp"

namespace neptune {

#define BOUNDARY_TOLERANCE 1e-12

/**
 * @brief A custom exception for when the next crossing cannot be found. This is
 * usually due to a particle being outside of the mesh or having a velocity that
 * does not point into any cell.
 */
class NoNextCrossingException : public std::exception {
   public:
    const char* what() const noexcept override
    {
        return "No next crossing could be found. This is usually due to a "
               "particle being outside of the mesh or having a velocity that "
               "does not point into any cell.";
    }
};

/**
 * @brief This is a very general class for a mesh
 *
 */
class Mesh {
   public:
    Mesh() = default;
    ~Mesh() = default;

    /**
     * @brief Return the number of cells in the mesh
     */
    virtual index_t n_cells() const = 0;

    /**
     * For a given particle, if it continues at the current velocity,
     * find the position of the crossing with a cell boundary
     *
     * @param p
     */
    virtual CellCrossing get_next_crossing(ParticleGeometry& p) const = 0;

    /**
     * Check if a given position is inside of the cell with the given index
     */
    virtual bool is_in_cell(Vec3 p, index_t cell_index) const = 0;

    /**
     * Find the cell index corresponding to a given position.
     * Note: this may be very slow for unstructured meshes.
     *
     * returns -1 if no cell was found.
     *
     * @deprecated You should use find_particle instead,
     * as it is more robust
     */
    virtual index_t get_cell_index(Vec3 position) const = 0;

    /**
     * @brief Modify the particle geometry to set the cell_id and boundary_id of
     * the particle. If a particle is on a boundary (or very close to a
     * boundary), the particle may be considered to be in either of the two
     * cells that share the boundary. In this case, the velocity of the particle
     * is used to determine which cell the particle is considered to be in.
     * (velocity pointing into the cell)
     */
    virtual void find_particle(ParticleGeometry& p) const = 0;

    /**
     * @brief Return the volume of a given cell
     */
    virtual real_t get_volume(index_t cell_index) const = 0;

    /**
     * @brief Return the normal of a given surface
     * This normal should point into the cell, i.e. in the direction of the cell
     * center.
     *
     * @param cell_index The index of the cell
     * @param surface_id The id of the surface within the cell
     * @param position The position on the surface where the normal is requested
     *                 (this is needed for curved surfaces)
     */
    virtual Vec3 get_surface_normal(index_t cell_index,
                                    index_t surface_id,
                                    Vec3 position) const = 0;

    /**
     * @brief Return the number of surfaces of a given cell
     */
    virtual index_t get_n_surfaces(index_t cell_index) const = 0;

    /**
     * @brief Return the index of the neighbour cell that shares a given surface
     */
    virtual index_t get_neighbour(index_t cell_index,
                                  index_t surface_id) const = 0;

    /**
     * @brief Return the boundary condition of a given surface
     * -1 if there is no boundary condition
     */
    virtual index_t get_boundary_condition(index_t cell_index,
                                           index_t surface_id) const = 0;

    /**
     * @brief Return the surface area of a given surface of a cell
     */
    virtual real_t get_surface_area(index_t cell_index,
                                    index_t surface_id) const = 0;

    /**
     * @brief Return the center of a given cell
     */
    virtual Vec3 get_center(index_t cell_index) const = 0;

    /**
     * @brief Return the center of a given surface of a cell
     */
    virtual Vec3 get_surface_center(index_t cell_index,
                                    index_t surface_id) const = 0;

    /**
     * @brief Return the first boundary crossing (i.e. an actual boundary, not a
     * cell boundary) That a particle will cross.
     */
    virtual std::pair<index_t, index_t> get_boundary_crossing(
        Vec3 position, Vec3 direction) const = 0;

    /**
     * @brief Generate a random position within a given cell
     */
    virtual Vec3 generate_random_position(index_t cell_index,
                                          rng_t& rng) const = 0;

    /**
     * @brief Generate a random position on a given surface of a cell
     */
    virtual Vec3 generate_random_surface_position(index_t cell_index,
                                                  index_t surface_id,
                                                  rng_t& rng) const = 0;

    /**
     * @brief Return whether the mesh uses toroidal coordinates or carthesian
     * coordinates
     */
    virtual bool is_toroidal() const { return false; }
};

}  // namespace neptune

#endif