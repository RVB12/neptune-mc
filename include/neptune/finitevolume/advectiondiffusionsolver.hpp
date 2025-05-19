#ifndef NEPTUNE_FINITE_VOLUME_SOLVER_HPP
#define NEPTUNE_FINITE_VOLUME_SOLVER_HPP

#include "neptune/common/utils.hpp"
#include "neptune/common/vector.hpp"
#include "neptune/geometry/mesh.hpp"
#include "neptune/geometry/trianglemesh.hpp"
#include "neptune/geometry/voxelmesh.hpp"
#include "neptune/physics/plasmabackground.hpp"
#include "neptune/physics/physicalconstants.hpp"

namespace neptune{

    /**
     * Solve the advection-diffusion equation using a finite volume method
     * 
     * @param mesh The mesh to use
     * @param advection_x The x component of the advection velocity in each cell (m/s)
     * @param advection_y The y component of the advection velocity in each cell (m/s)
     * @param advection_z The z component of the advection velocity in each cell (m/s)
     * @param diffusion The diffusion coefficient in each cell (m^2/s)
     * @param t_final The final time to simulate (s)
     * @param dt The time step to use (s)
     * @param density On input, the particle density in each cell (m^-3)
     *                On output, the updated density in each cell (m^-3)
     */
    void advection_diffusion_solver(
        const Mesh& mesh,
        const vector<real_t>& advection_x,
        const vector<real_t>& advection_y,
        const vector<real_t>& advection_z,
        const vector<real_t>& diffusion,
        real_t t_final,
        real_t dt,
        vector<real_t>& density);

}

#endif