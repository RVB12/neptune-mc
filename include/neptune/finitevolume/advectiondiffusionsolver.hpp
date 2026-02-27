#ifndef NEPTUNE_FINITE_VOLUME_SOLVER_HPP
#define NEPTUNE_FINITE_VOLUME_SOLVER_HPP

#include "neptune/common/utils.hpp"
#include "neptune/common/vector.hpp"
#include "neptune/geometry/mesh.hpp"
#include "neptune/geometry/trianglemesh.hpp"
#include "neptune/geometry/voxelmesh.hpp"
#include "neptune/physics/plasmabackground.hpp"
#include "neptune/physics/physicalconstants.hpp"
#include "neptune/finitevolume/fluidboundarycondition.hpp"

namespace neptune{

    /**
     * Solve the advection-diffusion equation using a finite volume method
     * 
     * @param mesh The mesh to use
     * @param advection_x The x component of the advection velocity in each cell (m/s)
     * @param advection_y The y component of the advection velocity in each cell (m/s)
     * @param advection_z The z component of the advection velocity in each cell (m/s)
     * @param diffusion_rate The diffusion rate in each cell (this is typically bgk_stddev^2)
     * @param diffusion_factor The diffusion factor in each cell (this is typically 1/bgk_rate)
     * @param ionisation_rate The ionisation rate in each cell (1/s)
     * @param t_final The final time to simulate (s)
     * @param dt The time step to use (s)
     * @param density On input, the particle density in each cell (m^-3)
     *                On output, the updated density in each cell (m^-3)
     * @param density_integral The integral of the density over time in each cell (m^-3 s)
     */
    void advection_diffusion_solver(
        const Mesh& mesh,
        const vector<std::unique_ptr<FluidBoundaryCondition>>& boundary_conditions,
        const vector<real_t>& advection_x,
        const vector<real_t>& advection_y,
        const vector<real_t>& advection_z,
        const vector<real_t>& diffusion_rate,
        const vector<real_t>& diffusion_factor,
        const vector<real_t>& ionisation_rate,
        real_t t_final,
        real_t dt,
        vector<real_t>& density,
        vector<real_t>& density_integral,
        vector<real_t>& momentum_x_integral,
        vector<real_t>& momentum_y_integral,
        vector<real_t>& momentum_z_integral,
        vector<real_t>& energy_integral,
        vector<real_t>& energy_x_integral);

}

#endif