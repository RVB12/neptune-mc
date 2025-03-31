#ifndef EIGRENE_FINITE_VOLUME_SOLVER_HPP
#define EIGRENE_FINITE_VOLUME_SOLVER_HPP

#include "eigrene/common/utils.hpp"
#include "eigrene/common/vector.hpp"
#include "eigrene/geometry/mesh.hpp"
#include "eigrene/geometry/trianglemesh.hpp"
#include "eigrene/geometry/voxelmesh.hpp"
#include "eigrene/physics/plasmabackground.hpp"
#include "eigrene/physics/physicalconstants.hpp"

namespace eigrene{

    /**
     * Solve the fluid model with as initial condition a dirac delta at the cell with the given id.
     */
    void solve_fluid_model(const VoxelMesh& mesh, const PlasmaBackground& background, index_t cell_id);

}

#endif