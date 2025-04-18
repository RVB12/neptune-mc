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
     * Solve the fluid model with as initial condition a dirac delta at the cell with the given id.
     */
    void solve_fluid_model(const VoxelMesh& mesh, const PlasmaBackground& background, index_t cell_id);

}

#endif