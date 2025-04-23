#ifndef NEPTUNE_COLLISION_PROCESSES_HPP
#define NEPTUNE_COLLISION_PROCESSES_HPP

#include "neptune/common/utils.hpp"
#include "neptune/common/vector.hpp"

namespace neptune {

class BGKScatteringProcess {
   public:
    BGKScatteringProcess() = default;

    BGKScatteringProcess(index_t n_cells)
        : rate(n_cells),
          bgk_stddev(n_cells),
          bgk_mean_x(n_cells),
          bgk_mean_y(n_cells),
          bgk_mean_z(n_cells)
    {}

    // The collision rate per cell
    vector<real_t> rate;
    // The stddev of the post-collisional velocity per cell
    vector<real_t> bgk_stddev;
    // The mean of the post-collisional velocity in x direction per cell
    vector<real_t> bgk_mean_x;
    // The mean of the post-collisional velocity in y direction per cell
    vector<real_t> bgk_mean_y;
    // The mean of the post-collisional velocity in z direction per cell
    vector<real_t> bgk_mean_z;
};

class IonisationProcess {
   public:
    IonisationProcess() = default;

    IonisationProcess(index_t n_cells) : rate(n_cells) {}

    // The collision rate per cell
    vector<real_t> rate;
};

}  // namespace neptune

#endif