#ifndef EIRON2_GEOMETRY_TRIANGLE_HPP
#define EIRON2_GEOMETRY_TRIANGLE_HPP

#include "eiron2/common/utils.hpp"

namespace eiron2 {

/**
 * @brief A triangle, defined by three indices into the vertex array
 *        The indices should be listed in counter clockwise order.
 */
class Triangle {
   public:
    index_t ind1;
    index_t ind2;
    index_t ind3;
};

}  // namespace eiron2

#endif