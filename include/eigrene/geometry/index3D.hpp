#ifndef EI_GRENE_GEOMETRY_INDEX3D_HPP
#define EI_GRENE_GEOMETRY_INDEX3D_HPP

#include "eigrene/common/utils.hpp"

namespace eigrene {

class Index3D {
   public:
    Index3D() {};
    Index3D(index_t i_, index_t j_, index_t k_) : i(i_), j(j_), k(k_) {}

    inline index_t& operator[](index_t ind)
    {
        assert(ind >= 0);
        assert(ind < 3);
        switch (ind) {
            case 0:
                return i;
            case 1:
                return j;
            case 2:
                return k;
        }
    }

    // Member variables
    index_t i = -1;
    index_t j = -1;
    index_t k = -1;
};

}  // namespace eigrene

#endif