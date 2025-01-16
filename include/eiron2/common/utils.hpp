#ifndef EIRON2_UTILS_HPP
#define EIRON2_UTILS_HPP

#include <cassert>
#include <random>

namespace eiron2 {

    // The floating point type used throughout the code
    typedef double real_t;

    // The index type used throughout the code
    typedef int64_t index_t;

    // The random number generator used throughout the code
    typedef std::mt19937 rng_t;

    inline real_t square(real_t x){
        return x * x;
    }

}  // namespace eiron2

#endif