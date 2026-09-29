#ifndef NEPTUNE_UTILS_HPP
#define NEPTUNE_UTILS_HPP

#include <cassert>
#include <stdexcept>
#include <cstdint>

namespace neptune {

// The floating point type used throughout the code
typedef double real_t;

// The index type used throughout the code
typedef int64_t index_t;


inline real_t square(real_t x) { return x * x; }

}  // namespace neptune

#endif