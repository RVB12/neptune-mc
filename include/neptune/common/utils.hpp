#ifndef NEPTUNE_UTILS_HPP
#define NEPTUNE_UTILS_HPP

#include <cassert>
#include <random>

namespace neptune {

// The floating point type used throughout the code
typedef double real_t;

// The index type used throughout the code
typedef int64_t index_t;

// The random number generator used throughout the code
typedef std::mt19937 rng_t;

inline real_t square(real_t x) { return x * x; }

/**
 * Sample the first time a brownian particle with initial position 0,
 * drift velocity v, and diffusion coefficient D exceeds alpha.
 */
inline real_t sample_first_exit_time(real_t v,
                                     real_t D,
                                     real_t alpha,
                                     rng_t& rng)
{
    std::uniform_real_distribution<real_t> r_unif(0, 1);
    std::normal_distribution<real_t> r_normal(0, 1);

    if (std::abs(v) < 1e-10) {
        // Sample from levy distribution with c = alpha * alpha / (2 * D) and mu
        // = 0 method from
        // https://stats.libretexts.org/Bookshelves/Probability_Theory/Probability_Mathematical_Statistics_and_Stochastic_Processes_(Siegrist)/05%3A_Special_Distributions/5.16%3A_The_Levy_Distribution
        real_t n = r_normal(rng);
        // x is now a sample from the standard levy distribution
        real_t x = 1. / (n * n);
        // scale x to get the scaled levy distribution
        return alpha * alpha / (2 * D) * x;
    }
    else {
        if (v > 0) {
            // Sample from the inverse gaussian distribution with mu = alpha / v
            // and lambda = alpha * alpha / (2 * D) method from
            // https://en.wikipedia.org/wiki/Inverse_Gaussian_distribution#Sampling_from_an_inverse-Gaussian_distribution
            real_t mu = alpha / v;
            real_t lambda = alpha * alpha / (2 * D);

            real_t nu = r_normal(rng);
            real_t y = nu * nu;
            real_t x = mu + (mu * mu * y) / (2 * lambda) -
                       (mu / (2 * lambda)) *
                           std::sqrt(4 * mu * lambda * y + mu * mu * y * y);
            real_t z = r_unif(rng);
            if (z * (mu + x) <= mu) {
                return x;
            }
            else {
                return mu * mu / x;
            }
        }
        else {
            // particle is drifting away from the boundary, the chance of
            // hitting the boundary at t = infinity is nonzero

            // TODO: correct this, we should sample sometimes
            return std::numeric_limits<real_t>::infinity();
        }
    }
}

}  // namespace neptune

#endif