#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include "neptune/common/utils.hpp"
#include "neptune/common/vector.hpp"

using namespace neptune;
TEST_CASE("FET: no drift works", "[FET]") {

    real_t D = 0.5;
    real_t v = 0;

    rng_ctr_t ctr{0, 0};
    rng_key_t key{42};
    rng_t rng(ctr, key);

    real_t alpha = 1;

    // bins for the histogram
    index_t n_bins = 100;
    real_t bin_width = 0.01;
    vector<real_t> bins(n_bins, 0);

    index_t n_samples = 10000;

    for (index_t i = 0; i < n_samples; i++) {
        real_t t = sample_first_exit_time(v, D, alpha, rng);
        index_t bin = t / bin_width;
        if(bin >= n_bins) {
            continue;
        }
        bins[bin]++;
    }

    // Normalize the histogram
    for (index_t i = 0; i < n_bins; i++) {
        bins[i] /= n_samples * bin_width;
    }

    // Now make a second histogram using the levy distribution directly
    vector<real_t> expected(n_bins, 0);
    real_t c = alpha * alpha / (2 * D);
    expected[0] = 0;
    for(index_t i = 1; i < n_bins; i++) {
        real_t t = bin_width * i;
        expected[i] = std::sqrt(c / (2 * M_PI)) * std::exp(-c / (2 * t)) / std::pow(t, 1.5);
    }

    // Check that the histograms are similar
    real_t max_error = 0;
    for (index_t i = 0; i < n_bins; i++) {
        max_error = std::max(max_error, std::abs(bins[i] - expected[i]));
    }
    REQUIRE(max_error < 0.2);
}

TEST_CASE("FET: positive drift works", "[FET]") {

    real_t D = 0.5;
    real_t v = 0.1;

    rng_ctr_t ctr{0, 0};
    rng_key_t key{42};
    rng_t rng(ctr, key);

    real_t alpha = 1;

    // bins for the histogram
    index_t n_bins = 100;
    real_t bin_width = 0.01;
    vector<real_t> bins(n_bins, 0);

    index_t n_samples = 10000;

    for (index_t i = 0; i < n_samples; i++) {
        real_t t = sample_first_exit_time(v, D, alpha, rng);
        index_t bin = t / bin_width;
        if(bin >= n_bins) {
            continue;
        }
        bins[bin]++;
    }

    // Normalize the histogram
    for (index_t i = 0; i < n_bins; i++) {
        bins[i] /= n_samples * bin_width;
    }

    // Now make a second histogram using the inverse gaussian distribution directly
    vector<real_t> expected(n_bins, 0);
    real_t mu = alpha / v;
    real_t lambda = alpha * alpha / (2 * D);
    expected[0] = 0;
    for(index_t i = 1; i < n_bins; i++) {
        real_t t = bin_width * i;
        expected[i] = std::sqrt(lambda / (2 * M_PI * t * t * t)) * std::exp(-lambda * (t - mu) * (t - mu) / (2 * mu * mu * t));
    }

    // Check that the histograms are similar
    real_t max_error = 0;
    for (index_t i = 0; i < n_bins; i++) {
        max_error = std::max(max_error, std::abs(bins[i] - expected[i]));
    }
    REQUIRE(max_error < 0.2);
    // Print the histogram
    for (index_t i = 0; i < n_bins; i++) {
        std::cout << "bin " << i*bin_width << ": " << bins[i] << " " << expected[i] << std::endl;
    }
}