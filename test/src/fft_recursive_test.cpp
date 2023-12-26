//
// Created by usatiynyan on 12/24/23.
//

#include "sl/calc/fourier/discrete.hpp"
#include "sl/calc/fourier/fast.hpp"

#include "fourier_fixtures.hpp"

#include <gtest/gtest.h>
#include <random>

namespace sl::calc::fourier {

constexpr std::size_t N = 64;
constexpr double ERR = 1e-12;

TEST(fftRecursive, sin) {
    const auto in = produce_wave_samples<double>([](double theta) { return std::complex{ std::sin(theta), 0.0 }; }, N);
    const auto out = fft_recursive<direction::time_to_freq>(std::span{ in });
    const auto dft_out = dft<direction::time_to_freq>(std::span{ in });
    for (std::size_t k = 0; k < N; ++k) {
        EXPECT_NEAR(out[k].real(), dft_out[k].real(), ERR);
        EXPECT_NEAR(out[k].imag(), dft_out[k].imag(), ERR);
    }
    const auto normalized_out = normalize<double>(out);
    write_test_data("fft_recursive_sin", std::span{ in }, std::span{ normalized_out });
}

TEST(fftRecursive, cos) {
    const auto in = produce_wave_samples<double>([](double theta) { return std::complex{ std::cos(theta), 0.0 }; }, N);
    const auto out = fft_recursive<direction::time_to_freq>(std::span{ in });
    const auto dft_out = dft<direction::time_to_freq>(std::span{ in });
    for (std::size_t k = 0; k < N; ++k) {
        EXPECT_NEAR(out[k].real(), dft_out[k].real(), ERR);
        EXPECT_NEAR(out[k].imag(), dft_out[k].imag(), ERR);
    }
    const auto normalized_out = normalize<double>(out);
    write_test_data("fft_recursive_cos", std::span{ in }, std::span{ normalized_out });
}

TEST(fftRecursive, identity) {
    const auto in = produce_wave_samples<double>([](double theta) { return std::polar(1.0, theta); }, N);
    const auto out = fft_recursive<direction::time_to_freq>(std::span{ in });
    const auto dft_out = dft<direction::time_to_freq>(std::span{ in });
    for (std::size_t k = 0; k < N; ++k) {
        EXPECT_NEAR(out[k].real(), dft_out[k].real(), ERR);
        EXPECT_NEAR(out[k].imag(), dft_out[k].imag(), ERR);
    }
    const auto normalized_out = normalize<double>(out);
    write_test_data("fft_recursive_identity", std::span{ in }, std::span{ normalized_out });
}

TEST(fftRecursive, harmonicSum) {
    const auto in = produce_wave_samples<double>(
        [](double theta) {
            return std::complex{ std::sin(theta) + std::cos(2 * theta), 0.0 }
                   + std::complex{ std::sin(3 * theta), 0.0 };
        },
        N
    );
    const auto out = fft_recursive<direction::time_to_freq>(std::span{ in });
    const auto dft_out = dft<direction::time_to_freq>(std::span{ in });
    for (std::size_t k = 0; k < N; ++k) {
        EXPECT_NEAR(out[k].real(), dft_out[k].real(), ERR);
        EXPECT_NEAR(out[k].imag(), dft_out[k].imag(), ERR);
    }
    const auto normalized_out = normalize<double>(out);
    write_test_data("fft_recursive_harmonic_sum", std::span{ in }, std::span{ normalized_out });
}

TEST(fftRecursive, random) {
    // Choose a random mean between 1 and 6
    std::default_random_engine re(std::random_device{}());
    std::uniform_real_distribution<double> uniform_dist(0.0, 2 * std::numbers::pi);
    const auto in =
        produce_wave_samples<double>([&uniform_dist, &re](double) { return std::polar(1.0, uniform_dist(re)); }, N);
    const auto out = fft_recursive<direction::time_to_freq>(std::span{ in });
    const auto dft_out = dft<direction::time_to_freq>(std::span{ in });
    for (std::size_t k = 0; k < N; ++k) {
        EXPECT_NEAR(out[k].real(), dft_out[k].real(), ERR);
        EXPECT_NEAR(out[k].imag(), dft_out[k].imag(), ERR);
    }
    const auto normalized_out = normalize<double>(out);
    write_test_data("fft_recursive_random", std::span{ in }, std::span{ normalized_out });
}

} // namespace sl::calc::fourier
