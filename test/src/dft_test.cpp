//
// Created by usatiynyan on 12/23/23.
//

#include "sl/calc/fourier/discrete.hpp"

#include "fourier_fixtures.hpp"

#include <gtest/gtest.h>
#include <random>

namespace sl::calc::fourier {

constexpr std::size_t N = 64;

TEST(dft, sin) {
    const auto in = produce_wave_samples<double>([](double theta) { return std::complex{ std::sin(theta), 0.0 }; }, N);
    const auto out = normalize<double>(dft<direction::time_to_freq>(std::span{ in }));
    write_test_data("dft_sin", std::span{ in }, std::span{ out });
}

TEST(dft, cos) {
    const auto in = produce_wave_samples<double>([](double theta) { return std::complex{ std::cos(theta), 0.0 }; }, N);
    const auto out = normalize<double>(dft<direction::time_to_freq>(std::span{ in }));
    write_test_data("dft_cos", std::span{ in }, std::span{ out });
}

TEST(dft, identity) {
    const auto in = produce_wave_samples<double>([](double theta) { return std::polar(1.0, theta); }, N);
    const auto out = normalize<double>(dft<direction::time_to_freq>(std::span{ in }));
    write_test_data("dft_identity", std::span{ in }, std::span{ out });
}

TEST(dft, harmonicSum) {
    const auto in = produce_wave_samples<double>(
        [](double theta) {
            return std::complex{ std::sin(theta) + std::cos(2 * theta), 0.0 }
                   + std::complex{ std::sin(3 * theta), 0.0 };
        },
        N
    );
    const auto out = normalize<double>(dft<direction::time_to_freq>(std::span{ in }));
    write_test_data("dft_harmonic_sum", std::span{ in }, std::span{ out });
}

TEST(dft, random) {
    // Choose a random mean between 1 and 6
    std::default_random_engine re(std::random_device{}());
    std::uniform_real_distribution<double> uniform_dist(0.0, 2 * std::numbers::pi);
    const auto in =
        produce_wave_samples<double>([&uniform_dist, &re](double) { return std::polar(1.0, uniform_dist(re)); }, N);
    const auto out = normalize<double>(dft<direction::time_to_freq>(std::span{ in }));
    write_test_data("dft_random", std::span{ in }, std::span{ out });
}

} // namespace sl::calc::fourier
