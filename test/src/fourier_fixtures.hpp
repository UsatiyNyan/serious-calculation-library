//
// Created by usatiynyan on 12/25/23.
//

#pragma once

#include <complex>
#include <fstream>
#include <numbers>
#include <type_traits>
#include <vector>

#include "sl/calc/fourier/detail.hpp"
#include <fmt/format.h>

template <typename T, typename Char>
struct fmt::formatter<std::complex<T>, Char> : public fmt::formatter<T, Char> {
    template <typename FormatCtx>
    auto format(const std::complex<T>& x, FormatCtx& ctx) const -> decltype(ctx.out()) {
        return format_to(ctx.out(), "[{}, {}]", x.real(), x.imag());
    }
};

namespace sl::calc::fourier {

template <typename F, typename T>
concept is_sample_producer = requires(F f, T x) {
    { f(x) } -> std::same_as<std::complex<T>>;
};

template <typename FloatT, is_sample_producer<FloatT> sample_producer>
    requires std::is_floating_point_v<FloatT>
std::vector<std::complex<FloatT>> produce_wave_samples(sample_producer f, std::size_t N) {
    std::vector<std::complex<FloatT>> wave_samples(N);
    for (std::size_t k = 0; k < N; ++k) {
        wave_samples[k] = f(2 * std::numbers::pi_v<FloatT> * static_cast<FloatT>(k) / static_cast<FloatT>(N));
    }
    return wave_samples;
}

template <typename FloatT, std::size_t extent_in_, std::size_t extent_out_>
void write_test_data(
    std::string_view name,
    std::span<const std::complex<FloatT>, extent_in_> in,
    std::span<const std::complex<FloatT>, extent_out_> out
) {
    std::ofstream{ fmt::format("{}.json", name) }
        << fmt::format(R"({{ "in": [{}], "out": [{}] }})", fmt::join(in, ", "), fmt::join(out, ", "));
}

template <typename FloatT>
auto normalize(auto out) {
    const std::size_t N = out.size();
    for (auto& x: out) {
        x /= static_cast<FloatT>(N);
    }
    return out;
}

} // namespace sl::calc::fourier
