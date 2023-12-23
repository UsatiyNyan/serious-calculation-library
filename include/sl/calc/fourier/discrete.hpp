//
// Created by usatiynyan on 12/22/23.
//

#pragma once

#include <complex>
#include <span>
#include <type_traits>
#include <vector>

#include "sl/calc/fourier/detail.hpp"

namespace sl::calc::fourier {

template <direction direction_, typename FloatT, std::size_t extent_>
    requires std::is_floating_point_v<FloatT>
std::vector<std::complex<FloatT>> dft(std::span<const std::complex<FloatT>, extent_> in) {
    const std::size_t N = in.size();
    std::vector<std::complex<FloatT>> out(N);

    for (std::size_t k = 0; k != N; ++k) {
        for (std::size_t n = 0; n != N; ++n) {
            // $$X_k = \sum_{n=0}^{N-1} x_n \cdot e^{-\frac{2\pi i}{N}kn}$$
            out[k] += in[n] * detail::polar(detail::theta<direction_, FloatT>(k * n, N));
        }
        if constexpr (direction_ == direction::freq_to_time) {
            out[k] /= static_cast<FloatT>(N);
        }
    }

    return out;
}

} // namespace sl::calc::fourier
