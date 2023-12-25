//
// Created by usatiynyan on 12/22/23.
//

#pragma once

#include <complex>
#include <span>
#include <type_traits>
#include <vector>

#include <assert.hpp>

#include "sl/calc/fourier/detail.hpp"
#include "sl/calc/fourier/discrete.hpp"

namespace sl::calc::fourier {

template <direction direction_, typename FloatT, std::size_t extent_>
    requires std::is_floating_point_v<FloatT> && detail::extent_is_power_of_2<extent_>
std::vector<std::complex<FloatT>>
    fft_recursive(std::span<const std::complex<FloatT>, extent_> in, std::size_t offset = 0, std::size_t stride = 1) {
    const std::size_t N = in.size() / stride;
    ASSERT(std::has_single_bit(N), "only accepting powers of 2");

    if (N == 1) {
        return { in[offset] };
    }

    const auto even_out = fft_recursive<direction_>(in, offset, stride * 2);
    const auto odd_out = fft_recursive<direction_>(in, offset + stride, stride * 2);

    std::vector<std::complex<FloatT>> out(N);

    for (std::size_t k = 0; k < N / 2; ++k) {
        // $$ e^{-i 2 \pi \frac{k}{N}} $$
        const auto twiddle_factor = detail::polar(detail::theta<direction_, FloatT>(k, N));
        // $$ e^{-i 2 \pi \frac{k}{N}} O_k $$
        const auto twiddle_factor_x_odd = twiddle_factor * odd_out[k];
        // $$ X_k         = E_k + e^{-i 2 \pi \frac{k}{N}} O_k $$
        out[k /*   */] = even_out[k] + twiddle_factor_x_odd;
        // $$ X_{k + N/2} = E_k - e^{-i 2 \pi \frac{k}{N}} O_k $$
        out[k + N / 2] = even_out[k] - twiddle_factor_x_odd;
    }

    return out;
}

// TODO: wrap impl
// TODO: for inverse transform divide output by N

} // namespace sl::calc::fourier
