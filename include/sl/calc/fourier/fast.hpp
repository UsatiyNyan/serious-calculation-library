//
// Created by usatiynyan on 12/22/23.
//

#pragma once

#include <complex>
#include <type_traits>
#include <vector>

#include "sl/calc/bits.hpp"
#include "sl/calc/fourier/detail.hpp"

#include <libassert/assert.hpp>

namespace sl::calc::fourier {
namespace detail {

template <direction direction_, typename FloatT, std::size_t extent_>
std::vector<std::complex<FloatT>>
    fft_recursive_impl(std::span<const std::complex<FloatT>, extent_> in, std::size_t offset, std::size_t stride) {
    const std::size_t N = in.size() / stride;

    if (N == 1) {
        return { in[offset] };
    }

    const auto even_out = fft_recursive_impl<direction_>(in, offset, stride * 2);
    const auto odd_out = fft_recursive_impl<direction_>(in, offset + stride, stride * 2);

    std::vector<std::complex<FloatT>> out(N);

    for (std::size_t k = 0; k < N / 2; ++k) {
        // $$ e^{-i 2 \pi \frac{k}{N}} $$
        const auto twiddle_factor = detail::polar(detail::theta<direction_, FloatT>(k, N));
        // $$ e^{-i 2 \pi \frac{k}{N}} O_k $$
        const auto twiddle_factor_x_odd = twiddle_factor * odd_out[k];
        // $$ X_k         = E_k + e^{-i 2 \pi \frac{k}{N}} O_k $$
        out[k] /*   */ = even_out[k] + twiddle_factor_x_odd;
        // $$ X_{k + N/2} = E_k - e^{-i 2 \pi \frac{k}{N}} O_k $$
        out[k + N / 2] = even_out[k] - twiddle_factor_x_odd;
    }

    return out;
}

// decimation-in-time (DIT)
template <direction direction_, typename FloatT, std::size_t extent_>
std::vector<std::complex<FloatT>> fft_impl(std::span<const std::complex<FloatT>, extent_> in) {
    const std::size_t N = in.size();
    const auto half_N_bit_width = static_cast<std::size_t>(std::bit_width(N >> 1));

    std::vector<std::complex<FloatT>> out(N);

    // step 1: bit-reversal permutation
    for (std::size_t k = 0; k < N; ++k) {
        const std::size_t k_bitswapped = bitswap(k, half_N_bit_width);
        out[k] = in[k_bitswapped];
    }

    // step 2: iterative computation
    for (std::size_t stride = 1; stride <= N; stride <<= 1) {
        // compute the fundamental frequency twiddle factor for this stride
        // $$ \omega = e^{-i 2 \pi \frac{1}{stride}} $$
        const auto fundamental_freq = detail::polar(detail::theta<direction_, FloatT>(1, stride));

        // perform FFT for each segment of the current stride
        for (size_t offset = 0; offset < N; offset += stride) {
            // initialize twiddle factor
            // $$ \omega^0 = 1 $$
            std::complex<FloatT> twiddle_factor{ 1 };

            for (size_t k = 0; k != stride / 2; ++k) {
                const auto even /*           */ = /*            */ out[offset + k];
                const auto twiddle_factor_x_odd = twiddle_factor * out[offset + k + stride / 2];

                // apply the butterfly operation
                out[offset + k] /*        */ = even + twiddle_factor_x_odd;
                out[offset + k + stride / 2] = even - twiddle_factor_x_odd;

                // update the twiddle factor for the next butterfly operation
                // $$ \omega^1 = e^{-i 2 \pi \frac{1}{stride}} $$
                // $$ \omega^2 = e^{-i 2 \pi \frac{2}{stride}} $$
                // ...
                twiddle_factor *= fundamental_freq;
            }
        }
    }

    return out;
}

} // namespace detail

template <direction direction_, typename FloatT, std::size_t extent_>
    requires std::is_floating_point_v<FloatT> && detail::extent_is_power_of_2<extent_>
std::vector<std::complex<FloatT>> fft_recursive(std::span<const std::complex<FloatT>, extent_> in) {
    const std::size_t N = in.size();
    ASSERT(std::has_single_bit(N), "only accepting powers of 2");

    constexpr std::size_t starting_offset = 0;
    constexpr std::size_t starting_stride = 1;
    auto out = detail::fft_recursive_impl<direction_>(in, starting_offset, starting_stride);

    if constexpr (direction_ == direction::freq_to_time) {
        for (auto& out_elem : out) {
            out_elem /= static_cast<FloatT>(N);
        }
    }

    return out;
}

template <direction direction_, typename FloatT, std::size_t extent_>
    requires std::is_floating_point_v<FloatT> && detail::extent_is_power_of_2<extent_>
std::vector<std::complex<FloatT>> fft(std::span<const std::complex<FloatT>, extent_> in) {
    const std::size_t N = in.size();
    ASSERT(std::has_single_bit(N), "only accepting powers of 2");

    auto out = detail::fft_impl<direction_>(in);

    if constexpr (direction_ == direction::freq_to_time) {
        for (auto& out_elem : out) {
            out_elem /= static_cast<FloatT>(N);
        }
    }

    return out;
}

} // namespace sl::calc::fourier
