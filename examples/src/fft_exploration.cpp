//
// Created by usatiynyan on 12/24/23.
//

#include "sl/calc/fourier.hpp"

#include <sl/meta/assert.hpp>
#include <range/v3/view.hpp>

namespace sl::calc::fourier::detail {

template <direction direction_, typename FloatT, std::size_t extent_>
    requires std::is_floating_point_v<FloatT>
std::vector<std::complex<FloatT>> dft(std::span<const std::complex<FloatT>, extent_> in) {
    const std::size_t N = in.size();
    std::vector<std::complex<FloatT>> out(N);

    for (std::size_t k = 0; k != N; ++k) {
        for (std::size_t n = 0; n != N; ++n) {
            out[k] += in[n] * detail::polar(detail::theta<direction_, FloatT>(k * n, N));
        }
        if constexpr (direction_ == direction::freq_to_time) {
            out[k] /= static_cast<FloatT>(N);
        }
    }

    return out;
}

template <direction direction_, typename FloatT, std::size_t extent_>
    requires std::is_floating_point_v<FloatT> && detail::extent_is_power_of_2<extent_>
std::vector<std::complex<FloatT>> fft_first_step(std::span<const std::complex<FloatT>, extent_> in) {
    const std::size_t N = in.size();
    ASSERT(std::has_single_bit(N), "only accepting powers of 2");

    const auto even_in = in | ranges::views::stride(2);
    const auto odd_in = in | ranges::views::drop_exactly(1) | ranges::views::stride(2);

    std::vector<std::complex<FloatT>> even_out = dft<direction_>(std::span{ even_in });
    std::vector<std::complex<FloatT>> odd_out = dft<direction_>(std::span{ odd_in });

    std::vector<std::complex<FloatT>> out(N);

    for (std::size_t k = 0; k < N / 2; ++k) {
        // $$ e^{-i 2 \pi \frac{k}{N}} $$
        const auto twiddle_factor = detail::polar(detail::theta<direction_, FloatT>(k, N));
        // $$ X_k         = E_k + e^{-i 2 \pi \frac{k}{N}} O_k $$
        out[k /*   */] = even_out[k] + twiddle_factor * odd_out[k];
        // $$ X_{k + N/2} = E_k - e^{-i 2 \pi \frac{k}{N}} O_k $$
        out[k + N / 2] = even_out[k] - twiddle_factor * odd_out[k];
    }

    return out;
}

template <direction direction_, typename FloatT, std::size_t extent_>
    requires std::is_floating_point_v<FloatT> && detail::extent_is_power_of_2<extent_>
std::vector<std::complex<FloatT>> fft_recursive(std::span<const std::complex<FloatT>, extent_> in) {
    const std::size_t N = in.size();
    ASSERT(std::has_single_bit(N), "only accepting powers of 2");

    if (N == 1) {
        return { in[0] };
    }

    const auto even_in = in | ranges::views::stride(2);
    const auto odd_in = in | ranges::views::drop_exactly(1) | ranges::views::stride(2);

    std::vector<std::complex<FloatT>> even_out = fft_recursive<direction_>(std::span{ even_in });
    std::vector<std::complex<FloatT>> odd_out = fft_recursive<direction_>(std::span{ odd_in });

    std::vector<std::complex<FloatT>> out(N);

    for (std::size_t k = 0; k < N / 2; ++k) {
        // $$ e^{-i 2 \pi \frac{k}{N}} $$
        const auto twiddle_factor = detail::polar(detail::theta<direction_, FloatT>(k, N));
        // $$ X_k         = E_k + e^{-i 2 \pi \frac{k}{N}} O_k $$
        out[k /*   */] = even_out[k] + twiddle_factor * odd_out[k];
        // $$ X_{k + N/2} = E_k - e^{-i 2 \pi \frac{k}{N}} O_k $$
        out[k + N / 2] = even_out[k] - twiddle_factor * odd_out[k];
    }

    return out;
}

template <direction direction_, typename FloatT, std::size_t extent_>
    requires std::is_floating_point_v<FloatT> && detail::extent_is_power_of_2<extent_>
std::vector<std::complex<FloatT>> fft_recursive_inplace(
    std::span<const std::complex<FloatT>, extent_> in,
    std::size_t offset = 0,
    std::size_t stride = 1
) {
    const std::size_t N = in.size() / stride;
    ASSERT(std::has_single_bit(N), "only accepting powers of 2");

    if (N == 1) {
        return { in[offset] };
    }

    std::vector<std::complex<FloatT>> even_out = fft_recursive_inplace<direction_>(in, offset, stride * 2);
    std::vector<std::complex<FloatT>> odd_out = fft_recursive_inplace<direction_>(in, offset + stride, stride * 2);

    std::vector<std::complex<FloatT>> out(N);

    for (std::size_t k = 0; k < N / 2; ++k) {
        // $$ e^{-i 2 \pi \frac{k}{N}} $$
        const auto twiddle_factor = detail::polar(detail::theta<direction_, FloatT>(k, N));
        // $$ X_k         = E_k + e^{-i 2 \pi \frac{k}{N}} O_k $$
        out[k /*   */] = even_out[k] + twiddle_factor * odd_out[k];
        // $$ X_{k + N/2} = E_k - e^{-i 2 \pi \frac{k}{N}} O_k $$
        out[k + N / 2] = even_out[k] - twiddle_factor * odd_out[k];
    }

    return out;
}

} // namespace sl::calc::fourier::detail


int main() {}
