//
// Created by usatiynyan on 12/22/23.
//

#pragma once

#include <bit>
#include <cstdint>
#include <numeric>
#include <type_traits>

namespace sl::calc {
namespace detail {

template <typename T>
constexpr std::size_t sizeof_bits() {
    constexpr std::size_t bits_in_byte = 8;
    return bits_in_byte * sizeof(T);
}

template <typename IntT>
concept bit_ops_supported = std::is_integral_v<IntT> && sizeof(IntT) <= sizeof(std::int64_t);

} // namespace detail

template <typename IntT>
    requires detail::bit_ops_supported<IntT>
constexpr IntT fill_ones(IntT n) {
    constexpr std::size_t byte_size = sizeof(IntT);
    n |= (n >> 1);
    n |= (n >> 2);
    n |= (n >> 4);
    if constexpr (byte_size >= 2) {
        n |= (n >> 8);
    }
    if constexpr (byte_size >= 4) {
        n |= (n >> 16);
    }
    if constexpr (byte_size >= 8) {
        n |= (n >> 32);
    }
    return n;
}

template <typename IntT>
    requires detail::bit_ops_supported<IntT>
constexpr IntT repeat_byte(std::byte byte) {
    constexpr std::size_t byte_size = sizeof(IntT);
    IntT result = static_cast<IntT>(byte);
    if constexpr (byte_size >= 2) {
        result |= result << 8;
    }
    if constexpr (byte_size >= 4) {
        result |= result << 16;
    }
    if constexpr (byte_size >= 8) {
        result |= result << 32;
    }
    return result;
}

template <typename IntT>
    requires detail::bit_ops_supported<IntT>
constexpr IntT bitswap(IntT n) {
    constexpr std::size_t byte_size = sizeof(IntT);
    constexpr auto cast = [](auto a_n) { return static_cast<IntT>(a_n); };
    n = cast((n & 0xAAAAAAAAAAAAAAAAu) >> 1 | (n & 0x5555555555555555u) << 1);
    n = cast((n & 0xCCCCCCCCCCCCCCCCu) >> 2 | (n & 0x3333333333333333u) << 2);
    n = cast((n & 0xF0F0F0F0F0F0F0F0u) >> 4 | (n & 0x0F0F0F0F0F0F0F0Fu) << 4);
    if constexpr (byte_size >= 2) {
        n = cast((n & 0xFF00FF00FF00FF00u) >> 8 | (n & 0x00FF00FF00FF00FFu) << 8);
    }
    if constexpr (byte_size >= 4) {
        n = cast((n & 0xFFFF0000FFFF0000u) >> 16 | (n & 0x0000FFFF0000FFFFu) << 16);
    }
    if constexpr (byte_size >= 8) {
        n = cast((n & 0xFFFFFFFF00000000u) >> 32 | (n & 0x00000000FFFFFFFFu) << 32);
    }
    return n;
}

template <typename IntT>
    requires detail::bit_ops_supported<IntT>
constexpr IntT bitswap(IntT n, std::size_t bit_width) {
    constexpr auto cast = [](auto a_n) { return static_cast<IntT>(a_n); };
    const IntT unchanged = cast(cast(n >> bit_width) << bit_width);
    const std::size_t shift = detail::sizeof_bits<IntT>() - bit_width;
    return unchanged | (bitswap(n) >> shift);
}
} // namespace sl::calc
