//
// Created by usatiynyan on 12/22/23.
//

#pragma once

#include <bit>
#include <complex>
#include <numbers>
#include <type_traits>

namespace sl::calc::fourier {
enum class direction {
    time_to_freq,
    freq_to_time,
};

namespace detail {

template <std::size_t extent_>
constexpr bool extent_is_power_of_2 = extent_ == std::dynamic_extent || std::has_single_bit(extent_);

template <direction direction_, typename FloatT>
    requires std::is_floating_point_v<FloatT>
constexpr FloatT direction_sign_v = static_cast<FloatT>(direction_ == direction::time_to_freq ? -1.0 : +1.0);

//      2pi * numerator
// sign ---------------
//        denominator
template <direction direction_, typename FloatT>
    requires std::is_floating_point_v<FloatT>
FloatT theta(size_t numerator, size_t denominator) {
    constexpr FloatT two_pi = 2 * std::numbers::pi_v<FloatT>;
    constexpr FloatT sign = direction_sign_v<direction_, FloatT>;
    return sign * two_pi * static_cast<FloatT>(numerator) / static_cast<FloatT>(denominator);
}

// slightly faster std::polar
template <typename FloatT>
    requires std::is_floating_point_v<FloatT>
constexpr auto polar(FloatT theta) {
    return std::complex<FloatT>{ std::cos(theta), std::sin(theta) };
}

} // namespace detail
} // namespace sl::calc::fourier
