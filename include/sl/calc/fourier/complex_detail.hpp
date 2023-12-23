//
// Created by usatiynyan on 12/22/23.
//

#pragma once

#include <complex>
#include <numbers>
#include <type_traits>

namespace sl::calc::fourier {
enum class direction {
    forward,
    backward,
};

template <direction direction_, typename FloatT>
    requires std::is_floating_point_v<FloatT>
constexpr FloatT direction_sign_v = static_cast<FloatT>(direction_ == direction::forward ? -1.0 : +1.0);

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
} // namespace sl::calc::fourier
