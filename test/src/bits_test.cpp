//
// Created by usatiynyan on 12/23/23.
// Some tests here written with the power of AI, can you spot them?
//

#include "sl/calc/bits.hpp"

#include <gtest/gtest.h>

namespace sl::calc {

TEST(Bits, fillOnes) {
    static_assert(fill_ones(0b0u) == 0b0u);
    static_assert(fill_ones(0b1u) == 0b1u);
    static_assert(fill_ones(0b10u) == 0b11u);
    static_assert(fill_ones(0b11u) == 0b11u);
    static_assert(fill_ones(0b100u) == 0b111u);
    static_assert(fill_ones(0b101u) == 0b111u);
    static_assert(fill_ones(0b110u) == 0b111u);
    static_assert(fill_ones(0b111u) == 0b111u);
    static_assert(fill_ones(0b1000u) == 0b1111u);
    static_assert(fill_ones(0b1001u) == 0b1111u);
    static_assert(fill_ones(0b1010u) == 0b1111u);
    static_assert(fill_ones(0b1011u) == 0b1111u);
    static_assert(fill_ones(0b1100u) == 0b1111u);
    static_assert(fill_ones(0b1101u) == 0b1111u);
    static_assert(fill_ones(0b1110u) == 0b1111u);
    static_assert(fill_ones(0b1111u) == 0b1111u);
    static_assert(fill_ones(0b10000u) == 0b11111u);
}

TEST(Bits, repeatByte) {
    static_assert(repeat_byte<uint8_t>(std::byte{ 0xAA }) == 0xAAu);
    static_assert(repeat_byte<uint16_t>(std::byte{ 0xAA }) == 0xAAAAu);
    static_assert(repeat_byte<uint32_t>(std::byte{ 0xAA }) == 0xAAAAAAAAu);
    static_assert(repeat_byte<uint64_t>(std::byte{ 0xAA }) == 0xAAAAAAAAAAAAAAAAu);

    static_assert(repeat_byte<uint8_t>(std::byte{ 0x01 }) == 0x01u);
    static_assert(repeat_byte<uint16_t>(std::byte{ 0x01 }) == 0x0101u);
    static_assert(repeat_byte<uint32_t>(std::byte{ 0x01 }) == 0x01010101u);
    static_assert(repeat_byte<uint64_t>(std::byte{ 0x01 }) == 0x0101010101010101u);
}

TEST(Bits, bitswap) {
    static_assert(bitswap<uint8_t>(0b00000000) == uint8_t{ 0b00000000 });
    static_assert(bitswap<uint8_t>(0b11111111) == uint8_t{ 0b11111111 });
    static_assert(bitswap<uint8_t>(0b00000001) == uint8_t{ 0b10000000 });
    static_assert(bitswap<uint8_t>(0b00000010) == uint8_t{ 0b01000000 });
    static_assert(bitswap<uint8_t>(0b00000100) == uint8_t{ 0b00100000 });
    static_assert(bitswap<uint8_t>(0b00001000) == uint8_t{ 0b00010000 });
    static_assert(bitswap<uint8_t>(0b00010000) == uint8_t{ 0b00001000 });
    static_assert(bitswap<uint8_t>(0b00100000) == uint8_t{ 0b00000100 });
    static_assert(bitswap<uint8_t>(0b01000000) == uint8_t{ 0b00000010 });
    static_assert(bitswap<uint8_t>(0b10000000) == uint8_t{ 0b00000001 });
    static_assert(bitswap<uint8_t>(0b10100000) == uint8_t{ 0b00000101 });
}

TEST(Bits, partialBitswap) {
    static_assert(bitswap<uint8_t>(0b10011101, std::bit_width(0b0u)) == uint8_t{ 0b10011101 });
    static_assert(bitswap<uint8_t>(0b10011101, std::bit_width(0b1u)) == uint8_t{ 0b10011101 });
    static_assert(bitswap<uint8_t>(0b10011101, std::bit_width(0b10000000u)) == bitswap<uint8_t>(0b10011101));
    
    static_assert(bitswap<uint8_t>(0b10000110, std::bit_width(0b10u)) == uint8_t{ 0b10000101 });
    static_assert(bitswap<uint8_t>(0b10001100, std::bit_width(0b100u)) == uint8_t{ 0b10001001 });
    static_assert(bitswap<uint8_t>(0b10011000, std::bit_width(0b1000u)) == uint8_t{ 0b10010001 });
    static_assert(bitswap<uint8_t>(0b10110000, std::bit_width(0b10000u)) == uint8_t{ 0b10100001 });
    static_assert(bitswap<uint8_t>(0b11100000, std::bit_width(0b100000u)) == uint8_t{ 0b11000001 });
    static_assert(bitswap<uint8_t>(0b11000000, std::bit_width(0b1000000u)) == uint8_t{ 0b10000001 });
    static_assert(bitswap<uint8_t>(0b10000000, std::bit_width(0b10000000u)) == uint8_t{ 0b00000001 });
}

// Tests covering stdlib <bit> functionality, purely for documenting purposes.
TEST(Bits, isPowerOf2) {
    static_assert(!std::has_single_bit(0u));
    static_assert(std::has_single_bit(1u));
    static_assert(std::has_single_bit(2u));
    static_assert(!std::has_single_bit(3u));
    static_assert(std::has_single_bit(4u));
    static_assert(!std::has_single_bit(5u));
    static_assert(!std::has_single_bit(6u));
    static_assert(!std::has_single_bit(7u));
    static_assert(std::has_single_bit(8u));
}

TEST(Bits, mostSignificantByte) {
    static_assert(std::bit_floor(0b0u) == 0u);
    static_assert(std::bit_floor(0b1u) == 0b1u);
    static_assert(std::bit_floor(0b10u) == 0b10u);
    static_assert(std::bit_floor(0b11u) == 0b10u);
    static_assert(std::bit_floor(0b100u) == 0b100u);
    static_assert(std::bit_floor(0b101u) == 0b100u);
    static_assert(std::bit_floor(0b110u) == 0b100u);
    static_assert(std::bit_floor(0b111u) == 0b100u);
    static_assert(std::bit_floor(0b1000u) == 0b1000u);
    static_assert(std::bit_floor(0b1001u) == 0b1000u);
    static_assert(std::bit_floor(0b1010u) == 0b1000u);
    static_assert(std::bit_floor(0b1011u) == 0b1000u);
    static_assert(std::bit_floor(0b1100u) == 0b1000u);
    static_assert(std::bit_floor(0b1101u) == 0b1000u);
    static_assert(std::bit_floor(0b1110u) == 0b1000u);
    static_assert(std::bit_floor(0b1111u) == 0b1000u);
    static_assert(std::bit_floor(0b10000u) == 0b10000u);
}

TEST(Bits, mostSignificantByteShift) {
    static_assert(std::bit_width(0b0u) == 0u);
    static_assert(std::bit_width(0b1u) == 1u);
    static_assert(std::bit_width(0b10u) == 2u);
    static_assert(std::bit_width(0b11u) == 2u);
    static_assert(std::bit_width(0b100u) == 3u);
    static_assert(std::bit_width(0b101u) == 3u);
    static_assert(std::bit_width(0b110u) == 3u);
    static_assert(std::bit_width(0b111u) == 3u);
    static_assert(std::bit_width(0b1000u) == 4u);
    static_assert(std::bit_width(0b1001u) == 4u);
    static_assert(std::bit_width(0b1010u) == 4u);
    static_assert(std::bit_width(0b1011u) == 4u);
    static_assert(std::bit_width(0b1100u) == 4u);
    static_assert(std::bit_width(0b1101u) == 4u);
    static_assert(std::bit_width(0b1110u) == 4u);
    static_assert(std::bit_width(0b1111u) == 4u);
    static_assert(std::bit_width(0b10000u) == 5u);
}

} // namespace sl::calc
