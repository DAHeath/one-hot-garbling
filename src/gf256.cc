#include "gf256.h"


#include <vector>


std::vector<std::uint8_t> init_mult_table() {
  std::vector<std::uint8_t> table(1 << 16);

  for (std::size_t x = 0; x < 256; ++x) {
    for (std::size_t y = 0; y < 256; ++y) {
      std::size_t a = x;
      std::size_t b = y;
      std::size_t c = 0;
      for (std::size_t i = 0; i < 8; ++i) {
        c ^= a & -(b & 1);
        b >>= 1;
        a <<= 1;
        a ^= (0x11B & -(a >> 8));
      }
      const std::size_t index = x << 8 | y;
      table[index] = c;
    }
  }
  return table;
}


std::vector<std::uint8_t> gf256_multiplication_table = init_mult_table();


std::uint8_t mul_gf256(std::uint8_t x, std::uint8_t y) {
  std::size_t index = x;
  index <<= 8;
  index |= y;
  return gf256_multiplication_table[index];
}


std::vector<std::uint8_t> init_invert_table() {
  std::vector<std::uint8_t> table(1 << 8);

  for (std::size_t x = 0; x < 256; ++x) {
    std::uint8_t z;

    z = x;
    for (std::size_t i = 0; i < 6; ++i) {
      z = mul_gf256(z, z);
      z = mul_gf256(z, x);
    }
    table[x] = mul_gf256(z, z);
  }
  return table;
}


std::vector<std::uint8_t> gf256_invert_table = init_invert_table();


std::uint8_t invert_gf256(std::uint8_t x) {
  return gf256_invert_table[x];
}
