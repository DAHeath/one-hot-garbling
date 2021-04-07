#ifndef GF256_H__
#define GF256_H__

#include <cstdint>


std::uint8_t mul_gf256(std::uint8_t, std::uint8_t);
std::uint8_t invert_gf256(std::uint8_t);


#endif
