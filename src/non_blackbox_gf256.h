#ifndef NON_BLACKBOX_GF256_H__
#define NON_BLACKBOX_GF256_H__

#include "share_matrix.h"
#include "gf256.h"


inline Matrix byte_to_vector(std::uint8_t x) {
  auto v = Matrix::vector(8);
  for (std::size_t i = 0; i < 8; ++i) {
    v[i] = x & 1;
    x >>= 1;
  }
  return v;
}


inline std::uint8_t vector_to_byte(const Matrix& v) {
  std::uint8_t x = 0;
  for (std::size_t i = 0; i < 8; ++i) {
    x |= (v[i] << i);
  }
  return x;
}


// The truth table that multiplies its argument by c in the field GF(256).
inline Matrix mul_by_constant_table(std::uint8_t c) {
  return truth_table([c](const auto& x) { 
    return byte_to_vector(mul_gf256(vector_to_byte(x), c));
  }, 8, 8);
}


// Multiply (x + color(x)) by y in GF(256)
template <Mode mode>
ShareMatrix<mode> half_mul_gf256(const ShareMatrix<mode>& x, const ShareMatrix<mode>& y);


// Multiply x by y in GF(256)
template <Mode mode>
ShareMatrix<mode> mul_gf256(const ShareMatrix<mode>& x, const ShareMatrix<mode>& y) {
  return
    half_mul_gf256(x, y) ^
    half_mul_gf256(y, ShareMatrix<mode>::constant(color<mode>(x))) ^
    ShareMatrix<mode>::constant(
        byte_to_vector(
          mul_gf256(
            vector_to_byte(color<mode>(x)),
            vector_to_byte(color<mode>(y)))));
}


// input x must be known to be non-zero
template <Mode mode>
ShareMatrix<mode> gf256_invert(const ShareMatrix<mode>& x);


template <Mode mode>
ShareMatrix<mode> aes_sbox(const ShareMatrix<mode>& x);


#endif
