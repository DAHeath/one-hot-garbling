#include "share_matrix.h"
#include "gf256.h"


Matrix byte_to_vector(std::uint8_t x) {
  auto v = Matrix::vector(8);
  for (std::size_t i = 0; i < 8; ++i) {
    v[i] = x & 1;
    x >>= 1;
  }
  return v;
}


std::uint8_t vector_to_byte(const Matrix& v) {
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


// The truth table that multiplies its argument by c in the field GF(256).
inline Matrix inverse_table() {
  return truth_table([](const auto& x) { 
    return byte_to_vector(invert_gf256(vector_to_byte(x)));
  }, 8, 8);
}


// Multiply (x + color(x)) by y in GF(256)
template <Mode mode>
ShareMatrix<mode> half_mul_gf256(const ShareMatrix<mode>& x, const ShareMatrix<mode>& y) {
  const auto xy_outer = half_outer_product<mode>(x, y);

  // The vector of coefficients of x*y without modular reduction.
  auto xy_unreduced = ShareMatrix<mode>::vector(15);
  for (std::size_t i = 0; i < 8; ++i) {
    for (std::size_t j = 0; j < 8; ++j) {
      xy_unreduced[i + j] ^= xy_outer(i, j);
    }
  }


  // Modular reduction in GF(256) can be computed by a matrix.
  // https://web.stanford.edu/class/ee387/handouts/notes18.pdf
  Matrix reduction_table(15, 8);
  // top of table is identity
  for (std::size_t i = 0; i < 8; ++i) {
    reduction_table(i, i) = 1;
  }
  reduction_table(8, 0) = 1;
  reduction_table(8, 1) = 1;
  reduction_table(8, 3) = 1;
  reduction_table(8, 4) = 1;
  reduction_table(9, 1) = 1;
  reduction_table(9, 2) = 1;
  reduction_table(9, 4) = 1;
  reduction_table(9, 5) = 1;
  reduction_table(10, 2) = 1;
  reduction_table(10, 3) = 1;
  reduction_table(10, 5) = 1;
  reduction_table(10, 6) = 1;
  reduction_table(11, 3) = 1;
  reduction_table(11, 4) = 1;
  reduction_table(11, 6) = 1;
  reduction_table(11, 7) = 1;
  reduction_table(12, 0) = 1;
  reduction_table(12, 1) = 1;
  reduction_table(12, 3) = 1;
  reduction_table(12, 5) = 1;
  reduction_table(12, 7) = 1;
  reduction_table(13, 0) = 1;
  reduction_table(13, 2) = 1;
  reduction_table(13, 3) = 1;
  reduction_table(13, 6) = 1;
  reduction_table(14, 1) = 1;
  reduction_table(14, 3) = 1;
  reduction_table(14, 4) = 1;
  reduction_table(14, 7) = 1;

  return (xy_unreduced.transposed() * reduction_table).transposed();
}


// Multiply x by y in GF(256)
template <Mode mode>
ShareMatrix<mode> full_mul_gf256(const ShareMatrix<mode>& x, const ShareMatrix<mode>& y) {
  return
    half_mul_gf256(x, y) ^
    half_mul_gf256(y, ShareMatrix<mode>::constant(x.color())) ^
    ShareMatrix<mode>::constant(
        byte_to_vector(mul_gf256(vector_to_byte(x.color()), vector_to_byte(y.color()))));
}




// input x must be known to be non-zero
template <Mode mode>
ShareMatrix<mode> gf256_invert(const ShareMatrix<mode>& x) {
  auto y = ShareMatrix<mode>::vector(8);

  // G draws a uniform, non-zero value.
  if constexpr (mode == Mode::G) {
    while (vector_to_byte(y.color()) == 0) {
      y = ShareMatrix<mode>::uniform(8, 1);
    }
  }

  auto xy = half_mul_gf256(x, y) ^
    ShareMatrix<mode>::constant(
      byte_to_vector(mul_gf256(vector_to_byte(x.color()), vector_to_byte(y.color()))));

  // It is secure to show x*y to E
  xy.reveal();

  return inverse_table() * xy.unary_outer_product(y);
}


template <Mode mode>
ShareMatrix<mode> aes_sbox(const ShareMatrix<mode>& x) {
}
