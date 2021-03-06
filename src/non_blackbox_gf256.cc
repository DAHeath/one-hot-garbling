#include "non_blackbox_gf256.h"
#include "table.h"


Matrix make_reduction_table() {
  // Modular reduction in GF(256) can be computed by a matrix.
  // https://web.stanford.edu/class/ee387/handouts/notes18.pdf

  Matrix table(15, 8);
  // top of table is identity
  for (std::size_t i = 0; i < 8; ++i) {
    table(i, i) = 1;
  }
  table(8, 0) = 1; table(8, 1) = 1; table(8, 3) = 1; table(8, 4) = 1;
  table(9, 1) = 1; table(9, 2) = 1; table(9, 4) = 1; table(9, 5) = 1;
  table(10, 2) = 1; table(10, 3) = 1; table(10, 5) = 1; table(10, 6) = 1;
  table(11, 3) = 1; table(11, 4) = 1; table(11, 6) = 1; table(11, 7) = 1;
  table(12, 0) = 1; table(12, 1) = 1; table(12, 3) = 1; table(12, 5) = 1; table(12, 7) = 1;
  table(13, 0) = 1; table(13, 2) = 1; table(13, 3) = 1; table(13, 6) = 1;
  table(14, 1) = 1; table(14, 3) = 1; table(14, 4) = 1; table(14, 7) = 1;
  table.transpose();
  return table;
}


struct InverseTable : public Table {
  std::size_t operator()(std::size_t i) const {
    return invert_gf256(i);
  }
};


Matrix reduction_table = make_reduction_table();



// Multiply (x + color(x)) by y in GF(256)
template <Mode mode>
ShareMatrix<mode> half_mul_gf256(const ShareMatrix<mode>& x, const ShareMatrix<mode>& y) {
  const auto xy_outer = half_outer_product<mode>(x, y);

  // The vector of coefficients of x*y without modular reduction.
  auto unreduced = ShareMatrix<mode>::vector(15);
  for (std::size_t i = 0; i < 8; ++i) {
    for (std::size_t j = 0; j < 8; ++j) {
      unreduced[i + j] ^= xy_outer(i, j);
    }
  }
  return reduction_table * unreduced;
}


template ShareMatrix<Mode::G> half_mul_gf256(const ShareMatrix<Mode::G>&, const ShareMatrix<Mode::G>&);
template ShareMatrix<Mode::E> half_mul_gf256(const ShareMatrix<Mode::E>&, const ShareMatrix<Mode::E>&);


template <Mode mode>
ShareMatrix<mode> gf256_invert(const ShareMatrix<mode>& x) {
  auto y = ShareMatrix<mode>::vector(8);

  // G draws a uniform, non-zero value.
  if constexpr (mode == Mode::G) {
    while (vector_to_byte(color<mode>(y)) == 0) {
      y = ShareMatrix<mode>::uniform(8, 1);
    }
  }

  auto cols = ShareMatrix<mode>::constant(
      byte_to_vector(mul_gf256(vector_to_byte(color<mode>(x)), vector_to_byte(color<mode>(y)))));
  MatrixView<const Share<mode>> colsv = cols;
  auto xy = half_mul_gf256(x, y);
  MatrixView<Share<mode>> xyv = xy;


  xyv ^= colsv;

  // It is secure to show x*y to E
  xy.reveal();

  InverseTable inv;
  ShareMatrix<mode> xy_outer(8, 8);
  unary_outer_product<mode>(inv, xy, y, xy_outer);

  auto unreduced = ShareMatrix<mode>::vector(15);
  for (std::size_t i = 0; i < 8; ++i) {
    for (std::size_t j = 0; j < 8; ++j) {
      unreduced[i + j] ^= xy_outer(i, j);
    }
  }
  return reduction_table * unreduced;
}


template ShareMatrix<Mode::G> gf256_invert(const ShareMatrix<Mode::G>&);
template ShareMatrix<Mode::E> gf256_invert(const ShareMatrix<Mode::E>&);


Matrix make_aes_linear_matrix() {
  Matrix m(8, 8);
  for (std::size_t i = 0; i < 8; ++i) {
    m(i, (0 + i)%8 ) = 1;
    m(i, (4 + i)%8) = 1;
    m(i, (5 + i)%8) = 1;
    m(i, (6 + i)%8) = 1;
    m(i, (7 + i)%8) = 1;
  }
  return m;
}

Matrix make_aes_linear_shift() {
  Matrix m(8, 1);
  m[0] = 1;
  m[1] = 1;
  m[5] = 1;
  m[6] = 1;
  return m;
}

Matrix aes_linear_matrix = make_aes_linear_matrix();
Matrix aes_linear_shift = make_aes_linear_shift();


template <Mode mode>
ShareMatrix<mode> aes_sbox(const ShareMatrix<mode>& x) {

  auto zero = ~x[0];
  for (std::size_t i = 1; i < 8; ++i) {
    zero &= x[i];
  }
  auto z = ShareMatrix<mode>::vector(8);
  z[0] = zero;

  return
    (aes_linear_matrix *
      (gf256_invert<mode>(x ^ z) ^ z))
    ^ ShareMatrix<mode>::constant(aes_linear_shift);
}


template ShareMatrix<Mode::G> aes_sbox(const ShareMatrix<Mode::G>&);
template ShareMatrix<Mode::E> aes_sbox(const ShareMatrix<Mode::E>&);
