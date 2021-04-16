#ifndef INTEGER_H__
#define INTEGER_H__


#include "share_matrix.h"

inline Matrix from_uint32(std::uint32_t x) {
  auto out = Matrix::vector(32);
  for (std::size_t i = 0; i < 32; ++i) {
    out[i] = x & 1;
    x >>= 1;
  }
  return out;
}


inline std::uint32_t to_uint32(const Matrix& m) {
  std::uint32_t out = 0;
  for (std::size_t i = 0; i < 32; ++i) {
    out |= (m[i] << i);
  }
  return out;
}


template <Mode mode>
ShareMatrix<mode> integer_add(
    std::size_t bits_to_add, const ShareMatrix<mode>& x, const ShareMatrix<mode>& y) {
  const auto n = bits_to_add;

  assert (x.cols() == 1);
  assert (y.cols() == 1);
  assert (x.rows() > n);
  assert (y.rows() > n);

  auto out = ShareMatrix<mode>::vector(n);

  auto carry = Share<mode>::bit(false);
  for (std::size_t i = 0; i < n-1; ++i) {
    const auto xc = x[i] ^ carry;
    out[i] = xc ^ y[i];
    carry ^= xc & (y[i] ^ carry);
  }
  out[n-1] = x[n-1] ^ y[n-1] ^ carry;
  return out;
}

template <Mode mode>
ShareMatrix<mode> integer_add(
    const ShareMatrix<mode>& x, const ShareMatrix<mode>& y) {
  return integer_add(x.rows(), x, y);
}


template <Mode mode>
ShareMatrix<mode> integer_multiply(const ShareMatrix<mode>& x, const ShareMatrix<mode>& y) {
  const auto n = x.rows();

  assert(x.cols() == 1);
  assert(y.cols() == 1);
  assert(y.rows() == n);

  const auto xy = outer_product<mode>(x, y);

  auto sum = ShareMatrix<mode>::vector(1);
  sum[0] = xy(n-1, n-1);

  for (std::size_t i = 1; i < n; ++i) {
    auto sum_swap = ShareMatrix<mode>::vector(sum.rows() + 1);
    for (std::size_t j = 0; j < sum.rows(); ++j) {
      sum_swap[j+1] = sum[j];
    }
    sum = sum_swap;
    sum = integer_add(i, sum, xy.row(n-i-1));
  }
  return sum;
}



#endif
