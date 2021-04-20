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
    std::size_t bits_to_add,
    const MatrixView<const Share<mode>>& x,
    const MatrixView<const Share<mode>>& y) {
  const auto n = bits_to_add;

  assert (x.cols() == 1);
  assert (y.cols() == 1);
  assert (x.rows() >= n);
  assert (y.rows() >= n);

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
ShareMatrix<mode> integer_multiply(
    const MatrixView<const Share<mode>>& x,
    const MatrixView<const Share<mode>>& y) {
  const auto n = x.rows();

  assert(x.cols() == 1);
  assert(y.cols() == 1);
  assert(y.rows() == n);

  const auto xy = outer_product<mode>(x, y);
  MatrixView<const Share<mode>> xy_view = xy;

  auto sum = ShareMatrix<mode>::vector(1);
  sum[0] = xy(n-1, n-1);

  for (std::size_t i = 1; i < n; ++i) {
    auto sum_swap = ShareMatrix<mode>::vector(sum.rows() + 1);
    for (std::size_t j = 0; j < sum.rows(); ++j) {
      sum_swap[j+1] = sum[j];
    }
    sum = sum_swap;
    MatrixView<const Share<mode>> sum_view = sum;
    sum = integer_add(i+1, sum_view, row(n-i-1, xy_view));
  }
  return sum;
}


template <Mode mode>
ShareMatrix<mode> naive_integer_multiply(
    const MatrixView<const Share<mode>>& x,
    const MatrixView<const Share<mode>>& y) {
  const auto n = x.rows();

  assert(x.cols() == 1);
  assert(y.cols() == 1);
  assert(y.rows() == n);

  const auto xy = outer_product<mode>(x, y);

  auto sum = ShareMatrix<mode>::vector(1);
  sum[0] = x[n-1] & y[0];

  for (std::size_t i = 1; i < n; ++i) {
    auto sum_swap = ShareMatrix<mode>::vector(sum.rows() + 1);
    for (std::size_t j = 0; j < sum.rows(); ++j) {
      sum_swap[j+1] = sum[j];
    }
    sum = sum_swap;

    auto prod = ShareMatrix<mode>::vector(i + 1);
    for (std::size_t j = 0; j < i+1; ++j) {
      prod[j] = x[n-i-1] & y[j];
    }

    MatrixView<const Share<mode>> sum_view = sum;
    MatrixView<const Share<mode>> prod_view = prod;
    sum = integer_add(i+1, sum_view, prod_view);
  }
  return sum;
}



#endif
