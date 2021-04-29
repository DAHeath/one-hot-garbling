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


inline Matrix from_uint64(std::uint64_t x) {
  auto out = Matrix::vector(64);
  for (std::size_t i = 0; i < 64; ++i) {
    out[i] = x & 1;
    x >>= 1;
  }
  return out;
}


inline std::uint64_t to_uint64(const Matrix& m) {
  std::uint64_t out = 0;
  for (std::size_t i = 0; i < 64; ++i) {
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
  MatrixView<const Share<mode>> xx = x;
  MatrixView<const Share<mode>> yy = y;
  return integer_add(x.rows(), xx, yy);
}


template <Mode mode>
ShareMatrix<mode> integer_sub(
    const MatrixView<const Share<mode>>& x,
    const MatrixView<const Share<mode>>& y) {
  const auto n = x.rows();
  assert (x.cols() == 1);
  assert (y.cols() == 1);
  assert (y.rows() == n);

  auto out = ShareMatrix<mode>::vector(n);

  auto borrow = Share<mode>::bit(false);
  for (std::size_t i = 0; i < n-1; ++i) {
    const auto xy = x[i] ^ y[i];
    const auto yc = borrow ^ y[i];
    out[i] = xy ^ borrow;
    borrow ^= xy & yc;
  }
  out[n-1] = x[n-1] ^ y[n-1] ^ borrow;
  return out;
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

  /* for (std::size_t i = 0; i < n; ++i) { */
  /* } */


  /* const auto xy = outer_product<mode>(x, y); */

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


template <Mode mode>
ShareMatrix<mode> naive_exponent(std::uint32_t x, const ShareMatrix<mode>& y) {
  auto out = ShareMatrix<mode>::constant(from_uint32(x));
  for (std::size_t i = 0; i < 32; ++i) {
    out[i] &= y[0];
  }

  for (std::size_t i = 1; i < 32; ++i) {
    auto mul = ShareMatrix<mode>::constant(from_uint32(x * (1 << i)));
    for (std::size_t i = 0; i < 32; ++i) {
      mul[i] &= y[i];
    }
    MatrixView<const Share<mode>> xx = out;
    MatrixView<const Share<mode>> yy = mul;
    out = naive_integer_multiply(xx, yy);
  }
  return out;
}


constexpr std::uint32_t pow32(std::uint32_t x, std::uint32_t p) {
  if (p == 0) return 1;
  if (p == 1) return x;
  
  std::uint32_t tmp = pow32(x, p/2);
  if (p % 2 == 0) return tmp * tmp;
  else return x * tmp * tmp;
}


struct ExpTable : public Table {
  ExpTable() { }
  ExpTable(std::uint32_t base, std::size_t shift) : base(base), shift(shift) { }

  std::size_t operator()(std::size_t i) const {
    return pow32(base, i << shift);
  }

  std::uint32_t base; 
  std::size_t shift;
};


template <Mode mode>
ShareMatrix<mode> exponent(std::uint32_t x, const ShareMatrix<mode>& y) {

  std::size_t n_chunks = (32 + chunking_factor() - 1) / chunking_factor();

  const auto mask = ShareMatrix<mode>::uniform(32, 1);
  auto masked = integer_sub<mode>(y, mask);
  masked.reveal();

  auto one = ShareMatrix<mode>(1, 1);
  one[0] = Share<mode>::bit(true);

  ShareMatrix<mode> result(32, 1);
  { // handle first chunk as a special case
    ShareMatrix<mode> chunk(chunking_factor(), 1);
    for (std::size_t i = 0; i < chunking_factor(); ++i) {
      chunk[i] = masked[i];
    }
    ExpTable etable { x, 0 };
    unary_outer_product<mode>(etable, chunk, one, result);
  }

  for (std::size_t i = 1; i < n_chunks; ++i) {
    // handle remaining chunks by multiplying them in
    std::size_t chunk_size = std::min(32 - i * chunking_factor(), chunking_factor());

    ShareMatrix<mode> chunk(chunk_size, 1);
    for (std::size_t j = 0; j < chunk_size; ++j) {
      chunk[j] = masked[j + i*chunking_factor()];
    }
    ExpTable etable { x, i*chunking_factor() };
    ShareMatrix<mode> mul(32, 1);
    unary_outer_product<mode>(etable, chunk, one, mul);
    result = integer_multiply<mode>(result, mul);
  }

  // strip off mask by multiplication
  const auto pow_mask = ShareMatrix<mode>::constant(from_uint32(pow32(x, to_uint32(color<mode>(mask)))));

  result = integer_multiply<mode>(result, pow_mask);

  return result;
}


// if s, x, else y
template <Mode mode>
ShareMatrix<mode> swap(
    const Share<mode>& s,
    const ShareMatrix<mode>& x,
    const ShareMatrix<mode>& y) {

  const auto n = x.rows();
  const auto m = x.cols();

  assert(y.rows() == n);
  assert(y.cols() == m);

  auto diff = x ^ y;

  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < m; ++j) {
      diff(i, j) &= s;
    }
  }

  return diff ^ y;
}



template <Mode mode>
ShareMatrix<mode> sub_if_greater(const ShareMatrix<mode>& x, const ShareMatrix<mode>& y) {
  assert(x.rows() == 32 && y.rows() == 32);
  const auto diff = integer_sub<mode>(x, y);

  const auto gt = ((x[31] == y[31]) & (x[31] != diff[31])) | (x[31] & ~y[31]);

  return swap<mode>(gt, x, diff);
}


constexpr std::uint32_t p = 65521;

template <Mode mode>
ShareMatrix<mode> naive_mod_p(ShareMatrix<mode> x) {
  for (std::uint32_t i = 0; i <= 16; ++i) {
    x = sub_if_greater(x, ShareMatrix<mode>::constant(from_uint32(p * (1 << (16-i)))));
  }
  return x;
}


struct ModpTable : public Table {
  ModpTable() { }
  ModpTable(std::size_t shift) : shift(shift) { }

  std::size_t operator()(std::size_t i) const {
    return (i * (1 << shift)) % p;
  }

  std::size_t shift;
};



template <Mode mode>
ShareMatrix<mode> mod_p(const ShareMatrix<mode>& x) {
  auto one = ShareMatrix<mode>(1, 1);
  one[0] = Share<mode>::bit(true);

  /* const auto mask = ShareMatrix<mode>::uniform(32, 1); */
  const auto mask = ShareMatrix<mode>(32, 1);
  auto masked = integer_sub<mode>(x, mask);
  masked.reveal();

  // split into three chunks
  ShareMatrix<mode> low(32, 1);
  ShareMatrix<mode> mid(32, 1);
  ShareMatrix<mode> high(32, 1);
  for (std::size_t i = 0; i < 16; ++i) { low[i] = masked[i]; }

  { // use hot to compute mod p for mid chunk
    ShareMatrix<mode> mid_chunk(8, 1);
    for (std::size_t i = 0; i < 8; ++i) { mid_chunk[i] = masked[i+16]; }

    ModpTable table(16);

    unary_outer_product<mode>(table, mid_chunk, one, mid);
  }

  { // use hot to compute mod p for high chunk
    ShareMatrix<mode> high_chunk(8, 1);
    for (std::size_t i = 0; i < 8; ++i) { high_chunk[i] = masked[i+24]; }

    ModpTable table(24);

    unary_outer_product<mode>(table, high_chunk, one, high);
  }

  auto out = integer_add<mode>(low, mid);
  out = integer_add<mode>(out, high);
  out = integer_add<mode>(out,
      ShareMatrix<mode>::constant(from_uint32(to_uint32(color<mode>(mask)) % p)));

  // sum cannot be more than 5p-1
  out = sub_if_greater(out, ShareMatrix<mode>::constant(from_uint32(p * 4)));
  out = sub_if_greater(out, ShareMatrix<mode>::constant(from_uint32(p * 2)));
  out = sub_if_greater(out, ShareMatrix<mode>::constant(from_uint32(p)));

  return out;
}


#endif
