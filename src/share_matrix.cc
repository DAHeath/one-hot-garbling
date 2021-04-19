#include "share_matrix.h"
#include "unary_outer_product.h"
#include "table.h"

#include <iostream>

constexpr std::size_t default_outer_product_slice_size = 5;


struct IdentityTable : public Table {
  std::size_t operator()(std::size_t i) const {
    return i;
  }
};

IdentityTable the_identity_table { };


template <Mode mode>
void partial_half_outer_product(
    const ShareCSpan<mode>& x,
    const ShareCSpan<mode>& y,
    const MatrixView<Share<mode>&>& out) {
  assert(x.cols() == 1);
  assert(y.cols() == 1);
  const auto n = x.rows();
  const auto m = y.rows();
  assert(out.rows() == n);
  assert(out.cols() == m);

  unary_outer_product<mode>(the_identity_table, x, y, out);
}


template <Mode mode>
void half_outer_product(
    const ShareCSpan<mode>& X,
    const ShareCSpan<mode>& Y,
    const MatrixView<Share<mode>&>& out) {
  assert(X.cols() == 1);
  assert(Y.cols() == 1);

  const auto n = X.rows();
  const auto m = Y.rows();

  assert(out.rows() == n);
  assert(out.cols() == m);

  const auto def = default_outer_product_slice_size;

  for (std::size_t s = 0; s < (n + def-1)/def; ++s) {
    std::size_t slice_size = def;
    if (slice_size*(s + 1) > n) { slice_size = n % slice_size; }

    ShareMatrix<mode> slice(slice_size, 1);
    for (std::size_t i = 0; i < slice_size; ++i) {
      slice[i] = X[i + s*def];
    }
    auto part = subrows(def*s, def*s + slice_size, out);
    partial_half_outer_product<mode>(slice, Y, part);
  }
}

template void half_outer_product(
    const ShareCSpan<Mode::G>&,
    const ShareCSpan<Mode::G>&,
    const MatrixView<Share<Mode::G>&>&);
template void half_outer_product(
    const ShareCSpan<Mode::E>&,
    const ShareCSpan<Mode::E>&,
    const MatrixView<Share<Mode::E>&>&);




template <Mode mode>
void outer_product(ShareCSpan<mode> X, ShareCSpan<mode> Y, ShareMatrix<mode>& out) {
  // An outer product can be computed by half outer products.
  // Unfortunately, due to exponential computation scaling, we need to "chunk" the outer product into slices.
  assert(X.cols() == 1);
  assert(Y.cols() == 1);

  const auto x = ShareMatrix<mode>::constant(X.color());
  const auto xy = ShareMatrix<mode>::constant(X.color().outer_product(Y.color()));

  half_outer_product<mode>(X, Y, out);
  out.transpose();
  half_outer_product<mode>(Y, x, out);
  out.transpose();
  out ^= xy;
}


template ShareMatrix<Mode::G> outer_product(const ShareCSpan<Mode::G>&, const ShareCSpan<Mode::G>&);
template ShareMatrix<Mode::E> outer_product(const ShareCSpan<Mode::E>&, const ShareCSpan<Mode::E>&);


template <Mode mode>
ShareMatrix<mode> ShareMatrix<mode>::operator*(const ShareMatrix<mode>& o) const {
  const auto l = rows();
  const auto n = o.rows();
  const auto m = o.cols();
  assert (cols() == n);

  ShareMatrix<mode> out(l, m);
  for (std::size_t i = 0; i < n; ++i) {
    outer_product(column(i), o.row(i), out);
  }
  return out;
}


template ShareMatrix<Mode::G> ShareMatrix<Mode::G>::operator*(const ShareMatrix<Mode::G>&) const;
template ShareMatrix<Mode::E> ShareMatrix<Mode::E>::operator*(const ShareMatrix<Mode::E>&) const;


template <Mode mode>
ShareMatrix<mode> naive_matrix_multiplication(const ShareMatrix<mode>& x, const ShareMatrix<mode>& y) {
  const auto l = x.rows();
  const auto n = y.rows();
  const auto m = y.cols();
  assert (x.cols() == n);

  ShareMatrix<mode> out(l, m);
  for (std::size_t i = 0; i < n; ++i) {
    naive_outer_product<mode>(x.column(i), y.row(i), out);
  }
  return out;
}


template ShareMatrix<Mode::G> naive_matrix_multiplication<Mode::G>(
    const ShareMatrix<Mode::G>&, const ShareMatrix<Mode::G>&);
template ShareMatrix<Mode::E> naive_matrix_multiplication<Mode::E>(
    const ShareMatrix<Mode::E>&, const ShareMatrix<Mode::E>&);
