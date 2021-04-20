#include "share_matrix.h"
#include "unary_outer_product.h"
#include "table.h"

#include <iostream>

constexpr std::size_t default_outer_product_slice_size = 8;


struct IdentityTable : public Table {
  std::size_t operator()(std::size_t i) const {
    return i;
  }
};

static IdentityTable the_identity_table { };


template <Mode mode>
void partial_half_outer_product(
    const MatrixView<const Share<mode>>& x,
    const MatrixView<const Share<mode>>& y,
    const MatrixView<Share<mode>>& out) {
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
    const MatrixView<const Share<mode>>& X,
    const MatrixView<const Share<mode>>& Y,
    const MatrixView<Share<mode>>& out) {
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
    auto part = subrows(def*s, slice_size, out);
    partial_half_outer_product<mode>(slice, Y, part);
  }
}

template void half_outer_product(
    const MatrixView<const Share<Mode::G>>&,
    const MatrixView<const Share<Mode::G>>&,
    const MatrixView<Share<Mode::G>>&);
template void half_outer_product(
    const MatrixView<const Share<Mode::E>>&,
    const MatrixView<const Share<Mode::E>>&,
    const MatrixView<Share<Mode::E>>&);




template <Mode mode>
void outer_product(
    const MatrixView<const Share<mode>>& X,
    const MatrixView<const Share<mode>>& Y,
    const MatrixView<Share<mode>>& out) {
  // An outer product can be computed by half outer products.
  // Unfortunately, due to exponential computation scaling, we need to "chunk" the outer product into slices.
  assert(X.cols() == 1);
  assert(Y.cols() == 1);

  const auto x = ShareMatrix<mode>::constant(color(X));
  const auto xy = ShareMatrix<mode>::constant(color(X).outer_product(color(Y)));


  half_outer_product<mode>(X, Y, out);
  half_outer_product<mode>(Y, x, transpose(out));
  out ^= xy;
}


/* template ShareMatrix<Mode::G> outer_product(const ShareCSpan<Mode::G>&, const ShareCSpan<Mode::G>&); */
/* template ShareMatrix<Mode::E> outer_product(const ShareCSpan<Mode::E>&, const ShareCSpan<Mode::E>&); */


template <Mode mode>
ShareMatrix<mode> operator*(
    const MatrixView<const Share<mode>>& x,
    const MatrixView<const Share<mode>>& y) {
  const auto l = x.rows();
  const auto n = y.rows();
  const auto m = y.cols();
  assert (x.cols() == n);

  ShareMatrix<mode> out(l, m);
  MatrixView<Share<mode>> out_view = out;
  for (std::size_t i = 0; i < n; ++i) {
    outer_product(column(i, x), row(i, y), out_view);
  }
  return out;
}


template ShareMatrix<Mode::G> operator*(
    const MatrixView<const Share<Mode::G>>&,
    const MatrixView<const Share<Mode::G>>&);
template ShareMatrix<Mode::E> operator*(
    const MatrixView<const Share<Mode::E>>&,
    const MatrixView<const Share<Mode::E>>&);


template <Mode mode>
ShareMatrix<mode> naive_matrix_multiplication(
    const MatrixView<const Share<mode>>& x,
    const MatrixView<const Share<mode>>& y) {
  const auto l = x.rows();
  const auto n = y.rows();
  const auto m = y.cols();
  assert (x.cols() == n);

  ShareMatrix<mode> out(l, m);
  for (std::size_t i = 0; i < n; ++i) {
    naive_outer_product<mode>(column(i, x), row(i, y), out);
  }
  return out;
}


template ShareMatrix<Mode::G> naive_matrix_multiplication<Mode::G>(
    const MatrixView<const Share<Mode::G>>&,
    const MatrixView<const Share<Mode::G>>&);
template ShareMatrix<Mode::E> naive_matrix_multiplication<Mode::E>(
    const MatrixView<const Share<Mode::E>>&,
    const MatrixView<const Share<Mode::E>>&);
