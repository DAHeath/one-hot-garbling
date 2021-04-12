#ifndef SHARE_MATRIX_H__
#define SHARE_MATRIX_H__


#include "share.h"
#include "matrix.h"
#include <vector>


template <Mode mode>
struct ShareMatrix {
public:
  ShareMatrix() { }
  ShareMatrix(std::size_t n, std::size_t m) :
    transposed(false), n(n), m(m), vals(n*m) { }

  static ShareMatrix uniform(std::size_t n, std::size_t m) {
    ShareMatrix out(n, m);
    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t j = 0; j < m; ++j) {
        out(i, j) = Share<mode>::uniform();
      }
    }
    return out;
  }

  static ShareMatrix constant(const Matrix& c) {
    const auto n = c.rows();
    const auto m = c.cols();
    ShareMatrix out(n, m);
    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t j = 0; j < m; ++j) {
        out(i, j) = Share<mode>::bit(c(i, j));
      }
    }
    return out;
  }

  static ShareMatrix vector(std::size_t n) {
    return { n, 1 };
  }

  // Column vector containing the `i`th row.
  ShareMatrix<mode> row(std::size_t i) const {
    auto out = ShareMatrix<mode>::vector(cols());
    for (std::size_t j = 0; j < cols(); ++j) {
      out[j] = (*this)(i, j);
    }
    return out;
  }

  // Column vector containing the `j`th column.
  ShareMatrix<mode> column(std::size_t j) const {
    auto out = ShareMatrix<mode>::vector(rows());
    for (std::size_t i = 0; i < rows(); ++i) {
      out[i] = (*this)(i, j);
    }
    return out;
  }

  Share<mode>& operator()(std::size_t i, std::size_t j) {
    if (transposed) { return get(j, i); } else { return get(i, j); }
  }

  const Share<mode>& operator()(std::size_t i, std::size_t j) const {
    if (transposed) { return get(j, i); } else { return get(i, j); }
  }

  // vector access
  Share<mode>& operator[](std::size_t i) {
    assert(cols() == 1);
    return (*this)(i, 0);
  }

  // vector access
  const Share<mode>& operator[](std::size_t i) const {
    assert(cols() == 1);
    return (*this)(i, 0);
  }

  ShareMatrix& operator^=(const ShareMatrix& o) {
    assert (rows() == o.rows());
    assert (cols() == o.cols());

    for (std::size_t i = 0; i < rows(); ++i) {
      for (std::size_t j = 0; j < cols(); ++j) {
        (*this)(i, j) ^= o(i, j);
      }
    }
    return *this;
  }

  ShareMatrix operator^(const ShareMatrix& o) const {
    ShareMatrix out = *this;
    out ^= o;
    return out;
  }

  void transpose() { transposed = !transposed; }

  const std::size_t rows() const { return transposed ? m : n; }
  const std::size_t cols() const { return transposed ? n : m; }

  friend std::ostream& operator<<(std::ostream& os, const ShareMatrix& m) {
    for (std::size_t i = 0; i < m.rows(); ++i) {
      for (std::size_t j = 0; j < m.cols(); ++j) {
        os << m(i, j) << ' ';
      }
      os << '\n';
    }
    return os;
  }

  Matrix color() const {
    Matrix out(rows(), cols());
    for (std::size_t i = 0; i < rows(); ++i) {
      for (std::size_t j = 0; j < cols(); ++j) {
        out(i, j) = (*this)(i, j).color();
      }
    }
    return out;
  }

  // Let a be this vector and b be the other vector.
  // Let c be the color of this vector.
  //
  // This computes [|U(a + c) x b |] where x denotes the vector outer product.
  ShareMatrix<mode> unary_outer_product(const ShareMatrix<mode>&) const;

  ShareMatrix<mode> operator*(const ShareMatrix<mode>&) const;

private:
  Share<mode>& get(std::size_t i, std::size_t j) {
    return vals[j*n + i];
  }

  const Share<mode>& get(std::size_t i, std::size_t j) const {
    return vals[j*n + i];
  }

  bool transposed;
  std::size_t n;
  std::size_t m;
  std::vector<Share<mode>> vals;
};


template <Mode mode>
ShareMatrix<mode> operator*(const Matrix& x, const ShareMatrix<mode>& y) {
  const auto l = x.rows();
  const auto n = y.rows();
  const auto m = y.cols();
  assert (x.cols() == n);

  ShareMatrix<mode> out(l, m);
  for (std::size_t i = 0; i < l; ++i) {
    for (std::size_t j = 0; j < m; ++j) {
      for (std::size_t k = 0; k < n; ++k) {
        if (x(i, k)) { out(i, j) ^= y(k, j); }
      }
    }
  }
  return out;
}

template <Mode mode>
ShareMatrix<mode> outer_product(const ShareMatrix<mode>&, const ShareMatrix<mode>&);


#endif
