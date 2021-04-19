#ifndef SHARE_MATRIX_H__
#define SHARE_MATRIX_H__


#include "share.h"
#include "matrix.h"
#include <vector>
#include <span>
#include <functional>


template <typename T>
struct MatrixView {
public:
  MatrixView(std::size_t n, std::size_t m, const std::function<T(std::size_t, std::size_t)>& access)
    : n(n), m(m), access(access) { }

  std::size_t rows() const { return n; }
  std::size_t cols() const { return m; }
  T operator()(std::size_t i, std::size_t j) const { return access(i, j); }
  T operator[](std::size_t i) const {
    assert(cols() == 1);
    return (*this)(i, 0);
  }

private:
  std::size_t n;
  std::size_t m;
  std::function<T(std::size_t, std::size_t)> access;
};


template <typename T>
MatrixView<T&> matrix_span(std::size_t n, std::size_t m, std::span<T> vals) {
  return {
    n, m,
    [=](std::size_t i, std::size_t j) -> T& { return vals[j*n + i]; }
  };
}


template <typename T>
MatrixView<const T&> matrix_const_span(std::size_t n, std::size_t m, std::span<const T> vals) {
  return {
    n, m,
    [=](std::size_t i, std::size_t j) -> const T& { return vals[j*n + i]; }
  };
}


template <typename T>
MatrixView<T> transpose(const MatrixView<T>& mat) {
  return {
    mat.cols(),
    mat.rows(),
    [=](std::size_t i, std::size_t j) -> T { return mat(j, i); }
  };
}


template <typename T>
MatrixView<T> subrows(std::size_t i0, std::size_t n, const MatrixView<T>& mat) {
  std::function<T(std::size_t, std::size_t)> f = [=](std::size_t i, std::size_t j) -> T { return mat(i + i0, j); };
  return {
    n,
    mat.cols(),
    f
  };
}


template <typename T>
MatrixView<T> column(std::size_t j0, const MatrixView<T>& mat) {
  return {
    mat.rows(),
    1,
    [=](std::size_t i, std::size_t j) -> T { return mat(i, j + j0); }
  };
}


template <typename T>
MatrixView<T> row(std::size_t i, const MatrixView<T>& mat) {
  return column(i, transpose(mat));
}


template <Mode mode>
struct ShareMatrix {
public:
  ShareMatrix() { }
  ShareMatrix(std::size_t n, std::size_t m) :
    n(n), m(m), vals(n*m) { }

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

  operator MatrixView<const Share<mode>&>() const {
    return matrix_const_span(n, m, std::span<const Share<mode>> { vals });
  }

  operator MatrixView<Share<mode>&>() {
    return matrix_span(n, m, std::span { vals });
  }

  static ShareMatrix vector(std::size_t n) { return { n, 1 }; }

/*   // Column vector containing the `i`th row. */
/*   MatrixView<const Share<mode>&> row(std::size_t i) const { */
/*     std::span<const Share<mode>> v = vals; */
/*     v = v.subspan(i); */
/*     return transpose(matrix_const_span(1, m, v)); */
/*   } */

/*   // Column vector containing the `j`th column. */
/*   MatrixView<const Share<mode>&> column(std::size_t j) const { */
/*     std::span v = vals; */
/*     v = v.subspan(j*n); */
/*     return matrix_const_span(n, 1, v); */
/*   } */

  Share<mode>& operator()(std::size_t i, std::size_t j) {
    return vals[j*n + i];
  }

  const Share<mode>& operator()(std::size_t i, std::size_t j) const {
    return vals[j*n + i];
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

  const std::size_t rows() const { return n; }
  const std::size_t cols() const { return m; }

  friend std::ostream& operator<<(std::ostream& os, const ShareMatrix& m) {
    for (std::size_t i = 0; i < m.rows(); ++i) {
      for (std::size_t j = 0; j < m.cols(); ++j) {
        os << m(i, j) << ' ';
      }
      os << '\n';
    }
    return os;
  }

  void reveal() {
    for (auto& s: vals) { s.reveal(); }
  }

private:
  Share<mode>& get(std::size_t i, std::size_t j) {
    return vals[j*n + i];
  }

  const Share<mode>& get(std::size_t i, std::size_t j) const {
    return vals[j*n + i];
  }

  std::size_t n;
  std::size_t m;
  std::vector<Share<mode>> vals;
};


template <Mode mode>
Matrix color(const MatrixView<const Share<mode>&>& x) {
  Matrix out(x.rows(), x.cols());
  for (std::size_t i = 0; i < x.rows(); ++i) {
    for (std::size_t j = 0; j < x.cols(); ++j) {
      out(i, j) = x(i, j).color();
    }
  }
  return out;
}


template <Mode mode>
const MatrixView<Share<mode>&>& operator^=(
    const MatrixView<Share<mode>&>& x,
    const MatrixView<const Share<mode>&>& y) {
  assert (x.rows() == y.rows());
  assert (x.cols() == y.cols());

  for (std::size_t i = 0; i < x.rows(); ++i) {
    for (std::size_t j = 0; j < x.cols(); ++j) {
      x(i, j) ^= y(i, j);
    }
  }
  return x;
}


template <Mode mode>
ShareMatrix<mode> operator*(
    const MatrixView<const Share<mode>&>&,
    const MatrixView<const Share<mode>&>&);


template <Mode mode>
const MatrixView<Share<mode>&>& operator^=(
    const MatrixView<Share<mode>&>& x,
    const ShareMatrix<mode>& y) {
  MatrixView<const Share<mode>&> yy = y;
  x ^= yy;
  return x;
}


template <Mode mode>
ShareMatrix<mode> operator*(const Matrix& x, const MatrixView<const Share<mode>&>& y) {
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
ShareMatrix<mode> operator*(const Matrix& x, const ShareMatrix<mode>& y) {
  MatrixView<const Share<mode>&> yy = y;
  return x * yy;
}


template <Mode mode>
ShareMatrix<mode> operator*(const MatrixView<Share<mode>&>& x, const Matrix& y) {
  const auto l = x.rows();
  const auto n = y.rows();
  const auto m = y.cols();
  assert (x.cols() == n);

  ShareMatrix<mode> out(l, m);
  for (std::size_t i = 0; i < l; ++i) {
    for (std::size_t j = 0; j < m; ++j) {
      for (std::size_t k = 0; k < n; ++k) {
        if (y(k, j)) { out(i, j) ^= x(i, k); }
      }
    }
  }
  return out;
}

template <Mode mode>
ShareMatrix<mode> operator*(const ShareMatrix<mode>& x, const Matrix& y) {
  MatrixView<const Share<mode>&> xx = x;
  return xx * y;
}

template <Mode mode>
void outer_product(
    const MatrixView<const Share<mode>&>&,
    const MatrixView<const Share<mode>&>&,
    const MatrixView<Share<mode>&>&);

template <Mode mode>
ShareMatrix<mode> outer_product(
    const MatrixView<const Share<mode>&>& x,
    const MatrixView<const Share<mode>&>& y) {
  ShareMatrix<mode> out(x.rows(), y.cols());
  outer_product(x, y, out);
  return out;
}


// Computes (a + color(a)) x b where x denotes the vector outer product.
template <Mode mode>
void half_outer_product(
    const MatrixView<const Share<mode>&>&,
    const MatrixView<const Share<mode>&>&,
    const MatrixView<Share<mode>&>&);

template <Mode mode>
ShareMatrix<mode> half_outer_product(
    const MatrixView<const Share<mode>&>& x,
    const MatrixView<const Share<mode>&>& y) {
  ShareMatrix<mode> out(x.rows(), y.cols());
  half_outer_product(x, y, out);
  return out;
}



inline Matrix decode(const ShareMatrix<Mode::G>& g, const ShareMatrix<Mode::E>& e) {
  const auto n = g.rows();
  const auto m = g.cols();

  Matrix out(n, m);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < m; ++j) {
      out(i, j) = decode(g(i, j), e(i, j));
    }
  }
  return out;
}



template <Mode mode>
void naive_outer_product(
    const MatrixView<const Share<mode>&>& x,
    const MatrixView<const Share<mode>&>& y,
    ShareMatrix<mode>& out) {
  for (std::size_t i = 0; i < x.rows(); ++i) {
    for (std::size_t j = 0; j < y.rows(); ++j) {
      out(i, j) = x[i] & y[j];
    }
  }
}


template <Mode mode>
ShareMatrix<mode> naive_matrix_multiplication(const ShareMatrix<mode>&, const ShareMatrix<mode>&);


#endif
