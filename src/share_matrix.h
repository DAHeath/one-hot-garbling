#ifndef SHARE_MATRIX_H__
#define SHARE_MATRIX_H__


#include "share.h"
#include "matrix.h"
#include <vector>
#include <span>


template <Mode mode>
struct ShareSpan {
public:
  ShareSpan() { }
  ShareSpan(
    bool t,
    std::size_t n,
    std::size_t m,
    std::size_t i_step,
    std::size_t j_step,
    std::span<const Share<mode>> vals)
    : t(t), n(n), m(m), i_step(i_step), j_step(j_step), vals(vals) { }

  const Share<mode>& operator()(std::size_t i, std::size_t j) const {
    if (t) { return get(j, i); } else { return get(i, j); }
  }

  // vector access
  const Share<mode>& operator[](std::size_t i) const {
    assert(cols() == 1);
    return (*this)(i, 0);
  }

  ShareSpan transposed() const {
    return { !t, n, m, i_step, j_step, vals };
  }

  const std::size_t rows() const { return t ? m : n; }
  const std::size_t cols() const { return t ? n : m; }

  Matrix color() const {
    Matrix out(rows(), cols());
    for (std::size_t i = 0; i < rows(); ++i) {
      for (std::size_t j = 0; j < cols(); ++j) {
        out(i, j) = (*this)(i, j).color();
      }
    }
    return out;
  }

private:
  const Share<mode>& get(std::size_t i, std::size_t j) const {
    return vals[j*j_step + i*i_step];
  }

  bool t;
  std::size_t n;
  std::size_t m;
  std::size_t i_step;
  std::size_t j_step;
  std::span<const Share<mode>> vals;
};


template struct ShareSpan<Mode::G>;

template <Mode mode>
struct ShareMatrix {
public:
  ShareMatrix() { }
  ShareMatrix(std::size_t n, std::size_t m) :
    t(false), n(n), m(m), vals(n*m) { }

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

  operator ShareSpan<mode>() const {
    return { t, n, m, 1, n, vals };
  }

  static ShareMatrix vector(std::size_t n) {
    return { n, 1 };
  }

  // Column vector containing the `i`th row.
  ShareSpan<mode> row(std::size_t i) const {
    // TODO ensure transposition works
    std::span v = vals;
    v = v.subspan(i);
    return { false, m, 1, n, 1, v };
  }

  // Column vector containing the `j`th column.
  ShareSpan<mode> column(std::size_t j) const {
    // TODO ensure transposition works
    std::span v = vals;
    v = v.subspan(j*n);
    return { false, n, 1, 1, 1, v };
  }

  Share<mode>& operator()(std::size_t i, std::size_t j) {
    if (t) { return get(j, i); } else { return get(i, j); }
  }

  const Share<mode>& operator()(std::size_t i, std::size_t j) const {
    if (t) { return get(j, i); } else { return get(i, j); }
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

  ShareMatrix& operator^=(ShareSpan<mode> o) {
    assert (rows() == o.rows());
    assert (cols() == o.cols());

    for (std::size_t i = 0; i < rows(); ++i) {
      for (std::size_t j = 0; j < cols(); ++j) {
        (*this)(i, j) ^= o(i, j);
      }
    }
    return *this;
  }

  ShareMatrix operator^(ShareSpan<mode> o) const {
    ShareMatrix out = *this;
    out ^= o;
    return out;
  }

  void transpose() { t = !t; }

  ShareMatrix transposed() const {
    ShareMatrix out = *this;
    out.transpose();
    return out;
  }

  const std::size_t rows() const { return t ? m : n; }
  const std::size_t cols() const { return t ? n : m; }

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

  ShareMatrix<mode> operator*(const ShareMatrix<mode>&) const;

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

  // is the matrix transposed?
  bool t;
  std::size_t n;
  std::size_t m;
  std::vector<Share<mode>> vals;
};


template <Mode mode>
ShareMatrix<mode> operator*(const Matrix& x, ShareSpan<mode> y) {
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
  const ShareSpan<mode> yy = y;
  return x * yy;
}


template <Mode mode>
ShareMatrix<mode> operator*(ShareSpan<mode> x, const Matrix& y) {
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
  const ShareSpan<mode> xx = x;
  return xx * y;
}

template <Mode mode>
void outer_product(ShareSpan<mode>, ShareSpan<mode>, ShareMatrix<mode>&);

template <Mode mode>
ShareMatrix<mode> outer_product(ShareSpan<mode> x, ShareSpan<mode> y) {
  ShareMatrix<mode> out(x.rows(), y.cols());
  outer_product(x, y, out);
  return out;
}


// Computes (a + color(a)) x b where x denotes the vector outer product.
template <Mode mode>
void half_outer_product(const ShareSpan<mode>&, const ShareSpan<mode>&, ShareMatrix<mode>&);

template <Mode mode>
ShareMatrix<mode> half_outer_product(ShareSpan<mode> x, ShareSpan<mode> y) {
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
void naive_outer_product(const ShareSpan<mode>& x, const ShareSpan<mode>& y, ShareMatrix<mode>& out) {
  for (std::size_t i = 0; i < x.rows(); ++i) {
    for (std::size_t j = 0; j < y.rows(); ++j) {
      out(i, j) = x[i] & y[j];
    }
  }
}


template <Mode mode>
ShareMatrix<mode> naive_matrix_multiplication(const ShareMatrix<mode>&, const ShareMatrix<mode>&);


#endif
