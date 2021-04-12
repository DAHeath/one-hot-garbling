#ifndef MATRIX_H__
#define MATRIX_H__

#include <boost/dynamic_bitset.hpp>

struct Matrix {
private:
  auto get(std::size_t i, std::size_t j) {
    return vals[j*n + i];
  }

  auto get(std::size_t i, std::size_t j) const {
    return vals[j*n + i];
  }

public:
  Matrix() { }
  Matrix(std::size_t n, std::size_t m) : transposed(false), n(n), m(m), vals(n*m) { }
  static Matrix vector(std::size_t n) { return { n, 1 }; }


  auto operator()(std::size_t i, std::size_t j) {
    if (transposed) { return get(j, i); } else { return get(i, j); }
  }

  auto operator()(std::size_t i, std::size_t j) const {
    if (transposed) { return get(j, i); } else { return get(i, j); }
  }

  // vector access
  auto operator[](std::size_t i) {
    assert(cols() == 1);
    return (*this)(i, 0);
  }

  // vector access
  auto operator[](std::size_t i) const {
    assert(cols() == 1);
    return (*this)(i, 0);
  }

  Matrix& operator^=(const Matrix& o) {
    assert (rows() == o.rows());
    assert (cols() == o.cols());

    if (transposed == o.transposed) {
      vals ^= o.vals;
    } else {
      for (std::size_t i = 0; i < rows(); ++i) {
        for (std::size_t j = 0; j < cols(); ++j) {
          (*this)(i, j) ^= o(i, j);
        }
      }
    }
    return *this;
  }

  Matrix& operator&=(const Matrix& o) {
    assert (rows() == o.rows());
    assert (cols() == o.cols());

    if (transposed == o.transposed) {
      vals &= o.vals;
    } else {
      for (std::size_t i = 0; i < rows(); ++i) {
        for (std::size_t j = 0; j < cols(); ++j) {
          (*this)(i, j) &= o(i, j);
        }
      }
    }
    return *this;
  }

  Matrix operator^(const Matrix& o) const {
    Matrix out = *this;
    out ^= o;
    return out;
  }

  Matrix operator&(const Matrix& o) const {
    Matrix out = *this;
    out &= o;
    return out;
  }

  void transpose() { transposed = !transposed; }

  Matrix outer_product(const Matrix& o) const {
    assert (cols() == 1);
    assert (o.cols() == 1);

    const auto n = rows();
    const auto m = o.rows();

    Matrix out(n, m);
    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t j = 0; j < m; ++j) {
        out(i, j) = (*this)[i] & o[j];
      }
    }
    return out;
  }

  const std::size_t rows() const { return transposed ? m : n; }
  const std::size_t cols() const { return transposed ? n : m; }

  friend std::ostream& operator<<(std::ostream& os, const Matrix& m) {
    for (std::size_t i = 0; i < m.rows(); ++i) {
      for (std::size_t j = 0; j < m.cols(); ++j) {
        os << m(i, j);
      }
      os << '\n';
    }
    return os;
  }

private:
  bool transposed;
  std::size_t n;
  std::size_t m;

  boost::dynamic_bitset<> vals;
};


template <typename F>
Matrix truth_table(F f, std::size_t n, std::size_t m) {
  Matrix table(m, 1 << n);
  for (std::size_t j = 0; j < (1 << n); ++j) {
    auto inp = Matrix::vector(n);
    for (std::size_t i = 0; i < n; ++i) {
      if ((1 << i) & j) { inp[i] = 1; }
    }

    const auto out = f(inp);

    for (std::size_t i = 0; i < m; ++i) {
      table(i, j) = out[i];
    }
  }
  return table;
}


inline Matrix identity_table(std::size_t n) {
  return truth_table([](const auto& x) { return x; }, n, n);
}


#endif
