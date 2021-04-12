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

  auto operator()(std::size_t i, std::size_t j) {
    if (transposed) { return get(j, i); } else { return get(i, j); }
  }

  auto operator()(std::size_t i, std::size_t j) const {
    if (transposed) { return get(j, i); } else { return get(i, j); }
  }

  void transpose() { transposed = !transposed; }

  const std::size_t rows() const { return transposed ? m : n; }
  const std::size_t cols() const { return transposed ? n : m; }

  friend std::ostream& operator<<(std::ostream& os, const Matrix& m) {
    for (std::size_t i = 0; i < m.rows(); ++i) {
      for (std::size_t j = 0; j < m.cols(); ++j) {
        os << m(i, j);
      }
      os << '\n';
    }
    os << '\n';
    return os;
  }

private:
  bool transposed;
  std::size_t n;
  std::size_t m;

  boost::dynamic_bitset<> vals;
};


#endif
