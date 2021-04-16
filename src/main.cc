#include "share_matrix.h"
#include "non_blackbox_gf256.h"
#include "integer.h"


template <Mode mode>
ShareMatrix<mode> test_integer() {
  return integer_multiply(
      ShareMatrix<mode>::constant(from_uint32(5555555)),
      ShareMatrix<mode>::constant(from_uint32(44445)));
}


template <Mode mode>
ShareMatrix<mode> test_matrix() {
  auto x = ShareMatrix<mode>(128, 128);
  auto y = ShareMatrix<mode>(128, 128);

  x(5, 1) = Share<mode>::ginput(true);
  y(10, 1) = Share<mode>::ginput(true);

  return x * y;
}


template <Mode mode>
ShareMatrix<mode> test_mul_gf256() {
  auto x = ShareMatrix<mode>(8, 1);
  auto y = ShareMatrix<mode>(8, 1);

  x[0] = Share<mode>::bit(true);
  x[1] = Share<mode>::bit(true);

  y[0] = Share<mode>::bit(true);
  y[7] = Share<mode>::bit(true);

  return aes_sbox<mode>(x);
}


Matrix decode(const ShareMatrix<Mode::G>& g, const ShareMatrix<Mode::E>& e) {
  const auto n = g.rows();
  const auto m = g.cols();

  Matrix out(n, m);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < m; ++j) {
      if (((*g(i, j)) ^ (*e(i, j))) == 0) {
        out(i, j) = 0;
      } else if ((*g(i, j) ^ *(e(i, j))) == Share<Mode::G>::delta()) {
        out(i, j) = 1;
      } else {
        std::cerr << "ERROR: BAD LABEL!\n";
        std::cerr << *g(i, j) << '\n';
        std::cerr << *e(i, j) << '\n';
        std::cerr << Share<Mode::G>::delta() << '\n';
        std::exit(1);
      }
    }
  }
  return out;
}


int main() {
  PRG prg;
  const auto key = prg();
  const auto seed = prg();

  Share<Mode::G>::initialize(key, seed);
  Share<Mode::E>::initialize(key, seed);

  const auto g = test_matrix<Mode::G>();
  const auto e = test_matrix<Mode::E>();

  std::cout << decode(g, e) << '\n';
  std::cout << std::dec << n_ciphertexts() << "\n";

  /* std::cout << byte_to_vector(invert_gf256(129)) << '\n'; */
}
