#include "share_matrix.h"
#include "non_blackbox_gf256.h"
#include "integer.h"
#include "unary_outer_product.h"

#include <iostream>


template <Mode mode>
ShareMatrix<mode> test_integer() {
  return integer_multiply(
      ShareMatrix<mode>::constant(from_uint32(5555555)),
      ShareMatrix<mode>::constant(from_uint32(44445)));
}


template <Mode mode>
ShareMatrix<mode> test_matrix() {
  constexpr std::size_t n = 4;

  auto x = ShareMatrix<mode>(n, n);
  for (std::size_t i = 0; i < n; ++i) {
    x(i, i) = Share<mode>::ginput(true);
  }

  auto y = ShareMatrix<mode>(n, n);
  y(3, 3) = Share<mode>::ginput(true);
  y(1, 3) = Share<mode>::ginput(true);

  return x*y;
  /* return x * y; */
  /* return naive_matrix_multiplication<mode>(x, y); */
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


int main() {
  PRG prg;
  const auto key = prg();
  const auto seed = prg();

  Share<Mode::G>::initialize(key, seed);
  Share<Mode::E>::initialize(key, seed);

  initialize_gjobs();
  const auto g = test_matrix<Mode::G>();
  finalize_jobs();

  initialize_ejobs();
  const auto e = test_matrix<Mode::E>();
  finalize_jobs();

  std::cout << decode(g, e) << '\n';
  std::cout << std::dec << n_ciphertexts() << "\n";

  /* std::cout << byte_to_vector(invert_gf256(129)) << '\n'; */
}
