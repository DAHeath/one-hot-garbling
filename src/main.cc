#include "share_matrix.h"
#include "non_blackbox_gf256.h"
#include "integer.h"
#include "unary_outer_product.h"
#include "ferret.h"
#include "net_link.h"

#include <thread>
#include <iostream>


template <Mode mode>
ShareMatrix<mode> test_integer() {
  return integer_multiply(
      ShareMatrix<mode>::constant(from_uint32(5555555)),
      ShareMatrix<mode>::constant(from_uint32(44445)));
}


template <Mode mode>
ShareMatrix<mode> test_matrix() {
  constexpr std::size_t n = 256;

  auto x = ShareMatrix<mode>(n, n);
  for (std::size_t i = 0; i < n; ++i) {
    x(i, i) = Share<mode>::ginput(true);
  }

  auto y = ShareMatrix<mode>(n, n);
  y(3, 3) = Share<mode>::ginput(true);
  y(1, 3) = Share<mode>::ginput(true);

  MatrixView<const Share<mode>> xx = x;
  MatrixView<const Share<mode>> yy = y;

  return xx*yy;
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


void protocol() {

  std::thread th { [] {
    // Generator
    NetLink link { nullptr, 11111 };
    *the_link() = &link;

    initialize_gjobs();
    test_matrix<Mode::G>();
    finalize_gjobs();

    link.flush();
  } };

  {
    // Evaluator
    NetLink link { "127.0.0.1", 11111 };
    *the_link() = &link;
    initialize_ejobs();
    test_matrix<Mode::E>();
    finalize_ejobs();
  }


  th.join();
  /* const char* m = (const char*)message.data(); */
  /* std::cout << m << '\n'; */
}


int main() {
  protocol();
  /* PRG prg; */
  /* const auto key = prg(); */
  /* const auto seed = prg(); */

  /* Share<Mode::G>::initialize(key, seed); */
  /* Share<Mode::E>::initialize(key, seed); */

  /* initialize_gjobs(); */
  /* const auto g = test_matrix<Mode::G>(); */
  /* finalize_jobs(); */

  /* initialize_ejobs(); */
  /* const auto e = test_matrix<Mode::E>(); */
  /* finalize_jobs(); */

  /* std::cout << decode(g, e) << '\n'; */
  /* std::cout << std::dec << n_ciphertexts() << "\n"; */

  /* std::cout << byte_to_vector(invert_gf256(129)) << '\n'; */
}
