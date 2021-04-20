#include "share_matrix.h"
#include "non_blackbox_gf256.h"
#include "integer.h"
#include "unary_outer_product.h"
#include "ferret.h"
#include "net_link.h"
#include "measure_link.h"

#include <thread>
#include <iostream>


template <Mode mode>
ShareMatrix<mode> test_integer() {
  const auto x = ShareMatrix<mode>::constant(from_uint32(5555555));
  const auto y = ShareMatrix<mode>::constant(from_uint32(44445));
  MatrixView<const Share<mode>> xx = x;
  MatrixView<const Share<mode>> yy = y;

  return naive_integer_multiply(xx, yy);
  /* return integer_multiply(xx, yy); */
}


template <Mode mode>
ShareMatrix<mode> test_outer_product() {
  constexpr std::size_t n = 128;
  auto x = ShareMatrix<mode>(n, 1);
  auto y = ShareMatrix<mode>(n, 1);


  return outer_product<mode>(x, y);
}


template <Mode mode>
ShareMatrix<mode> test_matrix_multiplication() {
  constexpr std::size_t n = 128;

  auto x = ShareMatrix<mode>(n, n);
  auto y = ShareMatrix<mode>(n, n);

  MatrixView<const Share<mode>> xx = x;
  MatrixView<const Share<mode>> yy = y;

  return xx*yy;
  /* return naive_matrix_multiplication(xx,yy); */
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


constexpr std::size_t experiment_repetitions = 200;


void protocol() {
  PRG prg;
  const auto key = prg();
  const auto seed = prg();

  ShareMatrix<Mode::G> g;
  ShareMatrix<Mode::E> e;

  std::thread th { [&] {
    // Generator
    NetLink link { nullptr, 11111 };
    MeasureLink<NetLink> mlink { &link };

    *the_link() = &mlink;

    Share<Mode::G>::initialize(key, seed);
    initialize_gjobs();
    for (std::size_t i = 0; i < experiment_repetitions; ++i) {
      g = test_integer<Mode::G>();
    }
    finalize_gjobs();

    mlink.flush();

    std::cout << mlink.count() << '\n';
  } };

  {

    // Evaluator
    NetLink link { "127.0.0.1", 11111 };
    MeasureLink<NetLink> mlink { &link };
    *the_link() = &mlink;
    Share<Mode::E>::initialize(key, seed);
    initialize_ejobs();
    for (std::size_t i = 0; i < experiment_repetitions; ++i) {
      e = test_integer<Mode::E>();
    }
    finalize_ejobs();

    std::cout << mlink.count() << '\n';
  }

  th.join();
  
  std::cout << to_uint32(decode(g, e)) << '\n';
}


int main() {
  protocol();
}
