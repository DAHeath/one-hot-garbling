#include "share_matrix.h"
#include "non_blackbox_gf256.h"
#include "integer.h"
#include "unary_outer_product.h"
#include "ferret.h"
#include "net_link.h"
#include "measure_link.h"
#include "standard_sbox.h"
#include "standard_mul_gf256.h"

#include <thread>
#include <iostream>
#include <chrono>


thread_local MeasureLink<GT::NetLink>* party_link;


std::size_t reps = 1000;
bool naive = false;


template <Mode mode>
ShareMatrix<mode> test_integer_exp() {

  std::uint32_t x = 13;
  const auto y = ShareMatrix<mode>::constant(from_uint32(15));

  for (std::size_t i = 0; i < reps; ++i) {
    if (naive) {
      naive_exponent<mode>(x, y);
    } else {
      exponent<mode>(x, y);
    }
  }

  return { };
}


template <Mode mode>
ShareMatrix<mode> test_integer_modp() {

  const auto x = ShareMatrix<mode>::constant(from_uint32(3215152));

  for (std::size_t i = 0; i < reps; ++i) {
    if (naive) {
      naive_mod_p<mode>(x);
    } else {
      mod_p<mode>(x);
    }
  }

  return { };
}


template <Mode mode>
ShareMatrix<mode> test_integer_mul() {
  /* const auto x = ShareMatrix<mode>::constant(from_uint64(0)); */
  /* const auto y = ShareMatrix<mode>::constant(from_uint64(0)); */
  const auto x = ShareMatrix<mode>::constant(from_uint32(0));
  const auto y = ShareMatrix<mode>::constant(from_uint32(0));
  MatrixView<const Share<mode>> xx = x;
  MatrixView<const Share<mode>> yy = y;

  for (std::size_t i = 0; i < reps; ++i) {
    if (naive) {
      naive_integer_multiply<mode>(x, y);
    } else {
      integer_multiply<mode>(xx, yy);
    }
  }


  return { };
}

template <Mode mode>
ShareMatrix<mode> test_aes() {
  const auto x = ShareMatrix<mode>::constant(byte_to_vector(1));

  for (std::size_t i = 0; i < reps; ++i) {
    if (naive) {
      standard_aes_sbox<mode>(x);
    } else {
      aes_sbox(x);
    }
  }

  return { };
}


template <Mode mode>
ShareMatrix<mode> test_mul_gf256() {
  const auto x = ShareMatrix<mode>::constant(byte_to_vector(1));
  const auto y = ShareMatrix<mode>::constant(byte_to_vector(1));

  for (std::size_t i = 0; i < reps; ++i) {
    if (naive) {
      standard_mul_gf256<mode>(x, y);
    } else {
      mul_gf256<mode>(x, y);
    }
  }

  return { };
}


template <typename F>
auto timed(F f) {
  auto start = std::chrono::high_resolution_clock::now();
  f();
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  return elapsed.count();
}


template <Mode mode>
ShareMatrix<mode> test_outer_product() {
  for (std::size_t i = 0; i < reps; ++i) {
    auto x = ShareMatrix<mode>(chunking_factor(), 1);
    auto y = ShareMatrix<mode>(chunking_factor(), 1);

    if (naive) {
      naive_outer_product<mode>(x, y);
    } else {
      outer_product<mode>(x, y);
    }
  }

  return { };
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


void protocol() {
  PRG prg;
  const auto key = prg();
  const auto seed = prg();

  ShareMatrix<Mode::G> g;
  ShareMatrix<Mode::E> e;

  std::thread th { [&] {
    // Generator
    GT::NetLink link { nullptr, 11111 };
    MeasureLink<GT::NetLink> mlink { &link };
    party_link = &mlink;

    *the_link() = &mlink;

    Share<Mode::G>::initialize(key, seed);
    initialize_gjobs();
    g = test_outer_product<Mode::G>();
    /* g = test_integer_mul<Mode::G>(); */
    /* g = test_integer_exp<Mode::G>(); */
    /* g = test_integer_modp<Mode::G>(); */
    /* g = test_mul_gf256<Mode::G>(); */
    finalize_gjobs();

    mlink.flush();

    /* std::cout << mlink.count() << '\n'; */
  } };

  {

    // Evaluator
    GT::NetLink link { "127.0.0.1", 11111 };
    MeasureLink<GT::NetLink> mlink { &link };
    party_link = &mlink;
    *the_link() = &mlink;
    Share<Mode::E>::initialize(key, seed);
    initialize_ejobs();
    e = test_outer_product<Mode::E>();
    /* e = test_integer_mul<Mode::E>(); */
    /* e = test_integer_exp<Mode::E>(); */
    /* e = test_integer_modp<Mode::E>(); */
    /* e = test_mul_gf256<Mode::E>(); */
    finalize_ejobs();

    std::cout << "GC size in bytes: " << mlink.count() << '\n';
  }

  th.join();
  
  /* std::cout << to_uint32(decode(g, e)) << '\n'; */
}


int main(int argc, char** argv) {

  if (argc < 3) {
    std::cerr << "usage: " << argv[0] << " <test repetitions> <naive{0,1}> <outer product size>\n";
    std::exit(1);
  }

  reps = atoi(argv[1]);
  naive = atoi(argv[2]);
  chunking_factor() = atoi(argv[3]);

  std::cout << naive << ' ' << chunking_factor() << '\n';

  protocol();
}
