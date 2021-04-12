#include "private_function.h"
#include "non_blackbox_prf.h"
#include "share_matrix.h"
#include "matrix.h"

template <Mode mode>
void test_eval() {
  constexpr std::size_t n = 16;
  constexpr std::size_t m = 16;
  constexpr std::size_t pack_factor = 2;
  constexpr std::size_t p = (m + pack_factor - 1) / pack_factor;

  std::vector<Share<mode>> x(n);
  std::vector<Share<mode>> packed_fx(p);

  for (auto& i: x) { i = Share<mode>::bit(false); }

  TruthTable f(n, m);

  private_function<mode>(f, x, packed_fx);

  std::vector<Share<mode>> fx(m);
  unpackage<mode>(packed_fx, fx);

  /* for (auto p: fx) { */
  /*   std::cout << p; */
  /*   if constexpr (mode == Mode::G) { */
  /*     std::cout << ' ' << (p ^ Share<mode>::bit(true)); */
  /*   } */
  /*   std::cout << '\n'; */
  /* } */
  /* std::cout << '\n'; */

  std::cout << std::dec << (n_ciphertexts() + n_table_ciphertexts()) << "\n";
}


/* template <Mode mode> */
/* void test_eval() { */
/*   constexpr std::size_t n = 4; */
/*   constexpr std::size_t m = 128; */
/*   constexpr std::size_t pack_factor = 2; */
/*   constexpr std::size_t p = (m + pack_factor - 1) / pack_factor; */

/*   PRG g; */
/*   const auto f = NonBlackboxPRF<mode>::from_seed(n, m, g); */

/*   std::vector<Share<mode>> x(n); */
/*   std::vector<Share<mode>> y(m); */

/*   f(x, y); */

/*   for (auto p: y) { */
/*     std::cout << p; */
/*     if constexpr (mode == Mode::G) { */
/*       std::cout << ' ' << (p ^ Share<mode>::bit(true)); */
/*     } */
/*     std::cout << '\n'; */
/*   } */
/*   std::cout << '\n'; */

/*   std::cout << std::dec << (n_ciphertexts() + n_table_ciphertexts()) << "\n"; */
/* } */

template <Mode mode>
void test_matrix() {

  const auto m = ShareMatrix<mode>::uniform(3, 2);

  std::cout << m << '\n';
  std::cout << m.color() << '\n';
}


int main() {
  PRG g;
  const auto key = g();
  const auto seed = g();

  Share<Mode::G>::initialize(key, seed);
  Share<Mode::E>::initialize(key, seed);

  test_matrix<Mode::G>();
  test_matrix<Mode::E>();

  /* test_eval<Mode::G>(); */
  /* test_eval<Mode::E>(); */


  /* Matrix m(3, 5); */

  /* m(2, 1) = true; */
  /* std::cout << m << '\n'; */
}
