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

/* } */

template <Mode mode>
ShareMatrix<mode> test_matrix() {

  auto x = ShareMatrix<mode>::vector(10);
  auto y = ShareMatrix<mode>::vector(10);

  x[0] = Share<mode>::ginput(true);
  /* x[1] = Share<mode>::ginput(true); */
  /* x[2] = Share<mode>::ginput(false); */
  /* x[3] = Share<mode>::ginput(true); */

  /* y[0] = Share<mode>::ginput(false); */
  y[1] = Share<mode>::ginput(true);
  /* y[2] = Share<mode>::ginput(false); */

  /* const auto xy = x.unary_outer_product(y); */
  return outer_product<mode>(x, y);

  /* std::cout << x.color().outer_product(y.color()) << "\n"; */

  /* return ShareMatrix<mode>::constant(x.color().outer_product(y.color())); */

  /* return xy; */
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
  std::cout << std::dec << (n_ciphertexts() + n_table_ciphertexts()) << "\n";

  /* std::cout << identity_table(3) << '\n'; */

  /* test_eval<Mode::G>(); */
  /* test_eval<Mode::E>(); */


  /* Matrix m(3, 5); */

  /* m(2, 1) = true; */
  /* std::cout << m << '\n'; */
}
