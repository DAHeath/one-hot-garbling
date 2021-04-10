#include "private_function.h"
#include "non_blackbox_prf.h"


template <Mode mode>
void test_eval() {
  constexpr std::size_t n = 8;
  constexpr std::size_t m = 8;
  constexpr std::size_t pack_factor = 8;
  constexpr std::size_t p = (m + pack_factor - 1) / pack_factor;

  /* std::cout << "HERE\n"; */
  /* PRG g; */
  /* const auto f = NonBlackboxPRF<mode>::from_seed(8, 8, g); */

  /* std::cout << "HERE\n"; */

  /* std::vector<Share<mode>> x(8); */
  /* std::vector<Share<mode>> y(8); */

  /* f(x, y); */

  std::vector<Share<mode>> x(n);
  std::vector<Share<mode>> packed_fx(p);

  for (auto& i: x) { i = Share<mode>::bit(false); }

  TruthTable f(n, m);

  private_function<mode>(f, x, packed_fx);

  std::vector<Share<mode>> fx(m);
  unpackage<mode>(packed_fx, fx);

  for (auto p: fx) {
    std::cout << p;
    if constexpr (mode == Mode::G) {
      std::cout << ' ' << (p ^ Share<mode>::bit(true));
    }
    std::cout << '\n';
  }
  std::cout << '\n';

  std::cout << std::dec << (n_ciphertexts()) << "\n";
  std::cout << std::dec << (n_ciphertexts() + n_table_ciphertexts()) << "\n";

}


int main() {
  PRG g;
  const auto key = g();
  const auto seed = g();

  Share<Mode::G>::initialize(key, seed);
  Share<Mode::E>::initialize(key, seed);

  test_eval<Mode::G>();
  test_eval<Mode::E>();


  std::cout << "HERE\n";
  std::cout << std::dec << (int)invert_gf256(16) << '\n';
  std::cout << std::dec << (int)invert_gf256(8) << '\n';
  std::cout << std::dec << (int)invert_gf256(4) << '\n';
  std::cout << std::dec << (int)invert_gf256(2) << '\n';
  std::cout << std::dec << (int)invert_gf256(1) << '\n';
  std::cout << std::dec << (int)invert_gf256(0) << '\n';
}
