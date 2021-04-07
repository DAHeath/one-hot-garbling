#include "privacy_free_point.h"
#include "private_function.h"
#include "random_function.h"
#include "prg.h"
#include <vector>
#include <iostream>
#include "truth_table.h"

#include "gf256.h"


template <Mode mode>
void test_eval() {
  std::vector<Share<mode>> x(3);
  std::vector<Share<mode>> fx(2);
  x[0] = Share<mode>::ginput(false);
  x[1] = Share<mode>::ginput(false);
  x[2] = Share<mode>::ginput(false);


  TruthTable f(3, 2);
  f.set(1, 0, true);
  f.set(0, 1, true);


  std::cout << f << '\n';


  private_function<mode>(f, x, fx);

  for (auto p: fx) {
    std::cout << p;
    if constexpr (mode == Mode::G) {
      std::cout << ' ' << (p ^ Share<mode>::bit(true));
    }
    std::cout << '\n';
  }
  std::cout << '\n';

}


std::bitset<2> gf4mult(std::bitset<2> a, std::bitset<2> b) {
  const std::bitset<2> mask = 1;
  const auto a0 = a & mask;
  const auto a1 = (a & ~mask) >> 1;
  const auto b0 = b & mask;
  const auto b1 = (b & ~mask) >> 1;
  const auto a1b1 = a1 & b1;

  return (((a1&b0) ^ (a0&b1) ^ a1b1) << 1) ^ (a1b1 ^ (a0&b0));
}


int main() {
  /* Share<Mode::G>::delta = PRG()() | std::bitset<128> { 1 }; */
  /* Share<Mode::G>::nonce = 0; */
  /* Share<Mode::G>::fixed_key = PRF(); */

  /* Share<Mode::E>::delta = 0; */
  /* Share<Mode::E>::nonce = 0; */
  /* Share<Mode::E>::fixed_key = Share<Mode::G>::fixed_key; */

  PRG g;
  const auto key = g();
  const auto seed = g();

  Share<Mode::G>::initialize(key, seed);
  Share<Mode::E>::initialize(key, seed);

  test_eval<Mode::G>();
  test_eval<Mode::E>();

  /* std::uint64_t mask = 0xaaaaaaaaaaaaaaaa; */
  /* std::bitset<128> onezero = mask; */
  /* onezero <<= 64; */
  /* onezero |= mask; */
  /* const auto zeroone = onezero >> 1; */

  /* std::cout << (int)mul_gf256(55, 55) << '\n'; */
  /* std::cout << (int)invert_gf256(98) << '\n'; */
  /* std::cout << (int)(mul_gf256(175, 98)) << '\n'; */
}
