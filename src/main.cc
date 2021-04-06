#include "privacy_free_point.h"
#include "prg.h"
#include <vector>
#include <iostream>


template <Mode mode>
void test() {

  std::vector<Share<mode>> x(3);
  std::vector<Share<mode>> point(8);
  x[0] = Share<mode>::ginput(false);
  x[1] = Share<mode>::ginput(true);
  x[2] = Share<mode>::ginput(true);

  std::cout << privacy_free_point<mode>(x, point) << '\n';

  for (auto p: point) {
    std::cout << p;
    if constexpr (mode == Mode::G) {
      std::cout << ' ' << (p ^ Share<mode>::constant(true));
    }
    std::cout << '\n';
  }
  std::cout << '\n';
}


int main() {
  Share<Mode::G>::delta = PRG()() | std::bitset<128> { 1 };
  Share<Mode::G>::nonce = 0;
  Share<Mode::G>::fixed_key = PRF();

  Share<Mode::E>::delta = 0;
  Share<Mode::E>::nonce = 0;
  Share<Mode::E>::fixed_key = Share<Mode::G>::fixed_key;

  test<Mode::G>();
  test<Mode::E>();
}
