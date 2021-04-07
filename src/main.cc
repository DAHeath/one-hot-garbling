#include "privacy_free_point.h"
#include "private_function.h"
#include "prg.h"
#include <vector>
#include <iostream>


template <Mode mode>
void test_eval() {
  std::vector<Share<mode>> x(3);
  std::vector<Share<mode>> fx(2);
  x[0] = Share<mode>::ginput(false);
  x[1] = Share<mode>::ginput(false);
  x[2] = Share<mode>::ginput(false);


  std::vector<std::bitset<128>> f(8);
  f[0][0] = false;
  f[0][1] = true;
  f[0][2] = false;
  f[0][3] = false;
  f[0][4] = false;
  f[0][5] = false;
  f[0][6] = false;
  f[0][7] = false;

  f[0][8] = true;
  f[0][9] = false;
  f[0][10] = false;
  f[0][11] = false;
  f[0][12] = false;
  f[0][13] = false;
  f[0][14] = false;
  f[0][15] = false;


  private_function<mode>(f, x, fx);

  for (auto p: fx) {
    std::cout << p;
    if constexpr (mode == Mode::G) {
      std::cout << ' ' << (p ^ Share<mode>::constant(true));
    }
    std::cout << '\n';
  }
  std::cout << '\n';

}


template <Mode mode>
void test_point() {

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

  test_eval<Mode::G>();
  test_eval<Mode::E>();
}
