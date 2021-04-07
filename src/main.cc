#include "privacy_free_point.h"
#include "private_function.h"
#include "random_function.h"
#include "prg.h"
#include <vector>
#include <iostream>
#include "truth_table.h"


void slice(std::size_t ix, std::size_t nn) {
  for (std::size_t i = 0; i < nn; ++i) {
    if ((i & (1 << ix)) > 0) { std::cout << 1; } else { std::cout << 0; }
  }
  std::cout << '\n';
}


template <Mode mode>
void test_eval() {
  std::vector<Share<mode>> x(3);
  std::vector<Share<mode>> fx(2);
  x[0] = Share<mode>::ginput(false);
  x[1] = Share<mode>::ginput(false);
  x[2] = Share<mode>::ginput(true);


  TruthTable f(3, 2);
  f.set(0, 0, false);
  f.set(1, 0, true);
  f.set(2, 0, false);
  f.set(3, 0, false);
  f.set(4, 0, false);
  f.set(5, 0, false);
  f.set(6, 0, false);
  f.set(7, 0, false);

  f.set(0, 1, true);
  f.set(1, 1, false);
  f.set(2, 1, false);
  f.set(3, 1, false);
  f.set(4, 1, false);
  f.set(5, 1, false);
  f.set(6, 1, false);
  f.set(7, 1, false);


  std::cout << f << '\n';


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


template <Mode mode>
void test_half_random() {
  TruthTable tt(3, 2);

  std::vector<Share<mode>> x(3);
  std::vector<Share<mode>> point(8);
  x[0] = Share<mode>::ginput(false);
  x[1] = Share<mode>::ginput(true);
  x[2] = Share<mode>::ginput(true);

  std::cout << privacy_free_point<mode>(x, point) << '\n';


  std::vector<Share<mode>> rx(2);

  /* half_random_function<mode>(0, x, point, tt, rx); */
  random_function<mode>(x, point, tt, rx);

  if constexpr (mode == Mode::G) {
    std::cout << tt << "\n";
  }

  for (auto p: rx) {
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

  /* test_half_random<Mode::G>(); */
  /* test_half_random<Mode::E>(); */
  test_eval<Mode::G>();
  test_eval<Mode::E>();


/*   for (int i = 0; i < 3; ++i) { */
/*     std::cout << TruthTable::input_column(3, 1, i) << '\n'; */
/*   } */

}
