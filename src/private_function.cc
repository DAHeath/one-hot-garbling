#include "private_function.h"
#include "privacy_free_point.h"
#include "uniform_function.h"
#include "util.h"

#include <vector>


template <Mode mode>
void private_function(
    const TruthTable& f,
    std::span<const Share<mode>> sx,
    std::span<Share<mode>> sfx) {

  const auto n = f.n_inp();
  const auto m = f.n_out();

  std::cout << "AFTER POINT: " << std::dec << n_ciphertexts() << "\n";

  std::vector<Share<mode>> sux(1 << n);
  const auto point = privacy_free_point<mode>(sx, std::span { sux });

  std::cout << "AFTER POINT: " << std::dec << n_ciphertexts() << "\n";

  TruthTable rf(n, m);
  if constexpr (mode == Mode::G) {
    rf = f;
    rf = rf.linear_shuffle(point);
  }

  uniform_function<mode>(sx, sux, rf, sfx);

  if constexpr (mode == Mode::G) {
    rf.send();
  } else {
    rf = TruthTable::recv(n, m);
  }

  std::vector<Share<mode>> unpacked_output(m);
  rf.apply<mode>(sux, unpacked_output);

  package<mode>(unpacked_output, sfx);
}


template void private_function<Mode::G>(
    const TruthTable&,
    std::span<const Share<Mode::G>>,
    std::span<Share<Mode::G>>);
template void private_function<Mode::E>(
    const TruthTable&,
    std::span<const Share<Mode::E>>,
    std::span<Share<Mode::E>>);
