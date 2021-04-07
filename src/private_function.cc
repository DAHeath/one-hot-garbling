#include "private_function.h"
#include "privacy_free_point.h"
#include "random_function.h"
#include "util.h"

#include <vector>


template <Mode mode>
void private_function(
    const TruthTable& f,
    std::span<const Share<mode>> sx,
    std::span<Share<mode>> sfx) {

  const auto n = sx.size();
  const auto m = sfx.size();


  std::vector<Share<mode>> sux(1 << n);
  const auto point = privacy_free_point<mode>(sx, std::span { sux });

  TruthTable rf(n, m);

  random_function<mode>(sx, sux, rf, sfx);
  if constexpr (mode == Mode::G) {
    rf = rf.linear_shuffle(point);
  }

  if constexpr (mode == Mode::G) {
    rf.send();
  } else {
    rf = TruthTable::recv(n, m);
  }

  rf.apply<mode>(sux, sfx);
}


template void private_function<Mode::G>(
    const TruthTable&,
    std::span<const Share<Mode::G>>,
    std::span<Share<Mode::G>>);
template void private_function<Mode::E>(
    const TruthTable&,
    std::span<const Share<Mode::E>>,
    std::span<Share<Mode::E>>);
