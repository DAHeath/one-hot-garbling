#include "private_function.h"
#include "privacy_free_point.h"
#include "random_function.h"
#include "util.h"

#include <vector>


void rearrange_table(
    std::size_t n,
    std::size_t m,
    std::size_t point,
    std::span<const std::bitset<128>> f,
    std::span<std::bitset<128>> out) {

  for (std::size_t j = 0; j < m; ++j) {
    for (std::size_t i = 0; i < (1 << n); ++i) {
      const auto i0 = j * (1<<n) + i;
      const auto i1 = j * (1<<n) + (i ^ point);
      out[i0/128][i0%128] = f[i1/128][i1%128];
    }
  }
}


template <Mode mode>
void private_function(
    std::span<const std::bitset<128>> f,
    std::span<const Share<mode>> sx,
    std::span<Share<mode>> sfx) {

  const auto n = sx.size();
  const auto m = sfx.size();


  std::vector<Share<mode>> sux(1 << n);
  const auto point = privacy_free_point<mode>(sx, std::span { sux });

  std::vector<std::bitset<128>> rf;
  const auto sz = ((1 << n)*m + 127)/128;
  if constexpr (mode == Mode::G) {
    assert(f.size() == sz);
  }
  rf.resize(sz);

  random_function<mode>(point, sx, sux, rf, sfx);
  if constexpr (mode == Mode::G) {
    rearrange_table(n, m, point, f, rf);
  }

  for (std::size_t i = 0; i < sz; ++i) {
    if constexpr (mode == Mode::G) {
      (Share<Mode::G> { rf[i] }).send();
    } else {
      const auto buf = Share<Mode::E>::recv();
      rf[i] = *buf;
    }
  }


  for (std::size_t j = 0; j < m; ++j) {
    for (std::size_t i = 0; i < (1 << n); ++i) {
      sfx[j] ^= index(rf, j*(1 << n)+i) ? sux[i] : Share<mode> { 0 };
    }
  }
}


template void private_function<Mode::G>(
    std::span<const std::bitset<128>>,
    std::span<const Share<Mode::G>>,
    std::span<Share<Mode::G>>);
template void private_function<Mode::E>(
    std::span<const std::bitset<128>>,
    std::span<const Share<Mode::E>>,
    std::span<Share<Mode::E>>);
