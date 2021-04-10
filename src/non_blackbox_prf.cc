#include "non_blackbox_prf.h"
#include "private_function.h"

constexpr std::size_t one_layer_bound = 8;

template <Mode mode>
NonBlackboxPRF<mode> NonBlackboxPRF<mode>::from_seed(
    std::size_t n, std::size_t m, PRG& g) {
  if constexpr (mode == Mode::G) {
    if (n <= one_layer_bound) {
      return { n, m, 1, { TruthTable::uniform(n, m, g) } };
    } else {
      std::vector<TruthTable> tables;
      tables.push_back(TruthTable::uniform(n/2, m, g));
      tables.push_back(TruthTable::uniform(n/2, m, g));

      for (std::size_t i = 0; i < m/(n/4); ++i) {
        tables.push_back(TruthTable::uniform(n/2, n/4, g));
      }
    }
  } else {
    return { n, m, 1, { } };
  }
}


template NonBlackboxPRF<Mode::G> NonBlackboxPRF<Mode::G>::from_seed(
    std::size_t, std::size_t, PRG&);
template NonBlackboxPRF<Mode::E> NonBlackboxPRF<Mode::E>::from_seed(
    std::size_t, std::size_t, PRG&);

template <Mode mode>
void NonBlackboxPRF<mode>::operator()(
    std::span<const Share<mode>> inp,
    std::span<Share<mode>> out) const {
  assert(inp.size() == n);
  assert(out.size() == m);

  if (n <= one_layer_bound) {
    private_function(tables[0], inp, out);
  } else {
    std::vector<Share<mode>> layer1_out(m);
    std::span<Share<mode>> layer1_out0 = layer1_out;
    std::span<Share<mode>> layer1_out1 = layer1_out0.subspan(m/2);
    layer1_out0 = layer1_out0.subspan(0, m/2);
  }
}

template void NonBlackboxPRF<Mode::G>::operator()(
    std::span<const Share<Mode::G>>,
    std::span<Share<Mode::G>>) const;
template void NonBlackboxPRF<Mode::E>::operator()(
    std::span<const Share<Mode::E>>,
    std::span<Share<Mode::E>>) const;
