#include "non_blackbox_prf.h"
#include "private_function.h"

constexpr std::size_t one_layer_bound = 8;
constexpr std::size_t block_size = 8;

template <Mode mode>
NonBlackboxPRF<mode> NonBlackboxPRF<mode>::from_seed(
    std::size_t n, std::size_t m, PRG& g) {
  if constexpr (mode == Mode::G) {
    if (n <= one_layer_bound) {
      return { n, m, { TruthTable::uniform(n, m, g) } };
    } else {
      std::vector<TruthTable> tables;
      tables.push_back(TruthTable::uniform(n/2, m, g));
      tables.push_back(TruthTable::uniform(n - n/2, m, g));

      for (std::size_t i = 0; i < 2*m/block_size; ++i) {
        tables.push_back(TruthTable::uniform(block_size, block_size/2, g));
      }
      return { n, m, tables };
    }
  } else {
    std::vector<TruthTable> tables;
    if (n <= one_layer_bound) {
      tables.push_back(TruthTable(n, m));
    } else {
      tables.push_back(TruthTable(n/2, m));
      tables.push_back(TruthTable(n/2, m));
      for (std::size_t i = 0; i < m/(n/4); ++i) {
        tables.push_back(TruthTable(n/2, n/4));
      }
    }
    return { n, m, tables };
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

  constexpr static std::size_t p = 2; // packing factor

  assert(inp.size() == n);
  assert(out.size() == m);

  if (n <= one_layer_bound) {
    private_function(tables[0], inp, out);
  } else {
    // pack the first layer's 2m outputs into 2m/pack labels
    std::vector<Share<mode>> layer1(2*m/p);
    std::span<Share<mode>> layer1_out0 = layer1;
    std::span<Share<mode>> layer1_out1 = layer1_out0.subspan(m/p);
    layer1_out0 = layer1_out0.subspan(0, m/p);

    private_function(tables[0], inp.subspan(0, n/2), layer1_out0);
    private_function(tables[1], inp.subspan(n/2), layer1_out1);

    std::vector<Share<mode>> layer1_unpacked(2*m);
    unpackage<mode>(layer1, layer1_unpacked);

    std::vector<Share<mode>> layer1_shuffled(2*m);
    for (std::size_t i = 0; i < m; ++i) {
      layer1_shuffled[2*i] = layer1_unpacked[i];
      layer1_shuffled[2*i+1] = layer1_unpacked[i + m];
    }

    std::cout << tables.size() << '\n';

    std::vector<Share<mode>> layer2(m/p);
    std::span<Share<mode>> layer2_inp = layer1_shuffled;
    std::span<Share<mode>> layer2_out = layer2;
    for (std::size_t i = 0; i < m/(n/4); ++i) {
      private_function<mode>(tables[2 + i], layer2_inp.subspan(i*n/2, n/2), layer2_out.subspan(i*n/(4*p), n/(4*p)));
    }
    unpackage<mode>(layer2_out, out);
  }
}

template void NonBlackboxPRF<Mode::G>::operator()(
    std::span<const Share<Mode::G>>,
    std::span<Share<Mode::G>>) const;
template void NonBlackboxPRF<Mode::E>::operator()(
    std::span<const Share<Mode::E>>,
    std::span<Share<Mode::E>>) const;
