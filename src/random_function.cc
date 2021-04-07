#include "random_function.h"


template <Mode mode>
void half_random_function(
    std::size_t ix,
    std::span<const Share<mode>> sx,
    std::span<const Share<mode>> sux,
    TruthTable& r,
    std::span<Share<mode>> out) {

  const auto n = sx.size();
  const auto m = out.size();

  const auto sxi = sx[ix];

  const auto mask = TruthTable::input_column(n, m, n - ix - 1);
  if constexpr (mode == Mode::G) {
    TruthTable table0, table1;
    if (sxi.color()) {
      table0 = TruthTable::random(n, m, *(sxi ^ Share<mode>::constant(true)));
      table0 &= ~mask;
      table1 = TruthTable::random(n, m, *sxi);
      table1 &= mask;
    } else {
      table0 = TruthTable::random(n, m, *sxi);
      table0 &= ~mask;
      table1 = TruthTable::random(n, m, *(sxi ^ Share<mode>::constant(true)));
      table1 &= mask;
    }

    r ^= table0;
    r ^= table1;


    if constexpr (mode == Mode::G) {
      TruthTable tt(n, m);
      tt ^= table0;
      tt ^= table1;
    }

    std::vector<Share<mode>> prod0(m);
    std::vector<Share<mode>> prod1(m);

    table0.template apply<mode>(sux, prod0);
    table1.template apply<mode>(sux, prod1);

    for (std::size_t j = 0; j < m; ++j) {
      out[j] ^= prod0[j] ^ prod1[j];

      if (sxi.color()) {
        const auto X = (sxi ^ Share<mode>::constant(true)).H() ^ prod1[j];
        const auto row = sxi.H() ^ prod0[j] ^ X;
        row.send();
        ++Share<mode>::nonce;
        out[j] ^= X;
      } else {
        const auto X = sxi.H() ^ prod1[j];
        const auto row = (sxi ^ Share<mode>::constant(true)).H() ^ prod0[j] ^ X;
        row.send();
        ++Share<mode>::nonce;
        out[j] ^= X;
      }
    }
  } else {
    TruthTable partial_table = TruthTable::random(n, m, *sxi);

    partial_table &= sxi.color() ? mask : ~mask;

    partial_table.apply(sux, out);


    for (std::size_t j = 0; j < m; ++j) {
      const auto row = Share<mode>::recv();
      out[j] ^= sxi.H();
      if (sxi.color()) { out[j] ^= row; }
      ++Share<mode>::nonce;
    }
  }
}


template <Mode mode>
void random_function(
    std::span<const Share<mode>> sx, // [|x|]
    std::span<const Share<mode>> sux, // [|U(x)|]
    TruthTable& r,
    std::span<Share<mode>> out) {

  const auto n = sx.size();
  const auto m = out.size();

  for (std::size_t i = 0; i < n; ++i) {
    half_random_function(i, sx, sux, r, out);
  }

  for (std::size_t j = 0; j < m; ++j) {
    const auto s = Share<mode>::random();
    out[j] ^= s;
    if (s.color()) {
      r.flip_column(j);
    }
  }
}

template void half_random_function<Mode::G>(
    std::size_t,
    std::span<const Share<Mode::G>>,
    std::span<const Share<Mode::G>>,
    TruthTable&,
    std::span<Share<Mode::G>>);
template void half_random_function<Mode::E>(
    std::size_t,
    std::span<const Share<Mode::E>>,
    std::span<const Share<Mode::E>>,
    TruthTable&,
    std::span<Share<Mode::E>>);

template void random_function<Mode::G>(
    std::span<const Share<Mode::G>>,
    std::span<const Share<Mode::G>>,
    TruthTable&,
    std::span<Share<Mode::G>>);
template void random_function<Mode::E>(
    std::span<const Share<Mode::E>>,
    std::span<const Share<Mode::E>>,
    TruthTable&,
    std::span<Share<Mode::E>>);
