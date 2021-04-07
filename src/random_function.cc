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

  const auto one = Share<mode>::bit(true);
  const auto zero = Share<mode>::bit(false);

  const auto mask = TruthTable::input_column(n, m, n - ix - 1);
  if constexpr (mode == Mode::G) {
    // Depending on the color, G should swap which key encrypts which secrets
    const auto key0 = sxi ^ (sxi.color() ? one : zero);
    const auto key1 = key0 ^ one;

    // Set up two halves of truth table according to possible labels of [|x_i|].
    const auto table0 = TruthTable::random(n, m, *key0) & ~mask;
    const auto table1 = TruthTable::random(n, m, *key1) & mask;
    r ^= table0 ^ table1;

    // Compute two possible half inner products.
    const auto prod0 = table0.template apply<mode>(sux);
    const auto prod1 = table1.template apply<mode>(sux);

    for (std::size_t j = 0; j < m; ++j) {
      // Now, encrypt the two inner products corresponding to the part of the
      // vector E does not hold.
      // Garbled row reduction technique.
      const auto X = key0.H() ^ prod1[j];
      const auto row = key1.H() ^ prod0[j] ^ X;
      row.send();
      ++Share<mode>::nonce;
      out[j] ^= prod0[j] ^ prod1[j] ^ X;
    }
  } else {
    // Take inner product of half of truth table
    (TruthTable::random(n, m, *sxi) & (sxi.color() ? mask : ~mask)).apply(sux, out);

    // Other half of inner product comes by decrypting the proper row
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
