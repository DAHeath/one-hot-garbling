#include "uniform_function.h"


template <Mode mode>
void half_uniform_function(
    std::size_t ix,
    std::span<const Share<mode>> sx,
    std::span<const Share<mode>> sux,
    TruthTable& r,
    std::span<Share<mode>> out) {

  const auto n = r.n_inp();
  const auto m = r.n_out();
  const auto p = out.size(); // number of packs

  const auto sxi = sx[ix];

  const auto one = Share<mode>::bit(true);
  const auto zero = Share<mode>::bit(false);

  const auto mask = TruthTable::input_column(n, m, n - ix - 1);
  if constexpr (mode == Mode::G) {
    // Depending on the color, G should swap which key encrypts which secrets
    const auto key0 = sxi ^ (sxi.color() ? one : zero);
    const auto key1 = key0 ^ one;

    // Set up two halves of truth table according to possible labels of [|x_i|].
    const auto table0 = TruthTable::uniform(n, m, *key0) & ~mask;
    const auto table1 = TruthTable::uniform(n, m, *key1) & mask;
    r ^= table0 ^ table1;

    // Compute two possible half inner products.
    const auto prod0 = table0.template apply<mode>(sux);
    const auto prod1 = table1.template apply<mode>(sux);

    // Package the m inner products into p shares.
    std::vector<Share<mode>> packed_prod0(p);
    std::vector<Share<mode>> packed_prod1(p);
    package<mode>(prod0, packed_prod0);
    package<mode>(prod1, packed_prod1);

    for (std::size_t i = 0; i < p; ++i) {
      // The last pack may not be full; calculate its size.
      // Now, encrypt the two inner products corresponding to the part of the
      // vector E does not hold.
      // Garbled row reduction technique.
      const auto X = key0.H() ^ packed_prod1[i];
      const auto row = key1.H() ^ packed_prod0[i] ^ X;
      row.send();
      ++Share<mode>::nonce;
      out[i] ^= packed_prod0[i] ^ packed_prod1[i] ^ X;
    }
  } else {
    // Take inner product of half of truth table
    std::vector<Share<mode>> prod(m);
    (TruthTable::uniform(n, m, *sxi) & (sxi.color() ? mask : ~mask)).template apply<mode>(sux, prod);

    // Package the m inner products into the p outputs.
    package<mode>(prod, out);

    // Other half of inner product comes by decrypting the proper row
    for (std::size_t i = 0; i < p; ++i) {
      const auto row = Share<mode>::recv();
      out[i] ^= sxi.H();
      if (sxi.color()) { out[i] ^= row; }
      ++Share<mode>::nonce;
    }
  }
}


template <Mode mode>
void uniform_function(
    std::span<const Share<mode>> sx, // [|x|]
    std::span<const Share<mode>> sux, // [|U(x)|]
    TruthTable& r,
    std::span<Share<mode>> out) {

  const auto n = r.n_inp();
  const auto m = r.n_out();

  for (std::size_t i = 0; i < n; ++i) {
    half_uniform_function(i, sx, sux, r, out);
  }

  std::vector<Share<mode>> masks(m);
  for (std::size_t j = 0; j < m; ++j) {
    const auto s = Share<mode>::uniform();
    masks[j] = s;
    if (s.color()) {
      r.flip_column(j);
    }
  }

  // Pack the m random masks on top of the p shares.
  package<mode>(masks, out);
}

template void half_uniform_function<Mode::G>(
    std::size_t,
    std::span<const Share<Mode::G>>,
    std::span<const Share<Mode::G>>,
    TruthTable&,
    std::span<Share<Mode::G>>);
template void half_uniform_function<Mode::E>(
    std::size_t,
    std::span<const Share<Mode::E>>,
    std::span<const Share<Mode::E>>,
    TruthTable&,
    std::span<Share<Mode::E>>);

template void uniform_function<Mode::G>(
    std::span<const Share<Mode::G>>,
    std::span<const Share<Mode::G>>,
    TruthTable&,
    std::span<Share<Mode::G>>);
template void uniform_function<Mode::E>(
    std::span<const Share<Mode::E>>,
    std::span<const Share<Mode::E>>,
    TruthTable&,
    std::span<Share<Mode::E>>);
