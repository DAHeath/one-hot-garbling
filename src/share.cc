#include "share.h"
#include "prg.h"

#include <vector>


std::vector<Share<Mode::G>> deltas;

PRG prg;

template<>
std::bitset<128> Share<Mode::G>::delta() { return *deltas[1]; }


template<> void Share<Mode::G>::initialize(std::bitset<128> fixed_key, std::bitset<128> seed) {
  Share<Mode::G>::nonce = 0;
  Share<Mode::G>::fixed_key = fixed_key;
  prg = seed;
  auto delta = prg();
  delta[0] = 1;
  for (std::size_t i = 1; i < 8; ++i) { delta[i] = 0; }

  // instantiate each multiple of delta for simplicity
  deltas.resize(1 << 8);
  std::uint8_t inp[16];
  memcpy(inp, &delta, 16);
  std::uint8_t out[16];
  for (std::size_t i = 0; i < (1 << 8); ++i) {
    for (std::size_t j = 0; j < 16; ++j) {
      out[j] = mul_gf256(inp[j], i);
    }
    memcpy(deltas.data() + i, out, 16);
  }
}


template<> void Share<Mode::E>::initialize(std::bitset<128> fixed_key, std::bitset<128> seed) {
  Share<Mode::E>::nonce = 0;
  Share<Mode::E>::fixed_key = fixed_key;
  prg = seed;
}


template<> Share<Mode::E> Share<Mode::E>::bit(bool b) {
  return { 0 };
}

template<> Share<Mode::G> Share<Mode::G>::bit(bool b) {
  if (b) {
    return { deltas[1] };
  } else {
    return { 0 };
  }
}


std::size_t ptr = 0;
std::vector<std::bitset<128>> messages;


std::size_t n_ciphertexts() {
  return messages.size();
}


template<> void Share<Mode::G>::send() const {
  messages.push_back(val);
}


template<> Share<Mode::E> Share<Mode::E>::recv() {
  return messages[ptr++];
}


template<> Share<Mode::G> Share<Mode::G>::ginput(bool b) {
  const auto out = prg();
  (Share<Mode::G>::bit(b) ^ out).send();
  return out;
}


template<> Share<Mode::E> Share<Mode::E>::ginput(bool b) {
  return Share<Mode::E>::recv();
}


template<> Share<Mode::G> Share<Mode::G>::uniform() {
  const auto b = prg()[0];
  return Share<Mode::G>::bit(b);
}


template<> Share<Mode::E> Share<Mode::E>::uniform() {
  return Share<Mode::E>::bit(false);
}


template <Mode mode>
void Share<mode>::unpack(std::span<Share<mode>> out) const {
  const auto n = out.size();
  if (n == 1) {
    // Base case: the word contains exactly one bit.
    out[0] ^= *this;
  } else {
    // In the general case, split the word into its low bits and its high bits
    // (using garbled rows) and then recursively unpack the two halves.
    const auto half = (n+1)/2;

    const std::size_t mask = (1 << n) - 1;
    const std::size_t half_mask = (1 << half) - 1;
    const auto col = color256() & mask;

    Share<mode> low;
    if constexpr (mode == Mode::G) {
      // G constructs 2^n-1 garbled rows such that each row encrypts the low
      // bits corresponding to the current word.
      // We use color to permute the rows, as is typical.
      const auto lastrow = col ^ mask;
      const auto X = ((*this) ^ deltas[lastrow]).H() ^ deltas[lastrow & half_mask];

      for (std::uint8_t i = 0; i < ((1 << n) - 1); ++i) {
        const auto row = i ^ col;

        (((*this) ^ deltas[row]).H() ^ X ^ deltas[row & half_mask]).send();
      }

      low = X;
    } else {
      std::vector<Share<Mode::E>> rows(1 << n);
      for (std::size_t i = 0; i < ((1 << n) - 1); ++i) {
        rows[i] = Share<Mode::E>::recv();
      }
      rows[(1 << n) - 1] = Share<Mode::E>::bit(false); // By GRR, last row is 0.
      low = (*this).H() ^ rows[col];
    }
    ++Share<mode>::nonce;
    // High bits can be computed by subtracting off low bits.
    // Now with the high bits in a separate word, we move the high bits to the
    // low bits by scaling by the appropriate magic number.
    const auto magic = invert_gf256(1 << (n/2));
    const auto high = ((*this) ^ low).scale(magic);

    // Recursively unpack the two words.
    auto out0 = out.subspan(0, half);
    auto out1 = out.subspan(half);
    low.unpack(out0);
    high.unpack(out1);
  }
}

template void Share<Mode::G>::unpack(std::span<Share<Mode::G>>) const;
template void Share<Mode::E>::unpack(std::span<Share<Mode::E>>) const;
