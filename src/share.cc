#include "share.h"
#include "prg.h"

#include <vector>


std::vector<Share<Mode::G>> deltas;

PRG prg;



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


template<> Share<Mode::G> Share<Mode::G>::random() {
  const auto b = prg()[0];
  return Share<Mode::G>::bit(b);
}


template<> Share<Mode::E> Share<Mode::E>::random() {
  return Share<Mode::E>::bit(false);
}


template <Mode mode>
void Share<mode>::unpack(std::span<Share<mode>> out) const {
  // Every field element of form abcd0000 maps to an element of form
  // 0000wxyz when multiplied by 116.
  constexpr static std::uint8_t magic8 = 116;
  // Every field element of form 00abc000 maps to an element of form
  // 00000xyz when multiplied by 232.
  constexpr static std::uint8_t magic6 = 232;
  // Every field element of form 0000ab00 maps to an element of form
  // 000000xy when multiplied by 203.
  constexpr static std::uint8_t magic4 = 203;
  // Every field element of form 000000a0 maps to an element of form
  // 0000000x when multiplied by 141.
  constexpr static std::uint8_t magic2 = 141;

  constexpr static std::uint8_t magic[9] =
    { 0, 1, magic2, magic4, magic4, magic6, magic6, magic8, magic8 };

  const auto n = out.size();
  if (n == 1) {
    out[0] = *this;
  } else {

    std::cout << n << '\n';

    const auto half = (n+1)/2;

    const std::size_t mask = (1 << n) - 1;
    const std::size_t half_mask = (1 << half) - 1;
    const auto col = color256() & mask;

    std::cout << col << '\n';

    Share<mode> low;

    if constexpr (mode == Mode::G) {
      const auto lastrow = col ^ mask;
      const auto X = ((*this) ^ deltas[lastrow]).H() ^ deltas[lastrow & half_mask];

      for (std::uint8_t i = 0; i < ((1 << n) - 1); ++i) {
        const auto row = i ^ col;

        (((*this) ^ deltas[row]).H() ^ X ^ deltas[row & half_mask]).send();
      }

      low = X;
    } else {
      std::vector<Share<Mode::E>> rows((1 << n) - 1);
      for (auto& row: rows) {
        row = Share<Mode::E>::recv();
      }
      low = (*this).H() ^ (col == ((1 << n) - 1) ? Share<Mode::E>(0) : rows[col]);
    }
    ++Share<mode>::nonce;
    const auto high = ((*this) ^ low).scale(magic[n]);
    auto out0 = out.subspan(0, half);
    auto out1 = out.subspan(half);

    low.unpack(out0);
    high.unpack(out1);
  }
}

template void Share<Mode::G>::unpack(std::span<Share<Mode::G>>) const;
template void Share<Mode::E>::unpack(std::span<Share<Mode::E>>) const;

