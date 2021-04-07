#include "share.h"
#include "prg.h"

#include <vector>


std::vector<std::bitset<128>> deltas;


template<> void Share<Mode::G>::initialize(std::bitset<128> fixed_key, std::bitset<128> seed) {
  Share<Mode::G>::nonce = 0;
  Share<Mode::G>::fixed_key = fixed_key;
  Share<Mode::G>::prg = seed;
  auto delta = Share<Mode::G>::prg();
  delta[0] = 1;
  for (std::size_t i = 1; i < 8; ++i) { delta[i] = 0; }

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
  std::cout << deltas[1] << '\n';
}


template<> void Share<Mode::E>::initialize(std::bitset<128> fixed_key, std::bitset<128> seed) {
  Share<Mode::E>::nonce = 0;
  Share<Mode::E>::fixed_key = fixed_key;
  Share<Mode::E>::prg = seed;
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
