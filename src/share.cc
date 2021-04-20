#include "share.h"
#include "prg.h"
#include "link.h"

#include <vector>
#include <iomanip>
#include <iostream>


std::bitset<128> delta;

PRG prg;


template<> void Share<Mode::G>::initialize(std::bitset<128> fixed_key, std::bitset<128> seed) {
  Share<Mode::G>::nonce = 0;
  Share<Mode::G>::fixed_key = fixed_key;
  prg = seed;
  delta = prg();
  delta[0] = 1;
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
    return { delta };
  } else {
    return { 0 };
  }
}



thread_local Link* link;

Link** the_link() { return &link; }


std::size_t ptr = 0;
std::vector<std::bitset<128>> messages;


std::size_t n_ciphertexts() {
  return 0;
  /* return messages.size(); */
}


template<> void Share<Mode::G>::send() const {
  static std::array<std::byte, 16> buffer;
  memcpy(buffer.data(), &val, 16);
  link->send(buffer);
  /* messages.push_back(val); */
}


template<> Share<Mode::E> Share<Mode::E>::recv() {
  static std::array<std::byte, 16> buffer;
  link->recv(buffer);
  Share<Mode::E> out;
  memcpy(&out.val, buffer.data(), 16);
  return out;
  /* return messages[ptr++]; */
}


std::size_t reveal_ptr = 0;
std::vector<bool> revelations;

template<Mode mode>
void Share<mode>::reveal() {
  if constexpr (mode == Mode::G) {
    revelations.push_back(color());
    val[0] = 0;
  } else {
    val[0] = val[0] ^ revelations[reveal_ptr++];
  }
}

template void Share<Mode::G>::reveal();
template void Share<Mode::E>::reveal();



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

template<> Share<Mode::E>& Share<Mode::E>::operator&=(const Share<Mode::E>& o) {
  const auto A = *this;
  const auto B = o;

  const auto a = A.color();
  const auto b = B.color();

  const auto zero = Share<Mode::E>::bit(0);

  // E gate
  const auto e_row = Share<Mode::E>::recv();
  const auto X = A.H() ^ (a ? (e_row ^ B) : zero);
  ++Share<Mode::E>::nonce;

  // G gate
  const auto g_row = Share<Mode::E>::recv();
  const auto Y = B.H() ^ (b ? g_row : zero);
  ++Share<Mode::E>::nonce;

  *this = X ^ Y;
  return *this;
}

template<> Share<Mode::G>& Share<Mode::G>::operator&=(const Share<Mode::G>& o) {
  const auto A = *this;
  const auto B = o;

  const auto a = A.color();
  const auto b = B.color();

  const auto zero = Share<Mode::G>::bit(0);
  const auto one = Share<Mode::G>::bit(1);

  // E gate
  const auto X = (A ^ (a ? one : zero)).H();
  const auto e_row = (A ^ (a ? zero : one)).H() ^ X ^ B;
  ++Share<Mode::G>::nonce;
  e_row.send();

  // G gate
  const auto Y = (B ^ (b ? one : zero)).H() ^ ((a && b) ? one : zero);
  const auto g_row = (B ^ (b ? zero : one)).H() ^ Y ^ ((a && !b) ? one : zero);
  ++Share<Mode::G>::nonce;
  g_row.send();

  *this = X ^ Y;
  return *this;
}


template <Mode mode>
std::ostream& operator<<(std::ostream& os, const Share<mode> s) {
  std::uint64_t xs[2];
  memcpy(xs, &(*s), 16);
  os << std::setfill('0') << std::setw(16) << std::right << std::hex << xs[0];
  os << std::setfill('0') << std::setw(16) << std::right << std::hex << xs[1];
  return os;
}

template std::ostream& operator<<(std::ostream&, const Share<Mode::G>);
template std::ostream& operator<<(std::ostream&, const Share<Mode::E>);


bool decode(const Share<Mode::G>& g, const Share<Mode::E>& e) {
  if ((*g ^ *e) == 0) {
    return false;
  } else if ((*g ^ *e) == delta) {
    return true;
  } else {
    std::cerr << "ERROR: BAD LABEL!\n";
    std::cerr << g << '\n';
    std::cerr << e << '\n';
    std::cerr << delta << '\n';
    std::exit(1);
  }
}
