#include "share.h"
#include "prg.h"

#include <vector>


std::size_t ptr = 0;
std::vector<std::bitset<128>> messages;


PRG g;


template<> void Share<Mode::G>::send() const {
  messages.push_back(val);
}


template<> Share<Mode::E> Share<Mode::E>::recv() {
  return messages[ptr++];
}


template<> Share<Mode::G> Share<Mode::G>::ginput(bool b) {
  const auto out = g();
  (Share<Mode::G>::constant(b) ^ out).send();
  return out;
}


template<> Share<Mode::E> Share<Mode::E>::ginput(bool b) {
  return Share<Mode::E>::recv();
}
