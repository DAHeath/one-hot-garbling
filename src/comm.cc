#include "share.h"
#include "truth_table.h"
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
  (Share<Mode::G>::bit(b) ^ out).send();
  return out;
}


template<> Share<Mode::E> Share<Mode::E>::ginput(bool b) {
  return Share<Mode::E>::recv();
}


template<> Share<Mode::G> Share<Mode::G>::random() {
  const auto b = g()[0];
  return Share<Mode::G>::bit(b);
}


template<> Share<Mode::E> Share<Mode::E>::random() {
  return Share<Mode::E>::bit(false);
}


std::size_t truth_ptr;
std::vector<std::uint8_t> table_messages;


TruthTable TruthTable::recv(std::size_t n, std::size_t m) {
  TruthTable out(n, m);
  memcpy(out.val.data(), table_messages.data() + truth_ptr, out.val.size());
  truth_ptr += out.val.size();
  return out;
}


void TruthTable::send() const {
  const auto n = table_messages.size();
  table_messages.resize(n + val.size());
  memcpy(table_messages.data() + n, val.data(), val.size());
}
