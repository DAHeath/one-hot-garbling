#include "truth_table.h"
#include "prg.h"

#include <vector>


std::size_t truth_ptr;
std::vector<std::uint8_t> table_messages;

std::size_t n_table_ciphertexts() {
  return (table_messages.size() + 15) / 16;
}


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
