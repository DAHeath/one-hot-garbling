#include "share_matrix.h"

std::size_t c_factor = 6;

std::size_t& chunking_factor() {
  return c_factor;
}
