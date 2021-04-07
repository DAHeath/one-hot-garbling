#ifndef UTIL_H__
#define UTIL_H__


#include "share.h"
#include <span>


inline auto index(std::span<const std::bitset<128>> x, std::size_t i) {
  return x[i/128][i%128];
}


#endif
