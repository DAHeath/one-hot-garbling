#ifndef SHARE_H__
#define SHARE_H__


#include "mode.h"
#include "prf.h"
#include <bitset>
#include <iostream>
#include <iomanip>


inline std::bitset<128> zeroone() {
  constexpr std::uint64_t half_zo = 0x5555555555555555;
  std::bitset<128> zo = half_zo;
  zo <<= 64;
  zo |= half_zo;
  return zo;
}


inline std::bitset<128> mulgf4(const std::bitset<128>& a, const std::bitset<128>& b) {
  const static std::bitset<128> mask = zeroone();
  const auto a0 = a & mask;
  const auto a1 = (a & ~mask) >> 1;
  const auto b0 = b & mask;
  const auto b1 = (b & ~mask) >> 1;
  const auto a1b1 = a1 & b1;

  return (((a1&b0) ^ (a0&b1) ^ a1b1) << 1) ^ (a1b1 ^ (a0&b0));
}



template <Mode mode>
struct Share {
public:
  static inline std::bitset<128> delta;
  static inline std::bitset<128> delta_gf4;
  static inline std::size_t nonce;
  static inline PRF fixed_key;

  constexpr Share() { }
  constexpr Share(std::bitset<128> val) : val(val) { }
  constexpr static Share bit(bool b) {
    Share out;
    if constexpr (mode == Mode::E) { out.val = 0; }
    else {
      if (b) { out.val = delta; } else { out.val = 0; }
    }
    return out;
  }

  static Share gf4(std::bitset<2> x) {
    Share out;
    if constexpr (mode == Mode::E) {
      out.val = 0;
    } else {
      static const auto zo = zeroone();
      static const auto oz = zo << 1;

      if (x[0]) { out.val ^= mulgf4(delta_gf4, zo); }
      if (x[1]) { out.val ^= mulgf4(delta_gf4, oz); }
    }
    return out;
  }

  static Share random();

  constexpr Share& operator^=(const Share& o) { val ^= o.val; return *this; }
  constexpr Share operator^(const Share& o) const {
    Share out = *this;
    out ^= o;
    return out;
  }

  constexpr Share operator~() const {
    return (*this) ^ bit(true);
  }

  Share H() const {
    return fixed_key(val ^ std::bitset<128> { nonce });
  }

  Share H(std::size_t nonce) const {
    return fixed_key(val ^ std::bitset<128> { nonce });
  }

  void send() const;
  static Share recv();

  static Share ginput(bool);

  constexpr bool color() const {
    return val[0];
  }


  friend std::ostream& operator<<(std::ostream& os, const Share<mode> s) {
    std::uint64_t xs[2];
    memcpy(xs, &s.val, 16);
    os << std::setfill('0') << std::setw(16) << std::right << std::hex << xs[0];
    os << std::setfill('0') << std::setw(16) << std::right << std::hex << xs[1];
    return os;
  }

  const std::bitset<128>& operator*() const { return val; }
  std::bitset<128>& operator*() { return val; }

private:
  std::bitset<128> val;
};


#endif
