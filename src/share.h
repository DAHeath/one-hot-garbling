#ifndef SHARE_H__
#define SHARE_H__


#include "mode.h"
#include "prf.h"
#include <bitset>
#include <iostream>
#include <iomanip>


template <Mode mode>
struct Share {
public:
  static inline std::bitset<128> delta;
  static inline std::size_t nonce;
  static inline PRF fixed_key;

  constexpr Share() { }
  constexpr Share(std::bitset<128> val) : val(val) { }
  constexpr static Share constant(bool b) {
    Share out;
    if constexpr (mode == Mode::E) { out.val = 0; }
    else {
      if (b) { out.val = delta; } else { out.val = 0; }
    }
    return out;
  }

  static Share bit(bool);

  constexpr Share& operator^=(const Share& o) { val ^= o.val; return *this; }
  constexpr Share operator^(const Share& o) const {
    Share out = *this;
    out ^= o;
    return out;
  }

  constexpr Share operator~() const {
    return (*this) ^ constant(true);
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
