#ifndef SHARE_H__
#define SHARE_H__


#include "mode.h"
#include "gf256.h"
#include "prg.h"

#include <bitset>
#include <span>
#include <iostream>
#include <iomanip>


std::size_t n_ciphertexts();


template <Mode mode>
struct Share {
public:
  static inline std::size_t nonce;
  static inline PRF fixed_key;

  static void initialize(std::bitset<128> fixed_key, std::bitset<128> seed);
  static std::bitset<128> delta();

  constexpr Share() { }
  constexpr Share(std::bitset<128> val) : val(val) { }
  static Share<mode> bit(bool b);

  // G selects a uniform input. E does not learn the input.
  static Share uniform();

  constexpr Share& operator^=(const Share& o) { val ^= o.val; return *this; }
  constexpr Share operator^(const Share& o) const {
    Share out = *this;
    out ^= o;
    return out;
  }

  Share& operator&=(const Share&);
  Share operator&(const Share& o) const {
    Share out = *this;
    out &= o;
    return out;
  }

  constexpr Share operator~() const {
    return (*this) ^ bit(true);
  }

  // Pack up to 8 bits into a single label.
  // Packing is homomorphic.
  static Share pack(std::span<const Share> args) {
    assert(args.size() <= 8 && args.size() > 0);

    Share out = args[0];
    for (std::size_t i = 1; i < args.size(); ++i) {
      out ^= args[i].scale(1 << i);
    }
    return out;
  }

  // Unpack a packed label into individual bits.
  // For a label that packs n bits, `unpack` consumes O(2^n) garbled rows.
  void unpack(std::span<Share> out) const;

  // Scale each byte of the label by a public scalar using multiplication in GF(256).
  Share scale(std::uint8_t scalar) const {
    Share scaled;
    std::uint8_t inp[16];
    std::uint8_t out[16];
    memcpy(inp, &val, 16);
    for (std::size_t i = 0; i < 16; ++i) {
      out[i] = mul_gf256(scalar, inp[i]);
    }
    memcpy(&scaled.val, out, 16);
    return scaled;
  }

  Share H() const {
    return fixed_key(val ^ std::bitset<128> { nonce });
  }

  Share H(std::size_t nonce) const {
    return fixed_key(val ^ std::bitset<128> { nonce });
  }

  void send() const;
  static Share recv();
  

  // G reveals the value to E.
  // After this call, the color will return the semantic value.
  void reveal();

  static Share ginput(bool);

  constexpr bool color() const {
    return val[0];
  }

  constexpr std::uint8_t color256() const {
    std::uint8_t out = 0;
    for (std::size_t i = 0; i < 8; ++i) {
      out |= (val[i] << i);
    }
    return out;
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


template <Mode mode>
void package(std::span<const Share<mode>> to_pack, std::span<Share<mode>> out) {
  const auto n = to_pack.size();
  const auto m = out.size();
  const auto p = (n + m - 1) / m;
  for (std::size_t i = 0; i < m; ++i) {
    out[i] ^= Share<mode>::pack(to_pack.subspan(i*p, p));
  }
}


template <Mode mode>
void unpackage(std::span<const Share<mode>> to_unpack, std::span<Share<mode>> out) {
  const auto n = to_unpack.size();
  const auto m = out.size();
  const auto p = (m + n - 1) / n;
  for (std::size_t i = 0; i < n; ++i) {
    to_unpack[i].unpack(out.subspan(i*p, p));
  }
}


#endif
