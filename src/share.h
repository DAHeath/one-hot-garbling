#ifndef SHARE_H__
#define SHARE_H__


#include "mode.h"
#include "prg.h"
#include "link.h"

#include <bitset>
#include <span>
#include <ostream>


std::size_t n_ciphertexts();


Link** the_link();


template <Mode mode>
struct Share {
public:
  static inline std::size_t nonce;
  static inline PRF fixed_key;

  static void initialize(std::bitset<128> fixed_key, std::bitset<128> seed);

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


  constexpr Share operator==(const Share& o) const { return ~((*this) ^ o); }
  constexpr Share operator!=(const Share& o) const { return (*this) ^ o; }
  Share operator|(const Share& o) const { return ~(~(*this) & ~o); }

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

  const std::bitset<128>& operator*() const { return val; }
  std::bitset<128>& operator*() { return val; }

private:
  std::bitset<128> val;
};


template <Mode mode>
std::ostream& operator<<(std::ostream&, const Share<mode>);

bool decode(const Share<Mode::G>&, const Share<Mode::E>&);


#endif
