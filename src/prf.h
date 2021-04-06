#ifndef PRF_H__
#define PRF_H__


#include <immintrin.h>
#include <bitset>
#include <array>


std::bitset<128> rand_key();


struct PRF {
public:
  static constexpr std::size_t nrounds = 10;

  PRF();
  PRF(std::bitset<128>);
  std::bitset<128> operator()(std::bitset<128> inp) const;

private:
  std::array<__m128i, nrounds+1> key;
};


#endif

