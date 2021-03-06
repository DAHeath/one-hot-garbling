#ifndef PRG_H__
#define PRG_H__


#include "prf.h"


struct PRG {
public:
  PRG() : nonce(0) { }
  PRG(PRF prf) : prf(std::move(prf)), nonce(0) { }
  PRG(std::bitset<128> seed) : prf(std::move(seed)), nonce(0) { }

  std::bitset<128> operator()() { return prf(nonce++); }

private:
  PRF prf;
  std::size_t nonce;
};


#endif
