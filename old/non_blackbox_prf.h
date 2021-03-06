#ifndef NON_BLACKBOX_PRF_H__
#define NON_BLACKBOX_PRF_H__


#include "share.h"
#include "truth_table.h"
#include "prg.h"


template <Mode mode>
struct NonBlackboxPRF {
public:
  NonBlackboxPRF(
      std::size_t n, std::size_t m, const std::vector<TruthTable>& tables)
    : n(n), m(m), tables(tables) { }

  static NonBlackboxPRF from_seed(std::size_t n, std::size_t m, PRG&);

  void operator()(std::span<const Share<mode>>, std::span<Share<mode>>) const;

private:
  std::size_t n;
  std::size_t m;

  std::vector<TruthTable> tables;
};


#endif
