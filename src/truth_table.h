#ifndef TRUTH_TABLE_H__
#define TRUTH_TABLE_H__

#include "prg.h"
#include "share.h"

#include <span>
#include <vector>
#include <iostream>


/**
 * Truth table representation.
 * Truth tables support arithmetic bitwise manipulations.
 * WARNING: we handle truth tables w/ minimum 3 inputs; hence 8 bit values.
 */
struct TruthTable {
  public:
    TruthTable() { }
    TruthTable(std::size_t n, std::size_t m) : n(n), m(m), nn(1 << n), val(nn/8*m) { }

    TruthTable& operator^=(const TruthTable& o) {
      assert(n == o.n);
      assert(m == o.m);

      for (std::size_t i = 0; i < val.size(); ++i) {
        val[i] ^= o.val[i];
      }
      return *this;
    }

    TruthTable& operator&=(const TruthTable& o) {
      assert(n == o.n);
      assert(m == o.m);

      for (std::size_t i = 0; i < val.size(); ++i) {
        val[i] &= o.val[i];
      }
      return *this;
    }

    TruthTable operator^(const TruthTable& o) const {
      TruthTable out = *this;
      out ^= o;
      return out;
    }

    TruthTable operator&(const TruthTable& o) const {
      TruthTable out = *this;
      out &= o;
      return out;
    }

    TruthTable operator~() const {
      TruthTable out(n, m);
      for (std::size_t i = 0; i < val.size(); ++i) {
        out.val[i] = ~val[i];
      }
      return out;
    }

    /**
     * Construct a truth table that models a particular truth table input column in every column.
     */
    static TruthTable input_column(std::size_t n, std::size_t m, std::size_t ix) {
      TruthTable tt(n, m);

      for (std::size_t i = 0; i < tt.nn; ++i) {
        tt.set(i, 0, (i & (1 << ix)) > 0);
      }

      std::size_t words_per_column = tt.nn/8;
      for (std::size_t j = 1; j < m; ++j) {
        for (std::size_t i = 0; i < tt.nn/8; ++i) {
          tt.val[j * words_per_column + i] = tt.val[i];
        }
      }
      return tt;
    }

    /**
     * Construct a uniform n input m output truth table from a PRG.
     */
    static TruthTable uniform(std::size_t n, std::size_t m, const std::bitset<128>& seed) {
      PRG g(seed);
      TruthTable tt(n, m);
      const std::size_t slices = (tt.val.size() + 15) / 16;
      for (std::size_t i = 0; i < slices-1; ++i) {
        const auto slice = g();
        memcpy(&tt.val[i*16], &slice, 16);
      }
      const auto slice = g();
      memcpy(&tt.val[(slices-1)*16], &slice, tt.val.size() % 16);
      return tt;
    }

    /**
     * Rearrange the rows of the truth table by applying a linear shift.
     * Namely, each row `i` is sent to row `i ^ delta`.
     */
    TruthTable linear_shuffle(std::size_t delta) const {
      TruthTable out(n, m);

      for (std::size_t j = 0; j < m; ++j) {
        for (std::size_t i = 0; i < nn; ++i) {
          out.set(i, j, (*this)(i ^ delta, j));
        }
      }

      return out;
    }

    // Flip all bits in column j
    void flip_column(std::size_t j) {
      const auto column_size = nn/8;
      for (std::size_t i = j*column_size; i < (j + 1)*column_size; ++i) {
        val[i] = ~val[i];
      }
    }

    /**
     * Apply the truth table to a *unary* encoding inp.
     * Result is in-place.
     * I.e., maps [|U(x)|] -> [|f(x)|]
     *
     * This uses the identity.
     * <T(f), U(x)> = f(x)
     * where `< , > denotes an inner product, `T` denotes truth table and `U` denotes unary.
     */
    template <Mode mode>
    void apply(
        std::span<const Share<mode>> inp,
        std::span<Share<mode>> out) const {
      assert(inp.size() == nn);
      assert(out.size() == m);

      for (std::size_t i = 0; i < nn; ++i) {
        for (std::size_t j = 0; j < m; ++j) {
          if ((*this)(i, j)) {
            out[j] ^= inp[i];
          }
        }
      }
    }

    /**
     * Apply truth table to unary encoding, result is out-of-place.
     */
    template <Mode mode>
    std::vector<Share<mode>> apply(std::span<const Share<mode>> inp) const {
      std::vector<Share<mode>> out(m);
      apply<mode>(inp, out);
      return out;
    }

    // Look up entry i, j
    bool operator()(std::size_t i, std::size_t j) const {
      const auto ix = (nn*j + i) / 8;
      const std::uint8_t slot = val[ix];
      const std::uint8_t bit = i % 8;
      const std::uint8_t mask = 1 << bit;
      return slot & mask;
    }

    // Set entry i, j to v
    void set(std::size_t i, std::size_t j, bool v) {
      const auto ix = (nn*j + i) / 8;
      std::uint8_t& slot = val[ix];
      const uint8_t bit = i % 8;
      const uint8_t diff = v << bit;
      slot |= diff;
    }

    friend std::ostream& operator<<(std::ostream& os, const TruthTable& tab) {
      for (std::size_t i = 0; i < tab.nn; ++i) {
        for (std::size_t j = 0; j < tab.m; ++j) {
          os << tab(i, j);
        }
        os << '\n';
      }
      return os;
    }

    void send() const;
    static TruthTable recv(std::size_t n, std::size_t m);

    constexpr std::size_t n_inp() const { return n; }
    constexpr std::size_t n_out() const { return m; }

  private:
    std::size_t n;
    std::size_t m;
    std::size_t nn;

    // the table is stored by column, then by row
    std::vector<std::uint8_t> val;
};


#endif
