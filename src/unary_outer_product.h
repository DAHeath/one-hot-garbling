#ifndef UNARY_OUTER_PRODUCT_H__
#define UNARY_OUTER_PRODUCT_H__


#include "share_matrix.h"

// Let c be the color of x.
// Let n be the length of x.
// Let m be the length of y.
// Let f be a function from n bits to l bits
//
// This computes [|T(f) * (U(x + c) & y) |] where & denotes the vector outer
// product and T(f) denotes the truth table of f.
// The resulting matrix is l x m, and is in-place added into `out`.
template <Mode mode, typename F>
void unary_outer_product(
    const F& f, const ShareSpan<mode>& x, const ShareSpan<mode>& y, ShareMatrix<mode>& out) {
  assert(x.rows() <= maximum_unary_outer_product_size);
  assert(x.cols() == 1);
  assert(y.cols() == 1);

  const auto n = x.rows();
  const auto m = y.rows();
  assert(out.cols() == m);

  std::vector<Share<mode>> seeds(1 << n);

  // We maintain the seed buffer by putting seeds into appropriate tree locations.
  // The buffer only has to be large enough for the final layer as we only
  // store intermediate seeds temporarily.

  // E keeps track of the missing tree node.
  std::size_t missing = 0;

  // As a base case, we can derive the first two seeds from the possible labels for x[0].
  if constexpr (mode == Mode::G) {
    if (x[n-1].color()) {
      seeds[0] = x[n-1].H(); // S00
      seeds[1] = (~x[n-1]).H(); // S01
    } else {
      seeds[1] = x[n-1].H(); // S00
      seeds[0] = (~x[n-1]).H(); // S01
    }
  } else {
    seeds[!x[n-1].color()] = x[n-1].H();
  }
  missing |= x[n-1].color();
  ++Share<mode>::nonce;

  const auto one = Share<mode>::bit(true);
  const auto zero = Share<mode>::bit(false);

  // Now, iterate over the levels of the tree.
  for (std::size_t i = 1; i < n; ++i) {

    const auto key0 = x[n-i-1] ^ (x[n-i-1].color() ? zero : one);
    const auto key1 = key0 ^ one;

    if constexpr (mode == Mode::G) {
      // Maintain xor sums of all odd seeds/all even seeds
      Share<mode> odds = std::bitset<128> { 0 };
      Share<mode> evens = std::bitset<128> { 0 };
      // Work backwards across the level so as to not overwrite the parent seed
      // until it is no longer needed.
      for (int j = (1 << (i-1)); j >= 0; --j) {
        seeds[j*2 + 1] = seeds[j].H(0);
        seeds[j*2] = seeds[j].H(1);
        evens ^= seeds[j*2];
        odds ^= seeds[j*2 + 1];
      }

      (evens ^ (key0.H())).send();
      (odds ^ (key1.H())).send();

      const auto bit = x[n-i-1].color();
      missing = (missing << 1) | bit;
    } else {
      const auto g_evens = Share<mode>::recv();
      const auto g_odds = Share<mode>::recv();

      Share<mode> e_evens = std::bitset<128> { 0 };
      Share<mode> e_odds = std::bitset<128> { 0 };

      for (int j = (1 << (i-1)); j >= 0; --j) {
        if (j != missing) {
          seeds[j*2 + 1] = seeds[j].H(0);
          seeds[j*2] = seeds[j].H(1);
          e_evens ^= seeds[j*2];
          e_odds ^= seeds[j*2 + 1];
        }
      }

      // use the color of the `i`th share to figure out which element is missing at the next level.
      const auto bit = x[n-i-1].color();
      missing = (missing << 1) | bit;

      // assign the sibling of the missing node by (1) decrypting the appropriate row given by
      seeds[missing ^ 1] = x[n-i-1].H() ^ (bit ? (g_evens ^ e_evens) : (g_odds ^ e_odds));
    }

    ++Share<mode>::nonce;
  }


  // Now we are ready to compute the outer product.
  // For each share (B, B + bDelta)
  // G sends the sum (XOR_i A_i) + B, which allows E to obtain A_{x + gamma} + bDelta
  if constexpr (mode == Mode::G) {
    for (std::size_t j = 0; j < m; ++j) {
      Share<mode> sum = std::bitset<128> { 0 };
      for (std::size_t i = 0; i < (1 << n); ++i) {
        const auto s = seeds[i].H();
        sum ^= s;
        ++Share<mode>::nonce;

        f(i, j, s, out);
      }
      sum ^= y[j];
      sum.send();
    }
  } else {
    for (std::size_t j = 0; j < m; ++j) {
      const Share<mode> g_sum = Share<mode>::recv();
      Share<mode> e_sum = std::bitset<128> { 0 };
      for (std::size_t i = 0; i < (1 << n); ++i) {
        if (i != missing) {
          const auto s = seeds[i].H();
          e_sum ^= s;
          f(i, j, s, out);
        }
        ++Share<mode>::nonce;
      }
      const auto s = e_sum ^ g_sum ^ y[j];
      f(missing, j, s, out);
    }
  }
}


#endif
