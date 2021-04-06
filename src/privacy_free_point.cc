#include "privacy_free_point.h"


/**
 * Parties map [|x|] to [|U(x + gamma)|] where gamma is the color of the input labels [|x|];
 * Notice, E knows x + gamma.
 */
template <Mode mode>
std::size_t privacy_free_point(std::span<const Share<mode>> point, std::span<Share<mode>> out) {
  const auto n = point.size();
  assert(out.size() == (1 << n));
  // We maintain the output buffer by putting seeds into appropriate tree locations.
  // The buffer only has to be large enough for the output, and we only store intermediate seeds.

  // E keeps track of the missing tree node.
  std::size_t missing = 0;

  // As a base case, we can derive the first two seeds from the possible labels for x[0].
  if constexpr (mode == Mode::G) {
    if (point[0].color()) {
      out[0] = point[0].H(); // S00
      out[1] = (~point[0]).H(); // S01
    } else {
      out[1] = point[0].H(); // S00
      out[0] = (~point[0]).H(); // S01
    }
  } else {
    out[!point[0].color()] = point[0].H();
    missing |= point[0].color();
  }
  ++Share<mode>::nonce;

  // Now, iterate over the levels of the tree.
  for (std::size_t i = 1; i < n; ++i) {
    if constexpr (mode == Mode::G) {
      // Maintain xor sums of all odd seeds/all even seeds
      Share<mode> odds = std::bitset<128> { 0 };
      Share<mode> evens = std::bitset<128> { 0 };
      // Work backwards across the level so as to not overwrite the parent seed
      // until it is no longer needed.
      for (int j = (1 << (i-1)); j >= 0; --j) {
        out[j*2 + 1] = out[j].H(0);
        out[j*2] = out[j].H(1);
        evens ^= out[j*2];
        odds ^= out[j*2 + 1];
      }

      if (point[i].color()) {
        (evens ^ (point[i].H())).send();
        (odds ^ ((~point[i]).H())).send();
      } else {
        (evens ^ ((~point[i]).H())).send();
        (odds ^ (point[i].H())).send();
      }
    } else {
      const auto g_evens = Share<mode>::recv();
      const auto g_odds = Share<mode>::recv();

      Share<mode> e_evens = std::bitset<128> { 0 };
      Share<mode> e_odds = std::bitset<128> { 0 };

      for (int j = (1 << (i-1)); j >= 0; --j) {
        if (j != missing) {
          out[j*2 + 1] = out[j].H(0);
          out[j*2] = out[j].H(1);
          e_evens ^= out[j*2];
          e_odds ^= out[j*2 + 1];
        }
      }

      // use the color of the `i`th share to figure out which element is missing at the next level.
      const auto bit = point[i].color();
      missing = (missing << 1) | bit;

      const auto dec = point[i].H() ^ (bit ? g_evens : g_odds);

      // assign the sibling of the missing node by (1) decrypting the appropriate row given by
      out[missing ^ 1] = point[i].H() ^ (bit ? (g_evens ^ e_evens) : (g_odds ^ e_odds));
    }

    ++Share<mode>::nonce;
  }

  // Now, each label A_i is derived from the leaf seed Sni.
  // G sends the sum (XOR_i A_i) + delta, which allows E to obtain A_{x + gamma} + delta
  if constexpr (mode == Mode::G) {
    Share<mode> sum = std::bitset<128> { 0 };
    for (std::size_t i = 0; i < (1 << n); ++i) {
      out[i] = out[i].H();
      sum ^= out[i];
      ++Share<mode>::nonce;
    }
    sum ^= Share<mode>::delta;
    sum.send();
  } else {
    const Share<mode> g_sum = Share<mode>::recv();
    Share<mode> e_sum = std::bitset<128> { 0 };
    for (std::size_t i = 0; i < (1 << n); ++i) {
      if (i != missing) {
        out[i] = out[i].H();
        e_sum ^= out[i];
      }
      ++Share<mode>::nonce;
    }
    out[missing] = e_sum ^ g_sum;
  }
  return missing;
}


template std::size_t privacy_free_point<Mode::G>(std::span<const Share<Mode::G>>, std::span<Share<Mode::G>>);
template std::size_t privacy_free_point<Mode::E>(std::span<const Share<Mode::E>>, std::span<Share<Mode::E>>);
