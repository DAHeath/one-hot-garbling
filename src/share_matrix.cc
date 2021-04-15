#include "share_matrix.h"

template <Mode mode>
ShareMatrix<mode> ShareMatrix<mode>::unary_outer_product(const ShareMatrix<mode>& o) const {
  assert(cols() == 1);
  assert(o.cols() == 1);

  const auto n = rows();
  const auto m = o.rows();

  std::vector<Share<mode>> seeds(1 << n);
  // We maintain the seed buffer by putting seeds into appropriate tree locations.
  // The buffer only has to be large enough for the final layer as we only
  // store intermediate seeds temporarily.

  // E keeps track of the missing tree node.
  std::size_t missing = 0;

  // As a base case, we can derive the first two seeds from the possible labels for x[0].
  if constexpr (mode == Mode::G) {
    if ((*this)[n-1].color()) {
      seeds[0] = (*this)[n-1].H(); // S00
      seeds[1] = (~(*this)[n-1]).H(); // S01
    } else {
      seeds[1] = (*this)[n-1].H(); // S00
      seeds[0] = (~(*this)[n-1]).H(); // S01
    }
  } else {
    seeds[!(*this)[n-1].color()] = (*this)[n-1].H();
  }
  missing |= (*this)[n-1].color();
  ++Share<mode>::nonce;

  const auto one = Share<mode>::bit(true);
  const auto zero = Share<mode>::bit(false);

  // Now, iterate over the levels of the tree.
  for (std::size_t i = 1; i < n; ++i) {

    const auto key0 = (*this)[n-i-1] ^ ((*this)[n-i-1].color() ? zero : one);
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

      const auto bit = (*this)[n-i-1].color();
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
      const auto bit = (*this)[n-i-1].color();
      missing = (missing << 1) | bit;

      // assign the sibling of the missing node by (1) decrypting the appropriate row given by
      seeds[missing ^ 1] = (*this)[n-i-1].H() ^ (bit ? (g_evens ^ e_evens) : (g_odds ^ e_odds));
    }

    ++Share<mode>::nonce;
  }


  ShareMatrix<mode> out(1 << n, m);

  // Now we are ready to compute the outer product.
  // For each share (B, B + bDelta)
  // G sends the sum (XOR_i A_i) + B, which allows E to obtain A_{x + gamma} + bDelta
  if constexpr (mode == Mode::G) {
    for (std::size_t j = 0; j < m; ++j) {
      Share<mode> sum = std::bitset<128> { 0 };
      for (std::size_t i = 0; i < (1 << n); ++i) {
        out(i, j) = seeds[i].H();
        sum ^= out(i, j);
        ++Share<mode>::nonce;
      }
      sum ^= o[j];
      sum.send();
    }
  } else {
    for (std::size_t j = 0; j < m; ++j) {
      const Share<mode> g_sum = Share<mode>::recv();
      Share<mode> e_sum = std::bitset<128> { 0 };
      for (std::size_t i = 0; i < (1 << n); ++i) {
        if (i != missing) {
          out(i, j) = seeds[i].H();
          e_sum ^= out(i, j);
        }
        ++Share<mode>::nonce;
      }
      out(missing, j) = e_sum ^ g_sum ^ o[j];
    }
  }
  return out;
}


template ShareMatrix<Mode::G>
ShareMatrix<Mode::G>::unary_outer_product(const ShareMatrix<Mode::G>&) const;
template ShareMatrix<Mode::E>
ShareMatrix<Mode::E>::unary_outer_product(const ShareMatrix<Mode::E>&) const;


template <Mode mode>
ShareMatrix<mode> half_outer_product(const ShareMatrix<mode>& x, const ShareMatrix<mode>& y) {
  assert(x.cols() == 1);
  assert(y.cols() == 1);

  auto id_n = identity_table(x.rows());
  return id_n * x.unary_outer_product(y);
}

template ShareMatrix<Mode::G> half_outer_product(const ShareMatrix<Mode::G>&, const ShareMatrix<Mode::G>&);
template ShareMatrix<Mode::E> half_outer_product(const ShareMatrix<Mode::E>&, const ShareMatrix<Mode::E>&);


constexpr std::size_t default_outer_product_slice_size = 8;


template <Mode mode>
ShareMatrix<mode> wide_half_outer_product(const ShareMatrix<mode>& X, const ShareMatrix<mode>& Y) {
  assert(X.cols() == 1);
  assert(Y.cols() == 1);

  const auto n = X.rows();
  const auto m = Y.rows();

  ShareMatrix<mode> out(n, m);

  const auto def = default_outer_product_slice_size;

  for (std::size_t s = 0; s < (n + def-1)/def; ++s) {
    std::size_t slice_size = def;
    if (slice_size*(s + 1) > n) { slice_size = n % slice_size; }

    ShareMatrix<mode> slice(slice_size, 1);
    for (std::size_t i = 0; i < slice_size; ++i) {
      slice[i] = X[i + s*def];
    }
    const auto half = half_outer_product<mode>(slice, Y);
    for (std::size_t i = 0; i < slice_size; ++i) {
      for (std::size_t j = 0; j < m; ++j) {
        out(def*s + i, j) = half(i, j);
      }
    }
  }
  return out;
}


template <Mode mode>
ShareMatrix<mode> outer_product(const ShareMatrix<mode>& X, const ShareMatrix<mode>& Y) {
  // An outer product can be computed by half outer products.
  // Unfortunately, due to exponential computation scaling, we need to "chunk" the outer product into slices.
  assert(X.cols() == 1);
  assert(Y.cols() == 1);

  const auto x = ShareMatrix<mode>::constant(X.color());
  const auto xy = ShareMatrix<mode>::constant(X.color().outer_product(Y.color()));

  return
    wide_half_outer_product<mode>(X, Y) ^
    wide_half_outer_product<mode>(Y, x).transposed() ^
    xy;
}

template ShareMatrix<Mode::G> outer_product(const ShareMatrix<Mode::G>&, const ShareMatrix<Mode::G>&);
template ShareMatrix<Mode::E> outer_product(const ShareMatrix<Mode::E>&, const ShareMatrix<Mode::E>&);


template <Mode mode>
ShareMatrix<mode> ShareMatrix<mode>::operator*(const ShareMatrix<mode>& o) const {
  const auto l = rows();
  const auto n = o.rows();
  const auto m = o.cols();
  assert (cols() == n);

  ShareMatrix<mode> out(l, m);
  for (std::size_t i = 0; i < n; ++i) {
    out ^= outer_product(column(i), o.row(i));
  }
  return out;
}


template ShareMatrix<Mode::G> ShareMatrix<Mode::G>::operator*(const ShareMatrix<Mode::G>&) const;
template ShareMatrix<Mode::E> ShareMatrix<Mode::E>::operator*(const ShareMatrix<Mode::E>&) const;
