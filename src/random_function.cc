#include "random_function.h"


template <Mode mode>
void half_random_function(
    std::size_t i,
    const Share<mode>& key,
    std::span<const Share<mode>> unary,
    std::span<std::bitset<128>> table,
    std::span<Share<mode>> out) {
  // TODO
}


template <Mode mode>
void random_function(
    std::size_t, // E inputs x
    std::span<const Share<mode>>, // [|x|]
    std::span<const Share<mode>>, // [|U(x)|]
    std::span<std::bitset<128>> table,
    std::span<Share<mode>> out) {
  // TODO
}

template void half_random_function<Mode::G>(
    std::size_t,
    const Share<Mode::G>&,
    std::span<const Share<Mode::G>>,
    std::span<std::bitset<128>>,
    std::span<Share<Mode::G>>);
template void half_random_function<Mode::E>(
    std::size_t,
    const Share<Mode::E>&,
    std::span<const Share<Mode::E>>,
    std::span<std::bitset<128>>,
    std::span<Share<Mode::E>>);

template void random_function<Mode::G>(
    std::size_t,
    std::span<const Share<Mode::G>>,
    std::span<const Share<Mode::G>>,
    std::span<std::bitset<128>>,
    std::span<Share<Mode::G>>);
template void random_function<Mode::E>(
    std::size_t,
    std::span<const Share<Mode::E>>,
    std::span<const Share<Mode::E>>,
    std::span<std::bitset<128>>,
    std::span<Share<Mode::E>>);
