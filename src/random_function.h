#ifndef HALF_RANDOM_FUNCTION_H__
#define HALF_RANDOM_FUNCTION_H__


#include "share.h"
#include <span>


template <Mode mode>
void half_random_function(
    std::size_t i,
    const Share<mode>& key,
    std::span<const Share<mode>> unary,
    std::span<std::bitset<128>> table,
    std::span<Share<mode>> out);


template <Mode mode>
void random_function(
    std::size_t, // E inputs x
    std::span<const Share<mode>>, // [|x|]
    std::span<const Share<mode>>, // [|U(x)|]
    std::span<std::bitset<128>> table,
    std::span<Share<mode>> out);


#endif
