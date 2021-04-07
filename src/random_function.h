#ifndef HALF_RANDOM_FUNCTION_H__
#define HALF_RANDOM_FUNCTION_H__


#include "share.h"
#include "truth_table.h"
#include <span>


template <Mode mode>
void half_random_function(
    std::size_t i,
    std::span<const Share<mode>> sx, // [| x |]
    std::span<const Share<mode>>, // [| U(x + delta) |]
    TruthTable&,
    std::span<Share<mode>> out);


template <Mode mode>
void random_function(
    std::span<const Share<mode>>, // [|x|]
    std::span<const Share<mode>>, // [|U(x)|]
    TruthTable&,
    std::span<Share<mode>> out);


#endif
