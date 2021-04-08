#ifndef HALF_UNIFORM_FUNCTION_H__
#define HALF_UNIFORM_FUNCTION_H__


#include "share.h"
#include "truth_table.h"
#include <span>


/**
 * Apply a uniformly chosen function `r` to the input using unary computation.
 * For G, the output truth table is XOR stacked into the argument `r`.
 * The output is XOR stacked into `out`.
 *
 * E learns "half" of `r` corresponding to truth table input column `i`.
 */
template <Mode mode>
void half_uniform_function(
    std::size_t i,
    std::span<const Share<mode>> sx, // [| x |]
    std::span<const Share<mode>>, // [| U(x + delta) |]
    TruthTable&,
    std::span<Share<mode>> out);


/**
 * Apply a uniformly chosen function `r` to the input using unary computation.
 * For G, the output truth table is XOR stacked into the argument `r`.
 * The output is XOR stacked into `out`.
 *
 * This is achieved by `n` applications of `half_uniform_function` and by
 * masking each column with a uniform bit.
 */
template <Mode mode>
void uniform_function(
    std::span<const Share<mode>>, // [|x|]
    std::span<const Share<mode>>, // [|U(x)|]
    TruthTable&,
    std::span<Share<mode>> out);


#endif
