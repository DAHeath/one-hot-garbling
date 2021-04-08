#ifndef PRIVATE_FUNCTION_H__
#define PRIVATE_FUNCTION_H__

#include "share.h"
#include "truth_table.h"

#include <span>


/**
 * G inputs a private function `f` described by a truth table and the parties
 * input a sharing [|x|].
 * The parties output [|f(x)|].
 */
template <Mode mode>
void private_function(
    const TruthTable&, // G inputs f
    std::span<const Share<mode>>, // [|x|]
    std::span<Share<mode>>); // [|f(x)|]


#endif
