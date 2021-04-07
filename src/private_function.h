#ifndef PRIVATE_FUNCTION_H__
#define PRIVATE_FUNCTION_H__

#include "share.h"
#include "truth_table.h"

#include <span>


template <Mode mode>
void private_function(
    const TruthTable&, // G inputs f
    std::span<const Share<mode>>, // [|x|]
    std::span<Share<mode>>); // [|f(x)|]


#endif
