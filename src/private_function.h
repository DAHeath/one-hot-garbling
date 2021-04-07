#ifndef PRIVATE_FUNCTION_H__
#define PRIVATE_FUNCTION_H__

#include "share.h"

#include <span>


template <Mode mode>
void private_function(
    std::span<const std::bitset<128>>, // G inputs f
    std::span<const Share<mode>>, // [|x|]
    std::span<Share<mode>>); // [|f(x)|]


#endif
