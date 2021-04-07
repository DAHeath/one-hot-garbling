#ifndef PRIVACY_FREE_POINT_H__
#define PRIVACY_FREE_POINT_H__


#include "share.h"
#include <span>


// Let delta be the color of the input labels
template <Mode mode>
std::size_t // x + delta
privacy_free_point(
    std::span<const Share<mode>>, //  [|x|]
    std::span<Share<mode>>); // [|U(x + delta)|]


#endif
