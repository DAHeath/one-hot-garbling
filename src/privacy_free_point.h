#ifndef PRIVACY_FREE_POINT_H__
#define PRIVACY_FREE_POINT_H__


#include "share.h"
#include <span>


template <Mode mode>
std::size_t privacy_free_point(std::span<const Share<mode>>, std::span<Share<mode>>);


#endif
