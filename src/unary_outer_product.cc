#include "unary_outer_product.h"


thread_local std::vector<Share<Mode::G>> g_seeds(1 << maximum_unary_outer_product_size);
thread_local std::vector<Share<Mode::E>> e_seeds(1 << maximum_unary_outer_product_size);

std::span<Share<Mode::G>> g_unary_seeds() { return g_seeds; }
std::span<Share<Mode::E>> e_unary_seeds() { return e_seeds; }
