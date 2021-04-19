#ifndef UNARY_OUTER_PRODUCT_H__
#define UNARY_OUTER_PRODUCT_H__


#include "share_matrix.h"
#include "table.h"

// Let c be the color of x.
// Let n be the length of x.
// Let m be the length of y.
// Let f be a function from n bits to l bits
//
// This computes [|T(f) * (U(x + c) & y) |] where & denotes the vector outer
// product and T(f) denotes the truth table of f.
// The resulting matrix is l x m, and is in-place added into `out`.
template <Mode mode>
void unary_outer_product(
    const Table&,
    const ShareCSpan<mode>&,
    const ShareCSpan<mode>&,
    const MatrixView<Share<mode>&>&);

void initialize_gjobs();
void initialize_ejobs();
void finalize_jobs();

#endif
