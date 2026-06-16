#ifndef _UTILS_H
#define _UTILS_H

#include "poly.h"


// Useful Subroutines
// Inner Product
Fr innerProd(FrVec &v1, FrVec &v2);

// Hadamard (i.e., Elementwise) Product
FrVec hadProduct(FrVec &v1, FrVec &v2);

// Sum all components
Fr sumAllComps(FrVec &v);


// Membership Polynomial
// Complexity: O(nlogn^2; stable algorithm)
FrVec constructMemPoly(
    NttCtx &nttctx,
    FrVec &S
);

// Membership Polynomial for M - N
FrVec constructMemPolyBatch (
    NttCtx &nttctx,
    FrVec &M,
    FrVec &N
);

FrVec constructMemPolyBatchSorted (
    NttCtx &nttctx,
    FrVec &M,
    FrVec &N
);

#endif