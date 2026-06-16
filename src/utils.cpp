#include "utils.h"

// Useful Subroutines
Fr innerProd(FrVec &v1, FrVec &v2) {
    uint32_t n = v1.size();
    Fr ret = Fr(0);
    for (uint32_t i = 0; i < n; i++) {
        ret += (v1[i] * v2[i]);
    }
    return ret;
}

// Hadamard (i.e., Elementwise) Product
FrVec hadProduct(FrVec &v1, FrVec &v2) {
    uint32_t n = v1.size();
    FrVec ret;
    for (uint32_t i = 0; i < n; i++) {
        ret.push_back(v1[i] * v2[i]);
    }
    return ret;
}

// Sum all components
Fr sumAllComps(FrVec &v) {
    uint32_t n = v.size();
    Fr ret = v[0];
    for (uint32_t i = 1; i < n; i++) {
        ret += v[i];
    }
    return ret;
}


// Membership Polynomial
// Complexity: O(nlogn^2; stable algorithm)
FrVec constructMemPoly(
    NttCtx &nttctx,
    FrVec &S
) {
    uint32_t numElts = S.size();
    if (numElts == 0) {
        return {1};
    }
    else if (numElts == 1) {
        return {S[0], 1};
    }
    FrVec left(S.begin(), S.begin() + numElts/2);
    FrVec right(S.begin() + numElts/2, S.end());

    FrVec leftRet = constructMemPoly(nttctx, left);
    FrVec rightRet = constructMemPoly(nttctx, right);
    FrVec ret = PolyMul(nttctx, leftRet, rightRet);
    return ret;
}

// Membership Polynomial for M - N
FrVec constructMemPolyBatch (
    NttCtx &nttctx,
    FrVec &M,
    FrVec &N
) {
    // Find Set Difference for sorted vectors
    std::sort(M.begin(), M.end());
    std::sort(N.begin(), N.end());
    return constructMemPolyBatchSorted(nttctx, M, N);
}

// Membership Polynomial for M - N
// For Sorted Vectors
FrVec constructMemPolyBatchSorted (
    NttCtx &nttctx,
    FrVec &M,
    FrVec &N
) {
    // Assume that the output is already sorted.
    FrVec ret;
    std::set_difference(
        M.begin(), M.end(), N.begin(), N.end(), 
        std::back_inserter(ret)
    );
    return constructMemPoly(nttctx, ret);
}