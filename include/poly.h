#ifndef _POLY_H
#define _POLY_H

#include <mcl/bls12_381.hpp>
#include <mcl/ntt.hpp>

#include <cstddef>
#include <tuple>
#include <vector>
#include <unordered_map>
#include <memory>
#include <stdexcept>

using namespace mcl;
using namespace mcl::bn;
using namespace std;

std::size_t nextPow2(std::size_t n);

// Coefficients are stored in little-endian order:
//   a[0] + a[1] X + ... + a[d] X^d.
// The zero polynomial is normalized as {0}; public functions also accept {}.
typedef std::vector<Fr> FrVec;
typedef std::tuple<FrVec, FrVec> FrvT_2;          // (quotient, remainder)
typedef std::tuple<FrVec, FrVec, FrVec> FrvT_3;   // (s, t, gcd), s*A + t*B = gcd

// NTT Context
struct NttCtx {
    std::unordered_map<std::size_t, std::unique_ptr<Ntt<Fr>>> cache;

    void init(std::size_t N) {
        for (std::size_t n = 2; n <= nextPow2(N); n <<= 1) build(n);
    }

    Ntt<Fr>& build(std::size_t n) {
        auto& slot = cache[n];
        if (!slot) {
            slot = std::make_unique<Ntt<Fr>>();
            if (!slot->init(n)) throw std::runtime_error("NTT init failed");
        }
        return *slot;
    }        

    Ntt<Fr>& getNtt(std::size_t n) const {
        auto it = cache.find(n);
        if (it == cache.end())
            throw std::runtime_error("NTT size not pre-warmed; call init() first");
        return *it->second;
    }
};


// Basic polynomial utilities.
bool IsPolyZero(const FrVec& a);
FrVec PolyCondense(const FrVec& a);
bool IsPolyEqual(const FrVec& a, const FrVec& b);

// Arithmetic.
FrVec PolyAdd(const FrVec& a, const FrVec& b);
FrVec PolySub(const FrVec& a, const FrVec& b);
FrVec PolyMul(NttCtx &nttctx, const FrVec& a, const FrVec& b);

// Division: A = quotient * B + remainder, deg(remainder) < deg(B).
FrVec PolyLongDiv(const FrVec& A, const FrVec& B);
FrvT_2 PolyDiv(const FrVec& A, const FrVec& B);
FrvT_2 PolyDiv_na(NttCtx &nttctx, const FrVec& A, const FrVec& B);

// Evaluation and calculus.
Fr PolyEvaluate(const FrVec& a, const Fr& x);
FrVec PolyDifferentiate(const FrVec& A);

// Extended Euclidean algorithm.
FrvT_3 xGCD(NttCtx &nttctx, const FrVec& a, const FrVec& b);

#endif
