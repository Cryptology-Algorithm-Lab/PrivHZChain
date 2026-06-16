#include "poly.h"

#include <algorithm>
#include <cstdint>
#include <stdexcept>

constexpr std::size_t kNaiveMulThreshold = 64;

Fr zeroFr() {
    Fr z;
    z.clear();
    return z;
}

Fr oneFr() {
    Fr o;
    o = 1;
    return o;
}

std::size_t nextPow2(std::size_t n) {
    if (n <= 1) {
        return 1;
    }
    --n;
    for (std::size_t shift = 1; shift < sizeof(std::size_t) * 8; shift <<= 1) {
        n |= (n >> shift);
    }
    return n + 1;
}

std::size_t degree(const FrVec& a) {
    FrVec c = PolyCondense(a);
    if (IsPolyZero(c)) {
        return 0;
    }
    return c.size() - 1;
}

void scaleInPlace(FrVec& a, const Fr& c) {
    for (Fr& x : a) {
        x *= c;
    }
    a = PolyCondense(a);
}

FrVec PolyMulNaive(const FrVec& a, const FrVec& b) {
    FrVec lhs = PolyCondense(a);
    FrVec rhs = PolyCondense(b);

    if (IsPolyZero(lhs) || IsPolyZero(rhs)) {
        return {zeroFr()};
    }

    FrVec out(lhs.size() + rhs.size() - 1, zeroFr());
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        for (std::size_t j = 0; j < rhs.size(); ++j) {
            out[i + j] += lhs[i] * rhs[j];
        }
    }
    return PolyCondense(out);
}

bool IsPolyZero(const FrVec& a) {
    for (const Fr& coeff : a) {
        if (!coeff.isZero()) {
            return false;
        }
    }
    return true;
}

FrVec PolyCondense(const FrVec& a) {
    if (a.empty()) {
        return {zeroFr()};
    }

    std::size_t n = a.size();
    while (n > 1 && a[n - 1].isZero()) {
        --n;
    }

    if (n == 0) {
        return {zeroFr()};
    }
    return FrVec(a.begin(), a.begin() + n);
}

bool IsPolyEqual(const FrVec& a, const FrVec& b) {
    FrVec lhs = PolyCondense(a);
    FrVec rhs = PolyCondense(b);

    if (lhs.size() != rhs.size()) {
        return false;
    }
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        if (lhs[i] != rhs[i]) {
            return false;
        }
    }
    return true;
}

FrVec PolyAdd(const FrVec& a, const FrVec& b) {
    FrVec lhs = PolyCondense(a);
    FrVec rhs = PolyCondense(b);

    const std::size_t n = std::max(lhs.size(), rhs.size());
    FrVec out(n, zeroFr());

    for (std::size_t i = 0; i < n; ++i) {
        if (i < lhs.size()) {
            out[i] += lhs[i];
        }
        if (i < rhs.size()) {
            out[i] += rhs[i];
        }
    }
    return PolyCondense(out);
}

FrVec PolySub(const FrVec& a, const FrVec& b) {
    FrVec lhs = PolyCondense(a);
    FrVec rhs = PolyCondense(b);

    const std::size_t n = std::max(lhs.size(), rhs.size());
    FrVec out(n, zeroFr());

    for (std::size_t i = 0; i < n; ++i) {
        if (i < lhs.size()) {
            out[i] += lhs[i];
        }
        if (i < rhs.size()) {
            out[i] -= rhs[i];
        }
    }
    return PolyCondense(out);
}

FrVec PolyMul(NttCtx &nttctx, const FrVec& a, const FrVec& b) {
    FrVec lhs = PolyCondense(a);
    FrVec rhs = PolyCondense(b);

    if (IsPolyZero(lhs) || IsPolyZero(rhs)) {
        return {zeroFr()};
    }

    const std::size_t outSize = lhs.size() + rhs.size() - 1;
    if (std::min(lhs.size(), rhs.size()) <= kNaiveMulThreshold) {
        return PolyMulNaive(lhs, rhs);
    }

    const std::size_t n = nextPow2(outSize);
    lhs.resize(n, zeroFr());
    rhs.resize(n, zeroFr());
    FrVec out(n, zeroFr());

    Ntt<Fr>& ntt = nttctx.getNtt(n);
    ntt.ntt(lhs.data());
    ntt.ntt(rhs.data());

    for (std::size_t i = 0; i < n; ++i) {
        out[i] = lhs[i] * rhs[i];
    }

    ntt.intt(out.data());
    out.resize(outSize);
    return PolyCondense(out);
}

FrVec PolyLongDiv(const FrVec& A, const FrVec& B) {
    return get<0>(PolyDiv(A, B));
}

FrvT_2 PolyDiv(const FrVec& A, const FrVec& B) {
    FrVec divisor = PolyCondense(B);
    if (IsPolyZero(divisor)) {
        throw std::runtime_error("PolyDiv: division by zero polynomial");
    }

    FrVec rem = PolyCondense(A);
    if (IsPolyZero(rem) || degree(rem) < degree(divisor)) {
        return FrvT_2({zeroFr()}, rem);
    }

    const std::size_t degA = rem.size() - 1;
    const std::size_t degB = divisor.size() - 1;
    FrVec quot(degA - degB + 1, zeroFr());

    const Fr lead = divisor.back();
    for (std::size_t k = degA - degB + 1; k > 0; --k) {
        const std::size_t offset = k - 1;
        Fr q = rem[degB + offset] / lead;
        quot[offset] = q;

        if (q.isZero()) {
            continue;
        }
        for (std::size_t j = 0; j <= degB; ++j) {
            rem[j + offset] -= q * divisor[j];
        }
    }

    return FrvT_2(PolyCondense(quot), PolyCondense(rem));
}

FrvT_2 PolyDiv_na(NttCtx &nttctx, const FrVec& A, const FrVec& B) {
    FrVec quotient = PolyLongDiv(A, B);
    FrVec remainder = PolySub(A, PolyMul(nttctx, quotient, B));
    return FrvT_2(PolyCondense(quotient), PolyCondense(remainder));
}

Fr PolyEvaluate(const FrVec& a, const Fr& x) {
    FrVec p = PolyCondense(a);
    Fr acc = zeroFr();

    for (std::size_t i = p.size(); i > 0; --i) {
        acc *= x;
        acc += p[i - 1];
    }
    return acc;
}

FrVec PolyDifferentiate(const FrVec& A) {
    FrVec p = PolyCondense(A);
    if (p.size() <= 1) {
        return {zeroFr()};
    }

    FrVec out(p.size() - 1, zeroFr());
    for (std::size_t i = 1; i < p.size(); ++i) {
        Fr scalar;
        scalar = static_cast<uint64_t>(i);
        out[i - 1] = p[i] * scalar;
    }
    return PolyCondense(out);
}

FrvT_3 xGCD(NttCtx &nttctx, const FrVec& a, const FrVec& b) {
    FrVec r0 = PolyCondense(a);
    FrVec r1 = PolyCondense(b);

    if (IsPolyZero(r0) && IsPolyZero(r1)) {
        throw std::runtime_error("xGCD: gcd(0, 0) is undefined");
    }

    FrVec s0 = {oneFr()};
    FrVec s1 = {zeroFr()};
    FrVec t0 = {zeroFr()};
    FrVec t1 = {oneFr()};

    while (!IsPolyZero(r1)) {
        FrvT_2 qr = PolyDiv(r0, r1);
        FrVec q = get<0>(qr);
        FrVec r2 = get<1>(qr);

        Ntt<Fr>& ntt = nttctx.getNtt(nextPow2(r0.size()));
        FrVec s2 = PolySub(s0, PolyMul(nttctx, q, s1));
        FrVec t2 = PolySub(t0, PolyMul(nttctx, q, t1));

        r0.swap(r1);
        r1.swap(r2);
        s0.swap(s1);
        s1.swap(s2);
        t0.swap(t1);
        t1.swap(t2);
    }

    Fr lead = r0.back();
    Fr invLead = oneFr() / lead;
    scaleInPlace(r0, invLead);
    scaleInPlace(s0, invLead);
    scaleInPlace(t0, invLead);

    return FrvT_3(PolyCondense(s0), PolyCondense(t0), PolyCondense(r0));
}