#include "core.h"

// Hidden ACC
HACCParams HACC_Setup (
    uint32_t maxDegree
) {
    // Select Trapdoor
    Fr s;
    s.setByCSPRNG();

    // Select base points
    G1 g1, gfrak, hfrak; G2 g2, h2;
    hashAndMapToG1(g1,"g1");
    hashAndMapToG1(gfrak,"gfrak");
    hashAndMapToG1(hfrak,"hfrak");
    hashAndMapToG2(g2,"g2");
    hashAndMapToG2(h2,"h2");    

    // Construct the base Vector
    std::vector<G1> g1base(maxDegree, g1);
    std::vector<G2> g2base(maxDegree, g2);

    for (uint32_t i = 1; i < maxDegree; i++) {
        G1::mul(g1base[i], g1base[i-1], s);
        G2::mul(g2base[i], g2base[i-1], s);
    }
    return HACCParams {
        g1, g2, h2, gfrak, hfrak,
        g1base, g2base
    };
}

// Commitment from Coefficients
G1 HACC_commitCoeffsG1(
    HACCParams pp, 
    std::vector<Fr> coeffs
) {
    uint32_t numCoeffs = coeffs.size();
    G1 ret;
    G1::mulVec(ret, &pp.g1base[0], &coeffs[0], numCoeffs);
    return ret;
}

G2 HACC_commitCoeffsG2(
    HACCParams pp, 
    std::vector<Fr> coeffs
) {
    uint32_t numCoeffs = coeffs.size();
    G2 ret;
    G2::mulVec(ret, &pp.g2base[0], &coeffs[0], numCoeffs);
    return ret;
}

G1 HACC_MemWitGen(
    HACCParams pp, 
    std::vector<Fr> A,
    std::vector<Fr> I
) {
    std::vector<Fr> div = PolyLongDiv(A, I);
    return HACC_commitCoeffsG1(pp, div);
}

tuple<G1, G1> HACC_NonMemWitGen(
    HACCParams pp, 
    std::vector<Fr> A,
    std::vector<Fr> I
) {
    FrvT_3 ret = xGCD(A, I);
    FrVec alpha = get<0>(ret);
    FrVec beta =  get<1>(ret);
    FrVec gcd =  get<2>(ret);

    if (!IsPolyEqual(gcd, FrVec({1}))) {
        throw std::runtime_error("GCD is NOT 1...");
    }  

    G1 w_A = HACC_commitCoeffsG1(pp, alpha);
    G1 w_I = HACC_commitCoeffsG1(pp, beta);
    return tuple<G1, G1>(w_A, w_I);
}

// ZKMP
ZKMPProof ZKMP_prove(
    HACCParams pp,
    G2 C_I,
    G2 C_A,
    G1 W_I,
    Fr r_I,
    Fr r_A
) {
    // Pick Randomness
    Fr t_1, t_2, r_rI, r_rA, r_t1, r_t2, r_d1, r_d2;
    t_1.setByCSPRNG(); t_2.setByCSPRNG();
    r_rI.setByCSPRNG(); r_rA.setByCSPRNG();
    r_t1.setByCSPRNG(); r_t2.setByCSPRNG();
    r_d1.setByCSPRNG(); r_d2.setByCSPRNG();

    // Compute Delta
    Fr d_1 = r_I * t_1;
    Fr d_2 = r_I * t_2;

    // Comms
    G1 P1 = pp.g1 * t_1 + pp.gfrak * t_2;
    G1 R1 = pp.g1 * r_t1 + pp.gfrak * r_t2;
    G1 P2 = W_I + pp.gfrak * t_1;
    G1 R2 = P1 * r_rI + pp.g1 * (-r_d1) + pp.gfrak * (-r_d2);

    // Compute R3
    G2 _op3, _op4;
    G2 first = C_I * r_t1 + pp.h2 * (-r_d1);
    G1 second = P2 * r_rI + pp.g1 * (-r_rA);

    G1 leftMP[2] = {pp.gfrak, second};
    G2 rightMP[2] = {first, pp.h2};
    GT R3;
    millerLoopVec(R3, leftMP, rightMP, 2);
    finalExp(R3, R3);

    // Get Challenge
    string buf = P1.getStr() + P2.getStr()+ R1.getStr() + R2.getStr() + R3.getStr();
    Fr c;
    c.setHashOf(buf);

    // Compute S Values
    Fr s_rI = r_rI + c * r_I;
    Fr s_rA = r_rA + c * r_A;
    Fr s_t1 = r_t1 + c * t_1;
    Fr s_t2 = r_t2 + c * t_2;
    Fr s_d1 = r_d1 + c * d_1;
    Fr s_d2 = r_d2 + c * d_2;

    // Done
    return ZKMPProof {
        P1, P2, R1, R2, R3,
        s_rI, s_rA, s_t1, s_t2, s_d1, s_d2
    };

}

bool ZKMP_verify(
    HACCParams pp,
    ZKMPProof proof,
    G2 C_I,
    G2 C_A
) {
    // Extract Challenge from proof
    string buf = proof.P1.getStr() + proof.P2.getStr() + 
                 proof.R1.getStr() + proof.R2.getStr() + proof.R3.getStr();
    Fr c;
    c.setHashOf(buf);

    // Useful Vector
    G1 base[3] = {proof.P1, pp.g1, pp.gfrak};
    // Check R1
    G1 right1;
    Fr expon1[3] = {-c, proof.s_t1, proof.s_t2};
    G1::mulVec(right1, base, expon1, 3);
    bool check1 = (proof.R1 == right1);

    // Check R2
    G1 right2;
    Fr expon2[3] = {proof.s_rI, -proof.s_d1, -proof.s_d2};
    G1::mulVec(right2, base, expon2, 3);
    bool check2 = (proof.R2 == right2);

    // Check R3
    G1 leftMP[3] = {pp.gfrak, proof.P2, pp.g1};

    // Prepare for the left term
    G2 firstTerm = C_I * proof.s_t1 + pp.h2 * (-proof.s_d1); 
    G2 secondTerm = pp.h2 * proof.s_rI + C_I * (-c);
    G2 thirdTerm = pp.h2 * (-proof.s_rA) + C_A * c;

    G2 rightMP[3] = {firstTerm, secondTerm, thirdTerm};

    // Do miller loop
    GT right3;
    millerLoopVec(right3,leftMP,rightMP,3);
    finalExp(right3,right3);
    bool check3 = (proof.R3 == right3);

    return check1 && check2 && check3;
}

int main() {


    initPairing(mcl::BLS12_381);

    // Setup
    HACCParams pp = HACC_Setup(8);
    FrVec A = {1, 3, 3, 1};
    FrVec I = {1,1};

    Fr r_A, r_I;
    r_A.setByCSPRNG(); r_I.setByCSPRNG();

    // Commit 
    G2 C_I = HACC_commitCoeffsG2(pp, I) + pp.h2 * r_I;
    G2 C_A = HACC_commitCoeffsG2(pp, A) + pp.h2 * r_A;    
    G1 W_I = HACC_MemWitGen(pp, A, I);

    double PTtot = 0;
    double VTtot = 0;

    for (int i = 0; i < 1000; i++) {
        // Proof
        auto t1p = chrono::steady_clock::now();
        ZKMPProof proof = ZKMP_prove(pp, C_I, C_A, W_I, r_I, r_A);
        auto t2p = chrono::steady_clock::now();
        double PT = chrono::duration_cast<chrono::microseconds>(t2p-t1p).count();

        // Verify 
        auto t1v = chrono::high_resolution_clock::now();
        bool ret = ZKMP_verify(pp, proof, C_I, C_A);
        auto t2v = chrono::high_resolution_clock::now();
        double VT = chrono::duration_cast<chrono::microseconds>(t2v-t1v).count();

        if (ret == 0) {
            throw std::runtime_error("Verification Failed...");
        }

        PTtot += PT;
        VTtot += VT;
    }

    cout << "Prover Time (ms): " << PTtot / 1000 / 1000 << endl;
    cout << "Verifier Time (ms): " << VTtot / 1000 / 1000 << endl;
}
