#include "ownpf.h"

// Witness Generation
G1 memWitGen(
    PCS &pcs,
    FrVec &S,
    Fr n
) {
    // Remove
    // S.erase(remove(S.begin(), S.end(), n), S.end());
    FrVec SPoly = constructMemPoly(pcs.pp.nttctx, S);
    FrVec nPoly = {n, 1};
    FrvT_2 divret = PolyDiv(SPoly, nPoly);
    FrVec retPoly = get<0>(divret);
    return pcs.commit_g1(retPoly);
}

tuple<G1, G1> nonMemWitGen(
    PCS &pcs,
    FrVec &S,
    Fr n    
) {
    // Make Polynomials
    FrVec SPoly = constructMemPoly(pcs.pp.nttctx, S);
    FrVec nPoly = {n, 1};

    FrvT_3 ret = xGCD(pcs.pp.nttctx, SPoly, nPoly);
    FrVec alpha = get<0>(ret);
    FrVec beta = get<1>(ret);
    FrVec gcd = get<2>(ret);

    if (gcd.size() > 1) {
        throw std::runtime_error("GCD is NOT correctly calculated");
    }
    G1 w_1 = pcs.commit_g1(beta);
    G1 w_2 = pcs.commit_g1(alpha);

    return tuple<G1, G1>(w_1, w_2);
} 

// Ownership Proof
void OwnPf::init(ZKAccSetup &setup) {
    this->setup = setup;
    this->piMP.init(setup);
    this->piPoK.init(setup);
}

void OwnPf::prove(
    G2 C_id, G2 C_A, Fr n, 
    FrVec &A, Fr id, Fr r, Fr delta
) {
    G2 C_pok = C_id - setup.g2si[1];

    // Witness Generation
    G1 w_n = memWitGen(*setup.pcs, A, n);

    // Compute C_I
    FrVec nPoly = {n, 1};
    G2 C_I = setup.pcs->commit_g2(nPoly);

    // Do Proof
    piPoK.prove(C_pok, id, r);
    piMP.prove(C_I, C_A, w_n, Fr(0), delta);
}



// Verify
bool OwnPf_verify(OwnPf *pi) {
    bool flag1 = PoK2_verify(&pi->piPoK);
    bool flag2 = ZKMP_verify(&pi->piMP);
    return flag1 && flag2;
}


// Non-Ownership Proof
void NonOwnPf::init(ZKAccSetup &setup) {
    this->setup = setup;
    this->piNMP.init(setup);
    this->piPoK.init(setup);
}

void NonOwnPf::prove(
    G2 C_id, G2 C_A, Fr n, 
    FrVec &A, Fr id, Fr r, Fr delta
) {
    G2 C_pok = C_id - setup.g2si[1];

    // Witness Generation
    tuple<G1,G1> wits = nonMemWitGen(*setup.pcs, A, n);
    G1 w_1 = get<0>(wits);
    G1 w_2 = get<1>(wits);

    // Compute C_I
    FrVec nPoly = {n, 1};
    G2 C_I = setup.pcs->commit_g2(nPoly);


    piPoK.prove(C_pok, id, r);
    piNMP.prove(C_I, C_A, w_1, w_2, Fr(0), delta);
}

// Verify
bool NonOwnPf_verify(NonOwnPf *pi) {
    bool flag1 = PoK2_verify(&pi->piPoK);
    bool flag2 = ZKNMP_verify(&pi->piNMP);

    return flag1 && flag2;
}

