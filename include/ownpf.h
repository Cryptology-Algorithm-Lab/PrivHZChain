#ifndef _OWNPF_H
#define _OWNPF_H

#include "zkacc.h"
#include "utils.h"

// Helpers
G1 memWitGen(
    PCS &pcs,
    FrVec &S,
    Fr n
);

tuple<G1, G1> nonMemWitGen(
    PCS &pcs,
    FrVec &S,
    Fr n    
);


// Ownership Proofs
struct OwnPf {
    ZKAccSetup &setup;
    ZKMP &piMP;
    PoK2 &piPoK;

    void init(ZKAccSetup &setup);

    void prove(
        G2 C_id, G2 C_A, Fr n, 
        FrVec &A, Fr id, Fr r, Fr delta
    );
};

struct NonOwnPf {
    ZKAccSetup &setup;
    ZKNMP &piNMP;
    PoK2 &piPoK;

    void init(ZKAccSetup &setup);

    void prove(
        G2 C_id, G2 C_A, Fr n, 
        FrVec &A, Fr id, Fr r, Fr delta
    );
};

bool OwnPf_verify(OwnPf *pi);
bool NonOwnPf_verify(NonOwnPf *pi);

#endif 