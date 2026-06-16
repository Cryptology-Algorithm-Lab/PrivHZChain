#ifndef _TXS_H
#define _TXS_H

#include "zkacc.h"
#include <chrono>

using namespace std;
using namespace std::chrono;

// Transaction Structures
struct Tx_Base {
    Fr Adr;
    G2 C;  
    G2 A;  
};

struct Tx_Entry {
    Tx_Base Tx_out;
    ZKMP PoA;
};

struct Tx_Exit {
    Tx_Base Tx_out1;
    Tx_Base Tx_out2;
    PoK2 PoID;
    ZKMP PoA;
    ZKSP PoQ;
};

struct Tx_Trans {
    Tx_Base Tx_out1;
    Tx_Base Tx_out2;
    PoK2 PoID;
    ZKMP PoA1;
    ZKMP PoA2;
    ZKSP PoQ;
};

struct Tx_Summary {
    Tx_Base Tx_out1;
    Tx_Base Tx_out2;
    PoKn PoID;
    AZKMP PoAgg;
    ZKMP PoA;
};

// Transcation Generation
Tx_Entry Tx_Entry_Gen (
    // Public Input
    PCS &pcs, G2 C_A,
    // Private Input
    FrVec &SNs, Fr id, G1 wit,
    // Provide Internal Randomnesses
    Fr gamma, Fr delta
);

Tx_Exit Tx_Exit_Gen (
    // Public Input
    Tx_Base &Tx_in, PCS &pcs, G2 C_A,
    // Private Input
    FrVec &S, FrVec &S_1, FrVec &S_2, Fr id, G1 wit, Fr gamma, Fr delta,
    // Provide Internal Randomnesses
    Fr gammap, Fr delta1, Fr delta2    
);

Tx_Trans Tx_Trans_Gen (
    // Public Inputs
    Tx_Base &Tx_in, PCS &pcs, G2 C_A,
    // Private Inputs
    FrVec &S, FrVec &S_1, FrVec &S_2, Fr id1, Fr id2, 
    G1 wit1, G1 wit2, Fr gamma, Fr delta,
    // Provide Internal Randomnesses
    Fr gamma1, Fr gamma2, Fr delta1, Fr delta2    
);

Tx_Summary Tx_Summary_Gen (
    // Public Inputs
    Tx_Base &Tx_in, PCS &pcs, FrVec &SN_loc, G2 C_A, vector<G2> &C_Is,
    // Private Inputs
    Fr gamma_loc, G1 wit_loc, vector<G1> &W_Is, FrVec &gamma_Is,
    // Provide Internal Randomnesses
    Fr gamma_locp
);

// Transaction Verification
bool Tx_Entry_Vrfy (
    Tx_Entry &tx, PCS &pcs
);

bool Tx_Exit_Vrfy (
    Tx_Exit &tx, PCS &pcs
);

bool Tx_Trans_Vrfy (
    Tx_Trans &tx, PCS &pcs
);

bool Tx_Summary_Vrfy (
    Tx_Summary &tx, PCS &pcs
);

#endif 