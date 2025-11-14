#ifndef TX_GEN_H
#define TX_GEN_H
#include <mcl/bls12_381.hpp>
#include "mainalgorithms.hpp"
#include "polyonfr.hpp"

using namespace chrono;

// Structs for Tx
struct Tx_Base;
struct Tx_Entry;
struct Tx_Exit;
struct Tx_Trans;
struct Tx_Summary;

// Entry Transaction
Tx_Entry Tx_Entry_Gen (
    // Public Input
    PCS &pcs, G2 C_A,
    // Private Input
    FrVec &SNs, Fr id, G1 wit,
    // Provide Internal Randomnesses
    Fr gamma, Fr delta
);

// Exit Transaction
Tx_Exit Tx_Exit_Gen (
    // Public Input
    Tx_Base &Tx_in, PCS &pcs, G2 C_A,
    // Private Input
    FrVec &S, FrVec &S_1, FrVec &S_2, Fr id, G1 wit, Fr gamma, Fr delta,
    // Provide Internal Randomnesses
    Fr gammap, Fr delta1, Fr delta2    
);

// Transfer Transaction
Tx_Trans Tx_Trans_Gen (
    // Public Inputs
    Tx_Base &Tx_in, PCS &pcs, G2 C_A,
    // Private Inputs
    FrVec &S, FrVec &S_1, FrVec &S_2, Fr id1, Fr id2, 
    G1 wit1, G1 wit2, Fr gamma, Fr delta,
    // Provide Internal Randomnesses
    Fr gamma1, Fr gamma2, Fr delta1, Fr delta2    
);

// Summary Transaction
Tx_Summary Tx_Summary_Gen (
    // Public Inputs
    Tx_Base &Tx_in, PCS &pcs, FrVec &SN_loc, G2 C_A, vector<G2> &C_Is,
    // Private Inputs
    Fr gamma_loc, G1 wit_loc, vector<G1> &W_Is, FrVec &gamma_Is,
    // Provide Internal Randomnesses
    Fr gamma_locp
);

// Verification Algorithms

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
