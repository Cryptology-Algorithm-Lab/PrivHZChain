#ifndef TX_GEN_H
#define TX_GEN_H
#include <mcl/bls12_381.hpp>
#include "mainalgorithms.hpp"

struct Tx_Entry;
struct Entry_SN;
struct output;
struct Transfer_input;
struct Tx_Transfer;
struct Tx_Exit;

Tx_Entry Tx_Entry_Gen(PCS pcs, ZKMP PI_PoA, Fr id, G1 wit, G2 C_A, FrVec SNs, Fr gamma, Fr delta, const string& file_name);

Tx_Transfer Tx_Transfer_Gen(
    PCS pcs, ZKMP PI_PoA_1, ZKMP PI_PoA_2, ZKSP PI_PoQ, Fr id1, Fr id2, G1 wit1, G1 wit2, G2 C_A, FrVec SNs, FrVec S_1, FrVec S_2, output output_old, 
    Fr gamma, Fr gamma1, Fr gamma2, Fr delta, Fr delta1, Fr delta2, const std::string& file_name, Fr Tid_old);


Tx_Exit Tx_Exit_Gen(PCS pcs, output output_old, ZKMP PI_PoA, ZKSP PI_PoQ, Fr id, G1 wit, G2 C_A, Fr gamma, Fr gamma1, Fr delta, Fr delta1, FrVec S, FrVec M, FrVec N, const std::string& file_name);

#endif