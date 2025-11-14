#include "TxGen.hpp"

#define DEBUG 0
#define VERBOSE 1

struct Tx_Base {
    Fr Adr;
    G2 C;  
    G2 A;  
};

// Entry Transction

struct Tx_Entry {
    Tx_Base Tx_out;
    ZKMP PoA;
};


Tx_Entry Tx_Entry_Gen (
    // Public Input
    PCS &pcs, G2 C_A,
    // Private Input
    FrVec &SNs, Fr id, G1 wit,
    // Provide Internal Randomnesses
    Fr gamma, Fr delta
) {
    // TxCom
    auto time1 = chrono::steady_clock::now();

    G2 C = pcs.commit({id, 1}) + pcs.pp.h2si[0] * gamma;
    G2 A = pcs.commit(constructMemPoly(SNs)) + pcs.pp.h2si[0] * delta;
    string buf = C.getStr() + A.getStr(); 
    Fr Adr; Adr.setHashOf(buf);
    Tx_Base Tx_out = {Adr, C, A};

    auto time2 = chrono::steady_clock::now();

    // TxPrv
    zkbpacc_setup setup;
    setup.setup_from_pcs(pcs);
    ZKMP PoA(setup); PoA.prove(C, C_A, wit, gamma, Fr(0));

    auto time3 = chrono::steady_clock::now();


#ifdef VERBOSE
    auto t_Com = duration_cast<microseconds>(time2-time1).count();
    auto t_Prv = duration_cast<microseconds>(time3-time2).count();

    cout << "TxEntryGEN - \t\t"
         << "TxCom: " << t_Com << "us \t\t"
         << "TxPrv: " << t_Prv << "us" << endl;
#endif 

    // Done
    return Tx_Entry {
        Tx_out, PoA
    };
}

bool Tx_Entry_Vrfy (
    Tx_Entry &tx, PCS &pcs
) {
    // Check the PoA
    bool checkPoA = ZKMP_verify(&tx.PoA);

#if DEBUG    
    cout << checkPoA << endl;
#endif

    return (checkPoA);
}

// Exit Transaction

struct Tx_Exit {
    Tx_Base Tx_out1;
    Tx_Base Tx_out2;
    PoK2 PoID;
    ZKMP PoA;
    ZKSP PoQ;
};

Tx_Exit Tx_Exit_Gen (
    // Public Input
    Tx_Base &Tx_in, PCS &pcs, G2 C_A,
    // Private Input
    FrVec &S, FrVec &S_1, FrVec &S_2, Fr id, G1 wit, Fr gamma, Fr delta,
    // Provide Internal Randomnesses
    Fr gammap, Fr delta1, Fr delta2    
) {

    auto time0 = chrono::steady_clock::now();

    // Line 1
    G1 w1 = pcs.commit_G1(constructMemPolyBatch(S, S_1));
    G1 w2 = pcs.commit_G1(constructMemPolyBatch(S, S_2));

    // TxCom
    auto time1 = chrono::steady_clock::now();

    G2 Cp = pcs.commit({id, 1}) + pcs.pp.h2si[0] * gammap;
    G2 A1 = pcs.commit(constructMemPoly(S_1)) + pcs.pp.h2si[0] * delta1;
    G2 A2 = pcs.commit(constructMemPoly(S_2)) + pcs.pp.h2si[0] * delta2;

    string buf = Cp.getStr() + A2.getStr();
    Fr Adr; Adr.setHashOf(buf);
    Tx_Base Tx_out1; Tx_out1.A = A1;
    Tx_Base Tx_out2 = {Adr, Cp, A2};

    auto time2 = chrono::steady_clock::now();

    // TxPrv
    zkbpacc_setup setup; setup.setup_from_pcs(pcs);
    PoK2 PoID;
    ZKMP PoA(setup); ZKSP PoQ(setup);
    PoID.prove(pcs, Tx_in.C - pcs.pp.g2si[1], id, gamma);
    PoA.prove(Cp, C_A, wit, gammap, Fr(0));
    PoQ.prove(A1, A2, Tx_in.A, w1, w2, delta1, delta2, delta);

    auto time3 = chrono::steady_clock::now();    


#ifdef VERBOSE
    auto t_MWG = duration_cast<microseconds>(time1-time0).count();
    auto t_Com = duration_cast<microseconds>(time2-time1).count();
    auto t_Prv = duration_cast<microseconds>(time3-time2).count();

    cout << "TxExitGEN - \t\t"
         << "MWG: " << t_MWG << "us \t\t"
         << "TxCom: " << t_Com << "us \t\t"
         << "TxPrv: " << t_Prv << "us" << endl;
#endif 

    return Tx_Exit {
        Tx_out1, Tx_out2,
        PoID, PoA, PoQ
    };
}

bool Tx_Exit_Vrfy (
    Tx_Exit &tx, PCS &pcs
) {
    // Check PoID
    bool checkPoID = Pok2_verify(pcs, &tx.PoID);

    // Check PoA
    bool checkPoA = ZKMP_verify(&tx.PoA);

    // Check PoQ
    bool checkPoQ = ZKSP_verify(&tx.PoQ);

#if DEBUG        
    cout << checkPoID << checkPoA << checkPoQ << endl;    
#endif    

    return (checkPoID & checkPoA & checkPoQ);
}


// Transfer Transaction
struct Tx_Trans {
    Tx_Base Tx_out1;
    Tx_Base Tx_out2;
    PoK2 PoID;
    ZKMP PoA1;
    ZKMP PoA2;
    ZKSP PoQ;
};

Tx_Trans Tx_Trans_Gen (
    // Public Inputs
    Tx_Base &Tx_in, PCS &pcs, G2 C_A,
    // Private Inputs
    FrVec &S, FrVec &S_1, FrVec &S_2, Fr id1, Fr id2, 
    G1 wit1, G1 wit2, Fr gamma, Fr delta,
    // Provide Internal Randomnesses
    Fr gamma1, Fr gamma2, Fr delta1, Fr delta2    
) {
    auto time0 = chrono::steady_clock::now();

    // Line 1
    G1 w1 = pcs.commit_G1(constructMemPolyBatch(S, S_1));
    G1 w2 = pcs.commit_G1(constructMemPolyBatch(S, S_2));    

    // TxCom
    auto time1 = chrono::steady_clock::now();

    G2 C1 = pcs.commit({id1, 1}) + pcs.pp.h2si[0] * gamma1;
    G2 C2 = pcs.commit({id2, 1}) + pcs.pp.h2si[0] * gamma2;
    G2 A1 = pcs.commit(constructMemPoly(S_1)) + pcs.pp.h2si[0] * delta1;
    G2 A2 = pcs.commit(constructMemPoly(S_2)) + pcs.pp.h2si[0] * delta2;
    string buf1 = C1.getStr() + A1.getStr();
    string buf2 = C2.getStr() + A2.getStr();
    Fr Adr1, Adr2; Adr1.setHashOf(buf1); Adr2.setHashOf(buf2);
    Tx_Base Tx_out1 = {Adr1, C1, A1}; Tx_Base Tx_out2 = {Adr2, C2, A2};

    auto time2 = chrono::steady_clock::now();

    // TxPrv
    zkbpacc_setup setup; setup.setup_from_pcs(pcs);
    PoK2 PoID;
    ZKMP PoA1(setup); ZKMP PoA2(setup); ZKSP PoQ(setup);
    PoID.prove(pcs, Tx_in.C - pcs.pp.g2si[1], id1, gamma);
    PoA1.prove(C1, C_A, wit1, gamma1, Fr(0));
    PoA2.prove(C2, C_A, wit2, gamma2, Fr(0));
    PoQ.prove(A1, A2, Tx_in.A, w1, w2, delta1, delta2, delta);
    
    auto time3 = chrono::steady_clock::now();

#ifdef VERBOSE
    auto t_MWG = duration_cast<microseconds>(time1-time0).count();
    auto t_Com = duration_cast<microseconds>(time2-time1).count();
    auto t_Prv = duration_cast<microseconds>(time3-time2).count();

    cout << "TxTransGEN - \t\t"
         << "MWG: " << t_MWG << "us \t\t"
         << "TxCom: " << t_Com << "us \t\t"
         << "TxPrv: " << t_Prv << "us" << endl;
#endif 

    return Tx_Trans {
        Tx_out1, Tx_out2,
        PoID, PoA1, PoA2, PoQ
    };
}

bool Tx_Trans_Vrfy (
    Tx_Trans &tx, PCS &pcs
) {
    // Check PoID
    bool checkPoID = Pok2_verify(pcs, &tx.PoID);

    // Check PoA
    bool checkPoA1 = ZKMP_verify(&tx.PoA1);
    bool checkPoA2 = ZKMP_verify(&tx.PoA2);

    // Check PoQ
    bool checkPoQ = ZKSP_verify(&tx.PoQ);

#if DEBUG        
    cout << checkPoID << checkPoA1 << checkPoA2 << checkPoQ << endl;
#endif

    return (checkPoID & checkPoA1 & checkPoA2 & checkPoQ);
}


// Summary Transaction

struct Tx_Summary {
    Tx_Base Tx_out1;
    Tx_Base Tx_out2;
    PoKn PoID;
    AZKMP PoAgg;
    ZKMP PoA;
};

Tx_Summary Tx_Summary_Gen (
    // Public Inputs
    Tx_Base &Tx_in, PCS &pcs, FrVec &SN_loc, G2 C_A, vector<G2> &C_Is,
    // Private Inputs
    Fr gamma_loc, G1 wit_loc, vector<G1> &W_Is, FrVec &gamma_Is,
    // Provide Internal Randomnesses
    Fr gamma_locp
) {
    uint32_t numComs = C_Is.size();

    // TxCom
    auto time1 = chrono::steady_clock::now();

    FrVec SN_memPoly = constructMemPoly(SN_loc);
    G2 C_locp = pcs.commit(SN_memPoly) + pcs.pp.h2si[0] * gamma_locp;
    string buf = C_locp.getStr();
    for (uint32_t i = 0; i < numComs; i++) {
        buf += C_Is[i].getStr();
    }
    buf += Tx_in.A.getStr();
    Fr Adr; Adr.setHashOf(buf);
    Tx_Base Tx_out1 = {Adr, C_locp, Tx_in.A}; Tx_Base Tx_out2;

    auto time2 = chrono::steady_clock::now();    

    // TxPrv
    zkbpacc_setup setup; setup.setup_from_pcs(pcs);    
    PoKn PoID;
    AZKMP PoAgg(setup);
    ZKMP PoA(setup);

    PoID.prove(pcs, Tx_in.C, SN_memPoly, gamma_loc);
    PoAgg.prove(C_Is, C_locp, W_Is, gamma_Is, gamma_locp);
    PoA.prove(C_locp, C_A, wit_loc, gamma_locp, Fr(0));

    auto time3 = chrono::steady_clock::now();    


#ifdef VERBOSE
    auto t_Com = duration_cast<microseconds>(time2-time1).count();
    auto t_Prv = duration_cast<microseconds>(time3-time2).count();

    cout << "TxSummary - \t\t"
         << "TxCom: " << t_Com << "us \t\t"
         << "TxPrv: " << t_Prv << "us" << endl;
#endif 

    return Tx_Summary {
        Tx_out1, Tx_out2,
        PoID, PoAgg, PoA
    };
}

bool Tx_Summary_Vrfy (
    Tx_Summary &tx, PCS &pcs
) {
    // Check PoID
    bool checkPoID = PoKn_verify(pcs, &tx.PoID);

    // Check PoA
    bool checkPoA = ZKMP_verify(&tx.PoA);

    // Check PoAgg
    bool checkAgg = AZKMP_verify(&tx.PoAgg);

#if DEBUG        
    cout << checkPoID << checkPoA << checkAgg << endl;
#endif    

    return (checkPoID & checkPoA & checkAgg);
}
