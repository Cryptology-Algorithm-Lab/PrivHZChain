#ifndef _ZKACC_H
#define _ZKACC_H

#include "poly.h"
#include "utils.h"

// PCS
struct PCS {
    struct pp {
        vector<G1> g1si;
        vector<G2> g2si;
        vector<G2> h2si;
        NttCtx nttctx;
    } pp;
    Fr s;
    void setup(uint32_t D);
    G1 commit_g1(const FrVec &p);
    G2 commit_g2(const FrVec &p);
};


// Core ZK proofs for Accumulators
struct ZKAccSetup {
    PCS* pcs;
    G1* g1si;
    G2* g2si;
    G2* h2si;
    G1 g1;
    G2 g2;
    G1 Frakg;
    G1 Frakh;
    G2 h2;
    
    void setup_from_pcs(PCS& pcs);
};


struct ZKMP {
    ZKAccSetup setup;
    G2 C_I;
    G2 C_A;

    ZKMP(ZKAccSetup &setup){
        this->setup=setup;
    }

    void init(ZKAccSetup &setup);

    // Message
    struct msg{
        G1 P_1;
        G1 P_2;
        G1 R_1;
        G1 R_2;
        GT R_3;
    } msg;

    // Response
    struct response{
        Fr s_r_I;
        Fr s_r_A;
        Fr s_tau_1;
        Fr s_tau_2;
        Fr s_delta_1;
        Fr s_delta_2;
    } response;

    // Membership Proof Generation
    void prove(G2 C_I, G2 C_A, G1 pi_I, Fr r_I, Fr r_A);
};

struct AZKMP {
    ZKAccSetup setup;
    vector<G2> C_I;
    G2 C_A;

    AZKMP(ZKAccSetup &setup) {
        this->setup = setup;
    }

    void init(ZKAccSetup &setup);

    struct msg {
        G1 P_1;
        vector<G1> P_2s;
        G1 R_1;
        G1 R_2;
        GT R_3;
    } msg;

    struct response {
        Fr s_r_A;
        Fr s_tau_1;
        Fr s_tau_2;
        FrVec s_r_Is;
        FrVec s_delta_1s;
        FrVec s_delta_2s;
    } response;

    void prove(vector<G2> &C_I, G2 C_A, vector<G1> &pi_I, FrVec r_I, Fr r_A);
};

struct ZKNMP {
    ZKAccSetup setup;
    G2 C_I;
    G2 C_A;

    ZKNMP(ZKAccSetup &setup){
        this->setup=setup;
    }

    void init(ZKAccSetup &setup);

    struct msg{

        G1 P_1;
        G1 P_2;
        G1 Q_1;
        G1 Q_2;
        G1 R_1;
        G1 R_2;
        G1 R_3;
        G1 R_4;
        GT R_5;

    } msg;

    struct response{

        Fr s_r_I;
        Fr s_r_A;
        Fr s_tau_i[4];
        Fr s_delta_i[4];

    }response;

    void prove(G2 C_I, G2 C_A, G1 w_I, G1 w_A, Fr r_I, Fr r_A);
};

struct ZKSP {

    ZKAccSetup setup;
    G2 C_I;
    G2 C_J;
    G2 C_A;

    ZKSP(ZKAccSetup &setup){
        this->setup=setup;
    }

    void init(ZKAccSetup &setup);

    struct msg{

        G1 P_1;
        G1 P_2;
        G1 P_3;
        G1 P_4;
        G1 R_1;
        G1 R_2;
        G1 R_3;
        G1 R_4;
        GT R_5;
        GT R_6;
        GT R_7;

    }msg;

    struct response{

        Fr s_r_I;
        Fr s_r_J;
        Fr s_r_A;
        Fr s_tau_i[4];
        Fr s_delta_i[4];

    }response;

    void prove(G2 C_I, G2 C_J, G2 C_A, G1 pi_I, G1 pi_J, Fr r_I, Fr r_J, Fr r_A);
};


// Verification
bool ZKMP_verify(ZKMP *pi);
bool AZKMP_verify(AZKMP *pi);
bool ZKNMP_verify(ZKNMP *pi);
bool ZKSP_verify(ZKSP *pi);

// Proof of Knowledge for g^x * h^r
struct PoK2 {
    ZKAccSetup setup;
    G2 C; G2 R; Fr z1; Fr z2;

    PoK2(ZKAccSetup &setup) {
        this->setup = setup;
    }

    void init(ZKAccSetup &setup);
    void prove(G2 C, Fr id, Fr r);
};

bool PoK2_verify(PoK2 *pi);

// proof of knowledge for n values (for aggregation)
struct PoKn {
    G2 P; G2 R; FrVec zs;
    void prove(PCS &pcs, G2 P, FrVec x, Fr r);
};

bool PoKn_verify(PCS &pcs, PoKn *pi);


#endif