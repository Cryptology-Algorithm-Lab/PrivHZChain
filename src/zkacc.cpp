#include "zkacc.h"

// PCS
void PCS::setup(uint32_t D) {
    pp.nttctx.init(nextPow2(D+1));

    Fr s;
    s.setByCSPRNG();
    this->s = s;

    G1 g1; G2 g2; G2 h2;
    hashAndMapToG1(g1,"PrivHZChain-g1");
    hashAndMapToG2(g2,"PrivHZChain-g2");
    hashAndMapToG2(h2,"PrivHZChain-h2");

    vector<G1> g1si(D+1);
    vector<G2> g2si(D+1);
    vector<G2> h2si(D+1);

    g1si[0] = g1;
    g2si[0] = g2;
    h2si[0] = h2;

    vector<Fr> spow(D+1);
    spow[0] = 1;
    for (uint32_t i = 1; i <= D; i++) {
        spow[i] = spow[i-1] * s; 
    }

    // Somewhat empirically found values
    int bitwidth;
    if (D < 32768) {
        bitwidth = 4;
    } else {
        bitwidth = 8;
    }

    // Optimization for fixed-base multiplication
    mcl::fp::WindowMethod<G1> wmG1; wmG1.init(g1, 256, bitwidth);
    mcl::fp::WindowMethod<G2> wmG2; wmG2.init(g2, 256, bitwidth);
    mcl::fp::WindowMethod<G2> wmH2; wmH2.init(h2, 256, bitwidth);    

    // Use Multithreading
    // This is done by the server.
    #pragma omp parallel for
    for (uint32_t i = 0; i <= D; i++) {
        wmG1.mul(g1si[i], spow[i]);
        wmG2.mul(g2si[i], spow[i]);
        wmH2.mul(h2si[i], spow[i]);
    }    
    this->pp.g1si = g1si;
    this->pp.g2si = g2si;
    this->pp.h2si = h2si;
}

G1 PCS::commit_g1(const FrVec &p) {
    G1 ret; ret.clear();
    G1::mulVec(ret, &pp.g1si[0], &p[0], p.size());
    return ret;
}

G2 PCS::commit_g2(const FrVec &p) {
    G2 ret; ret.clear();
    G2::mulVec(ret, &pp.g2si[0], &p[0], p.size());
    return ret;
}

// Setup 
void ZKAccSetup::setup_from_pcs(PCS& pcs) {
    this->pcs = &pcs;
    g1si = &pcs.pp.g1si[0];
    g2si = &pcs.pp.g2si[0];
    h2si = &pcs.pp.h2si[0];
    g1 = pcs.pp.g1si[0];
    g2 = pcs.pp.g2si[0];
    h2 = pcs.pp.h2si[0];
    hashAndMapToG1(Frakg,"PrivHZChain-Frakg");
    hashAndMapToG1(Frakh,"PrivHZChain-Frakh");
}

// ZKMP
void ZKMP::init(ZKAccSetup &setup) {
    this->setup = setup;
}

void ZKMP::prove(G2 C_I, G2 C_A, G1 pi_I, Fr r_I, Fr r_A) {

    this->C_I=C_I;
    this->C_A=C_A;

    // Assign Scalar Values
    Fr tau_1; Fr tau_2;
    Fr r_r_I; Fr r_r_A;
    Fr r_tau_1; Fr r_tau_2;
    Fr r_delta_1; Fr r_delta_2;

    tau_1.setByCSPRNG(); tau_2.setByCSPRNG();
    r_r_I.setByCSPRNG(); r_r_A.setByCSPRNG();
    r_tau_1.setByCSPRNG(); r_tau_2.setByCSPRNG();
    r_delta_1.setByCSPRNG(); r_delta_2.setByCSPRNG();

    Fr delta_1 = r_I*tau_1; Fr delta_2 = r_I*tau_2;

    // Read Setup Values
    G1 g1 = setup.g1;
    G1 Frakg = setup.Frakg;
    G2 h2 = this->setup.h2;

    // Compute P, R
    G1 P_1 = (g1 * tau_1) + (Frakg * tau_2);
    G1 P_2 = pi_I + (Frakg * tau_1);
    G1 R_1 = (g1 * r_tau_1)+(Frakg * r_tau_2);
    G1 R_2 = (P_1 * r_r_I) + (g1 * (-r_delta_1)) + (Frakg * (-r_delta_2));
    
    // Compute R_3
    G2 first = C_I * r_tau_1 + h2 * (-r_delta_1);
    G1 second = P_2 * r_r_I + g1 * (-r_r_A);
    G1 leftMP[2] = {Frakg, second};
    G2 rightMP[2] = {first, h2};
    GT R_3;
    millerLoopVec(R_3, leftMP, rightMP, 2);
    finalExp(R_3, R_3);

    // Make a challenge
    string buf=P_1.getStr()+P_2.getStr()+R_1.getStr()+R_2.getStr()+R_3.getStr();
    Fr c;
    c.setHashOf(buf);

    // Assign Values
    // Message
    this->msg.P_1 = P_1;
    this->msg.P_2 = P_2;
    this->msg.R_1 = R_1;
    this->msg.R_2 = R_2;
    this->msg.R_3 = R_3;

    // Response
    this->response.s_r_I=r_r_I+c*r_I;
    this->response.s_r_A=r_r_A+c*r_A;
    this->response.s_tau_1=r_tau_1+c*tau_1;
    this->response.s_tau_2=r_tau_2+c*tau_2;
    this->response.s_delta_1=r_delta_1+c*delta_1;
    this->response.s_delta_2=r_delta_2+c*delta_2;
}

bool ZKMP_verify(ZKMP *pi) {

    string buf=pi->msg.P_1.getStr()+pi->msg.P_2.getStr()+pi->msg.R_1.getStr()+pi->msg.R_2.getStr()+pi->msg.R_3.getStr();
    Fr c;
    c.setHashOf(buf);

    G1 base[3] = {pi->msg.P_1, pi->setup.g1, pi->setup.Frakg};

    // Check R1
    G1 right1;
    Fr expon1[3] = {-c, pi->response.s_tau_1, pi->response.s_tau_2};
    G1::mulVec(right1, base, expon1, 3);

    // Check R2
    G1 right2;
    Fr expon2[3] = {pi->response.s_r_I, -pi->response.s_delta_1, -pi->response.s_delta_2};
    G1::mulVec(right2, base, expon2, 3);

    // Check R3
    G1 leftMP[3] = {pi->setup.Frakg, pi->msg.P_2, pi->setup.g1};
    G2 firstTerm = pi->C_I * pi->response.s_tau_1 + pi->setup.h2 * (-pi->response.s_delta_1);
    G2 secondTerm = pi->setup.h2 * pi->response.s_r_I + pi->C_I * (-c);
    G2 thirdTerm = pi->setup.h2 * (-pi->response.s_r_A) + pi->C_A * c;
    G2 rightMP[3] = {firstTerm, secondTerm, thirdTerm};

    GT right3;
    millerLoopVec(right3, leftMP, rightMP, 3);
    finalExp(right3, right3);

    bool flag1 = (pi->msg.R_1 == right1);
    bool flag2 = (pi->msg.R_2 == right2);
    bool flag3 = (pi->msg.R_3 == right3);

    return flag1 && flag2 && flag3;
};


// AZKMP
void AZKMP::init(ZKAccSetup &setup) {
    this->setup = setup;
}

void AZKMP::prove(vector<G2> &C_I, G2 C_A, vector<G1> &pi_I, FrVec r_I, Fr r_A) {
    uint32_t n = C_I.size();

    this->C_I = C_I;
    this->C_A = C_A;

    // Sample Randoms
    Fr tau_1; Fr tau_2; Fr r_r_A; Fr r_tau_1; Fr r_tau_2;
    FrVec r_r_is; FrVec r_delta_1s; FrVec r_delta_2s;

    r_r_is.resize(n); r_delta_1s.resize(n); r_delta_2s.resize(n);

    for (uint32_t i = 0; i < n; i++) {
        r_r_is[i].setByCSPRNG();
        r_delta_1s[i].setByCSPRNG();
        r_delta_2s[i].setByCSPRNG();         
    }    

    tau_1.setByCSPRNG(); tau_2.setByCSPRNG();
    r_r_A.setByCSPRNG(); r_tau_1.setByCSPRNG(); r_tau_2.setByCSPRNG();


    // Assign Memory for Deltas 
    FrVec delta_1s; FrVec delta_2s;
    delta_1s.resize(n); delta_2s.resize(n);

    // Compute P1 & P2
    G1 Frakg_tau_1 = this->setup.Frakg * tau_1;
    this->msg.P_1 = (this->setup.g1 * tau_1) + (this->setup.Frakg * tau_2);        
    this->msg.P_2s.resize(n);

    // Merged Loops for Parallelization Efficiency
    for (uint32_t i = 0; i < n; i++) {
        delta_1s[i] = (r_I[i] * tau_1);
        delta_2s[i] = (r_I[i] * tau_2);        
        this->msg.P_2s[i] = (pi_I[i] + Frakg_tau_1);
    }

    // Obtain a challenge alpha
    string buf = this->msg.P_1.getStr();
    for (uint32_t i = 0; i < n; i++) {
        buf += this->msg.P_2s[i].getStr();
    }
    Fr alpha;
    alpha.setHashOf(buf);

    // Prepare alpha's powers
    FrVec alphaPows;
    alphaPows.push_back(Fr(1));
    for (uint32_t i = 1; i < n; i++) {
        alphaPows.push_back(alphaPows[i-1] * alpha);
    }

    // Compute R1
    this->msg.R_1 = (this->setup.g1 * r_tau_1) + (this->setup.Frakg * r_tau_2);

    // Compute R2
    // Step 1. Aggregate powers
    Fr agg_r_r_I = innerProd(alphaPows, r_r_is);
    Fr agg_r_delta_1s = innerProd(alphaPows, r_delta_1s);
    Fr agg_r_delta_2s = innerProd(alphaPows, r_delta_2s);
    
    // Step 2. Compute R2
    this->msg.R_2 = (this->msg.P_1 * agg_r_r_I) + (this->setup.g1 * (- agg_r_delta_1s)) + (this->setup.Frakg * (- agg_r_delta_2s));

    // Compute R3
    // IDEA: Merge MSM operations in a part
    // Prepare Bases
    G2 first; vector<G2> firstMSMBase; firstMSMBase.resize(n+1);
    G1 second; vector<G1> secondMSMBase; secondMSMBase.resize(n+1);
    for (uint32_t i = 0; i < n; i++) {
        firstMSMBase[i] = C_I[i];
    }    
    for (uint32_t i = 0; i < n; i++) {
        secondMSMBase[i] = (this->msg.P_2s[i]);
    }        
    firstMSMBase[n] = (this->setup.h2);
    secondMSMBase[n] = (this->setup.g1);

    // Prepare Exponents
    FrVec firstMSMExpon; firstMSMExpon.resize(n+1);
    FrVec secondMSMExpon; secondMSMExpon.resize(n+1);
    for (uint32_t i = 0; i < n; i++) {
        firstMSMExpon[i] = (alphaPows[i] * r_tau_1);
        secondMSMExpon[i] = (alphaPows[i] * r_r_is[i]);
    }
    firstMSMExpon[n] = (
        (- innerProd(alphaPows, r_delta_1s))
    );
    secondMSMExpon[n] = (
        (- sumAllComps(alphaPows)) * r_r_A
    );
    // Final MSM for each part
    G2::mulVecMT(first, &firstMSMBase[0], &firstMSMExpon[0], n+1);
    G1::mulVecMT(second, &secondMSMBase[0], &secondMSMExpon[0], n+1);
    
    // Do a final pairing loop
    G1 leftMP[2] = {this->setup.Frakg, second};
    G2 rightMP[2] = {first, this->setup.h2};
    GT R3; 
    millerLoopVec(R3, leftMP, rightMP, 2);
    finalExp(R3, R3);
    this->msg.R_3 = R3;

    buf += this->msg.R_1.getStr();
    buf += this->msg.R_2.getStr();
    buf += this->msg.R_3.getStr();

    Fr c;
    c.setHashOf(buf);

    // Prepare Responses
    this->response.s_r_A = r_r_A + c * r_A;
    this->response.s_tau_1 = r_tau_1 + c * tau_1;
    this->response.s_tau_2 = r_tau_2 + c * tau_2;

    this->response.s_r_Is.resize(n);
    this->response.s_delta_1s.resize(n);
    this->response.s_delta_2s.resize(n);

    for (uint32_t i = 0; i < n; i++) {
        this->response.s_r_Is[i] = (
            r_r_is[i] + c * r_I[i]
        );
        this->response.s_delta_1s[i] = (
            r_delta_1s[i] + c * delta_1s[i]
        );
        this->response.s_delta_2s[i] = (
            r_delta_2s[i] + c * delta_2s[i]
        );                        
    }
    // Done!
}



bool AZKMP_verify(AZKMP *pi) {
    uint32_t n = pi->msg.P_2s.size();

    // Retrieve Challenges
    string buf = pi->msg.P_1.getStr();
    for (uint32_t i = 0; i < n; i++) {
        buf += pi->msg.P_2s[i].getStr();
    }


    Fr alpha;
    alpha.setHashOf(buf);

    buf += pi->msg.R_1.getStr();
    buf += pi->msg.R_2.getStr();
    buf += pi->msg.R_3.getStr();

    Fr c;
    c.setHashOf(buf);

    // Prepare alpha's powers
    FrVec alphaPows;
    alphaPows.push_back(Fr(1));
    for (uint32_t i = 1; i < n; i++) {
        alphaPows.push_back(alphaPows[i-1] * alpha);
    }    

    // Check Relations
    G1 base[3] = {pi->msg.P_1, pi->setup.g1, pi->setup.Frakg};

    // First Relation
    G1 right1;
    Fr expon1[3] = {-c, pi->response.s_tau_1, pi->response.s_tau_2};
    G1::mulVec(right1, base, expon1, 3);
    
    // Second Relation
    G1 right2;
    Fr expon2[3] = {
        innerProd(alphaPows, pi->response.s_r_Is),
        -innerProd(alphaPows, pi->response.s_delta_1s),
        -innerProd(alphaPows, pi->response.s_delta_2s)
    };
    G1::mulVec(right2, base, expon2, 3);

    // Third Relation
    GT right3;
    vector<G1> leftMP; leftMP.resize(n+2);
    vector<G2> rightMP; rightMP.resize(n+2);
    for (uint32_t i = 0; i < n; i++) {
        leftMP[i] = (
            (pi->setup.Frakg * (pi->response.s_tau_1 * alphaPows[i])) + \
            (pi->msg.P_2s[i] * ((-c) * alphaPows[i]))
        );
        rightMP[i] = (pi->C_I[i]);
    }

    // First "Easy Term"
    Fr sumAlpha = sumAllComps(alphaPows);
    leftMP[n] = (pi->setup.g1);
    rightMP[n] = (
        (pi->C_A * (c * sumAlpha)) + \    
        (pi->setup.h2 * ((-pi->response.s_r_A) * sumAlpha ))
    );

    // Second "Easy Term"
    vector<G1> _tmp_base; _tmp_base.resize(n + 1);
    for (uint32_t i = 0; i < n; i++) {
        _tmp_base[i] = (pi->msg.P_2s[i]);
    }
    _tmp_base[n] = (pi->setup.Frakg);

    FrVec _tmp_expon = hadProduct(alphaPows, pi->response.s_r_Is);
    _tmp_expon.push_back(-innerProd(pi->response.s_delta_1s, alphaPows));

    G1 _tmp_MP_base;
    G1::mulVec(_tmp_MP_base, &_tmp_base[0], &_tmp_expon[0], n+1);

    leftMP[n+1] = (_tmp_MP_base);
    rightMP[n+1] = (pi->setup.h2);    

    // Do a final multipairing
    millerLoopVec(right3, &leftMP[0], &rightMP[0], n+2);
    finalExp(right3, right3);

    bool flag1 = (pi->msg.R_1 == right1);
    bool flag2 = (pi->msg.R_2 == right2);
    bool flag3 = (pi->msg.R_3 == right3);

    return flag1 && flag2 && flag3;
}

// Non-membership Proof
void ZKNMP::init(ZKAccSetup &setup) {
    this->setup = setup;
}

void ZKNMP::prove(G2 C_I, G2 C_A, G1 w_I, G1 w_A, Fr r_I, Fr r_A) {
    this->C_I = C_I;
    this->C_A = C_A;

    Fr r_r_I;
    Fr r_r_A;
    Fr tau_i[4];
    Fr r_tau_i[4];
    Fr r_delta_i[4];
    Fr delta_i[4];
    r_r_I.setByCSPRNG();
    r_r_A.setByCSPRNG();

    for(int i=0;i<4;i++){

        tau_i[i].setByCSPRNG();
        r_tau_i[i].setByCSPRNG();
        r_delta_i[i].setByCSPRNG();

    }

    delta_i[0]=r_I*tau_i[0];
    delta_i[1]=r_I*tau_i[1];
    delta_i[2]=r_A*tau_i[2];
    delta_i[3]=r_A*tau_i[3];

    this->msg.P_1 = this->setup.g1*tau_i[0]+this->setup.Frakg*tau_i[1];
    this->msg.P_2 = w_I+this->setup.Frakg*tau_i[0];
    this->msg.Q_1 = this->setup.g1*tau_i[2]+this->setup.Frakh*tau_i[3];
    this->msg.Q_2 = w_A+this->setup.Frakh*tau_i[2];
    this->msg.R_1 = this->setup.g1*r_tau_i[0]+this->setup.Frakg*r_tau_i[1];
    this->msg.R_2 = this->msg.P_1*r_r_I+this->setup.g1*(-r_delta_i[0])+this->setup.Frakg*(-r_delta_i[1]);
    this->msg.R_3 = this->setup.g1*r_tau_i[2]+this->setup.Frakh*r_tau_i[3];
    this->msg.R_4 = this->msg.Q_1*r_r_A+this->setup.g1*(-r_delta_i[2])+this->setup.Frakh*(-r_delta_i[3]);

    GT R_5;
    G1 first;
    G1 base[4] = {
        this->msg.P_2,
        this->msg.Q_2,
        this->setup.Frakg,
        this->setup.Frakh
    };
    Fr expon[4] = {r_r_I, r_r_A, -r_delta_i[0], -r_delta_i[2]};
    G1::mulVec(first, base, expon, 4);
    G1 leftMP[3] = {this->setup.Frakg * r_tau_i[0], this->setup.Frakh * r_tau_i[2], first};
    G2 rightMP[3] = {C_I,C_A,this->setup.h2};
    millerLoopVec(R_5, leftMP, rightMP, 3);
    finalExp(R_5, R_5);
    this->msg.R_5 = R_5;

    string buf=this->msg.P_1.getStr()+this->msg.P_2.getStr()+this->msg.Q_1.getStr()+this->msg.Q_2.getStr()+this->msg.R_1.getStr()+this->msg.R_2.getStr()
    +this->msg.R_3.getStr()+this->msg.R_4.getStr()+this->msg.R_5.getStr();
    Fr c;
    c.setHashOf(buf);      

    this->response.s_r_I=r_r_I+c*r_I;
    this->response.s_r_A=r_r_A+c*r_A;
    for(int i=0;i<4;i++){

        this->response.s_tau_i[i]=r_tau_i[i]+c*tau_i[i];
        this->response.s_delta_i[i]=r_delta_i[i]+c*delta_i[i];

    }
}


bool ZKNMP_verify(ZKNMP *pi) {
    string buf=pi->msg.P_1.getStr()+pi->msg.P_2.getStr()+pi->msg.Q_1.getStr()+pi->msg.Q_2.getStr()+pi->msg.R_1.getStr()+pi->msg.R_2.getStr()
    +pi->msg.R_3.getStr()+pi->msg.R_4.getStr()+pi->msg.R_5.getStr();
    Fr c;
    c.setHashOf(buf);

    G1 base1[3] = {pi->msg.P_1, pi->setup.g1, pi->setup.Frakg};
    G1 base2[3] = {pi->msg.Q_1, pi->setup.g1, pi->setup.Frakh};

    Fr expon1[3] = {-c, pi->response.s_tau_i[0], pi->response.s_tau_i[1]};
    Fr expon2[3] = {pi->response.s_r_I, -pi->response.s_delta_i[0], -pi->response.s_delta_i[1]};
    Fr expon3[3] = {-c, pi->response.s_tau_i[2], pi->response.s_tau_i[3]};
    Fr expon4[3] = {pi->response.s_r_A, -pi->response.s_delta_i[2], -pi->response.s_delta_i[3]};

    G1 R_1v, R_2v, R_3v, R_4v;
    G1::mulVec(R_1v, base1, expon1, 3);
    G1::mulVec(R_2v, base1, expon2, 3);
    G1::mulVec(R_3v, base2, expon3, 3);
    G1::mulVec(R_4v, base2, expon4, 3);

    GT R_5v;
    G1 first;
    G1 base[4] = {
        pi->msg.P_2, 
        pi->msg.Q_2,
        pi->setup.Frakg, 
        pi->setup.Frakh
    };
    Fr expon[4] = {
        pi->response.s_r_I,
        pi->response.s_r_A,
        -pi->response.s_delta_i[0],
        -pi->response.s_delta_i[2],
    };
    G1::mulVec(first, base, expon, 4);

    G1 leftMP[4] = {
        pi->setup.g1 * (c),
        pi->setup.Frakg * pi->response.s_tau_i[0] + pi->msg.P_2 * (-c),
        pi->setup.Frakh * pi->response.s_tau_i[2] + pi->msg.Q_2 * (-c),
        first,
    };
    G2 rightMP[4] = {
        pi->setup.g2,
        pi->C_I,
        pi->C_A,
        pi->setup.h2
    };

    millerLoopVec(R_5v, leftMP, rightMP, 4);
    finalExp(R_5v, R_5v);

    bool flag1 = (pi->msg.R_1 == R_1v);
    bool flag2 = (pi->msg.R_2 == R_2v);
    bool flag3 = (pi->msg.R_3 == R_3v);
    bool flag4 = (pi->msg.R_4 == R_4v);
    bool flag5 = (pi->msg.R_5 == R_5v);
    
    return flag1 && flag2 && flag3 && flag4 && flag5;
}

// Split Proof
void ZKSP::init(ZKAccSetup &setup) {
    this->setup = setup;
}

void ZKSP::prove(G2 C_I, G2 C_J, G2 C_A, G1 pi_I, G1 pi_J, Fr r_I, Fr r_J, Fr r_A){

    this->C_I=C_I;
    this->C_J=C_J;
    this->C_A=C_A;

    Fr r_r_I;
    Fr r_r_J;
    Fr r_r_A;
    Fr tau_i[4];
    Fr r_tau_i[4];
    Fr r_delta_i[4];
    Fr delta_i[4];
    r_r_I.setByCSPRNG();
    r_r_A.setByCSPRNG();
    r_r_J.setByCSPRNG();

    for(int i=0;i<4;i++){

        tau_i[i].setByCSPRNG();
        r_tau_i[i].setByCSPRNG();
        r_delta_i[i].setByCSPRNG();

    }

    delta_i[0]=r_I*tau_i[0];
    delta_i[1]=r_I*tau_i[1];
    delta_i[2]=r_J*tau_i[2];
    delta_i[3]=r_J*tau_i[3];

    this->msg.P_1 = ((this->setup.g1)*tau_i[0])+((this->setup.Frakg)*tau_i[1]);
    this->msg.P_2 = pi_I+((this->setup.Frakg)*tau_i[0]);
    this->msg.P_3 = ((this->setup.g1)*tau_i[2])+((this->setup.Frakg)*tau_i[3]);
    this->msg.P_4 = pi_J+((this->setup.Frakg)*tau_i[2]);

    this->msg.R_1 = ((this->setup.g1)*r_tau_i[0])+((this->setup.Frakg)*r_tau_i[1]);
    this->msg.R_2 = ((this->msg.P_1)*r_r_I)+((this->setup.g1)*(-r_delta_i[0]))+((this->setup.Frakg)*(-r_delta_i[1]));
    this->msg.R_3 = ((this->setup.g1)*r_tau_i[2])+((this->setup.Frakg)*r_tau_i[3]);
    this->msg.R_4 = ((this->msg.P_3)*r_r_J)+((this->setup.g1)*(-r_delta_i[2]))+((this->setup.Frakg)*(-r_delta_i[3]));

    G1 first;
    G1 base5[3] = {
        this->msg.P_2,
        this->setup.Frakg,
        this->setup.g1
    };
    Fr expon5[3] = {
        r_r_I, -r_delta_i[0], -r_r_A
    };
    G1::mulVec(first, base5, expon5, 3);
    G1 leftMP5[2] = {this->setup.Frakg * r_tau_i[0], first};
    G2 rightMP5[2] = {C_I, this->setup.h2};
    GT R_5;
    millerLoopVec(R_5, leftMP5, rightMP5, 2);
    finalExp(R_5, R_5);
    this->msg.R_5 = R_5;

    G1 base6[3] = {
        this->msg.P_4, this->setup.Frakg, this->setup.g1
    };
    Fr expon6[3] = {
        r_r_J, -r_delta_i[2], -r_r_A
    };
    G1::mulVec(first, base6, expon6, 3);
    G1 leftMP6[2] = {this->setup.Frakg * r_tau_i[2], first};
    G2 rightMP6[2] = {C_J, this->setup.h2};
    GT R_6;
    millerLoopVec(R_6, leftMP6, rightMP6, 2);
    finalExp(R_6, R_6);        
    this->msg.R_6 = R_6;

    G1 leftMP7[2] = {this->setup.g1 * r_r_I, this->setup.Frakg * (-r_tau_i[2])};
    G2 rightMP7[2] = {this->setup.h2, this->setup.g2};
    GT R_7;
    millerLoopVec(R_7, leftMP7, rightMP7, 2);
    finalExp(R_7, R_7);        
    this->msg.R_7 = R_7;

    string buf=this->msg.P_1.getStr()+this->msg.P_2.getStr()+this->msg.P_3.getStr()+this->msg.P_4.getStr()+this->msg.R_1.getStr()+this->msg.R_2.getStr()
    +this->msg.R_3.getStr()+this->msg.R_4.getStr()+this->msg.R_5.getStr()+this->msg.R_6.getStr()+this->msg.R_7.getStr();
    Fr c;
    c.setHashOf(buf);


    this->response.s_r_I=r_r_I+c*r_I;
    this->response.s_r_J=r_r_J+c*r_J;
    this->response.s_r_A=r_r_A+c*r_A;
    for(int i=0;i<4;i++){

        this->response.s_tau_i[i]=r_tau_i[i]+c*tau_i[i];
        this->response.s_delta_i[i]=r_delta_i[i]+c*delta_i[i];

    }
}


bool ZKSP_verify(ZKSP *pi) {

    string buf=pi->msg.P_1.getStr()+pi->msg.P_2.getStr()+pi->msg.P_3.getStr()+pi->msg.P_4.getStr()+pi->msg.R_1.getStr()+pi->msg.R_2.getStr()
        +pi->msg.R_3.getStr()+pi->msg.R_4.getStr()+pi->msg.R_5.getStr()+pi->msg.R_6.getStr()+pi->msg.R_7.getStr();
    Fr c;
    c.setHashOf(buf);

    // Prepare MSM bases/exponents for R1 -- R4
    G1 base1[3] = {pi->msg.P_1, pi->setup.g1, pi->setup.Frakg};
    G1 base2[3] = {pi->msg.P_3, pi->setup.g1, pi->setup.Frakg};

    Fr expon1[3] = {-c, pi->response.s_tau_i[0], pi->response.s_tau_i[1]};
    Fr expon2[3] = {pi->response.s_r_I, -pi->response.s_delta_i[0], -pi->response.s_delta_i[1]};
    Fr expon3[3] = {-c, pi->response.s_tau_i[2], pi->response.s_tau_i[3]};
    Fr expon4[3] = {pi->response.s_r_J, -pi->response.s_delta_i[2], -pi->response.s_delta_i[3]};

    // Compute R1 -- R4
    G1 R_1v, R_2v, R_3v, R_4v;
    G1::mulVec(R_1v, base1, expon1, 3);
    G1::mulVec(R_2v, base1, expon2, 3);
    G1::mulVec(R_3v, base2, expon3, 3);
    G1::mulVec(R_4v, base2, expon4, 3);

    // Compute R5 -- R7
    G1 base5[3] = {
        pi->msg.P_2, pi->setup.Frakg, pi->setup.g1
    };
    Fr expon5[3] = {
        pi->response.s_r_I, -pi->response.s_delta_i[0], -pi->response.s_r_A
    };
    G1 first5;
    G1::mulVec(first5, base5, expon5, 3);

    G1 leftMP5[3] = {
        pi->setup.Frakg * pi->response.s_tau_i[0] + pi->msg.P_2 * (-c),
        pi->setup.g1 * c,
        first5
    };
    G2 rightMP5[3] = {
        pi->C_I,
        pi->C_A,
        pi->setup.h2
    };
    GT R_5v;
    millerLoopVec(R_5v, leftMP5, rightMP5, 3);
    finalExp(R_5v, R_5v);    

    G1 base6[3] = {
        pi->msg.P_4, pi->setup.Frakg, pi->setup.g1
    };
    Fr expon6[3] = {
        pi->response.s_r_J, -pi->response.s_delta_i[2], -pi->response.s_r_A
    };
    G1 first6;
    G1::mulVec(first6, base6, expon6, 3);

    G1 leftMP6[3] = {
        pi->setup.Frakg * pi->response.s_tau_i[2] + pi->msg.P_4 * (-c),
        pi->setup.g1 * c,
        first6
    };

    G2 rightMP6[3] = {
        pi->C_J,
        pi->C_A,
        pi->setup.h2
    };
    GT R_6v;
    millerLoopVec(R_6v, leftMP6, rightMP6, 3);
    finalExp(R_6v, R_6v);


    G1 leftMP7[2] = {
        pi->setup.g1,
        pi->msg.P_4 * c + pi->setup.Frakg * (-pi->response.s_tau_i[2])
    };

    G2 rightMP7[2] = {
        pi->setup.h2 * pi->response.s_r_I + pi->C_I * (-c),        
        pi->setup.g2si[0]
    };
    GT R_7v;
    millerLoopVec(R_7v, leftMP7, rightMP7, 2);
    finalExp(R_7v, R_7v);

    // Check that all R's are correctly computed.
    bool flag=true;
    flag = flag && (pi->msg.R_1 == R_1v);
    flag = flag && (pi->msg.R_2 == R_2v);
    flag = flag && (pi->msg.R_3 == R_3v);
    flag = flag && (pi->msg.R_4 == R_4v);
    flag = flag && (pi->msg.R_5 == R_5v);
    flag = flag && (pi->msg.R_6 == R_6v);
    flag = flag && (pi->msg.R_7 == R_7v);

    return flag;
}


// Proof of Knowledge
// For g^x * h^r
void PoK2::init(ZKAccSetup &setup) {
    this->setup = setup;
}

void PoK2::prove(G2 C, Fr id, Fr r) {
    this->C = C;
    Fr tau1, tau2;
    tau1.setByCSPRNG(); tau2.setByCSPRNG();
    R = setup.g2si[0] * tau1 + setup.h2 * tau2;
    string buf = C.getStr() + R.getStr();
    Fr c;
    c.setHashOf(buf);
    z1 = tau1 + id * c;
    z2 = tau2 + r * c;
}


bool PoK2_verify(PoK2 *pi) {
    string buf = pi->C.getStr() + pi->R.getStr();
    Fr c;
    c.setHashOf(buf);
    G2 left = pi->R + pi->C * c;
    G2 right = pi->setup.g2 * pi->z1 + pi->setup.h2 * pi->z2;
    return (left == right);
}


// Proof of Knowledge for n elements
// It is only used for Tx Summary
// So we apply multithreading here.
void PoKn::prove(PCS &pcs, G2 P, FrVec x, Fr r) {
    uint32_t n = x.size();
    FrVec taus; taus.resize(n+1);
    for (uint32_t i = 0; i < n+1; i++) {
        taus[i].setByCSPRNG();
    }
    G2 R;
    G2::mulVecMT(R, &pcs.pp.g2si[0], &taus[0], n);
    R += pcs.pp.h2si[0] * taus[n];

    string buf = P.getStr() + R.getStr();  
    Fr c; c.setHashOf(buf);

    FrVec zs; zs.resize(n + 1);
    #pragma omp parallel for
    for (uint32_t i = 0; i < n; i++) {
        zs[i] = taus[i] + c * x[i];
    }
    zs[n] = taus[n] + c * r;

    this->P = P;
    this->R = R;
    this->zs = zs;
}

bool PoKn_verify(PCS &pcs, PoKn *pi) {
    string buf = pi->P.getStr() + pi->R.getStr();
    Fr c; c.setHashOf(buf);
    uint32_t n = pi->zs.size();
    G2 left = pi->R + pi->P * c;
    G2 right;
    G2::mulVec(right, &pcs.pp.g2si[0], &pi->zs[0], n-1);
    right += pcs.pp.h2si[0] * pi->zs[n-1];
    return (left == right);
}

