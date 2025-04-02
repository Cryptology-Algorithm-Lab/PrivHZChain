#include <mcl/bls12_381.hpp>
#include <iostream>
#include <limits>
#include <cmath>
#include <vector>
#include <malloc.h>
#include "primitiveroots.hpp"
#include "ntt.hpp"
#include "polyonfr.hpp"

using namespace mcl;
using namespace mcl::bn;
using namespace std;
struct pcs;



template <typename T>
void putG(vector<T> a){
    for(uint32_t i=0;i<a.size();i++){
        cout<<"<"<<a[i].getStr(10).c_str()<<">"<<endl;
    }
}

GT pairing_eq(G1 a, G2 b){

    GT res;
    pairing(res, a, b);
    return res;

}

struct IDset{
    Fr* idset;
    FrVec idVec;
    void init(uint32_t n)
    {
        FrVec ids;
        ids.resize(n,0);
   
        for(int i=0;i<n;i++){
            Fr s;
            s.setByCSPRNG();
            ids[i]=s;
        }
        idVec = ids;
        idset = &ids[0];
    }
    
};

struct PCS{

    
    struct pp{
        vector<G1> g1si;
        vector<G2> g2si;
        vector<G1> h1si;
        vector<G2> h2si;
    } pp;
    Fr s;

    void setup(uint32_t D){
        Fr s;
        s.setByCSPRNG();
        this->s = s;

        G1 g1;
        G2 g2;
        G1 h1;
        G2 h2;
        hashAndMapToG1(g1,"g1");
        hashAndMapToG2(g2,"g2");
        hashAndMapToG1(h1,"h1");
        hashAndMapToG2(h2,"h2");

        vector<G1> g1si(D+1);
        vector<G2> g2si(D+1);
        vector<G1> h1si(D+1);
        vector<G2> h2si(D+1);

        g1si[0]=g1;
        g2si[0]=g2;
        h1si[0]=h1;
        h2si[0]=h2;
    
        for(int i=1;i<=D;i++){
        
            g1si[i]=g1si[i-1]*s;
            g2si[i]=g2si[i-1]*s;
            h1si[i]=h1si[i-1]*s;
            h2si[i]=h2si[i-1]*s;

        }

        this->pp.g1si=g1si;
        this->pp.g2si=g2si;
        this->pp.h1si=h1si;
        this->pp.h2si=h2si;

    }

    G2 commit(FrVec p){

        G2 tmp;
        G2 res;
        res.clear();
        G2::mulVec(res, &pp.g2si[0], &p[0], p.size());

        return res;
    }
    G2 commit(Fr* p, uint32_t n){

        G2 tmp;
        G2 res;
        res.clear();
        G2::mulVec(res, &pp.g2si[0], p, n);

        return res;
    }

    G1 commit_G1(FrVec p){

        G1 tmp;
        G1 res;
        res.clear();
        G1::mulVec(res, &pp.g1si[0], &p[0], p.size());
        return res;
    }

    G1 commit_h1(FrVec p){

        G1 tmp;
        G1 res;
        res.clear();
        G1::mulVec(res, &pp.h1si[0], &p[0], p.size());
        return res;
    }

    G2 commit_h2(FrVec p){

        G2 tmp;
        G2 res;
        res.clear();
        G2::mulVec(res, &pp.h2si[0], &p[0], p.size());
        return res;
    }



    G2 Eval(G2 C, Fr z, Fr y, FrVec p){
        
        FrVec pz ={z};
        FrVec pmz=PolySub(p,pz);
        FrVec xmy={-y,1};
        FrVec w=PolyLongDiv(pmz,xmy);

        G2 pi_open=commit(w);  

        return pi_open;      
    }

    G1 Eval(G1 C, Fr z, Fr y, FrVec p){
        
        FrVec pz ={z};
        FrVec pmz=PolySub(p,pz);
        FrVec xmy={-y,1};
        FrVec w=PolyLongDiv(pmz,xmy);

        G1 pi_open=commit_G1(w);  

        return pi_open;      
    }



    bool verify(G2 C, Fr z, Fr y, G2 pi_open){
        bool flag=false;
        GT res1;
        GT res2;

        pairing(res1, this->pp.g1si[0],C+(this->pp.g2si[0]*(-z)));
        pairing(res2, (this->pp.g1si[1])+(this->pp.g1si[0]*(-y)),pi_open);

        if(res1==res2){
            flag = true;
        }

        return flag;       
    }

    bool verify(G1 C, Fr z, Fr y, G1 pi_open){
        bool flag=false;
        GT res1;
        GT res2;

        pairing(res1, C+(this->pp.g1si[0]*(-z)),this->pp.g2si[0]);
        pairing(res2, pi_open, (this->pp.g2si[1])+(this->pp.g2si[0]*(-y)));

        if(res1==res2){
            flag = true;
        }

        return flag;       
    }
};

struct PCS_h{

    
    struct pp{
        vector<G1> g1si;
        vector<G2> g2si;
    } pp;

    void setup(uint32_t D){
        Fr s;
        s.setByCSPRNG();

        G1 g1;
        G2 g2;
        hashAndMapToG1(g1,"h1");
        hashAndMapToG2(g2,"h2");

        vector<G1> g1si(D+1);
        vector<G2> g2si(D+1);

        g1si[0]=g1;
        g2si[0]=g2;
    
        for(int i=1;i<=D;i++){
        
            g1si[i]=g1si[i-1]*s;
            g2si[i]=g2si[i-1]*s;

        }

        this->pp.g1si=g1si;
        this->pp.g2si=g2si;

    }

    G2 commit(FrVec p){

        G2 tmp;
        G2 res;
        res.clear();
        G2::mulVec(res, &pp.g2si[0], &p[0], p.size());

        return res;
    }
    G2 commit(Fr* p, uint32_t n){

        G2 tmp;
        G2 res;
        res.clear();
        G2::mulVec(res, &pp.g2si[0], p, n);

        return res;
    }

    G1 commit_G1(FrVec p){

        G1 tmp;
        G1 res;
        res.clear();
        G1::mulVec(res, &pp.g1si[0], &p[0], p.size());
        return res;
    }

    G2 Eval(G2 C, Fr z, Fr y, FrVec p){
        
        FrVec pz ={z};
        FrVec pmz=PolySub(p,pz);
        FrVec xmy={-y,1};
        FrVec w=PolyLongDiv(pmz,xmy);

        G2 pi_open=commit(w);  

        return pi_open;      
    }

    G1 Eval(G1 C, Fr z, Fr y, FrVec p){
        
        FrVec pz ={z};
        FrVec pmz=PolySub(p,pz);
        FrVec xmy={-y,1};
        FrVec w=PolyLongDiv(pmz,xmy);

        G1 pi_open=commit_G1(w);  

        return pi_open;      
    }



    bool verify(G2 C, Fr z, Fr y, G2 pi_open){
        bool flag=false;
        GT res1;
        GT res2;

        pairing(res1, this->pp.g1si[0],C+(this->pp.g2si[0]*(-z)));
        pairing(res2, (this->pp.g1si[1])+(this->pp.g1si[0]*(-y)),pi_open);

        if(res1==res2){
            flag = true;
        }

        return flag;       
    }

    bool verify(G1 C, Fr z, Fr y, G1 pi_open){
        bool flag=false;
        GT res1;
        GT res2;

        pairing(res1, C+(this->pp.g1si[0]*(-z)),this->pp.g2si[0]);
        pairing(res2, pi_open, (this->pp.g2si[1])+(this->pp.g2si[0]*(-y)));

        if(res1==res2){
            flag = true;
        }

        return flag;       
    }
};


struct zkbpacc_setup {
    PCS pcs;
    G1* g1si;
    G2* g2si;
    G1 g1;
    G2 g2;
    G1  Frakg;
    G1  Frakh;
    G2  h2;

    void init(uint32_t n)
    {   
        G1 g1;
        G2 g2;

    /*Getting generator*/
        hashAndMapToG1(g1,"g1");
        hashAndMapToG2(g2,"g2");
        hashAndMapToG1(Frakg,"Frakg");
        hashAndMapToG1(Frakh,"Frakh");
        hashAndMapToG2(h2,"h2");


        Fr s;
        s.setByCSPRNG();

        g1si=(G1*) malloc(sizeof(G1)*(n+1));
        g2si=(G2*) malloc(sizeof(G2)*(n+1));

        g1si[0]=g1;
        g2si[0]=g2;  
    
        for(int i=1;i<=n;i++){
        
            g1si[i]=g1si[i-1]*s;
            g2si[i]=g2si[i-1]*s;

        }
    }

    void setup_from_pcs(PCS& pcs) {
        this->pcs = pcs;
        g1si = &pcs.pp.g1si[0];
        g2si = &pcs.pp.g2si[0];
        g1 = pcs.pp.g1si[0];
        g2 = pcs.pp.g2si[0];        
        hashAndMapToG1(Frakg,"Frakg");
        hashAndMapToG1(Frakh,"Frakh");
        hashAndMapToG2(h2,"h2");
    }
};

struct zkbpacc_setup_pcs {
    
    PCS& pcs;
    G1  Frakg;
    G1  Frakh;
    G2  h2;

    zkbpacc_setup_pcs(PCS& _pcs) : pcs(_pcs) {
        init();
    }

    void init()
    {   
        hashAndMapToG1(Frakg,"Frakg");
        hashAndMapToG1(Frakh,"Frakh");
        hashAndMapToG2(h2,"h2");
    }

};

template <typename T>
T* Vec2arr(vector<T> a){

    uint32_t len=a.size();
    T* res=(T*) malloc((sizeof(T))*len);

    for(uint32_t i=0;i<len;i++){
        res[i]=a[i];
    }

    return res;
}

template <typename T>
vector<T> arr2Vec(T* a, uint32_t n){

    vector<T> dest(n);
    memcpy(&dest[0],&a[0],n*sizeof(T));

    return dest;
}

template <typename T>
vector<T> convG(vector<T> a, vector<T> b){

    if(a.size()!=b.size()){
        throw runtime_error("Two vectors must have same length");
    }

    uint32_t n = a.size();
    vector<T> res;
    for(uint32_t i=0;i<n;i++){

        T element = a[i]+b[i];
        res.insert(res.begin(),element);

    }
    return res;
    
}

/*
a in GG1^n, b in GG2
return \prod_{i=0}^{n-1} e(a[i],b[i]);
*/
GT MultiPairing(vector<G1> a, vector<G2> b){
    if(a.size()!=b.size()){
        throw runtime_error("Two vectors must have same length");
    }

    uint32_t n=a.size();
    G1* a_arr=&a[0]; G2* b_arr=&b[0];
    GT c;
    millerLoopVec(c,a_arr,b_arr,n);
    finalExp(c,c);

    return c;
}

struct IPPproof {

        vector<G2> gg;
        GT P;
        GT Q;
        vector<GT> L;
        vector<GT> R;
        vector<Fr> x;
        G1 w;
        G2 g;

        IPPproof(vector<G2> gg_0, GT P){
                this->gg=gg_0;
                this->P=P;
                this->L.clear();
                this->R.clear();
                this->x.clear();
        }

        void Prove(vector<G2> gg, GT P, vector<G1> ww){
            if((IsPowerOfTwo(ww.size())==false) || (IsPowerOfTwo(gg.size())==false)){
                throw runtime_error("gg and ww need to be power of two");
            }

            uint32_t n = ww.size();
            GT l, r, Lx, Rinvx, P_hat;
            Fr c; string buf;

            while (n > 1) {
                n>>=1;
                millerLoopVec(l, &ww[n], &gg[0], n); millerLoopVec(r, &ww[0], &gg[n], n);
                finalExp(l, l); finalExp(r, r);
                L.push_back(l); R.push_back(r);

                buf = this->P.getStr()+l.getStr()+r.getStr();
                c.setHashOf(buf);
                x.push_back(c);

                for (uint32_t i = 0; i < n; i++) {
                    ww[i]=ww[i]+(ww[i + n]*(c));
                    gg[i]=gg[i]+(gg[i + n]*(1/c));                
                }
            }

            this->w = ww[0];
            this->g = gg[0];

        }
};

bool IPPverify(IPPproof pi){
    
    uint32_t lgn = pi.L.size();
    FrVec x_vec;

    // Compute Left
    FrVec expon;
    vector<GT> base;
    GT P_v;
    
    for(uint32_t i=0;i<lgn;i++){
        string buf = pi.P.getStr()+pi.L[i].getStr()+pi.R[i].getStr();
        Fr c;
        c.setHashOf(buf);
        expon.push_back(c);
        expon.push_back(1/c);
        x_vec.push_back(1/c);
        base.push_back(pi.L[i]);
        base.push_back(pi.R[i]);
    }
        
    GT::powVec(P_v, &base[0], &expon[0], 2*lgn);
    P_v = P_v * pi.P;
        
    // Compute Right
    Fr one = Fr(1);
    
    uint32_t n = pi.gg.size();
    expon.clear();
    expon.push_back(one);
       
    // Prepare Indices 
    for (uint32_t i = 1; i < n; i++) {
        uint32_t lgi = (int)(log2(i));
        uint32_t k = 1<<lgi;
        Fr u = x_vec[lgn - lgi - 1];                
        expon.push_back(expon[i-k] * u);        
    }
    
    // Do MSM & Pairing
    G2 g_v;
    G2::mulVec(g_v, &pi.gg[0], &expon[0], n);
    GT res;
    pairing(res, pi.w, g_v);

    // Check Left == Right
    bool flag=(P_v==res);
    return flag;
}

struct zkIPPproof{
    
    IPPproof pi;
    GT Q;
    GT R;

    void Prove(vector<G2> gg, GT P, vector<G1> ww){
        
        vector<G1> vv(gg.size());
        for(uint32_t i=0;i<ww.size();i++){

            Fr c;
            c.setByCSPRNG();
            G1 v_i;
            hashAndMapToG1(v_i,c.getStr());
            vv[i]=v_i;

        }

        millerLoopVec(Q, &vv[0], &gg[0], gg.size());
        finalExp(Q, Q);
        Fr chal;
        chal.setHashOf(Q.getStr());
        vector<G1> uu(ww.size());
        
        for(uint32_t i=0;i<ww.size();i++){
            uu[i]=(ww[i]*chal)+vv[i];
        }

        GT P_c;
        GT::pow(P_c,P,chal);
        this->R=P_c*(this->Q);
        
        IPPproof pi_zk(gg,P);
        pi_zk.Prove(gg,R,uu);
        this->pi=pi_zk;
        pi_zk.P=P;

    }

};

IPPproof zkIPPprove(vector<G2> gg, GT P, vector<G1> ww)
{
    vector<G1> vv(gg.size());
        for(uint32_t i=0;i<ww.size();i++){

            Fr c;
            c.setByCSPRNG();
            G1 v_i;
            hashAndMapToG1(v_i,c.getStr());
            vv[i]=v_i;

        }
        GT Q;
        millerLoopVec(Q, &vv[0], &gg[0], gg.size());
        finalExp(Q, Q);
        Fr chal;
        chal.setHashOf(Q.getStr());
        vector<G1> uu(ww.size());
        
        for(uint32_t i=0;i<ww.size();i++){
            uu[i]=(ww[i]*chal)+vv[i];
        }
        
        GT P_c;
        GT::pow(P_c,P,chal);
        GT R=P_c*Q;
        IPPproof pi_zk(gg,R);
        pi_zk.Q=Q;
        pi_zk.Prove(gg,R,uu);
        pi_zk.P=P;

        return pi_zk;
}

bool zkIPPverify(IPPproof pi){

    Fr chal;
    chal.setHashOf(pi.Q.getStr());
    GT P_c;
    GT::pow(P_c,pi.P,chal);
    GT R=P_c*pi.Q;
    pi.P=R;

    bool flag=IPPverify(pi);
    return flag;
}


struct ZKMP {

    zkbpacc_setup setup;
    G2 C_I;
    G2 C_A;

    ZKMP(zkbpacc_setup setup){
        this->setup=setup;
    }

    void init(zkbpacc_setup setup) {
        this->setup = setup;
    }

    struct msg{

        G1 P_1;
        G1 P_2;
        G1 R_1;
        G1 R_2;
        GT R_3;

    }msg;

    struct response{

        Fr s_r_I;
        Fr s_r_A;
        Fr s_tau_1;
        Fr s_tau_2;
        Fr s_delta_1;
        Fr s_delta_2;

    }response;

    void prove(G2 C_I, G2 C_A, G1 pi_I, Fr r_I, Fr r_A){

        this->C_I=C_I;
        this->C_A=C_A;

        Fr tau_1;
        Fr tau_2;
        Fr r_r_I;
        Fr r_r_A;
        Fr r_tau_1;
        Fr r_tau_2;
        Fr r_delta_1;
        Fr r_delta_2;

        tau_1.setByCSPRNG();
        tau_2.setByCSPRNG();
        r_r_I.setByCSPRNG();
        r_r_A.setByCSPRNG();
        r_tau_1.setByCSPRNG();
        r_tau_2.setByCSPRNG();
        r_delta_1.setByCSPRNG();
        r_delta_2.setByCSPRNG();

        Fr delta_1 = r_I*tau_1;
        Fr delta_2 = r_I*tau_2;

        this->msg.P_1 = ((this->setup.g1)*tau_1)+((this->setup.Frakg)*tau_2);
        this->msg.P_2 = pi_I+((this->setup.Frakg)*tau_1);
        this->msg.R_1 = ((this->setup.g1)*r_tau_1)+((this->setup.Frakg)*r_tau_2);
        this->msg.R_2 = ((this->msg.P_1)*r_r_I)+((this->setup.g1)*(-r_delta_1))+((this->setup.Frakg)*(-r_delta_2));
        
        // GT R_3_1;
        // GT R_3_2;
        // GT R_3_3;
        // GT R_3_4;

        // pairing(R_3_1, this->setup.Frakg*r_tau_1, C_I);
        // pairing(R_3_2, this->msg.P_2*r_r_I, this->setup.h2);
        // pairing(R_3_3, this->setup.Frakg*r_delta_1, this->setup.h2);

        G2 first = C_I * r_tau_1 + this->setup.h2 * (-r_delta_1);
        G1 second = this->msg.P_2 * r_r_I + this->setup.g1 * (-r_r_A);
        G1 leftMP[2] = {this->setup.Frakg, second};
        G2 rightMP[2] = {first, this->setup.h2};
        GT R3;
        millerLoopVec(R3, leftMP, rightMP, 2);
        finalExp(R3, R3);
        this->msg.R_3=R3;

        string buf=this->msg.P_1.getStr()+this->msg.P_2.getStr()+this->msg.R_1.getStr()+this->msg.R_2.getStr()+this->msg.R_3.getStr();
        Fr c;
        c.setHashOf(buf);


        this->response.s_r_I=r_r_I+c*r_I;
        this->response.s_r_A=r_r_A+c*r_A;
        this->response.s_tau_1=r_tau_1+c*tau_1;
        this->response.s_tau_2=r_tau_2+c*tau_2;
        this->response.s_delta_1=r_delta_1+c*delta_1;
        this->response.s_delta_2=r_delta_2+c*delta_2;

    }

};

bool ZKMP_verify(ZKMP pi){

    string buf=pi.msg.P_1.getStr()+pi.msg.P_2.getStr()+pi.msg.R_1.getStr()+pi.msg.R_2.getStr()+pi.msg.R_3.getStr();
    Fr c;
    c.setHashOf(buf);

    G1 base[3] = {pi.msg.P_1, pi.setup.g1, pi.setup.Frakg};

    // Check R1
    G1 right1;
    Fr expon1[3] = {-c, pi.response.s_tau_1, pi.response.s_tau_2};
    G1::mulVec(right1, base, expon1, 3);

    // Check R2
    G1 right2;
    Fr expon2[3] = {pi.response.s_r_I, -pi.response.s_delta_1, -pi.response.s_delta_2};
    G1::mulVec(right2, base, expon2, 3);

    // Check R3
    G1 leftMP[3] = {pi.setup.Frakg, pi.msg.P_2, pi.setup.g1};
    G2 firstTerm = pi.C_I * pi.response.s_tau_1 + pi.setup.h2 * (-pi.response.s_delta_1);
    G2 secondTerm = pi.setup.h2 * pi.response.s_r_I + pi.C_I * (-c);
    G2 thirdTerm = pi.setup.h2 * (-pi.response.s_r_A) + pi.C_A * c;
    G2 rightMP[3] = {firstTerm, secondTerm, thirdTerm};

    GT right3;
    millerLoopVec(right3, leftMP, rightMP, 3);
    finalExp(right3, right3);

    bool flag1 = (pi.msg.R_1 == right1);
    bool flag2 = (pi.msg.R_2 == right2);
    bool flag3 = (pi.msg.R_3 == right3);

    return flag1 && flag2 && flag3;
};

struct ZKNMP {
    zkbpacc_setup setup;
    G2 C_I;
    G2 C_A;

    ZKNMP(zkbpacc_setup setup){
        this->setup=setup;
    }

    void init(zkbpacc_setup setup) {
        this->setup = setup;
    }

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

    }msg;

    struct response{

        Fr s_r_I;
        Fr s_r_A;
        Fr s_tau_i[4];
        Fr s_delta_i[4];

    }response;

    void prove(G2 C_I, G2 C_A, G1 w_I, G1 w_A, Fr r_I, Fr r_A){
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

        // 
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
};

bool ZKNMP_verify(ZKNMP pi){
    string buf=pi.msg.P_1.getStr()+pi.msg.P_2.getStr()+pi.msg.Q_1.getStr()+pi.msg.Q_2.getStr()+pi.msg.R_1.getStr()+pi.msg.R_2.getStr()
    +pi.msg.R_3.getStr()+pi.msg.R_4.getStr()+pi.msg.R_5.getStr();
    Fr c;
    c.setHashOf(buf);

    G1 base1[3] = {pi.msg.P_1, pi.setup.g1, pi.setup.Frakg};
    G1 base2[3] = {pi.msg.Q_1, pi.setup.g1, pi.setup.Frakh};

    Fr expon1[3] = {-c, pi.response.s_tau_i[0], pi.response.s_tau_i[1]};
    Fr expon2[3] = {pi.response.s_r_I, -pi.response.s_delta_i[0], -pi.response.s_delta_i[1]};
    Fr expon3[3] = {-c, pi.response.s_tau_i[2], pi.response.s_tau_i[3]};
    Fr expon4[3] = {pi.response.s_r_A, -pi.response.s_delta_i[2], -pi.response.s_delta_i[3]};

    G1 R_1v, R_2v, R_3v, R_4v;
    G1::mulVec(R_1v, base1, expon1, 3);
    G1::mulVec(R_2v, base1, expon2, 3);
    G1::mulVec(R_3v, base2, expon3, 3);
    G1::mulVec(R_4v, base2, expon4, 3);

    GT R_5v;
    G1 first;
    G1 base[4] = {
        pi.msg.P_2, 
        pi.msg.Q_2,
        pi.setup.Frakg, 
        pi.setup.Frakh
    };
    Fr expon[4] = {
        pi.response.s_r_I,
        pi.response.s_r_A,
        -pi.response.s_delta_i[0],
        -pi.response.s_delta_i[2],
    };
    G1::mulVec(first, base, expon, 4);

    G1 leftMP[4] = {
        pi.setup.g1 * (c),
        pi.setup.Frakg * pi.response.s_tau_i[0] + pi.msg.P_2 * (-c),
        pi.setup.Frakh * pi.response.s_tau_i[2] + pi.msg.Q_2 * (-c),
        first,
    };
    G2 rightMP[4] = {
        pi.setup.g2,
        pi.C_I,
        pi.C_A,
        pi.setup.h2
    };

    millerLoopVec(R_5v, leftMP, rightMP, 4);
    finalExp(R_5v, R_5v);

    bool flag1 = (pi.msg.R_1 == R_1v);
    bool flag2 = (pi.msg.R_2 == R_2v);
    bool flag3 = (pi.msg.R_3 == R_3v);
    bool flag4 = (pi.msg.R_4 == R_4v);
    bool flag5 = (pi.msg.R_5 == R_5v);
    
    return flag1 && flag2 && flag3 && flag4 && flag5;
}

struct ZKSP {

    zkbpacc_setup setup;
    G2 C_I;
    G2 C_J;
    G2 C_A;

    ZKSP(zkbpacc_setup setup){
        this->setup=setup;
    }

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

    void prove(G2 C_I, G2 C_J, G2 C_A, G1 pi_I, G1 pi_J, Fr r_I, Fr r_J, Fr r_A){

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

};

bool ZKSP_verify(ZKSP pi){

    string buf=pi.msg.P_1.getStr()+pi.msg.P_2.getStr()+pi.msg.P_3.getStr()+pi.msg.P_4.getStr()+pi.msg.R_1.getStr()+pi.msg.R_2.getStr()
        +pi.msg.R_3.getStr()+pi.msg.R_4.getStr()+pi.msg.R_5.getStr()+pi.msg.R_6.getStr()+pi.msg.R_7.getStr();
    Fr c;
    c.setHashOf(buf);

    // Prepare MSM bases/exponents for R1 -- R4
    G1 base1[3] = {pi.msg.P_1, pi.setup.g1, pi.setup.Frakg};
    G1 base2[3] = {pi.msg.P_3, pi.setup.g1, pi.setup.Frakg};

    Fr expon1[3] = {-c, pi.response.s_tau_i[0], pi.response.s_tau_i[1]};
    Fr expon2[3] = {pi.response.s_r_I, -pi.response.s_delta_i[0], -pi.response.s_delta_i[1]};
    Fr expon3[3] = {-c, pi.response.s_tau_i[2], pi.response.s_tau_i[3]};
    Fr expon4[3] = {pi.response.s_r_J, -pi.response.s_delta_i[2], -pi.response.s_delta_i[3]};

    // Compute R1 -- R4
    G1 R_1v, R_2v, R_3v, R_4v;
    G1::mulVec(R_1v, base1, expon1, 3);
    G1::mulVec(R_2v, base1, expon2, 3);
    G1::mulVec(R_3v, base2, expon3, 3);
    G1::mulVec(R_4v, base2, expon4, 3);

    // Compute R5 -- R7
    G1 base5[3] = {
        pi.msg.P_2, pi.setup.Frakg, pi.setup.g1
    };
    Fr expon5[3] = {
        pi.response.s_r_I, -pi.response.s_delta_i[0], -pi.response.s_r_A
    };
    G1 first5;
    G1::mulVec(first5, base5, expon5, 3);

    G1 leftMP5[3] = {
        pi.setup.Frakg * pi.response.s_tau_i[0] + pi.msg.P_2 * (-c),
        pi.setup.g1 * c,
        first5
    };
    G2 rightMP5[3] = {
        pi.C_I,
        pi.C_A,
        pi.setup.h2
    };
    GT R_5v;
    millerLoopVec(R_5v, leftMP5, rightMP5, 3);
    finalExp(R_5v, R_5v);    

    G1 base6[3] = {
        pi.msg.P_4, pi.setup.Frakg, pi.setup.g1
    };
    Fr expon6[3] = {
        pi.response.s_r_J, -pi.response.s_delta_i[2], -pi.response.s_r_A
    };
    G1 first6;
    G1::mulVec(first6, base6, expon6, 3);

    G1 leftMP6[3] = {
        pi.setup.Frakg * pi.response.s_tau_i[2] + pi.msg.P_4 * (-c),
        pi.setup.g1 * c,
        first6
    };

    G2 rightMP6[3] = {
        pi.C_J,
        pi.C_A,
        pi.setup.h2
    };
    GT R_6v;
    millerLoopVec(R_6v, leftMP6, rightMP6, 3);
    finalExp(R_6v, R_6v);


    G1 leftMP7[2] = {
        pi.setup.g1,
        pi.msg.P_4 * c + pi.setup.Frakg * (-pi.response.s_tau_i[2])
    };

    G2 rightMP7[2] = {
        pi.setup.h2 * pi.response.s_r_I + pi.C_I * (-c),        
        pi.setup.g2si[0]
    };
    GT R_7v;
    millerLoopVec(R_7v, leftMP7, rightMP7, 2);
    finalExp(R_7v, R_7v);

    // Check that all R's are correctly computed.
    bool flag=true;
    flag = flag && (pi.msg.R_1 == R_1v);
    flag = flag && (pi.msg.R_2 == R_2v);
    flag = flag && (pi.msg.R_3 == R_3v);
    flag = flag && (pi.msg.R_4 == R_4v);
    flag = flag && (pi.msg.R_5 == R_5v);
    flag = flag && (pi.msg.R_6 == R_6v);
    flag = flag && (pi.msg.R_7 == R_7v);

    return flag;
}

struct PoK_proof_g2{

    G2 P;
    G2 R;
    G2 g2;
    Fr s;

    void prove(G2 P, Fr x){

        this->P = P;
        Fr tau;
        tau.setByCSPRNG();
        G2 g2;
        hashAndMapToG2(g2,"g2");
        this->g2 = g2;
        R = g2*tau;
        Fr c;
        c.setHashOf(R.getStr());
        this->s = tau + c*x;


    }

};

bool PoK_proof_g2_verify(PoK_proof_g2 pi) {
    Fr c; c.setHashOf(pi.R.getStr());
    G2 left = pi.g2 * pi.s;
    G2 right = pi.R + pi.P * c;
    return (left == right);
}

struct PoK_proof_g2_n{

    G2 P;
    G2 R;
    FrVec ss;
    vector<G2> gg;

    void prove(vector<G2> gg, G2 P, FrVec xx, uint32_t n){

        assert(gg.size()==n);
        this->P = P;

        FrVec tau; tau.resize(n);
        for(uint32_t i=0;i<n;i++){
            tau[i].setByCSPRNG();
        }

        G2::mulVec(this->R, &gg[0], &tau[0], n);
        Fr c; c.setHashOf(R.getStr());
        this->ss = PolyAdd(tau, PolyMul({c},xx));
        this->gg = gg;
    }

};

bool PoK_proof_g2_n_verify(PoK_proof_g2_n pi) {
    Fr c; c.setHashOf(pi.R.getStr());
    G2 left;
    G2::mulVec(left, &pi.gg[0], &pi.ss[0], pi.ss.size());
    G2 right = pi.R + pi.P * c;
    return (left == right);
}

struct PoK2_G2{
    zkbpacc_setup setup;
    G2 C; G2 R;
    Fr z1; Fr z2;

    PoK2_G2(zkbpacc_setup setup) {
        this->setup = setup;
    }

    void init(zkbpacc_setup setup) {
        this->setup = setup;
    }

    void prove(G2 C, Fr id, Fr r) {
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
};

bool PoK2_G2_verify (PoK2_G2 pi) {
    string buf = pi.C.getStr() + pi.R.getStr();
    Fr c;
    c.setHashOf(buf);
    G2 left = pi.R + pi.C * c;
    G2 right = pi.setup.g2 * pi.z1 + pi.setup.h2 * pi.z2;
    return (left == right);
}


// Polynomial Coefficient Calculation
// Complexity: O(nlogn^2; stable algorithm)
FrVec constructMemPoly(
    FrVec S
) {
    uint32_t numElts = S.size();
    if (numElts == 0) {
        return {1};
    }
    else if (numElts == 1) {
        return {S[0], 1};
    }
    FrVec left(S.begin(), S.begin() + numElts/2);
    FrVec right(S.begin() + numElts/2, S.end());

    FrVec leftRet = constructMemPoly(left);
    FrVec rightRet = constructMemPoly(right);
    FrVec ret = PolyMul(leftRet, rightRet);
    return ret;
}



// It constructs a vector that contains (I(X)/id_i+X) with the membership polynomial
// Fast but Memory Intensive
typedef struct _cachedAggAwResult {
    std::vector<FrVec> Payload;
    FrVec memPoly;
} cachedAggAwResult;

cachedAggAwResult cachedAggAWPolyGen(
    FrVec ID
) {
    uint32_t numElts = ID.size();
    if (numElts == 0) {
        // Do Nothing
        return cachedAggAwResult{
            {{1}}, {1}
        };
    }
    else if (numElts == 1) {
        return cachedAggAwResult{
            {{1}}, {ID[0], 1}
        };
    }
    FrVec left(ID.begin(), ID.begin() + numElts/2);
    FrVec right(ID.begin() + numElts/2, ID.end());

    cachedAggAwResult leftRet = cachedAggAWPolyGen(left);
    cachedAggAwResult rightRet = cachedAggAWPolyGen(right);

    // Merge
    std::vector<FrVec> newPayload;
    FrVec _tmp;

    for (FrVec val : leftRet.Payload) {
        _tmp = PolyMul(val, rightRet.memPoly);
        newPayload.push_back(_tmp);
    }

    for (FrVec val : rightRet.Payload) {
        _tmp = PolyMul(val, leftRet.memPoly);
        newPayload.push_back(_tmp);
    }

    FrVec newMemPoly = PolyMul(leftRet.memPoly, rightRet.memPoly);

    return cachedAggAwResult {
        newPayload, newMemPoly
    };
}


FrVec constructMemPolyBatch (
    FrVec M,
    FrVec N
) {
    // Find Set Difference for sorted vectors
    FrVec ret;
    std::sort(M.begin(), M.end());
    std::sort(N.begin(), N.end());
    std::set_difference(
        M.begin(), M.end(), N.begin(), N.end(), 
        std::back_inserter(ret)
    );
    return constructMemPoly(ret);
}

FrVec constructMemPolyBatchSorted (
    FrVec M,
    FrVec N
) {
    // Find Set Difference for sorted vectors
    FrVec ret;
    std::set_difference(
        M.begin(), M.end(), N.begin(), N.end(), 
        std::back_inserter(ret)
    );
    return constructMemPoly(ret);
}

G1 memWitGen(
    PCS pcs,
    FrVec S,
    Fr n
) {
    // Remove
    // S.erase(remove(S.begin(), S.end(), n), S.end());
    FrVec SPoly = constructMemPoly(S);
    FrVec nPoly = {n, 1};
    FrvT_2 divret = PolyDiv(SPoly, nPoly);
    FrVec retPoly = get<0>(divret);
    return pcs.commit_G1(retPoly);
}

tuple<G1, G1> nonMemWitGen(
    PCS pcs,
    FrVec S,
    Fr n    
) {
    // Make Polynomials
    FrVec SPoly = constructMemPoly(S);
    FrVec nPoly = {n, 1};

    FrvT_3 ret = xGCD(SPoly, nPoly);
    FrVec alpha = get<0>(ret);
    FrVec beta = get<1>(ret);
    FrVec gcd = get<2>(ret);

    if (gcd.size() > 1) {
        throw std::runtime_error("GCD is NOT correctly calculated");
    }
    G1 w_1 = pcs.commit_G1(beta);
    G1 w_2 = pcs.commit_G1(alpha);

    return tuple<G1, G1>(w_1, w_2);
} 

struct OwnPf {
    zkbpacc_setup setup;
    ZKMP piMP;
    PoK2_G2 piPoK;

    void init(zkbpacc_setup setup) {
        this->setup = setup;
        this->piMP.init(setup);
        this->piPoK.init(setup);
    }

    void prove(
        G2 C_id, G2 C_A, Fr n, 
        FrVec A, Fr id, Fr r, Fr delta
    ) {
        G2 C_pok = C_id - setup.g2si[1];

        // Witness Generation
        G1 w_n = memWitGen(setup.pcs, A, n);

        // Compute C_I
        FrVec nPoly = {n, 1};
        G2 C_I = setup.pcs.commit(nPoly);

        // Do Proof
        piPoK.prove(C_pok, id, r);
        piMP.prove(C_I, C_A, w_n, Fr(0), delta);
    }
};



bool OwnPf_verify(OwnPf pi) {
    bool flag1 = PoK2_G2_verify(pi.piPoK);
    bool flag2 = ZKMP_verify(pi.piMP);
    return flag1 && flag2;
}


struct NonOwnPf {
    zkbpacc_setup setup;
    ZKNMP piNMP;
    PoK2_G2 piPoK;

    void init(zkbpacc_setup setup) {
        this->setup = setup;
        this->piNMP.init(setup);
        this->piPoK.init(setup);
    }

    void prove(
        G2 C_id, G2 C_A, Fr n, 
        FrVec A, Fr id, Fr r, Fr delta
    ) {
        G2 C_pok = C_id - setup.g2si[1];

        // Witness Generation
        tuple<G1,G1> wits = nonMemWitGen(setup.pcs, A, n);
        G1 w_1 = get<0>(wits);
        G1 w_2 = get<1>(wits);

        // Compute C_I
        FrVec nPoly = {n, 1};
        G2 C_I = setup.pcs.commit(nPoly);


        piPoK.prove(C_pok, id, r);
        piNMP.prove(C_I, C_A, w_1, w_2, Fr(0), delta);
    }
};

bool NonOwnPf_verify(NonOwnPf pi) {
    bool flag1 = PoK2_G2_verify(pi.piPoK);
    bool flag2 = ZKNMP_verify(pi.piNMP);

    return flag1 && flag2;
}






// int main(){
    
//     initPairing(mcl::BLS12_381);


//     zkbpacc_setup setup1;
//     setup1.init(8);

//     FrVec poly1 = {6, 11, 6, 1};
//     FrVec poly2 = {1, 4};

//     FrvT_3 out=xGCD(poly1, poly2);

//     FrVec alpha = get<0>(out);
//     FrVec beta = get<1>(out);
//     FrVec gcd = get<2>(out);

//     PCS pcs1;
//     pcs1.setup(8);

//     G1 pi_1 = pcs1.commit_G1(alpha);
//     G1 pi_2 = pcs1.commit_G1(beta);

//     G2 C_A = pcs1.commit(poly1);
//     G2 C_I = pcs1.commit(poly2);
//     Fr r_I;
//     r_I.setByCSPRNG();
//     Fr r_A;
//     r_A.setByCSPRNG();
//     C_I = C_I+setup1.h2*r_I;
//     C_A = C_A+setup1.h2*r_A;

//     ZKNMP nmp(setup1);
//     nmp.prove(C_I, C_A, pi_1, pi_2, r_I, r_A);

//     bool flag = ZKNMP_verify(nmp);

//     std::cout<<"===================================="<<endl;

//     zkbpacc_setup setup2;
//     setup2.init(8);

//     FrVec Ax = {24, -50, 35, -10, 1};
//     FrVec Ix = {2, -3, 1};
//     FrVec Jx = {12, -7, 1};

//     PCS pcs2;
//     pcs2.setup(8);

//     C_A = pcs2.commit(Ax);
//     C_I = pcs2.commit(Ix);
//     G2 C_J = pcs2.commit(Jx);

//     G1 pi_I = pcs2.commit_G1(PolyLongDiv(Ax,Ix));
//     G1 pi_J = pcs2.commit_G1(PolyLongDiv(Ax,Jx));

//     r_A;
//     r_A.setByCSPRNG();
//     r_I;
//     r_I.setByCSPRNG();
//     Fr r_J;
//     r_J.setByCSPRNG();

//     C_A = C_A + setup2.h2*r_A;
//     C_I = C_I + setup2.h2*r_I;
//     C_J = C_J + setup2.h2*r_J;

//     ZKSP pi3(setup2);
//     pi3.prove(C_I, C_J, C_A, pi_I, pi_J, r_I, r_J, r_A);
//     ZKSP_verify(pi3);

//     cout<<"===================================="<<endl;

//     ZKMP pi4(setup2);
//     pi4.prove(C_I, C_A, pi_I, r_I, r_A);
//     bool fla = ZKMP_verify(pi4);
//     cout<<fla<<endl;

// }

