#include <mcl/bls12_381.hpp>
#include <iostream>
#include <limits>
#include <cmath>
#include <vector>
#include <malloc.h>
#include <mcl/ntt.hpp>
#include <chrono>
#include <fstream>
#include <iomanip> 

// Background Maths 
#include "primitiveroots.hpp"
#include "polyonfr.hpp"

// Main Logics
#include "mainalgorithms.hpp"
#include "TxGen.cpp"
#include "Aggtest.cpp"


using namespace mcl;
using namespace mcl::bn;
using namespace std;
using namespace chrono;

void ZKMP_test(){
    cout << "<<< ZKMP TEST START >>>" << endl;

    ofstream fout1;
    ofstream fout2;
    fout1.open("ZKMP_p.txt");
    fout2.open("ZKMP_v.txt");

    PCS pcs;
    pcs.setup(8);

    zkbpacc_setup setup;
    setup.setup_from_pcs(pcs);    

    FrVec A = {120, 274, 225, 85, 15, 1};
    FrVec I = {1,1};

    G2 C_A = pcs.commit(A);
    G2 C_I = pcs.commit(I);
    Fr r_I;
    r_I.setByCSPRNG();
    Fr r_A;
    r_A.setByCSPRNG();
    C_I = C_I+setup.h2*r_I;
    C_A = C_A+setup.h2*r_A;

    G1 pi_I = pcs.commit_G1(PolyLongDiv(A,I));

    ZKMP zkmp(setup);

    microseconds total_p(0);
    microseconds total_v(0);
    microseconds elaps;

    for(uint32_t i=0;i<1000;i++){

        auto start_1 = chrono::steady_clock::now();

        zkmp.prove(C_I, C_A, pi_I, r_I, r_A);

        auto end_1 = chrono::steady_clock::now();
        elaps = duration_cast<microseconds>(end_1-start_1);
        fout1<< elaps.count()<< ", ";
        total_p = total_p+elaps;

        auto start_2 = chrono::steady_clock::now();

        bool flag = ZKMP_verify(&zkmp);

        if (flag == 0) {
            throw std::runtime_error("Verification Failed on ZKMP_test");
        }

        auto end_2 = chrono::steady_clock::now();
        elaps = duration_cast<microseconds>(end_2-start_2);
        fout2<< elaps.count()<< ", ";
        total_v = total_v+elaps;

    }

    std::cout << "ZKMP_prove = : " << total_p.count()/1000 << "(μs)"<<std::endl;
    std::cout << "ZKMP_verifiy = : " << total_v.count()/1000 << "(μs)"<<std::endl;

}

void AZKMP_test(){
    cout << "<<< AZKMP TEST START >>>" << endl;

    ofstream fout1;
    ofstream fout2;
    fout1.open("AZKMP_p.txt");
    fout2.open("AZKMP_v.txt");

    uint32_t numItems = 1024;
    uint32_t numBatches = 1024;

    FrVec IDSet;
    IDSet.resize(numItems);

    for (uint32_t i = 0; i < numItems; i++) {
        IDSet[i].setByCSPRNG();
    }
    FrVec A = constructMemPoly(IDSet);

    vector<FrVec> Is;

    for (uint32_t i = 0; i < numBatches; i++) {
        Is.push_back({IDSet[i], 1});
    }


    PCS pcs;
    pcs.setup(numItems+1);

    zkbpacc_setup setup;
    setup.setup_from_pcs(pcs);    

    G2 C_A = pcs.commit(A);
    Fr r_A;
    r_A.setByCSPRNG();    
    C_A = C_A+setup.h2*r_A;
    vector<G2> C_Is; vector<G1> pi_Is; FrVec r_Is;

    G2 _tmp; Fr _tmp_Fr;

    for (uint32_t i = 0; i < numBatches; i++) {
        _tmp = pcs.commit(Is[i]);
        _tmp_Fr.setByCSPRNG();
        _tmp = _tmp + setup.h2 * _tmp_Fr;
        C_Is.push_back(_tmp); r_Is.push_back(_tmp_Fr);
        pi_Is.push_back(
            pcs.commit_G1(PolyLongDiv(A, Is[i]))
        );
    }
    cout << "Setup Done!!" << endl;
    

    microseconds total_p(0);
    microseconds total_v(0);
    microseconds elaps;

    uint64_t ps; uint64_t g1Size; uint64_t gtSize; uint64_t scSize;

    for(uint32_t i=0;i<100;i++){

        AZKMP azkmp(setup); ps = 0;

        auto start_1 = chrono::steady_clock::now();

        azkmp.prove(C_Is, C_A, pi_Is, r_Is, r_A);

        auto end_1 = chrono::steady_clock::now();
        elaps = duration_cast<microseconds>(end_1-start_1);
        fout1<< elaps.count()<< ", ";
        total_p = total_p+elaps;


        auto start_2 = chrono::steady_clock::now();

        bool flag = AZKMP_verify(&azkmp);


        if (flag == 0) {
            throw std::runtime_error("Verification Failed on AZKMP_test");
        }

        auto end_2 = chrono::steady_clock::now();
        elaps = duration_cast<microseconds>(end_2-start_2);
        fout2<< elaps.count()<< ", ";
        total_v = total_v+elaps;
        
        // Proof Size Calculator
        g1Size = sizeof(azkmp.msg.P_1);
        gtSize = sizeof(azkmp.msg.R_3);
        scSize = sizeof(azkmp.response.s_r_A);
        ps += g1Size * 3 + g1Size * azkmp.msg.P_2s.size() + gtSize;
        ps += scSize * (
            3 + azkmp.response.s_r_Is.size() + azkmp.response.s_delta_1s.size()
            + azkmp.response.s_delta_2s.size()
        );

        if (i == 0) {
            cout << g1Size << gtSize << scSize << endl;
        }
    }

    std::cout << "AZKMP_prove = : " << total_p.count()/100 << "(μs)"<<std::endl;
    std::cout << "AZKMP_verifiy = : " << total_v.count()/100 << "(μs)"<<std::endl;

    std::cout << "AZKMP_prove Amorized = : " << total_p.count()/100/numBatches << "(μs)"<<std::endl;
    std::cout << "AZKMP_verifiy Amorized = : " << total_v.count()/100/numBatches << "(μs)"<<std::endl;    

    std::cout << "Proof Size (bytes) = : " << ps << std::endl;
}

void ZKNMP_test(){
    cout << "<<< ZKNMP TEST START >>>" << endl;

    ofstream fout1;
    ofstream fout2;
    fout1.open("ZKNMP_p.txt");
    fout2.open("ZKNMP_v.txt");

    FrVec polyA = {6, 11, 6, 1};
    FrVec polyI = {4, 1};

    FrvT_3 out=xGCD(polyA, polyI);

    FrVec alpha = get<0>(out);
    FrVec beta = get<1>(out);
    FrVec gcd = get<2>(out);
    PCS pcs1;
    pcs1.setup(8);

    zkbpacc_setup setup1;
    setup1.setup_from_pcs(pcs1);    

    G1 w_A = pcs1.commit_G1(alpha);
    G1 w_I = pcs1.commit_G1(beta);

    G2 C_I = pcs1.commit(polyI);
    G2 C_A = pcs1.commit(polyA);

    Fr r_I;
    r_I.setByCSPRNG();
    Fr r_A;
    r_A.setByCSPRNG();
    C_I = C_I+setup1.h2*r_I;
    C_A = C_A+setup1.h2*r_A;

    ZKNMP nmp(setup1);

    microseconds total_p(0);
    microseconds total_v(0);
    microseconds elaps;

    for(uint32_t i=0;i<1000;i++){

        auto start_1 = chrono::steady_clock::now();

        nmp.prove(C_I, C_A, w_I, w_A, r_I, r_A);

        auto end_1 = chrono::steady_clock::now();
        elaps = duration_cast<microseconds>(end_1-start_1);
        fout1<< elaps.count()<< ", ";
        total_p = total_p+elaps;

        auto start_2 = chrono::steady_clock::now();

        bool flag = ZKNMP_verify(&nmp);

        if (flag == 0) {
            throw std::runtime_error("Verification Failed on ZKNMP_test");
        }        

        auto end_2 = chrono::steady_clock::now();
        elaps = duration_cast<microseconds>(end_2-start_2);
        fout2<< elaps.count()<< ", ";
        total_v = total_v+elaps;

    }

    std::cout << "ZKNMP_prove = : " << total_p.count()/1000 << "(μs)"<<std::endl;
    std::cout << "ZKNMP_verifiy = : " << total_v.count()/1000 << "(μs)"<<std::endl;

}



void ZKSP_test(){
    cout << "<<< ZKSP TEST START >>>" << endl;

    ofstream fout1;
    ofstream fout2;
    fout1.open("ZKSP_p.txt");
    fout2.open("ZKSP_v.txt");

    FrVec Ax = {24, -50, 35, -10, 1};
    FrVec Ix = {2, -3, 1};
    FrVec Jx = {12, -7, 1};

    PCS pcs2;
    pcs2.setup(8);
    zkbpacc_setup setup2;
    setup2.setup_from_pcs(pcs2);    

    G2 C_A = pcs2.commit(Ax);
    G2 C_I = pcs2.commit(Ix);
    G2 C_J = pcs2.commit(Jx);

    G1 pi_I = pcs2.commit_G1(PolyLongDiv(Ax,Ix));
    G1 pi_J = pcs2.commit_G1(PolyLongDiv(Ax,Jx));

    Fr r_A;
    r_A.setByCSPRNG();
    Fr r_I;
    r_I.setByCSPRNG();
    Fr r_J;
    r_J.setByCSPRNG();

    C_A = C_A + setup2.h2*r_A;
    C_I = C_I + setup2.h2*r_I;
    C_J = C_J + setup2.h2*r_J;

    ZKSP pi3(setup2);

    microseconds total_p(0);
    microseconds total_v(0);
    microseconds elaps;

    for(uint32_t i=0;i<1000;i++){

        auto start_1 = chrono::steady_clock::now();

        pi3.prove(C_I, C_J, C_A, pi_I, pi_J, r_I, r_J, r_A);

        auto end_1 = chrono::steady_clock::now();
        elaps = duration_cast<microseconds>(end_1-start_1);
        fout1<< elaps.count()<< ", ";
        total_p = total_p+elaps;

        auto start_2 = chrono::steady_clock::now();

        bool flag = ZKSP_verify(&pi3);

        if (flag == 0) {
            throw std::runtime_error("Verification Failed on ZKSP_test");
        }                

        auto end_2 = chrono::steady_clock::now();
        elaps = duration_cast<microseconds>(end_2-start_2);
        fout2<< elaps.count()<< ", ";
        total_v = total_v+elaps;

    }

    std::cout << "ZKSP_prove = : " << total_p.count()/1000 << "(μs)"<<std::endl;
    std::cout << "ZKSP_verifiy = : " << total_v.count()/1000 << "(μs)"<<std::endl;

}

void ZKIPP_test(const std::string& file_name){
    cout << "<<< ZKIPP TEST START >>>" << endl;

    vector<G2> gg;
    vector<G1> ww;

    std::ofstream fout1(file_name, std::ios::app);

    fout1<<"input size, "<<"Prove time"<<"verification time, "<<endl;

    //int values[] = {127, 255, 511, 1023, 2047, 4095, 8191, 16383};

    for(int i=128; i<=16384; i=i*2){

        
        i=i-1;
        microseconds elaps;

        i=nextPowOf2(i);
        fout1<< i << ", ";

        gg.resize(i);
        ww.resize(i);

        for(uint32_t i=0;i<ww.size();i++){

            hashAndMapToG1(ww[i],to_string(i));
            hashAndMapToG2(gg[i],to_string(i));

        }

        GT P =MultiPairing(ww,gg);

        auto start_1 = chrono::steady_clock::now();

        IPPproof ipp_pi = zkIPPprove(gg,P,ww);

        auto end_1 = chrono::steady_clock::now();
        elaps = duration_cast<microseconds>(end_1-start_1);
        fout1<< elaps.count()<< ", ";
        

        auto start_2 = chrono::steady_clock::now();

        bool flag = zkIPPverify(&ipp_pi);

        if (flag == 0) {
            throw std::runtime_error("Verification Failed on ZKIPP_test");
        }                

        auto end_2 = chrono::steady_clock::now();

        elaps = duration_cast<microseconds>(end_2-start_2);
        fout1<< elaps.count()<< ", ";

        fout1<<endl;
        


    }


}

void SetupTest(const std::string& file_name){
    cout << "<<< SETUP TEST START >>>" << endl;

    std::ofstream fout1(file_name, std::ios::app);
    microseconds elaps;

    uint32_t vals[4] = {5000, 20000, 80000, 320000};
        

        for(uint32_t i:vals){
            
            for(uint32_t j=0;j<100;j++){

                fout1<<i<<", ";
                

                FrVec IDset;
                IDset.resize(i);
                for(uint32_t k=0;k<i;k++){
                    IDset[k].setByCSPRNG();
                }

                auto start_1 = chrono::steady_clock::now();

                PCS pcs;
                pcs.setup(i+1);
                zkbpacc_setup setup;
                setup.setup_from_pcs(pcs);
                // setup.init(i+1);
                setup.g1si = &pcs.pp.g1si[0];
                setup.g2si = &pcs.pp.g2si[0];

                auto end_1 = chrono::steady_clock::now();
                elaps = duration_cast<microseconds>(end_1-start_1);
                fout1<< elaps.count()<< ", ";

                Fr s = pcs.s;

                start_1 = chrono::steady_clock::now();

                // FrVec Ax = Polytree(IDset);
                // G2 C_A = pcs.commit(Ax);

                Fr Ax=1;
                G2 C_A;

                for(uint32_t k=0;k<IDset.size();k++){

                    Ax = Ax*(s+IDset[k]);

                    // float progress = (float)k / (IDset.size() - 1);
                    // int barWidth = 70;
                    
                    // std::cout << "[";
                    // int pos = barWidth * progress;
                    // for (int i = 0; i < barWidth; ++i) {
                    //     if (i < pos) std::cout << "=";
                    //     else if (i == pos) std::cout << ">";
                    //     else std::cout << " ";
                    // }
                    // std::cout << "] " << int(progress * 100.0) << " %\r";
                    // std::cout.flush();

                }

                C_A = setup.g2si[0]*Ax;

                end_1 = chrono::steady_clock::now();
                elaps = duration_cast<microseconds>(end_1-start_1);
                fout1<< elaps.count()<< ", ";

                start_1 = chrono::steady_clock::now();

                Fr id = IDset[rand()%i];
                Fr Ix = s+id;
                G2 C_id = setup.g2si[0]*(s+id);

                end_1 = chrono::steady_clock::now();
                elaps = duration_cast<microseconds>(end_1-start_1);
                fout1<< elaps.count()<< ", ";
                

                start_1 = chrono::steady_clock::now();

                G1 wit = setup.g1si[0]*(Ax/Ix);

                end_1 = chrono::steady_clock::now();
                elaps = duration_cast<microseconds>(end_1-start_1);
                fout1<< elaps.count()<< ", ";
                fout1<<endl;

            }

        }
    
}

// #define PCS_BASE_LENGTH ( 2500 )
// #define SET_SIZE 2500


void OwnPf_test() {
    cout << "<<< OwnPf TEST START >>>" << endl;

    uint32_t vals[4] = {250, 500, 1000, 2000}; 

    // Initialize Stuffs
    for (auto val : vals) {
        cout << "Current Set Size: " << val << endl;


        PCS pcs;
        pcs.setup(val);
        zkbpacc_setup setup; setup.setup_from_pcs(pcs);
        ZKMP zkmp(setup);
        PoK2_G2 pipok(setup);
        OwnPf proof = {setup,zkmp,pipok};
    
        // Initialize Inputs 
        FrVec A(val);
        for (int i = 0; i < val; i++) {
            A[i] = Fr(i + 2);
        }    
        FrVec APoly = constructMemPoly(A);
    
        Fr n = 2; Fr id = 424242; FrVec idPoly = {id, 1};
        G2 C_id = pcs.commit(idPoly);
        G2 C_A = pcs.commit(APoly);
        Fr r, delta;
        r.setByCSPRNG();
        delta.setByCSPRNG();
        C_id = C_id + setup.h2 * r;
        C_A = C_A + setup.h2 * delta;
    
        microseconds total_p(0);
        microseconds total_v(0);
        microseconds elaps;
    
        for (int i = 0; i < 1000; i++) {
            auto start_1 = chrono::steady_clock::now();
            proof.prove(C_id, C_A, n, A, id, r, delta);
            auto end_1 = chrono::steady_clock::now();
            elaps = duration_cast<microseconds>(end_1-start_1);
            total_p = total_p+elaps;
    
            auto start_2 = chrono::steady_clock::now();
    
            bool flag = OwnPf_verify(&proof);
    
            if (flag == 0) {
                throw std::runtime_error("OwnPf Verification Failed.");
            }        
    
            auto end_2 = chrono::steady_clock::now();
            elaps = duration_cast<microseconds>(end_2-start_2);
            total_v = total_v+elaps;        
        }
    
        std::cout << "OwnPf_prove = : " << total_p.count()/1000 << "(μs)"<<std::endl;
        std::cout << "OwnPf_verifiy = : " << total_v.count()/1000 << "(μs)"<<std::endl;
    }
}


void NonOwnPf_test() {
    cout << "<<< NonOwnPf TEST START >>>" << endl;

    uint32_t vals[4] = {250, 500, 1000, 2000}; 

    for (auto val : vals) {
        cout << "Current Set Size: " << val << endl;

        // Initialize Stuffs
        PCS pcs;
        pcs.setup(val);
        zkbpacc_setup setup; setup.setup_from_pcs(pcs);
        ZKNMP zknmp(setup);
        PoK2_G2 pipok(setup);
        NonOwnPf proof = {setup,zknmp,pipok};

        // Initialize Inputs 
        FrVec A(val);
        for (int i = 0; i < val; i++) {
            A[i] = Fr(i + 2);
        }    
        FrVec APoly = constructMemPoly(A);

        Fr n = 1; Fr id = 424242; FrVec idPoly = {id, 1};
        G2 C_id = pcs.commit(idPoly);
        G2 C_A = pcs.commit(APoly);
        Fr r, delta;
        r.setByCSPRNG();
        delta.setByCSPRNG();
        C_id = C_id + setup.h2 * r;
        C_A = C_A + setup.h2 * delta;

        microseconds total_p(0);
        microseconds total_v(0);
        microseconds elaps;

        for (int i = 0; i < 1000; i++) {
            auto start_1 = chrono::steady_clock::now();
            proof.prove(C_id, C_A, n, A, id, r, delta);
            auto end_1 = chrono::steady_clock::now();
            elaps = duration_cast<microseconds>(end_1-start_1);
            total_p = total_p+elaps;

            auto start_2 = chrono::steady_clock::now();

            bool flag = NonOwnPf_verify(&proof);

            if (flag == 0) {
                throw std::runtime_error("NonOwnPf Verification Failed.");
            }        

            auto end_2 = chrono::steady_clock::now();
            elaps = duration_cast<microseconds>(end_2-start_2);
            total_v = total_v+elaps;        
        }

        std::cout << "NonOwnPf_prove = : " << total_p.count()/1000 << "(μs)"<<std::endl;
        std::cout << "NonOwnPf_verifiy = : " << total_v.count()/1000 << "(μs)"<<std::endl;
        
    }
}

void TxGenTest() {

    cout << "<<< TxGen Test >>>" << endl;

    vector<uint32_t> numTxs = {
        256, 512, 1024, 2048
    };



    for (auto val : numTxs) {
        cout << "Current # of Txs: " << val << endl;

        PCS pcs;
        pcs.setup(4243);

        // Setup by Authority
        FrVec IDSet; IDSet.resize(4242);
        for (uint32_t i = 0; i < 4242; i++) {
            IDSet[i].setByCSPRNG();
        }
        Fr id1 = IDSet[42];
        Fr id2 = IDSet[56];
        G1 wit1 = memWitGen(pcs, IDSet, id1);
        G1 wit2 = memWitGen(pcs, IDSet, id2);
        G2 C_A = pcs.commit(constructMemPoly(IDSet));
        
        // Prepare SNs
        FrVec SN; SN.resize(val);
        for (uint32_t i = 0; i < val; i++) {
            SN[i].setByCSPRNG();
        }

        FrVec S1 = FrVec(SN.begin(), SN.begin() + val/5);
        FrVec S2 = FrVec(SN.begin() + val/5, SN.end());

        // Test for Tx_Entry 
        // id1 enter with SN

        auto start_entry = chrono::steady_clock::now();

        Fr gamma; Fr delta;
        gamma.setByCSPRNG(); delta.setByCSPRNG();
        Tx_Entry TxEntry = Tx_Entry_Gen(
            pcs, C_A, SN, id1, wit1, gamma, delta
        );

        auto end_entry = chrono::steady_clock::now();

        auto start_entry_vrfy = chrono::steady_clock::now();

        bool entry_vrfy = Tx_Entry_Vrfy(TxEntry, pcs);

        auto end_entry_vrfy = chrono::steady_clock::now();


        // Test for Tx_Transfer
        // Transfer id1 -> id2 via S1

        auto start_trans = chrono::steady_clock::now();

        Fr gamma1; Fr gamma2; Fr delta1; Fr delta2;
        gamma1.setByCSPRNG(); gamma2.setByCSPRNG();
        delta1.setByCSPRNG(); delta2.setByCSPRNG();        
        Tx_Trans TxTrans = Tx_Trans_Gen(
            TxEntry.Tx_out, pcs, C_A,
            SN, S1, S2, id1, id2, wit1, wit2,
            gamma, delta,
            gamma1, gamma2, delta1, delta2
        );

        auto end_trans = chrono::steady_clock::now();

        auto start_trans_vrfy = chrono::steady_clock::now();

        bool trans_vrfy = Tx_Trans_Vrfy(TxTrans, pcs);

        auto end_trans_vrfy = chrono::steady_clock::now();        

        // Test for Tx_Exit
        /// id1 deletes S1 from S

        auto start_exit = chrono::steady_clock::now();

        Fr gammap; gammap.setByCSPRNG();
        Tx_Exit  TxExit = Tx_Exit_Gen(
            TxEntry.Tx_out, pcs, C_A, 
            SN, S1, S2, id1, wit1, gamma, delta,
            gammap, delta1, delta2
        );

        auto end_exit = chrono::steady_clock::now();        

        auto start_exit_vrfy = chrono::steady_clock::now();

        bool exit_vrfy = Tx_Exit_Vrfy(TxExit, pcs);

        auto end_exit_vrfy = chrono::steady_clock::now();        

        cout << "Entry: " 
             << duration_cast<microseconds>(end_entry - start_entry).count()
             << "us"
             << endl;
        cout << "Entry Verify: " 
             << duration_cast<microseconds>(end_entry_vrfy - start_entry_vrfy).count()
             << "us"
             << endl;             

        cout << "Transfer: " 
             << duration_cast<microseconds>(end_trans - start_trans).count()
             << "us"
             << endl;
        cout << "Transfer Verfy: " 
             << duration_cast<microseconds>(end_trans_vrfy - start_trans_vrfy).count()
             << "us"
             << endl;             
             
        cout << "Exit: " 
             << duration_cast<microseconds>(end_exit - start_exit).count()
             << "us"
             << endl;             
        cout << "Exit Verify: " 
             << duration_cast<microseconds>(end_exit_vrfy - start_exit_vrfy).count()
             << "us\n"
             << endl;                  
    }
}

void TxSummaryTest() {
    cout << "<<< TxSummary Test >>>" << endl;

    vector<uint32_t> numTxs = {
        128, 256, 512, 1024, 2048
    };

    for (auto val : numTxs) {
        cout << "Current # of Txs * 2: " << val << endl;
        PCS pcs;
        pcs.setup(4243);

        FrVec IDSet; IDSet.resize(4242);
        for (uint32_t i = 0; i < 4242; i++) {
            IDSet[i].setByCSPRNG();
        }        
        G2 C_A = pcs.commit(constructMemPoly(IDSet));

        // Simulate the aggregation scenario.        
        // Transactions that needs to be aggregated
        FrVec SN_loc; vector<G2> C_Is; vector<G1> W_Is; FrVec gamma_Is;
        SN_loc.resize(val); C_Is.resize(val);
        W_Is.resize(val); gamma_Is.resize(val);

        for (uint32_t i = 0; i < val; i++) {
            SN_loc[i] = IDSet[i];
            gamma_Is[i].setByCSPRNG();
            C_Is[i] = pcs.commit({SN_loc[i], 1}) + pcs.pp.h2si[0] * gamma_Is[i];
        }
        for (uint32_t i = 0; i < val; i++) {
            W_Is[i] = memWitGen(pcs, SN_loc, SN_loc[i]);
        }
        

        // Local Manager's Transcript
        // Fr ID_loc = IDSet[4200];
        Fr gamma_loc; gamma_loc.setByCSPRNG();
        // G1 W_loc = memWitGen(pcs, IDSet, ID_loc);
        G2 C_loc = pcs.commit(constructMemPoly(SN_loc)) + pcs.pp.h2si[0] * gamma_loc;
        G1 W_loc = pcs.commit_G1(constructMemPolyBatch (IDSet, SN_loc));
        
        // Some dummy Tx
        // This does not affect the correctness of the test
        Fr Adr(42);
        Tx_Base Tx_in = {Adr, C_loc, C_loc};        

        // Do Summary
        auto start_summary = chrono::steady_clock::now();

        Fr gamma_locp; gamma_locp.setByCSPRNG();
        Tx_Summary TxSum = Tx_Summary_Gen (
            Tx_in, pcs, SN_loc, C_A, C_Is,
            gamma_loc, W_loc, W_Is, gamma_Is,
            gamma_locp
        );

        auto end_summary = chrono::steady_clock::now();


        // Summary Verification 
        auto start_vrfy = chrono::steady_clock::now();

        bool ret_sum = Tx_Summary_Vrfy(TxSum, pcs);

        auto end_vrfy = chrono::steady_clock::now();

        // if (~ret_sum) {
        //     throw std::runtime_error("Verification Failed...");
        // }

        cout << "Summary: "
             << duration_cast<microseconds>(end_summary - start_summary).count()
             << "us"
             << endl;
        cout << "Summary Verify: "
             << duration_cast<microseconds>(end_vrfy - start_vrfy).count()
             << "us\n"
             << endl;             
    }
}





// Deprecated
// void TxGenTest() {
//     PCS pcs;
//     pcs.setup(2001);
//     zkbpacc_setup setup;
//     setup.setup_from_pcs(pcs);
//     setup.g1si = &pcs.pp.g1si[0];
//     setup.g2si = &pcs.pp.g2si[0];

//     ZKMP PI_PoA(setup);
//     FrVec IDset;
//     for(uint32_t i=0;i<10;i++){
//         Fr id;id.setByCSPRNG();
//         IDset.push_back(id);
//     }
//     srand(time(0));
//     uint32_t j = rand()%10;
//     Fr id = IDset[j];

//     const std::string& file_name = "EntryTxGen.csv"; 
//     ZKMP PI_PoA1(setup);

//     G1 wit = pcs.commit_G1(PolyLongDiv(Polytree(IDset),Polytree({id})));
//     G2 C_A = pcs.commit(Polytree(IDset));

//     FrVec SNs;
//     for(uint32_t i=0;i<125;i++){
//         Fr SN;SN.setByCSPRNG();
//         SNs.push_back(SN);
//     }

//     FrVec SNs2;
//     for(uint32_t i=0;i<250;i++){
//         Fr SN;SN.setByCSPRNG();
//         SNs2.push_back(SN);
//     }

//     FrVec SNs3;
//     for(uint32_t i=0;i<500;i++){
//         Fr SN;SN.setByCSPRNG();
//         SNs3.push_back(SN);
//     }

//     FrVec SNs4;
//     for(uint32_t i=0;i<1000;i++){
//         Fr SN;SN.setByCSPRNG();
//         SNs4.push_back(SN);
//     }

//     FrVec SNs5;
//     for(uint32_t i=0;i<2000;i++){
//         Fr SN;SN.setByCSPRNG();
//         SNs5.push_back(SN);
//     }

//     Fr gamma; gamma.setByCSPRNG();
//     Fr delta; delta.setByCSPRNG();

//     ZKMP PI_PoA_1(setup);
//     ZKMP PI_PoA_2(setup);
//     ZKSP PI_PoQ(setup);

//     uint32_t k = rand()%10;
//     Fr id2 = IDset[k];
//     G1 wit2 = pcs.commit_G1(PolyLongDiv(Polytree(IDset),Polytree({id2})));
//     Fr gamma1; gamma1.setByCSPRNG();
//     Fr gamma2; gamma2.setByCSPRNG();
//     Fr delta1; delta1.setByCSPRNG();
//     Fr delta2; delta2.setByCSPRNG();

//     FrVec S_2(SNs.begin(),SNs.begin()+SNs.size()/5);
//     FrVec S_1(SNs.begin()+SNs.size()/5,SNs.end());
    
//     FrVec S_4(SNs2.begin(),SNs2.begin()+SNs2.size()/5);
//     FrVec S_3(SNs2.begin()+SNs2.size()/5,SNs2.end());

//     FrVec S_6(SNs3.begin(),SNs3.begin()+SNs3.size()/5);
//     FrVec S_5(SNs3.begin()+SNs3.size()/5,SNs3.end());

//     FrVec S_8(SNs4.begin(),SNs4.begin()+SNs4.size()/5);
//     FrVec S_7(SNs4.begin()+SNs4.size()/5,SNs4.end());

//     FrVec S_10(SNs5.begin(),SNs5.begin()+SNs5.size()/5);
//     FrVec S_9(SNs5.begin()+SNs5.size()/5,SNs5.end());

//     for(int i=0;i<100;i++){

//         Tx_Entry Tx1 = Tx_Entry_Gen(pcs, PI_PoA1, id, wit, C_A, SNs, gamma, delta, "Tx_Entry_125.csv");
//         Tx_Entry Tx2 = Tx_Entry_Gen(pcs, PI_PoA1, id, wit, C_A, SNs2, gamma, delta, "Tx_Entry_250.csv");
//         Tx_Entry Tx3 = Tx_Entry_Gen(pcs, PI_PoA1, id, wit, C_A, SNs3, gamma, delta, "Tx_Entry_500.csv");
//         Tx_Entry Tx4 = Tx_Entry_Gen(pcs, PI_PoA1, id, wit, C_A, SNs4, gamma, delta, "Tx_Entry_1000.csv");   
//         Tx_Entry Tx5 = Tx_Entry_Gen(pcs, PI_PoA1, id, wit, C_A, SNs5, gamma, delta, "Tx_Entry_2000.csv");

//         Tx_Transfer Tx6 = Tx_Transfer_Gen(pcs, PI_PoA_1, PI_PoA_2, PI_PoQ, id, id2, wit, wit2, C_A, SNs, 
//         S_1, S_2, Tx1.output1, gamma, gamma1, gamma2, delta, delta1, delta2, "Transfer_125.csv", Tx1.Tid);

//         Tx_Transfer Tx7 = Tx_Transfer_Gen(pcs, PI_PoA_1, PI_PoA_2, PI_PoQ, id, id2, wit, wit2, C_A, SNs2, 
//         S_3, S_4, Tx2.output1, gamma, gamma1, gamma2, delta, delta1, delta2, "Transfer_250.csv", Tx2.Tid);

//         Tx_Transfer Tx8 = Tx_Transfer_Gen(pcs, PI_PoA_1, PI_PoA_2, PI_PoQ, id, id2, wit, wit2, C_A, SNs3, 
//         S_5, S_6, Tx3.output1, gamma, gamma1, gamma2, delta, delta1, delta2, "Transfer_500.csv", Tx3.Tid);

//         Tx_Transfer Tx9 = Tx_Transfer_Gen(pcs, PI_PoA_1, PI_PoA_2, PI_PoQ, id, id2, wit, wit2, C_A, SNs4, 
//         S_7, S_8, Tx4.output1, gamma, gamma1, gamma2, delta, delta1, delta2, "Transfer_1000.csv", Tx4.Tid);

//         Tx_Transfer Tx10 = Tx_Transfer_Gen(pcs, PI_PoA_1, PI_PoA_2, PI_PoQ, id, id2, wit, wit2, C_A, SNs5, 
//         S_9, S_10, Tx5.output1, gamma, gamma1, gamma2, delta, delta1, delta2, "Transfer_2000.csv", Tx5.Tid);

//         Tx_Exit Tx11 = Tx_Exit_Gen(pcs,Tx6.output1, PI_PoA, PI_PoQ, id2, wit2, C_A, gamma1, gamma2, delta1, delta2, SNs, {}, SNs, "Exit_125.csv");
//         Tx_Exit Tx12 = Tx_Exit_Gen(pcs,Tx7.output1, PI_PoA, PI_PoQ, id2, wit2, C_A, gamma1, gamma2, delta1, delta2, SNs2, {}, SNs2, "Exit_250.csv");
//         Tx_Exit Tx13 = Tx_Exit_Gen(pcs,Tx8.output1, PI_PoA, PI_PoQ, id2, wit2, C_A, gamma1, gamma2, delta1, delta2, SNs3, {}, SNs3, "Exit_500.csv");
//         Tx_Exit Tx14 = Tx_Exit_Gen(pcs,Tx9.output1, PI_PoA, PI_PoQ, id2, wit2, C_A, gamma1, gamma2, delta1, delta2, SNs4, {}, SNs4, "Exit_1000.csv");
//         Tx_Exit Tx15 = Tx_Exit_Gen(pcs,Tx10.output1, PI_PoA, PI_PoQ, id2, wit2, C_A, gamma1, gamma2, delta1, delta2, SNs5, {}, SNs5, "Exit_2000.csv");
//     }
// }

// void TxVrfyTest() {
//     cout << "<<< TxGEN Verify TEST >>>" << endl;

//     PCS pcs; 
//     pcs.setup(8);
//     zkbpacc_setup setup;
//     setup.setup_from_pcs(pcs);

//     FrVec IDset(8);
//     for (uint32_t i = 0; i < 8; i++) {
//         IDset[i].setByCSPRNG();
//     }
//     Fr id = IDset[4];
//     Fr id2 = IDset[5];
//     G1 wit = pcs.commit_G1(PolyLongDiv(Polytree(IDset),Polytree({id})));
//     G1 wit2 = pcs.commit_G1(PolyLongDiv(Polytree(IDset),Polytree({id2})));    
//     G2 C_A = pcs.commit(constructMemPoly(IDset));

//     FrVec SNs(100);
//     for (uint32_t i = 0; i < 100; i++) {
//         SNs[i].setByCSPRNG();
//     }
//     Fr gamma; gamma.setByCSPRNG();
//     Fr delta; delta.setByCSPRNG();
//     Fr gamma1; gamma1.setByCSPRNG();
//     Fr gamma2; gamma2.setByCSPRNG();
//     Fr delta1; delta1.setByCSPRNG();
//     Fr delta2; delta2.setByCSPRNG();



//     ZKMP PI_PoA(setup);
//     ZKMP PI_PoA1(setup);
//     ZKMP PI_PoA_1(setup);
//     ZKMP PI_PoA_2(setup);
//     ZKSP PI_PoQ(setup);


//     microseconds entry(0);
//     microseconds transfer(0);
//     microseconds exit(0);
//     microseconds elaps;

//     for (int i = 0; i < 1000; i++) {
//         Tx_Entry txe = Tx_Entry_Gen(pcs, PI_PoA1, id, wit, C_A, SNs, gamma, delta, "tt1");
//         Tx_Transfer txt = Tx_Transfer_Gen(pcs, PI_PoA_1, PI_PoA_2, PI_PoQ, id, id2, wit, wit2, C_A, SNs, 
//         {}, SNs, txe.output1, gamma, gamma1, gamma2, delta, delta1, delta2, "tt2", txe.Tid); 
//         Tx_Exit txo = Tx_Exit_Gen(pcs,txt.output1, PI_PoA, PI_PoQ, id2, wit2, C_A, gamma1, gamma2, delta1, delta2, SNs, {}, SNs, "tt3");       

//         // Entry Verification
//         auto start_1 = chrono::steady_clock::now();
//         PoK_proof_g2_verify(txe.SN.PI_PoK);
//         ZKMP_verify(&txe.output1.PI_PoA);
//         auto end_1 = chrono::steady_clock::now();
//         elaps = duration_cast<microseconds>(end_1 - start_1);
//         entry += elaps;
        

//         // Transfer Verification
//         auto start_2 = chrono::steady_clock::now();
//         PoK_proof_g2_n_verify(&txt.input.PI_PoKE);
//         ZKMP_verify(&txt.output1.PI_PoA);
//         ZKMP_verify(&txt.output2.PI_PoA);
//         ZKSP_verify(&txt.PI_PoQ);
//         auto end_2 = chrono::steady_clock::now();
//         elaps = duration_cast<microseconds>(end_2 - start_2);
//         transfer += elaps;
        
//         // Exit Verification
//         auto start_3 = chrono::steady_clock::now();
//         PoK_proof_g2_n_verify(&txo.input.PI_PoKE);
//         ZKMP_verify(&txo.output_.PI_PoA);
//         ZKSP_verify(&txo.PI_PoQ);
//         auto end_3 = chrono::steady_clock::now();        
//         elaps = duration_cast<microseconds>(end_3 - start_3);
//         exit += elaps;
//     }

//     std::cout << "TxEntry_Vrfy: " << entry.count() / 1000 << "(μs)"<<std::endl;
//     std::cout << "TxTransfer_Vrfy: " << transfer.count() / 1000 << "(μs)"<<std::endl;
//     std::cout << "TxExit_Vrfy: " << exit.count() / 1000 << "(μs)"<<std::endl;
// }



// Deprecated
// void TxAggTest() {
//     std::string file_name = "Agg_ZKMP.csv";
//     std::string file_name1 = "Agg_ZKMP_commit.csv";     
//     std::string file_name2 = "Agg_ZKMP_PME.csv";  
//     std::string file_name3 = "Agg_ZKMP_ZKMP.csv";  
//     std::string file_name4 = "Agg_ZKMP_ZKIPP.csv";  
//     std::string file_name5 = "Agg_ZKMP_PCopenP.csv";  
//     std::string file_name6 = "Agg_ZKMP_StructGen.csv";  
    
//     AggZKMP_test(128,file_name, file_name1, file_name2, file_name3, file_name4, file_name5, file_name6);
//     AggZKMP_test(256,file_name, file_name1, file_name2, file_name3, file_name4, file_name5, file_name6);
//     AggZKMP_test(512,file_name, file_name1, file_name2, file_name3, file_name4, file_name5, file_name6);
//     AggZKMP_test(1024,file_name, file_name1, file_name2, file_name3, file_name4, file_name5, file_name6);
//     AggZKMP_test(2048,file_name, file_name1, file_name2, file_name3, file_name4, file_name5, file_name6);
// }


int main(){
    initPairing(mcl::BLS12_381);
    // Core Protocols
    ZKMP_test();
    AZKMP_test();
    ZKSP_test();    
    ZKNMP_test();
    OwnPf_test();
    NonOwnPf_test();
    // Tx Generation & Verification
    TxGenTest();
    TxSummaryTest();
    // Setup Phase 
    SetupTest("setup_new.csv");
    return 0;
}
