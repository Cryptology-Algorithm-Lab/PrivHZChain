#include "test.h"

#include <algorithm>
#include <numeric>
#include <utility>
#include <string>
#include <vector>

using Clock = std::chrono::steady_clock;

static inline double us_since(Clock::time_point a, Clock::time_point b) {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count() / 1000.0;
}

static inline double mean_of(const std::vector<double>& v) {
    return v.empty() ? 0.0 : std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}

struct Bench {
    struct Row {
        std::string test;
        long        size;
        std::string metric;
        double mean, median, p95, min, max;
        long   n;
        std::string unit = "us";   // "us" for timings, "B" for sizes
    };
    std::vector<Row> rows;
    bool verbose = true;

    void echo(const Row& r) const {
        if (!verbose) return;
        bool isTime = (r.unit == "us");
        std::cout << "  [" << std::left << std::setw(10) << r.test << "]"
                  << " size=" << std::setw(7) << std::left << r.size
                  << std::setw(16) << std::left << r.metric << std::right
                  << std::fixed << std::setprecision(isTime ? 2 : 0);
        if (r.n > 1)
            std::cout << " mean=" << std::setw(11) << r.mean << r.unit
                      << "  (median " << r.median << ", p95 " << r.p95
                      << ", n=" << r.n << ")";
        else
            std::cout << " " << std::setw(11) << r.mean << " " << r.unit;
        std::cout << "\n";
    }

    void record(const std::string& test, long size,
                const std::string& metric, std::vector<double> s) {
        if (s.empty()) return;
        std::sort(s.begin(), s.end());
        Row r;
        r.test = test; r.size = size; r.metric = metric; r.n = (long)s.size();
        r.min = s.front(); r.max = s.back();
        r.mean = mean_of(s);
        r.median = s[s.size() / 2];
        r.p95 = s[std::min(s.size() - 1, (size_t)(s.size() * 0.95))];
        rows.push_back(r);
        echo(r);
    }

    void value(const std::string& test, long size, const std::string& metric,
               double v, const std::string& unit = "us") {
        Row r{test, size, metric, v, v, v, v, v, 1, unit};
        rows.push_back(r);
        echo(r);
    }

    void dumpConsole() const {
        std::cout << "\n========================= BENCHMARK SUMMARY =========================\n";
        std::cout << std::left
                  << std::setw(12) << "test"
                  << std::setw(9)  << "size"
                  << std::setw(18) << "metric"
                  << std::right
                  << std::setw(13) << "mean"
                  << std::setw(13) << "median"
                  << std::setw(13) << "p95"
                  << std::setw(7)  << "n"
                  << std::setw(7)  << "unit"
                  << "\n";
        std::cout << std::string(91, '-') << "\n";
        std::string cur;
        for (const auto& r : rows) {
            if (r.test != cur) {
                if (!cur.empty()) std::cout << "\n";
                cur = r.test;
            }
            bool isTime = (r.unit == "us");
            std::cout << std::left
                      << std::setw(12) << r.test
                      << std::setw(9)  << r.size
                      << std::setw(18) << r.metric
                      << std::right << std::fixed << std::setprecision(isTime ? 2 : 0)
                      << std::setw(13) << r.mean
                      << std::setw(13) << r.median
                      << std::setw(13) << r.p95
                      << std::setw(7)  << r.n
                      << std::setw(7)  << r.unit
                      << "\n";
        }
        std::cout << std::string(91, '=') << "\n";
    }

    void dumpCsv(const std::string& path) const {
        std::ofstream f(path);
        f << "test,size,metric,mean,median,p95,min,max,n,unit\n";
        f << std::fixed << std::setprecision(3);
        for (const auto& r : rows) {
            f << r.test << "," << r.size << "," << r.metric << ","
              << r.mean << "," << r.median << "," << r.p95 << ","
              << r.min << "," << r.max << "," << r.n << "," << r.unit << "\n";
        }
        std::cout << "[bench] wrote " << rows.size() << " rows to " << path << "\n";
    }
};

static Bench g_bench;

// throws if any verify returns false
template <class P, class V>
static std::pair<std::vector<double>, std::vector<double>>
benchProveVerify(int iters, P prove, V verify) {
    std::vector<double> ps, vs;
    ps.reserve(iters);
    vs.reserve(iters);
    for (int i = 0; i < iters; ++i) {
        auto t0 = Clock::now(); prove();
        auto t1 = Clock::now(); bool ok = verify();
        auto t2 = Clock::now();
        if (!ok) throw std::runtime_error("verification failed");
        ps.push_back(us_since(t0, t1));
        vs.push_back(us_since(t1, t2));
    }
    return {ps, vs};
}

// ===== Tests =====

void ZKMP_test() {
    cout << "<<< ZKMP TEST >>>" << endl;

    PCS pcs; pcs.setup(8);
    ZKAccSetup setup; setup.setup_from_pcs(pcs);

    FrVec A = {120, 274, 225, 85, 15, 1};
    FrVec I = {1, 1};

    G2 C_A = pcs.commit_g2(A);
    G2 C_I = pcs.commit_g2(I);
    Fr r_I; r_I.setByCSPRNG();
    Fr r_A; r_A.setByCSPRNG();
    C_I = C_I + setup.h2 * r_I;
    C_A = C_A + setup.h2 * r_A;

    G1 pi_I = pcs.commit_g1(PolyLongDiv(A, I));

    ZKMP zkmp(setup);

    auto [ps, vs] = benchProveVerify(1000,
        [&] { zkmp.prove(C_I, C_A, pi_I, r_I, r_A); },
        [&] { return ZKMP_verify(&zkmp); });

    g_bench.record("ZKMP", 8, "prove",  ps);
    g_bench.record("ZKMP", 8, "verify", vs);
}

void AZKMP_test() {
    cout << "<<< AZKMP TEST >>>" << endl;

    uint32_t numItems   = 1024;
    uint32_t numBatches = 1024;

    PCS pcs; pcs.setup(numItems + 1);
    ZKAccSetup setup; setup.setup_from_pcs(pcs);

    FrVec IDSet(numItems);
    for (uint32_t i = 0; i < numItems; i++) IDSet[i].setByCSPRNG();
    FrVec A = constructMemPoly(pcs.pp.nttctx, IDSet);

    vector<FrVec> Is;
    for (uint32_t i = 0; i < numBatches; i++) Is.push_back({IDSet[i], 1});

    G2 C_A = pcs.commit_g2(A);
    Fr r_A; r_A.setByCSPRNG();
    C_A = C_A + setup.h2 * r_A;

    vector<G2> C_Is; vector<G1> pi_Is; FrVec r_Is;
    G2 _tmp; Fr _tmp_Fr;
    for (uint32_t i = 0; i < numBatches; i++) {
        _tmp = pcs.commit_g2(Is[i]);
        _tmp_Fr.setByCSPRNG();
        _tmp = _tmp + setup.h2 * _tmp_Fr;
        C_Is.push_back(_tmp); r_Is.push_back(_tmp_Fr);
        pi_Is.push_back(pcs.commit_g1(PolyLongDiv(A, Is[i])));
    }
    cout << "Setup Done!!" << endl;

    AZKMP azkmp(setup);

    auto [ps, vs] = benchProveVerify(100,
        [&] { azkmp.prove(C_Is, C_A, pi_Is, r_Is, r_A); },
        [&] { return AZKMP_verify(&azkmp); });

    g_bench.record("AZKMP", numBatches, "prove",  ps);
    g_bench.record("AZKMP", numBatches, "verify", vs);
    g_bench.value("AZKMP", numBatches, "prove_amort",  mean_of(ps) / numBatches);
    g_bench.value("AZKMP", numBatches, "verify_amort", mean_of(vs) / numBatches);

    uint64_t g1Size = sizeof(azkmp.msg.P_1);
    uint64_t gtSize = sizeof(azkmp.msg.R_3);
    uint64_t scSize = sizeof(azkmp.response.s_r_A);
    uint64_t pbytes = g1Size * 3 + g1Size * azkmp.msg.P_2s.size() + gtSize;
    pbytes += scSize * (3 + azkmp.response.s_r_Is.size()
                          + azkmp.response.s_delta_1s.size()
                          + azkmp.response.s_delta_2s.size());
    g_bench.value("AZKMP", numBatches, "proof_bytes", (double)pbytes, "B");
}

void ZKNMP_test() {
    cout << "<<< ZKNMP TEST >>>" << endl;

    PCS pcs; pcs.setup(8);
    ZKAccSetup setup; setup.setup_from_pcs(pcs);

    FrVec polyA = {6, 11, 6, 1};
    FrVec polyI = {4, 1};

    FrvT_3 out = xGCD(pcs.pp.nttctx, polyA, polyI);
    FrVec alpha = get<0>(out);
    FrVec beta  = get<1>(out);

    G1 w_A = pcs.commit_g1(alpha);
    G1 w_I = pcs.commit_g1(beta);

    G2 C_I = pcs.commit_g2(polyI);
    G2 C_A = pcs.commit_g2(polyA);

    Fr r_I; r_I.setByCSPRNG();
    Fr r_A; r_A.setByCSPRNG();
    C_I = C_I + setup.h2 * r_I;
    C_A = C_A + setup.h2 * r_A;

    ZKNMP nmp(setup);

    auto [ps, vs] = benchProveVerify(1000,
        [&] { nmp.prove(C_I, C_A, w_I, w_A, r_I, r_A); },
        [&] { return ZKNMP_verify(&nmp); });

    g_bench.record("ZKNMP", 8, "prove",  ps);
    g_bench.record("ZKNMP", 8, "verify", vs);
}

void ZKSP_test() {
    cout << "<<< ZKSP TEST >>>" << endl;

    PCS pcs; pcs.setup(8);
    ZKAccSetup setup; setup.setup_from_pcs(pcs);

    FrVec Ax = {24, -50, 35, -10, 1};
    FrVec Ix = {2, -3, 1};
    FrVec Jx = {12, -7, 1};

    G2 C_A = pcs.commit_g2(Ax);
    G2 C_I = pcs.commit_g2(Ix);
    G2 C_J = pcs.commit_g2(Jx);

    G1 pi_I = pcs.commit_g1(PolyLongDiv(Ax, Ix));
    G1 pi_J = pcs.commit_g1(PolyLongDiv(Ax, Jx));

    Fr r_A; r_A.setByCSPRNG();
    Fr r_I; r_I.setByCSPRNG();
    Fr r_J; r_J.setByCSPRNG();

    C_A = C_A + setup.h2 * r_A;
    C_I = C_I + setup.h2 * r_I;
    C_J = C_J + setup.h2 * r_J;

    ZKSP pi3(setup);

    auto [ps, vs] = benchProveVerify(1000,
        [&] { pi3.prove(C_I, C_J, C_A, pi_I, pi_J, r_I, r_J, r_A); },
        [&] { return ZKSP_verify(&pi3); });

    g_bench.record("ZKSP", 8, "prove",  ps);
    g_bench.record("ZKSP", 8, "verify", vs);
}

void SetupTest() {
    cout << "<<< SETUP TEST >>>" << endl;

    uint32_t vals[4] = {5000, 20000, 80000, 320000};
    const int reps = 100;

    for (uint32_t n : vals) {
        std::vector<double> t_srs, t_acc, t_idc, t_wit;
        t_srs.reserve(reps); t_acc.reserve(reps);
        t_idc.reserve(reps); t_wit.reserve(reps);

        for (int j = 0; j < reps; j++) {
            FrVec IDset(n);
            for (uint32_t k = 0; k < n; k++) IDset[k].setByCSPRNG();

            auto a0 = Clock::now();
            PCS pcs; pcs.setup(n + 1);
            ZKAccSetup setup; setup.setup_from_pcs(pcs);
            auto a1 = Clock::now();
            t_srs.push_back(us_since(a0, a1));

            Fr s = pcs.s;

            a0 = Clock::now();
            Fr Ax = 1;
            for (uint32_t k = 0; k < IDset.size(); k++) Ax = Ax * (s + IDset[k]);
            G2 C_A = setup.g2si[0] * Ax;
            a1 = Clock::now();
            t_acc.push_back(us_since(a0, a1));

            a0 = Clock::now();
            Fr id = IDset[rand() % n];
            Fr Ix = s + id;
            G2 C_id = setup.g2 * (s + id);
            a1 = Clock::now();
            t_idc.push_back(us_since(a0, a1));

            a0 = Clock::now();
            G1 wit = setup.g1 * (Ax / Ix);
            a1 = Clock::now();
            t_wit.push_back(us_since(a0, a1));
        }

        g_bench.record("Setup", n, "srs_gen",    t_srs);
        g_bench.record("Setup", n, "acc_commit", t_acc);
        g_bench.record("Setup", n, "id_commit",  t_idc);
        g_bench.record("Setup", n, "wit_gen",    t_wit);
    }
}

void OwnPf_test() {
    cout << "<<< OwnPf TEST >>>" << endl;

    uint32_t vals[4] = {250, 500, 1000, 2000};
    for (auto val : vals) {
        cout << "  set size: " << val << endl;

        PCS pcs; pcs.setup(val);
        ZKAccSetup setup; setup.setup_from_pcs(pcs);
        ZKMP zkmp(setup);
        PoK2 pipok(setup);
        OwnPf proof = {setup, zkmp, pipok};

        FrVec A(val);
        for (uint32_t i = 0; i < val; i++) A[i] = Fr(i + 2);
        FrVec APoly = constructMemPoly(pcs.pp.nttctx, A);

        Fr n = 2; Fr id = 424242; FrVec idPoly = {id, 1};
        G2 C_id = pcs.commit_g2(idPoly);
        G2 C_A  = pcs.commit_g2(APoly);
        Fr r, delta; r.setByCSPRNG(); delta.setByCSPRNG();
        C_id = C_id + setup.h2 * r;
        C_A  = C_A  + setup.h2 * delta;

        auto [ps, vs] = benchProveVerify(1000,
            [&] { proof.prove(C_id, C_A, n, A, id, r, delta); },
            [&] { return OwnPf_verify(&proof); });

        g_bench.record("OwnPf", val, "prove",  ps);
        g_bench.record("OwnPf", val, "verify", vs);
    }
}

void NonOwnPf_test() {
    cout << "<<< NonOwnPf TEST >>>" << endl;

    uint32_t vals[4] = {250, 500, 1000, 2000};
    for (auto val : vals) {
        cout << "  set size: " << val << endl;

        PCS pcs; pcs.setup(val);
        ZKAccSetup setup; setup.setup_from_pcs(pcs);
        ZKNMP zknmp(setup);
        PoK2 pipok(setup);
        NonOwnPf proof = {setup, zknmp, pipok};

        FrVec A(val);
        for (uint32_t i = 0; i < val; i++) A[i] = Fr(i + 2);
        FrVec APoly = constructMemPoly(pcs.pp.nttctx, A);

        Fr n = 1; Fr id = 424242; FrVec idPoly = {id, 1};
        G2 C_id = pcs.commit_g2(idPoly);
        G2 C_A  = pcs.commit_g2(APoly);
        Fr r, delta; r.setByCSPRNG(); delta.setByCSPRNG();
        C_id = C_id + setup.h2 * r;
        C_A  = C_A  + setup.h2 * delta;

        auto [ps, vs] = benchProveVerify(1000,
            [&] { proof.prove(C_id, C_A, n, A, id, r, delta); },
            [&] { return NonOwnPf_verify(&proof); });

        g_bench.record("NonOwnPf", val, "prove",  ps);
        g_bench.record("NonOwnPf", val, "verify", vs);
    }
}

void TxGenTest() {
    cout << "<<< TxGen TEST >>>" << endl;

    vector<uint32_t> numTxs = {256, 512, 1024, 2048};
    for (auto val : numTxs) {
        cout << "  # Txs: " << val << endl;

        PCS pcs; pcs.setup(4243);

        FrVec IDSet(4242);
        for (uint32_t i = 0; i < 4242; i++) IDSet[i].setByCSPRNG();
        Fr id1 = IDSet[42];
        Fr id2 = IDSet[56];
        G1 wit1 = memWitGen(pcs, IDSet, id1);
        G1 wit2 = memWitGen(pcs, IDSet, id2);
        G2 C_A  = pcs.commit_g2(constructMemPoly(pcs.pp.nttctx, IDSet));

        FrVec SN(val);
        for (uint32_t i = 0; i < val; i++) SN[i].setByCSPRNG();
        FrVec S1 = FrVec(SN.begin(), SN.begin() + val / 5);
        FrVec S2 = FrVec(SN.begin() + val / 5, SN.end());

        Fr gamma, delta; gamma.setByCSPRNG(); delta.setByCSPRNG();
        auto a0 = Clock::now();
        Tx_Entry TxEntry = Tx_Entry_Gen(pcs, C_A, SN, id1, wit1, gamma, delta);
        auto a1 = Clock::now();
        bool entry_vrfy = Tx_Entry_Vrfy(TxEntry, pcs);
        auto a2 = Clock::now();
        if (!entry_vrfy) throw std::runtime_error("Tx_Entry verify failed");
        g_bench.value("TxGen", val, "entry_gen",    us_since(a0, a1));
        g_bench.value("TxGen", val, "entry_verify", us_since(a1, a2));

        Fr gamma1, gamma2, delta1, delta2;
        gamma1.setByCSPRNG(); gamma2.setByCSPRNG();
        delta1.setByCSPRNG(); delta2.setByCSPRNG();
        a0 = Clock::now();
        Tx_Trans TxTrans = Tx_Trans_Gen(
            TxEntry.Tx_out, pcs, C_A, SN, S1, S2, id1, id2, wit1, wit2,
            gamma, delta, gamma1, gamma2, delta1, delta2);
        a1 = Clock::now();
        bool trans_vrfy = Tx_Trans_Vrfy(TxTrans, pcs);
        a2 = Clock::now();
        if (!trans_vrfy) throw std::runtime_error("Tx_Trans verify failed");
        g_bench.value("TxGen", val, "transfer_gen",    us_since(a0, a1));
        g_bench.value("TxGen", val, "transfer_verify", us_since(a1, a2));

        Fr gammap; gammap.setByCSPRNG();
        a0 = Clock::now();
        Tx_Exit TxExit = Tx_Exit_Gen(
            TxEntry.Tx_out, pcs, C_A, SN, S1, S2, id1, wit1, gamma, delta,
            gammap, delta1, delta2);
        a1 = Clock::now();
        bool exit_vrfy = Tx_Exit_Vrfy(TxExit, pcs);
        a2 = Clock::now();
        if (!exit_vrfy) throw std::runtime_error("Tx_Exit verify failed");
        g_bench.value("TxGen", val, "exit_gen",    us_since(a0, a1));
        g_bench.value("TxGen", val, "exit_verify", us_since(a1, a2));
    }
}

void TxSummaryTest() {
    cout << "<<< TxSummary TEST >>>" << endl;

    vector<uint32_t> numTxs = {128, 256, 512, 1024, 2048};
    for (auto val : numTxs) {
        cout << "  # Txs*2: " << val << endl;

        PCS pcs; pcs.setup(4243);

        FrVec IDSet(4242);
        for (uint32_t i = 0; i < 4242; i++) IDSet[i].setByCSPRNG();
        G2 C_A = pcs.commit_g2(constructMemPoly(pcs.pp.nttctx, IDSet));

        FrVec SN_loc(val); vector<G2> C_Is(val);
        vector<G1> W_Is(val); FrVec gamma_Is(val);

        for (uint32_t i = 0; i < val; i++) {
            SN_loc[i] = IDSet[i];
            gamma_Is[i].setByCSPRNG();
            C_Is[i] = pcs.commit_g2({SN_loc[i], 1}) + pcs.pp.h2si[0] * gamma_Is[i];
        }
        for (uint32_t i = 0; i < val; i++)
            W_Is[i] = memWitGen(pcs, SN_loc, SN_loc[i]);

        Fr gamma_loc; gamma_loc.setByCSPRNG();
        G2 C_loc = pcs.commit_g2(constructMemPoly(pcs.pp.nttctx, SN_loc))
                   + pcs.pp.h2si[0] * gamma_loc;
        G1 W_loc = pcs.commit_g1(constructMemPolyBatch(pcs.pp.nttctx, IDSet, SN_loc));

        Fr Adr(42);
        Tx_Base Tx_in = {Adr, C_loc, C_loc};

        Fr gamma_locp; gamma_locp.setByCSPRNG();
        auto a0 = Clock::now();
        Tx_Summary TxSum = Tx_Summary_Gen(
            Tx_in, pcs, SN_loc, C_A, C_Is,
            gamma_loc, W_loc, W_Is, gamma_Is, gamma_locp);
        auto a1 = Clock::now();
        bool ret_sum = Tx_Summary_Vrfy(TxSum, pcs);
        auto a2 = Clock::now();
        if (!ret_sum) throw std::runtime_error("Tx_Summary verify failed");

        g_bench.value("TxSummary", val, "summary_gen",    us_since(a0, a1));
        g_bench.value("TxSummary", val, "summary_verify", us_since(a1, a2));
    }
}

static void usage(const char* prog) {
    std::cerr << "usage: " << prog
              << " [--mode local|server|all] [--verbose|--no-verbose]\n"
              << "  --mode local   ZKMP, AZKMP, ZKSP, ZKNMP, OwnPf, NonOwnPf, TxGen\n"
              << "  --mode server  TxSummary, Setup\n"
              << "  --mode all     everything (default)\n";
}

int main(int argc, char** argv) {
    initPairing(mcl::BLS12_381);

    std::string mode = "all";
    bool verbose = true;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--mode" && i + 1 < argc) {
            mode = argv[++i];
        } else if (a.rfind("--mode=", 0) == 0) {
            mode = a.substr(7);
        } else if (a == "--verbose") {
            verbose = true;
        } else if (a == "--no-verbose") {
            verbose = false;
        } else if (a.rfind("--verbose=", 0) == 0) {
            std::string v = a.substr(10);
            verbose = (v == "1" || v == "true" || v == "on" || v == "yes");
        } else if (a == "-h" || a == "--help") {
            usage(argv[0]);
            return 0;
        } else {
            std::cerr << "unknown arg: " << a << "\n";
            usage(argv[0]);
            return 1;
        }
    }

    g_bench.verbose = verbose;

    auto run_local = [] {
        ZKMP_test();
        ZKSP_test();
        ZKNMP_test();
        OwnPf_test();
        NonOwnPf_test();
        TxGenTest();
    };
    auto run_server = [] {
        AZKMP_test();
        TxSummaryTest();
        SetupTest();
    };

    if (mode == "local") {
        run_local();
    } else if (mode == "server") {
        run_server();
    } else if (mode == "all") {
        run_local();
        run_server();
    } else {
        std::cerr << "unknown mode: " << mode << " (expected local|server|all)\n";
        usage(argv[0]);
        return 1;
    }

    g_bench.dumpConsole();
    g_bench.dumpCsv("bench_" + mode + ".csv");
    return 0;
}