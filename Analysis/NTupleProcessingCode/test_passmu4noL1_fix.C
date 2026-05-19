// Test: verify passmu4noL1 fallback logic fix.
// Run with: root -l -b -q test_passmu4noL1_fix.C
#include "DataAnalysisClasses.h"

void test_passmu4noL1_fix() {
    // PbPb 2023, batch 1, 500 events
    PbPbAnalysis pbpb(23, 1);
    pbpb.nevents_max = 500;
    pbpb.is_test_run = true;
    pbpb.Run();

    std::string outpath = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_2023/muon_pairs_pbpb_2023_part1_single_mu4_mindR_0_02_res_cut_v2_test.root";
    std::cout << "\n=== TEST: reading back " << outpath << " ===" << std::endl;

    auto f = TFile::Open(outpath.c_str());
    if (!f || f->IsZombie()) { std::cerr << "Cannot open output!" << std::endl; return; }

    for (const char* tname : {"muon_pair_tree_sign1", "muon_pair_tree_sign2"}) {
        auto t = (TTree*)f->Get(tname);
        if (!t) { std::cout << tname << ": not found" << std::endl; continue; }

        MuonPairPbPb* p = nullptr;
        t->SetBranchAddress("MuonPairObj", &p);

        int N = t->GetEntries();
        int n_m1mu4=0, n_m2mu4=0, n_m1noL1=0, n_m2noL1=0;
        int n_passmu4mu4noL1=0, n_passmu4noL1=0, n_pass2mu4=0;
        int n_both_mu4_and_noL1 = 0;  // m1 passes both mu4 AND noL1 (should not happen frequently)

        for (int i = 0; i < N; i++) {
            t->GetEntry(i);
            if (p->m1.passmu4) n_m1mu4++;
            if (p->m2.passmu4) n_m2mu4++;
            if (p->m1.passmu4noL1) n_m1noL1++;
            if (p->m2.passmu4noL1) n_m2noL1++;
            if (p->passmu4mu4noL1) n_passmu4mu4noL1++;
            if (p->passmu4noL1) n_passmu4noL1++;
            if (p->pass2mu4) n_pass2mu4++;

            // Consistency: if m1.passmu4noL1 is true, the pair should pass mu4_mu4noL1
            if (p->m1.passmu4noL1 && !p->passmu4mu4noL1)
                std::cerr << "INCONSISTENCY: m1.passmu4noL1=true but passmu4mu4noL1=false at entry " << i << std::endl;
            if (p->m2.passmu4noL1 && !p->passmu4mu4noL1)
                std::cerr << "INCONSISTENCY: m2.passmu4noL1=true but passmu4mu4noL1=false at entry " << i << std::endl;

            // Consistency: if m1.passmu4noL1 is true, m2 should pass mu4
            if (p->m1.passmu4noL1 && !p->m2.passmu4)
                std::cerr << "INCONSISTENCY: m1.passmu4noL1=true but m2.passmu4=false at entry " << i << std::endl;
            if (p->m2.passmu4noL1 && !p->m1.passmu4)
                std::cerr << "INCONSISTENCY: m2.passmu4noL1=true but m1.passmu4=false at entry " << i << std::endl;
        }

        std::cout << tname << " (N=" << N << "):" << std::endl;
        std::cout << "  m1.passmu4=" << n_m1mu4 << " (" << (N>0 ? 100.*n_m1mu4/N : 0) << "%)" << std::endl;
        std::cout << "  m2.passmu4=" << n_m2mu4 << " (" << (N>0 ? 100.*n_m2mu4/N : 0) << "%)" << std::endl;
        std::cout << "  m1.passmu4noL1=" << n_m1noL1 << " (" << (N>0 ? 100.*n_m1noL1/N : 0) << "%)" << std::endl;
        std::cout << "  m2.passmu4noL1=" << n_m2noL1 << " (" << (N>0 ? 100.*n_m2noL1/N : 0) << "%)" << std::endl;
        std::cout << "  passmu4mu4noL1=" << n_passmu4mu4noL1 << " (" << (N>0 ? 100.*n_passmu4mu4noL1/N : 0) << "%)" << std::endl;
        std::cout << "  passmu4noL1=" << n_passmu4noL1 << " (" << (N>0 ? 100.*n_passmu4noL1/N : 0) << "%)" << std::endl;
        std::cout << "  pass2mu4=" << n_pass2mu4 << " (" << (N>0 ? 100.*n_pass2mu4/N : 0) << "%)" << std::endl;

        // KEY CHECK: passmu4noL1 should no longer be all-zero
        if (N > 0 && n_m1noL1 == 0 && n_m2noL1 == 0)
            std::cerr << "  *** FAIL: passmu4noL1 is still all-zero! ***" << std::endl;
        else if (N > 0)
            std::cout << "  PASS: passmu4noL1 has non-zero entries" << std::endl;
    }

    f->Close();
}
