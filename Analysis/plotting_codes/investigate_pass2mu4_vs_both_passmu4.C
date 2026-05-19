// H6 test: compare pass2mu4 (pair-level) vs m1.passmu4 && m2.passmu4 (individual)
// If pass2mu4 is more restrictive, it explains why eps_nc < P(probe fires mu4)
// and why the inverse-weighted ratio > 1.0

#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TROOT.h"
#include "/usatlas/u/yuhanguo/workarea/dimuon_codes/Analysis/MuonObjectsParamsAndHelpers/MuonPairPbPb.h"
#include <cstdio>
#include <algorithm>

R__LOAD_LIBRARY(/usatlas/u/yuhanguo/workarea/dimuon_codes/Analysis/MuonObjectsParamsAndHelpers/MuonPairPbPb_h.so)

void investigate_pass2mu4_vs_both_passmu4(int year = 25) {
    std::string yr = std::to_string(year);
    std::string base = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_20" + yr + "/";

    std::vector<std::string> candidates = {
        base + "muon_pairs_pbpb_20" + yr + "_single_mu4_mindR_0_02_res_cut_v2.root",
        base + "muon_pairs_pbpb_20" + yr + "_single_mu4_mindR_0_02.root",
        base + "muon_pairs_pbpb_20" + yr + "_single_mu4.root"
    };

    TFile* f = nullptr;
    for (auto& c : candidates) {
        if (!gSystem->AccessPathName(c.c_str())) {
            f = TFile::Open(c.c_str());
            if (f && !f->IsZombie()) { printf("Opened: %s\n", c.c_str()); break; }
        }
    }
    if (!f) { printf("Cannot open ntuple for year %d\n", year); return; }

    for (auto treename : {"muon_pair_tree_sign2", "muon_pair_tree_sign1"}) {
        auto t = (TTree*)f->Get(treename);
        if (!t) { printf("%s: not found\n", treename); continue; }

        MuonPairPbPb* mp = nullptr;
        t->SetBranchAddress("mpair", &mp);

        long long n_total = 0, n_both = 0, n_2mu4 = 0, n_both_not_2mu4 = 0, n_2mu4_not_both = 0;
        long long n_tag1_total = 0, n_tag1_probe_pass = 0, n_tag1_2mu4 = 0;
        long long nev = t->GetEntries();

        for (long long i = 0; i < nev; i++) {
            t->GetEntry(i);
            if (!mp->passSeparated) continue;
            n_total++;

            bool b1 = mp->m1.passmu4;
            bool b2 = mp->m2.passmu4;
            bool both = b1 && b2;
            bool pair = mp->pass2mu4;

            if (both) n_both++;
            if (pair) n_2mu4++;
            if (both && !pair) n_both_not_2mu4++;
            if (pair && !both) n_2mu4_not_both++;

            // Tag = m1, probe = m2
            if (b1) {
                n_tag1_total++;
                if (b2) n_tag1_probe_pass++;
                if (pair) n_tag1_2mu4++;
            }
        }

        printf("\n=== %s (dR > 0.8, %lld entries) ===\n", treename, nev);
        printf("Pairs with dR > 0.8: %lld\n", n_total);
        printf("m1.passmu4 && m2.passmu4: %lld (%.2f%%)\n", n_both, 100.0*n_both/n_total);
        printf("pass2mu4:                 %lld (%.2f%%)\n", n_2mu4, 100.0*n_2mu4/n_total);
        printf("Both BUT NOT pass2mu4:    %lld (%.2f%% of those with both)\n", n_both_not_2mu4, n_both>0 ? 100.0*n_both_not_2mu4/n_both : 0);
        printf("pass2mu4 BUT NOT both:    %lld (%.2f%% of pass2mu4)\n", n_2mu4_not_both, n_2mu4>0 ? 100.0*n_2mu4_not_both/n_2mu4 : 0);
        printf("\nTag-probe (m1=tag):\n");
        printf("  Tag fires mu4:         %lld\n", n_tag1_total);
        printf("  Probe fires mu4:       %lld (%.4f)\n", n_tag1_probe_pass, (double)n_tag1_probe_pass/n_tag1_total);
        printf("  Pair fires 2mu4:       %lld (%.4f)\n", n_tag1_2mu4, (double)n_tag1_2mu4/n_tag1_total);
        printf("  Ratio probe/2mu4:      %.4f\n", (double)n_tag1_probe_pass/n_tag1_2mu4);
        printf("  --> If this ratio > 1, pass2mu4 is more restrictive than individual matching\n");
    }
    f->Close();
}
