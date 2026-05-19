// Compare mu4/2mu4 flags between old and new outputs to verify the passmu4noL1 fix
// doesn't affect mu4 or 2mu4 assignments.
// Run with: root -l -b -q compare_mu4_flags.C
#include "DataAnalysisClasses.h"

void compare_mu4_flags() {
    auto f_old = TFile::Open("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_2023/muon_pairs_pbpb_2023_part1_single_mu4_mindR_0_02_res_cut_v2.root");
    auto f_new = TFile::Open("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_2023/muon_pairs_pbpb_2023_part1_single_mu4_mindR_0_02_res_cut_v2_test.root");
    if (!f_old) { printf("Old file not found\n"); return; }
    if (!f_new) { printf("New file not found\n"); return; }

    for (const char* tname : {"muon_pair_tree_sign1", "muon_pair_tree_sign2"}) {
        auto t_old = (TTree*)f_old->Get(tname);
        auto t_new = (TTree*)f_new->Get(tname);
        if (!t_old || !t_new) { printf("%s: missing in one file\n", tname); continue; }

        MuonPairPbPb *p_old = nullptr, *p_new = nullptr;
        t_old->SetBranchAddress("MuonPairObj", &p_old);
        t_new->SetBranchAddress("MuonPairObj", &p_new);

        int N_old = t_old->GetEntries();
        int N_new = t_new->GetEntries();
        printf("%s: old=%d new=%d\n", tname, N_old, N_new);

        int N = std::min(N_old, N_new);
        int n_diff_mu4_m1=0, n_diff_mu4_m2=0, n_diff_2mu4=0, n_diff_mu4mu4noL1=0;
        for (int i = 0; i < N; i++) {
            t_old->GetEntry(i);
            t_new->GetEntry(i);
            if (p_old->m1.passmu4 != p_new->m1.passmu4) n_diff_mu4_m1++;
            if (p_old->m2.passmu4 != p_new->m2.passmu4) n_diff_mu4_m2++;
            if (p_old->pass2mu4 != p_new->pass2mu4) n_diff_2mu4++;
            if (p_old->passmu4mu4noL1 != p_new->passmu4mu4noL1) n_diff_mu4mu4noL1++;
        }
        printf("  diff m1.passmu4: %d\n", n_diff_mu4_m1);
        printf("  diff m2.passmu4: %d\n", n_diff_mu4_m2);
        printf("  diff pass2mu4: %d\n", n_diff_2mu4);
        printf("  diff passmu4mu4noL1: %d\n", n_diff_mu4mu4noL1);

        int old_noL1 = 0;
        for (int i = 0; i < N; i++) {
            t_old->GetEntry(i);
            if (p_old->m1.passmu4noL1 || p_old->m2.passmu4noL1) old_noL1++;
        }
        printf("  old passmu4noL1 non-zero entries: %d (confirms bug was present)\n", old_noL1);
    }
    f_old->Close();
    f_new->Close();
}
