void check_test_ntuple() {
    const char* f_path = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_2024/muon_pairs_pbpb_2024_part1_single_mu4_mindR_0_02_test.root";
    TFile* f = TFile::Open(f_path, "READ");
    if (!f || f->IsZombie()) { std::cout << "Cannot open: " << f_path << std::endl; return; }

    for (const char* tname : {"muon_pair_tree_sign1", "muon_pair_tree_sign2"}) {
        TTree* t = (TTree*)f->Get(tname);
        if (!t) { std::cout << tname << ": not found" << std::endl; continue; }
        std::cout << tname << " entries=" << t->GetEntries() << std::endl;
        if (t->GetEntries() == 0) continue;

        t->SetMakeClass(1);
        Float_t fcal_tot=-999, fcal_a=-999, fcal_c=-999, zdc_tot=-999;
        Int_t avg_cen=-999, ntrk_loose=-999, ntrk_tight=-999;
        t->SetBranchAddress("MuonPairObj.FCal_Et",       &fcal_tot);
        t->SetBranchAddress("MuonPairObj.FCal_Et_A",     &fcal_a);
        t->SetBranchAddress("MuonPairObj.FCal_Et_C",     &fcal_c);
        t->SetBranchAddress("MuonPairObj.ZDC_E_tot",     &zdc_tot);
        t->SetBranchAddress("MuonPairObj.avg_centrality",&avg_cen);
        t->SetBranchAddress("MuonPairObj.ntrk_HIloose",  &ntrk_loose);
        t->SetBranchAddress("MuonPairObj.ntrk_HItight",  &ntrk_tight);

        int nprint = std::min((Long64_t)5, t->GetEntries());
        std::cout << "  [entry] FCal_Et_TeV  FCal_A_TeV  FCal_C_TeV  ZDC_GeV   avg_cen  ntrk_loose  ntrk_tight" << std::endl;
        for (int i = 0; i < nprint; i++) {
            t->GetEntry(i);
            std::cout << "  [" << i << "]"
                      << "  " << fcal_tot   // stored in TeV after fix
                      << "  " << fcal_a
                      << "  " << fcal_c
                      << "  " << zdc_tot    // stored in GeV
                      << "  " << avg_cen
                      << "  " << ntrk_loose
                      << "  " << ntrk_tight << std::endl;
        }
    }
    f->Close();
}
