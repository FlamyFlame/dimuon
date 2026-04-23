// test_check_skim_branches.cxx
// Diagnostic: check branch types in source skim and values in pair ntuple.
// Run: root -b -q test_check_skim_branches.cxx

void test_check_skim_branches() {
    const char* skim = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_2024/data_pbpb24_part1.root";
    const char* pairs = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_2024/muon_pairs_pbpb_2024_part1_single_mu4_mindR_0_02_res_cut_v2.root";

    // ---- 1. Print branch types for ZDC/FCal in source skim ----
    TFile* fs = TFile::Open(skim, "READ");
    TTree* ts = (TTree*)fs->Get("HeavyIonD3PD");
    std::cout << "\n=== Source skim: TTree::Print() for ZDC/FCal/trk branches ===" << std::endl;
    ts->Print("zdc_ZdcEnergy");
    ts->Print("zdc_ZdcTime");
    ts->Print("FCal_Et_P");
    ts->Print("FCal_Et_N");
    ts->Print("FCal_Et");
    ts->Print("trk_numqual");

    std::cout << "\n=== Leaf type summary ===" << std::endl;
    const char* branches[] = {
        "FCal_Et", "FCal_Et_P", "FCal_Et_N",
        "zdc_ZdcEnergy", "zdc_ZdcTime", "zdc_ZdcModulePreSampleAmp",
        "trk_numqual", nullptr
    };
    for (int i = 0; branches[i] != nullptr; i++) {
        TBranch* b = ts->GetBranch(branches[i]);
        if (!b) { std::cout << branches[i] << ": NOT FOUND" << std::endl; continue; }
        // Print leaf type info
        TLeaf* leaf = (TLeaf*)b->GetListOfLeaves()->At(0);
        if (leaf)
            std::cout << branches[i] << ": " << b->GetClassName()
                      << "  leaf type=" << leaf->GetTypeName()
                      << "  len=" << leaf->GetLen()
                      << "  IsRange=" << leaf->IsRange() << std::endl;
        else
            std::cout << branches[i] << ": " << b->GetClassName() << " (no leaves)" << std::endl;
    }

    // ---- 2. Try reading zdc_ZdcEnergy both ways and compare ----
    std::cout << "\n=== Try reading zdc_ZdcEnergy as Float_t[2] ===" << std::endl;
    Float_t zdc_arr[2] = {-999.f, -999.f};
    int ret1 = ts->SetBranchAddress("zdc_ZdcEnergy", zdc_arr);
    std::cout << "SetBranchAddress(Float_t[2]) return code = " << ret1 << "  (-1=mismatch, 0=ok, >0=ok)" << std::endl;
    ts->GetEntry(0);
    std::cout << "  event 0: zdc_arr[0]=" << zdc_arr[0] << " zdc_arr[1]=" << zdc_arr[1] << std::endl;

    std::cout << "\n=== Try reading zdc_ZdcEnergy as vector<float>* ===" << std::endl;
    ts->ResetBranchAddresses();
    std::vector<float>* zdc_vec = nullptr;
    int ret2 = ts->SetBranchAddress("zdc_ZdcEnergy", &zdc_vec);
    std::cout << "SetBranchAddress(vector<float>*) return code = " << ret2 << std::endl;
    ts->GetEntry(0);
    if (zdc_vec) {
        std::cout << "  event 0: zdc_vec.size()=" << zdc_vec->size();
        for (size_t k = 0; k < zdc_vec->size() && k < 4; k++)
            std::cout << "  [" << k << "]=" << zdc_vec->at(k);
        std::cout << std::endl;
    } else std::cout << "  zdc_vec is null" << std::endl;

    // ---- 3. Check FCal_Et_P as Float_t scalar ----
    ts->ResetBranchAddresses();
    Float_t fcal_p = -999.f, fcal_n = -999.f;
    int r3 = ts->SetBranchAddress("FCal_Et_P", &fcal_p);
    int r4 = ts->SetBranchAddress("FCal_Et_N", &fcal_n);
    ts->GetEntry(0);
    std::cout << "\n=== FCal_Et_P/N (ret=" << r3 << "/" << r4 << ") ===" << std::endl;
    std::cout << "  event 0: FCal_Et_P=" << fcal_p << " FCal_Et_N=" << fcal_n << " TeV" << std::endl;

    fs->Close();

    // ---- 4. Check values in pair ntuple ----
    std::cout << "\n=== Pair ntuple: first 5 events ===" << std::endl;
    TFile* fp = TFile::Open(pairs, "READ");
    // sign1 = first tree
    TTree* tp = nullptr;
    fp->GetObject("muon_pair_tree_sign1", tp);
    if (!tp) fp->GetObject("muon_pair_tree_sign2", tp);
    if (!tp) { std::cout << "No muon_pair_tree found!" << std::endl; fp->Close(); return; }
    std::cout << "Tree: " << tp->GetName() << "  entries=" << tp->GetEntries() << std::endl;

    // Print branch list for MuonPairObj
    std::cout << "MuonPairObj sub-branches:" << std::endl;
    TBranch* bpair = tp->GetBranch("MuonPairObj");
    if (bpair) {
        TObjArray* subbranches = bpair->GetListOfBranches();
        for (int k = 0; k < subbranches->GetEntries(); k++) {
            TBranch* sub = (TBranch*)subbranches->At(k);
            std::cout << "  " << sub->GetName() << std::endl;
        }
    }

    // Scan the columns via TTree::Scan (no dictionary needed)
    std::cout << "\nScan FCal_Et:FCal_Et_A:FCal_Et_C:ZDC_E_tot:ntrk_HIloose (first 5):" << std::endl;
    tp->Scan("FCal_Et:FCal_Et_A:FCal_Et_C:ZDC_E_tot:ntrk_HIloose", "", "", 5);

    fp->Close();
}
