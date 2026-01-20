#pragma once
#include "../MuonObjectsParamsAndHelpers/TruthQQPair.h"
#include "../MuonObjectsParamsAndHelpers/struct_particle.h"
#include "../MuonObjectsParamsAndHelpers/muon_pair_enums_MC.h"

template <class Derived>
class PowhegTruthExtras{
protected:
    Derived& self() { return static_cast<Derived&>(*this); }
    const Derived& self() const { return static_cast<const Derived&>(*this); }

    auto& mc_mode() { return self().mc_mode; }
    auto& pms()     { return self().pms; }
    auto& fChain()   { return self().fChain; }
    auto& mpair()   { return self().mpair; }
    auto& mcdir()   { return self().mcdir; }
    auto& EventWeights()   { return self().EventWeights; }

// --------------------- input files & trees & data for setting branches---------------------------

    std::vector<int>* truth_id = nullptr;
    std::vector<int>* truth_barcode = nullptr;
    std::vector<int>* truth_qual = nullptr;
    std::vector<float>* truth_pt = nullptr;
    std::vector<float>* truth_eta = nullptr;
    std::vector<float>* truth_phi = nullptr;
    std::vector<float>* truth_m = nullptr;
    std::vector<std::vector<int>>* truth_parents = nullptr;
    std::vector<std::vector<int>>* truth_children = nullptr;
    // std::vector<std::vector<int>>* truth_muon_parent_ids = nullptr;
    // std::vector<std::vector<int>>* truth_muon_parent_bars = nullptr;

// --------------------- truth origin parameters ---------------------------

    static const int nAncestorGroups = 4;

    // configuration
    std::string sign_labels[ParamsSet::nSigns] = {"same_sign", "op_sign"};
    std::string signs[ParamsSet::nSigns] = {"_ss", "_op"};
    std::string dphis[2] = {"_near", "_away"};
    std::string ancestor_grps[nAncestorGroups] = {"_gg", "_qg","_single_g","_qq"};

    // std::vector<std::string> parentGroupLabels = {"direct b","b to c","direct c","s/light","single photon", "Drell-Yan"};
    std::vector<std::string> parentGroupLabels = {"direct b","b to c","direct c","s/light","single photon"};
    int nParentGroups = parentGroupLabels.size();    
    
    std::vector<std::string> ancestor_labels = {"gg", "gq", "single g", "q qbar"};
    
    std::vector<std::string> samePrtsLabels = {"Same Parents", "Different Parents"};
    std::vector<std::string> bb_op_one_b_one_btoc_labels = {"Same b", "Involve osc(s)", "From different ancestors", "Others"};

    // std::vector<std::string> osc_labels = {"0 osc, one b one b-to-c", "0 osc, others", "1 osc, one b one b-to-c", "1 osc, regular", "2 oscs, one b one b-to-c", "2 oscs, others"};
    std::vector<std::string> osc_labels = {"0 osc, one c-tag", "0 osc, others(*)", "1 osc, one c-tag(*)", "1 osc, regular", "2 oscs, one c-tag", "2 oscs, others(*)"};
    std::vector<std::string> num_hard_scatt_out_labels = {"2","3","more"};

// --------------------- temporary variables (muon, muonpair objects, vectors, etc.) ---------------------------

    // std::vector<std::shared_ptr<MuonPairPowheg>> muon_pair_list_cur_event;
    TruthQQPair* qqpair = nullptr;

    int cur_m1_earliest_parent_barcode;
    int cur_m2_earliest_parent_barcode;

    // bool same_ancestors;
    bool m1_ancestor_is_incoming;
    bool m2_ancestor_is_incoming;

    bool m1_from_J_psi;
    bool m2_from_J_psi;

    bool m1_from_Upsilon;
    bool m2_from_Upsilon;

    float crossx_2_to_2 = 0.;
    float crossx_2_to_3 = 0.;
    float crossx_relevant_hard_isnt_hardest = 0.;

    double skipped_event_crossx = 0;
    bool skip_event = false;
    // boolean to tag whether 
    bool pre_return = false;
    
    std::vector<std::vector<int>>* single_gluon_history;
    std::vector<std::vector<int>>* m1_history;
    std::vector<std::vector<int>>* m2_history;
    std::vector<std::vector<Particle>>* m1_history_particle;
    std::vector<std::vector<Particle>>* m2_history_particle;


    std::vector<int> cur_m1_ancestor_ids;
    std::vector<int> cur_m2_ancestor_ids;
    std::vector<int> cur_m1_ancestor_bars;
    std::vector<int> cur_m2_ancestor_bars;

    std::vector<int> m1_multi_hf_quark_ids;
    std::vector<int> m2_multi_hf_quark_ids;

    bool m1_c_tag;
    bool m2_c_tag;
    bool m1_osc;
    bool m2_osc;
    int m1_earliest_parent_id;
    int m2_earliest_parent_id;

    // std::vector<float>* m1_last_b_hadron_prt_pt_eta_phi_m;
    // std::vector<float>* m2_last_b_hadron_prt_pt_eta_phi_m;
    // std::vector<float>* m1_last_hf_hadron_prt_pt_eta_phi_m;
    // std::vector<float>* m2_last_hf_hadron_prt_pt_eta_phi_m;
    std::vector<float>* m1_first_hf_hadron_prt_pt_eta_phi_m;
    std::vector<float>* m2_first_hf_hadron_prt_pt_eta_phi_m;
    // std::vector<float>* m1_hq_ancestor_pt_eta_phi_m;
    // std::vector<float>* m2_hq_ancestor_pt_eta_phi_m;

// --------------------- output file, histograms & trees ---------------------------
  
    std::ofstream* m_unspecified_parent_file = nullptr;
    std::ofstream* m_b_parent_file[ParamsSet::nSigns][2];
    std::ofstream* m_c_parent_file[ParamsSet::nSigns][2];
    std::ofstream* m_cc_ss_small_dphi_file = nullptr;
    std::ofstream* m_bb_ss_near_file = nullptr;
    std::ofstream* m_bb_ss_away_file = nullptr;
    std::ofstream* m_bb_op_near_one_b_one_btoc_others_file = nullptr;

    TH2D* h_parent_groups[ParamsSet::nSigns][2];

    TH1D* h_QQ_DR[ParamsSet::nSigns][2][4];
    TH1D* h_QQ_Dphi[ParamsSet::nSigns][2][4];
    TH1D* h_QQ_minv[ParamsSet::nSigns][2][4];
    TH1D* h_QQ_pair_pt_ptlead_ratio[ParamsSet::nSigns][2][4];
    TH1D* h_QQ_pt_avg[ParamsSet::nSigns][2][4];
    TH1D* h_QQ_asym[ParamsSet::nSigns][2][4];
    TH1D* h_QQ_minv_mHard_ratio[ParamsSet::nSigns][2][4];

    TH2D* h_QQ_ptlead_pair_pt[ParamsSet::nSigns][2][4];
    TH2D* h_QQ_pt1_pt2[ParamsSet::nSigns][2][4];
    TH2D* h_QQ_Deta_Dphi[ParamsSet::nSigns][2][4];
    TH2D* h_QQ_eta1_eta2[ParamsSet::nSigns][2][4];
    TH2D* h_QQ_minv_pair_pt[ParamsSet::nSigns][2][4];
    TH2D* h_QQ_minv_Dphi[ParamsSet::nSigns][2][4];

    static const int npTbins = 3;
    TH1D* h_crossx;
    TH1D* h_crossx_pt_binned[npTbins];

    // TH1D* h_bb_op_one_b_one_btoc[2];

    TH1D* h_QQ_both_from_Q_same_ancestors[ParamsSet::nSigns][2];
    TH1D* h_QQ_both_from_Q_ancestor_sp[ParamsSet::nSigns][2];
    TH2D* h_QQ_both_from_Q_ancestor_dp[ParamsSet::nSigns][2];

    TH1D* h_num_hard_scatt_out[ParamsSet::nSigns][2];
    TH1D* h_pt_muon_pt_closest_hadr_ratio[ParamsSet::nSigns][2];
    TH1D* h_pt_closest_hadr_pt_furthest_hadr_ratio[ParamsSet::nSigns][2];
    TH1D* h_pt_hadr_hq_ratio[ParamsSet::nSigns][2];
    TH1D* h_dphi_muon_closest_hadr[ParamsSet::nSigns][2];

    TTree* QQPairOutTree[ParamsSet::nSigns][2][nAncestorGroups];

    void MuonPairAncestorTracing();
    void SingleMuonAncestorTracing(bool isMuon1);
    int  ParentGrouping(std::vector<int>& parent_ids, bool c_tag, bool prev_is_lepton);
    void GetPtEtaPhiMFromBarcode(int barcode, std::vector<float>* pt_eta_phi_m);
    int  UpdateCurParents(bool isMuon1, std::vector<int>& cur_prt_bars, std::vector<int>& cur_prt_ids, int hf_quark_index = -1);
    int  FindHeavyQuarks(std::vector<int>& cur_prt_ids, int quark_type, bool isMuon1, int hadron_child_id = 0);
    int  AncestorGrouping(std::vector<int>& ancestor_ids);
    void HardScatteringAnalysis(std::vector<int>& ancestor_bars, std::vector<int>& ancestor_ids, int sign_dphi_mode, int ancestor_grp);
    void PrintHistory(std::ostream* f, bool print_single, bool muon1_sameancestor);
    void KinematicCorrPlots(int isign, int jdphi);
    void MuonPairTagsReinit();
    void CheckIfFromSameB();

    void InitializeExtra();
    void InitParamsExtra(){
        self().perform_truth = true;
    }
    void InitInputExtra();
    void InitTempVariablesExtra();

    void InitOutStreamFiles();
    void InitOutputExtra();
    void InitOutputTreesExtra();
    void InitOutputHistsExtra();
    void HistAdjustExtra();
    void TruthPairAnalysis();
    
    void FinalizeExtra();
public: 
    bool output_QQpair_tree = true;
    bool output_truth_hists = true;

    bool print_prt_history = false;
    bool print_specific_prt_history = false;

    ~PowhegTruthExtras(){}
};
