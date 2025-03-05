#ifndef __PowhegNTupleFirstPass_h__
#define __PowhegNTupleFirstPass_h__

#include "../MuonObjectsParamsAndHelpers/MuonPairPowheg.h"
#include "../MuonObjectsParamsAndHelpers/TruthQQPair.h"
#include "../MuonObjectsParamsAndHelpers/muon_pair_enums_MC.h"
#include "../MuonObjectsParamsAndHelpers/struct_particle.h"
#include "MuonNTupleFirstPassBaseClass.c"


class PowhegNTupleFirstPass : public MuonNTupleFirstPassBaseClass{
    // Read through the N-tuple, apply appropriate cuts
    // Then fill in histograms and/or output trees
    // mode = 1: output single-muon information into a TTree
    // mode = 2: output muon-pair information into a TTree
    // NOW ONLY HAVE MODE = 1, 2


private:
// --------------------- general settings ---------------------------

    double crossx_cut;
    double filter_effcy;
    double filter_effcy_bb = 0.003;
    double filter_effcy_cc = 0.001108;
    // static const int nMCmodes = 2;

    ParamsSet pms;

    std::string mcdir;

// --------------------- input files & trees & data for setting branches---------------------------
  
    TChain          *fChain;   //!pointer to the analyzed TTree or TChain
  
    // Declaration of leaf types
    // UInt_t          RunNumber;
    float                 Q;
    std::vector<float>   *EventWeights           =nullptr;

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

    std::vector<float>   *muon_pair_muon1_pt           =nullptr;
    std::vector<float>   *muon_pair_muon1_eta       =nullptr;
    std::vector<float>   *muon_pair_muon1_phi          =nullptr;
    std::vector<int>     *muon_pair_muon1_ch         =nullptr;
    std::vector<int>     *muon_pair_muon1_bar         =nullptr;

    std::vector<float>   *muon_pair_muon2_pt       =nullptr;
    std::vector<float>   *muon_pair_muon2_eta          =nullptr;
    std::vector<float>   *muon_pair_muon2_phi         =nullptr;
    std::vector<int>     *muon_pair_muon2_ch         =nullptr;
    std::vector<int>     *muon_pair_muon2_bar         =nullptr;

    std::vector<float>   *truth_mupair_asym         =nullptr;
    std::vector<float>   *truth_mupair_acop         =nullptr;
    std::vector<float>   *truth_mupair_pt         =nullptr;
    std::vector<float>   *truth_mupair_y         =nullptr;
    // std::vector<float>   *truth_mupair_phi         =nullptr;
    std::vector<float>   *truth_mupair_m         =nullptr;
    
// --------------------- temporary variables (muon, muonpair objects, vectors, etc.) ---------------------------
  
    Long64_t nentries;
    std::string sub_dir;

    Muon* tempmuon = nullptr;
    std::shared_ptr<MuonPairPowheg> mpair;
    MuonPairPowheg* mpair_raw_ptr = nullptr;

    // std::vector<std::shared_ptr<MuonPairPowheg>> muon_pair_list_cur_event;
    std::vector<std::shared_ptr<MuonPairPowheg>> muon_pair_list_cur_event_pre_resonance_cut;
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

    int nmuonpairs;

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

    TTree* meta_tree;
    TTree* muonPairOutTree[ParamsSet::nSigns];
    static const int nAncestorGroups = 4;
    TTree* QQPairOutTree[ParamsSet::nSigns][2][nAncestorGroups];
    TTree* muonPairOutTreeBinned[ParamsSet::ndRselcs][ParamsSet::nSigns];

    // TH1D* h_numParents;
    TH1D* h_numMuonPairs;
    TH2D* h_ptlead_pair_pt[ParamsSet::nSigns][2];
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

    // std::vector<std::string> parentGroupLabels = {"direct b","b to c","direct c","s/light","single photon", "Drell-Yan"};
    std::vector<std::string> parentGroupLabels = {"direct b","b to c","direct c","s/light","single photon"};
    int nParentGroups = parentGroupLabels.size();    
    
    std::vector<std::string> ancestor_labels = {"gg", "gq", "single g", "q qbar"};
    
    std::vector<std::string> samePrtsLabels = {"Same Parents", "Different Parents"};
    std::vector<std::string> bb_op_one_b_one_btoc_labels = {"Same b", "Involve osc(s)", "From different ancestors", "Others"};

    // std::vector<std::string> osc_labels = {"0 osc, one b one b-to-c", "0 osc, others", "1 osc, one b one b-to-c", "1 osc, regular", "2 oscs, one b one b-to-c", "2 oscs, others"};
    std::vector<std::string> osc_labels = {"0 osc, one c-tag", "0 osc, others(*)", "1 osc, one c-tag(*)", "1 osc, regular", "2 oscs, one c-tag", "2 oscs, others(*)"};
    std::vector<std::string> num_hard_scatt_out_labels = {"2","3","more"};

// --------------------- class methods ---------------------------
  
    void InitInput() override;
    void InitTempVariables() override;
    void ProcessData() override;
    bool PassCuts(const std::shared_ptr<MuonPair>& mpair) override;

    void InitOutput() override;
    void HistAdjust() override;

    void FillMuonPairTree();
    void FillMuonPair(int pair_ind, std::shared_ptr<MuonPairPowheg>& mpair);
    int  ParentGrouping(std::vector<int>& parent_ids, bool c_tag, bool prev_is_lepton);
    void GetPtEtaPhiMFromBarcode(int barcode, std::vector<float>* pt_eta_phi_m);
    int  UpdateCurParents(bool isMuon1, std::vector<int>& cur_prt_bars, std::vector<int>& cur_prt_ids, int hf_quark_index = -1);
    int  FindHeavyQuarks(std::vector<int>& cur_prt_ids, int quark_type, bool isMuon1, int hadron_child_id = 0);
    void SingleMuonAncestorTracing(bool isMuon1);
    int  AncestorGrouping(std::vector<int>& ancestor_ids);
    void HardScatteringAnalysis(std::vector<int>& ancestor_bars, std::vector<int>& ancestor_ids, int sign_dphi_mode, int ancestor_grp);
    void PrintHistory(std::ostream* f, bool print_single, bool muon1_sameancestor);
    // void SameSignSameAncestorsAnalysis(bool near_side, bool one_b_one_btoc, bool print_history = false);
    void KinematicCorrPlots(int isign, int jdphi);
    void MuonPairTagsReinit();
    void CheckIfFromSameB();
    void MuonPairAncestorTracing();
    void Finalize();
  

public :
    // int mode = 2;
    bool is_full_sample = true;
    int full_sample_batch_num;
    std::string mc_mode = "mc_truth_cc";
    bool print_prt_history = false;
    bool print_specific_prt_history = false;
    PowhegNTupleFirstPass(){
        numCuts = numCuts_MC;
        cutLabels = cutLabels_MC;
        crossx_cut = 5 * pow(10,8);
        std::cout << "Powheg Ntuple processing script:" << std::endl;
        std::cout << "The following public variable(s) MUST(????? UNSURE???) be checked:" << std::endl;
        std::cout << "mc_mode: string that takes value mc_truth_cc or mc_truth_bb" << std::endl;
        std::cout << "is_full_sample" << std::endl;
        std::cout << "full_sample_batch_num" << std::endl;
    }
    ~PowhegNTupleFirstPass(){}
    void Run() override;
};


#endif


