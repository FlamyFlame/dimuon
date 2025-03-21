#ifndef __PythiaNTupleFirstPass_h__
#define __PythiaNTupleFirstPass_h__

#include "../MuonObjectsParamsAndHelpers/MuonPairPythia.h"
#include "../MuonObjectsParamsAndHelpers/TruthQQPair.h"
#include "../MuonObjectsParamsAndHelpers/muon_pair_enums_MC.h"
#include "../MuonObjectsParamsAndHelpers/struct_particle.h"
#include "MuonNTupleFirstPassBaseClass.c"
#include "time.h"

class PythiaNTupleFirstPass : public MuonNTupleFirstPassBaseClass{


private:

// --------------------- general settings ---------------------------

    ParamsSet pms;

    static const int nBeamTypes = 4;
	static const int nKinRanges = 5;

    std::string kin_dirs[nKinRanges] = {"k0/", "k1/", "k2/", "k3/", "k4/"};
    std::string beam_dirs[nBeamTypes] = {"pp/", "pn/", "np/", "nn/"};
	int nfiles_base[nBeamTypes] = {4, 6, 6, 9};
    // int nfiles_base[nBeamTypes] = {1,1,1,1};


	Long64_t nevents[nKinRanges] = {0,0,0,0,0};
    Long64_t nevents_accum[nKinRanges] = {0,0,0,0,0};
    Long64_t njobs_accum[nKinRanges] = {0,0,0,0,0};
	Long64_t njobs[nKinRanges] = {0,0,0,0,0};
	int nevents_per_file[nKinRanges] = {10,100,5000,20000,20000};

    std::string py_dir = "/usatlas/u/yuhanguo/usatlasdata/pythia/";

    std::vector<float> kinRanges = {5., 10., 25., 60., 120., 3200.};
    
    bool new_run;
    std::string batch_suffix = "";
    std::string outfile_name;
    std::string outhistfile_name;
    std::vector<std::string> job_dirs;
    std::vector<std::vector<bool>> kn_in_job;
    std::vector<std::vector<int>> nfiles_factor;

// --------------------- resonance tagging & tracking ---------------------------

    std::vector<int> resonance_ids;
    std::map<int, std::pair<std::string,double>> resonance_id_to_name_and_crossx_map;

// --------------------- crossx tracking ---------------------------
    int     pair_counter = 0;
    double  total_crossx = 0.;
    double  from_resonance_total_crossx = 0.;
    double  from_same_b_total_crossx = 0.;
    double  either_from_tau_total_crossx = 0.;
    double  both_incoming_total_crossx = 0.;
    double  both_incoming_FE_QQ_from_same_b_total_crossx = 0.;
    double  FE_total_crossx = 0.;
    double  from_same_gluon_spitting_total_crossx = 0.;
    double  hard_QQ_scatt_total_crossx = 0.;
    double  FE_from_same_b_total_crossx = 0.;
    double  FE_from_same_GS_total_crossx = 0.;
    double  FE_from_diff_ancestors_total_crossx = 0.;
    double  FE_from_same_ancestors_not_same_b_or_gs_total_crossx = 0.;
// --------------------- input files & trees & data for setting branches ---------------------------

    // std::vector <TChain*> evChain;	// each element: pointer to a collection of the ttrees named PyTree (recording event-level info)
    // std::vector <TChain*> metaChain;	// each element: pointer to a collection of the ttrees named meta_tree (recording job-level info)
    TChain* evChain;
    TChain* metaChain;
    
    std::shared_ptr<MuonPairPythia> mpair;
    MuonPairPythia* mpair_raw_ptr = nullptr;

    // for meta tree
    double efficiency = 1.;
 	double ev_weight = 0.;

 	std::vector<double> eff_list{};
 	std::vector<double> sigma_list{};

    // for PyTree

	double QHard;
	double pTHat;
	double mHat;

    std::vector<int>* truth_id = nullptr;
    std::vector<int>* truth_barcode = nullptr;
    std::vector<int>* truth_status = nullptr;
    // std::vector<int>* truth_qual = nullptr;
    std::vector<double>* truth_m = nullptr;
    std::vector<double>* truth_pt = nullptr;
    std::vector<double>* truth_eta = nullptr;
    std::vector<double>* truth_phi = nullptr;
    std::vector<int>* truth_mother1 = nullptr;
    std::vector<int>* truth_mother2 = nullptr;
    // std::vector<int>* truth_daughter1 = nullptr;
    // std::vector<int>* truth_daughter2 = nullptr;

    std::vector<double>   *muon_pair_muon1_pt           =nullptr;
    std::vector<double>   *muon_pair_muon1_eta       =nullptr;
    std::vector<double>   *muon_pair_muon1_phi          =nullptr;
    std::vector<int>     *muon_pair_muon1_ch         =nullptr;
    std::vector<int>     *muon_pair_muon1_bar         =nullptr;

    std::vector<double>   *muon_pair_muon2_pt       =nullptr;
    std::vector<double>   *muon_pair_muon2_eta          =nullptr;
    std::vector<double>   *muon_pair_muon2_phi         =nullptr;
    std::vector<int>     *muon_pair_muon2_ch         =nullptr;
    std::vector<int>     *muon_pair_muon2_bar         =nullptr;

// --------------------- temporary variables (muon, muonpair objects, vectors, etc.) ---------------------------
  
    long nentries_k0;
    long nentries_k1;
    long nentries_k2;
    long nentries_k3;
    long nentries_k4;

    Muon* tempmuon = nullptr;
    // std::vector<std::shared_ptr<MuonPairPythia>> muon_pair_list_cur_event;
    std::vector<std::shared_ptr<MuonPairPythia>> muon_pair_list_cur_event_pre_resonance_cut;

    // TruthQQPair* qqpair = nullptr;

    int m1_youngest_non_chadron_parent_barcode;
    int m2_youngest_non_chadron_parent_barcode;

    int m1_eldest_bhadron_barcode;
    int m2_eldest_bhadron_barcode;

    // bool same_ancestors;
    bool m1_ancestor_is_incoming = false;
    bool m2_ancestor_is_incoming = false;

    bool m1_from_tau;
    bool m2_from_tau;

    int m1_resonance_barcode; // check if from same resonance
    int m2_resonance_barcode; // check if from same resonance

    int nmuonpairs;

    TLorentzVector vg, vQ1, vQ2, vm1out1, vm1out2, vm2out1, vm2out2;
    // bool gs_4vec_correctly_set;

    std::vector<std::vector<int>>* m1_history;
    std::vector<std::vector<int>>* m2_history;
    std::vector<std::vector<Particle>>* m1_history_particle;
    std::vector<std::vector<Particle>>* m2_history_particle;
    // std::vector<std::vector<int>>* m1_single_gluon_history;
    // std::vector<std::vector<int>>* m2_single_gluon_history;
    // std::vector<std::vector<Particle>>* m1_single_gluon_history_particle;
    // std::vector<std::vector<Particle>>* m2_single_gluon_history_particle;

    std::vector<int> m1_parton_ancestor_ids;
    std::vector<int> m2_parton_ancestor_ids;
    std::vector<int> m1_parton_ancestor_bars;
    std::vector<int> m2_parton_ancestor_bars;

    std::vector<int> m1_multi_hf_quark_ids;
    std::vector<int> m2_multi_hf_quark_ids;

    bool m1_c_tag;
    bool m2_c_tag;
    bool m1_osc;
    bool m2_osc;

    bool from_same_gluon_splitting;

    double skipped_event_crossx = 0;
    bool skip_event; // this should only be used rarely
                     // for now, this is used for the "old" runs where the truth_daughters are not recorded
                     // so that if two heavy quarks appear in the same string
                     // we "throw away" the event (still store it as an entry in the muon-pair TTree)
                     // but do not fill it in the HF-muon-pair-ancestor-category histograms (since the category info may not be accurate)
                     // we'll evaluate the contributions of such events to the total cross section
                     // to double check if skipping them has a negligible effect

    bool m1_from_hard_scatt_before_gs;
    bool m1_from_hard_scatt_after_gs;
    int  m1_hard_scatt_in_bar1;
    bool m2_from_hard_scatt_before_gs;
    bool m2_from_hard_scatt_after_gs;
    int  m2_hard_scatt_in_bar1;

    double m1_hard_scatt_Q;
    double m2_hard_scatt_Q;

    // int m1_hard_scatt_category;
    // int m2_hard_scatt_category;
    // int muon_pair_origin_category;

    int m1_earliest_parent_id;
    int m2_earliest_parent_id;

    // std::vector<float>* m1_last_hf_hadron_prt_pt_eta_phi_m;
    // std::vector<float>* m2_last_hf_hadron_prt_pt_eta_phi_m;
    std::vector<float>* m1_first_hf_hadron_prt_pt_eta_phi_m;
    std::vector<float>* m2_first_hf_hadron_prt_pt_eta_phi_m;
    // std::vector<float>* m1_hq_ancestor_pt_eta_phi_m;
    // std::vector<float>* m2_hq_ancestor_pt_eta_phi_m;

// --------------------- output file, histograms & trees ---------------------------
  
    std::ofstream* m_crossx_summary_file = nullptr;
    std::ofstream* m_unspecified_parent_file = nullptr;
    std::ofstream* m_FE_file = nullptr;
    std::ofstream* m_very_bad_warning_file = nullptr;
    std::ofstream* m_very_low_minv_resonance_file = nullptr;
    std::ofstream* m_hard_scattering_warning_file = nullptr;
    std::ofstream* m_others_category_file = nullptr;
    std::ofstream* m_b_parent_file[ParamsSet::nSigns][2];
    std::ofstream* m_c_parent_file[ParamsSet::nSigns][2];

    TTree* meta_tree;
    TTree* muonPairOutTree[ParamsSet::nSigns];
    TTree* muonPairOutTreeKinRange[nKinRanges][ParamsSet::nSigns];

    static const int nAncestorGroups = 4;
    // TTree* QQPairOutTree[ParamsSet::nSigns][2][nAncestorGroups];

    TH1D* h_numMuonPairs;

    TH2D* h_ptlead_pair_pt[ParamsSet::nSigns][2];
    TH2D* h_parent_groups[ParamsSet::nSigns][2];

    TH1D* h_Qsplit_gs_ISR_one_hard_scatt[ParamsSet::nSigns][2];
    TH1D* h_Qsplit_gs_FSR[ParamsSet::nSigns][2];
    // TH1D* h_mHard_FC;

    TH1D* h_Qsplit_to_mHard_gs_ISR_one_hard_scatt[ParamsSet::nSigns][2];
    TH1D* h_Qsplit_to_mHard_gs_FSR[ParamsSet::nSigns][2];
    // TH1D* h_mQQ_to_mHard_gs_ISR_one_hard_scatt;
    // TH1D* h_mQQ_to_Qsplit_gs_ISR_one_hard_scatt;
    
    TH2D* h_hard_scatt_categr_vs_origin_categr[ParamsSet::nSigns][2];
    TH2D* h_both_from_b_hard_scatt_categr_vs_origin_categr[ParamsSet::nSigns][2];
    TH2D* h_both_from_c_hard_scatt_categr_vs_origin_categr[ParamsSet::nSigns][2];
    TH2D* h_one_from_b_one_from_c_hard_scatt_categr_vs_origin_categr[ParamsSet::nSigns][2];
    
    TH1D* h_muon_pair_origin_categr[ParamsSet::nSigns][2];
    TH1D* h_both_from_b_muon_pair_origin_categr[ParamsSet::nSigns][2];
    TH1D* h_both_from_c_muon_pair_origin_categr[ParamsSet::nSigns][2];
    TH1D* h_one_from_b_one_from_c_muon_pair_origin_categr[ParamsSet::nSigns][2];

    TH1D* h_QQ_DR[ParamsSet::nSigns][2][4];
    TH1D* h_QQ_Dphi[ParamsSet::nSigns][2][4];
    TH1D* h_QQ_minv[ParamsSet::nSigns][2][4];
    TH1D* h_QQ_pair_pt_ptlead_ratio[ParamsSet::nSigns][2][4];
    TH1D* h_QQ_pt_avg[ParamsSet::nSigns][2][4];
    TH1D* h_QQ_asym[ParamsSet::nSigns][2][4];
    TH1D* h_QQ_minv_s_cm_ratio[ParamsSet::nSigns][2][4];

    TH2D* h_QQ_ptlead_pair_pt[ParamsSet::nSigns][2][4];
    TH2D* h_QQ_pt1_pt2[ParamsSet::nSigns][2][4];
    TH2D* h_QQ_Deta_Dphi[ParamsSet::nSigns][2][4];
    TH2D* h_QQ_eta1_eta2[ParamsSet::nSigns][2][4];
    TH2D* h_QQ_minv_pair_pt[ParamsSet::nSigns][2][4];
    TH2D* h_QQ_minv_Dphi[ParamsSet::nSigns][2][4];

    static const int npTbins = 3;
    TH1D* h_crossx;
    TH1D* h_crossx_pt_binned[npTbins];
    TH1D* h_op_one_b_one_btoc[2];
    TH1D* h_both_from_b_same_ancestors[ParamsSet::nSigns][2];
    TH1D* h_both_from_c_same_ancestors[ParamsSet::nSigns][2];
    TH1D* h_both_from_b_ancestor_sp[ParamsSet::nSigns][2];
    TH2D* h_both_from_b_ancestor_dp[ParamsSet::nSigns][2];
    TH1D* h_both_from_c_ancestor_sp[ParamsSet::nSigns][2];
    TH2D* h_both_from_c_ancestor_dp[ParamsSet::nSigns][2];

    // TH1D* h_dphi_bb_op_near_both_from_b;
    // TH1D* h_dphi_bb_op_near_one_b_one_btoc;
    // TH1D* h_dphi_bb_ss_near;
    // TH1D* h_dphi_bb_op_near;
    // TH1D* h_dphi_bb_op_near_from_same_b;
    // TH1D* h_bb_ss_near_involv_osc;
    // TH1D* h_bb_ss_away_involv_osc;

    // TH2D* h_cc_ss_small_dphi_prt_gps;
    // TH1D* h_cc_ss_small_dphi_same_ancestors;
    // TH1D* h_cc_ss_small_dphi_sp;
    // TH2D* h_cc_ss_small_dphi_dp;
    // TH2D* h_cc_ss_plateau_prt_gps;
    // TH1D* h_cc_ss_plateau_same_ancestors;
    // TH1D* h_cc_ss_plateau_sp;
    // TH2D* h_cc_ss_plateau_dp;

    // TH1D* h_num_hard_scatt_out[ParamsSet::nSigns][2];
    // TH1D* h_pt_muon_pt_closest_hadr_ratio[ParamsSet::nSigns][2];
    // TH1D* h_pt_closest_hadr_pt_furthest_hadr_ratio[ParamsSet::nSigns][2];
    // TH1D* h_pt_hadr_hq_ratio[ParamsSet::nSigns][2];
    // TH1D* h_dphi_muon_closest_hadr[ParamsSet::nSigns][2];

    std::vector<std::string> parentGroupLabels = {"direct b","b to c","direct c","s/light","direct photon", "Drell-Yan"};
    // std::vector<std::string> parentGroupLabels = {"direct b","b to c","direct c","s/light","direct photon"};
    int nParentGroups = parentGroupLabels.size();
    
    // std::vector<std::string> HardScattCategoryLabels = {"gg->gg",   "qqbar->gg",  "gq->gq",  "flavor creation", "flavor excitation", "QQ'->QQ'",  "QQbar->gg",  "gg->qqbar", "Drell-Yan", "qq'->qq'"};
    std::vector<std::string> HardScattCategoryLabels = {"gg->gg",   "qqbar->gg",  "gq->gq",  "flavor creation", "flavor excitation", "QQ'->QQ'",  "gg->qqbar", "qq'->qq'"};
    int nHardScattCategories = HardScattCategoryLabels.size();
    
    // std::vector<std::string> HFMuonPairOriginCategoryLabels = {"flavor creation", "same G(/#gamma)S FSR", "same GS ISR 0 hard scatt", "same GS ISR 1 hard scatt", "same GS ISR both hard scatt", "diff GS same hard scatt", "Drell-Yan", "photon to QQbar", "others"};
    std::vector<std::string> HFMuonPairOriginCategoryLabels = {"flavor creation", "same G(/#gamma)S FSR", "same GS ISR 0 hard scatt", "same GS ISR 1 hard scatt", "same GS ISR both hard scatt", "diff GS same hard scatt", "others"};
    int nHFMuonPairOriginCategories = HFMuonPairOriginCategoryLabels.size();
    // cout << "Number of muon pair origin categories (should be 8): " << nHFMuonPairOriginCategories << std::endl;
    // std::vector<std::string> HFMuonPairOriginCategoryLabels = {"flavor creation", "same GS FSR", "same GS ISR 0 HS", "same GS ISR 1 HS", "same GS ISR both HS", "diff GS same HS", "others"};

    // std::vector<std::string> ancestor_labels = {"gg", "gq", "single g", "q qbar"};
    
    std::vector<std::string> samePrtsLabels = {"Same Ancestors", "Different Ancestors"};
    std::vector<std::string> bb_op_one_b_one_btoc_labels = {"Same b", "Involve osc(s)", "From different ancestors", "Others"};

    // std::vector<std::string> osc_labels = {"0 osc, one b one b-to-c", "0 osc, others", "1 osc, one b one b-to-c", "1 osc, regular", "2 oscs, one b one b-to-c", "2 oscs, others"};
    std::vector<std::string> osc_labels = {"0 osc, one c-tag", "0 osc, others(*)", "1 osc, one c-tag(*)", "1 osc, regular", "2 oscs, one c-tag", "2 oscs, others(*)"};
    std::vector<std::string> num_hard_scatt_out_labels = {"2","3","more"};


// --------------------- class methods ---------------------------
  
    void InitInput() override;
    void InitOutput() override;
    void InitTempVariables() override;
    void ProcessData() override;
    bool PassCuts(const std::shared_ptr<MuonPair>& mpair) override;

    void FillMuonPairTreePythia(int nkin);
    void HistAdjust() override;

    void FillMuonPair(int pair_ind, std::shared_ptr<MuonPairPythia>& mpair);
    void SetInputOutputFilesFromBatch();
    void InputSanityCheck();
    void ResonanceNameMap();
    int ParentGrouping(std::vector<int>& parent_ids, bool c_tag, bool prev_is_lepton);
    void GetPtEtaPhiMFromBarcode(int barcode, std::vector<float>* pt_eta_phi_m);
    std::pair<int,int>  UpdateCurParents(bool isMuon1, std::vector<int>& cur_prt_bars, std::vector<int>& cur_prt_ids, bool before_gs, bool& prev_out_hard_scatt, int hf_quark_index = -1000000);    
    int  FindHeavyQuarks(std::vector<int>& cur_prt_ids, std::vector<int>& cur_prt_bars, int quark_type, bool isMuon1, int prev_hq_bar, int hadron_child_id = 0);
    void PrintHistory(std::ostream* f, bool print_single, bool muon1_sameancestor, bool print_if_same_ancestor = false);
    void SingleMuonAncestorTracing(bool isMuon1);
    // int  AncestorGrouping(std::vector<int>& ancestor_ids, bool sameprts, std::string mc_mode);
    int  HardAnalysisCategr(int in_bar1, int in_bar2);
    // void HardScatteringAnalysis(std::vector<int>& ancestor_bars, std::vector<int>& ancestor_ids, int sign_dphi_mode, int ancestor_grp);
    void PrintHistory(std::ostream* f, int mode_single_both, int mode_info);
    // void SameSignSameAncestorsAnalysis(bool near_side, bool one_b_one_btoc, bool print_history = false);
    // void KinematicCorrPlots(int isign, int jdphi);
    int  GluonHistoryTracking(int gluon_bar, bool isMuon1);
    void MuonPairTagsReinit();
    void FillCategoryHistograms(TH1D* horig, TH2D* hhard_vs_orig);
    void HFMuonPairAnalysis();
    void MuonPairAncestorTracing();
    void CrossxClear();
    void PerPairCrossxUpdate();
    void WriteCrossxSummary();
    void Finalize();

public :

    bool print_prt_history = false;
    bool print_others_history = true;
    bool print_unspecified_parent = false;
    bool print_bad_warnings = false;
    bool print_very_low_minv = false;
    bool print_FE = false;

    bool turn_data_resonance_cuts_on = false; // turn on minv-based data cuts for resonances: by default false
    int batch_num = 0; // between 1-6: process the privately-generated pythia sample in batch

    PythiaNTupleFirstPass();
    ~PythiaNTupleFirstPass(){}
    void Run() override;
    // float filter_effcy;

};


#endif