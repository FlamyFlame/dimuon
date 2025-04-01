// Assumes that we exclude resonances

#ifndef __MuonPairPlottingPythia_h__
#define __MuonPairPlottingPythia_h__

#include <TROOT.h>
// #include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <string.h>
#include  <stdlib.h>
#include <fstream>
#include "../MuonObjectsParamsAndHelpers/ParamsSet.h"
#include "../MuonObjectsParamsAndHelpers/muon_pair_enums_MC_utils.h"

#include "vector"
#include "TH1D.h"
#include "TH2D.h"

class MuonPairPlottingPythia{
  	//updates histograms - centrality binned
    //if need histograms with all centrality bins added together
    //can add at the histograms level


private:
// --------------------- general settings ---------------------------

    // int mode = 1;
    // bool isScram = true;
    // bool isScram = false;

  	ParamsSet pms;
    // static const int nFiles = 3;
    static const int nFiles = 2;
    static const int nKinRanges = 5;
    static const int nDphis = 2;
    // static const int pair_flavor_index::nFlavors = ; // both from b, both from c, one b one c
    std::vector<float> kinRanges = {5., 10., 25., 60., 120., 3200.};

    std::string with_data_resonance_cuts_suffix;

// --------------------- histogram settings needed for multiple class methods ---------------------------
    int npair_eta_bins_coarse = 24;
    int npair_pT_bins_coarse = 30;
    float pair_eta_min = -2.4;
    float pair_eta_max = 2.4;
    float pair_pt_min = 0;
    float pair_pt_max = 30;

// --------------------- output file & histograms ---------------------------

  	TFile *outFile = nullptr;
    std::string output_file_path;

    std::string dphi_regions[nDphis] = {"near", "away"};

    TH2D* h_crossx_truth_from_single_b_vs_pair_pt_pair_eta;

    TH1D* h_minv_sub_GeV_signal_no_res_cut;
    TH1D* h_minv_sub_GeV_signal_old_res_cut;
    TH1D* h_minv_sub_GeV_signal_new_res_cut;
    TH1D* h_minv_sub_GeV_jacobian_corrected_signal_no_res_cut;
    TH1D* h_minv_sub_GeV_jacobian_corrected_signal_old_res_cut;
    TH1D* h_minv_sub_GeV_jacobian_corrected_signal_new_res_cut;
    TH1D* h_minv_sub_GeV_resonance_and_res_contam_bkg_no_res_cut;
    TH1D* h_minv_sub_GeV_resonance_and_res_contam_bkg_old_res_cut;
    TH1D* h_minv_sub_GeV_resonance_and_res_contam_bkg_new_res_cut;
    TH1D* h_minv_sub_GeV_jacobian_corrected_resonance_and_res_contam_bkg_no_res_cut;
    TH1D* h_minv_sub_GeV_jacobian_corrected_resonance_and_res_contam_bkg_old_res_cut;
    TH1D* h_minv_sub_GeV_jacobian_corrected_resonance_and_res_contam_bkg_new_res_cut;
    TH1D* h_minv_single_b_region_signal_no_res_cut;
    TH1D* h_minv_single_b_region_signal_old_res_cut;
    TH1D* h_minv_single_b_region_signal_new_res_cut;
    TH1D* h_minv_single_b_region_jacobian_corrected_signal_no_res_cut;
    TH1D* h_minv_single_b_region_jacobian_corrected_signal_old_res_cut;
    TH1D* h_minv_single_b_region_jacobian_corrected_signal_new_res_cut;
    TH1D* h_minv_single_b_region_resonance_and_res_contam_bkg_no_res_cut;
    TH1D* h_minv_single_b_region_resonance_and_res_contam_bkg_old_res_cut;
    TH1D* h_minv_single_b_region_resonance_and_res_contam_bkg_new_res_cut;
    TH1D* h_minv_single_b_region_jacobian_corrected_resonance_and_res_contam_bkg_no_res_cut;
    TH1D* h_minv_single_b_region_jacobian_corrected_resonance_and_res_contam_bkg_old_res_cut;
    TH1D* h_minv_single_b_region_jacobian_corrected_resonance_and_res_contam_bkg_new_res_cut;

    TH1D* h_pair_dP_overP[nDphis][ParamsSet::nSigns];
    TH1D* h_Dphi[nDphis][ParamsSet::nSigns];
    TH1D* h_DR[nDphis][ParamsSet::nSigns];
    TH1D* h_DR_zoomin[nDphis][ParamsSet::nSigns];
    TH1D* h_Deta_zoomin[nDphis][ParamsSet::nSigns];
    TH1D* h_Dphi_zoomin[nDphis][ParamsSet::nSigns];
    TH1D* h_DR_jacobian_corrected[nDphis][ParamsSet::nSigns];
    TH1D* h_DR_zoomin_jacobian_corrected[nDphis][ParamsSet::nSigns];
    TH1D* h_pair_y[nDphis][ParamsSet::nSigns];
    TH1D* h_pt_asym[nDphis][ParamsSet::nSigns];
    TH1D* h_psrapidity_ordered_pt_asym[nDphis][ParamsSet::nSigns];
    TH1D* h_pair_pt_ptlead_ratio[nDphis][ParamsSet::nSigns];
    TH1D* h_pair_pt_jacobian_corrected[nDphis][ParamsSet::nSigns];
    TH2D* h_eta_avg_Dphi[nDphis][ParamsSet::nSigns];
    TH2D* h_Deta_Dphi[nDphis][ParamsSet::nSigns];
    TH2D* h_eta1_eta2[nDphis][ParamsSet::nSigns];
    TH2D* h_eta_avg_Deta[nDphis][ParamsSet::nSigns];
    TH2D* h_pt1_pt2[nDphis][ParamsSet::nSigns];
    TH2D* h_ptlead_pair_pt[nDphis][ParamsSet::nSigns];
    TH2D* h_ptlead_pair_pt_zoomin[nDphis][ParamsSet::nSigns];
    TH2D* h_ptlead_pair_pt_log[nDphis][ParamsSet::nSigns];
    TH2D* h_minv_pair_pt[nDphis][ParamsSet::nSigns];
    TH2D* h_minv_pair_pt_zoomin[nDphis][ParamsSet::nSigns];
    TH2D* h_minv_pair_pt_jacobian_corrected[nDphis][ParamsSet::nSigns];
    TH2D* h_minv_pair_pt_zoomin_jacobian_corrected[nDphis][ParamsSet::nSigns];
    TH2D* h_minv_pair_pt_log[nDphis][ParamsSet::nSigns];

    TH2D* h_eta1_eta2_dphicut[ParamsSet::ndphiselcs][ParamsSet::nSigns];

    std::map<int,std::string> origin_grp_map;
    std::map<int,std::string> flavor_grp_map;

    TH1D* h_DR_origin_binned[ParamsSet::nSigns][muon_pair_both_from_open_HF_origin_catgr::nOrigins];
    TH1D* h_DR_zoomin_origin_binned[ParamsSet::nSigns][muon_pair_both_from_open_HF_origin_catgr::nOrigins];
    TH1D* h_Deta_zoomin_origin_binned[ParamsSet::nSigns][muon_pair_both_from_open_HF_origin_catgr::nOrigins];
    TH1D* h_Dphi_zoomin_origin_binned[ParamsSet::nSigns][muon_pair_both_from_open_HF_origin_catgr::nOrigins];
    TH1D* h_DR_jacobian_corrected_origin_binned[ParamsSet::nSigns][muon_pair_both_from_open_HF_origin_catgr::nOrigins];
    TH1D* h_DR_zoomin_jacobian_corrected_origin_binned[ParamsSet::nSigns][muon_pair_both_from_open_HF_origin_catgr::nOrigins];
    TH1D* h_psrapidity_ordered_pt_asym_origin_binned[ParamsSet::nSigns][muon_pair_both_from_open_HF_origin_catgr::nOrigins];
    TH1D* h_pt_asym_origin_binned[ParamsSet::nSigns][muon_pair_both_from_open_HF_origin_catgr::nOrigins];
    TH1D* h_pair_pt_ptlead_ratio_origin_binned[ParamsSet::nSigns][muon_pair_both_from_open_HF_origin_catgr::nOrigins];
    TH1D* h_pair_pt_jacobian_corrected_origin_binned[ParamsSet::nSigns][muon_pair_both_from_open_HF_origin_catgr::nOrigins];

    TH1D* h_Qsplit_origin_binned[ParamsSet::nSigns][muon_pair_both_from_open_HF_origin_catgr::nOrigins];
    TH1D* h_Qsplit_pTHat_ratio_origin_binned[ParamsSet::nSigns][muon_pair_both_from_open_HF_origin_catgr::nOrigins];
    TH1D* h_Qsplit_mHat_ratio_origin_binned[ParamsSet::nSigns][muon_pair_both_from_open_HF_origin_catgr::nOrigins];
    
    TH2D* h_Deta_Dphi_origin_binned[ParamsSet::nSigns][muon_pair_both_from_open_HF_origin_catgr::nOrigins];
    TH2D* h_minv_pair_pt_origin_binned[ParamsSet::nSigns][muon_pair_both_from_open_HF_origin_catgr::nOrigins];
    TH2D* h_minv_pair_pt_zoomin_origin_binned[ParamsSet::nSigns][muon_pair_both_from_open_HF_origin_catgr::nOrigins];
    TH2D* h_minv_pair_pt_jacobian_corrected_origin_binned[ParamsSet::nSigns][muon_pair_both_from_open_HF_origin_catgr::nOrigins];
    TH2D* h_minv_pair_pt_zoomin_jacobian_corrected_origin_binned[ParamsSet::nSigns][muon_pair_both_from_open_HF_origin_catgr::nOrigins];
    TH2D* h_minv_pair_pt_log_origin_binned[ParamsSet::nSigns][muon_pair_both_from_open_HF_origin_catgr::nOrigins];
    TH2D* h_ptlead_pair_pt_origin_binned[ParamsSet::nSigns][muon_pair_both_from_open_HF_origin_catgr::nOrigins];
    TH2D* h_ptlead_pair_pt_zoomin_origin_binned[ParamsSet::nSigns][muon_pair_both_from_open_HF_origin_catgr::nOrigins];
    TH2D* h_ptlead_pair_pt_log_origin_binned[ParamsSet::nSigns][muon_pair_both_from_open_HF_origin_catgr::nOrigins];

    TH1D* h_DR_flavor_binned[ParamsSet::nSigns][pair_flavor_index::nFlavors];
    TH1D* h_DR_zoomin_flavor_binned[ParamsSet::nSigns][pair_flavor_index::nFlavors];
    TH1D* h_Deta_zoomin_flavor_binned[ParamsSet::nSigns][pair_flavor_index::nFlavors];
    TH1D* h_Dphi_zoomin_flavor_binned[ParamsSet::nSigns][pair_flavor_index::nFlavors];
    TH1D* h_DR_jacobian_corrected_flavor_binned[ParamsSet::nSigns][pair_flavor_index::nFlavors];
    TH1D* h_DR_zoomin_jacobian_corrected_flavor_binned[ParamsSet::nSigns][pair_flavor_index::nFlavors];
    TH1D* h_psrapidity_ordered_pt_asym_flavor_binned[ParamsSet::nSigns][pair_flavor_index::nFlavors];
    TH1D* h_pt_asym_flavor_binned[ParamsSet::nSigns][pair_flavor_index::nFlavors];
    TH1D* h_pair_pt_ptlead_ratio_flavor_binned[ParamsSet::nSigns][pair_flavor_index::nFlavors];
    TH1D* h_pair_pt_jacobian_corrected_flavor_binned[ParamsSet::nSigns][pair_flavor_index::nFlavors];
    
    TH2D* h_Deta_Dphi_flavor_binned[ParamsSet::nSigns][pair_flavor_index::nFlavors];
    TH2D* h_minv_pair_pt_flavor_binned[ParamsSet::nSigns][pair_flavor_index::nFlavors];
    TH2D* h_minv_pair_pt_zoomin_flavor_binned[ParamsSet::nSigns][pair_flavor_index::nFlavors];
    TH2D* h_minv_pair_pt_jacobian_corrected_flavor_binned[ParamsSet::nSigns][pair_flavor_index::nFlavors];
    TH2D* h_minv_pair_pt_zoomin_jacobian_corrected_flavor_binned[ParamsSet::nSigns][pair_flavor_index::nFlavors];
    TH2D* h_minv_pair_pt_log_flavor_binned[ParamsSet::nSigns][pair_flavor_index::nFlavors];
    TH2D* h_ptlead_pair_pt_flavor_binned[ParamsSet::nSigns][pair_flavor_index::nFlavors];
    TH2D* h_ptlead_pair_pt_zoomin_flavor_binned[ParamsSet::nSigns][pair_flavor_index::nFlavors];
    TH2D* h_ptlead_pair_pt_log_flavor_binned[ParamsSet::nSigns][pair_flavor_index::nFlavors];


    TH1D* h_kinbin_pair_dP_overP[nKinRanges][nDphis][ParamsSet::nSigns];
    TH1D* h_kinbin_Dphi[nKinRanges][nDphis][ParamsSet::nSigns];
    TH1D* h_kinbin_DR[nKinRanges][nDphis][ParamsSet::nSigns];
    TH1D* h_kinbin_pair_y[nKinRanges][nDphis][ParamsSet::nSigns];
    TH1D* h_kinbin_pt_asym[nKinRanges][nDphis][ParamsSet::nSigns];
    TH1D* h_kinbin_pair_pt_ptlead_ratio[nKinRanges][nDphis][ParamsSet::nSigns];
    TH2D* h_kinbin_eta_avg_Dphi[nKinRanges][nDphis][ParamsSet::nSigns];
    TH2D* h_kinbin_Deta_Dphi[nKinRanges][nDphis][ParamsSet::nSigns];
    TH2D* h_kinbin_eta1_eta2[nKinRanges][nDphis][ParamsSet::nSigns];
    TH2D* h_kinbin_eta_avg_Deta[nKinRanges][nDphis][ParamsSet::nSigns];
    TH2D* h_kinbin_pt1_pt2[nKinRanges][nDphis][ParamsSet::nSigns];
    TH2D* h_kinbin_ptlead_pair_pt[nKinRanges][nDphis][ParamsSet::nSigns];
    TH2D* h_kinbin_minv_pair_pt[nKinRanges][nDphis][ParamsSet::nSigns];
    TH2D* h_kinbin_minv_pair_pt_log[nKinRanges][nDphis][ParamsSet::nSigns];


// --------------------- input files & trees & data for setting branches---------------------------

   	TFile* inFile[nFiles];
    // NEED TO CHECK IF NEW AFTER0322 CONTAINS 0429 + IF THE nevents_before_cuts ARE CORRECT
    std::vector<std::string> file_names = {"/usatlas/u/yuhanguo/usatlasdata/pythia/muon_pairs_pythia_allto0318", "/usatlas/u/yuhanguo/usatlasdata/pythia/muon_pairs_pythia_after0322"};
    
    // we need to multiply the previous weight for each of the five different kinematic ranges (k_i)
    // in each of the two files (file A, B) by the corresponding reweight factor
    // N^i_A / (N^i_A + N^i_B) or N^i_B / (N^i_A + N^i_B), respectively
    std::vector<std::vector<int>> nevents_before_cuts = {{20000, 50000, 1000000, 1000000, 1000000}, {99910, 400000, 1500000, 1500000, 1500000}};


    std::vector<int> nevents_before_cuts_combined {};

    // TTree *inTree_ctrbin[nKinRanges][ParamsSet::nSigns];
    // TTree *inTree[ParamsSet::nSigns];
    TTree *inTree[nFiles][nKinRanges][ParamsSet::nSigns];

    double weight[nFiles][nKinRanges][ParamsSet::nSigns];
    double pTHat[nFiles][nKinRanges][ParamsSet::nSigns];
    double Qsplit[nFiles][nKinRanges][ParamsSet::nSigns];
    double mHard_relevant[nFiles][nKinRanges][ParamsSet::nSigns];

    int muon_pair_flavor_category[nFiles][nKinRanges][ParamsSet::nSigns];
    int muon_pair_origin_category[nFiles][nKinRanges][ParamsSet::nSigns];

    bool from_same_resonance[nFiles][nKinRanges][ParamsSet::nSigns];
    bool resonance_contaminated[nFiles][nKinRanges][ParamsSet::nSigns];
    bool from_same_b[nFiles][nKinRanges][ParamsSet::nSigns];
    bool both_from_b[nFiles][nKinRanges][ParamsSet::nSigns];
    bool both_from_c[nFiles][nKinRanges][ParamsSet::nSigns];
    bool data_resonance_or_reso_contam_tagged_old[nFiles][nKinRanges][ParamsSet::nSigns];
    bool data_resonance_or_reso_contam_tagged_new[nFiles][nKinRanges][ParamsSet::nSigns];

    float pair_dPoverP[nFiles][nKinRanges][ParamsSet::nSigns];
    float pt_lead[nFiles][nKinRanges][ParamsSet::nSigns];
    float pair_pt[nFiles][nKinRanges][ParamsSet::nSigns];
    float pair_eta[nFiles][nKinRanges][ParamsSet::nSigns];
    float pair_y[nFiles][nKinRanges][ParamsSet::nSigns];
    float asym[nFiles][nKinRanges][ParamsSet::nSigns];
    float dpt[nFiles][nKinRanges][ParamsSet::nSigns];
    float deta[nFiles][nKinRanges][ParamsSet::nSigns];
    float etaavg[nFiles][nKinRanges][ParamsSet::nSigns];
    float phiavg[nFiles][nKinRanges][ParamsSet::nSigns];
    float dphi[nFiles][nKinRanges][ParamsSet::nSigns];
    float dr[nFiles][nKinRanges][ParamsSet::nSigns];
    float minv[nFiles][nKinRanges][ParamsSet::nSigns];
    float m1pt[nFiles][nKinRanges][ParamsSet::nSigns];
    float m2pt[nFiles][nKinRanges][ParamsSet::nSigns];
    float m1eta[nFiles][nKinRanges][ParamsSet::nSigns];
    float m2eta[nFiles][nKinRanges][ParamsSet::nSigns];
    float m1phi[nFiles][nKinRanges][ParamsSet::nSigns];
    float m2phi[nFiles][nKinRanges][ParamsSet::nSigns];
    // int m1charge[nFiles][nKinRanges][ParamsSet::nSigns];
    // int m2charge[nFiles][nKinRanges][ParamsSet::nSigns];

// --------------------- class methods ---------------------------

   	void InitInput();
    void InitOutput();
   	void InitHists();
   	void ProcessData();
    void WriteOutput();
    bool PassSingleMuonGapCut(float meta, float mpt, int mcharge);
   	void FillHistograms(int nfile, int nkin, int nsign);
    // void FillCrossxHistogramsTruthSingleB(int nfile, int nkin, int nsign);    
   	void FillPtBinnedHistograms(int nkin, int npt, int nsign);

public:
    int mode = 1;
    bool with_data_resonance_cuts = false;
    bool plot_kin_binned_histograms = false; // if true, plot kinematic-range-binned histograms (default false)
    // std::string file_suffix = "_allto0318";
    // std::string file_suffix = "_after0322";
  	MuonPairPlottingPythia();
  	~MuonPairPlottingPythia(){}
  	void Run();

};

#endif