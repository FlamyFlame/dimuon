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
    static const int nDphis = 2;

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


// --------------------- input files & trees & data for setting branches---------------------------

   	TFile* inFile;
    std::string file_name = "/usatlas/u/yuhanguo/usatlasdata/pythia/muon_pairs_pythia_combined";
    
    TTree *inTree[ParamsSet::nSigns];

    double weight[ParamsSet::nSigns];
    double pTHat[ParamsSet::nSigns];
    double Qsplit[ParamsSet::nSigns];
    double mHard_relevant[ParamsSet::nSigns];

    int muon_pair_flavor_category[ParamsSet::nSigns];
    int muon_pair_origin_category[ParamsSet::nSigns];

    bool from_same_resonance[ParamsSet::nSigns];
    bool resonance_contaminated[ParamsSet::nSigns];
    bool from_same_b[ParamsSet::nSigns];
    bool both_from_b[ParamsSet::nSigns];
    bool both_from_c[ParamsSet::nSigns];
    bool data_resonance_or_reso_contam_tagged_old[ParamsSet::nSigns];
    bool data_resonance_or_reso_contam_tagged_new[ParamsSet::nSigns];

    float pair_dPoverP[ParamsSet::nSigns];
    float pt_lead[ParamsSet::nSigns];
    float pair_pt[ParamsSet::nSigns];
    float pair_eta[ParamsSet::nSigns];
    float pair_y[ParamsSet::nSigns];
    float asym[ParamsSet::nSigns];
    float dpt[ParamsSet::nSigns];
    float deta[ParamsSet::nSigns];
    float etaavg[ParamsSet::nSigns];
    float phiavg[ParamsSet::nSigns];
    float dphi[ParamsSet::nSigns];
    float dr[ParamsSet::nSigns];
    float minv[ParamsSet::nSigns];
    float m1pt[ParamsSet::nSigns];
    float m2pt[ParamsSet::nSigns];
    float m1eta[ParamsSet::nSigns];
    float m2eta[ParamsSet::nSigns];
    float m1phi[ParamsSet::nSigns];
    float m2phi[ParamsSet::nSigns];
    // int m1charge[ParamsSet::nSigns];
    // int m2charge[ParamsSet::nSigns];

// --------------------- class methods ---------------------------

   	void InitInput();
    void InitOutput();
   	void InitHists();
   	void ProcessData();
    void WriteOutput();
    bool PassSingleMuonGapCut(float meta, float mpt, int mcharge);
   	void FillHistograms(int nsign);
    // void FillCrossxHistogramsTruthSingleB(int nsign);    
   	void FillPtBinnedHistograms(int npt, int nsign);

public:
    int mode = 1;
    bool with_data_resonance_cuts = false;
    // std::string file_suffix = "_allto0318";
    // std::string file_suffix = "_after0322";
  	MuonPairPlottingPythia();
  	~MuonPairPlottingPythia(){}
  	void Run();

};

#endif