// Assumes that we exclude resonances

#ifndef __MuonPairPlottingPP_h__
#define __MuonPairPlottingPP_h__

#include <TROOT.h>
// #include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <string.h>
#include  <stdlib.h>
#include <fstream>
#include "../MuonObjectsParamsAndHelpers/ParamsSet.h"
#include "../MuonObjectsParamsAndHelpers/MuonPair.h"
#include "vector"
#include "TH1D.h"
#include "TH2D.h"

class MuonPairPlottingPP{
  	//updates histograms - centrality binned
    //if need histograms with all centrality bins added together
    //can add at the histograms level


private:
  	// --------------------- general settings ---------------------------

    // bool isScram = true;
    // bool isScram = false;

  	ParamsSet pms;

// --------------------- histogram settings needed for multiple class methods ---------------------------
    int npair_eta_bins_coarse = 24;
    int npair_pT_bins_coarse = 30;
    float pair_eta_min = -2.4;
    float pair_eta_max = 2.4;
    float pair_pt_min = 0;
    float pair_pt_max = 30;

  	// --------------------- output file & histograms ---------------------------

  	TFile *outFile = nullptr;
    std::string trig_suffix;
    
    std::string dphi_regions[2] = {"near", "away"};

    TH2D* h_crossx_dR_cut_vs_pair_pt_pair_eta;

    TH1D* h_pair_dP_overP[2][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH1D* h_Dphi[2][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH1D* h_DR[2][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH1D* h_DR_zoomin[2][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH1D* h_DR_jacobian_corrected[2][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH1D* h_DR_zoomin_jacobian_corrected[2][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH1D* h_pair_y[2][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH1D* h_pt_asym[2][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH1D* h_pair_pt_ptlead_ratio[2][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_eta_avg_Dphi[2][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_Deta_Dphi[2][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_eta1_eta2[2][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_eta_avg_Deta[2][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_pt1_pt2[2][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_ptlead_pair_pt[2][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_ptlead_pair_pt_zoomin[2][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_ptlead_pair_pt_log[2][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_minv_pair_pt[2][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_minv_pair_pt_zoomin[2][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_minv_pair_pt_jacobian_corrected[2][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_minv_pair_pt_zoomin_jacobian_corrected[2][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_minv_pair_pt_log[2][ParamsSet::nSigns][ParamsSet::nGapCuts];

    TH1D* h_Deta_mu4[ParamsSet::nSigns];
    TH1D* h_Deta_zoomin_mu4[ParamsSet::nSigns];
    TH1D* h_Dphi_mu4[ParamsSet::nSigns];
    TH1D* h_Dphi_zoomin_mu4[ParamsSet::nSigns];
    TH1D* h_DR_mu4[ParamsSet::nSigns];
    TH1D* h_DR_zoomin_mu4[ParamsSet::nSigns];
    TH1D* h_pt2nd_mu4[ParamsSet::nSigns];
    TH2D* h_pt2nd_vs_q_eta_2nd_mu4[ParamsSet::nSigns];
    TH2D* h_pair_eta_vs_pair_pT_mu4[ParamsSet::nSigns];
    TH2D* h_Deta_Dphi_mu4[ParamsSet::nSigns];
    TH2D* h_eta1_eta2_mu4[ParamsSet::nSigns];
    TH2D* h_eta_avg_Deta_mu4[ParamsSet::nSigns];
    TH2D* h_eta_avg_Dphi_mu4[ParamsSet::nSigns];
    TH1D* h_minv_zoomin_mu4[ParamsSet::nSigns];
    TH1D* h_pair_pt_log_mu4[ParamsSet::nSigns];
    TH2D* h_minv_pair_pt_log_mu4[ParamsSet::nSigns];

    TH1D* h_Deta_mu4_mu4noL1[ParamsSet::nSigns];
    TH1D* h_Deta_zoomin_mu4_mu4noL1[ParamsSet::nSigns];
    TH1D* h_Dphi_mu4_mu4noL1[ParamsSet::nSigns];
    TH1D* h_Dphi_zoomin_mu4_mu4noL1[ParamsSet::nSigns];
    TH1D* h_DR_mu4_mu4noL1[ParamsSet::nSigns];
    TH1D* h_DR_zoomin_mu4_mu4noL1[ParamsSet::nSigns];
    TH1D* h_pt2nd_mu4_mu4noL1[ParamsSet::nSigns];
    TH2D* h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1[ParamsSet::nSigns];
    TH2D* h_pair_eta_vs_pair_pT_mu4_mu4noL1[ParamsSet::nSigns];
    TH2D* h_Deta_Dphi_mu4_mu4noL1[ParamsSet::nSigns];
    TH2D* h_eta1_eta2_mu4_mu4noL1[ParamsSet::nSigns];
    TH2D* h_eta_avg_Deta_mu4_mu4noL1[ParamsSet::nSigns];
    TH2D* h_eta_avg_Dphi_mu4_mu4noL1[ParamsSet::nSigns];
    TH1D* h_minv_zoomin_mu4_mu4noL1[ParamsSet::nSigns];
    TH1D* h_pair_pt_log_mu4_mu4noL1[ParamsSet::nSigns];
    TH2D* h_minv_pair_pt_log_mu4_mu4noL1[ParamsSet::nSigns];

    TH1D* h_Deta_mu4_mu4noL1_excl[ParamsSet::nSigns];
    TH1D* h_Deta_zoomin_mu4_mu4noL1_excl[ParamsSet::nSigns];
    TH1D* h_Dphi_mu4_mu4noL1_excl[ParamsSet::nSigns];
    TH1D* h_Dphi_zoomin_mu4_mu4noL1_excl[ParamsSet::nSigns];
    TH1D* h_DR_mu4_mu4noL1_excl[ParamsSet::nSigns];
    TH1D* h_DR_zoomin_mu4_mu4noL1_excl[ParamsSet::nSigns];
    TH1D* h_pt2nd_mu4_mu4noL1_excl[ParamsSet::nSigns];
    TH2D* h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_excl[ParamsSet::nSigns];
    TH2D* h_pair_eta_vs_pair_pT_mu4_mu4noL1_excl[ParamsSet::nSigns];
    TH2D* h_Deta_Dphi_mu4_mu4noL1_excl[ParamsSet::nSigns];
    TH2D* h_eta1_eta2_mu4_mu4noL1_excl[ParamsSet::nSigns];
    TH2D* h_eta_avg_Deta_mu4_mu4noL1_excl[ParamsSet::nSigns];
    TH2D* h_eta_avg_Dphi_mu4_mu4noL1_excl[ParamsSet::nSigns];
    TH1D* h_minv_zoomin_mu4_mu4noL1_excl[ParamsSet::nSigns];
    TH1D* h_pair_pt_log_mu4_mu4noL1_excl[ParamsSet::nSigns];
    TH2D* h_minv_pair_pt_log_mu4_mu4noL1_excl[ParamsSet::nSigns];

    TH1D* h_Deta_2mu4[ParamsSet::nSigns];
    TH1D* h_Deta_zoomin_2mu4[ParamsSet::nSigns];
    TH1D* h_Dphi_2mu4[ParamsSet::nSigns];
    TH1D* h_Dphi_zoomin_2mu4[ParamsSet::nSigns];
    TH1D* h_DR_2mu4[ParamsSet::nSigns];
    TH1D* h_DR_zoomin_2mu4[ParamsSet::nSigns];
    TH1D* h_pt2nd_2mu4[ParamsSet::nSigns];
    TH2D* h_pt2nd_vs_q_eta_2nd_2mu4[ParamsSet::nSigns];
    TH2D* h_pair_eta_vs_pair_pT_2mu4[ParamsSet::nSigns];
    TH2D* h_Deta_Dphi_2mu4[ParamsSet::nSigns];
    TH2D* h_eta1_eta2_2mu4[ParamsSet::nSigns];
    TH2D* h_eta_avg_Deta_2mu4[ParamsSet::nSigns];
    TH2D* h_eta_avg_Dphi_2mu4[ParamsSet::nSigns];
    TH1D* h_minv_zoomin_2mu4[ParamsSet::nSigns];
    TH1D* h_pair_pt_log_2mu4[ParamsSet::nSigns];
    TH2D* h_minv_pair_pt_log_2mu4[ParamsSet::nSigns];

    // separated (dR > 0.8)
    TH1D* h_Deta_mu4_sepr[ParamsSet::nSigns];
    TH1D* h_Deta_zoomin_mu4_sepr[ParamsSet::nSigns];
    TH1D* h_Dphi_mu4_sepr[ParamsSet::nSigns];
    TH1D* h_Dphi_zoomin_mu4_sepr[ParamsSet::nSigns];
    TH1D* h_DR_mu4_sepr[ParamsSet::nSigns];
    TH1D* h_DR_zoomin_mu4_sepr[ParamsSet::nSigns];
    TH1D* h_pt2nd_mu4_sepr[ParamsSet::nSigns];
    TH2D* h_pt2nd_vs_q_eta_2nd_mu4_sepr[ParamsSet::nSigns];
    TH2D* h_pair_eta_vs_pair_pT_mu4_sepr[ParamsSet::nSigns];
    TH2D* h_Deta_Dphi_mu4_sepr[ParamsSet::nSigns];
    TH2D* h_eta1_eta2_mu4_sepr[ParamsSet::nSigns];
    TH2D* h_eta_avg_Deta_mu4_sepr[ParamsSet::nSigns];
    TH2D* h_eta_avg_Dphi_mu4_sepr[ParamsSet::nSigns];
    TH1D* h_minv_zoomin_mu4_sepr[ParamsSet::nSigns];
    TH1D* h_pair_pt_log_mu4_sepr[ParamsSet::nSigns];
    TH2D* h_minv_pair_pt_log_mu4_sepr[ParamsSet::nSigns];

    TH1D* h_Deta_mu4_mu4noL1_sepr[ParamsSet::nSigns];
    TH1D* h_Deta_zoomin_mu4_mu4noL1_sepr[ParamsSet::nSigns];
    TH1D* h_Dphi_mu4_mu4noL1_sepr[ParamsSet::nSigns];
    TH1D* h_Dphi_zoomin_mu4_mu4noL1_sepr[ParamsSet::nSigns];
    TH1D* h_DR_mu4_mu4noL1_sepr[ParamsSet::nSigns];
    TH1D* h_DR_zoomin_mu4_mu4noL1_sepr[ParamsSet::nSigns];
    TH1D* h_pt2nd_mu4_mu4noL1_sepr[ParamsSet::nSigns];
    TH2D* h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sepr[ParamsSet::nSigns];
    TH2D* h_pair_eta_vs_pair_pT_mu4_mu4noL1_sepr[ParamsSet::nSigns];
    TH2D* h_Deta_Dphi_mu4_mu4noL1_sepr[ParamsSet::nSigns];
    TH2D* h_eta1_eta2_mu4_mu4noL1_sepr[ParamsSet::nSigns];
    TH2D* h_eta_avg_Deta_mu4_mu4noL1_sepr[ParamsSet::nSigns];
    TH2D* h_eta_avg_Dphi_mu4_mu4noL1_sepr[ParamsSet::nSigns];
    TH1D* h_minv_zoomin_mu4_mu4noL1_sepr[ParamsSet::nSigns];
    TH1D* h_pair_pt_log_mu4_mu4noL1_sepr[ParamsSet::nSigns];
    TH2D* h_minv_pair_pt_log_mu4_mu4noL1_sepr[ParamsSet::nSigns];

    TH1D* h_Deta_mu4_mu4noL1_excl_sepr[ParamsSet::nSigns];
    TH1D* h_Deta_zoomin_mu4_mu4noL1_excl_sepr[ParamsSet::nSigns];
    TH1D* h_Dphi_mu4_mu4noL1_excl_sepr[ParamsSet::nSigns];
    TH1D* h_Dphi_zoomin_mu4_mu4noL1_excl_sepr[ParamsSet::nSigns];
    TH1D* h_DR_mu4_mu4noL1_excl_sepr[ParamsSet::nSigns];
    TH1D* h_DR_zoomin_mu4_mu4noL1_excl_sepr[ParamsSet::nSigns];
    TH1D* h_pt2nd_mu4_mu4noL1_excl_sepr[ParamsSet::nSigns];
    TH2D* h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_excl_sepr[ParamsSet::nSigns];
    TH2D* h_pair_eta_vs_pair_pT_mu4_mu4noL1_excl_sepr[ParamsSet::nSigns];
    TH2D* h_Deta_Dphi_mu4_mu4noL1_excl_sepr[ParamsSet::nSigns];
    TH2D* h_eta1_eta2_mu4_mu4noL1_excl_sepr[ParamsSet::nSigns];
    TH2D* h_eta_avg_Deta_mu4_mu4noL1_excl_sepr[ParamsSet::nSigns];
    TH2D* h_eta_avg_Dphi_mu4_mu4noL1_excl_sepr[ParamsSet::nSigns];
    TH1D* h_minv_zoomin_mu4_mu4noL1_excl_sepr[ParamsSet::nSigns];
    TH1D* h_pair_pt_log_mu4_mu4noL1_excl_sepr[ParamsSet::nSigns];
    TH2D* h_minv_pair_pt_log_mu4_mu4noL1_excl_sepr[ParamsSet::nSigns];

    TH1D* h_Deta_2mu4_sepr[ParamsSet::nSigns];
    TH1D* h_Deta_zoomin_2mu4_sepr[ParamsSet::nSigns];
    TH1D* h_Dphi_2mu4_sepr[ParamsSet::nSigns];
    TH1D* h_Dphi_zoomin_2mu4_sepr[ParamsSet::nSigns];
    TH1D* h_DR_2mu4_sepr[ParamsSet::nSigns];
    TH1D* h_DR_zoomin_2mu4_sepr[ParamsSet::nSigns];
    TH1D* h_pt2nd_2mu4_sepr[ParamsSet::nSigns];
    TH2D* h_pt2nd_vs_q_eta_2nd_2mu4_sepr[ParamsSet::nSigns];
    TH2D* h_pair_eta_vs_pair_pT_2mu4_sepr[ParamsSet::nSigns];
    TH2D* h_Deta_Dphi_2mu4_sepr[ParamsSet::nSigns];
    TH2D* h_eta1_eta2_2mu4_sepr[ParamsSet::nSigns];
    TH2D* h_eta_avg_Deta_2mu4_sepr[ParamsSet::nSigns];
    TH2D* h_eta_avg_Dphi_2mu4_sepr[ParamsSet::nSigns];
    TH1D* h_minv_zoomin_2mu4_sepr[ParamsSet::nSigns];
    TH1D* h_pair_pt_log_2mu4_sepr[ParamsSet::nSigns];
    TH2D* h_minv_pair_pt_log_2mu4_sepr[ParamsSet::nSigns];

    // with signal selection 

    TH1D* h_Deta_mu4_w_single_b_sig_sel;
    TH1D* h_Deta_zoomin_mu4_w_single_b_sig_sel;
    TH1D* h_Dphi_mu4_w_single_b_sig_sel;
    TH1D* h_Dphi_zoomin_mu4_w_single_b_sig_sel;
    TH1D* h_DR_mu4_w_single_b_sig_sel;
    TH1D* h_DR_zoomin_mu4_w_single_b_sig_sel;
    TH1D* h_pt2nd_mu4_w_single_b_sig_sel;
    TH2D* h_pt2nd_vs_q_eta_2nd_mu4_w_single_b_sig_sel;
    TH2D* h_pair_eta_vs_pair_pT_mu4_w_single_b_sig_sel;
    TH2D* h_Deta_Dphi_mu4_w_single_b_sig_sel;
    TH2D* h_eta1_eta2_mu4_w_single_b_sig_sel;
    TH2D* h_eta_avg_Deta_mu4_w_single_b_sig_sel;
    TH2D* h_eta_avg_Dphi_mu4_w_single_b_sig_sel;
    TH1D* h_minv_zoomin_mu4_w_single_b_sig_sel;
    TH1D* h_pair_pt_log_mu4_w_single_b_sig_sel;
    TH2D* h_minv_pair_pt_log_mu4_w_single_b_sig_sel;

    TH1D* h_Deta_mu4_mu4noL1_w_single_b_sig_sel;
    TH1D* h_Deta_zoomin_mu4_mu4noL1_w_single_b_sig_sel;
    TH1D* h_Dphi_mu4_mu4noL1_w_single_b_sig_sel;
    TH1D* h_Dphi_zoomin_mu4_mu4noL1_w_single_b_sig_sel;
    TH1D* h_DR_mu4_mu4noL1_w_single_b_sig_sel;
    TH1D* h_DR_zoomin_mu4_mu4noL1_w_single_b_sig_sel;
    TH1D* h_pt2nd_mu4_mu4noL1_w_single_b_sig_sel;
    TH2D* h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_w_single_b_sig_sel;
    TH2D* h_pair_eta_vs_pair_pT_mu4_mu4noL1_w_single_b_sig_sel;
    TH2D* h_Deta_Dphi_mu4_mu4noL1_w_single_b_sig_sel;
    TH2D* h_eta1_eta2_mu4_mu4noL1_w_single_b_sig_sel;
    TH2D* h_eta_avg_Deta_mu4_mu4noL1_w_single_b_sig_sel;
    TH2D* h_eta_avg_Dphi_mu4_mu4noL1_w_single_b_sig_sel;
    TH1D* h_minv_zoomin_mu4_mu4noL1_w_single_b_sig_sel;
    TH1D* h_pair_pt_log_mu4_mu4noL1_w_single_b_sig_sel;
    TH2D* h_minv_pair_pt_log_mu4_mu4noL1_w_single_b_sig_sel;

    TH1D* h_Deta_mu4_mu4noL1_excl_w_single_b_sig_sel;
    TH1D* h_Deta_zoomin_mu4_mu4noL1_excl_w_single_b_sig_sel;
    TH1D* h_Dphi_mu4_mu4noL1_excl_w_single_b_sig_sel;
    TH1D* h_Dphi_zoomin_mu4_mu4noL1_excl_w_single_b_sig_sel;
    TH1D* h_DR_mu4_mu4noL1_excl_w_single_b_sig_sel;
    TH1D* h_DR_zoomin_mu4_mu4noL1_excl_w_single_b_sig_sel;
    TH1D* h_pt2nd_mu4_mu4noL1_excl_w_single_b_sig_sel;
    TH2D* h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_excl_w_single_b_sig_sel;
    TH2D* h_pair_eta_vs_pair_pT_mu4_mu4noL1_excl_w_single_b_sig_sel;
    TH2D* h_Deta_Dphi_mu4_mu4noL1_excl_w_single_b_sig_sel;
    TH2D* h_eta1_eta2_mu4_mu4noL1_excl_w_single_b_sig_sel;
    TH2D* h_eta_avg_Deta_mu4_mu4noL1_excl_w_single_b_sig_sel;
    TH2D* h_eta_avg_Dphi_mu4_mu4noL1_excl_w_single_b_sig_sel;
    TH1D* h_minv_zoomin_mu4_mu4noL1_excl_w_single_b_sig_sel;
    TH1D* h_pair_pt_log_mu4_mu4noL1_excl_w_single_b_sig_sel;
    TH2D* h_minv_pair_pt_log_mu4_mu4noL1_excl_w_single_b_sig_sel;

    TH1D* h_Deta_2mu4_w_single_b_sig_sel;
    TH1D* h_Deta_zoomin_2mu4_w_single_b_sig_sel;
    TH1D* h_Dphi_2mu4_w_single_b_sig_sel;
    TH1D* h_Dphi_zoomin_2mu4_w_single_b_sig_sel;
    TH1D* h_DR_2mu4_w_single_b_sig_sel;
    TH1D* h_DR_zoomin_2mu4_w_single_b_sig_sel;
    TH1D* h_pt2nd_2mu4_w_single_b_sig_sel;
    TH2D* h_pt2nd_vs_q_eta_2nd_2mu4_w_single_b_sig_sel;
    TH2D* h_pair_eta_vs_pair_pT_2mu4_w_single_b_sig_sel;
    TH2D* h_Deta_Dphi_2mu4_w_single_b_sig_sel;
    TH2D* h_eta1_eta2_2mu4_w_single_b_sig_sel;
    TH2D* h_eta_avg_Deta_2mu4_w_single_b_sig_sel;
    TH2D* h_eta_avg_Dphi_2mu4_w_single_b_sig_sel;
    TH1D* h_minv_zoomin_2mu4_w_single_b_sig_sel;
    TH1D* h_pair_pt_log_2mu4_w_single_b_sig_sel;
    TH2D* h_minv_pair_pt_log_2mu4_w_single_b_sig_sel;

    TH2D* h_DR_zoomin_vs_pair_pT_mu4_w_single_b_sig_sel;

    // ------------------- trigger efficiency vs pair pT or pT1st histograms (full range) -------------------
        
    // pair kinematics requiring only single mu4 trigger
    TH2D* h_Deta_vs_pT_1st_mu4[ParamsSet::nSigns];
    TH2D* h_Deta_zoomin_vs_pT_1st_mu4[ParamsSet::nSigns];
    TH2D* h_Dphi_vs_pT_1st_mu4[ParamsSet::nSigns];
    TH2D* h_Dphi_zoomin_vs_pT_1st_mu4[ParamsSet::nSigns];
    TH2D* h_DR_vs_pT_1st_mu4[ParamsSet::nSigns];
    TH2D* h_DR_zoomin_vs_pT_1st_mu4[ParamsSet::nSigns];
    TH2D* h_minv_zoomin_vs_pT_1st_mu4[ParamsSet::nSigns];
    TH2D* h_pair_pt_log_vs_pT_1st_mu4[ParamsSet::nSigns];
    TH2D* h_pt2nd_vs_pT_1st_mu4[ParamsSet::nSigns];

    TH2D* h_Deta_vs_pair_pT_mu4[ParamsSet::nSigns];
    TH2D* h_Deta_zoomin_vs_pair_pT_mu4[ParamsSet::nSigns];
    TH2D* h_Dphi_vs_pair_pT_mu4[ParamsSet::nSigns];
    TH2D* h_Dphi_zoomin_vs_pair_pT_mu4[ParamsSet::nSigns];
    TH2D* h_DR_vs_pair_pT_mu4[ParamsSet::nSigns];
    TH2D* h_DR_zoomin_vs_pair_pT_mu4[ParamsSet::nSigns];
    TH2D* h_minv_zoomin_vs_pair_pT_mu4[ParamsSet::nSigns];
    TH2D* h_pair_pt_log_vs_pair_pT_mu4[ParamsSet::nSigns];
    TH2D* h_pt2nd_vs_pair_pT_mu4[ParamsSet::nSigns];

    // pair kinematics requiring mu4_mu4noL1 trigger
    TH2D* h_Deta_vs_pT_1st_mu4_mu4noL1[ParamsSet::nSigns];
    TH2D* h_Deta_zoomin_vs_pT_1st_mu4_mu4noL1[ParamsSet::nSigns];
    TH2D* h_Dphi_vs_pT_1st_mu4_mu4noL1[ParamsSet::nSigns];
    TH2D* h_Dphi_zoomin_vs_pT_1st_mu4_mu4noL1[ParamsSet::nSigns];
    TH2D* h_DR_vs_pT_1st_mu4_mu4noL1[ParamsSet::nSigns];
    TH2D* h_DR_zoomin_vs_pT_1st_mu4_mu4noL1[ParamsSet::nSigns];
    TH2D* h_minv_zoomin_vs_pT_1st_mu4_mu4noL1[ParamsSet::nSigns];
    TH2D* h_pair_pt_log_vs_pT_1st_mu4_mu4noL1[ParamsSet::nSigns];
    TH2D* h_pt2nd_vs_pT_1st_mu4_mu4noL1[ParamsSet::nSigns];

    TH2D* h_Deta_vs_pair_pT_mu4_mu4noL1[ParamsSet::nSigns];
    TH2D* h_Deta_zoomin_vs_pair_pT_mu4_mu4noL1[ParamsSet::nSigns];
    TH2D* h_Dphi_vs_pair_pT_mu4_mu4noL1[ParamsSet::nSigns];
    TH2D* h_Dphi_zoomin_vs_pair_pT_mu4_mu4noL1[ParamsSet::nSigns];
    TH2D* h_DR_vs_pair_pT_mu4_mu4noL1[ParamsSet::nSigns];
    TH2D* h_DR_zoomin_vs_pair_pT_mu4_mu4noL1[ParamsSet::nSigns];
    TH2D* h_minv_zoomin_vs_pair_pT_mu4_mu4noL1[ParamsSet::nSigns];
    TH2D* h_pair_pt_log_vs_pair_pT_mu4_mu4noL1[ParamsSet::nSigns];
    TH2D* h_pt2nd_vs_pair_pT_mu4_mu4noL1[ParamsSet::nSigns];

    TH2D* h_Deta_vs_pT_1st_mu4_mu4noL1_excl[ParamsSet::nSigns];
    TH2D* h_Deta_zoomin_vs_pT_1st_mu4_mu4noL1_excl[ParamsSet::nSigns];
    TH2D* h_Dphi_vs_pT_1st_mu4_mu4noL1_excl[ParamsSet::nSigns];
    TH2D* h_Dphi_zoomin_vs_pT_1st_mu4_mu4noL1_excl[ParamsSet::nSigns];
    TH2D* h_DR_vs_pT_1st_mu4_mu4noL1_excl[ParamsSet::nSigns];
    TH2D* h_DR_zoomin_vs_pT_1st_mu4_mu4noL1_excl[ParamsSet::nSigns];
    TH2D* h_minv_zoomin_vs_pT_1st_mu4_mu4noL1_excl[ParamsSet::nSigns];
    TH2D* h_pair_pt_log_vs_pT_1st_mu4_mu4noL1_excl[ParamsSet::nSigns];
    TH2D* h_pt2nd_vs_pT_1st_mu4_mu4noL1_excl[ParamsSet::nSigns];

    TH2D* h_Deta_vs_pair_pT_mu4_mu4noL1_excl[ParamsSet::nSigns];
    TH2D* h_Deta_zoomin_vs_pair_pT_mu4_mu4noL1_excl[ParamsSet::nSigns];
    TH2D* h_Dphi_vs_pair_pT_mu4_mu4noL1_excl[ParamsSet::nSigns];
    TH2D* h_Dphi_zoomin_vs_pair_pT_mu4_mu4noL1_excl[ParamsSet::nSigns];
    TH2D* h_DR_vs_pair_pT_mu4_mu4noL1_excl[ParamsSet::nSigns];
    TH2D* h_DR_zoomin_vs_pair_pT_mu4_mu4noL1_excl[ParamsSet::nSigns];
    TH2D* h_minv_zoomin_vs_pair_pT_mu4_mu4noL1_excl[ParamsSet::nSigns];
    TH2D* h_pair_pt_log_vs_pair_pT_mu4_mu4noL1_excl[ParamsSet::nSigns];
    TH2D* h_pt2nd_vs_pair_pT_mu4_mu4noL1_excl[ParamsSet::nSigns];

    // pair kinematics requiring 2mu4 trigger
    TH2D* h_Deta_vs_pT_1st_2mu4[ParamsSet::nSigns];
    TH2D* h_Deta_zoomin_vs_pT_1st_2mu4[ParamsSet::nSigns];
    TH2D* h_Dphi_vs_pT_1st_2mu4[ParamsSet::nSigns];
    TH2D* h_Dphi_zoomin_vs_pT_1st_2mu4[ParamsSet::nSigns];
    TH2D* h_DR_vs_pT_1st_2mu4[ParamsSet::nSigns];
    TH2D* h_DR_zoomin_vs_pT_1st_2mu4[ParamsSet::nSigns];
    TH2D* h_minv_zoomin_vs_pT_1st_2mu4[ParamsSet::nSigns];
    TH2D* h_pair_pt_log_vs_pT_1st_2mu4[ParamsSet::nSigns];
    TH2D* h_pt2nd_vs_pT_1st_2mu4[ParamsSet::nSigns];

    TH2D* h_Deta_vs_pair_pT_2mu4[ParamsSet::nSigns];
    TH2D* h_Deta_zoomin_vs_pair_pT_2mu4[ParamsSet::nSigns];
    TH2D* h_Dphi_vs_pair_pT_2mu4[ParamsSet::nSigns];
    TH2D* h_Dphi_zoomin_vs_pair_pT_2mu4[ParamsSet::nSigns];
    TH2D* h_DR_vs_pair_pT_2mu4[ParamsSet::nSigns];
    TH2D* h_DR_zoomin_vs_pair_pT_2mu4[ParamsSet::nSigns];
    TH2D* h_minv_zoomin_vs_pair_pT_2mu4[ParamsSet::nSigns];
    TH2D* h_pair_pt_log_vs_pair_pT_2mu4[ParamsSet::nSigns];
    TH2D* h_pt2nd_vs_pair_pT_2mu4[ParamsSet::nSigns];

    // --------------------- input files & trees & data for setting branches---------------------------

   	TFile *inFile;
    // TTree *inTree_ctrbin[ParamsSet::nSigns];
    // TTree *inTree[ParamsSet::nSigns];
    TTree *inTree[ParamsSet::nSigns];

    double weight[ParamsSet::nSigns];
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
    int m1charge[ParamsSet::nSigns];
    int m2charge[ParamsSet::nSigns];
    bool passmu4mu4noL1[ParamsSet::nSigns];
    bool pass2mu4[ParamsSet::nSigns];
    bool mu1PassSingle[ParamsSet::nSigns];
    bool mu2PassSingle[ParamsSet::nSigns];
    bool passSeparated[ParamsSet::nSigns];
    bool passSeparatedDeta[ParamsSet::nSigns];


    // --------------------- class methods ---------------------------

   	void InitInput();
    void InitOutput();
   	void InitHists();
   	void ProcessData();
    void WriteOutput();
    bool PassSingleMuonGapCut(float meta, float mpt, int mcharge);
   	void FillHistograms(int nsign);
    void CalculateSingleTrigEffcyRatio(TH1* h1, TH1* h2, TH1* h3);
    void CalculateTrigEffcyRatio();

public:
    bool output_non_trig_effcy_hists;
    bool isScram;
    bool isTight;
    bool isRun3;
    // bool isMu4mu4noL1;
    int trigger_mode = 1;
    // bool require_exclusive_mu4_for_mu4_mu4noL1 = false; // if true, require that only one muon fire mu4, to ensure unambiguous knowledge of which muon fires mu4noL1
    bool doTrigEffcy = true;
    bool filter_out_photo_resn_for_trig_effcy = true;
  	MuonPairPlottingPP();
  	~MuonPairPlottingPP(){}
  	void Run();

};



#endif