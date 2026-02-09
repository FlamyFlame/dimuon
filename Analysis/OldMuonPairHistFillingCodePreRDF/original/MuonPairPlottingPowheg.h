// Assumes that we exclude resonances

#ifndef __MuonPairPlottingPowheg_h__
#define __MuonPairPlottingPowheg_h__

#include <TROOT.h>
// #include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <string.h>
#include  <stdlib.h>
#include <fstream>
#include " ../MuonObjectsParamsAndHelpers/ParamsSet.h"
#include " ../MuonObjectsParamsAndHelpers/MuonPair.h"
#include "vector"
#include "TH1D.h"
#include "TH2D.h"

class MuonPairPlottingPowheg{
  	//updates histograms - centrality binned
    //if need histograms with all centrality bins added together
    //can add at the histograms level


private:
// --------------------- general settings ---------------------------

    // int mode = 1;

  	ParamsSet pms;

    double crossx_cut;
    double filter_effcy;
    double filter_effcy_bb = 0.003;
    double filter_effcy_cc = 0.001108;

    static const int nBatches = 6;
    std::string mcdir = "/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample/";
    std::string sub_dir;
    std::string mc_mode;

    // int dr_bin_start;

// --------------------- output file & histograms ---------------------------

  	TFile *outFile = nullptr;
    int nevents_before_cuts_total = 0;
    std::string dphi_regions[2] = {"near", "away"};

    // TH1D* h_pair_dP_overP_dr_binned[ParamsSet::ndRselcs][ParamsSet::nSigns];
    // TH1D* h_Dphi_dr_binned[ParamsSet::ndRselcs][ParamsSet::nSigns];
    // TH1D* h_DR_dr_binned[ParamsSet::ndRselcs][ParamsSet::nSigns];
    // TH1D* h_pair_y_dr_binned[ParamsSet::ndRselcs][ParamsSet::nSigns];
    // TH1D* h_pt_asym_dr_binned[ParamsSet::ndRselcs][ParamsSet::nSigns];
    // TH1D* h_pair_pt_ptlead_ratio_dr_binned[ParamsSet::ndRselcs][ParamsSet::nSigns];

    // TH2D* h_eta_avg_Dphi_dr_binned[ParamsSet::ndRselcs][ParamsSet::nSigns];
    // TH2D* h_Deta_Dphi_dr_binned[ParamsSet::ndRselcs][ParamsSet::nSigns];
    // TH2D* h_eta1_eta2_dr_binned[ParamsSet::ndRselcs][ParamsSet::nSigns];
    // TH2D* h_eta_avg_Deta_dr_binned[ParamsSet::ndRselcs][ParamsSet::nSigns];
    // TH2D* h_pt1_pt2_dr_binned[ParamsSet::ndRselcs][ParamsSet::nSigns];
    // TH2D* h_ptlead_pair_pt_dr_binned[ParamsSet::ndRselcs][ParamsSet::nSigns];
    // TH2D* h_minv_pair_pt_dr_binned[ParamsSet::ndRselcs][ParamsSet::nSigns];

    // TH2D* h_unweighted_Deta_Dphi_dr_binned[ParamsSet::ndRselcs][ParamsSet::nSigns];
    
    TH1D* h_mQQ[2][ParamsSet::nSigns];
    TH1D* h_mQQ_Q_ratio[2][ParamsSet::nSigns];
    TH1D* h_mQQ_mHard_ratio[2][ParamsSet::nSigns];

    TH1D* h_pair_dP_overP[2][ParamsSet::nSigns];
    TH1D* h_Dphi[2][ParamsSet::nSigns];
    TH1D* h_DR[2][ParamsSet::nSigns];
    TH1D* h_DR_zoomin[2][ParamsSet::nSigns];
    TH1D* h_DR_jacobian_corrected[2][ParamsSet::nSigns];
    TH1D* h_DR_zoomin_jacobian_corrected[2][ParamsSet::nSigns];
    TH1D* h_pair_y[2][ParamsSet::nSigns];
    TH1D* h_pt_asym[2][ParamsSet::nSigns];
    TH1D* h_pair_pt_ptlead_ratio[2][ParamsSet::nSigns];
    TH2D* h_eta_avg_Dphi[2][ParamsSet::nSigns];
    TH2D* h_Deta_Dphi[2][ParamsSet::nSigns];
    TH2D* h_eta1_eta2[2][ParamsSet::nSigns];
    TH2D* h_eta_avg_Deta[2][ParamsSet::nSigns];
    TH2D* h_pt1_pt2[2][ParamsSet::nSigns];
    TH2D* h_ptlead_pair_pt[2][ParamsSet::nSigns];
    TH2D* h_ptlead_pair_pt_zoomin[2][ParamsSet::nSigns];
    TH2D* h_ptlead_pair_pt_log[2][ParamsSet::nSigns];
    TH2D* h_minv_pair_pt[2][ParamsSet::nSigns];
    TH2D* h_minv_pair_pt_zoomin[2][ParamsSet::nSigns];
    TH2D* h_minv_pair_pt_log[2][ParamsSet::nSigns];

    static const int nFlavors = 4; // single b, bb, cc, others
    static const int nAncestorGroups = 4;

    TH1D* h_mQQ_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 2];
    TH1D* h_mQQ_Q_ratio_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 2];
    TH1D* h_mQQ_mHard_ratio_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 2];

    TH1D* h_DR_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 2];
    TH1D* h_DR_zoomin_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 2];
    TH1D* h_DR_jacobian_corrected_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 2];
    TH1D* h_DR_zoomin_jacobian_corrected_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 2];
    TH1D* h_pt_asym_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 2];
    // TH1D* h_psrapidity_ordered_pt_asym_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 2];
    TH1D* h_pair_pt_ptlead_ratio_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 2];
    TH2D* h_ptlead_pair_pt_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 2];
    TH2D* h_Deta_Dphi_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 2];
    TH2D* h_minv_pair_pt_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 2];
    TH2D* h_minv_pair_pt_zoomin_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 2];
    TH2D* h_minv_pair_pt_jacobian_corrected_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 2];
    TH2D* h_minv_pair_pt_zoomin_jacobian_corrected_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 2];

    TH1D* h_DR_flavor_binned[ParamsSet::nSigns][nFlavors];
    TH1D* h_DR_zoomin_flavor_binned[ParamsSet::nSigns][nFlavors];
    TH1D* h_DR_jacobian_corrected_flavor_binned[ParamsSet::nSigns][nFlavors];
    TH1D* h_DR_zoomin_jacobian_corrected_flavor_binned[ParamsSet::nSigns][nFlavors];
    TH1D* h_pt_asym_flavor_binned[ParamsSet::nSigns][nFlavors];
    // TH1D* h_psrapidity_ordered_pt_asym_flavor_binned[ParamsSet::nSigns][nFlavors];
    TH1D* h_pair_pt_ptlead_ratio_flavor_binned[ParamsSet::nSigns][nFlavors];
    TH2D* h_ptlead_pair_pt_flavor_binned[ParamsSet::nSigns][nFlavors];
    TH2D* h_Deta_Dphi_flavor_binned[ParamsSet::nSigns][nFlavors];
    TH2D* h_minv_pair_pt_flavor_binned[ParamsSet::nSigns][nFlavors];
    TH2D* h_minv_pair_pt_zoomin_flavor_binned[ParamsSet::nSigns][nFlavors];
    TH2D* h_minv_pair_pt_jacobian_corrected_flavor_binned[ParamsSet::nSigns][nFlavors];
    TH2D* h_minv_pair_pt_zoomin_jacobian_corrected_flavor_binned[ParamsSet::nSigns][nFlavors];

    // TH1D* h_mQQ_flavor_binned[ParamsSet::nSigns][nFlavors];
    // TH1D* h_mQQ_Q_ratio_flavor_binned[ParamsSet::nSigns][nFlavors];
    // TH1D* h_mQQ_mHard_ratio_flavor_binned[ParamsSet::nSigns][nFlavors];


    // TH2D* h_unweighted_eta1_eta2_dphicut[ParamsSet::ndphiselcs][ParamsSet::nSigns];

// --------------------- input files & trees & data for setting branches---------------------------

   	TFile *inFile[nBatches];
    // TTree *inTree_ctrbin[nBatches][ParamsSet::nSigns];
    // TTree *inTree[ParamsSet::nSigns];
    TTree *meta_tree[nBatches];
    TTree *inTree[nBatches][ParamsSet::nSigns];

    Long64_t nevents_before_cuts[nBatches];
    double weight[nBatches][ParamsSet::nSigns];
    float Q[nBatches][ParamsSet::nSigns];
    float mQQ[nBatches][ParamsSet::nSigns];
    float mHard_relevant[nBatches][ParamsSet::nSigns];
    bool from_same_b[nBatches][ParamsSet::nSigns];
    bool both_from_b[nBatches][ParamsSet::nSigns];
    bool both_from_c[nBatches][ParamsSet::nSigns];
    bool from_same_ancestors[nBatches][ParamsSet::nSigns];
    int m1_ancestor_category[nBatches][ParamsSet::nSigns];
    int m2_ancestor_category[nBatches][ParamsSet::nSigns];
  	float pair_dPoverP[nBatches][ParamsSet::nSigns];
  	float pt_lead[nBatches][ParamsSet::nSigns];
  	float pair_pt[nBatches][ParamsSet::nSigns];
  	// float pair_eta[nBatches][ParamsSet::nSigns];
  	float pair_y[nBatches][ParamsSet::nSigns];
    float asym[nBatches][ParamsSet::nSigns];
  	float dpt[nBatches][ParamsSet::nSigns];
  	float deta[nBatches][ParamsSet::nSigns];
  	float etaavg[nBatches][ParamsSet::nSigns];
  	float phiavg[nBatches][ParamsSet::nSigns];
  	float dphi[nBatches][ParamsSet::nSigns];
  	float dr[nBatches][ParamsSet::nSigns];
  	float minv[nBatches][ParamsSet::nSigns];
  	float m1pt[nBatches][ParamsSet::nSigns];
  	float m2pt[nBatches][ParamsSet::nSigns];
  	float m1eta[nBatches][ParamsSet::nSigns];
  	float m2eta[nBatches][ParamsSet::nSigns];
  	float m1phi[nBatches][ParamsSet::nSigns];
  	float m2phi[nBatches][ParamsSet::nSigns];
    int m1charge[nBatches][ParamsSet::nSigns];
    int m2charge[nBatches][ParamsSet::nSigns];


// --------------------- class methods ---------------------------

   	void InitInput();
    void InitOutput();
   	void InitHists();
   	void ProcessData();
    void WriteOutput();
    bool PassSingleMuonGapCut(float meta, float mpt, int mcharge);
   	void FillHistograms(int ndr, int nsign);
   	void FillPtBinnedHistograms(int ndr, int npt, int nsign);

public:
    int mode = 1;
    bool isScram;
    bool isTight;
    bool isMCTruthBB;
    bool isMCTruthCC;
    // bool saveDRbinned;
  	MuonPairPlottingPowheg(){
        crossx_cut = 5 * pow(10,8);
    }
  	~MuonPairPlottingPowheg(){}
  	void Run();

};

void MuonPairPlottingPowheg::InitInput(){

    // if (saveDRbinned) dr_bin_start = 0;
    
    for (int ibatch = 0; ibatch < nBatches; ibatch++){
        char infilename[100];

        std::sprintf(infilename, "%s%smuon_pairs_%s_%d-%d.root", mcdir.c_str(), sub_dir.c_str(), mc_mode.c_str(), 5 * ibatch + 1, 5 * ibatch + 5);
        // std::cout << infilename << std::endl;
        std::ifstream in_file(infilename);
        
        if (!in_file.good()){
            std::cout << "Warning: File " << infilename << " not found. Skip.\n";
            continue; // skip this file
        }
        inFile[ibatch] = new TFile(infilename,"read");
        if (!inFile[ibatch] || inFile[ibatch]->IsZombie()) {
            std::cout << "the root file at " << infilename << " does not give a valid TFile" << std::endl;
            throw std::exception();
        }
    }

    for (int ibatch = 0; ibatch < nBatches; ibatch++){
        meta_tree[ibatch] = (TTree*) inFile[ibatch]->Get("meta_tree");  
        meta_tree[ibatch] ->SetBranchAddress("nentries_before_cuts"          , &nevents_before_cuts[ibatch]);
        for (int isign = 0; isign < ParamsSet::nSigns; isign++){
            inTree[ibatch][isign]  = (TTree*) inFile[ibatch]->Get(Form("muon_pair_tree_dr3_sign%d",isign+1));  
            inTree[ibatch][isign] ->SetBranchAddress("weight"          , &weight[ibatch][isign]);
            inTree[ibatch][isign] ->SetBranchAddress("Q"          , &Q[ibatch][isign]);
            inTree[ibatch][isign] ->SetBranchAddress("mQQ"          , &mQQ[ibatch][isign]);
            inTree[ibatch][isign] ->SetBranchAddress("mHard_relevant"          , &mHard_relevant[ibatch][isign]);
            inTree[ibatch][isign] ->SetBranchAddress("from_same_b"          , &from_same_b[ibatch][isign]);
            inTree[ibatch][isign] ->SetBranchAddress("both_from_b"          , &both_from_b[ibatch][isign]);
            inTree[ibatch][isign] ->SetBranchAddress("both_from_c"          , &both_from_c[ibatch][isign]);
            inTree[ibatch][isign] ->SetBranchAddress("from_same_ancestors"          , &from_same_ancestors[ibatch][isign]);
            inTree[ibatch][isign] ->SetBranchAddress("m1_ancestor_category"          , &m1_ancestor_category[ibatch][isign]);
            inTree[ibatch][isign] ->SetBranchAddress("m2_ancestor_category"          , &m2_ancestor_category[ibatch][isign]);
        	inTree[ibatch][isign] ->SetBranchAddress("pair_dPoverP"    , &pair_dPoverP[ibatch][isign]);
        	inTree[ibatch][isign] ->SetBranchAddress("pt_lead"          , &pt_lead[ibatch][isign]);
        	inTree[ibatch][isign] ->SetBranchAddress("pair_pt"          , &pair_pt[ibatch][isign]);
        	// inTree[ibatch][isign] ->SetBranchAddress("pair_eta"     , &pair_eta[ibatch][isign]);
        	inTree[ibatch][isign] ->SetBranchAddress("pair_y"           , &pair_y[ibatch][isign]);
            inTree[ibatch][isign] ->SetBranchAddress("asym"             , &asym[ibatch][isign]);
        	inTree[ibatch][isign] ->SetBranchAddress("dpt"           , &dpt[ibatch][isign]);
        	inTree[ibatch][isign] ->SetBranchAddress("deta"       , &deta[ibatch][isign]);
        	inTree[ibatch][isign] ->SetBranchAddress("etaavg"      , &etaavg[ibatch][isign]);
        	inTree[ibatch][isign] ->SetBranchAddress("phiavg"            , &phiavg[ibatch][isign]);
        	inTree[ibatch][isign] ->SetBranchAddress("dphi"     , &dphi[ibatch][isign]);
        	inTree[ibatch][isign] ->SetBranchAddress("dr"        , &dr[ibatch][isign]);
        	inTree[ibatch][isign] ->SetBranchAddress("minv"        , &minv[ibatch][isign]);
        	inTree[ibatch][isign] ->SetBranchAddress("m1.pt"           , &m1pt[ibatch][isign]);
        	inTree[ibatch][isign] ->SetBranchAddress("m2.pt"           , &m2pt[ibatch][isign]);
        	inTree[ibatch][isign] ->SetBranchAddress("m1.eta"       , &m1eta[ibatch][isign]);
        	inTree[ibatch][isign] ->SetBranchAddress("m2.eta"       , &m2eta[ibatch][isign]);
        	inTree[ibatch][isign] ->SetBranchAddress("m1.phi"     	, &m1phi[ibatch][isign]);
        	inTree[ibatch][isign] ->SetBranchAddress("m2.phi"     	, &m2phi[ibatch][isign]);
            inTree[ibatch][isign] ->SetBranchAddress("m1.charge"           , &m1charge[ibatch][isign]);
            inTree[ibatch][isign] ->SetBranchAddress("m2.charge"           , &m2charge[ibatch][isign]);

            inTree[ibatch][isign] ->SetBranchStatus("*"                  ,0);//switch off all branches, then enable just the ones that we need
            inTree[ibatch][isign] ->SetBranchStatus("weight"           ,1);
            inTree[ibatch][isign] ->SetBranchStatus("Q"           ,1);
            inTree[ibatch][isign] ->SetBranchStatus("mQQ"           ,1);
            inTree[ibatch][isign] ->SetBranchStatus("mHard_relevant"           ,1);
            inTree[ibatch][isign] ->SetBranchStatus("from_same_b"           ,1);
            inTree[ibatch][isign] ->SetBranchStatus("both_from_b"           ,1);
            inTree[ibatch][isign] ->SetBranchStatus("both_from_c"           ,1);
            inTree[ibatch][isign] ->SetBranchStatus("from_same_ancestors"           ,1);
            inTree[ibatch][isign] ->SetBranchStatus("m1_ancestor_category"           ,1);
            inTree[ibatch][isign] ->SetBranchStatus("m2_ancestor_category"           ,1);
            inTree[ibatch][isign] ->SetBranchStatus("pair_dPoverP"           ,1);
            inTree[ibatch][isign] ->SetBranchStatus("pt_lead"           ,1);
            inTree[ibatch][isign] ->SetBranchStatus("pair_pt"           ,1);
            inTree[ibatch][isign] ->SetBranchStatus("pair_y"           ,1);
            inTree[ibatch][isign] ->SetBranchStatus("asym"           ,1);
            inTree[ibatch][isign] ->SetBranchStatus("dpt"           ,1);
            inTree[ibatch][isign] ->SetBranchStatus("deta"           ,1);
            inTree[ibatch][isign] ->SetBranchStatus("etaavg"           ,1);
            inTree[ibatch][isign] ->SetBranchStatus("phiavg"           ,1);
            inTree[ibatch][isign] ->SetBranchStatus("dphi"           ,1);
            inTree[ibatch][isign] ->SetBranchStatus("dr"           ,1);
            inTree[ibatch][isign] ->SetBranchStatus("minv"           ,1);
            inTree[ibatch][isign] ->SetBranchStatus("m1.pt"           ,1);
            inTree[ibatch][isign] ->SetBranchStatus("m2.pt"           ,1);
            inTree[ibatch][isign] ->SetBranchStatus("m1.eta"           ,1);
            inTree[ibatch][isign] ->SetBranchStatus("m2.eta"           ,1);
            inTree[ibatch][isign] ->SetBranchStatus("m1.phi"           ,1);
            inTree[ibatch][isign] ->SetBranchStatus("m2.phi"           ,1);
            inTree[ibatch][isign] ->SetBranchStatus("m1.charge"           ,1);
            inTree[ibatch][isign] ->SetBranchStatus("m2.charge"           ,1);

        }
    }
}

void MuonPairPlottingPowheg::InitOutput(){
    // sub_dir = (isMCTruthBB)? "bb_full_sample/" : "cc_full_sample/";
    outFile = new TFile(Form("%s%shistograms_%s_combined.root", mcdir.c_str(), sub_dir.c_str(), mc_mode.c_str()),"recreate");
}

void MuonPairPlottingPowheg::InitHists(){
    int nDR_bins = (isMCTruthBB || isMCTruthCC)? 100 : 200;
    int nDR_zoomin_bins = 100;
    int nDphi_bins = (isMCTruthBB || isMCTruthCC)? 64 : 128;
    int neta_bins = (isMCTruthBB || isMCTruthCC)? 50 : 100;
    int nDeta_bins = (isMCTruthBB || isMCTruthCC)? 100 : 200;
    int npair_y_bins = (isMCTruthBB || isMCTruthCC)? 45 : 90;
    int npt_asym_bins = (isMCTruthBB || isMCTruthCC)? 50 : 100;
    int npair_pt_ptlead_ratio_bins = (isMCTruthBB || isMCTruthCC)? 50 : 100;
    int nminv_bins_linear = (isMCTruthBB || isMCTruthCC)? 50 : 100;
    int npair_pT_bins_linear = (isMCTruthBB || isMCTruthCC)? 50 : 100;
    int npT_lead_bins_linear = (isMCTruthBB || isMCTruthCC)? 50 : 100;

    static const int nminv_bins_log = 40;
    float minv_logpow[ParamsSet::nSigns];
    minv_logpow[0] = 0.062;
    minv_logpow[1] = 0.045;
    float minv_max = 60;
    double minv_bins_log[ParamsSet::nSigns][nminv_bins_log+1];

    for (unsigned int isign = 0; isign < ParamsSet::nSigns; isign++){
        for(int iminv = 0; iminv <= nminv_bins_log; iminv++){
            minv_bins_log[isign][iminv] = minv_max * pow(10.0, (static_cast<float>(iminv - nminv_bins_log)) * minv_logpow[isign]);
        }
    }

    std::string ancestor_grps[nAncestorGroups + 2] = {"_gg", "_qg","_single_g","_qq", "_from_same_b", "_others"};
    std::string flavor_grps[nFlavors] = {"_single_b", "_bb", "_cc", "_other_flavors"};

   	if (mode == 1){
        for (unsigned int isign = 0; isign < ParamsSet::nSigns; isign++){
            for (int igrp = 0; igrp < nAncestorGroups + 2; igrp++){
                h_DR_ancestor_binned[isign][igrp] = new TH1D(Form("h_DR_sign%d%s",isign+1,ancestor_grps[igrp].c_str()),";#Delta R;1/N_{evt} dN/d#Delta R", nDR_bins,0,pms.deltaR_thrsh[2]);
                h_DR_zoomin_ancestor_binned[isign][igrp] = new TH1D(Form("h_DR_zoomin_sign%d%s",isign+1,ancestor_grps[igrp].c_str()),";#Delta R;1/N_{evt} dN/d#Delta R", nDR_zoomin_bins,0,pms.deltaR_thrsh_zoomin);
                h_DR_jacobian_corrected_ancestor_binned[isign][igrp] = new TH1D(Form("h_DR_jacobian_corrected_sign%d%s",isign+1,ancestor_grps[igrp].c_str()),";#Delta R;#frac{1}{N_{evt}} #frac{1}{#Delta R} #frac{dN}{d#Delta R}", nDR_bins,0,pms.deltaR_thrsh[2]);
                h_DR_zoomin_jacobian_corrected_ancestor_binned[isign][igrp] = new TH1D(Form("h_DR_zoomin_jacobian_corrected_sign%d%s",isign+1,ancestor_grps[igrp].c_str()),";#Delta R;#frac{1}{N_{evt}} #frac{1}{#Delta R} #frac{dN}{d#Delta R}", nDR_zoomin_bins,0,pms.deltaR_thrsh_zoomin);
                h_pt_asym_ancestor_binned[isign][igrp] = new TH1D(Form("h_pt_asym_sign%d%s",isign+1,ancestor_grps[igrp].c_str()),";A = (pT1 - pT2)/(pT1 + pT2);d#sigma/dA", npt_asym_bins,0,1.);
                // h_psrapidity_ordered_pt_asym_ancestor_binned
                h_pair_pt_ptlead_ratio_ancestor_binned[isign][igrp] = new TH1D(Form("h_pair_pt_ptlead_ratio_sign%d%s",isign+1,ancestor_grps[igrp].c_str()),";#frac{p_{T}^{pair}}{p_{T}^{lead}};d#sigma/d#frac{p_{T}^{pair}}{p_{T}^{lead}}", npair_pt_ptlead_ratio_bins,0,2.);
                h_ptlead_pair_pt_ancestor_binned[isign][igrp] = new TH2D(Form("h_ptlead_pair_pt_sign%d%s",isign+1,ancestor_grps[igrp].c_str()),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",npair_pT_bins_linear,0,30,npT_lead_bins_linear,0,30);
                h_Deta_Dphi_ancestor_binned[isign][igrp] = new TH2D(Form("h_Deta_Dphi_sign%d%s",isign+1,ancestor_grps[igrp].c_str()),";#Delta#phi;#Delta#eta", nDphi_bins,-pms.PI,pms.PI,nDeta_bins,-4.8,4.8);
                h_minv_pair_pt_ancestor_binned[isign][igrp] = new TH2D(Form("h_minv_pair_pt_sign%d%s",isign+1,ancestor_grps[igrp].c_str()),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",npair_pT_bins_linear,0,30,nminv_bins_linear,0,30);
                h_minv_pair_pt_zoomin_ancestor_binned[isign][igrp] = new TH2D(Form("h_minv_pair_pt_sign%d%s",isign+1,ancestor_grps[igrp].c_str()),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",npair_pT_bins_linear,0,20,nminv_bins_linear,0,4);
                h_minv_pair_pt_jacobian_corrected_ancestor_binned[isign][igrp] = new TH2D(Form("h_minv_pair_pt_sign%d%s",isign+1,ancestor_grps[igrp].c_str()),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",npair_pT_bins_linear,0,30,nminv_bins_linear,0,30);
                h_minv_pair_pt_zoomin_jacobian_corrected_ancestor_binned[isign][igrp] = new TH2D(Form("h_minv_pair_pt_sign%d%s",isign+1,ancestor_grps[igrp].c_str()),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",npair_pT_bins_linear,0,20,nminv_bins_linear,0,4);
                                
                h_mQQ_ancestor_binned[isign][igrp] = new TH1D(Form("h_mQQ_sign%d%s",isign+1,ancestor_grps[igrp].c_str()),";m_{QQ}", 200, 0, 200.);
                h_mQQ_Q_ratio_ancestor_binned[isign][igrp] = new TH1D(Form("h_mQQ_Q_ratio_sign%d%s",isign+1,ancestor_grps[igrp].c_str()),";#frac{m_{QQ}}{Q};d#sigma/d#frac{m_{QQ}}{Q}", 100,0,5.);
                h_mQQ_mHard_ratio_ancestor_binned[isign][igrp] = new TH1D(Form("h_mQQ_mHard_ratio_sign%d%s",isign+1,ancestor_grps[igrp].c_str()),";#frac{m_{QQ}}{#sqrt{#hat{s}}};d#sigma/d#frac{m_{QQ}}{#sqrt{#hat{s}}}", 50,0,1.);
                
                h_DR_ancestor_binned[isign][igrp]->Sumw2();
                h_DR_zoomin_ancestor_binned[isign][igrp]->Sumw2();
                h_DR_jacobian_corrected_ancestor_binned[isign][igrp]->Sumw2();
                h_DR_zoomin_jacobian_corrected_ancestor_binned[isign][igrp]->Sumw2();
                h_pt_asym_ancestor_binned[isign][igrp]->Sumw2();
                h_pair_pt_ptlead_ratio_ancestor_binned[isign][igrp]->Sumw2();
                h_ptlead_pair_pt_ancestor_binned[isign][igrp]->Sumw2();
                h_Deta_Dphi_ancestor_binned[isign][igrp]->Sumw2();
                h_minv_pair_pt_ancestor_binned[isign][igrp]->Sumw2();
                h_minv_pair_pt_zoomin_ancestor_binned[isign][igrp]->Sumw2();
                h_minv_pair_pt_jacobian_corrected_ancestor_binned[isign][igrp]->Sumw2();
                h_minv_pair_pt_zoomin_jacobian_corrected_ancestor_binned[isign][igrp]->Sumw2();
                h_mQQ_ancestor_binned[isign][igrp]->Sumw2();
                h_mQQ_Q_ratio_ancestor_binned[isign][igrp]->Sumw2();
                h_mQQ_mHard_ratio_ancestor_binned[isign][igrp]->Sumw2();
            }

            for (int kflav = 0; kflav < nFlavors; kflav++){
                h_DR_flavor_binned[isign][kflav] = new TH1D(Form("h_DR_sign%d%s",isign+1,flavor_grps[kflav].c_str()),";#Delta R;1/N_{evt} dN/d#Delta R", nDR_bins,0,pms.deltaR_thrsh[2]);
                h_DR_zoomin_flavor_binned[isign][kflav] = new TH1D(Form("h_DR_zoomin_sign%d%s",isign+1,flavor_grps[kflav].c_str()),";#Delta R;1/N_{evt} dN/d#Delta R", nDR_zoomin_bins,0,pms.deltaR_thrsh_zoomin);
                h_DR_jacobian_corrected_flavor_binned[isign][kflav] = new TH1D(Form("h_DR_jacobian_corrected_sign%d%s",isign+1,flavor_grps[kflav].c_str()),";#Delta R;#frac{1}{N_{evt}} #frac{1}{#Delta R} #frac{dN}{d#Delta R}", nDR_bins,0,pms.deltaR_thrsh[2]);
                h_DR_zoomin_jacobian_corrected_flavor_binned[isign][kflav] = new TH1D(Form("h_DR_zoomin_jacobian_corrected_sign%d%s",isign+1,flavor_grps[kflav].c_str()),";#Delta R;#frac{1}{N_{evt}} #frac{1}{#Delta R} #frac{dN}{d#Delta R}", nDR_zoomin_bins,0,pms.deltaR_thrsh_zoomin);
                h_pt_asym_flavor_binned[isign][kflav] = new TH1D(Form("h_pt_asym_sign%d%s",isign+1,flavor_grps[kflav].c_str()),";A = (pT1 - pT2)/(pT1 + pT2);d#sigma/dA", npt_asym_bins,0,1.);
                // h_psrapidity_ordered_pt_asym_flavor_binned
                h_pair_pt_ptlead_ratio_flavor_binned[isign][kflav] = new TH1D(Form("h_pair_pt_ptlead_ratio_sign%d%s",isign+1,flavor_grps[kflav].c_str()),";#frac{p_{T}^{pair}}{p_{T}^{lead}};d#sigma/d#frac{p_{T}^{pair}}{p_{T}^{lead}}", npair_pt_ptlead_ratio_bins,0,2.);
                h_ptlead_pair_pt_flavor_binned[isign][kflav] = new TH2D(Form("h_ptlead_pair_pt_sign%d%s",isign+1,flavor_grps[kflav].c_str()),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",npair_pT_bins_linear,0,30,npT_lead_bins_linear,0,30);
                h_Deta_Dphi_flavor_binned[isign][kflav] = new TH2D(Form("h_Deta_Dphi_sign%d%s",isign+1,flavor_grps[kflav].c_str()),";#Delta#phi;#Delta#eta", nDphi_bins,-pms.PI,pms.PI,nDeta_bins,-4.8,4.8);
                h_minv_pair_pt_flavor_binned[isign][kflav] = new TH2D(Form("h_minv_pair_pt_sign%d%s",isign+1,flavor_grps[kflav].c_str()),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",npair_pT_bins_linear,0,30,nminv_bins_linear,0,30);
                h_minv_pair_pt_zoomin_flavor_binned[isign][kflav] = new TH2D(Form("h_minv_pair_pt_sign%d%s",isign+1,flavor_grps[kflav].c_str()),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",npair_pT_bins_linear,0,20,nminv_bins_linear,0,4);
                h_minv_pair_pt_jacobian_corrected_flavor_binned[isign][kflav] = new TH2D(Form("h_minv_pair_pt_sign%d%s",isign+1,flavor_grps[kflav].c_str()),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",npair_pT_bins_linear,0,30,nminv_bins_linear,0,30);
                h_minv_pair_pt_zoomin_jacobian_corrected_flavor_binned[isign][kflav] = new TH2D(Form("h_minv_pair_pt_sign%d%s",isign+1,flavor_grps[kflav].c_str()),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",npair_pT_bins_linear,0,20,nminv_bins_linear,0,4);
                                
                // h_mQQ_flavor_binned[isign][kflav] = new TH1D(Form("h_mQQ_sign%d%s",isign+1,flavor_grps[kflav].c_str()),";m_{QQ}", 200, 0, 200.);
                // h_mQQ_Q_ratio_flavor_binned[isign][kflav] = new TH1D(Form("h_mQQ_Q_ratio_sign%d%s",isign+1,flavor_grps[kflav].c_str()),";#frac{m_{QQ}}{Q};d#sigma/d#frac{m_{QQ}}{Q}", 100,0,5.);
                // h_mQQ_mHard_ratio_flavor_binned[isign][kflav] = new TH1D(Form("h_mQQ_mHard_ratio_sign%d%s",isign+1,flavor_grps[kflav].c_str()),";#frac{m_{QQ}}{#sqrt{#hat{s}}};d#sigma/d#frac{m_{QQ}}{#sqrt{#hat{s}}}", 50,0,1.);
                
                h_DR_flavor_binned[isign][kflav]->Sumw2();
                h_pt_asym_flavor_binned[isign][kflav]->Sumw2();
                h_pair_pt_ptlead_ratio_flavor_binned[isign][kflav]->Sumw2();
                h_ptlead_pair_pt_flavor_binned[isign][kflav]->Sumw2();
                h_Deta_Dphi_flavor_binned[isign][kflav]->Sumw2();
                h_minv_pair_pt_flavor_binned[isign][kflav]->Sumw2();
                h_minv_pair_pt_zoomin_flavor_binned[isign][kflav]->Sumw2();
                h_minv_pair_pt_jacobian_corrected_flavor_binned[isign][kflav]->Sumw2();
                h_minv_pair_pt_zoomin_jacobian_corrected_flavor_binned[isign][kflav]->Sumw2();
                // h_mQQ_flavor_binned[isign][kflav]->Sumw2();
                // h_mQQ_Q_ratio_flavor_binned[isign][kflav]->Sumw2();
                // h_mQQ_mHard_ratio_flavor_binned[isign][kflav]->Sumw2();
            }

            for (unsigned int jdphi = 0; jdphi < 2; jdphi++){
                h_pair_dP_overP[jdphi][isign] = new TH1D(Form("h_pair_dP_overP_%s_sign%d",dphi_regions[jdphi].c_str(),isign+1),";(#Delta p / p)_{pair};1/N_{evt} dN/d(#Delta p / p)_{pair}" ,pms.deltaP_overP_nbins,0,pms.deltaP_overP_max);
                h_DR[jdphi][isign] = new TH1D(Form("h_DR_%s_sign%d",dphi_regions[jdphi].c_str(),isign+1),";#Delta R;1/N_{evt} dN/d#Delta R", nDR_bins,0,pms.deltaR_thrsh[2]);
                h_DR_zoomin[jdphi][isign] = new TH1D(Form("h_DR_zoomin_%s_sign%d",dphi_regions[jdphi].c_str(),isign+1),";#Delta R;1/N_{evt} dN/d#Delta R", nDR_zoomin_bins,0,pms.deltaR_thrsh_zoomin);
                h_DR_jacobian_corrected[jdphi][isign] = new TH1D(Form("h_DR_jacobian_corrected_%s_sign%d",dphi_regions[jdphi].c_str(),isign+1),";#Delta R;#frac{1}{N_{evt}} #frac{1}{#Delta R} #frac{dN}{d#Delta R}", nDR_bins,0,pms.deltaR_thrsh[2]);
                h_DR_zoomin_jacobian_corrected[jdphi][isign] = new TH1D(Form("h_DR_zoomin_jacobian_corrected_%s_sign%d",dphi_regions[jdphi].c_str(),isign+1),";#Delta R;#frac{1}{N_{evt}} #frac{1}{#Delta R} #frac{dN}{d#Delta R}", nDR_zoomin_bins,0,pms.deltaR_thrsh_zoomin);
                h_Dphi[jdphi][isign] = new TH1D(Form("h_Dphi_%s_sign%d",dphi_regions[jdphi].c_str(),isign+1),";#Delta#phi;1/N_{evt} dN/d#Delta#phi", nDphi_bins,-pms.PI,pms.PI);
                h_pair_y[jdphi][isign] = new TH1D(Form("h_pair_y_%s_sign%d",dphi_regions[jdphi].c_str(),isign+1),";y_{pair};1/N_{evt} dN/dy_{pair}" ,npair_y_bins,-3,3);
                h_pt_asym[jdphi][isign] = new TH1D(Form("h_pt_asym_%s_sign%d",dphi_regions[jdphi].c_str(),isign+1),";A = (pT1 - pT2)/(pT1 + pT2);d#sigma/dA", npt_asym_bins,0,1.);
                h_pair_pt_ptlead_ratio[jdphi][isign] = new TH1D(Form("h_pair_pt_ptlead_ratio_%s_sign%d",dphi_regions[jdphi].c_str(),isign+1),";#frac{p_{T}^{pair}}{p_{T}^{lead}};d#sigma/d#frac{p_{T}^{pair}}{p_{T}^{lead}}", npair_pt_ptlead_ratio_bins,0,2.);
                h_mQQ[jdphi][isign] = new TH1D(Form("h_mQQ_%s_sign%d",dphi_regions[jdphi].c_str(),isign+1),";m_{Q#Bar{Q}}", 200, 0, 200.);
                h_mQQ_Q_ratio[jdphi][isign] = new TH1D(Form("h_mQQ_Q_ratio_%s_sign%d",dphi_regions[jdphi].c_str(),isign+1),";#frac{m_{Q#Bar{Q}}}{Q};d#sigma/d#frac{m_{Q#Bar{Q}}}{Q}", 100,0,5.);
                h_mQQ_mHard_ratio[jdphi][isign] = new TH1D(Form("h_mQQ_mHard_ratio_%s_sign%d",dphi_regions[jdphi].c_str(),isign+1),";#frac{m_{Q#Bar{Q}}}{#sqrt{#hat{s}};d#sigma/d#frac{m_{Q#Bar{Q}}}{#sqrt{#hat{s}}", 50,0,1.);

                h_eta_avg_Dphi[jdphi][isign] = new TH2D(Form("h_eta_avg_Dphi_%s_sign%d",dphi_regions[jdphi].c_str(),isign+1),";#Delta#phi;#bar{#eta}", nDphi_bins,-pms.PI,pms.PI,neta_bins,-2.4,2.4);
                h_Deta_Dphi[jdphi][isign] = new TH2D(Form("h_Deta_Dphi_%s_sign%d",dphi_regions[jdphi].c_str(),isign+1),";#Delta#phi;#Delta#eta", nDphi_bins,-pms.PI,pms.PI,nDeta_bins,-4.8,4.8);
                h_eta1_eta2[jdphi][isign] = new TH2D(Form("h_eta1_eta2_%s_sign%d",dphi_regions[jdphi].c_str(),isign+1),";#eta_{sublead};#eta_{lead}",neta_bins,-2.4,2.4, neta_bins,-2.4,2.4);
                h_eta_avg_Deta[jdphi][isign] = new TH2D(Form("h_eta_avg_Deta_%s_sign%d",dphi_regions[jdphi].c_str(),isign+1),";#Delta#eta;#bar{#eta}",nDeta_bins,-4.8,4.8,neta_bins,-2.4,2.4);
                h_pt1_pt2[jdphi][isign] = new TH2D(Form("h_pt1_pt2_%s_sign%d",dphi_regions[jdphi].c_str(),isign+1),";p_{T}^{sublead} [GeV];p_{T}^{lead} [GeV]",pms.npt_bins,pms.pTBins,pms.npt_bins,pms.pTBins);
                h_ptlead_pair_pt[jdphi][isign] = new TH2D(Form("h_ptlead_pair_pt_%s_sign%d",dphi_regions[jdphi].c_str(),isign+1),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",npair_pT_bins_linear,0,30,npT_lead_bins_linear,0,30);
                h_ptlead_pair_pt_zoomin[jdphi][isign] = new TH2D(Form("h_ptlead_pair_pt_zoomin_%s_sign%d",dphi_regions[jdphi].c_str(),isign+1),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",npair_pT_bins_linear,0,20,npT_lead_bins_linear,0,15);
                h_ptlead_pair_pt_log[jdphi][isign] = new TH2D(Form("h_ptlead_pair_pt_log_%s_sign%d",dphi_regions[jdphi].c_str(),isign+1),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",pms.npairPT_bins,pms.pairPTBins[isign][2],pms.npt_bins,pms.pTBins);
                h_minv_pair_pt[jdphi][isign] = new TH2D(Form("h_minv_pair_pt_%s_sign%d",dphi_regions[jdphi].c_str(),isign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",npair_pT_bins_linear,0,30,nminv_bins_linear,0,30);
                h_minv_pair_pt_zoomin[jdphi][isign] = new TH2D(Form("h_minv_pair_pt_zoomin_%s_sign%d",dphi_regions[jdphi].c_str(),isign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",npair_pT_bins_linear,0,20,nminv_bins_linear,0,4);
                h_minv_pair_pt_log[jdphi][isign] = new TH2D(Form("h_minv_pair_pt_log_%s_sign%d",dphi_regions[jdphi].c_str(),isign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",pms.npairPT_bins,pms.pairPTBins[isign][2],nminv_bins_log,minv_bins_log[isign]);
                
                h_pair_dP_overP[jdphi][isign]->Sumw2();
                h_DR[jdphi][isign]->Sumw2();
                h_DR_zoomin[jdphi][isign]->Sumw2();
                h_DR_jacobian_corrected[jdphi][isign]->Sumw2();
                h_DR_zoomin_jacobian_corrected[jdphi][isign]->Sumw2();
                h_Dphi[jdphi][isign]->Sumw2();
                h_pair_y[jdphi][isign]->Sumw2();
                h_pt_asym[jdphi][isign]->Sumw2();
                h_pair_pt_ptlead_ratio[jdphi][isign]->Sumw2();
                h_mQQ[jdphi][isign]->Sumw2();
                h_mQQ_Q_ratio[jdphi][isign]->Sumw2();
                h_mQQ_mHard_ratio[jdphi][isign]->Sumw2();

                h_eta_avg_Dphi[jdphi][isign]->Sumw2();
                h_Deta_Dphi[jdphi][isign]->Sumw2();
                h_eta1_eta2[jdphi][isign]->Sumw2();
                h_eta_avg_Deta[jdphi][isign]->Sumw2();
                h_pt1_pt2[jdphi][isign]->Sumw2();
                h_ptlead_pair_pt[jdphi][isign]->Sumw2();
                h_ptlead_pair_pt_zoomin[jdphi][isign]->Sumw2();
                h_ptlead_pair_pt_log[jdphi][isign]->Sumw2();
                h_minv_pair_pt[jdphi][isign]->Sumw2();
                h_minv_pair_pt_zoomin[jdphi][isign]->Sumw2();
                h_minv_pair_pt_log[jdphi][isign]->Sumw2();
            }
        }
    }
}


#endif