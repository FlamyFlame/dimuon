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
#include "ParamsSet.h"
#include "MuonPair.h"
#include "vector"
#include "TH1D.h"
#include "TH2D.h"

class MuonPairPlottingPP{
  	//updates histograms - centrality binned
    //if need histograms with all centrality bins added together
    //can add at the histograms level


private:
  	// --------------------- general settings ---------------------------

    // int mode = 1;
    // bool isScram = true;
    // bool isScram = false;

  	ParamsSet pms;
    int dr_bin_start;

  	// --------------------- output file & histograms ---------------------------

  	TFile *outFile = nullptr;

    TH1D* h_pair_dP_overP_dr_binned[ParamsSet::ndRselcs][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH1D* h_Dphi_dr_binned[ParamsSet::ndRselcs][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH1D* h_DR_dr_binned[ParamsSet::ndRselcs][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH1D* h_pair_y_dr_binned[ParamsSet::ndRselcs][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH1D* h_pt_asym_dr_binned[ParamsSet::ndRselcs][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH1D* h_pair_pt_ptlead_ratio_dr_binned[ParamsSet::ndRselcs][ParamsSet::nSigns][ParamsSet::nGapCuts];

    TH2D* h_eta_avg_Dphi_dr_binned[ParamsSet::ndRselcs][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_Deta_Dphi_dr_binned[ParamsSet::ndRselcs][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_eta1_eta2_dr_binned[ParamsSet::ndRselcs][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_eta_avg_Deta_dr_binned[ParamsSet::ndRselcs][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_pt1_pt2_dr_binned[ParamsSet::ndRselcs][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_ptlead_pair_pt_dr_binned[ParamsSet::ndRselcs][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_minv_pair_pt_dr_binned[ParamsSet::ndRselcs][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_minv_pair_pt_zoomin_dr_binned[ParamsSet::ndRselcs][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_minv_pair_pt_log_dr_binned[ParamsSet::ndRselcs][ParamsSet::nSigns][ParamsSet::nGapCuts];

    // TH2D* h_unweighted_Deta_Dphi_dr_binned[ParamsSet::ndRselcs][ParamsSet::nSigns][ParamsSet::nGapCuts];
    
    std::string dphi_regions[2] = {"near", "away"};

    TH1D* h_pair_dP_overP[2][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH1D* h_Dphi[2][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH1D* h_DR[2][ParamsSet::nSigns][ParamsSet::nGapCuts];
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
    TH2D* h_minv_pair_pt_log[2][ParamsSet::nSigns][ParamsSet::nGapCuts];

    static const int nAncestorGroups = 4;

    TH1D* h_DR_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 1];
    TH1D* h_pt_asym_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 1];
    // TH1D* h_psrapidity_ordered_pt_asym_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 1];
    TH1D* h_pair_pt_ptlead_ratio_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 1];
    TH2D* h_ptlead_pair_pt_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 1];
    TH2D* h_Deta_Dphi_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 1];
    TH2D* h_minv_pair_pt_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 1];

    // TH2D* h_unweighted_eta1_eta2_dphicut[ParamsSet::ndphiselcs][ParamsSet::nSigns][ParamsSet::nGapCuts];

    // --------------------- input files & trees & data for setting branches---------------------------

   	TFile *inFile;
    // TTree *inTree_ctrbin[ParamsSet::ndRselcs][ParamsSet::nSigns];
    // TTree *inTree[ParamsSet::nSigns];
    TTree *inTree[ParamsSet::ndRselcs][ParamsSet::nSigns];

    double weight[ParamsSet::ndRselcs][ParamsSet::nSigns];
  	float pair_dPoverP[ParamsSet::ndRselcs][ParamsSet::nSigns];
  	float pt_lead[ParamsSet::ndRselcs][ParamsSet::nSigns];
  	float pair_pt[ParamsSet::ndRselcs][ParamsSet::nSigns];
  	// float pair_eta[ParamsSet::ndRselcs][ParamsSet::nSigns];
  	float pair_y[ParamsSet::ndRselcs][ParamsSet::nSigns];
    float asym[ParamsSet::ndRselcs][ParamsSet::nSigns];
  	float dpt[ParamsSet::ndRselcs][ParamsSet::nSigns];
  	float deta[ParamsSet::ndRselcs][ParamsSet::nSigns];
  	float etaavg[ParamsSet::ndRselcs][ParamsSet::nSigns];
  	float phiavg[ParamsSet::ndRselcs][ParamsSet::nSigns];
  	float dphi[ParamsSet::ndRselcs][ParamsSet::nSigns];
  	float dr[ParamsSet::ndRselcs][ParamsSet::nSigns];
  	float minv[ParamsSet::ndRselcs][ParamsSet::nSigns];
  	float m1pt[ParamsSet::ndRselcs][ParamsSet::nSigns];
  	float m2pt[ParamsSet::ndRselcs][ParamsSet::nSigns];
  	float m1eta[ParamsSet::ndRselcs][ParamsSet::nSigns];
  	float m2eta[ParamsSet::ndRselcs][ParamsSet::nSigns];
  	float m1phi[ParamsSet::ndRselcs][ParamsSet::nSigns];
  	float m2phi[ParamsSet::ndRselcs][ParamsSet::nSigns];
    int m1charge[ParamsSet::ndRselcs][ParamsSet::nSigns];
    int m2charge[ParamsSet::ndRselcs][ParamsSet::nSigns];


    // --------------------- class methods ---------------------------

   	void InitInput();
   	void InitHists();
   	void ProcessData();
    // void InitOutput();
    void WriteOutput();
    bool PassSingleMuonGapCut(float meta, float mpt, int mcharge);
   	void FillHistograms(int ndr, int nsign);
   	void FillPtBinnedHistograms(int ndr, int npt, int nsign);

public:
    int mode = 1;
    bool isScram;
    bool isTight;
    // bool isMCTruthBB;
    // bool isMCTruthCC;
    bool saveDRbinned;
    // assert(isScram & isMCTruthBB & isMCTruthCC == 1 || isScram & isMCTruthBB & isMCTruthCC == 1); // at most one can be true (if all false: real data)
  	MuonPairPlottingPP();
  	~MuonPairPlottingPP(){}
  	void Run();

};

MuonPairPlottingPP::MuonPairPlottingPP(){
    // if (mode != 1 && mode != 3){
    //     std::cout<<"Error:: Mode has to be 1 (no binning) or 3 (binning by pT),  quitting"<<std::endl;
    //     throw std::exception();
    // }
    isScram = false;
    isTight = false;
    // isMCTruthBB = false;
    // isMCTruthCC = false;
    saveDRbinned = true;
    dr_bin_start = 2;
}

void MuonPairPlottingPP::InitInput(){

    if (saveDRbinned) dr_bin_start = 0;

    if (isScram){
        inFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/scrambled_muon_pairs_pp.root","read");
    }else if (isTight){
        inFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/muon_pairs_pp_tight.root","read");
    // }else if(isMCTruthBB){
    //     inFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/athena/runMCV2/muon_pairs_mc_truth_bb.root","read");
    // }else if(isMCTruthCC){
    //     inFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/athena/runMCV2/muon_pairs_mc_truth_cc.root","read");
    }else{
   	    // inFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/muon_pairs_pp_no_resn_cut.root","read");
        // inFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/muon_pairs_pp_no_resn_no_photo.root","read");
        inFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/muon_pairs_pp.root","read");
        // inFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/muon_pairs_pp_1s_only.root","read");
    }

    for (int idr = 0; idr < ParamsSet::ndRselcs; idr++){
        for (int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
            if (isScram){
                inTree[idr][ksign] = (TTree*) inFile->Get(Form("scramb_muon_pair_tree_dr%d_sign%d",idr+1,ksign+1));  
            }else{
        	    inTree[idr][ksign] = (TTree*) inFile->Get(Form("muon_pair_tree_dr%d_sign%d",idr+1,ksign+1));  
            }
            inTree[idr][ksign]->SetBranchAddress("weight"          , &weight[idr][ksign]);
        	inTree[idr][ksign]->SetBranchAddress("pair_dPoverP"    , &pair_dPoverP[idr][ksign]);
        	inTree[idr][ksign]->SetBranchAddress("pt_lead"          , &pt_lead[idr][ksign]);
        	inTree[idr][ksign]->SetBranchAddress("pair_pt"          , &pair_pt[idr][ksign]);
        	// inTree[idr][ksign]->SetBranchAddress("pair_eta"     , &pair_eta[idr][ksign]);
        	inTree[idr][ksign]->SetBranchAddress("pair_y"           , &pair_y[idr][ksign]);
            inTree[idr][ksign]->SetBranchAddress("asym"             , &asym[idr][ksign]);
        	inTree[idr][ksign]->SetBranchAddress("dpt"           , &dpt[idr][ksign]);
        	inTree[idr][ksign]->SetBranchAddress("deta"       , &deta[idr][ksign]);
        	inTree[idr][ksign]->SetBranchAddress("etaavg"      , &etaavg[idr][ksign]);
        	inTree[idr][ksign]->SetBranchAddress("phiavg"            , &phiavg[idr][ksign]);
        	inTree[idr][ksign]->SetBranchAddress("dphi"     , &dphi[idr][ksign]);
        	inTree[idr][ksign]->SetBranchAddress("dr"        , &dr[idr][ksign]);
        	inTree[idr][ksign]->SetBranchAddress("minv"        , &minv[idr][ksign]);
        	inTree[idr][ksign]->SetBranchAddress("m1.pt"           , &m1pt[idr][ksign]);
        	inTree[idr][ksign]->SetBranchAddress("m2.pt"           , &m2pt[idr][ksign]);
        	inTree[idr][ksign]->SetBranchAddress("m1.eta"       , &m1eta[idr][ksign]);
        	inTree[idr][ksign]->SetBranchAddress("m2.eta"       , &m2eta[idr][ksign]);
        	inTree[idr][ksign]->SetBranchAddress("m1.phi"     	, &m1phi[idr][ksign]);
        	inTree[idr][ksign]->SetBranchAddress("m2.phi"     	, &m2phi[idr][ksign]);
            inTree[idr][ksign]->SetBranchAddress("m1.charge"           , &m1charge[idr][ksign]);
            inTree[idr][ksign]->SetBranchAddress("m2.charge"           , &m2charge[idr][ksign]);

            inTree[idr][ksign]->SetBranchStatus("*"                  ,0);//switch off all branches, then enable just the ones that we need
            inTree[idr][ksign]->SetBranchStatus("weight"           ,1);
            inTree[idr][ksign]->SetBranchStatus("pair_dPoverP"           ,1);
            inTree[idr][ksign]->SetBranchStatus("pt_lead"           ,1);
            inTree[idr][ksign]->SetBranchStatus("pair_pt"           ,1);
            inTree[idr][ksign]->SetBranchStatus("pair_y"           ,1);
            inTree[idr][ksign]->SetBranchStatus("asym"           ,1);
            inTree[idr][ksign]->SetBranchStatus("dpt"           ,1);
            inTree[idr][ksign]->SetBranchStatus("deta"           ,1);
            inTree[idr][ksign]->SetBranchStatus("etaavg"           ,1);
            inTree[idr][ksign]->SetBranchStatus("phiavg"           ,1);
            inTree[idr][ksign]->SetBranchStatus("dphi"           ,1);
            inTree[idr][ksign]->SetBranchStatus("dr"           ,1);
            inTree[idr][ksign]->SetBranchStatus("minv"           ,1);
            inTree[idr][ksign]->SetBranchStatus("m1.pt"           ,1);
            inTree[idr][ksign]->SetBranchStatus("m2.pt"           ,1);
            inTree[idr][ksign]->SetBranchStatus("m1.eta"           ,1);
            inTree[idr][ksign]->SetBranchStatus("m2.eta"           ,1);
            inTree[idr][ksign]->SetBranchStatus("m1.phi"           ,1);
            inTree[idr][ksign]->SetBranchStatus("m2.phi"           ,1);
            inTree[idr][ksign]->SetBranchStatus("m1.charge"           ,1);
            inTree[idr][ksign]->SetBranchStatus("m2.charge"           ,1);

        }
    }
}

void MuonPairPlottingPP::InitHists(){
    // int nDR_bins = (isMCTruthBB || isMCTruthCC)? 80 : 200;
    // int nDphi_bins = (isMCTruthBB || isMCTruthCC)? 64 : 128;
    // int neta_bins = (isMCTruthBB || isMCTruthCC)? 50 : 100;
    // int nDeta_bins = (isMCTruthBB || isMCTruthCC)? 100 : 200;
    // int npair_y_bins = (isMCTruthBB || isMCTruthCC)? 45 : 90;
    // int npt_asym_bins = (isMCTruthBB || isMCTruthCC)? 50 : 100;
    // int npair_pt_ptlead_ratio_bins = (isMCTruthBB || isMCTruthCC)? 50 : 100;
    // int nminv_bins_linear = (isMCTruthBB || isMCTruthCC)? 50 : 100;
    // int npair_pT_bins_linear = (isMCTruthBB || isMCTruthCC)? 50 : 100;
    // int npT_lead_bins_linear = (isMCTruthBB || isMCTruthCC)? 50 : 100;
    int nDR_bins = 200;
    int nDphi_bins = 128;
    int neta_bins = 100;
    int nDeta_bins =  200;
    int npair_y_bins = 90;
    int npt_asym_bins = 100;
    int npair_pt_ptlead_ratio_bins = 100;
    int nminv_bins_linear = 100;
    int npair_pT_bins_linear = 100;
    int npT_lead_bins_linear = 100;

    static const int nminv_bins_log = 40;
    float minv_logpow[ParamsSet::nSigns];
    minv_logpow[0] = 0.062;
    minv_logpow[1] = 0.045;
    float minv_max = 60;
    double minv_bins_log[ParamsSet::nSigns][nminv_bins_log+1];

    for (unsigned int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
        for(int iminv = 0; iminv <= nminv_bins_log; iminv++){
            minv_bins_log[ksign][iminv] = minv_max * pow(10.0, (static_cast<float>(iminv - nminv_bins_log)) * minv_logpow[ksign]);
        }
    }

    std::string ancestor_grps[nAncestorGroups + 1] = {"_from_same_b", "_gg", "_qg","_single_g","_qq"};

   	if (mode == 1){
        for (unsigned int ksign = 0; ksign < ParamsSet::nSigns; ksign++){

            for (int kgrp = 0; kgrp < nAncestorGroups + 1; kgrp++){
                h_DR_ancestor_binned[ksign][kgrp] = new TH1D(Form("h_DR_sign%d%s",ksign+1,ancestor_grps[kgrp].c_str()),";#Delta R;1/N_{evt} dN/d#Delta R", nDR_bins,0,pms.deltaR_thrsh[2]);
                h_pt_asym_ancestor_binned[ksign][kgrp] = new TH1D(Form("h_pt_asym_sign%d%s",ksign+1,ancestor_grps[kgrp].c_str()),";A = (pT1 - pT2)/(pT1 + pT2);d#sigma/dA", npt_asym_bins,0,1.);
                // h_psrapidity_ordered_pt_asym_ancestor_binned
                h_pair_pt_ptlead_ratio_ancestor_binned[ksign][kgrp] = new TH1D(Form("h_pair_pt_ptlead_ratio_sign%d%s",ksign+1,ancestor_grps[kgrp].c_str()),";#frac{p_{T}^{pair}}{p_{T}^{lead}};d#sigma/d#frac{p_{T}^{pair}}{p_{T}^{lead}}", npair_pt_ptlead_ratio_bins,0,2.);
                h_ptlead_pair_pt_ancestor_binned[ksign][kgrp] = new TH2D(Form("h_ptlead_pair_pt_sign%d%s",ksign+1,ancestor_grps[kgrp].c_str()),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",npair_pT_bins_linear,0,30,npT_lead_bins_linear,0,30);
                h_Deta_Dphi_ancestor_binned[ksign][kgrp] = new TH2D(Form("h_Deta_Dphi_sign%d%s",ksign+1,ancestor_grps[kgrp].c_str()),";#Delta#phi;#Delta#eta", nDphi_bins,-pms.PI,pms.PI,nDeta_bins,-4.8,4.8);
                h_minv_pair_pt_ancestor_binned[ksign][kgrp] = new TH2D(Form("h_minv_pair_pt_sign%d%s",ksign+1,ancestor_grps[kgrp].c_str()),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",npair_pT_bins_linear,0,30,nminv_bins_linear,0,30);
            }

                

            for (unsigned int lgapcut = 0; lgapcut < ParamsSet::nGapCuts; lgapcut++){

                for (unsigned int jdphi = 0; jdphi < 2; jdphi++){
                    h_pair_dP_overP[jdphi][ksign][lgapcut] = new TH1D(Form("h_pair_dP_overP_%s_sign%d_gapcut%d",dphi_regions[jdphi].c_str(),ksign+1,lgapcut+1),";(#Delta p / p)_{pair};1/N_{evt} dN/d(#Delta p / p)_{pair}" ,pms.deltaP_overP_nbins,0,pms.deltaP_overP_max);
                    h_DR[jdphi][ksign][lgapcut] = new TH1D(Form("h_DR_%s_sign%d_gapcut%d",dphi_regions[jdphi].c_str(),ksign+1,lgapcut+1),";#Delta R;1/N_{evt} dN/d#Delta R", nDR_bins,0,pms.deltaR_thrsh[2]);
                    h_Dphi[jdphi][ksign][lgapcut] = new TH1D(Form("h_Dphi_%s_sign%d_gapcut%d",dphi_regions[jdphi].c_str(),ksign+1,lgapcut+1),";#Delta#phi;1/N_{evt} dN/d#Delta#phi", nDphi_bins,-pms.PI,pms.PI);
                    h_pair_y[jdphi][ksign][lgapcut] = new TH1D(Form("h_pair_y_%s_sign%d_gapcut%d",dphi_regions[jdphi].c_str(),ksign+1,lgapcut+1),";y_{pair};1/N_{evt} dN/dy_{pair}" ,npair_y_bins,-3,3);
                    h_pt_asym[jdphi][ksign][lgapcut] = new TH1D(Form("h_pt_asym_%s_sign%d_gapcut%d",dphi_regions[jdphi].c_str(),ksign+1,lgapcut+1),";A = (pT1 - pT2)/(pT1 + pT2);d#sigma/dA", npt_asym_bins,0,1.);
                    h_pair_pt_ptlead_ratio[jdphi][ksign][lgapcut] = new TH1D(Form("h_pair_pt_ptlead_ratio_%s_sign%d_gapcut%d",dphi_regions[jdphi].c_str(),ksign+1,lgapcut+1),";#frac{p_{T}^{pair}}{p_{T}^{lead}};d#sigma/d#frac{p_{T}^{pair}}{p_{T}^{lead}}", npair_pt_ptlead_ratio_bins,0,2.);

                    h_eta_avg_Dphi[jdphi][ksign][lgapcut] = new TH2D(Form("h_eta_avg_Dphi_%s_sign%d_gapcut%d",dphi_regions[jdphi].c_str(),ksign+1,lgapcut+1),";#Delta#phi;#bar{#eta}", nDphi_bins,-pms.PI,pms.PI,neta_bins,-2.4,2.4);
                    h_Deta_Dphi[jdphi][ksign][lgapcut] = new TH2D(Form("h_Deta_Dphi_%s_sign%d_gapcut%d",dphi_regions[jdphi].c_str(),ksign+1,lgapcut+1),";#Delta#phi;#Delta#eta", nDphi_bins,-pms.PI,pms.PI,nDeta_bins,-4.8,4.8);
                    h_eta1_eta2[jdphi][ksign][lgapcut] = new TH2D(Form("h_eta1_eta2_%s_sign%d_gapcut%d",dphi_regions[jdphi].c_str(),ksign+1,lgapcut+1),";#eta_{sublead};#eta_{lead}",neta_bins,-2.4,2.4, neta_bins,-2.4,2.4);
                    h_eta_avg_Deta[jdphi][ksign][lgapcut] = new TH2D(Form("h_eta_avg_Deta_%s_sign%d_gapcut%d",dphi_regions[jdphi].c_str(),ksign+1,lgapcut+1),";#Delta#eta;#bar{#eta}",nDeta_bins,-4.8,4.8,neta_bins,-2.4,2.4);
                    h_pt1_pt2[jdphi][ksign][lgapcut] = new TH2D(Form("h_pt1_pt2_%s_sign%d_gapcut%d",dphi_regions[jdphi].c_str(),ksign+1,lgapcut+1),";p_{T}^{sublead} [GeV];p_{T}^{lead} [GeV]",pms.npt_bins,pms.pTBins,pms.npt_bins,pms.pTBins);
                    h_ptlead_pair_pt[jdphi][ksign][lgapcut] = new TH2D(Form("h_ptlead_pair_pt_%s_sign%d_gapcut%d",dphi_regions[jdphi].c_str(),ksign+1,lgapcut+1),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",npair_pT_bins_linear,0,30,npT_lead_bins_linear,0,30);
                    h_ptlead_pair_pt_zoomin[jdphi][ksign][lgapcut] = new TH2D(Form("h_ptlead_pair_pt_zoomin_%s_sign%d_gapcut%d",dphi_regions[jdphi].c_str(),ksign+1,lgapcut+1),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",npair_pT_bins_linear,0,20,npT_lead_bins_linear,0,15);
                    h_ptlead_pair_pt_log[jdphi][ksign][lgapcut] = new TH2D(Form("h_ptlead_pair_pt_log_%s_sign%d_gapcut%d",dphi_regions[jdphi].c_str(),ksign+1,lgapcut+1),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",pms.npairPT_bins,pms.pairPTBins[ksign][2],pms.npt_bins,pms.pTBins);
                    h_minv_pair_pt[jdphi][ksign][lgapcut] = new TH2D(Form("h_minv_pair_pt_%s_sign%d_gapcut%d",dphi_regions[jdphi].c_str(),ksign+1,lgapcut+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",npair_pT_bins_linear,0,30,nminv_bins_linear,0,30);
                    h_minv_pair_pt_zoomin[jdphi][ksign][lgapcut] = new TH2D(Form("h_minv_pair_pt_zoomin_%s_sign%d_gapcut%d",dphi_regions[jdphi].c_str(),ksign+1,lgapcut+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",npair_pT_bins_linear,0,20,nminv_bins_linear,5,15);
                    h_minv_pair_pt_log[jdphi][ksign][lgapcut] = new TH2D(Form("h_minv_pair_pt_log_%s_sign%d_gapcut%d",dphi_regions[jdphi].c_str(),ksign+1,lgapcut+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",pms.npairPT_bins,pms.pairPTBins[ksign][2],nminv_bins_log,minv_bins_log[ksign]);
                }

                if (saveDRbinned){
       	            for (unsigned int idr = 0; idr < ParamsSet::ndRselcs; idr++){
                        h_pair_dP_overP_dr_binned[idr][ksign][lgapcut] = new TH1D(Form("h_pair_dP_overP_dr%d_sign%d_gapcut%d",idr+1,ksign+1,lgapcut+1),";(#Delta p / p)_{pair};1/N_{evt} dN/d(#Delta p / p)_{pair}" ,pms.deltaP_overP_nbins,0,pms.deltaP_overP_max);
                        h_DR_dr_binned[idr][ksign][lgapcut] = new TH1D(Form("h_DR_dr%d_sign%d_gapcut%d",idr+1,ksign+1,lgapcut+1),";#Delta R;1/N_{evt} dN/d#Delta R", static_cast<int>(pms.deltaR_nbins[idr]/2.),0,pms.deltaR_thrsh[idr]);
                        h_Dphi_dr_binned[idr][ksign][lgapcut] = new TH1D(Form("h_Dphi_dr%d_sign%d_gapcut%d",idr+1,ksign+1,lgapcut+1),";#Delta#phi;1/N_{evt} dN/d#Delta#phi", nDphi_bins,-pms.PI,pms.PI);
                        h_pair_y_dr_binned[idr][ksign][lgapcut] = new TH1D(Form("h_pair_y_dr%d_sign%d_gapcut%d",idr+1,ksign+1,lgapcut+1),";y_{pair};1/N_{evt} dN/dy_{pair}" ,npair_y_bins,-3,3);
                        h_pt_asym_dr_binned[idr][ksign][lgapcut] = new TH1D(Form("h_pt_asym_dr%d_sign%d_gapcut%d",idr+1,ksign+1,lgapcut+1),";A = (pT1 - pT2)/(pT1 + pT2);d#sigma/dA", npt_asym_bins,0,1.);
                        h_pair_pt_ptlead_ratio_dr_binned[idr][ksign][lgapcut] = new TH1D(Form("h_pair_pt_ptlead_ratio_dr%d_sign%d_gapcut%d",idr+1,ksign+1,lgapcut+1),";#frac{p_{T}^{pair}}{p_{T}^{lead}};d#sigma/d#frac{p_{T}^{pair}}{p_{T}^{lead}}", npair_pt_ptlead_ratio_bins,0,2.);

                        h_eta_avg_Dphi_dr_binned[idr][ksign][lgapcut] = new TH2D(Form("h_eta_avg_Dphi_dr%d_sign%d_gapcut%d",idr+1,ksign+1,lgapcut+1),";#Delta#phi;#bar{#eta}", nDphi_bins,-pms.PI,pms.PI,neta_bins,-2.4,2.4);
                        h_Deta_Dphi_dr_binned[idr][ksign][lgapcut] = new TH2D(Form("h_Deta_Dphi_dr%d_sign%d_gapcut%d",idr+1,ksign+1,lgapcut+1),";#Delta#phi;#Delta#eta", nDphi_bins,-pms.PI,pms.PI,nDeta_bins,-4.8,4.8);
                        h_eta1_eta2_dr_binned[idr][ksign][lgapcut] = new TH2D(Form("h_eta1_eta2_dr%d_sign%d_gapcut%d",idr+1,ksign+1,lgapcut+1),";#eta_{sublead};#eta_{lead}",neta_bins,-2.4,2.4, neta_bins,-2.4,2.4);
                        h_eta_avg_Deta_dr_binned[idr][ksign][lgapcut] = new TH2D(Form("h_eta_avg_Deta_dr%d_sign%d_gapcut%d",idr+1,ksign+1,lgapcut+1),";#Delta#eta;#bar{#eta}",nDeta_bins,-4.8,4.8,neta_bins,-2.4,2.4);
                        h_pt1_pt2_dr_binned[idr][ksign][lgapcut] = new TH2D(Form("h_pt1_pt2_dr%d_sign%d_gapcut%d",idr+1,ksign+1,lgapcut+1),";p_{T}^{sublead} [GeV];p_{T}^{lead} [GeV]",pms.npt_bins,pms.pTBins,pms.npt_bins,pms.pTBins);
                        h_ptlead_pair_pt_dr_binned[idr][ksign][lgapcut] = new TH2D(Form("h_ptlead_pair_pt_dr%d_sign%d_gapcut%d",idr+1,ksign+1,lgapcut+1),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",pms.npairPT_bins,pms.pairPTBins[ksign][idr],pms.npt_bins,pms.pTBins);
                        h_minv_pair_pt_dr_binned[idr][ksign][lgapcut] = new TH2D(Form("h_minv_pair_pt_dr%d_sign%d_gapcut%d",idr+1,ksign+1,lgapcut+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",npair_pT_bins_linear,0,30,nminv_bins_linear,0,30);
                        h_minv_pair_pt_zoomin_dr_binned[idr][ksign][lgapcut] = new TH2D(Form("h_minv_pair_pt_zoomin_dr%d_sign%d_gapcut%d",idr+1,ksign+1,lgapcut+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",npair_pT_bins_linear,0,20,nminv_bins_linear,5,15);
                        h_minv_pair_pt_log_dr_binned[idr][ksign][lgapcut] = new TH2D(Form("h_minv_pair_pt_log_dr%d_sign%d_gapcut%d",idr+1,ksign+1,lgapcut+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",pms.npairPT_bins,pms.pairPTBins[ksign][2],nminv_bins_log,minv_bins_log[ksign]);

                        // h_minv_pair_pt_dr_binned[idr][ksign][lgapcut] = new TH2D(Form("h_minv_pair_pt_dr%d_sign%d_gapcut%d",idr+1,ksign+1,lgapcut+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",pms.npairPT_bins,pms.pairPTBins[ksign][idr],pms.minv_nbins[idr],0,pms.minv_max[idr]);
                        // h_minv_pair_pt_dr_binned[idr][ksign][lgapcut] = new TH2D(Form("h_minv_pair_pt_dr%d_sign%d_gapcut%d",idr+1,ksign+1,lgapcut+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",pms.npairPT_bins,pms.pairPTBins[ksign][idr],nminv_bins,minv_bins[ksign]);
                        
                        // if (isMCTruthBB || isMCTruthCC){
                        //     h_unweighted_Deta_Dphi_dr_binned[idr][ksign][lgapcut] = new TH2D(Form("h_unweighted_Deta_Dphi_dr%d_sign%d_gapcut%d",idr+1,ksign+1,lgapcut+1),";#Delta#phi;#Delta#eta", nDphi_bins,-pms.PI,pms.PI,nDeta_bins,-4.8,4.8);
                        // }
                    }
                }



                // for (unsigned int idphi= 0; idphi < ParamsSet::ndphiselcs; idphi++){
                //     // h_eta1_eta2_dphicut[idphi][ksign][lgapcut] = new TH2D(Form("h_eta1_eta2_dphi%d_sign%d_gapcut%d",idphi+1,ksign+1,lgapcut+1),";#eta_{sublead};#eta_{lead}",100,-2.4,2.4, 100,-2.4,2.4);
                //     h_eta1_eta2_dphicut[idphi][ksign][lgapcut] = new TH2D(Form("h_eta1_eta2_dphi%d_sign%d_gapcut%d",idphi+1,ksign+1,lgapcut+1),";#eta_{sublead};#eta_{lead}",50,-2.4,2.4, 50,-2.4,2.4);
                //     // if (isMCTruthBB || isMCTruthCC) h_unweighted_eta1_eta2_dphicut[idphi][ksign][lgapcut] = new TH2D(Form("h_unweighted_eta1_eta2_dphi%d_sign%d_gapcut%d",idphi+1,ksign+1,lgapcut+1),";#eta_{sublead};#eta_{lead}",50,-2.4,2.4, 50,-2.4,2.4);
                // }
            }
        }
    }
}


#endif