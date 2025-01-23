// Assumes that we exclude resonances

#ifndef __MuonPairPlotting_h__
#define __MuonPairPlotting_h__

#include <TROOT.h>
// #include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <string.h>
#include  <stdlib.h>
#include <assert.h>
#include <fstream>
#include "ParamsSet.h"
#include "MuonPair.h"
#include "vector"
#include "TH1D.h"
#include "TH2D.h"

class MuonPairPlotting{
  	//updates histograms - centrality binned
    //if need histograms with all centrality bins added together
    //can add at the histograms level


private:
  	// --------------------- general settings ---------------------------

    // bool isScram = false;

  	ParamsSet pms;

    static const int nCtrIntvls = 16; // number of size-5% small centrality intervals
    static const int nCtrBins = 6;
    static const int ndphis = 2;

    // 0-5, 5-10, 10-20, 20-30, 30-50, 50-80
    std::vector<std::vector<int>> ctrBins = {{0},{1},{2,3},{4,5},{6,7,8,9},{10,11,12,13,14,15}}; // centrality bins as collections of centrality intervals
    int scaleFactorCtrs[nCtrBins] = {1,1,1,1,1,1}; // no rebinning
    // int scaleFactorCtrs = {1,1,1,1,1,2}; // rebinning for the last centrality bin (50-80)
    // assert (nCtrBins == ctrBins.size()); // number of centrality bins

  	// --------------------- output file & histograms ---------------------------

  	TFile *outFile = nullptr;

    // TH1D* h_FCal_Et;
    // TH1D* h_deltaP_overP;
    TH1D* h_pair_dP_overP[ParamsSet::nEffCorr][ParamsSet::nPtBins][ParamsSet::nCtrBins][ParamsSet::nSigns][ParamsSet::ndphiRegions][ParamsSet::ndetaRegions][ParamsSet::nPhotoProdCuts][ParamsSet::nGapCuts];
    TH1D* h_Dphi[ParamsSet::nEffCorr][ParamsSet::nPtBins][ParamsSet::nCtrBins][ParamsSet::nSigns][ParamsSet::ndphiRegions][ParamsSet::ndetaRegions][ParamsSet::nPhotoProdCuts][ParamsSet::nGapCuts];
    TH1D* h_DR[ParamsSet[ParamsSet::nEffCorr]::nPtBins][ParamsSet::nCtrBins][ParamsSet::nSigns][ParamsSet::ndphiRegions][ParamsSet::ndetaRegions][ParamsSet::nPhotoProdCuts][ParamsSet::nGapCuts];
    TH1D* h_pt_asym[ParamsSet::nEffCorr][ParamsSet::nPtBins][ParamsSet::nCtrBins][ParamsSet::nSigns][ParamsSet::ndphiRegions][ParamsSet::ndetaRegions][ParamsSet::nPhotoProdCuts][ParamsSet::nGapCuts];
    TH1D* h_pair_pt_ptlead_ratio[ParamsSet::nEffCorr][ParamsSet::nPtBins][ParamsSet::nCtrBins][ParamsSet::nSigns][ParamsSet::ndphiRegions][ParamsSet::ndetaRegions][ParamsSet::nPhotoProdCuts][ParamsSet::nGapCuts];
    TH1D* h_pair_y[ParamsSet::nEffCorr][ParamsSet::nPtBins][ParamsSet::nCtrBins][ParamsSet::nSigns][ParamsSet::ndphiRegions][ParamsSet::ndetaRegions][ParamsSet::nPhotoProdCuts][ParamsSet::nGapCuts];
    TH2D* h_eta_avg_Dphi[ParamsSet::nEffCorr][ParamsSet::nPtBins][ParamsSet::nCtrBins][ParamsSet::nSigns][ParamsSet::ndphiRegions][ParamsSet::ndetaRegions][ParamsSet::nPhotoProdCuts][ParamsSet::nGapCuts];
    TH2D* h_Deta_Dphi[ParamsSet::nEffCorr][ParamsSet::nPtBins][ParamsSet::nCtrBins][ParamsSet::nSigns][ParamsSet::ndphiRegions][ParamsSet::ndetaRegions][ParamsSet::nPhotoProdCuts][ParamsSet::nGapCuts];
    TH2D* h_eta1_eta2[ParamsSet::nEffCorr][ParamsSet::nPtBins][ParamsSet::nCtrBins][ParamsSet::nSigns][ParamsSet::ndphiRegions][ParamsSet::ndetaRegions][ParamsSet::nPhotoProdCuts][ParamsSet::nGapCuts];
    TH2D* h_eta_avg_Deta[ParamsSet::nEffCorr][ParamsSet::nPtBins][ParamsSet::nCtrBins][ParamsSet::nSigns][ParamsSet::ndphiRegions][ParamsSet::ndetaRegions][ParamsSet::nPhotoProdCuts][ParamsSet::nGapCuts];
    TH2D* h_eta_avg_pair_eta[ParamsSet::nEffCorr][ParamsSet::nPtBins][ParamsSet::nCtrBins][ParamsSet::nSigns][ParamsSet::ndphiRegions][ParamsSet::ndetaRegions][ParamsSet::nPhotoProdCuts][ParamsSet::nGapCuts];
    TH2D* h_pt1_pt2[ParamsSet::nEffCorr][ParamsSet::nPtBins][ParamsSet::nCtrBins][ParamsSet::nSigns][ParamsSet::ndphiRegions][ParamsSet::ndetaRegions][ParamsSet::nPhotoProdCuts][ParamsSet::nGapCuts];
    TH2D* h_ptlead_pair_pt[ParamsSet::nEffCorr][ParamsSet::nPtBins][ParamsSet::nCtrBins][ParamsSet::nSigns][ParamsSet::ndphiRegions][ParamsSet::ndetaRegions][ParamsSet::nPhotoProdCuts][ParamsSet::nGapCuts];
    TH2D* h_minv_pair_pt[ParamsSet::nEffCorr][ParamsSet::nPtBins][ParamsSet::nCtrBins][ParamsSet::nSigns][ParamsSet::ndphiRegions][ParamsSet::ndetaRegions][ParamsSet::nPhotoProdCuts][ParamsSet::nGapCuts];
    
    // Trigger & reconstruction efficiency corrections
    TH2D* h_trig_reco_eff_pt_avg[ParamsSet::nPtBins][ParamsSet::nCtrBins][ParamsSet::nSigns][ParamsSet::ndphiRegions][ParamsSet::ndetaRegions][ParamsSet::nPhotoProdCuts][ParamsSet::nGapCuts];
    TH2D* h_trig_reco_eff_eta_avg[ParamsSet::nPtBins][ParamsSet::nCtrBins][ParamsSet::nSigns][ParamsSet::ndphiRegions][ParamsSet::ndetaRegions][ParamsSet::nPhotoProdCuts][ParamsSet::nGapCuts];
    TH2D* h_trig_reco_eff_deta[ParamsSet::nPtBins][ParamsSet::nCtrBins][ParamsSet::nSigns][ParamsSet::ndphiRegions][ParamsSet::ndetaRegions][ParamsSet::nPhotoProdCuts][ParamsSet::nGapCuts];
    TH2D* h_trig_reco_eff_dphi[ParamsSet::nPtBins][ParamsSet::nCtrBins][ParamsSet::nSigns][ParamsSet::ndphiRegions][ParamsSet::ndetaRegions][ParamsSet::nPhotoProdCuts][ParamsSet::nGapCuts];
    TH2D* h_trig_reco_eff_centrality[ParamsSet::nPtBins][ParamsSet::nCtrBins][ParamsSet::nSigns][ParamsSet::ndphiRegions][ParamsSet::ndetaRegions][ParamsSet::nPhotoProdCuts][ParamsSet::nGapCuts];

    // TH2D* h_ctrbin_eta1_eta2_dphicut[ParamsSet::ndphiselcs][nCtrBins][ParamsSet::nSigns][ParamsSet::nGapCuts];


    // TH1D* h_ctrbin_Dphi_failgapcut[nCtrBins][ParamsSet::nSigns];
    // TH2D* h_ctrbin_pair_eta_Dphi_failgapcut[nCtrBins][ParamsSet::nSigns];



    // --------------------- input files & trees & data for setting branches---------------------------

   	TFile *inFile;
    // TTree *inTree_ctrbin[nCtrIntvls][ParamsSet::nSigns];
    // TTree *inTree[ParamsSet::nSigns];
    TTree *inTree[nCtrIntvls][ParamsSet::nSigns];

    int avg_centrality; // event centrality for real pairs; average centrality for scrambled pairs (different events)

  	float pair_dPoverP[nCtrIntvls][ParamsSet::nSigns];
  	float pt_lead[nCtrIntvls][ParamsSet::nSigns];
  	float pair_pt[nCtrIntvls][ParamsSet::nSigns];
  	float pair_eta[nCtrIntvls][ParamsSet::nSigns];
  	float pair_y[nCtrIntvls][ParamsSet::nSigns];
    float asym[nCtrIntvls][ParamsSet::nSigns];
  	float dpt[nCtrIntvls][ParamsSet::nSigns];
  	float deta[nCtrIntvls][ParamsSet::nSigns];
  	float etaavg[nCtrIntvls][ParamsSet::nSigns];
  	float phiavg[nCtrIntvls][ParamsSet::nSigns];
  	float dphi[nCtrIntvls][ParamsSet::nSigns];
  	float dr[nCtrIntvls][ParamsSet::nSigns];
  	float minv[nCtrIntvls][ParamsSet::nSigns];
    float m1pt[nCtrIntvls][ParamsSet::nSigns];
  	float m2pt[nCtrIntvls][ParamsSet::nSigns];
  	float m1eta[nCtrIntvls][ParamsSet::nSigns];
  	float m2eta[nCtrIntvls][ParamsSet::nSigns];
  	float m1phi[nCtrIntvls][ParamsSet::nSigns];
  	float m2phi[nCtrIntvls][ParamsSet::nSigns];
    float m1charge[nCtrIntvls][ParamsSet::nSigns];
    float m2charge[nCtrIntvls][ParamsSet::nSigns];

    // --------------------- intermediate variables ---------------------------

    float m1eff;
    float m2eff;
    float pair_eff;

    // --------------------- class methods ---------------------------

   	void InitInput();
   	void InitHists();
   	void ProcessData();
    void WriteOutput();
    bool PassSingleMuonGapCut(float meta, float mpt, int mcharge);
   	void FillHistograms(int nctr_bin, int nctr_intvl, int nsign);

public:
    bool isScram = false;
    bool rebinning = false;
    std::string postfix = "_5_10_20_30_50_80";

  	// MuonPairPlotting(){}
    MuonPairPlotting();
  	~MuonPairPlotting(){}
  	void Run();

};


MuonPairPlotting::MuonPairPlotting(){}


void MuonPairPlotting::InitInput(){

    if (isScram){
        inFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/scrambled_muon_pairs_small_ctr_intvls.root","read");
    }else{
   	    inFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/muon_pairs_small_ctr_intvls.root","read");
    }

    for (int jjctr = 0; jjctr < nCtrIntvls; jjctr++){
        for (int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
            if (isScram){
                inTree[jjctr][ksign] = (TTree*) inFile->Get(Form("scramb_muon_pair_tree_dr%d_ctr%d_sign%d",idr+1,jjctr+1,ksign+1));  
            }else{
        	    inTree[jjctr][ksign] = (TTree*) inFile->Get(Form("muon_pair_tree_dr%d_ctr%d_sign%d",idr+1,jjctr+1,ksign+1));  
            }
        	inTree[jjctr][ksign]->SetBranchAddress("pair_dPoverP"           , &pair_dPoverP[jjctr][ksign]);
        	inTree[jjctr][ksign]->SetBranchAddress("pt_lead"          , &pt_lead[jjctr][ksign]);
        	inTree[jjctr][ksign]->SetBranchAddress("pair_pt"          , &pair_pt[jjctr][ksign]);
        	inTree[jjctr][ksign]->SetBranchAddress("pair_eta"     , &pair_eta[jjctr][ksign]);
        	inTree[jjctr][ksign]->SetBranchAddress("pair_y"           , &pair_y[jjctr][ksign]);
            inTree[jjctr][ksign]->SetBranchAddress("asym"           , &asym[jjctr][ksign]);
        	inTree[jjctr][ksign]->SetBranchAddress("dpt"           , &dpt[jjctr][ksign]);
        	inTree[jjctr][ksign]->SetBranchAddress("deta"       , &deta[jjctr][ksign]);
        	inTree[jjctr][ksign]->SetBranchAddress("etaavg"      , &etaavg[jjctr][ksign]);
        	inTree[jjctr][ksign]->SetBranchAddress("phiavg"            , &phiavg[jjctr][ksign]);
        	inTree[jjctr][ksign]->SetBranchAddress("dphi"     , &dphi[jjctr][ksign]);
        	inTree[jjctr][ksign]->SetBranchAddress("dr"        , &dr[jjctr][ksign]);
        	inTree[jjctr][ksign]->SetBranchAddress("minv"        , &minv[jjctr][ksign]);
        	inTree[jjctr][ksign]->SetBranchAddress("m1.pt"           , &m1pt[jjctr][ksign]);
        	inTree[jjctr][ksign]->SetBranchAddress("m2.pt"           , &m2pt[jjctr][ksign]);
        	inTree[jjctr][ksign]->SetBranchAddress("m1.eta"       , &m1eta[jjctr][ksign]);
        	inTree[jjctr][ksign]->SetBranchAddress("m2.eta"       , &m2eta[jjctr][ksign]);
        	inTree[jjctr][ksign]->SetBranchAddress("m1.phi"     	, &m1phi[jjctr][ksign]);
        	inTree[jjctr][ksign]->SetBranchAddress("m2.phi"     	, &m2phi[jjctr][ksign]);
            inTree[jjctr][ksign]->SetBranchAddress("m1.charge"        , &m1charge[jjctr][ksign]);
            inTree[jjctr][ksign]->SetBranchAddress("m2.charge"        , &m2charge[jjctr][ksign]);
            inTree[jjctr][ksign]->SetBranchAddress("avg_centrality"        , &avg_centrality);

            inTree[jjctr][ksign]->SetBranchStatus("*"                  ,0);//switch off all branches, then enable just the ones that we need
            inTree[jjctr][ksign]->SetBranchStatus("pair_dPoverP"           ,1);
            inTree[jjctr][ksign]->SetBranchStatus("pt_lead"           ,1);
            inTree[jjctr][ksign]->SetBranchStatus("pair_pt"           ,1);
            inTree[jjctr][ksign]->SetBranchStatus("pair_eta"           ,1);
            inTree[jjctr][ksign]->SetBranchStatus("pair_y"           ,1);
            inTree[jjctr][ksign]->SetBranchStatus("asym"           ,1);
            inTree[jjctr][ksign]->SetBranchStatus("dpt"           ,1);
            inTree[jjctr][ksign]->SetBranchStatus("deta"           ,1);
            inTree[jjctr][ksign]->SetBranchStatus("etaavg"           ,1);
            inTree[jjctr][ksign]->SetBranchStatus("phiavg"           ,1);
            inTree[jjctr][ksign]->SetBranchStatus("dphi"           ,1);
            inTree[jjctr][ksign]->SetBranchStatus("dr"           ,1);
            inTree[jjctr][ksign]->SetBranchStatus("minv"           ,1);
            inTree[jjctr][ksign]->SetBranchStatus("m1.pt"           ,1);
            inTree[jjctr][ksign]->SetBranchStatus("m2.pt"           ,1);
            inTree[jjctr][ksign]->SetBranchStatus("m1.eta"           ,1);
            inTree[jjctr][ksign]->SetBranchStatus("m2.eta"           ,1);
            inTree[jjctr][ksign]->SetBranchStatus("m1.phi"           ,1);
            inTree[jjctr][ksign]->SetBranchStatus("m2.phi"           ,1);
            inTree[jjctr][ksign]->SetBranchStatus("m1.charge"           ,1);
            inTree[jjctr][ksign]->SetBranchStatus("m2.charge"           ,1);
            inTree[jjctr][ksign]->SetBranchStatus("avg_centrality"           ,1);
       	}
    }
}

void MuonPairPlotting::InitHists(){

    std::string dphis[ndphis] = {"_near", "_away"};

    // int nDR_bins = 200;
    // int nDphi_bins = 128;
    // int neta_bins = 100;
    // int nDeta_bins =  200;
    // int npair_y_bins = 90;
    // int npt_asym_bins = 100;
    // int npair_pt_ptlead_ratio_bins = 100;
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

    for (unsigned int icorr = 0; icorr < ParamsSet::nPtBins; icorr++){
        for (unsigned int ipt = 0; ipt < ParamsSet::nPtBins; ipt++){
            for (unsigned int ictr = 0; ictr < ParamsSet::nCtrBins; ictr++){
                for (unsigned int isign = 0; isign < ParamsSet::nSigns; isign++){
                    for (unsigned int idphi = 0; idphi < ParamsSet::ndphiRegions; idphi++){
                        for (unsigned int ideta = 0; ideta < ParamsSet::ndetaRegions; ideta++){
                            for (unsigned int igap = 0; igap < ParamsSet::nGapCuts; igap++){

                                h_pair_dP_overP[ipt][ictr][isign][idphi][ideta][igap] = new TH1D(Form("h_pair_dP_overP_%s_%s_%s_%s_%s_%s",pms.pt_labels[ipt], pms.ctr_labels[ictr], pms.signs[isign], pms.dphi_regions[idphi], pms.deta_regions[ideta], pms.gapcut_labels[igap]),";(#Delta p / p)_{pair};1/N_{evt} dN/d" ,pms.deltaP_overP_nbins,0,pms.deltaP_overP_max);
                                h_Dphi[ipt][ictr][isign][idphi][ideta][igap] = new TH1D(Form("h_Dphi_%s_%s_%s_%s_%s_%s",pms.pt_labels[ipt], pms.ctr_labels[ictr], pms.signs[isign], pms.dphi_regions[idphi], pms.deta_regions[ideta], pms.gapcut_labels[igap]),";#Delta#phi;1/N_{evt} dN/d#Delta#phi" ,128,-pms.PI ,pms.PI );
                                h_DR[ipt][ictr][isign][idphi][ideta][igap] = new TH1D(Form("h_DR_%s_%s_%s_%s_%s_%s",pms.pt_labels[ipt], pms.ctr_labels[ictr], pms.signs[isign], pms.dphi_regions[idphi], pms.deta_regions[ideta], pms.gapcut_labels[igap]),";#Delta R;1/N_{evt} dN/d#Delta R", pms.deltaR_nbins,0,pms.deltaR_thrsh);
                                h_pt_asym[ipt][ictr][isign][idphi][ideta][igap] = new TH1D(Form("h_pt_asym_%s_%s_%s_%s_%s_%s",pms.pt_labels[ipt], pms.ctr_labels[ictr], pms.signs[isign], pms.dphi_regions[idphi], pms.deta_regions[ideta], pms.gapcut_labels[igap]),";A = (pT1 - pT2)/(pT1 + pT2);d#sigma/dA", 100,0,1.);
                                h_pair_pt_ptlead_ratio[ipt][ictr][isign][idphi][ideta][igap] = new TH1D(Form("h_pair_pt_ptlead_ratio_%s_%s_%s_%s_%s_%s",pms.pt_labels[ipt], pms.ctr_labels[ictr], pms.signs[isign], pms.dphi_regions[idphi], pms.deta_regions[ideta], pms.gapcut_labels[igap]),";#frac{p_{T}^{pair}}{p_{T}^{lead}};d#sigma/d#frac{p_{T}^{pair}}{p_{T}^{lead}}", 100,0,2.);
                                // h_Minv[ipt][ictr][isign][idphi][ideta][igap] = new TH1D(Form("h_Minv_%s_%s_%s_%s_%s_%s",pms.pt_labels[ipt], pms.ctr_labels[ictr], pms.signs[isign], pms.dphi_regions[idphi], pms.deta_regions[ideta], pms.gapcut_labels[igap]),";m_{#mu#mu} [GeV];1/N_{evt} dN/dm_{#mu#mu}", pms.minv_nbins,0   ,pms.minv_max);
                                // h_pt_lead[ipt][ictr][isign][idphi][ideta][igap] = new TH1D(Form("h_pt_lead_%s_%s_%s_%s_%s_%s",pms.pt_labels[ipt], pms.ctr_labels[ictr], pms.signs[isign], pms.dphi_regions[idphi], pms.deta_regions[ideta], pms.gapcut_labels[igap]),";p_{T}^{lead} [GeV];1/N_{evt} dN/dp_{T}^{lead}" ,pms.npt_bins,pms.pTBins );
                                // h_eta_avg[ipt][ictr][isign][idphi][ideta][igap] = new TH1D(Form("h_eta_avg_%s_%s_%s_%s_%s_%s",pms.pt_labels[ipt], pms.ctr_labels[ictr], pms.signs[isign], pms.dphi_regions[idphi], pms.deta_regions[ideta], pms.gapcut_labels[igap]),";#bar{#eta};1/N_{evt} dN/d#bar{#eta}",100,-2.4,2.4);
                                // h_pair_pt[ipt][ictr][isign][idphi][ideta][igap] = new TH1D(Form("h_pair_pt_%s_%s_%s_%s_%s_%s",pms.pt_labels[ipt], pms.ctr_labels[ictr], pms.signs[isign], pms.dphi_regions[idphi], pms.deta_regions[ideta], pms.gapcut_labels[igap]),";p_{T}^{pair} [GeV];1/N_{evt} dN/dp_{T}^{pair}" ,pms.npairPT_bins,pms.pairPTBins[ksign] );
                                // h_pair_eta[ipt][ictr][isign][idphi][ideta][igap] = new TH1D(Form("h_pair_eta_%s_%s_%s_%s_%s_%s",pms.pt_labels[ipt], pms.ctr_labels[ictr], pms.signs[isign], pms.dphi_regions[idphi], pms.deta_regions[ideta], pms.gapcut_labels[igap]),";#eta_{pair};1/N_{evt} dN/d#eta_{pair}" ,100,-2.4,2.4);
                                h_pair_y[ipt][ictr][isign][idphi][ideta][igap] = new TH1D(Form("h_pair_y_%s_%s_%s_%s_%s_%s",pms.pt_labels[ipt], pms.ctr_labels[ictr], pms.signs[isign], pms.dphi_regions[idphi], pms.deta_regions[ideta], pms.gapcut_labels[igap]),";y_{pair};1/N_{evt} dN/dy_{pair}" ,90,-3,3);
                                h_eta_avg_Dphi[ipt][ictr][isign][idphi][ideta][igap] = new TH2D(Form("h_eta_avg_Dphi_%s_%s_%s_%s_%s_%s",pms.pt_labels[ipt], pms.ctr_labels[ictr], pms.signs[isign], pms.dphi_regions[idphi], pms.deta_regions[ideta], pms.gapcut_labels[igap]),";#Delta#phi;#bar{#eta}",128,-pms.PI,pms.PI,100,-2.4,2.4);
                                h_Deta_Dphi[ipt][ictr][isign][idphi][ideta][igap] = new TH2D(Form("h_Deta_Dphi_%s_%s_%s_%s_%s_%s",pms.pt_labels[ipt], pms.ctr_labels[ictr], pms.signs[isign], pms.dphi_regions[idphi], pms.deta_regions[ideta], pms.gapcut_labels[igap]),";#Delta#phi;#Delta#eta",128,-pms.PI,pms.PI,200,-4.8,4.8);
                                h_eta1_eta2[ipt][ictr][isign][idphi][ideta][igap] = new TH2D(Form("h_eta1_eta2_%s_%s_%s_%s_%s_%s",pms.pt_labels[ipt], pms.ctr_labels[ictr], pms.signs[isign], pms.dphi_regions[idphi], pms.deta_regions[ideta], pms.gapcut_labels[igap]),";#eta_{sublead};#eta_{lead}",100,-2.4,2.4, 100,-2.4,2.4);
                                h_eta_avg_Deta[ipt][ictr][isign][idphi][ideta][igap] = new TH2D(Form("h_eta_avg_Deta_%s_%s_%s_%s_%s_%s",pms.pt_labels[ipt], pms.ctr_labels[ictr], pms.signs[isign], pms.dphi_regions[idphi], pms.deta_regions[ideta], pms.gapcut_labels[igap]),";#Delta#eta;#bar{#eta}",200,-4.8,4.8, 100,-2.4,2.4);
                                h_eta_avg_pair_eta[ipt][ictr][isign][idphi][ideta][igap] = new TH2D(Form("h_eta_avg_pair_eta_%s_%s_%s_%s_%s_%s",pms.pt_labels[ipt], pms.ctr_labels[ictr], pms.signs[isign], pms.dphi_regions[idphi], pms.deta_regions[ideta], pms.gapcut_labels[igap]),";#eta_{pair};#bar{#eta}",100,-2.4,2.4, 100,-2.4,2.4);
                                h_pt1_pt2[ipt][ictr][isign][idphi][ideta][igap] = new TH2D(Form("h_pt1_pt2_%s_%s_%s_%s_%s_%s",pms.pt_labels[ipt], pms.ctr_labels[ictr], pms.signs[isign], pms.dphi_regions[idphi], pms.deta_regions[ideta], pms.gapcut_labels[igap]),";p_{T}^{sublead} [GeV];p_{T}^{lead} [GeV]",pms.npt_bins,pms.pTBins,pms.npt_bins,pms.pTBins);
                                h_ptlead_pair_pt[ipt][ictr][isign][idphi][ideta][igap] = new TH2D(Form("h_ptlead_pair_pt_%s_%s_%s_%s_%s_%s",pms.pt_labels[ipt], pms.ctr_labels[ictr], pms.signs[isign], pms.dphi_regions[idphi], pms.deta_regions[ideta], pms.gapcut_labels[igap]),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",pms.npairPT_bins,pms.pairPTBins[ksign],pms.npt_bins,pms.pTBins);
                                h_minv_pair_pt[ipt][ictr][isign][idphi][ideta][igap] = new TH2D(Form("h_minv_pair_pt_%s_%s_%s_%s_%s_%s",pms.pt_labels[ipt], pms.ctr_labels[ictr], pms.signs[isign], pms.dphi_regions[idphi], pms.deta_regions[ideta], pms.gapcut_labels[igap]),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",pms.npairPT_bins,pms.pairPTBins[ksign],pms.minv_nbins,0,pms.minv_max);

                                // int nminvbins = pms.minv_bins[isign].size(); // size of a vector
                                // double minv_arr[pms.minv_nbins+1];
                                // std::copy(pms.minv_bins[ksign].begin(),pms.minv_bins[ksign].end(),minv_arr);
                                // h_minv_pair_pt[ipt][ictr][isign][idphi][ideta][igap] = new TH2D(Form("h_minv_pair_pt_dr%u__sign%u",idr+1,ksign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",pms.npairPT_bins,pms.pairPTBins[ksign],pms.minv_nbins,minv_arr);
        	  	            }
                        }
                    }
                }
            }
        }
    }
}




#endif