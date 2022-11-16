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

    // int mode = 2;
    // bool isScram = false;

  	ParamsSet pms;
    static const int nCtrIntvls = 16; // number of size-5% small centrality intervals
    static const int nCtrBins = 6;
    // 0-5, 5-10, 10-20, 20-30, 30-50, 50-80
    std::vector<std::vector<int>> ctrBins = {{0},{1},{2,3},{4,5},{6,7,8,9},{10,11,12,13,14,15}}; // centrality bins as collections of centrality intervals
    int scaleFactorCtrs[nCtrBins] = {1,1,1,1,1,1}; // no rebinning
    // int scaleFactorCtrs = {1,1,1,1,1,2}; // rebinning for the last centrality bin (50-80)
    // assert (nCtrBins == ctrBins.size()); // number of centrality bins

  	// --------------------- output file & histograms ---------------------------

  	TFile *outFile = nullptr;

    TH1D* h_FCal_Et[ParamsSet::ndRselcs];
    // TH1D* h_deltaP_overP[ParamsSet::ndRselcs];
    TH1D* h_pair_dP_overP[ParamsSet::ndRselcs][ParamsSet::nSigns];
    TH1D* h_Dphi[ParamsSet::ndRselcs][ParamsSet::nSigns];
    TH1D* h_DR[ParamsSet::ndRselcs][ParamsSet::nSigns];
    // TH1D* h_Minv[ParamsSet::ndRselcs][ParamsSet::nSigns];
    // TH1D* h_pt_lead[ParamsSet::ndRselcs][ParamsSet::nSigns];
    // TH1D* h_eta_avg[ParamsSet::ndRselcs][ParamsSet::nSigns];
    // TH1D* h_pair_pt[ParamsSet::ndRselcs][ParamsSet::nSigns];
    // TH1D* h_pair_eta[ParamsSet::ndRselcs][ParamsSet::nSigns];
    TH1D* h_pair_y[ParamsSet::ndRselcs][ParamsSet::nSigns];
    TH2D* h_eta_avg_Dphi[ParamsSet::ndRselcs][ParamsSet::nSigns];
    TH2D* h_Deta_Dphi[ParamsSet::ndRselcs][ParamsSet::nSigns];
    TH2D* h_eta1_eta2[ParamsSet::ndRselcs][ParamsSet::nSigns];
    TH2D* h_eta_avg_Deta[ParamsSet::ndRselcs][ParamsSet::nSigns];
    TH2D* h_eta_avg_pair_eta[ParamsSet::ndRselcs][ParamsSet::nSigns];
    TH2D* h_pt1_pt2[ParamsSet::ndRselcs][ParamsSet::nSigns];
    TH2D* h_ptlead_pair_pt[ParamsSet::ndRselcs][ParamsSet::nSigns];
    TH2D* h_minv_pair_pt[ParamsSet::ndRselcs][ParamsSet::nSigns];

    TH1D* h_ctrbin_pair_dP_overP[ParamsSet::ndRselcs][nCtrBins][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH1D* h_ctrbin_Dphi[ParamsSet::ndRselcs][nCtrBins][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH1D* h_ctrbin_DR[ParamsSet::ndRselcs][nCtrBins][ParamsSet::nSigns][ParamsSet::nGapCuts];
    // TH1D* h_ctrbin_Minv[ParamsSet::ndRselcs][nCtrBins][ParamsSet::nSigns][ParamsSet::nGapCuts];
    // TH1D* h_ctrbin_pt_lead[ParamsSet::ndRselcs][nCtrBins][ParamsSet::nSigns][ParamsSet::nGapCuts];
    // TH1D* h_ctrbin_eta_avg[ParamsSet::ndRselcs][nCtrBins][ParamsSet::nSigns][ParamsSet::nGapCuts];
    // TH1D* h_ctrbin_pair_pt[ParamsSet::ndRselcs][nCtrBins][ParamsSet::nSigns][ParamsSet::nGapCuts];
    // TH1D* h_ctrbin_pair_eta[ParamsSet::ndRselcs][nCtrBins][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH1D* h_ctrbin_pair_y[ParamsSet::ndRselcs][nCtrBins][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_ctrbin_eta_avg_Dphi[ParamsSet::ndRselcs][nCtrBins][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_ctrbin_Deta_Dphi[ParamsSet::ndRselcs][nCtrBins][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_ctrbin_eta1_eta2[ParamsSet::ndRselcs][nCtrBins][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_ctrbin_eta_avg_Deta[ParamsSet::ndRselcs][nCtrBins][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_ctrbin_eta_avg_pair_eta[ParamsSet::ndRselcs][nCtrBins][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_ctrbin_pt1_pt2[ParamsSet::ndRselcs][nCtrBins][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_ctrbin_ptlead_pair_pt[ParamsSet::ndRselcs][nCtrBins][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_ctrbin_minv_pair_pt[ParamsSet::ndRselcs][nCtrBins][ParamsSet::nSigns][ParamsSet::nGapCuts];

    TH2D* h_ctrbin_eta1_eta2_dphicut[ParamsSet::ndphiselcs][nCtrBins][ParamsSet::nSigns][ParamsSet::nGapCuts];


    // TH1D* h_ctrbin_Dphi_failgapcut[ParamsSet::ndRselcs][nCtrBins][ParamsSet::nSigns];
    // TH2D* h_ctrbin_pair_eta_Dphi_failgapcut[ParamsSet::ndRselcs][nCtrBins][ParamsSet::nSigns];



    // --------------------- input files & trees & data for setting branches---------------------------

   	TFile *inFile;
    // TTree *inTree_ctrbin[ParamsSet::ndRselcs][nCtrIntvls][ParamsSet::nSigns];
    // TTree *inTree[ParamsSet::nSigns];
    TTree *inTree[ParamsSet::ndRselcs][nCtrIntvls][ParamsSet::nSigns];

  	float pair_dPoverP[ParamsSet::ndRselcs][nCtrIntvls][ParamsSet::nSigns];
  	float pt_lead[ParamsSet::ndRselcs][nCtrIntvls][ParamsSet::nSigns];
  	float pair_pt[ParamsSet::ndRselcs][nCtrIntvls][ParamsSet::nSigns];
  	float pair_eta[ParamsSet::ndRselcs][nCtrIntvls][ParamsSet::nSigns];
  	float pair_y[ParamsSet::ndRselcs][nCtrIntvls][ParamsSet::nSigns];
  	float dpt[ParamsSet::ndRselcs][nCtrIntvls][ParamsSet::nSigns];
  	float deta[ParamsSet::ndRselcs][nCtrIntvls][ParamsSet::nSigns];
  	float etaavg[ParamsSet::ndRselcs][nCtrIntvls][ParamsSet::nSigns];
  	float phiavg[ParamsSet::ndRselcs][nCtrIntvls][ParamsSet::nSigns];
  	float dphi[ParamsSet::ndRselcs][nCtrIntvls][ParamsSet::nSigns];
  	float dr[ParamsSet::ndRselcs][nCtrIntvls][ParamsSet::nSigns];
  	float minv[ParamsSet::ndRselcs][nCtrIntvls][ParamsSet::nSigns];
  	float m1pt[ParamsSet::ndRselcs][nCtrIntvls][ParamsSet::nSigns];
  	float m2pt[ParamsSet::ndRselcs][nCtrIntvls][ParamsSet::nSigns];
  	float m1eta[ParamsSet::ndRselcs][nCtrIntvls][ParamsSet::nSigns];
  	float m2eta[ParamsSet::ndRselcs][nCtrIntvls][ParamsSet::nSigns];
  	float m1phi[ParamsSet::ndRselcs][nCtrIntvls][ParamsSet::nSigns];
  	float m2phi[ParamsSet::ndRselcs][nCtrIntvls][ParamsSet::nSigns];
    float m1charge[ParamsSet::ndRselcs][nCtrIntvls][ParamsSet::nSigns];
    float m2charge[ParamsSet::ndRselcs][nCtrIntvls][ParamsSet::nSigns];

    // --------------------- class methods ---------------------------

   	void InitInput();
   	void InitHists();
   	void ProcessData();
    void WriteOutput();
    bool PassSingleMuonGapCut(float meta, float mpt);
   	void FillUnbinnedHistograms(int ndr, int nctr_intvl, int nsign);
   	void FillCtrBinnedHistograms(int ndr, int nctr_bin, int nctr_intvl, int nsign);
   	void FillPtBinnedHistograms(int ndr, int nctr_intvl, int nsign);

public:
    int mode = 2;
    bool isScram = false;
    bool rebinning = false;
    std::string postfix = "_5_10_20_30_50_80";

  	// MuonPairPlotting(){}
    MuonPairPlotting();
  	~MuonPairPlotting(){}
  	void Run();

};


MuonPairPlotting::MuonPairPlotting(){
    if (isScram){
        std::cout << "Cannot do scrambled right now. Haven't adjusted the codes for the new centrality bins." << std::endl;
        throw std::exception();
    }
}


void MuonPairPlotting::InitInput(){

    if (isScram){
        inFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/scrambled_muon_pairs_small_ctr_intvls.root","read");
    }else{
   	    inFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/muon_pairs_small_ctr_intvls.root","read");
    }

    for (int idr = 0; idr < ParamsSet::ndRselcs; idr++){
       	for (int jjctr = 0; jjctr < nCtrIntvls; jjctr++){
            for (int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
                if (isScram){
                    inTree[idr][jjctr][ksign] = (TTree*) inFile->Get(Form("scramb_muon_pair_tree_dr%d_ctr%d_sign%d",idr+1,jjctr+1,ksign+1));  
                }else{
            	    inTree[idr][jjctr][ksign] = (TTree*) inFile->Get(Form("muon_pair_tree_dr%d_ctr%d_sign%d",idr+1,jjctr+1,ksign+1));  
                }
            	inTree[idr][jjctr][ksign]->SetBranchAddress("pair_dPoverP"           , &pair_dPoverP[idr][jjctr][ksign]);
            	inTree[idr][jjctr][ksign]->SetBranchAddress("pt_lead"          , &pt_lead[idr][jjctr][ksign]);
            	inTree[idr][jjctr][ksign]->SetBranchAddress("pair_pt"          , &pair_pt[idr][jjctr][ksign]);
            	inTree[idr][jjctr][ksign]->SetBranchAddress("pair_eta"     , &pair_eta[idr][jjctr][ksign]);
            	inTree[idr][jjctr][ksign]->SetBranchAddress("pair_y"           , &pair_y[idr][jjctr][ksign]);
            	inTree[idr][jjctr][ksign]->SetBranchAddress("dpt"           , &dpt[idr][jjctr][ksign]);
            	inTree[idr][jjctr][ksign]->SetBranchAddress("deta"       , &deta[idr][jjctr][ksign]);
            	inTree[idr][jjctr][ksign]->SetBranchAddress("etaavg"      , &etaavg[idr][jjctr][ksign]);
            	inTree[idr][jjctr][ksign]->SetBranchAddress("phiavg"            , &phiavg[idr][jjctr][ksign]);
            	inTree[idr][jjctr][ksign]->SetBranchAddress("dphi"     , &dphi[idr][jjctr][ksign]);
            	inTree[idr][jjctr][ksign]->SetBranchAddress("dr"        , &dr[idr][jjctr][ksign]);
            	inTree[idr][jjctr][ksign]->SetBranchAddress("minv"        , &minv[idr][jjctr][ksign]);
            	inTree[idr][jjctr][ksign]->SetBranchAddress("m1.pt"           , &m1pt[idr][jjctr][ksign]);
            	inTree[idr][jjctr][ksign]->SetBranchAddress("m2.pt"           , &m2pt[idr][jjctr][ksign]);
            	inTree[idr][jjctr][ksign]->SetBranchAddress("m1.eta"       , &m1eta[idr][jjctr][ksign]);
            	inTree[idr][jjctr][ksign]->SetBranchAddress("m2.eta"       , &m2eta[idr][jjctr][ksign]);
            	inTree[idr][jjctr][ksign]->SetBranchAddress("m1.phi"     	, &m1phi[idr][jjctr][ksign]);
            	inTree[idr][jjctr][ksign]->SetBranchAddress("m2.phi"     	, &m2phi[idr][jjctr][ksign]);
                inTree[idr][jjctr][ksign]->SetBranchAddress("m1.charge"        , &m1charge[idr][jjctr][ksign]);
                inTree[idr][jjctr][ksign]->SetBranchAddress("m2.charge"        , &m2charge[idr][jjctr][ksign]);
           	}
        }
    }
}

void MuonPairPlotting::InitHists(){

   	if (mode == 1){
        
    	for (unsigned int idr = 0; idr < ParamsSet::ndRselcs; idr++){
    	    h_FCal_Et[idr]=new TH1D(Form("h_FCal_Et_dr%u",idr+1),";FCal #it{E}_{T} [TeV];1/N_{evt} dN/dFCal #it{E}_{T}",550,-0.5,5.0);
    	    // h_deltaP_overP[idr] =new TH1D(Form("h_deltaP_overP_dr%u",idr+1),";(p_{ID}-p_{MS}-#Delta p_{calo}) / p_{ID};1/N_{evt} dN/d" ,deltaP_overP_nbins*2,-deltaP_overP_thrsh,deltaP_overP_thrsh);
    	
    	    for (unsigned int ksign = 0; ksign < ParamsSet::nSigns; ksign++){

            /////////// --------------------- to be updated (copy) -----------------------------
          	
                h_pair_dP_overP[idr][ksign] = new TH1D(Form("h_pair_dP_overP_dr%u_sign%u",idr+1,ksign+1),";(#Delta p / p)_{pair};1/N_{evt} dN/d" ,pms.deltaP_overP_nbins,0,pms.deltaP_overP_max);
                h_Dphi[idr][ksign] = new TH1D(Form("h_Dphi_dr%u_sign%u",idr+1,ksign+1),";#Delta#phi;1/N_{evt} dN/d#Delta#phi" ,128,-pms.PI ,pms.PI );
                h_DR[idr][ksign] = new TH1D(Form("h_DR_dr%u_sign%u",idr+1,ksign+1),";#Delta R;1/N_{evt} dN/d#Delta R", pms.deltaR_nbins[idr],0,pms.deltaR_thrsh[idr]);
                // h_Minv[idr][ksign] = new TH1D(Form("h_Minv_dr%u_sign%u",idr+1,ksign+1),";m_{#mu#mu} [GeV];1/N_{evt} dN/dm_{#mu#mu}", pms.minv_nbins[idr],0   ,pms.minv_max[idr]);
                // h_pt_lead[idr][ksign] = new TH1D(Form("h_pt_lead_dr%u_sign%u",idr+1,ksign+1),";p_{T}^{lead} [GeV];1/N_{evt} dN/dp_{T}^{lead}" ,pms.npt_bins,pms.pTBins );
                // h_eta_avg[idr][ksign] = new TH1D(Form("h_eta_avg_dr%u_sign%u",idr+1,ksign+1),";#bar{#eta};1/N_{evt} dN/d#bar{#eta}",100,-2.4,2.4);
                // h_pair_pt[idr][ksign] = new TH1D(Form("h_pair_pt_dr%u_sign%u",idr+1,ksign+1),";p_{T}^{pair} [GeV];1/N_{evt} dN/dp_{T}^{pair}" ,pms.npairPT_bins,pms.pairPTBins[ksign][idr] );
                // h_pair_eta[idr][ksign] = new TH1D(Form("h_pair_eta_dr%u_sign%u",idr+1,ksign+1),";#eta_{pair};1/N_{evt} dN/d#eta_{pair}" ,100,-2.4,2.4);
                h_pair_y[idr][ksign] = new TH1D(Form("h_pair_y_dr%u_sign%u",idr+1,ksign+1),";y_{pair};1/N_{evt} dN/dy_{pair}" ,90,-3,3);
                h_eta_avg_Dphi[idr][ksign] = new TH2D(Form("h_eta_avg_Dphi_dr%u_sign%u",idr+1,ksign+1),";#Delta#phi;#bar{#eta}",128,-pms.PI,pms.PI,100,-2.4,2.4);
                h_Deta_Dphi[idr][ksign] = new TH2D(Form("h_Deta_Dphi_dr%u_sign%u",idr+1,ksign+1),";#Delta#phi;#Delta#eta",128,-pms.PI,pms.PI,200,-4.8,4.8);
                h_eta1_eta2[idr][ksign] = new TH2D(Form("h_eta1_eta2_dr%u_sign%u",idr+1,ksign+1),";#eta_{sublead};#eta_{lead}",100,-2.4,2.4, 100,-2.4,2.4);
                h_eta_avg_Deta[idr][ksign] = new TH2D(Form("h_eta_avg_Deta_dr%u_sign%u",idr+1,ksign+1),";#Delta#eta;#bar{#eta}",200,-4.8,4.8, 100,-2.4,2.4);
                h_eta_avg_pair_eta[idr][ksign] = new TH2D(Form("h_eta_avg_pair_eta_dr%u_sign%u",idr+1,ksign+1),";#eta_{pair};#bar{#eta}",100,-2.4,2.4, 100,-2.4,2.4);
                h_pt1_pt2[idr][ksign] = new TH2D(Form("h_pt1_pt2_dr%u_sign%u",idr+1,ksign+1),";p_{T}^{sublead} [GeV];p_{T}^{lead} [GeV]",pms.npt_bins,pms.pTBins,pms.npt_bins,pms.pTBins);
                h_ptlead_pair_pt[idr][ksign] = new TH2D(Form("h_ptlead_pair_pt_dr%u_sign%u",idr+1,ksign+1),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",pms.npairPT_bins,pms.pairPTBins[ksign][idr],pms.npt_bins,pms.pTBins);
                h_minv_pair_pt[idr][ksign] = new TH2D(Form("h_minv_pair_pt_dr%u_sign%u",idr+1,ksign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",pms.npairPT_bins,pms.pairPTBins[ksign][idr],pms.minv_nbins[idr],0,pms.minv_max[idr]);

                // int nminvbins = pms.minv_bins[isign][idr].size(); // size of a vector
                // double minv_arr[pms.minv_nbins[idr]+1];
                // std::copy(pms.minv_bins[ksign][idr].begin(),pms.minv_bins[ksign][idr].end(),minv_arr);
                // h_minv_pair_pt[idr][ksign] = new TH2D(Form("h_minv_pair_pt_dr%u_sign%u",idr+1,ksign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",pms.npairPT_bins,pms.pairPTBins[ksign][idr],pms.minv_nbins[idr],minv_arr);
    	  	}
        }
    }
    else if (mode == 2){
     	for (unsigned int jctr = 0; jctr < nCtrBins; jctr++){
            for (unsigned int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
   	            for (unsigned int idr = 0; idr < ParamsSet::ndRselcs; idr++){
                    // h_ctrbin_Dphi_failgapcut[idr][jctr][ksign] = new TH1D(Form("h_Dphi_dr%d_ctr%d_sign%d_failgapcut",idr+1,jctr+1,ksign+1),";#Delta#phi;1/N_{evt} dN/d#Delta#phi", static_cast<int>(128/scaleFactorCtrs[jctr]),-pms.PI,pms.PI);
                    // h_ctrbin_pair_eta_Dphi_failgapcut[idr][jctr][ksign] = new TH2D(Form("h_pair_eta_Dphi_dr%d_ctr%d_sign%d_failgapcut",idr+1,jctr+1,ksign+1),";#Delta#phi;#eta_{pair}", static_cast<int>(128/scaleFactorCtrs[jctr]),-pms.PI,pms.PI,20,-1.1 * pms.eta_gap_cut,1.1 * pms.eta_gap_cut);

                    for (unsigned int lgapcut = 0; lgapcut < ParamsSet::nGapCuts; lgapcut++){
                        h_ctrbin_pair_dP_overP[idr][jctr][ksign][lgapcut] = new TH1D(Form("h_pair_dP_overP_dr%d_ctr%d_sign%d_gapcut%d",idr+1,jctr+1,ksign+1,lgapcut+1),";(#Delta p / p)_{pair};1/N_{evt} dN/d(#Delta p / p)_{pair}" ,static_cast<int>(pms.deltaP_overP_nbins/scaleFactorCtrs[jctr]),0,pms.deltaP_overP_max);
                        h_ctrbin_DR[idr][jctr][ksign][lgapcut] = new TH1D(Form("h_DR_dr%d_ctr%d_sign%d_gapcut%d",idr+1,jctr+1,ksign+1,lgapcut+1),";#Delta R;1/N_{evt} dN/d#Delta R", static_cast<int>(pms.deltaR_nbins[idr]/scaleFactorCtrs[jctr]),0,pms.deltaR_thrsh[idr]);
                        h_ctrbin_Dphi[idr][jctr][ksign][lgapcut] = new TH1D(Form("h_Dphi_dr%d_ctr%d_sign%d_gapcut%d",idr+1,jctr+1,ksign+1,lgapcut+1),";#Delta#phi;1/N_{evt} dN/d#Delta#phi", static_cast<int>(128/scaleFactorCtrs[jctr]),-pms.PI,pms.PI);
                        h_ctrbin_pair_y[idr][jctr][ksign][lgapcut] = new TH1D(Form("h_pair_y_dr%d_ctr%d_sign%d_gapcut%d",idr+1,jctr+1,ksign+1,lgapcut+1),";y_{pair};1/N_{evt} dN/dy_{pair}" ,static_cast<int>(90/scaleFactorCtrs[jctr]),-3,3);
                        h_ctrbin_eta_avg_Dphi[idr][jctr][ksign][lgapcut] = new TH2D(Form("h_eta_avg_Dphi_dr%d_ctr%d_sign%d_gapcut%d",idr+1,jctr+1,ksign+1,lgapcut+1),";#Delta#phi;#bar{#eta}", static_cast<int>(128/scaleFactorCtrs[jctr]),-pms.PI,pms.PI,static_cast<int>(100/scaleFactorCtrs[jctr]),-2.4,2.4);
                        h_ctrbin_Deta_Dphi[idr][jctr][ksign][lgapcut] = new TH2D(Form("h_Deta_Dphi_dr%d_ctr%d_sign%d_gapcut%d",idr+1,jctr+1,ksign+1,lgapcut+1),";#Delta#phi;#Delta#eta", static_cast<int>(128/scaleFactorCtrs[jctr]),-pms.PI,pms.PI,static_cast<int>(200/scaleFactorCtrs[jctr]),-4.8,4.8);
                        h_ctrbin_eta1_eta2[idr][jctr][ksign][lgapcut] = new TH2D(Form("h_eta1_eta2_dr%d_ctr%d_sign%d_gapcut%d",idr+1,jctr+1,ksign+1,lgapcut+1),";#eta_{sublead};#eta_{lead}",static_cast<int>(100/scaleFactorCtrs[jctr]),-2.4,2.4, static_cast<int>(100/scaleFactorCtrs[jctr]),-2.4,2.4);
                        h_ctrbin_eta_avg_Deta[idr][jctr][ksign][lgapcut] = new TH2D(Form("h_eta_avg_Deta_dr%d_ctr%d_sign%d_gapcut%d",idr+1,jctr+1,ksign+1,lgapcut+1),";#Delta#eta;#bar{#eta}",static_cast<int>(200/scaleFactorCtrs[jctr]),-4.8,4.8, static_cast<int>(100/scaleFactorCtrs[jctr]),-2.4,2.4);
                        h_ctrbin_eta_avg_pair_eta[idr][jctr][ksign][lgapcut] = new TH2D(Form("h_eta_avg_pair_eta_dr%d_ctr%d_sign%d_gapcut%d",idr+1,jctr+1,ksign+1,lgapcut+1),";#eta_{pair};#bar{#eta}",static_cast<int>(100/scaleFactorCtrs[jctr]),-2.4,2.4, static_cast<int>(100/scaleFactorCtrs[jctr]),-2.4,2.4);
                        h_ctrbin_pt1_pt2[idr][jctr][ksign][lgapcut] = new TH2D(Form("h_pt1_pt2_dr%d_ctr%d_sign%d_gapcut%d",idr+1,jctr+1,ksign+1,lgapcut+1),";p_{T}^{sublead} [GeV];p_{T}^{lead} [GeV]",pms.npt_bins,pms.pTBins,pms.npt_bins,pms.pTBins);
                        h_ctrbin_ptlead_pair_pt[idr][jctr][ksign][lgapcut] = new TH2D(Form("h_ptlead_pair_pt_dr%d_ctr%d_sign%d_gapcut%d",idr+1,jctr+1,ksign+1,lgapcut+1),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",pms.npairPT_bins,pms.pairPTBins[ksign][idr],pms.npt_bins,pms.pTBins);
                        h_ctrbin_minv_pair_pt[idr][jctr][ksign][lgapcut] = new TH2D(Form("h_minv_pair_pt_dr%d_ctr%d_sign%d_gapcut%d",idr+1,jctr+1,ksign+1,lgapcut+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",pms.npairPT_bins,pms.pairPTBins[ksign][idr],pms.minv_nbins[idr],0,pms.minv_max[idr]);
                        
                        // double minv_arr[pms.minv_nbins[idr]+1];
                        // std::copy(pms.minv_bins[ksign][idr].begin(),pms.minv_bins[ksign][idr].end(),minv_arr);
                        // h_ctrbin_minv_pair_pt[idr][jctr][ksign][lgapcut] = new TH2D(Form("h_minv_pair_pt_dr%d_ctr%d_sign%d_gapcut%d",idr+1,jctr+1,ksign+1,lgapcut+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",pms.npairPT_bins,pms.pairPTBins[ksign][idr],pms.minv_nbins[idr],minv_arr);
                    }
                }

                for (unsigned int idphi= 0; idphi < ParamsSet::ndphiselcs; idphi++){
                    for (unsigned int lgapcut = 0; lgapcut < ParamsSet::nGapCuts; lgapcut++){
                        h_ctrbin_eta1_eta2_dphicut[idphi][jctr][ksign][lgapcut] = new TH2D(Form("h_eta1_eta2_dphi%d_ctr%d_sign%d_gapcut%d",idphi+1,jctr+1,ksign+1,lgapcut+1),";#eta_{sublead};#eta_{lead}",static_cast<int>(100/scaleFactorCtrs[jctr]),-2.4,2.4, static_cast<int>(100/scaleFactorCtrs[jctr]),-2.4,2.4);
                    }
                }
            }
        }
    }
}




#endif