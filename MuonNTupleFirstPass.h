#ifndef __MuonNTupleFirstPass_h__
#define __MuonNTupleFirstPass_h__

#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <string.h>
#include <fstream>
// #include <cmath>  
#include "MuonPair.h"
#include "ParamsSet.h"
#include "vector"
#include "TH1D.h"
#include "TH2D.h"

class MuonNTupleFirstPass{
  // Read through the N-tuple, apply appropriate cuts
  // Then fill in histograms and/or output trees
  // mode = 1: raw histograms
  // mode = 2: kill resonances & make no-binning histograms
  // mode = 3: kill resonances & output single-muon information into a TTree
  // mode = 4: kill resonances & output muon-pair information into a TTree
  // mode = 5: kill resonances & make centrality-binned histograms
  // mode = 6: kill resonances & make pT-binned histograms


private:
  // --------------------- general settings ---------------------------

  int mode = 3;
  ParamsSet pms;
  // const ParamsSet pms;
  // int scaleFactorCtrs[ParamsSet::nCtrBins] = {1,1,3,5,5};
  int scaleFactorPts[ParamsSet::nPtBins] = {1,1,1,1,1};

    // --------------------- input files & trees & data for setting branches---------------------------

  TChain          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   UInt_t          RunNumber;
   Float_t         ZdcEtA;
   Float_t         ZdcEtC;
   Float_t         FCal_Et;
   Int_t           centrality;
   std::vector<float>   *muon_deltaP_overP           =nullptr;

   std::vector<int>     *muon_pair_muon1_index       =nullptr;
   std::vector<float>   *muon_pair_muon1_pt          =nullptr;
   std::vector<float>   *muon_pair_muon1_eta         =nullptr;
   std::vector<float>   *muon_pair_muon1_phi         =nullptr;
   std::vector<int>     *muon_pair_muon1_quality     =nullptr;
   std::vector<float>   *muon_pair_muon1_d0          =nullptr;
   std::vector<float>   *muon_pair_muon1_z0          =nullptr;

   std::vector<int>     *muon_pair_muon2_index       =nullptr;
   std::vector<float>   *muon_pair_muon2_pt          =nullptr;
   std::vector<float>   *muon_pair_muon2_eta         =nullptr;
   std::vector<float>   *muon_pair_muon2_phi         =nullptr;
   std::vector<int>     *muon_pair_muon2_quality     =nullptr;
   std::vector<float>   *muon_pair_muon2_d0          =nullptr;
   std::vector<float>   *muon_pair_muon2_z0          =nullptr;

   Bool_t               b_HLT_mu4_mu4noL1;
   std::vector<bool>    *dimuon_b_HLT_mu4_mu4noL1    =nullptr;

     // --------------------- temporary muon and muonpair objects ---------------------------

   Muon* tempmuon = nullptr;
   MuonPair* mpair = nullptr;

    // --------------------- class methods ---------------------------

   void Init();
   void InitOutput();
   void InitHists();
   void ProcessData();
   bool PassCuts();
   void FillUnbinnedHistograms();
   void FillCtrBinnedHistograms();
   void FillPtBinnedHistograms();
   void FillSingleMuonTree();
   void FillMuonPairTree();
   void WriteOutput();


  // --------------------- output file, histograms & trees ---------------------------

   TFile *m_outfile = nullptr;

   TH1D* h_FCal_Et[ParamsSet::ndRselcs];
   // TH1D* h_deltaP_overP[ParamsSet::ndRselcs];
   TH1D* h_pair_dP_overP[ParamsSet::ndRselcs][ParamsSet::nSigns];
   TH1D* h_Dphi[ParamsSet::ndRselcs][ParamsSet::nSigns];
   TH1D* h_Deta[ParamsSet::ndRselcs][ParamsSet::nSigns];
   TH1D* h_DR[ParamsSet::ndRselcs][ParamsSet::nSigns];
   TH1D* h_Minv[ParamsSet::ndRselcs][ParamsSet::nSigns];
   TH1D* h_pt_lead[ParamsSet::ndRselcs][ParamsSet::nSigns];
   TH1D* h_eta_lead[ParamsSet::ndRselcs][ParamsSet::nSigns];
   TH1D* h_eta_sublead[ParamsSet::ndRselcs][ParamsSet::nSigns];
   TH1D* h_eta_avg[ParamsSet::ndRselcs][ParamsSet::nSigns];
   TH1D* h_pair_pt[ParamsSet::ndRselcs][ParamsSet::nSigns];
   TH1D* h_pair_eta[ParamsSet::ndRselcs][ParamsSet::nSigns];
   TH1D* h_pair_y[ParamsSet::ndRselcs][ParamsSet::nSigns];
   TH2D* h_eta_phi[ParamsSet::ndRselcs][ParamsSet::nSigns];
   TH2D* h_eta1_eta2[ParamsSet::ndRselcs][ParamsSet::nSigns];
   TH2D* h_eta_avg_Deta[ParamsSet::ndRselcs][ParamsSet::nSigns];
   TH2D* h_eta_avg_pair_eta[ParamsSet::ndRselcs][ParamsSet::nSigns];
   TH2D* h_pt1_pt2[ParamsSet::ndRselcs][ParamsSet::nSigns];
   TH2D* h_ptlead_pair_pt[ParamsSet::ndRselcs][ParamsSet::nSigns];
   TH2D* h_minv_pair_pt[ParamsSet::ndRselcs][ParamsSet::nSigns];


   TH1D* h_ctrbin_pair_dP_overP[ParamsSet::ndRselcs][ParamsSet::nCtrBins][ParamsSet::nSigns];
   TH1D* h_ctrbin_Dphi[ParamsSet::ndRselcs][ParamsSet::nCtrBins][ParamsSet::nSigns];
   TH1D* h_ctrbin_Deta[ParamsSet::ndRselcs][ParamsSet::nCtrBins][ParamsSet::nSigns];
   TH1D* h_ctrbin_DR[ParamsSet::ndRselcs][ParamsSet::nCtrBins][ParamsSet::nSigns];
   TH1D* h_ctrbin_Minv[ParamsSet::ndRselcs][ParamsSet::nCtrBins][ParamsSet::nSigns];
   TH1D* h_ctrbin_pt_lead[ParamsSet::ndRselcs][ParamsSet::nCtrBins][ParamsSet::nSigns];
   TH1D* h_ctrbin_eta_lead[ParamsSet::ndRselcs][ParamsSet::nCtrBins][ParamsSet::nSigns];
   TH1D* h_ctrbin_eta_sublead[ParamsSet::ndRselcs][ParamsSet::nCtrBins][ParamsSet::nSigns];
   TH1D* h_ctrbin_eta_avg[ParamsSet::ndRselcs][ParamsSet::nCtrBins][ParamsSet::nSigns];
   TH1D* h_ctrbin_pair_pt[ParamsSet::ndRselcs][ParamsSet::nCtrBins][ParamsSet::nSigns];
   TH1D* h_ctrbin_pair_eta[ParamsSet::ndRselcs][ParamsSet::nCtrBins][ParamsSet::nSigns];
   TH1D* h_ctrbin_pair_y[ParamsSet::ndRselcs][ParamsSet::nCtrBins][ParamsSet::nSigns];
   TH2D* h_ctrbin_eta_phi[ParamsSet::ndRselcs][ParamsSet::nCtrBins][ParamsSet::nSigns];
   TH2D* h_ctrbin_eta1_eta2[ParamsSet::ndRselcs][ParamsSet::nCtrBins][ParamsSet::nSigns];
   TH2D* h_ctrbin_eta_avg_Deta[ParamsSet::ndRselcs][ParamsSet::nCtrBins][ParamsSet::nSigns];
   TH2D* h_ctrbin_eta_avg_pair_eta[ParamsSet::ndRselcs][ParamsSet::nCtrBins][ParamsSet::nSigns];
   TH2D* h_ctrbin_pt1_pt2[ParamsSet::ndRselcs][ParamsSet::nCtrBins][ParamsSet::nSigns];
   TH2D* h_ctrbin_ptlead_pair_pt[ParamsSet::ndRselcs][ParamsSet::nCtrBins][ParamsSet::nSigns];
   TH2D* h_ctrbin_minv_pair_pt[ParamsSet::ndRselcs][ParamsSet::nCtrBins][ParamsSet::nSigns];


   TH1D* h_ptbin_pair_dP_overP[ParamsSet::ndRselcs][ParamsSet::nPtBins][ParamsSet::nSigns];
   TH1D* h_ptbin_Dphi[ParamsSet::ndRselcs][ParamsSet::nPtBins][ParamsSet::nSigns];
   TH1D* h_ptbin_Deta[ParamsSet::ndRselcs][ParamsSet::nPtBins][ParamsSet::nSigns];
   TH1D* h_ptbin_DR[ParamsSet::ndRselcs][ParamsSet::nPtBins][ParamsSet::nSigns];
   TH1D* h_ptbin_Minv[ParamsSet::ndRselcs][ParamsSet::nPtBins][ParamsSet::nSigns];
   TH1D* h_ptbin_pt_lead[ParamsSet::ndRselcs][ParamsSet::nPtBins][ParamsSet::nSigns];
   TH1D* h_ptbin_eta_lead[ParamsSet::ndRselcs][ParamsSet::nPtBins][ParamsSet::nSigns];
   TH1D* h_ptbin_eta_sublead[ParamsSet::ndRselcs][ParamsSet::nPtBins][ParamsSet::nSigns];
   TH1D* h_ptbin_eta_avg[ParamsSet::ndRselcs][ParamsSet::nPtBins][ParamsSet::nSigns];
   TH1D* h_ptbin_pair_pt[ParamsSet::ndRselcs][ParamsSet::nPtBins][ParamsSet::nSigns];
   TH1D* h_ptbin_pair_eta[ParamsSet::ndRselcs][ParamsSet::nPtBins][ParamsSet::nSigns];
   TH1D* h_ptbin_pair_y[ParamsSet::ndRselcs][ParamsSet::nPtBins][ParamsSet::nSigns];
   TH2D* h_ptbin_eta_phi[ParamsSet::ndRselcs][ParamsSet::nPtBins][ParamsSet::nSigns];
   TH2D* h_ptbin_eta1_eta2[ParamsSet::ndRselcs][ParamsSet::nPtBins][ParamsSet::nSigns];
   TH2D* h_ptbin_eta_avg_Deta[ParamsSet::ndRselcs][ParamsSet::nPtBins][ParamsSet::nSigns];
   TH2D* h_ptbin_eta_avg_pair_eta[ParamsSet::ndRselcs][ParamsSet::nPtBins][ParamsSet::nSigns];
   TH2D* h_ptbin_pt1_pt2[ParamsSet::ndRselcs][ParamsSet::nPtBins][ParamsSet::nSigns];
   TH2D* h_ptbin_ptlead_pair_pt[ParamsSet::ndRselcs][ParamsSet::nPtBins][ParamsSet::nSigns];
   TH2D* h_ptbin_minv_pair_pt[ParamsSet::ndRselcs][ParamsSet::nPtBins][ParamsSet::nSigns];


   TTree* muonOutTree;
   TTree* muonOutTreeBinned[ParamsSet::nCtrBins];
   TTree* muonPairOutTree[ParamsSet::nSigns];
   TTree* muonPairOutTreeBinned[ParamsSet::ndRselcs][ParamsSet::nCtrBins][ParamsSet::nSigns];


public :
   MuonNTupleFirstPass();
   ~MuonNTupleFirstPass(){}
   void Run();
};


MuonNTupleFirstPass::MuonNTupleFirstPass(){
  std::cout << "Sanity test ensuring the constructor is called properly." << std::endl;
  std::cout << "Before initiation; mode = " << mode << std::endl;
}

//initialize the TChain
void MuonNTupleFirstPass::Init(){

   fChain = new TChain("HeavyIonD3PD","HeavyIonD3PD");
   fChain->SetMakeClass(1);
   fChain->Add("all_Data2018_12March2022.root");
   fChain->Add("all_Data2015_12March2022.root"); 

   fChain->SetBranchAddress("RunNumber"                   , &RunNumber);
   fChain->SetBranchAddress("ZdcEtA"                      , &ZdcEtA              );
   fChain->SetBranchAddress("ZdcEtC"                      , &ZdcEtC              );
   fChain->SetBranchAddress("FCal_Et"                     , &FCal_Et);
   fChain->SetBranchAddress("centrality"                  , &centrality);
   fChain->SetBranchAddress("muon_deltaP_overP"           , &muon_deltaP_overP);

   fChain->SetBranchAddress("muon_pair_muon1_index"       , &muon_pair_muon1_index);
   fChain->SetBranchAddress("muon_pair_muon1_pt"          , &muon_pair_muon1_pt);
   fChain->SetBranchAddress("muon_pair_muon1_eta"         , &muon_pair_muon1_eta);
   fChain->SetBranchAddress("muon_pair_muon1_phi"         , &muon_pair_muon1_phi);
   fChain->SetBranchAddress("muon_pair_muon1_quality"     , &muon_pair_muon1_quality);
   fChain->SetBranchAddress("muon_pair_muon1_d0"          , &muon_pair_muon1_d0);
   fChain->SetBranchAddress("muon_pair_muon1_z0"          , &muon_pair_muon1_z0);

   fChain->SetBranchAddress("muon_pair_muon2_index"       , &muon_pair_muon2_index);
   fChain->SetBranchAddress("muon_pair_muon2_pt"          , &muon_pair_muon2_pt);
   fChain->SetBranchAddress("muon_pair_muon2_eta"         , &muon_pair_muon2_eta);
   fChain->SetBranchAddress("muon_pair_muon2_phi"         , &muon_pair_muon2_phi);
   fChain->SetBranchAddress("muon_pair_muon2_quality"     , &muon_pair_muon2_quality);
   fChain->SetBranchAddress("muon_pair_muon2_d0"          , &muon_pair_muon2_d0);
   fChain->SetBranchAddress("muon_pair_muon2_z0"          , &muon_pair_muon2_z0);
   fChain->SetBranchAddress("b_HLT_mu4_mu4noL1"           , &b_HLT_mu4_mu4noL1);
   fChain->SetBranchAddress("dimuon_b_HLT_mu4_mu4noL1"    , &dimuon_b_HLT_mu4_mu4noL1);


   //SetBranch Status
   fChain->SetBranchStatus("*"                               ,0);//switch off all branches, then enable just the ones that we need
   fChain->SetBranchStatus("RunNumber"                       ,1);
   fChain->SetBranchStatus("ZdcEtA"                          ,1);
   fChain->SetBranchStatus("ZdcEtC"                          ,1);
   fChain->SetBranchStatus("FCal_Et"                         ,1);
   fChain->SetBranchStatus("centrality"                      ,1);
   fChain->SetBranchStatus("muon_deltaP_overP"               ,1);

   fChain->SetBranchStatus("muon_pair_muon1_index"           ,1);
   fChain->SetBranchStatus("muon_pair_muon1_pt"              ,1);
   fChain->SetBranchStatus("muon_pair_muon1_eta"             ,1);
   fChain->SetBranchStatus("muon_pair_muon1_phi"             ,1);
   fChain->SetBranchStatus("muon_pair_muon1_quality"         ,1);
   fChain->SetBranchStatus("muon_pair_muon1_d0"              ,1);
   fChain->SetBranchStatus("muon_pair_muon1_z0"              ,1);

   fChain->SetBranchStatus("muon_pair_muon2_index"           ,1);
   fChain->SetBranchStatus("muon_pair_muon2_pt"              ,1);
   fChain->SetBranchStatus("muon_pair_muon2_eta"             ,1);
   fChain->SetBranchStatus("muon_pair_muon2_phi"             ,1);
   fChain->SetBranchStatus("muon_pair_muon2_quality"         ,1);
   fChain->SetBranchStatus("muon_pair_muon2_d0"              ,1);
   fChain->SetBranchStatus("muon_pair_muon2_z0"              ,1);

   fChain->SetBranchStatus("b_HLT_mu4_mu4noL1"               ,1);
   fChain->SetBranchStatus("dimuon_b_HLT_mu4_mu4noL1"        ,1);
}


void MuonNTupleFirstPass::InitHists(){
  if (mode == 3 || mode == 4) return; // output tree only (no histogram)

  if (mode == 1 || mode == 2){
    for (unsigned int idr = 0; idr < ParamsSet::ndRselcs; idr++){
      h_FCal_Et[idr]=new TH1D(Form("h_FCal_Et_dr%u",idr+1),";FCal #it{E}_{T} [TeV];1/N_{evt} dN/dFCal #it{E}_{T}",550,-0.5,5.0);
      // h_deltaP_overP[idr] =new TH1D(Form("h_deltaP_overP_dr%u",idr+1),";(p_{ID}-p_{MS}-#Delta p_{calo}) / p_{ID};1/N_{evt} dN/d" ,pms.deltaP_overP_nbins*2,-pms.deltaP_overP_thrsh,pms.deltaP_overP_thrsh);
      
      for (unsigned int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
        h_pair_dP_overP[idr][ksign] = new TH1D(Form("h_pair_dP_overP_dr%u_sign%u",idr+1,ksign+1),";(#Delta p / p)_{pair};1/N_{evt} dN/d" ,pms.deltaP_overP_nbins,0,pms.deltaP_overP_max);
        h_Dphi[idr][ksign] = new TH1D(Form("h_Dphi_dr%u_sign%u",idr+1,ksign+1),";#Delta#phi;1/N_{evt} dN/d#Delta#phi" ,128,-pms.PI ,pms.PI );
        h_Deta[idr][ksign] = new TH1D(Form("h_Deta_dr%u_sign%u",idr+1,ksign+1),";#Delta#eta;1/N_{evt} dN/d#Delta#eta" ,100,-2.4,2.4);
        h_DR[idr][ksign] = new TH1D(Form("h_DR_dr%u_sign%u",idr+1,ksign+1),";#Delta R;1/N_{evt} dN/d#Delta R", pms.deltaR_nbins[idr],0,pms.deltaR_thrsh[idr]);
        h_Minv[idr][ksign] = new TH1D(Form("h_Minv_dr%u_sign%u",idr+1,ksign+1),";m_{#mu#mu} [GeV];1/N_{evt} dN/dm_{#mu#mu}", 300, 0, 15 );
        h_pt_lead[idr][ksign] = new TH1D(Form("h_pt_lead_dr%u_sign%u",idr+1,ksign+1),";p_{T}^{lead} [GeV];1/N_{evt} dN/dp_{T}^{lead}" ,pms.npt_bins,pms.pTBins );
        h_eta_lead[idr][ksign] = new TH1D(Form("h_eta_lead_dr%u_sign%u",idr+1,ksign+1),";#eta_{lead};1/N_{evt} dN/d#eta_{lead}" ,100,-2.4,2.4);
        h_eta_sublead[idr][ksign] = new TH1D(Form("h_eta_sublead_dr%u_sign%u",idr+1,ksign+1),";#eta_{sublead};1/N_{evt} dN/d#eta_{sublead}",100,-2.4,2.4);
        h_eta_avg[idr][ksign] = new TH1D(Form("h_eta_avg_dr%u_sign%u",idr+1,ksign+1),";#bar{#eta};1/N_{evt} dN/d#bar{#eta}",100,-2.4,2.4);
        h_pair_pt[idr][ksign] = new TH1D(Form("h_pair_pt_dr%u_sign%u",idr+1,ksign+1),";p_{T}^{pair} [GeV];1/N_{evt} dN/dp_{T}^{pair}" ,pms.npairPT_bins,pms.pairPTBins );
        h_pair_eta[idr][ksign] = new TH1D(Form("h_pair_eta_dr%u_sign%u",idr+1,ksign+1),";#eta_{pair};1/N_{evt} dN/d#eta_{pair}" ,100,-2.4,2.4);
        h_pair_y[idr][ksign] = new TH1D(Form("h_pair_y_dr%u_sign%u",idr+1,ksign+1),";y_{pair};1/N_{evt} dN/dy_{pair}" ,90,-3,3);
        h_eta_phi[idr][ksign] = new TH2D(Form("h_eta_phi_dr%u_sign%u",idr+1,ksign+1),";#Delta#phi;#bar{#eta}",128,-pms.PI,pms.PI,100,-2.4,2.4);
        h_eta1_eta2[idr][ksign] = new TH2D(Form("h_eta1_eta2_dr%u_sign%u",idr+1,ksign+1),";#eta_{sublead};#eta_{lead}",100,-2.4,2.4, 100,-2.4,2.4);
        h_eta_avg_Deta[idr][ksign] = new TH2D(Form("h_eta_avg_Deta_dr%u_sign%u",idr+1,ksign+1),";#Delta#eta;#bar{#eta}",100,-2.4,2.4, 100,-2.4,2.4);
        h_eta_avg_pair_eta[idr][ksign] = new TH2D(Form("h_eta_avg_pair_eta_dr%u_sign%u",idr+1,ksign+1),";#eta_{pair};#bar{#eta}",100,-2.4,2.4, 100,-2.4,2.4);
        h_pt1_pt2[idr][ksign] = new TH2D(Form("h_pt1_pt2_dr%u_sign%u",idr+1,ksign+1),";p_{T}^{sublead} [GeV];p_{T}^{lead} [GeV]",pms.npt_bins,pms.pTBins,pms.npt_bins,pms.pTBins);
        h_ptlead_pair_pt[idr][ksign] = new TH2D(Form("h_ptlead_pair_pt_dr%u_sign%u",idr+1,ksign+1),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",pms.npairPT_bins,pms.pairPTBins,pms.npt_bins,pms.pTBins);
        h_minv_pair_pt[idr][ksign] = new TH2D(Form("h_minv_pair_pt_dr%u_sign%u",idr+1,ksign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",pms.npairPT_bins,pms.pairPTBins,2000,0   ,100);
      }
    }
  }
  else if (mode == 5){
    for (unsigned int idr = 0; idr < ParamsSet::ndRselcs; idr++){
      for (unsigned int jctr = 0; jctr < ParamsSet::nCtrBins; jctr++){
        for (unsigned int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
          h_ctrbin_pair_dP_overP[idr][jctr][ksign] = new TH1D(Form("h_pair_dP_overP_dr%d_ctr%d_sign%d",idr+1,jctr+1,ksign+1),";(#Delta p / p)_{pair};1/N_{evt} dN/d(#Delta p / p)_{pair}" ,static_cast<int>(pms.deltaP_overP_nbins/pms.scaleFactorCtrs[jctr]),0,pms.deltaP_overP_max);
          h_ctrbin_Dphi[idr][jctr][ksign] = new TH1D(Form("h_Dphi_dr%d_ctr%d_sign%d",idr+1,jctr+1,ksign+1),";#Delta#phi;1/N_{evt} dN/d#Delta#phi" ,static_cast<int>(128/pms.scaleFactorCtrs[jctr]) ,-pms.PI ,pms.PI );
          h_ctrbin_Deta[idr][jctr][ksign] = new TH1D(Form("h_Deta_dr%d_ctr%d_sign%d",idr+1,jctr+1,ksign+1),";#Delta#eta;1/N_{evt} dN/d#Delta#eta" ,static_cast<int>(100/pms.scaleFactorCtrs[jctr]),-2.4,2.4);
          h_ctrbin_DR[idr][jctr][ksign] = new TH1D(Form("h_DR_dr%d_ctr%d_sign%d",idr+1,jctr+1,ksign+1),";#Delta R;1/N_{evt} dN/d#Delta R", static_cast<int>(pms.deltaR_nbins[idr]/pms.scaleFactorCtrs[jctr]),0,pms.deltaR_thrsh[idr]);
          h_ctrbin_Minv[idr][jctr][ksign] = new TH1D(Form("h_Minv_dr%d_ctr%d_sign%d",idr+1,jctr+1,ksign+1),";m_{#mu#mu} [GeV];1/N_{evt} dN/dm_{#mu#mu}", static_cast<int>(pms.minv_nbins[idr]/pms.scaleFactorCtrs[jctr]),0   ,pms.minv_max[idr]);
          h_ctrbin_pt_lead[idr][jctr][ksign] = new TH1D(Form("h_pt_lead_dr%d_ctr%d_sign%d",idr+1,jctr+1,ksign+1),";p_{T}^{lead} [GeV];1/N_{evt} dN/dp_{T}^{lead}" ,pms.npt_bins,pms.pTBins );
          // h_ctrbin_eta_lead[idr][jctr][ksign] = new TH1D(Form("h_eta_lead_dr%d_ctr%d_sign%d",idr+1,jctr+1,ksign+1),";#eta_{lead};1/N_{evt} dN/d#eta_{lead}" ,static_cast<int>(100/pms.scaleFactorCtrs[jctr]),-2.4,2.4);
          // h_ctrbin_eta_sublead[idr][jctr][ksign] = new TH1D(Form("h_eta_sublead_dr%d_ctr%d_sign%d",idr+1,jctr+1,ksign+1),";#eta_{sublead};1/N_{evt} dN/d#eta_{sublead}",static_cast<int>(100/pms.scaleFactorCtrs[jctr]),-2.4,2.4);
          h_ctrbin_eta_avg[idr][jctr][ksign] = new TH1D(Form("h_eta_avg_dr%d_ctr%d_sign%d",idr+1,jctr+1,ksign+1),";#bar{#eta};1/N_{evt} dN/d#bar{#eta}",static_cast<int>(100/pms.scaleFactorCtrs[jctr]),-2.4,2.4);
          h_ctrbin_pair_pt[idr][jctr][ksign] = new TH1D(Form("h_pair_pt_dr%d_ctr%d_sign%d",idr+1,jctr+1,ksign+1),";p_{T}^{pair} [GeV];1/N_{evt} dN/dp_{T}^{pair}" ,pms.npairPT_bins,pms.pairPTBins );
          h_ctrbin_pair_eta[idr][jctr][ksign] = new TH1D(Form("h_pair_eta_dr%d_ctr%d_sign%d",idr+1,jctr+1,ksign+1),";#eta_{pair};1/N_{evt} dN/d#eta_{pair}" ,static_cast<int>(100/pms.scaleFactorCtrs[jctr]),-2.4,2.4);
          h_ctrbin_pair_y[idr][jctr][ksign] = new TH1D(Form("h_pair_y_dr%d_ctr%d_sign%d",idr+1,jctr+1,ksign+1),";y_{pair};1/N_{evt} dN/dy_{pair}" ,static_cast<int>(90/pms.scaleFactorCtrs[jctr]),-3,3);
          h_ctrbin_eta_phi[idr][jctr][ksign] = new TH2D(Form("h_eta_phi_dr%d_ctr%d_sign%d",idr+1,jctr+1,ksign+1),";#Delta#phi;#bar{#eta}", static_cast<int>(128/pms.scaleFactorCtrs[jctr]),-pms.PI,pms.PI,static_cast<int>(100/pms.scaleFactorCtrs[jctr]),-2.4,2.4);
          h_ctrbin_eta1_eta2[idr][jctr][ksign] = new TH2D(Form("h_eta1_eta2_dr%d_ctr%d_sign%d",idr+1,jctr+1,ksign+1),";#eta_{sublead};#eta_{lead}",static_cast<int>(100/pms.scaleFactorCtrs[jctr]),-2.4,2.4, static_cast<int>(100/pms.scaleFactorCtrs[jctr]),-2.4,2.4);
          h_ctrbin_eta_avg_Deta[idr][jctr][ksign] = new TH2D(Form("h_eta_avg_Deta_dr%d_ctr%d_sign%d",idr+1,jctr+1,ksign+1),";#Delta#eta;#bar{#eta}",static_cast<int>(100/pms.scaleFactorCtrs[jctr]),-2.4,2.4, static_cast<int>(100/pms.scaleFactorCtrs[jctr]),-2.4,2.4);
          h_ctrbin_eta_avg_pair_eta[idr][jctr][ksign] = new TH2D(Form("h_eta_avg_pair_eta_dr%d_ctr%d_sign%d",idr+1,jctr+1,ksign+1),";#eta_{pair};#bar{#eta}",static_cast<int>(100/pms.scaleFactorCtrs[jctr]),-2.4,2.4, static_cast<int>(100/pms.scaleFactorCtrs[jctr]),-2.4,2.4);
          h_ctrbin_pt1_pt2[idr][jctr][ksign] = new TH2D(Form("h_pt1_pt2_dr%d_ctr%d_sign%d",idr+1,jctr+1,ksign+1),";p_{T}^{sublead} [GeV];p_{T}^{lead} [GeV]",static_cast<int>(pms.npt_bins/pms.scaleFactorCtrs[jctr]),pms.pTBins,static_cast<int>(pms.npt_bins/pms.scaleFactorCtrs[jctr]),pms.pTBins);
          h_ctrbin_ptlead_pair_pt[idr][jctr][ksign] = new TH2D(Form("h_ptlead_pair_pt_dr%d_ctr%d_sign%d",idr+1,jctr+1,ksign+1),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",static_cast<int>(pms.npt_bins/pms.scaleFactorCtrs[jctr]),pms.pTBins,static_cast<int>(pms.npt_bins/pms.scaleFactorCtrs[jctr]),pms.pTBins);
          h_ctrbin_minv_pair_pt[idr][jctr][ksign] = new TH2D(Form("h_minv_pair_pt_dr%d_ctr%d_sign%d",idr+1,jctr+1,ksign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",static_cast<int>(pms.npairPT_bins/pms.scaleFactorCtrs[jctr]),pms.pairPTBins,static_cast<int>(pms.minv_nbins[idr]/pms.scaleFactorCtrs[jctr]),0,pms.minv_max[idr]);
        }
      }
    }
  }
  else if (mode == 6){
    for (unsigned int idr = 0; idr < ParamsSet::ndRselcs; idr++){
      for (unsigned int jpt = 0; jpt < ParamsSet::nPtBins; jpt++){
        for (unsigned int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
          h_ptbin_pair_dP_overP[idr][jpt][ksign] = new TH1D(Form("h_pair_dP_overP_dr%d_pt%d_sign%d",idr+1,jpt+1,ksign+1),";(p_{ID}-p_{MS}-#Delta p_{calo}) / p_{ID};1/N_{evt} dN/d" ,pms.deltaP_overP_nbins,0,pms.deltaP_overP_max);
          h_ptbin_Dphi[idr][jpt][ksign] = new TH1D(Form("h_Dphi_dr%d_pt%d_sign%d",idr+1,jpt+1,ksign+1),";#Delta#phi;1/N_{evt} dN/d#Delta#phi" ,128 ,-pms.PI ,pms.PI );
          h_ptbin_Deta[idr][jpt][ksign] = new TH1D(Form("h_Deta_dr%d_pt%d_sign%d",idr+1,jpt+1,ksign+1),";#Delta#eta;1/N_{evt} dN/d#Delta#eta" ,100,-2.4,2.4);
          h_ptbin_DR[idr][jpt][ksign] = new TH1D(Form("h_DR_dr%d_pt%d_sign%d",idr+1,jpt+1,ksign+1),";#Delta R;1/N_{evt} dN/d#Delta R", pms.deltaR_nbins[idr],0,pms.deltaR_thrsh[idr]);
          h_ptbin_Minv[idr][jpt][ksign] = new TH1D(Form("h_Minv_dr%d_pt%d_sign%d",idr+1,jpt+1,ksign+1),";m_{#mu#mu} [GeV];1/N_{evt} dN/dm_{#mu#mu}" ,300,0   ,15 );
          h_ptbin_pt_lead[idr][jpt][ksign] = new TH1D(Form("h_pt_lead_dr%d_pt%d_sign%d",idr+1,jpt+1,ksign+1),";p_{T}^{lead} [GeV];1/N_{evt} dN/dp_{T}^{lead}" ,ParamsSet::nPtBins,pms.pTBins );
          h_ptbin_eta_lead[idr][jpt][ksign] = new TH1D(Form("h_eta_lead_dr%d_pt%d_sign%u",idr+1,jpt+1,ksign+1),";#eta_{lead};1/N_{evt} dN/d#eta_{lead}" ,100,-2.4,2.4);
          h_ptbin_eta_sublead[idr][jpt][ksign] = new TH1D(Form("h_eta_sublead_dr%d_pt%d_sign%u",idr+1,jpt+1,ksign+1),";#eta_{sublead};1/N_{evt} dN/d#eta_{sublead}",100,-2.4,2.4);
          h_ptbin_eta_avg[idr][jpt][ksign] = new TH1D(Form("h_eta_avg_dr%d_pt%d_sign%d",idr+1,jpt+1,ksign+1),";#bar{#eta};1/N_{evt} dN/d#bar{#eta}",100,-2.4,2.4);
          h_ptbin_pair_pt[idr][jpt][ksign] = new TH1D(Form("h_pair_pt_dr%d_pt%d_sign%d",idr+1,jpt+1,ksign+1),";p_{T}^{pair} [GeV];1/N_{evt} dN/dp_{T}^{pair}" ,pms.npairPT_bins,pms.pairPTBins );
          h_ptbin_pair_eta[idr][jpt][ksign] = new TH1D(Form("h_pair_eta_dr%d_pt%d_sign%d",idr+1,jpt+1,ksign+1),";#eta_{pair};1/N_{evt} dN/d#eta_{pair}" ,100,-2.4,2.4);
          h_ptbin_pair_y[idr][jpt][ksign] = new TH1D(Form("h_pair_y_dr%d_pt%d_sign%d",idr+1,jpt+1,ksign+1),";y_{pair};1/N_{evt} dN/dy_{pair}" ,90,-3,3);
          h_ptbin_eta_phi[idr][jpt][ksign] = new TH2D(Form("h_eta_phi_dr%d_pt%d_sign%d",idr+1,jpt+1,ksign+1),";#Delta#phi;#bar{#eta}", 128,-pms.PI,pms.PI,100,-2.4,2.4);
          h_ptbin_eta1_eta2[idr][jpt][ksign] = new TH2D(Form("h_eta1_eta2_dr%d_pt%d_sign%d",idr+1,jpt+1,ksign+1),";#eta_{sublead};#eta_{lead}",100,-2.4,2.4, 100,-2.4,2.4);
          h_ptbin_eta_avg_Deta[idr][jpt][ksign] = new TH2D(Form("h_eta_avg_Deta_dr%d_pt%d_sign%d",idr+1,jpt+1,ksign+1),";#Delta#eta;#bar{#eta}",100,-2.4,2.4, 100,-2.4,2.4);
          h_ptbin_eta_avg_pair_eta[idr][jpt][ksign] = new TH2D(Form("h_eta_avg_pair_eta_dr%d_pt%d_sign%d",idr+1,jpt+1,ksign+1),";#eta_{pair};#bar{#eta}",100,-2.4,2.4, 100,-2.4,2.4);
          h_ptbin_pt1_pt2[idr][jpt][ksign] = new TH2D(Form("h_pt1_pt2_dr%d_pt%d_sign%d",idr+1,jpt+1,ksign+1),";p_{T}^{sublead} [GeV];p_{T}^{lead} [GeV]",ParamsSet::nPtBins,pms.pTBins,ParamsSet::nPtBins,pms.pTBins);
          h_ptbin_ptlead_pair_pt[idr][jpt][ksign] = new TH2D(Form("h_ptlead_pair_pt_dr%d_pt%d_sign%d",idr+1,jpt+1,ksign+1),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",pms.npairPT_bins,pms.pairPTBins,ParamsSet::nPtBins,pms.pTBins);
          h_ptbin_minv_pair_pt[idr][jpt][ksign] = new TH2D(Form("h_minv_pair_pt_dr%d_pt%d_sign%d",idr+1,jpt+1,ksign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",pms.npairPT_bins,pms.pairPTBins,static_cast<int>(2000/scaleFactorPts[jpt]),0   ,100);
        }
      }
    }
  }
}

void MuonNTupleFirstPass::InitOutput(){
  // m_outfile=new TFile("analysis_kill_resonances.root","recreate");

  if(mode == 1){
    m_outfile=new TFile("histograms_raw.root","recreate");
  }
  else if (mode == 2 || mode == 5 || mode == 6){
    // m_outfile=new TFile("histograms_kill_resonances.root","update");
    m_outfile=new TFile("histograms_kill_resonances_TRYPMS.root","update");
  }
  else if (mode == 3){
    // m_outfile=new TFile("single_muon_trees.root","recreate");
    m_outfile=new TFile("single_muon_trees_CHECK.root","recreate");
    muonOutTree = new TTree("muon_tree","all single muons");
    muonOutTree->Branch("MuonObj",&tempmuon);

    for (int jctr = 0; jctr < ParamsSet::nCtrBins; jctr++){
      muonOutTreeBinned[jctr] = new TTree(Form("muon_tree_ctr%d",jctr+1),Form("all muons, centrality bin %d",jctr+1));
      muonOutTreeBinned[jctr]->Branch("MuonObj",&tempmuon);
    }
  }
  else if (mode == 4){ //output muon pair trees
    m_outfile=new TFile("muon_pairs.root","recreate");

    for (unsigned int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
      muonPairOutTree[ksign] = new TTree(Form("muon_pair_tree_sign%u",ksign+1),Form("all muon pairs, sign%u",ksign+1));
      muonPairOutTree[ksign]->Branch("MuonPairObj",&mpair);
      for (unsigned int idr = 0; idr < ParamsSet::ndRselcs; idr++){
        for (unsigned int jctr = 0; jctr < ParamsSet::nCtrBins; jctr++){
          muonPairOutTreeBinned[idr][jctr][ksign] = new TTree(Form("muon_pair_tree_dr%u_ctr%u_sign%u",idr+1,jctr+1,ksign+1),Form("all muon pairs, dr%u, ctr%u, sign%u",idr+1,jctr+1,ksign+1));
          muonPairOutTreeBinned[idr][jctr][ksign]->Branch("MuonPairObj",&mpair);
        }
      }
    }
  }
}

#endif


