#ifndef ScrambSampleGen_h
#define ScrambSampleGen_h
// #define nSigns 2
// #define nScramb 10000000 //Scramb sample size for each centrality bin
// #define nScramb 50000 //Scramb sample size for each centrality bin

//#include "struct_muon.h" 
#include "MuonPair.h"
#include "ParamsSet.h"
#include <TROOT.h>
//#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
//#include <string.h>
//#include <fstream>
#include "vector"
#include "TH1D.h"
#include "TH2D.h"
#include <stdlib.h>

class ScrambSampleGen{
private:
	//static const unsigned int ParamsSet::ndRselcs = 3;
   // static const unsigned int nCtrBins = 5;

   //TChain          *fChain[ParamsSet::ndRselcs];


   // int nScramb[nCtrBins];
   int nScramb[nCtrBins][ParamsSet::nSigns] = {{6417170, 7166880}, {2041080, 2584180}, {431020, 693860}, {61550, 151490}, {7120, 488690}};

   ParamsSet pms;
   // static const unsigned int ctr_step = 5;
   static const unsigned int nctr_intvls = 20;
   static const int nCtrBins = 6;
   // 0-5, 5-10, 10-20, 20-30, 30-50, 50-80
   std::vector<std::vector<int>> ctrBins = {{0},{1},{2,3},{4,5},{6,7,8,9},{10,11,12,13,14,15}}; // centrality bins as collections of centrality intervals

   TFile* inFile = nullptr;
   TFile* outFile = nullptr;
   TTree *inTree[nCtrBins];
   TTree *outTree[ParamsSet::ndRselcs][nCtrBins][ParamsSet::nSigns];

   //these variables are for setting branch address (per entry) for the input trees
   float pt[nCtrBins];
   float eta[nCtrBins];
   float phi[nCtrBins];
   float dP_overP[nCtrBins];
   float d0[nCtrBins];
   float z0[nCtrBins];
   int charge[nCtrBins];
   int quality[nCtrBins];
   int event_num[nCtrBins];
   int centrality[nCtrBins];
   float FCal_Et[nCtrBins];

   //these vectors are for storing data from the input trees
   std::vector<float>   *muon_pt[nCtrBins];
   std::vector<float>   *muon_eta[nCtrBins];
   std::vector<float>   *muon_phi[nCtrBins];
   std::vector<float>   *muon_dP_overP[nCtrBins];
   std::vector<float>   *muon_d0[nCtrBins];
   std::vector<float>   *muon_z0[nCtrBins];
   std::vector<int>     *muon_charge[nCtrBins];
   std::vector<int>     *muon_quality[nCtrBins];
   std::vector<int>     *ev_num[nCtrBins];
   std::vector<int>     *ev_centrality[nCtrBins];
   std::vector<float>   *ev_FCal_Et[nCtrBins];

   //these Muon instances & variables are for setting branch address (per entry) for the input trees
   // MuonPair *mpair = nullptr;
   // MuonPair *mpair;
   MuonPair *mpair = new MuonPair();

   int n_ss_scr_pairs = 0;
   int n_op_scr_pairs = 0;


   void InitInput();
   void InitOutput();
   // void InitHists();
   void ReadData();
   void GenerateRandPair(int num_muon, int nctr, bool opsign_only);
   bool IsResonance();
   bool IsPhotoProduction();
   // void FillHistograms(unsigned int ndr, unsigned int nctr, unsigned int nsign);
   void ImplementOneScramPair(int num_muon, int nctr, bool opsign_only = false);


   // TH1D* h_pair_dP_overP[ParamsSet::ndRselcs][nCtrBins][ParamsSet::nSigns] = {nullptr};
   // TH1D* h_Dphi[ParamsSet::ndRselcs][nCtrBins][ParamsSet::nSigns] = {nullptr};
   // TH1D* h_Deta[ParamsSet::ndRselcs][nCtrBins][ParamsSet::nSigns] = {nullptr};
   // TH1D* h_DR[ParamsSet::ndRselcs][nCtrBins][ParamsSet::nSigns] = {nullptr};
   // TH1D* h_Minv[ParamsSet::ndRselcs][nCtrBins][ParamsSet::nSigns] = {nullptr};
   // TH1D* h_pt_lead[ParamsSet::ndRselcs][nCtrBins][ParamsSet::nSigns] = {nullptr};
   // TH1D* h_eta_lead[ParamsSet::ndRselcs][nCtrBins][ParamsSet::nSigns] = {nullptr};
   // TH1D* h_eta_sublead[ParamsSet::ndRselcs][nCtrBins][ParamsSet::nSigns] = {nullptr};
   // TH1D* h_eta_avg[ParamsSet::ndRselcs][nCtrBins][ParamsSet::nSigns] = {nullptr};
   // TH1D* h_pair_pt[ParamsSet::ndRselcs][nCtrBins][ParamsSet::nSigns] = {nullptr};
   // TH1D* h_pair_eta[ParamsSet::ndRselcs][nCtrBins][ParamsSet::nSigns] = {nullptr};
   // TH1D* h_pair_y[ParamsSet::ndRselcs][nCtrBins][ParamsSet::nSigns] = {nullptr};
   // TH2D* h_eta_phi[ParamsSet::ndRselcs][nCtrBins][ParamsSet::nSigns] = {nullptr};
   // TH2D* h_eta1_eta2[ParamsSet::ndRselcs][nCtrBins][ParamsSet::nSigns] = {nullptr};
   // TH2D* h_eta_avg_Deta[ParamsSet::ndRselcs][nCtrBins][ParamsSet::nSigns] = {nullptr};
   // TH2D* h_eta_avg_pair_eta[ParamsSet::ndRselcs][nCtrBins][ParamsSet::nSigns] = {nullptr};
   // TH2D* h_pt1_pt2[ParamsSet::ndRselcs][nCtrBins][ParamsSet::nSigns] = {nullptr};
   // TH2D* h_ptlead_pair_pt[ParamsSet::ndRselcs][nCtrBins][ParamsSet::nSigns] = {nullptr};

public:
   ScrambSampleGen();
   //void Loop();
   ~ScrambSampleGen(){}
   void Run();
};

ScrambSampleGen::ScrambSampleGen(){
}

void ScrambSampleGen::InitInput(){

   inFile = new TFile("single_muon_trees_coarse_ctr_bins.root","read");

   for (int jctr = 0; jctr < nCtrBins; jctr++){
      // TTree* tempTree = (TTree*) inFile->Get(Form("muon_tree%d",jctr+1));
      // inTree[jctr] = tempTree->CloneTree();
      inTree[jctr] = (TTree*) inFile->Get(Form("muon_tree_ctr%d",jctr+1));
      inTree[jctr]->SetBranchAddress("pt"           , &pt[jctr]);
      inTree[jctr]->SetBranchAddress("eta"          , &eta[jctr]);
      inTree[jctr]->SetBranchAddress("phi"          , &phi[jctr]);
      inTree[jctr]->SetBranchAddress("dP_overP"     , &dP_overP[jctr]);
      inTree[jctr]->SetBranchAddress("d0"           , &d0[jctr]);
      inTree[jctr]->SetBranchAddress("z0"           , &z0[jctr]);
      inTree[jctr]->SetBranchAddress("charge"       , &charge[jctr]);
      inTree[jctr]->SetBranchAddress("quality"      , &quality[jctr]);
      inTree[jctr]->SetBranchAddress("ev_num"            , &event_num[jctr]);
      inTree[jctr]->SetBranchAddress("ev_centrality"     , &centrality[jctr]);
      inTree[jctr]->SetBranchAddress("ev_FCal_Et"        , &FCal_Et[jctr]);
   }
}

void ScrambSampleGen::InitOutput(){
   outFile = new TFile("scrambled_muon_pairs_coarse_ctr_bins.root","recreate");

   for (unsigned int idr = 0; idr < ParamsSet::ndRselcs; idr++){
      for (unsigned int jctr = 0; jctr < nCtrBins; jctr++){
         for (unsigned int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
            outTree[idr][jctr][ksign] = new TTree(Form("scramb_muon_pair_tree_dr%u_ctr%u_sign%u",idr+1,jctr+1,ksign+1),Form("scramb_muon_pair_tree_dr%u_ctr%u_sign%u",idr+1,jctr+1,ksign+1));
            outTree[idr][jctr][ksign]->Branch("MuonPairObj",&mpair);
         }
      }
   }
}

//initialize the histograms
// void ScrambSampleGen::InitHists(){
   
//    for (unsigned int idr = 0; idr < ParamsSet::ndRselcs; idr++){
//       for (unsigned int jctr = 0; jctr < nCtrBins; jctr++){
//          for (unsigned int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
//             // h_FCal_Et[idr]=new TH1D(Form("h_FCal_Et%d",idr),";FCal #it{E}_{T} [TeV];1/N_{evt} dN/dFCal #it{E}_{T}",550,-0.5,5.0);
//             // h_deltaP_overP[idr][jctr] =new TH1D(Form("h_deltaP_overP%d_ctr%d",idr,jctr),";(p_{ID}-p_{MS}-#Delta p_{calo}) / p_{ID};1/N_{evt} dN/d" ,pms.deltaP_overP_nbins*2,-pms.deltaP_overP_thrsh,pms.deltaP_overP_thrsh);
//             h_pair_dP_overP[idr][jctr][ksign] = new TH1D(Form("h_pair_dP_overP_dr%d_ctr%d_sign%d",idr+1,jctr+1,ksign+1),";(p_{ID}-p_{MS}-#Delta p_{calo}) / p_{ID};1/N_{evt} dN/d" ,pms.deltaP_overP_nbins,0,pms.deltaP_overP_max);
//             h_Dphi[idr][jctr][ksign] = new TH1D(Form("h_Dphi_dr%d_ctr%d_sign%d",idr+1,jctr+1,ksign+1),";#Delta#phi;1/N_{evt} dN/d#Delta#phi" ,128 ,-pms.PI ,pms.PI );
//             h_Deta[idr][jctr][ksign] = new TH1D(Form("h_Deta_dr%d_ctr%d_sign%d",idr+1,jctr+1,ksign+1),";#Delta#eta;1/N_{evt} dN/d#Delta#eta" ,100,-2.4,2.4);
//             h_DR[idr][jctr][ksign] = new TH1D(Form("h_DR_dr%d_ctr%d_sign%d",idr+1,jctr+1,ksign+1),";#Delta R;1/N_{evt} dN/d#Delta R", pms.deltaR_nbins[idr], 0, pms.deltaR_thrsh[idr]);
//             h_Minv[idr][jctr][ksign] = new TH1D(Form("h_Minv_dr%d_ctr%d_sign%d",idr+1,jctr+1,ksign+1),";m_{#mu#mu} [GeV];1/N_{evt} dN/dm_{#mu#mu}" ,300,0   ,15 );
//             h_Minv[idr][jctr][ksign] = new TH1D(Form("h_Minv_dr%d_ctr%d_sign%d",idr+1,jctr+1,ksign+1),";m_{#mu#mu} [GeV];1/N_{evt} dN/dm_{#mu#mu}" ,pms.minv_nbins[idr],0   ,pms.minv_max[idr] );
//             h_pt_lead[idr][jctr][ksign] = new TH1D(Form("h_pt_lead_dr%d_ctr%d_sign%d",idr+1,jctr+1,ksign+1),";p_{T}^{lead} [GeV];1/N_{evt} dN/dp_{T}^{lead}" ,pms.npt_bins, pms.pTBins );
//             h_eta_lead[idr][jctr][ksign] = new TH1D(Form("h_pt_lead_dr%d_ctr%d_sign%d",idr+1,jctr+1,ksign+1),";p_{T}^{lead} [GeV];1/N_{evt} dN/dp_{T}^{lead}" ,pms.npt_bins, pms.pTBins );
//             h_eta_sublead[idr][jctr][ksign] = new TH1D(Form("h_pt_lead_dr%d_ctr%d_sign%d",idr+1,jctr+1,ksign+1),";p_{T}^{lead} [GeV];1/N_{evt} dN/dp_{T}^{lead}" ,pms.npt_bins, pms.pTBins );
//             h_eta_avg[idr][jctr][ksign] = new TH1D(Form("h_eta_avg_dr%d_ctr%d_sign%d",idr+1,jctr+1,ksign+1),";#bar{#eta};1/N_{evt} dN/d#bar{#eta}",100,-2.4,2.4);
//             h_pair_pt[idr][jctr][ksign] = new TH1D(Form("h_pair_pt_dr%d_ctr%d_sign%d",idr+1,jctr+1,ksign+1),";p_{T}^{pair} [GeV];1/N_{evt} dN/dp_{T}^{pair}" ,pms.npairPT_bins, pms.pairPTBins );
//             h_pair_eta[idr][jctr][ksign] = new TH1D(Form("h_pair_eta_dr%d_ctr%d_sign%d",idr+1,jctr+1,ksign+1),";#eta_{pair};1/N_{evt} dN/d#eta_{pair}" ,100,-2.4,2.4);
//             h_pair_y[idr][jctr][ksign] = new TH1D(Form("h_pair_y_dr%d_ctr%d_sign%d",idr+1,jctr+1,ksign+1),";y_{pair};1/N_{evt} dN/dy_{pair}" ,90,-3,3);
//             h_eta_phi[idr][jctr][ksign] = new TH2D(Form("h_eta_phi_dr%d_ctr%d_sign%d",idr+1,jctr+1,ksign+1),";#Delta#phi;#bar{#eta}", 128,-pms.PI,pms.PI,100,-2.4,2.4);
//             h_eta1_eta2[idr][jctr][ksign] = new TH2D(Form("h_eta1_eta2_dr%d_ctr%d_sign%d",idr+1,jctr+1,ksign+1),";#eta_{sublead};#eta_{lead}",100,-2.4,2.4, 100,-2.4,2.4);
//             h_eta_avg_Deta[idr][jctr][ksign] = new TH2D(Form("h_eta_avg_Deta_dr%d_ctr%d_sign%d",idr+1,jctr+1,ksign+1),";#Delta#eta;#bar{#eta}",100,-2.4,2.4, 100,-2.4,2.4);
//             h_eta_avg_pair_eta[idr][jctr][ksign] = new TH2D(Form("h_eta_avg_pair_eta_dr%d_ctr%d_sign%d",idr+1,jctr+1,ksign+1),";#eta_{pair};#bar{#eta}",100,-2.4,2.4, 100,-2.4,2.4);
//             h_pt1_pt2[idr][jctr][ksign] = new TH2D(Form("h_pt1_pt2_dr%d_ctr%d_sign%d",idr+1,jctr+1,ksign+1),";p_{T}^{sublead} [GeV];p_{T}^{lead} [GeV]",pms.npt_bins, pms.pTBins, pms.npt_bins, pms.pTBins);
//             h_ptlead_pair_pt[idr][jctr][ksign] = new TH2D(Form("h_ptlead_pair_pt_dr%d_ctr%d_sign%d",idr+1,jctr+1,ksign+1),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",pms.npairPT_bins, pms.pairPTBins, pms.npt_bins, pms.pTBins);
//          }
//       }
//    }
// }

#endif

