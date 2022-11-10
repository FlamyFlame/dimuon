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

  int mode = 4;
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
   bool IsResonance();
   bool IsPhotoProduction();
   void FillSingleMuonTree();
   void FillMuonPairTree();
   void WriteOutput();


  // --------------------- output file, histograms & trees ---------------------------

   TFile *m_outfile = nullptr;

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
    if(mode != 3 && mode != 4){
        std::cout<<"Error:: Mode has to equal 3 or 4; code is used for outputting muon / muon-pair trees only."<<std::endl;
        throw std::exception();
    }
}

//initialize the TChain
void MuonNTupleFirstPass::Init(){

   fChain = new TChain("HeavyIonD3PD","HeavyIonD3PD");
   fChain->SetMakeClass(1);
   fChain->Add("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/All_Data2018_12March2022.root");
   fChain->Add("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/All_Data2015_12March2022.root"); 

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


void MuonNTupleFirstPass::InitOutput(){

  if (mode == 3){
    m_outfile=new TFile("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/single_muon_trees_coarse_ctr_bins.root","recreate");
    muonOutTree = new TTree("muon_tree","all single muons");
    muonOutTree->Branch("MuonObj",&tempmuon);

    for (int jctr = 0; jctr < ParamsSet::nCtrBins; jctr++){
      muonOutTreeBinned[jctr] = new TTree(Form("muon_tree_ctr%d",jctr+1),Form("all muons, centrality bin %d",jctr+1));
      muonOutTreeBinned[jctr]->Branch("MuonObj",&tempmuon);
    }
  }
  else{ // mode == 4
    m_outfile=new TFile("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/muon_pairs_coarse_ctr_bins.root","recreate");

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


