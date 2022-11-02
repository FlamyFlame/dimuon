#ifndef ScrambSampleGenPP_h
#define ScrambSampleGenPP_h

#include "MuonPair.h"
#include "ParamsSet.h"
#include <TROOT.h>
//#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
//#include <string.h>
#include "vector"
#include "TH1D.h"
#include "TH2D.h"
#include <stdlib.h>

class ScrambSampleGenPP{
private:

   int nScramb[ParamsSet::nSigns] = {8213832,21092256};

   ParamsSet pms;

   TFile* inFile = nullptr;
   TFile* outFile = nullptr;
   TTree *inTree;
   TTree *outTree[ParamsSet::ndRselcs][ParamsSet::nSigns];

   //these variables are for setting branch address (per entry) for the input trees
   float pt;
   float eta;
   float phi;
   float dP_overP;
   float d0;
   float z0;
   int charge;
   int quality;
   int event_num;

   //these vectors are for storing data from the input trees
   std::vector<float>   *muon_pt;
   std::vector<float>   *muon_eta;
   std::vector<float>   *muon_phi;
   std::vector<float>   *muon_dP_overP;
   std::vector<float>   *muon_d0;
   std::vector<float>   *muon_z0;
   std::vector<int>     *muon_charge;
   std::vector<int>     *muon_quality;
   std::vector<int>     *ev_num;

   //for setting branch address (per entry) for the output trees
   MuonPair *mpair = new MuonPair();

   int n_ss_scr_pairs = 0;
   int n_op_scr_pairs = 0;


   void InitInput();
   void InitOutput();
   void ReadData();
   void GenerateRandPair(int num_muon, bool opsign_only);
   bool CheckResonance();
   void ImplementOneScramPair(int num_muon, bool opsign_only = false);

public:
   ScrambSampleGenPP();
   ~ScrambSampleGenPP(){}
   void Run();
};

ScrambSampleGenPP::ScrambSampleGenPP(){
}

void ScrambSampleGenPP::InitInput(){

   inFile = new TFile("single_muon_trees_pp.root","read");
   inTree = (TTree*) inFile->Get("muon_tree");
   inTree->SetBranchAddress("pt"           , &pt);
   inTree->SetBranchAddress("eta"          , &eta);
   inTree->SetBranchAddress("phi"          , &phi);
   inTree->SetBranchAddress("dP_overP"     , &dP_overP);
   inTree->SetBranchAddress("d0"           , &d0);
   inTree->SetBranchAddress("z0"           , &z0);
   inTree->SetBranchAddress("charge"       , &charge);
   inTree->SetBranchAddress("quality"      , &quality);
   inTree->SetBranchAddress("ev_num"            , &event_num);
}

void ScrambSampleGenPP::InitOutput(){
   outFile = new TFile("scrambled_muon_pairs_pp.root","recreate");

   for (unsigned int idr = 0; idr < ParamsSet::ndRselcs; idr++){
      for (unsigned int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
         outTree[idr][ksign] = new TTree(Form("scramb_muon_pair_tree_dr%u_sign%u",idr+1,ksign+1),Form("scramb_muon_pair_tree_dr%u_sign%u",idr+1,ksign+1));
         outTree[idr][ksign]->Branch("MuonPairObj",&mpair);
      }
   }
}

#endif

