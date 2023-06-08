#ifndef ScrambSampleGen_h
#define ScrambSampleGen_h
// #define nSigns 2
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

   ParamsSet pms;
   static constexpr float ctr_step = 5.;
   static const unsigned int nctr_intvls = 16;
   int nmodes = 3;
   // 0-5, 5-10, 10-20, 20-30, 30-50, 50-80
   // std::vector<std::vector<int>> ctrBins = {{0},{1},{2,3},{4,5},{6,7,8,9},{10,11,12,13,14,15}}; // centrality bins as collections of centrality intervals


   // std::vector<int> nScramb_op = {2216608, 1800972, 1442700, 1137412, 880732, 668284, 495680, 358545, 252502, 173171, 116172, 77127, 51656, 35380, 23920, 12898};
   // std::vector<int> nScrambLeftSpill_op = {483028, 388815, 308215, 240132, 183470, 137135, 100033, 71068, 49146, 33172, 22052, 14690, 9991, 6862, 4206};
   // std::vector<int> nScrambRightSpill_op = {517644, 418652, 333638, 261503, 201156, 151500, 111442, 79887, 55739, 37904, 25287, 16794, 11329, 7798, 5107};
   // std::vector<int> nScramb_ss= {2175940, 1740376, 1368424, 1054992, 794980, 583300, 414856, 284559, 187316, 118032, 71615, 42972, 27012, 18640, 12764, 4292};
   // std::vector<int> nScrambLeftSpill_ss= {469336, 371204, 288128, 218833, 162044, 116490, 80896, 53990, 34498, 21147, 12664, 7776, 5210, 3692, 1950};
   // std::vector<int> nScrambRightSpill_ss= {505608, 402180, 314227, 240479, 179663, 130505, 91732, 62072, 40250, 24994, 15030, 9085, 5887, 4161, 2635};
   
   // // the 3 lines are: within one centrality bin, "left spill" and "right spill"
   std::vector<std::vector<int>> nScramb_op = 
      {{2216608, 1800972, 1442700, 1137412, 880732, 668284, 495680, 358545, 252502, 173171, 116172, 77127, 51656, 35380, 23920, 12898},
       {0, 483028, 388815, 308215, 240132, 183470, 137135, 100033, 71068, 49146, 33172, 22052, 14690, 9991, 6862, 4206},
       {517644, 418652, 333638, 261503, 201156, 151500, 111442, 79887, 55739, 37904, 25287, 16794, 11329, 7798, 5107, 0}};
   std::vector<std::vector<int>> nScramb_ss= 
      {{2175940, 1740376, 1368424, 1054992, 794980, 583300, 414856, 284559, 187316, 118032, 71615, 42972, 27012, 18640, 12764, 4292},
       {0, 469336, 371204, 288128, 218833, 162044, 116490, 80896, 53990, 34498, 21147, 12664, 7776, 5210, 3692, 1950},
       {505608, 402180, 314227, 240479, 179663, 130505, 91732, 62072, 40250, 24994, 15030, 9085, 5887, 4161, 2635, 0}};

   // std::vector<std::vector<int>> nScramb_op = 
   //    {{20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20},
   //     {15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15},
   //     {15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15}};
   // std::vector<std::vector<int>> nScramb_ss= 
   //    {{15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15},
   //     {10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10},
   //     {10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10}};

   TFile* inFile = nullptr;
   TFile* outFile = nullptr;
   TTree *inTree[nctr_intvls];
   TTree *outTree[ParamsSet::ndRselcs][nctr_intvls][ParamsSet::nSigns];

   //these variables are for setting branch address (per entry) for the input trees
   float pt[nctr_intvls];
   float eta[nctr_intvls];
   float phi[nctr_intvls];
   float dP_overP[nctr_intvls];
   float d0[nctr_intvls];
   float z0[nctr_intvls];
   int charge[nctr_intvls];
   int quality[nctr_intvls];
   int event_num[nctr_intvls];
   int centrality[nctr_intvls];
   float FCal_Et[nctr_intvls];

   //these vectors are for storing data from the input trees
   std::vector<float>   *muon_pt[nctr_intvls];
   std::vector<float>   *muon_eta[nctr_intvls];
   std::vector<float>   *muon_phi[nctr_intvls];
   std::vector<float>   *muon_dP_overP[nctr_intvls];
   std::vector<float>   *muon_d0[nctr_intvls];
   std::vector<float>   *muon_z0[nctr_intvls];
   std::vector<int>     *muon_charge[nctr_intvls];
   std::vector<int>     *muon_quality[nctr_intvls];
   std::vector<int>     *ev_num[nctr_intvls];
   std::vector<int>     *ev_centrality[nctr_intvls];
   std::vector<float>   *ev_FCal_Et[nctr_intvls];

   //these Muon instances & variables are for setting branch address (per entry) for the input trees
   // MuonPair *mpair = nullptr;
   // MuonPair *mpair;
   MuonPair *mpair = new MuonPair();

   int n_ss_scr_pairs = 0;
   int n_op_scr_pairs = 0;


   void InitInput();
   void InitOutput();
   // void InitHists();
   void InitTempVariables();
   void GenerateRandPair(int nctr, int mode, bool opsign_only);
   bool ResonanceCut();
   bool PhotoProductionCut();
   // void FillHistograms(unsigned int ndr, unsigned int nctr, unsigned int nsign);
   void ImplementOneScramPair(int nctr, int mode, bool opsign_only = false);


public:
   bool keep_spill = true;
   ScrambSampleGen();
   //void Loop();
   ~ScrambSampleGen(){}
   void Run();
};

ScrambSampleGen::ScrambSampleGen(){
   if (nScramb_op[0].size() != nctr_intvls || nScramb_op[1].size() != nctr_intvls || nScramb_op[2].size() != nctr_intvls || nScramb_ss[0].size() != nctr_intvls || nScramb_ss[1].size() != nctr_intvls || nScramb_ss[2].size() != nctr_intvls){
      std::cout << "The numbers of elements in the scrambled-pair-count vectors do not match number of centrality bins." << std::endl;
      throw std::exception();
   }
   if (nScramb_op[1][0] != 0 || nScramb_op[2][nctr_intvls-1] != 0 || nScramb_ss[1][0] != 0 || nScramb_ss[2][nctr_intvls-1] != 0){
      std::cout << "The first element in the left-spill mode and last element in the right-spill mode in nScramb must be zero." << std::endl;
      throw std::exception();
   } 
}

void ScrambSampleGen::InitInput(){

   inFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/single_muon_trees_small_ctr_intvls.root","read");

   for (int jctr = 0; jctr < nctr_intvls; jctr++){
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
      inTree[jctr]->SetBranchAddress("ev_num"            , &event_num[jctr]); // change the names of these values for reading so that they don't coincide with the vector names
      inTree[jctr]->SetBranchAddress("ev_centrality"     , &centrality[jctr]);
      inTree[jctr]->SetBranchAddress("ev_FCal_Et"        , &FCal_Et[jctr]);
   }
}

void ScrambSampleGen::InitOutput(){
   if (keep_spill){
      outFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/scrambled_muon_pairs_small_ctr_intvls.root","recreate");
   }else{
      outFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/scrambled_muon_pairs_small_ctr_intvls_nospill.root","recreate");   
   }
   

   for (unsigned int idr = 0; idr < ParamsSet::ndRselcs; idr++){
      for (unsigned int jctr = 0; jctr < nctr_intvls; jctr++){
         for (unsigned int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
            outTree[idr][jctr][ksign] = new TTree(Form("scramb_muon_pair_tree_dr%u_ctr%u_sign%u",idr+1,jctr+1,ksign+1),Form("scramb_muon_pair_tree_dr%u_ctr%u_sign%u",idr+1,jctr+1,ksign+1));
            outTree[idr][jctr][ksign]->Branch("MuonPairObj",&mpair);
         }
      }
   }
}

void ScrambSampleGen::InitTempVariables(){
   for (int jctr = 0; jctr < nctr_intvls; jctr++){
      muon_pt[jctr] =  new std::vector<float>();
      muon_eta[jctr] =  new std::vector<float>();
      muon_phi[jctr] =  new std::vector<float>();
      muon_dP_overP[jctr] =  new std::vector<float>();
      muon_d0[jctr] =  new std::vector<float>();
      muon_z0[jctr] =  new std::vector<float>();
      muon_charge[jctr] =  new std::vector<int>();
      muon_quality[jctr] =  new std::vector<int>();
      ev_num[jctr] =  new std::vector<int>();
      ev_centrality[jctr] =  new std::vector<int>();
      ev_FCal_Et[jctr] =  new std::vector<float>();
    }
}

   
#endif
