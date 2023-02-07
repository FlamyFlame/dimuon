#ifndef __SimpleOutput_h__
#define __SimpleOutput_h__

#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <string.h>
#include "vector"
// #include <fstream>
// #include <cmath>  
// #include "MuonPairMC.h"
// #include "ParamsSet.h"
// #include "TH1D.h"
// #include "TH2D.h"

class SimpleOutput{

private:
// --------------------- general settings ---------------------------

    float filter_effcy_bb = 0.0003774031;
    float filter_effcy_cc = 0.000005964574;

    // ParamsSet pms;

    std::string mcdir = "/usatlas/u/yuhanguo/usatlasdata/athena/runMCV2/";
    
// --------------------- input files & trees & data for setting branches---------------------------
  
    TChain          *fChain;   //!pointer to the analyzed TTree or TChain
  
    // Declaration of leaf types
    // UInt_t          RunNumber;

    std::vector<float>   *EventWeights           =nullptr;
    std::vector<int>* truth_id = nullptr;
    std::vector<int>* truth_barcode = nullptr;
    std::vector<float>* truth_pt = nullptr;
    std::vector<float>* truth_eta = nullptr;
    std::vector<std::vector<int>>* truth_parents = nullptr;
    std::vector<std::vector<int>>* truth_children = nullptr;
      


// --------------------- events of interest ---------------------------

    // bb: 3 from HF-diquark, 1 from 1103, 1 from 203?
    std::vector<int> bb_events_of_interest = {2853, 4414, 33844, 72637, 96789};
    // cc: 4 from [21 21], 3 from [3], 2 from [21 2]
    // std::vector<int> cc_events_of_interest = {164, 3198, 9872, 15617, 23391, 32701, 48796, 24135, 49037};
    std::vector<int> cc_events_of_interest = {25279, 27004, 36467, 40757, 14995};
    std::vector<int> events_of_interest;

// --------------------- output file, histograms & trees ---------------------------
  
    // ofstream outfile;
    FILE * outfile;

// --------------------- class methods ---------------------------
  
    void InitInput();
    void ProcessData();  

public :
    bool isTruth = true;
    std::string mc_mode;
    // std::string mc_mode = "mc_truth_cc";
    SimpleOutput();
    ~SimpleOutput(){}
    void Run();
    float filter_effcy;
};


SimpleOutput::SimpleOutput(){
}

//initialize the TChain
void SimpleOutput::InitInput(){

    fChain = new TChain("HeavyIonD3PD","HeavyIonD3PD");
    fChain->SetMakeClass(1);
    fChain->Add(Form("%s%s.root", mcdir.c_str(), mc_mode.c_str()));
    
    // MC muons have no quality, d0 or z0 recorded
    fChain->SetBranchAddress("truth_id"                   , &truth_id);
    fChain->SetBranchAddress("truth_barcode"              , &truth_barcode);
    fChain->SetBranchAddress("truth_parents"              , &truth_parents);
    fChain->SetBranchAddress("truth_children"             , &truth_children);
    fChain->SetBranchAddress("truth_pt"                   , &truth_pt);
    fChain->SetBranchAddress("truth_eta"                  , &truth_eta);
    fChain->SetBranchAddress("EventWeights"               , &EventWeights);

    //SetBranch Status
    fChain->SetBranchStatus("*"                         ,0);//switch off all branches, then enable just the ones that we need
    fChain->SetBranchStatus("truth_id"                  ,1);
    fChain->SetBranchStatus("truth_barcode"             ,1);
    fChain->SetBranchStatus("truth_parents"             ,1);
    fChain->SetBranchStatus("truth_children"            ,1);
    fChain->SetBranchStatus("truth_pt"                  ,1);
    fChain->SetBranchStatus("truth_eta"                 ,1);
    fChain->SetBranchStatus("EventWeights"              ,1);
}

    
#endif
