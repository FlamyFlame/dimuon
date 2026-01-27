#pragma once

#include "../MuonObjectsParamsAndHelpers/MuonPairPbPb.h"
#include "DimuonDataAnalysisBaseClass.c"

class PbPbExtras{
    // Read through the N-tuple, apply appropriate cuts, including killing pairs from resonances or photoproduction
    // output muon-pair or single-muon information into a TTree
protected:
    // --------------------- general settings ---------------------------

    // const ParamsSet pms;
    // int scaleFactorCtrs[ParamsSet::nCtrIntvls] = {1,1,3,5,5};
    // int scaleFactorPts[ParamsSet::nPtBins] = {1,1,1,1,1};

      // --------------------- input files & trees & data for setting branches---------------------------

    // Declaration of additional leaf types
    Float_t         FCal_Et;
    Int_t           centrality;

    // --------------------- class methods ---------------------------

    void InitParamsExtra();
    void InitInputExtra();
    void PerformTChainFill();
    void InitOutputSettingsExtra();

    void FillSingleMuonTreeExtra();
    void FillMuonPairTreeExtra();

    void FillMuonPairExtra(int pair_ind);
    // --------------------- output file, histograms & trees ---------------------------

    TTree* muonOutTreeBinned[ParamsSet::nCtrIntvls];
    TTree* muonPairOutTreeBinned[ParamsSet::nCtrIntvls][ParamsSet::nSigns];


public :
    bool turn_on_ctr_binned_tree_writing = false;
    
    explicit PbPbExtras(int run_year_input, int file_batch_input);
    ~PbPbExtras(){}
};

