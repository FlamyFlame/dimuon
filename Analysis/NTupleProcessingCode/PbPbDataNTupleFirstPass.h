#pragma once

#include "../MuonObjectsParamsAndHelpers/MuonPairPbPb.h"
#include "DimuonDataAnalysisBaseClass.c"

class PbPbDataNTupleFirstPass : public DimuonDataAnalysisBaseClass{
    // Read through the N-tuple, apply appropriate cuts, including killing pairs from resonances or photoproduction
    // output muon-pair or single-muon information into a TTree
private:
    // --------------------- general settings ---------------------------

    // const ParamsSet pms;
    // int scaleFactorCtrs[ParamsSet::nCtrIntvls] = {1,1,3,5,5};
    // int scaleFactorPts[ParamsSet::nPtBins] = {1,1,1,1,1};

      // --------------------- input files & trees & data for setting branches---------------------------

    // Declaration of additional leaf types
    Float_t         FCal_Et;
    Int_t           centrality;

    // --------------------- class methods ---------------------------

    virtual PairPtr MakeMuonPair() const override {
        return std::make_shared<MuonPairPbPb>();
    }

    void ParamCheck() override;
    void InitInput() override;
    void TChainFill() override;
    void InitOutput() override;

    void FillSingleMuonTree() override;
    void FillMuonPairTree() override;

    void FillMuonPair    (int pair_ind, std::shared_ptr<MuonPairData> const& mpair) override;
    void FillMuonPairPbPb(int pair_ind, std::shared_ptr<MuonPairPbPb> const& mpair);
    // --------------------- output file, histograms & trees ---------------------------

    TTree* muonOutTreeBinned[ParamsSet::nCtrIntvls];
    TTree* muonPairOutTreeBinned[ParamsSet::nCtrIntvls][ParamsSet::nSigns];


public :
    bool turn_on_ctr_binned_tree_writing = false;
    
    PbPbDataNTupleFirstPass();
    ~PbPbDataNTupleFirstPass(){}
    void Run() override;
};

