#pragma once

#include "../MuonObjectsParamsAndHelpers/MuonPairPbPb.h"
#include "TTree.h"

template <class PairT, class MuonT, class Derived, class... Extras>
class DimuonDataAlgCoreT;

template <class Derived>
class PbPbExtras{
    template <class, class, class, class...> friend class DimuonDataAlgCoreT;

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

    Derived& self() { return static_cast<Derived&>(*this); }
    const Derived& self() const { return static_cast<const Derived&>(*this); }

    ~PbPbExtras(){}
};

