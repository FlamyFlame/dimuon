#pragma once

#include <unordered_map>
#include <utility>
#include "../MuonObjectsParamsAndHelpers/MuonPairPbPb.h"
#include "TTree.h"
#include "TGraph.h"
#include "PbPbEventSelConfig.h"

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
    Float_t         FCal_Et_P, FCal_Et_N;               // side A, side C
    Int_t           centrality;
    // Run 3 only
    Float_t         zdc_ZdcEnergy[2]{};                 // [0]=A, [1]=C
    Float_t         zdc_ZdcTime[2]{};                   // [0]=C, [1]=A
    Float_t         zdc_ZdcModulePreSampleAmp[2][4]{};  // [side][module], [0]=C, [1]=A
    std::vector<int>* trk_numqual{nullptr};              // 8-element; indices: [2]=HIloose+pT, [3]=HItight+pT, [6]=HIloose, [7]=HItight

    // --------------------- class methods ---------------------------

    void InitParamsExtra();
    void InitInputExtra();
    void PerformTChainFill();
    void InitOutputSettingsExtra();

    void FillSingleMuonTreeExtra();
    void FillMuonPairTreeExtra();

    void FillMuonPairExtra(int pair_ind);

    // --------------------- event-level selection (Run 3 PbPb only) ---------------
    void InitEventSel();
    bool PassEventSel() const;
    bool PassEventSelExtra() { return PassEventSel(); }

    // Loaded by InitEventSel() from event_sel_cuts_pbpb_20YY.root
    TGraph* g_evsel_cut1_{nullptr};       // ZDC-FCal banana upper bound
    double  evsel_cut2_ns_{1.8};          // ZDC time window [ns]
    float   evsel_cut3_A_{385.f};         // ZDC preamp A upper cut [ADC] (scalar fallback / 23+24)
    float   evsel_cut3_C_{385.f};         // ZDC preamp C upper cut [ADC] (scalar fallback / 23+24)
    // PbPb25 only: per-run mu+7sigma preamp cuts, keyed by run number
    std::unordered_map<int, std::pair<float,float>> evsel_cut3_per_run_;
    // Runs excluded from all event selection and analysis (loaded by InitEventSel)
    std::unordered_set<int> evsel_bad_runs_;
    TGraph* g_evsel_cut4_{nullptr};       // nTrk frac lower bound
    TGraph* g_evsel_cut5_lo_{nullptr};    // nTrk-FCal lower bound
    TGraph* g_evsel_cut5_hi_{nullptr};    // nTrk-FCal upper bound

    // --------------------- output file, histograms & trees ---------------------------

    TTree* muonOutTreeBinned[ParamsSet::nCtrIntvls];
    TTree* muonPairOutTreeBinned[ParamsSet::nCtrIntvls][ParamsSet::nSigns];


public :
    bool turn_on_ctr_binned_tree_writing = false;

    Derived& self() { return static_cast<Derived&>(*this); }
    const Derived& self() const { return static_cast<const Derived&>(*this); }

    ~PbPbExtras(){}
};

