#pragma once

#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <string.h>
#include <fstream>
#include <map>
#include <memory>
#include <cmath>  
#include "vector"
#include "TH1D.h"
#include "TH2D.h"
#include "../MuonObjectsParamsAndHelpers/MuonPairBase.h"
#include "../MuonObjectsParamsAndHelpers/ParamsSet.h"


template <class PairT, class MuonT, class Derived>
class DimuonAlgCoreT{
public:
  	using pair_t = PairT;
  	using muon_t = MuonT;

protected:
    Derived& self() { return static_cast<Derived&>(*this); }
    const Derived& self() const { return static_cast<const Derived&>(*this); }

// --------------------- input parameter ---------------------------
    ParamsSet pms;

// --------------------- input files & trees & data for setting branches ---------------------------
    TChain          *fChain;   //!pointer to the analyzed TTree or TChain

// --------------------- temporary muon and muonpair objects ---------------------------
    Long64_t nentries;

    muon_t* muon_raw_ptr{nullptr};

    using PairPtr = std::shared_ptr<pair_t>;
	 
	PairPtr mpair;
  	pair_t* mpair_raw_ptr{nullptr};

    std::vector<PairPtr> muon_pair_list_cur_event_pre_resonance_cut;

	std::vector<int> resonance_tagged_muon_index_list_v2 {};
	std::vector<int> resonance_tagged_muon_index_list {};

// --------------------- output file, histograms & trees ---------------------------
    TFile *m_outfile = nullptr;
    TFile *m_outHistFile = nullptr;
    TTree* muonOutTree;
    TTree* muonPairOutTree[ParamsSet::nSigns];

    std::string output_file_path;
    std::string output_hist_file_path;

    int numCuts = 0;
	std::vector<std::string> cutLabels = {};
	TH1D* h_cutAcceptance[ParamsSet::nSigns];

// --------------------- protected class methods ---------------------------

	void Initialize();
  	
  	void InitOutput();
	void InitOutputTrees();    
    void InitOutputHists();

    void FillSingleMuonTree();
	void FillMuonPairTree();

	bool PassCuts(){return self().PassCutsHook();} // required

    void FillMuonPair(int pair_ind){self().FillMuonPairHook(pair_ind);}
	void ResonanceTagging();
	void ResonanceTaggingV2();
	bool IsPhotoProduction();
	
	void HistAdjust();

    void Finalize();

    void InitializeExtraHook() {}
    void InitParamsHook(){}
    void InitInputHook() {}
    void InitTempVariablesHook() {}
    void InitOutputSettingsHook() {}

    void OutputTreePath() {return self().OutputTreePathHook();} // required
    void OutputHistPath() {return self().OutputHistPathHook();} // required
    void InitOutputTreesExtraHook(){}
    void InitOutputHistsExtraHook(){}

	void ProcessDataHook(){}
	void FillMuonPairTreeHook(){}
    void FillSingleMuonTreeHook(){}

    void HistAdjustHook() {}
    void FinalizeHook() {}


public:
    bool debug_mode = false;
	bool output_single_muon_tree = false;

// --------------------- public class methods ---------------------------
    void Run();
};


template <class Derived>
struct PowhegHooksDefault {
};
