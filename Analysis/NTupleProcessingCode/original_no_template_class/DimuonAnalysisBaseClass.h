#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <string.h>
#include <fstream>
#include <map>
#include <cmath>  
#include "vector"
#include "TH1D.h"
#include "TH2D.h"
#include "../MuonObjectsParamsAndHelpers/MuonPair.h"
#include "../MuonObjectsParamsAndHelpers/ParamsSet.h"


#ifndef __DimuonAnalysisBaseClass_h__
#define __DimuonAnalysisBaseClass_h__


class DimuonAnalysisBaseClass{


protected:

// --------------------- input parameter ---------------------------
    ParamsSet pms;

// --------------------- input files & trees & data for setting branches ---------------------------
    TChain          *fChain;   //!pointer to the analyzed TTree or TChain
// --------------------- temporary muon and muonpair objects ---------------------------
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

	virtual void InitInput(){}
  	virtual void InitTempVariables(){}
  	virtual void InitOutput();
	
	virtual void ProcessData(){}

    virtual void FillSingleMuonTree();
	virtual void FillMuonPairTree() = 0;

	virtual bool PassCuts(std::shared_ptr<MuonPair> const& mpair) = 0;
    // virtual void FillMuonPair(int pair_ind, std::shared_ptr<MuonPair> const& mpair){}
	virtual void ResonanceTagging(std::shared_ptr<MuonPair> const& mpair);
	virtual void ResonanceTaggingV2(std::shared_ptr<MuonPair> const& mpair);
	virtual bool IsPhotoProduction(std::shared_ptr<MuonPair> const& mpair);
	
	virtual void HistAdjust();

public:
    bool debug_mode = false;

// --------------------- public class methods ---------------------------
	~DimuonAnalysisBaseClass(){}
	virtual void Run(){}
};

#endif