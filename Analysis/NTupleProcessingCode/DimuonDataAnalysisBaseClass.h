#pragma once

#include "../MuonObjectsParamsAndHelpers/MuonPairData.h"
#include "../MuonObjectsParamsAndHelpers/muon_pair_enums_data.h"
#include "DimuonAnalysisBaseClass.c"


class DimuonDataAnalysisBaseClass : public DimuonAnalysisBaseClass{
protected:
    DimuonDataAnalysisBaseClass(int run_year_input, int file_batch_input)
        : run_year (run_year_input), file_batch (file_batch_input){}

// --------------------- general settings ---------------------------
    std::string trig_suffix = "";
    ParamsSet pms;
    std::string data_dir;

    int run_year; // used (only) by PbPb
    int file_batch;
    bool isRun3;

    bool isPbPb; // true for PbPb, false for pp

    bool trigger_effcy_calc; // if true, not care about physics origin, e.g, skip photoproduction + resonance cuts to gain more statistics
    bool use_mu6_for_trg_eff; // if true (now only for Pb+Pb 23), use mu6 as a supporting trigger for P[dimuon trg | supporting] efficiency analysis
    bool use_mu8_for_trg_eff; // if true (now only for Pb+Pb 23), use mu6 as a supporting trigger for P[dimuon trg | supporting] efficiency analysis
// --------------------- input files & trees & data for setting branches---------------------------

    TChain          *fChain;   //!pointer to the analyzed TTree or TChain

    // Declaration of leaf types
    UInt_t          RunNumber;
    UInt_t          lbn;
    UInt_t          bcid;

    Bool_t               b_HLT_mu4; // single muon 4GeV trigger
    std::vector<bool>    *muon_b_HLT_mu4     =nullptr; // single muon 4GeV trigger

    Bool_t               b_HLT_mu6_L1MU3V; // single muon 6GeV trigger (same L1 trigger as mu4)
    std::vector<bool>    *muon_b_HLT_mu6_L1MU3V     =nullptr; // single muon 6GeV trigger (same L1 trigger as mu4)

    Bool_t               b_HLT_mu8_L1MU5VF; // single muon 6GeV trigger (same L1 trigger as mu4)
    std::vector<bool>    *muon_b_HLT_mu8_L1MU5VF     =nullptr; // single muon 6GeV trigger (same L1 trigger as mu4)

    Bool_t               b_HLT_2mu4;
    std::vector<bool>    *dimuon_b_HLT_2mu4    =nullptr;

    Bool_t               b_HLT_mu4_mu4noL1;
    std::vector<bool>    *dimuon_b_HLT_mu4_mu4noL1    =nullptr;

    std::vector<float>   *muon_deltaP_overP           =nullptr;

    std::vector<float>   *muon_pt          =nullptr; // only used for mu4 trigger-efficiency study
    std::vector<float>   *muon_eta         =nullptr; // only used for mu4 trigger-efficiency study
    std::vector<float>   *muon_phi         =nullptr; // only used for mu4 trigger-efficiency study
    std::vector<int>     *muon_quality     =nullptr; // only used for mu4 trigger-efficiency study
    std::vector<float>   *muon_d0          =nullptr; // only used for mu4 trigger-efficiency study
    std::vector<float>   *muon_z0          =nullptr; // only used for mu4 trigger-efficiency study
  
    std::vector<int>     *muon_pair_muon1_index       =nullptr;
    std::vector<float>   *muon_pair_muon1_pt          =nullptr;
    std::vector<float>   *muon_pair_muon1_eta         =nullptr;
    std::vector<float>   *muon_pair_muon1_phi         =nullptr;
    std::vector<int>     *muon_pair_muon1_quality     =nullptr;
    std::vector<float>   *muon_pair_muon1_d0          =nullptr;
    std::vector<float>   *muon_pair_muon1_z0          =nullptr;
    std::vector<float>   *muon_pair_muon1_trk_pt          =nullptr;
    std::vector<float>   *muon_pair_muon1_trk_eta         =nullptr;
    std::vector<float>   *muon_pair_muon1_trk_phi         =nullptr;
  
    std::vector<int>     *muon_pair_muon2_index       =nullptr;
    std::vector<float>   *muon_pair_muon2_pt          =nullptr;
    std::vector<float>   *muon_pair_muon2_eta         =nullptr;
    std::vector<float>   *muon_pair_muon2_phi         =nullptr;
    std::vector<int>     *muon_pair_muon2_quality     =nullptr;
    std::vector<float>   *muon_pair_muon2_d0          =nullptr;
    std::vector<float>   *muon_pair_muon2_z0          =nullptr;
    std::vector<float>   *muon_pair_muon2_trk_pt          =nullptr;
    std::vector<float>   *muon_pair_muon2_trk_eta         =nullptr;
    std::vector<float>   *muon_pair_muon2_trk_phi         =nullptr;

    // --------------------- temporary muon and muonpair objects ---------------------------

    using PairPtr = std::shared_ptr<MuonPairData>;

    // one copy, visible to everyone
    PairPtr mpair;
    MuonPairData* mpair_raw_ptr = nullptr;
    Muon* tempmuon = nullptr;
    std::vector<PairPtr> muon_pair_list_cur_event_pre_resonance_cut;

// --------------------- output file, histograms & trees ---------------------------

// --------------------- private class methods ---------------------------

    void TrigModeToSuffixMap();
    bool TrigMatch(int pair_ind, int m1_ind, int m2_ind);

    virtual void InitInput() override;
    virtual void InitInputBranchesDimuonAnalysis(); // set branches for muon-pair variables
    virtual void InitInputBranchesSingleMuonAnalysis(); // set branches for single-muon variables (for MB-data mu4-trigger-efficiency study)
    virtual void InitOutput() override;
    virtual void TChainFill() = 0; // purely virtual method

    virtual PairPtr MakeMuonPair() const {
        return std::make_shared<MuonPairData>();
    }

    virtual void FillMuonPair (int pair_ind, std::shared_ptr<MuonPairData> const& mpair);
    bool PassCuts(const std::shared_ptr<MuonPairData>& mpair, bool requireTight);
    bool PassCuts(const std::shared_ptr<MuonPair>& mpair) override;

    virtual void FillMuonPairTree() override;

    virtual void InitParams();

    virtual void PerformAdditionalPairAnalysis(){}
    void ProcessData() override;

public:
    int trigger_mode = 1; // default to single-muon trigger

    bool isMinBias = false; // whether use MinBias data for single-muon trigger efficiency
    bool output_single_muon_tree = false;

    bool requireTight = false;
    int resonance_cut_mode = 1; // only applies to nominal (not-trigger-efficiency) analysis

    bool pbpb24_mu4_NO_trig_calc = false; // only affects Pb+Pb24: set to true if for Pb+Pb 24, single mu4, do NOT want to perform trigger efficiency study

    bool turn_on_track_charge = false; // turn on if track charge is stored

    bool filter_out_photo_resn_for_trig_effcy = true; // if true: filter out pairs from photoproduction or resonance decay for trigger efficiency study as well
// --------------------- public class methods ---------------------------
	~DimuonDataAnalysisBaseClass(){}
    virtual void PrintInstructions();
    void Run() override;

};
