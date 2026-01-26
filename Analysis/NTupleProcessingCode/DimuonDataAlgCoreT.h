#pragma once

#include "../MuonObjectsParamsAndHelpers/MuonPairReco.h"
#include "../MuonObjectsParamsAndHelpers/muon_pair_enums_data.h"
#include "DimuonAlgCoreT.c"


template <class PairT, class MuonT, class Derived, class... Extras>
class DimuonDataAlgCoreT
    : public DimuonAlgCoreT<PairT, MuonT, Derived>
{
public:
    using pair_t = PairT;

protected:
    TChain*&                    fChainRef()            { return this->fChain; }
    std::shared_ptr<pair_t>&    mpairRef()             { return this->mpair; }
    TH1D* (&h_cutAcceptanceRef())[ParamsSet::nSigns] {
        return this->h_cutAcceptance;
    }

    DimuonDataAlgCoreT(int run_year_input, int file_batch_input)
        : run_year (run_year_input), file_batch (file_batch_input){
            PrintInstructions();
        }

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

// --------------------- output file, histograms & trees ---------------------------

// --------------------- private class methods ---------------------------

    void TrigModeToSuffixMap();
    bool TrigMatch(int pair_ind, int m1_ind, int m2_ind);

    void InitInputBranchesDimuonAnalysis(); // set branches for muon-pair variables
    void InitInputBranchesSingleMuonAnalysis(); // set branches for single-muon variables (for MB-data mu4-trigger-efficiency study)
    void TChainFill(){ return self().TChainFill(); }

    bool IsPhotoProduction();


    // --------------- InitInputImpl ---------------
    template <class E>
    void CallInitInput() {
        if constexpr (requires(E& e){ e.InitInputExtra(); }) {
            static_cast<E&>(*this).InitInputExtra();
        }
    }

    void InitInput_DataCore();

    // --------------- InitInputBranchesSingleMuonAnalysisImpl ---------------
    template <class E>
    void CallInitInputBranchesSingleMuonAnalysis() {
        if constexpr (requires(E& e){ e.InitInputBranchesSingleMuonAnalysisExtra(); }) {
            static_cast<E&>(*this).InitInputBranchesSingleMuonAnalysisExtra();
        }
    }

    void InitInputBranchesSingleMuonAnalysis_DataCore();

    // --------------- InitInputBranchesDimuonAnalysisImpl ---------------
    template <class E>
    void CallInitInputBranchesDimuonAnalysis() {
        if constexpr (requires(E& e){ e.InitInputBranchesDimuonAnalysisExtra(); }) {
            static_cast<E&>(*this).InitInputBranchesDimuonAnalysisExtra();
        }
    }

    void InitInputBranchesDimuonAnalysis_DataCore();

    // --------------- InitTempVariablesImpl ---------------
    template <class E>
    void CallInitTempVariables() {
        if constexpr (requires(E& e){ e.InitTempVariablesExtra(); }) {
            static_cast<E&>(*this).InitTempVariablesExtra();
        }
    }

    void InitTempVariables_DataCore(){}

    // --------------- InitParamsImpl ---------------
    template <class E>
    void CallInitParams() {
        if constexpr (requires(E& e){ e.InitParamsExtra(); }) {
            static_cast<E&>(*this).InitParamsExtra();
        }
    }

    void InitParams_DataCore();

    // --------------- InitOutputSettingsImpl ---------------
    template <class E>
    void CallInitOutputSettings() {
        if constexpr (requires(E& e){ e.InitOutputSettingsExtra(); }) {
            static_cast<E&>(*this).InitOutputSettingsExtra();
        }
    }

    void InitOutputSettings_DataCore(){}

    // --------------- InitOutputTreesExtraImpl ---------------
    template <class E>
    void CallInitOutputTreesExtra() {
        if constexpr (requires(E& e){ e.InitOutputTreesExtra(); }) {
            static_cast<E&>(*this).InitOutputTreesExtra();
        }
    }

    // --------------- InitOutputHistsExtraImpl ---------------
    template <class E>
    void CallInitOutputHistsExtra() {
        if constexpr (requires(E& e){ e.InitOutputHistsExtra(); }) {
            static_cast<E&>(*this).InitOutputHistsExtra();
        }
    }

    // --------------- InitializeExtraImpl ---------------
    template <class E>
    void CallInitializeExtra() {
        if constexpr (requires(E& e){ e.InitializeExtra(); }) {
            static_cast<E&>(*this).InitializeExtra();
        }
    }

    void InitializeExtra_DataCore(){}

    // --------------- PassCutsImpl ---------------
    template <class E>
    bool CallPassCuts() {
        if constexpr (requires(E& e){ e.PassCutsExtra(); }) {
            return static_cast<E&>(*this).PassCutsExtra();
        } else {
            return true;
        }
    }

    bool PassCuts_DataCore(bool requireTight);

    // --------------- FillMuonPairTreeImpl ---------------
    template <class E>
    void CallFillMuonPairTree() {
        if constexpr (requires(E& e){ e.FillMuonPairTreeExtra(); }) {
            static_cast<E&>(*this).FillMuonPairTreeExtra();
        }
    }

    // --------------- FillSingleMuonTreeImpl ---------------
    template <class E>
    void CallFillSingleMuonTree() {
        if constexpr (requires(E& e){ e.FillSingleMuonTreeExtra(); }) {
            static_cast<E&>(*this).FillSingleMuonTreeExtra();
        }
    }

    // --------------- FillMuonPairImpl ---------------
    template <class E>
    void CallFillMuonPair(int pair_ind) {
        if constexpr (requires(E& e, int ind) { e.FillMuonPairExtra(ind); }) {
            static_cast<E&>(*this).FillMuonPairExtra(pair_ind);
        }
    }

    void FillMuonPair_DataCore(int pair_ind);

    // --------------- HistAdjustImpl ---------------
    template <class E>
    void CallHistAdjust() {
        if constexpr (requires(E& e){ e.HistAdjustExtra(); }) {
            static_cast<E&>(*this).HistAdjustExtra();
        }
    }

    void HistAdjust_DataCore(){}

   // --------------- PrintInstructionsImpl ---------------
    template <class E>
    void CallPrintInstructions() {
        if constexpr (requires(E& e){ e.PrintInstructionsExtra(); }) {
            static_cast<E&>(*this).PrintInstructionsExtra();
        }
    }

    void PrintInstructions_DataCore();

    // --------------- PerformAdditionalPairAnalysisImpl ---------------
    template <class E>
    void CallAdditionalPairAnalysis() {
        if constexpr (requires(E& e){ e.AdditionalPairAnalysis(); }) {
            static_cast<E&>(*this).AdditionalPairAnalysis();
        }
    }

    // --------------- FinalizeImpl ---------------
    template <class E>
    void CallFinalize() {
        if constexpr (requires(E& e){ e.FinalizeExtra(); }) {
            static_cast<E&>(*this).FinalizeExtra();
        }
    }

    void Finalize_DataCore(){}

public:
    int trigger_mode = 1; // default to single-muon trigger

    bool isMinBias = false; // whether use MinBias data for single-muon trigger efficiency

    bool requireTight = false;
    int resonance_cut_mode = 1; // only applies to nominal (not-trigger-efficiency) analysis

    bool pbpb24_mu4_NO_trig_calc = false; // only affects Pb+Pb24: set to true if for Pb+Pb 24, single mu4, do NOT want to perform trigger efficiency study

    bool turn_on_track_charge = false; // turn on if track charge is stored

    bool filter_out_photo_resn_for_trig_effcy = true; // if true: filter out pairs from photoproduction or resonance decay for trigger efficiency study as well
// --------------------- public class methods ---------------------------
	~DimuonDataAlgCoreT(){}
    void PrintInstructions();
    void ProcessDataHook();


    // --------------- InitInputHook ---------------
    void InitInputHook(){
        InitInput_DataCore();
        (CallInitInput<Extras>(), ...);
    }

    // --------------- InitInputBranchesSingleMuonAnalysisHook ---------------
    void InitInputBranchesSingleMuonAnalysisHook(){
        InitInputBranchesSingleMuonAnalysis_DataCore();
        (CallInitInputBranchesSingleMuonAnalysis<Extras>(), ...);
    }

    // --------------- InitInputBranchesDimuonAnalysisHook ---------------
    void InitInputBranchesDimuonAnalysisHook(){
        InitInputBranchesDimuonAnalysis_DataCore();
        (CallInitInputBranchesDimuonAnalysis<Extras>(), ...);
    }

    // --------------- InitTempVariablesHook ---------------
    void InitTempVariablesHook(){
        InitTempVariables_DataCore();
        (CallInitTempVariables<Extras>(), ...);
    }

    // --------------- InitParamsHook ---------------
    void InitParamsHook(){
        InitParams_DataCore();
        (CallInitParams<Extras>(), ...);
    }

    // --------------- InitOutputSettingsHookHook ---------------
    void InitOutputSettings(){
        InitOutputSettings_DataCore();
        (CallInitOutputSettings<Extras>(), ...);
    }

    // --------------- InitOutputTreesExtraHook ---------------
    void InitOutputTreesExtraHook(){
        (CallInitOutputTreesExtra<Extras>(), ...);
    }

    // --------------- InitOutputHistsExtraHook ---------------
    void InitOutputHistsExtraHook(){
        (CallInitOutputHistsExtra<Extras>(), ...);
    }

    // --------------- InitializeExtraHook ---------------
    void InitializeExtraHook(){
        InitializeExtra_DataCore();
        (CallInitializeExtra<Extras>(), ...);
    }

    // --------------- PassCutsHook ---------------
    bool PassCutsHook(){
        return (PassCuts_DataCore(bool requireTight) && (CallPassCuts<Extras>(), ...));
    }

    // --------------- FillMuonPairTreeHook ---------------
    void FillMuonPairTreeHook(){
        (CallFillMuonPairTree<Extras>(), ...);
    }

    // --------------- FillSingleMuonTreeHook ---------------
    void FillSingleMuonTreeHook(){
        (CallFillSingleMuonTree<Extras>(), ...);
    }

    // --------------- FillMuonPairHook ---------------
    void FillMuonPairHook(int pair_ind){
        FillMuonPair_DataCore(pair_ind);
        (CallFillMuonPair<Extras>(pair_ind), ...);
    }

    // --------------- HistAdjustHook ---------------
    void HistAdjustHook(){
        HistAdjust_DataCore();
        (CallHistAdjust<Extras>(), ...);
    }

   // --------------- PrintInstructionsHook ---------------
    void PrintInstructionsHook(){
        PrintInstructions_DataCore();
        (CallPrintInstructions<Extras>(), ...);
    }

    // --------------- PerformAdditionalPairAnalysisHook ---------------
    void PerformAdditionalPairAnalysisHook(){
        (CallAdditionalPairAnalysis<Extras>(), ...);
    }

    // --------------- FinalizeHook ---------------
    void FinalizeHook(){
        Finalize_DataCore();
        (CallFinalize<Extras>(), ...);
    }

};
