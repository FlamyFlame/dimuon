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
    using Base = DimuonAlgCoreT<PairT, MuonT, Derived>;
    using Base::self;
    using Base::fChainRef;
    using Base::mpairRef;
    using Base::pmsRef;
    using Base::h_cutAcceptanceRef;

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
    void TChainFill(){ return self().TChainFillHook(); }

    bool IsPhotoProduction();


    // --------------- TChainFillImpl ---------------
        template <class E>
    void CallTChainFill() {
        // Evaluate access in Derived context so Extras may keep TChainFill(int) protected
        if constexpr (requires(Derived& d){ static_cast<E&>(d).PerformTChainFill(); }) {
            static_cast<E&>(self()).PerformTChainFill();
        }
    }

    // --------------- InitInputImpl ---------------
    template <class E>
    void CallInitInput() {
        if constexpr (requires(Derived& d){ static_cast<E&>(d).InitInputExtra(); }) {
            static_cast<E&>(self()).InitInputExtra();
        }
    }

    void InitInput_DataCore();

    // --------------- InitInputBranchesSingleMuonAnalysisImpl ---------------
    template <class E>
    void CallInitInputBranchesSingleMuonAnalysis() {
        if constexpr (requires(Derived& d){ static_cast<E&>(d).InitInputBranchesSingleMuonAnalysisExtra(); }) {
            static_cast<E&>(self()).InitInputBranchesSingleMuonAnalysisExtra();
        }
    }

    void InitInputBranchesSingleMuonAnalysis_DataCore();

    // --------------- InitInputBranchesDimuonAnalysisImpl ---------------
    template <class E>
    void CallInitInputBranchesDimuonAnalysis() {
        if constexpr (requires(Derived& d){ static_cast<E&>(d).InitInputBranchesDimuonAnalysisExtra(); }) {
            static_cast<E&>(self()).InitInputBranchesDimuonAnalysisExtra();
        }
    }

    void InitInputBranchesDimuonAnalysis_DataCore();

    // --------------- InitTempVariablesImpl ---------------
    template <class E>
    void CallInitTempVariables() {
        if constexpr (requires(Derived& d){ static_cast<E&>(d).InitTempVariablesExtra(); }) {
            static_cast<E&>(self()).InitTempVariablesExtra();
        }
    }

    void InitTempVariables_DataCore(){}

    // --------------- InitParamsImpl ---------------
    template <class E>
    void CallInitParams() {
        if constexpr (requires(Derived& d){ static_cast<E&>(d).InitParamsExtra(); }) {
            static_cast<E&>(self()).InitParamsExtra();
        }
    }

    void InitParams_DataCore();

    // --------------- InitOutputSettingsImpl ---------------
    template <class E>
    void CallInitOutputSettings() {
        if constexpr (requires(Derived& d){ static_cast<E&>(d).InitOutputSettingsExtra(); }) {
            static_cast<E&>(self()).InitOutputSettingsExtra();
        }
    }

    void InitOutputSettings_DataCore(){}

    // --------------- InitOutputTreesExtraImpl ---------------
    template <class E>
    void CallInitOutputTreesExtra() {
        if constexpr (requires(Derived& d){ static_cast<E&>(d).InitOutputTreesExtra(); }) {
            static_cast<E&>(self()).InitOutputTreesExtra();
        }
    }

    // --------------- InitOutputHistsExtraImpl ---------------
    template <class E>
    void CallInitOutputHistsExtra() {
        if constexpr (requires(Derived& d){ static_cast<E&>(d).InitOutputHistsExtra(); }) {
            static_cast<E&>(self()).InitOutputHistsExtra();
        }
    }

    // --------------- InitializeExtraImpl ---------------
    template <class E>
    void CallInitializeExtra() {
        if constexpr (requires(Derived& d){ static_cast<E&>(d).InitializeExtra(); }) {
            static_cast<E&>(self()).InitializeExtra();
        }
    }

    void InitializeExtra_DataCore(){}

    // --------------- PassCutsImpl ---------------
    template <class E>
    bool CallPassCuts() {
        if constexpr (requires(Derived& d){ static_cast<E&>(d).PassCutsExtra(); }) {
            return static_cast<E&>(self()).PassCutsExtra();
        } else {
            return true;
        }
    }

    bool PassCuts_DataCore(bool requireTight);

    // --------------- FillMuonPairTreeImpl ---------------
    template <class E>
    void CallFillMuonPairTree() {
        if constexpr (requires(Derived& d){ static_cast<E&>(d).FillMuonPairTreeExtra(); }) {
            static_cast<E&>(self()).FillMuonPairTreeExtra();
        }
    }

    // --------------- FillSingleMuonTreeImpl ---------------
    template <class E>
    void CallFillSingleMuonTree() {
        if constexpr (requires(Derived& d){ static_cast<E&>(d).FillSingleMuonTreeExtra(); }) {
            static_cast<E&>(self()).FillSingleMuonTreeExtra();
        }
    }

    // --------------- FillMuonPairImpl ---------------
    template <class E>
    void CallFillMuonPair(int pair_ind) {
        if constexpr (requires(E& e, int ind) { e.FillMuonPairExtra(ind); }) {
            static_cast<E&>(self()).FillMuonPairExtra(pair_ind);
        }
    }

    void FillMuonPair_DataCore(int pair_ind);

    // --------------- HistAdjustImpl ---------------
    template <class E>
    void CallHistAdjust() {
        if constexpr (requires(Derived& d){ static_cast<E&>(d).HistAdjustExtra(); }) {
            static_cast<E&>(self()).HistAdjustExtra();
        }
    }

    void HistAdjust_DataCore(){}

   // --------------- PrintInstructionsImpl ---------------
    template <class E>
    void CallPrintInstructions() {
        if constexpr (requires(Derived& d){ static_cast<E&>(d).PrintInstructionsExtra(); }) {
            static_cast<E&>(self()).PrintInstructionsExtra();
        }
    }

    void PrintInstructions_DataCore();

    // --------------- PerformAdditionalPairAnalysisImpl ---------------
    template <class E>
    void CallAdditionalPairAnalysis() {
        if constexpr (requires(Derived& d){ static_cast<E&>(d).AdditionalPairAnalysis(); }) {
            static_cast<E&>(self()).AdditionalPairAnalysis();
        }
    }

    // --------------- FinalizeImpl ---------------
    template <class E>
    void CallFinalize() {
        if constexpr (requires(Derived& d){ static_cast<E&>(d).FinalizeExtra(); }) {
            static_cast<E&>(self()).FinalizeExtra();
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

    // --------------- TChainFillHook ---------------
    void TChainFillHook(){
        (CallTChainFill<Extras>(), ...);
    }

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
