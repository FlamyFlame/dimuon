#pragma once

#include "../MuonObjectsParamsAndHelpers/Muon.h"
#include "../MuonObjectsParamsAndHelpers/MuonPairPowheg.h"
#include "../MuonObjectsParamsAndHelpers/muon_pair_enums_MC.h"
#include "DimuonAlgCoreT.c"

template <class PairT, class MuonT, class Derived, class... Extras>
class PowhegAlgCoreT
    : public DimuonAlgCoreT<PairT, MuonT, Derived>
{
protected:
    auto& fChain()   { return this->fChain; }
    auto& mpair()   { return this->mpair; }
    auto& h_cutAcceptance()   { return this->h_cutAcceptance; }

// --------------------- general settings ---------------------------

    int file_batch;
    std::string mc_mode;
    bool    is_fullsim{false}; // if is fullsim, load reco quantities
    bool    is_fullsim_overlay{false}; // if is fullsim overlay, load centrality info
    bool    perform_truth{false};
    int     fullsim_run_yr {-1};
    double crossx_cut;
    double filter_effcy;
    double filter_effcy_bb = 0.003;
    double filter_effcy_cc = 0.001108;
    // static const int nMCmodes = 2;

    ParamsSet pms;

    std::string mcdir;
    std::string data_subdir;
    std::string dt_suffix;
    std::string truth_suffix;

// --------------------- input files & trees & data for setting branches---------------------------
  
    // Declaration of leaf types
    float                 Q;
    std::vector<float>   *EventWeights           =nullptr;
  
    std::vector<float>   *truth_mupair_pt1           =nullptr;
    std::vector<float>   *truth_mupair_eta1       =nullptr;
    std::vector<float>   *truth_mupair_phi1          =nullptr;
    std::vector<int>     *truth_mupair_ch1         =nullptr;
    std::vector<int>     *truth_mupair_bar1         =nullptr;

    std::vector<float>   *truth_mupair_pt2       =nullptr;
    std::vector<float>   *truth_mupair_eta2          =nullptr;
    std::vector<float>   *truth_mupair_phi2         =nullptr;
    std::vector<int>     *truth_mupair_ch2         =nullptr;
    std::vector<int>     *truth_mupair_bar2         =nullptr;

    std::vector<float>   *truth_mupair_asym         =nullptr;
    std::vector<float>   *truth_mupair_acop         =nullptr;
    std::vector<float>   *truth_mupair_pt         =nullptr;
    std::vector<float>   *truth_mupair_y         =nullptr;
    std::vector<float>   *truth_mupair_m         =nullptr;
    
  // --------------------- output file, histograms & trees ---------------------------
  
    TTree* meta_tree;

// --------------------- class methods ---------------------------
  
    PowhegAlgCoreT(int file_batch_input, std::string mc_mode_input)
    : file_batch (file_batch_input)
    , mc_mode (mc_mode_input){
        
        crossx_cut = 5 * pow(10,8);
        std::cout << "Powheg Ntuple processing script:" << std::endl;
        std::cout << "mc_mode: " << mc_mode << std::endl;
    }

    void ProcessDataHook();
    void OutputTreePathHook();
    void OutputHistPathHook();

    // --------------- InitInputHook ---------------
    template <class E>
    void CallInitInput() {
        if constexpr (requires(E& e){ e.InitInputExtra(); }) {
            static_cast<E&>(*this).InitInputExtra();
        }
    }

    void InitInput_PowhegCore();

    void InitInputHook(){
        InitInput_PowhegCore();
        (CallInitInput<Extras>(), ...);
    }

    // --------------- InitParamsHook ---------------
    template <class E>
    void CallInitParams() {
        if constexpr (requires(E& e){ e.InitParamsExtra(); }) {
            static_cast<E&>(*this).InitParamsExtra();
        }
    }

    void InitParams_PowhegCore();

    void InitParamsHook(){
        (CallInitParams<Extras>(), ...); // first init parameters from the subanalyses
        InitParams_PowhegCore();
    }

    // --------------- InitTempVariablesHook ---------------
    template <class E>
    void CallInitTempVariables() {
        if constexpr (requires(E& e){ e.InitTempVariablesExtra(); }) {
            static_cast<E&>(*this).InitTempVariablesExtra();
        }
    }

    void InitTempVariables_PowhegCore(){}

    void InitTempVariablesHook(){
        InitTempVariables_PowhegCore();
        (CallInitTempVariables<Extras>(), ...);
    }

    // --------------- InitOutputSettingsHookHook ---------------
    template <class E>
    void CallInitOutputSettings() {
        if constexpr (requires(E& e){ e.InitOutputSettingsExtra(); }) {
            static_cast<E&>(*this).InitOutputSettingsExtra();
        }
    }

    void InitOutputSettings_PowhegCore(){}

    void InitOutputSettings(){
        InitOutputSettings_PowhegCore();
        (CallInitOutputSettings<Extras>(), ...);
    }

    // --------------- InitOutputTreesExtraHook ---------------
    template <class E>
    void CallInitOutputTreesExtra() {
        if constexpr (requires(E& e){ e.InitOutputTreesExtra(); }) {
            static_cast<E&>(*this).InitOutputTreesExtra();
        }
    }

    void InitOutputTreesExtra_PowhegCore();

    void InitOutputTreesExtraHook(){
        InitOutputTreesExtra_PowhegCore();
        (CallInitOutputTreesExtra<Extras>(), ...);
    }

    // --------------- InitializeExtraHook ---------------
    template <class E>
    void CallInitializeExtra() {
        if constexpr (requires(E& e){ e.InitializeExtra(); }) {
            static_cast<E&>(*this).InitializeExtra();
        }
    }

    void InitializeExtra_PowhegCore(){}

    void InitializeExtraHook(){
        InitializeExtra_PowhegCore();
        (CallInitializeExtra<Extras>(), ...);
    }

    // --------------- PassCutsHook ---------------
    template <class E>
    bool CallPassCuts() {
        if constexpr (requires(E& e){ e.PassCutsExtra(); }) {
            return static_cast<E&>(*this).PassCutsExtra();
        } else {
            return true;
        }
    }

    bool PassCuts_PowhegCore();

    bool PassCutsHook(){
        return (PassCuts_PowhegCore() && (CallPassCuts<Extras>(), ...));
    }

    // --------------- FillMuonPairTreeHook ---------------
    template <class E>
    void CallFillMuonPairTree() {
        if constexpr (requires(E& e){ e.FillMuonPairTreeExtra(); }) {
            static_cast<E&>(*this).FillMuonPairTreeExtra();
        }
    }

    void FillMuonPairTreeHook(){
        (CallFillMuonPairTree<Extras>(), ...);
    }

    // --------------- FillMuonPairHook ---------------
    template <class E>
    void CallFillMuonPair(int pair_ind) {
        if constexpr (requires(E& e, int ind) { e.FillMuonPairExtra(ind); }) {
            static_cast<E&>(*this).FillMuonPairExtra(pair_ind);
        }
    }

    void FillMuonPair_PowhegCore(int pair_ind);

    void FillMuonPairHook(int pair_ind){
        FillMuonPair_PowhegCore(pair_ind);
        (CallFillMuonPair<Extras>(pair_ind), ...);
    }

    // --------------- HistAdjustHook ---------------
    template <class E>
    void CallHistAdjust() {
        if constexpr (requires(E& e){ e.HistAdjustExtra(); }) {
            static_cast<E&>(*this).HistAdjustExtra();
        }
    }

    void HistAdjust_PowhegCore(){}

    void HistAdjustHook(){
        HistAdjust_PowhegCore();
        (CallHistAdjust<Extras>(), ...);
    }

    // --------------- PerformTruthPairAnalysisHook ---------------
    template <class E>
    void CallTruthPairAnalysis() {
        if constexpr (requires(E& e){ e.TruthPairAnalysis(); }) {
            static_cast<E&>(*this).TruthPairAnalysis();
        }
    }

    void PerformTruthPairAnalysisHook(){
        (CallTruthPairAnalysis<Extras>(), ...);
    }

    // --------------- FinalizeHook ---------------
    template <class E>
    void CallFinalize() {
        if constexpr (requires(E& e){ e.FinalizeExtra(); }) {
            static_cast<E&>(*this).FinalizeExtra();
        }
    }

    void Finalize_PowhegCore(){}

    void FinalizeHook(){
        Finalize_PowhegCore();
        (CallFinalize<Extras>(), ...);
    }

public :
    ~PowhegAlgCoreT(){}
};
