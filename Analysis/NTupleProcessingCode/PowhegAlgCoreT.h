#pragma once

#include "../MuonObjectsParamsAndHelpers/Muon.h"
#include "../MuonObjectsParamsAndHelpers/MuonPairPowheg.h"
#include "../MuonObjectsParamsAndHelpers/muon_pair_enums_MC.h"
#include "DimuonAlgCoreT.c"

template <class PairT, class MuonT, class Derived, class... Extras>
class PowhegAlgCoreT
    : public DimuonAlgCoreT<PairT, MuonT, Derived>
{
    template<class, class> friend class PowhegTruthExtras;  // match your template params
    template<class, class, class> friend class PowhegFullSimExtras;  // match your template params
    template<class> friend class PowhegFullSimOverlayExtras;  // match your template params


public:
    using pair_t = PairT;

protected:    
    using Base = DimuonAlgCoreT<PairT, MuonT, Derived>;
    using Base::self;
    using Base::fChainRef;
    using Base::mpairRef;
    using Base::pmsRef;
    using Base::h_cutAcceptanceRef;

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

    std::string mcdir;
    std::string data_subdir;
    std::string dt_suffix;
    std::string truth_suffix;

// --------------------- input files & trees & data for setting branches---------------------------
  
    // Declaration of leaf types
    std::vector<float>   *EventWeights           =nullptr;
  
    std::vector<int>     *truth_muon_barcode        =nullptr;

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

  // --------------------- reference API, setters & getters ---------------------------

    std::string&                mcModeRef() { return mc_mode; }
    std::string&                mcdirRef()   { return mcdir; }
    std::vector<float> *&       EventWeightsRef()   { return EventWeights; }
    
    bool                        getIsFullsim(){ return is_fullsim; }
    bool                        getIsFullsimOverlay(){ return is_fullsim_overlay; }
    bool                        getPerformTruth(){ return perform_truth; }

    void                        setIsFullsim(bool is_fullsim_input){ is_fullsim = is_fullsim_input; }
    void                        setIsFullsimOverlay(bool is_fullsim_overlay_input){ is_fullsim_overlay = is_fullsim_overlay_input; }
    void                        setPerformTruth(bool perform_truth_input){ perform_truth = perform_truth_input; }

// --------------------- class methods ---------------------------
  
    PowhegAlgCoreT(int file_batch_input, std::string mc_mode_input)
    : file_batch (file_batch_input)
    , mc_mode (mc_mode_input){
        
        crossx_cut = 5 * pow(10,8);
        std::cout << "Powheg Ntuple processing script:" << std::endl;
        std::cout << "mc_mode: " << mc_mode << std::endl;
    }

    // --------------- CheckBranchPtrs ---------------
        template <class E>
    void CallCheckBranchPtrs() {
        if constexpr (requires(Derived& d){ static_cast<E&>(d).CheckBranchPtrsExtra(); }) {
            static_cast<E&>(self()).CheckBranchPtrsExtra();
        }
    }

    void CheckBranchPtrs_PowhegCore();

    void CheckBranchPtrs(){
        CheckBranchPtrs_PowhegCore();
        (CallCheckBranchPtrs<Extras>(), ...);
    }

    // --------------- InitInputImpl ---------------
        template <class E>
    void CallInitInput() {
        if constexpr (requires(Derived& d){ static_cast<E&>(d).InitInputExtra(); }) {
            static_cast<E&>(self()).InitInputExtra();
        }
    }

    void InitInput_PowhegCore();

    // --------------- InitParamsImpl ---------------
        template <class E>
    void CallInitParams() {
        if constexpr (requires(Derived& d){ static_cast<E&>(d).InitParamsExtra(); }) {
            static_cast<E&>(self()).InitParamsExtra();
        }
    }

    void InitParams_PowhegCore();

    // --------------- InitTempVariablesImpl ---------------
        template <class E>
    void CallInitTempVariables() {
        if constexpr (requires(Derived& d){ static_cast<E&>(d).InitTempVariablesExtra(); }) {
            static_cast<E&>(self()).InitTempVariablesExtra();
        }
    }

    void InitTempVariables_PowhegCore(){}

    // --------------- InitOutputSettingsHookImpl ---------------
        template <class E>
    void CallInitOutputSettings() {
        if constexpr (requires(Derived& d){ static_cast<E&>(d).InitOutputSettingsExtra(); }) {
            static_cast<E&>(self()).InitOutputSettingsExtra();
        }
    }

    void InitOutputSettings_PowhegCore(){}

    // --------------- InitOutputTreesExtraImpl ---------------
        template <class E>
    void CallInitOutputTreesExtra() {
        if constexpr (requires(Derived& d){ static_cast<E&>(d).InitOutputTreesExtra(); }) {
            static_cast<E&>(self()).InitOutputTreesExtra();
        }
    }

    void InitOutputTreesExtra_PowhegCore();

    // --------------- InitOutputHistsExtraImpl ---------------
        template <class E>
    void CallInitOutputHistsExtra() {
        if constexpr (requires(Derived& d){ static_cast<E&>(d).InitOutputHistsExtra(); }) {
            static_cast<E&>(self()).InitOutputHistsExtra();
        }
    }

    // --------------- InitOutputExtraImpl ---------------
        template <class E>
    void CallInitOutputExtra() {
        if constexpr (requires(Derived& d){ static_cast<E&>(d).InitOutputExtra(); }) {
            static_cast<E&>(self()).InitOutputExtra();
        }
    }

    // --------------- InitializeExtraImpl ---------------
        template <class E>
    void CallInitializeExtra() {
        if constexpr (requires(Derived& d){ static_cast<E&>(d).InitializeExtra(); }) {
            static_cast<E&>(self()).InitializeExtra();
        }
    }

    void InitializeExtra_PowhegCore(){}

    // --------------- PassCutsImpl ---------------
        template <class E>
    bool CallPassCuts() {
        if constexpr (requires(Derived& d){ static_cast<E&>(d).PassCutsExtra(); }) {
            return static_cast<E&>(self()).PassCutsExtra();
        } else {
            return true;
        }
    }

    bool PassCuts_PowhegCore();

    // --------------- FillMuonPairTreeImpl ---------------
        template <class E>
    void CallFillMuonPairTree() {
        if constexpr (requires(Derived& d){ static_cast<E&>(d).FillMuonPairTreeExtra(); }) {
            static_cast<E&>(self()).FillMuonPairTreeExtra();
        }
    }

    // --------------- FillMuonPairImpl ---------------
        template <class E>
    void CallFillMuonPair(int pair_ind) {
        if constexpr (requires(Derived& d, int ind) { static_cast<E&>(d).FillMuonPairExtra(ind); }) {
            static_cast<E&>(self()).FillMuonPairExtra(pair_ind);
        }
    }

    void FillMuonPair_PowhegCore(int pair_ind);

    // --------------- HistAdjustImpl ---------------
        template <class E>
    void CallHistAdjust() {
        if constexpr (requires(Derived& d){ static_cast<E&>(d).HistAdjustExtra(); }) {
            static_cast<E&>(self()).HistAdjustExtra();
        }
    }

    void HistAdjust_PowhegCore(){}

    // --------------- PerformTruthPairAnalysisImpl ---------------
        template <class E>
    void CallTruthPairAnalysis() {
        if constexpr (requires(Derived& d){ static_cast<E&>(d).PerformTruthPairAnalysis(); }) {
            static_cast<E&>(self()).PerformTruthPairAnalysis();
        }
    }

    // --------------- ProcessEventFullsimImpl ---------------
        template <class E>
    void CallProcessEventFullsim(int ev_num) {
        if constexpr (requires(Derived& d, int num){ static_cast<E&>(d).ProcessEventFullsim(num); }) {
            static_cast<E&>(self()).ProcessEventFullsim(ev_num);
        }
    }

    // --------------- FinalizeImpl ---------------
        template <class E>
    void CallFinalize() {
        if constexpr (requires(Derived& d){ static_cast<E&>(d).FinalizeExtra(); }) {
            static_cast<E&>(self()).FinalizeExtra();
        }
    }

    void Finalize_PowhegCore(){}

public:
    void ProcessDataHook();
    void ProcessEventTruthOnly(int ev_num);
    void OutputTreePathHook();
    void OutputHistPathHook();

    // --------------- InitInputHook ---------------
    void InitInputHook(){
        InitInput_PowhegCore();
        (CallInitInput<Extras>(), ...);
    }

    // --------------- InitParamsHook ---------------
    void InitParamsHook(){
        (CallInitParams<Extras>(), ...); // first init parameters from the subanalyses
        InitParams_PowhegCore();
    }

    // --------------- InitTempVariablesHook ---------------
    void InitTempVariablesHook(){
        InitTempVariables_PowhegCore();
        (CallInitTempVariables<Extras>(), ...);
    }

    // --------------- InitOutputSettingsHookHook ---------------
    void InitOutputSettings(){
        InitOutputSettings_PowhegCore();
        (CallInitOutputSettings<Extras>(), ...);
    }

    // --------------- InitOutputTreesExtraHook ---------------
    void InitOutputTreesExtraHook(){
        InitOutputTreesExtra_PowhegCore();
        (CallInitOutputTreesExtra<Extras>(), ...);
    }

    // --------------- InitOutputHistsExtraHook ---------------
    void InitOutputHistsExtraHook(){
        (CallInitOutputHistsExtra<Extras>(), ...);
    }

    // --------------- InitOutputExtraHook ---------------
    void InitOutputExtraHook(){
        (CallInitOutputExtra<Extras>(), ...);
    }

    // --------------- InitializeExtraHook ---------------
    void InitializeExtraHook(){
        InitializeExtra_PowhegCore();
        (CallInitializeExtra<Extras>(), ...);
    }

    // --------------- PassCutsHook ---------------
    bool PassCutsHook(){
        return (PassCuts_PowhegCore() && (CallPassCuts<Extras>(), ...));
    }

    // --------------- FillMuonPairTreeHook ---------------
    void FillMuonPairTreeHook(){
        (CallFillMuonPairTree<Extras>(), ...);
    }

    // --------------- FillMuonPairHook ---------------
    void FillMuonPairHook(int pair_ind){
        FillMuonPair_PowhegCore(pair_ind);
        (CallFillMuonPair<Extras>(pair_ind), ...);
    }

    // --------------- HistAdjustHook ---------------
    void HistAdjustHook(){
        HistAdjust_PowhegCore();
        (CallHistAdjust<Extras>(), ...);
    }

    // --------------- PerformTruthPairAnalysisHook ---------------
    void PerformTruthPairAnalysisHook(){
        (CallTruthPairAnalysis<Extras>(), ...);
    }

    // --------------- ProcessEventFullsimHook ---------------
    void ProcessEventFullsimHook(int ev_num){
        (CallProcessEventFullsim<Extras>(ev_num), ...);
    }

    // --------------- FinalizeHook ---------------
    void FinalizeHook(){
        Finalize_PowhegCore();
        (CallFinalize<Extras>(), ...);
    }

public :
    ~PowhegAlgCoreT(){}
};
