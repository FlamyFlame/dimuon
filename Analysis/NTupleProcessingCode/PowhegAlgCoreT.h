#pragma once

#include "../MuonObjectsParamsAndHelpers/MuonPairPowheg.h"
#include "../MuonObjectsParamsAndHelpers/muon_pair_enums_MC.h"
#include "DimuonAlgCoreT.c"

template <class PairT, class MuonT, class Derived, class... Extras>
class PowhegAlgCoreT
    : DimuonAlgCoreT<PairT, MuonT, Derived>
{
protected:
// --------------------- general settings ---------------------------

    int file_batch;
    std::string mc_mode;
    bool   is_fullsim; // if is fullsim, load reco quantities
    
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
  
    TChain          *fChain;   //!pointer to the analyzed TTree or TChain
  
    // Declaration of leaf types
    float                 Q;
    std::vector<float>   *EventWeights           =nullptr;
  
  // --------------------- output file, histograms & trees ---------------------------
  
    TTree* meta_tree;
    static const int nAncestorGroups = 4;

    // configuration
    std::string sign_labels[ParamsSet::nSigns] = {"same_sign", "op_sign"};
    std::string signs[ParamsSet::nSigns] = {"_ss", "_op"};
    std::string dphis[2] = {"_near", "_away"};
    std::string ancestor_grps[nAncestorGroups] = {"_gg", "_qg","_single_g","_qq"};



    // std::vector<std::string> parentGroupLabels = {"direct b","b to c","direct c","s/light","single photon", "Drell-Yan"};
    std::vector<std::string> parentGroupLabels = {"direct b","b to c","direct c","s/light","single photon"};
    int nParentGroups = parentGroupLabels.size();    
    
    std::vector<std::string> ancestor_labels = {"gg", "gq", "single g", "q qbar"};
    
    std::vector<std::string> samePrtsLabels = {"Same Parents", "Different Parents"};
    std::vector<std::string> bb_op_one_b_one_btoc_labels = {"Same b", "Involve osc(s)", "From different ancestors", "Others"};

    // std::vector<std::string> osc_labels = {"0 osc, one b one b-to-c", "0 osc, others", "1 osc, one b one b-to-c", "1 osc, regular", "2 oscs, one b one b-to-c", "2 oscs, others"};
    std::vector<std::string> osc_labels = {"0 osc, one c-tag", "0 osc, others(*)", "1 osc, one c-tag(*)", "1 osc, regular", "2 oscs, one c-tag", "2 oscs, others(*)"};
    std::vector<std::string> num_hard_scatt_out_labels = {"2","3","more"};

// --------------------- class methods ---------------------------
  
    PowhegAlgCoreT(int file_batch_input, std::string mc_mode_input, bool is_fullsim_input, bool is_fullsim_overlay_input = false)
    : is_fullsim (is_fullsim_input)
    : is_fullsim_overlay (is_fullsim_overlay_input)
    , file_batch (file_batch_input) {
        mc_mode = mc_mode_input;

        is_fullsim |= is_fullsim_overlay;
        
        cutLabels = cutLabels_MC;
        numCuts = cutLabels_MC.size();
        crossx_cut = 5 * pow(10,8);
        std::cout << "Powheg Ntuple processing script:" << std::endl;
        std::cout << "mc_mode: " << mc_mode << std::endl;
        std::cout << "file_batch: " << file_batch << std::endl;
        std::cout << "is_fullsim? " << is_fullsim << std::endl;
        std::cout << "is_fullsim_overlay? " << is_fullsim_overlay << std::endl;
    }

    
    
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
        InitParams_PowhegCore();
        (CallInitParams<Extras>(), ...);
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

    void InitOutputTreesExtra_PowhegCore(){
        meta_tree = new TTree("meta_tree","meta_tree");
        meta_tree->Branch("nentries_before_cuts", &nentries, "nentries_before_cuts/L");
    }

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
    bool perform_truth {false}; // if true, perform truth analysis even if fullsim

    ~PowhegAlgCoreT(){}
};
