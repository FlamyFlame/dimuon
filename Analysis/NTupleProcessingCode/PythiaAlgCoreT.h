#pragma once

#include <map>
#include <string>
#include <vector>
#include "../MuonObjectsParamsAndHelpers/Muon.h"
#include "../MuonObjectsParamsAndHelpers/MuonPairPythia.h"
#include "../MuonObjectsParamsAndHelpers/muon_pair_enums_MC.h"
#include "DimuonAlgCoreT.c"

template <class PairT, class MuonT, class Derived, class... Extras>
class PythiaAlgCoreT
    : public DimuonAlgCoreT<PairT, MuonT, Derived>
{
    template <class, class> friend class PythiaTruthExtras;
    template <class, class, class> friend class PythiaFullSimExtras;
    template <class> friend class PythiaFullSimOverlayExtras;

public:
    using pair_t = PairT;

protected:
    using Base = DimuonAlgCoreT<PairT, MuonT, Derived>;
    using Base::self;
    using Base::fChainRef;
    using Base::mpairRef;
    using Base::pmsRef;
    using Base::h_cutAcceptanceRef;

    static const int nBeamTypes = 4;

// --------------------- general settings ---------------------------

    int run_year = 23;
    bool isRun3 = true;
    bool is_fullsim = false;
    bool is_fullsim_overlay = false;
    bool perform_truth = true;
    bool useLocal = false;

    int batch_num = 0;
    int kn_batch = 0;

// --------------------- kin ranges & file structure ---------------------------

    int nKinRanges = 5;
    std::vector<std::string> kin_dirs;
    std::vector<float> kinRanges;
    std::vector<Long64_t> nevents;
    std::vector<Long64_t> nevents_accum;
    std::vector<Long64_t> njobs_accum;
    std::vector<Long64_t> njobs;
    std::vector<int> nevents_per_file;
    std::vector<Long64_t> njobs_all_files_combined;
    std::vector<std::string> job_dirs;
    std::vector<std::string> beam_dirs = {"pp/", "pn/", "np/", "nn/"};
    int nfiles_base[4] = {4, 6, 6, 9};

    std::vector<std::vector<bool>> kn_in_job;
    std::vector<std::vector<int>> nfiles_factor;

    bool new_run = false;
    std::string batch_suffix;
    std::string outfile_name;
    std::string outhistfile_name;
    std::string py_dir;
    std::string fullsim_input_dir;  // directory holding Pythia fullsim NTUP files

// --------------------- input files & trees & data for setting branches ---------------------------

    TChain* evChain = nullptr;
    TChain* metaChain = nullptr;
    std::vector<std::vector<TChain*>> evChains_kn_beam;
    std::vector<std::vector<Long64_t>> nentries_kn_beam;
    std::vector<Long64_t> nentries_kn_sum;
    std::vector<std::vector<double>> ami_weight_kn_beam;
    std::map<std::string, double> nominal_beam_ratio;

    double efficiency = 1.;
    double ev_weight = 0.;

    // Event-level kinematic branches
    // private: double; Q_float used for non-private "Q" (float)
    double QHard = 0.;
    float  Q_float = 0.f;
    double pTHat = -1000.;
    double mHat = -1000.;

    // Muon-pair branches — named after the ROOT branch name
    // private sample: vector<double>*; non-private: vector<float>* (_f suffix)
    std::vector<double>* truth_mupair_pt1  = nullptr;
    std::vector<float>*  truth_mupair_pt1_f = nullptr;
    std::vector<double>* truth_mupair_eta1 = nullptr;
    std::vector<float>*  truth_mupair_eta1_f = nullptr;
    std::vector<double>* truth_mupair_phi1 = nullptr;
    std::vector<float>*  truth_mupair_phi1_f = nullptr;
    std::vector<int>*    truth_mupair_ch1  = nullptr;
    std::vector<int>*    truth_mupair_bar1 = nullptr;

    std::vector<double>* truth_mupair_pt2  = nullptr;
    std::vector<float>*  truth_mupair_pt2_f = nullptr;
    std::vector<double>* truth_mupair_eta2 = nullptr;
    std::vector<float>*  truth_mupair_eta2_f = nullptr;
    std::vector<double>* truth_mupair_phi2 = nullptr;
    std::vector<float>*  truth_mupair_phi2_f = nullptr;
    std::vector<int>*    truth_mupair_ch2  = nullptr;
    std::vector<int>*    truth_mupair_bar2 = nullptr;

    // Fullsim-specific branches (bound by InitInputFullsim_PythiaCore)
    std::vector<int>*    truth_muon_barcode = nullptr; // barcodes of truth muons (indexed by truth muon index)

    // Fullsim per-beam/kn normalization factor: ami_weight * isospin_ratio / N_beam
    double fullsim_weight_factor = 1.0;

// --------------------- output trees ---------------------------

    std::vector<long> nentries_per_kin;
    TTree* meta_tree_out = nullptr;
    std::vector<std::vector<TTree*>> muonPairOutTreeKinRange;

    int current_ikin = 0;  // set during ProcessData for FillMuonPairTreePythia

// --------------------- getters ---------------------------

    bool getIsPrivate() const { return self().isPrivate; }
    bool getIsFullsim() const { return is_fullsim; }
    bool getIsFullsimOverlay() const { return is_fullsim_overlay; }
    bool getPerformTruth() const { return perform_truth; }
    bool getUseLocal() const { return useLocal; }
    void setIsFullsim(bool v) { is_fullsim = v; }
    void setIsFullsimOverlay(bool v) { is_fullsim_overlay = v; }
    void setPerformTruth(bool v) { perform_truth = v; }
    void setUseLocal(bool v) { useLocal = v; }

    std::vector<int>*& TruthMuonBarcodeRef() { return truth_muon_barcode; }

    TChain* GetChainForBranchSetup() const {
        if (getIsPrivate() && evChain) return evChain;
        if (!getIsPrivate() && !evChains_kn_beam.empty() && !evChains_kn_beam[0].empty())
            return evChains_kn_beam[0][0];
        return nullptr;
    }

// --------------------- core method declarations ---------------------------

    void InitInput_PythiaCore();
    void InitInputPrivate_PythiaCore();
    void InitInputCentrProd_PythiaCore();
    void InitInputFullsim_PythiaCore();
    void SetInputOutputFilesFromBatch_PythiaCore();
    void InputSanityCheck_PythiaCore();

    void InitParams_PythiaCore();
    void InitTempVariables_PythiaCore() {}
    void InitOutputTreesExtra_PythiaCore();
    void InitializeExtra_PythiaCore();

    void FillMuonPair_PythiaCore(int pair_ind);
    bool PassCuts_PythiaCore();
    void FillMuonPairTreePythia(int nkin);
    void HistAdjust_PythiaCore() {}
    void Finalize_PythiaCore();

// --------------------- Extras call helpers ---------------------------

    template <class E>
    void CallInitInput() {
        if constexpr (requires(Derived& d){ static_cast<E&>(d).InitInputExtra(); }) {
            static_cast<E&>(self()).InitInputExtra();
        }
    }
    template <class E>
    void CallInitParams() {
        if constexpr (requires(Derived& d){ static_cast<E&>(d).InitParamsExtra(); }) {
            static_cast<E&>(self()).InitParamsExtra();
        }
    }
    template <class E>
    void CallInitTempVariables() {
        if constexpr (requires(Derived& d){ static_cast<E&>(d).InitTempVariablesExtra(); }) {
            static_cast<E&>(self()).InitTempVariablesExtra();
        }
    }
    template <class E>
    void CallInitOutputTreesExtra() {
        if constexpr (requires(Derived& d){ static_cast<E&>(d).InitOutputTreesExtra(); }) {
            static_cast<E&>(self()).InitOutputTreesExtra();
        }
    }
    template <class E>
    void CallInitOutputHistsExtra() {
        if constexpr (requires(Derived& d){ static_cast<E&>(d).InitOutputHistsExtra(); }) {
            static_cast<E&>(self()).InitOutputHistsExtra();
        }
    }
    template <class E>
    void CallInitOutputExtra() {
        if constexpr (requires(Derived& d){ static_cast<E&>(d).InitOutputExtra(); }) {
            static_cast<E&>(self()).InitOutputExtra();
        }
    }
    template <class E>
    void CallInitializeExtra() {
        if constexpr (requires(Derived& d){ static_cast<E&>(d).InitializeExtra(); }) {
            static_cast<E&>(self()).InitializeExtra();
        }
    }
    template <class E>
    void CallPerformTruthPairAnalysis() {
        if constexpr (requires(Derived& d){ static_cast<E&>(d).PerformTruthPairAnalysis(); }) {
            static_cast<E&>(self()).PerformTruthPairAnalysis();
        }
    }
    template <class E>
    void CallPerPairCrossxUpdate() {
        if constexpr (requires(Derived& d){ static_cast<E&>(d).PerPairCrossxUpdate(); }) {
            static_cast<E&>(self()).PerPairCrossxUpdate();
        }
    }
    template <class E>
    void CallFillNumMuonPairsHist(int nAfter, double w) {
        if constexpr (requires(Derived& d){ static_cast<E&>(d).FillNumMuonPairsHist(nAfter, w); }) {
            static_cast<E&>(self()).FillNumMuonPairsHist(nAfter, w);
        }
    }
    template <class E>
    void CallProcessEventFullsim(int ev_num) {
        if constexpr (requires(Derived& d){ static_cast<E&>(d).ProcessEventFullsim(ev_num); }) {
            static_cast<E&>(self()).ProcessEventFullsim(ev_num);
        }
    }
    template <class E>
    void CallHistAdjust() {
        if constexpr (requires(Derived& d){ static_cast<E&>(d).HistAdjustExtra(); }) {
            static_cast<E&>(self()).HistAdjustExtra();
        }
    }
    template <class E>
    void CallFinalize() {
        if constexpr (requires(Derived& d){ static_cast<E&>(d).FinalizeExtra(); }) {
            static_cast<E&>(self()).FinalizeExtra();
        }
    }

public:
    // public: let user set batch and kn index before Run()
    int GetKnBatch() const { return kn_batch; }
    void SetKnBatch(int kn) { kn_batch = kn; }

    bool turn_data_resonance_cuts_on = false;
    bool fill_kn_trees_fullsim = false;  // set true to bin fullsim pairs into per-kn trees

    explicit PythiaAlgCoreT(int batch_num_input, bool useLocal_input = false)
        : batch_num(batch_num_input)
        , useLocal(useLocal_input)
    {}

    void ProcessDataHook();
    void PerformTruthPairAnalysisHook() {
        (CallPerformTruthPairAnalysis<Extras>(), ...);
    }
    void ProcessEventFullsimHook(int ev_num) {
        (CallProcessEventFullsim<Extras>(ev_num), ...);
    }
    void OutputTreePathHook();
    void OutputHistPathHook();

    void InitInputHook() {
        InitInput_PythiaCore();
        (CallInitInput<Extras>(), ...);
    }
    void InitParamsHook() {
        (CallInitParams<Extras>(), ...);
        InitParams_PythiaCore();
    }
    void InitTempVariablesHook() {
        InitTempVariables_PythiaCore();
        (CallInitTempVariables<Extras>(), ...);
    }
    void InitOutputTreesExtraHook() {
        InitOutputTreesExtra_PythiaCore();
        (CallInitOutputTreesExtra<Extras>(), ...);
    }
    void InitOutputHistsExtraHook() {
        (CallInitOutputHistsExtra<Extras>(), ...);
    }
    void InitOutputExtraHook() {
        (CallInitOutputExtra<Extras>(), ...);
    }
    void InitializeExtraHook() {
        InitializeExtra_PythiaCore();
        (CallInitializeExtra<Extras>(), ...);
    }
    bool PassCutsHook() {
        return PassCuts_PythiaCore();
    }
    void FillMuonPairHook(int pair_ind) {
        FillMuonPair_PythiaCore(pair_ind);
    }
    void FillMuonPairTreeHook() {
        FillMuonPairTreePythia(current_ikin);
    }
    void HistAdjustHook() {
        HistAdjust_PythiaCore();
        (CallHistAdjust<Extras>(), ...);
    }
    void FinalizeHook() {
        Finalize_PythiaCore();
        (CallFinalize<Extras>(), ...);
    }

    ~PythiaAlgCoreT() {}
};
