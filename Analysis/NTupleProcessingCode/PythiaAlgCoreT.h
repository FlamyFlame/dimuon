#pragma once

#include <map>
#include <string>
#include <vector>
#include "../MuonObjectsParamsAndHelpers/Muon.h"
#include "../MuonObjectsParamsAndHelpers/MuonPairPythia.h"
#include "../MuonObjectsParamsAndHelpers/muon_pair_enums_MC.h"
#include "../MuonObjectsParamsAndHelpers/FullSimSampleType.h"
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
    FullSimSampleType fullsim_sample_type = FullSimSampleType::pp;
    bool perform_truth = true;
    bool useLocal = false;
    bool only_pp_isospin = false;   // TRUTH path only (InitInputCentrProd)

    // ---- Isospin content of the SIMULATED COLLISION SYSTEM (fullsim paths) ----
    // The beam content must match the system the sample simulates:
    //   pp CONDITIONS fullsim (FullSimSampleType::pp) simulates pp collisions
    //     -> ONE beam (pp), and the isospin weight is 1 (there is nothing to average).
    //   PbPb CONDITIONS HIJING overlay simulates Pb+Pb, whose nucleons are a p/n mix
    //     -> FOUR beams {pp,pn,np,nn}, combined with the Pb ratio 4:6:6:9.
    // The TEST samples are exceptions in BOTH directions -- they are what happens to
    // exist on disk, not what the physics wants:
    //   pp24 test sample        = 4 beams (a production mistake; kept so the existing
    //                             test-sample results stay reproducible)
    //   HIJING overlay test smp = pp beam only (the 4-beam overlay full sample is in
    //                             production)
    // so each run must be able to override the sample-type default.
    // -1 = use the sample-type default | 0 = force pp-beam only | 1 = force 4 beams
    int isospin_beams_override = -1;

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
    FullSimSampleType getFullSimSampleType() const { return fullsim_sample_type; }
    void setIsFullsim(bool v) { is_fullsim = v; }
    void setIsFullsimOverlay(bool v) { is_fullsim_overlay = v; }
    void setFullSimSampleType(FullSimSampleType t) { fullsim_sample_type = t; }

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
    void FillMuonPairTreePythia(int nkin);        // global pair tree + kn tree
    void FillMuonPairTreeKinRangePythia(int nkin); // kn tree only
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

    // Beam content of a fullsim run (see isospin_beams_override above).
    // Default: the HIJING overlay simulates Pb+Pb -> 4 isospin beams (4:6:6:9);
    // pp-conditions fullsim simulates pp -> the pp beam alone, isospin weight 1.
    // Public so a run script can override it for the two test samples, which are
    // exceptions in both directions.
    bool UseFourIsospinBeams() const {
        if (isospin_beams_override >= 0) return isospin_beams_override == 1;
        return FullSimSampleIsOverlay(fullsim_sample_type);
    }
    void setIsospinBeams(bool four_beams) { isospin_beams_override = four_beams ? 1 : 0; }

    bool turn_data_resonance_cuts_on = false;
    bool fill_kn_trees_fullsim = false;  // set true to bin fullsim pairs into per-kn trees
    std::string fullsim_input_dir_override;  // if non-empty, replaces computed fullsim_input_dir

    // DIAGNOSTIC ONLY (default false = strict). When true, a missing pT-hat slice is a
    // warning instead of a fatal error. Required for single-slice studies (e.g. the r17662
    // signal-only-truth sample, which exists ONLY for pTH8_14). NEVER set this for a
    // cross-section-weighted production run: a missing slice biases the sigma-weighted
    // combination (which is exactly what the strict check exists to prevent).
    bool allow_missing_slices = false;

    // ---- AMI provenance (BLOCKING; see InitInputFullsim) ----
    // AMI files are named by BEAM+SLICE only, so they do NOT identify the production. The pp24
    // TEST sample (802758-802781) and the pp24 FULL "_pdf" sample (803015-803020) have different,
    // slice-dependent cross-sections, so reading the wrong one silently corrupts every
    // sigma-weighted quantity (and does NOT cancel in ratios).
    //   ami_info_dir_override : if non-empty, replaces <py_dir>/ami_info/
    //   expected_ami_dsids    : if non-empty, the datasetNumber in each AMI file MUST be in this
    //                           list, else InitInputFullsim throws.
    std::string ami_info_dir_override;
    std::vector<int> expected_ami_dsids;
    // Propagate trigger decisions/matching from the trigger-enabled MC skims (_July2026)
    // into the output trees (m1/m2.passmu4, pair pass2mu4). Adds "_mc_trig" to the output
    // file name so nominal outputs are never clobbered. Input files that lack the trigger
    // branches (old trigger-off skims) are skipped entirely — never default-filled, which
    // would bias any efficiency built from the output. See
    // docs/tracking/mc_trigger_efficiency.md (Physics Procedure + R1/R2).
    bool store_mc_trigger = false;

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
    // Called from DimuonAlgCoreT::FillMuonPairTree(), which has ALREADY filled the
    // global pair tree -- so only the kinematic-range tree is filled here.  (The
    // Pythia truth path bypasses this hook and calls FillMuonPairTreePythia directly.)
    void FillMuonPairTreeHook() {
        FillMuonPairTreeKinRangePythia(current_ikin);
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
