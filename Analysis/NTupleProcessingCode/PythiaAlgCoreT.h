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

    // Escape hatch only. -1 = derive from (sample type, isTestSample) -- the correct behaviour.
    // 0 = force pp-beam only | 1 = force 4 beams. Do NOT use to paper over a wrong isTestSample.
    // (isTestSample itself is PUBLIC -- see the config section below.)
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
    // ---- FULLSIM event bookkeeping written to meta_tree_out (added 2026-09-09) ----
    // WHY: fullsim_weight_factor is normalised to N_beam (ProcessDataHook), but the event
    // loop runs over N_proc = min(N_beam, nevents_max). A run truncated by nevents_max
    // therefore writes N_proc events' worth of pairs carrying an N_beam normalisation, i.e.
    // an absolute cross-section low by exactly N_proc/N_beam. Every other handle on the
    // denominator -- the weight itself, the NTUP chain entry count, and AMI totalEvents --
    // returns N_beam and they all AGREE with each other, so nothing could detect it. This is
    // reachable: the pipeline smoke test passes NEVENTS_MAX through the same run script while
    // the output suffix stays "_full", so a truncated file can land on the nominal path.
    // Recording both numbers makes the discrepancy visible to any consumer.
    // Flat index [ikin * nBeamTypes + ibeam]; sized once in InitOutputTreesExtra_PythiaCore
    // and never resized after, so the Branch() addresses stay valid.
    std::vector<Long64_t> meta_nproc_kn_beam;   // events actually looped over
    std::vector<Long64_t> meta_nbeam_kn_beam;   // events the weight is normalised to
    bool meta_fullsim_truncated = false;        // true if any (kn,beam) had N_proc < N_beam
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

    // Beam content of a fullsim run. Derived from (sample type, isTestSample) -- ONE switch,
    // so the input path and the isospin weight cannot drift apart. See FullSimSampleType.h.
    bool UseFourIsospinBeams() const {
        if (isospin_beams_override >= 0) return isospin_beams_override == 1;  // escape hatch
        return FullSimSampleUsesFourBeams(fullsim_sample_type, isTestSample);
    }
    void setIsospinBeams(bool four_beams) { isospin_beams_override = four_beams ? 1 : 0; }

    bool turn_data_resonance_cuts_on = false;
    bool fill_kn_trees_fullsim = false;  // set true to bin fullsim pairs into per-kn trees
    std::string fullsim_input_dir_override;  // if non-empty, replaces computed fullsim_input_dir

    // ---- THE fullsim sample switch: TEST sample vs FULL production ----
    // ONE flag drives the input directory, the AMI cross-section directory, AND the isospin
    // treatment, so they can never disagree. The authoritative rule lives in FullSimSampleType.h
    // (FullSimSampleInputDir / FullSimSampleUsesFourBeams):
    //   isTestSample = false (DEFAULT) = the FULL production -- the physics sample
    //       pp conditions  -> pp full sample,      pp beam only, isospin weight 1
    //       HIJING overlay -> overlay full sample, 4 beams,      Pb ratio 4:6:6:9
    //   isTestSample = true            = the small TEST sample
    //       pp conditions  -> pp24 test sample,    4 beams (produced that way BY MISTAKE)
    //       HIJING overlay -> overlay test sample, pp beam only (only beam produced)
    // Default false because the full production is the physics sample; a test-sample run must
    // declare itself. NOTE: a cross-section from the 4-beam pp24 TEST sample carries the Pb
    // isospin average and is therefore NOT a physical pp cross-section -- label it honestly.
    bool isTestSample = false;

    // HIJING-overlay CONDITIONS YEAR (FullSimSampleType.h): 24 (DEFAULT) = Pb+Pb 2024 conditions,
    // pythia_fullsim_hijing_overlay_test_sample/ (r17864) and the coming full production;
    // 23 = Pb+Pb 2023 conditions, ..._test_sample_pbpb23/ (r17618 / r17662). Drives the input
    // directory AND the output label (hijing_overlay_pbpb<yy>) together, so an output can never
    // carry the wrong year. Ignored for every non-overlay sample type.
    int overlay_pbpb_year = 24;

    // DIAGNOSTIC ONLY (default false = strict). When true, a missing pT-hat slice is a
    // warning instead of a fatal error. Required for single-slice studies (e.g. the r17662
    // signal-only-truth sample, which exists ONLY for pTH8_14). NEVER set this for a
    // cross-section-weighted production run: a missing slice biases the sigma-weighted
    // combination (which is exactly what the strict check exists to prevent).
    bool allow_missing_slices = false;

    // ---- AMI provenance (BLOCKING; see InitInputFullsim and FullSimSampleType.h) ----
    // The AMI weight is that of the PYTHIA EVGEN a sample was simulated from (never of its AOD
    // chain, never of a copy inside the sample directory). `ami_evgen` is DERIVED from
    // (fullsim_sample_type, isTestSample) -- the same switch as the input directory and the
    // isospin treatment: pp24 FULL -> PDF (803015-803020, pp only, proton PDF); everything else
    // -> nPDF (802758-802781, 4 isospins, nuclear PDF). AMI files are named by BEAM+SLICE only,
    // and the two evgens' cross-sections differ SLICE-DEPENDENTLY, so a wrong directory would
    // silently corrupt every sigma-weighted quantity (and cancels in no ratio).
    //   ami_info_dir_override : diagnostic escape hatch; if non-empty, replaces PythiaEvgenAmiDir(ami_evgen)
    //   expected_ami_dsids    : the datasetNumber in each AMI file MUST be in this list, else
    //                           InitInputFullsim throws. Defaults to PythiaEvgenDsids(ami_evgen);
    //                           a run script may narrow it (e.g. the single-DSID noovl sample).
    PythiaEvgen ami_evgen = PythiaEvgen::nPDF;
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
