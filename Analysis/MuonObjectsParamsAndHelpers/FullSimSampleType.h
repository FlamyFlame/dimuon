#pragma once
#include <string>
#include <stdexcept>

// `noovl` = the r17663 NO-OVERLAY diagnostic sample: Pythia pp collisions reconstructed
// with the SAME PbPb23-conditions pass as the HIJING overlay (Athena 24.0.58,
// OFLCOND-MC23-SDR-RUN3-05, ConditionsRunNumber=460000, L1 MC_HI_run3_v1) but with
// `Digitization.PileUp=False` and no HIJING input. Single pT-hat slice (pTH8_14, DSID
// 802781). It exists ONLY to separate "conditions/L1 configuration" from "HIJING
// occupancy" as the cause of the forward-endcap L1 anomaly
// (mc_trigger_efficiency.md R8; round-4 outcome tree). It is NOT an overlay: no HIJING,
// no centrality, no overlay Extras.
enum class FullSimSampleType { pp, hijing, zmumu, data, noovl };

inline std::string FullSimSampleSuffix(FullSimSampleType t) {
    switch (t) {
    case FullSimSampleType::pp:     return "";
    case FullSimSampleType::hijing: return "_hijing";
    case FullSimSampleType::zmumu:  return "_zmumu";
    case FullSimSampleType::data:   return "_data";
    case FullSimSampleType::noovl:  return "_noovl";
    }
    throw std::runtime_error("FullSimSampleSuffix: unknown type");
}

// HIJING-OVERLAY CONDITIONS YEAR. Two overlay TEST productions exist on disk, in two directories:
//   24 (DEFAULT) -> pythia_fullsim_hijing_overlay_test_sample/         Pb+Pb 2024 conditions
//                   (r17864; docs/tracking/hijing_overlay_pbpb24_test_sample_skim.md). The
//                   forward-looking default: the FULL overlay production will be 2024 conditions.
//   23           -> pythia_fullsim_hijing_overlay_test_sample_pbpb23/  Pb+Pb 2023 conditions
//                   (r17618 / r17662; the sample every overlay result up to 2026-09 was made on).
// The year is an explicit knob, never an implicit default that changed underneath the code: every
// overlay product carries it in its label (hijing_overlay_pbpb<yy>) and the MC trigger-efficiency
// chain exposes it as `overlay_year` (dr_correction_sample_cfg.h). Sanity-check codes that study
// the r17618 / r17662 / r17663 productions name the _pbpb23 directory explicitly.
// Ignored for every other sample type.
inline void FullSimCheckPbPbYear(int pbpb_year) {
    if (pbpb_year != 23 && pbpb_year != 24)
        throw std::runtime_error("FullSimSampleType: HIJING-overlay conditions year must be 23 or 24, got "
                                 + std::to_string(pbpb_year));
}

// Input directory. `is_test_sample` selects the small TEST sample vs the FULL production.
// It is the SAME switch that selects the isospin treatment (FullSimSampleUsesFourBeams below)
// -- deliberately one flag, so the input path and the isospin weight can never disagree.
inline std::string FullSimSampleInputDir(FullSimSampleType t, bool is_test_sample, int pbpb_year = 24) {
    switch (t) {
    case FullSimSampleType::pp:
        return is_test_sample
            ? "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_test_sample/"
            : "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_full_sample/";
    case FullSimSampleType::hijing:
        FullSimCheckPbPbYear(pbpb_year);
        if (!is_test_sample)
            return "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_hijing_overlay_full_sample/";
        return pbpb_year == 24
            ? "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_hijing_overlay_test_sample/"
            : "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_hijing_overlay_test_sample_pbpb23/";
    case FullSimSampleType::zmumu:
        return "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_zmumu_overlay_test_sample/";
    case FullSimSampleType::data:
        return "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_data_overlay_test_sample/";
    case FullSimSampleType::noovl:
        // Only one production exists (10 k events, single slice) -- the same directory
        // either way, so the isTestSample switch cannot point it anywhere wrong.
        return "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_no_overlay_test_sample/";
    }
    throw std::runtime_error("FullSimSampleInputDir: unknown type");
}

// "Overlay" = a minimum-bias/underlying event was overlaid on the Pythia signal, so the
// event carries a heavy-ion environment (centrality, FCal, HIJING truth). `noovl` is
// reconstructed with PbPb conditions but has NO overlaid event -> NOT an overlay.
inline bool FullSimSampleIsOverlay(FullSimSampleType t) {
    return t != FullSimSampleType::pp && t != FullSimSampleType::noovl;
}

// ISOSPIN CONTENT OF THE SIMULATED COLLISION SYSTEM -- the single rule, in one place.
//
//   pp CONDITIONS fullsim simulates pp COLLISIONS
//       -> ONE beam (pp), isospin weight 1. There is nothing to isospin-average.
//   PbPb CONDITIONS HIJING overlay simulates Pb+Pb, whose nucleons are a p/n mix
//       -> FOUR beams {pp,pn,np,nn}, combined with the Pb ratio 4:6:6:9 (Z=82, N=126).
//
// The two TEST samples are exceptions in OPPOSITE directions -- they are what happens to
// exist on disk, not what the physics wants:
//   pp24 TEST sample        : produced with 4 isospin beams BY MISTAKE  -> 4 beams
//   HIJING overlay TEST smp : only the pp beam was produced             -> 1 beam
// so the exception is exactly "is this a test sample?", and the rule collapses to a XOR:
inline bool FullSimSampleUsesFourBeams(FullSimSampleType t, bool is_test_sample) {
    // r17663: pp COLLISIONS, and only the pp-beam DSID (802781) was ever produced ->
    // one beam, isospin weight 1, regardless of the test/full switch. (It is a
    // single-slice diagnostic anyway: one global weight cancels in every ratio.)
    if (t == FullSimSampleType::noovl) return false;
    return FullSimSampleIsOverlay(t) != is_test_sample;
}

// NTUP file tag.  Frozen: it is baked into the skimmed NTUP file names on disk and
// into the grid output-dataset names (SkimCode/run_pythia_fullsim_HIJING_overlay/
// grid_sub*.sh).  "PP24" here is a legacy misnomer for the HIJING overlay -- the
// overlay is Pb+Pb (see FullSimSampleLabel) -- but renaming it would orphan the
// existing NTUPs and grid datasets.
inline std::string FullSimSampleFileTag(FullSimSampleType t) {
    switch (t) {
    case FullSimSampleType::pp:     return "FullSimPP24";
    case FullSimSampleType::hijing: return "FullSimHIJINGOverlayPP24";
    case FullSimSampleType::zmumu:  return "FullSimZmumuOverlayPP24";
    case FullSimSampleType::data:   return "FullSimDataOverlayPP24";
    // r17663: skimmed with the overlay run mode, so it inherits the overlay basename,
    // distinguished by the _r17663 suffix (grid_sub_r17663_nooverlay.sh).
    case FullSimSampleType::noovl:  return "FullSimHIJINGOverlayPP24_r17663";
    }
    throw std::runtime_error("FullSimSampleFileTag: unknown type");
}

// Output-file / plot-directory label.  The HIJING overlay simulates Pb+Pb collisions,
// never pp: the label carries the Pb+Pb CONDITIONS year (see FullSimSampleInputDir) --
// "hijing_overlay_pbpb23" for the r17618 / r17662 test sample (ConditionsRunNumber=460000),
// "hijing_overlay_pbpb24" for the r17864 test sample and the coming full production.
inline std::string FullSimSampleLabel(FullSimSampleType t, int pbpb_year = 24) {
    switch (t) {
    case FullSimSampleType::pp:     return "pp24";
    case FullSimSampleType::hijing:
        FullSimCheckPbPbYear(pbpb_year);
        return "hijing_overlay_pbpb" + std::to_string(pbpb_year);
    case FullSimSampleType::zmumu:  return "zmumu_overlay_pp24";
    case FullSimSampleType::data:   return "data_overlay_pp24";
    case FullSimSampleType::noovl:  return "r17663_no_overlay";
    }
    throw std::runtime_error("FullSimSampleLabel: unknown type");
}

inline std::string FullSimSamplePlotDir(FullSimSampleType t, int pbpb_year = 24) {
    // Identical to the label today; kept as its own function because plot directories and file
    // labels are allowed to diverge (they did for the pp samples' _full suffix).
    return FullSimSampleLabel(t, pbpb_year);
}

// =============================================================================================
// PER-SAMPLE DIRECTORY LAYOUT (docs/tracking/fullsim_sample_dir_layout.md, user decision
// 2026-09-16). THE single source of truth for where a product lives inside a sample directory;
// nothing composes one of these subdirectory names anywhere else (the shell twin is
// Analysis/pipelines/fullsim_sample_layout.sh -- keep the two in step).
//
//   <sample>/                          raw NTUP, ami_info/, merging-record.txt, r-tag records,
//                                      grid-monitor state -- the SAMPLE itself (SkimCode-owned)
//     muon_pairs_*.root                ntuple-processing output        } FLAT, exactly like the
//     hists_pythia_ntuple_processing_* ntuple-processing histograms    } data directories
//     histograms_pythia_fullsim_*.root RDF hist-filling output         } (dimuon_data/pp_2024/)
//     mc_trig_eff/hists/               mc_trig_eff_hists_<label>*.root      FillMCTrigEffHists
//     mc_trig_eff/singles_fits/        single_mu_effcy_pT_fit_mc_<label>*   FitMCSinglesEffcy
//     mc_trig_eff/dr_correction/       dr_correction_plateaus_* / _fits_*   plot_mc_trig_eff,
//                                                                           fit_dr_corrections
//     mc_trig_eff/pair_eff/            pair_trig_eff_<label>*.root          FillMCTrigEffPairEff
//     mc_trig_eff/closure/             mc_trig_eff_closure_<label>*.root    FillMCTrigEffClosure
//     reco_eff/                        pair_reco_eff_<label>.root           build_pp24_fullsim_pair_reco_eff
//     plots/                           per-sample plots (reco-eff, det-response, kn tables, ...)
//     backup/                          *.bak_<date> copies -- NEVER read by code
//     logs/                            pipeline / farm logs
//
// MC-ONLY derived products go into subtrees; the ntuple-processing and hist-filling outputs stay
// flat because that is the layout of the DATA directories and one convention is easier to keep
// than two. Producers create their subdirectory (mkdir -p / gSystem->mkdir(..., kTRUE)).
inline std::string FullSimMCTrigEffDir(const std::string& sample_dir)           { return sample_dir + "mc_trig_eff/"; }
inline std::string FullSimMCTrigEffHistsDir(const std::string& sample_dir)      { return sample_dir + "mc_trig_eff/hists/"; }
inline std::string FullSimMCTrigEffSinglesFitDir(const std::string& sample_dir) { return sample_dir + "mc_trig_eff/singles_fits/"; }
inline std::string FullSimMCTrigEffDrCorrDir(const std::string& sample_dir)     { return sample_dir + "mc_trig_eff/dr_correction/"; }
inline std::string FullSimMCTrigEffPairEffDir(const std::string& sample_dir)    { return sample_dir + "mc_trig_eff/pair_eff/"; }
inline std::string FullSimMCTrigEffClosureDir(const std::string& sample_dir)    { return sample_dir + "mc_trig_eff/closure/"; }
inline std::string FullSimRecoEffDir(const std::string& sample_dir)             { return sample_dir + "reco_eff/"; }
inline std::string FullSimPlotsDir(const std::string& sample_dir)               { return sample_dir + "plots/"; }
inline std::string FullSimBackupDir(const std::string& sample_dir)              { return sample_dir + "backup/"; }
inline std::string FullSimLogsDir(const std::string& sample_dir)                { return sample_dir + "logs/"; }

// The MC single-muon mu4 turn-on fit file (FitMCSinglesEffcy). The sample LABEL is part of the
// basename (mc_trigger_efficiency.md RW 3c): before 2026-09-16 the basename was identical across
// pp24 / overlay / noovl and only the directory told them apart, so a wrong directory would have
// silently overwritten a sibling. `corr_suffix` = "" | "_corrected" | "_corrected_sfclosure";
// `wp_suffix` = "" (Tight) | "_medium_wp".
inline std::string FullSimMCSinglesFitFile(const std::string& sample_dir, const std::string& label,
                                           const std::string& corr_suffix, const std::string& wp_suffix) {
    return FullSimMCTrigEffSinglesFitDir(sample_dir) + "single_mu_effcy_pT_fit_mc_" + label
         + corr_suffix + wp_suffix + ".root";
}
