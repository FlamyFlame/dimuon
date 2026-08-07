// dr_correction_sample_cfg.h
//
// SINGLE SOURCE OF TRUTH for the MC trigger-efficiency sample identities used by the whole
// DeltaR-correction chain (mc_trigger_efficiency.md §3.3 Step 3 / §3.4 Step 4):
//
//   plot_mc_trig_eff.cxx        measures eps_dR / eps_single and WRITES the per-cell large-dR
//                               plateau ROOT file
//   fit_dr_corrections.cxx      READS that plateau file, guards it, normalizes and FITS
//   plot_dr_correction_fits.cxx re-opens the fit file in a FRESH process and draws it
//
// The three stages must agree BYTE-FOR-BYTE on the input directory, the file label and the plot
// root, or the fit stage silently fits one sample's histograms against another sample's plateaus.
// Hence one table, included by all three.
//
// `is_full_sample` is the sample IDENTITY flag that drives the plateau guard (a FULL sample must
// have |plateau-1| <= 0.1 in every (pair pT, pair eta) cell; a 10 000-event TEST sample is exempt
// because its high-pair-pT cells are noise-dominated). It is a property of the PRODUCTION, never
// a hand-maintained list of measured numbers.

#ifndef DR_CORRECTION_SAMPLE_CFG_H
#define DR_CORRECTION_SAMPLE_CFG_H

#include <stdexcept>
#include <string>
#include "../../../Utilities/MCTrigEffPairPtBinning.h"

struct DrCorrSample {
    std::string key;             // sample token used on the command line
    std::string mc_dir;          // directory holding mc_trig_eff_hists_* and the fit outputs
    std::string mc_label;        // file label inside those names
    std::string out_base;        // plots root for this sample
    std::string sample_text;     // canvas headline (without the working point)
    std::string eps_dr_text;     // Step-3 symbol (2mu4 product vs mu4 cross term)
    bool        is_full_sample;  // true = FULL production -> plateau guard is ENFORCED
};

// UNMEASURABLE-CELL SCREEN -- shared by the fit stage (which refuses to fit such a cell) and the
// plot stage (which refuses to draw one). eps_dR is a ratio that must tend to 1 at large dR, so a
// plateau far below 1 is not a normalization offset but an empty cell; dividing by 0.01 inflates
// the curve 100x rather than normalizing it. Likewise a plateau whose error exceeds itself carries
// no information. User policy: such cells are NOT plotted and their count IS reported.
inline bool DrCorrPlateauUsable(double plateau, double err) {
    constexpr double kPlateauMinUsable = 0.5;
    return plateau > 0. && plateau >= kPlateauMinUsable && !(err >= plateau);
}

inline std::string DrCorrOutTag(bool use_tight_wp) {
    return MCTrigEffPairPt::FileSuffix() + std::string(use_tight_wp ? "" : "_medium");
}

inline DrCorrSample GetDrCorrSample(const std::string& key, bool use_tight_wp)
{
    DrCorrSample s;
    s.key = key;
    if (key == "pp") {
        // pp24 fullsim TEST sample (24 x 10k). Superseded by "pp_full"; kept so the older
        // outputs remain reproducible.
        s.mc_dir         = "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_test_sample/";
        s.mc_label       = "pp24";
        s.out_base       = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/"
                           "pp_trigger_efficiency/mc_based" + DrCorrOutTag(use_tight_wp) + "/";
        // Canvas headlines carry the PHYSICAL sample identity only (beam, energy, run
        // conditions) -- production bookkeeping ("fullsim", "FULL sample", r-tags) means
        // nothing to a physics audience and is kept in this doc / the file names instead.
        // The only exception is the r17663 diagnostic below, whose whole point IS the tag.
        s.sample_text    = "Pythia8 pp, #sqrt{s} = 5.36 TeV (2024 conditions)";
        s.eps_dr_text    = "#varepsilon_{#DeltaR}^{2mu4}";
        s.is_full_sample = false;
    } else if (key == "pp_full") {
        s.mc_dir         = "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_full_sample/";
        s.mc_label       = "pp24_full";
        s.out_base       = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/"
                           "pp_trigger_efficiency/mc_based" + DrCorrOutTag(use_tight_wp) + "/";
        // Same PHYSICS as the "pp" test sample above, hence the same headline: the test/full
        // distinction is a production fact, not a property the reader of the figure can act on.
        s.sample_text    = "Pythia8 pp, #sqrt{s} = 5.36 TeV (2024 conditions)";
        s.eps_dr_text    = "#varepsilon_{#DeltaR}^{2mu4}";
        s.is_full_sample = true;    // the only FULL production so far
    } else if (key == "overlay") {
        s.mc_dir         = "/usatlas/u/yuhanguo/usatlasdata/"
                           "pythia_fullsim_hijing_overlay_test_sample/";
        s.mc_label       = "hijing_overlay_pbpb23";
        s.out_base       = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/"
                           "pbpb_trigger_efficiency/mc_based" + DrCorrOutTag(use_tight_wp) + "/";
        s.sample_text    = "Pythia8 + HIJING overlay, Pb+Pb #sqrt{s_{NN}} = 5.36 TeV, "
                           "0-5% (2023 conditions)";
        s.eps_dr_text    = "#varepsilon_{#DeltaR}^{cross}";
        s.is_full_sample = false;   // 10 000-event TEST sample
    } else if (key == "noovl") {
        s.mc_dir         = "/usatlas/u/yuhanguo/usatlasdata/"
                           "pythia_fullsim_no_overlay_test_sample/";
        s.mc_label       = "r17663_no_overlay";
        s.out_base       = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/"
                           "r17663_no_overlay_trigger_efficiency/mc_based" + DrCorrOutTag(use_tight_wp) + "/";
        // The ONE headline that keeps a production tag (user decision 2026-08-04): this plot set
        // exists to compare reconstruction CONFIGURATIONS, so the tag is the subject of the
        // figure rather than bookkeeping. The series labels are still physical.
        s.sample_text    = "Pythia8 pp, #sqrt{s} = 5.36 TeV, Pb+Pb reconstruction conditions, "
                           "no overlay (r17663)";
        s.eps_dr_text    = "#varepsilon_{#DeltaR}^{2mu4}";
        s.is_full_sample = false;   // 10 000-event TEST sample
    } else {
        throw std::runtime_error("GetDrCorrSample: sample must be 'pp', 'pp_full', 'overlay' or "
                                 "'noovl', got '" + key + "'");
    }
    return s;
}

// WP-dependent file suffix (Tight nominal = unsuffixed; Medium = _medium_wp), used for every
// file name in the chain. Registry: Analysis/docs/muon_wp_registry.md.
inline std::string DrCorrWpSuffix(bool use_tight_wp) { return use_tight_wp ? "" : "_medium_wp"; }

// VARIANT LAYOUT (user, 2026-08-05): every non-nominal variant lives in its OWN TOP-LEVEL plot
// directory -- `mc_based`, `mc_based_medium`, `mc_based_pt4bin`, `mc_based_pt4bin_medium` --
// instead of a `medium/` or `pt4bin/` subdirectory sprinkled through every step directory. One
// variant = one self-contained tree, so a whole plot set can be diffed, copied or deleted as a
// unit and two variants can never interleave inside one directory (which is exactly how the
// 4-bin fit plots ended up mixed with the 8-bin ones on 2026-08-04).
// Retired: the per-directory WP subdirectory. The working point is now in the top-level base
// (DrCorrOutTag), so this must return "" -- keeping it non-empty would nest medium/ INSIDE
// mc_based_medium/.
inline std::string DrCorrWpDir(bool /*use_tight_wp*/) { return ""; }

// The per-(pair pT, pair eta) large-dR plateau file: WRITTEN by plot_mc_trig_eff.cxx (the step
// that measures it) and READ by fit_dr_corrections.cxx. It is the machine-readable source --
// the .txt tables beside the plots are for humans and must never be parsed by code.
inline std::string DrCorrPlateauFile(const DrCorrSample& s, bool use_tight_wp)
{
    // The pair-pT-binning token keeps a 4-bin comparison run from overwriting the nominal
    // plateau map -- and, more importantly, keeps the fit stage's cell-count consistency guard
    // from silently pairing an 8-bin plateau file with 4-bin histograms.
    return s.mc_dir + "dr_correction_plateaus_" + s.mc_label + DrCorrWpSuffix(use_tight_wp)
         + MCTrigEffPairPt::FileSuffix() + ".root";
}

// Fit output, one file per (sample, WP, step, method).
inline std::string DrCorrFitFile(const DrCorrSample& s, bool use_tight_wp, int step,
                                 const std::string& method)
{
    return s.mc_dir + "dr_correction_fits_" + s.mc_label + DrCorrWpSuffix(use_tight_wp)
         + MCTrigEffPairPt::FileSuffix()
         + "_step" + std::to_string(step) + "_" + method + ".root";
}

#endif // DR_CORRECTION_SAMPLE_CFG_H
