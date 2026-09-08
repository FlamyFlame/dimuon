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

// ---------------------------------------------------------------------------------------------
// SIGN SERIES of the dR correction (added 2026-08-11).
//   ""   sign-integrated -- the NOMINAL correction the analysis applies
//   "ss" same sign        (repo convention: muon_pair_tree_sign1 = same sign)
//   "os" opposite sign    (muon_pair_tree_sign2 = opposite sign)
// The token enters (a) the histogram-name prefix, (b) the plateau-file key tag, (c) the fit-FILE
// name and (d) the report-file names. It deliberately does NOT enter the keys INSIDE a fit file:
// the sign is already in the file name, so all three series carry identical internal keys and
// every consumer reads them the same way.
// Canvas/report wording is always the SPELLED-OUT physics ("same sign" / "opposite sign"); bare
// SS/OS and sign1/sign2 are forbidden on a figure (.claude/conventions/atlas-plotting.md).
inline std::string DrCorrSignText(const std::string& sign)
{
    if (sign.empty())  return "";
    if (sign == "ss")  return "same sign";
    if (sign == "os")  return "opposite sign";
    throw std::runtime_error("DrCorrSignText: sign must be \"\" (sign-integrated), \"ss\" or "
                             "\"os\", got '" + sign + "'");
}

// Report-file token, e.g. fit_report.txt -> fit_report_same_sign.txt. Physically named so a
// same-sign report can never be mistaken for -- or clobber -- the opposite-sign or the nominal one.
inline std::string DrCorrSignFileTag(const std::string& sign)
{
    if (sign.empty()) return "";
    std::string t = DrCorrSignText(sign);
    for (char& c : t) if (c == ' ') c = '_';
    return "_" + t;
}

// ---------------------------------------------------------------------------------------------
// PLATEAU MODE of the dR fit (added 2026-08-11).
//   "corr"   (or "") -- NOMINAL: each cell's eps_dR is divided by ITS OWN large-dR plateau
//                       (dR in [2, 3.5], MCTrigEffPlateauWindow.h) and the normalized curve is
//                       fitted with a shape that tends to 1.
//   "nocorr"         -- the RAW, un-normalized eps_dR is fitted with a FREE additive baseline C.
//                       The baseline is then determined by the dR < 1 data itself and the
//                       [2, 3.5] window NEVER enters the fit. Motivation (user, 2026-08-11): the
//                       inverse-weighted dR distribution shows structure out to large dR --
//                       worst in the pair-eta bins that enclose the detector gap -- so a plateau
//                       measured far away may not be the right baseline for the small-dR region
//                       the correction is actually about.
// The mode is the TOP level of the plot tree (step<N>_dr_fit/<mode dir>/<method>/<sign mode>/)
// and a token in the fit file name; both are built HERE so the fit stage, the plot stage and the
// driver script cannot disagree about them.
//   "nocorr_ptmerge" -- the SAME raw fit as "nocorr", with the LAST TWO pair-pT BINS MERGED into
//                       one cell (added 2026-08-17, user request). It exists because the top two
//                       cells of the 8-bin log axis, p_T^pair in [72.1,104) and [104,150) GeV, run
//                       past where pp Pythia has yield: they hold the plateau-guard failures and
//                       their dR fits are noise-dominated. Merging them buys back statistics in
//                       exactly one place and leaves every other cell untouched.
//                       IT IS NOT A NEW BINNING (.claude/CLAUDE.md 'Binnings'): the FILLED
//                       histograms keep ParamsSet::pair_pt_coarse_bins, the two source bins are
//                       simply PROJECTED TOGETHER at the fit/plot stage -- numerically identical
//                       to having filled a 7-bin axis -- and the variant is opt-in and suffixed so
//                       it can neither overwrite nor be mistaken for the other two modes.
//                       Defined for the 8-bin nominal axis ONLY; asking for it with
//                       MCTRIGEFF_PAIRPT_4BIN set THROWS (dr_correction_cell_groups.h).
//   "nocorr_etamerge" -- the SAME raw fit as "nocorr", with the 9 pair-eta bins MERGED into THREE
//                       SIGN-INDEPENDENT |eta^pair| BINS: |eta| < 1.0 (barrel), 1.0 <= |eta| < 2.0
//                       and 2.0 <= |eta| < 2.4 (added 2026-08-24, user request; SUPERSEDED
//                       2026-09-03 -- the ORIGINAL grouping was signed: negative-eta endcap
//                       (-2.4,-1.0) / barrel (-1.0,1.0) / positive-eta endcap (1.0,2.4). Replaced
//                       because the dR correlation was found to barely depend on the SIGN of pair
//                       eta, while the |eta|>2-vs-<2 split inside the endcap is a much bigger
//                       effect than any negative/positive asymmetry --
//                       docs/tracking/mc_trigeff_dr_binning_approaches.md). It is the pair-ETA
//                       analogue of "nocorr_ptmerge": the 9-bin pair-eta grid is the
//                       CROSS-SECTION's presentation binning and was never chosen for eps_dR's
//                       statistics, so grouping it triples the pairs per fit while keeping the one
//                       distinction that is now physically motivated (forward vs less-forward, not
//                       positive vs negative). The barrel group is one contiguous source range;
//                       each forward group FOLDS its negative- and positive-eta source bins
//                       together (dr_correction_cell_groups.h).
//   "nocorr_etamerge_ptmerge" -- both merges at once (7 pair-pT x 3 pair-eta = 21 cells).
// Like the pair-pT merge, NEITHER is a new binning: the filled histograms are untouched and the
// source bins are PROJECTED TOGETHER at the fit stage (dr_correction_cell_groups.h).
inline std::string DrCorrPlateauModeDir(const std::string& mode)
{
    if (mode.empty() || mode == "corr")  return "plateau_corrected/";
    if (mode == "nocorr")                return "no_plateau_correction/";
    if (mode == "nocorr_ptmerge")        return "no_plateau_correction_last2ptbins_merged/";
    if (mode == "nocorr_etamerge")       return "no_plateau_correction_paireta_merged/";
    if (mode == "nocorr_etamerge_ptmerge")
        return "no_plateau_correction_paireta_merged_last2ptbins_merged/";
    throw std::runtime_error("DrCorrPlateauModeDir: plateau mode must be \"corr\", \"nocorr\", "
                             "\"nocorr_ptmerge\", \"nocorr_etamerge\" or "
                             "\"nocorr_etamerge_ptmerge\", got '" + mode + "'");
}

// THE TWO QUESTIONS every stage actually asks about a plateau mode. Ask these -- never compare the
// token to a literal: with three modes, `mode == "nocorr"` silently answers NO for the merged mode,
// which is also a no-plateau-correction one, and the plateau would be applied where it must not be.
inline bool DrCorrModeNoPlateau(const std::string& mode)
{
    DrCorrPlateauModeDir(mode);                       // validates the token
    return mode.rfind("nocorr", 0) == 0;              // every nocorr* mode
}

inline bool DrCorrModeMergeLastTwoPt(const std::string& mode)
{
    DrCorrPlateauModeDir(mode);                       // validates the token
    return mode == "nocorr_ptmerge" || mode == "nocorr_etamerge_ptmerge";
}

// Does this mode group the 9 filled pair-eta bins into the 3 sign-independent |eta^pair| bins?
inline bool DrCorrModeMergeEta(const std::string& mode)
{
    DrCorrPlateauModeDir(mode);                       // validates the token
    return mode == "nocorr_etamerge" || mode == "nocorr_etamerge_ptmerge";
}

// File-name token. EMPTY for the nominal mode, so every pre-existing fit file keeps its current
// name byte-for-byte and no consumer of the nominal correction has to be touched.
inline std::string DrCorrPlateauModeTag(const std::string& mode)
{
    if (mode.empty() || mode == "corr")       return "";
    if (mode == "nocorr")                     return "_nocorr";
    if (mode == "nocorr_ptmerge")             return "_nocorr_ptmerge";
    if (mode == "nocorr_etamerge")            return "_nocorr_etamerge";
    if (mode == "nocorr_etamerge_ptmerge")    return "_nocorr_etamerge_ptmerge";
    throw std::runtime_error("DrCorrPlateauModeTag: plateau mode must be \"corr\", \"nocorr\", "
                             "\"nocorr_ptmerge\", \"nocorr_etamerge\" or "
                             "\"nocorr_etamerge_ptmerge\", got '" + mode + "'");
}

// Fit output, one file per (sample, WP, step, method, sign, plateau mode). `sign` and
// `plateau_mode` both default to the NOMINAL choice so every pre-existing call site keeps its
// current file name byte-for-byte.
inline std::string DrCorrFitFile(const DrCorrSample& s, bool use_tight_wp, int step,
                                 const std::string& method, const std::string& sign = "",
                                 const std::string& plateau_mode = "")
{
    return s.mc_dir + "dr_correction_fits_" + s.mc_label + DrCorrWpSuffix(use_tight_wp)
         + MCTrigEffPairPt::FileSuffix()
         + "_step" + std::to_string(step) + "_" + method
         + (sign.empty() ? "" : "_" + sign) + DrCorrPlateauModeTag(plateau_mode) + ".root";
}

// ---------------------------------------------------------------------------------------------
// THE VARIANT THE pp24 CROSSX APPLICATION USES (user, 2026-08-17) -- TEMPORARY.
// Named here so the choice is made in ONE place and no consumer retypes it. It is the
// no-plateau-correction fit with the last two pair-pT bins merged, opposite-sign pairs, `expo`
// form; the correction is applied as eps_dR = f(dR)/C for dR < 1 and 1 above
// (dr_correction_apply.h, mc_trig_eff_closure.md 3.2). The other methods and sign series are still
// produced -- they are the comparison, not the deliverable.
// TEMPORARY because mc_trigger_efficiency.md R26 (both parametric forms fail in a minority of
// cells, and `usable` carries no chi2 term) is still OPEN.
inline const char* DrCorrCrossxMethod() { return "expo"; }
inline const char* DrCorrCrossxSign()   { return "os"; }
inline const char* DrCorrCrossxMode()   { return "nocorr_ptmerge"; }

#endif // DR_CORRECTION_SAMPLE_CFG_H
