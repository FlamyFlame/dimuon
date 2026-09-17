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
#include <vector>
#include "../../../Utilities/MCTrigEffPairPtBinning.h"
#include "../../../MuonObjectsParamsAndHelpers/FullSimSampleType.h"

struct DrCorrSample {
    std::string key;             // sample token used on the command line
    std::string sample_dir;      // the sample's ROOT directory (raw NTUP level). Products live in
                                 // its subtrees -- compose with the FullSimMCTrigEff*Dir /
                                 // DrCorr*File helpers, NEVER by appending a basename to this.
    std::string mc_label;        // file label inside the product names
    int         overlay_year;    // HIJING-overlay conditions year (24 default / 23); see below
    std::string plot_root;       // ".../<beam>_trigger_efficiency/" -- the sample's plot root
    std::string plot_leaf;       // "" or "pbpb<yy>/" -- the conditions-year leaf (overlay only)
    std::string out_base;        // plot_root + "mc_based<tag>/" + plot_leaf: the variant tree
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

// RUNAWAY-BASELINE SCREEN (user decision 2026-09-10). `DrCorrPlateauUsable` above bounds the
// baseline only FROM BELOW, because it was written against a COLLAPSING one ("dividing by 0.01
// inflates the curve 100x"). It is structurally blind to the opposite failure, which is real:
// in the crossx-consumed expo/OS/Tight series the cell pair-pT [52.2,74.2) x pair-eta [1.0,1.5)
// had lambda railed at its 3.0 limit, and MINUIT bought that nearly-straight-line fit by driving
// the free baseline to C = 3.9755 +- 1.1774 against a plateau the SAME cell measured as 0.8225.
// Since eps_dR = f/C, the delivered correction was 0.159 at dR=0 and 0.323 at dR=1 -- it never
// approached 1, which it must by construction -- and the pp24 weight 1/(eps1*eps2*eps_dR) was
// inflated 3.1-6.3x for every affected pair.
//
// C is the baseline the curve is normalized BY, so it must be consistent with the plateau the
// cell actually MEASURED. Measured over the 63 accepted cells of that series, C/plateau has
// median 1.003 and a 5th-95th percentile of 0.950-1.234; the pathological cell sits at 4.834 and
// the next-worst at 1.991 (a cell whose PLATEAU is a marginal 0.539, not a runaway C). A
// two-sided factor-1.5 window is therefore generous by a wide margin and still isolates both.
// A cell rejected here falls through to the polyu / interp / raw-bin tiers, which exist for
// exactly this and were being reached in 0 of 63 cells because the expo screen accepted almost
// everything.
//
// NOT a bound on the delivered correction f/C. That is a DIFFERENT question and remains open:
// 26 of those 63 cells deliver f/C < 0.5 somewhere in dR < 1 (down to 0.036), which the C screen
// cannot see because their C is perfectly normal. See the parent tracking doc.
inline constexpr double kDrCorrBaselineMaxRatio = 1.5;
inline bool DrCorrBaselineConsistent(double C, double measured_plateau) {
    if (!(measured_plateau > 0.)) return true;   // no measured plateau -> screen inapplicable, not failed
    const double r = C / measured_plateau;
    return r <= kDrCorrBaselineMaxRatio && r >= 1.0 / kDrCorrBaselineMaxRatio;
}

inline std::string DrCorrOutTag(bool use_tight_wp) {
    return MCTrigEffPairPt::FileSuffix() + std::string(use_tight_wp ? "" : "_medium");
}

// `overlay_year` (user, 2026-09-16): which HIJING-overlay TEST production the "overlay" key
// means -- 24 (DEFAULT, Pb+Pb 2024 conditions, pythia_fullsim_hijing_overlay_test_sample/) or
// 23 (Pb+Pb 2023 conditions, ..._pbpb23/, the sample of every overlay result up to 2026-09).
// Directory, label, headline and plot root all derive from it (FullSimSampleType.h), so the
// two productions can never share a file or a figure. Ignored for the pp / noovl keys.
inline std::string DrCorrSiblingTree(const DrCorrSample& s, const std::string& name);

inline DrCorrSample GetDrCorrSample(const std::string& key, bool use_tight_wp, int overlay_year = 24)
{
    DrCorrSample s;
    s.key = key;
    s.overlay_year = overlay_year;
    if (key == "pp") {
        // pp24 fullsim TEST sample (24 x 10k). Superseded by "pp_full"; kept so the older
        // outputs remain reproducible.
        s.sample_dir     = FullSimSampleInputDir(FullSimSampleType::pp, /*is_test_sample=*/true);
        s.mc_label       = "pp24";
        s.plot_root      = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/pp_trigger_efficiency/";
        // Canvas headlines carry the PHYSICAL sample identity only (beam, energy, run
        // conditions) -- production bookkeeping ("fullsim", "FULL sample", r-tags) means
        // nothing to a physics audience and is kept in this doc / the file names instead.
        // The only exception is the r17663 diagnostic below, whose whole point IS the tag.
        s.sample_text    = "Pythia8 pp, #sqrt{s} = 5.36 TeV (2024 conditions)";
        s.eps_dr_text    = "#varepsilon_{#DeltaR}^{2mu4}";
        s.is_full_sample = false;
    } else if (key == "pp_full") {
        s.sample_dir     = FullSimSampleInputDir(FullSimSampleType::pp, /*is_test_sample=*/false);
        s.mc_label       = "pp24_full";
        s.plot_root      = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/pp_trigger_efficiency/";
        // Same PHYSICS as the "pp" test sample above, hence the same headline: the test/full
        // distinction is a production fact, not a property the reader of the figure can act on.
        s.sample_text    = "Pythia8 pp, #sqrt{s} = 5.36 TeV (2024 conditions)";
        s.eps_dr_text    = "#varepsilon_{#DeltaR}^{2mu4}";
        s.is_full_sample = true;    // the only FULL production so far
    } else if (key == "overlay") {
        FullSimCheckPbPbYear(overlay_year);
        s.sample_dir     = FullSimSampleInputDir(FullSimSampleType::hijing, /*is_test_sample=*/true,
                                                 overlay_year);
        s.mc_label       = FullSimSampleLabel(FullSimSampleType::hijing, overlay_year);
        // The conditions year is the LEAF of the plot root, exactly as the data-side trees put it
        // (pbpb_trigger_efficiency/mu4/no_corr/pbpb<yy>), so the two overlay productions' plot
        // sets sit side by side under one variant tree and never interleave.
        s.plot_root      = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/pbpb_trigger_efficiency/";
        s.plot_leaf      = "pbpb" + std::to_string(overlay_year) + "/";
        s.sample_text    = "Pythia8 + HIJING overlay, Pb+Pb #sqrt{s_{NN}} = 5.36 TeV, "
                           "0-5% (20" + std::to_string(overlay_year) + " conditions)";
        s.eps_dr_text    = "#varepsilon_{#DeltaR}^{cross}";
        s.is_full_sample = false;   // 10 000-event TEST sample
    } else if (key == "noovl") {
        s.sample_dir     = FullSimSampleInputDir(FullSimSampleType::noovl, /*is_test_sample=*/true);
        s.mc_label       = FullSimSampleLabel(FullSimSampleType::noovl);
        s.plot_root      = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/r17663_no_overlay_trigger_efficiency/";
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
    // VARIANT LAYOUT (below): the WP / binning variant is a top-level tree under the plot root;
    // the overlay's conditions year is the leaf INSIDE it (D8 of fullsim_sample_dir_layout.md).
    s.out_base = DrCorrSiblingTree(s, "mc_based" + DrCorrOutTag(use_tight_wp));
    return s;
}

// A SIBLING top-level tree of the variant tree, e.g. "mc_statistics_pt8bins" next to
// "mc_based": plot_root + "<name>/" + plot_leaf. Every product set that sits BESIDE the mc_based
// tree composes its directory here, so the year leaf can never be dropped or doubled.
inline std::string DrCorrSiblingTree(const DrCorrSample& s, const std::string& name)
{
    return s.plot_root + name + "/" + s.plot_leaf;
}

// WP-dependent file suffix (Tight nominal = unsuffixed; Medium = _medium_wp), used for every
// file name in the chain. Registry: Analysis/docs/muon_wp_registry.md.
inline std::string DrCorrWpSuffix(bool use_tight_wp) { return use_tight_wp ? "" : "_medium_wp"; }

// ---------------------------------------------------------------------------------------------
// PRODUCT FILES of the chain, each in its subtree of the sample directory (FullSimSampleType.h
// "PER-SAMPLE DIRECTORY LAYOUT"). Every stage composes its input and output names HERE, so a
// producer and its consumers cannot disagree about where a file lives.
//
// The Step-1..4 histogram file written by FillMCTrigEffHists. `variant` is everything that
// follows the WP token in the basename, e.g. "" (Steps 1-2), "_step3", "_step4", "_sanity",
// "_corrected", "_corrected_sfclosure_step3", "_nogapcut", or the 4-bin token followed by one
// of those -- the caller keeps composing it, only the directory and the prefix are fixed here.
inline std::string DrCorrHistFile(const DrCorrSample& s, bool use_tight_wp,
                                  const std::string& variant = "")
{
    return FullSimMCTrigEffHistsDir(s.sample_dir) + "mc_trig_eff_hists_" + s.mc_label
         + DrCorrWpSuffix(use_tight_wp) + variant + ".root";
}

// The MC single-muon turn-on fit file (FitMCSinglesEffcy), `corr_suffix` = "" | "_corrected" |
// "_corrected_sfclosure". Carries the sample label since 2026-09-16 (RW 3c).
inline std::string DrCorrSinglesFitFile(const DrCorrSample& s, bool use_tight_wp,
                                        const std::string& corr_suffix = "")
{
    return FullSimMCSinglesFitFile(s.sample_dir, s.mc_label, corr_suffix, DrCorrWpSuffix(use_tight_wp));
}

// The MC closure output (FillMCTrigEffClosure); `variant` = the plateau-mode / binning tokens
// that follow the WP token, composed by the caller exactly as before.
inline std::string DrCorrClosureFile(const DrCorrSample& s, bool use_tight_wp,
                                     const std::string& variant = "")
{
    return FullSimMCTrigEffClosureDir(s.sample_dir) + "mc_trig_eff_closure_" + s.mc_label
         + DrCorrWpSuffix(use_tight_wp) + variant + ".root";
}

// The single-value pair trigger efficiency map (FillMCTrigEffPairEff) is named by
// PairTrigEff::FileName(s.sample_dir, s.mc_label, wp) in Utilities/PairTrigEffEvaluator.h.

// DATA tag-and-probe REFERENCE of a sample (Step-1 comparison, corrected-MC study): the data
// taken under the SAME conditions the sample simulates -- pp24 data for every pp-collision sample
// (pp, pp_full, and the r17663 no-overlay diagnostic, which simulates pp collisions), Pb+Pb data
// of the overlay's CONDITIONS YEAR for the HIJING overlay (23 -> pbpb_2023, the D2 choice for the
// r17618 sample; 24 -> pbpb_2024, like-for-like for the r17864 sample). Returned as the data
// directory + the year text, so a caller composes the file it needs without retyping the year.
inline std::string DrCorrDataRefDir(const DrCorrSample& s)
{
    if (s.key == "overlay")
        return "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_20" + std::to_string(s.overlay_year) + "/";
    return "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/";
}
inline std::string DrCorrDataRefTag(const DrCorrSample& s)     // file-name stem token
{
    if (s.key == "overlay") return "pbpb_20" + std::to_string(s.overlay_year);
    return "pp_2024";
}
inline std::string DrCorrDataRefText(const DrCorrSample& s)    // legend / headline
{
    if (s.key == "overlay") return "Pb+Pb 20" + std::to_string(s.overlay_year) + " data, 0-5%";
    return "pp 2024 data";
}

// The pp24-fullsim PAIR reconstruction efficiency map (build_pp24_fullsim_pair_reco_eff.C),
// consumed by the pp24 crossx RDF stage. Not a trigger product, but it lives in the same sample
// directory and is named from the same label.
inline std::string DrCorrPairRecoEffFile(const DrCorrSample& s)
{
    return FullSimRecoEffDir(s.sample_dir) + "pair_reco_eff_" + s.mc_label + ".root";
}

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
    return FullSimMCTrigEffDrCorrDir(s.sample_dir) + "dr_correction_plateaus_" + s.mc_label
         + DrCorrWpSuffix(use_tight_wp) + MCTrigEffPairPt::FileSuffix() + ".root";
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
//                       and 2.0 <= |eta| < 2.2 (added 2026-08-24, user request; SUPERSEDED
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
    return FullSimMCTrigEffDrCorrDir(s.sample_dir) + "dr_correction_fits_" + s.mc_label
         + DrCorrWpSuffix(use_tight_wp) + MCTrigEffPairPt::FileSuffix()
         + "_step" + std::to_string(step) + "_" + method
         + (sign.empty() ? "" : "_" + sign) + DrCorrPlateauModeTag(plateau_mode) + ".root";
}

// ---------------------------------------------------------------------------------------------
// THE pp24 CROSSX APPLICATION OF THE PAIR TRIGGER EFFICIENCY (user, 2026-09-17;
// docs/tracking/pp24_trig_eff_hybrid_application.md Physics Procedure §2-3). Named HERE so the
// choice is made in ONE place and no consumer retypes it. Every edge value quoted in the comments
// of this block and of the two evaluator headers is ILLUSTRATIVE of the 2026-09 axis
// (ParamsSet::pair_pt_coarse_bins, 8 log bins 9 -> 150 GeV): the code reads the axes and prints
// the ranges it resolved at load time; the comments are not a second source. Consumed by
// Utilities/DrCorrectionCrossxEvaluator.h (region A) and Utilities/PairTrigEffCrossxEvaluator.h
// (the hybrid). The other methods, sign series and plateau modes are still produced by the fit
// stage -- they are the comparison, not the deliverable.
//
// The weight is a HYBRID over the canonical coarse pair-pT axis (ParamsSet::pair_pt_coarse_bins,
// N = 8 bins):
//   REGION A, bins 1..N-2 ([9, 74.24) GeV):  eps^nc_data(1) * eps^nc_data(2) * eps_dR(dR; cell)
//       eps_dR from the NO-PLATEAU-CORRECTION fit (free baseline C, eps_dR = f/C for dR < 1),
//       OPPOSITE-sign series, on the 3-group |eta^pair| fold (DrCorrCrossxMode below). The
//       primary fit form per cell is DrCorrCrossxMethod() ("expo"), EXCEPT the cells named by
//       DrCorrCrossxPolyPtBins() x the LAST |eta| group, whose primary is the constrained
//       polynomial (user choice on the 2026-09-10 fit figures: the polynomial describes the
//       points where expo rails at p = 8). The ONLY fallback is DrCorrCrossxFallbackMethod()
//       ("interp"); there is NO raw-bin tier any more, and a cell with neither its primary nor
//       the interpolation accepted THROWS at load time (user, 2026-09-17).
//   REGION B, bins N-1..N ([74.24, 150) GeV): eps^pair_MC(cell) * SF(1) * SF(2)
//       the SINGLE-VALUE pair 2mu4 efficiency (Utilities/PairTrigEffEvaluator.h, PURE form) in
//       the SIGNAL mass window, opposite sign, on the UN-MERGED cells, times the product of the
//       two single-muon DATA/MC scale factors SF = eps^nc_data / eps_MC. It is a PAIR-level
//       efficiency, NOT a dR correction: it must never be multiplied by the single-muon
//       efficiencies themselves (docs/tracking/mc_trigeff_single_value_pair_eff.md D4).
//
// Previous application (2026-08-17 .. 2026-09-17): `nocorr_ptmerge` in every bin, expo -> polyu
// -> raw-bin placeholder. Retired because the raw tier delivered eps_dR down to 0.044 and the
// polynomial tier was never reached (mc_trigger_efficiency.md R35).
inline const char* DrCorrCrossxMethod()         { return "expo"; }            // region-A primary
inline const char* DrCorrCrossxPolyMethod()     { return "polyu_fixedRp"; }   // primary in the poly cells
inline const char* DrCorrCrossxFallbackMethod() { return "interp"; }          // the ONLY fallback
inline const char* DrCorrCrossxSign()           { return "os"; }
inline const char* DrCorrCrossxMode()           { return "nocorr_etamerge"; }
// The region-A cells whose PRIMARY is the polynomial: pair-pT bins 2, 3, 4 (1-based on
// ParamsSet::pair_pt_coarse_bins -- [12.8,18.2), [18.2,25.8), [25.8,36.7) GeV on the 9 -> 150
// axis) x the LAST |eta^pair| group ([2.0, 2.2), the most forward). Named by INDEX so the
// choice follows the canonical axes; the physical ranges are printed from the axes at load time.
inline const std::vector<int>& DrCorrCrossxPolyPtBins()
{
    static const std::vector<int> b = {2, 3, 4};
    return b;
}
// Region B: the single-value pair efficiency's (sign, mass window, cell mode).
inline const char* PairTrigEffCrossxSign()   { return "os"; }
inline const char* PairTrigEffCrossxWindow() { return "sig"; }
inline const char* PairTrigEffCrossxMode()   { return "nomerge"; }
// ⚠ TEMPORARY (user, 2026-09-17; D1 of pp24_trig_eff_hybrid_application.md). ONE region-B cell,
// [105.53,150) x |eta| [2.0,2.2) opposite sign, holds 3 raw MC pairs and is refused by the
// delivery gate. Until the additional high-pT MC statistics being requested arrive, a REFUSED
// un-merged cell is served from the pT-MERGED cell of the same |eta| group ([74.24,150) x
// [2.0,2.2), 62 raw pairs). Once that sample is in, set this to "" (no fallback) so the |eta|
// group goes back to the nominal un-merged cells -- do NOT widen it.
inline const char* PairTrigEffCrossxRefusedCellFallbackMode() { return "ptmerge"; }

#endif // DR_CORRECTION_SAMPLE_CFG_H
