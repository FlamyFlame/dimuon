// fit_dr_corrections.cxx
//
// FIT STAGE of the DeltaR-correction chain (mc_trigger_efficiency.md §3.3 Step 3 / §3.4 Step 4,
// round-7 Autonomy-Contract item 6).
//
//   plot_mc_trig_eff.cxx  -> measures eps_dR / eps_single AND writes the per-(pair pT, pair eta)
//                            large-dR plateau to  dr_correction_plateaus_<label><wp>.root
//   THIS MACRO            -> reads that plateau ROOT file, GUARDS it, divides each cell's dR
//                            curve by its own plateau, and fits the normalized curve
//   plot_dr_correction_fits.cxx -> re-opens the fit file in a FRESH process and draws it
//
// WHAT IS FITTED
//   Step 3: eps_dR^cross(dR)  (PbPb, the cross term dressing eps1*eps2 in the union weight) or
//           eps_dR^2mu4(dR)   (pp, the 2mu4 product correction).
//   Step 4: eps_dR^single(dR) (the single-leg correction dressing the PbPb union LINEAR terms).
//   Per (pair pT, pair eta) cell from the 3D histograms, plus the fully inclusive curve
//   (cell index 0,0). Errors are the round-7 CONDITIONAL ones (dr_correction_ratio.h): with the
//   pre-round-7 errors every chi2/ndf was ~0.25 and every function looked good (R12).
//
// PLATEAU NORMALIZATION -- and why it is read from a ROOT file
//   The correction is only defined up to the large-dR normalization: well-separated muons must
//   decorrelate, so eps -> plateau there, and the physical content is the small-dR shape
//   RELATIVE to that plateau. Each cell is therefore divided by ITS OWN plateau before fitting.
//   The plateau is measured by the step above and travels as DATA. Nothing here parses a .txt
//   or a .md table: a number retyped by hand is wrong the moment the chain is re-run.
//
// PLATEAU MODE (`plateau_mode` argument, added 2026-08-11)
//   "corr"   NOMINAL, and everything above describes it: the curve is divided by the plateau
//            and the fitted shape tends to 1.
//   "nocorr" The plateau is NOT applied. The RAW, un-normalized eps_dR is fitted with a FREE
//            additive baseline C -- expo becomes C + A exp(-(dR/lambda)^p), polyu becomes
//            C + u^2(a2 + a3 u + a4 u^2), and the interpolation pins its flat branch to the LAST
//            MEASURED knot instead of to 1. The baseline is then determined by the dR < 1 data
//            ITSELF and the [2, 3.5] window never enters the fit.
//            WHY (user, 2026-08-11): the inverse-weighted dR distribution carries structure out
//            to large dR -- worst in the pair-eta bins that enclose the detector gap -- so a
//            plateau measured far from the small-dR region may not be the right baseline for it.
//            CONSEQUENCE FOR THE GUARD, stated here because a silently vanishing guard is worse
//            than no guard: in "nocorr" NOTHING is normalized by the plateau, so |plateau-1| can
//            no longer disqualify a cell and a cell whose far-dR plateau is unmeasurable is still
//            perfectly fittable. The plateau is still MEASURED and REPORTED per cell (clearly
//            labelled as not used); only its CONSEQUENCES are switched off. A cell is skipped
//            in this mode for exactly two reasons: too few points, or a failed/unphysical fit.
//   The mode is the TOP level of the plot tree and a file-name token, both built in
//   dr_correction_sample_cfg.h (DrCorrPlateauModeDir / DrCorrPlateauModeTag).
//
// FULL-SAMPLE GUARD
//   For a FULL production (`DrCorrSample::is_full_sample`, currently only pp_full) the
//   inverse-weighting closure must hold cell by cell. TWO TIERS: |plateau - 1| > 0.10 is
//   FLAGGED (reported loudly, not fatal); |plateau - 1| > 0.15 is a FAILURE and this macro
//   lists every offending cell and THROWS. A 10 000-event TEST sample
//   (HIJING overlay, r17663) is exempt -- its high-pair-pT cells are noise-dominated -- but its
//   offending cells are still REPORTED. FULL-vs-TEST comes from the sample identity in
//   dr_correction_sample_cfg.h, never from a list of measured numbers.
//   `allow_plateau_violation=true` downgrades the throw to a loud warning; it exists only so a
//   human can look at the plots of a failing sample, and the driver script never uses it on the
//   nominal path.
//
// SIGN SERIES (`sign` argument; "" = sign-integrated, "ss" = same sign, "os" = opposite sign)
//   The same measurement restricted to one dimuon charge combination. Everything is identical
//   except the input histogram prefix, the plateau key tag, the output file name and the report
//   names (dr_correction_sample_cfg.h).
//   GUARD POLICY -- deliberately asymmetric:
//     * sign-INTEGRATED: unchanged. A FULL production that fails |plateau-1| <= 0.15 THROWS. That
//       series is the nominal deliverable, so a broken inverse-weighting closure must stop it.
//     * sign-SEPARATED: measured, reported and listed cell by cell exactly as above, but NEVER
//       fatal. Each sign carries a FRACTION of the statistics (same sign is ~13% of the selected
//       pairs), so cells that fail there are a STATEMENT ABOUT THE STATISTICS of that subsample,
//       not a defect of the nominal correction.
//     * unchanged for all three: the unmeasurable-cell screen (DrCorrPlateauUsable) and the
//       `fit_ok` flag, so nothing unusable is ever published no matter which series produced it.
//   A sample whose per-sign histograms/plateaus do not exist yet (e.g. an overlay production
//   filled before the per-sign booking) SKIPS the sign-separated series with a printed note --
//   never a throw, so its sign-integrated nominal series stays runnable.
//
// SHAPE CONSTRAINT
//   Step 4 must be flat for dR >~ 0.3, Step 3 for dR >~ 0.5 (kFlatOnsetStep{3,4}).
//
// METHODS (one output file and one plot subdirectory each)
//   powerlaw_fixedRp : f = 1 + A*max(0, 1 - dR/Rp)^n     Rp FIXED at the flat onset.
//                      Exactly 1 and continuous for dR >= Rp -- the constraint is built in.
//   powerlaw_floatRp : same, Rp free within [0.6, 1.6] x the nominal onset.
//   expo             : f = 1 + A*exp(-(dR/lambda)^p)     smooth, approaches 1 asymptotically.
//   interp           : linear interpolation through the measured points below Rp, hard 1 above
//                      (stored as a TGraph of knots -- no free parameters, no chi2).
//
// PERSISTENCE (two traps this repo has already been bitten by)
//   1. A fit is a CONTINUOUS function: evaluate it at the exact dR, never resample-to-nearest.
//   2. Compiled/lambda-based TF1s return 0 outside their range on READ-BACK -- that silently
//      floored a single-muon efficiency and produced a spurious pp cross-section jump
//      (pp_trig_eff_highpt_jump.md). Every TF1 written here is TFORMULA-STRING based and is
//      given a generous range [0, kTF1RangeHi]; the consumer must still CLAMP at the point of
//      use rather than trust range behaviour. plot_dr_correction_fits.cxx re-opens these files
//      in a separate process and re-evaluates them, inside and outside the fit range, as a
//      standing read-back test.
//
// Compile/run (ACLiC, from this directory):
//   root -l -b -q -e '.L fit_dr_corrections.cxx+'                         // compile only
//   root -l -b -q 'fit_dr_corrections.cxx+("pp_full", true, 3, "powerlaw_fixedRp")'
//   root -l -b -q 'fit_dr_corrections.cxx+("pp_full", true, 3, "expo", false, "ss")'   // same sign
//   root -l -b -q 'fit_dr_corrections_all.cxx...'  -> see fit_dr_corrections_all() below

#include <TAxis.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TNamed.h>
#include <TROOT.h>
#include <TSystem.h>

#include "dr_correction_sample_cfg.h"
#include "dr_correction_ratio.h"
#include "../../../Utilities/MCTrigEffPlateauWindow.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

// ---------------------------------------------------------------- shared constants

// Flat-onset radius Rp: beyond it the correction is 1 by construction (user requirement).
// Step 4 (single leg) recovers sooner than Step 3 (both legs), because the cross term needs BOTH
// muons out of the shared L1 RoI / MS sector.
constexpr double kFlatOnsetStep3 = 0.5;
constexpr double kFlatOnsetStep4 = 0.3;

// The zoom histogram spans dR in [0,1] with 20 bins; that is the fit domain. The bins inside
// the PLATEAU WINDOW (MCTrigEffPlateauWindow.h) are what DEFINED the plateau, so re-fitting
// them would just re-fit the normalization.
constexpr double kFitLo = 0.0;
constexpr double kFitHi = 1.0;

// Generous stored range (see the read-back trap in the header comment).
constexpr double kTF1RangeHi = 10.0;

// Guard thresholds on |plateau - 1| for a FULL sample (user decision 2026-08-04, two-tier):
//   > kPlateauFlagTol  -> FLAGGED: listed loudly in the report and on stdout, but NOT fatal.
//   > kPlateauGuardTol -> FAILURE: the macro throws.
// The 0.10 tier is deliberately loose already; 0.15 is the "this cell is not usable" line. The
// gap between them is the band where the inverse-weighting closure is imperfect but the cell is
// still normalizable -- those cells carry |plateau-1| as a systematic
// (docs/systematic_uncertainties.md §1a).
constexpr double kPlateauGuardTol = 0.15;
constexpr double kPlateauFlagTol  = 0.10;


namespace {

struct StepCfg {
    int         step;
    std::string file_suffix;   // input hist file suffix
    std::string h_prefix;      // histogram name prefix
    bool        has_cov;       // Step 4 carries the leg-leg covariance terms
    double      flat_onset;    // Rp
    std::string quantity;      // short name used in the plateau file / provenance
};

// THE ONE PLACE where the sign token enters the histogram names. FillMCTrigEffHists.cxx books the
// per-sign copies under the sign-integrated name with "<sign>_" inserted directly after the common
// prefix: h_mc_dr_ -> h_mc_dr_ss_ , h_mc_single_dr_ -> h_mc_single_dr_os_ .
StepCfg MakeStepCfg(int step, const std::string& sign = "")
{
    const std::string sg = sign.empty() ? "" : sign + "_";
    if (step == 3) return {3, "_step3.root", "h_mc_dr_" + sg,        false, kFlatOnsetStep3, "eps_dR"};
    if (step == 4) return {4, "_step4.root", "h_mc_single_dr_" + sg, true,  kFlatOnsetStep4, "eps_single"};
    throw std::runtime_error("fit_dr_corrections: step must be 3 or 4, got "
                             + std::to_string(step));
}

struct MethodCfg {
    std::string name;
    std::string formula;   // empty => interpolation, no TF1
    int         npar;      // total TF1 parameters (fixed ones included)
    int         nfree;     // free parameters -- drives the "too few points" guard
    bool        float_rp;
};

// `u` below is the reduced separation u = max(0, 1 - dR/Rp): u = 1 at dR = 0, u = 0 at and
// beyond the flat onset Rp. Writing every shape in u makes "exactly 1 beyond Rp" automatic
// instead of a piecewise `if` that TFormula would have to carry.

// `nocorr` = the plateau is NOT applied, so the fitted shape must supply its own asymptote:
// the leading constant 1 of every formula becomes a FREE parameter C, appended as the LAST
// parameter so the meaning of [0], [1], ... is the same in both modes.
MethodCfg MakeMethodCfg(const std::string& m, bool nocorr = false)
{
    // TFormula strings ONLY (never a C++ lambda) so the TF1s survive write/read -- see the
    // read-back trap in the header comment.
    // ---- REJECTED (user 2026-08-04): the two power laws are NOT SMOOTH ----------------------
    // f = 1 + A u^n with u = max(0, 1 - dR/Rp) is continuous in VALUE at Rp but its slope is
    // -A n u^(n-1)/Rp, which DIVERGES for n < 1 -- a cusp. The fitted n rails to its lower limit
    // 0.2 in a large fraction of cells (pp step3 11/37, step4 9/37; overlay step3 19/37,
    // step4 17/37), so this is the typical case, not an edge case. Retained only so an old
    // output can be reproduced; NOT in the driver's default METHODS list.
    if (m == "powerlaw_fixedRp")   // f = 1 + A u^n ,          Rp fixed        -- REJECTED
        return nocorr
            ? MethodCfg{m, "[3]+[0]*TMath::Power(TMath::Max(0.,1.-x/[2]),[1])", 4, 3, false}
            : MethodCfg{m, "1+[0]*TMath::Power(TMath::Max(0.,1.-x/[2]),[1])",   3, 2, false};
    if (m == "powerlaw_floatRp")   // f = 1 + A u^n ,          Rp free         -- REJECTED
        return nocorr
            ? MethodCfg{m, "[3]+[0]*TMath::Power(TMath::Max(0.,1.-x/[2]),[1])", 4, 4, true}
            : MethodCfg{m, "1+[0]*TMath::Power(TMath::Max(0.,1.-x/[2]),[1])",   3, 3, true};
    // ---- NOMINAL (user 2026-08-04) ----------------------------------------------------------
    // Smooth everywhere and -> 1 as dR -> infinity. It approaches 1 ASYMPTOTICALLY rather than
    // reaching it exactly at Rp; that is fine -- "flat for dR >~ 0.5" is an ESTIMATE, not a
    // strict bound. Measured residual |f-1| at dR = 1 over pp cells: median 0.0000, 90th pct
    // 0.0143 (step3) / 0.0000 (step4); the large overlay residuals sit in cells already marked
    // fit_ok = 0. Preferred over polyu because it has no high-order polynomial terms that could
    // fit procedure artefacts rather than real shape.
    if (m == "expo")               // f = 1 + A exp(-(dR/lambda)^p) -- NOMINAL
                                   // nocorr: f = C + A exp(-(dR/lambda)^p), C free
        return nocorr
            ? MethodCfg{m, "[3]+[0]*TMath::Exp(-TMath::Power(x/[1],[2]))", 4, 4, false}
            : MethodCfg{m, "1+[0]*TMath::Exp(-TMath::Power(x/[1],[2]))",   3, 3, false};
    // BACKUP: C^1 at Rp by construction and flexible enough for a non-monotonic small-dR shape,
    // but the higher-order terms can also absorb shapes that are procedure artefacts rather than
    // physics -- which is why `expo` is nominal and this is the cross-check.
    if (m == "polyu_fixedRp")      // f = 1 + a2 u^2 + a3 u^3 + a4 u^4 -- C^1 at Rp (value AND
                                   // slope -> 0), and flexible enough for a non-monotonic
                                   // small-dR shape, which the single power law cannot do
        return nocorr
            ? MethodCfg{m, "[4]+TMath::Power(TMath::Max(0.,1.-x/[3]),2)*([0]+[1]*TMath::Max(0.,1.-x/[3])"
                           "+[2]*TMath::Power(TMath::Max(0.,1.-x/[3]),2))", 5, 4, false}
            : MethodCfg{m, "1+TMath::Power(TMath::Max(0.,1.-x/[3]),2)*([0]+[1]*TMath::Max(0.,1.-x/[3])"
                           "+[2]*TMath::Power(TMath::Max(0.,1.-x/[3]),2))", 4, 3, false};
    if (m == "interp")             // linear interpolation of the measured points, 1 above Rp
        return {m, "", 0, 0, false};
    throw std::runtime_error("fit_dr_corrections: unknown method '" + m +
                             "' (powerlaw_fixedRp | powerlaw_floatRp | expo | polyu_fixedRp | "
                             "interp)");
}

template <typename T>
T* GetObj(TFile* f, const std::string& name)
{
    T* o = dynamic_cast<T*>(f->Get(name.c_str()));
    if (!o) throw std::runtime_error("fit_dr_corrections: missing object '" + name + "' in "
                                     + f->GetName());
    return o;
}

TFile* OpenRead(const std::string& path)
{
    TFile* f = TFile::Open(path.c_str(), "READ");
    if (!f || f->IsZombie())
        throw std::runtime_error("fit_dr_corrections: cannot open " + path);
    return f;
}

TH2D* BookLike(const TH2D* like, const std::string& name, const std::string& ztitle)
{
    auto* h = (TH2D*)like->Clone(name.c_str());
    h->SetDirectory(nullptr);
    h->Reset();
    h->SetTitle((";p_{T}^{pair} [GeV];#eta^{pair};" + ztitle).c_str());
    return h;
}

// A free parameter that came out sitting ON one of its fit limits is NOT a measurement: MINUIT
// parks it on the boundary and still returns a parabolic error, which then spans the limit
// ("#lambda = 0.02 +- 0.058"). Such a parameter is flagged here and on the canvas
// (plot_dr_correction_fits.cxx uses the identical test on the PERSISTED TF1, whose limits survive
// the write/read). The limits themselves are never touched -- changing them would move the
// nominal correction.
bool ParAtLimit(TF1* f, int ip)
{
    double lo = 0., hi = 0.;
    f->GetParLimits(ip, lo, hi);
    if (!(lo < hi)) return false;          // unbounded, or FixParameter (which sets lo == hi)
    const double tol = 1e-3 * (hi - lo);
    const double v   = f->GetParameter(ip);
    return std::fabs(v - lo) <= tol || std::fabs(v - hi) <= tol;
}

// "#lambda" -> "lambda", "R_{p}" -> "R_p": the reports are plain text, the canvas is LaTeX.
std::string ParPlainName(const char* n)
{
    std::string s;
    for (const char* p = n; *p; ++p)
        if (*p != '#' && *p != '{' && *p != '}') s += *p;
    return s;
}

std::string CellName(const std::string& base, int step, int iy, int iz)
{
    // iy = 0 and iz = 0 is the fully inclusive cell.
    if (iy == 0 && iz == 0) return Form("%s_step%d_incl", base.c_str(), step);
    return Form("%s_step%d_pt%d_eta%d", base.c_str(), step, iy, iz);
}

}  // namespace

// ================================================================= main

// step   : 3 (cross / 2mu4 term) or 4 (single leg)
// method : powerlaw_fixedRp | powerlaw_floatRp | expo | interp
// sign   : "" (sign-integrated, the NOMINAL series) | "ss" (same sign) | "os" (opposite sign)
// plateau_mode : "corr" (NOMINAL, divide by the large-dR plateau) | "nocorr" (fit the raw
//                eps_dR with a free baseline C -- see the PLATEAU MODE block in the header)
void fit_dr_corrections(const std::string& sample = "pp_full", bool use_tight_wp = true,
                        int step = 3, const std::string& method = "powerlaw_fixedRp",
                        bool allow_plateau_violation = false, const std::string& sign = "",
                        const std::string& plateau_mode = "corr")
{
    gROOT->SetBatch(kTRUE);

    const DrCorrSample cfg = GetDrCorrSample(sample, use_tight_wp);
    const StepCfg      S   = MakeStepCfg(step, sign);
    // THE mode switch. `mode_dir` also validates the token (it throws on anything else).
    const bool         nocorr   = (plateau_mode == "nocorr");
    const std::string  mode_dir = DrCorrPlateauModeDir(plateau_mode);
    const std::string  mode_text = nocorr
        ? "NO plateau correction (raw eps, free baseline C)" : "plateau-normalized";
    const MethodCfg    M   = MakeMethodCfg(method, nocorr);
    const std::string  wp_suf  = DrCorrWpSuffix(use_tight_wp);
    const std::string  wp_text = use_tight_wp ? "Tight muons" : "Medium muons";
    const std::string  sign_text = DrCorrSignText(sign);          // throws on an unknown token
    const std::string  sign_ftag = DrCorrSignFileTag(sign);
    const std::string  series    = sign.empty() ? "sign-integrated" : sign_text;

    std::cout << "\n================ fit_dr_corrections: " << sample << " / step " << step
              << " / " << method << " / " << wp_text << " / " << series << " / " << mode_text
              << " ================\n";

    // The ONE place the plateau mode enters the output PATH: everything this job writes below
    // (the guard report and the per-method report directory) hangs off `sdir`.
    const std::string sdir = cfg.out_base + "step" + std::to_string(step) + "_dr_fit/" + mode_dir;

    // The pair-pT-binning token must be in BOTH input names: without it a 4-bin run reads the
    // NOMINAL 8-bin histograms while its plateau map is 4x9 (the cell-count guard below catches
    // that, but only after the fact).
    const std::string plateau_path = DrCorrPlateauFile(cfg, use_tight_wp);
    const std::string hist_path = cfg.mc_dir + "mc_trig_eff_hists_" + cfg.mc_label + wp_suf
                                + MCTrigEffPairPt::FileSuffix()
                                + S.file_suffix;
    // Plateau-file key tag. The SIGN is part of it (h_step3_ss_plateau / prov_step3_ss), because
    // the plateau file holds all three series side by side. The keys WRITTEN below keep the plain
    // "step<N>" tag -- the sign is in the fit FILE name, so every consumer reads one set of names.
    const std::string ptag = "step" + std::to_string(step) + (sign.empty() ? "" : "_" + sign);
    const std::string tag  = "step" + std::to_string(step);

    // ---- availability probe for a SIGN-SEPARATED series (non-fatal by design) -----------------
    // A production filled before the per-sign booking (e.g. the HIJING overlay) has neither the
    // per-sign histograms nor the per-sign plateaus. Skipping with a printed note keeps that
    // sample's sign-integrated NOMINAL series runnable; throwing here would take it down too.
    if (!sign.empty()) {
        bool have_plateau = false, have_hists = false;
        if (TFile* p1 = TFile::Open(plateau_path.c_str(), "READ")) {
            if (!p1->IsZombie()) have_plateau = p1->Get(("h_" + ptag + "_plateau").c_str());
            p1->Close();
        }
        if (TFile* p2 = TFile::Open(hist_path.c_str(), "READ")) {
            if (!p2->IsZombie())
                have_hists = p2->Get((S.h_prefix + "zoom_vs_pt_eta_num").c_str());
            p2->Close();
        }
        if (!have_plateau || !have_hists) {
            std::cout << "  SKIPPED: sample '" << sample << "' has no " << sign_text
                      << " inputs yet (" << (have_plateau ? "" : "no h_" + ptag + "_plateau in "
                                                                 + plateau_path + "  ")
                      << (have_hists ? "" : "no " + S.h_prefix + "zoom_vs_pt_eta_num in "
                                            + hist_path)
                      << ").\n  Re-fill the histograms with the per-sign booking to produce it."
                         " The sign-integrated series is unaffected.\n";
            return;
        }
    }

    // ---------------------------------------------------------------- 1. plateaus (ROOT file)
    TFile* fpl = OpenRead(plateau_path);
    TH2D* hplat = GetObj<TH2D>(fpl, "h_" + ptag + "_plateau");
    TH2D* hpnb  = GetObj<TH2D>(fpl, "h_" + ptag + "_plateau_nbins");
    TH1D* hpinc = GetObj<TH1D>(fpl, "h_" + ptag + "_plateau_inclusive");
    // PROVENANCE. Every sample shares the same (pair pT, pair eta) binning, so an edge check
    // could never catch a cross-sample mix-up; the only thing that can is the stamp the producer
    // wrote. Refuse to normalize one sample's curves by another sample's (or another WP's)
    // plateaus.
    {
        TNamed* prov = GetObj<TNamed>(fpl, "prov_" + ptag);
        const std::string title = prov->GetTitle();
        if (title.find("sample=" + sample + ";") == std::string::npos ||
            title.find("WP=" + wp_text + ";") == std::string::npos)
            throw std::runtime_error("fit_dr_corrections: plateau file " + plateau_path
                + " was written for a DIFFERENT sample/WP -- its stamp is '" + title
                + "', this job is sample=" + sample + ", WP=" + wp_text);
    }
    hplat = (TH2D*)hplat->Clone(("plat_" + ptag).c_str()); hplat->SetDirectory(nullptr);
    hpnb  = (TH2D*)hpnb ->Clone(("pnb_"  + ptag).c_str()); hpnb ->SetDirectory(nullptr);
    // Plateau-window normalization systematic |p[kLo,kHi] - p[kSystLo,kSystHi]| (nominal
    // [2, 3.5] against the retired [1, 4]; MCTrigEffPlateauWindow.h), added 2026-08-04. OPTIONAL:
    // a plateau file written before that date has no such key, and the fit is unaffected by it
    // (it is reported, never applied), so an older file must keep working rather than throw.
    TH2D* hpsys = (TH2D*)fpl->Get(("h_" + ptag + "_plateau_syst").c_str());
    if (hpsys) { hpsys = (TH2D*)hpsys->Clone(("psys_" + ptag).c_str()); hpsys->SetDirectory(nullptr); }
    const double plat_incl     = hpinc->GetBinContent(1);
    const double plat_incl_err = hpinc->GetBinError(1);
    fpl->Close();
    std::cout << "  plateaus read from " << plateau_path << "\n"
              << "  inclusive plateau = " << Form("%.4f +- %.4f", plat_incl, plat_incl_err) << "\n";

    const int npt  = hplat->GetNbinsX();
    const int neta = hplat->GetNbinsY();

    // ---------------------------------------------------------------- 2. the guard
    // GUARD POLICY (user, 2026-08-11). The measurement, the reporting and the per-cell lists are
    // IDENTICAL for all three series. Only the CONSEQUENCE differs: a FULL-sample failure is fatal
    // for the sign-INTEGRATED series (the nominal deliverable) and is reported-but-not-fatal for a
    // sign-separated one, because each sign holds only a fraction of the pairs (same sign ~13% of
    // the selected sample) -- a per-sign failure is a statistics statement about that subsample,
    // not a defect in the nominal correction. The unmeasurable-cell screen and the fit_ok flag are
    // untouched by this, so an unusable cell is still never published in ANY series.
    // ... and NEVER in "nocorr": nothing there is normalized by the plateau, so |plateau-1|
    // cannot disqualify a cell. The measurement and the per-cell lists below are produced in
    // both modes; only the consequence differs.
    const bool guard_is_fatal = cfg.is_full_sample && sign.empty() && !nocorr;
    std::vector<std::string> violations, flagged, unmeasured, wsyst;
    double wsyst_max = 0.;
    for (int iy = 1; iy <= npt; ++iy) {
        for (int iz = 1; iz <= neta; ++iz) {
            const double v = hplat->GetBinContent(iy, iz);
            const double e = hplat->GetBinError(iy, iz);
            const char* cell = Form("pT_pair[%.1f,%.1f) x eta_pair[%.1f,%.1f)",
                                    hplat->GetXaxis()->GetBinLowEdge(iy),
                                    hplat->GetXaxis()->GetBinUpEdge(iy),
                                    hplat->GetYaxis()->GetBinLowEdge(iz),
                                    hplat->GetYaxis()->GetBinUpEdge(iz));
            // UNMEASURABLE SCREEN (user policy: such cells are NOT plotted, and their count IS
            // reported; no silent fallback). `v <= 0.` alone was far too weak -- it let through
            // overlay cells with a plateau of 0.0078-0.27, which are then used as DIVISORS, so
            // eps/plateau exploded to O(15) and the "fit" was parameters pinned at their bounds
            // through panels containing no visible data at all. Two additional objective screens:
            //   (a) plateau too far below 1 to be a normalization at all. The plateau is a RATIO
            //       that must tend to 1; a value below kPlateauMinUsable is not a small offset,
            //       it is an empty cell, and dividing by it manufactures a huge correction.
            //   (b) relative error >= 100%: the cell carries no information (e.g. the pp cell
            //       1.0879 +- 1.0879 built from a SINGLE dR bin, which the old guard called "ok"
            //       because it only ever tested |p-1| and never the error).
            const double e_cell = hplat->GetBinError(iy, iz);
            const bool   degenerate    = (v > 0. && v < 0.5);
            const bool   uninformative = (v > 0. && e_cell >= v);
            if (hpnb->GetBinContent(iy, iz) <= 0 || !DrCorrPlateauUsable(v, e_cell)) {
                unmeasured.push_back(std::string(cell) +
                    (degenerate    ? Form("  [plateau %.4f < 0.5: not a normalization]", v)
                   : uninformative ? Form("  [rel. error %.0f%%: no information]", 100. * e_cell / v)
                   : ""));
                continue;
            }
            if (hpsys) {
                const double s = hpsys->GetBinContent(iy, iz);
                if (s >= 0.) {
                    wsyst.push_back(Form("%s : %.4f", cell, s));
                    wsyst_max = std::max(wsyst_max, s);
                }
            }
            const double dev = std::fabs(v - 1.0);
            if (dev > kPlateauGuardTol)
                violations.push_back(Form("%s : plateau = %.4f +- %.4f  (|plateau-1| = %.4f)",
                                          cell, v, e, dev));
            else if (dev > kPlateauFlagTol)
                flagged.push_back(Form("%s : plateau = %.4f +- %.4f  (|plateau-1| = %.4f)",
                                       cell, v, e, dev));
        }
    }
    {
        const std::string gdir = sdir + DrCorrWpDir(use_tight_wp);
        gSystem->mkdir(gdir.c_str(), kTRUE);
        // The tier rule and WHO enforces it. In "nocorr" the tiers are still COMPUTED and the
        // same cells are still listed -- so the two modes' reports can be compared line by line --
        // but nothing enforces them, and the report must say so where it states the rule rather
        // than leave a "the fatal tier is ENFORCED" sentence standing next to a mode that applies
        // no plateau at all. The nominal branch below is byte-for-byte the text this report has
        // always carried.
        std::string rule_block;
        if (nocorr) {
            rule_block =
                std::string("# the two |plateau-1| tiers (FLAGGED > ") + Form("%g", kPlateauFlagTol)
                + ", FAILING > " + Form("%g", kPlateauGuardTol) + ") are computed and listed below"
                  " exactly as in the\n"
                  "#       plateau-corrected mode, so the same cells can be compared line by line."
                  " NEITHER TIER IS ENFORCED HERE, for any series\n"
                  "#       and any production: no curve is divided by the plateau, so no value of"
                  " it can make a cell unusable. The screens that\n"
                  "#       DO apply in this mode are the fit's own: too few points, a failed fit,"
                  " or a correction that is not > 0 over [0, Rp].\n";
        } else {
            rule_block =
                std::string("# rule: a FULL production must satisfy |plateau-1| <= ")
                + Form("%g", kPlateauGuardTol)
                + " in EVERY (pair pT, pair eta) cell (FATAL above that); cells with |plateau-1| > "
                + Form("%g", kPlateauFlagTol) + " are FLAGGED but allowed, and carry |plateau-1| as a\n"
                  "#       systematic (docs/systematic_uncertainties.md 1a). A TEST sample is exempt"
                  " from the fatal tier but still reported.\n"
                + (sign.empty()
                   ? std::string("# this is the SIGN-INTEGRATED series: the fatal tier is ENFORCED"
                                 " on a FULL production.\n")
                   : "# this is the " + series + " series: every flagged/failing/unmeasurable cell is"
                     " measured and listed below exactly as for the\n"
                     "#       sign-integrated series, but the fatal tier is NOT enforced. One sign"
                     " carries only a fraction of the pairs\n"
                     "#       (same sign ~13% of the selected sample), so a failing cell here is a"
                     " statement about the statistics of this\n"
                     "#       subsample, not a defect of the nominal correction. The unmeasurable-cell"
                     " screen and the fit_ok flag are\n"
                     "#       unchanged, so nothing unusable is published from this series either.\n");
        }
        std::ofstream os(gdir + "plateau_guard_report" + sign_ftag + ".txt");
        os << "# Large-dR plateau guard, " << S.quantity << " (Step " << step << ")\n"
           << "# sample=" << sample << " (" << (cfg.is_full_sample ? "FULL" : "TEST")
           << " production)  WP=" << wp_text << "  series=" << series << "\n"
           << "# source: " << plateau_path << "  (keys h_" << ptag << "_*)\n"
           << (nocorr
               ? "# PLATEAU MODE: no plateau correction. NOTHING below is applied to the fits --"
                 " each cell's raw eps is fitted with a FREE\n"
                 "#       baseline C determined from the dR < 1 data alone, so the plateau"
                 " measured here is REPORTED FOR INFORMATION\n"
                 "#       ONLY: it disqualifies no cell, and a cell whose plateau is unmeasurable"
                 " is still fitted. In this mode a cell is\n"
                 "#       skipped only for too few points or a failed/unphysical fit.\n"
               : "")   // the nominal mode's report is left BYTE-IDENTICAL to the pre-2026-08-11
                       // one -- it only moved into plateau_corrected/, and the directory says so

           << rule_block
           << "# plateau window: dR in [" << MCTrigEffPlateau::kLo << ","
           << MCTrigEffPlateau::kHi << "]  (systematic variation: dR in ["
           << MCTrigEffPlateau::kSystLo << "," << MCTrigEffPlateau::kSystHi << "])\n"
           << "# inclusive plateau = " << Form("%.4f +- %.4f", plat_incl, plat_incl_err) << "\n\n";
        os << "unmeasurable cells (no dR bin in the plateau window): " << unmeasured.size() << "\n";
        for (const auto& u : unmeasured) os << "  " << u << "\n";
        os << "\nFLAGGED cells (" << kPlateauFlagTol << " < |plateau-1| <= " << kPlateauGuardTol
           << ", allowed, carry a systematic) : " << flagged.size() << "\n";
        for (const auto& f : flagged) os << "  " << f << "\n";
        os << "\nFAILING cells (|plateau-1| > " << kPlateauGuardTol << ") : " << violations.size()
           << "\n";
        for (const auto& v : violations) os << "  " << v << "\n";
        if (hpsys) {
            os << "\nplateau-WINDOW systematic |plateau[nominal] - plateau[retired [1,4]]| per"
                  " cell (reported, NOT applied here; it is an uncertainty on the normalization,"
                  " docs/systematic_uncertainties.md 1a) -- max = "
               << Form("%.4f", wsyst_max) << "\n";
            for (const auto& w : wsyst) os << "  " << w << "\n";
        }
        os << "\nverdict: " << (nocorr
                              ? "not applicable (no plateau correction: the plateau is measured"
                                " and reported, nothing is normalized by it)"
                              : violations.empty()
                                  ? "PASS"
                                  : (guard_is_fatal
                                        ? "FAIL (FULL sample)"
                                        : (cfg.is_full_sample
                                               ? "reported, not enforced (" + series + " series)"
                                               : "reported, exempt (TEST sample)")))
           << "\n";
        std::cout << "  wrote " << gdir << "plateau_guard_report" << sign_ftag << ".txt\n";
    }
    if (!flagged.empty()) {
        std::cout << "  ~~ FLAGGED: " << flagged.size() << " cell(s) with "
                  << kPlateauFlagTol << " < |plateau-1| <= " << kPlateauGuardTol
                  << " (allowed; each carries |plateau-1| as a systematic):\n";
        for (const auto& f : flagged) std::cout << "     " << f << "\n";
    }
    if (!violations.empty()) {
        std::cout << "  !! " << violations.size() << " cell(s) with |plateau-1| > "
                  << kPlateauGuardTol << ":\n";
        for (const auto& v : violations) std::cout << "     " << v << "\n";
    }
    if (guard_is_fatal && !violations.empty()) {
        const std::string msg =
            "fit_dr_corrections: PLATEAU GUARD FAILED for FULL sample '" + sample + "' (step "
            + std::to_string(step) + ", " + wp_text + "): " + std::to_string(violations.size())
            + " (pair pT, pair eta) cell(s) have |plateau-1| > "
            + std::to_string(kPlateauGuardTol) + " -- see the list above and "
            + cfg.out_base + "step" + std::to_string(step) + "_dr_fit/"
            + DrCorrWpDir(use_tight_wp) + "plateau_guard_report" + sign_ftag + ".txt";
        // Flush BEFORE throwing: an uncaught exception out of a ROOT macro aborts the process,
        // and abort() does not flush stdout -- the violation list printed above would be lost.
        std::cout << std::flush;
        if (!allow_plateau_violation) throw std::runtime_error(msg);
        std::cout << "  ###############################################################\n"
                  << "  ## OVERRIDE: " << msg << "\n"
                  << "  ## allow_plateau_violation=true -- continuing ON PURPOSE.\n"
                  << "  ###############################################################\n";
    } else if (!violations.empty()) {
        std::cout << "  (" << (nocorr ? std::string("no plateau correction")
                             : cfg.is_full_sample ? series + " series"
                                                  : std::string("TEST sample"))
                  << " -> guard NOT enforced; the cells above are reported only.)\n";
    } else {
        std::cout << "  plateau guard PASSED (all |plateau-1| <= " << kPlateauGuardTol
                  << (flagged.empty() ? "" : "; see the FLAGGED cells above") << ").\n";
    }

    // ---------------------------------------------------------------- 3. input histograms
    TFile* fh = OpenRead(hist_path);
    // Staleness: the plateau file is produced BY the histograms, so it must be at least as new.
    // (A plot_mc_trig_eff run that predates a hist refill would normalize the new curves by the
    // old plateaus -- the exact silent-wrong-result this ROOT-file hand-off exists to prevent.)
    {
        Long_t id, sz, fl, mt_h, mt_p;
        gSystem->GetPathInfo(hist_path.c_str(),    &id, &sz, &fl, &mt_h);
        gSystem->GetPathInfo(plateau_path.c_str(), &id, &sz, &fl, &mt_p);
        if (mt_p < mt_h)
            std::cout << "  ** WARNING: the plateau file is OLDER than the histogram file\n"
                      << "     " << plateau_path << "\n     is older than\n     " << hist_path
                      << "\n     -> re-run plot_mc_trig_eff() for this sample/WP before trusting "
                         "these fits.\n";
    }
    TH3D* h3n = GetObj<TH3D>(fh, S.h_prefix + "zoom_vs_pt_eta_num");
    TH3D* h3d = GetObj<TH3D>(fh, S.h_prefix + "zoom_vs_pt_eta_denom");
    TH3D* h3a = GetObj<TH3D>(fh, S.h_prefix + "zoom_vs_pt_eta_errA");
    TH3D* h3b = GetObj<TH3D>(fh, S.h_prefix + "zoom_vs_pt_eta_errB");
    TH3D* h3p = S.has_cov ? GetObj<TH3D>(fh, S.h_prefix + "zoom_vs_pt_eta_covP") : nullptr;
    TH3D* h3q = S.has_cov ? GetObj<TH3D>(fh, S.h_prefix + "zoom_vs_pt_eta_covQ") : nullptr;

    // The plateau map and the histograms MUST describe the same cells, or a cell would be
    // normalized by another cell's plateau. Check the binning, do not assume it.
    if (h3n->GetYaxis()->GetNbins() != npt || h3n->GetZaxis()->GetNbins() != neta)
        throw std::runtime_error("fit_dr_corrections: plateau map (" + std::to_string(npt) + "x"
            + std::to_string(neta) + ") does not match the histogram cells ("
            + std::to_string(h3n->GetYaxis()->GetNbins()) + "x"
            + std::to_string(h3n->GetZaxis()->GetNbins()) + ") -- stale plateau file?");
    for (int iy = 1; iy <= npt + 1; ++iy)
        if (std::fabs(h3n->GetYaxis()->GetBinLowEdge(iy) - hplat->GetXaxis()->GetBinLowEdge(iy))
            > 1e-6)
            throw std::runtime_error("fit_dr_corrections: pair-pT edges differ between the "
                                     "plateau file and the histograms -- stale plateau file?");
    for (int iz = 1; iz <= neta + 1; ++iz)
        if (std::fabs(h3n->GetZaxis()->GetBinLowEdge(iz) - hplat->GetYaxis()->GetBinLowEdge(iz))
            > 1e-6)
            throw std::runtime_error("fit_dr_corrections: pair-eta edges differ between the "
                                     "plateau file and the histograms -- stale plateau file?");

    // ---------------------------------------------------------------- 4. fit every cell
    TH2D* hchi  = BookLike(hplat, "h_" + tag + "_chi2ndf", "#chi^{2}/ndf");
    std::vector<TH2D*> hpar;
    // 4 in the nominal mode (unchanged); one more only when the formula really has a 5th
    // parameter (nocorr polyu, whose free baseline C is [4]).
    for (int ip = 0; ip < std::max(4, M.npar); ++ip)
        hpar.push_back(BookLike(hplat, "h_" + tag + "_par" + std::to_string(ip),
                                "fit parameter " + std::to_string(ip)));
    TH2D* hf0   = BookLike(hplat, "h_" + tag + "_f_at_0",  "fitted correction at #DeltaR = 0");
    TH2D* hstat = BookLike(hplat, "h_" + tag + "_fit_ok",  "1 = fit converged, 0 = no fit");
    TH2D* hknot = BookLike(hplat, "h_" + tag + "_knot_rel_err",
                           "mean relative stat. error of the interpolated points");

    const std::string out_path = DrCorrFitFile(cfg, use_tight_wp, step, method, sign,
                                              plateau_mode);
    TFile* fout = TFile::Open(out_path.c_str(), "RECREATE");
    if (!fout || fout->IsZombie())
        throw std::runtime_error("fit_dr_corrections: cannot write " + out_path);
    fout->cd();

    // one plot/report subdirectory per method (and medium/ inside it for the Medium WP)
    const std::string mdir = sdir + method + "/" + DrCorrWpDir(use_tight_wp);
    gSystem->mkdir(mdir.c_str(), kTRUE);
    std::ofstream rep(mdir + "fit_report" + sign_ftag + ".txt");
    rep << "# " << S.quantity << " (Step " << step << ") "
        << (nocorr ? "RAW (no plateau correction)" : "plateau-normalized")
        << " fit, method = " << method << "\n"
        << "# sample=" << sample << "  WP=" << wp_text << "  series=" << series << "  "
        << cfg.sample_text << "\n"
        << "# formula: "
        << (M.formula.empty()
                ? (nocorr ? "linear interpolation of the measured points, flat at the LAST "
                            "measured knot value above Rp"
                          : "linear interpolation of the measured points, 1 above Rp")
                : M.formula) << "\n"
        << "# flat onset Rp = " << S.flat_onset << " ; fit range dR in [" << kFitLo << ","
        << kFitHi << "]\n"
        << (nocorr
            ? "# NO PLATEAU CORRECTION: nothing is divided by the plateau. The asymptote is the "
              "FREE parameter C, determined by the\n"
              "#   dR < 1 data alone -- the fit sees ONLY the zoom histogram (dR in [0,1], 20 "
              "bins); no bin of the plateau window\n"
              "#   is ever fitted, in this mode or the nominal one. The plateau column below is MEASURED AND REPORTED FOR "
              "INFORMATION ONLY -- it is NOT used, it\n"
              "#   normalizes nothing, and it disqualifies no cell: the |plateau-1| guard and the "
              "unmeasurable-plateau screen do NOT\n"
              "#   apply in this mode (a cell whose far-dR plateau cannot be measured is still "
              "perfectly fittable). A cell is skipped\n"
              "#   here only for too few points or a failed/unphysical fit, and fit_ok = (fit "
              "valid) AND (f > 0 over [0, Rp]).\n"
              "#   Source of the reported plateau: "
            : "# each cell's curve is divided by ITS OWN large-dR plateau (from ")
        << plateau_path << ", keys h_" << ptag << "_*"
        << (nocorr ? "\n" : ") before fitting\n")
        << "# measured from " << hist_path << " (" << S.h_prefix << "zoom_vs_pt_eta_*)\n\n"
        << std::left << std::setw(22) << "pT_pair" << std::setw(16) << "eta_pair"
        << std::setw(12) << "plateau" << std::setw(10) << "npts"
        << std::setw(14) << "p0" << std::setw(14) << "p1" << std::setw(14) << "p2"
        << std::setw(12) << "chi2/ndf" << std::setw(12) << "f(0)" << "\n";

    // Fit-quality bookkeeping. A PLAIN MEAN over all cells is useless for a TEST sample: the
    // 10k-event overlay has cells whose plateau is 0.04 or 16 (too few pairs in the plateau
    // window to
    // define one at all), and normalizing by such a plateau blows the curve up and drags the
    // mean to O(100). Report, in addition: the MEDIAN, the mean over cells whose plateau is
    // sane (|plateau-1| <= the guard tolerance), and the INCLUSIVE cell -- which is the only
    // statistically meaningful curve for a 10k-event sample anyway.
    std::vector<double> chi2_all, chi2_sane;
    int n_unusable = 0;   // EVERY cell persisted with fit_ok = 0, by whichever route: no
                          // measurable plateau / too few points to fit / fit failed / plateau
                          // outside the guard tolerance / correction not > 0 over [0, Rp]
    double chi2_incl = -1.;
    // The inclusive cell has no bin in the per-cell TH2Ds, so its usability flag has to travel
    // separately -- exactly as its plateau does (h_stepN_plateau_inclusive). Without it the plot
    // stage had no way to know that the inclusive fit was rejected and drew it like any other.
    double incl_fit_ok = 0.;
    // F6 bookkeeping: which cells ended with a free parameter pinned on a fit limit.
    std::vector<std::string> atlimit_cells;
    std::map<std::string, int> atlimit_par_count;

    // iy/iz = 0 is the inclusive cell, fitted first.
    for (int iy = 0; iy <= npt; ++iy) {
        for (int iz = 0; iz <= neta; ++iz) {
            const bool inclusive = (iy == 0 && iz == 0);
            if (!inclusive && (iy == 0 || iz == 0)) continue;   // only the full inclusive cell

            const double plateau     = inclusive ? plat_incl : hplat->GetBinContent(iy, iz);
            const double plateau_err = inclusive ? plat_incl_err : hplat->GetBinError(iy, iz);
            const int    pnb     = inclusive ? 1         : (int)hpnb->GetBinContent(iy, iz);
            const std::string nm = CellName("r", step, iy, iz);

            // ONE DECIMAL, as the canvases and the PNG file names already use: at "%.0f" the
            // log edges 11.54 / 16.65 / 24.01 / 34.64 / 49.97 / 72.08 / 103.98 printed as
            // 12 / 17 / 24 / 35 / 50 / 72 / 104, so the human-readable report named cells the
            // figures next to it call 11.5 / 16.6 / 24.0 ... -- two namings of one binning,
            // exactly the drift .claude/CLAUDE.md 'Binnings' exists to stop. Label only; no
            // fitted quantity depends on this string.
            const char* ptlab  = inclusive ? "inclusive"
                : Form("[%.1f,%.1f)", hplat->GetXaxis()->GetBinLowEdge(iy),
                                      hplat->GetXaxis()->GetBinUpEdge(iy));
            const char* etalab = inclusive ? "inclusive"
                : Form("[%.1f,%.1f)", hplat->GetYaxis()->GetBinLowEdge(iz),
                                      hplat->GetYaxis()->GetBinUpEdge(iz));

            // Do not even attempt a fit on an unmeasurable cell -- the same screen the guard
            // and the plot stage use. Fitting one wastes the fit and writes chi2/parameters
            // for a curve divided by a near-zero plateau.
            // "nocorr" DOES NOT DIVIDE, so this screen must not run there: the danger it guards
            // against (a near-zero divisor inflating the curve) does not exist, and a cell whose
            // far-dR plateau is unmeasurable is still perfectly fittable from its dR < 1 points.
            if (!nocorr && (pnb <= 0 || !DrCorrPlateauUsable(plateau, plateau_err))) {
                // n_unusable must count every cell PERSISTED with fit_ok = 0, not only the ones
                // that reach the end of the loop: it is the caption of the fit_ok map, and the
                // two early exits used to write the flag without counting it, so the report
                // under-stated its own map (same-sign expo said 25 for 32 zeroed cells).
                if (!inclusive) ++n_unusable;
                if (!inclusive) hstat->SetBinContent(iy, iz, 0.);
                rep << std::left << std::setw(22) << ptlab << std::setw(16) << etalab
                    << std::setw(12) << "--" << std::setw(10) << 0
                    << "  no plateau -> cell skipped\n";
                continue;
            }

            // measured curve: normalized to 1 at large dR (corr) or RAW (nocorr, where the
            // asymptote is the free parameter C instead of an external divisor)
            TH1D* r = DrCellRatio(h3n, h3d, h3a, h3b, h3p, h3q, iy, iz, nm.c_str());
            if (!nocorr) r->Scale(1.0 / plateau);   // scales contents AND errors

            // points that carry information (a zero-denominator bin has error 0)
            auto* g = new TGraphErrors();
            int k = 0;
            double knot_rel = 0.;
            for (int i = 1; i <= r->GetNbinsX(); ++i) {
                const double x = r->GetBinCenter(i);
                if (x < kFitLo || x > kFitHi) continue;
                const double y = r->GetBinContent(i), e = r->GetBinError(i);
                if (e <= 0.) continue;
                g->SetPoint(k, x, y);
                g->SetPointError(k, 0., e);
                if (y > 0.) knot_rel += e / y;
                ++k;
            }
            delete r;
            g->SetName(CellName("g", step, iy, iz).c_str());
            if (k > 0) knot_rel /= k;

            if (k < M.nfree + 2) {
                // Persisted with fit_ok = 0 -> counted (see the "no plateau" branch above).
                if (!inclusive) ++n_unusable;
                if (!inclusive) hstat->SetBinContent(iy, iz, 0.);
                rep << std::left << std::setw(22) << ptlab << std::setw(16) << etalab
                    << std::setw(12) << Form("%.4f", plateau) << std::setw(10) << k
                    << "  too few points -> no fit\n";
                delete g;
                continue;
            }

            // Starting value (and scale for the limits) of the FREE baseline C in "nocorr":
            // the mean of the measured points in the UPPER part of the FIT DOMAIN, dR in
            // [Rp, 1] -- the flattest part of the data the fit actually sees. It is a starting
            // value only; nothing outside dR < 1 is used, so the [2, 3.5] plateau still plays no
            // role. Falls back to the mean over all fitted points if no point sits above Rp.
            double C0 = 1.;
            if (nocorr) {
                double sum = 0.; int nc = 0;
                for (int i = 0; i < g->GetN(); ++i) {
                    double x, y;
                    g->GetPoint(i, x, y);
                    if (x >= S.flat_onset) { sum += y; ++nc; }
                }
                if (nc == 0)
                    for (int i = 0; i < g->GetN(); ++i) {
                        double x, y;
                        g->GetPoint(i, x, y);
                        sum += y; ++nc;
                    }
                C0 = (nc > 0 && sum > 0.) ? sum / nc : 1.;
            }
            // Generous but FINITE limits, scaled to the data: an unbounded baseline lets MINUIT
            // trade C against A without ever converging.
            const double C_lo = 0., C_hi = std::max(5.0 * C0, 2.0);

            if (M.formula.empty()) {
                // ---- interpolation: knots below Rp, then a hard 1 --------------------------
                // No free parameters and no chi2: it passes through every point, so it also
                // transmits every point's statistical fluctuation straight into the correction.
                // That transmitted noise is what `knot_rel_err` records.
                auto* gk = new TGraph();
                int kk = 0;
                bool first_knot = true;
                double last_knot = C0;    // nocorr fallback if no point sits below Rp
                for (int i = 0; i < g->GetN(); ++i) {
                    double x, y;
                    g->GetPoint(i, x, y);
                    if (x >= S.flat_onset) continue;
                    // Pin dR = 0 to the first measured value: TGraph::Eval extrapolates
                    // LINEARLY below the first knot, which would invent a dR -> 0 trend the
                    // data does not contain (the first bin centre is 0.025).
                    if (first_knot) { gk->SetPoint(kk++, 0.0, y); first_knot = false; }
                    gk->SetPoint(kk++, x, y);
                    last_knot = y;
                }
                // Pin the flat branch out to the generous range end. corr: the value is 1 by
                // construction (that is what the plateau normalization bought). nocorr: there is
                // no external normalization to pin to, so it is pinned to the LAST MEASURED KNOT
                // -- the baseline the data itself supplies just below Rp.
                const double flat_val = nocorr ? last_knot : 1.0;
                for (double x : {S.flat_onset, 1.0, 2.0, 5.0, kTF1RangeHi})
                    if (x >= S.flat_onset) gk->SetPoint(kk++, x, flat_val);
                gk->SetName(CellName("gknots", step, iy, iz).c_str());
                gk->SetTitle(Form("linear-interpolation knots; #DeltaR; %s (%s)",
                                  S.quantity.c_str(),
                                  nocorr ? "no plateau correction" : "plateau-normalized"));
                gk->Write();
                // USABILITY. corr: the plateau screens (the curve is only meaningful relative
                // to a sane plateau). nocorr: no plateau is involved, so the only requirement is
                // that the interpolated correction is physical -- positive over [0, Rp].
                auto interp_physical = [&]() {
                    for (int t = 0; t <= 200; ++t)
                        if (gk->Eval(kFitLo + (S.flat_onset - kFitLo) * t / 200.0) <= 0.)
                            return false;
                    return true;
                };
                const bool interp_ok = nocorr
                    ? interp_physical()
                    : (DrCorrPlateauUsable(plateau, plateau_err)
                       && std::fabs(plateau - 1.0) <= kPlateauGuardTol);
                if (inclusive) incl_fit_ok = interp_ok ? 1. : 0.;
                if (!inclusive) {
                    // An interpolation always "succeeds" numerically, but that says nothing about
                    // whether the CELL is usable. This branch used to write fit_ok = 1
                    // unconditionally, publishing 20 overlay cells (e.g. plateau 0.0078 +- 0.0051,
                    // a x128 inflation) that the guard listed as unmeasurable and the plot stage
                    // refused to draw. Same screen as the parametric branch.
                    if (!interp_ok) ++n_unusable;
                    hstat->SetBinContent(iy, iz, interp_ok ? 1. : 0.);
                    hchi ->SetBinContent(iy, iz, -1.);     // n/a for an interpolation
                    hknot->SetBinContent(iy, iz, knot_rel);
                    hf0  ->SetBinContent(iy, iz, gk->Eval(0.0));
                }
                rep << std::left << std::setw(22) << ptlab << std::setw(16) << etalab
                    << std::setw(12) << Form("%.4f", plateau) << std::setw(10) << k
                    << std::setw(14) << "--" << std::setw(14) << "--" << std::setw(14) << "--"
                    << std::setw(12) << "n/a"
                    << std::setw(12) << Form("%.4f", gk->Eval(0.0))
                    << (nocorr ? Form("  C(flat branch)=%.4f", flat_val) : "")
                    << "  mean rel. stat. err of knots = " << Form("%.4f", knot_rel) << "\n";
                g->Write();
                delete g;
                continue;
            }

            // ---- parametric fit ------------------------------------------------------------
            auto* f = new TF1(CellName("f", step, iy, iz).c_str(), M.formula.c_str(),
                              kFitLo, kFitHi);
            double y0 = 0., x0 = 0.;
            g->GetPoint(0, x0, y0);
            // The amplitude is measured FROM THE BASELINE: 1 when the curve was normalized,
            // the fitted-baseline estimate C0 when it was not.
            const double A0 = y0 - (nocorr ? C0 : 1.0);
            if (method == "expo") {
                if (nocorr) {
                    f->SetParNames("A", "#lambda", "p", "C");
                    f->SetParameters(A0 != 0. ? A0 : 0.2, 0.25, 1.5, C0);
                } else {
                    f->SetParNames("A", "#lambda", "p");
                    f->SetParameters(A0 != 0. ? A0 : 0.2, 0.25, 1.5);
                }
                f->SetParLimits(0, -5.0, 20.0);
                f->SetParLimits(1, 0.02, 3.0);
                f->SetParLimits(2, 0.3, 8.0);
                if (nocorr) f->SetParLimits(3, C_lo, C_hi);
            } else if (method == "polyu_fixedRp") {
                if (nocorr) {
                    f->SetParNames("a_{2}", "a_{3}", "a_{4}", "R_{p}", "C");
                    f->SetParameters(A0 != 0. ? A0 : 0.2, 0.0, 0.0, S.flat_onset, C0);
                } else {
                    f->SetParNames("a_{2}", "a_{3}", "a_{4}", "R_{p}");
                    f->SetParameters(A0 != 0. ? A0 : 0.2, 0.0, 0.0, S.flat_onset);
                }
                f->SetParLimits(0, -50.0, 50.0);
                f->SetParLimits(1, -100.0, 100.0);
                f->SetParLimits(2, -100.0, 100.0);
                f->FixParameter(3, S.flat_onset);
                if (nocorr) f->SetParLimits(4, C_lo, C_hi);
            } else {
                if (nocorr) {
                    f->SetParNames("A", "n", "R_{p}", "C");
                    f->SetParameters(A0 != 0. ? A0 : 0.2, 2.0, S.flat_onset, C0);
                } else {
                    f->SetParNames("A", "n", "R_{p}");
                    f->SetParameters(A0 != 0. ? A0 : 0.2, 2.0, S.flat_onset);
                }
                f->SetParLimits(0, -5.0, 20.0);
                f->SetParLimits(1, 0.2, 20.0);
                if (M.float_rp) f->SetParLimits(2, 0.6 * S.flat_onset, 1.6 * S.flat_onset);
                else            f->FixParameter(2, S.flat_onset);
                if (nocorr) f->SetParLimits(3, C_lo, C_hi);
            }

            TFitResultPtr fr = g->Fit(f, "QRNS");
            const bool ok = (fr.Get() != nullptr) && fr->IsValid() && fr->Ndf() > 0;
            const double chi2ndf = ok ? fr->Chi2() / fr->Ndf() : -1.;

            // Store with a GENEROUS range: a TFormula TF1 evaluates its formula everywhere, but
            // a wide range keeps Draw()/Integral() honest too. The CONSUMER still clamps.
            f->SetRange(kFitLo, kTF1RangeHi);
            f->Write();
            g->Write();

            // USABILITY, not just "the fitter returned": a cell is usable only if the fit
            // converged AND its plateau is inside the guard tolerance AND the resulting
            // correction is physical (a trigger-probability correction can never be <= 0) over
            // the whole fitted range. Without the last two conditions the overlay persisted
            // cells with f(0) = -1.005 / -0.559 carrying fit_ok = 1, i.e. a NEGATIVE trigger
            // correction advertised as good -- the class of trap pp_trig_eff_highpt_jump.md
            // records. Consumers must require h_stepN_fit_ok == 1.
            bool physical = true;
            for (int k = 0; k <= 200; ++k) {
                const double x = kFitLo + (S.flat_onset - kFitLo) * k / 200.0;
                if (f->Eval(x) <= 0.) { physical = false; break; }
            }
            // fit_ok is the flag CONSUMERS gate on, so it must apply the SAME screen as the
            // guard and the plot stage. It previously tested only |plateau-1|, never the plateau
            // ERROR, so a cell like 1.0879 +- 1.0879 (a single dR bin, 100% relative error) was
            // published as fit_ok = 1 while the guard called it unmeasurable and the plot drew
            // "no fit" -- the artefact the analysis consumes disagreed with both.
            // In "nocorr" the two plateau conditions are DROPPED, and deliberately so: nothing
            // is normalized by the plateau there, so it cannot make a fit unusable. What remains
            // is the part that is about the fit itself -- it converged, and the correction is
            // positive over [0, Rp].
            const bool usable = ok && physical
                             && (nocorr || (DrCorrPlateauUsable(plateau, plateau_err)
                                            && std::fabs(plateau - 1.0) <= kPlateauGuardTol));
            if (!usable && !inclusive) ++n_unusable;
            if (inclusive) incl_fit_ok = usable ? 1. : 0.;
            if (!inclusive) {
                hstat->SetBinContent(iy, iz, usable ? 1. : 0.);
                hchi ->SetBinContent(iy, iz, chi2ndf);
                for (int ip = 0; ip < M.npar && ip < (int)hpar.size(); ++ip) {
                    hpar[ip]->SetBinContent(iy, iz, f->GetParameter(ip));
                    hpar[ip]->SetBinError  (iy, iz, f->GetParError(ip));
                }
                hf0  ->SetBinContent(iy, iz, f->Eval(0.0));
                hknot->SetBinContent(iy, iz, knot_rel);
            }
            if (ok) {
                if (inclusive) chi2_incl = chi2ndf;
                else {
                    chi2_all.push_back(chi2ndf);
                    // The "sane plateau" subset exists because a curve divided by a bad plateau
                    // has an inflated chi2. Nothing is divided in "nocorr", so the subset has no
                    // meaning there and is left empty (the report says so).
                    if (!nocorr && std::fabs(plateau - 1.0) <= kPlateauFlagTol)
                        chi2_sane.push_back(chi2ndf);
                }
            }

            // F6: a parameter pinned on its limit is reported as such, here and on the canvas.
            std::string atlim;
            for (int ip = 0; ip < M.npar; ++ip) {
                if (!ParAtLimit(f, ip)) continue;
                const std::string pn = ParPlainName(f->GetParName(ip));
                atlim += (atlim.empty() ? "" : ",") + pn;
                ++atlimit_par_count[pn];
            }
            if (!atlim.empty())
                atlimit_cells.push_back(std::string(ptlab) + " x " + etalab + " : " + atlim);

            rep << std::left << std::setw(22) << ptlab << std::setw(16) << etalab
                << std::setw(12) << Form("%.4f", plateau) << std::setw(10) << k;
            for (int ip = 0; ip < 3; ++ip)
                rep << std::setw(14) << (ip < M.npar ? Form("%.4f", f->GetParameter(ip)) : "--");
            rep << std::setw(12) << (ok ? Form("%.3f", chi2ndf) : "FAILED")
                << std::setw(12) << Form("%.4f", f->Eval(0.0))
                ;
            for (int ip = 3; ip < M.npar; ++ip)
                rep << (nocorr && ip == M.npar - 1
                            ? Form("  C=%.4f+-%.4f", f->GetParameter(ip), f->GetParError(ip))
                            : Form("  p%d=%.4f", ip, f->GetParameter(ip)));
            rep << (atlim.empty() ? "" : "  AT LIMIT: " + atlim) << "\n";
            delete f;
            delete g;
        }
    }

    // plateau maps travel WITH the fits, so a consumer of the fit file needs nothing else
    hplat->SetName(("h_" + tag + "_plateau").c_str());
    hplat->Write();
    hpnb ->SetName(("h_" + tag + "_plateau_nbins").c_str());
    hpnb ->Write();
    hchi->Write();
    for (TH2D* h : hpar) h->Write();
    hf0->Write(); hstat->Write(); hknot->Write();
    {
        auto* hi = new TH1D(("h_" + tag + "_plateau_inclusive").c_str(), "", 1, 0., 1.);
        hi->SetDirectory(nullptr);
        hi->SetBinContent(1, plat_incl);
        hi->SetBinError(1, plat_incl_err);
        hi->Write();
        delete hi;
        // Consumers must gate on fit_ok; the inclusive cell needs its own carrier.
        auto* hk = new TH1D(("h_" + tag + "_fit_ok_inclusive").c_str(),
                            ";;1 = fit usable, 0 = rejected", 1, 0., 1.);
        hk->SetDirectory(nullptr);
        hk->SetBinContent(1, incl_fit_ok);
        hk->Write();
        delete hk;
    }
    TNamed("provenance",
           Form("sample=%s (%s); WP=%s; series=%s; plateau mode=%s; step=%d (%s); method=%s; "
                "formula=%s; Rp=%.2f; "
                "fit range dR=[%.2f,%.2f]; stored TF1 range=[0,%.1f]; plateau source=%s (keys "
                "h_%s_*); histograms=%s (%szoom_vs_pt_eta_*); guard=%s; "
                "producer=fit_dr_corrections.cxx",
                sample.c_str(), cfg.is_full_sample ? "FULL" : "TEST", wp_text.c_str(),
                series.c_str(),
                nocorr ? "nocorr (raw eps, free baseline C -- the plateau is NOT applied)"
                       : "corr (each cell divided by its own large-dR plateau)",
                step,
                S.quantity.c_str(), method.c_str(),
                M.formula.empty() ? "linear interpolation (TGraph knots)" : M.formula.c_str(),
                S.flat_onset, kFitLo, kFitHi, kTF1RangeHi, plateau_path.c_str(), ptag.c_str(),
                hist_path.c_str(), S.h_prefix.c_str(),
                violations.empty() ? "PASS"
                                   : (guard_is_fatal ? "FAIL(overridden)"
                                                     : "reported, not enforced")))
        .Write();
    fout->Close();
    fh->Close();

    auto stats = [](std::vector<double> v) {
        if (v.empty()) return std::string("n/a");
        std::sort(v.begin(), v.end());
        double s = 0.;
        for (double x : v) s += x;
        const double med = (v.size() % 2) ? v[v.size() / 2]
                                          : 0.5 * (v[v.size() / 2 - 1] + v[v.size() / 2]);
        return std::string(Form("mean %.3f, median %.3f  (n=%zu)", s / v.size(), med, v.size()));
    };
    const std::string line_all  = stats(chi2_all);
    const std::string line_sane = stats(chi2_sane);
    const std::string line_incl = (chi2_incl >= 0.) ? Form("%.3f", chi2_incl) : "n/a";
    rep << "\n# chi2/ndf, INCLUSIVE cell (the only statistically meaningful curve for a 10k-event"
           " TEST sample) = " << line_incl << "\n"
        << "# chi2/ndf over all converged cells:            " << line_all << "\n"
        << (nocorr
            ? Form("# cells marked UNUSABLE (h_stepN_fit_ok = 0: too few points to fit, fit"
                   " failed, or the correction is not > 0 over [0, Rp] -- the plateau screens do"
                   " NOT apply in this mode): %d  -- consumers MUST require fit_ok == 1\n",
                   n_unusable)
            : Form("# cells marked UNUSABLE (h_stepN_fit_ok = 0: no measurable plateau, too few"
                   " points to fit, fit failed, |plateau-1| > %g, or the correction is not > 0"
                   " over [0, Rp]): %d  -- consumers MUST require fit_ok == 1\n",
                   kPlateauGuardTol, n_unusable))
        << (nocorr
            ? std::string("# chi2/ndf over cells with a near-unity plateau: n/a -- no curve is"
                          " divided by the plateau in this mode, so that subset carries no"
                          " information here.\n")
            : std::string("# chi2/ndf over cells with |plateau-1| <= ")
              + Form("%g", kPlateauFlagTol) + ": " + line_sane + "\n"
              + "# (cells whose plateau is far from 1 have too few pairs inside the plateau window"
              + Form(" dR in [%g,%g] to define one;", MCTrigEffPlateau::kLo, MCTrigEffPlateau::kHi)
              + " normalizing by such a plateau inflates chi2 without saying anything about the"
                " fit function.)\n")
        // kept for backwards compatibility with the driver's summary parser
        << "# parameters pinned ON a fit limit (MINUIT parks the value on the boundary and still"
           " returns an error that spans it -- such a value is a constraint, not a measurement;"
           " marked '(at limit)' on the canvas): " << atlimit_cells.size() << " cell(s)\n";
    for (const auto& kv : atlimit_par_count)
        rep << "#   " << kv.first << ": " << kv.second << " cell(s)\n";
    for (const auto& c : atlimit_cells) rep << "#     " << c << "\n";
    rep << "# mean chi2/ndf over " << chi2_all.size() << " converged cells = "
        << (chi2_all.empty() ? "n/a"
                             : Form("%.3f", std::accumulate(chi2_all.begin(), chi2_all.end(), 0.)
                                            / chi2_all.size()))
        << "\n";
    rep.close();

    std::cout << "  chi2/ndf inclusive = " << line_incl << "\n"
              << "  chi2/ndf all cells: " << line_all << "\n"
              << "  chi2/ndf sane-plateau cells: " << line_sane << "\n"
              << "  wrote " << out_path << "\n"
              << "  wrote " << mdir << "fit_report" << sign_ftag << ".txt\n";
}

// Convenience: every method for one (sample, WP, step, sign).
void fit_dr_corrections_all(const std::string& sample = "pp_full", bool use_tight_wp = true,
                            int step = 3, bool allow_plateau_violation = false,
                            const std::string& sign = "",
                            const std::string& plateau_mode = "corr")
{
    for (const std::string& m : {"powerlaw_fixedRp", "powerlaw_floatRp", "expo",
                                 "polyu_fixedRp", "interp"})
        fit_dr_corrections(sample, use_tight_wp, step, m, allow_plateau_violation, sign,
                           plateau_mode);
}
