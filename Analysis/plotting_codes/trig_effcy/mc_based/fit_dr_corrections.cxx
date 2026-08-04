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
// FULL-SAMPLE GUARD
//   For a FULL production (`DrCorrSample::is_full_sample`, currently only pp_full) the
//   inverse-weighting closure must hold cell by cell, so |plateau - 1| > 0.1 anywhere is a
//   FAILURE: this macro lists every offending cell and THROWS. A 10 000-event TEST sample
//   (HIJING overlay, r17663) is exempt -- its high-pair-pT cells are noise-dominated -- but its
//   offending cells are still REPORTED. FULL-vs-TEST comes from the sample identity in
//   dr_correction_sample_cfg.h, never from a list of measured numbers.
//   `allow_plateau_violation=true` downgrades the throw to a loud warning; it exists only so a
//   human can look at the plots of a failing sample, and the driver script never uses it on the
//   nominal path.
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

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
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

// The zoom histogram spans dR in [0,1] with 20 bins; that is the fit domain. The bins in
// [1,4] are what DEFINED the plateau, so re-fitting them would just re-fit the normalization.
constexpr double kFitLo = 0.0;
constexpr double kFitHi = 1.0;

// Generous stored range (see the read-back trap in the header comment).
constexpr double kTF1RangeHi = 10.0;

// Guard threshold on |plateau - 1| for a FULL sample.
constexpr double kPlateauGuardTol = 0.1;

namespace {

struct StepCfg {
    int         step;
    std::string file_suffix;   // input hist file suffix
    std::string h_prefix;      // histogram name prefix
    bool        has_cov;       // Step 4 carries the leg-leg covariance terms
    double      flat_onset;    // Rp
    std::string quantity;      // short name used in the plateau file / provenance
};

StepCfg MakeStepCfg(int step)
{
    if (step == 3) return {3, "_step3.root", "h_mc_dr_",        false, kFlatOnsetStep3, "eps_dR"};
    if (step == 4) return {4, "_step4.root", "h_mc_single_dr_", true,  kFlatOnsetStep4, "eps_single"};
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

MethodCfg MakeMethodCfg(const std::string& m)
{
    // TFormula strings ONLY (never a C++ lambda) so the TF1s survive write/read -- see the
    // read-back trap in the header comment.
    if (m == "powerlaw_fixedRp")   // f = 1 + A u^n ,          Rp fixed
        return {m, "1+[0]*TMath::Power(TMath::Max(0.,1.-x/[2]),[1])", 3, 2, false};
    if (m == "powerlaw_floatRp")   // f = 1 + A u^n ,          Rp free in a small window
        return {m, "1+[0]*TMath::Power(TMath::Max(0.,1.-x/[2]),[1])", 3, 3, true};
    if (m == "expo")               // f = 1 + A exp(-(dR/lambda)^p) -- smooth, 1 asymptotically
        return {m, "1+[0]*TMath::Exp(-TMath::Power(x/[1],[2]))", 3, 3, false};
    if (m == "polyu_fixedRp")      // f = 1 + a2 u^2 + a3 u^3 + a4 u^4 -- C^1 at Rp (value AND
                                   // slope -> 0), and flexible enough for a non-monotonic
                                   // small-dR shape, which the single power law cannot do
        return {m, "1+TMath::Power(TMath::Max(0.,1.-x/[3]),2)*([0]+[1]*TMath::Max(0.,1.-x/[3])"
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
void fit_dr_corrections(const std::string& sample = "pp_full", bool use_tight_wp = true,
                        int step = 3, const std::string& method = "powerlaw_fixedRp",
                        bool allow_plateau_violation = false)
{
    gROOT->SetBatch(kTRUE);

    const DrCorrSample cfg = GetDrCorrSample(sample);
    const StepCfg      S   = MakeStepCfg(step);
    const MethodCfg    M   = MakeMethodCfg(method);
    const std::string  wp_suf  = DrCorrWpSuffix(use_tight_wp);
    const std::string  wp_text = use_tight_wp ? "Tight muons" : "Medium muons";

    std::cout << "\n================ fit_dr_corrections: " << sample << " / step " << step
              << " / " << method << " / " << wp_text << " ================\n";

    // ---------------------------------------------------------------- 1. plateaus (ROOT file)
    const std::string plateau_path = DrCorrPlateauFile(cfg, use_tight_wp);
    TFile* fpl = OpenRead(plateau_path);
    const std::string tag = "step" + std::to_string(step);
    TH2D* hplat = GetObj<TH2D>(fpl, "h_" + tag + "_plateau");
    TH2D* hpnb  = GetObj<TH2D>(fpl, "h_" + tag + "_plateau_nbins");
    TH1D* hpinc = GetObj<TH1D>(fpl, "h_" + tag + "_plateau_inclusive");
    // PROVENANCE. Every sample shares the same (pair pT, pair eta) binning, so an edge check
    // could never catch a cross-sample mix-up; the only thing that can is the stamp the producer
    // wrote. Refuse to normalize one sample's curves by another sample's (or another WP's)
    // plateaus.
    {
        TNamed* prov = GetObj<TNamed>(fpl, "prov_" + tag);
        const std::string title = prov->GetTitle();
        if (title.find("sample=" + sample + ";") == std::string::npos ||
            title.find("WP=" + wp_text + ";") == std::string::npos)
            throw std::runtime_error("fit_dr_corrections: plateau file " + plateau_path
                + " was written for a DIFFERENT sample/WP -- its stamp is '" + title
                + "', this job is sample=" + sample + ", WP=" + wp_text);
    }
    hplat = (TH2D*)hplat->Clone(("plat_" + tag).c_str()); hplat->SetDirectory(nullptr);
    hpnb  = (TH2D*)hpnb ->Clone(("pnb_"  + tag).c_str()); hpnb ->SetDirectory(nullptr);
    const double plat_incl     = hpinc->GetBinContent(1);
    const double plat_incl_err = hpinc->GetBinError(1);
    fpl->Close();
    std::cout << "  plateaus read from " << plateau_path << "\n"
              << "  inclusive plateau = " << Form("%.4f +- %.4f", plat_incl, plat_incl_err) << "\n";

    const int npt  = hplat->GetNbinsX();
    const int neta = hplat->GetNbinsY();

    // ---------------------------------------------------------------- 2. the guard
    std::vector<std::string> violations, unmeasured;
    for (int iy = 1; iy <= npt; ++iy) {
        for (int iz = 1; iz <= neta; ++iz) {
            const double v = hplat->GetBinContent(iy, iz);
            const double e = hplat->GetBinError(iy, iz);
            const char* cell = Form("pT_pair[%.0f,%.0f) x eta_pair[%.1f,%.1f)",
                                    hplat->GetXaxis()->GetBinLowEdge(iy),
                                    hplat->GetXaxis()->GetBinUpEdge(iy),
                                    hplat->GetYaxis()->GetBinLowEdge(iz),
                                    hplat->GetYaxis()->GetBinUpEdge(iz));
            if (hpnb->GetBinContent(iy, iz) <= 0 || v <= 0.) {
                unmeasured.push_back(cell);
                continue;
            }
            if (std::fabs(v - 1.0) > kPlateauGuardTol)
                violations.push_back(Form("%s : plateau = %.4f +- %.4f  (|plateau-1| = %.4f)",
                                          cell, v, e, std::fabs(v - 1.0)));
        }
    }
    {
        const std::string gdir = cfg.out_base + "step" + std::to_string(step) + "_dr_fit/"
                               + DrCorrWpDir(use_tight_wp);
        gSystem->mkdir(gdir.c_str(), kTRUE);
        std::ofstream os(gdir + "plateau_guard_report.txt");
        os << "# Large-dR plateau guard, " << S.quantity << " (Step " << step << ")\n"
           << "# sample=" << sample << " (" << (cfg.is_full_sample ? "FULL" : "TEST")
           << " production)  WP=" << wp_text << "\n"
           << "# source: " << plateau_path << "\n"
           << "# rule: a FULL production must satisfy |plateau-1| <= " << kPlateauGuardTol
           << " in EVERY (pair pT, pair eta) cell; a TEST sample is exempt but still reported.\n"
           << "# inclusive plateau = " << Form("%.4f +- %.4f", plat_incl, plat_incl_err) << "\n\n";
        os << "unmeasurable cells (no dR bin in the plateau window): " << unmeasured.size() << "\n";
        for (const auto& u : unmeasured) os << "  " << u << "\n";
        os << "\ncells with |plateau-1| > " << kPlateauGuardTol << " : " << violations.size()
           << "\n";
        for (const auto& v : violations) os << "  " << v << "\n";
        os << "\nverdict: " << (violations.empty() ? "PASS"
                                                   : (cfg.is_full_sample ? "FAIL (FULL sample)"
                                                                         : "reported, exempt (TEST sample)"))
           << "\n";
        std::cout << "  wrote " << gdir << "plateau_guard_report.txt\n";
    }
    if (!violations.empty()) {
        std::cout << "  !! " << violations.size() << " cell(s) with |plateau-1| > "
                  << kPlateauGuardTol << ":\n";
        for (const auto& v : violations) std::cout << "     " << v << "\n";
    }
    if (cfg.is_full_sample && !violations.empty()) {
        const std::string msg =
            "fit_dr_corrections: PLATEAU GUARD FAILED for FULL sample '" + sample + "' (step "
            + std::to_string(step) + ", " + wp_text + "): " + std::to_string(violations.size())
            + " (pair pT, pair eta) cell(s) have |plateau-1| > "
            + std::to_string(kPlateauGuardTol) + " -- see the list above and "
            + cfg.out_base + "step" + std::to_string(step) + "_dr_fit/"
            + DrCorrWpDir(use_tight_wp) + "plateau_guard_report.txt";
        // Flush BEFORE throwing: an uncaught exception out of a ROOT macro aborts the process,
        // and abort() does not flush stdout -- the violation list printed above would be lost.
        std::cout << std::flush;
        if (!allow_plateau_violation) throw std::runtime_error(msg);
        std::cout << "  ###############################################################\n"
                  << "  ## OVERRIDE: " << msg << "\n"
                  << "  ## allow_plateau_violation=true -- continuing ON PURPOSE.\n"
                  << "  ###############################################################\n";
    } else if (!cfg.is_full_sample && !violations.empty()) {
        std::cout << "  (TEST sample -> guard NOT enforced; the cells above are reported only.)\n";
    } else {
        std::cout << "  plateau guard PASSED (all |plateau-1| <= " << kPlateauGuardTol << ").\n";
    }

    // ---------------------------------------------------------------- 3. input histograms
    const std::string hist_path = cfg.mc_dir + "mc_trig_eff_hists_" + cfg.mc_label + wp_suf
                                + S.file_suffix;
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
    for (int ip = 0; ip < 4; ++ip)
        hpar.push_back(BookLike(hplat, "h_" + tag + "_par" + std::to_string(ip),
                                "fit parameter " + std::to_string(ip)));
    TH2D* hf0   = BookLike(hplat, "h_" + tag + "_f_at_0",  "fitted correction at #DeltaR = 0");
    TH2D* hstat = BookLike(hplat, "h_" + tag + "_fit_ok",  "1 = fit converged, 0 = no fit");
    TH2D* hknot = BookLike(hplat, "h_" + tag + "_knot_rel_err",
                           "mean relative stat. error of the interpolated points");

    const std::string out_path = DrCorrFitFile(cfg, use_tight_wp, step, method);
    TFile* fout = TFile::Open(out_path.c_str(), "RECREATE");
    if (!fout || fout->IsZombie())
        throw std::runtime_error("fit_dr_corrections: cannot write " + out_path);
    fout->cd();

    // one plot/report subdirectory per method (and medium/ inside it for the Medium WP)
    const std::string mdir = cfg.out_base + "step" + std::to_string(step) + "_dr_fit/" + method
                           + "/" + DrCorrWpDir(use_tight_wp);
    gSystem->mkdir(mdir.c_str(), kTRUE);
    std::ofstream rep(mdir + "fit_report.txt");
    rep << "# " << S.quantity << " (Step " << step << ") plateau-normalized fit, method = "
        << method << "\n"
        << "# sample=" << sample << "  WP=" << wp_text << "  " << cfg.sample_text << "\n"
        << "# formula: " << (M.formula.empty() ? "linear interpolation of the measured points, "
                                                 "1 above Rp" : M.formula) << "\n"
        << "# flat onset Rp = " << S.flat_onset << " ; fit range dR in [" << kFitLo << ","
        << kFitHi << "]\n"
        << "# each cell's curve is divided by ITS OWN large-dR plateau (from " << plateau_path
        << ") before fitting\n\n"
        << std::left << std::setw(22) << "pT_pair" << std::setw(16) << "eta_pair"
        << std::setw(12) << "plateau" << std::setw(10) << "npts"
        << std::setw(14) << "p0" << std::setw(14) << "p1" << std::setw(14) << "p2"
        << std::setw(12) << "chi2/ndf" << std::setw(12) << "f(0)" << "\n";

    // Fit-quality bookkeeping. A PLAIN MEAN over all cells is useless for a TEST sample: the
    // 10k-event overlay has cells whose plateau is 0.04 or 16 (too few pairs in dR in [1,4] to
    // define one at all), and normalizing by such a plateau blows the curve up and drags the
    // mean to O(100). Report, in addition: the MEDIAN, the mean over cells whose plateau is
    // sane (|plateau-1| <= the guard tolerance), and the INCLUSIVE cell -- which is the only
    // statistically meaningful curve for a 10k-event sample anyway.
    std::vector<double> chi2_all, chi2_sane;
    double chi2_incl = -1.;

    // iy/iz = 0 is the inclusive cell, fitted first.
    for (int iy = 0; iy <= npt; ++iy) {
        for (int iz = 0; iz <= neta; ++iz) {
            const bool inclusive = (iy == 0 && iz == 0);
            if (!inclusive && (iy == 0 || iz == 0)) continue;   // only the full inclusive cell

            const double plateau = inclusive ? plat_incl : hplat->GetBinContent(iy, iz);
            const int    pnb     = inclusive ? 1         : (int)hpnb->GetBinContent(iy, iz);
            const std::string nm = CellName("r", step, iy, iz);

            const char* ptlab  = inclusive ? "inclusive"
                : Form("[%.0f,%.0f)", hplat->GetXaxis()->GetBinLowEdge(iy),
                                      hplat->GetXaxis()->GetBinUpEdge(iy));
            const char* etalab = inclusive ? "inclusive"
                : Form("[%.1f,%.1f)", hplat->GetYaxis()->GetBinLowEdge(iz),
                                      hplat->GetYaxis()->GetBinUpEdge(iz));

            if (pnb <= 0 || plateau <= 0.) {
                if (!inclusive) hstat->SetBinContent(iy, iz, 0.);
                rep << std::left << std::setw(22) << ptlab << std::setw(16) << etalab
                    << std::setw(12) << "--" << std::setw(10) << 0
                    << "  no plateau -> cell skipped\n";
                continue;
            }

            // measured curve, normalized to 1 at large dR
            TH1D* r = DrCellRatio(h3n, h3d, h3a, h3b, h3p, h3q, iy, iz, nm.c_str());
            r->Scale(1.0 / plateau);   // scales contents AND errors

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
                if (!inclusive) hstat->SetBinContent(iy, iz, 0.);
                rep << std::left << std::setw(22) << ptlab << std::setw(16) << etalab
                    << std::setw(12) << Form("%.4f", plateau) << std::setw(10) << k
                    << "  too few points -> no fit\n";
                delete g;
                continue;
            }

            if (M.formula.empty()) {
                // ---- interpolation: knots below Rp, then a hard 1 --------------------------
                // No free parameters and no chi2: it passes through every point, so it also
                // transmits every point's statistical fluctuation straight into the correction.
                // That transmitted noise is what `knot_rel_err` records.
                auto* gk = new TGraph();
                int kk = 0;
                bool first_knot = true;
                for (int i = 0; i < g->GetN(); ++i) {
                    double x, y;
                    g->GetPoint(i, x, y);
                    if (x >= S.flat_onset) continue;
                    // Pin dR = 0 to the first measured value: TGraph::Eval extrapolates
                    // LINEARLY below the first knot, which would invent a dR -> 0 trend the
                    // data does not contain (the first bin centre is 0.025).
                    if (first_knot) { gk->SetPoint(kk++, 0.0, y); first_knot = false; }
                    gk->SetPoint(kk++, x, y);
                }
                // pin the flat branch: value 1 from Rp out to the generous range end
                for (double x : {S.flat_onset, 1.0, 2.0, 5.0, kTF1RangeHi})
                    if (x >= S.flat_onset) gk->SetPoint(kk++, x, 1.0);
                gk->SetName(CellName("gknots", step, iy, iz).c_str());
                gk->SetTitle(Form("linear-interpolation knots; #DeltaR; %s (plateau-normalized)",
                                  S.quantity.c_str()));
                gk->Write();
                if (!inclusive) {
                    hstat->SetBinContent(iy, iz, 1.);
                    hchi ->SetBinContent(iy, iz, -1.);     // n/a for an interpolation
                    hknot->SetBinContent(iy, iz, knot_rel);
                    hf0  ->SetBinContent(iy, iz, gk->Eval(0.0));
                }
                rep << std::left << std::setw(22) << ptlab << std::setw(16) << etalab
                    << std::setw(12) << Form("%.4f", plateau) << std::setw(10) << k
                    << std::setw(14) << "--" << std::setw(14) << "--" << std::setw(14) << "--"
                    << std::setw(12) << "n/a"
                    << std::setw(12) << Form("%.4f", gk->Eval(0.0))
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
            const double A0 = y0 - 1.0;
            if (method == "expo") {
                f->SetParNames("A", "#lambda", "p");
                f->SetParameters(A0 != 0. ? A0 : 0.2, 0.25, 1.5);
                f->SetParLimits(0, -5.0, 20.0);
                f->SetParLimits(1, 0.02, 3.0);
                f->SetParLimits(2, 0.3, 8.0);
            } else if (method == "polyu_fixedRp") {
                f->SetParNames("a_{2}", "a_{3}", "a_{4}", "R_{p}");
                f->SetParameters(A0 != 0. ? A0 : 0.2, 0.0, 0.0, S.flat_onset);
                f->SetParLimits(0, -50.0, 50.0);
                f->SetParLimits(1, -100.0, 100.0);
                f->SetParLimits(2, -100.0, 100.0);
                f->FixParameter(3, S.flat_onset);
            } else {
                f->SetParNames("A", "n", "R_{p}");
                f->SetParameters(A0 != 0. ? A0 : 0.2, 2.0, S.flat_onset);
                f->SetParLimits(0, -5.0, 20.0);
                f->SetParLimits(1, 0.2, 20.0);
                if (M.float_rp) f->SetParLimits(2, 0.6 * S.flat_onset, 1.6 * S.flat_onset);
                else            f->FixParameter(2, S.flat_onset);
            }

            TFitResultPtr fr = g->Fit(f, "QRNS");
            const bool ok = (fr.Get() != nullptr) && fr->IsValid() && fr->Ndf() > 0;
            const double chi2ndf = ok ? fr->Chi2() / fr->Ndf() : -1.;

            // Store with a GENEROUS range: a TFormula TF1 evaluates its formula everywhere, but
            // a wide range keeps Draw()/Integral() honest too. The CONSUMER still clamps.
            f->SetRange(kFitLo, kTF1RangeHi);
            f->Write();
            g->Write();

            if (!inclusive) {
                hstat->SetBinContent(iy, iz, ok ? 1. : 0.);
                hchi ->SetBinContent(iy, iz, chi2ndf);
                for (int ip = 0; ip < M.npar && ip < 4; ++ip) {
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
                    if (std::fabs(plateau - 1.0) <= kPlateauGuardTol) chi2_sane.push_back(chi2ndf);
                }
            }

            rep << std::left << std::setw(22) << ptlab << std::setw(16) << etalab
                << std::setw(12) << Form("%.4f", plateau) << std::setw(10) << k;
            for (int ip = 0; ip < 3; ++ip)
                rep << std::setw(14) << (ip < M.npar ? Form("%.4f", f->GetParameter(ip)) : "--");
            rep << std::setw(12) << (ok ? Form("%.3f", chi2ndf) : "FAILED")
                << std::setw(12) << Form("%.4f", f->Eval(0.0))
                << (M.npar > 3 ? Form("  p3=%.4f", f->GetParameter(3)) : "") << "\n";
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
    }
    TNamed("provenance",
           Form("sample=%s (%s); WP=%s; step=%d (%s); method=%s; formula=%s; Rp=%.2f; "
                "fit range dR=[%.2f,%.2f]; stored TF1 range=[0,%.1f]; plateau source=%s; "
                "guard=%s; producer=fit_dr_corrections.cxx",
                sample.c_str(), cfg.is_full_sample ? "FULL" : "TEST", wp_text.c_str(), step,
                S.quantity.c_str(), method.c_str(),
                M.formula.empty() ? "linear interpolation (TGraph knots)" : M.formula.c_str(),
                S.flat_onset, kFitLo, kFitHi, kTF1RangeHi, plateau_path.c_str(),
                violations.empty() ? "PASS"
                                   : (cfg.is_full_sample ? "FAIL(overridden)" : "TEST-exempt")))
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
        << "# chi2/ndf over cells with |plateau-1| <= " << kPlateauGuardTol << ": "
        << line_sane << "\n"
        << "# (cells whose plateau is far from 1 have too few dR-in-[1,4] pairs to define one;"
           " normalizing by such a plateau inflates chi2 without saying anything about the fit"
           " function.)\n"
        // kept for backwards compatibility with the driver's summary parser
        << "# mean chi2/ndf over " << chi2_all.size() << " converged cells = "
        << (chi2_all.empty() ? "n/a"
                             : Form("%.3f", std::accumulate(chi2_all.begin(), chi2_all.end(), 0.)
                                            / chi2_all.size()))
        << "\n";
    rep.close();

    std::cout << "  chi2/ndf inclusive = " << line_incl << "\n"
              << "  chi2/ndf all cells: " << line_all << "\n"
              << "  chi2/ndf sane-plateau cells: " << line_sane << "\n"
              << "  wrote " << out_path << "\n"
              << "  wrote " << mdir << "fit_report.txt\n";
}

// Convenience: every method for one (sample, WP, step).
void fit_dr_corrections_all(const std::string& sample = "pp_full", bool use_tight_wp = true,
                            int step = 3, bool allow_plateau_violation = false)
{
    for (const std::string& m : {"powerlaw_fixedRp", "powerlaw_floatRp", "expo",
                                 "polyu_fixedRp", "interp"})
        fit_dr_corrections(sample, use_tight_wp, step, m, allow_plateau_violation);
}
