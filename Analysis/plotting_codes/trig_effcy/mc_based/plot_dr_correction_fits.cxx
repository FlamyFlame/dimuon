// plot_dr_correction_fits.cxx
//
// PLOT + READ-BACK-VERIFICATION stage of the DeltaR-correction chain
// (mc_trigger_efficiency.md §3.3 Step 3 / §3.4 Step 4, round-7 Autonomy-Contract item 6).
//
//   plot_mc_trig_eff.cxx        measures eps_dR / eps_single, writes the plateau ROOT file
//   fit_dr_corrections.cxx      guards + normalizes + fits, writes the fit ROOT file
//   THIS MACRO                  re-opens the FIT FILE and draws it over the measured points
//
// WHY IT IS A SEPARATE MACRO (and a separate process). The persisted TF1s are what the analysis
// will actually apply, so they must be verified AS READ BACK, not as they were in the fitting
// job's memory. This repo has already been burned by exactly that: compiled/lambda-based TF1s
// return 0 outside their range on read-back, which floored a single-muon efficiency and produced
// a spurious pp cross-section jump (pp_trig_eff_highpt_jump.md). Every fitted curve drawn here
// comes from TF1::Eval / TGraph::Eval of the object read out of the file, evaluated at the EXACT
// DeltaR (never resampled to the nearest bin), and the read-back audit below additionally probes
// each function well outside its fit range and fails loudly if it collapses to 0 or drifts off 1.
//
// LAYOUT (user requirement): 1 PNG per pair-pT bin, 1 subplot per pair-eta bin (9 bins -> 3x3,
// nrows >= ncols ~ sqrt(N)); plus one PNG for the inclusive cell.
//
// TWO dR VIEWS of the same cells (user, 2026-08-12), STEP 3 ONLY -- see `dr_view` and kXhiFull:
//   <sign mode>/          DEFAULT. dR in [0, 1] = the FIT DOMAIN, the window the correction is
//                         applied in. No point outside it, so no open markers and no
//                         "outside the fit range" legend entry.
//   <sign mode>/dR0_2/    REFERENCE. dR in [0, 2]: the fit domain plus the large-dR points that
//                         DEFINED the plateau (drawn open), kept so that check stays available.
// Step 4 was NOT restructured (user asked for step3_dr_fit only): it keeps one flat directory
// holding the dR in [0, 2] view, and the pipeline calls it with dr_view = "dr0_2".
// In sign_sepr, the same-sign/opposite-sign OVERLAY is the main figure and stays at the top of
// its view directory; the same-sign/opposite-sign RATIO canvases go one level down in ratio/.
//
// TWO SERIES MODES (user, 2026-08-11) -- one subdirectory per mode under <method>/:
//   sign_intgr/  the sign-INTEGRATED correction, i.e. the nominal one the analysis applies:
//                measured values black with error bars, fitted function red.
//   sign_sepr/   same cells, but the SAME-SIGN and the OPPOSITE-SIGN series and both of their
//                fitted curves overlaid on every subplot -- blue for same sign, red for opposite
//                sign; the fitted curves in kBlue/kRed and the markers + error bars one shade
//                darker (kBlue+1/kRed+1), with TRIANGLES for same sign against CIRCLES for
//                opposite sign so the two series stay separable in print and in greyscale.
//                This is a physics comparison
//                (do the two charge combinations need the same dR correction?), so both series
//                must be on the SAME axes, in the same cells, with both sets of fitted parameters
//                readable on the panel.
// TWO PLATEAU MODES (user, 2026-08-11) -- the TOP level of the tree, above <method>/:
//   plateau_corrected/       NOMINAL. Every measured point is divided by its cell's large-dR
//                            plateau and the fitted shape tends to 1; the y axis is the RATIO
//                            of the efficiency to that plateau.
//   no_plateau_correction/   The RAW, un-normalized efficiency, fitted with a FREE baseline C.
//                            The y axis is the efficiency itself -- deliberately a DIFFERENT
//                            axis title, because it is a different quantity, and the [2, 3.5]
//                            plateau is not used anywhere on these canvases. C is drawn per
//                            panel with the other fitted parameters (for the interpolation,
//                            which has no fitted parameters, the pinned flat-branch value is
//                            drawn instead), so the reader can see what the fit decided the
//                            baseline is in every cell.
//
// Only PNGs live in those two subdirectories; the text artefacts (fit_report*.txt from the fit
// stage, readback_check*.txt from this one) stay directly in <method>/ where they have always been.
//
// The measured points are the SAME plateau-normalized values the fit saw: fine zoom bins for
// dR < 1 (the fit domain) and, beyond dR = 1, the coarse full-range bins -- those are the ones
// that DEFINED the plateau, so they must sit at 1 by construction and are shown as a check, not
// as fit input. Filled markers = inside the fit domain, open markers = outside it.
//
// Compile/run (ACLiC, from this directory):
//   root -l -b -q -e '.L plot_dr_correction_fits.cxx+'                       // compile only
//   root -l -b -q 'plot_dr_correction_fits.cxx+("pp_full", true, 3, "polyu_fixedRp")'
//   root -l -b -q 'plot_dr_correction_fits.cxx+("pp_full", true, 3, "expo", "sign_sepr")'
//   root -l -b -q 'plot_dr_correction_fits.cxx+("pp_full", true, 3, "expo", "sign_sepr", "corr", "dr0_2")'

#include <TArrow.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>

#include "dr_correction_sample_cfg.h"
#include "dr_correction_ratio.h"
#include "../../../Utilities/MCTrigEffPlateauWindow.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double kFitHi   = 1.0;    // fit domain upper edge (must match fit_dr_corrections.cxx)
// DRAWN dR WINDOW -- two views of the SAME fitted cells (user, 2026-08-12), selected by `dr_view`:
//   "dr0_1"  DEFAULT. The FIT DOMAIN itself, dR in [0, kFitHi]. This is the region the correction
//            is actually applied in, so it is the main figure. It contains no point the fit did
//            not see -- hence no open markers and no "outside the fit range" legend entry.
//   "dr0_2"  REFERENCE, kept one level down in dR0_2/. The fit domain PLUS the large-dR points
//            that DEFINED the plateau, so the check that they sit at 1 stays available.
// The two views share every number: same histograms, same fits, same acceptance screens. Only the
// drawn x window, and with it the y range computed from the drawn points, differ.
constexpr double kXhiFull = 2.0;
// The interpolation's pinned flat branch is a property of the FIT, not of the drawn window, so it
// is probed at a FIXED dR inside that branch. Probing it at the view's upper edge would make the
// printed baseline -- and the read-back audit's asymptote -- change between two views of one fit.
constexpr double kFlatProbeX = 2.0;
constexpr double kYcapHi  = 3.0;    // hard cap on the auto y range; off-scale points are arrowed
constexpr double kYcapLo  = 0.0;

// Human-readable LaTeX form of each fitted function, drawn ABOVE the parameter values so the
// reader never has to know what "polyu_fixedRp" or a bare "A"/"n"/"R_{p}" means. These MUST stay
// in step with the TFormula strings in fit_dr_corrections.cxx::MakeMethodCfg.
//   powerlaw_fixedRp / powerlaw_floatRp : 1+[0]*pow(max(0,1-x/[2]),[1])
//   expo                                : 1+[0]*exp(-pow(x/[1],[2]))
//   polyu_fixedRp                       : 1+u^2*([0]+[1]*u+[2]*u^2), u = max(0,1-x/[3])
//   interp                              : linear interpolation, 1 above Rp
// Returns {formula line, definition line} -- the second may be empty.
// `nocorr` = the no-plateau-correction mode, where the leading 1 of every shape is the FREE
// fitted baseline C (fit_dr_corrections.cxx). The equation MUST show it: it is the parameter that
// replaces the external normalization, and it is drawn per panel with the other parameters.
inline std::pair<std::string, std::string> MethodFormulaTex(const std::string& method,
                                                            bool nocorr = false)
{
    const std::string base = nocorr ? "C" : "1";
    if (method == "powerlaw_fixedRp" || method == "powerlaw_floatRp")
        return {"f(#DeltaR) = " + base + " + A #upoint u^{n}",
                "u #equiv max(0, 1 - #DeltaR/R_{p})"};
    if (method == "expo")
        return {"f(#DeltaR) = " + base + " + A #upoint exp[-(#DeltaR/#lambda)^{p}]", ""};
    if (method == "polyu_fixedRp")
        return {"f(#DeltaR) = " + base + " + u^{2}(a_{2} + a_{3}u + a_{4}u^{2})",
                "u #equiv max(0, 1 - #DeltaR/R_{p})"};
    if (method == "interp")
        // P2: a real piecewise DEFINITION with R_p defined and its value drawn, matching the
        // treatment polyu_fixedRp already gets. "linear interpolation of the points" is prose,
        // not an equation, and R_p was never defined anywhere on the canvas.
        return {"f(#DeltaR) = piecewise-linear through the measured points for #DeltaR < R_{p}",
                "f(#DeltaR) = " + base + " for #DeltaR #geq R_{p}"};
    return {"", ""};
}

template <typename T>
T* GetObj(TFile* f, const std::string& name)
{
    T* o = dynamic_cast<T*>(f->Get(name.c_str()));
    if (!o) throw std::runtime_error("plot_dr_correction_fits: missing '" + name + "' in "
                                     + f->GetName());
    return o;
}

TFile* OpenRead(const std::string& path)
{
    TFile* f = TFile::Open(path.c_str(), "READ");
    if (!f || f->IsZombie())
        throw std::runtime_error("plot_dr_correction_fits: cannot open " + path);
    return f;
}

std::string CellSuffix(int iy, int iz)
{
    return (iy == 0 && iz == 0) ? "_incl" : Form("_pt%d_eta%d", iy, iz);
}

// One evaluator for both the parametric methods (TF1) and the interpolation (TGraph knots), so
// the drawing code -- and any future consumer -- never branches on the method.
struct FittedFunc {
    TF1*    f = nullptr;
    TGraph* g = nullptr;
    bool valid() const { return f || g; }
    double Eval(double x) const { return f ? f->Eval(x) : (g ? g->Eval(x) : 1.0); }
};

FittedFunc LoadFunc(TFile* fit, int step, int iy, int iz)
{
    FittedFunc F;
    const std::string suf = CellSuffix(iy, iz);
    F.f = dynamic_cast<TF1*>(fit->Get(Form("f_step%d%s", step, suf.c_str())));
    if (!F.f) F.g = dynamic_cast<TGraph*>(fit->Get(Form("gknots_step%d%s", step, suf.c_str())));
    return F;
}

// ONE drawn series = one charge combination (or the sign-integrated sum). It owns its measured
// 3D histograms, its plateau map and its fit file; the keys INSIDE a fit file are identical for
// every series (the sign lives in the file NAME), so nothing below branches on the sign.
struct Series {
    std::string sign;      // "" | "ss" | "os"
    std::string legend;    // spelled-out physics wording drawn on the canvas
    Color_t     c_curve;   // fitted function
    Color_t     c_mark;    // measured points inside the fit domain
    Color_t     c_exc;     // measured points outside the fit domain
    Style_t     m_fit;     // filled marker, inside the fit domain
    Style_t     m_exc;     // open marker, outside it
    TFile* ffit = nullptr;
    TH3D  *zn = nullptr, *zd = nullptr, *za = nullptr, *zb = nullptr, *zp = nullptr, *zq = nullptr;
    TH3D  *fn = nullptr, *fd = nullptr, *fa = nullptr, *fb = nullptr, *fp = nullptr, *fq = nullptr;
    TH2D  *hplat = nullptr, *hchi = nullptr, *hok = nullptr;
    TH1D  *hpinc = nullptr, *hokinc = nullptr;
};

void DrawUnity(double xlo, double xhi)
{
    auto* l = new TLine(xlo, 1.0, xhi, 1.0);
    l->SetLineStyle(2);
    l->SetLineColor(kGray + 2);
    l->Draw("same");
}

// A parameter the fitter held FIXED is a constraint, not a measurement: no uncertainty, and the
// same value in every series. Detected exactly as the annotation has always detected it.
bool ParIsFixed(TF1* f, int ip)
{
    double plo = 0., phi = 0.;
    f->GetParLimits(ip, plo, phi);
    return (f->GetParError(ip) == 0.) || (plo == phi && plo != 0.);
}

// A FREE parameter that came out sitting ON one of its limits is not a measurement: MINUIT parks
// it on the boundary and still returns a parabolic error, which then spans the limit
// ("#lambda = 0.02 #pm 0.058" reads exactly like a measured value). The limits survive the
// write/read, so the test is made on the PERSISTED function -- identical to the one
// fit_dr_corrections.cxx uses when it counts them in the fit report.
bool ParAtLimit(TF1* f, int ip)
{
    double lo = 0., hi = 0.;
    f->GetParLimits(ip, lo, hi);
    if (!(lo < hi)) return false;          // unbounded, or FixParameter (which sets lo == hi)
    const double tol = 1e-3 * (hi - lo);
    const double v   = f->GetParameter(ip);
    return std::fabs(v - lo) <= tol || std::fabs(v - hi) <= tol;
}

}  // namespace

// step   : 3 (cross / 2mu4 term) or 4 (single leg)
// method : powerlaw_fixedRp | powerlaw_floatRp | expo | polyu_fixedRp | interp
// mode   : "sign_intgr" (the nominal, sign-integrated correction) | "sign_sepr" (same sign and
//          opposite sign overlaid). One output subdirectory each, under <method>/.
// plateau_mode : "corr" (NOMINAL: points divided by the plateau) | "nocorr" (raw efficiency,
//          free fitted baseline C). One output subdirectory each, ABOVE <method>/.
// dr_view : "dr0_1" (DEFAULT, the fit domain -- the main figure) | "dr0_2" (the reference view
//          out to dR = 2, written to dR0_2/). STEP 3 ONLY: step 4 was not restructured (user), so
//          it keeps its flat layout and must be called with "dr0_2" to reproduce its plots.
void plot_dr_correction_fits(const std::string& sample = "pp_full", bool use_tight_wp = true,
                             int step = 3, const std::string& method = "polyu_fixedRp",
                             const std::string& mode = "sign_intgr",
                             const std::string& plateau_mode = "corr",
                             const std::string& dr_view = "dr0_1")
{
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    if (mode != "sign_intgr" && mode != "sign_sepr")
        throw std::runtime_error("plot_dr_correction_fits: mode must be 'sign_intgr' or "
                                 "'sign_sepr', got '" + mode + "'");
    if (dr_view != "dr0_1" && dr_view != "dr0_2")
        throw std::runtime_error("plot_dr_correction_fits: dr_view must be 'dr0_1' or 'dr0_2', "
                                 "got '" + dr_view + "'");
    // The drawn window. Everything downstream -- the point loading, the frame, the sampled curve,
    // the y range, the ratio canvases -- reads it from here, so the two views cannot drift apart.
    const double kXhi = (dr_view == "dr0_2") ? kXhiFull : kFitHi;
    // Points beyond the fit domain exist only in the reference view.
    const bool draw_check_region = (kXhi > kFitHi);
    const bool sepr = (mode == "sign_sepr");
    // THE mode switch. DrCorrPlateauModeDir validates the token (it throws on anything else).
    const bool        nocorr   = (plateau_mode == "nocorr");
    const std::string mode_dir = DrCorrPlateauModeDir(plateau_mode);

    const DrCorrSample cfg = GetDrCorrSample(sample, use_tight_wp);
    const std::string wp_suf  = DrCorrWpSuffix(use_tight_wp);
    const std::string wp_text = use_tight_wp ? "Tight muons" : "Medium muons";
    const std::string wp_dir  = DrCorrWpDir(use_tight_wp);
    const std::string h_base  = (step == 3) ? "h_mc_dr_" : "h_mc_single_dr_";
    const bool has_cov = (step == 4);
    const std::string quantity_tex = (step == 3)
        ? cfg.eps_dr_text                                   // eps_dR^2mu4 (pp) / ^cross (PbPb)
        : std::string("#varepsilon_{#DeltaR}^{single}");
    const std::string tag = "step" + std::to_string(step);

    // out_base carries the variant tag (mc_based / _medium / _pt4bin / _pt4bin_medium); the mode
    // subdirectory separates the sign-integrated plots from the sign-separated ones so the two
    // can never interleave under identical file names. Text reports stay in <method>/ (rdir).
    const std::string rdir = cfg.out_base + tag + "_dr_fit/" + mode_dir + method + "/" + wp_dir;
    // VIEW LEVEL. The fit-domain view is the MAIN figure and stays at the top of <sign mode>/;
    // the reference view goes one level down. Step 4 was deliberately left on its old flat layout
    // (user restructured step 3 only), so its dR0_2/ level is suppressed.
    const std::string view_dir = (step == 3 && dr_view == "dr0_2") ? "dR0_2/" : "";
    const std::string odir = rdir + mode + "/" + view_dir;
    // The same-sign/opposite-sign RATIO canvases get their own subdirectory (user, 2026-08-12):
    // the overlay of the two charge combinations is the main figure of sign_sepr, and one flat
    // directory of 18 files interleaved the two families under near-identical names.
    // GATED ON STEP 3, exactly like view_dir above: the user restructured step3_dr_fit only, and
    // an ungated ratio/ left step 4 with the new subdirectory AND its 9 top-level ratio PNGs from
    // the previous run -- two copies, identical today, silently diverging from the next run on.
    const std::string rodir = (step == 3) ? odir + "ratio/" : odir;
    // The directory is created only once the series have actually loaded -- a mode that is
    // skipped for lack of per-sign inputs must not leave an empty directory behind.

    // ---- inputs: measured histograms and the FITS ---------------------------------------------
    // The pair-pT-binning token belongs here exactly as it does in plot_mc_trig_eff.cxx and
    // fit_dr_corrections.cxx. It was MISSING until 2026-08-11, so a MCTRIGEFF_PAIRPT_4BIN=1 run
    // read the NOMINAL 8-bin histograms and drew them against the 4-bin fit file.
    const std::string hist_path = cfg.mc_dir + "mc_trig_eff_hists_" + cfg.mc_label + wp_suf
                                + MCTrigEffPairPt::FileSuffix() + "_" + tag + ".root";
    TFile* fh = OpenRead(hist_path);

    // Series definition. Sign-integrated: the historical black/blue/red scheme, unchanged.
    // Sign-separated (user, 2026-08-11, marker/colour revision 2026-08-12): opposite sign = kRed
    // curve with kRed+1 markers AND error bars, same sign = kBlue curve with kBlue+1 markers and
    // error bars, so a curve and its measurement read as one pair. The MARKER SHAPE, not only the
    // hue, carries the sign: same sign gets triangles (22 filled / 26 open), opposite sign keeps
    // circles (20 / 24) -- the earlier kBlue+3 vs kRed+3 pair was too dark to tell apart where the
    // two series overlap at small dR, which is exactly the region the comparison is about.
    std::vector<Series> series;
    if (!sepr) {
        series.push_back({"", "measurement", kRed + 1, kBlack, kBlue + 1, 20, 24});
    } else {
        series.push_back({"ss", DrCorrSignText("ss"), kBlue, kBlue + 1, kBlue + 1, 22, 26});
        series.push_back({"os", DrCorrSignText("os"), kRed,  kRed  + 1, kRed  + 1, 20, 24});
    }

    for (auto& s : series) {
        const std::string pref    = h_base + (s.sign.empty() ? "" : s.sign + "_");
        const std::string fitpath = DrCorrFitFile(cfg, use_tight_wp, step, method, s.sign,
                                                 plateau_mode);
        // GRACEFUL SKIP (non-fatal, by design): a production filled before the per-sign booking
        // has neither the per-sign histograms nor the per-sign fits. Say so and return, so the
        // sign-integrated plots of that same sample keep being produced.
        const bool have_hists = fh->Get((pref + "zoom_vs_pt_eta_num").c_str()) != nullptr;
        const bool have_fit   = (gSystem->AccessPathName(fitpath.c_str()) == false);
        if (!have_hists || !have_fit) {
            std::cout << "  SKIPPED " << mode << " for sample '" << sample << "': "
                      << (have_hists ? "" : "no " + pref + "zoom_vs_pt_eta_num in " + hist_path)
                      << (have_fit   ? "" : (have_hists ? "" : "; ") + std::string("no fit file ")
                                            + fitpath)
                      << ".\n  Re-fill the histograms with the per-sign booking and re-run the fit"
                         " stage for this series.\n";
            fh->Close();
            return;
        }
        s.ffit = OpenRead(fitpath);
        s.zn = GetObj<TH3D>(fh, pref + "zoom_vs_pt_eta_num");
        s.zd = GetObj<TH3D>(fh, pref + "zoom_vs_pt_eta_denom");
        s.za = GetObj<TH3D>(fh, pref + "zoom_vs_pt_eta_errA");
        s.zb = GetObj<TH3D>(fh, pref + "zoom_vs_pt_eta_errB");
        s.zp = has_cov ? GetObj<TH3D>(fh, pref + "zoom_vs_pt_eta_covP") : nullptr;
        s.zq = has_cov ? GetObj<TH3D>(fh, pref + "zoom_vs_pt_eta_covQ") : nullptr;
        s.fn = GetObj<TH3D>(fh, pref + "full_vs_pt_eta_num");
        s.fd = GetObj<TH3D>(fh, pref + "full_vs_pt_eta_denom");
        s.fa = GetObj<TH3D>(fh, pref + "full_vs_pt_eta_errA");
        s.fb = GetObj<TH3D>(fh, pref + "full_vs_pt_eta_errB");
        s.fp = has_cov ? GetObj<TH3D>(fh, pref + "full_vs_pt_eta_covP") : nullptr;
        s.fq = has_cov ? GetObj<TH3D>(fh, pref + "full_vs_pt_eta_covQ") : nullptr;
        // Keys inside a fit file are the same for every series (the sign is in the file name).
        s.hplat = GetObj<TH2D>(s.ffit, "h_" + tag + "_plateau");
        s.hchi  = GetObj<TH2D>(s.ffit, "h_" + tag + "_chi2ndf");
        s.hpinc = GetObj<TH1D>(s.ffit, "h_" + tag + "_plateau_inclusive");
        // THE SCREEN THE ANALYSIS ITSELF APPLIES. fit_dr_corrections.cxx writes fit_ok = 0 for a
        // cell whose fit did not converge, whose plateau is outside the guard tolerance, or whose
        // fitted correction is not > 0 over [0, Rp] -- a NEGATIVE trigger correction. Consumers
        // must require fit_ok == 1, so a canvas that draws such a curve with its parameters and
        // chi2/ndf is showing numbers the analysis will never use. The inclusive cell has no bin
        // in the TH2D and carries its own flag (added 2026-08-11); a fit file written before that
        // has none, in which case the inclusive fit is drawn unscreened and it is said so on
        // stdout rather than silently assumed good.
        s.hok    = GetObj<TH2D>(s.ffit, "h_" + tag + "_fit_ok");
        s.hokinc = dynamic_cast<TH1D*>(s.ffit->Get(("h_" + tag + "_fit_ok_inclusive").c_str()));
        if (!s.hokinc)
            std::cout << "  NOTE: no h_" << tag << "_fit_ok_inclusive in "
                      << s.ffit->GetName() << " -- the inclusive fit is drawn UNSCREENED."
                         " Re-run the fit stage.\n";
    }

    gSystem->mkdir(odir.c_str(), kTRUE);
    if (sepr) gSystem->mkdir(rodir.c_str(), kTRUE);

    const int npt  = series[0].hplat->GetNbinsX();
    const int neta = series[0].hplat->GetNbinsY();
    for (const auto& s : series)
        if (s.hplat->GetNbinsX() != npt || s.hplat->GetNbinsY() != neta)
            throw std::runtime_error("plot_dr_correction_fits: the overlaid series describe "
                                     "different (pair pT, pair eta) cells -- stale fit file?");

    // Flat onset R_p, READ FROM THE FIT FILE's provenance stamp (fit_dr_corrections.cxx writes it
    // from kFlatOnsetStep{3,4}) rather than retyped here. Needed on the canvas for the methods
    // whose equation mentions R_p but which have no R_p among the drawn parameters -- the
    // interpolation, whose piecewise definition is meaningless without the value of R_p.
    double rp_prov = -1.;
    if (auto* pr = dynamic_cast<TNamed*>(series[0].ffit->Get("provenance"))) {
        const TString t = pr->GetTitle();
        const Ssiz_t i = t.Index("Rp=");
        if (i >= 0) rp_prov = TString(t(i + 3, 8)).Atof();
    }

    auto plateau_of = [&](const Series& s, int iy, int iz) {
        return (iy == 0 && iz == 0) ? s.hpinc->GetBinContent(1) : s.hplat->GetBinContent(iy, iz);
    };
    auto plateau_err_of = [&](const Series& s, int iy, int iz) {
        return (iy == 0 && iz == 0) ? s.hpinc->GetBinError(1) : s.hplat->GetBinError(iy, iz);
    };
    // 1 = the analysis accepts this fit, 0 = it rejects it, -1 = the file predates the flag.
    auto fit_ok_of = [&](const Series& s, int iy, int iz) -> int {
        if (iy == 0 && iz == 0)
            return s.hokinc ? (s.hokinc->GetBinContent(1) > 0.5 ? 1 : 0) : -1;
        return s.hok->GetBinContent(iy, iz) > 0.5 ? 1 : 0;
    };

    // The DIVISOR applied to every measured point: the cell's plateau in the nominal mode, and
    // exactly 1 in the no-correction mode, where the baseline is a fitted parameter instead.
    auto norm_of = [&](const Series& s, int iy, int iz) {
        return nocorr ? 1.0 : plateau_of(s, iy, iz);
    };
    // The unmeasurable-plateau screen is a statement about a DIVISOR. Nothing is divided in the
    // no-correction mode, so no cell is screened out there -- it is drawn and fitted from its
    // own dR < 1 points, which is the entire point of the mode.
    auto cell_drawable = [&](const Series& s, int iy, int iz) {
        return nocorr || DrCorrPlateauUsable(plateau_of(s, iy, iz), plateau_err_of(s, iy, iz));
    };

    // ---- the measured points of one cell (fit input + large-dR check region) ------------------
    auto cell_points = [&](const Series& s, int iy, int iz, double plateau) -> TGraphErrors* {
        auto* g = new TGraphErrors();
        int k = 0;
        TH1D* rz = DrCellRatio(s.zn, s.zd, s.za, s.zb, s.zp, s.zq, iy, iz,
                               Form("pz_%s_%d_%d", s.sign.c_str(), iy, iz));
        for (int i = 1; i <= rz->GetNbinsX(); ++i) {
            const double e = rz->GetBinError(i);
            if (e <= 0.) continue;
            g->SetPoint(k, rz->GetBinCenter(i), rz->GetBinContent(i) / plateau);
            g->SetPointError(k, 0., e / plateau);
            ++k;
        }
        delete rz;
        TH1D* rf = DrCellRatio(s.fn, s.fd, s.fa, s.fb, s.fp, s.fq, iy, iz,
                               Form("pf_%s_%d_%d", s.sign.c_str(), iy, iz));
        for (int i = 1; i <= rf->GetNbinsX(); ++i) {
            const double x = rf->GetBinCenter(i), e = rf->GetBinError(i);
            if (x <= kFitHi || x > kXhi || e <= 0.) continue;   // only the dR > 1 check region
            g->SetPoint(k, x, rf->GetBinContent(i) / plateau);
            g->SetPointError(k, 0., e / plateau);
            ++k;
        }
        delete rf;
        g->SetMarkerStyle(s.m_fit);
        g->SetMarkerSize(0.7);
        g->SetMarkerColor(s.c_mark);
        g->SetLineColor(s.c_mark);
        g->SetLineWidth(1);
        return g;
    };

    // ---- WHAT IS IDENTICAL IN EVERY PANEL (drawn ONCE, in the canvas header strip) ---------
    // The fitted equation, its auxiliary definition, the parameters the method HOLDS FIXED and
    // R_p do not change from cell to cell. Repeating all of it inside all nine panels is what
    // forced the annotation block down into the data; the convention's answer to "no quadrant is
    // free" is to give the text its own space, so the shared part moves to the header strip and
    // each panel keeps only its own numbers.
    const auto ftex = MethodFormulaTex(method, nocorr);
    std::vector<int> free_par;             // drawn per cell, per series
    std::vector<std::string> fixed_line;   // drawn once, shared
    bool has_tf1 = false;
    {
        std::vector<FittedFunc> Fr(series.size());
        for (size_t is = 0; is < series.size(); ++is) Fr[is] = LoadFunc(series[is].ffit, step, 0, 0);
        has_tf1 = (Fr[0].f != nullptr);
        if (Fr[0].f) {
            for (int ip = 0; ip < Fr[0].f->GetNpar(); ++ip) {
                bool shared_fixed = ParIsFixed(Fr[0].f, ip);
                for (size_t is = 1; is < series.size() && shared_fixed; ++is)
                    shared_fixed = Fr[is].f && ParIsFixed(Fr[is].f, ip)
                                && Fr[is].f->GetParameter(ip) == Fr[0].f->GetParameter(ip);
                if (shared_fixed)
                    fixed_line.push_back(Form("%s = %.3g (fixed)", Fr[0].f->GetParName(ip),
                                              Fr[0].f->GetParameter(ip)));
                else
                    free_par.push_back(ip);
            }
        }
    }
    // The equation quotes R_p; if no drawn parameter carries it (the interpolation has no TF1 at
    // all), its VALUE comes from the fit file's provenance stamp, never a number retyped here.
    bool rp_drawn = false;
    for (const auto& fl : fixed_line) if (fl.find("R_{p}") != std::string::npos) rp_drawn = true;
    {
        const FittedFunc Fr0 = LoadFunc(series[0].ffit, step, 0, 0);
        if (Fr0.f) for (int ip : free_par)
            if (std::string(Fr0.f->GetParName(ip)) == "R_{p}") rp_drawn = true;
    }
    const bool rp_ext_line = !rp_drawn && rp_prov > 0.
                          && (ftex.first.find("R_{p}") != std::string::npos
                           || ftex.second.find("R_{p}") != std::string::npos);

    // Geometry of the RESERVED annotation band. Its height follows from the number of lines the
    // columns actually carry, and the shared y range below is stretched until no measured point
    // can reach into it -- reserving space, not painting a box over the data.
    // RESERVED STRIP ABOVE THE FRAME -- not a band stolen from the top of the y range.
    // History, because the two wrong answers were both tried here: (1) drawing the per-cell
    // numbers INSIDE the frame put "plateau = 1.0043" through the fitted curve, the unity line
    // and two measured points; (2) reserving the space by STRETCHING the y range until no point
    // could reach the text worked only while the drawn window was dR in [0, 2], where the text
    // sat over the flat right half. In the dR in [0, 1] view the same columns cover dR > 0.34,
    // i.e. real structure in every cell, and the exact reservation inflated the shared range to
    // [0, 4.88] -- the figure became 76% empty to keep text off the data. Clamping that stretch
    // at the y cap then put error bars back through the text. The strip removes the trade-off:
    // the text can never collide with data because it is not in the frame, and the y range is
    // free to follow the points (capped) in BOTH views.
    constexpr double kLabRow  = 0.045;             // the pair-eta label, at the top of the strip
    constexpr double kStripTop = 0.965;            // NDC y of the label baseline
    const double kColStep = sepr ? 0.040 : 0.045;
    // rows per column: series name (sign_sepr only) + the plateau line (nominal mode only) +
    // one per free parameter + chi2/ndf -- or, for the interpolation in the no-correction mode,
    // the single pinned-baseline line that takes the place of both.
    const int n_ann_rows = std::max((sepr ? 1 : 0) + (nocorr ? 0 : 1) + (int)free_par.size()
                                    + (has_tf1 ? 1 : (nocorr ? 1 : 0)),
                                    (sepr ? 1 : 0) + 2);
    const double h_ann   = n_ann_rows * kColStep + 0.02;
    // Pad top margin = label + annotation rows + a little air. The frame starts below it.
    const double kTopMarg = 0.035 + kLabRow + h_ann;
    const double kAnnTop  = kStripTop - kLabRow;   // first annotation row, under the label
    // Single-series canvases start their one column further left in the no-correction mode: the
    // baseline line ("C (fitted plateau) = 0.979 #pm 0.011") is the longest string any panel
    // carries and at 0.55 its tail was clipped by the frame edge. The nominal mode keeps 0.55.
    const double xcol[2] = {sepr ? 0.42 : (nocorr ? 0.45 : 0.55), 0.71};

    // ---- one subplot ------------------------------------------------------------------------
    // Returns the list of points pushed off the (capped) frame, so the caller can list them on
    // the canvas: a hidden point must never look like a missing one.
    auto draw_cell = [&](int iy, int iz, double ylo, double yhi,
                         std::vector<std::string>& offscale) {
        gPad->SetLeftMargin(0.14);
        gPad->SetBottomMargin(0.14);
        // The strip: the pair-eta label AND this cell's fitted numbers live above the frame.
        gPad->SetTopMargin(kTopMarg);
        TH1* fr = gPad->DrawFrame(0.0, ylo, kXhi, yhi);
        fr->GetXaxis()->SetTitle("#DeltaR(#mu_{1}, #mu_{2})");
        // The two modes plot DIFFERENT quantities and must never carry the same axis title:
        // a ratio to the large-dR plateau, or the efficiency ratio itself.
        fr->GetYaxis()->SetTitle((nocorr ? quantity_tex : quantity_tex + " / plateau").c_str());
        fr->GetXaxis()->SetTitleSize(0.050);
        fr->GetYaxis()->SetTitleSize(0.050);
        fr->GetXaxis()->SetLabelSize(0.042);
        fr->GetYaxis()->SetLabelSize(0.042);
        fr->GetYaxis()->SetTitleOffset(1.25);
        DrawUnity(0.0, kXhi);

        // THE CELL LABEL COMES FIRST. It used to be drawn after the "no usable plateau" early
        // return, so a cell with no measurement came out as an unlabelled empty frame and the
        // reader could not tell WHICH cell it was.
        // WHICH panel a point belongs to must travel WITH the point: the off-scale record is
        // drawn once for the whole nine-panel canvas, so "#DeltaR=0.03: 4.92 #pm 4.92" on its own
        // does not say which pair-eta cell it came from (nor, in the sign-separated mode, which
        // charge combination). Empty on the single-panel inclusive canvas, where there is only
        // one cell and the panel already names it.
        const std::string cell_id = (iy == 0 && iz == 0)
            ? std::string()
            : std::string(Form("#eta^{pair} #in [%.1f,%.1f), ",
                               series[0].hplat->GetYaxis()->GetBinLowEdge(iz),
                               series[0].hplat->GetYaxis()->GetBinUpEdge(iz)));
        TLatex t;
        t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.045);
        if (iy == 0 && iz == 0) {
            t.DrawLatex(0.18, kStripTop, "inclusive (all p_{T}^{pair}, all #eta^{pair})");
        } else {
            t.DrawLatex(0.18, kStripTop, Form("%.1f < #eta^{pair} < %.1f",
                                         series[0].hplat->GetYaxis()->GetBinLowEdge(iz),
                                         series[0].hplat->GetYaxis()->GetBinUpEdge(iz)));
        }

        // Same screen as the fit stage (shared helper): an unmeasurable cell is NOT drawn --
        // drawing one means dividing the curve by a near-zero plateau and presenting the
        // resulting O(10) excursion as a correction.
        std::vector<bool> usable(series.size(), false);
        int n_usable = 0;
        for (size_t is = 0; is < series.size(); ++is) {
            usable[is] = cell_drawable(series[is], iy, iz);
            if (usable[is]) ++n_usable;
        }
        if (n_usable == 0) {
            t.SetTextSize(0.05);
            t.DrawLatex(0.20, 0.55, "no fit");
            return;
        }
        // ... and the screen the ANALYSIS applies to the fit itself. A rejected fit keeps its
        // measured points (the measurement is real) but loses its curve and its parameters.
        std::vector<bool> accepted(series.size(), false);
        for (size_t is = 0; is < series.size(); ++is)
            accepted[is] = usable[is] && fit_ok_of(series[is], iy, iz) != 0;

        // The fitted curves are drawn BEFORE the points: a 2px line over 0.7-size markers hides
        // the measurement, and the measurement is what the reader has to judge the fit against.
        std::vector<FittedFunc> F(series.size());
        for (size_t is = 0; is < series.size(); ++is) {
            if (!usable[is]) continue;
            F[is] = LoadFunc(series[is].ffit, step, iy, iz);
            if (!F[is].valid() || !accepted[is]) continue;
            auto* gf = new TGraph();
            const int NP = 500;
            for (int i = 0; i <= NP; ++i) {
                const double x = kXhi * i / NP;
                gf->SetPoint(i, x, F[is].Eval(x));   // continuous function at the EXACT dR
            }
            gf->SetLineColor(series[is].c_curve);
            gf->SetLineWidth(2);
            gf->Draw("L same");
        }

        // POINTS OUTSIDE THE FIT DOMAIN ARE DRAWN OPEN / IN A SECOND COLOUR (user, 2026-08-06),
        // so a human can see at a glance which points the fit was actually constrained by. Only
        // dR < kFitHi enters the fit; the points beyond it are the plateau-window bins, shown as
        // a check that the normalized curve really does sit at 1 there.
        for (size_t is = 0; is < series.size(); ++is) {
            if (!usable[is]) continue;
            const Series& s = series[is];
            const std::string sign_id = sepr ? s.legend + ", " : std::string();
            auto* g = cell_points(s, iy, iz, norm_of(s, iy, iz));
            auto* g_fit = new TGraphErrors();
            auto* g_exc = new TGraphErrors();
            for (int i = 0, kf = 0, ke = 0; i < g->GetN(); ++i) {
                double x, y;
                g->GetPoint(i, x, y);
                auto* dst = (x < kFitHi) ? g_fit : g_exc;
                int&  k   = (x < kFitHi) ? kf : ke;
                dst->SetPoint(k, x, y);
                dst->SetPointError(k, 0., g->GetErrorY(i));
                ++k;
            }
            g_fit->SetMarkerStyle(s.m_fit); g_fit->SetMarkerSize(0.7);
            g_fit->SetMarkerColor(s.c_mark); g_fit->SetLineColor(s.c_mark);
            g_exc->SetMarkerStyle(s.m_exc); g_exc->SetMarkerSize(0.7);
            g_exc->SetMarkerColor(s.c_exc);  g_exc->SetLineColor(s.c_exc);
            if (g_fit->GetN()) g_fit->Draw("PZ same");
            if (g_exc->GetN()) g_exc->Draw("PZ same");
            for (int i = 0; i < g->GetN(); ++i) {
                double x, y;
                g->GetPoint(i, x, y);
                if (y <= yhi && y >= ylo) continue;
                const bool up = (y > yhi);
                // An up-arrow under the reserved annotation band is drawn BELOW the band: the
                // band is text-only space, and an arrowhead behind a parameter value is exactly
                // the collision the reservation exists to prevent.
                // The whole frame is data space now (the numbers are in the strip above it), so
                // an up-arrow goes to the frame top wherever the point is.
                constexpr double ftop = 0.98;
                auto* ar = new TArrow(x, ylo + (up ? ftop - 0.10 : 0.12) * (yhi - ylo),
                                      x, ylo + (up ? ftop        : 0.02) * (yhi - ylo),
                                      0.008, "|>");
                const Color_t acol = (x < kFitHi) ? s.c_mark : s.c_exc;
                ar->SetLineColor(acol); ar->SetFillColor(acol); ar->SetLineWidth(2);
                ar->Draw();
                offscale.push_back(Form("%s%s#DeltaR=%.2f: %.2f #pm %.2f", cell_id.c_str(),
                                        sign_id.c_str(), x, y, g->GetErrorY(i)));
            }
            // `g` is intentionally NOT deleted: ~TGraph removes the object from every pad that
            // draws it, so deleting it here silently erased all the points from the canvas.
        }

        auto chi2_of = [&](size_t is) {
            if (!F[is].f) return -1.;
            return (iy == 0 && iz == 0)
                ? (F[is].f->GetNDF() > 0 ? F[is].f->GetChisquare() / F[is].f->GetNDF() : -1.)
                : series[is].hchi->GetBinContent(iy, iz);
        };

        // ---- annotation: ONLY the numbers that belong to THIS cell --------------------------
        // One column per series, inside the reserved band at the top of the frame. A fit the
        // analysis rejects prints "fit rejected" INSTEAD of its parameters and chi2/ndf: quoting
        // a value that will never be applied presents it as a measurement.
        TLatex q; q.SetNDC(); q.SetTextFont(42); q.SetTextSize(sepr ? 0.023 : 0.030);
        for (size_t is = 0; is < series.size() && is < 2; ++is) {
            double ty = kAnnTop;
            const double xc = xcol[is];
            q.SetTextColor(sepr ? series[is].c_mark : kBlack);
            if (sepr) { q.DrawLatex(xc, ty, series[is].legend.c_str()); ty -= kColStep; }
            if (!usable[is]) { q.DrawLatex(xc, ty, "no fit"); q.SetTextColor(kBlack); continue; }
            // The [2, 3.5] plateau is drawn ONLY where it is actually applied. In the
            // no-correction mode it is measured but unused, and a second number called "plateau"
            // sitting next to the fitted baseline C would be read as the thing the curve tends
            // to. It stays in the fit report instead.
            if (!nocorr) {
                q.DrawLatex(xc, ty, Form("plateau = %.4f", plateau_of(series[is], iy, iz)));
                ty -= kColStep;
            }
            // A cell can have a usable plateau and still carry no fitted function for ONE sign
            // (too few points with a non-zero error in that charge combination).
            if (!F[is].valid()) { q.DrawLatex(xc, ty, "no fit"); q.SetTextColor(kBlack); continue; }
            if (!accepted[is]) {
                q.DrawLatex(xc, ty, "fit rejected");
                q.SetTextColor(kBlack);
                continue;
            }
            // The fitted baseline gets a NAMED label wherever the column is wide enough for it
            // (one series per panel); with two colour-coded columns side by side there is room
            // for the symbol only, and the header line next to the equation says what it is.
            // ("fitted" only where it really is fitted: the interpolation PINS its baseline to
            // the last measured point rather than fitting it.)
            const char* c_label = sepr ? "C"
                                : (method == "interp" ? "C (plateau)" : "C (fitted plateau)");
            if (F[is].f) {
                for (int ip : free_par) {
                    // A parameter pinned ON its limit is a constraint, not a measurement, and its
                    // error spans the limit -- say so instead of printing "0.02 +- 0.058".
                    const char* pn = F[is].f->GetParName(ip);
                    if (nocorr && std::string(pn) == "C") pn = c_label;
                    q.DrawLatex(xc, ty, ParAtLimit(F[is].f, ip)
                        ? Form("%s = %.3g (at limit)", pn, F[is].f->GetParameter(ip))
                        : Form("%s = %.3g #pm %.2g", pn, F[is].f->GetParameter(ip),
                               F[is].f->GetParError(ip)));
                    ty -= kColStep;
                }
                const double c = chi2_of(is);
                if (c >= 0.) q.DrawLatex(xc, ty, Form("#chi^{2}/ndf = %.2f", c));
            } else if (nocorr && F[is].g) {
                // The interpolation has no fitted parameters, but it does have a BASELINE: the
                // value its flat branch was pinned to (the last measured point below R_p). The
                // asymptote must be readable off every panel here exactly as it is for the
                // parametric fits, since nothing external normalizes these curves.
                q.DrawLatex(xc, ty, Form("%s = %.4f", c_label, F[is].Eval(kFlatProbeX)));
            }
            q.SetTextColor(kBlack);
        }
    };

    // ---- legend entries, built once and reused on every canvas --------------------------------
    // Markers and fitted curve of one series are separate entries so both the marker shape and
    // the line colour are defined. Wording is the spelled-out physics -- bare SS/OS and
    // sign1/sign2 are forbidden on a canvas (.claude/conventions/atlas-plotting.md).
    const int n_leg = (int)series.size() * 2 + (draw_check_region ? 1 : 0);
    auto fill_legend = [&](TLegend* leg) {
        for (const auto& s : series) {
            auto* gd = new TGraphErrors();
            gd->SetMarkerStyle(s.m_fit); gd->SetMarkerColor(s.c_mark); gd->SetLineColor(s.c_mark);
            leg->AddEntry(gd, s.legend.c_str(), "lp");
        }
        for (const auto& s : series) {
            auto* gf = new TGraph();
            gf->SetLineColor(s.c_curve); gf->SetLineWidth(2);
            leg->AddEntry(gf, sepr ? ("fit, " + s.legend).c_str() : "fit", "l");
        }
        // The open markers were drawn but never identified. They are the bins BEYOND the fit
        // domain -- the large-dR check region the fit was not constrained by. In the fit-domain
        // view there are none, and a key for a marker the reader cannot find is worse than no key.
        if (!draw_check_region) return;
        auto* ge = new TGraphErrors();
        ge->SetMarkerStyle(24);
        ge->SetMarkerColor(sepr ? kGray + 2 : series[0].c_exc);
        ge->SetLineColor  (sepr ? kGray + 2 : series[0].c_exc);
        leg->AddEntry(ge, "outside the fit range", "lp");
    };

    // ---- the DEFINING equation of the plotted quantity (drawn once per canvas) ---------------
    // The symbol eps_dR^{2mu4} / ^{cross} / ^{single} is invented by this analysis, so the canvas
    // has to say what it MEANS -- the drawn f(dR) is the FIT function, not the definition. Wording
    // follows the Physics Procedure (mc_trigger_efficiency.md §3.3 for Step 3, §3.4 for Step 4).
    // NO NEGATIVE KERN around the conditioning bar: #kern[-0.20]{#DeltaR} pulled the #Delta on top
    // of the '|' on the INCLUSIVE canvas (900 px, smaller absolute font), rendering 'P(pair passes
    // 2mu4 |#DeltaR)' as one garbled glyph while the 1556 px 3x3 canvases looked fine. A kern is an
    // absolute-size-dependent nudge; the canvases here come in two sizes, so it cannot be tuned once.
    const bool is_2mu4 = cfg.eps_dr_text.find("2mu4") != std::string::npos;
    std::string def_line;
    if (step == 3)
        def_line = quantity_tex + "(#DeltaR) #equiv P("
                 + (is_2mu4 ? "pair passes 2mu4" : "both muons mu4-matched")
                 + " | #DeltaR) / (#varepsilon_{1}#varepsilon_{2}),   "
                   "#varepsilon_{i} = #varepsilon^{mu4}(p_{T}^{i}, q^{i}#eta^{i})";
    else
        def_line = quantity_tex + "(#DeltaR) #equiv P(leg fires the full mu4 chain "
                   "| leg of a reco pair at #DeltaR) "
                   "/ #varepsilon^{mu4}(p_{T}, q#eta)";
    // The large-dR plateau is DEFINED here only where the figure uses it. In the no-correction
    // mode nothing on the canvas is divided by it, so defining it would introduce a quantity the
    // reader then has to look for and never finds.
    if (!nocorr)
        def_line += Form(";   plateau = #LT%s#GT over #DeltaR #in [%.1f, %.1f]",
                         quantity_tex.c_str(), MCTrigEffPlateau::kLo, MCTrigEffPlateau::kHi);

    // The fit function and everything about it that is the same in all nine panels.
    const std::string eq_line1 = ftex.first;
    std::string eq_line2 = ftex.second;
    // C replaces the external normalization in the no-correction mode, so the canvas says what it
    // is once, next to the equation; its per-cell VALUE is drawn in every panel.
    if (nocorr)
        eq_line2 += (eq_line2.empty() ? "" : ",   ")
                  + std::string(method == "interp"
                        ? "C = the plateau, pinned to the last measured point"
                        : "C = the fitted plateau, free in the fit");
    for (const auto& fl : fixed_line) eq_line2 += (eq_line2.empty() ? "" : ",   ") + fl;
    if (rp_ext_line) eq_line2 += (eq_line2.empty() ? "" : std::string(",   "))
                               + Form("R_{p} = %.2f", rp_prov);
    eq_line2 += std::string(eq_line2.empty() ? "" : ",   ")
              + Form("fit range: #DeltaR < %.1f", kFitHi);

    // ---- y range: ONE range for the WHOLE METHOD DIRECTORY (user, 2026-08-05) ---------------
    // Common not just across the eta panels of one canvas, but across EVERY pair-pT bin -- i.e.
    // every PNG in this <step>_dr_fit/<method>/<mode>/ directory shares one y axis. Otherwise each
    // file silently rescales and the pair-pT dependence, which is the whole point of binning in
    // pair pT, cannot be read by flipping between the files. With two overlaid series the range
    // spans BOTH, so the two signs stay directly comparable panel by panel.
    //
    // The range is set from the CENTRAL VALUES only, NOT value +- error (user): in the sparse
    // high-pair-pT cells the error bars are several times the correction itself, and including
    // them let one noisy cell dictate the axis for every plot, squashing all the real structure
    // into a sliver. Points whose central value still falls outside the capped range are marked
    // with an arrow by draw_cell(), so nothing is dropped silently.
    // TWO ranges, deliberately different, because the two figure types answer different
    // questions:
    //   * the NINE-PANEL pair-pT files share ONE range across the whole method/mode directory,
    //     so the pair-pT dependence can be read by flipping between the files (user, 2026-08-05);
    //   * the single-panel INCLUSIVE canvas is a figure in its own right and is not compared to
    //     anything by flipping, so it gets a range from its OWN points. Inheriting the
    //     directory-wide range put the Step-4 inclusive points (0.99-1.16) inside a 0-2.4 frame,
    //     i.e. ~7% of the height, and turned the ~15% small-dR rise that IS the content of the
    //     figure into a barely visible kink.
    // Both are set from the CENTRAL VALUES only, NOT value +- error (user): in the sparse
    // high-pair-pT cells the error bars are several times the correction itself, and including
    // them let one noisy cell dictate the axis, squashing all the real structure into a sliver.
    // Points whose central value still falls outside the capped range are marked with an arrow
    // by draw_cell(), so nothing is dropped silently.
    auto compute_range = [&](bool incl_only) -> std::pair<double, double> {
        std::vector<std::pair<int, int>> cells;
        if (incl_only) cells.push_back({0, 0});
        else for (int iy = 1; iy <= npt; ++iy)
                 for (int iz = 1; iz <= neta; ++iz) cells.push_back({iy, iz});

        double ylo = 1., yhi = 1.;
        for (const auto& s : series) {
            for (const auto& c : cells) {
                // Screen with the SAME test the drawing code uses: an unmeasurable cell is never
                // drawn, so letting its eps/plateau values (O(10-130) when the plateau is ~0.008)
                // set the range pushed every PNG to the 3.0 cap and squashed the structure the
                // panels exist to show.
                if (!cell_drawable(s, c.first, c.second)) continue;
                auto* g = cell_points(s, c.first, c.second, norm_of(s, c.first, c.second));
                for (int i = 0; i < g->GetN(); ++i) {
                    double x, y;
                    g->GetPoint(i, x, y);
                    ylo = std::min(ylo, y);
                    yhi = std::max(yhi, y);
                }
                delete g;
            }
        }
        {
            const double span = yhi - ylo;
            ylo = std::max(kYcapLo, ylo - 0.05 * span);
            yhi = std::min(kYcapHi, yhi + 0.05 * span);
            // Minimum span, so a cell whose points all sit at ~1 does not get an absurdly zoomed
            // axis. It widens AROUND THE DATA, never to a fixed window: the fixed [0.9, 1.12] it
            // used to jump to sat below the Step-4 inclusive maximum (1.16) and would have
            // arrowed off scale the very rise the figure exists to show.
            if (yhi - ylo < 0.2) {
                const double mid = 0.5 * (ylo + yhi);
                ylo = std::max(kYcapLo, mid - 0.10);
                yhi = std::min(kYcapHi, mid + 0.10);
            }
        }
        // NO annotation-clearance stretch any more: the numbers sit in the strip above the frame
        // (see kTopMarg), so the frame is data space and the range is set by the data alone.
        return {ylo, yhi};
    };

    const std::pair<double, double> grid_rng = compute_range(false);
    const std::pair<double, double> incl_rng = compute_range(true);
    const double g_ylo = grid_rng.first, g_yhi = grid_rng.second;
    std::cout << "  common y range for every pair-pT PNG in this method/mode dir: ["
              << g_ylo << ", " << g_yhi << "]  (the fit annotation is in the strip above the "
                 "frame, so none of this range is reserved for it);  inclusive canvas: ["
              << incl_rng.first << ", " << incl_rng.second << "]\n";

    // ---- the canvas HEADER STRIP -------------------------------------------------------------
    // Everything the nine panels share: the sample/working point, the DEFINITION of the plotted
    // quantity, the fitted equation with its fixed parameters and fit range, any point pushed off
    // scale, and the legend. Sizes are given in PIXELS of this canvas and converted, so the strip
    // looks the same on the 3x3 grid and on the single-panel inclusive canvas.
    auto header_px = [&](int leg_ncol) {
        return 124 + ((n_leg + leg_ncol - 1) / leg_ncol) * 22 + 8;
    };
    // n_show = how many off-scale entries are spelled out before the running total takes over.
    // Each entry now names its pair-eta cell and its charge combination, so it is ~3x longer than
    // the bare "#DeltaR=..: .." it replaced and the 900 px inclusive canvas fits fewer of them
    // than the 1560 px grid; the count is therefore the caller's (it knows the canvas width).
    auto draw_header = [&](double H, const std::string& cell_text,
                           const std::vector<std::string>& offscale, int leg_ncol, size_t n_show) {
        auto ny = [&](double px) { return 1.0 - px / H; };
        TLatex t; t.SetNDC(); t.SetTextFont(42);
        t.SetTextSize(22.0 / H);
        t.DrawLatex(0.02, ny(26), Form("%s, %s,  %s%s%s", cfg.sample_text.c_str(),
                                       wp_text.c_str(), quantity_tex.c_str(),
                                       nocorr ? "" : " / plateau", cell_text.c_str()));
        t.SetTextSize(17.0 / H);
        t.DrawLatex(0.02, ny(56), def_line.c_str());
        if (!eq_line1.empty()) t.DrawLatex(0.02, ny(78), eq_line1.c_str());
        if (!eq_line2.empty()) t.DrawLatex(0.02, ny(96), eq_line2.c_str());
        if (!offscale.empty()) {
            std::string sub = "off scale: ";
            for (size_t i = 0; i < offscale.size() && i < n_show; ++i)
                sub += offscale[i]
                     + std::string(i + 1 < offscale.size() && i + 1 < n_show ? ";  " : "");
            if (offscale.size() > n_show) sub += Form(";  ... (%zu total)", offscale.size());
            t.SetTextColor(kGray + 3);
            t.SetTextSize(15.0 / H);
            t.DrawLatex(0.02, ny(114), sub.c_str());
            t.SetTextColor(kBlack);
        }
        auto* leg = new TLegend(0.02, ny(header_px(leg_ncol) - 8.), 0.98, ny(126));
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(16.0 / H);
        leg->SetNColumns(leg_ncol);
        fill_legend(leg);
        leg->Draw();
    };

    // ---- one canvas per pair-pT bin ----------------------------------------------------------
    // subplot grid: nrows >= ncols, nrows ~ sqrt(N)  (feedback_subplot_layout)
    const int ncol = (int)std::ceil(std::sqrt((double)neta));
    const int nrow = (int)std::ceil((double)neta / ncol);
    // The pair-pT edges in the FILE NAME must be the same numbers the canvas prints. "%.0f" wrote
    // pairpt_8_12 for the bin the canvas labels 8.0-11.5 GeV: two namings of one binning, which
    // is exactly the drift .claude/CLAUDE.md §Binnings exists to stop.
    auto pt_png_tag = [&](int iy) {
        return std::string(Form("pairpt_%.1f_%.1f", series[0].hplat->GetXaxis()->GetBinLowEdge(iy),
                                                    series[0].hplat->GetXaxis()->GetBinUpEdge(iy)));
    };
    auto pt_label = [&](int iy) {
        return std::string(Form(",  %.1f < p_{T}^{pair} < %.1f GeV",
                                series[0].hplat->GetXaxis()->GetBinLowEdge(iy),
                                series[0].hplat->GetXaxis()->GetBinUpEdge(iy)));
    };
    const int kHeaderPx = header_px(n_leg);
    for (int iy = 1; iy <= npt; ++iy) {
        const int    canv_h = 470 * nrow + kHeaderPx;
        const double hfrac  = double(kHeaderPx) / canv_h;
        TCanvas c(Form("c_%s_%s_%s_pt%d", tag.c_str(), method.c_str(), mode.c_str(), iy), "",
                  520 * ncol, canv_h);
        auto* grid = new TPad(Form("grid_%s_%s_%s_pt%d", tag.c_str(), method.c_str(),
                                   mode.c_str(), iy), "", 0., 0., 1., 1. - hfrac);
        grid->Draw();
        grid->Divide(ncol, nrow);
        std::vector<std::string> offscale;
        for (int iz = 1; iz <= neta; ++iz) {
            grid->cd(iz);
            draw_cell(iy, iz, g_ylo, g_yhi, offscale);
        }
        c.cd(0);
        draw_header(canv_h, pt_label(iy), offscale, n_leg, 3);
        const std::string png = odir + tag + "_dr_fit_" + method + "_" + pt_png_tag(iy) + ".png";
        c.SaveAs(png.c_str());
        std::cout << "  wrote " << png << "  (" << offscale.size() << " points off scale)\n";
    }

    // ---- the inclusive cell, on its own canvas ------------------------------------------------
    {
        const int    kInclHeaderPx = header_px(3);
        const int    incl_h = 700 + kInclHeaderPx;
        const double ihfrac = double(kInclHeaderPx) / incl_h;
        TCanvas c(Form("c_%s_%s_%s_incl", tag.c_str(), method.c_str(), mode.c_str()), "",
                  900, incl_h);
        auto* ipad = new TPad("ipad", "", 0., 0., 1., 1. - ihfrac);
        ipad->SetBottomMargin(0.13); ipad->Draw(); ipad->cd();
        std::vector<std::string> offscale;
        // The inclusive canvas uses the range built from the INCLUSIVE points (compute_range(true)
        // above), not the directory-wide one: it is a stand-alone figure, nobody flips between it
        // and the pair-pT files, and the shared range compressed its data into <10% of the frame.
        // Still central values only, never y +- error.
        draw_cell(0, 0, incl_rng.first, incl_rng.second, offscale);
        c.cd(0);
        draw_header(incl_h, "", offscale, 3, 2);   // 900 px wide: two entries fit on the line
        const std::string png = odir + tag + "_dr_fit_" + method + "_inclusive.png";
        c.SaveAs(png.c_str());
        std::cout << "  wrote " << png << "  (" << offscale.size() << " points off scale)\n";
    }

    // ---- SAME SIGN / OPPOSITE SIGN RATIO (sign_sepr only) -------------------------------------
    // The sign split exists to answer ONE question -- does the dR correction depend on the charge
    // combination? -- so the two series are expected to agree and the figure owes the reader their
    // ratio (R3). It gets its OWN canvas rather than a sub-pad under each of the nine panels: at
    // 520x470 px a panel already carries two parameter columns, and splitting it 70/30 leaves a
    // ratio pad ~140 px tall whose tick labels are unreadable at this canvas size.
    // The ratio is formed from the MEASURED, plateau-normalized points of the two signs. Same-sign
    // and opposite-sign pairs are DISJOINT samples, so the two errors are independent and add in
    // quadrature.
    if (sepr) {
        auto cell_ratio = [&](int iy, int iz) -> TGraphErrors* {
            auto* gr = new TGraphErrors();
            if (!cell_drawable(series[0], iy, iz) || !cell_drawable(series[1], iy, iz))
                return gr;
            auto* gs = cell_points(series[0], iy, iz, norm_of(series[0], iy, iz));
            auto* go = cell_points(series[1], iy, iz, norm_of(series[1], iy, iz));
            int k = 0;
            for (int i = 0; i < gs->GetN(); ++i) {
                double x1, y1;
                gs->GetPoint(i, x1, y1);
                for (int j2 = 0; j2 < go->GetN(); ++j2) {
                    double x2, y2;
                    go->GetPoint(j2, x2, y2);
                    if (std::fabs(x1 - x2) > 1e-6 || y2 == 0.) continue;
                    const double r  = y1 / y2;
                    const double e1 = gs->GetErrorY(i), e2 = go->GetErrorY(j2);
                    const double er = (y1 != 0.)
                        ? std::fabs(r) * std::sqrt((e1 / y1) * (e1 / y1) + (e2 / y2) * (e2 / y2))
                        : e1 / std::fabs(y2);
                    gr->SetPoint(k, x1, r);
                    gr->SetPointError(k, 0., er);
                    ++k;
                    break;
                }
            }
            delete gs; delete go;
            return gr;
        };

        // TWO ranges, for the same reason as the main panels above: the nine-panel pair-pT
        // files share ONE range across the directory so the pair-pT dependence can be read by
        // flipping between them, while the single-panel inclusive canvas is a stand-alone figure
        // and is scaled to its OWN points -- the shared range put the inclusive ratio points
        // (0.42-1.10) inside a 0-3.0 frame. Central values only in both cases.
        auto ratio_range = [&](bool incl_only) -> std::pair<double, double> {
            std::vector<std::pair<int, int>> cells;
            if (incl_only) cells.push_back({0, 0});
            else for (int iy = 1; iy <= npt; ++iy)
                     for (int iz = 1; iz <= neta; ++iz) cells.push_back({iy, iz});
            double lo = 1., hi = 1.;
            for (const auto& c : cells) {
                auto* g = cell_ratio(c.first, c.second);
                for (int i = 0; i < g->GetN(); ++i) {
                    double x, y;
                    g->GetPoint(i, x, y);
                    lo = std::min(lo, y);
                    hi = std::max(hi, y);
                }
                delete g;
            }
            const double span = hi - lo;
            lo = std::max(0.0,     lo - 0.10 * span);
            hi = std::min(kYcapHi, hi + 0.10 * span);
            if (hi - lo < 0.2) {   // widen AROUND the data, never to a fixed window -- see above
                const double mid = 0.5 * (lo + hi);
                lo = std::max(0.0,     mid - 0.10);
                hi = std::min(kYcapHi, mid + 0.10);
            }
            return {lo, hi};
        };
        const std::pair<double, double> rgrid = ratio_range(false);
        const std::pair<double, double> rincl = ratio_range(true);
        double r_lo = rgrid.first, r_hi = rgrid.second;   // draw_ratio_cell captures these

        auto draw_ratio_cell = [&](int iy, int iz, std::vector<std::string>& offscale) {
            gPad->SetLeftMargin(0.14);
            gPad->SetBottomMargin(0.14);
            gPad->SetTopMargin(0.10);
            TH1* fr = gPad->DrawFrame(0.0, r_lo, kXhi, r_hi);
            fr->GetXaxis()->SetTitle("#DeltaR(#mu_{1}, #mu_{2})");
            fr->GetYaxis()->SetTitle("same sign / opposite sign");
            fr->GetXaxis()->SetTitleSize(0.050);
            fr->GetYaxis()->SetTitleSize(0.050);
            fr->GetXaxis()->SetLabelSize(0.042);
            fr->GetYaxis()->SetLabelSize(0.042);
            fr->GetYaxis()->SetTitleOffset(1.25);
            DrawUnity(0.0, kXhi);
            // Same reason as in draw_cell(): the off-scale record is canvas-level, so the point
            // has to carry its own pair-eta cell. The two charge combinations need no tag here --
            // this figure IS their ratio and every y axis on it says so.
            const std::string cell_id = (iy == 0 && iz == 0)
                ? std::string()
                : std::string(Form("#eta^{pair} #in [%.1f,%.1f), ",
                                   series[0].hplat->GetYaxis()->GetBinLowEdge(iz),
                                   series[0].hplat->GetYaxis()->GetBinUpEdge(iz)));
            TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.045);
            if (iy == 0 && iz == 0)
                t.DrawLatex(0.18, 0.925, "inclusive (all p_{T}^{pair}, all #eta^{pair})");
            else
                t.DrawLatex(0.18, 0.925, Form("%.1f < #eta^{pair} < %.1f",
                                              series[0].hplat->GetYaxis()->GetBinLowEdge(iz),
                                              series[0].hplat->GetYaxis()->GetBinUpEdge(iz)));
            auto* g = cell_ratio(iy, iz);
            if (g->GetN() == 0) {
                t.SetTextSize(0.05);
                t.DrawLatex(0.20, 0.55, "no measurement");
                return;
            }
            g->SetMarkerStyle(20); g->SetMarkerSize(0.7);
            g->SetMarkerColor(kBlack); g->SetLineColor(kBlack); g->SetLineWidth(1);
            g->Draw("PZ same");
            for (int i = 0; i < g->GetN(); ++i) {
                double x, y;
                g->GetPoint(i, x, y);
                if (y <= r_hi && y >= r_lo) continue;
                const bool up = (y > r_hi);
                auto* ar = new TArrow(x, r_lo + (up ? 0.88 : 0.12) * (r_hi - r_lo),
                                      x, r_lo + (up ? 0.98 : 0.02) * (r_hi - r_lo), 0.008, "|>");
                ar->SetLineColor(kBlack); ar->SetFillColor(kBlack); ar->SetLineWidth(2);
                ar->Draw();
                offscale.push_back(Form("%s#DeltaR=%.2f: %.2f #pm %.2f", cell_id.c_str(),
                                        x, y, g->GetErrorY(i)));
            }
        };

        auto ratio_header = [&](double H, const std::string& cell_text,
                                const std::vector<std::string>& offscale, size_t n_show) {
            auto ny = [&](double px) { return 1.0 - px / H; };
            TLatex t; t.SetNDC(); t.SetTextFont(42);
            t.SetTextSize(22.0 / H);
            // The two charge combinations are named by the y-axis title of every panel, so the
            // headline does NOT repeat them: with them it overran the right edge of the 900 px
            // inclusive canvas and was silently cut off mid-word.
            t.DrawLatex(0.02, ny(26), Form("%s, %s,  %s%s%s",
                                           cfg.sample_text.c_str(), wp_text.c_str(),
                                           quantity_tex.c_str(), nocorr ? "" : " / plateau",
                                           cell_text.c_str()));
            t.SetTextSize(17.0 / H);
            t.DrawLatex(0.02, ny(56), def_line.c_str());
            if (!offscale.empty()) {
                std::string sub = "off scale: ";
                for (size_t i = 0; i < offscale.size() && i < n_show; ++i)
                    sub += offscale[i]
                         + std::string(i + 1 < offscale.size() && i + 1 < n_show ? ";  " : "");
                if (offscale.size() > n_show) sub += Form(";  ... (%zu total)", offscale.size());
                t.SetTextColor(kGray + 3);
                t.SetTextSize(15.0 / H);
                // 30 px below the definition line, not 22: both lines carry super/subscripts
                // (eps^{2mu4}_{dR} above, eta^{pair} below) and at 22 px the eta superscript was
                // drawn into the epsilon subscript.
                t.DrawLatex(0.02, ny(86), sub.c_str());
                t.SetTextColor(kBlack);
            }
        };

        const int kRatioHeaderPx = 104;   // title 26 + definition 56 + off-scale record 86 + margin
        for (int iy = 1; iy <= npt; ++iy) {
            const int    canv_h = 470 * nrow + kRatioHeaderPx;
            const double hfrac  = double(kRatioHeaderPx) / canv_h;
            TCanvas c(Form("cr_%s_%s_pt%d", tag.c_str(), method.c_str(), iy), "",
                      520 * ncol, canv_h);
            auto* grid = new TPad(Form("rgrid_%s_%s_pt%d", tag.c_str(), method.c_str(), iy), "",
                                  0., 0., 1., 1. - hfrac);
            grid->Draw();
            grid->Divide(ncol, nrow);
            std::vector<std::string> offscale;
            for (int iz = 1; iz <= neta; ++iz) {
                grid->cd(iz);
                draw_ratio_cell(iy, iz, offscale);
            }
            c.cd(0);
            ratio_header(canv_h, pt_label(iy), offscale, 3);
            const std::string png = rodir + tag + "_dr_fit_" + method + "_" + pt_png_tag(iy)
                                  + "_ratio.png";
            c.SaveAs(png.c_str());
            std::cout << "  wrote " << png << "  (" << offscale.size() << " points off scale)\n";
        }
        {
            const int    incl_h = 700 + kRatioHeaderPx;
            const double ihfrac = double(kRatioHeaderPx) / incl_h;
            TCanvas c(Form("cr_%s_%s_incl", tag.c_str(), method.c_str()), "", 900, incl_h);
            auto* ipad = new TPad("ripad", "", 0., 0., 1., 1. - ihfrac);
            ipad->SetBottomMargin(0.13); ipad->Draw(); ipad->cd();
            std::vector<std::string> offscale;
            r_lo = rincl.first; r_hi = rincl.second;   // this canvas is scaled to its own points
            draw_ratio_cell(0, 0, offscale);
            c.cd(0);
            ratio_header(incl_h, "", offscale, 2);
            const std::string png = rodir + tag + "_dr_fit_" + method + "_inclusive_ratio.png";
            c.SaveAs(png.c_str());
            std::cout << "  wrote " << png << "  (" << offscale.size() << " points off scale)\n";
        }
    }

    // ---- READ-BACK VERIFICATION + FLATNESS AUDIT ----------------------------------------------
    // TWO INDEPENDENT criteria, deliberately reported separately -- conflating them would have
    // made `expo` look like a persistence bug when it is really a shape-constraint failure:
    //
    //   PERSISTENCE : the object read back out of the file evaluates to a FINITE number
    //                 everywhere, AND -- for the COMPACT-SUPPORT shapes only -- reproduces the
    //                 exact 1 they are built to return beyond Rp. Those shapes (powerlaw*,
    //                 polyu*, interp: everything written in u = max(0, 1 - dR/Rp)) are
    //                 analytically 1 for every dR >= Rp, so ANY deviation on read-back -- above
    //                 all the 0 that a compiled/lambda TF1 returns outside its stored range
    //                 (pp_trig_eff_highpt_jump.md) -- is a persistence failure and nothing else.
    //                 `expo` has no compact support, so for it only finiteness is a persistence
    //                 statement; how close it gets to 1 is a FIT property, judged below.
    //   FLATNESS    : |f(dR) - 1| <= kFlatTol for every dR >= Rp -- the user's requirement that
    //                 Step 4 be flat beyond ~0.3 and Step 3 beyond ~0.5. Compact support gives
    //                 it for free; `expo` has to earn it and, in the overlay's unmeasurable
    //                 cells, does not (its A runs to the parameter limit and the "correction"
    //                 never returns to 1 even at dR = 20). That is a reason to reject `expo`,
    //                 not a bug in the persistence layer.
    //
    // One report per SERIES (the fits are per series), named after the physics:
    // readback_check.txt / readback_check_same_sign.txt / readback_check_opposite_sign.txt.
    for (const auto& s : series) {
        const double Rp = (step == 3) ? 0.5 : 0.3;
        constexpr double kFlatTol = 1e-3;   // shape: allowed deviation from 1 beyond Rp
        // Compact support is a property of the stored formula, read back from the file itself
        // rather than assumed from the method name.
        bool compact = (method == "interp");
        {
            const FittedFunc Fp = LoadFunc(s.ffit, step, 0, 0);
            if (Fp.f) compact = TString(Fp.f->GetExpFormula()).Contains("TMath::Max");
        }
        const std::string rpath = rdir + "readback_check" + DrCorrSignFileTag(s.sign) + ".txt";
        std::ofstream os(rpath);
        os << "# Read-back verification of "
           << DrCorrFitFile(cfg, use_tight_wp, step, method, s.sign)
           << "\n# performed in a SEPARATE ROOT process from the fit (that is the point).\n"
           << "# series: " << (s.sign.empty() ? "sign-integrated" : s.legend) << "\n"
           << "# shape has compact support (exactly 1 for dR >= Rp by construction): "
           << (compact ? "YES" : "NO") << "\n"
           << "# PERSISTENCE: every value finite"
           << (compact ? (nocorr
                          ? ", and exactly the fitted baseline C for dR >= Rp (a read-back"
                            " returning 0 outside the stored range -- the compiled-TF1 bug --"
                            " fails here)"
                          : ", and exactly 1 for dR >= Rp (a read-back returning 0 outside the"
                            " stored range -- the compiled-TF1 bug -- fails here)") : "") << "\n"
           << (nocorr
               ? "# FLATNESS   : |f(dR) - C| <= "     // no external normalization: the asymptote
               : "# FLATNESS   : |f(dR) - 1| <= ")    // is the fitted baseline C, not 1
           << kFlatTol << " for dR >= Rp = " << Rp
           << (nocorr ? "  (C = the fitted/pinned baseline of that cell)" : "") << "\n\n";
        int nchecked = 0, nfail_persist = 0, nfail_flat = 0;
        double worst_flat = 0.;
        for (int iy = 0; iy <= npt; ++iy) {
            for (int iz = 0; iz <= neta; ++iz) {
                if (!(iy == 0 && iz == 0) && (iy == 0 || iz == 0)) continue;
                const FittedFunc F = LoadFunc(s.ffit, step, iy, iz);
                if (!F.valid()) continue;
                ++nchecked;
                bool bad_persist = false, bad_flat = false;
                std::string line = (iy == 0 && iz == 0) ? "inclusive"
                                                        : Form("pt%d_eta%d", iy, iz);
                // PERSISTENCE uses the CELL's OWN support radius: `powerlaw_floatRp` fits R_p
                // per cell and it may land above the nominal onset, so the function is
                // legitimately != 1 just above the nominal Rp. FLATNESS keeps using the NOMINAL
                // Rp, because "flat beyond 0.3 / 0.5" is the requirement, and a cell whose
                // fitted R_p drifted past it genuinely fails that requirement.
                double rp_cell = Rp;
                if (F.f) {
                    const int ip = F.f->GetParNumber("R_{p}");
                    if (ip >= 0) rp_cell = std::max(rp_cell, F.f->GetParameter(ip));
                }
                // THE ASYMPTOTE this cell is supposed to reach. Exactly 1 when the curve was
                // normalized by the plateau; the cell's own fitted baseline C when it was not
                // (for the interpolation, the constant its flat branch was pinned to). Testing
                // a no-correction curve against 1 would call every cell non-flat and would say
                // nothing about either persistence or flatness.
                double flat_ref = 1.0;
                if (nocorr) {
                    if (F.f) {
                        const int ic = F.f->GetParNumber("C");
                        flat_ref = (ic >= 0) ? F.f->GetParameter(ic) : F.Eval(20.0);
                    } else {
                        flat_ref = F.Eval(kFlatProbeX);   // pinned flat branch of the interpolation
                    }
                }
                for (double x : {0.0, 0.5 * Rp, Rp, 0.9, 1.5, 3.0, 8.0, 20.0}) {
                    const double v = F.Eval(x);
                    line += Form("  f(%.2f)=%.5f", x, v);
                    if (!std::isfinite(v)) bad_persist = true;
                    if (compact && x >= rp_cell && std::fabs(v - flat_ref) > 1e-9)
                        bad_persist = true;
                    if (x >= Rp) {
                        worst_flat = std::max(worst_flat, std::fabs(v - flat_ref));
                        if (std::fabs(v - flat_ref) > kFlatTol) bad_flat = true;
                    }
                }
                if (bad_persist) { ++nfail_persist; line += "   <-- PERSISTENCE FAIL"; }
                if (bad_flat)    { ++nfail_flat;    line += "   <-- not flat beyond Rp"; }
                os << line << "\n";
            }
        }
        // The deviation is measured from the asymptote the mode actually has: 1 when the curves
        // were normalized by the plateau, the cell's own baseline C when they were not.
        const char* dev_sym = nocorr ? "|f-C|" : "|f-1|";
        os << "\nchecked " << nchecked << " functions, " << nfail_persist << " FAILED persistence"
           << "\nflatness beyond Rp: " << nfail_flat << " function(s) exceed " << kFlatTol
           << ", worst " << dev_sym << " = " << Form("%.2e", worst_flat) << "\n";
        std::cout << "  read-back check ["
                  << (s.sign.empty() ? "sign-integrated" : s.legend) << "]: " << nchecked
                  << " functions, " << nfail_persist << " FAILED persistence; " << nfail_flat
                  << " not flat beyond Rp (worst " << dev_sym << " = "
                  << Form("%.2e", worst_flat) << ")\n"
                  << "  wrote " << rpath << "\n";
        if (nfail_persist) std::cout << "  ** PERSISTENCE FAILURE -- do NOT use these fits **\n";
    }

    for (auto& s : series) s.ffit->Close();
    fh->Close();
}

// Convenience: every method for one (sample, WP, step, mode).
void plot_dr_correction_fits_all(const std::string& sample = "pp_full", bool use_tight_wp = true,
                                 int step = 3, const std::string& mode = "sign_intgr",
                                 const std::string& plateau_mode = "corr")
{
    for (const std::string& m : {"powerlaw_fixedRp", "powerlaw_floatRp", "expo",
                                 "polyu_fixedRp", "interp"})
        plot_dr_correction_fits(sample, use_tight_wp, step, m, mode, plateau_mode);
}
