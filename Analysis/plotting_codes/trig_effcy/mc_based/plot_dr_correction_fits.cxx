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
// a spurious pp cross-section jump (pp_trig_eff_highpt_jump.md). Everything drawn in red here
// comes from TF1::Eval / TGraph::Eval of the object read out of the file, evaluated at the EXACT
// DeltaR (never resampled to the nearest bin), and CheckReadBack() below additionally probes each
// function well outside its fit range and fails loudly if it collapses to 0 or drifts off 1.
//
// LAYOUT (user requirement): 1 PNG per pair-pT bin, 1 subplot per pair-eta bin (9 bins -> 3x3,
// nrows >= ncols ~ sqrt(N)); measured values BLACK with error bars, fitted function RED; plus one
// PNG for the inclusive cell.
//
// The black points are the SAME plateau-normalized values the fit saw: fine zoom bins for
// dR < 1 (the fit domain) and, beyond dR = 1, the coarse full-range bins -- those are the ones
// that DEFINED the plateau, so they must sit at 1 by construction and are shown as a check, not
// as fit input.
//
// Compile/run (ACLiC, from this directory):
//   root -l -b -q -e '.L plot_dr_correction_fits.cxx+'                       // compile only
//   root -l -b -q 'plot_dr_correction_fits.cxx+("pp_full", true, 3, "polyu_fixedRp")'

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

#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double kXhi     = 2.0;    // drawn dR range (fit domain 0-1 + plateau check to 2)
constexpr double kFitHi   = 1.0;    // fit domain upper edge (must match fit_dr_corrections.cxx)
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
inline std::pair<std::string, std::string> MethodFormulaTex(const std::string& method)
{
    if (method == "powerlaw_fixedRp" || method == "powerlaw_floatRp")
        return {"f(#DeltaR) = 1 + A #upoint u^{n}",
                "u #equiv max(0, 1 - #DeltaR/R_{p})"};
    if (method == "expo")
        return {"f(#DeltaR) = 1 + A #upoint exp[-(#DeltaR/#lambda)^{p}]", ""};
    if (method == "polyu_fixedRp")
        return {"f(#DeltaR) = 1 + u^{2}(a_{2} + a_{3}u + a_{4}u^{2})",
                "u #equiv max(0, 1 - #DeltaR/R_{p})"};
    if (method == "interp")
        // P2: a real piecewise DEFINITION with R_p defined and its value drawn, matching the
        // treatment polyu_fixedRp already gets. "linear interpolation of the points" is prose,
        // not an equation, and R_p was never defined anywhere on the canvas.
        return {"f(#DeltaR) = piecewise-linear through the measured points for #DeltaR < R_{p}",
                "f(#DeltaR) = 1 for #DeltaR #geq R_{p}"};
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

void DrawUnity(double xlo, double xhi)
{
    auto* l = new TLine(xlo, 1.0, xhi, 1.0);
    l->SetLineStyle(2);
    l->SetLineColor(kGray + 2);
    l->Draw("same");
}

}  // namespace

// step   : 3 (cross / 2mu4 term) or 4 (single leg)
// method : powerlaw_fixedRp | powerlaw_floatRp | expo | polyu_fixedRp | interp
void plot_dr_correction_fits(const std::string& sample = "pp_full", bool use_tight_wp = true,
                             int step = 3, const std::string& method = "polyu_fixedRp")
{
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    const DrCorrSample cfg = GetDrCorrSample(sample, use_tight_wp);
    const std::string wp_suf  = DrCorrWpSuffix(use_tight_wp);
    const std::string wp_text = use_tight_wp ? "Tight muons" : "Medium muons";
    const std::string wp_dir  = DrCorrWpDir(use_tight_wp);
    const std::string h_pref  = (step == 3) ? "h_mc_dr_" : "h_mc_single_dr_";
    const bool has_cov = (step == 4);
    const std::string quantity_tex = (step == 3)
        ? cfg.eps_dr_text                                   // eps_dR^2mu4 (pp) / ^cross (PbPb)
        : std::string("#varepsilon_{#DeltaR}^{single}");
    const std::string tag = "step" + std::to_string(step);

    // out_base now carries the variant tag (mc_based / _medium / _pt4bin / _pt4bin_medium),
    // so this needs no per-directory token. Before 2026-08-05 it had neither, which put the
    // 4-bin fit plots into the 8-bin directories under identical names.
    const std::string odir = cfg.out_base + tag + "_dr_fit/" + method + "/" + wp_dir;
    gSystem->mkdir(odir.c_str(), kTRUE);

    // ---- inputs: measured histograms, plateaus (inside the fit file) and the FITS ------------
    const std::string hist_path = cfg.mc_dir + "mc_trig_eff_hists_" + cfg.mc_label + wp_suf
                                + "_" + tag + ".root";
    TFile* fh  = OpenRead(hist_path);
    TFile* ffit = OpenRead(DrCorrFitFile(cfg, use_tight_wp, step, method));

    TH3D* zn = GetObj<TH3D>(fh, h_pref + "zoom_vs_pt_eta_num");
    TH3D* zd = GetObj<TH3D>(fh, h_pref + "zoom_vs_pt_eta_denom");
    TH3D* za = GetObj<TH3D>(fh, h_pref + "zoom_vs_pt_eta_errA");
    TH3D* zb = GetObj<TH3D>(fh, h_pref + "zoom_vs_pt_eta_errB");
    TH3D* zp = has_cov ? GetObj<TH3D>(fh, h_pref + "zoom_vs_pt_eta_covP") : nullptr;
    TH3D* zq = has_cov ? GetObj<TH3D>(fh, h_pref + "zoom_vs_pt_eta_covQ") : nullptr;
    TH3D* fn = GetObj<TH3D>(fh, h_pref + "full_vs_pt_eta_num");
    TH3D* fd = GetObj<TH3D>(fh, h_pref + "full_vs_pt_eta_denom");
    TH3D* fa = GetObj<TH3D>(fh, h_pref + "full_vs_pt_eta_errA");
    TH3D* fb = GetObj<TH3D>(fh, h_pref + "full_vs_pt_eta_errB");
    TH3D* fp = has_cov ? GetObj<TH3D>(fh, h_pref + "full_vs_pt_eta_covP") : nullptr;
    TH3D* fq = has_cov ? GetObj<TH3D>(fh, h_pref + "full_vs_pt_eta_covQ") : nullptr;

    TH2D* hplat = GetObj<TH2D>(ffit, "h_" + tag + "_plateau");
    TH2D* hchi  = GetObj<TH2D>(ffit, "h_" + tag + "_chi2ndf");
    TH1D* hpinc = GetObj<TH1D>(ffit, "h_" + tag + "_plateau_inclusive");

    const int npt  = hplat->GetNbinsX();
    const int neta = hplat->GetNbinsY();

    // ---- the measured, plateau-normalized points of one cell (fit input + plateau check) -----
    auto cell_points = [&](int iy, int iz, double plateau) -> TGraphErrors* {
        auto* g = new TGraphErrors();
        int k = 0;
        TH1D* rz = DrCellRatio(zn, zd, za, zb, zp, zq, iy, iz, Form("pz_%d_%d", iy, iz));
        for (int i = 1; i <= rz->GetNbinsX(); ++i) {
            const double e = rz->GetBinError(i);
            if (e <= 0.) continue;
            g->SetPoint(k, rz->GetBinCenter(i), rz->GetBinContent(i) / plateau);
            g->SetPointError(k, 0., e / plateau);
            ++k;
        }
        delete rz;
        TH1D* rf = DrCellRatio(fn, fd, fa, fb, fp, fq, iy, iz, Form("pf_%d_%d", iy, iz));
        for (int i = 1; i <= rf->GetNbinsX(); ++i) {
            const double x = rf->GetBinCenter(i), e = rf->GetBinError(i);
            if (x <= kFitHi || x > kXhi || e <= 0.) continue;   // only the dR > 1 check region
            g->SetPoint(k, x, rf->GetBinContent(i) / plateau);
            g->SetPointError(k, 0., e / plateau);
            ++k;
        }
        delete rf;
        g->SetMarkerStyle(20);
        g->SetMarkerSize(0.7);
        g->SetMarkerColor(kBlack);
        g->SetLineColor(kBlack);
        g->SetLineWidth(1);
        return g;
    };

    // ---- one subplot ------------------------------------------------------------------------
    // Returns the list of points pushed off the (capped) frame, so the caller can list them on
    // the canvas: a hidden point must never look like a missing one.
    auto draw_cell = [&](int iy, int iz, double ylo, double yhi,
                         std::vector<std::string>& offscale) {
        const double plateau = (iy == 0 && iz == 0) ? hpinc->GetBinContent(1)
                                                    : hplat->GetBinContent(iy, iz);
        gPad->SetLeftMargin(0.14);
        gPad->SetBottomMargin(0.14);
        // Reserved strip above the frame for the pair-eta label: when the y range is clamped at
        // the cap, the off-scale arrows sit just under the frame top and ran through the label.
        gPad->SetTopMargin(0.10);
        TH1* fr = gPad->DrawFrame(0.0, ylo, kXhi, yhi);
        fr->GetXaxis()->SetTitle("#DeltaR(#mu_{1}, #mu_{2})");
        fr->GetYaxis()->SetTitle((quantity_tex + " / plateau").c_str());
        fr->GetXaxis()->SetTitleSize(0.050);
        fr->GetYaxis()->SetTitleSize(0.050);
        fr->GetXaxis()->SetLabelSize(0.042);
        fr->GetYaxis()->SetLabelSize(0.042);
        fr->GetYaxis()->SetTitleOffset(1.25);
        DrawUnity(0.0, kXhi);

        // Same screen as the fit stage (shared helper): an unmeasurable cell is NOT drawn --
        // drawing one means dividing the curve by a near-zero plateau and presenting the
        // resulting O(10) excursion as a correction.
        const double plat_err = (iy == 0 && iz == 0) ? hpinc->GetBinError(1)
                                                     : hplat->GetBinError(iy, iz);
        if (!DrCorrPlateauUsable(plateau, plat_err)) {
            TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.05);
            t.DrawLatex(0.20, 0.55, "no fit");
            return;
        }

        auto* g = cell_points(iy, iz, plateau);
        const FittedFunc Fpre = LoadFunc(ffit, step, iy, iz);

        // RED = the function AS READ BACK FROM THE FILE, evaluated at the exact dR.
        // Drawn BEFORE the points: a 2px line over 0.7-size markers hides the measurement, and
        // the measurement is what the reader has to be able to judge the fit against.
        const FittedFunc F = Fpre;
        if (F.valid()) {
            auto* gf = new TGraph();
            const int NP = 500;
            for (int i = 0; i <= NP; ++i) {
                const double x = kXhi * i / NP;
                gf->SetPoint(i, x, F.Eval(x));   // continuous function at the EXACT dR
            }
            gf->SetLineColor(kRed + 1);
            gf->SetLineWidth(2);
            gf->Draw("L same");
        }

        // POINTS OUTSIDE THE FIT DOMAIN ARE DRAWN IN BLUE (user, 2026-08-06), so a human can see
        // at a glance which points the fit was actually constrained by. Only dR < kFitHi enters
        // the fit; the points beyond it are the plateau-window bins, shown as a check that the
        // normalized curve really does sit at 1 there. Black = fitted, blue = excluded.
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
        g_fit->SetMarkerStyle(20); g_fit->SetMarkerSize(0.7);
        g_fit->SetMarkerColor(kBlack);   g_fit->SetLineColor(kBlack);
        g_exc->SetMarkerStyle(24); g_exc->SetMarkerSize(0.7);
        g_exc->SetMarkerColor(kBlue + 1); g_exc->SetLineColor(kBlue + 1);
        if (g_fit->GetN()) g_fit->Draw("PZ same");
        if (g_exc->GetN()) g_exc->Draw("PZ same");
        for (int i = 0; i < g->GetN(); ++i) {
            double x, y;
            g->GetPoint(i, x, y);
            if (y <= yhi && y >= ylo) continue;
            const bool up = (y > yhi);
            auto* ar = new TArrow(x, ylo + (up ? 0.88 : 0.12) * (yhi - ylo),
                                  x, ylo + (up ? 0.98 : 0.02) * (yhi - ylo), 0.008, "|>");
            const Color_t acol = (x < kFitHi) ? kBlack : kBlue + 1;
            ar->SetLineColor(acol); ar->SetFillColor(acol); ar->SetLineWidth(2);
            ar->Draw();
            offscale.push_back(Form("#DeltaR=%.2f: %.2f #pm %.2f", x, y, g->GetErrorY(i)));
        }

        TLatex t;
        t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.045);
        if (iy == 0 && iz == 0) {
            t.DrawLatex(0.18, 0.925, "inclusive (all p_{T}^{pair}, all #eta^{pair})");
        } else {
            t.DrawLatex(0.18, 0.925, Form("%.1f < #eta^{pair} < %.1f",
                                         hplat->GetYaxis()->GetBinLowEdge(iz),
                                         hplat->GetYaxis()->GetBinUpEdge(iz)));
        }
        // Annotation block in the lower-RIGHT quadrant: the correction lies at 1 for dR > Rp, so
        // everything right of the turn-on and below unity is empty. (It used to start at NDC
        // y=0.32 on the left and ran off the bottom of the frame, through the axis labels.)
        t.SetTextSize(0.033);
        double ty = 0.50;
        // The FORMULA first: parameter names alone ("A", "n", "R_{p}") are meaningless to a
        // reader who has not read the fitter source.
        {
            const auto ftex = MethodFormulaTex(method);
            if (!ftex.first.empty()) {
                t.SetTextSize(0.029);
                t.DrawLatex(0.52, ty, ftex.first.c_str());
                ty -= 0.045;
                if (!ftex.second.empty()) {
                    t.DrawLatex(0.52, ty, ftex.second.c_str());
                    ty -= 0.045;
                }
                t.SetTextSize(0.033);
            }
        }
        t.DrawLatex(0.52, ty, Form("plateau = %.4f", plateau));
        ty -= 0.052;
        if (F.f) {
            for (int ip = 0; ip < F.f->GetNpar(); ++ip) {
                // A parameter that was FIXED in the fit has zero error; printing "0.5 +- 0"
                // reads as an infinitely precise measurement. Mark it as fixed instead.
                double plo = 0., phi = 0.;
                F.f->GetParLimits(ip, plo, phi);
                const bool fixed = (F.f->GetParError(ip) == 0.) || (plo == phi && plo != 0.);
                t.DrawLatex(0.52, ty, fixed
                    ? Form("%s = %.3g (fixed)", F.f->GetParName(ip), F.f->GetParameter(ip))
                    : Form("%s = %.3g #pm %.2g", F.f->GetParName(ip),
                           F.f->GetParameter(ip), F.f->GetParError(ip)));
                ty -= 0.052;
            }
        }
        // chi2/ndf ALWAYS as a number. The inclusive cell has no entry in the per-cell TH2D, so
        // take it from the TF1 itself (TF1 persists its chi2 and ndf), never a pointer to a file:
        // the audience sees only the plot.
        if (F.f) {
            const double c = (iy == 0 && iz == 0)
                ? (F.f->GetNDF() > 0 ? F.f->GetChisquare() / F.f->GetNDF() : -1.)
                : hchi->GetBinContent(iy, iz);
            if (c >= 0.) t.DrawLatex(0.52, ty, Form("#chi^{2}/ndf = %.2f", c));
        }
        // `g` is intentionally NOT deleted: ~TGraph removes the object from every pad that
        // draws it, so deleting it here silently erased all the black points from the canvas.
    };

    // ---- y range: ONE range for the WHOLE METHOD DIRECTORY (user, 2026-08-05) ---------------
    // Common not just across the eta panels of one canvas, but across EVERY pair-pT bin -- i.e.
    // every PNG in this <step>_dr_fit/<method>/ directory shares one y axis. Otherwise each file
    // silently rescales and the pair-pT dependence, which is the whole point of binning in
    // pair pT, cannot be read by flipping between the files.
    //
    // The range is set from the CENTRAL VALUES only, NOT value +- error (user): in the sparse
    // high-pair-pT cells the error bars are several times the correction itself, and including
    // them let one noisy cell dictate the axis for every plot, squashing all the real structure
    // into a sliver. Points whose central value still falls outside the capped range are marked
    // with an arrow by draw_cell(), so nothing is dropped silently.
    double g_ylo = 1., g_yhi = 1.;
    for (int iy = 1; iy <= npt; ++iy) {
        for (int iz = 1; iz <= neta; ++iz) {
            const double plateau = hplat->GetBinContent(iy, iz);
            // Screen with the SAME test the drawing code uses: an unmeasurable cell is never
            // drawn, so letting its eps/plateau values (O(10-130) when the plateau is ~0.008)
            // set the shared range pushed every PNG in the directory to the 3.0 cap and squashed
            // the structure the panels exist to show.
            if (!DrCorrPlateauUsable(plateau, hplat->GetBinError(iy, iz))) continue;
            auto* g = cell_points(iy, iz, plateau);
            for (int i = 0; i < g->GetN(); ++i) {
                double x, y;
                g->GetPoint(i, x, y);
                g_ylo = std::min(g_ylo, y);
                g_yhi = std::max(g_yhi, y);
            }
            delete g;
        }
    }
    {
        // 25% headroom at the top for the legend; a legend over data points hides them just as
        // badly as an axis that crops them.
        const double span = g_yhi - g_ylo;
        g_ylo = std::max(kYcapLo, g_ylo - 0.05 * span);
        g_yhi = std::min(kYcapHi, g_yhi + 0.25 * span);
        if (g_yhi - g_ylo < 0.2) { g_ylo = std::max(kYcapLo, 0.9); g_yhi = 1.12; }
    }
    std::cout << "  common y range for every PNG in this method dir: ["
              << g_ylo << ", " << g_yhi << "]\n";
    auto y_range = [&](int /*iy*/, double& ylo, double& yhi) { ylo = g_ylo; yhi = g_yhi; };

    // ---- one canvas per pair-pT bin ----------------------------------------------------------
    // subplot grid: nrows >= ncols, nrows ~ sqrt(N)  (feedback_subplot_layout)
    const int ncol = (int)std::ceil(std::sqrt((double)neta));
    const int nrow = (int)std::ceil((double)neta / ncol);
    for (int iy = 1; iy <= npt; ++iy) {
        double ylo, yhi;
        y_range(iy, ylo, yhi);
        // A dedicated HEADER STRIP at the top of the canvas: the divided pads otherwise reach
        // the very top edge, and the title + fit-range line were drawn on top of the first
        // panel's pair-eta label. The panel grid lives in its own pad below the strip.
        const int   kHeaderPx = 70;
        const int   canv_h    = 470 * nrow + kHeaderPx;
        const double hfrac    = double(kHeaderPx) / canv_h;
        TCanvas c(Form("c_%s_%s_pt%d", tag.c_str(), method.c_str(), iy), "",
                  520 * ncol, canv_h);
        auto* grid = new TPad(Form("grid_%s_%s_pt%d", tag.c_str(), method.c_str(), iy), "",
                              0., 0., 1., 1. - hfrac);
        grid->Draw();
        grid->Divide(ncol, nrow);
        std::vector<std::string> offscale;
        for (int iz = 1; iz <= neta; ++iz) {
            grid->cd(iz);
            draw_cell(iy, iz, ylo, yhi, offscale);
            if (iz == 1) {   // legend once
                auto* gd = new TGraphErrors(); gd->SetMarkerStyle(20); gd->SetMarkerColor(kBlack);
                gd->SetLineColor(kBlack);
                auto* gf = new TGraph(); gf->SetLineColor(kRed + 1); gf->SetLineWidth(2);
                auto* leg = new TLegend(0.45, 0.78, 0.96, 0.90);
                leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.036);
                leg->AddEntry(gd, "measurement", "lp");
                leg->AddEntry(gf, "fit", "l");
                leg->Draw();
            }
        }
        c.cd(0);
        TLatex st; st.SetNDC(); st.SetTextFont(42); st.SetTextSize(0.019);
        // No "Step %d": the step number is an internal pipeline index. The quantity itself
        // (eps_dR^2mu4 / eps_dR^single) already identifies what is drawn.
        st.DrawLatex(0.02, 1. - 0.40 * hfrac,
            Form("%s, %s,  %s / plateau,  %.1f < p_{T}^{pair} < %.1f GeV",
                 cfg.sample_text.c_str(), wp_text.c_str(), quantity_tex.c_str(),
                 hplat->GetXaxis()->GetBinLowEdge(iy),
                 hplat->GetXaxis()->GetBinUpEdge(iy)));
        TLatex n; n.SetNDC(); n.SetTextFont(42); n.SetTextSize(0.0135);
        n.SetTextColor(kGray + 3);
        // Formal statement of the fit range only -- no explanatory prose on a physics figure.
        std::string sub = Form("fit range: #DeltaR < %.1f", kFitHi);
        if (!offscale.empty()) {
            sub += "   |   off scale: ";
            for (size_t i = 0; i < offscale.size() && i < 6; ++i)
                sub += offscale[i] + std::string(i + 1 < offscale.size() && i < 5 ? ", " : "");
            if (offscale.size() > 6) sub += Form(", ... (%zu total)", offscale.size());
        }
        n.DrawLatex(0.02, 1. - 0.80 * hfrac, sub.c_str());
        const std::string png = odir + Form("%s_dr_fit_%s_pairpt_%.0f_%.0f.png", tag.c_str(),
                                            method.c_str(),
                                            hplat->GetXaxis()->GetBinLowEdge(iy),
                                            hplat->GetXaxis()->GetBinUpEdge(iy));
        c.SaveAs(png.c_str());
        std::cout << "  wrote " << png << "\n";
    }

    // ---- the inclusive cell, on its own canvas ------------------------------------------------
    {
        // Reserved header strip, as on the grid canvases: without it the headline was drawn at
        // y = 0.96 straight on top of the panel label "inclusive (all p_T^pair, all eta^pair)".
        const int kInclHeaderPx = 46;
        const int incl_h = 700 + kInclHeaderPx;
        const double ihfrac = double(kInclHeaderPx) / incl_h;
        TCanvas c(Form("c_%s_%s_incl", tag.c_str(), method.c_str()), "", 900, incl_h);
        auto* ipad = new TPad("ipad", "", 0., 0., 1., 1. - ihfrac);
        ipad->SetBottomMargin(0.13); ipad->Draw(); ipad->cd();
        std::vector<std::string> offscale;
        // SAME y range as every other PNG in this method directory (user): the inclusive canvas
        // used to compute its own from y +- error, which broke both the common-range rule and the
        // explicit "central values only, not error bars" instruction.
        draw_cell(0, 0, g_ylo, g_yhi, offscale);
        auto* gd = new TGraphErrors(); gd->SetMarkerStyle(20); gd->SetMarkerColor(kBlack);
        gd->SetLineColor(kBlack);
        auto* gf = new TGraph(); gf->SetLineColor(kRed + 1); gf->SetLineWidth(2);
        // Below the panel label (which is long on this canvas -- "inclusive (all pT, all eta)")
        // but still above the plateau points, which sit at y = 1.
        auto* leg = new TLegend(0.45, 0.69, 0.95, 0.81);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.032);
        leg->AddEntry(gd, "measurement", "lp");
        leg->AddEntry(gf, "fit", "l");
        leg->Draw();
        // 0.024, not 0.030: the longest headline (overlay + Medium + `powerlaw_fixedRp`) ran off
        // the right edge of the 900 px canvas and lost the method name.
        c.cd(0);   // headline goes on the CANVAS, in the reserved strip -- not inside the pad
        TLatex st; st.SetNDC(); st.SetTextFont(42); st.SetTextSize(0.024);
        // No internal method name on the canvas -- the equation in the annotation box identifies
        // the fitted function to the audience (the method name only tags the output directory).
        st.DrawLatex(0.06, 1. - 0.62 * ihfrac, Form("%s, %s,  %s / plateau",
                                       cfg.sample_text.c_str(), wp_text.c_str(),
                                       quantity_tex.c_str()));
        const std::string png = odir + tag + "_dr_fit_" + method + "_inclusive.png";
        c.SaveAs(png.c_str());
        std::cout << "  wrote " << png << "\n";
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
    {
        const double Rp = (step == 3) ? 0.5 : 0.3;
        constexpr double kFlatTol = 1e-3;   // shape: allowed deviation from 1 beyond Rp
        // Compact support is a property of the stored formula, read back from the file itself
        // rather than assumed from the method name.
        bool compact = (method == "interp");
        {
            const FittedFunc Fp = LoadFunc(ffit, step, 0, 0);
            if (Fp.f) compact = TString(Fp.f->GetExpFormula()).Contains("TMath::Max");
        }
        std::ofstream os(odir + "readback_check.txt");
        os << "# Read-back verification of " << DrCorrFitFile(cfg, use_tight_wp, step, method)
           << "\n# performed in a SEPARATE ROOT process from the fit (that is the point).\n"
           << "# shape has compact support (exactly 1 for dR >= Rp by construction): "
           << (compact ? "YES" : "NO") << "\n"
           << "# PERSISTENCE: every value finite"
           << (compact ? ", and exactly 1 for dR >= Rp (a read-back returning 0 outside the"
                         " stored range -- the compiled-TF1 bug -- fails here)" : "") << "\n"
           << "# FLATNESS   : |f(dR) - 1| <= " << kFlatTol << " for dR >= Rp = " << Rp << "\n\n";
        int nchecked = 0, nfail_persist = 0, nfail_flat = 0;
        double worst_flat = 0.;
        for (int iy = 0; iy <= npt; ++iy) {
            for (int iz = 0; iz <= neta; ++iz) {
                if (!(iy == 0 && iz == 0) && (iy == 0 || iz == 0)) continue;
                const FittedFunc F = LoadFunc(ffit, step, iy, iz);
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
                for (double x : {0.0, 0.5 * Rp, Rp, 0.9, 1.5, 3.0, 8.0, 20.0}) {
                    const double v = F.Eval(x);
                    line += Form("  f(%.2f)=%.5f", x, v);
                    if (!std::isfinite(v)) bad_persist = true;
                    if (compact && x >= rp_cell && std::fabs(v - 1.0) > 1e-9) bad_persist = true;
                    if (x >= Rp) {
                        worst_flat = std::max(worst_flat, std::fabs(v - 1.0));
                        if (std::fabs(v - 1.0) > kFlatTol) bad_flat = true;
                    }
                }
                if (bad_persist) { ++nfail_persist; line += "   <-- PERSISTENCE FAIL"; }
                if (bad_flat)    { ++nfail_flat;    line += "   <-- not flat beyond Rp"; }
                os << line << "\n";
            }
        }
        os << "\nchecked " << nchecked << " functions, " << nfail_persist << " FAILED persistence"
           << "\nflatness beyond Rp: " << nfail_flat << " function(s) exceed " << kFlatTol
           << ", worst |f-1| = " << Form("%.2e", worst_flat) << "\n";
        std::cout << "  read-back check: " << nchecked << " functions, " << nfail_persist
                  << " FAILED persistence; " << nfail_flat << " not flat beyond Rp (worst |f-1| = "
                  << Form("%.2e", worst_flat) << ")\n";
        if (nfail_persist) std::cout << "  ** PERSISTENCE FAILURE -- do NOT use these fits **\n";
    }

    ffit->Close();
    fh->Close();
}

// Convenience: every method for one (sample, WP, step).
void plot_dr_correction_fits_all(const std::string& sample = "pp_full", bool use_tight_wp = true,
                                 int step = 3)
{
    for (const std::string& m : {"powerlaw_fixedRp", "powerlaw_floatRp", "expo",
                                 "polyu_fixedRp", "interp"})
        plot_dr_correction_fits(sample, use_tight_wp, step, m);
}
