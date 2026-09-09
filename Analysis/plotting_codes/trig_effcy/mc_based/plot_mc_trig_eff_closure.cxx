// =============================================================================
// plot_mc_trig_eff_closure.cxx
//
// The MC-closure figure of the pp24 2mu4 trigger correction
// (docs/tracking/mc_trigeff_dr_binning_approaches.md §PP-4; machinery and physics:
//  mc_trig_eff_closure.md; inputs from FillMCTrigEffClosure.cxx).
//
// ONE PNG per (sample version, figure set). Subplots = pair-eta bins
// (CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap, 9 bins -> 3x3). Each subplot:
//
//   upper pad   dsigma/dpT^pair -- the no-trigger-requirement series, plus one corrected series
//               per entry of the figure set.
//   lower pad   the closure ratio of each corrected series over the no-trigger one; it must be 1
//               everywhere.
//
// THE FOUR CELL-GROUPING APPROACHES (`plateau_mode`; user 2026-08-24). `nocorr` (8 pair-pT x 9
// pair-eta), `nocorr_ptmerge` (7x9), `nocorr_etamerge` (8x3) and `nocorr_etamerge_ptmerge` (7x3)
// each go to their OWN subdirectory, named by the repo's own mode helper (DrCorrPlateauModeDir).
// The DELIVERED series is the same cross-method cascade in all four -- expo fit -> polyu fit where
// the exponential is rejected -> the interpolation where both are -> the measured values where
// every tier is -- so what differs between the figures is the CELL GROUPING and not the cascade.
//
// THE UN-MERGED REFERENCE PRODUCES A SECOND FIGURE SET, in a `separate_fit_forms/` subdirectory of
// its own (user, 2026-08-24): the `expo` and `polyu_fixedRp` corrections drawn as two SEPARATE
// series. That is a different question from the delivered correction -- how far apart the two
// parametric forms are, cell by cell (mc_trigger_efficiency.md R26, OPEN) -- and mixing it into the
// deliverable's directory is what the separation prevents.
//
// THE X AXIS IS THE pp24 CROSS-SECTION's, NOT THE CORRECTION's (doc D8, superseding
// mc_trig_eff_closure.md D2): ParamsSet::pT_bins_150, 16 log bins over 9-150 GeV, IN EVERY
// APPROACH -- because the question the comparison answers is which approach corrects the
// cross-section most accurately, and four figures on four different x-axes could not be compared.
// The axis comes from the filled histogram itself here; the fill stage is the one place that reads
// ParamsSet.
//
// THE APPLIED SINGLE-MUON EFFICIENCY IS eps_MC (user, 2026-08-18; closure doc D4): the closure is
// self-contained, so what it tests is the dR correction and its fit alone. The fill stage also
// books an eps^nc_data numerator, which is the printed diagnostic, NOT drawn here.
//
// SCALES (user, 2026-08-24). The RATIO pad range is SHARED BY THE TWO PNGs of a figure set -- the
// closure is the same dimensionless test on two samples, so a y-position must mean the same number
// in both files. The SPECTRUM pad is NOT shared: those are cross sections of two different samples
// and differ in normalization by physics, so each PNG keeps its own range (still one range for all
// nine panels within it, ApplyCommonLogYRange property 2).
// The ratio frame CONTAINS EVERY DRAWN POINT (central values; error bars may clip), the rule
// ApplyCommonLogYRange already states for the spectra: a cell drawn in the upper pad and missing
// from the ratio pad below it reads as absent data, and the cells that sit far from 1 are this
// figure's findings.
//
// The ratio carries the CONDITIONAL (binomial-correct) error, because the numerator is a
// re-weighted SUBSET of the denominator (mc_trigger_efficiency.md R12). It is computed by the
// repo's ONE implementation, SetConditionalRatioErrors in dr_correction_ratio.h -- Var = A - R*B
// with the k=n boundary fallback. The spectra carry their own sqrt(sum w^2).
//
// Usage (from Analysis/plotting_codes/trig_effcy/mc_based/):
//   root -l -b -q 'plot_mc_trig_eff_closure.cxx+("pp_full", true)'
//   root -l -b -q 'plot_mc_trig_eff_closure.cxx+("pp_full", true, "nocorr_ptmerge")'
//   root -l -b -q 'plot_mc_trig_eff_closure.cxx+("pp_full", true, "nocorr_etamerge")'
//   MCTRIGEFF_PAIRPT_4BIN=1 root -l -b -q 'plot_mc_trig_eff_closure.cxx+("pp_full", true)'
// =============================================================================

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TNamed.h>
#include <TSystem.h>

#include "dr_correction_sample_cfg.h"
#include "dr_correction_apply.h"   // the applied dR ceiling, so the drawn value cannot drift
#include "dr_correction_ratio.h"   // THE conditional-error implementation -- never a second copy
#include "../../../Utilities/CommonLogYRange.h"
#include "../../../Utilities/MCTrigEffPairPtBinning.h"
#include "../../../RDFBasedHistFilling/CommonEffcyConfig.h"
#include "../../../Utilities/PairEtaPanelBins.h"

namespace {

// Series identity. NO internal method name ever reaches the canvas
// (.claude/conventions/atlas-plotting.md P1a) -- the corrected series are named by the FUNCTIONAL
// FORM of the dR correction, whose equations are drawn in the header strip.
struct Series {
    std::string key;        // histogram token (fill-stage series key); never drawn
    std::string legend;     // what the reader sees
    Color_t     colour;
    Style_t     marker;
};

// THE DELIVERED correction -- one series, in every approach. Naming it after a fit form would
// misdescribe the cells routed to another form or to the interpolation.
const std::vector<Series> kCorrectedCascade = {
    {"cascade", "2mu4, corrected (#varepsilon_{#DeltaR})", kRed + 1, 20},
};
// The un-merged reference's SECOND figure set: the two parametric forms side by side.
const std::vector<Series> kCorrectedPerForm = {
    {"expo",          "2mu4, corrected (exponential #varepsilon_{#DeltaR})", kRed + 1,  20},
    {"polyu_fixedRp", "2mu4, corrected (polynomial #varepsilon_{#DeltaR})",  kBlue + 1, 22},
};
const Color_t kUncorrColour = kBlack;
const Style_t kUncorrMarker = 21;

// The subdirectory the per-form set goes to, relative to the approach's own directory.
const char* kPerFormSubdir = "separate_fit_forms/";

struct VersionCfg {
    std::string key;        // histogram token
    std::string file;       // PNG basename
    std::string headline;   // what the sample IS, in physics words
};
const std::vector<VersionCfg> kVersions = {
    {"all_os", "closure_pair_pt_all_opposite_sign",
     "all opposite-sign muon pairs"},
    {"signal", "closure_pair_pt_single_b_signal_cuts",
     "opposite-sign pairs in the single-b signal region"},
};

std::string EtaLabel(const std::pair<float, float>& r)
{
    return Form("%.1f < #eta^{pair} < %.1f", r.first, r.second);
}

template <typename T>
T* Get(TFile* f, const std::string& n)
{
    T* o = dynamic_cast<T*>(f->Get(n.c_str()));
    if (!o) throw std::runtime_error("plot_mc_trig_eff_closure: missing '" + n + "' in "
                                     + f->GetName());
    return o;
}

// ProjectionX of one pair-eta row, detached from the file.
TH1D* Row(TH2D* h, int iz, const std::string& nm)
{
    TH1D* p = h->ProjectionX(nm.c_str(), iz, iz, "e");
    p->SetDirectory(nullptr);
    return p;
}

// =============================================================================================
// ONE FIGURE SET: the same series list drawn for both sample versions, sharing one ratio range.
// The un-merged approach calls this twice (the cascade set, then the per-form set); every other
// approach calls it once.
// =============================================================================================
void DrawClosureSet(TFile* fin, const DrCorrSample& cfg, const std::string& wp_text,
                    const std::string& plateau_mode, const std::vector<Series>& kCorrected,
                    const std::string& outdir, bool per_form_set,
                    const std::string& set_headline, const std::string& tag)
{
    static const CommonEffcyConfig ecfg{};
    const auto& eta_ranges = ecfg.pair_eta_proj_ranges_coarse_incl_gap;
    const int neta = static_cast<int>(eta_ranges.size());
    // 3x3 for the 9 pair-eta bins; nrows >= ncols, nrows ~ sqrt(N) (memory feedback_subplot_layout)
    const int ncol = static_cast<int>(std::ceil(std::sqrt((double)neta)));
    const int nrow = static_cast<int>(std::ceil((double)neta / ncol));

    gSystem->mkdir(outdir.c_str(), kTRUE);

    // ==========================================================================================
    // THE RATIO RANGE IS SHARED BY THE TWO FIGURES; THE SPECTRUM RANGE IS NOT (user, 2026-08-24).
    //
    // The closure ratio is the same dimensionless test on two samples (all opposite-sign pairs /
    // the single-b signal region), and the reader compares the two files against each other, so a
    // ratio-pad y-position must mean the same number in both -- the failure ApplyCommonLogYRange's
    // property 2 removes WITHIN a figure, reappearing BETWEEN them. The spectra are the opposite
    // case: they are cross sections of two different samples whose normalizations differ by
    // physics, so forcing one range on both would waste most of each frame; each PNG keeps its own
    // (still common to its nine panels).
    // Hence every histogram of both versions is built FIRST -- the shared ratio range needs them
    // all -- and only then is anything drawn.
    // ==========================================================================================
    struct VersionData {
        TH2D*                                     h_den = nullptr;
        std::map<std::string, TH2D*>              h_num;
        std::vector<TH1D*>                        spec_unc;   // [pair-eta bin]
        std::vector<std::map<std::string, TH1D*>> spec_cor, ratio;
    };
    std::vector<VersionData> vdata(kVersions.size());
    std::vector<double> ratio_pts;    // every ratio point drawn in EITHER figure

    // ---------------- build every panel's histograms first (the shared ranges need them all) ---
    for (size_t iv = 0; iv < kVersions.size(); ++iv) {
        const VersionCfg& V = kVersions[iv];
        VersionData&      D = vdata[iv];

        const std::string pre = "h_closure_" + V.key + "_";
        D.h_den = Get<TH2D>(fin, pre + "den");
        std::map<std::string, TH2D*> h_A, h_B;
        // `_epsmc_`: the numerator weighted with the APPLIED eps_MC. The file also holds an
        // `_epsdata_` numerator -- the diagnostic (closure doc D4) -- which is deliberately NOT
        // drawn.
        for (const auto& S : kCorrected) {
            D.h_num[S.key] = Get<TH2D>(fin, pre + "num_epsmc_"  + S.key);
            h_A[S.key]     = Get<TH2D>(fin, pre + "numA_epsmc_" + S.key);
            h_B[S.key]     = Get<TH2D>(fin, pre + "numB_epsmc_" + S.key);
        }
        if (D.h_den->GetNbinsY() != neta)
            throw std::runtime_error("plot_mc_trig_eff_closure: the histograms carry "
                + std::to_string(D.h_den->GetNbinsY()) + " pair-eta bins but CommonEffcyConfig has "
                + std::to_string(neta) + " -- stale input?");
        // BIN COUNT IS NOT ENOUGH. Panel `iz` is drawn with the label EtaLabel(eta_ranges[iz-1])
        // read LIVE from CommonEffcyConfig, while its content is Y-bin `iz` of a file written
        // earlier. A count check passes even when the two describe different ranges -- which is
        // exactly what happened on 2026-09-07, when the coarse pair-eta outer bins moved
        // +-(2.0,2.4) -> +-(2.0,2.2) to follow the pair-level gap cut: every closure file on disk
        // still had 9 bins over +-2.4 and would have been relabelled silently. Check the EDGES,
        // per panel, through the shared guard (muon_gap_cuts_acceptance.md F18).
        PairEtaPanels::CheckPanelsMatchFiducialCut();
        for (int iz = 1; iz <= neta; ++iz)
            PairEtaPanels::CheckAxisAligned(D.h_den->GetYaxis(), eta_ranges[iz - 1], iz, iz,
                                            "plot_mc_trig_eff_closure");

        D.spec_unc.assign(neta, nullptr);
        D.spec_cor.resize(neta);
        D.ratio.resize(neta);
        std::vector<TH1*> for_range;   // every spectrum of THIS figure -- its own range (above)

        for (int iz = 1; iz <= neta; ++iz) {
            TH1D* den = Row(D.h_den, iz, Form("%s%s_den_eta%d", tag.c_str(), V.key.c_str(), iz));
            for (const auto& S : kCorrected) {
                TH1D* num = Row(D.h_num[S.key], iz, Form("%s%s_num_%s_eta%d", tag.c_str(),
                                                         V.key.c_str(), S.key.c_str(), iz));
                TH1D* A   = Row(h_A[S.key], iz, Form("%s%s_A_%s_eta%d", tag.c_str(),
                                                     V.key.c_str(), S.key.c_str(), iz));
                TH1D* B   = Row(h_B[S.key], iz, Form("%s%s_B_%s_eta%d", tag.c_str(),
                                                     V.key.c_str(), S.key.c_str(), iz));
                // CLOSURE RATIO with the conditional (binomial-correct) error. Built BEFORE the
                // width scaling: dividing two spectra scaled by the same widths is identical, but
                // A and B are sums of squared weights and do NOT scale the same way, so the error
                // must be formed on the unscaled histograms.
                //
                // Var = A - R*B, NOT A - B. B = sum_fired a^2 * p carries the PREDICTED per-pair
                // probability p; the variance needs the TRUE one, whose first-order estimate is
                // R*p (R = this bin's closure ratio). A - B is the special case R = 1 -- i.e. it
                // assumes the closure already holds, which is the very thing being tested.
                // SetConditionalRatioErrors is the repo's single implementation and carries the
                // k=n boundary fallback.
                TH1D* r = static_cast<TH1D*>(num->Clone(Form("%s%s_ratio_%s_eta%d", tag.c_str(),
                                                             V.key.c_str(), S.key.c_str(), iz)));
                r->SetDirectory(nullptr);
                r->Divide(den);
                SetConditionalRatioErrors(r, den, A, B);
                delete A; delete B;
                D.ratio[iz - 1][S.key] = r;
                for (int b = 1; b <= r->GetNbinsX(); ++b)
                    if (r->GetBinContent(b) > 0.) ratio_pts.push_back(r->GetBinContent(b));
                num->Scale(1.0, "width");
                D.spec_cor[iz - 1][S.key] = num;
                for_range.push_back(num);
            }
            den->Scale(1.0, "width");
            D.spec_unc[iz - 1] = den;
            for_range.push_back(den);
        }
        ApplyCommonLogYRange(for_range);   // per FIGURE, not across the two (see the block above)
    }

    // ---------------- the shared ratio range --------------------------------------------------
    // RATIO PAD: EVERY DRAWN POINT IS INSIDE THE FRAME (user, 2026-08-24) -- the same rule
    // ApplyCommonLogYRange states for the spectra (property 1), now applied to the ratio.
    // Central values only: in the statistics-poor top pair-pT cells the conditional error reaches
    // ~0.6, and sizing the frame to contain every error BAR would compress the band the figure is
    // about for no gain -- a marker inside the frame with its bar clipped at the edge is the same
    // convention the spectra use.
    double rmin = 1.0, rmax = 1.0;   // seeded at the closure reference
    if (!ratio_pts.empty()) {
        rmin = *std::min_element(ratio_pts.begin(), ratio_pts.end());
        rmax = *std::max_element(ratio_pts.begin(), ratio_pts.end());
    }
    // The range always CONTAINS 1: perfect closure is the reference the reader measures the
    // result against, so an axis that excludes it would hide how far away the result is. A
    // minimum width keeps a panel from being blown up by rounding when everything coincides.
    rmin = std::min(rmin, 1.0); rmax = std::max(rmax, 1.0);
    if (rmax - rmin < 0.10) { const double m = 0.5 * (rmin + rmax); rmin = m - 0.05; rmax = m + 0.05; }
    const double rpad = 0.08 * (rmax - rmin);
    rmin -= rpad; rmax += rpad;

    // ---------------- draw --------------------------------------------------------------------
    for (size_t iv = 0; iv < kVersions.size(); ++iv) {
        const VersionCfg& V = kVersions[iv];
        VersionData&      D = vdata[iv];

        // Guard, not a range policy: with the rule above nothing can fall outside the frame, so a
        // non-zero count here means the range and the drawing have drifted apart.
        int n_offscale = 0;
        for (int iz = 0; iz < neta; ++iz)
            for (const auto& S : kCorrected) {
                TH1D* r = D.ratio[iz][S.key];
                for (int b = 1; b <= r->GetNbinsX(); ++b) {
                    const double y = r->GetBinContent(b);
                    if (y > 0. && (y < rmin || y > rmax)) ++n_offscale;
                }
            }

        // Header strip reserved at the TOP of the canvas: with nine dense two-pad panels no
        // quadrant is free, so the legend and the defining equations get their own space rather
        // than being laid over data (.claude/conventions/atlas-plotting.md, legend placement).
        const double head = 0.17;
        TCanvas c((tag + "c_closure_" + V.key).c_str(), "", 1800, 1750);

        for (int iz = 1; iz <= neta; ++iz) {
            // c.cd() EVERY iteration: a new TPad is adopted by gPad, and gPad is whatever pad was
            // last cd()'d into. Without this the second panel is created INSIDE the first panel's
            // ratio pad and only one panel ever appears on the canvas.
            c.cd();
            const int col = (iz - 1) % ncol, row = (iz - 1) / ncol;
            const double x1 = col / (double)ncol, x2 = (col + 1) / (double)ncol;
            const double ytop = 1.0 - head - row * (1.0 - head) / nrow;
            const double ybot = 1.0 - head - (row + 1) * (1.0 - head) / nrow;
            const double ysplit = ybot + 0.36 * (ytop - ybot);   // 64 % spectrum / 36 % ratio

            auto* pu = new TPad(Form("%spu_%s_%d", tag.c_str(), V.key.c_str(), iz), "",
                                x1, ysplit, x2, ytop);
            auto* pl = new TPad(Form("%spl_%s_%d", tag.c_str(), V.key.c_str(), iz), "",
                                x1, ybot,   x2, ysplit);
            for (TPad* p : {pu, pl}) {
                p->SetLeftMargin(0.235); p->SetRightMargin(0.02);
                p->SetLogx(1);
                p->SetTicks(1, 1);
            }
            pu->SetTopMargin(0.04); pu->SetBottomMargin(0.02);
            pl->SetTopMargin(0.02); pl->SetBottomMargin(0.32);
            pu->Draw(); pl->Draw();

            // ---- spectra ----
            pu->cd();
            pu->SetLogy(1);
            TH1D* du = D.spec_unc[iz - 1];
            du->SetTitle("");
            du->GetYaxis()->SetTitle("d#sigma/dp_{T}^{pair} [nb/GeV]");
            du->GetYaxis()->SetTitleSize(0.070); du->GetYaxis()->SetLabelSize(0.060);
            du->GetYaxis()->SetTitleOffset(1.55);
            du->GetXaxis()->SetLabelSize(0.0);
            du->SetMarkerColor(kUncorrColour); du->SetLineColor(kUncorrColour);
            du->SetMarkerStyle(kUncorrMarker); du->SetMarkerSize(1.1);
            du->Draw("PE");
            for (const auto& S : kCorrected) {
                TH1D* h = D.spec_cor[iz - 1][S.key];
                h->SetMarkerColor(S.colour); h->SetLineColor(S.colour);
                h->SetMarkerStyle(S.marker); h->SetMarkerSize(1.1);
                h->Draw("PE SAME");
            }
            auto* tl = new TLatex();
            tl->SetNDC(); tl->SetTextFont(42); tl->SetTextSize(0.068);
            tl->DrawLatex(0.28, 0.09, EtaLabel(eta_ranges[iz - 1]).c_str());

            // ---- closure ratio ----
            pl->cd();
            bool first = true;
            for (const auto& S : kCorrected) {
                TH1D* r = D.ratio[iz - 1][S.key];
                r->SetTitle("");
                r->SetMarkerColor(S.colour); r->SetLineColor(S.colour);
                r->SetMarkerStyle(S.marker); r->SetMarkerSize(1.1);
                r->GetYaxis()->SetRangeUser(rmin, rmax);
                r->GetYaxis()->SetTitle("corrected / no trigger");
                r->GetYaxis()->SetNdivisions(505);
                r->GetYaxis()->SetTitleSize(0.100); r->GetYaxis()->SetLabelSize(0.095);
                r->GetYaxis()->SetTitleOffset(1.08);
                r->GetXaxis()->SetTitle("p_{T}^{pair} [GeV]");
                r->GetXaxis()->SetTitleSize(0.125); r->GetXaxis()->SetLabelSize(0.105);
                r->GetXaxis()->SetTitleOffset(1.00);
                // NO SetMoreLogLabels(): the axis starts at 8 GeV, so it draws "9" and "10"
                // adjacent on a log scale and at this panel width they touch and read as "910".
                r->Draw(first ? "PE" : "PE SAME");
                first = false;
            }
            // Perfect closure, the reference the whole figure is measured against.
            auto* ln = new TLine(D.spec_unc[iz - 1]->GetXaxis()->GetXmin(), 1.0,
                                 D.spec_unc[iz - 1]->GetXaxis()->GetXmax(), 1.0);
            ln->SetLineStyle(2); ln->SetLineColor(kGray + 2);
            ln->Draw();
        }

        // ---- header strip: identity, legend, and the defining equations --------------------
        c.cd();
        auto* hd = new TLatex();
        hd->SetNDC(); hd->SetTextFont(42);
        hd->SetTextSize(0.0180);
        hd->DrawLatex(0.035, 0.982,
                      (cfg.sample_text + ",  " + wp_text + ",  " + V.headline + ",  "
                       + set_headline).c_str());
        hd->SetTextSize(0.0150);
        hd->DrawLatex(0.035, 0.918,
                      "w = 1 / [#varepsilon(p_{T,1},q#eta_{1}) "
                      "#varepsilon(p_{T,2},q#eta_{2}) #varepsilon_{#DeltaR}(#DeltaR)],"
                      "   #varepsilon = single-muon mu4 efficiency from the MC turn-on");
        hd->DrawLatex(0.035, 0.897,
                      Form("#varepsilon_{#DeltaR}(#DeltaR) = f(#DeltaR)/C for #DeltaR < %g and 1 "
                           "above;   exponential  f = C + A e^{-(#DeltaR/#lambda)^{p}}",
                           DrCorrectionEvaluator::kDrMax));
        hd->DrawLatex(0.035, 0.876,
                      "polynomial  f = C + u^{2}[A + a_{3}(u - 1) + a_{4}(u^{2} - 1)],"
                      "   u #equiv max(0, 1 - #DeltaR/R_{p}),   A #equiv f(0) - C");
        // Which f a cell uses is not a code detail in the cascade set -- it is part of what
        // eps_dR MEANS there, so the canvas has to define it.
        if (!per_form_set)
            hd->DrawLatex(0.035, 0.855,
                          "f = the exponential where its fit is accepted, the polynomial where it "
                          "is not, the linear interpolation of the measured points where neither "
                          "is, and those measured values themselves where none is");
        if (n_offscale > 0)
            hd->DrawLatex(0.035, per_form_set ? 0.855 : 0.834,
                          Form("%d ratio point(s) outside the lower-pad range", n_offscale));

        // One canvas-level legend, laid out as a single row in the reserved strip.
        auto* leg = new TLegend(0.035, 0.943, 0.990, 0.974);
        leg->SetNColumns(1 + (int)kCorrected.size());
        leg->SetBorderSize(0); leg->SetFillStyle(0);
        leg->SetTextFont(42); leg->SetTextSize(0.0150);
        leg->AddEntry(D.spec_unc[0], "no trigger requirement", "PE");
        for (const auto& S : kCorrected)
            leg->AddEntry(D.spec_cor[0][S.key], S.legend.c_str(), "PE");
        leg->Draw();

        const std::string png = outdir + V.file + ".png";
        c.SaveAs(png.c_str());
        std::cout << "  wrote " << png << std::endl;

        // ---- the number the figure is about, in text -----------------------------------------
        std::cout << "  inclusive closure [" << plateau_mode << "] (" << V.key << "):";
        for (const auto& S : kCorrected) {
            double n = 0., d = 0.;
            for (int iz = 1; iz <= neta; ++iz) {
                n += D.h_num[S.key]->Integral(1, D.h_num[S.key]->GetNbinsX(), iz, iz);
                d += D.h_den->Integral(1, D.h_den->GetNbinsX(), iz, iz);
            }
            std::cout << "  " << S.key << " = " << (d > 0 ? n / d : -1.);
        }
        std::cout << std::endl;
    }
    std::cout << "  shared ratio-pad range [" << rmin << ", " << rmax << "]  (" << outdir << ")"
              << std::endl;
}

}  // namespace

// =============================================================================
void plot_mc_trig_eff_closure(const std::string& sample = "pp_full", bool use_tight_wp = true,
                              const std::string& plateau_mode = "nocorr")
{
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);

    // Validates the token; DrCorrPlateauModeDir throws on anything else, so a typo can never
    // silently land in the wrong directory.
    const std::string mode_dir = DrCorrPlateauModeDir(plateau_mode);
    if (!DrCorrModeNoPlateau(plateau_mode))
        throw std::runtime_error("plot_mc_trig_eff_closure: the closure is defined for the "
                                 "no-plateau-correction approaches only, got '" + plateau_mode
                                 + "'");

    const DrCorrSample cfg = GetDrCorrSample(sample, use_tight_wp);
    const std::string wp_suf  = DrCorrWpSuffix(use_tight_wp);
    const std::string wp_text = use_tight_wp ? "Tight muons" : "Medium muons";

    const std::string in_path = cfg.mc_dir + "mc_trig_eff_closure_" + cfg.mc_label + wp_suf
                              + MCTrigEffPairPt::FileSuffix()
                              + DrCorrPlateauModeTag(plateau_mode) + ".root";
    TFile* fin = TFile::Open(in_path.c_str(), "READ");
    if (!fin || fin->IsZombie())
        throw std::runtime_error("plot_mc_trig_eff_closure: cannot open " + in_path
                                 + " -- run FillMCTrigEffClosure first");

    // What the cells of THIS approach are, in physics words, for the canvas headline. Read from
    // the fill stage's own provenance stamp rather than re-derived, so the figure cannot claim a
    // grouping the file was not built with.
    std::string cell_text = "#varepsilon_{#DeltaR} cells: see provenance";
    if (auto* prov = dynamic_cast<TNamed*>(fin->Get("provenance"))) {
        const TString t = prov->GetTitle();
        const Ssiz_t i = t.Index("eps_dR CELLS: ");
        if (i >= 0) {
            const Ssiz_t j = t.Index(" (", i);
            if (j > i) {
                TString cells = t(i + 14, j - i - 14);        // "8 pair pT x 9 pair eta"
                cells.ReplaceAll(" pair pT x ", " #times ");
                cells.ReplaceAll(" pair eta", "");
                cell_text = std::string(cells.Data())
                          + " #varepsilon_{#DeltaR} cells (p_{T}^{pair} #times #eta^{pair})";
            }
        }
    }

    const std::string outdir = cfg.out_base + "closure/" + mode_dir;

    // The DELIVERED correction: one series, in every approach.
    DrawClosureSet(fin, cfg, wp_text, plateau_mode, kCorrectedCascade, outdir,
                   /*per_form_set=*/false, cell_text, "casc_");

    // The un-merged reference's second figure set. Its series exist in the file only for that
    // approach (FillMCTrigEffClosure.cxx::WantsSeparateFormSeries), so the presence of the
    // histograms is what decides -- not a mode literal repeated here.
    const bool have_per_form =
        fin->Get("h_closure_all_os_num_epsmc_expo") &&
        fin->Get("h_closure_all_os_num_epsmc_polyu_fixedRp");
    if (have_per_form)
        DrawClosureSet(fin, cfg, wp_text, plateau_mode, kCorrectedPerForm,
                       outdir + kPerFormSubdir, /*per_form_set=*/true,
                       cell_text + ", the two parametric forms separately", "form_");

    fin->Close();
    std::cout << "done." << std::endl;
}
