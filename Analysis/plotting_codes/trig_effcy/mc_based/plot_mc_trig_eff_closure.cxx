// =============================================================================
// plot_mc_trig_eff_closure.cxx
//
// The MC-closure figure of the pp24 2mu4 trigger correction
// (docs/tracking/mc_trig_eff_closure.md §3.4; inputs from FillMCTrigEffClosure.cxx).
//
// ONE PNG per sample version. Subplots = pair-eta bins
// (CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap, 9 bins -> 3x3). Each subplot:
//
//   upper pad   dsigma/dpT^pair, three series:
//                 (1) all pairs, NO trigger requirement
//                 (2) pairs passing 2mu4, corrected with the exponential dR-correction form
//                 (3) pairs passing 2mu4, corrected with the polynomial dR-correction form
//   lower pad   the closure ratio (2)/(1) and (3)/(1); it must be 1 everywhere.
//
// The ratio carries the CONDITIONAL (binomial-correct) error, because the numerator is a
// re-weighted SUBSET of the denominator (mc_trigger_efficiency.md R12). It is computed by the
// repo's ONE implementation, SetConditionalRatioErrors in dr_correction_ratio.h -- Var = A - R*B
// with the k=n boundary fallback. The spectra carry their own sqrt(sum w^2).
//
// The pair-pT axis IS the binning the dR correction's cells are defined on (doc D2) -- 8 log bins
// by default, the 4-bin variant under MCTRIGEFF_PAIRPT_4BIN. Both are read from
// MCTrigEffPairPt::Edges, never retyped.
//
// Usage (from Analysis/plotting_codes/trig_effcy/mc_based/):
//   root -l -b -q 'plot_mc_trig_eff_closure.cxx+("pp_full", true)'
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
#include <TSystem.h>

#include "dr_correction_sample_cfg.h"
#include "dr_correction_apply.h"   // the applied dR ceiling, so the drawn value cannot drift
#include "dr_correction_ratio.h"   // THE conditional-error implementation -- never a second copy
#include "../../../Utilities/CommonLogYRange.h"
#include "../../../Utilities/MCTrigEffPairPtBinning.h"
#include "../../../RDFBasedHistFilling/CommonEffcyConfig.h"

namespace {

// Series identity. NO internal method name ever reaches the canvas
// (.claude/conventions/atlas-plotting.md P1a) -- the two corrected series are named by the
// FUNCTIONAL FORM of the dR correction, whose equations are drawn in the header strip.
struct Series {
    std::string key;        // histogram token (fill-stage method name); never drawn
    std::string legend;     // what the reader sees
    Color_t     colour;
    Style_t     marker;
};
const std::vector<Series> kCorrected = {
    {"expo",          "2mu4, corrected (exponential #varepsilon_{#DeltaR})", kRed + 1,  20},
    {"polyu_fixedRp", "2mu4, corrected (polynomial #varepsilon_{#DeltaR})",  kBlue + 1, 22},
};
const Color_t kUncorrColour = kBlack;
const Style_t kUncorrMarker = 21;

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

}  // namespace

// =============================================================================
void plot_mc_trig_eff_closure(const std::string& sample = "pp_full", bool use_tight_wp = true)
{
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);

    const DrCorrSample cfg = GetDrCorrSample(sample, use_tight_wp);
    const std::string wp_suf  = DrCorrWpSuffix(use_tight_wp);
    const std::string wp_text = use_tight_wp ? "Tight muons" : "Medium muons";

    const std::string in_path = cfg.mc_dir + "mc_trig_eff_closure_" + cfg.mc_label + wp_suf
                              + MCTrigEffPairPt::FileSuffix() + ".root";
    TFile* fin = TFile::Open(in_path.c_str(), "READ");
    if (!fin || fin->IsZombie())
        throw std::runtime_error("plot_mc_trig_eff_closure: cannot open " + in_path
                                 + " -- run FillMCTrigEffClosure first");

    const std::string outdir = cfg.out_base + "closure/";
    gSystem->mkdir(outdir.c_str(), kTRUE);

    static const CommonEffcyConfig ecfg{};
    const auto& eta_ranges = ecfg.pair_eta_proj_ranges_coarse_incl_gap;
    const int neta = static_cast<int>(eta_ranges.size());
    // 3x3 for the 9 pair-eta bins; nrows >= ncols, nrows ~ sqrt(N) (memory feedback_subplot_layout)
    const int ncol = static_cast<int>(std::ceil(std::sqrt((double)neta)));
    const int nrow = static_cast<int>(std::ceil((double)neta / ncol));

    for (const auto& V : kVersions) {
        const std::string pre = "h_closure_" + V.key + "_";
        TH2D* h_den = Get<TH2D>(fin, pre + "den");
        std::map<std::string, TH2D*> h_num, h_A, h_B;
        for (const auto& S : kCorrected) {
            h_num[S.key] = Get<TH2D>(fin, pre + "num_"  + S.key);
            h_A  [S.key] = Get<TH2D>(fin, pre + "numA_" + S.key);
            h_B  [S.key] = Get<TH2D>(fin, pre + "numB_" + S.key);
        }
        if (h_den->GetNbinsY() != neta)
            throw std::runtime_error("plot_mc_trig_eff_closure: the histograms carry "
                + std::to_string(h_den->GetNbinsY()) + " pair-eta bins but CommonEffcyConfig has "
                + std::to_string(neta) + " -- stale input?");

        // ---------------- build every panel's histograms first (shared ranges need them all) ----
        std::vector<TH1D*> spec_unc(neta, nullptr);
        std::vector<std::map<std::string, TH1D*>> spec_cor(neta), ratio(neta);
        std::vector<TH1*> for_range;
        double rmin = 1.0, rmax = 1.0;   // seeded at the closure reference
        int n_offscale = 0;

        for (int iz = 1; iz <= neta; ++iz) {
            TH1D* den = Row(h_den, iz, Form("%s_den_eta%d", V.key.c_str(), iz));
            for (const auto& S : kCorrected) {
                TH1D* num = Row(h_num[S.key], iz, Form("%s_num_%s_eta%d", V.key.c_str(),
                                                       S.key.c_str(), iz));
                TH1D* A   = Row(h_A  [S.key], iz, Form("%s_A_%s_eta%d", V.key.c_str(),
                                                       S.key.c_str(), iz));
                TH1D* B   = Row(h_B  [S.key], iz, Form("%s_B_%s_eta%d", V.key.c_str(),
                                                       S.key.c_str(), iz));
                // CLOSURE RATIO with the conditional (binomial-correct) error. Built BEFORE the
                // width scaling: dividing two spectra scaled by the same widths is identical, but
                // A and B are sums of squared weights and do NOT scale the same way, so the error
                // must be formed on the unscaled histograms.
                //
                // Var = A - R*B, NOT A - B. B = sum_fired a^2 * p carries the PREDICTED per-pair
                // probability p; the variance needs the TRUE one, whose first-order estimate is
                // R*p (R = this bin's closure ratio). A - B is the special case R = 1 -- i.e. it
                // assumes the closure already holds, which is the very thing being tested and
                // which misses by ~27 % here. It also goes NEGATIVE in cells where R > 1, and a
                // negative variance silently became a zero error, which both consumers drop.
                // SetConditionalRatioErrors is the repo's single implementation and carries the
                // k=n boundary fallback for exactly that case.
                TH1D* r = static_cast<TH1D*>(num->Clone(Form("%s_ratio_%s_eta%d", V.key.c_str(),
                                                             S.key.c_str(), iz)));
                r->SetDirectory(nullptr);
                r->Divide(den);
                SetConditionalRatioErrors(r, den, A, B);
                delete A; delete B;
                ratio[iz - 1][S.key] = r;
                num->Scale(1.0, "width");
                spec_cor[iz - 1][S.key] = num;
                for_range.push_back(num);
            }
            den->Scale(1.0, "width");
            spec_unc[iz - 1] = den;
            for_range.push_back(den);
        }

        // Ratio-pad range: ONE range for all nine panels, derived from the data by PERCENTILE
        // rather than by min/max.
        //
        // Why not min/max: a single statistics-poor high-pair-pT cell that is nonetheless
        // "informative" by any error-based gate (measured: 3.007 +- 0.621, 21 % relative) stretched
        // the signal-region figure to [0.13, 3.89] and squeezed the band the figure exists to show
        // -- the ~1.27 offset and the barrel/endcap step -- into ~16 % of the pad height. Gating on
        // the error cannot remove such a point, because its error is genuinely small.
        //
        // Percentiles keep the frame on the bulk while the excursions remain VISIBLE (drawn, just
        // off the top or bottom of the frame) and, crucially, COUNTED: the count is printed on the
        // canvas below, so nothing is silently dropped.
        {
            std::vector<double> vals;
            for (int iz = 0; iz < neta; ++iz)
                for (const auto& S : kCorrected) {
                    TH1D* r = ratio[iz][S.key];
                    for (int b = 1; b <= r->GetNbinsX(); ++b) {
                        const double y = r->GetBinContent(b), e = r->GetBinError(b);
                        if (y <= 0. || e <= 0.) continue;
                        if (e > 0.5 * y) continue;   // no information at all
                        vals.push_back(y - e);
                        vals.push_back(y + e);
                    }
                }
            if (!vals.empty()) {
                std::sort(vals.begin(), vals.end());
                const size_t n = vals.size();
                rmin = vals[(size_t)(0.02 * (n - 1))];
                rmax = vals[(size_t)(0.98 * (n - 1))];
            }
        }
        // The range always CONTAINS 1: perfect closure is the reference the reader measures the
        // result against, so an axis that excludes it would hide how far away the result is. A
        // minimum width keeps a panel from being blown up by rounding when everything coincides.
        rmin = std::min(rmin, 1.0); rmax = std::max(rmax, 1.0);
        if (rmax - rmin < 0.10) { const double m = 0.5 * (rmin + rmax); rmin = m - 0.05; rmax = m + 0.05; }
        const double pad = 0.08 * (rmax - rmin);
        rmin -= pad; rmax += pad;
        for (int iz = 0; iz < neta; ++iz)
            for (const auto& S : kCorrected) {
                TH1D* r = ratio[iz][S.key];
                for (int b = 1; b <= r->GetNbinsX(); ++b) {
                    const double y = r->GetBinContent(b);
                    if (y > 0. && (y < rmin || y > rmax)) ++n_offscale;
                }
            }

        ApplyCommonLogYRange(for_range);

        // ---------------- draw ----------------------------------------------------------------
        // Header strip reserved at the TOP of the canvas: with nine dense two-pad panels no
        // quadrant is free, so the legend and the defining equations get their own space rather
        // than being laid over data (.claude/conventions/atlas-plotting.md, legend placement).
        const double head = 0.17;
        TCanvas c(("c_closure_" + V.key).c_str(), "", 1800, 1750);

        std::vector<TPad*> upper, lower;
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

            auto* pu = new TPad(Form("pu_%s_%d", V.key.c_str(), iz), "", x1, ysplit, x2, ytop);
            auto* pl = new TPad(Form("pl_%s_%d", V.key.c_str(), iz), "", x1, ybot,   x2, ysplit);
            for (TPad* p : {pu, pl}) {
                p->SetLeftMargin(0.235); p->SetRightMargin(0.02);
                p->SetLogx(1);
                p->SetTicks(1, 1);
            }
            pu->SetTopMargin(0.04); pu->SetBottomMargin(0.02);
            pl->SetTopMargin(0.02); pl->SetBottomMargin(0.32);
            pu->Draw(); pl->Draw();
            upper.push_back(pu); lower.push_back(pl);

            // ---- spectra ----
            pu->cd();
            pu->SetLogy(1);
            TH1D* du = spec_unc[iz - 1];
            du->SetTitle("");
            du->GetYaxis()->SetTitle("d#sigma/dp_{T}^{pair} [nb/GeV]");
            du->GetYaxis()->SetTitleSize(0.070); du->GetYaxis()->SetLabelSize(0.060);
            du->GetYaxis()->SetTitleOffset(1.55);
            du->GetXaxis()->SetLabelSize(0.0);
            du->SetMarkerColor(kUncorrColour); du->SetLineColor(kUncorrColour);
            du->SetMarkerStyle(kUncorrMarker); du->SetMarkerSize(1.1);
            du->Draw("PE");
            for (const auto& S : kCorrected) {
                TH1D* h = spec_cor[iz - 1][S.key];
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
                TH1D* r = ratio[iz - 1][S.key];
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
                // The decade labels plus minor ticks are unambiguous on their own.
                r->Draw(first ? "PE" : "PE SAME");
                first = false;
            }
            // Perfect closure, the reference the whole figure is measured against.
            auto* ln = new TLine(spec_unc[iz - 1]->GetXaxis()->GetXmin(), 1.0,
                                 spec_unc[iz - 1]->GetXaxis()->GetXmax(), 1.0);
            ln->SetLineStyle(2); ln->SetLineColor(kGray + 2);
            ln->Draw();
        }

        // ---- header strip: identity, legend, and the defining equations --------------------
        c.cd();
        auto* hd = new TLatex();
        hd->SetNDC(); hd->SetTextFont(42);
        hd->SetTextSize(0.0180);
        hd->DrawLatex(0.035, 0.982,
                      (cfg.sample_text + ",  " + wp_text + ",  " + V.headline).c_str());
        hd->SetTextSize(0.0150);
        hd->DrawLatex(0.035, 0.918,
                      "w = 1 / [#varepsilon(p_{T,1},q#eta_{1}) "
                      "#varepsilon(p_{T,2},q#eta_{2}) #varepsilon_{#DeltaR}(#DeltaR)],"
                      "   #varepsilon = single-muon mu4 efficiency from the data tag-and-probe "
                      "turn-on");
        hd->DrawLatex(0.035, 0.897,
                      Form("#varepsilon_{#DeltaR}(#DeltaR) = f(#DeltaR)/C for #DeltaR < %g and 1 "
                           "above;   exponential  f = C + A e^{-(#DeltaR/#lambda)^{p}}", 
                           DrCorrectionEvaluator::kDrMax));
        hd->DrawLatex(0.035, 0.876,
                      "polynomial  f = C + u^{2}(a_{2} + a_{3}u + a_{4}u^{2}),"
                      "   u #equiv max(0, 1 - #DeltaR/R_{p})");
        if (n_offscale > 0)
            hd->DrawLatex(0.035, 0.855,
                          Form("%d ratio point(s) outside the lower-pad range", n_offscale));

        // One canvas-level legend, laid out as a single row in the reserved strip.
        auto* leg = new TLegend(0.035, 0.943, 0.990, 0.974);
        leg->SetNColumns(3);
        leg->SetBorderSize(0); leg->SetFillStyle(0);
        leg->SetTextFont(42); leg->SetTextSize(0.0150);
        leg->AddEntry(spec_unc[0], "no trigger requirement", "PE");
        for (const auto& S : kCorrected)
            leg->AddEntry(spec_cor[0][S.key], S.legend.c_str(), "PE");
        leg->Draw();

        const std::string png = outdir + V.file + ".png";
        c.SaveAs(png.c_str());
        std::cout << "  wrote " << png << std::endl;

        // ---- the number the figure is about, in text -----------------------------------------
        std::cout << "  inclusive closure (" << V.key << "):";
        for (const auto& S : kCorrected) {
            double n = 0., d = 0.;
            for (int iz = 1; iz <= neta; ++iz) {
                n += h_num[S.key]->Integral(1, h_num[S.key]->GetNbinsX(), iz, iz);
                d += h_den->Integral(1, h_den->GetNbinsX(), iz, iz);
            }
            std::cout << "  " << S.key << " = " << (d > 0 ? n / d : -1.);
        }
        std::cout << std::endl;
    }

    fin->Close();
    std::cout << "done." << std::endl;
}
