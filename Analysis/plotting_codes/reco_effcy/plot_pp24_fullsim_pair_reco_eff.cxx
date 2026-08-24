// =================================================================================================
// plot_pp24_fullsim_pair_reco_eff.cxx
//
// Draw the pair reconstruction efficiency THAT THE CROSS-SECTION ACTUALLY APPLIES, i.e. the
// contents of `pair_reco_eff_pp24_full.root` written by build_pp24_fullsim_pair_reco_eff.C and
// read at fill time by Utilities/PairRecoEffEvaluator.h.
//
// WHY THIS FILE EXISTS. build_pp24_fullsim_pair_reco_eff.C persists the efficiency but draws
// nothing, and PythiaFullsimRecoEffPlotter draws only the FINE diagnostic axes (pT_bins_80 x 48
// pair-eta x 20 dR). So until now the applied 8 x 9 x 4 map had no picture anywhere -- the object
// entering every corrected yield was the one object nobody could look at. These figures are what
// the internal note shows for the reconstruction efficiency.
//
// WHAT IS DRAWN (one set per working point; NOMINAL WP = Tight):
//   1. pair_reco_eff_vs_dr.png            eps_reco vs DeltaR, one curve per pair-pT bin, plus the
//                                         inclusive curve. THE figure for the pair treatment: if
//                                         eps_reco were eps_1 x eps_2 with no close-by effect it
//                                         would not need a DeltaR axis at all.
//   2. pair_reco_eff_vs_pair_pt.png       eps_reco vs pair pT, one curve per DeltaR bin, plus the
//                                         DeltaR-integrated curve. Log x -- the axis is log-binned.
//   3. pair_reco_eff_vs_pair_eta.png      eps_reco vs pair eta, one curve per DeltaR bin, plus the
//                                         DeltaR-integrated curve.
//   4. pair_reco_eff_map_dr_integrated.png  the DeltaR-integrated 2D map eps_reco(pair pT, pair eta),
//                                         COLZ with the cell value printed -- this is the FIRST
//                                         fallback level of the evaluator.
//   5. pair_reco_eff_3d_cells.png         the applied 3D map in full: one panel per pair-eta bin
//                                         (9 -> 3x3), eps_reco vs pair pT, one curve per DeltaR bin.
//                                         Cells with no measure are absent, which is exactly the
//                                         set that routes to the fallback.
// A companion `pair_reco_eff_values.txt` carries every drawn number, so a value quoted in the note
// can be traced without reopening the ROOT file.
//
// EVERY CURVE IS BUILT FROM THE STORED NUMERATOR AND DENOMINATOR, never from the ratio histogram:
// a projection of a ratio is not the ratio of the projections. The numerator and denominator 3D
// histograms live in the same file for this reason (build macro, "keep the raw ingredients").
// Errors are binomial (TH1::Divide option "B"), matching how the applied map itself was built;
// the histograms are MC-weighted, so these are the usual weighted-binomial approximation.
//
// The pair-pT and pair-eta axes are verified against ParamsSet::pair_pt_coarse_bins and
// CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap and the macro THROWS on a mismatch --
// the same guard PairRecoEffEvaluator applies, so a stale efficiency file cannot be plotted as if
// it were the applied one. Axis edges are never retyped here.
//
// Usage (from Analysis/plotting_codes/reco_effcy/):
//   root -l -b -q 'plot_pp24_fullsim_pair_reco_eff.cxx+()'                 // FULL sample, Tight
//   root -l -b -q 'plot_pp24_fullsim_pair_reco_eff.cxx+(true, false)'      // FULL sample, Medium
//   root -l -b -q 'plot_pp24_fullsim_pair_reco_eff.cxx+(false)'            // TEST sample, Tight
// =================================================================================================

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include "../../MuonObjectsParamsAndHelpers/ParamsSet.h"
#include "../../RDFBasedHistFilling/CommonEffcyConfig.h"

namespace {

// One colour/marker per curve, reused by all four 1D figures. Kept short and high-contrast;
// eight pair-pT bins is the longest series.
const std::vector<Color_t> kCols = {kBlack,   kRed + 1,  kBlue + 1,   kGreen + 2, kMagenta + 1,
                                    kOrange + 7, kCyan + 2, kViolet + 1, kSpring - 6};
const std::vector<Style_t> kMrks = {20, 21, 22, 23, 33, 34, 29, 47, 24};

template <typename T>
T* Get(TFile* f, const std::string& n)
{
    T* o = dynamic_cast<T*>(f->Get(n.c_str()));
    if (!o) throw std::runtime_error("plot_pp24_fullsim_pair_reco_eff: missing '" + n + "' in "
                                     + f->GetName()
                                     + " -- rerun build_pp24_fullsim_pair_reco_eff.C");
    return o;
}

// The evaluator's guard, repeated here so a stale file cannot be drawn as the applied one.
void CheckCanonicalBinning(const TH3D* h, const std::string& path)
{
    static const ParamsSet pms{};
    static const CommonEffcyConfig cfg{};

    std::vector<double> eta;
    eta.push_back(cfg.pair_eta_proj_ranges_coarse_incl_gap.front().first);
    for (const auto& r : cfg.pair_eta_proj_ranges_coarse_incl_gap) eta.push_back(r.second);

    auto check = [&path](const TAxis* ax, const std::vector<double>& edges, const char* what) {
        if (ax->GetNbins() != (int)edges.size() - 1)
            throw std::runtime_error("plot_pp24_fullsim_pair_reco_eff: " + std::string(what)
                + " axis of " + path + " has " + std::to_string(ax->GetNbins())
                + " bins but the canonical binning has " + std::to_string(edges.size() - 1)
                + " -- stale efficiency file");
        for (size_t i = 0; i < edges.size(); ++i) {
            const double got = (i + 1 <= (size_t)ax->GetNbins())
                             ? ax->GetBinLowEdge(i + 1) : ax->GetBinUpEdge(ax->GetNbins());
            if (std::fabs(got - edges[i]) > 1e-6)
                throw std::runtime_error("plot_pp24_fullsim_pair_reco_eff: " + std::string(what)
                    + " edge " + std::to_string(i) + " is " + std::to_string(got) + " but "
                    + std::to_string(edges[i]) + " canonically -- stale efficiency file");
        }
    };
    check(h->GetXaxis(), pms.pair_pt_coarse_bins, "pair-pT");
    check(h->GetYaxis(), eta,                     "pair-eta");
    std::vector<double> dr(cfg.dr_bins_edges_for_reco_effcy.begin(),
                           cfg.dr_bins_edges_for_reco_effcy.end());
    check(h->GetZaxis(), dr, "dR");
}

// A bin with an empty denominator has NO measured efficiency. Drawing it at 0 would read as
// "the efficiency is zero here", which is a different statement, so such bins are pushed far
// below the frame and simply do not appear. This is the same set of cells the evaluator routes
// to its fallback.
const double kNoMeasure = -999.0;

// A bin in which every generated pair was reconstructed has eff = 1 exactly, and one in which
// none was has eff = 0; the binomial formula returns an error of ZERO in both cases. Printed or
// drawn as "1.00 +- 0.00" on a cell holding 1e-7 of the simulated weight, that reads as a precise
// measurement, which is the opposite of the truth. For those two degenerate cases the drawn
// uncertainty is instead the width of the 1-sigma WILSON interval, which at p = 1 (or 0) runs
// from N/(N+1) to 1 and is therefore 1/(N+1) wide, with the effective entry count
// N = (content/error)^2 of the denominator. It is a one-sided excursion, not a symmetric error;
// only the DISPLAY is affected, never the applied efficiency file.
void DegenerateBinError(TH1* h, int b, double den_content, double den_error)
{
    const double v = h->GetBinContent(b);
    if (h->GetBinError(b) > 0.) return;
    if (v > 0. && v < 1.) return;            // a genuine zero-error interior value: leave alone
    if (den_content <= 0. || den_error <= 0.) return;
    const double neff = (den_content / den_error) * (den_content / den_error);
    h->SetBinError(b, 1.0 / (neff + 1.0));
}

// eps = num/den for one projection, with binomial errors. `axis` selects which axis survives;
// the other two are restricted to the given INCLUSIVE bin range (0 for that axis = all bins).
// The bin ranges are passed to TH3::Projection* explicitly rather than through TAxis::SetRange,
// so the caller's histograms are never modified and there is no dependence on ROOT's
// "default arguments mean the currently set range" behaviour.
TH1D* Efficiency(TH3D* num, TH3D* den, char axis, const std::string& name,
                 int x1 = 0, int x2 = 0, int y1 = 0, int y2 = 0, int z1 = 0, int z2 = 0)
{
    const int nx = num->GetNbinsX(), ny = num->GetNbinsY(), nz = num->GetNbinsZ();
    if (!x1) { x1 = 1; x2 = nx; }
    if (!y1) { y1 = 1; y2 = ny; }
    if (!z1) { z1 = 1; z2 = nz; }

    auto project = [&](TH3D* h, const std::string& nm) {
        TH1D* p = (axis == 'x') ? h->ProjectionX(nm.c_str(), y1, y2, z1, z2, "e")
                : (axis == 'y') ? h->ProjectionY(nm.c_str(), x1, x2, z1, z2, "e")
                                : h->ProjectionZ(nm.c_str(), x1, x2, y1, y2, "e");
        p->SetDirectory(nullptr);
        return p;
    };

    TH1D* n = project(num, name + "_num");
    TH1D* d = project(den, name + "_den");
    auto* e = static_cast<TH1D*>(n->Clone(name.c_str()));
    e->SetDirectory(nullptr);
    e->Divide(n, d, 1.0, 1.0, "B");
    for (int b = 1; b <= e->GetNbinsX(); ++b) {
        if (d->GetBinContent(b) <= 0.) { e->SetBinContent(b, kNoMeasure); e->SetBinError(b, 0.); }
        else DegenerateBinError(e, b, d->GetBinContent(b), d->GetBinError(b));
    }
    delete n;
    delete d;
    return e;
}

void StyleCurve(TH1* h, int i)
{
    const Color_t c = kCols[i % kCols.size()];
    h->SetLineColor(c);
    h->SetMarkerColor(c);
    h->SetMarkerStyle(kMrks[i % kMrks.size()]);
    h->SetMarkerSize(1.1);
    h->SetLineWidth(2);
    h->SetStats(0);
}

// One value for the companion table. A cell with NO generated pair prints "-"; a cell that has
// a generated pair but no reconstructed one is a genuine measured zero and prints 0.0000 --
// conflating the two would hide the second, which is the one that costs signal.
std::string Val(const TH1* h, int b)
{
    return (h->GetBinContent(b) <= kNoMeasure / 2.0)
         ? std::string("-")
         : std::string(Form("%.4f +- %.4f", h->GetBinContent(b), h->GetBinError(b)));
}

// The companion table is plain text, so the labels must be too -- ROOT-LaTeX markup belongs on
// the canvas, not in a file meant to be read with `cat`.
std::string PtLabelAscii(const TAxis* a, int i)
{
    return Form("%.4g < pT_pair < %.4g GeV", a->GetBinLowEdge(i), a->GetBinUpEdge(i));
}
std::string DrLabelAscii(const TAxis* a, int i)
{
    return Form("%.1f < dR < %.1f", a->GetBinLowEdge(i), a->GetBinUpEdge(i));
}
std::string EtaLabelAscii(const TAxis* a, int j)
{
    return Form("%.1f < eta_pair < %.1f", a->GetBinLowEdge(j), a->GetBinUpEdge(j));
}

std::string PtLabel(const TAxis* a, int i)
{
    return Form("%.4g < p_{T}^{pair} < %.4g GeV", a->GetBinLowEdge(i), a->GetBinUpEdge(i));
}
std::string DrLabel(const TAxis* a, int i)
{
    return Form("%.1f < #DeltaR < %.1f", a->GetBinLowEdge(i), a->GetBinUpEdge(i));
}
std::string EtaLabel(const TAxis* a, int i)
{
    return Form("%.1f < #eta^{pair} < %.1f", a->GetBinLowEdge(i), a->GetBinUpEdge(i));
}

// Sample identity, drawn once per canvas. No method names, no file paths (ATLAS plotting P1a).
void DrawHeader(const std::string& wp_text, const std::string& sample_text,
                double y = 0.955, double size = 0.030, bool wp_here = true)
{
    TLatex t;
    t.SetNDC();
    t.SetTextFont(42);
    t.SetTextSize(size);
    t.DrawLatex(0.13, y, sample_text.c_str());
    if (!wp_here) return;
    t.SetTextAlign(31);
    t.DrawLatex(0.97, y, wp_text.c_str());
}

// The defining equation of the plotted symbol. eps_reco^pair is specific to this analysis -- it
// is a PAIR efficiency, not eps_1 x eps_2 -- so the canvas must say what it is.
void DrawDefinition(double y, const std::string& wp_text = "", double size = 0.026)
{
    TLatex t;
    t.SetNDC();
    t.SetTextFont(42);
    t.SetTextSize(size);
    t.DrawLatex(0.13, y,
        "#varepsilon_{reco}^{pair} = N_{reco} / N_{gen}, both counted in the single-b signal "
        "region");
    if (wp_text.empty()) return;
    t.SetTextAlign(31);
    t.DrawLatex(0.97, y, wp_text.c_str());
}

}  // namespace

// =================================================================================================
void plot_pp24_fullsim_pair_reco_eff(bool use_full_sample = true, bool tight_WP = true)
{
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat("4.2f");

    const std::string dir = use_full_sample
        ? "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_full_sample/"
        : "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_test_sample/";
    const std::string sfx     = use_full_sample ? "_full" : "";
    const std::string in_path = dir + "pair_reco_eff_pp24" + sfx + ".root";
    const std::string wp      = tight_WP ? "tight" : "medium";
    const std::string wp_text = tight_WP ? "Tight muons" : "Medium muons";
    const std::string sample_text = std::string("Pythia8 full simulation, ")
                                  + (use_full_sample ? "pp 2024 conditions"
                                                     : "pp 2024 conditions (test sample)")
                                  + ", single-b opposite-sign pairs";

    std::unique_ptr<TFile> fin(TFile::Open(in_path.c_str(), "READ"));
    if (!fin || fin->IsZombie())
        throw std::runtime_error("plot_pp24_fullsim_pair_reco_eff: cannot open " + in_path
                                 + " -- run build_pp24_fullsim_pair_reco_eff.C first");

    TH3D* num = Get<TH3D>(fin.get(), "h_pair_reco_eff_num_" + wp);
    TH3D* den = Get<TH3D>(fin.get(), "h_pair_reco_eff_denom");
    TH2D* eff2d = Get<TH2D>(fin.get(), "h_pair_reco_eff_" + wp + "_dr_integrated");
    CheckCanonicalBinning(num, in_path);

    const TAxis* apt  = num->GetXaxis();
    const TAxis* aeta = num->GetYaxis();
    const TAxis* adr  = num->GetZaxis();
    const int npt = apt->GetNbins(), neta = aeta->GetNbins(), ndr = adr->GetNbins();

    const std::string outdir = dir + "plots/pp24_reco_effcy_plots/" + wp + "/applied/";
    gSystem->mkdir(outdir.c_str(), kTRUE);
    std::ofstream tab(outdir + "pair_reco_eff_values.txt");
    tab << "# Applied single-b pair reconstruction efficiency, " << wp << " WP\n"
        << "# source " << in_path << "\n"
        << "# eps = N(both muons reconstructed, WP, reco pair in the signal region)\n"
        << "#     / N(truth pair in the signal region);  binomial errors\n";

    // ---------------------------------------------------------------------------------------
    // 1. eps vs DeltaR, one curve per pair-pT bin (+ inclusive)
    // ---------------------------------------------------------------------------------------
    {
        TCanvas c("c_dr", "", 900, 760);
        c.SetTicks(1, 1);
        // The legend gets its OWN space above the frame. A transparent legend placed inside the
        // frame was drawn straight through the lowest curve (atlas-plotting.md: a legend backing
        // does not solve an overlap -- reserve space instead).
        c.SetTopMargin(0.26);
        TLegend leg(0.13, 0.755, 0.97, 0.925);
        leg.SetNColumns(3);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        leg.SetTextSize(0.025);

        TH1D* incl = Efficiency(num, den, 'z', "eff_dr_incl");
        StyleCurve(incl, 0);
        incl->SetTitle("");
        incl->GetXaxis()->SetTitle("#DeltaR(#mu,#mu)");
        incl->GetYaxis()->SetTitle("#varepsilon_{reco}^{pair}");
        incl->GetYaxis()->SetRangeUser(-0.04, 1.05);   // a measured zero sits ON the axis otherwise
        incl->SetLineWidth(3);
        incl->Draw("PE");
        leg.AddEntry(incl, "all p_{T}^{pair}", "lp");
        tab << "\n[eps vs dR, integrated over pair pT and pair eta]\n";
        for (int k = 1; k <= ndr; ++k)
            tab << "  " << DrLabelAscii(adr, k) << " : " << Val(incl, k) << "\n";

        tab << "\n[eps vs dR, per pair-pT bin]\n";
        for (int i = 1; i <= npt; ++i) {
            TH1D* h = Efficiency(num, den, 'z', Form("eff_dr_pt%d", i), i, i);
            StyleCurve(h, i);
            h->Draw("PE SAME");
            leg.AddEntry(h, PtLabel(apt, i).c_str(), "lp");
            tab << "  " << PtLabelAscii(apt, i) << "\n";
            for (int k = 1; k <= ndr; ++k)
                tab << "     " << DrLabelAscii(adr, k) << " : " << Val(h, k) << "\n";
        }
        leg.Draw();
        DrawHeader(wp_text, sample_text, 0.975, 0.030, false);
        DrawDefinition(0.945, wp_text);
        c.SaveAs((outdir + "pair_reco_eff_vs_dr.png").c_str());
    }

    // ---------------------------------------------------------------------------------------
    // 2. eps vs pair pT, one curve per DeltaR bin (+ DeltaR-integrated)
    // ---------------------------------------------------------------------------------------
    {
        TCanvas c("c_pt", "", 900, 760);
        c.SetTicks(1, 1);
        c.SetLogx();   // the pair-pT axis is log-binned (ParamsSet::pair_pt_coarse_bins)
        c.SetTopMargin(0.26);
        TLegend leg(0.13, 0.755, 0.97, 0.925);
        leg.SetNColumns(3);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        leg.SetTextSize(0.025);

        TH1D* incl = Efficiency(num, den, 'x', "eff_pt_incl");
        StyleCurve(incl, 0);
        incl->SetTitle("");
        incl->GetXaxis()->SetTitle("p_{T}^{pair} [GeV]");
        // No SetMoreLogLabels(): on an axis starting at 8 GeV it prints "9" and "10" on top of
        // each other. The decade labels are enough.
        incl->GetXaxis()->SetNoExponent();
        incl->GetYaxis()->SetTitle("#varepsilon_{reco}^{pair}");
        incl->GetYaxis()->SetRangeUser(-0.04, 1.05);   // a measured zero sits ON the axis otherwise
        incl->SetLineWidth(3);
        incl->Draw("PE");
        leg.AddEntry(incl, "all #DeltaR", "lp");
        tab << "\n[eps vs pair pT, integrated over pair eta and dR]\n";
        for (int i = 1; i <= npt; ++i)
            tab << "  " << PtLabelAscii(apt, i) << " : " << Val(incl, i) << "\n";

        tab << "\n[eps vs pair pT, per dR bin]\n";
        for (int k = 1; k <= ndr; ++k) {
            TH1D* h = Efficiency(num, den, 'x', Form("eff_pt_dr%d", k), 0, 0, 0, 0, k, k);
            StyleCurve(h, k);
            h->Draw("PE SAME");
            leg.AddEntry(h, DrLabel(adr, k).c_str(), "lp");
            tab << "  " << DrLabelAscii(adr, k) << "\n";
            for (int i = 1; i <= npt; ++i)
                tab << "     " << PtLabelAscii(apt, i) << " : " << Val(h, i) << "\n";
        }
        leg.Draw();
        DrawHeader(wp_text, sample_text, 0.975, 0.030, false);
        DrawDefinition(0.945, wp_text);
        c.SaveAs((outdir + "pair_reco_eff_vs_pair_pt.png").c_str());
    }

    // ---------------------------------------------------------------------------------------
    // 3. eps vs pair eta, one curve per DeltaR bin (+ DeltaR-integrated)
    // ---------------------------------------------------------------------------------------
    {
        TCanvas c("c_eta", "", 900, 760);
        c.SetTicks(1, 1);
        c.SetTopMargin(0.26);
        TLegend leg(0.13, 0.755, 0.97, 0.925);
        leg.SetNColumns(3);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        leg.SetTextSize(0.025);

        TH1D* incl = Efficiency(num, den, 'y', "eff_eta_incl");
        StyleCurve(incl, 0);
        incl->SetTitle("");
        incl->GetXaxis()->SetTitle("#eta^{pair}");
        incl->GetYaxis()->SetTitle("#varepsilon_{reco}^{pair}");
        incl->GetYaxis()->SetRangeUser(-0.04, 1.05);   // a measured zero sits ON the axis otherwise
        incl->SetLineWidth(3);
        incl->Draw("PE");
        leg.AddEntry(incl, "all #DeltaR", "lp");
        tab << "\n[eps vs pair eta, integrated over pair pT and dR]\n";
        for (int j = 1; j <= neta; ++j)
            tab << "  " << EtaLabelAscii(aeta, j) << " : " << Val(incl, j) << "\n";

        tab << "\n[eps vs pair eta, per dR bin]\n";
        for (int k = 1; k <= ndr; ++k) {
            TH1D* h = Efficiency(num, den, 'y', Form("eff_eta_dr%d", k), 0, 0, 0, 0, k, k);
            StyleCurve(h, k);
            h->Draw("PE SAME");
            leg.AddEntry(h, DrLabel(adr, k).c_str(), "lp");
            tab << "  " << DrLabelAscii(adr, k) << "\n";
            for (int j = 1; j <= neta; ++j)
                tab << "     " << EtaLabelAscii(aeta, j) << " : " << Val(h, j) << "\n";
        }
        leg.Draw();
        DrawHeader(wp_text, sample_text, 0.975, 0.030, false);
        DrawDefinition(0.945, wp_text);
        c.SaveAs((outdir + "pair_reco_eff_vs_pair_eta.png").c_str());
    }

    // ---------------------------------------------------------------------------------------
    // 4. the DeltaR-integrated 2D map -- the evaluator's first fallback level
    // ---------------------------------------------------------------------------------------
    {
        TCanvas c("c_map", "", 1000, 800);
        c.SetTicks(1, 1);
        c.SetLogx();
        c.SetRightMargin(0.14);
        c.SetTopMargin(0.14);
        TH2D* den2d = Get<TH2D>(fin.get(), "h_pair_reco_eff_denom2d");
        auto* h = static_cast<TH2D*>(eff2d->Clone("map_dr_integrated"));
        h->SetDirectory(nullptr);
        h->SetTitle("");
        h->GetXaxis()->SetTitle("p_{T}^{pair} [GeV]");
        h->GetXaxis()->SetNoExponent();
        h->GetYaxis()->SetTitle("#eta^{pair}");
        h->GetZaxis()->SetTitle("#varepsilon_{reco}^{pair}");
        // The colour axis spans the POPULATED range, not [0,1]: with every cell between ~0.5 and
        // 1 a fixed [0,1] palette leaves half its dynamic range dead and renders the whole map
        // in one band. The uncertainty is drawn with each value ("TEXTE"), so a cell carrying
        // almost no simulated weight -- the top pair-pT row holds 1e-5 to 1e-7 of it -- cannot be
        // read as a precise measurement.
        // The span is taken over BOTH working points, which live in the same file, so the Tight
        // and Medium maps share one colour scale: a per-file range would make the working-point
        // comparison a different scale on each figure.
        double zlo = 1.0, zhi = 0.0;
        for (const std::string w : {std::string("tight"), std::string("medium")}) {
            auto* hw = Get<TH2D>(fin.get(), "h_pair_reco_eff_" + w + "_dr_integrated");
            for (int i = 1; i <= hw->GetNbinsX(); ++i)
                for (int j = 1; j <= hw->GetNbinsY(); ++j) {
                    const double v = hw->GetBinContent(i, j);
                    if (v <= 0.) continue;
                    zlo = std::min(zlo, v);
                    zhi = std::max(zhi, v);
                }
        }
        h->SetMinimum(std::max(0.0, std::floor(zlo * 20.0) / 20.0));
        h->SetMaximum(std::min(1.0, std::ceil(zhi * 20.0) / 20.0));
        // Same degenerate-cell treatment as the 1D projections (see DegenerateBinError).
        TH2D* den2 = den2d;
        int n_degenerate = 0;
        for (int i = 1; i <= h->GetNbinsX(); ++i)
            for (int j = 1; j <= h->GetNbinsY(); ++j) {
                if (h->GetBinError(i, j) > 0.) continue;
                const double v = h->GetBinContent(i, j);
                if (v > 0. && v < 1.) continue;
                const double dc = den2->GetBinContent(i, j), de = den2->GetBinError(i, j);
                if (dc <= 0. || de <= 0.) continue;
                const double neff = (dc / de) * (dc / de);
                h->SetBinError(i, j, 1.0 / (neff + 1.0));
                ++n_degenerate;
            }
        std::cout << "  2D map cells at eff = 0 or 1 given a Wilson-interval width : "
                  << n_degenerate << std::endl;
        h->SetMarkerSize(0.85);
        h->Draw("COLZ TEXTE");
        DrawHeader(wp_text, sample_text, 0.955, 0.030, false);
        DrawDefinition(0.925, wp_text);
        c.SaveAs((outdir + "pair_reco_eff_map_dr_integrated.png").c_str());

        tab << "\n[eps, dR-integrated 2D map (rows = pair pT, cols = pair eta)]\n";
        for (int i = 1; i <= npt; ++i)
            for (int j = 1; j <= neta; ++j)
                tab << "  " << PtLabelAscii(apt, i) << " x " << EtaLabelAscii(aeta, j) << " : "
                    << (den2d->GetBinContent(i, j) <= 0.
                            ? std::string("-")
                            : std::string(Form("%.4f +- %.4f", h->GetBinContent(i, j),
                                               h->GetBinError(i, j))))
                    << "\n";
    }

    // ---------------------------------------------------------------------------------------
    // 5. the applied 3D map in full: one panel per pair-eta bin, eps vs pair pT per DeltaR bin
    // ---------------------------------------------------------------------------------------
    {
        // nrows >= ncols, nrows ~ sqrt(N)  (subplot-layout convention)
        const int ncol = (int)std::ceil(std::sqrt((double)neta));
        const int nrow = (int)std::ceil((double)neta / ncol);
        // A reserved header strip, so the sample line cannot land on top of the first row of
        // panels -- the panels of a Divide()d canvas run right up to y = 1.
        const int hstrip = 108;
        TCanvas c("c_cells", "", 420 * ncol, 340 * nrow + hstrip);
        const double fh = (double)hstrip / (340.0 * nrow + hstrip);
        TPad* head = new TPad("head", "", 0.0, 1.0 - fh, 1.0, 1.0);
        TPad* body = new TPad("body", "", 0.0, 0.0, 1.0, 1.0 - fh);
        head->SetFillStyle(0);
        body->SetFillStyle(0);
        head->Draw();
        body->Draw();
        head->cd();
        DrawHeader(wp_text, sample_text, 0.80, 0.26, false);
        DrawDefinition(0.58, wp_text, 0.22);
        body->cd();
        body->Divide(ncol, nrow);

        int n_measured = 0, n_zero = 0, n_total = 0;
        tab << "\n[applied 3D cells: eps per (pair pT, pair eta, dR); '-' = no measure]\n";
        // ONE legend for the whole figure, in the header strip. A per-panel legend has nowhere
        // safe to sit: the lower-left corner it used is exactly where the measured-zero cell and
        // the low-efficiency points land, and it was drawn straight through them.
        TLegend* head_leg = nullptr;
        for (int j = 1; j <= neta; ++j) {
            body->cd(j);
            gPad->SetTicks(1, 1);
            gPad->SetLogx();
            gPad->SetLeftMargin(0.13);
            gPad->SetBottomMargin(0.13);

            tab << "  " << EtaLabelAscii(aeta, j) << "\n";
            for (int k = 1; k <= ndr; ++k) {
                TH1D* h = Efficiency(num, den, 'x', Form("cell_eta%d_dr%d", j, k),
                                     0, 0, j, j, k, k);
                StyleCurve(h, k);
                h->SetTitle("");
                h->GetXaxis()->SetTitle("p_{T}^{pair} [GeV]");
                h->GetXaxis()->SetNoExponent();
                h->GetYaxis()->SetTitle("#varepsilon_{reco}^{pair}");
                h->GetYaxis()->SetRangeUser(-0.04, 1.19);   // a measured zero sits ON the axis otherwise
                h->GetXaxis()->SetTitleSize(0.050);
                h->GetXaxis()->SetLabelSize(0.042);
                h->GetYaxis()->SetTitleSize(0.050);
                h->GetYaxis()->SetLabelSize(0.042);
                h->GetYaxis()->SetTitleOffset(1.15);
                h->Draw(k == 1 ? "PE" : "PE SAME");
                if (j == 1) {
                    if (!head_leg) {
                        head_leg = new TLegend(0.13, 0.06, 0.97, 0.40);
                        head_leg->SetNColumns(ndr);
                        head_leg->SetBorderSize(0);
                        head_leg->SetFillStyle(0);
                        head_leg->SetTextSize(0.24);
                    }
                    head_leg->AddEntry(h, DrLabel(adr, k).c_str(), "lp");
                }

                tab << "     " << DrLabelAscii(adr, k) << " :";
                for (int i = 1; i <= npt; ++i) {
                    ++n_total;
                    const double e = h->GetBinContent(i);
                    if (e <= kNoMeasure / 2.0) { tab << " -"; }          // no generated pair
                    else if (e > 0.)           { ++n_measured; tab << " " << Form("%.4f", e); }
                    else                       { ++n_zero; tab << " 0.0000"; }  // measured zero
                }
                tab << "\n";
            }
            TLatex t;
            t.SetNDC();
            t.SetTextFont(42);
            t.SetTextSize(0.048);
            t.DrawLatex(0.17, 0.90, EtaLabel(aeta, j).c_str());
        }
        head->cd();
        if (head_leg) head_leg->Draw();
        c.SaveAs((outdir + "pair_reco_eff_3d_cells.png").c_str());

        std::cout << "  3D cells delivering a non-zero efficiency : " << n_measured << " of "
                  << n_total << "  (" << (n_total ? 100.0 * n_measured / n_total : 0.0) << " %)\n"
                  << "  of the rest, " << n_zero << (n_zero == 1 ? " has" : " have")
                  << " a generated pair but no reconstructed one (a measured zero) and "
                  << (n_total - n_measured - n_zero)
                  << " have no generated pair at all; both route to the dR-integrated fallback"
                  << std::endl;
        tab << "\n[cell coverage] " << n_measured << " of " << n_total
            << " 3D cells deliver a non-zero efficiency ("
            << (n_total ? 100.0 * n_measured / n_total : 0.0) << " %); " << n_zero
            << (n_zero == 1 ? " is a measured zero, " : " are a measured zero, ") << (n_total - n_measured - n_zero)
            << " have no generated pair\n";
    }

    tab.close();
    std::cout << "plot_pp24_fullsim_pair_reco_eff: wrote 5 PNGs + pair_reco_eff_values.txt to\n  "
              << outdir << std::endl;
}
