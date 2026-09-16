// =================================================================================================
// plot_mc_singles_2d_effcy.cxx
//
// MC single-muon mu4 trigger efficiency as a 2D map in (q*eta, pT) -- the MC counterpart of the
// DATA tag-and-probe map produced by
//   Analysis/plotting_codes/trig_effcy/TrigEffPlotterPP.cxx :: plot2D_SingleMuonEffcy()  (line 494)
//   -> .../plots/pp_trigger_efficiency/mu4/no_corr/pp24/single_muon_effcy/
//      pt2nd_vs_q_eta2nd_2D_trig_effcy_sepr.png
// and drawn with the SAME conventions as that reference: x = q*eta (linear), y = pT (LOG), COLZ,
// z in [0, 1], mu+ on the LEFT pad and mu- on the RIGHT pad.
//
// PHYSICS (mc_trigger_efficiency.md Step 1 / §3.1): per reconstructed muon passing the singles
// selection,
//      eps_MC(q*eta, pT) = N(mu4 fired) / N(all)      -- weighted, per 2D cell,
// with NO trigger requirement on the denominator and no tag-and-probe pairing. Numerator and
// denominator are read straight from the RDF stage output (FillMCTrigEffHists.cxx:958-964), so
// every selection cut and every MC weight is exactly the nominal one -- this macro never touches
// the ntuples.
//
// BINNING: taken from the input histograms themselves, never retyped. Both axes are bin-for-bin
// IDENTICAL to the data map (verified: 184 x 41 bins, max |edge difference| = 0 on both axes),
// because MakeBinnings() in FillMCTrigEffHists.cxx copies the data construction:
//   x : ParamsSet::makeEtaTrigEffcyBinning(1)                 (= data "eta_bins_trig_effcy")
//   y : ParamsSet::pT_bins_8 ++ ParamsSet::pT_bins_60         (= data "pT_bins_single_muon")
// The y construction concatenates the two log-binned vectors WITHOUT dropping their shared 8.0 GeV
// edge, so bin 21 is a ZERO-WIDTH bin [8, 8]. That is deliberate and inherited from data
// (FillMCTrigEffHists.cxx:234-237, RDFBasedHistFillingData.cxx:315-317): it keeps the MC and data
// maps bin-for-bin comparable. No muon can ever land in it (no pT satisfies 8 <= pT < 8), so it is
// permanently empty and, being zero-width, invisible. It is handled here like any other cell with
// an empty denominator: left UNFILLED, never painted as efficiency 0.
//
// Charge combination (the "comb" canvas) sums the NUMERATORS and sums the DENOMINATORS and then
// divides -- the yield-weighted efficiency of the combined sample. It is NOT the average of the
// two per-charge efficiencies, which would weight a sparsely populated charge as heavily as the
// other one.
//
// Run (from this directory):
//   root -l -b -q 'plot_mc_singles_2d_effcy.cxx+("pp_full", true)'    // Tight  -> mc_based/
//   root -l -b -q 'plot_mc_singles_2d_effcy.cxx+("pp_full", false)'   // Medium -> mc_based_medium/
//   root -l -b -q 'plot_mc_singles_2d_effcy.cxx+("overlay", true)'    // Pb+Pb HIJING overlay
// =================================================================================================

#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

#include "dr_correction_sample_cfg.h"
#include "../../../AtlasStyle.C"

namespace {

// Sentinel for a cell whose denominator is empty: the efficiency is UNDEFINED there, not zero.
// It sits below the palette minimum (0), so COLZ skips the cell and leaves it white -- exactly
// how the data map renders its unpopulated high-pT corners. A cell with a populated denominator
// and an empty numerator keeps a genuine 0 and IS painted.
constexpr double kUndefined = -1.0;

template <typename T>
T* GetObj(TFile* f, const std::string& name) {
    auto* o = dynamic_cast<T*>(f->Get(name.c_str()));
    if (!o) throw std::runtime_error("plot_mc_singles_2d_effcy: missing '" + name + "' in " +
                                     f->GetName());
    return o;
}

// eff = num/den per cell; undefined where den <= 0.
TH2D* MakeEffMap(const TH2D* num, const TH2D* den, const std::string& name) {
    if (num->GetNbinsX() != den->GetNbinsX() || num->GetNbinsY() != den->GetNbinsY())
        throw std::runtime_error("plot_mc_singles_2d_effcy: numerator/denominator binning differs");

    auto* h = static_cast<TH2D*>(num->Clone(name.c_str()));
    h->SetDirectory(nullptr);
    h->Reset();

    for (int ix = 1; ix <= den->GetNbinsX(); ++ix) {
        for (int iy = 1; iy <= den->GetNbinsY(); ++iy) {
            const double d = den->GetBinContent(ix, iy);
            if (d <= 0.) {                       // includes the zero-width pT = 8 GeV bin
                h->SetBinContent(ix, iy, kUndefined);
                h->SetBinError(ix, iy, 0.);
                continue;
            }
            const double n = num->GetBinContent(ix, iy);
            const double eff = n / d;
            // Binomial error on a weighted ratio of a subset to its parent sample.
            const double err = std::sqrt(std::max(0., eff * (1. - eff)) / d);
            h->SetBinContent(ix, iy, eff);
            h->SetBinError(ix, iy, err);
        }
    }
    return h;
}

// One COLZ panel. Axis titles and ranges come from the histogram, so the panel can never
// disagree with the binning it draws.
void DrawEffPanel(TH2D* h, const std::string& charge_label) {
    gPad->SetLeftMargin(0.14);
    gPad->SetRightMargin(0.17);   // room for the palette AND its axis title
    gPad->SetBottomMargin(0.14);
    gPad->SetTopMargin(0.08);
    gPad->SetLogy();              // log-binned pT axis => log scale (as in the data map)

    h->SetMinimum(0.);
    h->SetMaximum(1.);
    h->GetXaxis()->SetTitle("q #upoint #eta");
    h->GetYaxis()->SetTitle("p_{T} [GeV]");
    h->GetZaxis()->SetTitle("#varepsilon(mu4)");
    // The pT axis spans less than one and a half decades (4-60 GeV), so the default log
    // labelling prints only "10". Label the intermediate values too, in plain notation.
    h->GetYaxis()->SetMoreLogLabels();
    h->GetYaxis()->SetNoExponent();
    h->GetXaxis()->SetTitleOffset(1.1);
    h->GetYaxis()->SetTitleOffset(1.1);
    h->GetZaxis()->SetTitleOffset(1.1);
    h->GetXaxis()->SetTitleSize(0.050);
    h->GetYaxis()->SetTitleSize(0.050);
    h->GetZaxis()->SetTitleSize(0.045);
    h->GetXaxis()->SetLabelSize(0.042);
    h->GetYaxis()->SetLabelSize(0.042);
    h->GetZaxis()->SetLabelSize(0.038);
    h->Draw("COLZ");

    TLatex t;
    t.SetNDC();
    t.SetTextFont(42);
    t.SetTextSize(0.055);
    t.DrawLatex(0.17, 0.935, charge_label.c_str());
}

void DrawHeadline(const std::string& text, double size) {
    TLatex t;
    t.SetNDC();
    t.SetTextFont(42);
    t.SetTextAlign(22);
    t.SetTextSize(size);
    t.DrawLatex(0.5, 0.965, text.c_str());
}

}  // namespace

void plot_mc_singles_2d_effcy(const std::string& sample = "pp_full", bool use_tight_wp = true,
                              int overlay_year = 24)
{
    gROOT->SetBatch(kTRUE);
    SetAtlasStyle();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetNumberContours(100);   // smooth palette for a 0-1 efficiency map

    // Sample identity (input dir, file label, output root, canvas headline) -- single source of
    // truth, shared with the whole dR-correction chain. Never hardcode a path here.
    const DrCorrSample cfg = GetDrCorrSample(sample, use_tight_wp, overlay_year);
    const std::string wp_suf  = DrCorrWpSuffix(use_tight_wp);
    const std::string wp_text = use_tight_wp ? "Tight muons" : "Medium muons";
    const std::string headline = cfg.sample_text + ", " + wp_text;

    const std::string in_path = DrCorrHistFile(cfg, use_tight_wp);
    TFile* fin = TFile::Open(in_path.c_str(), "READ");
    if (!fin || fin->IsZombie())
        throw std::runtime_error("plot_mc_singles_2d_effcy: cannot open " + in_path);

    TH2D* num_p = GetObj<TH2D>(fin, "h_mc_pt_vs_q_eta_num_muplus");
    TH2D* den_p = GetObj<TH2D>(fin, "h_mc_pt_vs_q_eta_denom_muplus");
    TH2D* num_m = GetObj<TH2D>(fin, "h_mc_pt_vs_q_eta_num_muminus");
    TH2D* den_m = GetObj<TH2D>(fin, "h_mc_pt_vs_q_eta_denom_muminus");

    TH2D* eff_p = MakeEffMap(num_p, den_p, "eff_muplus");
    TH2D* eff_m = MakeEffMap(num_m, den_m, "eff_muminus");

    // Charge-combined: sum numerators, sum denominators, THEN divide.
    auto* num_c = static_cast<TH2D*>(num_p->Clone("num_comb"));
    num_c->SetDirectory(nullptr);
    num_c->Add(num_m);
    auto* den_c = static_cast<TH2D*>(den_p->Clone("den_comb"));
    den_c->SetDirectory(nullptr);
    den_c->Add(den_m);
    TH2D* eff_c = MakeEffMap(num_c, den_c, "eff_comb");

    const std::string out_dir = cfg.out_base + "step1_singles_data_mc/";
    gSystem->mkdir(out_dir.c_str(), kTRUE);

    // ---------------- charge-separated: mu+ LEFT, mu- RIGHT ----------------
    TCanvas c_sepr("c_mc_singles_2d_sepr", "", 1300, 560);
    // Pads stop at y = 0.93 so the canvas headline has its own strip and cannot land inside a pad.
    TPad p_left("p_left", "", 0.00, 0.00, 0.50, 0.93);
    TPad p_right("p_right", "", 0.50, 0.00, 1.00, 0.93);
    p_left.Draw();
    p_right.Draw();
    p_left.cd();  DrawEffPanel(eff_p, "#mu^{+}");
    p_right.cd(); DrawEffPanel(eff_m, "#mu^{-}");
    c_sepr.cd();
    DrawHeadline(headline, 0.040);
    const std::string png_sepr = out_dir + "step1_eff_2d_pt_vs_q_eta_charge_sepr.png";
    c_sepr.SaveAs(png_sepr.c_str());

    // ---------------- charge-combined: one pad ----------------
    TCanvas c_comb("c_mc_singles_2d_comb", "", 750, 620);
    TPad p_comb("p_comb", "", 0.00, 0.00, 1.00, 0.93);
    p_comb.Draw();
    p_comb.cd();
    DrawEffPanel(eff_c, "#mu^{+} and #mu^{-}");
    c_comb.cd();
    DrawHeadline(headline, 0.035);
    const std::string png_comb = out_dir + "step1_eff_2d_pt_vs_q_eta_charge_comb.png";
    c_comb.SaveAs(png_comb.c_str());

    // Numbers behind the figures, so the maps can be checked against values rather than eyeballed.
    auto report = [](const char* tag, const TH2D* num, const TH2D* den, const TH2D* eff) {
        const int ny = den->GetNbinsY();
        const int b8 = den->GetYaxis()->FindBin(8.0 + 1e-6);   // first bin above the 8 GeV edge
        const double N = num->Integral(1, den->GetNbinsX(), b8, ny);
        const double D = den->Integral(1, den->GetNbinsX(), b8, ny);
        int filled = 0, undef = 0;
        for (int ix = 1; ix <= eff->GetNbinsX(); ++ix)
            for (int iy = 1; iy <= eff->GetNbinsY(); ++iy)
                (eff->GetBinContent(ix, iy) < 0. ? undef : filled)++;
        std::cout << "    " << tag << ": plateau (p_{T} > 8 GeV) eff = "
                  << (D > 0 ? N / D : 0.) << "  (denominator weight " << D << "),  "
                  << filled << " drawn cells, " << undef << " undefined cells\n";
    };
    std::cout << "\n===== MC single-muon 2D mu4 efficiency (" << sample << ", " << wp_text
              << ") =====\n";
    report("mu+ ", num_p, den_p, eff_p);
    report("mu- ", num_m, den_m, eff_m);
    report("comb", num_c, den_c, eff_c);
    std::cout << "  wrote " << png_sepr << "\n  wrote " << png_comb << std::endl;

    fin->Close();
}
