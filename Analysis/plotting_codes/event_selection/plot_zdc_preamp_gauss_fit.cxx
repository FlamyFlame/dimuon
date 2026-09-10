// plot_zdc_preamp_gauss_fit.cxx
//
// 2×N canvas of ZDC presample amplitude distributions with Gaussian fits.
//   Rows: Side A (top), Side C (bottom)
//   Columns: one per AVAILABLE PbPb year (subset of kYearsCfgPP = 2023, 2024, 2025, 2026)
//
// Data: events passing trigger + Cut1 (ZDC-FCal banana) + Cut2 (ZDC time).
// Fit range matches the per-year zoom window used in the standalone cut-3 plot.
// Each subplot legend reports the fitted Gaussian μ and σ.
//
// Usage (run from Analysis/):
//   source /usatlas/u/yuhanguo/setup.sh
//   root -l -b -q 'plotting_codes/event_selection/plot_zdc_preamp_gauss_fit.cxx+'

#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include "TChain.h"
#include "TFile.h"
#include "TParameter.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TSystem.h"
#include "TStyle.h"
#include "../../NTupleProcessingCode/PbPbEventSelConfig.h"

static const std::string kBase    = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/";
static const std::string kPlotDir = kBase + "plots/single_b_analysis/event_selection/";

// The PbPb years this macro KNOWS ABOUT (configured superset), in column order.
// It is NOT the list that gets drawn: AvailableYears() filters it at run time down to
// the years whose cuts file and skim are actually on disk, and the column count, pad
// numbering, canvas width and output filename are all sized from that FILTERED list.
static const int kYearsCfgPP[] = {23, 24, 25, 26};
static const int kNYearsCfg    = (int)(sizeof(kYearsCfgPP) / sizeof(kYearsCfgPP[0]));

static std::vector<std::string> FilesForYear(int yr) {
    std::string base = kBase + "pbpb_20" + std::to_string(yr) + "/";
    // NOTE: 2023 lists 3 parts here while plot_pbpb_event_sel_cuts.cxx lists 4.
    // PRE-EXISTING inconsistency, left as found — do not "fix" it silently.
    if (yr == 23) return { base+"data_pbpb23_part1.root",
                           base+"data_pbpb23_part2.root",
                           base+"data_pbpb23_part3.root" };
    if (yr == 24) return { base+"data_pbpb24_part1.root",
                           base+"data_pbpb24_part2.root" };
    if (yr == 25) return { base+"data_pbpb25_part1.root",
                           base+"data_pbpb25_part2.root",
                           base+"data_pbpb25_part3.root",
                           base+"data_pbpb25_part4.root",
                           base+"data_pbpb25_part5.root",
                           base+"data_pbpb25_part6.root" };
    // PLACEHOLDER (2026): 5 grid tasks were submitted
    // (SkimCode/run_26hi/InDstxt_PbPb2026_5p36TeV_part1..5.txt) so there are AT
    // LEAST 5 parts; grid_monitor's chunked-hadd fallback can produce more.
    // CONFIRM against pbpb_2026/ and extend — an UNDER-count silently drops data.
    // See docs/tracking/pbpb2026_analysis_support.md.
    if (yr == 26) return { base+"data_pbpb26_part1.root",
                           base+"data_pbpb26_part2.root",
                           base+"data_pbpb26_part3.root",
                           base+"data_pbpb26_part4.root",
                           base+"data_pbpb26_part5.root" };
    // Reachable only if a year is added to kYearsCfgPP without being added here.
    // No silent empty list: it would build an empty TChain and produce an empty
    // histogram that still gets fitted and drawn as if it were data.
    throw std::runtime_error("FilesForYear: no input files configured for PbPb year "
                             + std::to_string(yr));
}

// The years that can actually be drawn: cuts file present AND at least one skim part
// present.  A configured-but-not-yet-produced year (2026 while the grid skim runs) is
// SKIPPED with an [INFO] line instead of aborting the whole figure.  The canvas grid
// and the output filename are sized from this list, so a 2×N canvas can never carry a
// column for a year that contributed nothing.
static std::vector<int> AvailableYears() {
    std::vector<int> years;
    for (int i = 0; i < kNYearsCfg; ++i) {
        const int yr = kYearsCfgPP[i];
        const std::string cpath = PbPbEvSelCutsPath(2000 + yr);
        if (gSystem->AccessPathName(cpath.c_str())) {
            std::cout << "[INFO] PbPb 20" << yr << ": skipping — cuts file not found: "
                      << cpath << std::endl;
            continue;
        }
        bool any_part = false;
        for (const auto& f : FilesForYear(yr))
            if (!gSystem->AccessPathName(f.c_str())) { any_part = true; break; }
        if (!any_part) {
            std::cout << "[INFO] PbPb 20" << yr << ": skipping — no skim part file found under "
                      << kBase << "pbpb_20" << yr << "/" << std::endl;
            continue;
        }
        std::cout << "[INFO] PbPb 20" << yr << ": available." << std::endl;
        years.push_back(yr);
    }
    if (years.empty())
        throw std::runtime_error("AvailableYears: no PbPb year has both a cuts file and a skim "
                                 "input on disk — nothing to plot.");
    return years;
}

// Per-year zoom range (mirrors DrawStandaloneCut3 in plot_pbpb_event_sel_cuts.cxx)
static void ZoomRange(int yr, double& lo, double& hi) {
    if      (yr == 23) { lo = -500.; hi =  700.; }
    else if (yr == 24) { lo = -800.; hi =  800.; }
    else if (yr == 25) { lo = -800.; hi = 1500.; }
    // PLACEHOLDER (2026): copied from 2025; the fit range is also the Gaussian
    // seed window, so re-check it against the 2026 preamp distribution.
    // See docs/tracking/pbpb2026_analysis_support.md.
    else if (yr == 26) { lo = -800.; hi = 1500.; }
    // No silent default window: the zoom range IS the fit range here, so an
    // unlisted year would be fitted over a window nobody chose.
    else throw std::runtime_error("ZoomRange: no preamp zoom/fit range defined for PbPb year "
                                  + std::to_string(yr));
}

struct YearHists {
    TH1D* hA{nullptr};
    TH1D* hC{nullptr};
    double cut_A{0.}, cut_C{0.};
};

static YearHists ProcessYear(int yr) {
    const std::string cpath = PbPbEvSelCutsPath(2000 + yr);
    TFile* fc = TFile::Open(cpath.c_str(), "READ");
    if (!fc || fc->IsZombie())
        throw std::runtime_error("Cuts file not found: " + cpath);
    // Throwing accessors: an unchecked fc->Get() segfaults on a file that opens but
    // is missing a key (e.g. an interrupted cut-derivation run).
    auto loadG = [&](const char* key) -> TGraph* {
        TGraph* g = (TGraph*)fc->Get(key);
        if (!g) { fc->Close(); throw std::runtime_error(
            std::string("Missing TGraph '") + key + "' in " + cpath); }
        return (TGraph*)g->Clone();
    };
    auto loadP = [&](const char* key) -> double {
        TParameter<double>* p = (TParameter<double>*)fc->Get(key);
        if (!p) { fc->Close(); throw std::runtime_error(
            std::string("Missing TParameter '") + key + "' in " + cpath); }
        return p->GetVal();
    };
    TGraph* g_cut1  = loadG(PbPbEvSelKey::kZDCFCalCut);
    double cut2_ns  = loadP(PbPbEvSelKey::kZDCTimeCutNs);
    double cut_A    = loadP(PbPbEvSelKey::kPreampACutADC);
    double cut_C    = loadP(PbPbEvSelKey::kPreampCCutADC);
    fc->Close();

    // 198 bins, −1000 to 2960 → ~20 ADC/bin; rebin×3 for display → 60 ADC/bin
    TH1D* hA = new TH1D(Form("hA_%d", yr), ";ZDC preamp sum [ADC];Events / 60 ADC", 198, -1000., 2960.);
    TH1D* hC = new TH1D(Form("hC_%d", yr), ";ZDC preamp sum [ADC];Events / 60 ADC", 198, -1000., 2960.);
    hA->SetDirectory(nullptr);
    hC->SetDirectory(nullptr);

    TChain chain("HeavyIonD3PD", "HeavyIonD3PD");
    for (const auto& f : FilesForYear(yr)) {
        if (gSystem->AccessPathName(f.c_str()))
            { std::cerr << "Skipping: " << f << std::endl; continue; }
        chain.Add(f.c_str());
    }
    // An empty chain would leave hA/hC empty, and the empty histogram would still be
    // fitted and drawn as if it were data.  Fail loudly instead.
    if (chain.GetEntries() == 0) {
        delete g_cut1; delete hA; delete hC;
        throw std::runtime_error("TChain empty for PbPb 20" + std::to_string(yr));
    }
    chain.SetMakeClass(1);
    chain.SetBranchStatus("*", 0);
    Int_t   b_HLT = 0;
    Float_t FCal_P = 0.f, FCal_N = 0.f;
    Float_t zdcE[2]{}, zdcT[2]{}, preamp[2][4]{};
    chain.SetBranchStatus("b_HLT_mu4_L1MU3V",         1);
    chain.SetBranchStatus("FCal_Et_P",                 1);
    chain.SetBranchStatus("FCal_Et_N",                 1);
    chain.SetBranchStatus("zdc_ZdcEnergy",             1);
    chain.SetBranchStatus("zdc_ZdcTime",               1);
    chain.SetBranchStatus("zdc_ZdcModulePreSampleAmp", 1);
    chain.SetBranchAddress("b_HLT_mu4_L1MU3V",         &b_HLT);
    chain.SetBranchAddress("FCal_Et_P",                 &FCal_P);
    chain.SetBranchAddress("FCal_Et_N",                 &FCal_N);
    chain.SetBranchAddress("zdc_ZdcEnergy",             zdcE);
    chain.SetBranchAddress("zdc_ZdcTime",               zdcT);
    chain.SetBranchAddress("zdc_ZdcModulePreSampleAmp", preamp);

    long long n_sel = 0;
    const Long64_t n = chain.GetEntries();
    std::cout << "20" << yr << ": " << n << " events..." << std::flush;
    for (Long64_t i = 0; i < n; ++i) {
        chain.GetEntry(i);
        if (!b_HLT) continue;
        const float fcal_AC = (FCal_P + FCal_N) * 1e-6f;
        const float zdcTot  = (zdcE[0] + zdcE[1]) / 1000.f;
        if (zdcTot > (float)PbPbEvSelEvalCut(g_cut1, fcal_AC)) continue; // Cut1
        if (std::abs(zdcT[1]) >= (float)cut2_ns) continue;               // Cut2 A
        if (std::abs(zdcT[0]) >= (float)cut2_ns) continue;               // Cut2 C
        float pA = 0.f, pC = 0.f;
        for (int k = 0; k < 4; ++k) pA += preamp[1][k]; // [1]=A
        for (int k = 0; k < 4; ++k) pC += preamp[0][k]; // [0]=C
        hA->Fill(pA);
        hC->Fill(pC);
        ++n_sel;
    }
    delete g_cut1;
    std::cout << "  " << n_sel << " pass trigger+Cut1+Cut2\n";
    if (n_sel == 0) {
        delete hA; delete hC;
        throw std::runtime_error("No event passes trigger+Cut1+Cut2 for PbPb 20"
                                 + std::to_string(yr) + " — refusing to fit an empty histogram.");
    }

    YearHists res;
    res.hA    = hA;    res.hC    = hC;
    res.cut_A = cut_A; res.cut_C = cut_C;
    return res;
}

// Fit Gaussian to histogram peak. Returns (mu, sigma) from final fit.
static std::pair<double,double> FitGauss(TH1D* h, double zoom_lo, double zoom_hi,
                                          const char* fname) {
    // Pass 1: rough Gaussian seed over the zoom range (left 60% to avoid tail)
    double seed_hi = zoom_lo + 0.60 * (zoom_hi - zoom_lo);
    TF1 g1(Form("%s_p1", fname), "gaus", zoom_lo, seed_hi);
    double amp0 = h->GetMaximum();
    g1.SetParameters(amp0, (zoom_lo + seed_hi) * 0.5, (seed_hi - zoom_lo) * 0.25);
    h->Fit(&g1, "RQN");
    double mu1  = g1.GetParameter(1);
    double sig1 = std::abs(g1.GetParameter(2));

    // Pass 2: fit Gaussian in [mu-3σ, mu+1.5σ] to focus on peak core
    double fit_lo = mu1 - 3. * sig1;
    double fit_hi = mu1 + 1.5 * sig1;
    TF1* gf = new TF1(fname, "gaus", fit_lo, fit_hi);
    gf->SetParameters(g1.GetParameter(0), mu1, sig1);
    h->Fit(gf, "RQN");
    double mu_f  = gf->GetParameter(1);
    double sig_f = std::abs(gf->GetParameter(2));
    delete gf;
    return {mu_f, sig_f};
}

void plot_zdc_preamp_gauss_fit() {
    gStyle->SetOptStat(0);
    gSystem->mkdir(kPlotDir.c_str(), true);

    const std::vector<int> years = AvailableYears();
    const int nY = (int)years.size();

    std::vector<YearHists> Y(nY);
    for (int i = 0; i < nY; ++i) Y[i] = ProcessYear(years[i]);

    // Rebin histograms for display (×3 → 60 ADC/bin)
    const int rebin = 3;
    for (int i = 0; i < nY; ++i) {
        Y[i].hA->Rebin(rebin); Y[i].hA->Scale(1.0 / rebin);
        Y[i].hC->Rebin(rebin); Y[i].hC->Scale(1.0 / rebin);
    }

    // Canvas: nY columns × 2 rows; pad index = col + row*nY + 1 (ROOT numbering).
    // Width scales with the column count actually drawn.
    TCanvas* cv = new TCanvas("cv_gauss", "", 500 * nY, 900);
    cv->Divide(nY, 2, 0.003, 0.003);

    const char* sideLabel[2] = {"A", "C"};

    for (int row = 0; row < 2; ++row) {      // row 0 = Side A, row 1 = Side C
        for (int col = 0; col < nY; ++col) { // col = available-year index
            int yr    = years[col];
            int padNo = col + 1 + row * nY;  // ROOT pad number (1-based)
            cv->cd(padNo);
            gPad->SetLogy();
            gPad->SetLeftMargin(0.15);
            gPad->SetRightMargin(0.04);
            gPad->SetTopMargin(0.08);
            gPad->SetBottomMargin(0.14);

            TH1D* h    = (row == 0) ? Y[col].hA : Y[col].hC;
            double cutV = (row == 0) ? Y[col].cut_A : Y[col].cut_C;

            double zoom_lo, zoom_hi;
            ZoomRange(yr, zoom_lo, zoom_hi);

            // Gaussian fit
            TString fname = Form("gfit_%s_%d", sideLabel[row], yr);
            auto [mu_f, sig_f] = FitGauss(h, zoom_lo, zoom_hi, fname.Data());

            // Attach fit TF1 for drawing (fit again to get the pointer)
            double draw_lo = mu_f - 3.*sig_f;
            double draw_hi = mu_f + 1.5*sig_f;
            TF1* gDraw = new TF1(Form("gdraw_%s_%d", sideLabel[row], yr),
                                 "gaus", draw_lo, draw_hi);
            gDraw->SetParameters(h->GetBinContent(h->FindBin(mu_f)), mu_f, sig_f);
            h->Fit(gDraw, "RQN");
            gDraw->SetLineColor(kRed+1); gDraw->SetLineWidth(2);

            // X range and Y range from local max within zoom
            h->GetXaxis()->SetRangeUser(zoom_lo, zoom_hi);
            int b_lo = h->FindBin(zoom_lo + 0.5);
            int b_hi = h->FindBin(zoom_hi - 0.5);
            double lmax = 0.;
            for (int b = b_lo; b <= b_hi; ++b)
                lmax = std::max(lmax, h->GetBinContent(b));
            h->GetYaxis()->SetRangeUser(0.5, lmax * 4.0);

            h->SetMarkerStyle(20); h->SetMarkerSize(0.5);
            h->SetMarkerColor(kBlack); h->SetLineColor(kBlack);
            h->GetXaxis()->SetTitle(Form("ZDC preamp sum side %s [ADC]", sideLabel[row]));
            h->GetXaxis()->SetTitleSize(0.050); h->GetXaxis()->SetLabelSize(0.043);
            h->GetYaxis()->SetTitleSize(0.047); h->GetYaxis()->SetLabelSize(0.040);
            h->GetYaxis()->SetTitleOffset(1.35);
            h->Draw("E");
            gDraw->Draw("same");

            // Cut line
            TLine lcut(cutV, 0.5, cutV, lmax * 4.0);
            lcut.SetLineColor(kBlue+1); lcut.SetLineWidth(2); lcut.SetLineStyle(2);
            lcut.DrawClone();

            // Legend
            TLegend* leg = new TLegend(0.48, 0.62, 0.95, 0.90);
            leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.042);
            leg->AddEntry(h,     "Data (after Cut1+2)", "lp");
            leg->AddEntry(gDraw, Form("Gauss fit  #mu=%.0f ADC", mu_f), "l");
            leg->AddEntry((TObject*)nullptr, Form("                     #sigma=%.0f ADC", sig_f), "");
            leg->Draw();

            // Panel label
            TLatex tl; tl.SetNDC(); tl.SetTextSize(0.048);
            tl.DrawLatex(0.18, 0.91, Form("Pb+Pb %d  Side %s", 2000+yr, sideLabel[row]));
            tl.SetTextSize(0.040); tl.SetTextColor(kBlue+1);
            tl.DrawLatex(0.18, 0.84, Form("cut = %.0f ADC", cutV));
            tl.SetTextColor(kBlack);
        }
    }

    // Filename carries the grid ACTUALLY drawn, so it can never claim a layout it does
    // not have: 3 available years -> zdc_preamp_gauss_fit_2x3.png (the file already on
    // disk), 4 -> zdc_preamp_gauss_fit_2x4.png.
    const std::string out = kPlotDir + Form("zdc_preamp_gauss_fit_2x%d.png", nY);
    cv->SaveAs(out.c_str());
    std::cout << "Saved: " << out << std::endl;

    for (int i = 0; i < nY; ++i) { delete Y[i].hA; delete Y[i].hC; }
}
