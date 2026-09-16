// plot_zdc_preamp_cut_over_mean.cxx
//
// Single-panel plot: (hard cut − mean preamp sum) per year and side.
// Baseline-subtracts the preamp cut to expose whether the threshold above
// the ADC pedestal is stable across run years.
//
// Method: apply trigger + Cut1 (ZDC-FCal banana) + Cut2 (ZDC time), then
// compute the event-by-event preamp sum for side A ([1][0..3]) and C ([0][0..3]).
// Mean of that distribution is subtracted from each year's hard cut.
//
// Usage (run from Analysis/):
//   source /usatlas/u/yuhanguo/setup.sh
//   root -l -b -q 'plotting_codes/event_selection/plot_zdc_preamp_cut_over_mean.cxx+'

#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include "TChain.h"
#include "TFile.h"
#include "TParameter.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TAxis.h"
#include "../../NTupleProcessingCode/PbPbEventSelConfig.h"

static const std::string kBase    = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/";
static const std::string kPlotDir = kBase + "plots/single_b_analysis/event_selection/";

// The PbPb years this macro KNOWS ABOUT (configured superset), in x-axis order.
// It is NOT the list that gets drawn: AvailableYears() filters it at run time down
// to the years whose cuts file and skim are actually on disk, and every array, loop,
// axis limit and label below is sized from that FILTERED list.
static const int kYearsCfgPP[] = {23, 24, 25, 26};
static const int kNYearsCfg    = (int)(sizeof(kYearsCfgPP) / sizeof(kYearsCfgPP[0]));

static std::vector<std::string> FilesForYear(int yr) {
    std::string base = kBase + "pbpb_20" + std::to_string(yr) + "/";
    if (yr == 23) return { base+"data_pbpb23_part1.root",
                           base+"data_pbpb23_part2.root",
                           base+"data_pbpb23_part3.root",
                           base+"data_pbpb23_part4.root",
                           base+"data_pbpb23_part5.root" };
    if (yr == 24) return { base+"data_pbpb24_part1.root",
                           base+"data_pbpb24_part2.root" };
    if (yr == 25) return { base+"data_pbpb25_part1.root",
                           base+"data_pbpb25_part2.root",
                           base+"data_pbpb25_part3.root",
                           base+"data_pbpb25_part4.root",
                           base+"data_pbpb25_part5.root",
                           base+"data_pbpb25_part6.root",
                           base+"data_pbpb25_part7.root" };
    // Part lists = merged files on disk (2023 part5 / 2025 part7 = Sep-2026 recovery skims;
    // 2026 = 7 parts).  Keep equal to file_batch_max in PbPbExtras.c;
    // pipelines/preflight_pbpb_year.sh <yr> cross-checks.  An UNDER-count silently drops data.
    if (yr == 26) return { base+"data_pbpb26_part1.root",
                           base+"data_pbpb26_part2.root",
                           base+"data_pbpb26_part3.root",
                           base+"data_pbpb26_part4.root",
                           base+"data_pbpb26_part5.root",
                           base+"data_pbpb26_part6.root",
                           base+"data_pbpb26_part7.root" };
    // Reachable only if a year is added to kYearsCfgPP without being added here.
    // No silent empty list: it would build an empty TChain and give a mean of 0,
    // i.e. a plotted (cut - mean) that is pure fiction.
    throw std::runtime_error("FilesForYear: no input files configured for PbPb year "
                             + std::to_string(yr));
}

// The years that can actually be drawn: cuts file present AND at least one skim part
// present.  A configured-but-not-yet-produced year (2026 while the grid skim runs) is
// SKIPPED with an [INFO] line instead of aborting the whole figure.  Everything
// downstream is sized from this list, so the figure can never show N-1 years under an
// N-year legend/axis.
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

struct YearResult {
    double mean_A{0.}, mean_C{0.};
    double cut_A{0.},  cut_C{0.};
};

static YearResult ProcessYear(int yr) {
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
    TGraph* g_cut1 = loadG(PbPbEvSelKey::kZDCFCalCut);
    double cut2_ns = loadP(PbPbEvSelKey::kZDCTimeCutNs);
    double cut_A   = loadP(PbPbEvSelKey::kPreampACutADC);
    double cut_C   = loadP(PbPbEvSelKey::kPreampCCutADC);
    fc->Close();

    TChain chain("HeavyIonD3PD", "HeavyIonD3PD");
    for (const auto& f : FilesForYear(yr)) {
        if (gSystem->AccessPathName(f.c_str()))
            { std::cerr << "Skipping: " << f << std::endl; continue; }
        chain.Add(f.c_str());
    }
    // A present-but-empty (or entirely absent) skim would give n_sel = 0 and hence a
    // mean of 0 below — a (cut - 0) plotted as if it were measured.  Fail loudly.
    if (chain.GetEntries() == 0) {
        delete g_cut1;
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

    double sumA = 0., sumC = 0.;
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
        sumA += pA; sumC += pC;
        ++n_sel;
    }
    delete g_cut1;
    std::cout << "  " << n_sel << " pass trigger+Cut1+Cut2\n";
    if (n_sel == 0)
        throw std::runtime_error("No event passes trigger+Cut1+Cut2 for PbPb 20"
                                 + std::to_string(yr) + " — refusing to plot (cut - mean) "
                                 "against a placeholder mean.");

    YearResult res;
    res.cut_A  = cut_A;  res.cut_C  = cut_C;
    res.mean_A = sumA / n_sel;
    res.mean_C = sumC / n_sel;
    std::cout << "  mean_A=" << res.mean_A << "  mean_C=" << res.mean_C
              << "  cut_A=" << cut_A       << "  cut_C=" << cut_C
              << "  (cut-mean)_A=" << cut_A - res.mean_A
              << "  (cut-mean)_C=" << cut_C - res.mean_C << "\n";
    return res;
}

void plot_zdc_preamp_cut_over_mean() {
    gStyle->SetOptStat(0);
    gSystem->mkdir(kPlotDir.c_str(), true);

    const std::vector<int> years = AvailableYears();
    const int nY = (int)years.size();

    std::vector<YearResult> R(nY);
    for (int i = 0; i < nY; ++i) R[i] = ProcessYear(years[i]);

    std::vector<double> diffA(nY), diffC(nY);
    for (int i = 0; i < nY; ++i) {
        diffA[i] = R[i].cut_A - R[i].mean_A;
        diffC[i] = R[i].cut_C - R[i].mean_C;
    }

    // y-axis range: span both A and C, always include zero
    double lo = *std::min_element(diffA.begin(), diffA.end());
    lo = std::min(lo, *std::min_element(diffC.begin(), diffC.end()));
    double hi = *std::max_element(diffA.begin(), diffA.end());
    hi = std::max(hi, *std::max_element(diffC.begin(), diffC.end()));
    double pad = 0.22 * (hi - lo);
    lo = std::min(lo - pad - 40., 0.);
    hi = hi + pad + 40.;

    TCanvas* c = new TCanvas("c", "ZDC preamp threshold above baseline", 580, 520);
    c->SetLeftMargin(0.18);
    c->SetRightMargin(0.05);
    c->SetBottomMargin(0.18);
    c->SetTopMargin(0.07);

    std::vector<double> xs(nY);
    for (int i = 0; i < nY; ++i) xs[i] = (double)i;
    const double x_lo = -0.5, x_hi = nY - 0.5;
    TGraph* gA = new TGraph(nY, xs.data(), diffA.data());
    TGraph* gC = new TGraph(nY, xs.data(), diffC.data());

    gA->SetMarkerStyle(20); gA->SetMarkerSize(1.8);
    gA->SetMarkerColor(kBlack); gA->SetLineColor(kBlack); gA->SetLineWidth(2);
    gC->SetMarkerStyle(21); gC->SetMarkerSize(1.8);
    gC->SetMarkerColor(kBlue+1); gC->SetLineColor(kBlue+1); gC->SetLineWidth(2);

    TGraph* fr = new TGraph(nY, xs.data(), diffA.data());
    fr->SetTitle(";;Cut #minus mean [ADC]");
    fr->GetXaxis()->SetLimits(x_lo, x_hi);
    fr->GetYaxis()->SetRangeUser(lo, hi);
    fr->GetXaxis()->SetNdivisions(nY);
    fr->SetMarkerStyle(1); fr->SetMarkerColor(0); fr->SetLineColor(0);
    fr->Draw("AP");
    fr->GetXaxis()->SetLimits(x_lo, x_hi);
    fr->GetYaxis()->SetRangeUser(lo, hi);
    fr->GetXaxis()->SetLabelOffset(999);
    fr->GetYaxis()->SetTitleOffset(1.5);
    fr->GetYaxis()->SetTitleSize(0.052);
    fr->GetXaxis()->SetTitleSize(0.052);

    if (lo < 0. && hi > 0.) {
        TLine* zl = new TLine(x_lo, 0., x_hi, 0.);
        zl->SetLineStyle(2); zl->SetLineColor(kGray+1);
        zl->Draw();
    }

    gA->Draw("LP same");
    gC->Draw("LP same");

    TLatex lab;
    lab.SetTextAlign(22); lab.SetTextSize(0.055);
    const double laby = lo - 0.09*(hi - lo);
    for (int i = 0; i < nY; ++i)
        lab.DrawLatex(xs[i], laby, Form("20%d", years[i]));

    TLatex val;
    val.SetTextSize(0.042); val.SetTextAlign(21);
    const double off = 0.04 * (hi - lo);
    for (int i = 0; i < nY; ++i) {
        val.SetTextColor(kBlack);
        val.DrawLatex(xs[i] - 0.09, diffA[i] + off, Form("%.0f", diffA[i]));
        val.SetTextColor(kBlue+1);
        val.DrawLatex(xs[i] + 0.09, diffC[i] + off, Form("%.0f", diffC[i]));
    }

    TLegend* leg = new TLegend(0.22, 0.73, 0.54, 0.91);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.050);
    leg->AddEntry(gA, "Side A", "lp");
    leg->AddEntry(gC, "Side C", "lp");
    leg->Draw();

    const std::string out = kPlotDir + "zdc_preamp_cut_minus_mean.png";
    c->SaveAs(out.c_str());
    std::cout << "Saved: " << out << std::endl;
}
