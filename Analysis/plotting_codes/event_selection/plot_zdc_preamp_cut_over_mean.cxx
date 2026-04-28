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

static std::vector<std::string> FilesForYear(int yr) {
    std::string base = kBase + "pbpb_20" + std::to_string(yr) + "/";
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
    return {};
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
    TGraph* g_cut1 = (TGraph*)((TGraph*)fc->Get(PbPbEvSelKey::kZDCFCalCut))->Clone();
    double cut2_ns = ((TParameter<double>*)fc->Get(PbPbEvSelKey::kZDCTimeCutNs))->GetVal();
    double cut_A   = ((TParameter<double>*)fc->Get(PbPbEvSelKey::kPreampACutADC))->GetVal();
    double cut_C   = ((TParameter<double>*)fc->Get(PbPbEvSelKey::kPreampCCutADC))->GetVal();
    fc->Close();

    TChain chain("HeavyIonD3PD", "HeavyIonD3PD");
    for (const auto& f : FilesForYear(yr)) {
        if (gSystem->AccessPathName(f.c_str()))
            { std::cerr << "Skipping: " << f << std::endl; continue; }
        chain.Add(f.c_str());
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

    YearResult res;
    res.cut_A  = cut_A;  res.cut_C  = cut_C;
    res.mean_A = (n_sel > 0) ? sumA / n_sel : 0.;
    res.mean_C = (n_sel > 0) ? sumC / n_sel : 0.;
    std::cout << "  mean_A=" << res.mean_A << "  mean_C=" << res.mean_C
              << "  cut_A=" << cut_A       << "  cut_C=" << cut_C
              << "  (cut-mean)_A=" << cut_A - res.mean_A
              << "  (cut-mean)_C=" << cut_C - res.mean_C << "\n";
    return res;
}

void plot_zdc_preamp_cut_over_mean() {
    gStyle->SetOptStat(0);
    gSystem->mkdir(kPlotDir.c_str(), true);

    const int years[3] = {23, 24, 25};
    YearResult R[3];
    for (int i = 0; i < 3; ++i) R[i] = ProcessYear(years[i]);

    double diffA[3], diffC[3];
    for (int i = 0; i < 3; ++i) {
        diffA[i] = R[i].cut_A - R[i].mean_A;
        diffC[i] = R[i].cut_C - R[i].mean_C;
    }

    // y-axis range: span both A and C, always include zero
    double lo = *std::min_element(diffA, diffA+3);
    lo = std::min(lo, *std::min_element(diffC, diffC+3));
    double hi = *std::max_element(diffA, diffA+3);
    hi = std::max(hi, *std::max_element(diffC, diffC+3));
    double pad = 0.22 * (hi - lo);
    lo = std::min(lo - pad - 40., 0.);
    hi = hi + pad + 40.;

    TCanvas* c = new TCanvas("c", "ZDC preamp threshold above baseline", 580, 520);
    c->SetLeftMargin(0.18);
    c->SetRightMargin(0.05);
    c->SetBottomMargin(0.18);
    c->SetTopMargin(0.07);

    const double xs[3] = {0., 1., 2.};
    TGraph* gA = new TGraph(3, xs, diffA);
    TGraph* gC = new TGraph(3, xs, diffC);

    gA->SetMarkerStyle(20); gA->SetMarkerSize(1.8);
    gA->SetMarkerColor(kBlack); gA->SetLineColor(kBlack); gA->SetLineWidth(2);
    gC->SetMarkerStyle(21); gC->SetMarkerSize(1.8);
    gC->SetMarkerColor(kBlue+1); gC->SetLineColor(kBlue+1); gC->SetLineWidth(2);

    TGraph* fr = new TGraph(3, xs, diffA);
    fr->SetTitle(";;Cut #minus mean [ADC]");
    fr->GetXaxis()->SetLimits(-0.5, 2.5);
    fr->GetYaxis()->SetRangeUser(lo, hi);
    fr->GetXaxis()->SetNdivisions(3);
    fr->SetMarkerStyle(1); fr->SetMarkerColor(0); fr->SetLineColor(0);
    fr->Draw("AP");
    fr->GetXaxis()->SetLimits(-0.5, 2.5);
    fr->GetYaxis()->SetRangeUser(lo, hi);
    fr->GetXaxis()->SetLabelOffset(999);
    fr->GetYaxis()->SetTitleOffset(1.5);
    fr->GetYaxis()->SetTitleSize(0.052);
    fr->GetXaxis()->SetTitleSize(0.052);

    if (lo < 0. && hi > 0.) {
        TLine* zl = new TLine(-0.5, 0., 2.5, 0.);
        zl->SetLineStyle(2); zl->SetLineColor(kGray+1);
        zl->Draw();
    }

    gA->Draw("LP same");
    gC->Draw("LP same");

    TLatex lab;
    lab.SetTextAlign(22); lab.SetTextSize(0.055);
    const double laby = lo - 0.09*(hi - lo);
    lab.DrawLatex(0., laby, "2023");
    lab.DrawLatex(1., laby, "2024");
    lab.DrawLatex(2., laby, "2025");

    TLatex val;
    val.SetTextSize(0.042); val.SetTextAlign(21);
    const double off = 0.04 * (hi - lo);
    for (int i = 0; i < 3; ++i) {
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
