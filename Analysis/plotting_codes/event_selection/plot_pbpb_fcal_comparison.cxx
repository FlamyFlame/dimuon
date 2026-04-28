// plot_pbpb_fcal_comparison.cxx
// FCal ET shape comparison: PbPb 2024 and 2025 vs 2023 after full 5-cut
// event selection (nominal two-band banana OR alternative quadratic banana).
// Each output is a single canvas with two ratio sub-plots (24/23, 25/23).
//
// Usage:
//   .L plot_pbpb_fcal_comparison.cxx+
//   plot_pbpb_fcal_comparison()       // nominal 2-band banana cuts
//   plot_pbpb_fcal_comparison_alt()   // alternative quadratic banana cuts

#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <iostream>
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TParameter.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "../../NTupleProcessingCode/PbPbEventSelConfig.h"

// ---- Input file map ---------------------------------------------------------
static std::map<int, std::vector<std::string>> BuildFilesFC() {
    const std::string base = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/";
    return {
        {23, {base+"pbpb_2023/data_pbpb23_part1.root",
              base+"pbpb_2023/data_pbpb23_part2.root",
              base+"pbpb_2023/data_pbpb23_part3.root"}},
        {24, {base+"pbpb_2024/data_pbpb24_part1.root",
              base+"pbpb_2024/data_pbpb24_part2.root"}},
        {25, {base+"pbpb_2025/data_pbpb25_part1.root",
              base+"pbpb_2025/data_pbpb25_part2.root",
              base+"pbpb_2025/data_pbpb25_part3.root",
              base+"pbpb_2025/data_pbpb25_part4.root",
              base+"pbpb_2025/data_pbpb25_part5.root",
              base+"pbpb_2025/data_pbpb25_part6.root"}},
    };
}

// ---- Per-year cut container -------------------------------------------------
struct CutSetFC {
    TGraph* g_cut1    = nullptr;  // ZDC-FCal banana
    double  t_ns      = 1.5;
    double  preamp_A  = 385.;
    double  preamp_C  = 385.;
    TGraph* g_cut4    = nullptr;  // nTrk frac lower bound
    TGraph* g_cut5_lo = nullptr;  // nTrk-FCal lower bound
    TGraph* g_cut5_hi = nullptr;  // nTrk-FCal upper bound
};

// Load all cuts from the per-year ROOT file (nominal or alt).
static CutSetFC LoadCutsFC(int yr, bool is_alt) {
    const std::string ys  = std::to_string(yr);
    const std::string sfx = is_alt ? "_alt" : "";
    const std::string path = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_20" +
                             ys + "/event_sel_cuts_pbpb_20" + ys + sfx + ".root";
    TFile* f = TFile::Open(path.c_str(), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Cannot open cuts file: " << path << std::endl;
        return {};
    }
    CutSetFC cs;
    auto loadG = [&](const char* key) -> TGraph* {
        TGraph* g = (TGraph*)f->Get(key);
        return g ? (TGraph*)g->Clone() : nullptr;
    };
    auto loadP = [&](const char* key, double def) -> double {
        TParameter<double>* p = (TParameter<double>*)f->Get(key);
        return p ? p->GetVal() : def;
    };
    cs.g_cut1    = loadG(PbPbEvSelKey::kZDCFCalCut);
    cs.t_ns      = loadP(PbPbEvSelKey::kZDCTimeCutNs,  1.5);
    cs.preamp_A  = loadP(PbPbEvSelKey::kPreampACutADC, 385.);
    cs.preamp_C  = loadP(PbPbEvSelKey::kPreampCCutADC, 385.);
    cs.g_cut4    = loadG(PbPbEvSelKey::kNTrkFracCutLo);
    cs.g_cut5_lo = loadG(PbPbEvSelKey::kNTrkFCalCutLo);
    cs.g_cut5_hi = loadG(PbPbEvSelKey::kNTrkFCalCutHi);
    f->Close();
    std::cout << Form("Loaded cuts for 20%d (%s): cut1=%s cut4=%s cut5lo=%s\n",
                      yr, is_alt ? "alt" : "nominal",
                      cs.g_cut1    ? "ok" : "MISSING",
                      cs.g_cut4    ? "ok" : "MISSING",
                      cs.g_cut5_lo ? "ok" : "MISSING");
    return cs;
}

// ---- Event loop: fill FCal ET histogram after all 5 cuts --------------------
static TH1D* FillFCalHist(int yr, const CutSetFC& cs) {
    const std::string ys = std::to_string(yr);
    auto fm = BuildFilesFC();
    auto it = fm.find(yr);
    if (it == fm.end()) { std::cerr << "No files for year " << yr << std::endl; return nullptr; }

    static int s_uid = 0;
    TH1D* h = new TH1D(Form("h_fcal_yr%d_%d", yr, s_uid++),
                        ";FCal E_{T}^{A+C} [TeV];Events (normalized to unit area)",
                        120, -0.5, 5.5);
    h->SetDirectory(nullptr);
    h->Sumw2();

    TChain chain("HeavyIonD3PD", "HeavyIonD3PD");
    for (const auto& fpath : it->second)
        if (!gSystem->AccessPathName(fpath.c_str())) chain.Add(fpath.c_str());
    if (chain.GetEntries() == 0) {
        std::cerr << "Empty TChain for year " << yr << std::endl; return h;
    }

    chain.SetMakeClass(1);
    chain.SetBranchStatus("*", 0);

    Int_t   b_HLT = 0;
    Float_t FCal_Et_P = 0.f, FCal_Et_N = 0.f;
    Float_t zdc_E[2] = {}, zdc_t[2] = {};
    Float_t preamp[2][4] = {};
    std::vector<int>* trk_numqual = nullptr;

    for (const char* br : {"b_HLT_mu4_L1MU3V","FCal_Et_P","FCal_Et_N",
                            "zdc_ZdcEnergy","zdc_ZdcTime","zdc_ZdcModulePreSampleAmp","trk_numqual"})
        chain.SetBranchStatus(br, 1);

    chain.SetBranchAddress("b_HLT_mu4_L1MU3V",        &b_HLT);
    chain.SetBranchAddress("FCal_Et_P",                 &FCal_Et_P);
    chain.SetBranchAddress("FCal_Et_N",                 &FCal_Et_N);
    chain.SetBranchAddress("zdc_ZdcEnergy",              zdc_E);
    chain.SetBranchAddress("zdc_ZdcTime",                zdc_t);
    chain.SetBranchAddress("zdc_ZdcModulePreSampleAmp",  preamp);
    chain.SetBranchAddress("trk_numqual",               &trk_numqual);

    const Long64_t ntot = chain.GetEntries();
    std::cout << Form("PbPb 20%d: processing %lld events...\n", yr, ntot);
    Long64_t n_pass = 0;

    for (Long64_t i = 0; i < ntot; ++i) {
        chain.GetEntry(i);
        if (!b_HLT) continue;

        const float fcal_AC = (FCal_Et_P + FCal_Et_N) * 1e-6f;
        const float zdc_tot = (zdc_E[0] + zdc_E[1]) / 1000.f;
        const float tA      = zdc_t[1];   // [1]=A
        const float tC      = zdc_t[0];   // [0]=C
        const float pA      = preamp[1][0]+preamp[1][1]+preamp[1][2]+preamp[1][3];
        const float pC      = preamp[0][0]+preamp[0][1]+preamp[0][2]+preamp[0][3];

        int ntrk_total = 0, ntrk_tight = 0;
        if (trk_numqual && (int)trk_numqual->size() >= 4) {
            ntrk_total = (*trk_numqual)[0];
            ntrk_tight = (*trk_numqual)[3];
        }

        // Cut 1: ZDC-FCal banana
        const double c1 = PbPbEvSelEvalCut(cs.g_cut1, fcal_AC);
        if (c1 > 0. && zdc_tot > c1) continue;
        // Cut 2: ZDC time box
        if (std::fabs(tA) >= (float)cs.t_ns || std::fabs(tC) >= (float)cs.t_ns) continue;
        // Cut 3: ZDC preamp
        if (pA >= (float)cs.preamp_A || pC >= (float)cs.preamp_C) continue;
        // Cuts 4 & 5: track-based — skip events with no tracks (consistent with cuts plotter)
        if (ntrk_total <= 0) continue;
        // Cut 4: nTrk frac lower bound
        const double frac = (double)ntrk_tight / ntrk_total;
        if (frac < PbPbEvSelEvalCut(cs.g_cut4, ntrk_total)) continue;
        // Cut 5: nTrk-FCal band
        if ((double)ntrk_tight < PbPbEvSelEvalCut(cs.g_cut5_lo, fcal_AC) ||
            (double)ntrk_tight > PbPbEvSelEvalCut(cs.g_cut5_hi, fcal_AC)) continue;

        h->Fill(fcal_AC);
        ++n_pass;
    }
    std::cout << Form("PbPb 20%d: %lld events passed all cuts\n", yr, n_pass);
    return h;
}

// ---- Ratio canvas -----------------------------------------------------------
static void DrawRatioColumn(TPad* p_top, TPad* p_bot,
                            TH1D* h_ref, TH1D* h_yr,
                            int col_ref, int col_yr,
                            const char* lbl_ref, const char* lbl_yr,
                            const char* proc_label)
{
    const double split   = p_bot->GetAbsHNDC() / (p_top->GetAbsHNDC() + p_bot->GetAbsHNDC());
    const double scale   = (1. - split) / split;  // for scaling axis labels in ratio pad

    // ----- top pad -----
    p_top->cd();
    p_top->SetBottomMargin(0.015);
    p_top->SetTopMargin(0.12);
    p_top->SetLeftMargin(0.18);
    p_top->SetRightMargin(0.04);

    TH1D* ht_ref = (TH1D*)h_ref->Clone();
    TH1D* ht_yr  = (TH1D*)h_yr ->Clone();
    ht_ref->SetDirectory(nullptr);
    ht_yr ->SetDirectory(nullptr);

    ht_ref->SetLineColor(col_ref);  ht_ref->SetMarkerColor(col_ref);
    ht_yr ->SetLineColor(col_yr);   ht_yr ->SetMarkerColor(col_yr);
    ht_ref->SetMarkerStyle(20);  ht_ref->SetMarkerSize(0.4);
    ht_yr ->SetMarkerStyle(24);  ht_yr ->SetMarkerSize(0.4);

    ht_ref->GetXaxis()->SetLabelSize(0.);
    ht_ref->GetYaxis()->SetTitle("Normalized");
    ht_ref->GetYaxis()->SetTitleSize(0.060);
    ht_ref->GetYaxis()->SetLabelSize(0.052);
    ht_ref->GetYaxis()->SetTitleOffset(1.30);

    const double ymax = std::max(ht_ref->GetMaximum(), ht_yr->GetMaximum()) * 1.20;
    ht_ref->GetYaxis()->SetRangeUser(0., ymax);
    ht_ref->Draw("E HIST");
    ht_yr ->Draw("E HIST SAME");

    TLegend* leg = new TLegend(0.55, 0.67, 0.93, 0.85);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.052);
    leg->AddEntry(ht_ref, lbl_ref, "lp");
    leg->AddEntry(ht_yr,  lbl_yr,  "lp");
    leg->Draw();

    TLatex tl; tl.SetNDC(); tl.SetTextSize(0.052); tl.SetTextColor(kGray+2);
    tl.DrawLatex(0.20, 0.86, proc_label);

    // ----- bottom pad -----
    p_bot->cd();
    p_bot->SetTopMargin(0.00);
    p_bot->SetBottomMargin(0.42);
    p_bot->SetLeftMargin(0.18);
    p_bot->SetRightMargin(0.04);

    TH1D* ratio = (TH1D*)h_yr->Clone();
    ratio->SetDirectory(nullptr);
    ratio->Divide(h_ref);
    ratio->SetLineColor(col_yr);
    ratio->SetMarkerColor(col_yr);
    ratio->SetMarkerStyle(24);
    ratio->SetMarkerSize(0.4 * scale);

    ratio->GetXaxis()->SetTitle("FCal E_{T}^{A+C} [TeV]");
    ratio->GetXaxis()->SetTitleSize(0.060 * scale);
    ratio->GetXaxis()->SetLabelSize(0.052 * scale);
    ratio->GetXaxis()->SetTitleOffset(1.0);
    ratio->GetYaxis()->SetTitle(Form("%s / %s", lbl_yr, lbl_ref));
    ratio->GetYaxis()->SetTitleSize(0.055 * scale);
    ratio->GetYaxis()->SetLabelSize(0.048 * scale);
    ratio->GetYaxis()->SetTitleOffset(1.30 / scale);
    ratio->GetYaxis()->SetNdivisions(505);
    ratio->GetYaxis()->SetRangeUser(0.50, 1.50);

    ratio->Draw("E");
    TLine line(h_ref->GetXaxis()->GetXmin(), 1., h_ref->GetXaxis()->GetXmax(), 1.);
    line.SetLineColor(kGray+1); line.SetLineStyle(2); line.SetLineWidth(1);
    line.DrawClone();
}

// ---- Main comparison function -----------------------------------------------
static void MakeComparisonPlot(bool is_alt) {
    CutSetFC cs23 = LoadCutsFC(23, is_alt);
    CutSetFC cs24 = LoadCutsFC(24, is_alt);
    CutSetFC cs25 = LoadCutsFC(25, is_alt);

    TH1D* h23 = FillFCalHist(23, cs23);
    TH1D* h24 = FillFCalHist(24, cs24);
    TH1D* h25 = FillFCalHist(25, cs25);

    if (!h23 || !h24 || !h25) {
        std::cerr << "Missing histogram — aborting." << std::endl; return;
    }
    // Normalize to unit area within 0-80% centrality (FCal_ET > PbPb2023 80% boundary).
    static const double kFCal80pct = 0.063208;  // TeV — FCal_ET_Bins_PbPb2023[79]
    for (TH1D* h : {h23, h24, h25}) {
        const int bin_lo = h->FindBin(kFCal80pct);
        const double integ = h->Integral(bin_lo, h->GetNbinsX());
        if (integ > 0.) h->Scale(1. / integ);
    }

    const double split = 0.30;   // fraction of canvas height for ratio pads

    TCanvas* c = new TCanvas("c_fcal_comp", "", 1400, 750);
    c->SetMargin(0, 0, 0, 0);

    // Four pads: left top, left bottom, right top, right bottom
    TPad* pLt = new TPad("pLt", "", 0.00, split, 0.50, 1.00);
    TPad* pLb = new TPad("pLb", "", 0.00, 0.00, 0.50, split);
    TPad* pRt = new TPad("pRt", "", 0.50, split, 1.00, 1.00);
    TPad* pRb = new TPad("pRb", "", 0.50, 0.00, 1.00, split);
    for (TPad* p : {pLt, pLb, pRt, pRb}) p->Draw();

    const char* proc_label = is_alt ? "Alt. banana cut (quadratic)" : "Nominal 2-band banana cut";

    DrawRatioColumn(pLt, pLb, h23, h24, kBlack, kRed+1,
                    "PbPb 2023", "PbPb 2024", proc_label);
    DrawRatioColumn(pRt, pRb, h23, h25, kBlack, kBlue+1,
                    "PbPb 2023", "PbPb 2025", proc_label);

    // Save
    const std::string base = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/single_b_analysis/";
    const std::string out_dir = is_alt
        ? base + "event_selection_alternative_banana/"
        : base + "event_selection/";
    gSystem->mkdir(out_dir.c_str(), true);
    const std::string outname = out_dir +
        (is_alt ? "fcal_comparison_pbpb_alt_banana.png"
                : "fcal_comparison_pbpb_nominal.png");
    c->SaveAs(outname.c_str());
    std::cout << "Saved: " << outname << std::endl;

    delete c;
    delete h23; delete h24; delete h25;
    for (auto* g : {cs23.g_cut1, cs23.g_cut4, cs23.g_cut5_lo, cs23.g_cut5_hi,
                    cs24.g_cut1, cs24.g_cut4, cs24.g_cut5_lo, cs24.g_cut5_hi,
                    cs25.g_cut1, cs25.g_cut4, cs25.g_cut5_lo, cs25.g_cut5_hi})
        delete g;
}

// =============================================================================
void plot_pbpb_fcal_comparison() {
    gStyle->SetOptStat(0);
    MakeComparisonPlot(false);
}

void plot_pbpb_fcal_comparison_alt() {
    gStyle->SetOptStat(0);
    MakeComparisonPlot(true);
}
