// derive_fcal_scaling.cxx
// Derives FCal-ET-dependent event weights for PbPb 2024/2025 alignment to 2023.
//
// Procedure:
//   1. Apply 5-cut event selection to raw skim for each year (23, 24, 25).
//   2. Fill FCal_Et = (A+C)*1e-6 [TeV] histograms.
//   3. For centrality 0-80% (FCal > kFCal_80_ref = 0.063208 TeV):
//      normalise per-year histograms to area=1; weight = h23_norm / hyr_norm per bin.
//   4. Save TGraph("g_fcal_weight") to fcal_weight_pbpb_20YY.root per year 24/25.
//   5. Output:
//      - fcal_before_scaling.png : raw 2-panel comparison (2024 left, 2025 right vs 2023)
//      - fcal_after_scaling_24.png, fcal_after_scaling_25.png : per-year ratio plots
//        (upper = weighted vs 2023 distribution; lower = weight(FCal) with w=1 reference)
//
// Usage (run from Analysis/):
//   source /usatlas/u/yuhanguo/setup.sh
//   root -l -b -q 'plotting_codes/fcal_scaling/derive_fcal_scaling.cxx+'

#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TParameter.h"
#include "../../NTupleProcessingCode/PbPbEventSelConfig.h"

// FCal_ET_Bins_PbPb2023[79]: boundary between 79% and 80% centrality [TeV].
// Events with FCal_Et > kFCal_80_ref are in centrality 0-80%.
static const double kFCal_80_ref = 0.063208;

static const std::string kBase    = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/";
static const std::string kPlotDir = kBase + "plots/fcal_scaling/";

// ---- input file lists -------------------------------------------------------
static std::map<int, std::vector<std::string>> BuildFileMap() {
    return {
        {23, {kBase+"pbpb_2023/data_pbpb23_part1.root",
              kBase+"pbpb_2023/data_pbpb23_part2.root",
              kBase+"pbpb_2023/data_pbpb23_part3.root"}},
        {24, {kBase+"pbpb_2024/data_pbpb24_part1.root",
              kBase+"pbpb_2024/data_pbpb24_part2.root"}},
        {25, {kBase+"pbpb_2025/data_pbpb25_part1.root",
              kBase+"pbpb_2025/data_pbpb25_part2.root",
              kBase+"pbpb_2025/data_pbpb25_part3.root",
              kBase+"pbpb_2025/data_pbpb25_part4.root",
              kBase+"pbpb_2025/data_pbpb25_part5.root",
              kBase+"pbpb_2025/data_pbpb25_part6.root"}},
    };
}

// ---- event sel cuts per year ------------------------------------------------
struct EvSelCuts {
    TGraph* g_cut1{nullptr}; double cut2_ns{1.8};
    float cut3_A{385.f}, cut3_C{385.f};
    TGraph* g_cut4{nullptr}, *g_cut5_lo{nullptr}, *g_cut5_hi{nullptr};
};

static EvSelCuts LoadEvSelCuts(int yr) {
    std::string path = PbPbEvSelCutsPath(2000 + yr);
    TFile* f = TFile::Open(path.c_str(), "READ");
    if (!f || f->IsZombie()) throw std::runtime_error("Cuts file not found: " + path);
    auto get_graph = [&](const char* k) -> TGraph* {
        TGraph* g = (TGraph*)f->Get(k);
        if (!g) throw std::runtime_error(std::string("Key not found: ") + k);
        return (TGraph*)g->Clone();
    };
    auto get_param = [&](const char* k) -> double {
        TParameter<double>* p = (TParameter<double>*)f->Get(k);
        if (!p) throw std::runtime_error(std::string("Key not found: ") + k);
        return p->GetVal();
    };
    EvSelCuts c;
    c.g_cut1    = get_graph(PbPbEvSelKey::kZDCFCalCut);
    c.cut2_ns   = get_param(PbPbEvSelKey::kZDCTimeCutNs);
    c.cut3_A    = (float)get_param(PbPbEvSelKey::kPreampACutADC);
    c.cut3_C    = (float)get_param(PbPbEvSelKey::kPreampCCutADC);
    c.g_cut4    = get_graph(PbPbEvSelKey::kNTrkFracCutLo);
    c.g_cut5_lo = get_graph(PbPbEvSelKey::kNTrkFCalCutLo);
    c.g_cut5_hi = get_graph(PbPbEvSelKey::kNTrkFCalCutHi);
    f->Close();
    std::cout << "Loaded cuts for 20" << yr << std::endl;
    return c;
}

// ---- fill FCal histogram after full 5-cut event selection -------------------
static TH1D* FillFCalHist(int yr, const EvSelCuts& cuts,
                           const std::vector<std::string>& files) {
    std::string name = "h_fcal_" + std::to_string(yr);
    TH1D* h = new TH1D(name.c_str(), ";FCal E_{T}^{A+C} [TeV];Events (normalised)",
                        600, 0., 6.);
    h->SetDirectory(nullptr);
    TChain chain("HeavyIonD3PD", "HeavyIonD3PD");
    for (const auto& f : files) {
        if (gSystem->AccessPathName(f.c_str())) { std::cerr << "Skipping: " << f << std::endl; continue; }
        chain.Add(f.c_str());
    }
    if (chain.GetEntries() == 0)
        throw std::runtime_error("Empty TChain for 20" + std::to_string(yr));
    chain.SetMakeClass(1);
    chain.SetBranchStatus("*", 0);
    Int_t b_HLT = 0; Float_t FCal_P = 0.f, FCal_N = 0.f;
    Float_t zdcE[2]{}, zdcT[2]{}, preamp[2][4]{};
    std::vector<int>* trk_numqual = nullptr;
    chain.SetBranchStatus("b_HLT_mu4_L1MU3V",         1);
    chain.SetBranchStatus("FCal_Et_P",                 1);
    chain.SetBranchStatus("FCal_Et_N",                 1);
    chain.SetBranchStatus("zdc_ZdcEnergy",             1);
    chain.SetBranchStatus("zdc_ZdcTime",               1);
    chain.SetBranchStatus("zdc_ZdcModulePreSampleAmp", 1);
    chain.SetBranchStatus("trk_numqual",               1);
    chain.SetBranchAddress("b_HLT_mu4_L1MU3V",         &b_HLT);
    chain.SetBranchAddress("FCal_Et_P",                 &FCal_P);
    chain.SetBranchAddress("FCal_Et_N",                 &FCal_N);
    chain.SetBranchAddress("zdc_ZdcEnergy",             zdcE);
    chain.SetBranchAddress("zdc_ZdcTime",               zdcT);
    chain.SetBranchAddress("zdc_ZdcModulePreSampleAmp", preamp);
    chain.SetBranchAddress("trk_numqual",               &trk_numqual);
    const Long64_t n = chain.GetEntries();
    std::cout << "20" << yr << ": processing " << n << " events..." << std::endl;
    Long64_t n_pass = 0;
    for (Long64_t i = 0; i < n; ++i) {
        chain.GetEntry(i);
        if (!b_HLT) continue;
        const float fcal = (FCal_P + FCal_N) * 1e-6f;
        const float zdcT_tot = (zdcE[0] + zdcE[1]) / 1000.f;
        if (zdcT_tot > (float)PbPbEvSelEvalCut(cuts.g_cut1, fcal)) continue;
        if (std::abs(zdcT[1]) >= (float)cuts.cut2_ns) continue;
        if (std::abs(zdcT[0]) >= (float)cuts.cut2_ns) continue;
        float pA = 0.f, pC = 0.f;
        for (int k = 0; k < 4; ++k) pA += preamp[1][k];
        for (int k = 0; k < 4; ++k) pC += preamp[0][k];
        if (pA >= cuts.cut3_A || pC >= cuts.cut3_C) continue;
        if (!trk_numqual || (int)trk_numqual->size() < 4) continue;
        const int nt = (*trk_numqual)[0], ntq = (*trk_numqual)[3];
        if (nt > 0 && (double)ntq/nt < PbPbEvSelEvalCut(cuts.g_cut4, nt)) continue;
        if ((double)ntq < PbPbEvSelEvalCut(cuts.g_cut5_lo, fcal)) continue;
        if ((double)ntq > PbPbEvSelEvalCut(cuts.g_cut5_hi, fcal)) continue;
        h->Fill(fcal);
        ++n_pass;
    }
    std::cout << "20" << yr << ": " << n_pass << " / " << n << " pass event selection" << std::endl;
    return h;
}

// ---- derive FCal-dependent weight for FCal > boundary ----------------------
// w(FCal) = h23_norm(FCal) / hyr_norm(FCal), normalised independently over FCal > boundary.
// Bins where hyr = 0 are skipped (TGraph::Eval interpolates across them).
// Weights are capped at 5 to suppress statistical noise at the high-FCal tail.
static TGraph* DeriveWeight(TH1D* h23, TH1D* hyr, double boundary) {
    const int b_lo = h23->FindBin(boundary + 1e-9);
    const double int23 = h23->Integral(b_lo, h23->GetNbinsX());
    const double intyr = hyr->Integral(b_lo, hyr->GetNbinsX());
    if (int23 == 0. || intyr == 0.) return nullptr;
    std::vector<double> xs, ys;
    for (int ib = b_lo; ib <= h23->GetNbinsX(); ++ib) {
        const double c23 = h23->GetBinContent(ib) / int23;
        const double cyr = hyr->GetBinContent(ib) / intyr;
        if (cyr <= 0.) continue;
        const double w = (c23 > 0.) ? std::min(c23 / cyr, 5.) : 0.;
        xs.push_back(h23->GetBinCenter(ib));
        ys.push_back(w);
    }
    if (xs.empty()) return nullptr;
    TGraph* g = new TGraph((int)xs.size(), xs.data(), ys.data());
    g->SetName("g_fcal_weight");
    g->SetTitle(";FCal E_{T}^{A+C} [TeV];Weight");
    return g;
}

// ---- build weighted FCal histogram ------------------------------------------
static TH1D* MakeWeightedHist(const TH1D* hyr, const TGraph* g_weight,
                               double boundary, const char* name) {
    TH1D* hw = (TH1D*)hyr->Clone(name);
    hw->SetDirectory(nullptr);
    for (int ib = 1; ib <= hyr->GetNbinsX(); ++ib) {
        const double fcal = hyr->GetBinCenter(ib);
        const double cnt  = hyr->GetBinContent(ib);
        if (fcal > boundary && cnt > 0.)
            hw->SetBinContent(ib, cnt * g_weight->Eval(fcal));
        else
            hw->SetBinContent(ib, 0.);
    }
    return hw;
}

// ---- before-scaling 2-panel comparison (raw) --------------------------------
static void SaveBeforePlot(TH1D* h23, TH1D* h24, TH1D* h25) {
    auto norm = [](TH1D* h) {
        TH1D* c = (TH1D*)h->Clone(); c->SetDirectory(nullptr);
        c->Scale(1./c->Integral()); return c;
    };
    TCanvas* cv = new TCanvas("c_before", "", 1400, 600);
    cv->Divide(2, 1);
    auto draw_panel = [&](int pad, TH1D* href, TH1D* hyr, int col, const char* lab) {
        cv->cd(pad); gPad->SetLogy();
        TH1D* r = norm(href); r->SetLineColor(kBlack); r->SetLineWidth(2);
        TH1D* y = norm(hyr);  y->SetLineColor(col);    y->SetLineWidth(2);
        r->GetXaxis()->SetRangeUser(0., 5.5); r->SetMinimum(1e-6);
        r->Draw("hist"); y->Draw("hist same");
        TLegend* leg = new TLegend(0.52, 0.72, 0.88, 0.88);
        leg->SetBorderSize(0); leg->SetFillStyle(0);
        leg->AddEntry(r, "PbPb 2023 (ref)", "l");
        leg->AddEntry(y, lab, "l");
        leg->Draw();
    };
    draw_panel(1, h23, h24, kRed+1,   "PbPb 2024 (raw)");
    draw_panel(2, h23, h25, kAzure+2, "PbPb 2025 (raw)");
    std::string out = kPlotDir + "fcal_before_scaling.png";
    cv->SaveAs(out.c_str()); std::cout << "Saved: " << out << std::endl;
    delete cv;
}

// ---- draw one ratio column into a given upper+lower pad pair ----------------
static void DrawRatioColumn(TPad* pad_top, TPad* pad_bot,
                             TH1D* h23, TH1D* hyr_weighted,
                             TGraph* g_weight, int yr, double boundary,
                             bool left_col) {
    const int yr_col = (yr == 24) ? kRed+1 : kAzure+2;
    const int b_lo   = h23->FindBin(boundary + 1e-9);
    const char* suf  = left_col ? "L" : "R";

    // --- upper pad ---
    pad_top->cd();
    TH1D* href = (TH1D*)h23->Clone(Form("href_%s", suf));          href->SetDirectory(nullptr);
    TH1D* hcmp = (TH1D*)hyr_weighted->Clone(Form("hcmp_%s", suf)); hcmp->SetDirectory(nullptr);
    double i23 = href->Integral(b_lo, href->GetNbinsX());
    double iyr = hcmp->Integral(b_lo, hcmp->GetNbinsX());
    if (i23 > 0.) href->Scale(1./i23);
    if (iyr > 0.) hcmp->Scale(1./iyr);

    href->SetLineColor(kBlack); href->SetLineWidth(2);
    hcmp->SetLineColor(yr_col); hcmp->SetLineWidth(2);
    href->GetXaxis()->SetRangeUser(0.04, 5.5);
    href->GetXaxis()->SetTitle(""); href->GetXaxis()->SetLabelSize(0.);
    if (left_col) {
        href->GetYaxis()->SetTitle("Events (norm., FCal > 80% bdy.)");
        href->GetYaxis()->SetTitleSize(0.06); href->GetYaxis()->SetTitleOffset(0.95);
        href->GetYaxis()->SetLabelSize(0.055);
    } else {
        href->GetYaxis()->SetTitle(""); href->GetYaxis()->SetLabelSize(0.);
    }
    href->SetMinimum(1e-7);
    href->Draw("hist"); hcmp->Draw("hist same");

    TLegend* leg = new TLegend(0.38, 0.72, 0.92, 0.88);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.055);
    leg->AddEntry(href, "PbPb 2023 (ref)", "l");
    leg->AddEntry(hcmp, Form("PbPb 20%d (weighted)", yr), "l");
    leg->Draw();

    TLatex tl; tl.SetNDC(); tl.SetTextSize(0.058);
    tl.DrawLatex(0.12, 0.88, Form("Pb+Pb 20%d / 2023  |  0-80%%", yr));

    // --- lower pad ---
    pad_bot->cd();
    TGraph* gw = (TGraph*)g_weight->Clone(Form("gw_%s", suf));
    gw->SetLineColor(yr_col); gw->SetLineWidth(2); gw->SetMarkerStyle(0);

    TH1D* hfr = new TH1D(Form("hframe_%s", suf), "", 600, 0.04, 5.5);
    hfr->SetDirectory(nullptr);
    hfr->GetXaxis()->SetTitle("FCal E_{T}^{A+C} [TeV]");
    hfr->GetXaxis()->SetTitleSize(0.13); hfr->GetXaxis()->SetLabelSize(0.11);
    if (left_col) {
        hfr->GetYaxis()->SetTitle("Weight");
        hfr->GetYaxis()->SetTitleSize(0.12); hfr->GetYaxis()->SetTitleOffset(0.45);
        hfr->GetYaxis()->SetLabelSize(0.10);
    } else {
        hfr->GetYaxis()->SetTitle(""); hfr->GetYaxis()->SetLabelSize(0.);
    }
    hfr->GetYaxis()->SetNdivisions(504);
    hfr->SetMinimum(0.); hfr->SetMaximum(2.5);
    hfr->Draw("axis");

    TLine* l1 = new TLine(0.04, 1., 5.5, 1.);
    l1->SetLineStyle(2); l1->SetLineColor(kBlack); l1->Draw();
    gw->Draw("L same");
}

// ---- combined after plot: 2-column ratio canvas (2024 left, 2025 right) -----
static void SaveCombinedAfterPlot(TH1D* h23,
                                   TH1D* h24w, TGraph* g24,
                                   TH1D* h25w, TGraph* g25,
                                   double boundary) {
    TCanvas* cv = new TCanvas("c_after", "", 1400, 700);

    // upper row: y in [0.35, 1.0]; lower row: y in [0.0, 0.35]
    TPad* pTL = new TPad("pTL","", 0.0, 0.35, 0.5, 1.0);
    TPad* pTR = new TPad("pTR","", 0.5, 0.35, 1.0, 1.0);
    TPad* pBL = new TPad("pBL","", 0.0, 0.0,  0.5, 0.35);
    TPad* pBR = new TPad("pBR","", 0.5, 0.0,  1.0, 0.35);

    pTL->SetBottomMargin(0.015); pTL->SetTopMargin(0.09); pTL->SetRightMargin(0.02); pTL->SetLogy();
    pTR->SetBottomMargin(0.015); pTR->SetTopMargin(0.09); pTR->SetLeftMargin(0.02);  pTR->SetLogy();
    pBL->SetTopMargin(0.015);    pBL->SetBottomMargin(0.35); pBL->SetRightMargin(0.02);
    pBR->SetTopMargin(0.015);    pBR->SetBottomMargin(0.35); pBR->SetLeftMargin(0.02);

    pTL->Draw(); pTR->Draw(); pBL->Draw(); pBR->Draw();

    DrawRatioColumn(pTL, pBL, h23, h24w, g24, 24, boundary, true);
    DrawRatioColumn(pTR, pBR, h23, h25w, g25, 25, boundary, false);

    std::string out = kPlotDir + "fcal_after_scaling.png";
    cv->SaveAs(out.c_str()); std::cout << "Saved: " << out << std::endl;
    delete cv;
}

// ---- main -------------------------------------------------------------------
void derive_fcal_scaling() {
    gStyle->SetOptStat(0);
    gSystem->mkdir(kPlotDir.c_str(), true);

    const auto fmap = BuildFileMap();

    // fill FCal histograms
    std::map<int, EvSelCuts> cuts;
    std::map<int, TH1D*>     hists;
    for (int yr : {23, 24, 25}) {
        cuts[yr]  = LoadEvSelCuts(yr);
        hists[yr] = FillFCalHist(yr, cuts[yr], fmap.at(yr));
    }

    // derive weights and save
    std::map<int, TGraph*> weights;
    for (int yr : {24, 25}) {
        weights[yr] = DeriveWeight(hists[23], hists[yr], kFCal_80_ref);
        if (!weights[yr]) throw std::runtime_error("Weight derivation failed for 20" + std::to_string(yr));
        std::string path = PbPbFCalWeightPath(2000 + yr);
        TFile* fout = TFile::Open(path.c_str(), "RECREATE");
        if (!fout || fout->IsZombie()) throw std::runtime_error("Cannot create: " + path);
        weights[yr]->Write("g_fcal_weight");
        fout->Close();
        std::cout << "Saved weight for 20" << yr << " → " << path << std::endl;
    }

    // before plot (raw distributions)
    SaveBeforePlot(hists[23], hists[24], hists[25]);

    // combined after ratio plot (2024 left, 2025 right)
    TH1D* h24w = MakeWeightedHist(hists[24], weights[24], kFCal_80_ref, "h_fcal_24_weighted");
    TH1D* h25w = MakeWeightedHist(hists[25], weights[25], kFCal_80_ref, "h_fcal_25_weighted");
    SaveCombinedAfterPlot(hists[23], h24w, weights[24], h25w, weights[25], kFCal_80_ref);
    delete h24w; delete h25w;
}
