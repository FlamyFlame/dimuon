// plot_zdc_run_dependence_pbpb23.cxx
//
// Per-run ZDC variable study for pbpb2023.
// Selection: b_HLT_mu4_L1MU3V only; no event-selection cuts.
//
// Four sets of outputs (saved to plots/.../pbpb_2023/run_dependence/):
//   zdc_time_run_dep_mu_sigma.png          — ZDC time Gaussian fit mu/sigma vs run
//   zdc_preamp_run_dep_mu_sigma.png        — ZDC preamp Gaussian fit mu/sigma vs run
//   zdc_energy_1n_run_dep_mu_sigma.png     — ZDC energy 1n-peak fit mu/sigma vs run
//   zdc_energy_mean_run_dep.png            — ZDC energy arithmetic mean vs run (no fit)
//
// Usage (run from Analysis/):
//   source /usatlas/u/yuhanguo/setup.sh
//   root -l -b -q 'plotting_codes/event_selection/plot_zdc_run_dependence_pbpb23.cxx+'

#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <numeric>
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TSystem.h"
#include "TStyle.h"

static const std::string kBase    = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/";
static const std::string kPlotDir = kBase + "plots/single_b_analysis/event_selection/pbpb_2023/run_dependence/";

// ---- Gaussian fit result ----
struct GaussFit { double mu, sig, amp; bool valid; };

// Generic Gaussian fit: seed with [seedLo, seedHi], then refine within mu ± [nLo, nHi]*sigma
static GaussFit DoGaussFit(TH1D* h, double seedLo, double seedHi,
                            double nLo, double nHi, const char* fname) {
    if (h->GetEntries() < 30) return {0., 0., 0., false};
    // Use local max within seed range as amplitude seed
    int b1 = h->FindBin(seedLo + 0.5*(h->GetBinWidth(1)));
    int b2 = h->FindBin(seedHi - 0.5*(h->GetBinWidth(1)));
    double localMax = 0.;
    int    bMax     = b1;
    for (int b = b1; b <= b2; ++b)
        if (h->GetBinContent(b) > localMax) { localMax = h->GetBinContent(b); bMax = b; }
    if (localMax < 1.) return {0., 0., 0., false};

    TF1 g1(Form("%s_seed", fname), "gaus", seedLo, seedHi);
    g1.SetParameters(localMax, h->GetBinCenter(bMax), 0.25*(seedHi-seedLo));
    h->Fit(&g1, "RQN");
    double mu1  = g1.GetParameter(1);
    double sig1 = std::abs(g1.GetParameter(2));
    // Reject if mu drifted outside seed window or sigma is degenerate/too large
    if (mu1 < seedLo || mu1 > seedHi) return {0., 0., 0., false};
    if (sig1 < 1e-9 || sig1 > (seedHi - seedLo)) return {0., 0., 0., false};

    TF1* gf = new TF1(fname, "gaus", mu1 - nLo*sig1, mu1 + nHi*sig1);
    gf->SetParameters(g1.GetParameter(0), mu1, sig1);
    int st = h->Fit(gf, "RQN");
    double mu_f  = gf->GetParameter(1);
    double sig_f = std::abs(gf->GetParameter(2));
    double amp_f = gf->GetParameter(0);
    delete gf;
    if (st != 0 || sig_f < 1e-9) return {0., 0., 0., false};
    // Final sanity: mu must stay within seed window; sigma must be narrower than it
    if (mu_f < seedLo || mu_f > seedHi) return {0., 0., 0., false};
    if (sig_f > (seedHi - seedLo))      return {0., 0., 0., false};
    return {mu_f, sig_f, amp_f, true};
}

// ---- Per-run histograms for one variable ----
struct VarHists {
    TH1D* hA = nullptr;
    TH1D* hC = nullptr;
    double sumA = 0., sumC = 0.;
    long long nA = 0, nC = 0;
};

// ---- Draw mu/sigma two-subplot summary ----
static void DrawMuSigma(const std::vector<int>& runs,
                         const std::vector<GaussFit>& fitA,
                         const std::vector<GaussFit>& fitC,
                         const char* ytitle_mu,
                         const char* ytitle_sig,
                         const char* header,
                         const char* outname) {
    const int nRuns = (int)runs.size();

    std::vector<double> xs(nRuns), muA_v(nRuns), muC_v(nRuns), sigA_v(nRuns), sigC_v(nRuns);
    for (int i = 0; i < nRuns; ++i) {
        xs[i]     = i + 1.;
        muA_v[i]  = fitA[i].valid ? fitA[i].mu  : 0.;
        muC_v[i]  = fitC[i].valid ? fitC[i].mu  : 0.;
        sigA_v[i] = fitA[i].valid ? fitA[i].sig : 0.;
        sigC_v[i] = fitC[i].valid ? fitC[i].sig : 0.;
    }

    auto makeGraph = [&](const std::vector<double>& y, int color, int marker) -> TGraph* {
        TGraph* g = new TGraph(nRuns, xs.data(), y.data());
        g->SetMarkerStyle(marker); g->SetMarkerSize(0.9);
        g->SetMarkerColor(color); g->SetLineColor(color); g->SetLineWidth(1);
        return g;
    };

    TGraph* gMuA  = makeGraph(muA_v,  kBlack,  20);
    TGraph* gMuC  = makeGraph(muC_v,  kRed+1,  21);
    TGraph* gSigA = makeGraph(sigA_v, kBlack,  20);
    TGraph* gSigC = makeGraph(sigC_v, kRed+1,  21);

    auto makeFrame = [&](TPad* pad,
                         const char* ytit,
                         const std::vector<double>& vA,
                         const std::vector<double>& vC) -> TH1D* {
        pad->cd();
        pad->SetLeftMargin(0.13); pad->SetRightMargin(0.03);
        pad->SetBottomMargin(0.22); pad->SetTopMargin(0.07);
        double lo = 1e38, hi = -1e38;
        for (int i = 0; i < nRuns; ++i) {
            if (vA[i] != 0.) { lo = std::min(lo, vA[i]); hi = std::max(hi, vA[i]); }
            if (vC[i] != 0.) { lo = std::min(lo, vC[i]); hi = std::max(hi, vC[i]); }
        }
        if (lo > hi) { lo = 0.; hi = 1.; }
        double pad_y = 0.20 * (hi - lo);
        lo -= pad_y; hi += pad_y;
        TH1D* fr = new TH1D(Form("fr_%s_%s", ytit, outname), Form(";Run;%s", ytit),
                             nRuns, 0.5, nRuns + 0.5);
        for (int i = 0; i < nRuns; ++i)
            fr->GetXaxis()->SetBinLabel(i+1, Form("%d", runs[i]));
        fr->GetYaxis()->SetRangeUser(lo, hi);
        fr->GetXaxis()->SetLabelSize(0.038);
        fr->GetXaxis()->LabelsOption("v");
        fr->GetYaxis()->SetTitleSize(0.052); fr->GetYaxis()->SetLabelSize(0.042);
        fr->GetYaxis()->SetTitleOffset(1.20);
        fr->GetXaxis()->SetTitleSize(0.);
        fr->SetLineColor(0); fr->SetMarkerColor(0);
        fr->Draw("AXIS");
        return fr;
    };

    TCanvas* cv = new TCanvas(Form("cv_%s", outname), "", 1600, 600);
    TPad* pL = new TPad("pL_mu",  "", 0.00, 0.00, 0.50, 1.00);
    TPad* pR = new TPad("pR_sig", "", 0.50, 0.00, 1.00, 1.00);
    pL->Draw(); pR->Draw();

    TH1D* frMu = makeFrame(pL, ytitle_mu, muA_v, muC_v);
    gMuA->Draw("P same"); gMuC->Draw("P same");
    pL->cd();
    TLegend* legL = new TLegend(0.15, 0.78, 0.42, 0.91);
    legL->SetBorderSize(0); legL->SetFillStyle(0); legL->SetTextSize(0.045);
    legL->AddEntry(gMuA, "Side A", "p"); legL->AddEntry(gMuC, "Side C", "p");
    legL->Draw();
    TLatex tlL; tlL.SetNDC(); tlL.SetTextSize(0.042); tlL.SetTextAlign(11);
    tlL.DrawLatex(0.15, 0.93, Form("Pb+Pb 2023  %s  #mu per run", header));

    TH1D* frSig = makeFrame(pR, ytitle_sig, sigA_v, sigC_v);
    gSigA->Draw("P same"); gSigC->Draw("P same");
    pR->cd();
    TLegend* legR = new TLegend(0.15, 0.78, 0.42, 0.91);
    legR->SetBorderSize(0); legR->SetFillStyle(0); legR->SetTextSize(0.045);
    legR->AddEntry(gSigA, "Side A", "p"); legR->AddEntry(gSigC, "Side C", "p");
    legR->Draw();
    TLatex tlR; tlR.SetNDC(); tlR.SetTextSize(0.042); tlR.SetTextAlign(11);
    tlR.DrawLatex(0.15, 0.93, Form("Pb+Pb 2023  %s  #sigma per run", header));

    std::string out = kPlotDir + outname + ".png";
    cv->SaveAs(out.c_str());
    std::cout << "Saved: " << out << std::endl;
    delete frMu; delete frSig;
    delete cv;
}

// ---- Draw mean-only two-subplot summary ----
static void DrawMeanPlot(const std::vector<int>& runs,
                          const std::vector<double>& meanA,
                          const std::vector<double>& meanC,
                          const char* ytitle,
                          const char* header,
                          const char* outname) {
    const int nRuns = (int)runs.size();
    std::vector<double> xs(nRuns);
    for (int i = 0; i < nRuns; ++i) xs[i] = i + 1.;

    auto makeGraph = [&](const std::vector<double>& y, int color, int marker) -> TGraph* {
        TGraph* g = new TGraph(nRuns, xs.data(), y.data());
        g->SetMarkerStyle(marker); g->SetMarkerSize(0.9);
        g->SetMarkerColor(color); g->SetLineColor(color); g->SetLineWidth(1);
        return g;
    };
    TGraph* gA = makeGraph(meanA, kBlack, 20);
    TGraph* gC = makeGraph(meanC, kRed+1, 21);

    TCanvas* cv = new TCanvas(Form("cv_%s", outname), "", 900, 600);
    cv->SetLeftMargin(0.12); cv->SetRightMargin(0.03);
    cv->SetBottomMargin(0.22); cv->SetTopMargin(0.07);

    double lo = 1e38, hi = -1e38;
    for (int i = 0; i < nRuns; ++i) {
        if (meanA[i] != 0.) { lo = std::min(lo, meanA[i]); hi = std::max(hi, meanA[i]); }
        if (meanC[i] != 0.) { lo = std::min(lo, meanC[i]); hi = std::max(hi, meanC[i]); }
    }
    if (lo > hi) { lo = 0.; hi = 1.; }
    double pad_y = 0.20 * (hi - lo);
    lo -= pad_y; hi += pad_y;

    TH1D* fr = new TH1D(Form("fr_%s", outname), Form(";Run;%s", ytitle),
                         nRuns, 0.5, nRuns + 0.5);
    for (int i = 0; i < nRuns; ++i)
        fr->GetXaxis()->SetBinLabel(i+1, Form("%d", runs[i]));
    fr->GetYaxis()->SetRangeUser(lo, hi);
    fr->GetXaxis()->SetLabelSize(0.038);
    fr->GetXaxis()->LabelsOption("v");
    fr->GetYaxis()->SetTitleSize(0.050); fr->GetYaxis()->SetLabelSize(0.040);
    fr->GetYaxis()->SetTitleOffset(1.25);
    fr->GetXaxis()->SetTitleSize(0.);
    fr->SetLineColor(0); fr->SetMarkerColor(0);
    fr->Draw("AXIS");
    gA->Draw("P same"); gC->Draw("P same");

    TLegend* leg = new TLegend(0.14, 0.78, 0.38, 0.91);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.043);
    leg->AddEntry(gA, "Side A", "p"); leg->AddEntry(gC, "Side C", "p");
    leg->Draw();
    TLatex tl; tl.SetNDC(); tl.SetTextSize(0.040); tl.SetTextAlign(11);
    tl.DrawLatex(0.14, 0.93, Form("Pb+Pb 2023  %s  mean per run", header));

    std::string out = kPlotDir + outname + ".png";
    cv->SaveAs(out.c_str());
    std::cout << "Saved: " << out << std::endl;
    delete fr; delete cv;
}

// ---- Main ----
void plot_zdc_run_dependence_pbpb23() {
    gStyle->SetOptStat(0);
    gSystem->mkdir(kPlotDir.c_str(), true);

    // Build TChain
    TChain chain("HeavyIonD3PD", "HeavyIonD3PD");
    for (int p = 1; p <= 4; ++p) {
        std::string f = kBase + "pbpb_2023/data_pbpb23_part" + std::to_string(p) + ".root";
        if (!gSystem->AccessPathName(f.c_str())) chain.Add(f.c_str());
        else std::cerr << "Skipping: " << f << std::endl;
    }
    chain.SetMakeClass(1);
    chain.SetBranchStatus("*", 0);

    Int_t   b_RunNum = 0;
    Bool_t  b_HLT   = false;
    Float_t zdcE[2]{}, zdcT[2]{}, preamp[2][4]{};

    chain.SetBranchStatus("RunNumber",                 1);
    chain.SetBranchStatus("b_HLT_mu4_L1MU3V",          1);
    chain.SetBranchStatus("zdc_ZdcEnergy",             1);
    chain.SetBranchStatus("zdc_ZdcTime",               1);
    chain.SetBranchStatus("zdc_ZdcModulePreSampleAmp", 1);
    chain.SetBranchAddress("RunNumber",                 &b_RunNum);
    chain.SetBranchAddress("b_HLT_mu4_L1MU3V",          &b_HLT);
    chain.SetBranchAddress("zdc_ZdcEnergy",             zdcE);
    chain.SetBranchAddress("zdc_ZdcTime",               zdcT);
    chain.SetBranchAddress("zdc_ZdcModulePreSampleAmp", preamp);

    // Per-run histogram maps
    std::map<int, TH1D*> htA_map, htC_map;     // ZDC time
    std::map<int, TH1D*> hpA_map, hpC_map;     // ZDC preamp
    std::map<int, TH1D*> heA_map, heC_map;     // ZDC energy
    std::map<int, double> sumEA_map, sumEC_map; // for mean calculation
    std::map<int, long long> nEA_map, nEC_map;

    const Long64_t n = chain.GetEntries();
    std::cout << "pbpb2023: " << n << " events, scanning..." << std::flush;
    long long n_sel = 0;
    for (Long64_t i = 0; i < n; ++i) {
        chain.GetEntry(i);
        if (!b_HLT) continue;
        const int run = b_RunNum;

        // Book histograms on first encounter of this run
        if (htA_map.find(run) == htA_map.end()) {
            htA_map[run] = new TH1D(Form("htA_%d", run), "", 200, -10., 10.);
            htC_map[run] = new TH1D(Form("htC_%d", run), "", 200, -10., 10.);
            hpA_map[run] = new TH1D(Form("hpA_%d", run), "", 200, -1000., 3000.);
            hpC_map[run] = new TH1D(Form("hpC_%d", run), "", 200, -1000., 3000.);
            heA_map[run] = new TH1D(Form("heA_%d", run), "", 200, 0., 8000.);
            heC_map[run] = new TH1D(Form("heC_%d", run), "", 200, 0., 8000.);
            htA_map[run]->SetDirectory(nullptr);
            htC_map[run]->SetDirectory(nullptr);
            hpA_map[run]->SetDirectory(nullptr);
            hpC_map[run]->SetDirectory(nullptr);
            heA_map[run]->SetDirectory(nullptr);
            heC_map[run]->SetDirectory(nullptr);
            sumEA_map[run] = 0.; sumEC_map[run] = 0.;
            nEA_map[run]   = 0;  nEC_map[run]   = 0;
        }

        // ZDC time: index 1=A, 0=C
        htA_map[run]->Fill(zdcT[1]);
        htC_map[run]->Fill(zdcT[0]);

        // ZDC preamp sum: preamp[1][k]=A, preamp[0][k]=C
        float pA = 0.f, pC = 0.f;
        for (int k = 0; k < 4; ++k) pA += preamp[1][k];
        for (int k = 0; k < 4; ++k) pC += preamp[0][k];
        hpA_map[run]->Fill(pA);
        hpC_map[run]->Fill(pC);

        // ZDC energy: index 1=A, 0=C (already in GeV from branch)
        heA_map[run]->Fill(zdcE[1]);
        heC_map[run]->Fill(zdcE[0]);
        sumEA_map[run] += zdcE[1];  nEA_map[run]++;
        sumEC_map[run] += zdcE[0];  nEC_map[run]++;

        ++n_sel;
    }
    std::cout << "  " << n_sel << " pass b_HLT_mu4_L1MU3V\n";

    // Sort runs
    std::vector<int> runs;
    for (auto& kv : htA_map) runs.push_back(kv.first);
    std::sort(runs.begin(), runs.end());
    std::cout << "Found " << runs.size() << " unique runs\n";

    const int nRuns = (int)runs.size();

    // ---- Fit ZDC time ----
    std::vector<GaussFit> ftA(nRuns), ftC(nRuns);
    for (int i = 0; i < nRuns; ++i) {
        int run = runs[i];
        ftA[i] = DoGaussFit(htA_map[run], -5., 5., 2.5, 2.5, Form("gtA_%d", run));
        ftC[i] = DoGaussFit(htC_map[run], -5., 5., 2.5, 2.5, Form("gtC_%d", run));
    }
    DrawMuSigma(runs, ftA, ftC,
                "#mu_{t} [ns]", "#sigma_{t} [ns]",
                "ZDC time",
                "zdc_time_run_dep_mu_sigma");

    // ---- Fit ZDC preamp ----
    std::vector<GaussFit> fpA(nRuns), fpC(nRuns);
    for (int i = 0; i < nRuns; ++i) {
        int run = runs[i];
        hpA_map[run]->Rebin(2); hpA_map[run]->Scale(0.5);
        hpC_map[run]->Rebin(2); hpC_map[run]->Scale(0.5);
        fpA[i] = DoGaussFit(hpA_map[run], -800., 1500., 3., 1.5, Form("gpA_%d", run));
        fpC[i] = DoGaussFit(hpC_map[run], -800., 1500., 3., 1.5, Form("gpC_%d", run));
    }
    DrawMuSigma(runs, fpA, fpC,
                "#mu_{p} [ADC]", "#sigma_{p} [ADC]",
                "ZDC PreSampleAmp",
                "zdc_preamp_run_dep_mu_sigma");

    // ---- Fit ZDC energy 1n peak ----
    std::vector<GaussFit> feA(nRuns), feC(nRuns);
    for (int i = 0; i < nRuns; ++i) {
        int run = runs[i];
        // Seed: 1n peak near 2680 GeV; narrow window to avoid multi-n contamination
        feA[i] = DoGaussFit(heA_map[run], 1500., 3500., 2., 1.5, Form("geA_%d", run));
        feC[i] = DoGaussFit(heC_map[run], 1500., 3500., 2., 1.5, Form("geC_%d", run));
    }
    DrawMuSigma(runs, feA, feC,
                "#mu_{E} [GeV]", "#sigma_{E} [GeV]",
                "ZDC energy 1n peak",
                "zdc_energy_1n_run_dep_mu_sigma");

    // ---- Mean ZDC energy (no fit) ----
    std::vector<double> meanA(nRuns, 0.), meanC(nRuns, 0.);
    for (int i = 0; i < nRuns; ++i) {
        int run = runs[i];
        if (nEA_map[run] > 0) meanA[i] = sumEA_map[run] / nEA_map[run];
        if (nEC_map[run] > 0) meanC[i] = sumEC_map[run] / nEC_map[run];
    }
    DrawMeanPlot(runs, meanA, meanC,
                 "#LT E_{ZDC} #GT [GeV]",
                 "ZDC energy mean",
                 "zdc_energy_mean_run_dep");

    // Cleanup
    for (auto& kv : htA_map) delete kv.second;
    for (auto& kv : htC_map) delete kv.second;
    for (auto& kv : hpA_map) delete kv.second;
    for (auto& kv : hpC_map) delete kv.second;
    for (auto& kv : heA_map) delete kv.second;
    for (auto& kv : heC_map) delete kv.second;
}
