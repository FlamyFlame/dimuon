// plot_zdc_preamp_per_run_pbpb25.cxx
//
// Per-run ZDC preamp study for pbpb2025.
// Investigates whether run-to-run variations explain the wide Gaussian peak
// seen in the combined 2025 sample.
//
// Selection: trigger + Cut1 (ZDC-FCal banana) + Cut2 (ZDC time), same as the
// nominal event selection up to (but not including) the preamp cut itself.
//
// Output (saved to plots/single_b_analysis/event_selection/pbpb_2025/):
//   zdc_preamp_per_run_dists_p{1,2,3}.png         — hard-cut vertical lines
//   zdc_preamp_per_run_dists_mu7sig_p{1,2,3}.png  — mu+7sigma vertical lines
//   zdc_preamp_per_run_mu_sigma.png               — fitted mu and sigma vs run
//
// Usage (run from Analysis/):
//   source /usatlas/u/yuhanguo/setup.sh
//   root -l -b -q 'plotting_codes/event_selection/plot_zdc_preamp_per_run_pbpb25.cxx+'

#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <iostream>
#include <cmath>
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
#include "TAxis.h"
#include "TSystem.h"
#include "TStyle.h"
#include "../../NTupleProcessingCode/PbPbEventSelConfig.h"

static const std::string kBase    = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/";
static const std::string kPlotDir = kBase + "plots/single_b_analysis/event_selection/pbpb_2025/";

static const double kZoomLo = -800., kZoomHi = 1500.;

// ---- Gaussian fit result ----
struct GaussFit { double mu, sig, amp; bool valid; };

static GaussFit DoGaussFit(TH1D* h, const char* fname) {
    if (h->GetEntries() < 50) return {0., 0., 0., false};
    double seed_hi = kZoomLo + 0.60 * (kZoomHi - kZoomLo);
    TF1 g1(Form("%s_p1", fname), "gaus", kZoomLo, seed_hi);
    g1.SetParameters(h->GetMaximum(), (kZoomLo + seed_hi) * 0.5, (seed_hi - kZoomLo) * 0.25);
    h->Fit(&g1, "RQN");
    double mu1  = g1.GetParameter(1);
    double sig1 = std::abs(g1.GetParameter(2));
    if (sig1 < 1.) return {0., 0., 0., false};

    TF1* gf = new TF1(fname, "gaus", mu1 - 3.*sig1, mu1 + 1.5*sig1);
    gf->SetParameters(g1.GetParameter(0), mu1, sig1);
    int st = h->Fit(gf, "RQN");
    double mu_f  = gf->GetParameter(1);
    double sig_f = std::abs(gf->GetParameter(2));
    double amp_f = gf->GetParameter(0);
    delete gf;
    if (st != 0 || sig_f < 1.) return {0., 0., 0., false};
    return {mu_f, sig_f, amp_f, true};
}

// ---- Per-run result ----
struct RunResult {
    int    run;
    GaussFit fitA, fitC;
    TH1D*  hA;
    TH1D*  hC;
};

// ---- Draw mode ----
enum DrawMode { kHardCut, kMu7Sig };

// ---- Distribution canvas pages ----
static void DrawDistPages(const std::vector<RunResult>& results,
                          double cut_A, double cut_C,
                          DrawMode mode, const char* suffix) {
    const int nRuns  = (int)results.size();
    const int nCols  = 5, nRows = 3, perPage = nCols * nRows;
    const int nPages = (nRuns + perPage - 1) / perPage;

    for (int pg = 0; pg < nPages; ++pg) {
        TCanvas* cv = new TCanvas(Form("cv_%s_p%d", suffix, pg+1), "", 1800, 1080);
        cv->Divide(nCols, nRows, 0.002, 0.002);

        for (int slot = 0; slot < perPage; ++slot) {
            const int idx = pg * perPage + slot;
            cv->cd(slot + 1);
            if (idx >= nRuns) continue;

            gPad->SetLogy();
            gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.03);
            gPad->SetTopMargin(0.10);  gPad->SetBottomMargin(0.14);

            const RunResult& R = results[idx];
            TH1D* hA = R.hA;
            TH1D* hC = R.hC;

            auto localMax = [](TH1D* h) {
                int b1 = h->FindBin(kZoomLo + 0.5), b2 = h->FindBin(kZoomHi - 0.5);
                double m = 0.;
                for (int b = b1; b <= b2; ++b) m = std::max(m, h->GetBinContent(b));
                return m;
            };
            double yhi = std::max(localMax(hA), localMax(hC)) * 5.;
            if (yhi < 1.) yhi = 10.;

            hA->GetXaxis()->SetRangeUser(kZoomLo, kZoomHi);
            hA->GetYaxis()->SetRangeUser(0.5, yhi);
            hA->SetMarkerStyle(20); hA->SetMarkerSize(0.3);
            hA->SetMarkerColor(kBlack); hA->SetLineColor(kBlack); hA->SetLineWidth(1);
            hA->GetXaxis()->SetTitle("Preamp sum [ADC]");
            hA->GetXaxis()->SetTitleSize(0.060); hA->GetXaxis()->SetLabelSize(0.053);
            hA->GetYaxis()->SetTitleSize(0.055); hA->GetYaxis()->SetLabelSize(0.050);
            hA->GetYaxis()->SetTitle("Events / 60 ADC");
            hA->GetYaxis()->SetTitleOffset(1.20);
            hA->Draw("E");

            hC->GetXaxis()->SetRangeUser(kZoomLo, kZoomHi);
            hC->SetMarkerStyle(20); hC->SetMarkerSize(0.3);
            hC->SetMarkerColor(kRed+1); hC->SetLineColor(kRed+1); hC->SetLineWidth(1);
            hC->Draw("E SAME");

            // Gaussian fit curves
            auto drawFitCurve = [&](const GaussFit& gf, TH1D* h, int color, const char* tag) {
                if (!gf.valid) return;
                TF1* gd = new TF1(Form("gd_%s_%s_%d", tag, suffix, R.run),
                                  "gaus", gf.mu - 3.*gf.sig, gf.mu + 1.5*gf.sig);
                gd->SetParameters(gf.amp, gf.mu, gf.sig);
                h->Fit(gd, "RQN");
                gd->SetLineColor(color); gd->SetLineWidth(2);
                gd->DrawClone("SAME");
                delete gd;
            };
            drawFitCurve(R.fitA, hA, kBlue+1,  "A");
            drawFitCurve(R.fitC, hC, kGreen+2, "C");

            // Vertical lines and value labels depending on mode
            double lineA = 0., lineC = 0.;
            if (mode == kHardCut) {
                lineA = cut_A; lineC = cut_C;
            } else { // kMu7Sig
                lineA = R.fitA.valid ? R.fitA.mu + 7.*R.fitA.sig : cut_A;
                lineC = R.fitC.valid ? R.fitC.mu + 7.*R.fitC.sig : cut_C;
            }

            double cutLineHi = yhi * 0.9;
            TLine lA(lineA, 0.5, lineA, cutLineHi);
            lA.SetLineColor(kBlack); lA.SetLineWidth(1); lA.SetLineStyle(2);
            lA.DrawClone();
            TLine lC(lineC, 0.5, lineC, cutLineHi);
            lC.SetLineColor(kRed+1); lC.SetLineWidth(1); lC.SetLineStyle(2);
            lC.DrawClone();

            // Text annotations
            TLatex tl; tl.SetNDC(); tl.SetTextSize(0.067); tl.SetTextAlign(11);
            tl.DrawLatex(0.16, 0.90, Form("Run %d", R.run));
            tl.SetTextSize(0.051);
            if (R.fitA.valid)
                tl.DrawLatex(0.16, 0.81,
                    Form("A: #mu=%.0f, #sigma=%.0f", R.fitA.mu, R.fitA.sig));
            if (R.fitC.valid)
                tl.DrawLatex(0.16, 0.72,
                    Form("#color[2]{C: #mu=%.0f, #sigma=%.0f}", R.fitC.mu, R.fitC.sig));
            // Show the mu+7sigma value for kMu7Sig mode
            if (mode == kMu7Sig) {
                tl.SetTextSize(0.046);
                if (R.fitA.valid)
                    tl.DrawLatex(0.16, 0.63, Form("#mu+7#sigma_{A}=%.0f", lineA));
                if (R.fitC.valid)
                    tl.DrawLatex(0.16, 0.55,
                        Form("#color[2]{#mu+7#sigma_{C}=%.0f}", lineC));
            }
        }

        std::string outpath = kPlotDir + Form("zdc_preamp_per_run_%s_p%d.png", suffix, pg+1);
        cv->SaveAs(outpath.c_str());
        std::cout << "Saved: " << outpath << std::endl;
        delete cv;
    }
}

// ---- mu/sigma summary ----
static void DrawMuSigma(const std::vector<RunResult>& results) {
    const int nRuns = (int)results.size();
    std::vector<double> xs(nRuns), muA_v(nRuns), muC_v(nRuns), sigA_v(nRuns), sigC_v(nRuns);
    for (int i = 0; i < nRuns; ++i) {
        xs[i]     = i + 1.;
        muA_v[i]  = results[i].fitA.valid ? results[i].fitA.mu  : 0.;
        muC_v[i]  = results[i].fitC.valid ? results[i].fitC.mu  : 0.;
        sigA_v[i] = results[i].fitA.valid ? results[i].fitA.sig : 0.;
        sigC_v[i] = results[i].fitC.valid ? results[i].fitC.sig : 0.;
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

    auto makeFrame = [&](TPad* pad, const char* ytitle,
                         const std::vector<double>& vA,
                         const std::vector<double>& vC) -> TH1D* {
        pad->cd();
        pad->SetLeftMargin(0.13); pad->SetRightMargin(0.03);
        pad->SetBottomMargin(0.22); pad->SetTopMargin(0.07);
        double lo = 1e9, hi = -1e9;
        for (int i = 0; i < nRuns; ++i) {
            if (vA[i] != 0.) { lo = std::min(lo, vA[i]); hi = std::max(hi, vA[i]); }
            if (vC[i] != 0.) { lo = std::min(lo, vC[i]); hi = std::max(hi, vC[i]); }
        }
        double pad_y = 0.20 * (hi - lo);
        lo -= pad_y; hi += pad_y;
        TH1D* fr = new TH1D(Form("fr_%s", ytitle), Form(";Run;%s", ytitle),
                             nRuns, 0.5, nRuns + 0.5);
        for (int i = 0; i < nRuns; ++i)
            fr->GetXaxis()->SetBinLabel(i+1, Form("%d", results[i].run));
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

    TCanvas* cvS = new TCanvas("cv_mu_sigma", "", 1600, 600);
    TPad* pL = new TPad("pL", "", 0.00, 0.00, 0.50, 1.00);
    TPad* pR = new TPad("pR", "", 0.50, 0.00, 1.00, 1.00);
    pL->Draw(); pR->Draw();

    TH1D* frMu  = makeFrame(pL, "#mu [ADC]",   muA_v,  muC_v);
    gMuA->Draw("P same"); gMuC->Draw("P same");
    pL->cd();
    TLegend* legL = new TLegend(0.15, 0.78, 0.42, 0.91);
    legL->SetBorderSize(0); legL->SetFillStyle(0); legL->SetTextSize(0.045);
    legL->AddEntry(gMuA, "Side A", "p"); legL->AddEntry(gMuC, "Side C", "p");
    legL->Draw();
    TLatex tlL; tlL.SetNDC(); tlL.SetTextSize(0.045); tlL.SetTextAlign(11);
    tlL.DrawLatex(0.15, 0.93, "Pb+Pb 2025  ZDC preamp Gaussian fit  #mu per run");

    TH1D* frSig = makeFrame(pR, "#sigma [ADC]", sigA_v, sigC_v);
    gSigA->Draw("P same"); gSigC->Draw("P same");
    pR->cd();
    TLegend* legR = new TLegend(0.15, 0.78, 0.42, 0.91);
    legR->SetBorderSize(0); legR->SetFillStyle(0); legR->SetTextSize(0.045);
    legR->AddEntry(gSigA, "Side A", "p"); legR->AddEntry(gSigC, "Side C", "p");
    legR->Draw();
    TLatex tlR; tlR.SetNDC(); tlR.SetTextSize(0.045); tlR.SetTextAlign(11);
    tlR.DrawLatex(0.15, 0.93, "Pb+Pb 2025  ZDC preamp Gaussian fit  #sigma per run");

    std::string outS = kPlotDir + "zdc_preamp_per_run_mu_sigma.png";
    cvS->SaveAs(outS.c_str());
    std::cout << "Saved: " << outS << std::endl;
    delete frMu; delete frSig;
}

// ---- Main ----
void plot_zdc_preamp_per_run_pbpb25() {
    gStyle->SetOptStat(0);
    gSystem->mkdir(kPlotDir.c_str(), true);

    // Load cuts
    const std::string cpath = PbPbEvSelCutsPath(2025);
    TFile* fc = TFile::Open(cpath.c_str(), "READ");
    if (!fc || fc->IsZombie()) throw std::runtime_error("Cuts file not found: " + cpath);
    TGraph* g_cut1 = (TGraph*)((TGraph*)fc->Get(PbPbEvSelKey::kZDCFCalCut))->Clone();
    double cut2_ns = ((TParameter<double>*)fc->Get(PbPbEvSelKey::kZDCTimeCutNs))->GetVal();
    double cut_A   = ((TParameter<double>*)fc->Get(PbPbEvSelKey::kPreampACutADC))->GetVal();
    double cut_C   = ((TParameter<double>*)fc->Get(PbPbEvSelKey::kPreampCCutADC))->GetVal();
    fc->Close();

    // Build TChain
    TChain chain("HeavyIonD3PD", "HeavyIonD3PD");
    for (int p = 1; p <= 6; ++p) {
        std::string f = kBase + "pbpb_2025/data_pbpb25_part" + std::to_string(p) + ".root";
        if (!gSystem->AccessPathName(f.c_str())) chain.Add(f.c_str());
        else std::cerr << "Skipping: " << f << std::endl;
    }
    chain.SetMakeClass(1);
    chain.SetBranchStatus("*", 0);
    Int_t   b_RunNum = 0, b_HLT = 0;
    Float_t FCal_P = 0.f, FCal_N = 0.f;
    Float_t zdcE[2]{}, zdcT[2]{}, preamp[2][4]{};
    chain.SetBranchStatus("RunNumber",                1);
    chain.SetBranchStatus("b_HLT_mu4_L1MU3V",         1);
    chain.SetBranchStatus("FCal_Et_P",                1);
    chain.SetBranchStatus("FCal_Et_N",                1);
    chain.SetBranchStatus("zdc_ZdcEnergy",            1);
    chain.SetBranchStatus("zdc_ZdcTime",              1);
    chain.SetBranchStatus("zdc_ZdcModulePreSampleAmp",1);
    chain.SetBranchAddress("RunNumber",                &b_RunNum);
    chain.SetBranchAddress("b_HLT_mu4_L1MU3V",         &b_HLT);
    chain.SetBranchAddress("FCal_Et_P",                &FCal_P);
    chain.SetBranchAddress("FCal_Et_N",                &FCal_N);
    chain.SetBranchAddress("zdc_ZdcEnergy",            zdcE);
    chain.SetBranchAddress("zdc_ZdcTime",              zdcT);
    chain.SetBranchAddress("zdc_ZdcModulePreSampleAmp", preamp);

    // Event loop
    std::map<int, TH1D*> hA_map, hC_map;
    const Long64_t n = chain.GetEntries();
    std::cout << "pbpb2025: " << n << " events, scanning..." << std::flush;
    long long n_sel = 0;
    for (Long64_t i = 0; i < n; ++i) {
        chain.GetEntry(i);
        if (!b_HLT) continue;
        const float fcal_AC = (FCal_P + FCal_N) * 1e-6f;
        const float zdcTot  = (zdcE[0] + zdcE[1]) / 1000.f;
        if (zdcTot > (float)PbPbEvSelEvalCut(g_cut1, fcal_AC)) continue;
        if (std::abs(zdcT[1]) >= (float)cut2_ns) continue;
        if (std::abs(zdcT[0]) >= (float)cut2_ns) continue;
        const int run = b_RunNum;
        if (hA_map.find(run) == hA_map.end()) {
            hA_map[run] = new TH1D(Form("hA_%d", run), "", 198, -1000., 2960.);
            hC_map[run] = new TH1D(Form("hC_%d", run), "", 198, -1000., 2960.);
            hA_map[run]->SetDirectory(nullptr);
            hC_map[run]->SetDirectory(nullptr);
        }
        float pA = 0.f, pC = 0.f;
        for (int k = 0; k < 4; ++k) pA += preamp[1][k];
        for (int k = 0; k < 4; ++k) pC += preamp[0][k];
        hA_map[run]->Fill(pA);
        hC_map[run]->Fill(pC);
        ++n_sel;
    }
    delete g_cut1;
    std::cout << "  " << n_sel << " pass trigger+Cut1+Cut2\n";

    // Sort runs
    std::vector<int> runs;
    for (auto& kv : hA_map) runs.push_back(kv.first);
    std::sort(runs.begin(), runs.end());
    std::cout << "Found " << runs.size() << " unique runs\n";

    // Rebin and fit
    std::vector<RunResult> results;
    results.reserve(runs.size());
    for (int run : runs) {
        TH1D* hA = hA_map[run];
        TH1D* hC = hC_map[run];
        hA->Rebin(3); hA->Scale(1.0/3.);
        hC->Rebin(3); hC->Scale(1.0/3.);
        GaussFit gfA = DoGaussFit(hA, Form("gA_%d", run));
        GaussFit gfC = DoGaussFit(hC, Form("gC_%d", run));
        results.push_back({run, gfA, gfC, hA, hC});
    }

    // Draw distribution versions
    DrawDistPages(results, cut_A, cut_C, kHardCut, "dists");
    DrawDistPages(results, cut_A, cut_C, kMu7Sig,  "dists_mu7sig");

    // Draw mu/sigma summary
    DrawMuSigma(results);

    for (auto& kv : hA_map) delete kv.second;
    for (auto& kv : hC_map) delete kv.second;
}
