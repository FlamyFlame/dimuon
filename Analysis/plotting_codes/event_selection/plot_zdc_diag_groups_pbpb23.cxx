// plot_zdc_diag_groups_pbpb23.cxx
//
// Diagnostic ZDC plots for 4 run groups in pbpb2023.
// Groups: Run 461674 / Run 462964 / Run 461641 / all other runs.
// Selection: b_HLT_mu4_L1MU3V only; no event-selection cuts.
//
// Five canvases (one per variable / debug):
//   zdc_diag_zdce_logx.png          — ZDC energy E/E_1N^truth log-x A/C overlay
//   zdc_diag_time_2d.png            — 2D ZDC time A vs C
//   zdc_diag_preamp.png             — ZDC preamp A/C overlay (x: 0–5000 ADC)
//   zdc_diag_zdce_vs_fcal_2d.png    — 2D ZDC energy vs FCal ET (linear xy, log z)
//   zdc_diag_time_461674_debug.png  — 2×2 ZDC time A/C with Gaussian fits:
//                                     top = run 461674, bottom = all other runs
//
// Usage (run from Analysis/):
//   source /usatlas/u/yuhanguo/setup.sh
//   root -l -b -q 'plotting_codes/event_selection/plot_zdc_diag_groups_pbpb23.cxx+'

#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TSystem.h"
#include "TStyle.h"

static const std::string kBase    = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/";
static const std::string kPlotDir = kBase + "plots/single_b_analysis/event_selection/pbpb_2023/run_dependence/";

static const int  kNGrp       = 4;
static const int  kSpecRuns[] = {461674, 462964, 461641};
static const char* kGrpLabel[] = {"Run 461674", "Run 462964", "Run 461641", "Other runs"};
static const int  kRefGrp     = 3;  // "all other runs" used for range determination

// 1-neutron truth energy at 5.36 TeV/nucleon beam
static const double kE1N = 2680.;   // GeV

// ---- Gaussian fit result ----
struct GaussFit { double mu, sig, amp; bool valid; };

// Two-pass Gaussian fit (seed window then refine within mu ± nLo/nHi*sigma)
static GaussFit DoGaussFit(TH1D* h, double seedLo, double seedHi,
                            double nLo, double nHi, const char* fname) {
    if (h->GetEntries() < 30) return {0., 0., 0., false};
    int b1 = h->FindBin(seedLo + 0.5*h->GetBinWidth(1));
    int b2 = h->FindBin(seedHi - 0.5*h->GetBinWidth(1));
    double localMax = 0.; int bMax = b1;
    for (int b = b1; b <= b2; ++b)
        if (h->GetBinContent(b) > localMax) { localMax = h->GetBinContent(b); bMax = b; }
    if (localMax < 1.) return {0., 0., 0., false};

    TF1 g1(Form("%s_seed", fname), "gaus", seedLo, seedHi);
    g1.SetParameters(localMax, h->GetBinCenter(bMax), 0.25*(seedHi-seedLo));
    h->Fit(&g1, "RQN");
    double mu1  = g1.GetParameter(1);
    double sig1 = std::abs(g1.GetParameter(2));
    if (mu1 < seedLo || mu1 > seedHi)         return {0., 0., 0., false};
    if (sig1 < 1e-9 || sig1 > (seedHi-seedLo)) return {0., 0., 0., false};

    TF1* gf = new TF1(fname, "gaus", mu1 - nLo*sig1, mu1 + nHi*sig1);
    gf->SetParameters(g1.GetParameter(0), mu1, sig1);
    int st = h->Fit(gf, "RQN");
    double mu_f  = gf->GetParameter(1);
    double sig_f = std::abs(gf->GetParameter(2));
    double amp_f = gf->GetParameter(0);
    delete gf;
    if (st != 0 || sig_f < 1e-9)              return {0., 0., 0., false};
    if (mu_f < seedLo || mu_f > seedHi)        return {0., 0., 0., false};
    if (sig_f > (seedHi - seedLo))             return {0., 0., 0., false};
    return {mu_f, sig_f, amp_f, true};
}

// Draw Gaussian fit on top of histogram and print mu/sigma via TLatex
static void FitAndDrawTime(TH1D* h, const char* fname, int fitColor = kBlue+1) {
    GaussFit gf = DoGaussFit(h, -5., 5., 2.5, 2.5, fname);
    if (!gf.valid) return;
    TF1* f = new TF1(Form("draw_%s", fname), "gaus",
                     gf.mu - 2.5*gf.sig, gf.mu + 2.5*gf.sig);
    f->SetParameters(gf.amp, gf.mu, gf.sig);
    f->SetLineColor(fitColor); f->SetLineWidth(2);
    f->Draw("SAME");
    TLatex tl; tl.SetNDC(); tl.SetTextSize(0.054); tl.SetTextAlign(11);
    tl.DrawLatex(0.56, 0.80, Form("#mu = %.3f ns",  gf.mu));
    tl.DrawLatex(0.56, 0.72, Form("#sigma = %.3f ns", gf.sig));
}

// ---- Log-spaced bin edges helper ----
static void LogBins(int n, double lo, double hi, double* edges) {
    double logLo = std::log10(lo), logHi = std::log10(hi);
    for (int i = 0; i <= n; ++i)
        edges[i] = std::pow(10., logLo + i*(logHi-logLo)/n);
}

// ---- Effective 1D range from first/last non-empty bin, clamped to histogram limits ----
static std::pair<double,double> EffRange1D(TH1D* h, bool logScale = false, double pad = 0.05) {
    int first = h->FindFirstBinAbove(0);
    int last  = h->FindLastBinAbove(0);
    const double hmin = h->GetXaxis()->GetXmin(), hmax = h->GetXaxis()->GetXmax();
    if (first <= 0) return {hmin, hmax};
    double lo = h->GetXaxis()->GetBinLowEdge(first);
    double hi = h->GetXaxis()->GetBinUpEdge(last);
    if (logScale) {
        if (lo <= 0.) lo = hmin;
        double logLo = std::log10(lo), logHi = std::log10(hi);
        double span  = logHi - logLo;
        lo = std::max(hmin, std::pow(10., logLo - pad*span));
        hi = std::min(hmax, std::pow(10., logHi + pad*span));
    } else {
        double span = hi - lo;
        lo = std::max(hmin, lo - pad*span);
        hi = std::min(hmax, hi + pad*span);
    }
    return {lo, hi};
}

// ---- Effective 2D axis ranges from projection ----
static std::pair<double,double> EffRange2D(TH2D* h, bool isX, double pad = 0.05) {
    TH1D* proj = isX ? h->ProjectionX("__px_tmp") : h->ProjectionY("__py_tmp");
    auto r = EffRange1D(proj, false, pad);
    delete proj;
    return r;
}

// ---- Set up a 2×2 canvas and return the 4 pads ----
static TCanvas* Make2x2Canvas(const char* name, TPad** pads) {
    TCanvas* cv = new TCanvas(name, "", 1600, 1200);
    pads[0] = new TPad(Form("%s_p0",name), "", 0.00, 0.50, 0.50, 1.00);
    pads[1] = new TPad(Form("%s_p1",name), "", 0.50, 0.50, 1.00, 1.00);
    pads[2] = new TPad(Form("%s_p2",name), "", 0.00, 0.00, 0.50, 0.50);
    pads[3] = new TPad(Form("%s_p3",name), "", 0.50, 0.00, 1.00, 0.50);
    for (int i = 0; i < 4; ++i) pads[i]->Draw();
    return cv;
}

static void SetPadStyle(TPad* p) {
    p->SetLeftMargin(0.14); p->SetRightMargin(0.04);
    p->SetTopMargin(0.08);  p->SetBottomMargin(0.14);
}

static void LabelPad(int grp, long long n) {
    TLatex tl; tl.SetNDC(); tl.SetTextSize(0.055); tl.SetTextAlign(11);
    tl.DrawLatex(0.16, 0.91, kGrpLabel[grp]);
    tl.SetTextSize(0.042);
    tl.DrawLatex(0.16, 0.84, Form("N = %lld", n));
}

// ============================================================
void plot_zdc_diag_groups_pbpb23() {
    gStyle->SetOptStat(0);
    gSystem->mkdir(kPlotDir.c_str(), true);

    // ---- Log-spaced bins for ZDC energy ratio E/E_1N (0.3 – 100) ----
    const int kNE = 80;
    double eBins[kNE+1];
    LogBins(kNE, 0.3, 100., eBins);

    // ---- Histograms [group] ----
    TH1D* hEA[kNGrp], *hEC[kNGrp];
    TH2D* hTAC[kNGrp];
    TH1D* hTA[kNGrp], *hTC[kNGrp];   // 1D ZDC time for Gaussian fit debug
    TH1D* hPA[kNGrp], *hPC[kNGrp];
    TH2D* hEF[kNGrp];
    long long nEvt[kNGrp] = {};

    for (int g = 0; g < kNGrp; ++g) {
        hEA[g]  = new TH1D(Form("hEA%d",g),
                            ";E_{ZDC}^{A}/E_{1N}^{truth};Events", kNE, eBins);
        hEC[g]  = new TH1D(Form("hEC%d",g),
                            ";E_{ZDC}^{C}/E_{1N}^{truth};Events", kNE, eBins);
        hTAC[g] = new TH2D(Form("hTAC%d",g),";ZDC t^{A} [ns];ZDC t^{C} [ns]",
                            100,-10.,10., 100,-10.,10.);
        hTA[g]  = new TH1D(Form("hTA%d",g),  ";ZDC t^{A} [ns];Events", 200,-10.,10.);
        hTC[g]  = new TH1D(Form("hTC%d",g),  ";ZDC t^{C} [ns];Events", 200,-10.,10.);
        // preamp: fixed range 0–5000 ADC
        hPA[g]  = new TH1D(Form("hPA%d",g),  ";PreSampleAmp sum [ADC];Events", 100, 0., 5000.);
        hPC[g]  = new TH1D(Form("hPC%d",g),  ";PreSampleAmp sum [ADC];Events", 100, 0., 5000.);
        // ZDC vs FCal 2D: linear y-axis, 100 linear bins 0–400 TeV ZDC total
        hEF[g]  = new TH2D(Form("hEF%d",g),
                            ";FCal E_{T}^{A+C} [TeV];ZDC E_{total} [TeV]",
                            120,-0.5,5.5, 100,0.,400.);
        hEA[g]->SetDirectory(nullptr); hEC[g]->SetDirectory(nullptr);
        hTAC[g]->SetDirectory(nullptr);
        hTA[g]->SetDirectory(nullptr);  hTC[g]->SetDirectory(nullptr);
        hPA[g]->SetDirectory(nullptr);  hPC[g]->SetDirectory(nullptr);
        hEF[g]->SetDirectory(nullptr);
    }

    // ---- Build TChain ----
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
    Float_t FCal_P  = 0.f, FCal_N = 0.f;
    Float_t zdcE[2]{}, zdcT[2]{}, preamp[2][4]{};

    chain.SetBranchStatus("RunNumber",                 1);
    chain.SetBranchStatus("b_HLT_mu4_L1MU3V",          1);
    chain.SetBranchStatus("FCal_Et_P",                 1);
    chain.SetBranchStatus("FCal_Et_N",                 1);
    chain.SetBranchStatus("zdc_ZdcEnergy",             1);
    chain.SetBranchStatus("zdc_ZdcTime",               1);
    chain.SetBranchStatus("zdc_ZdcModulePreSampleAmp", 1);
    chain.SetBranchAddress("RunNumber",                 &b_RunNum);
    chain.SetBranchAddress("b_HLT_mu4_L1MU3V",          &b_HLT);
    chain.SetBranchAddress("FCal_Et_P",                 &FCal_P);
    chain.SetBranchAddress("FCal_Et_N",                 &FCal_N);
    chain.SetBranchAddress("zdc_ZdcEnergy",             zdcE);
    chain.SetBranchAddress("zdc_ZdcTime",               zdcT);
    chain.SetBranchAddress("zdc_ZdcModulePreSampleAmp", preamp);

    // ---- Event loop ----
    const Long64_t nTot = chain.GetEntries();
    std::cout << "Scanning " << nTot << " events..." << std::flush;
    for (Long64_t i = 0; i < nTot; ++i) {
        chain.GetEntry(i);
        if (!b_HLT) continue;

        int grp = kRefGrp;
        for (int k = 0; k < 3; ++k)
            if (b_RunNum == kSpecRuns[k]) { grp = k; break; }

        ++nEvt[grp];
        const float fcal    = (FCal_P + FCal_N) * 1e-6f;        // TeV
        const float zdcTot  = (zdcE[0] + zdcE[1]) / 1000.f;    // TeV
        float pA = 0.f, pC = 0.f;
        for (int k = 0; k < 4; ++k) pA += preamp[1][k];
        for (int k = 0; k < 4; ++k) pC += preamp[0][k];

        // zdcE[1]=A, zdcE[0]=C; fill energy as ratio to 1N truth energy
        hEA[grp]->Fill(zdcE[1] / kE1N);
        hEC[grp]->Fill(zdcE[0] / kE1N);
        hTAC[grp]->Fill(zdcT[1], zdcT[0]);
        hTA[grp]->Fill(zdcT[1]);
        hTC[grp]->Fill(zdcT[0]);
        hPA[grp]->Fill(pA);
        hPC[grp]->Fill(pC);
        hEF[grp]->Fill(fcal, zdcTot);
    }
    std::cout << "  done.\n";
    for (int g = 0; g < kNGrp; ++g)
        std::cout << "  " << kGrpLabel[g] << ": " << nEvt[g] << " events\n";

    // ---- Merged ZDC time histograms for "all other runs" (grp 1 + 2 + 3) ----
    TH1D* hTA_other = (TH1D*)hTA[1]->Clone("hTA_other");
    TH1D* hTC_other = (TH1D*)hTC[1]->Clone("hTC_other");
    hTA_other->Add(hTA[2]); hTA_other->Add(hTA[3]);
    hTC_other->Add(hTC[2]); hTC_other->Add(hTC[3]);
    hTA_other->SetDirectory(nullptr);
    hTC_other->SetDirectory(nullptr);
    long long nOther = nEvt[1] + nEvt[2] + nEvt[3];

    // ---- Determine display ranges ----
    // Energy: fixed ratio range 0.3–100
    const double eLo = 0.3, eHi = 100.;

    // Time 2D
    auto [txLo, txHi] = EffRange2D(hTAC[kRefGrp], true);
    auto [tyLo, tyHi] = EffRange2D(hTAC[kRefGrp], false);

    // FCal 2D: fixed ranges matching event selection plots
    // x: [-0.5, 5.5] TeV (histogram range); y: [0, 400] TeV ZDC total

    // ---- Helper lambdas ----
    auto Style2D = [](TH2D* h) {
        h->SetContour(99);
        h->GetXaxis()->SetTitleSize(0.055); h->GetXaxis()->SetLabelSize(0.045);
        h->GetYaxis()->SetTitleSize(0.055); h->GetYaxis()->SetLabelSize(0.045);
        h->GetYaxis()->SetTitleOffset(1.20);
    };
    auto Style1D = [](TH1D* h, int color) {
        h->SetLineColor(color); h->SetLineWidth(1);
        h->SetMarkerColor(color); h->SetMarkerStyle(20); h->SetMarkerSize(0.3);
        h->GetXaxis()->SetTitleSize(0.055); h->GetXaxis()->SetLabelSize(0.045);
        h->GetYaxis()->SetTitleSize(0.050); h->GetYaxis()->SetLabelSize(0.045);
        h->GetYaxis()->SetTitleOffset(1.25);
    };

    // ============================================================
    // Canvas 1: ZDC energy E/E_1N^truth log-x, A/C overlay
    {
        TPad* pads[4];
        TCanvas* cv = Make2x2Canvas("cv_zdce", pads);
        for (int g = 0; g < kNGrp; ++g) {
            pads[g]->cd();
            SetPadStyle(pads[g]);
            pads[g]->SetLogx(); pads[g]->SetLogy();
            Style1D(hEA[g], kBlack); Style1D(hEC[g], kRed+1);
            hEA[g]->GetXaxis()->SetRangeUser(eLo, eHi);
            double yhi = std::max(hEA[kRefGrp]->GetMaximum(), hEC[kRefGrp]->GetMaximum()) * 5.;
            if (yhi < 1.) yhi = 10.;
            hEA[g]->SetMinimum(0.5); hEA[g]->SetMaximum(yhi);
            hEA[g]->Draw("E"); hEC[g]->Draw("E SAME");
            if (g == 0) {
                TLegend* leg = new TLegend(0.55, 0.78, 0.82, 0.91);
                leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.048);
                leg->AddEntry(hEA[g], "Side A", "l");
                leg->AddEntry(hEC[g], "Side C", "l");
                leg->Draw();
            }
            LabelPad(g, nEvt[g]);
        }
        std::string out = kPlotDir + "zdc_diag_zdce_logx.png";
        cv->SaveAs(out.c_str());
        std::cout << "Saved: " << out << std::endl;
        delete cv;
    }

    // ============================================================
    // Canvas 2: ZDC time 2D (A vs C)
    {
        TPad* pads[4];
        TCanvas* cv = Make2x2Canvas("cv_zdct", pads);
        for (int g = 0; g < kNGrp; ++g) {
            pads[g]->cd();
            SetPadStyle(pads[g]);
            pads[g]->SetRightMargin(0.14);
            pads[g]->SetLogz();
            Style2D(hTAC[g]);
            hTAC[g]->GetXaxis()->SetRangeUser(txLo, txHi);
            hTAC[g]->GetYaxis()->SetRangeUser(tyLo, tyHi);
            hTAC[g]->SetMinimum(0.5);
            hTAC[g]->Draw("COLZ");
            LabelPad(g, nEvt[g]);
        }
        std::string out = kPlotDir + "zdc_diag_time_2d.png";
        cv->SaveAs(out.c_str());
        std::cout << "Saved: " << out << std::endl;
        delete cv;
    }

    // ============================================================
    // Canvas 3: ZDC preamp A/C overlay (x fixed 0–5000 ADC)
    {
        TPad* pads[4];
        TCanvas* cv = Make2x2Canvas("cv_preamp", pads);
        for (int g = 0; g < kNGrp; ++g) {
            pads[g]->cd();
            SetPadStyle(pads[g]);
            pads[g]->SetLogy();
            Style1D(hPA[g], kBlack); Style1D(hPC[g], kRed+1);
            double yhi = std::max(hPA[kRefGrp]->GetMaximum(), hPC[kRefGrp]->GetMaximum()) * 5.;
            if (yhi < 1.) yhi = 10.;
            hPA[g]->SetMinimum(0.5); hPA[g]->SetMaximum(yhi);
            hPA[g]->Draw("E"); hPC[g]->Draw("E SAME");
            if (g == 0) {
                TLegend* leg = new TLegend(0.55, 0.78, 0.82, 0.91);
                leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.048);
                leg->AddEntry(hPA[g], "Side A", "l");
                leg->AddEntry(hPC[g], "Side C", "l");
                leg->Draw();
            }
            LabelPad(g, nEvt[g]);
        }
        std::string out = kPlotDir + "zdc_diag_preamp.png";
        cv->SaveAs(out.c_str());
        std::cout << "Saved: " << out << std::endl;
        delete cv;
    }

    // ============================================================
    // Canvas 4: ZDC energy total vs FCal ET 2D (linear x, y, z)
    {
        TPad* pads[4];
        TCanvas* cv = Make2x2Canvas("cv_efcal", pads);
        for (int g = 0; g < kNGrp; ++g) {
            pads[g]->cd();
            SetPadStyle(pads[g]);
            pads[g]->SetRightMargin(0.14);
            pads[g]->SetLogz();
            Style2D(hEF[g]);
            hEF[g]->SetMinimum(0.5);
            hEF[g]->Draw("COLZ");
            LabelPad(g, nEvt[g]);
        }
        std::string out = kPlotDir + "zdc_diag_zdce_vs_fcal_2d.png";
        cv->SaveAs(out.c_str());
        std::cout << "Saved: " << out << std::endl;
        delete cv;
    }

    // ============================================================
    // Canvas 5: ZDC time debug — 2×2
    //   [0] top-left  : run 461674, side A, with Gaussian fit + mu/sigma
    //   [1] top-right : run 461674, side C, with Gaussian fit + mu/sigma
    //   [2] bot-left  : all other runs, side A, with Gaussian fit + mu/sigma
    //   [3] bot-right : all other runs, side C, with Gaussian fit + mu/sigma
    {
        TPad* pads[4];
        TCanvas* cv = Make2x2Canvas("cv_tdbg", pads);

        // Helper to draw one ZDC time pad — y scale per-histogram, x scale fixed by histogram range
        auto DrawTimePad = [&](TPad* pad, TH1D* h, const char* label,
                                long long n, const char* fitName, int grpIdx) {
            pad->cd();
            SetPadStyle(pad);
            Style1D(h, kBlack);
            double yhi = h->GetMaximum() * 1.5;
            if (yhi < 1.) yhi = 10.;
            h->SetMaximum(yhi);
            h->SetMinimum(0.);
            h->Draw("E");
            FitAndDrawTime(h, fitName);
            TLatex tl; tl.SetNDC(); tl.SetTextSize(0.055); tl.SetTextAlign(11);
            tl.DrawLatex(0.16, 0.91, label);
            tl.SetTextSize(0.042);
            tl.DrawLatex(0.16, 0.84, Form("N = %lld", n));
        };

        DrawTimePad(pads[0], hTA[0], "Run 461674  Side A", nEvt[0], "gtA_461674_A", 0);
        DrawTimePad(pads[1], hTC[0], "Run 461674  Side C", nEvt[0], "gtA_461674_C", 0);
        DrawTimePad(pads[2], hTA_other, "Other runs  Side A",  nOther, "gtA_other_A", 3);
        DrawTimePad(pads[3], hTC_other, "Other runs  Side C",  nOther, "gtA_other_C", 3);

        std::string out = kPlotDir + "zdc_diag_time_461674_debug.png";
        cv->SaveAs(out.c_str());
        std::cout << "Saved: " << out << std::endl;
        delete cv;
    }

    for (int g = 0; g < kNGrp; ++g) {
        delete hEA[g]; delete hEC[g];
        delete hTAC[g];
        delete hTA[g];  delete hTC[g];
        delete hPA[g];  delete hPC[g];
        delete hEF[g];
    }
    delete hTA_other; delete hTC_other;
}
