// plot_pbpb_ZDC_preamp.cxx
// Overlay ZDC per-arm presample amplitude (sum of 4 modules) for A and C sides.
// Reads directly from the pair ntuple (sign1 + sign2 trees).
//
// Usage:
//   .L plot_pbpb_ZDC_preamp.cxx+
//   plot_pbpb_ZDC_preamp(24)

#include <string>
#include <stdexcept>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TF1.h"
#include "TMath.h"
#include "TPaveText.h"

void plot_pbpb_ZDC_preamp(int run_year = 24) {
    run_year %= 2000;
    const std::string yr = std::to_string(run_year);

    const std::string infile =
        "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_20" + yr +
        "/muon_pairs_pbpb_20" + yr + "_single_mu4_mindR_0_02_res_cut_v2.root";
    const std::string out_dir =
        "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/single_b_analysis/event_selection";

    if (gSystem->AccessPathName(infile.c_str()))
        throw std::runtime_error("Input file not found: " + infile);
    TFile* f = TFile::Open(infile.c_str(), "READ");
    if (!f || f->IsZombie())
        throw std::runtime_error("Cannot open: " + infile);

    gSystem->mkdir(out_dir.c_str(), true);

    // 150 bins, -1000 to 6500
    TH1D* hA = new TH1D("hA", ";ZDC preamp sum [ADC counts];Pairs", 150, -1000, 6500);
    TH1D* hC = new TH1D("hC", ";ZDC preamp sum [ADC counts];Pairs", 150, -1000, 6500);
    hA->SetDirectory(nullptr);
    hC->SetDirectory(nullptr);

    for (const char* tname : {"muon_pair_tree_sign1", "muon_pair_tree_sign2"}) {
        TTree* t = (TTree*)f->Get(tname);
        if (!t) { std::cerr << tname << " not found, skipping" << std::endl; continue; }
        t->SetMakeClass(1);
        Float_t pa = -999.f, pc = -999.f;
        t->SetBranchAddress("MuonPairObj.ZDC_preamp_A", &pa);
        t->SetBranchAddress("MuonPairObj.ZDC_preamp_C", &pc);
        const Long64_t n = t->GetEntries();
        for (Long64_t i = 0; i < n; ++i) {
            t->GetEntry(i);
            hA->Fill(pa);
            hC->Fill(pc);
        }
    }
    f->Close();

    std::cout << "preamp_A entries=" << hA->GetEntries()
              << "  preamp_C entries=" << hC->GetEntries() << std::endl;

    // No normalisation — raw pair counts

    hA->SetLineColor(kBlue);
    hA->SetLineWidth(2);
    hA->SetMarkerColor(kBlue);
    hA->SetMarkerStyle(20);
    hA->SetMarkerSize(0.5);

    hC->SetLineColor(kRed);
    hC->SetLineWidth(2);
    hC->SetMarkerColor(kRed);
    hC->SetMarkerStyle(20);
    hC->SetMarkerSize(0.5);

    gStyle->SetOptStat(0);

    TCanvas* c = new TCanvas("c_preamp", "", 800, 620);
    c->SetLeftMargin(0.13);
    c->SetBottomMargin(0.13);
    c->SetLogy();

    double ymax = 3.0 * std::max(hA->GetMaximum(), hC->GetMaximum());
    hA->GetYaxis()->SetRangeUser(0.5, ymax);
    hA->Draw("E");
    hC->Draw("E SAME");

    TLegend* leg = new TLegend(0.55, 0.72, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.036);
    leg->AddEntry(hA, "ZDC preamp A-side", "lp");
    leg->AddEntry(hC, "ZDC preamp C-side", "lp");
    leg->Draw();

    TLatex tl; tl.SetNDC(); tl.SetTextSize(0.036);
    tl.DrawLatex(0.15, 0.87, Form("Pb+Pb 20%s data", yr.c_str()));

    const std::string outpath = out_dir + "/ZDC_preamp_A_vs_C_pbpb_20" + yr + ".png";
    c->SaveAs(outpath.c_str());
    std::cout << "Saved: " << outpath << std::endl;
    delete c;
    delete hA;
    delete hC;
}

// Two-pass Gaussian fit of the ZDC preamp peak near -200 ADC.
// Saves full-range and zoomed (-1000 to 1000) canvases.
void plot_pbpb_ZDC_preamp_fit(int run_year = 24) {
    run_year %= 2000;
    const std::string yr = std::to_string(run_year);

    const std::string infile =
        "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_20" + yr +
        "/muon_pairs_pbpb_20" + yr + "_single_mu4_mindR_0_02_res_cut_v2.root";
    const std::string out_dir =
        "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/single_b_analysis/event_selection";

    if (gSystem->AccessPathName(infile.c_str()))
        throw std::runtime_error("Input file not found: " + infile);
    TFile* f = TFile::Open(infile.c_str(), "READ");
    if (!f || f->IsZombie())
        throw std::runtime_error("Cannot open: " + infile);

    gSystem->mkdir(out_dir.c_str(), true);

    // 300 bins from -1000 to 6500 = 25 ADC/bin for better peak resolution
    TH1D* hA = new TH1D("hAf", ";ZDC preamp sum [ADC counts];Pairs", 300, -1000, 6500);
    TH1D* hC = new TH1D("hCf", ";ZDC preamp sum [ADC counts];Pairs", 300, -1000, 6500);
    hA->SetDirectory(nullptr);
    hC->SetDirectory(nullptr);

    for (const char* tname : {"muon_pair_tree_sign1", "muon_pair_tree_sign2"}) {
        TTree* t = (TTree*)f->Get(tname);
        if (!t) { std::cerr << tname << " not found, skipping" << std::endl; continue; }
        t->SetMakeClass(1);
        Float_t pa = -999.f, pc = -999.f;
        t->SetBranchAddress("MuonPairObj.ZDC_preamp_A", &pa);
        t->SetBranchAddress("MuonPairObj.ZDC_preamp_C", &pc);
        const Long64_t n = t->GetEntries();
        for (Long64_t i = 0; i < n; ++i) {
            t->GetEntry(i);
            hA->Fill(pa);
            hC->Fill(pc);
        }
    }
    f->Close();

    for (TH1D* h : {hA, hC}) {
        h->SetMarkerStyle(20);
        h->SetMarkerSize(0.6);
        h->SetMarkerColor(kBlack);
        h->SetLineColor(kBlack);
        h->SetLineWidth(1);
    }

    gStyle->SetOptStat(0);

    // Pass 1: Gaussian seed; Pass 2: EMG fit. Returns heap-allocated TF1 (caller owns).
    // EMG = Gaussian convolved with one-sided exponential (right tail).
    //   f(x) = A * exp(lam*(mu + 0.5*lam*sig^2 - x)) * erfc((mu + lam*sig^2 - x)/(sig*sqrt2))
    // Also prints chi2/NDF comparison vs pure Gaussian in window [-400, 400].
    auto doFit = [](TH1D* h, const char* name) -> TF1* {
        // Pass 1: Gaussian seed, range [-700, +350], initial mean -150
        std::string p1name = std::string(name) + "_p1";
        TF1 g1(p1name.c_str(), "gaus", -700., 350.);
        double amp0 = h->GetBinContent(h->FindBin(-150.));
        if (amp0 <= 0.) amp0 = h->GetMaximum();
        g1.SetParameters(amp0, -150., 200.);
        h->Fit(&g1, "RQN");
        double mu1  = g1.GetParameter(1);
        double sig1 = std::abs(g1.GetParameter(2));

        double fit_lo = mu1 - 3.*sig1;
        double fit_hi = 400.;

        // Gaussian comparison in same range (for chi2 benchmark)
        std::string gcname = std::string(name) + "_gc";
        TF1 gc(gcname.c_str(), "gaus", fit_lo, fit_hi);
        gc.SetParameters(g1.GetParameter(0), mu1, sig1);
        h->Fit(&gc, "RQN");

        // Pass 2: EMG fit
        TF1* emg = new TF1(name, [](double* x, double* p) -> double {
            double A   = p[0];
            double mu  = p[1];
            double sig = std::abs(p[2]);
            double lam = std::abs(p[3]);
            return A * TMath::Exp(lam*(mu + 0.5*lam*sig*sig - x[0]))
                     * TMath::Erfc((mu + lam*sig*sig - x[0]) / (sig * TMath::Sqrt2()));
        }, fit_lo, fit_hi, 4);

        emg->SetParameters(amp0, mu1, sig1, 0.004);
        emg->SetParLimits(2, 5., 500.);
        emg->SetParLimits(3, 1e-5, 0.05);
        h->Fit(emg, "RQN");
        emg->SetLineColor(kRed);
        emg->SetLineWidth(2);

        // Chi2 comparison in [-400, 400]
        auto chi2_in_window = [&](TF1* f) {
            int b1 = h->FindBin(-399.5), b2 = h->FindBin(399.5);
            double chi2 = 0.; int ndof = 0;
            for (int b = b1; b <= b2; ++b) {
                double e = h->GetBinError(b);
                if (e <= 0.) continue;
                double res = h->GetBinContent(b) - f->Eval(h->GetBinCenter(b));
                chi2 += res*res / (e*e);
                ++ndof;
            }
            return std::make_pair(chi2, ndof);
        };

        auto rg  = chi2_in_window(&gc);
        auto re  = chi2_in_window(emg);
        int ndf_g   = rg.second - 3;
        int ndf_emg = re.second - 4;
        std::cout << Form("[%s] chi2/NDF in [-400,400]:  Gaussian = %.0f/%d = %.2f  |  EMG = %.0f/%d = %.2f\n",
                          name,
                          rg.first,  ndf_g,   rg.first  / ndf_g,
                          re.first,  ndf_emg, re.first  / ndf_emg);
        return emg;
    };

    TF1* fitA = doFit(hA, "fitA_preamp");
    TF1* fitC = doFit(hC, "fitC_preamp");

    // manual_y=true: set y-range from local max in [xlo,xhi]; false: use global max.
    // rebin: rebin factor applied to a clone before drawing (1 = no rebin).
    auto makeCanvas = [&](const char* cname, double xlo, double xhi,
                          bool manual_y, int rebin, double pt_x0,
                          const std::string& outpath) {
        TCanvas* c = new TCanvas(cname, "", 1200, 600);
        c->Divide(2, 1, 0.005, 0.005);

        TH1D* hArr[2] = {hA, hC};
        TF1*  fArr[2] = {fitA, fitC};
        const char* sides[2] = {"A-side", "C-side"};
        TH1D* clones[2] = {nullptr, nullptr};

        for (int i = 0; i < 2; ++i) {
            c->cd(i + 1);
            gPad->SetLeftMargin(0.14);
            gPad->SetBottomMargin(0.13);
            gPad->SetLogy();

            TH1D* hd = (rebin > 1)
                ? (TH1D*)hArr[i]->Clone(Form("hd_%s_%d", cname, i))
                : hArr[i];
            if (rebin > 1) {
                hd->Rebin(rebin);
                hd->Scale(1.0 / rebin);  // restore to counts per original bin width
                clones[i] = hd;
            }

            hd->GetYaxis()->UnZoom();
            hd->GetXaxis()->SetRangeUser(xlo, xhi);
            hd->GetXaxis()->SetTitle(
                Form("ZDC preamp sum side %s [ADC counts]", (i == 0 ? "A" : "C")));

            if (manual_y) {
                int b1 = hd->FindBin(xlo + 0.5);
                int b2 = hd->FindBin(xhi - 0.5);
                double local_max = 0.;
                for (int b = b1; b <= b2; ++b)
                    local_max = std::max(local_max, hd->GetBinContent(b));
                hd->GetYaxis()->SetRangeUser(0.5, local_max * 3.0);
            } else {
                hd->GetYaxis()->SetRangeUser(0.5, hd->GetMaximum() * 3.0);
            }

            // Equation box drawn first → behind histogram
            double mu_f  = fArr[i]->GetParameter(1);
            double sig_f = std::abs(fArr[i]->GetParameter(2));
            double lam_f = std::abs(fArr[i]->GetParameter(3));
            double c1    = mu_f + 0.5*lam_f*sig_f*sig_f;
            double c2    = mu_f + lam_f*sig_f*sig_f;
            double c3    = sig_f * TMath::Sqrt2();

            TPaveText* pt = new TPaveText(pt_x0, 0.63, pt_x0 + 0.47, 0.88, "NDC");
            pt->SetBorderSize(0);
            pt->SetFillStyle(0);   // transparent: histogram shows through
            pt->SetTextSize(0.031);
            pt->SetTextAlign(12);
            pt->AddText("f(x) = A e^{#lambda(c_{1}-x)} erfc((c_{2}-x)/c_{3})");
            pt->AddText(Form("c_{1} #equiv #mu+#lambda#sigma^{2}/2 = %.1f ADC", c1));
            pt->AddText(Form("c_{2} #equiv #mu+#lambda#sigma^{2}   = %.1f ADC", c2));
            pt->AddText(Form("c_{3} #equiv #sigma#sqrt{2}          = %.1f ADC", c3));
            pt->AddText(Form("#lambda = %.5f ADC^{-1}", lam_f));
            hd->Draw("E");
            fArr[i]->Draw("SAME");
            pt->Draw();

            TLatex tl; tl.SetNDC(); tl.SetTextSize(0.040);
            tl.DrawLatex(0.17, 0.88, Form("ZDC preamp %s", sides[i]));
            tl.SetTextSize(0.035);
            tl.DrawLatex(0.17, 0.82, Form("Pb+Pb 20%s data", yr.c_str()));
            tl.DrawLatex(0.17, 0.76, Form("#mu = %.1f ADC", fArr[i]->GetParameter(1)));
            tl.DrawLatex(0.17, 0.70, Form("#sigma = %.1f ADC",
                                           std::abs(fArr[i]->GetParameter(2))));
            tl.DrawLatex(0.17, 0.64, Form("1/#lambda = %.0f ADC",
                                           1./std::abs(fArr[i]->GetParameter(3))));

        }

        c->SaveAs(outpath.c_str());
        std::cout << "Saved: " << outpath << std::endl;
        delete c;
        for (TH1D* cl : clones) delete cl;
    };

    makeCanvas("c_preamp_fit_full", -1000., 6500., false, 1, 0.50,
               out_dir + "/ZDC_preamp_fit_pbpb_20" + yr + ".png");
    makeCanvas("c_preamp_fit_zoom", -800., 800., true, 3, 0.60,
               out_dir + "/ZDC_preamp_fit_zoom_pbpb_20" + yr + ".png");

    delete fitA;
    delete fitC;
    delete hA;
    delete hC;
}
