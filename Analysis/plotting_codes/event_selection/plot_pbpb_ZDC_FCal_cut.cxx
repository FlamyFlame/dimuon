// plot_pbpb_ZDC_FCal_cut.cxx
//
// Two-pass per-slice Gaussian fit on ZDC E_total vs FCal E_T^{A+C}.
// For each FCal Et slice:
//   Pass 1 — initial mean = Y bin with max content; fit in [mu0 ± 3*RMS]
//   Pass 2 — fit in [mu1 - 3*sigma1, mu1 + 3*sigma1]
// Cut per slice: mu2 + N_SIGMA_CUT * sigma2
//
// Canvas 1: 2D COLZ + red markers at cut value per slice
// Canvas 2: passing (ZDC <= cut) | failing (ZDC > cut) as COLZ subplots
//
// Usage:
//   .L plot_pbpb_ZDC_FCal_cut.cxx+
//   plot_pbpb_ZDC_FCal_cut(24)

#include <string>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLatex.h"

// ---- tunable parameters -----------------------------------------------------
static const int    N_BINS_PER_SLICE = 5;     // FCal Et x-bins grouped per slice
static const int    MIN_ENTRIES      = 100;    // skip slice if fewer entries
static const double N_SIGMA_CUT      = 5.0;   // cut threshold: mu + N*sigma

// ---- helpers ----------------------------------------------------------------
static TH2D* GetH2(TFile* f, const std::string& name) {
    TH2D* h = dynamic_cast<TH2D*>(f->Get(name.c_str()));
    if (!h) throw std::runtime_error("Histogram not found: " + name);
    h->SetDirectory(nullptr);
    return h;
}

// ---- main -------------------------------------------------------------------
void plot_pbpb_ZDC_FCal_cut(int run_year = 24) {
    run_year %= 2000;
    const std::string yr = std::to_string(run_year);

    const std::string infile =
        "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_20" + yr +
        "/histograms_real_pairs_pbpb_20" + yr + "_single_mu4_no_trg_plots_coarse_q_eta_bin.root";
    const std::string out_dir =
        "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/single_b_analysis/event_selection";

    if (gSystem->AccessPathName(infile.c_str()))
        throw std::runtime_error("Input file not found: " + infile);
    TFile* f = TFile::Open(infile.c_str(), "READ");
    if (!f || f->IsZombie())
        throw std::runtime_error("Cannot open: " + infile);

    gSystem->mkdir(out_dir.c_str(), true);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    TH2D* h2 = GetH2(f, "h2d_evsel_ZDC_E_tot_vs_FCal_Et_AC");
    const int nx = h2->GetNbinsX();
    const int ny = h2->GetNbinsY();
    const int n_slices = nx / N_BINS_PER_SLICE;

    // cut_per_xbin[ix] = cut value for x-bin ix (1-based); -1 = no valid fit
    std::vector<double> cut_per_xbin(nx + 2, -1.0);

    // for Canvas 1 overlay: (x_center, cut)
    std::vector<double> gr_x, gr_cut;
    // for debug canvas: mu, sigma, and slice half-width
    std::vector<double> gr_mu, gr_sigma, gr_ex;

    // first-converged-slice record for 1D debug canvas
    struct FirstSlice {
        bool   found = false;
        int    ix_lo, ix_hi;
        double xlo, xhi;
        double amp2, mu2, sig2, lo2, hi2;   // pass-2 fit result + range
    } first_sl;

    TF1 gfunc("gfit", "gaus");

    // ---- two-pass per-slice fitting -----------------------------------------
    for (int isl = 0; isl < n_slices; ++isl) {
        const int ix_lo = isl * N_BINS_PER_SLICE + 1;
        const int ix_hi = std::min((isl + 1) * N_BINS_PER_SLICE, nx);

        TH1D* hpy = (TH1D*)h2->ProjectionY(Form("__hpy_%d", isl), ix_lo, ix_hi);
        hpy->SetDirectory(nullptr);

        if (hpy->GetEntries() < MIN_ENTRIES) { delete hpy; continue; }

        // Initial guess: bin with max content
        const double mu0  = hpy->GetBinCenter(hpy->GetMaximumBin());
        const double rms  = (hpy->GetRMS() > 0.) ? hpy->GetRMS() : 1.0;
        const double ymin = hpy->GetXaxis()->GetXmin();
        const double ymax = hpy->GetXaxis()->GetXmax();

        // Pass 1
        const double lo1 = std::max(ymin, mu0 - 3.0 * rms);
        const double hi1 = std::min(ymax, mu0 + 3.0 * rms);
        gfunc.SetRange(lo1, hi1);
        gfunc.SetParameters(hpy->GetMaximum(), mu0, rms);
        if (hpy->Fit(&gfunc, "RNQ") != 0 || gfunc.GetParameter(2) <= 0.) {
            delete hpy; continue;
        }
        const double mu1 = gfunc.GetParameter(1), sig1 = gfunc.GetParameter(2);

        // Pass 2
        const double lo2 = std::max(ymin, mu1 - 3.0 * sig1);
        const double hi2 = std::min(ymax, mu1 + 3.0 * sig1);
        gfunc.SetRange(lo2, hi2);
        gfunc.SetParameters(gfunc.GetParameter(0), mu1, sig1);
        if (hpy->Fit(&gfunc, "RNQ") != 0 || gfunc.GetParameter(2) <= 0.) {
            delete hpy; continue;
        }

        const double mu2  = gfunc.GetParameter(1);
        const double sig2 = gfunc.GetParameter(2);
        const double cut  = mu2 + N_SIGMA_CUT * sig2;

        const double xlo = h2->GetXaxis()->GetBinLowEdge(ix_lo);
        const double xhi = h2->GetXaxis()->GetBinUpEdge(ix_hi);

        if (!first_sl.found) {
            first_sl = {true, ix_lo, ix_hi, xlo, xhi,
                        gfunc.GetParameter(0), mu2, sig2, lo2, hi2};
        }
        delete hpy;
        gr_x.push_back(0.5 * (xlo + xhi));
        gr_cut.push_back(cut);
        gr_mu.push_back(mu2);
        gr_sigma.push_back(sig2);
        gr_ex.push_back(0.5 * (xhi - xlo));

        for (int ix = ix_lo; ix <= ix_hi; ++ix)
            cut_per_xbin[ix] = cut;

        std::cout << Form("  slice %2d [%.2f-%.2f TeV]  mu=%.4f  sig=%.4f  cut=%.4f TeV\n",
                          isl, xlo, xhi, mu2, sig2, cut);
    }
    std::cout << "Fitting done: " << gr_x.size() << " / " << n_slices
              << " slices converged." << std::endl;

    // ---- Canvas 1: 2D + cut markers -----------------------------------------
    {
        TH2D* hd = (TH2D*)h2->Clone("__h2_c1");
        hd->SetDirectory(nullptr);

        TCanvas* c1 = new TCanvas("c_ZDC_FCal_cut", "", 800, 650);
        c1->SetLogz();
        c1->SetRightMargin(0.13);
        hd->SetContour(99);
        hd->Draw("COLZ");

        if (!gr_x.empty()) {
            TGraph* g = new TGraph((int)gr_x.size(), gr_x.data(), gr_cut.data());
            g->SetMarkerStyle(20);
            g->SetMarkerSize(0.9);
            g->SetMarkerColor(kRed);
            g->SetLineColor(kRed);
            g->Draw("P SAME");
        }

        TLatex tl; tl.SetNDC(); tl.SetTextSize(0.036);
        tl.DrawLatex(0.15, 0.87, Form("Pb+Pb 20%s data", yr.c_str()));
        tl.DrawLatex(0.15, 0.81, Form("red: #mu_{2} + %.0f#sigma_{2} per slice", N_SIGMA_CUT));

        c1->SaveAs((out_dir + "/ZDC_E_tot_vs_FCal_Et_AC_with_cut_pbpb_20" + yr + ".png").c_str());
        delete c1;
        delete hd;
    }

    // ---- Canvas 2: passing | failing ----------------------------------------
    {
        TH2D* h_pass = (TH2D*)h2->Clone("__h2_pass");
        TH2D* h_fail = (TH2D*)h2->Clone("__h2_fail");
        h_pass->Reset(); h_pass->SetDirectory(nullptr);
        h_fail->Reset(); h_fail->SetDirectory(nullptr);
        h_pass->SetTitle(Form(";FCal E_{T}^{A+C} [TeV];ZDC E_{total} [TeV]"));
        h_fail->SetTitle(Form(";FCal E_{T}^{A+C} [TeV];ZDC E_{total} [TeV]"));

        for (int ix = 1; ix <= nx; ++ix) {
            const double cut = cut_per_xbin[ix];
            for (int iy = 1; iy <= ny; ++iy) {
                const double w = h2->GetBinContent(ix, iy);
                if (w == 0.) continue;
                const double xc = h2->GetXaxis()->GetBinCenter(ix);
                const double yc = h2->GetYaxis()->GetBinCenter(iy);
                // Bins with no valid fit go to pass by default
                if (cut >= 0. && yc > cut)
                    h_fail->Fill(xc, yc, w);
                else
                    h_pass->Fill(xc, yc, w);
            }
        }

        TCanvas* c2 = new TCanvas("c_ZDC_FCal_pass_fail", "", 1400, 620);
        c2->Divide(2, 1);

        c2->cd(1);
        gPad->SetLogz(); gPad->SetRightMargin(0.15);
        h_pass->SetContour(99);
        h_pass->Draw("COLZ");
        TLatex tl1; tl1.SetNDC(); tl1.SetTextSize(0.036);
        tl1.DrawLatex(0.15, 0.87, Form("Pb+Pb 20%s data", yr.c_str()));
        tl1.DrawLatex(0.15, 0.81, Form("Passing: ZDC E_{total} #leq #mu_{2} + %.0f#sigma_{2}", N_SIGMA_CUT));

        c2->cd(2);
        gPad->SetLogz(); gPad->SetRightMargin(0.15);
        h_fail->SetContour(99);
        h_fail->Draw("COLZ");
        TLatex tl2; tl2.SetNDC(); tl2.SetTextSize(0.036);
        tl2.DrawLatex(0.15, 0.87, Form("Pb+Pb 20%s data", yr.c_str()));
        tl2.DrawLatex(0.15, 0.81, Form("Failing: ZDC E_{total} > #mu_{2} + %.0f#sigma_{2}", N_SIGMA_CUT));

        c2->SaveAs((out_dir + "/ZDC_E_tot_vs_FCal_Et_AC_pass_fail_pbpb_20" + yr + ".png").c_str());
        delete c2;
        delete h_pass;
        delete h_fail;
    }

    // ---- Canvas 3 (debug): 2D + mu+5sigma dots + mu with x/y=sigma error bars --
    {
        TH2D* hd = (TH2D*)h2->Clone("__h2_c3");
        hd->SetDirectory(nullptr);

        TCanvas* c3 = new TCanvas("c_ZDC_FCal_debug", "", 800, 650);
        c3->SetLogz();
        c3->SetRightMargin(0.13);
        hd->SetContour(99);
        hd->Draw("COLZ");

        const int n = (int)gr_x.size();
        if (n > 0) {
            // mu+5sigma dots (filled circle, red)
            TGraph* g_cut = new TGraph(n, gr_x.data(), gr_cut.data());
            g_cut->SetMarkerStyle(20);
            g_cut->SetMarkerSize(0.9);
            g_cut->SetMarkerColor(kRed);
            g_cut->SetLineColor(kRed);
            g_cut->Draw("P SAME");

            // mu with x=slice half-width, y=sigma error bars (open circle, red)
            TGraphErrors* g_mu = new TGraphErrors(n,
                gr_x.data(), gr_mu.data(),
                gr_ex.data(), gr_sigma.data());
            g_mu->SetMarkerStyle(24);   // open circle
            g_mu->SetMarkerSize(0.9);
            g_mu->SetMarkerColor(kRed);
            g_mu->SetLineColor(kRed);
            g_mu->Draw("P SAME");
        }

        TLatex tl; tl.SetNDC(); tl.SetTextSize(0.036);
        tl.DrawLatex(0.15, 0.87, Form("Pb+Pb 20%s data", yr.c_str()));
        tl.DrawLatex(0.15, 0.81, Form("#bullet  #mu_{2} + %.0f#sigma_{2}  (filled)", N_SIGMA_CUT));
        tl.DrawLatex(0.15, 0.75, "#circ  #mu_{2} #pm #sigma_{2}  (open, x = slice width)");

        c3->SaveAs((out_dir + "/ZDC_E_tot_vs_FCal_Et_AC_with_cut_pbpb_20" + yr + "_debug.png").c_str());
        delete c3;
        delete hd;
    }

    // ---- Canvas 4: 1D ZDC projection of first converged slice + pass-2 fit ----
    if (first_sl.found) {
        TH1D* h1d = (TH1D*)h2->ProjectionY("__h1d_first_slice",
                                            first_sl.ix_lo, first_sl.ix_hi);
        h1d->SetDirectory(nullptr);
        h1d->SetTitle(";ZDC E_{total} [TeV];Counts");

        // Re-run pass-2 fit so TF1 carries the correct range for drawing
        TF1 gfit2("gfit2_1d", "gaus");
        gfit2.SetRange(first_sl.lo2, first_sl.hi2);
        gfit2.SetParameters(first_sl.amp2, first_sl.mu2, first_sl.sig2);
        h1d->Fit(&gfit2, "RNQ");
        gfit2.SetLineColor(kRed);
        gfit2.SetLineWidth(2);

        TCanvas* c4 = new TCanvas("c_1d_first_slice", "", 750, 600);
        c4->SetLeftMargin(0.13);
        c4->SetBottomMargin(0.13);

        h1d->SetMarkerStyle(20);
        h1d->SetMarkerSize(0.7);
        h1d->Draw("E");                   // markers + Poisson error bars
        gfit2.Draw("SAME");

        // Extend drawn range of fit to full ±5 sigma for visual context
        TF1 gfit2_ext("gfit2_ext", "gaus",
                      first_sl.mu2 - 5.*first_sl.sig2,
                      first_sl.mu2 + 5.*first_sl.sig2);
        gfit2_ext.SetParameters(gfit2.GetParameter(0),
                                gfit2.GetParameter(1),
                                gfit2.GetParameter(2));
        gfit2_ext.SetLineColor(kRed);
        gfit2_ext.SetLineWidth(2);
        gfit2_ext.SetLineStyle(2);        // dashed outside fit range
        gfit2_ext.Draw("SAME");

        // Info box: FCal slice, mu, sigma, mu+5sigma
        TLatex tl; tl.SetNDC(); tl.SetTextSize(0.038);
        tl.DrawLatex(0.55, 0.82,
            Form("FCal E_{T}^{A+C} #in [%.2f, %.2f] TeV",
                 first_sl.xlo, first_sl.xhi));
        tl.DrawLatex(0.55, 0.75,
            Form("#mu_{2} = %.1f TeV", first_sl.mu2));
        tl.DrawLatex(0.55, 0.68,
            Form("#sigma_{2} = %.1f TeV", first_sl.sig2));
        tl.DrawLatex(0.55, 0.61,
            Form("#mu_{2} + 5#sigma_{2} = %.1f TeV",
                 first_sl.mu2 + N_SIGMA_CUT * first_sl.sig2));
        // mark fit range
        tl.SetTextSize(0.032);
        tl.DrawLatex(0.55, 0.53,
            Form("Fit range: [%.1f, %.1f] TeV", first_sl.lo2, first_sl.hi2));

        c4->SaveAs((out_dir + "/ZDC_1D_first_slice_fit_pbpb_20" + yr + "_debug.png").c_str());
        delete c4;
        delete h1d;
    }

    f->Close();
    std::cout << "Saved to " << out_dir << std::endl;
}
