// plot_pbpb_event_selection.cxx
// Draws event-selection diagnostic plots for PbPb data:
//   - FCal Et A vs C [TeV]
//   - ZDC E_total (Y) vs FCal Et A+C (X) [TeV]
//   - nTrk HIloose/HItight (Y) vs FCal Et A+C (X) [TeV] (2 subplots)
//   - ZDC time A vs C [ns]
//
// Usage in ROOT:
//   .L plot_pbpb_event_selection.cxx+
//   plot_pbpb_event_selection(24)   // run_year = 23/24/25

#include <iostream>
#include <string>
#include <stdexcept>
#include "TFile.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLatex.h"

// ---- helpers ----------------------------------------------------------------

static TH2D* GetH2(TFile* f, const std::string& name) {
    TH2D* h = dynamic_cast<TH2D*>(f->Get(name.c_str()));
    if (!h) throw std::runtime_error("Histogram not found: " + name);
    h->SetDirectory(nullptr);
    return h;
}

static void DrawH2(TH2D* h, const std::string& opt = "COLZ") {
    h->SetContour(99);
    h->Draw(opt.c_str());
}

static void SaveCanvas(TCanvas* c, const std::string& out_dir, const std::string& name) {
    c->SaveAs((out_dir + "/" + name + ".png").c_str());
}

// ---- main -------------------------------------------------------------------

void plot_pbpb_event_selection(int run_year = 24) {
    run_year %= 2000;
    const std::string yr = std::to_string(run_year);

    // Input file: RDFBasedHistFillingPbPb output for single_mu4 trigger mode
    const std::string infile =
        "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_20" + yr +
        "/histograms_real_pairs_pbpb_20" + yr + "_single_mu4_no_trg_plots_coarse_q_eta_bin.root";

    const std::string out_dir =
        "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/single_b_analysis/event_selection";

    // Open file
    if (gSystem->AccessPathName(infile.c_str()))
        throw std::runtime_error("Input file not found: " + infile);
    TFile* f = TFile::Open(infile.c_str(), "READ");
    if (!f || f->IsZombie())
        throw std::runtime_error("Cannot open: " + infile);

    gSystem->mkdir(out_dir.c_str(), true);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    const std::string label = "Pb+Pb 20" + yr + " data";

    auto AddLabel = [&](double x, double y) {
        TLatex* t = new TLatex(x, y, label.c_str());
        t->SetNDC();
        t->SetTextSize(0.038);
        t->DrawLatex(x, y, label.c_str());
    };

    // ---- 1. FCal Et A vs C ----
    {
        TH2D* h = GetH2(f, "h2d_evsel_FCal_Et_C_vs_A");
        TCanvas* c = new TCanvas("c_FCal_AC", "", 700, 600);
        c->SetLogz();
        c->SetRightMargin(0.13);
        DrawH2(h);
        AddLabel(0.15, 0.87);
        SaveCanvas(c, out_dir, "FCal_Et_A_vs_C_pbpb_20" + yr);
        delete c;
    }

    // ---- 2. ZDC total energy (Y) vs FCal Et A+C (X) [TeV] ----
    {
        TH2D* h = GetH2(f, "h2d_evsel_ZDC_E_tot_vs_FCal_Et_AC");
        TCanvas* c = new TCanvas("c_ZDC_FCal", "", 700, 600);
        c->SetLogz();
        c->SetRightMargin(0.13);
        DrawH2(h);
        AddLabel(0.15, 0.87);
        SaveCanvas(c, out_dir, "ZDC_E_tot_vs_FCal_Et_AC_pbpb_20" + yr);
        delete c;
    }

    // ---- 3. FCal Et (A+C) vs nTrk HIloose + HItight (2 subplots) ----
    {
        TH2D* h_loose = GetH2(f, "h2d_evsel_ntrk_HIloose_vs_FCal_Et_AC");
        TH2D* h_tight = GetH2(f, "h2d_evsel_ntrk_HItight_vs_FCal_Et_AC");

        TCanvas* c = new TCanvas("c_FCal_ntrk", "", 1300, 580);
        c->Divide(2, 1);

        c->cd(1);
        gPad->SetLogz();
        gPad->SetRightMargin(0.14);
        DrawH2(h_loose);
        TLatex tl;
        tl.SetNDC(); tl.SetTextSize(0.038);
        tl.DrawLatex(0.15, 0.87, label.c_str());
        tl.DrawLatex(0.15, 0.81, "N_{trk}^{HIloose}");

        c->cd(2);
        gPad->SetLogz();
        gPad->SetRightMargin(0.14);
        DrawH2(h_tight);
        TLatex tl2;
        tl2.SetNDC(); tl2.SetTextSize(0.038);
        tl2.DrawLatex(0.15, 0.87, label.c_str());
        tl2.DrawLatex(0.15, 0.81, "N_{trk}^{HItight}");

        SaveCanvas(c, out_dir, "FCal_Et_AC_vs_ntrk_pbpb_20" + yr);
        delete c;
    }

    // ---- 4. ZDC time A vs C ----
    {
        TH2D* h = GetH2(f, "h2d_evsel_ZDC_t_C_vs_A");
        TCanvas* c = new TCanvas("c_ZDC_t", "", 700, 600);
        c->SetLogz();
        c->SetRightMargin(0.13);
        DrawH2(h);
        AddLabel(0.15, 0.87);
        SaveCanvas(c, out_dir, "ZDC_t_A_vs_C_pbpb_20" + yr);
        delete c;
    }

    f->Close();
    std::cout << "Saved plots to " << out_dir << std::endl;
}
