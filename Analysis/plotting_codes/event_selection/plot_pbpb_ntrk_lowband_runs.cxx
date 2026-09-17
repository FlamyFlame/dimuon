// Evidence record for a run block whose N_trk^HItight-vs-FCal correlation sits on a SECOND band
// (first seen in Pb+Pb 2026: runs 522336 / 522355 / 522384 / 522408, ~30 % below the main band).
//
// Whether such runs are excluded from data AND luminosity is a decision for the human experts
// (DQ / detector groups) -- this macro only PRODUCES THE RECORD they need:
//   (1) <out>/ntrk_lowband_ntrk_1d_pbpb_20YY.png  -- N_trk^HItight distribution, one panel per
//       problematic run + one panel for all other runs, identical axes;
//   (2) <out>/ntrk_lowband_ntrk_vs_fcal_pbpb_20YY.png -- N_trk^HItight vs FCal E_T^{A+C}, the SAME
//       binning and range as the cut-5 figures of plot_pbpb_event_sel_cuts.cxx
//       ({120,-0.5,5.5} x {100,0,3500}), per problematic run + all other runs, with the cut-5
//       band (g_ntrk_fcal_cut_lo/hi from the year's cuts file) overlaid so one can read off
//       whether cut 5 already removes the population;
//   (3) <out>/ntrk_lowband_ratio_1d_pbpb_20YY.png -- N_trk^HItight / FCal E_T^{A+C} [1/TeV] for
//       FCal > 1 TeV, the variable that actually separates the two bands, same panel layout;
//   (4) <out>/ntrk_lowband_runs_pbpb_20YY.txt -- per-run table over ALL runs (events, fraction in
//       the low band, LB range of the low-band events, cut-5 rejection fraction), plus a per-LB
//       listing for the problematic runs, so the identification is self-contained.
//
// Event reading mirrors plot_pbpb_event_sel_cuts.cxx exactly: HLT_mu4_L1MU3V required,
// N_trk^HItight = trk_numqual[3], FCal E_T^{A+C} = (FCal_Et_P + FCal_Et_N) x 1e-6 TeV, all merged
// part files of the year (part list = file_batch_max, see pbpb2026_analysis_support.md).  No
// event-selection cut is applied to the distributions (the point is to see the population
// before cut 5); the cut-5 band is drawn, not applied.
//
// "Low band" for the table: N_trk^HItight below the cut-5 LOWER edge at FCal E_T^{A+C} > 1.5 TeV
// (below that the two bands merge and the classification is meaningless).
//
// Usage:  root -l -b -q 'plot_pbpb_ntrk_lowband_runs.cxx+(26, {522336, 522355, 522384, 522408})'
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"

#include "../../NTupleProcessingCode/PbPbEventSelConfig.h"

void plot_pbpb_ntrk_lowband_runs(int run_year, std::vector<int> problematic_runs, int nparts = 7)
{
    const int yr = run_year % 2000;
    const std::string yrs  = std::to_string(yr);
    const std::string base = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/";
    const std::string out  = base + "plots/single_b_analysis/event_selection/pbpb_20" + yrs + "/ntrk_lowband_runs";
    gSystem->mkdir(out.c_str(), true);

    // ---- cut-5 band from the year's cuts file (drawn, not applied) ----
    const std::string cuts_path = base + "pbpb_20" + yrs + "/event_sel_cuts_pbpb_20" + yrs + ".root";
    TFile* fc = TFile::Open(cuts_path.c_str(), "READ");
    if (!fc || fc->IsZombie()) throw std::runtime_error("cannot open " + cuts_path);
    TGraph* g_lo = dynamic_cast<TGraph*>(fc->Get(PbPbEvSelKey::kNTrkFCalCutLo));
    TGraph* g_hi = dynamic_cast<TGraph*>(fc->Get(PbPbEvSelKey::kNTrkFCalCutHi));
    if (!g_lo || !g_hi) throw std::runtime_error("cuts file lacks the cut-5 band graphs: " + cuts_path);

    // ---- chain: every merged part of the year ----
    TChain chain("HeavyIonD3PD");
    for (int p = 1; p <= nparts; ++p) {
        const std::string f = base + "pbpb_20" + yrs + "/data_pbpb" + yrs + "_part" + std::to_string(p) + ".root";
        if (gSystem->AccessPathName(f.c_str())) throw std::runtime_error("missing part file " + f);
        chain.Add(f.c_str());
    }
    chain.SetMakeClass(1);
    chain.SetBranchStatus("*", 0);
    Int_t b_HLT = 0, run_num = 0, lbn = 0;
    Float_t FCal_Et_P = 0.f, FCal_Et_N = 0.f;
    std::vector<int>* trk_numqual = nullptr;
    for (const char* b : {"b_HLT_mu4_L1MU3V", "RunNumber", "lbn", "FCal_Et_P", "FCal_Et_N", "trk_numqual"})
        chain.SetBranchStatus(b, 1);
    chain.SetBranchAddress("b_HLT_mu4_L1MU3V", &b_HLT);
    chain.SetBranchAddress("RunNumber",         &run_num);
    chain.SetBranchAddress("lbn",               &lbn);
    chain.SetBranchAddress("FCal_Et_P",         &FCal_Et_P);
    chain.SetBranchAddress("FCal_Et_N",         &FCal_Et_N);
    chain.SetBranchAddress("trk_numqual",       &trk_numqual);

    // ---- histograms: one set per problematic run + "all other runs" ----
    const std::set<int> prob(problematic_runs.begin(), problematic_runs.end());
    const int npan = static_cast<int>(problematic_runs.size()) + 1;
    auto panel_of = [&](int run) -> int {
        auto it = std::find(problematic_runs.begin(), problematic_runs.end(), run);
        return it == problematic_runs.end() ? npan - 1 : static_cast<int>(it - problematic_runs.begin());
    };
    auto panel_label = [&](int ip) -> std::string {
        return ip == npan - 1 ? "all other runs" : "run " + std::to_string(problematic_runs[ip]);
    };
    std::vector<TH1D*> h1(npan), hr(npan);
    std::vector<TH2D*> h2(npan);
    for (int ip = 0; ip < npan; ++ip) {
        h1[ip] = new TH1D(Form("h_ntrk_%d", ip), "", 100, 0., 3500.);
        h2[ip] = new TH2D(Form("h_ntrk_fcal_%d", ip), "", 120, -0.5, 5.5, 100, 0., 3500.);  // = cut-5 figure axes
        hr[ip] = new TH1D(Form("h_ratio_%d", ip), "", 100, 0., 1000.);
        for (TH1* h : {static_cast<TH1*>(h1[ip]), static_cast<TH1*>(h2[ip]), static_cast<TH1*>(hr[ip])})
            h->SetDirectory(nullptr);
    }

    // ---- per-run / per-LB bookkeeping ----
    struct RunRec { long n = 0, n_hi = 0, n_low = 0, n_c5fail = 0; int lb_min = 1 << 30, lb_max = -1, lb_low_min = 1 << 30, lb_low_max = -1;
                    std::map<int, std::pair<long,long>> lb; };  // lb -> (n at FCal>1.5, n low)
    std::map<int, RunRec> runs;

    const double FCAL_CLASSIFY_MIN = 1.5;  // TeV: below this the two bands merge
    const Long64_t n_tot = chain.GetEntries();
    std::cout << "PbPb 20" << yr << ": " << n_tot << " events in " << nparts << " parts" << std::endl;
    for (Long64_t i = 0; i < n_tot; ++i) {
        chain.GetEntry(i);
        if (!b_HLT) continue;
        if (!trk_numqual || trk_numqual->size() < 4) continue;
        const int    ntrk = (*trk_numqual)[3];
        const double fcal = (FCal_Et_P + FCal_Et_N) * 1e-6;
        RunRec& r = runs[run_num];
        ++r.n; r.lb_min = std::min(r.lb_min, lbn); r.lb_max = std::max(r.lb_max, lbn);
        const int ip = panel_of(run_num);
        h1[ip]->Fill(ntrk);
        h2[ip]->Fill(fcal, ntrk);
        if (fcal > 1.0) hr[ip]->Fill(ntrk / fcal);
        const double lo = PbPbEvSelEvalCut(g_lo, fcal), hi = PbPbEvSelEvalCut(g_hi, fcal);
        if (ntrk < lo || ntrk > hi) ++r.n_c5fail;
        if (fcal > FCAL_CLASSIFY_MIN) {
            ++r.n_hi; ++r.lb[lbn].first;
            if (ntrk < lo) { ++r.n_low; ++r.lb[lbn].second; r.lb_low_min = std::min(r.lb_low_min, lbn); r.lb_low_max = std::max(r.lb_low_max, lbn); }
        }
        if (i % 20000000 == 0 && i) std::cout << "  " << i << " / " << n_tot << std::endl;
    }
    fc->Close();

    // ---- layout: nrows >= ncols, nrows ~ sqrt(N) ----
    const int ncols = std::max(1, static_cast<int>(std::floor(std::sqrt(npan))));
    const int nrows = (npan + ncols - 1) / ncols;
    gStyle->SetOptStat(0);
    auto draw_set = [&](const char* stem, auto& hs, bool is2d, const char* xt, const char* yt, bool logy) {
        TCanvas c(Form("c_%s", stem), "", 600 * ncols, 500 * nrows);
        c.Divide(ncols, nrows, 0.003, 0.003);
        for (int ip = 0; ip < npan; ++ip) {
            c.cd(ip + 1);
            gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13); gPad->SetRightMargin(is2d ? 0.14 : 0.05);
            if (is2d) gPad->SetLogz(); else if (logy) gPad->SetLogy();
            hs[ip]->SetTitle(Form("Pb+Pb 20%d, %s;%s;%s", yr, panel_label(ip).c_str(), xt, yt));
            hs[ip]->Draw(is2d ? "COLZ" : "HIST");
            if (is2d) { g_lo->SetLineColor(kRed); g_hi->SetLineColor(kRed); g_lo->Draw("L same"); g_hi->Draw("L same"); }
        }
        const std::string png = out + "/ntrk_lowband_" + stem + "_pbpb_20" + yrs + ".png";
        c.SaveAs(png.c_str());
        std::cout << "[INFO] Saved: " << png << std::endl;
    };
    draw_set("ntrk_1d",       h1, false, "N_{trk}^{HItight}", "events", true);
    draw_set("ntrk_vs_fcal",  h2, true,  "FCal E_{T}^{A+C} [TeV]", "N_{trk}^{HItight}", false);
    draw_set("ratio_1d",      hr, false, "N_{trk}^{HItight} / FCal E_{T}^{A+C}  [TeV^{-1}]  (FCal > 1 TeV)", "events", false);

    // ---- table ----
    const std::string txt = out + "/ntrk_lowband_runs_pbpb_20" + yrs + ".txt";
    std::ofstream o(txt);
    o << "Pb+Pb 20" << yr << " -- N_trk^HItight low-band record (HLT_mu4_L1MU3V events, no event-selection cut applied)\n"
      << "low band := N_trk^HItight below the cut-5 LOWER band edge (" << cuts_path << ") at FCal E_T^{A+C} > "
      << FCAL_CLASSIFY_MIN << " TeV\n"
      << "cut5_fail := outside [lo, hi] of the cut-5 band at any FCal (what cut 5 would reject)\n\n";
    o << Form("%-8s %10s %8s %8s %10s %8s %9s %14s %14s\n", "run", "events", "LB_min", "LB_max", "N(FCal>1.5)", "N_low", "low_frac", "LB_low[min,max]", "cut5_fail_frac");
    for (auto& kv : runs) {
        const RunRec& r = kv.second;
        o << Form("%-8d %10ld %8d %8d %10ld %8ld %8.3f%% %6d,%-7d %13.3f%%\n", kv.first, r.n, r.lb_min, r.lb_max, r.n_hi, r.n_low,
                  r.n_hi ? 100. * r.n_low / r.n_hi : 0., r.n_low ? r.lb_low_min : -1, r.n_low ? r.lb_low_max : -1,
                  r.n ? 100. * r.n_c5fail / r.n : 0.);
    }
    o << "\nPer-LB low-band fraction (FCal > " << FCAL_CLASSIFY_MIN << " TeV) for the problematic runs:\n";
    for (int run : problematic_runs) {
        auto it = runs.find(run);
        if (it == runs.end()) { o << "run " << run << ": NOT IN THE SKIM\n"; continue; }
        o << "run " << run << ":\n";
        for (auto& lbkv : it->second.lb)
            o << Form("  LB %5d  N=%7ld  low=%7ld  (%6.2f%%)\n", lbkv.first, lbkv.second.first, lbkv.second.second,
                      lbkv.second.first ? 100. * lbkv.second.second / lbkv.second.first : 0.);
    }
    o.close();
    std::cout << "[INFO] Saved: " << txt << std::endl;
    // echo the per-run summary for the log
    std::ifstream in(txt); std::string line; int n = 0;
    while (std::getline(in, line) && n++ < 60) std::cout << line << std::endl;
}
