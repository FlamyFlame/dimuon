// Abnormal-run investigation (docs/tracking/data_per_run_signal_yield_qa.md R2d): where does a per-run N_signal/L
// deviation come from?  N_sig/L = (N_HLT/L) x (N_evsel/N_HLT) x (N_sig/N_evsel) -- one panel per factor, all runs,
// target runs highlighted:
//   N_HLT   = HLT_mu4 events (raw NTUP, good runs/LBs)                     [evsel_single_run_fill output]
//   N_evsel = of those, passing all 5 event-selection cuts, centrality 0-80 % [same]
//   N_sig   = signal-region OS Tight pairs (nominal muon-pair NTuple, centrality 0-80 %), counted as distinct events
//   L       = Prescale-Corrected HLT_mu4 lumi (ub^-1 -> nb^-1)
// A lumi/DAQ-type problem moves N_HLT/L; a physics/reco problem moves N_sig/N_evsel.
// Per target run: N_sig/N_evsel and N_evsel per LB (lumi-free, within-run time dependence).
// Also the two factors N_sig/N_evsel and N_evsel/L as stand-alone per-run PNGs (one per year).
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include "ROOT/RDataFrame.hxx"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include "MuonObjectsParamsAndHelpers/ParamsSet.h"
#include "NTupleProcessingCode/PbPbEventSelConfig.h"

namespace rnd {
inline std::map<int,double> ReadLumi(const std::string& path) {   // Prescale Corrected = column 7
    std::map<int,double> m; std::ifstream f(path); std::string line; std::getline(f, line);
    while (std::getline(f, line)) {
        std::stringstream ss(line); std::string c; std::vector<std::string> v;
        while (std::getline(ss, c, ',')) v.push_back(c);
        if (v.size() < 7 || v[0].find_first_not_of("0123456789 ") != std::string::npos) continue;
        m[std::stoi(v[0])] = std::stod(v[6]);
    }
    return m;
}
}

void run_norm_diagnostics(int yr, const std::string& evsel_root, const std::string& pair_file, const std::string& lumi_csv,
                          std::vector<int> targets, const std::string& plotbase, int lb_group = 10) {
    using namespace rnd;
    yr %= 2000;
    ROOT::EnableImplicitMT(); gStyle->SetOptStat(0);
    const auto lumi = ReadLumi(lumi_csv);
    const auto bad_runs = PbPbBadRuns(yr); const auto bad_lbs = PbPbBadLBRanges(yr);

    // signal events per (run, LB), same selection as plot_signal_yield_per_lumi_per_run
    const std::string sig =
        ParamsSet::SignalMinvCutExpr("minv") + " && " + ParamsSet::SignalPairPtCutExpr("pair_pt") + " && "
        + ParamsSet::FiducialGapCutExpr("m1.charge * m1.eta") + " && " + ParamsSet::FiducialGapCutExpr("m2.charge * m2.eta")
        + " && " + ParamsSet::PairFiducialEtaCutExpr("pair_eta");
    const std::string ctr = "avg_centrality >= " + std::to_string(ParamsSet::ctrbins.front())
                          + " && avg_centrality < " + std::to_string(ParamsSet::ctrbins.back());
    auto good = [bad_runs, bad_lbs](UInt_t r, UInt_t lb) { return !PbPbIsBadRunLB(bad_runs, bad_lbs, (int)r, lb); };
    auto n = ROOT::RDataFrame("muon_pair_tree_sign2", pair_file).Filter(good, {"run_number", "lb"})
                 .Filter("pair_pass_tight").Filter(ctr).Filter(sig);
    auto vr = n.Take<UInt_t>("run_number"); auto vl = n.Take<UInt_t>("lb");
    auto vb = n.Take<UInt_t>("bcid");       auto ve = n.Take<Int_t>("m1.ev_num");
    std::map<int, std::set<std::tuple<unsigned,unsigned,int>>> evs;       // distinct events per run
    std::map<int, std::map<unsigned, std::set<std::pair<unsigned,int>>>> evs_lb;
    for (size_t i = 0; i < vr->size(); ++i) {
        evs[(*vr)[i]].insert({(*vl)[i], (*vb)[i], (*ve)[i]});
        evs_lb[(*vr)[i]][(*vl)[i]].insert({(*vb)[i], (*ve)[i]});
    }

    TFile* f = TFile::Open(evsel_root.c_str());
    auto* h_all = (TH2D*)f->Get("h_stage_all"); auto* h_c80 = (TH2D*)f->Get("h_stage_c80");
    struct R { int run; double L, nhlt, nev, nsig; };
    std::vector<R> rows;
    for (int b = 1; b <= h_all->GetNbinsX(); ++b) {
        const double nhlt = h_all->GetBinContent(b, 1);
        if (nhlt <= 0) continue;
        const int r = (int)h_all->GetXaxis()->GetBinLowEdge(b);
        auto it = lumi.find(r);
        if (it == lumi.end() || it->second <= 0) { std::cout << "run " << r << ": no lumi, skipped\n"; continue; }
        // Partially-bad runs (PbPbBadLBRanges): the per-run lumi row covers the bad LBs too -> L unknown (NaN);
        // the run still enters the lumi-free factors.
        double L = it->second / 1000.;
        for (auto& br : bad_lbs) if (br.run == r) { L = std::nan(""); std::cout << "run " << r << ": good-LB lumi pending (L-dependent factors not drawn)\n"; }
        rows.push_back({r, L, nhlt, h_c80->GetBinContent(b, 6), (double)evs[r].size()});
    }

    // three factors vs run index
    struct Fac { const char* name; const char* ytitle; };
    const Fac fac[3] = {{"hlt_per_lumi", "N_{HLT_mu4} / L [nb]"},
                        {"evsel_frac", "N_{evsel, 0-80%} / N_{HLT_mu4}"},
                        {"sig_per_event", "N_{signal} / N_{evsel, 0-80%}"}};
    const std::string dir = plotbase + "/signal_yield_factored_per_run";   // the original N/L plots stay in signal_yield_per_lumi_per_run/
    gSystem->mkdir(dir.c_str(), true);
    std::ofstream o(dir + Form("/pbpb%d_run_norm_factors.csv", yr));
    o << "run,L_nb,n_hlt,n_evsel_c80,n_sig,hlt_per_L,evsel_frac,sig_per_evsel,pull_hlt_per_L,pull_sig_per_evsel\n";
    std::vector<double> val[3], err[3];
    for (auto& r : rows) {
        val[0].push_back(r.nhlt / r.L);  err[0].push_back(std::sqrt(r.nhlt) / r.L);
        val[1].push_back(r.nev / r.nhlt); err[1].push_back(std::sqrt(r.nev / r.nhlt * (1 - r.nev / r.nhlt) / r.nhlt));
        val[2].push_back(r.nsig / r.nev); err[2].push_back(std::sqrt(std::max(r.nsig, 1.)) / r.nev);
    }
    double avg[3];
    for (int k = 0; k < 3; ++k) {   // dataset average = ratio of sums (runs with unknown L excluded from k = 0)
        double num = 0, den = 0;
        for (auto& r : rows) { if (k == 0 && std::isnan(r.L)) continue;
            num += (k == 0 ? r.nhlt : k == 1 ? r.nev : r.nsig); den += (k == 0 ? r.L : k == 1 ? r.nhlt : r.nev); }
        avg[k] = num / den;
    }
    for (size_t i = 0; i < rows.size(); ++i)
        o << rows[i].run << "," << rows[i].L << "," << rows[i].nhlt << "," << rows[i].nev << "," << rows[i].nsig << ","
          << val[0][i] << "," << val[1][i] << "," << val[2][i] << ","
          << (val[0][i] - avg[0]) / err[0][i] << "," << (val[2][i] - avg[2]) / err[2][i] << "\n";

    TCanvas c("c_fac", "", 1100, 1000); c.Divide(1, 3, 0, 0);
    for (int k = 0; k < 3; ++k) {
        c.cd(k + 1); gPad->SetLeftMargin(0.10); gPad->SetRightMargin(0.03); gPad->SetGridx();
        if (k == 2) gPad->SetBottomMargin(0.22); if (k == 0) gPad->SetTopMargin(0.08);
        const int N = rows.size();
        TH1D* fr = new TH1D(Form("fr%d", k), "", N, 0, N);
        for (int i = 0; i < N; ++i) fr->GetXaxis()->SetBinLabel(i + 1, Form("%d", rows[i].run));
        fr->GetXaxis()->LabelsOption("v"); fr->GetXaxis()->SetLabelSize(k == 2 ? 0.06 : 0);
        double lo = 1e30, hi = -1e30;
        for (int i = 0; i < N; ++i) { if (std::isnan(val[k][i])) continue;
            lo = std::min(lo, val[k][i] - err[k][i]); hi = std::max(hi, val[k][i] + err[k][i]); }
        fr->SetMinimum(lo - 0.1 * (hi - lo)); fr->SetMaximum(hi + 0.25 * (hi - lo));
        fr->GetYaxis()->SetTitle(fac[k].ytitle); fr->GetYaxis()->SetTitleSize(0.07); fr->GetYaxis()->SetTitleOffset(0.65);
        fr->GetYaxis()->SetLabelSize(0.06); fr->Draw("AXIS");
        TGraphErrors* g = new TGraphErrors(); TGraphErrors* gt = new TGraphErrors();
        for (int i = 0; i < N; ++i) {
            if (std::isnan(val[k][i])) continue;
            bool is_t = std::find(targets.begin(), targets.end(), rows[i].run) != targets.end();
            TGraphErrors* gg = is_t ? gt : g; int p = gg->GetN();
            gg->SetPoint(p, i + 0.5, val[k][i]); gg->SetPointError(p, 0, err[k][i]);
        }
        g->SetMarkerStyle(20); g->SetMarkerSize(0.8); g->Draw("P SAME");
        gt->SetMarkerStyle(21); gt->SetMarkerColor(kRed + 1); gt->SetLineColor(kRed + 1); gt->SetMarkerSize(1.1); gt->Draw("P SAME");
        TLine l(0, avg[k], N, avg[k]); l.SetLineStyle(2); l.SetLineColor(kBlue + 1); l.DrawClone();
        if (k == 0) { TLatex tl; tl.SetNDC(); tl.SetTextSize(0.075);
            tl.DrawLatex(0.12, 0.80, Form("Pb+Pb 20%d data, HLT_mu4; red = investigated runs; dashed = dataset average", yr)); }
    }
    c.SaveAs((dir + Form("/pbpb%d_norm_factors_per_run.png", yr)).c_str());

    // The two factors of N_sig/L as stand-alone per-run plots, same style as plot_signal_yield_per_lumi_per_run:
    //   N_sig / N_evsel  (lumi-free; partial runs included)   and   N_evsel / L  (partial runs not drawn)
    auto per_run_plot = [&](const char* stem, const char* ytitle, auto num, auto den, double mean, const char* unit) {
        std::vector<int> idx;
        for (size_t i = 0; i < rows.size(); ++i) if (!std::isnan(den(rows[i]))) idx.push_back(i);
        const int nr = idx.size();
        const int pw = (int)std::floor(std::log10(mean));   // power of ten moved into the y title (no axis exponent label)
        const double sc = std::pow(10., -pw);
        TH1D* fr = new TH1D(Form("frp_%s", stem), "", nr, 0, nr);
        TGraphErrors* g = new TGraphErrors(nr);
        double ymax = 0;
        for (int j = 0; j < nr; ++j) {
            const R& r = rows[idx[j]];
            const double y = sc * num(r) / den(r), ey = sc * std::sqrt(std::max(num(r), 1.)) / den(r);
            g->SetPoint(j, j + 0.5, y); g->SetPointError(j, 0, ey); ymax = std::max(ymax, y + ey);
            fr->GetXaxis()->SetBinLabel(j + 1, std::to_string(r.run).c_str());
        }
        TCanvas* cc = new TCanvas(Form("cp_%s", stem), "", 1400, 750);
        cc->SetBottomMargin(0.16); cc->SetLeftMargin(0.09); cc->SetRightMargin(0.03); cc->SetTopMargin(0.07);
        fr->SetMinimum(0); fr->SetMaximum(1.35 * ymax);
        fr->GetXaxis()->SetLabelSize(nr > 40 ? 0.022 : 0.03); fr->GetXaxis()->LabelsOption("v"); fr->GetXaxis()->SetTickLength(0);
        fr->GetXaxis()->SetTitle("Run number"); fr->GetXaxis()->SetTitleOffset(2.2); fr->GetXaxis()->SetTitleSize(0.04);
        fr->GetYaxis()->SetTitle(Form("%s  [10^{%d}%s]", ytitle, pw, unit)); fr->GetYaxis()->SetTitleSize(0.04); fr->GetYaxis()->SetTitleOffset(1.1);
        fr->Draw("AXIS");
        g->SetMarkerStyle(20); g->SetMarkerSize(1.0); g->SetLineWidth(2);
        TLine* ln = new TLine(0, sc * mean, nr, sc * mean); ln->SetLineStyle(2); ln->SetLineColor(kRed + 1); ln->SetLineWidth(2);
        ln->Draw(); g->Draw("P SAME");
        TLegend* lg = new TLegend(0.60, 0.78, 0.96, 0.90); lg->SetBorderSize(0); lg->SetFillStyle(0);
        lg->AddEntry(g, "per run (stat. error)", "pe");
        lg->AddEntry(ln, Form("dataset average %.4g%s", mean, unit), "l"); lg->Draw();
        TLatex t; t.SetNDC(); t.SetTextSize(0.04); t.SetTextFont(72); t.DrawLatex(0.10, 0.94, "ATLAS");
        t.SetTextFont(42); t.DrawLatex(0.17, 0.94, Form("Internal   Pb+Pb 20%d, HLT_mu4, #sqrt{s_{NN}} = 5.36 TeV", yr));
        cc->SaveAs((dir + Form("/pbpb%d_%s_per_run.png", yr, stem)).c_str());
    };
    per_run_plot("signal_per_evsel_event", "N_{events}^{signal} / N_{events}^{evsel, 0-80%}",
                 [](const R& r) { return r.nsig; }, [](const R& r) { return r.nev; }, avg[2], "");
    double sev = 0, sL = 0;
    for (auto& r : rows) if (!std::isnan(r.L)) { sev += r.nev; sL += r.L; }
    per_run_plot("evsel_events_per_lumi", "N_{events}^{evsel, 0-80%} / L",
                 [](const R& r) { return r.nev; }, [](const R& r) { return r.L; }, sev / sL, " nb");

    // per-LB, target runs
    for (int r : targets) {
        auto* hh = (TH1D*)f->Get(Form("h_lb_hlt_run%d", r)); auto* hp = (TH1D*)f->Get(Form("h_lb_pass_run%d", r));
        if (!hh || !hp) { std::cout << "no LB hists for " << r << "\n"; continue; }
        int lb_lo = hh->FindFirstBinAbove(0) - 1, lb_hi = hh->FindLastBinAbove(0);
        TGraphErrors gs, gn;
        double tot_s = 0, tot_e = 0;
        for (int a = lb_lo; a < lb_hi; a += lb_group) {
            double ne = hp->Integral(a + 1, std::min(a + lb_group, lb_hi)), ns = 0;
            for (int lb = a; lb < a + lb_group; ++lb) { auto it = evs_lb[r].find(lb); if (it != evs_lb[r].end()) ns += it->second.size(); }
            if (ne <= 0) continue;
            const double x = a + 0.5 * lb_group; int p = gs.GetN();
            gs.SetPoint(p, x, ns / ne); gs.SetPointError(p, 0.5 * lb_group, std::sqrt(std::max(ns, 1.)) / ne);
            gn.SetPoint(p, x, ne / lb_group); gn.SetPointError(p, 0.5 * lb_group, 0);
            tot_s += ns; tot_e += ne;
        }
        TCanvas cl(Form("cl%d", r), "", 900, 800); cl.Divide(1, 2, 0, 0);
        cl.cd(1); gPad->SetLeftMargin(0.13); gPad->SetTopMargin(0.08);
        gn.SetMarkerStyle(20); gn.SetTitle(Form(";lumi block;N_{evsel, 0-80%%} / LB"));
        gn.GetYaxis()->SetTitleSize(0.06); gn.GetYaxis()->SetTitleOffset(1.0); gn.Draw("AP");
        TLatex tl; tl.SetNDC(); tl.SetTextSize(0.05);
        tl.DrawLatex(0.16, 0.84, Form("Pb+Pb 20%d data, run %d, groups of %d LB", yr, r, lb_group));
        cl.cd(2); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.14);
        gs.SetMarkerStyle(20); gs.SetTitle(";lumi block;N_{signal} / N_{evsel, 0-80%}");
        gs.GetYaxis()->SetTitleSize(0.06); gs.GetYaxis()->SetTitleOffset(1.0); gs.GetXaxis()->SetTitleSize(0.06);
        gs.SetMinimum(0); gs.SetMaximum(2.2 * avg[2]); gs.Draw("AP");
        TLine la(gs.GetXaxis()->GetXmin(), avg[2], gs.GetXaxis()->GetXmax(), avg[2]); la.SetLineStyle(2); la.SetLineColor(kBlue + 1); la.Draw();
        TLine lr(gs.GetXaxis()->GetXmin(), tot_s / tot_e, gs.GetXaxis()->GetXmax(), tot_s / tot_e); lr.SetLineColor(kRed + 1); lr.Draw();
        tl.DrawLatex(0.16, 0.88, "dashed blue = dataset average, red = this run");
        const std::string d = plotbase + Form("/pbpb%d_run%d", yr, r);
        gSystem->mkdir(d.c_str(), true);
        cl.SaveAs((d + Form("/signal_per_event_vs_lb_run%d.png", r)).c_str());
    }
}
