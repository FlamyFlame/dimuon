// Per-run Pb+Pb event-selection QA (abnormal-run investigation, docs/tracking/data_per_run_signal_yield_qa.md R2a).
//
// Applies the year's DERIVED cuts (event_sel_cuts_pbpb_20YY.root) exactly as the analysis does in
// PbPbExtras::PassEventSel (NTuple stage): bad run/LB first, then cuts 1-5 sequentially, cut 3 with the
// per-run mu+7sigma thresholds when the cuts file carries them.  Nothing is re-derived here.
// Event population = b_HLT_mu4_L1MU3V events (same as the cut derivation in plot_pbpb_event_sel_cuts.cxx).
//
//   evsel_single_run_fill(yr, infile, outfile, runs)  -- one raw-NTUP part -> per-run stage counts (ALL runs)
//                                                         + the cut-variable 2D histograms of the target runs
//   evsel_single_run_plot(yr, merged, runs, plotbase) -- hadd'ed fill outputs -> standalone PNG per cut per
//                                                         target run (the stage BEFORE that cut, cut drawn) +
//                                                         per-run pass-fraction CSV + summary printout
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLine.h"
#include "TParameter.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "../../NTupleProcessingCode/PbPbEventSelConfig.h"
#include "../event_selection/PbPbCentralityFCalMirror.h"

namespace evsr {

// Run-number axis of the per-run counters (one bin per run).  No silent default for a new year.
inline std::pair<int,int> RunRange(int yr) {
    if (yr == 23) return {460000, 464000};
    if (yr == 24) return {488000, 491000};
    if (yr == 25) return {508000, 514000};   // 2025 lumitable runs 510493-512049
    if (yr == 26) return {522000, 524000};
    throw std::runtime_error("evsel_single_run: no run range for year " + std::to_string(yr));
}
// Same rule as plot_pbpb_event_sel_cuts.cxx::UsesFCalCentralityRecompute and MuonPairPbPb.h (25, 26).
inline bool RecomputeCentrality(int yr) { return yr == 25 || yr == 26; }

// Stage index: 0 = HLT (good run/LB), k = passes cuts 1..k
constexpr int kNStage = 6;
// Views: the variables each cut acts on, filled at the stage BEFORE that cut
enum View { kZdcFcal = 0, kZdcTime, kPreampAC, kFrac, kNtrkFcal, kNView };
const int  kViewStage[kNView] = {0, 1, 2, 3, 4};  // view v filled for events passing cuts 1..v
struct Spec { int nx; double xlo, xhi; int ny; double ylo, yhi; const char* xt; const char* yt; };
// Binning = BookHists() of plot_pbpb_event_sel_cuts.cxx
const Spec kSpec[kNView] = {
    {120, -0.5, 5.5,   100,     0.,  400., "FCal E_{T}^{A+C} [TeV]",          "ZDC E_{total} [TeV]"},
    {100, -10., 10.,   100,   -10.,   10., "ZDC t^{A} [ns]",                  "ZDC t^{C} [ns]"},
    {200,-1000.,3000., 200, -1000., 3000., "ZDC preamp sum side A [ADC]",     "ZDC preamp sum side C [ADC]"},
    {100,   0., 5000., 110,     0.,   1.1, "N_{trk}^{total} (p_{T}>400 MeV)", "N_{trk}^{HItight} / N_{trk}^{total}"},
    {120, -0.5, 5.5,   100,     0., 3500., "FCal E_{T}^{A+C} [TeV]",          "N_{trk}^{HItight}"},
};
const char* kViewName[kNView] = {"zdc_fcal", "zdc_time", "preamp_AC", "ntrk_frac", "ntrk_fcal"};

struct Cuts {
    TGraph *c1 = nullptr, *c4 = nullptr, *c5lo = nullptr, *c5hi = nullptr;
    PbPbZdcTimeWindow c2{}; float c3A = 0.f, c3C = 0.f;
    std::map<int, std::pair<float,float>> c3_per_run;
    std::pair<float,float> Cut3(int run) const {
        auto it = c3_per_run.find(run);
        return it != c3_per_run.end() ? it->second : std::make_pair(c3A, c3C);
    }
};

// Same keys / same presence-driven per-run cut-3 rule as PbPbExtras::InitEventSel
inline Cuts LoadCuts(int yr) {
    const std::string path = PbPbEvSelCutsPath(yr);
    TFile* f = TFile::Open(path.c_str(), "READ");
    if (!f || f->IsZombie()) throw std::runtime_error("cannot open cuts file " + path);
    auto g = [&](const char* k) { auto* o = (TGraph*)f->Get(k); if (!o) throw std::runtime_error(std::string("missing ") + k); return (TGraph*)o->Clone(); };
    auto p = [&](const char* k) { auto* o = (TParameter<double>*)f->Get(k); if (!o) throw std::runtime_error(std::string("missing ") + k); return o->GetVal(); };
    Cuts c;
    c.c1 = g(PbPbEvSelKey::kZDCFCalCut);  c.c2 = PbPbReadZdcTimeWindow(f, yr);
    c.c3A = (float)p(PbPbEvSelKey::kPreampACutADC); c.c3C = (float)p(PbPbEvSelKey::kPreampCCutADC);
    c.c4 = g(PbPbEvSelKey::kNTrkFracCutLo);
    c.c5lo = g(PbPbEvSelKey::kNTrkFCalCutLo); c.c5hi = g(PbPbEvSelKey::kNTrkFCalCutHi);
    if (auto* t = dynamic_cast<TTree*>(f->Get(PbPbEvSelKey::kPreampPerRunTree))) {
        Int_t r = 0; Double_t a = 0., cc = 0.;
        t->SetBranchAddress("run_number", &r); t->SetBranchAddress("cut_A_ADC", &a); t->SetBranchAddress("cut_C_ADC", &cc);
        for (Long64_t i = 0; i < t->GetEntries(); ++i) { t->GetEntry(i); c.c3_per_run[r] = {(float)a, (float)cc}; }
    }
    f->Close();
    std::cout << "cuts: " << path << "  (per-run cut-3 thresholds for " << c.c3_per_run.size() << " runs)\n";
    return c;
}
} // namespace evsr

// -----------------------------------------------------------------------------------------------------------
void evsel_single_run_fill(int yr, const std::string& infile, const std::string& outfile, std::vector<int> runs,
                           Long64_t nmax = -1) {
    using namespace evsr;
    yr %= 2000;
    const Cuts cuts = LoadCuts(yr);
    const auto bad_runs = PbPbBadRuns(yr);
    const auto bad_lbs  = PbPbBadLBRanges(yr);
    const auto [r_lo, r_hi] = RunRange(yr);

    TFile fo(outfile.c_str(), "RECREATE");
    TH2D h_all("h_stage_all", ";run;stage", r_hi - r_lo, r_lo, r_hi, kNStage, 0, kNStage);
    TH2D h_c80("h_stage_c80", ";run;stage", r_hi - r_lo, r_lo, r_hi, kNStage, 0, kNStage);
    std::map<int, std::vector<TH2D*>> hv;
    std::map<int, std::pair<TH1D*,TH1D*>> hlb;   // target runs: HLT_mu4 events / events passing all cuts, per LB
    for (int r : runs) hlb[r] = {new TH1D(Form("h_lb_hlt_run%d", r), ";lumi block;events", 3000, 0, 3000),
                                 new TH1D(Form("h_lb_pass_run%d", r), ";lumi block;events", 3000, 0, 3000)};
    for (int r : runs) for (int v = 0; v < kNView; ++v) {
        const Spec& s = kSpec[v];
        hv[r].push_back(new TH2D(Form("h_%s_run%d", kViewName[v], r), Form(";%s;%s", s.xt, s.yt),
                                 s.nx, s.xlo, s.xhi, s.ny, s.ylo, s.yhi));
    }

    TChain ch("HeavyIonD3PD"); ch.Add(infile.c_str());
    ch.SetBranchStatus("*", 0);
    UInt_t run = 0, lbn = 0; Bool_t hlt = false; Int_t ctr_br = 0;
    Float_t fP = 0, fN = 0, zE[2] = {}, zT[2] = {}, pre[2][4] = {};
    std::vector<int>* nq = nullptr;
    for (auto b : {"RunNumber", "lbn", "b_HLT_mu4_L1MU3V", "centrality", "FCal_Et_P", "FCal_Et_N", "zdc_ZdcEnergy",
                   "zdc_ZdcTime", "zdc_ZdcModulePreSampleAmp", "trk_numqual"}) ch.SetBranchStatus(b, 1);
    ch.SetBranchAddress("RunNumber", &run);   ch.SetBranchAddress("lbn", &lbn);
    ch.SetBranchAddress("b_HLT_mu4_L1MU3V", &hlt); ch.SetBranchAddress("centrality", &ctr_br);
    ch.SetBranchAddress("FCal_Et_P", &fP);    ch.SetBranchAddress("FCal_Et_N", &fN);
    ch.SetBranchAddress("zdc_ZdcEnergy", zE); ch.SetBranchAddress("zdc_ZdcTime", zT);
    ch.SetBranchAddress("zdc_ZdcModulePreSampleAmp", pre); ch.SetBranchAddress("trk_numqual", &nq);

    const Long64_t n = nmax > 0 ? std::min(nmax, ch.GetEntries()) : ch.GetEntries();
    std::cout << infile << ": " << n << " entries\n";
    for (Long64_t i = 0; i < n; ++i) {
        if (i % 5000000 == 0) std::cout << "  " << i << "\n" << std::flush;
        ch.GetEntry(i);
        if (!hlt) continue;
        if (PbPbIsBadRunLB(bad_runs, bad_lbs, (int)run, lbn)) continue;   // PassEventSel: before any cut
        if ((int)run < r_lo || (int)run >= r_hi) throw std::runtime_error("run outside RunRange: " + std::to_string(run));

        const float fcal = (fP + fN) * 1e-6f;                    // MeV -> TeV
        const float ztot = (zE[0] + zE[1]) / 1000.f;             // GeV -> TeV
        const float tA = zT[1], tC = zT[0];                      // [1]=A, [0]=C
        float pA = 0.f, pC = 0.f;
        for (int k = 0; k < 4; ++k) { pA += pre[1][k]; pC += pre[0][k]; }
        const bool has_trk = nq && (int)nq->size() >= 4;
        const int ntot = has_trk ? (*nq)[0] : 0, ntight = has_trk ? (*nq)[3] : 0;
        const int ctr = RecomputeCentrality(yr) ? CentralityFromFCal2023(fcal) : ctr_br;
        const bool c80 = ctr >= 0 && ctr < 80;

        // sequential decision, identical to PbPbExtras::PassEventSel
        int stage = 0;
        auto it = hv.find((int)run);
        auto fillv = [&](int v) {
            if (it == hv.end()) return;
            TH2D* h = it->second[v];
            if (v == kZdcFcal) h->Fill(fcal, ztot);
            else if (v == kZdcTime) h->Fill(tA, tC);
            else if (v == kPreampAC) h->Fill(pA, pC);
            else if (v == kFrac) { if (ntot > 0) h->Fill(ntot, (double)ntight / ntot); }
            else if (v == kNtrkFcal) h->Fill(fcal, ntight);
        };
        do {
            fillv(kZdcFcal);
            if (ztot > (float)PbPbEvSelEvalCut(cuts.c1, fcal)) break;                                  stage = 1;
            fillv(kZdcTime);
            if (!PbPbPassZdcTime(cuts.c2, tA, tC)) break;                                              stage = 2;
            fillv(kPreampAC);
            const auto c3 = cuts.Cut3((int)run);
            if (pA >= c3.first && pC >= c3.second) break;                                              stage = 3;
            fillv(kFrac);
            if (!has_trk) break;
            if (ntot > 0 && (double)ntight / ntot < PbPbEvSelEvalCut(cuts.c4, ntot)) break;            stage = 4;
            fillv(kNtrkFcal);
            if ((double)ntight < PbPbEvSelEvalCut(cuts.c5lo, fcal) ||
                (double)ntight > PbPbEvSelEvalCut(cuts.c5hi, fcal)) break;                             stage = 5;
        } while (false);
        for (int s = 0; s <= stage; ++s) { h_all.Fill(run + 0.5, s + 0.5); if (c80) h_c80.Fill(run + 0.5, s + 0.5); }
        if (it != hv.end()) { hlb[(int)run].first->Fill(lbn + 0.5); if (stage == 5) hlb[(int)run].second->Fill(lbn + 0.5); }
    }
    fo.cd(); h_all.Write(); h_c80.Write();
    for (auto& kv : hv) for (auto* h : kv.second) h->Write();
    for (auto& kv : hlb) { kv.second.first->Write(); kv.second.second->Write(); }
    fo.Close();
    std::cout << "wrote " << outfile << "\n";
}

// -----------------------------------------------------------------------------------------------------------
void evsel_single_run_plot(int yr, const std::string& merged, std::vector<int> runs, const std::string& plotbase,
                           const std::string& csv) {
    using namespace evsr;
    yr %= 2000;
    gStyle->SetOptStat(0); gStyle->SetPalette(kBird);
    const Cuts cuts = LoadCuts(yr);
    TFile* f = TFile::Open(merged.c_str());
    auto* h_all = (TH2D*)f->Get("h_stage_all");
    auto* h_c80 = (TH2D*)f->Get("h_stage_c80");

    // per-run pass fractions, ALL runs (context for the outliers)
    std::ofstream o(csv);
    o << "run,n_hlt,n_c1,n_c2,n_c3,n_c4,n_c5,frac_all,n_hlt_c80,n_c5_c80,frac_c80\n";
    double S0 = 0, S5 = 0, S0c = 0, S5c = 0;
    for (int b = 1; b <= h_all->GetNbinsX(); ++b) {
        const double n0 = h_all->GetBinContent(b, 1);
        if (n0 <= 0) continue;
        const int r = (int)h_all->GetXaxis()->GetBinLowEdge(b);
        o << r;
        for (int s = 1; s <= kNStage; ++s) o << "," << (long)h_all->GetBinContent(b, s);
        const double n5 = h_all->GetBinContent(b, kNStage), c0 = h_c80->GetBinContent(b, 1), c5 = h_c80->GetBinContent(b, kNStage);
        o << "," << n5 / n0 << "," << (long)c0 << "," << (long)c5 << "," << (c0 > 0 ? c5 / c0 : 0.) << "\n";
        S0 += n0; S5 += n5; S0c += c0; S5c += c5;
    }
    std::cout << Form("Pb+Pb 20%d ALL runs: pass all evsel = %.3f %% (HLT_mu4), %.3f %% (0-80 %%)\n", yr, 100 * S5 / S0, 100 * S5c / S0c);

    for (int r : runs) {
        const int b = h_all->GetXaxis()->FindBin(r + 0.5);
        const double n0 = h_all->GetBinContent(b, 1);
        std::cout << Form("run %d: HLT_mu4 events %.0f;", r, n0);
        for (int s = 2; s <= kNStage; ++s) std::cout << Form(" c%d %.3f%%", s - 1, 100 * h_all->GetBinContent(b, s) / n0);
        std::cout << Form("  | 0-80%%: %.3f%%\n", 100 * h_c80->GetBinContent(b, kNStage) / h_c80->GetBinContent(b, 1));

        const std::string dir = plotbase + Form("/pbpb%d_run%d/evsel", yr, r);
        gSystem->mkdir(dir.c_str(), true);
        auto label = [&](TLatex& tl, int v) {
            tl.SetNDC(); tl.SetTextSize(0.034);
            tl.DrawLatex(0.15, 0.92, Form("Pb+Pb 20%d data, run %d, HLT_mu4", yr, r));
            const double frac = h_all->GetBinContent(b, kViewStage[v] + 1) / n0;
            tl.DrawLatex(0.15, 0.87, v == 0 ? "before event selection"
                                            : Form("after cuts 1-%d (%.2f %% of events)", v, 100 * frac));
        };
        for (int v = 0; v < kNView; ++v) {
            auto* h = (TH2D*)f->Get(Form("h_%s_run%d", kViewName[v], r));
            const int cut = v + 1;
            if (v == kPreampAC) {   // cut 3: 1D per side, both-sides-fail rule, run's own threshold
                const auto c3 = cuts.Cut3(r);
                TCanvas c(Form("c3_%d", r), "", 1200, 600); c.Divide(2, 1, 0.005, 0.005);
                TH1D* p[2] = {h->ProjectionX(Form("pA_%d", r)), h->ProjectionY(Form("pC_%d", r))};
                const float cv[2] = {c3.first, c3.second}; const char* side[2] = {"A", "C"};
                for (int i = 0; i < 2; ++i) {
                    c.cd(i + 1); gPad->SetLogy(); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
                    p[i]->Rebin(3); p[i]->SetMarkerStyle(20); p[i]->SetMarkerSize(0.5);
                    p[i]->GetYaxis()->SetTitle("Events / 60 ADC");
                    const double ymax = p[i]->GetBinContent(p[i]->GetMaximumBin()) * 5.;
                    p[i]->GetYaxis()->SetRangeUser(0.5, ymax); p[i]->Draw("E");
                    TLine l(cv[i], 0.5, cv[i], ymax); l.SetLineColor(kRed); l.SetLineWidth(2); l.SetLineStyle(2); l.DrawClone();
                    TLatex tl; label(tl, v);
                    tl.SetTextColor(kRed); tl.DrawLatex(0.15, 0.82, Form("side %s cut = %.0f ADC (event fails if A AND C above)", side[i], cv[i]));
                }
                c.SaveAs((dir + Form("/evsel_cut%d_%s_run%d.png", cut, kPbPbEvSelCutLabel(cut), r)).c_str());
                continue;
            }
            TCanvas c(Form("c%d_%d", v, r), "", 800, 650);
            c.SetLogz(); c.SetRightMargin(0.14); c.SetTopMargin(0.15);
            h->SetContour(99); h->Draw("COLZ");
            TLatex tl; label(tl, v); tl.SetTextColor(kRed);
            auto drawg = [](TGraph* g) { auto* d = (TGraph*)g->Clone(); d->SetLineColor(kRed); d->SetLineWidth(2); d->Draw("L SAME"); };
            if (v == kZdcFcal) { drawg(cuts.c1); tl.DrawLatex(0.15, 0.80, "ZDC-FCal cut (year)"); }
            else if (v == kZdcTime) {
                TLine l; l.SetLineColor(kRed); l.SetLineWidth(2); l.SetLineStyle(2); const auto& w = cuts.c2;
                l.DrawLine(w.a_lo, w.c_lo, w.a_lo, w.c_hi); l.DrawLine(w.a_hi, w.c_lo, w.a_hi, w.c_hi);
                l.DrawLine(w.a_lo, w.c_lo, w.a_hi, w.c_lo); l.DrawLine(w.a_lo, w.c_hi, w.a_hi, w.c_hi);
                tl.DrawLatex(0.15, 0.80, Form("%.1f < t^{A} < %.1f ns, %.1f < t^{C} < %.1f ns", w.a_lo, w.a_hi, w.c_lo, w.c_hi));
            }
            else if (v == kFrac) { drawg(cuts.c4); tl.DrawLatex(0.15, 0.80, "N_{trk} fraction lower cut (year)"); }
            else if (v == kNtrkFcal) { drawg(cuts.c5lo); drawg(cuts.c5hi); tl.DrawLatex(0.15, 0.80, "N_{trk}-FCal band (year)"); }
            c.SaveAs((dir + Form("/evsel_cut%d_%s_run%d.png", cut, kPbPbEvSelCutLabel(cut), r)).c_str());
        }
    }
}
