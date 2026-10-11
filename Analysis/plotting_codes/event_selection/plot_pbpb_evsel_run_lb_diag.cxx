// Draw stage of the Pb+Pb event-selection run / LB / ZDC-timing diagnostics.  Reads the hadd-ed
// output of fill_pbpb_evsel_run_lb_diag.cxx (see that file for the event selection, which mirrors
// PbPbExtras::PassEventSel).  Evidence for the DQ experts only -- no cut or run list is changed.
//
//   (A) cut-4 bad-run hunt (year 26): per-run cut-4 failure fraction (sample = after cuts 1-3, as
//       event_sel_cut4_nTrk_frac_standalone) next to the per-run cut-5 failure fraction (after
//       cuts 1-4) and the low N_trk/FCal fraction; per-run panels of the cut-4 figure.
//   (B) LB scan of the runs stored in t_lbscan (year 26): see DrawLBScan.
//   (C) ZDC t_A vs t_C and the PV0 z of its clusters (any year).
//
// Usage:  root -l -b -q 'plot_pbpb_evsel_run_lb_diag.cxx+(26)'
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <functional>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TKey.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TBox.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include "../../NTupleProcessingCode/PbPbEventSelConfig.h"

namespace {
const std::string kBase = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/";

// Low N_trk^HItight/FCal classification used for the cut-5 low-band runs
// (plot_pbpb_ntrk_lowband_runs.cxx: low-band runs peak at ~470 TeV^-1, normal runs at ~590).
const double kRatioLow = 530.;   // TeV^-1, events with FCal > 1 TeV

std::vector<int> RunsInFile(TFile* f) {
    std::vector<int> runs;
    for (TObject* o : *f->GetListOfKeys()) {
        std::string n = o->GetName();
        if (n.rfind("h_cnt_", 0) == 0) runs.push_back(std::stoi(n.substr(6)));
    }
    std::sort(runs.begin(), runs.end());
    return runs;
}
template <class T> T* Need(TFile* f, const std::string& k) {
    T* o = dynamic_cast<T*>(f->Get(k.c_str()));
    if (!o) throw std::runtime_error("missing " + k + " in " + f->GetName());
    return o;
}
// layout rule: nrows >= ncols, nrows ~ sqrt(N); N <= 3 -> one row
void Grid(int n, int& nc, int& nr) {
    if (n <= 3) { nc = n; nr = 1; return; }
    nc = std::max(1, (int)std::floor(std::sqrt(n)));
    nr = (n + nc - 1) / nc;
}
}  // namespace

// ============================================================================ (A)
void DrawCut4Hunt(TFile* f, TGraph* g4, const std::string& out, int yr)
{
    const std::vector<int> runs = RunsInFile(f);
    struct Rec { int run; double n13, f4, n14, f5, nhi, flow; };
    std::vector<Rec> recs;
    for (int r : runs) {
        TH1D* c = Need<TH1D>(f, Form("h_cnt_%d", r));
        TH1D* h = Need<TH1D>(f, Form("h_ratio_nocut_%d", r));
        Rec x{r, c->GetBinContent(2), 0., c->GetBinContent(4), 0., h->Integral(), 0.};
        x.f4   = x.n13 > 0 ? c->GetBinContent(3) / x.n13 : 0.;
        x.f5   = x.n14 > 0 ? c->GetBinContent(5) / x.n14 : 0.;
        x.flow = x.nhi > 0 ? h->Integral(1, h->FindBin(kRatioLow) - 1) / x.nhi : 0.;
        recs.push_back(x);
    }

    // per-run fractions, run index on x
    const int nr = recs.size();
    TH1D h4("h_f4", "", nr, 0, nr), h5("h_f5", "", nr, 0, nr), hl("h_fl", "", nr, 0, nr);
    for (int i = 0; i < nr; ++i) {
        for (TH1D* h : {&h4, &h5, &hl}) h->GetXaxis()->SetBinLabel(i + 1, Form("%d", recs[i].run));
        h4.SetBinContent(i + 1, 100. * recs[i].f4);
        h4.SetBinError(i + 1, 100. * std::sqrt(recs[i].f4 * (1 - recs[i].f4) / std::max(1., recs[i].n13)));
        h5.SetBinContent(i + 1, 100. * recs[i].f5);
        h5.SetBinError(i + 1, 100. * std::sqrt(recs[i].f5 * (1 - recs[i].f5) / std::max(1., recs[i].n14)));
        hl.SetBinContent(i + 1, 100. * recs[i].flow);
        hl.SetBinError(i + 1, 100. * std::sqrt(recs[i].flow * (1 - recs[i].flow) / std::max(1., recs[i].nhi)));
    }
    {
        TCanvas c("c_frac_run", "", 1400, 1500);
        c.Divide(1, 3, 0.001, 0.001);
        struct P { TH1D* h; const char* yt; int col; };
        const P ps[3] = {{&h4, "cut 4 failure [%]", kBlue + 1}, {&h5, "cut 5 failure [%]", kRed + 1}, {&hl, "f_{low} [%]", kGreen + 2}};
        const char* defs[3] = {"cut 4: N_{trk}^{HItight}/N_{trk}^{total} < #mu#minus5#sigma(N_{trk}^{total}); events after cuts 1-3",
                               "cut 5: N_{trk}^{HItight} outside #mu#pm5#sigma band vs FCal E_{T}^{A+C}; events after cuts 1-4",
                               Form("f_{low}: N_{trk}^{HItight}/FCal E_{T}^{A+C} < %.0f TeV^{-1}, FCal > 1 TeV; no cut", kRatioLow)};
        TLatex td; td.SetNDC(); td.SetTextSize(0.05);
        for (int ip = 0; ip < 3; ++ip) {
            c.cd(ip + 1);
            gPad->SetLogy();  // log y: per-run fractions span 0.05 % (bulk) to ~90 % (522384/522408)
            gPad->SetBottomMargin(0.2); gPad->SetLeftMargin(0.08); gPad->SetRightMargin(0.02);
            gPad->SetGridy();
            TH1D* h = ps[ip].h;
            h->SetTitle(Form("Pb+Pb 20%d;;%s", yr, ps[ip].yt));
            h->SetMarkerStyle(20); h->SetMarkerColor(ps[ip].col); h->SetLineColor(ps[ip].col);
            h->SetMinimum(1e-3); h->SetMaximum(1e4);  // headroom above 100 % keeps the definition text clear of the points
            h->GetXaxis()->LabelsOption("v"); h->GetXaxis()->SetLabelSize(0.06);
            h->GetYaxis()->SetTitleSize(0.05); h->GetYaxis()->SetTitleOffset(0.7); h->GetYaxis()->SetLabelSize(0.05);
            h->Draw("E1P");
            td.DrawLatex(0.1, 0.83, defs[ip]);
        }
        const std::string png = out + "/cut4_cut5_failure_per_run_pbpb_20" + std::to_string(yr) + ".png";
        c.SaveAs(png.c_str());
    }

    // panels: runs with cut-4 failure > 3x the median run, plus the cut-5 low-band runs, + all others
    std::vector<double> f4s; for (auto& x : recs) f4s.push_back(x.f4);
    std::nth_element(f4s.begin(), f4s.begin() + f4s.size() / 2, f4s.end());
    const double med = f4s[f4s.size() / 2];
    std::vector<int> pan;
    for (auto& x : recs) if (x.f4 > 3. * med || x.flow > 0.2) pan.push_back(x.run);
    if (pan.size() > 11) {  // keep the 11 largest cut-4 failure fractions
        std::sort(pan.begin(), pan.end(), [&](int a, int b) {
            auto fa = std::find_if(recs.begin(), recs.end(), [&](const Rec& r) { return r.run == a; })->f4;
            auto fb = std::find_if(recs.begin(), recs.end(), [&](const Rec& r) { return r.run == b; })->f4;
            return fa > fb; });
        pan.resize(11);
        std::sort(pan.begin(), pan.end());
    }
    TH2D* h_other = nullptr;
    for (int r : runs) {
        if (std::find(pan.begin(), pan.end(), r) != pan.end()) continue;
        TH2D* h = Need<TH2D>(f, Form("h_frac_c13_%d", r));
        if (!h_other) { h_other = (TH2D*)h->Clone("h_frac_c13_other"); h_other->SetDirectory(nullptr); }
        else h_other->Add(h);
    }
    std::vector<std::pair<TH2D*, std::string>> hs;
    for (int r : pan) {
        auto it = std::find_if(recs.begin(), recs.end(), [&](const Rec& x) { return x.run == r; });
        hs.push_back({Need<TH2D>(f, Form("h_frac_c13_%d", r)),
                      Form("run %d: cut 4 fail %.2f%%", r, 100. * it->f4)});
    }
    double o13 = 0, o4 = 0;
    for (auto& x : recs) if (std::find(pan.begin(), pan.end(), x.run) == pan.end()) { o13 += x.n13; o4 += x.f4 * x.n13; }
    hs.push_back({h_other, Form("all other runs: cut 4 fail %.2f%%", o13 > 0 ? 100. * o4 / o13 : 0.)});
    int nc, nrw; Grid(hs.size(), nc, nrw);
    TCanvas c("c_c4_runs", "", 560 * nc, 460 * nrw);
    c.Divide(nc, nrw, 0.002, 0.002);
    g4->SetLineColor(kRed); g4->SetLineWidth(2);
    TLatex tl; tl.SetNDC(); tl.SetTextSize(0.05);
    for (size_t i = 0; i < hs.size(); ++i) {
        c.cd(i + 1);
        gPad->SetLogz(); gPad->SetRightMargin(0.14); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);  // log z as cut-4 figure
        hs[i].first->SetTitle("");
        hs[i].first->SetContour(99);
        hs[i].first->Draw("COLZ");
        g4->Draw("L SAME");
        tl.DrawLatex(0.15, 0.92, Form("Pb+Pb 20%d, %s", yr, hs[i].second.c_str()));
        tl.SetTextAlign(31); tl.DrawLatex(0.84, 0.84, "red: cut 4 (#mu#minus5#sigma); events after cuts 1-3"); tl.SetTextAlign(11);
    }
    const std::string png = out + "/cut4_nTrk_frac_per_run_pbpb_20" + std::to_string(yr) + ".png";
    c.SaveAs(png.c_str());

    // table
    std::ofstream o(out + "/cut4_cut5_per_run_pbpb_20" + std::to_string(yr) + ".txt");
    o << "Pb+Pb 20" << yr << ", HLT_mu4_L1MU3V events. Cuts as PbPbExtras::PassEventSel (cuts file " << PbPbEvSelCutsPath(yr) << ")\n"
      << "N13 = pass cuts 1-3; f4 = fraction of N13 failing cut 4; N14 = pass cuts 1-4; f5 = fraction of N14 failing cut 5;\n"
      << "f_low = fraction of all HLT events with FCal > 1 TeV having N_trk^HItight/FCal < " << kRatioLow << " TeV^-1 (no cut)\n\n"
      << "N_ovf = events of N13 with N_trk^total > 5000 (beyond the cut-4 figure axis; cut-4 curve defined for 100-4900,\n"
      << "        clamped to its last point above); N_ovf_lo = of those, N_trk^HItight/N_trk^total < 0.45\n\n"
      << Form("%-8s %11s %9s %11s %9s %9s %9s %9s\n", "run", "N13", "f4[%]", "N14", "f5[%]", "f_low[%]", "N_ovf", "N_ovf_lo");
    for (auto& x : recs)
    {
        TH2D* hf = Need<TH2D>(f, Form("h_frac_c13_%d", x.run));
        const int ox = hf->GetNbinsX() + 1;
        const double novf = hf->Integral(ox, ox, 0, hf->GetNbinsY() + 1);
        const double nlo  = hf->Integral(ox, ox, 0, hf->GetYaxis()->FindBin(0.45 - 1e-6));
        o << Form("%-8d %11.0f %9.4f %11.0f %9.4f %9.3f %9.0f %9.0f\n", x.run, x.n13, 100 * x.f4, x.n14, 100 * x.f5, 100 * x.flow, novf, nlo);
    }
    o << "\nmedian per-run f4 = " << 100 * med << " %; panel runs (f4 > 3x median or f_low > 20%):";
    for (int r : pan) o << " " << r;
    o << "\n";
}

// ============================================================================ (C)
// ZDC t_A vs t_C (HLT_mu4, no event-selection cut, axes as ZDC_t_A_vs_C) and PV0 z per cluster.
// The analysis' ONLY vertex-z requirement is the skim's |z(PV0)| < 250 mm (TrigRates.cxx
// CleaningCuts, MaxZvtx = 250 in TrigRates_CA.py) -- every event here has already passed it, so
// a "before" sample does not exist on disk.  The right pad shows instead where each cluster sits in
// z(PV0) relative to that +-250 mm boundary: a cluster the cut could remove must populate |z| ~ 250.
// NOTE: z(PV0) in the skim is truncated UPSTREAM at |z| ~ 205 mm (2026: 7 events in [-220,-210) vs
// 2822 in [-210,-200), 1 in [200,210) vs 2658 in [190,200), plus 1 lone event in [240,250); same in 2025) -- a reconstruction-level
// z acceptance of the PV0 actually stored, not an analysis cut (source not identifiable from the skim
// content).  The skim's 250 mm requirement therefore never acts, and any population at |z| > ~205 mm
// is absent before it; a true "before" needs a re-skim from AOD with MaxZvtx < 0.
struct TBoxDef { const char* name; double a0, a1, c0, c1; int col; };
void DrawZdcVtx(TFile* f, const std::string& out, int yr, const std::vector<TBoxDef>& boxes)
{
    TH3F* h3 = Need<TH3F>(f, "h3_zdct_vtxz");
    TH2D* h2 = (TH2D*)h3->Project3D("yx");   // x = t_A, y = t_C
    h2->SetDirectory(nullptr);
    h2->SetTitle(";ZDC t^{A} [ns];ZDC t^{C} [ns]");
    TCanvas c("c_zdc_vtx", "", 1400, 600);
    c.Divide(2, 1, 0.002, 0.002);
    c.cd(1);
    gPad->SetLogz(); gPad->SetRightMargin(0.13); gPad->SetLeftMargin(0.12); gPad->SetTopMargin(0.08);  // log z as ZDC_t_A_vs_C (~7 decades)
    h2->SetContour(99); h2->Draw("COLZ");
    TLatex tl; tl.SetNDC(); tl.SetTextSize(0.038);
    tl.SetTextSize(0.034);
    tl.DrawLatex(0.12, 0.94, Form("Pb+Pb 20%d, HLT_mu4, no event selection", yr));  // in the top margin: no empty region inside the frame
    std::vector<TH1D*> hz;
    for (auto& b : boxes) {
        TBox* bx = new TBox(b.a0, b.c0, b.a1, b.c1);
        bx->SetFillStyle(0); bx->SetLineColor(b.col); bx->SetLineWidth(2); bx->Draw("same");
        auto ax = h3->GetXaxis(); auto ay = h3->GetYaxis();
        TH1D* h = h3->ProjectionZ(Form("hz_%s", b.name), ax->FindBin(b.a0 + 1e-4), ax->FindBin(b.a1 - 1e-4),
                                  ay->FindBin(b.c0 + 1e-4), ay->FindBin(b.c1 - 1e-4));
        h->SetDirectory(nullptr);
        hz.push_back(h);
    }
    c.cd(2);
    gPad->SetLogy(); gPad->SetLeftMargin(0.12); gPad->SetRightMargin(0.04);  // log y: shows the tails toward +-250 mm
    TLegend leg(0.14, 0.68, 0.7, 0.9); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.028);
    double ymax = 0, ymin = 1e9;
    for (size_t i = 0; i < hz.size(); ++i) {
        const double n = hz[i]->Integral();
        if (n > 0) hz[i]->Scale(1. / n);
        for (int k = 1; k <= hz[i]->GetNbinsX(); ++k) if (hz[i]->GetBinContent(k) > 0) ymin = std::min(ymin, hz[i]->GetBinContent(k));
        ymax = std::max(ymax, hz[i]->GetBinContent(hz[i]->GetMaximumBin()));
        hz[i]->SetLineColor(boxes[i].col); hz[i]->SetMarkerColor(boxes[i].col); hz[i]->SetLineWidth(2);
        leg.AddEntry(hz[i], Form("%s: N = %.3g, #LTz#GT = %.0f mm, RMS %.0f mm", boxes[i].name, n,
                                 hz[i]->GetMean(), hz[i]->GetRMS()), "l");
    }
    for (size_t i = 0; i < hz.size(); ++i) {
        hz[i]->SetTitle(";z_{vtx}(PV0) [mm];normalised events");
        hz[i]->SetMinimum(0.5 * ymin); hz[i]->SetMaximum(2000 * ymax);
        hz[i]->Draw(i ? "HIST SAME" : "HIST");
    }
    TLine l; l.SetLineStyle(2); l.SetLineColor(kBlack); l.SetLineWidth(2);
    l.DrawLine(-250, 0.5 * ymin, -250, 10 * ymax); l.DrawLine(250, 0.5 * ymin, 250, 10 * ymax);  // stop below the legend
    leg.AddEntry((TObject*)nullptr, "dashed: skim requirement |z_{vtx}(PV0)| < 250 mm", "");
    leg.Draw();
    const std::string png = out + "/ZDC_t_A_vs_C_and_vtx_z_pbpb_20" + std::to_string(yr) + ".png";
    c.SaveAs(png.c_str());
    for (size_t i = 0; i < hz.size(); ++i)
        std::cout << Form("[ZDC] %-40s N=%10.0f  <z>=%6.1f mm  RMS=%5.1f\n", boxes[i].name, hz[i]->GetEntries(), hz[i]->GetMean(),
                          hz[i]->GetRMS());
}

// ============================================================================ (B)
// LB localisation of the low-N_trk^HItight/FCal shoulder in one run (t_lbscan, all HLT_mu4 events
// of the run, no cut).  Per LB: f_low = N(ratio < kRatioLow)/N(FCal > 1 TeV).  Reference f_ref =
// the same fraction in all runs except the scanned ones and the low-band runs (f_low > 20 %).
// Bisection: a LB range is CLEAN if its f_low is within kZ binomial sigma of f_ref (or below);
// otherwise it is split at its median LB and both halves are tested again, down to single LBs; a
// single significant LB (> kZ sigma) is "second band" if f_low > kSecondBand, "single-band drift"
// if f_low > kDrift x f_ref, else counted good.  The bisection tree is the printed search record; the
// classification itself is then applied to EVERY LB (a pruned clean range can hide a weak LB).  Same-category LBs separated by <= kMaxGap other
// LBs are merged into one block (the gap LBs join it).
// Contiguous BAD single LBs (no clean LB between them) form a block.
// Cut-4 cross-check: per-LB cut-4 failure fraction among events passing cuts 1-3.
namespace {
const double kZ = 5.;
const double kSecondBand = 0.5;  // f_low > 50 %: the low band holds the majority of the LB's events
const double kDrift = 2.0;       // elevated single band: f_low > 2 x f_ref (with ~3e4 events/LB, 5 sigma
                                 // alone flags 0.4 % absolute drifts, f_low 2.3 vs 1.9 %)
const int    kMaxGap = 3;        // <= 3 unclassified LBs between two LBs of one category are absorbed
struct LBC { double n = 0, nlow = 0, n13 = 0, n4f = 0; };
}

void DrawLBScan(TFile* f, int run, const std::vector<int>& ref_excl, TGraph* g4, const std::string& out, int yr,
                const std::map<int, std::pair<int, std::string>>& lb_marks,
                std::ofstream& o)
{
    // ---- reference f_low ----
    double rn = 0, rl = 0;
    for (int r : RunsInFile(f)) {
        if (std::find(ref_excl.begin(), ref_excl.end(), r) != ref_excl.end()) continue;
        TH1D* h = Need<TH1D>(f, Form("h_ratio_nocut_%d", r));
        const double n = h->Integral(), l = h->Integral(1, h->FindBin(kRatioLow) - 1);
        if (n > 0 && l / n > 0.2) continue;  // low-band run (522384/522408 type)
        rn += n; rl += l;
    }
    const double fref = rl / rn;

    TTree* t = Need<TTree>(f, "t_lbscan");
    Int_t b_run, b_lbn, b_ntot, b_ntight, b_pass; Float_t b_fcal;
    t->SetBranchAddress("run", &b_run); t->SetBranchAddress("lbn", &b_lbn);
    t->SetBranchAddress("ntrk_total", &b_ntot); t->SetBranchAddress("ntrk_tight", &b_ntight);
    t->SetBranchAddress("pass_bits", &b_pass); t->SetBranchAddress("fcal", &b_fcal);
    auto is_low = [&]() { return b_fcal > 1.f && b_ntight / b_fcal < kRatioLow; };
    auto c4fail = [&]() { return (b_pass & 7) == 7 && !(b_pass & 8); };

    std::map<int, LBC> lb;
    const Long64_t ne = t->GetEntries();
    for (Long64_t i = 0; i < ne; ++i) {
        t->GetEntry(i);
        if (b_run != run) continue;
        LBC& x = lb[b_lbn];
        if (b_fcal > 1.f) { x.n += 1; if (is_low()) x.nlow += 1; }
        if ((b_pass & 7) == 7) { x.n13 += 1; if (c4fail()) x.n4f += 1; }
    }
    if (lb.empty()) throw std::runtime_error(Form("run %d not in t_lbscan", run));
    std::vector<int> lbs; for (auto& kv : lb) lbs.push_back(kv.first);

    // ---- bisection ----
    std::map<int, int> cat;  // LB -> 1 second band, 2 single-band drift (absent = good)
    o << "\n==== run " << run << ": " << lbs.size() << " LBs [" << lbs.front() << ", " << lbs.back() << "], f_ref = "
      << 100 * fref << " %  (f_low: ratio < " << kRatioLow << " TeV^-1, FCal > 1 TeV)\nBisection tree (a range is split while "
      << "f_low > f_ref + " << kZ << " sigma; leaves: SECOND-BAND f_low > " << 100 * kSecondBand << " %, DRIFT f_low > "
      << kDrift << " x f_ref):\n";
    std::function<void(size_t, size_t, int)> bis = [&](size_t i0, size_t i1, int depth) {
        double n = 0, l = 0;
        for (size_t i = i0; i < i1; ++i) { n += lb[lbs[i]].n; l += lb[lbs[i]].nlow; }
        const double fl = n > 0 ? l / n : 0., z = n > 0 ? (fl - fref) / std::sqrt(fref * (1 - fref) / n) : 0.;
        const bool sig = z > kZ, leaf = (i1 - i0 == 1);
        const int c = !(sig && leaf) ? 0 : (fl > kSecondBand ? 1 : (fl > kDrift * fref ? 2 : 0));
        // print the tree down to ranges of >= 8 LBs, plus every classified LB (all LBs: per-LB plot)
        if (i1 - i0 >= 8 || c)
            o << std::string(2 * depth, ' ') << Form("LB %4d-%-4d  N=%8.0f  f_low=%6.2f%%  z=%7.1f  %s\n", lbs[i0], lbs[i1 - 1], n,
                 100 * fl, z, !sig ? "clean" : (!leaf ? "-> split" : (c == 1 ? "SECOND-BAND" : (c == 2 ? "DRIFT" : "minor"))));
        if (!sig) return;
        if (leaf) { if (c) cat[lbs[i0]] = c; return; }
        const size_t m = (i0 + i1) / 2;
        bis(i0, m, depth + 1); bis(m, i1, depth + 1);
    };
    bis(0, lbs.size(), 0);
    for (int l : lbs) {
        const double n = lb[l].n, fl = n > 0 ? lb[l].nlow / n : 0.;
        if (n > 0 && (fl - fref) / std::sqrt(fref * (1 - fref) / n) > kZ)
            if (const int c = fl > kSecondBand ? 1 : (fl > kDrift * fref ? 2 : 0)) cat[l] = c;
    }
    // absorb short gaps between same-category LBs, then form blocks over the run's LB list
    std::vector<int> cv(lbs.size(), 0);
    for (size_t i = 0; i < lbs.size(); ++i) if (cat.count(lbs[i])) cv[i] = cat[lbs[i]];
    for (size_t i = 0; i < lbs.size(); ++i) {
        if (!cv[i]) continue;
        for (size_t j = i + 1; j < lbs.size() && j <= i + kMaxGap + 1; ++j) {
            if (cv[j] == cv[i]) { for (size_t k = i + 1; k < j; ++k) if (!cv[k]) cv[k] = cv[i]; break; }
            if (cv[j]) break;
        }
    }
    struct Block { int first, last, cat; };
    std::vector<Block> blocks;
    for (size_t i = 0; i < lbs.size(); ++i) {
        if (!cv[i]) continue;
        if (!blocks.empty() && i > 0 && cv[i - 1] == cv[i]) blocks.back().last = lbs[i];
        else blocks.push_back({lbs[i], lbs[i], cv[i]});
    }
    auto block_of = [&](int l) -> int {
        for (size_t k = 0; k < blocks.size(); ++k) if (l >= blocks[k].first && l <= blocks[k].last) return k;
        return -1;
    };
    auto lbrange = [](int a, int b) { return a == b ? std::string(Form("LB %d", a)) : std::string(Form("LB %d-%d", a, b)); };
    auto catname = [](int c) { return c == 1 ? "second band" : "drifting"; };

    // ---- summary per block ----
    o << "Blocks (" << blocks.size() << "):\n";
    double gn = 0, gl = 0, g13 = 0, g4f = 0;
    for (int l : lbs) if (block_of(l) < 0) { gn += lb[l].n; gl += lb[l].nlow; g13 += lb[l].n13; g4f += lb[l].n4f; }
    for (auto& b : blocks) {
        double n = 0, l = 0, n13 = 0, n4 = 0;
        for (int x : lbs) if (x >= b.first && x <= b.last) { n += lb[x].n; l += lb[x].nlow; n13 += lb[x].n13; n4 += lb[x].n4f; }
        o << Form("  %-12s %-17s N(FCal>1)=%8.0f  f_low=%6.2f%%   cut-4 fail (after 1-3) = %.4f%%\n", lbrange(b.first, b.last).c_str(),
                  catname(b.cat), n, 100 * l / n, n13 > 0 ? 100 * n4 / n13 : 0.);
    }
    o << Form("  %-12s %-17s N(FCal>1)=%8.0f  f_low=%6.2f%%   cut-4 fail (after 1-3) = %.4f%%\n", "other LBs", "good", gn,
              100 * gl / gn, g13 > 0 ? 100 * g4f / g13 : 0.);
    {
        const double r4 = g13 > 0 ? g4f / g13 : 0.;
        o << Form("LBs with cut-4 failure > %.0f sigma above the run's other-LB rate (%.4f %%), and LBs with < 500 events after cuts 1-3 next to a block:\n", kZ, 100 * r4);
        for (size_t i = 0; i < lbs.size(); ++i) {
            const LBC& x = lb[lbs[i]];
            if (x.n13 <= 0) continue;
            const double f4 = x.n4f / x.n13, z4 = (f4 - r4) / std::sqrt(std::max(r4 * (1 - r4), 1e-12) / x.n13);
            const bool edge = x.n13 < 500 && ((i > 0 && block_of(lbs[i - 1]) != block_of(lbs[i])) || block_of(lbs[i]) >= 0 ||
                                              (i + 1 < lbs.size() && block_of(lbs[i + 1]) >= 0));
            if (z4 > kZ || edge)
                o << Form("  LB %4d  N13=%7.0f  cut-4 fail=%7.0f (%7.3f%%)  f_low=%6.2f%%  %s\n", lbs[i], x.n13, x.n4f, 100 * f4,
                          x.n > 0 ? 100 * x.nlow / x.n : 0., block_of(lbs[i]) >= 0 ? catname(blocks[block_of(lbs[i])].cat) : "other");
        }
    }
    o << "  (block edges of the DRIFT category depend on the kDrift threshold; the SECOND-BAND edges are sharp steps)\n";

    // ---- histograms: bad (all blocks) vs good ----
    const int nb = blocks.size();
    std::vector<TH1D*> hr(nb + 1); std::vector<TH2D*> hnf(2), hfr(2);
    for (int k = 0; k <= nb; ++k) {
        hr[k] = new TH1D(Form("hr_%d_%d", run, k), "", 100, 0., 1000.); hr[k]->SetDirectory(nullptr);
    }
    for (int k = 0; k < 2; ++k) {
        hnf[k] = new TH2D(Form("hnf_%d_%d", run, k), ";FCal E_{T}^{A+C} [TeV];N_{trk}^{HItight}", 120, -0.5, 5.5, 100, 0., 3500.);
        hfr[k] = new TH2D(Form("hfr_%d_%d", run, k), ";N_{trk}^{total} (p_{T}>400 MeV);N_{trk}^{HItight} / N_{trk}^{total}",
                          100, 0., 5000., 110, 0., 1.1);
        hnf[k]->SetDirectory(nullptr); hfr[k]->SetDirectory(nullptr);
    }
    for (Long64_t i = 0; i < ne; ++i) {
        t->GetEntry(i);
        if (b_run != run) continue;
        const int k = block_of(b_lbn);
        const int g = k < 0 ? 0 : (blocks[k].cat == 1 ? 1 : -1);  // 2D: 0 = good LBs, 1 = second band
        if (b_fcal > 1.f) hr[k < 0 ? nb : k]->Fill(b_ntight / b_fcal);
        if (g < 0) continue;
        hnf[g]->Fill(b_fcal, b_ntight);
        if ((b_pass & 7) == 7 && b_ntot > 0) hfr[g]->Fill(b_ntot, (double)b_ntight / b_ntot);
    }

    // ---- page 1: per-LB scan + ratio 1D ----
    {
        TGraphErrors gl_lb;  // binomial errors: edge LBs hold only tens of events
        for (int l : lbs) {
            if (lb[l].n <= 0) continue;
            const double n = lb[l].n, fl = lb[l].nlow / n; const int i = gl_lb.GetN();
            gl_lb.SetPoint(i, l, 100 * fl); gl_lb.SetPointError(i, 0, 100 * std::sqrt(std::max(fl * (1 - fl), 1. / n) / n));
        }
        TCanvas c(Form("c_lb_%d", run), "", 1300, 1000);
        c.Divide(1, 2, 0.002, 0.002);
        TLatex tl; tl.SetNDC(); tl.SetTextSize(0.045);
        // (1) per LB.  log y: f_low runs from ~1-3 % (normal LBs) to ~70-100 % (second band)
        c.cd(1); gPad->SetLeftMargin(0.08); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.25); gPad->SetLogy();
        const double ylo = 0.5, yhi = 200.;
        TH1F* fr = gPad->DrawFrame(lbs.front() - 2, ylo, lbs.back() + 2, yhi, ";lumi block;f_{low} [%]");
        fr->GetYaxis()->SetTitleOffset(0.8);
        const int bcol[3] = {0, kOrange, kYellow - 9};
        for (auto& b : blocks) {
            TBox* bx = new TBox(b.first - 0.5, ylo, b.last + 0.5, yhi);
            bx->SetFillColorAlpha(bcol[b.cat], b.cat == 1 ? 0.35 : 0.6); bx->SetLineWidth(0); bx->Draw();
        }
        gl_lb.SetMarkerStyle(20); gl_lb.SetMarkerSize(0.6); gl_lb.SetMarkerColor(kGreen + 2); gl_lb.Draw("P");
        TLine lr; lr.SetLineStyle(2); lr.SetLineColor(kGreen + 2); lr.DrawLine(lbs.front() - 2, 100 * fref, lbs.back() + 2, 100 * fref);
        TLine lm; lm.SetLineStyle(2); lm.SetLineColor(kBlack); lm.SetLineWidth(2);
        const auto mk = lb_marks.find(run);
        if (mk != lb_marks.end()) lm.DrawLine(mk->second.first, ylo, mk->second.first, yhi);
        TLegend lg(0.08, 0.76, 0.97, 0.94); lg.SetBorderSize(0); lg.SetFillStyle(0); lg.SetTextSize(0.042); lg.SetNColumns(2);
        lg.AddEntry(&gl_lb, Form("f_{low}: N_{trk}^{HItight}/FCal E_{T}^{A+C} < %.0f TeV^{-1} (FCal > 1 TeV, no cut)", kRatioLow), "p");
        lg.AddEntry(&lr, "f_{low}, reference runs", "l");
        TH1D* lbx1 = new TH1D(Form("lbx1_%d", run), "", 1, 0, 1); lbx1->SetDirectory(nullptr);
        lbx1->SetFillColorAlpha(bcol[1], 0.35); lbx1->SetLineWidth(0);
        TH1D* lbx2 = new TH1D(Form("lbx2_%d", run), "", 1, 0, 1); lbx2->SetDirectory(nullptr);
        lbx2->SetFillColorAlpha(bcol[2], 0.6); lbx2->SetLineWidth(0);
        lg.AddEntry(lbx1, Form("second band: f_{low} > %.0f%% (> %.0f#sigma)", 100 * kSecondBand, kZ), "f");
        if (std::any_of(blocks.begin(), blocks.end(), [](const Block& b) { return b.cat == 2; }))
            lg.AddEntry(lbx2, Form("drifting: f_{low} > %.1f#times reference (> %.0f#sigma)", kDrift, kZ), "f");
        if (mk != lb_marks.end()) lg.AddEntry(&lm, mk->second.second.c_str(), "l");
        lg.Draw();
        tl.DrawLatex(0.08, 0.955, Form("Pb+Pb 20%d, run %d, per lumi block", yr, run));
        // (2) ratio 1D, normalised
        c.cd(2); gPad->SetLeftMargin(0.08); gPad->SetRightMargin(0.03);
        TLegend lg3(0.62, 0.5, 0.97, 0.9); lg3.SetBorderSize(0); lg3.SetFillStyle(0); lg3.SetTextSize(0.04);
        const int ccol[3] = {kBlack, kOrange + 7, kRed + 2};
        std::map<int, int> nseen;
        double hmax = 0;
        for (int k = 0; k <= nb; ++k) { if (hr[k]->Integral() > 0) hr[k]->Scale(1. / hr[k]->Integral());
                                        hmax = std::max(hmax, hr[k]->GetBinContent(hr[k]->GetMaximumBin())); }
        for (int k = nb; k >= 0; --k) {
            hr[k]->SetLineColor(k == nb ? kBlack : ccol[blocks[k].cat]); hr[k]->SetLineWidth(2);
            if (k < nb) hr[k]->SetLineStyle(1 + nseen[blocks[k].cat]++ % 4);  // several blocks of one category
            hr[k]->SetTitle(";N_{trk}^{HItight} / FCal E_{T}^{A+C} [TeV^{-1}] (FCal > 1 TeV);normalised events");
            hr[k]->SetMaximum(1.15 * hmax);
            hr[k]->GetXaxis()->SetRangeUser(300, 800);
            hr[k]->Draw(k == nb ? "HIST" : "HIST SAME");
            lg3.AddEntry(hr[k], k == nb ? "other LBs" : Form("%s (%s)", lbrange(blocks[k].first, blocks[k].last).c_str(), catname(blocks[k].cat)), "l");
        }
        TLine lt; lt.SetLineStyle(2); lt.SetLineColor(kGreen + 2); lt.SetLineWidth(2);
        lt.DrawLine(kRatioLow, 0., kRatioLow, 1.15 * hmax);
        lg3.AddEntry(&lt, Form("f_{low} threshold: %.0f TeV^{-1}", kRatioLow), "l");
        lg3.Draw();
        c.SaveAs((out + Form("/lbscan_run%d_summary_pbpb_20%d.png", run, yr)).c_str());
    }
    // ---- page 2: 2D bad vs good ----
    {
        TCanvas c(Form("c_lb2_%d", run), "", 1200, 1040);
        TPad* pm = new TPad(Form("pm_%d", run), "", 0., 0.04, 1., 1.);  // bottom strip keeps the definition clear of the axes
        pm->Draw(); pm->cd();
        pm->Divide(2, 2, 0.002, 0.002);
        TLatex tl; tl.SetNDC(); tl.SetTextSize(0.036);
        std::string badlab = "second-band";
        for (auto& b : blocks) if (b.cat == 1) badlab += " " + lbrange(b.first, b.last);
        const std::string labs[2] = {"LBs without second band or drift", badlab};
        for (int g = 0; g < 2; ++g) {
            pm->cd(1 + g);
            gPad->SetLogz();  // log z as the event-selection 2D figures (populations span ~5 decades)
            gPad->SetRightMargin(0.14); gPad->SetLeftMargin(0.14); gPad->SetTopMargin(0.09);
            hnf[g]->SetContour(99); hnf[g]->Draw("COLZ");
            tl.DrawLatex(0.13, 0.93, Form("run %d, %s, no cut", run, labs[g].c_str()));
            pm->cd(3 + g); gPad->SetLogz(); gPad->SetRightMargin(0.14); gPad->SetLeftMargin(0.14); gPad->SetTopMargin(0.09);
            hfr[g]->SetContour(99); hfr[g]->Draw("COLZ");
            g4->Draw("L SAME");
            tl.DrawLatex(0.13, 0.93, Form("run %d, %s, after cuts 1-3", run, labs[g].c_str()));
            tl.SetTextAlign(31); tl.DrawLatex(0.84, 0.84, "red: cut 4 (#mu#minus5#sigma)"); tl.SetTextAlign(11);
        }
        c.cd(0);
        TLatex tdef; tdef.SetNDC(); tdef.SetTextSize(0.019);
        tdef.DrawLatex(0.01, 0.012, Form("second band: LBs with f_{low} > %.0f%%, f_{low} = fraction with N_{trk}^{HItight}/FCal E_{T}^{A+C} < %.0f TeV^{-1} at FCal > 1 TeV", 100 * kSecondBand, kRatioLow));
        c.SaveAs((out + Form("/lbscan_run%d_2d_pbpb_20%d.png", run, yr)).c_str());
    }
}

// ============================================================================ driver
void plot_pbpb_evsel_run_lb_diag(int run_year, std::vector<int> lbscan_runs = {}, std::string in_file = "",
                                 std::string out_dir = "")
{
    const int yr = run_year % 2000;
    const std::string in  = !in_file.empty() ? in_file : kBase + "pbpb_20" + std::to_string(yr) + "/evsel_run_lb_diag_pbpb_20" + std::to_string(yr) + ".root";
    const std::string out = !out_dir.empty() ? out_dir : kBase + "plots/single_b_analysis/event_selection/pbpb_20" + std::to_string(yr) + "/evsel_run_lb_diag";
    gSystem->mkdir(out.c_str(), true);
    gStyle->SetOptStat(0);
    TFile* f = TFile::Open(in.c_str(), "READ");
    if (!f || f->IsZombie()) throw std::runtime_error("cannot open " + in);
    TFile* fc = TFile::Open(PbPbEvSelCutsPath(run_year).c_str(), "READ");
    TGraph* g4 = Need<TGraph>(fc, PbPbEvSelKey::kNTrkFracCutLo);
    g4->SetLineColor(kRed); g4->SetLineWidth(2);

    // ZDC clusters: boxes read off the year's own t_A projection at |t_C| < 0.5 ns (checked on the
    // data 2026-10-08).  2026: horizontal secondary at t_A ~ +1.3 ns.  2025: no +1.3 ns cluster;
    // instead a shoulder at t_A ~ -1.3 ns and a blob at t_A ~ -3 ns.  Box edges sit on the 0.2 ns bin edges.
    std::vector<TBoxDef> boxes = {{"main", -0.6, 0.6, -0.6, 0.6, kBlack}};
    if (yr == 26) boxes.push_back({"t^{A} #in [1, 2], t^{C} #in [-0.6, 1]", 1.0, 2.0, -0.6, 1.0, kRed + 1});
    if (yr == 25) {
        boxes.push_back({"t^{A} #in [-1.8, -1.0], |t^{C}| < 0.6", -1.8, -1.0, -0.6, 0.6, kRed + 1});
        boxes.push_back({"t^{A} #in [-3.6, -2.4], |t^{C}| < 0.8", -3.6, -2.4, -0.8, 0.8, kMagenta + 1});
    }
    boxes.push_back({"satellite (+5, +5)", 4.4, 5.4, 4.4, 5.4, kBlue + 1});
    boxes.push_back({"satellite (#minus5, #minus5)", -5.4, -4.4, -5.4, -4.4, kGreen + 2});
    DrawZdcVtx(f, out, yr, boxes);
    if (yr == 26) {
        DrawCut4Hunt(f, g4, out, yr);
        std::ofstream o(out + "/lbscan_pbpb_20" + std::to_string(yr) + ".txt");
        // External LB markers drawn on the per-LB pad (from the ATLAS eLog, supplied by the user)
        const std::map<int, std::pair<int, std::string>> lb_marks = {
            {522355, {440, "LB 440: first IBL error (eLog)"}},
        };
        for (int r : lbscan_runs) DrawLBScan(f, r, lbscan_runs, g4, out, yr, lb_marks, o);
    }
    std::cout << "[INFO] outputs in " << out << std::endl;
}
