// ZDC t^C vs t^A, one panel per run (docs/tracking/pbpb_zdc_time_cut2_window.md): are the off-centre t^A clusters a
// few-run feature (-> drop runs, keep one common cut-2 window for all years) or present in every run?
// Population = the cut-2 standalone figure of plot_pbpb_event_sel_cuts.cxx: b_HLT_mu4_L1MU3V, good run/LB
// (PbPbBadRuns/PbPbBadLBRanges), passing cut 1 (g_ZDC_FCal_cut of the year's cuts file; pass if cut <= 0 or E_ZDC <= cut).
// Each panel normalised to unit area, common log z range. Red dashed box = +-1.5 ns (common window); blue line = the year's
// PbPbZdcTimeCut() edges where they differ from +-1.5.
// Also <outdir>/zdc_time_regions_per_run_pbpb_20YY.csv: per run, with t^C inside +-1.5 ns, the fractions in
//   core (-0.6 < t^A < 1.5), shoulder (-1.5 < t^A < -0.6), blob (-4 < t^A < -2), and the counts.
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
#include "ROOT/RDataFrame.hxx"
#include "TCanvas.h"
#include "TChain.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include "../../NTupleProcessingCode/PbPbEventSelConfig.h"

void plot_pbpb_zdc_time_per_run(int yr, int nparts, const std::string& outdir) {
    yr %= 2000;
    ROOT::EnableImplicitMT();
    gStyle->SetOptStat(0); gStyle->SetPalette(kBird);
    TFile* fc = TFile::Open(PbPbEvSelCutsPath(yr).c_str());
    TGraph* g1 = (TGraph*)fc->Get(PbPbEvSelKey::kZDCFCalCut)->Clone("g1");
    const auto w = PbPbZdcTimeCut(yr);
    const auto bad_runs = PbPbBadRuns(yr); const auto bad_lbs = PbPbBadLBRanges(yr);
    TChain ch("HeavyIonD3PD");
    for (int p = 1; p <= nparts; ++p)
        ch.Add(PbPbRawInputPath(Form("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_20%d/data_pbpb%d_part%d.root", yr, yr, p)).c_str());
    ROOT::RDataFrame df(ch);
    auto base = df.Filter([](bool h) { return h; }, {"b_HLT_mu4_L1MU3V"})
                  .Filter([bad_runs, bad_lbs](UInt_t r, UInt_t lb) { return !PbPbIsBadRunLB(bad_runs, bad_lbs, (int)r, lb); }, {"RunNumber", "lbn"})
                  .Filter([g1](float fP, float fN, const ROOT::RVecF& e) {
                      const double c = PbPbEvSelEvalCut(g1, (fP + fN) * 1e-6); return c <= 0. || (e[0] + e[1]) / 1000. <= c; },
                      {"FCal_Et_P", "FCal_Et_N", "zdc_ZdcEnergy"})
                  .Define("tA", [](const ROOT::RVecF& t) { return t[1]; }, {"zdc_ZdcTime"})   // [1]=A
                  .Define("tC", [](const ROOT::RVecF& t) { return t[0]; }, {"zdc_ZdcTime"});  // [0]=C
    // run list = the year's HLT_mu4 lumitable (same run list as the data; checked in data_per_run_signal_yield_qa.md),
    // minus PbPbBadRuns; a run with no events after the selection is skipped
    std::vector<int> rl;
    {
        std::ifstream lf(Form("/gpfs/mnt/atlasgpfs01/usatlas/workarea/yuhanguo/dimuon_codes/IntNotes/data/luminosity/pbpb_20%d/lumitable_pbpb_%d_HLT_mu4.csv", yr, yr));
        std::string line; std::getline(lf, line);
        while (std::getline(lf, line)) {
            const std::string f0 = line.substr(0, line.find(','));
            if (f0.empty() || f0.find_first_not_of("0123456789 ") != std::string::npos) continue;
            const int r = std::stoi(f0);
            if (!bad_runs.count(r)) rl.push_back(r);
        }
        if (rl.empty()) throw std::runtime_error("empty run list for year " + std::to_string(yr));
    }
    std::map<int, ROOT::RDF::RResultPtr<TH2D>> h;
    for (int r : rl)
        h[r] = base.Filter([r](UInt_t x) { return (int)x == r; }, {"RunNumber"})
                   .Histo2D({Form("h%d", r), ";ZDC t^{A} [ns];ZDC t^{C} [ns]", 100, -10., 10., 100, -10., 10.}, "tA", "tC");

    gSystem->mkdir(outdir.c_str(), true);
    std::ofstream o(outdir + Form("/zdc_time_regions_per_run_pbpb_20%d.csv", yr));
    o << "run,n_after_cut1,n_tC_in,frac_core,frac_shoulder,frac_blob\n";
    std::vector<int> rk;
    for (int r : rl) if (h[r]->Integral() > 0) rk.push_back(r); else std::cout << "run " << r << ": no events, skipped\n";
    rl = rk;
    const int N = rl.size(), ncol = (int)std::floor(std::sqrt((double)N)), nrow = (N + ncol - 1) / ncol;   // nrows >= ncols
    TCanvas c("c", "", 360 * ncol, 300 * nrow);
    c.Divide(ncol, nrow, 0.001, 0.001);
    std::vector<TH2D*> keep;
    for (int i = 0; i < N; ++i) {
        const int r = rl[i];
        TH2D* hh = (TH2D*)h[r]->Clone(); keep.push_back(hh);
        const double ntot = hh->Integral(0, 101, 0, 101);
        auto ax = hh->GetXaxis(); auto ay = hh->GetYaxis();
        const int cy0 = ay->FindBin(-1.5 + 1e-4), cy1 = ay->FindBin(1.5 - 1e-4);
        auto reg = [&](double a0, double a1) { return hh->Integral(ax->FindBin(a0 + 1e-4), ax->FindBin(a1 - 1e-4), cy0, cy1); };
        const double nC = hh->Integral(0, 101, cy0, cy1);
        const double core = reg(-0.6, 1.5), sh = reg(-1.5, -0.6), blob = reg(-4., -2.);
        o << r << "," << (long)ntot << "," << (long)nC << "," << core / nC << "," << sh / nC << "," << blob / nC << "\n";
        hh->Scale(1. / ntot);
        hh->SetMinimum(1e-7); hh->SetMaximum(0.2);
        c.cd(i + 1); gPad->SetLogz(); gPad->SetRightMargin(0.03); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13); gPad->SetTopMargin(0.03);
        hh->GetXaxis()->SetTitleSize(0.06); hh->GetYaxis()->SetTitleSize(0.06); hh->GetXaxis()->SetLabelSize(0.05); hh->GetYaxis()->SetLabelSize(0.05);
        hh->GetYaxis()->SetTitleOffset(0.95);
        hh->SetContour(99); hh->Draw("COL");
        TLine l; l.SetLineColor(kRed); l.SetLineWidth(1); l.SetLineStyle(2);
        l.DrawLine(-1.5, -1.5, -1.5, 1.5); l.DrawLine(1.5, -1.5, 1.5, 1.5); l.DrawLine(-1.5, -1.5, 1.5, -1.5); l.DrawLine(-1.5, 1.5, 1.5, 1.5);
        TLine lb; lb.SetLineColor(kBlue + 1); lb.SetLineWidth(1);
        if (std::abs(w.a_lo + 1.5) > 1e-9) lb.DrawLine(w.a_lo, -1.5, w.a_lo, 1.5);
        if (std::abs(w.a_hi - 1.5) > 1e-9) lb.DrawLine(w.a_hi, -1.5, w.a_hi, 1.5);
        TLatex t; t.SetNDC(); t.SetTextSize(0.065);
        t.DrawLatex(0.17, 0.89, Form("run %d", r));
        t.SetTextSize(0.055);
        t.DrawLatex(0.17, 0.81, Form("shoulder %.1f%%, blob %.2f%%", 100 * sh / nC, 100 * blob / nC));
    }
    c.SaveAs((outdir + Form("/zdc_time_A_vs_C_per_run_pbpb_20%d.png", yr)).c_str());
}
