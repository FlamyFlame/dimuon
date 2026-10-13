// Abnormal-run / cut-2 study (docs/tracking/pbpb_zdc_time_cut2_window.md): per-run effect of the cut-2 window change.
// Events passing cut 1 (HLT_mu4, good run/LB; cut 1 from <cut1_file>), t^C inside its window; counts in the OLD |t^A| < 1.5 ns
// window and in the NEW PbPbZdcTimeCut(yr) window, per run -> <csv>.
#include "../../NTupleProcessingCode/PbPbEventSelConfig.h"
#include "ROOT/RDataFrame.hxx"
#include "TChain.h"
#include <fstream>
#include <map>
void cut2_window_per_run(int yr, int nparts, const std::string& cut1_file, const std::string& csv) {
    ROOT::EnableImplicitMT();
    TFile* fc = TFile::Open(cut1_file.c_str());
    TGraph* g1 = (TGraph*)fc->Get(PbPbEvSelKey::kZDCFCalCut)->Clone("g1");
    const auto w = PbPbZdcTimeCut(yr);
    const auto bad_runs = PbPbBadRuns(yr); const auto bad_lbs = PbPbBadLBRanges(yr);
    TChain ch("HeavyIonD3PD");
    for (int p = 1; p <= nparts; ++p)
        ch.Add(PbPbRawInputPath(Form("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_20%d/data_pbpb%d_part%d.root", yr, yr, p)).c_str());
    ROOT::RDataFrame df(ch);
    auto n = df.Filter([](bool h) { return h; }, {"b_HLT_mu4_L1MU3V"})
               .Filter([bad_runs, bad_lbs](UInt_t r, UInt_t lb) { return !PbPbIsBadRunLB(bad_runs, bad_lbs, (int)r, lb); }, {"RunNumber", "lbn"})
               .Filter([g1](float fP, float fN, const ROOT::RVecF& e) {
                   const double c = PbPbEvSelEvalCut(g1, (fP + fN) * 1e-6); return c <= 0. || (e[0] + e[1]) / 1000. <= c; },
                   {"FCal_Et_P", "FCal_Et_N", "zdc_ZdcEnergy"})
               .Filter([w](const ROOT::RVecF& t) { return t[0] > (float)w.c_lo && t[0] < (float)w.c_hi; }, {"zdc_ZdcTime"})
               .Define("cat", [w](const ROOT::RVecF& t) {   // 0: outside old, 1: in old only, 2: in new (new is inside old here)
                   const bool old = std::abs(t[1]) < 1.5f, nw = t[1] > (float)w.a_lo && t[1] < (float)w.a_hi;
                   return nw ? 2 : (old ? 1 : 0); }, {"zdc_ZdcTime"});
    const int rlo = yr == 25 ? 508000 : 522000, rhi = yr == 25 ? 514000 : 524000;
    auto h = n.Define("runD", [](UInt_t r) { return (double)r + 0.5; }, {"RunNumber"}).Define("catD", [](int c) { return c + 0.5; }, {"cat"})
              .Histo2D({"h", "", rhi - rlo, (double)rlo, (double)rhi, 3, 0., 3.}, "runD", "catD");
    std::map<unsigned, std::array<long,3>> m;
    for (int b = 1; b <= h->GetNbinsX(); ++b) for (int c = 0; c < 3; ++c) {
        const double v = h->GetBinContent(b, c + 1); if (v > 0) m[(unsigned)h->GetXaxis()->GetBinLowEdge(b)][c] = (long)v; }
    std::ofstream o(csv); o << "run,n_outside_old,n_old_only,n_new,new_over_old\n";
    for (auto& kv : m) { const double old = kv.second[1] + kv.second[2];
        o << kv.first << "," << kv.second[0] << "," << kv.second[1] << "," << kv.second[2] << "," << (old > 0 ? kv.second[2] / old : 0) << "\n"; }
}
