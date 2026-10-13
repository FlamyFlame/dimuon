// Is the pbpb24-vs-pbpb23 LOW-SIDE SHIFT of the HIJING overlay also present in Pb+Pb DATA?
//
// R12/R13 of docs/tracking/hijing_overlay_pbpb24_test_sample_skim.md found the 2024-conditions
// overlay (r17913) sitting ~2 % below the 2023-conditions one (r17618) in FCal and L1 sum-ET and
// -6 % in HITight tracks, and attributed it to the cross-year / cross-vertex comparison.  The
// alternative is that it is a REAL detector effect (e.g. radiation damage), in which case Pb+Pb
// 2023 vs 2024 DATA must show the same low-side shift.  This macro puts the two side by side:
// DATA on the top row, HIJING overlay on the bottom row, 2023 vs 2024 overlaid in every panel
// exactly as in plots/r17913_rtag_sanity/.
//
//   root -l -b -q 'plot_overlay_rtag_compr_pbpb_data.cxx+("fill")'       data event loop -> cache
//   root -l -b -q 'plot_overlay_rtag_compr_pbpb_data.cxx+("fcal_l1te")'  PNG 1
//   root -l -b -q 'plot_overlay_rtag_compr_pbpb_data.cxx+("ntrk")'       PNG 2
//   root -l -b -q 'plot_overlay_rtag_compr_pbpb_data.cxx+("all")'        fill if needed, then both
//
// DATA (top row) -- the established good-event definition, USER DECISION 2026-10-06:
//   HLT_mu4_L1MU3V + PbPbBadRuns(yr) + the full 5-cut Pb+Pb event selection, cuts read from that
//   year's own event_sel_cuts_pbpb_20yy.root.  This is the recipe of
//   plotting_codes/event_selection/plot_pbpb_fcal_comparison.cxx::FillFCalHist(), reproduced cut
//   for cut.  Without the selection the data tails carry out-of-time pileup and non-collision
//   background whose YEAR-DEPENDENT rate would itself mimic a calorimeter shift.
//   NO centrality cut: centrality is derived FROM FCal with the fixed 2023 calibration, so cutting
//   on it would partly define away the very shift under test.  Same for the overlay row, which has
//   no selection at all -- as in the r17913_rtag_sanity event-level panels.
//
// NTuple-processing provenance: FCal_Et_P/N, L1TE and trk_numqual are per-event skim quantities
// with no selection, weight or truth-navigation procedure behind them.  The only processing output
// that carries them (the muon-pair trees) is per-PAIR and post-selection, so it cannot give a
// per-event distribution.  The event-level selection itself is NOT re-derived here: the cut curves
// are read from the per-year cuts file produced by the event-selection pipeline.
//
// L1TE IN 2024 DATA DOES NOT EXIST.  Measured 2026-10-06 on data_pbpb24_part1/2.root at several
// offsets: L1TE and L1TE24 are identically zero, while 2023 is filled (mean 7.8 TeV).  The 2024
// skim carries no L1 sum-ET VALUE branch at all, only L1jTE* trigger-decision flags: TrigRates
// reads the legacy LVL1EnergySumRoI (SkimCode .../TrigRates.cxx:2099-2104), which the 2024 HI menu
// no longer fills after the move to Phase-I jFEX.  Per user decision the panel is kept with the
// 2023 curve alone.
#include <TFile.h>
#include <TH1D.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TParameter.h>
#include <TGraph.h>
#include <ROOT/RDataFrame.hxx>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <cstdio>
#include "../../MuonObjectsParamsAndHelpers/FullSimSampleType.h"
#include "../../NTupleProcessingCode/PbPbEventSelConfig.h"

namespace {
    // ---- the two conditions years ------------------------------------------------------------
    // r-tags of the overlay samples, as in plot_overlay_rtag_event_level.cxx.
    const char* NEW_RTAG = "r17913";   // Pb+Pb 2024 conditions
    const char* OLD_RTAG = "r17618";   // Pb+Pb 2023 conditions
    const int kC24 = kRed + 1, kC23 = kBlue + 1;

    // Resolved at first use, never at namespace scope: a std::string built by concatenation at
    // library-load time breaks cling's later RDF JIT symbol resolution (trap recorded in
    // plot_overlay_rtag_event_level.cxx).
    const std::string& MC24_LABEL() { static const std::string s = std::string("2024 conditions (") + NEW_RTAG + ")"; return s; }
    const std::string& MC23_LABEL() { static const std::string s = std::string("2023 conditions (") + OLD_RTAG + ")"; return s; }

    // pTH125_300 overlay NTUP of each conditions year (dir + file tag from FullSimSampleType.h).
    std::string ntup(int pbpb_year) {
        return FullSimSampleInputDir(FullSimSampleType::hijing, true, pbpb_year)
             + "Pythia_5p36TeV_pp_hQCD_DiMu_pTH125_300." + FullSimSampleFileTag(FullSimSampleType::hijing, pbpb_year)
             + ".NTUP.root";
    }
    const std::string& MC24_FILE() { static const std::string f = ntup(24); return f; }
    const std::string& MC23_FILE() { static const std::string f = ntup(23); return f; }

    const std::string& OUTDIR() {
        static const std::string d = FullSimPlotsDir(FullSimSampleInputDir(FullSimSampleType::hijing, true, 24))
                                   + NEW_RTAG + "_rtag_compr_PbPb_data/";
        return d;
    }
    // TEST mode: a prefix of each data chain, written to and read from a SEPARATE cache and
    // PNG name, so a quick validation run can neither overwrite nor be mistaken for the
    // full-statistics result.  One flag drives both names.
    bool& USE_TEST() { static bool b = false; return b; }
    // Cache of the DATA histograms: the loop below reads ~220 M events, so re-styling a plot must
    // never require a second pass.
    std::string CACHE() { return OUTDIR() + (USE_TEST() ? "pbpb_data_event_level_hists_TEST.root"
                                                        : "pbpb_data_event_level_hists.root"); }

    std::string outfile(const char* base) {
        return OUTDIR() + base + "_data_vs_hijing_pbpb23_vs_pbpb24" + (USE_TEST() ? "_TEST" : "") + ".png";
    }

    // ---- ONE binning per quantity, used by BOTH rows ------------------------------------------
    // Data (all centralities, mu4-triggered) and the overlay (b < 5 fm, ultra-central) populate
    // very different parts of each axis, but a shape comparison across rows is only meaningful on
    // identical bins, so the binning is defined once here.  Only the DRAWN x-range differs by row.
    struct Bins { int n; double lo, hi; };
    const Bins kFCal {320,   0.,  8.};      // TeV, 0.025 TeV bins
    const Bins kL1TE {400,   0., 40.};      // TeV, 0.1   TeV bins
    const Bins kNtrk0{320,   0., 16000.};   // tracks, 50 /bin
    const Bins kNtrk3{300,   0.,  6000.};   // tracks, 20 /bin

    void style() {
        gStyle->SetOptStat(0);
        gStyle->SetOptTitle(0);
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);
        gSystem->mkdir(OUTDIR().c_str(), true);
    }
    void axes(TH1* h, const char* xt, const char* yt) {
        h->GetXaxis()->SetTitle(xt); h->GetYaxis()->SetTitle(yt);
        h->GetXaxis()->SetTitleSize(0.052); h->GetYaxis()->SetTitleSize(0.052);
        h->GetXaxis()->SetLabelSize(0.046); h->GetYaxis()->SetLabelSize(0.046);
        h->GetXaxis()->SetTitleOffset(1.05); h->GetYaxis()->SetTitleOffset(1.30);
    }
    double median(TH1D* h) { double q, p = 0.5; h->GetQuantiles(1, &q, &p); return q; }

    // ---- data skim parts ----------------------------------------------------------------------
    // Equal to the part lists of plot_pbpb_fcal_comparison.cxx / PbPbExtras.c file_batch_max.
    // An UNDER-count silently drops data, so a missing part is fatal, not skipped.
    // Some skim parts are symlinks into /pnfs.  Reading them through the POSIX dCache mount can
    // wedge the process in uninterruptible disk sleep (seen 2026-10-08 on data_pbpb23_part2.root:
    // 0 bytes read in 20 min).  With PBPB_XROOTD=1 the /pnfs path is served through the BNL
    // xrootd door instead -- same bytes, different transport.  The door is port 1096 and the path
    // keeps its leading slash, hence the double slash.
    std::string MaybeXrootd(const std::string& path) {
        if (!gSystem->Getenv("PBPB_XROOTD")) return path;
        char* real = gSystem->ExpandPathName(path.c_str());
        TString resolved = gSystem->GetFromPipe(Form("readlink -f %s", real));
        delete[] real;
        if (resolved.BeginsWith("/pnfs/"))
            return std::string("root://dcgftp.usatlas.bnl.gov:1096/") + resolved.Data();
        return path;
    }
    std::vector<std::string> DataParts(int yr) {
        const std::string base = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/";
        if (yr == 23) return {base+"pbpb_2023/data_pbpb23_part1.root", base+"pbpb_2023/data_pbpb23_part2.root",
                              base+"pbpb_2023/data_pbpb23_part3.root", base+"pbpb_2023/data_pbpb23_part4.root",
                              base+"pbpb_2023/data_pbpb23_part5.root"};
        if (yr == 24) return {base+"pbpb_2024/data_pbpb24_part1.root", base+"pbpb_2024/data_pbpb24_part2.root"};
        throw std::runtime_error("DataParts: only Pb+Pb 2023 and 2024 are configured, got " + std::to_string(yr));
    }

    // ---- the 5 event-selection cuts of this year -----------------------------------------------
    struct CutSet {
        TGraph* g_cut1 = nullptr;   // ZDC-FCal banana
        PbPbZdcTimeWindow t{};   // cut 2, per side
        double  preamp_A = 0., preamp_C = 0.;
        TGraph* g_cut4 = nullptr;                      // nTrk HITight-fraction lower bound
        TGraph* g_cut5_lo = nullptr, *g_cut5_hi = nullptr;  // nTrk-FCal band
    };
    // No defaults anywhere: a missing graph or parameter throws.  A silent fallback would apply a
    // cut nobody chose (or, with null graphs, fail every event) and still produce a full-looking plot.
    CutSet LoadCuts(int yr) {
        const std::string path = PbPbEvSelCutsPath(yr);
        TFile* f = TFile::Open(path.c_str(), "READ");
        if (!f || f->IsZombie()) throw std::runtime_error("LoadCuts: cannot open " + path);
        auto loadG = [&](const char* k) {
            TGraph* g = (TGraph*)f->Get(k);
            if (!g) { f->Close(); throw std::runtime_error(std::string("LoadCuts: missing TGraph '") + k + "' in " + path); }
            return (TGraph*)g->Clone();
        };
        auto loadP = [&](const char* k) {
            auto* p = (TParameter<double>*)f->Get(k);
            if (!p) { f->Close(); throw std::runtime_error(std::string("LoadCuts: missing TParameter '") + k + "' in " + path); }
            return p->GetVal();
        };
        CutSet cs;
        cs.g_cut1    = loadG(PbPbEvSelKey::kZDCFCalCut);
        cs.t         = PbPbReadZdcTimeWindow(f, yr);
        cs.preamp_A  = loadP(PbPbEvSelKey::kPreampACutADC);
        cs.preamp_C  = loadP(PbPbEvSelKey::kPreampCCutADC);
        cs.g_cut4    = loadG(PbPbEvSelKey::kNTrkFracCutLo);
        cs.g_cut5_lo = loadG(PbPbEvSelKey::kNTrkFCalCutLo);
        cs.g_cut5_hi = loadG(PbPbEvSelKey::kNTrkFCalCutHi);
        f->Close();
        std::printf("Pb+Pb 20%d cuts: %.2f < tA < %.2f, %.2f < tC < %.2f ns, preamp A = %.0f, C = %.0f ADC\n", yr, cs.t.a_lo, cs.t.a_hi, cs.t.c_lo, cs.t.c_hi, cs.preamp_A, cs.preamp_C);
        return cs;
    }
}

// ============================ DATA: one pass, all four quantities ============================
// yr_sel/part_sel >= 0 process ONE skim part into its own file under cache_parts/, so the seven
// parts (five 2023 + two 2024) can run as seven concurrent processes and be hadd-ed afterwards;
// the counters are stored as a histogram precisely so that hadd sums them instead of keeping the
// first file's copy.  yr_sel < 0 keeps the original single-process behaviour.
void fill_data(Long64_t max_per_year = -1, int yr_sel = -1, int part_sel = -1) {
    style();
    std::string out = CACHE();
    if (yr_sel > 0) {
        gSystem->mkdir((OUTDIR() + "cache_parts").c_str(), true);
        out = OUTDIR() + "cache_parts/pbpb_data_hists_" + std::to_string(yr_sel) + "_p" +
              std::to_string(part_sel) + ".root";
    }
    TFile fout(out.c_str(), "RECREATE");
    for (int yr : {23, 24}) {
        if (yr_sel > 0 && yr != yr_sel) continue;
        const CutSet cs = LoadCuts(yr);
        TChain chain("HeavyIonD3PD");
        const auto parts = DataParts(yr);
        for (size_t ip = 0; ip < parts.size(); ++ip) {
            if (part_sel >= 0 && (int)ip != part_sel) continue;
            if (gSystem->AccessPathName(parts[ip].c_str()))
                throw std::runtime_error("fill_data: skim part missing, the year would be incomplete: " + parts[ip]);
            const std::string use = MaybeXrootd(parts[ip]);
            std::printf("  part %zu: %s\n", ip, use.c_str());
            chain.Add(use.c_str());
        }
        if (chain.GetNtrees() == 0)
            throw std::runtime_error("fill_data: no skim part selected for Pb+Pb 20" + std::to_string(yr));
        TH1D h_fcal(Form("h_data%d_fcal",  yr), "", kFCal.n,  kFCal.lo,  kFCal.hi);
        TH1D h_l1te(Form("h_data%d_l1te",  yr), "", kL1TE.n,  kL1TE.lo,  kL1TE.hi);
        TH1D h_nt0 (Form("h_data%d_ntrk0", yr), "", kNtrk0.n, kNtrk0.lo, kNtrk0.hi);
        TH1D h_nt3 (Form("h_data%d_ntrk3", yr), "", kNtrk3.n, kNtrk3.lo, kNtrk3.hi);
        for (TH1D* h : {&h_fcal, &h_l1te, &h_nt0, &h_nt3}) { h->SetDirectory(nullptr); h->Sumw2(); }

        chain.SetMakeClass(1);
        chain.SetBranchStatus("*", 0);
        Int_t   b_HLT = 0, run_num = 0;
        Float_t FCal_Et_P = 0.f, FCal_Et_N = 0.f, L1TE = 0.f;
        Float_t zdc_E[2] = {}, zdc_t[2] = {}, preamp[2][4] = {};
        std::vector<int>* trk_numqual = nullptr;
        for (const char* br : {"b_HLT_mu4_L1MU3V", "RunNumber", "FCal_Et_P", "FCal_Et_N", "L1TE",
                               "zdc_ZdcEnergy", "zdc_ZdcTime", "zdc_ZdcModulePreSampleAmp", "trk_numqual"})
            chain.SetBranchStatus(br, 1);
        chain.SetBranchAddress("b_HLT_mu4_L1MU3V", &b_HLT);
        chain.SetBranchAddress("RunNumber",        &run_num);
        chain.SetBranchAddress("FCal_Et_P",        &FCal_Et_P);
        chain.SetBranchAddress("FCal_Et_N",        &FCal_Et_N);
        chain.SetBranchAddress("L1TE",             &L1TE);
        chain.SetBranchAddress("zdc_ZdcEnergy",     zdc_E);
        chain.SetBranchAddress("zdc_ZdcTime",       zdc_t);
        chain.SetBranchAddress("zdc_ZdcModulePreSampleAmp", preamp);
        chain.SetBranchAddress("trk_numqual",      &trk_numqual);

        const auto bad_runs = PbPbBadRuns(yr);
        const Long64_t ntot = (max_per_year > 0) ? std::min<Long64_t>(max_per_year, chain.GetEntries()) : chain.GetEntries();
        Long64_t n_trig = 0, n_pass = 0, n_l1te_nonzero = 0;
        std::printf("Pb+Pb 20%d data: %lld events in the chain...\n", yr, ntot);
        for (Long64_t i = 0; i < ntot; ++i) {
            chain.GetEntry(i);
            if (!b_HLT) continue;
            if (bad_runs.count(run_num)) continue;
            ++n_trig;
            const float fcal_AC = (FCal_Et_P + FCal_Et_N) * 1e-6f;   // MeV -> TeV
            const float zdc_tot = (zdc_E[0] + zdc_E[1]) / 1000.f;
            const float tA = zdc_t[1], tC = zdc_t[0];                // [1] = A, [0] = C
            const float pA = preamp[1][0]+preamp[1][1]+preamp[1][2]+preamp[1][3];
            const float pC = preamp[0][0]+preamp[0][1]+preamp[0][2]+preamp[0][3];
            int ntrk_total = 0, ntrk_tight = 0;
            if (trk_numqual && (int)trk_numqual->size() >= 4) { ntrk_total = (*trk_numqual)[0]; ntrk_tight = (*trk_numqual)[3]; }
            // Cut 1: ZDC-FCal banana
            const double c1 = PbPbEvSelEvalCut(cs.g_cut1, fcal_AC);
            if (c1 > 0. && zdc_tot > c1) continue;
            // Cut 2: ZDC time box
            if (!PbPbPassZdcTime(cs.t, tA, tC)) continue;
            // Cut 3: ZDC preamp -- fail only if BOTH sides exceed the threshold
            if (pA >= (float)cs.preamp_A && pC >= (float)cs.preamp_C) continue;
            // Cuts 4 & 5 are track-based: an event with no track cannot be judged
            if (ntrk_total <= 0) continue;
            // Cut 4: HITight fraction lower bound
            if ((double)ntrk_tight / ntrk_total < PbPbEvSelEvalCut(cs.g_cut4, ntrk_total)) continue;
            // Cut 5: nTrk-FCal band
            if ((double)ntrk_tight < PbPbEvSelEvalCut(cs.g_cut5_lo, fcal_AC) ||
                (double)ntrk_tight > PbPbEvSelEvalCut(cs.g_cut5_hi, fcal_AC)) continue;
            ++n_pass;
            h_fcal.Fill(fcal_AC);
            h_l1te.Fill(L1TE / 1e3f);          // GeV -> TeV
            h_nt0 .Fill(ntrk_total);
            h_nt3 .Fill(ntrk_tight);
            if (L1TE != 0.f) ++n_l1te_nonzero;
        }
        std::printf("Pb+Pb 20%d data: %lld trigger+run, %lld pass all 5 cuts; L1TE non-zero in %lld of them\n",
                    yr, n_trig, n_pass, n_l1te_nonzero);
        fout.cd();
        for (TH1D* h : {&h_fcal, &h_l1te, &h_nt0, &h_nt3}) h->Write();
        // Counters live in a histogram, not a TParameter: hadd SUMS histograms but would keep
        // only the first part's TParameter, which would under-report every count after a merge.
        // The non-zero L1TE count is what makes "2024 has no L1TE" a measured statement rather
        // than an assumption, so it travels with the histograms.
        TH1D h_cnt(Form("h_data%d_counts", yr), "", 3, 0.5, 3.5);
        h_cnt.SetDirectory(nullptr);
        h_cnt.GetXaxis()->SetBinLabel(1, "trigger+run");
        h_cnt.GetXaxis()->SetBinLabel(2, "pass 5 cuts");
        h_cnt.GetXaxis()->SetBinLabel(3, "L1TE non-zero");
        h_cnt.SetBinContent(1, (double)n_trig);
        h_cnt.SetBinContent(2, (double)n_pass);
        h_cnt.SetBinContent(3, (double)n_l1te_nonzero);
        h_cnt.Write();
    }
    fout.Close();
    std::printf("wrote %s\n", out.c_str());
}

namespace {
    // ---- MC: the two overlay NTUPs, no selection, unweighted (one pT-hat slice each) -----------
    struct McHists { TH1D *fcal, *l1te, *nt0, *nt3; double n; };
    McHists mc_hists(const std::string& file, const char* tag) {
        ROOT::RDataFrame d("HeavyIonD3PD", file);
        auto n = d.Define("fcal_tev", "(FCal_Et_P + FCal_Et_N) / 1e6f")
                  .Define("l1te_tev", "L1TE / 1e3f")
                  .Define("ntrk0",    "(double)trk_numqual[0]")
                  .Define("ntrk3",    "(double)trk_numqual[3]");
        auto hf = n.Histo1D({Form("h_mc%s_fcal",  tag), "", kFCal.n,  kFCal.lo,  kFCal.hi},  "fcal_tev");
        auto hl = n.Histo1D({Form("h_mc%s_l1te",  tag), "", kL1TE.n,  kL1TE.lo,  kL1TE.hi},  "l1te_tev");
        auto h0 = n.Histo1D({Form("h_mc%s_ntrk0", tag), "", kNtrk0.n, kNtrk0.lo, kNtrk0.hi}, "ntrk0");
        auto h3 = n.Histo1D({Form("h_mc%s_ntrk3", tag), "", kNtrk3.n, kNtrk3.lo, kNtrk3.hi}, "ntrk3");
        McHists m;
        m.n    = *n.Count();
        m.fcal = (TH1D*)hf->Clone(); m.l1te = (TH1D*)hl->Clone();
        m.nt0  = (TH1D*)h0->Clone(); m.nt3  = (TH1D*)h3->Clone();
        for (TH1D* h : {m.fcal, m.l1te, m.nt0, m.nt3}) h->SetDirectory(nullptr);
        return m;
    }
    TH1D* get_cached(TFile& f, const char* name) {
        TH1D* h = (TH1D*)f.Get(name);
        if (!h) throw std::runtime_error(std::string("missing histogram '") + name + "' in " + CACHE() +
                                         " -- run the \"fill\" mode first");
        h = (TH1D*)h->Clone(); h->SetDirectory(nullptr);
        return h;
    }
    // One panel: 2023 vs 2024 of the same quantity, each normalised to unit area so the two
    // (very different) sample sizes overlay.  h23 may be empty-by-construction (2024 data L1TE),
    // in which case it is simply not drawn and the legend says so in three words.
    // rebin > 1 widens the common binning for a sparse row; both curves always share it.
    void panel(TH1D* h23, TH1D* h24, const char* xt, double xlo, double xhi,
               const char* header, const char* lbl23, const char* lbl24,
               const char* unit, int rebin = 1, const char* missing24 = nullptr)
    {
        gPad->SetLeftMargin(0.145); gPad->SetBottomMargin(0.135); gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.03);
        if (rebin > 1) { h23->Rebin(rebin); if (h24) h24->Rebin(rebin); }
        auto norm = [](TH1D* h) { const double s = h->Integral(0, h->GetNbinsX() + 1); if (s > 0) h->Scale(1. / s); };
        norm(h23); if (h24) norm(h24);
        const double bw = h23->GetXaxis()->GetBinWidth(1);
        axes(h23, xt, Form("fraction of events / %g %s", bw, unit));
        h23->SetLineColor(kC23); h23->SetLineWidth(2);
        if (h24) { h24->SetLineColor(kC24); h24->SetLineWidth(2); }
        h23->GetXaxis()->SetRangeUser(xlo, xhi);
        // GetMaximum() would return the STORED fMaximum once SetMaximum has been called, so the
        // peak is always taken from the bin contents.
        double ymax = 0.;
        for (TH1D* h : {h23, h24}) {
            if (!h) continue;
            for (int b = h->FindBin(xlo); b <= h->FindBin(xhi); ++b) ymax = std::max(ymax, h->GetBinContent(b));
        }
        h23->SetMaximum(1.55 * ymax); h23->SetMinimum(0.);
        h23->Draw("hist");
        if (h24) h24->Draw("hist same");
        TLegend l(0.40, 0.66, 0.955, 0.955);
        l.SetBorderSize(0); l.SetFillStyle(0); l.SetTextSize(0.038);
        l.SetHeader(header);
        const char* fmt = (h23->GetMean() > 100.) ? "median %.0f, mean %.0f" : "median %.3g, mean %.3g";
        l.AddEntry(h23, lbl23, "l");
        l.AddEntry((TObject*)nullptr, Form(fmt, median(h23), h23->GetMean()), "");
        if (h24) {
            l.AddEntry(h24, lbl24, "l");
            l.AddEntry((TObject*)nullptr, Form(fmt, median(h24), h24->GetMean()), "");
        } else if (missing24) {
            l.AddEntry((TObject*)nullptr, missing24, "");
        }
        l.DrawClone();
    }
}

// ============================ PNG 1: FCal sum-ET and L1 sum-ET ============================
void plot_fcal_l1te() {
    style();
    TFile fc(CACHE().c_str(), "READ");
    if (fc.IsZombie()) throw std::runtime_error("plot_fcal_l1te: no cache, run the \"fill\" mode first");
    TH1D *d23_f = get_cached(fc, "h_data23_fcal"), *d24_f = get_cached(fc, "h_data24_fcal");
    TH1D *d23_l = get_cached(fc, "h_data23_l1te"), *d24_l = get_cached(fc, "h_data24_l1te");
    const double n24_l1te = get_cached(fc, "h_data24_counts")->GetBinContent(3);
    fc.Close();
    const McHists m24 = mc_hists(MC24_FILE(), "24"), m23 = mc_hists(MC23_FILE(), "23");

    // Measured, not assumed: the 2024 data L1TE curve is dropped only because the loop counted
    // zero non-zero entries.  If a future re-skim fills the branch it is drawn automatically.
    const bool have_d24_l1te = (n24_l1te > 0.);
    std::printf("2024 data L1TE: %.0f non-zero events -> panel %s\n", n24_l1te, have_d24_l1te ? "overlaid" : "2023 only");

    TCanvas c("c", "", 1700, 1300);
    c.Divide(2, 2);
    c.cd(1); panel(d23_f, d24_f, "FCal #Sigma E_{T} [TeV]", 0., 6.,
                   "Pb+Pb data, HLT_mu4, full event selection",
                   "Pb+Pb 2023 data", "Pb+Pb 2024 data", "TeV", 8);
    c.cd(2); panel(d23_l, have_d24_l1te ? d24_l : nullptr, "L1 #Sigma E_{T} [TeV]", 0., 25.,
                   "Pb+Pb data, HLT_mu4, full event selection",
                   "Pb+Pb 2023 data", "Pb+Pb 2024 data", "TeV", 8,
                   "Pb+Pb 2024 data: L1TE not stored");
    c.cd(3); panel(m23.fcal, m24.fcal, "FCal #Sigma E_{T} [TeV]", 2., 6.5,
                   "HIJING overlay, b < 5 fm, #hat{p}_{T} 125-300 GeV",
                   MC23_LABEL().c_str(), MC24_LABEL().c_str(), "TeV", 8);
    c.cd(4); panel(m23.l1te, m24.l1te, "L1 #Sigma E_{T} [TeV]", 6., 26.,
                   "HIJING overlay, b < 5 fm, #hat{p}_{T} 125-300 GeV",
                   MC23_LABEL().c_str(), MC24_LABEL().c_str(), "TeV", 8);
    c.SaveAs(outfile("fcal_l1te").c_str());
}

// ============================ PNG 2: track multiplicity ============================
void plot_ntrk() {
    style();
    TFile fc(CACHE().c_str(), "READ");
    if (fc.IsZombie()) throw std::runtime_error("plot_ntrk: no cache, run the \"fill\" mode first");
    TH1D *d23_0 = get_cached(fc, "h_data23_ntrk0"), *d24_0 = get_cached(fc, "h_data24_ntrk0");
    TH1D *d23_3 = get_cached(fc, "h_data23_ntrk3"), *d24_3 = get_cached(fc, "h_data24_ntrk3");
    fc.Close();
    const McHists m24 = mc_hists(MC24_FILE(), "24"), m23 = mc_hists(MC23_FILE(), "23");

    // trk_numqual[0] = every InDet track with p_T > 400 MeV (no quality requirement);
    // trk_numqual[3] = the same with HITight.  [4-7] repeat the four without the p_T cut and are
    // bit-equal to [0-3] here, because the HI track reconstruction stores nothing below 400 MeV.
    TCanvas c("c", "", 1700, 1300);
    c.Divide(2, 2);
    c.cd(1); panel(d23_0, d24_0, "N_{trk} (p_{T} > 400 MeV)", 0., 6000.,
                   "Pb+Pb data, HLT_mu4, full event selection",
                   "Pb+Pb 2023 data", "Pb+Pb 2024 data", "tracks", 5);
    c.cd(2); panel(d23_3, d24_3, "N_{trk} (p_{T} > 400 MeV, HITight)", 0., 3500.,
                   "Pb+Pb data, HLT_mu4, full event selection",
                   "Pb+Pb 2023 data", "Pb+Pb 2024 data", "tracks", 5);
    c.cd(3); panel(m23.nt0, m24.nt0, "N_{trk} (p_{T} > 400 MeV)", 2000., 12000.,
                   "HIJING overlay, b < 5 fm, #hat{p}_{T} 125-300 GeV",
                   MC23_LABEL().c_str(), MC24_LABEL().c_str(), "tracks", 5);
    c.cd(4); panel(m23.nt3, m24.nt3, "N_{trk} (p_{T} > 400 MeV, HITight)", 1000., 4000.,
                   "HIJING overlay, b < 5 fm, #hat{p}_{T} 125-300 GeV",
                   MC23_LABEL().c_str(), MC24_LABEL().c_str(), "tracks", 5);
    c.SaveAs(outfile("trk_multiplicity").c_str());
}

void plot_overlay_rtag_compr_pbpb_data(const char* what = "all") {
    const std::string w = what;
    USE_TEST() = (w.size() > 5 && w.compare(w.size() - 5, 5, "_test") == 0);
    const std::string base = USE_TEST() ? w.substr(0, w.size() - 5) : w;
    if (base.rfind("fill:", 0) == 0) {           // fill:<year>:<part index>, one part per process
        const size_t c = base.find(':', 5);
        fill_data(-1, std::stoi(base.substr(5, c - 5)), std::stoi(base.substr(c + 1)));
    }
    else if (base == "fill")      fill_data(USE_TEST() ? 2000000 : -1);
    else if (base == "fcal_l1te") plot_fcal_l1te();
    else if (base == "ntrk")      plot_ntrk();
    else if (base == "all") {
        if (gSystem->AccessPathName(CACHE().c_str())) fill_data(USE_TEST() ? 2000000 : -1);
        plot_fcal_l1te();
        plot_ntrk();
    }
    else std::printf("unknown mode '%s' -- fill|fcal_l1te|ntrk|all, each with an optional "
                     "'_test' suffix (2 M events/year, separate cache and PNG names)\n", what);
}
