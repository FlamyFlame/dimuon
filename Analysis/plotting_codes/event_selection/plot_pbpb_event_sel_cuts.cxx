// plot_pbpb_event_sel_cuts.cxx
// Sequential 5-cut event selection pipeline for PbPb data.
//   Cut 1: ZDC vs FCal banana (loaded from event_sel_cuts_pbpb_YYYY.root)
//   Cut 2: ZDC time |tA| < 1.5 ns AND |tC| < 1.5 ns
//   Cut 3: ZDC preamp sum: fail only if BOTH A and C exceed threshold (hard cut per year; per-run mu+7sigma for 2025)
//   Cut 4: nTrk HItight fraction vs total nTrk — lower bound at mu-5sigma (per-slice Gauss)
//   Cut 5: nTrk HItight vs FCal ET — band [mu-5sigma, mu+5sigma] (per-slice Gauss)
// All cuts saved to event_sel_cuts_pbpb_YYYY.root.
//
// Output per cut: standalone canvas, 5-panel pass canvas, 5-panel fail canvas.
// Additionally: 5-panel no-cuts canvas, per-cut survival rate plot (centrality 0-79).
//
// Usage:
//   .L plot_pbpb_event_sel_cuts.cxx+
//   plot_pbpb_event_sel_cuts(24)

#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLatex.h"
#include "TLine.h"
#include "TBox.h"
#include "TParameter.h"
#include "TMath.h"
#include "../../NTupleProcessingCode/PbPbEventSelConfig.h"

// ---- Cut tuning -------------------------------------------------------------
static const double CUT2_T_NS_ALT    = 1.5;    // ZDC time window [ns]
// Per-year ZDC preamp upper cut [ADC] — evaluate from Gaussian-to-tail turnover in 1D preamp histogram.
// Values are {side_A, side_C}. pbpb2023 not yet evaluated (uses 2024 value as placeholder).
static std::pair<float,float> GetPreampCuts(int yr) {
    static const std::map<int, std::pair<float,float>> kMap = {
        {23, {300.f, 300.f}},  // side A ~300 ADC, side C ~300 ADC
        {24, {420.f, 420.f}},  // symmetric; both sides well above 420
        {25, {850.f, 700.f}},  // side A ~850 ADC, side C ~700 ADC
    };
    auto it = kMap.find(yr);
    return (it != kMap.end()) ? it->second : std::make_pair(385.f, 385.f);
}
static const double CUT4_N_SIGMA_ALT = 5.0;    // nTrk frac lower cut
static const double CUT5_N_SIGMA_ALT = 5.0;    // nTrk-FCal band half-width

static const int    N_BINS_FRAC_SL_ALT  = 4;
static const int    N_BINS_FCAL_SL_ALT  = 2;
static const int    MIN_ENTRIES_FIT_ALT = 100;
static const double FCAL_CUT5_LIN_MIN_ALT = 0.2;  // FCal lower bound for cut-5 linear fit [TeV]

// ---- Centrality recalculation from FCal ET ----------------------------------
// PbPb2023 Glauber v3.2 table: FCal ET lower bounds (TeV) per integer percentile bin (0-84).
// Mirrors PairPbPbExtras::FCal_ET_Bins_PbPb2023 in MuonPairPbPb.h.
// Used for years where the 'centrality' branch is unfilled (e.g. pbpb2025).
static const float kFCalBinsPbPb2023[85] = {
  4.51272f, 4.32043f, 4.15372f, 3.99602f, 3.84498f, 3.69944f, 3.55802f, 3.42045f, 3.28744f, 3.15972f, // 0-10
  3.03748f, 2.92012f, 2.80723f, 2.69878f, 2.59464f, 2.49406f, 2.39646f, 2.30180f, 2.21028f, 2.12188f, // 10-20
  2.03659f, 1.95428f, 1.87489f, 1.79842f, 1.72484f, 1.65387f, 1.58516f, 1.51853f, 1.45406f, 1.39178f, // 20-30
  1.33168f, 1.27363f, 1.21752f, 1.16336f, 1.11112f, 1.06069f, 1.01210f, 0.96518f, 0.91991f, 0.87632f, // 30-40
  0.83431f, 0.79384f, 0.75493f, 0.71753f, 0.68162f, 0.64714f, 0.61393f, 0.58195f, 0.55126f, 0.52171f, // 40-50
  0.49348f, 0.46640f, 0.44043f, 0.41569f, 0.39200f, 0.36933f, 0.34768f, 0.32701f, 0.30731f, 0.28853f, // 50-60
  0.27064f, 0.25357f, 0.23726f, 0.22178f, 0.20720f, 0.19325f, 0.17983f, 0.16748f, 0.15557f, 0.14435f, // 60-70
  0.13388f, 0.12389f, 0.11468f, 0.10598f, 0.09765f, 0.09015f, 0.08264f, 0.07589f, 0.06955f, 0.06321f, // 70-80
  0.05760f, 0.05273f, 0.04787f, 0.04300f, 0.03885f                                                     // 80-85
};
// Returns integer centrality percentile [0,84], or -1 if FCal ET below the 85th-percentile boundary.
// Identical logic to PairPbPbExtras::GetCentrality() with std::greater binary search.
static int CentralityFromFCal2023(float fcal_AC) {
    const int n = 85;
    // lower_bound with greater<float>: first position where bins[k] <= fcal_AC
    int lo = 0, hi = n;
    while (lo < hi) {
        int mid = (lo + hi) / 2;
        if (kFCalBinsPbPb2023[mid] > fcal_AC) lo = mid + 1;
        else                                   hi = mid;
    }
    return (lo >= n) ? -1 : lo;
}

// ---- Input file map ---------------------------------------------------------
static std::map<int, std::vector<std::string>> BuildFileMapAlt() {
    const std::string base = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/";
    return {
        // pbpb2023 should have part1..part4, but currently only part1..part3 finished skimming.
        // This needs to be rerun after part4 is available.
        {23, {
            base + "pbpb_2023/data_pbpb23_part1.root",
            base + "pbpb_2023/data_pbpb23_part2.root",
            base + "pbpb_2023/data_pbpb23_part3.root",
        }},
        {24, {
            base + "pbpb_2024/data_pbpb24_part1.root",
            base + "pbpb_2024/data_pbpb24_part2.root",
        }},
        {25, {
            base + "pbpb_2025/data_pbpb25_part1.root",
            base + "pbpb_2025/data_pbpb25_part2.root",
            base + "pbpb_2025/data_pbpb25_part3.root",
            base + "pbpb_2025/data_pbpb25_part4.root",
            base + "pbpb_2025/data_pbpb25_part5.root",
            base + "pbpb_2025/data_pbpb25_part6.root",
        }},
    };
}

// ---- Histogram stage/type indices ------------------------------------------
enum StageAlt {
    kNoCut_A  = 0,
    kC1Pass_A = 1, kC1Fail_A = 2,
    kC2Pass_A = 3, kC2Fail_A = 4,
    kC3Pass_A = 5, kC3Fail_A = 6,
    kC3FailOne_A = 11, kC3FailBoth_A = 12,  // fail exactly one side / both sides
    kC4Pass_A = 7, kC4Fail_A = 8,
    kC5Pass_A = 9, kC5Fail_A = 10,
    kNStages_A = 13
};
enum HTypeAlt {
    kZdcFcal_A  = 0,
    kZdcTime_A  = 1,
    kFCalAC_A   = 2,
    kFrac_A     = 3,
    kNtrkFcal_A = 4,
    kPreampAC_A = 5,
    kNHTypes_A  = 6
};
static const int kNPanelTypes_A = 5;

// Cut labels provided by kPbPbEvSelCutLabel() from PbPbEventSelConfig.h.

// ---- Per-event data cache ---------------------------------------------------
struct EvDataAlt {
    float fcal_A, fcal_C, fcal_AC, zdc_tot, tA, tC, preamp_A, preamp_C;
    int   ntrk_total, ntrk_tight;
    int   centrality;
};

// ---- Helpers ----------------------------------------------------------------
// EvalCutAlt: clamped TGraph evaluator — use shared PbPbEvSelEvalCut from PbPbEventSelConfig.h

struct SliceFitAlt { double x, mu, sig, ex; };

static std::vector<SliceFitAlt> FitSlices2DAlt(TH2D* h2, int n_bins_per_slice) {
    std::vector<SliceFitAlt> result;
    const int nx = h2->GetNbinsX();
    const int n_slices = nx / n_bins_per_slice;
    TF1 gf("__gf_slice_alt", "gaus");
    for (int isl = 0; isl < n_slices; ++isl) {
        const int ix_lo = isl * n_bins_per_slice + 1;
        const int ix_hi = std::min((isl + 1) * n_bins_per_slice, nx);
        TH1D* hp = (TH1D*)h2->ProjectionY(Form("__hp_alt_%s_%d", h2->GetName(), isl), ix_lo, ix_hi);
        hp->SetDirectory(nullptr);
        if (hp->GetEntries() < MIN_ENTRIES_FIT_ALT) { delete hp; continue; }
        const double mu0  = hp->GetBinCenter(hp->GetMaximumBin());
        const double rms  = (hp->GetRMS() > 0.) ? hp->GetRMS() : 1.0;
        const double ymin = hp->GetXaxis()->GetXmin();
        const double ymax = hp->GetXaxis()->GetXmax();
        gf.SetRange(std::max(ymin, mu0 - 3.*rms), std::min(ymax, mu0 + 3.*rms));
        gf.SetParameters(hp->GetMaximum(), mu0, rms);
        if (hp->Fit(&gf, "RNQ") != 0 || gf.GetParameter(2) <= 0.) { delete hp; continue; }
        double mu1 = gf.GetParameter(1), sig1 = gf.GetParameter(2);
        gf.SetRange(std::max(ymin, mu1 - 3.*sig1), std::min(ymax, mu1 + 3.*sig1));
        gf.SetParameters(gf.GetParameter(0), mu1, sig1);
        if (hp->Fit(&gf, "RNQ") != 0 || gf.GetParameter(2) <= 0.) { delete hp; continue; }
        const double xlo = h2->GetXaxis()->GetBinLowEdge(ix_lo);
        const double xhi = h2->GetXaxis()->GetBinUpEdge(ix_hi);
        result.push_back({0.5*(xlo+xhi), gf.GetParameter(1), gf.GetParameter(2), 0.5*(xhi-xlo)});
        delete hp;
    }
    return result;
}

// =============================================================================
class PbPbEventSelCutsAlt {
public:
    explicit PbPbEventSelCutsAlt(int run_year, bool is_alt = false)
        : run_year_(run_year % 2000), is_alt_(is_alt) {
        yr_       = std::to_string(run_year_);
        out_dir_  = is_alt_
            ? "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/"
              "single_b_analysis/event_selection_alternative_banana/pbpb_20" + yr_
            : "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/"
              "single_b_analysis/event_selection/pbpb_20" + yr_;
        cuts_dir_ = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_20" + yr_;
        auto fm = BuildFileMapAlt();
        if (fm.find(run_year_) == fm.end())
            throw std::runtime_error("No files for year " + yr_);
        infiles_ = fm.at(run_year_);
        auto [pa, pc] = GetPreampCuts(run_year_);
        cut3_preamp_A_ = pa;
        cut3_preamp_C_ = pc;
    }

    void Run() {
        gStyle->SetOptStat(0);
        gStyle->SetPalette(kBird);
        gSystem->mkdir(out_dir_.c_str(), true);
        bad_runs_ = PbPbBadRuns(run_year_);
        if (!bad_runs_.empty())
            std::cout << "Excluding " << bad_runs_.size()
                      << " bad run(s) from PbPb 20" << yr_ << std::endl;
        BookHists();
        LoadCut1();
        DerivePerRunPreampCuts25();  // no-op for 23/24; populates per_run_preamp_cuts_ for 25
        FillHists();
        DeriveCut4();
        FillCuts45();
        SaveCuts();
        SavePlots();
    }

private:
    int  run_year_;
    bool is_alt_ = false;
    std::string yr_, out_dir_, cuts_dir_;
    std::vector<std::string> infiles_;
    float cut3_preamp_A_ = 385.f;
    float cut3_preamp_C_ = 385.f;
    std::unordered_set<int> bad_runs_;

    TH2D*   hh_[kNStages_A][kNHTypes_A] = {};
    TH2D*   h_zdctime_ctr_[6] = {};   // ZDC time AC corr: [0-4] = top 1-5%, [5] = 5-10%
    TH1D*   h_centr_before_cuts_ = nullptr; // centrality distribution before any cuts (post-HLT)
    TH1D*   h_centr_after_cuts_  = nullptr; // centrality distribution after all 5 cuts
    // FCal ET lower bounds (TeV) per centrality bin — from PbPb2023 Glauber table (MuonPairPbPb.h).
    // Applied to all Run3 years until per-year tables are derived.
    //   [0]=0-1%, [1]=1-2%, [2]=2-3%, [3]=3-4%, [4]=4-5%, [5]=5-10% (lower bound = 9th-percentile entry)
    static constexpr float kCtrThresh[6] = {4.51272f, 4.32043f, 4.15372f, 3.99602f, 3.84498f, 3.15972f};
    TGraph* g_cut1_    = nullptr;
    TGraph* g_cut4_    = nullptr;
    TGraph* g_cut5_lo_ = nullptr;
    TGraph* g_cut5_hi_ = nullptr;

    Long64_t n_HLT_             = 0;
    Long64_t n_ctr80_[kNStages_A] = {};
    std::vector<EvDataAlt> stored_;  // events passing cuts 1-3

    // Per-run preamp cuts for PbPb25 (derived by DerivePerRunPreampCuts25)
    std::map<int, std::pair<float,float>> per_run_preamp_cuts_;

    // Part III: preamp AC correlation split by Cut-3 pass/fail status
    // Groups: [0]=both pass, [1]=exactly one side fails, [2]=both fail
    TH2D*     h_preamp_grp_[3] = {};
    long long n_after_c12_   = 0;
    long long n_preamp_grp_[3] = {};

    // -------------------------------------------------------------------------
    void BookHists() {
        struct Spec { int nx; double xlo, xhi; int ny; double ylo, yhi;
                      const char* xtitle; const char* ytitle; };
        Spec specs[kNHTypes_A] = {
            {120, -0.5, 5.5,    100,     0.,  400., "FCal E_{T}^{A+C} [TeV]",          "ZDC E_{total} [TeV]"},
            {100, -10., 10.,    100,   -10.,   10., "ZDC t^{A} [ns]",                  "ZDC t^{C} [ns]"},
            {110, -0.5, 3.0,    110,   -0.5,   3.0, "FCal E_{T}^{A} [TeV]",            "FCal E_{T}^{C} [TeV]"},
            {100,   0., 5000.,  110,     0.,   1.1, "N_{trk}^{total} (p_{T}>400 MeV)", "N_{trk}^{HItight} / N_{trk}^{total}"},
            {120, -0.5, 5.5,    100,     0., 3500., "FCal E_{T}^{A+C} [TeV]",          "N_{trk}^{HItight}"},
            {200,-1000., 3000., 200, -1000., 3000., "ZDC preamp sum side A [ADC]",      "ZDC preamp sum side C [ADC]"},
        };
        for (int s = 0; s < kNStages_A; ++s)
            for (int t = 0; t < kNHTypes_A; ++t) {
                const Spec& sp = specs[t];
                hh_[s][t] = new TH2D(Form("hh_alt_%d_%d", s, t),
                    Form(";%s;%s", sp.xtitle, sp.ytitle),
                    sp.nx, sp.xlo, sp.xhi, sp.ny, sp.ylo, sp.yhi);
                hh_[s][t]->SetDirectory(nullptr);
            }
        // Per-centrality ZDC time AC correlation (no cuts, post-HLT)
        // [0-4] = 0-1%,1-2%,...,4-5% (centrality==k); [5] = 5-10% (centrality 5-9)
        for (int k = 0; k < 6; ++k) {
            h_zdctime_ctr_[k] = new TH2D(Form("h_zdctime_ctr_%d", k),
                Form(";ZDC t^{A} [ns];ZDC t^{C} [ns]"),
                100, -10., 10., 100, -10., 10.);
            h_zdctime_ctr_[k]->SetDirectory(nullptr);
        }
        h_centr_before_cuts_ = new TH1D("h_centr_before_cuts",
            ";Centrality percentile;Events (before event selection)",
            80, -0.5, 79.5);
        h_centr_before_cuts_->SetDirectory(nullptr);
        h_centr_after_cuts_ = new TH1D("h_centr_after_cuts",
            ";Centrality percentile;Events (after all 5 cuts)",
            80, -0.5, 79.5);
        h_centr_after_cuts_->SetDirectory(nullptr);

        // Part III: preamp A vs C correlation for three Cut-3 outcome groups
        static const char* grp_tag[3] = {"both_pass", "one_fail", "both_fail"};
        for (int g = 0; g < 3; ++g) {
            h_preamp_grp_[g] = new TH2D(Form("h_pgrp_%s_%d", grp_tag[g], run_year_),
                ";ZDC preamp sum side A [ADC];ZDC preamp sum side C [ADC]",
                200, -1000., 3000., 200, -1000., 3000.);
            h_preamp_grp_[g]->SetDirectory(nullptr);
        }
    }

    // -------------------------------------------------------------------------
    void LoadCut1() {
        const std::string suffix = is_alt_ ? "_alt" : "";
        const std::string path = cuts_dir_ + "/event_sel_cuts_pbpb_20" + yr_ + suffix + ".root";
        TFile* f = TFile::Open(path.c_str(), "READ");
        if (!f || f->IsZombie())
            throw std::runtime_error("Cannot open cuts file: " + path);
        TGraph* g = (TGraph*)f->Get(PbPbEvSelKey::kZDCFCalCut);
        if (!g) { f->Close(); throw std::runtime_error(std::string(PbPbEvSelKey::kZDCFCalCut) + " not found in " + path); }
        g_cut1_ = (TGraph*)g->Clone("g_cut1_local_alt");
        f->Close();
        std::cout << "Loaded cut 1 TGraph: " << g_cut1_->GetN() << " points." << std::endl;
    }

    // -------------------------------------------------------------------------
    // PbPb25 only: derive per-run mu+7sigma ZDC preamp cuts from a first pass over
    // the raw skim (trigger + Cut1 + Cut2), populating per_run_preamp_cuts_.
    // No-op for 23/24.  Must be called after LoadCut1().
    void DerivePerRunPreampCuts25() {
        if (run_year_ != 25) return;

        static const double kFitLo  = -800., kFitHi = 1500.;
        static const int    kNSig   = 7;      // mu + 7*sigma
        static const int    kMinEnt = 50;

        TChain chain("HeavyIonD3PD", "HeavyIonD3PD");
        for (const auto& f : infiles_) {
            if (!gSystem->AccessPathName(f.c_str())) chain.Add(f.c_str());
        }
        chain.SetMakeClass(1);
        chain.SetBranchStatus("*", 0);

        Int_t   b_HLT = 0, run_num = 0;
        Float_t FCal_Et_P = 0.f, FCal_Et_N = 0.f;
        Float_t zdc_E[2] = {}, zdc_t[2] = {};
        Float_t preamp[2][4] = {};

        chain.SetBranchStatus("b_HLT_mu4_L1MU3V",         1);
        chain.SetBranchStatus("RunNumber",                  1);
        chain.SetBranchStatus("FCal_Et_P",                  1);
        chain.SetBranchStatus("FCal_Et_N",                  1);
        chain.SetBranchStatus("zdc_ZdcEnergy",              1);
        chain.SetBranchStatus("zdc_ZdcTime",                1);
        chain.SetBranchStatus("zdc_ZdcModulePreSampleAmp",  1);

        chain.SetBranchAddress("b_HLT_mu4_L1MU3V",         &b_HLT);
        chain.SetBranchAddress("RunNumber",                  &run_num);
        chain.SetBranchAddress("FCal_Et_P",                  &FCal_Et_P);
        chain.SetBranchAddress("FCal_Et_N",                  &FCal_Et_N);
        chain.SetBranchAddress("zdc_ZdcEnergy",              zdc_E);
        chain.SetBranchAddress("zdc_ZdcTime",                zdc_t);
        chain.SetBranchAddress("zdc_ZdcModulePreSampleAmp",  preamp);

        std::map<int, TH1D*> hA_map, hC_map;
        const Long64_t n_tot = chain.GetEntries();
        std::cout << "DerivePerRunPreampCuts25: " << n_tot << " events, scanning..." << std::flush;
        long long n_sel = 0;

        for (Long64_t i = 0; i < n_tot; ++i) {
            chain.GetEntry(i);
            if (!b_HLT) continue;
            if (bad_runs_.count(run_num)) continue;
            const float fcal_AC = (FCal_Et_P + FCal_Et_N) * 1e-6f;
            const float zdc_tot = (zdc_E[0] + zdc_E[1]) / 1000.f;
            if (zdc_tot > (float)PbPbEvSelEvalCut(g_cut1_, fcal_AC)) continue;  // Cut 1
            if (std::abs(zdc_t[1]) >= CUT2_T_NS_ALT) continue;                 // Cut 2 A
            if (std::abs(zdc_t[0]) >= CUT2_T_NS_ALT) continue;                 // Cut 2 C
            ++n_sel;

            if (hA_map.find(run_num) == hA_map.end()) {
                hA_map[run_num] = new TH1D(Form("prc25_hA_%d", run_num), "", 198, -1000., 2960.);
                hC_map[run_num] = new TH1D(Form("prc25_hC_%d", run_num), "", 198, -1000., 2960.);
                hA_map[run_num]->SetDirectory(nullptr);
                hC_map[run_num]->SetDirectory(nullptr);
            }
            float pA = 0.f, pC = 0.f;
            for (int k = 0; k < 4; ++k) pA += preamp[1][k];
            for (int k = 0; k < 4; ++k) pC += preamp[0][k];
            hA_map[run_num]->Fill(pA);
            hC_map[run_num]->Fill(pC);
        }
        std::cout << "  " << n_sel << " pass Cut1+Cut2; " << hA_map.size() << " runs." << std::endl;

        // Two-pass Gaussian fit per run per side; returns {mu+kNSig*sigma, ok}
        auto fitCut = [&](TH1D* h, const char* name) -> std::pair<double, bool> {
            if (h->GetEntries() < kMinEnt) return {0., false};
            const double seed_hi = kFitLo + 0.6 * (kFitHi - kFitLo);
            TF1 g1(Form("%s_p1", name), "gaus", kFitLo, seed_hi);
            g1.SetParameters(h->GetMaximum(), (kFitLo + seed_hi) * 0.5,
                             (seed_hi - kFitLo) * 0.25);
            h->Fit(&g1, "RQN");
            const double mu1 = g1.GetParameter(1), sig1 = std::abs(g1.GetParameter(2));
            if (sig1 < 1.) return {0., false};
            TF1 g2(Form("%s_p2", name), "gaus", mu1 - 3.*sig1, mu1 + 1.5*sig1);
            g2.SetParameters(g1.GetParameter(0), mu1, sig1);
            h->Fit(&g2, "RQN");
            const double mu2 = g2.GetParameter(1), sig2 = std::abs(g2.GetParameter(2));
            if (sig2 < 1.) return {0., false};
            return {mu2 + kNSig * sig2, true};
        };

        const auto [hard_A, hard_C] = GetPreampCuts(25);
        int n_fit_ok = 0;
        for (auto& [run, hA] : hA_map) {
            TH1D* hC = hC_map[run];
            auto [cutA, okA] = fitCut(hA, Form("prc25_fA_%d", run));
            auto [cutC, okC] = fitCut(hC, Form("prc25_fC_%d", run));
            per_run_preamp_cuts_[run] = {
                okA ? (float)cutA : hard_A,
                okC ? (float)cutC : hard_C
            };
            if (okA && okC) ++n_fit_ok;
            std::cout << Form("  Run %d  cutA=%.0f(%s)  cutC=%.0f(%s)\n",
                run, (double)per_run_preamp_cuts_[run].first,  (okA ? "fit" : "hard"),
                     (double)per_run_preamp_cuts_[run].second, (okC ? "fit" : "hard"));
        }
        std::cout << "DerivePerRunPreampCuts25: " << n_fit_ok << "/" << hA_map.size()
                  << " runs with successful fits on both sides." << std::endl;

        for (auto& [r, h] : hA_map) delete h;
        for (auto& [r, h] : hC_map) delete h;
    }

    // -------------------------------------------------------------------------
    void FillOneStage(StageAlt s, const EvDataAlt& ev) {
        hh_[s][kZdcFcal_A] ->Fill(ev.fcal_AC, ev.zdc_tot);
        hh_[s][kZdcTime_A] ->Fill(ev.tA,      ev.tC);
        hh_[s][kFCalAC_A]  ->Fill(ev.fcal_A,  ev.fcal_C);
        hh_[s][kPreampAC_A]->Fill(ev.preamp_A, ev.preamp_C);
        if (ev.ntrk_total > 0) {
            double frac = (double)ev.ntrk_tight / ev.ntrk_total;
            hh_[s][kFrac_A]    ->Fill(ev.ntrk_total, frac);
            hh_[s][kNtrkFcal_A]->Fill(ev.fcal_AC, ev.ntrk_tight);
        }
    }

    // -------------------------------------------------------------------------
    void FillHists() {
        TChain chain("HeavyIonD3PD", "HeavyIonD3PD");
        for (const auto& f : infiles_) {
            if (gSystem->AccessPathName(f.c_str())) {
                std::cerr << "File not found (skipping): " << f << std::endl; continue;
            }
            chain.Add(f.c_str());
        }
        if (chain.GetEntries() == 0)
            throw std::runtime_error("TChain empty for PbPb 20" + yr_);

        chain.SetMakeClass(1);
        chain.SetBranchStatus("*", 0);

        Int_t   b_HLT = 0, run_num_ev = 0;
        Float_t FCal_Et_P = 0.f, FCal_Et_N = 0.f;
        Float_t zdc_E[2] = {}, zdc_t[2] = {};
        Float_t preamp[2][4] = {};
        Int_t   centrality_i = 0;
        std::vector<int>* trk_numqual = nullptr;

        chain.SetBranchStatus("b_HLT_mu4_L1MU3V",         1);
        chain.SetBranchStatus("RunNumber",                  1);
        chain.SetBranchStatus("FCal_Et_P",                  1);
        chain.SetBranchStatus("FCal_Et_N",                  1);
        chain.SetBranchStatus("zdc_ZdcEnergy",              1);
        chain.SetBranchStatus("zdc_ZdcTime",                1);
        chain.SetBranchStatus("zdc_ZdcModulePreSampleAmp",  1);
        chain.SetBranchStatus("trk_numqual",                1);
        chain.SetBranchStatus("centrality",                 1);

        chain.SetBranchAddress("b_HLT_mu4_L1MU3V",         &b_HLT);
        chain.SetBranchAddress("RunNumber",                  &run_num_ev);
        chain.SetBranchAddress("FCal_Et_P",                 &FCal_Et_P);
        chain.SetBranchAddress("FCal_Et_N",                 &FCal_Et_N);
        chain.SetBranchAddress("zdc_ZdcEnergy",             zdc_E);
        chain.SetBranchAddress("zdc_ZdcTime",               zdc_t);
        chain.SetBranchAddress("zdc_ZdcModulePreSampleAmp", preamp);
        chain.SetBranchAddress("trk_numqual",               &trk_numqual);
        chain.SetBranchAddress("centrality",                &centrality_i);

        const Long64_t n_total = chain.GetEntries();
        std::cout << "FillHists: " << n_total << " total events (PbPb 20" << yr_ << ")" << std::endl;

        for (Long64_t i = 0; i < n_total; ++i) {
            chain.GetEntry(i);
            if (!b_HLT) continue;
            if (bad_runs_.count(run_num_ev)) continue;
            ++n_HLT_;

            EvDataAlt ev;
            ev.fcal_A   = FCal_Et_P * 1e-6f;
            ev.fcal_C   = FCal_Et_N * 1e-6f;
            ev.fcal_AC  = ev.fcal_A + ev.fcal_C;
            ev.zdc_tot  = (zdc_E[0] + zdc_E[1]) / 1000.f;
            ev.tA       = zdc_t[1];   // [1]=A
            ev.tC       = zdc_t[0];   // [0]=C
            ev.preamp_A = preamp[1][0] + preamp[1][1] + preamp[1][2] + preamp[1][3];  // [1]=A
            ev.preamp_C = preamp[0][0] + preamp[0][1] + preamp[0][2] + preamp[0][3];  // [0]=C
            ev.ntrk_total = 0; ev.ntrk_tight = 0;
            if (trk_numqual && (int)trk_numqual->size() >= 4) {
                ev.ntrk_total = (*trk_numqual)[0];
                ev.ntrk_tight = (*trk_numqual)[3];
            }
            // For years where 'centrality' branch is unfilled, recalculate from FCal ET
            // using the same Glauber table as MuonPairPbPb::UpdateCentrality().
            ev.centrality = (run_year_ == 25)
                            ? CentralityFromFCal2023(ev.fcal_AC)
                            : centrality_i;
            const bool is_ctr80 = (ev.centrality >= 0 && ev.centrality < 80);

            FillOneStage(kNoCut_A, ev);
            if (is_ctr80) {
                ++n_ctr80_[kNoCut_A];
                h_centr_before_cuts_->Fill(ev.centrality);
            }

            // Per-centrality ZDC time AC correlation (no cuts, PbPb2023 Glauber FCal thresholds)
            {
                const float f = ev.fcal_AC;
                if      (f >= kCtrThresh[0]) h_zdctime_ctr_[0]->Fill(ev.tA, ev.tC);
                else if (f >= kCtrThresh[1]) h_zdctime_ctr_[1]->Fill(ev.tA, ev.tC);
                else if (f >= kCtrThresh[2]) h_zdctime_ctr_[2]->Fill(ev.tA, ev.tC);
                else if (f >= kCtrThresh[3]) h_zdctime_ctr_[3]->Fill(ev.tA, ev.tC);
                else if (f >= kCtrThresh[4]) h_zdctime_ctr_[4]->Fill(ev.tA, ev.tC);
                else if (f >= kCtrThresh[5]) h_zdctime_ctr_[5]->Fill(ev.tA, ev.tC);
            }

            // Cut 1: ZDC banana
            const double cut1_val = PbPbEvSelEvalCut(g_cut1_, ev.fcal_AC);
            const bool   pass_c1  = (cut1_val <= 0. || ev.zdc_tot <= cut1_val);
            FillOneStage(pass_c1 ? kC1Pass_A : kC1Fail_A, ev);
            if (is_ctr80) ++n_ctr80_[pass_c1 ? kC1Pass_A : kC1Fail_A];
            if (!pass_c1) continue;

            // Cut 2: ZDC time box
            const bool pass_c2 = (std::abs(ev.tA) < CUT2_T_NS_ALT && std::abs(ev.tC) < CUT2_T_NS_ALT);
            FillOneStage(pass_c2 ? kC2Pass_A : kC2Fail_A, ev);
            if (is_ctr80) ++n_ctr80_[pass_c2 ? kC2Pass_A : kC2Fail_A];
            if (!pass_c2) continue;

            // Determine Cut 3 threshold: per-run mu+7sigma for yr25, hard cut for 23/24
            float c3_A = cut3_preamp_A_, c3_C = cut3_preamp_C_;
            if (run_year_ == 25 && !per_run_preamp_cuts_.empty()) {
                auto it = per_run_preamp_cuts_.find(run_num_ev);
                if (it != per_run_preamp_cuts_.end()) {
                    c3_A = it->second.first;
                    c3_C = it->second.second;
                }
            }

            // Part III: classify event by Cut-3 outcome on each side
            ++n_after_c12_;
            const bool fail_A = (ev.preamp_A >= c3_A);
            const bool fail_C = (ev.preamp_C >= c3_C);
            const int grp = (fail_A && fail_C) ? 2 : ((fail_A || fail_C) ? 1 : 0);
            h_preamp_grp_[grp]->Fill(ev.preamp_A, ev.preamp_C);
            ++n_preamp_grp_[grp];

            // Cut 3 decision: fail only if both sides exceed threshold
            const bool pass_c3 = !(fail_A && fail_C);
            FillOneStage(pass_c3 ? kC3Pass_A : kC3Fail_A, ev);
            if (is_ctr80) ++n_ctr80_[pass_c3 ? kC3Pass_A : kC3Fail_A];
            if (!pass_c3) {
                const StageAlt sub = (fail_A && fail_C) ? kC3FailBoth_A : kC3FailOne_A;
                FillOneStage(sub, ev);
                if (is_ctr80) ++n_ctr80_[sub];
                continue;
            }

            stored_.push_back(ev);
        }
        std::cout << "HLT events: " << n_HLT_ << ", passing cuts 1-3: " << stored_.size() << std::endl;
    }

    // -------------------------------------------------------------------------
    void DeriveCut4() {
        auto fits = FitSlices2DAlt(hh_[kC3Pass_A][kFrac_A], N_BINS_FRAC_SL_ALT);
        if (fits.empty()) { std::cerr << "DeriveCut4: no converged slices!" << std::endl; return; }
        std::vector<double> gx, gy_lo;
        for (const auto& f : fits) {
            gx.push_back(f.x);
            gy_lo.push_back(f.mu - CUT4_N_SIGMA_ALT * f.sig);
            std::cout << Form("  cut4 slice x=%.0f  mu=%.4f  sig=%.4f  cut_lo=%.4f\n",
                              f.x, f.mu, f.sig, f.mu - CUT4_N_SIGMA_ALT * f.sig);
        }
        g_cut4_ = new TGraph((int)gx.size(), gx.data(), gy_lo.data());
        g_cut4_->SetTitle(Form("nTrk HItight frac lower cut (mu-%.0fsigma);N_{{trk}}^{{total}};frac cut", CUT4_N_SIGMA_ALT));
        std::cout << "Cut 4 derived: " << g_cut4_->GetN() << " slices." << std::endl;
    }

    // -------------------------------------------------------------------------
    void FillCuts45() {
        if (!g_cut4_) { std::cerr << "FillCuts45: cut 4 not derived!" << std::endl; return; }

        for (const auto& ev : stored_) {
            if (ev.ntrk_total <= 0) continue;
            const double frac  = (double)ev.ntrk_tight / ev.ntrk_total;
            const double cut4  = PbPbEvSelEvalCut(g_cut4_, ev.ntrk_total);
            const bool pass_c4 = (frac >= cut4);
            FillOneStage(pass_c4 ? kC4Pass_A : kC4Fail_A, ev);
            if (ev.centrality >= 0 && ev.centrality < 80)
                ++n_ctr80_[pass_c4 ? kC4Pass_A : kC4Fail_A];
        }
        std::cout << "Cut 4 filled: pass=" << (Long64_t)hh_[kC4Pass_A][kZdcFcal_A]->Integral()
                  << " fail=" << (Long64_t)hh_[kC4Fail_A][kZdcFcal_A]->Integral() << std::endl;

        auto fits5 = FitSlices2DAlt(hh_[kC4Pass_A][kNtrkFcal_A], N_BINS_FCAL_SL_ALT);
        if (fits5.empty()) { std::cerr << "DeriveCut5: no converged slices!" << std::endl; return; }
        std::vector<double> gx, gy_lo, gy_hi, lin_x5, lin_y5_lo, lin_y5_hi;
        for (const auto& f : fits5) {
            gx.push_back(f.x);
            gy_lo.push_back(f.mu - CUT5_N_SIGMA_ALT * f.sig);
            gy_hi.push_back(f.mu + CUT5_N_SIGMA_ALT * f.sig);
            std::cout << Form("  cut5 slice x=%.3f  mu=%.1f  sig=%.1f  lo=%.1f  hi=%.1f\n",
                              f.x, f.mu, f.sig,
                              f.mu - CUT5_N_SIGMA_ALT * f.sig, f.mu + CUT5_N_SIGMA_ALT * f.sig);
            if (f.x > FCAL_CUT5_LIN_MIN_ALT) {
                lin_x5.push_back(f.x);
                lin_y5_lo.push_back(f.mu - CUT5_N_SIGMA_ALT * f.sig);
                lin_y5_hi.push_back(f.mu + CUT5_N_SIGMA_ALT * f.sig);
            }
        }
        if ((int)lin_x5.size() < 2) { std::cerr << "DeriveCut5: too few points for linear fit!" << std::endl; return; }
        TGraph g5_lo_raw((int)lin_x5.size(), lin_x5.data(), lin_y5_lo.data());
        TGraph g5_hi_raw((int)lin_x5.size(), lin_x5.data(), lin_y5_hi.data());
        TF1 f5a_lo("__f5a_lo", "pol1", lin_x5.front(), lin_x5.back());
        TF1 f5a_hi("__f5a_hi", "pol1", lin_x5.front(), lin_x5.back());
        g5_lo_raw.Fit(&f5a_lo, "RNQ");
        g5_hi_raw.Fit(&f5a_hi, "RNQ");
        std::cout << Form("  cut5 lin-fit lo: intercept=%.1f  slope=%.2f\n",
                          f5a_lo.GetParameter(0), f5a_lo.GetParameter(1));
        std::cout << Form("  cut5 lin-fit hi: intercept=%.1f  slope=%.2f\n",
                          f5a_hi.GetParameter(0), f5a_hi.GetParameter(1));
        const double fcal_lo_g5 = -0.5, fcal_hi_g5 = 5.5;
        const int n_lin5_pts = 12;
        std::vector<double> lg5x(n_lin5_pts), lg5y_lo(n_lin5_pts), lg5y_hi(n_lin5_pts);
        for (int k = 0; k < n_lin5_pts; ++k) {
            double x = fcal_lo_g5 + k * (fcal_hi_g5 - fcal_lo_g5) / (n_lin5_pts - 1);
            lg5x[k]    = x;
            lg5y_lo[k] = f5a_lo.Eval(x);
            lg5y_hi[k] = f5a_hi.Eval(x);
        }
        g_cut5_lo_ = new TGraph(n_lin5_pts, lg5x.data(), lg5y_lo.data());
        g_cut5_hi_ = new TGraph(n_lin5_pts, lg5x.data(), lg5y_hi.data());
        g_cut5_lo_->SetTitle(Form("nTrk HItight lower cut linear fit (FCal>%.1f TeV);FCal E_{{T}} [TeV];nTrk cut", FCAL_CUT5_LIN_MIN_ALT));
        g_cut5_hi_->SetTitle(Form("nTrk HItight upper cut linear fit (FCal>%.1f TeV);FCal E_{{T}} [TeV];nTrk cut", FCAL_CUT5_LIN_MIN_ALT));
        std::cout << "Cut 5 linear fit done: " << (int)lin_x5.size() << " points used." << std::endl;

        for (const auto& ev : stored_) {
            if (ev.ntrk_total <= 0) continue;
            const double frac    = (double)ev.ntrk_tight / ev.ntrk_total;
            const double cut4    = PbPbEvSelEvalCut(g_cut4_, ev.ntrk_total);
            if (frac < cut4) continue;
            const double cut5_lo = PbPbEvSelEvalCut(g_cut5_lo_, ev.fcal_AC);
            const double cut5_hi = PbPbEvSelEvalCut(g_cut5_hi_, ev.fcal_AC);
            const bool pass_c5   = (ev.ntrk_tight >= cut5_lo && ev.ntrk_tight <= cut5_hi);
            FillOneStage(pass_c5 ? kC5Pass_A : kC5Fail_A, ev);
            if (ev.centrality >= 0 && ev.centrality < 80) {
                ++n_ctr80_[pass_c5 ? kC5Pass_A : kC5Fail_A];
                if (pass_c5) h_centr_after_cuts_->Fill(ev.centrality);
            }
        }
        std::cout << "Cut 5 filled: pass=" << (Long64_t)hh_[kC5Pass_A][kZdcFcal_A]->Integral()
                  << " fail=" << (Long64_t)hh_[kC5Fail_A][kZdcFcal_A]->Integral() << std::endl;
    }

    // -------------------------------------------------------------------------
    void SaveCuts() {
        const std::string suffix = is_alt_ ? "_alt" : "";
        const std::string path = cuts_dir_ + "/event_sel_cuts_pbpb_20" + yr_ + suffix + ".root";
        // Alt mode: UPDATE to preserve alt_poly_a_b_c_x1 and other keys from FitAndPlotZDCAlt.
        // Nominal mode: RECREATE (clean slate each run).
        const char* open_mode = is_alt_ ? "UPDATE" : "RECREATE";
        TFile* f = TFile::Open(path.c_str(), open_mode);
        if (!f || f->IsZombie()) { std::cerr << "Cannot open " << path << std::endl; return; }
        auto wobj = [](TObject* o, const char* key) { o->Write(key, TObject::kOverwrite); };
        if (g_cut1_) wobj(g_cut1_->Clone(PbPbEvSelKey::kZDCFCalCut),    PbPbEvSelKey::kZDCFCalCut);
        if (g_cut4_) wobj(g_cut4_->Clone(PbPbEvSelKey::kNTrkFracCutLo), PbPbEvSelKey::kNTrkFracCutLo);
        if (g_cut5_lo_ && g_cut5_hi_) {
            wobj(g_cut5_lo_->Clone(PbPbEvSelKey::kNTrkFCalCutLo), PbPbEvSelKey::kNTrkFCalCutLo);
            wobj(g_cut5_hi_->Clone(PbPbEvSelKey::kNTrkFCalCutHi), PbPbEvSelKey::kNTrkFCalCutHi);
        }
        TParameter<double> p1(PbPbEvSelKey::kZDCTimeCutNs,  CUT2_T_NS_ALT);   p1.Write(PbPbEvSelKey::kZDCTimeCutNs,  TObject::kOverwrite);
        TParameter<double> p2(PbPbEvSelKey::kPreampACutADC, cut3_preamp_A_);   p2.Write(PbPbEvSelKey::kPreampACutADC, TObject::kOverwrite);
        TParameter<double> p3(PbPbEvSelKey::kPreampCCutADC, cut3_preamp_C_);   p3.Write(PbPbEvSelKey::kPreampCCutADC, TObject::kOverwrite);
        TParameter<double> p4("nTrk_frac_n_sigma",      CUT4_N_SIGMA_ALT);     p4.Write("nTrk_frac_n_sigma",          TObject::kOverwrite);
        TParameter<double> p5("nTrk_FCal_band_n_sigma", CUT5_N_SIGMA_ALT);     p5.Write("nTrk_FCal_band_n_sigma",     TObject::kOverwrite);

        // PbPb25: write per-run preamp cut TTree (branches: run_number, cut_A_ADC, cut_C_ADC)
        if (run_year_ == 25 && !per_run_preamp_cuts_.empty()) {
            f->cd();
            std::vector<int> runs;
            for (auto& kv : per_run_preamp_cuts_) runs.push_back(kv.first);
            std::sort(runs.begin(), runs.end());
            Int_t    run_num_w = 0;
            Double_t cut_a_w = 0., cut_c_w = 0.;
            TTree* t = new TTree(PbPbEvSelKey::kPreampPerRunTree,
                                 "Per-run ZDC preamp mu+7sigma cuts (PbPb25)");
            t->SetDirectory(nullptr);
            t->Branch("run_number", &run_num_w, "run_number/I");
            t->Branch("cut_A_ADC",  &cut_a_w,  "cut_A_ADC/D");
            t->Branch("cut_C_ADC",  &cut_c_w,  "cut_C_ADC/D");
            for (int r : runs) {
                run_num_w = r;
                cut_a_w   = per_run_preamp_cuts_.at(r).first;
                cut_c_w   = per_run_preamp_cuts_.at(r).second;
                t->Fill();
            }
            t->Write(PbPbEvSelKey::kPreampPerRunTree, TObject::kOverwrite);
            delete t;
            std::cout << "Per-run preamp cut TTree saved (" << runs.size() << " runs)." << std::endl;
        }

        f->Close();
        std::cout << "All cuts saved to: " << path << std::endl;
    }

    // -------------------------------------------------------------------------
    std::string OutPath(const std::string& stem) const {
        return out_dir_ + "/" + stem + "_pbpb_20" + yr_ + ".png";
    }

    void AddLabel(TLatex& tl, double x, double y) const {
        tl.DrawLatex(x, y, ("Pb+Pb 20" + yr_ + " data").c_str());
    }

    void DrawPad(int pad_idx, StageAlt s) {
        gPad->SetLogz();
        gPad->SetRightMargin(0.16);
        gPad->SetLeftMargin(0.13);
        gPad->SetBottomMargin(0.13);
        TH2D* h = hh_[s][pad_idx];
        h->SetContour(99);
        h->Draw("COLZ");
    }

    void DrawFivePanelCanvas(StageAlt s, const std::string& outpath) {
        TCanvas* c = new TCanvas(Form("c5a_%d", s), "", 1500, 900);
        c->Divide(3, 2, 0.003, 0.003);
        for (int t = 0; t < kNPanelTypes_A; ++t) {
            c->cd(t + 1);
            DrawPad(t, s);
        }
        c->SaveAs(outpath.c_str());
        delete c;
    }

    // Standalone cut 3: 1D preamp distributions with fixed cut line, no fitting
    void DrawStandaloneCut3() {
        TH2D* h2 = hh_[kC2Pass_A][kPreampAC_A];
        TH1D* hA = (TH1D*)h2->ProjectionX("__alt_preamp_A");
        TH1D* hC = (TH1D*)h2->ProjectionY("__alt_preamp_C");
        hA->SetDirectory(nullptr); hC->SetDirectory(nullptr);

        TH1D*       hArr[2]  = {hA,               hC};
        const char* sides[2] = {"A",               "C"};
        const int rebin = 3;
        TH1D* hd[2] = {};

        const float cutV[2] = {cut3_preamp_A_, cut3_preamp_C_};  // per-side cut values
        // Per-year zoom range for the standalone preamp plot
        double zoom_lo = -800., zoom_hi = 1500.;
        if      (run_year_ == 23) { zoom_lo = -500.; zoom_hi =  700.; }
        else if (run_year_ == 24) { zoom_lo = -800.; zoom_hi =  800.; }
        // year 25: keep [-800, 1500]

        TCanvas* c = new TCanvas("c_alt_cut3_sa", "", 1200, 600);
        c->Divide(2, 1, 0.005, 0.005);
        for (int i = 0; i < 2; ++i) {
            c->cd(i + 1);
            gPad->SetLogy();
            gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);

            hd[i] = (TH1D*)hArr[i]->Clone(Form("__alt_hd3_%d", i));
            hd[i]->Rebin(rebin);
            hd[i]->Scale(1.0 / rebin);
            hd[i]->GetXaxis()->SetRangeUser(zoom_lo, zoom_hi);

            // Find peak max within displayed window for y-axis scaling
            double lmax = 0.;
            const int b_lo = hd[i]->FindBin(zoom_lo);
            const int b_hi = hd[i]->FindBin(zoom_hi);
            for (int b = b_lo; b <= b_hi; ++b)
                lmax = std::max(lmax, hd[i]->GetBinContent(b));
            hd[i]->GetYaxis()->SetRangeUser(0.5, lmax * 3.0);
            hd[i]->GetXaxis()->SetTitle(Form("ZDC preamp sum side %s [ADC counts]", sides[i]));
            hd[i]->GetYaxis()->SetTitle("Events / 60 ADC");
            hd[i]->SetMarkerStyle(20); hd[i]->SetMarkerSize(0.5);
            hd[i]->Draw("E");

            TLine lcut(cutV[i], 0.5, cutV[i], lmax * 3.0);
            lcut.SetLineColor(kBlue+1); lcut.SetLineWidth(2); lcut.SetLineStyle(2);
            lcut.DrawClone();

            TLatex tl; tl.SetNDC(); tl.SetTextSize(0.038);
            tl.DrawLatex(0.17, 0.88, Form("Pb+Pb 20%s  side %s", yr_.c_str(), sides[i]));
            tl.SetTextColor(kBlue+1); tl.SetTextSize(0.033);
            tl.DrawLatex(0.17, 0.82, Form("cut = %.0f ADC", (double)cutV[i]));
            tl.SetTextColor(kBlack);
        }
        c->SaveAs(OutPath(Form("event_sel_cut3_%s_standalone", kPbPbEvSelCutLabel(3))).c_str());
        delete c;
        for (int i = 0; i < 2; ++i) delete hd[i];

        // Full-range preamp plot (no x-axis restriction)
        TH1D* hdf[2] = {};
        TCanvas* cf = new TCanvas("c_alt_cut3_sa_full", "", 1200, 600);
        cf->Divide(2, 1, 0.005, 0.005);
        for (int i = 0; i < 2; ++i) {
            cf->cd(i + 1);
            gPad->SetLogy();
            gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
            hdf[i] = (TH1D*)hArr[i]->Clone(Form("__alt_hdf3_%d", i));
            hdf[i]->Rebin(rebin);
            hdf[i]->Scale(1.0 / rebin);
            hdf[i]->GetXaxis()->SetTitle(Form("ZDC preamp sum side %s [ADC counts]", sides[i]));
            hdf[i]->GetYaxis()->SetTitle("Events / 60 ADC");
            hdf[i]->SetMarkerStyle(20); hdf[i]->SetMarkerSize(0.5);
            hdf[i]->Draw("E");

            TLine lfcut(cutV[i], hdf[i]->GetMinimum(1e-3), cutV[i], hdf[i]->GetMaximum() * 3.0);
            lfcut.SetLineColor(kBlue+1); lfcut.SetLineWidth(2); lfcut.SetLineStyle(2);
            lfcut.DrawClone();

            TLatex tl; tl.SetNDC(); tl.SetTextSize(0.038);
            tl.DrawLatex(0.17, 0.88, Form("Pb+Pb 20%s  side %s  (full range)", yr_.c_str(), sides[i]));
            tl.SetTextColor(kBlue+1); tl.SetTextSize(0.033);
            tl.DrawLatex(0.17, 0.82, Form("cut = %.0f ADC", (double)cutV[i]));
            tl.SetTextColor(kBlack);
        }
        cf->SaveAs(OutPath(Form("event_sel_cut3_%s_standalone_fullrange", kPbPbEvSelCutLabel(3))).c_str());
        delete cf;
        for (int i = 0; i < 2; ++i) delete hdf[i];

        delete hA; delete hC;
    }

    // Standalone 2D cut canvas for cut index 1-5
    void DrawStandaloneCut(int cut_idx) {
        if (cut_idx == 3) { DrawStandaloneCut3(); return; }

        static const HTypeAlt cut_htype[6] = {
            kZdcFcal_A, kZdcFcal_A, kZdcTime_A, kPreampAC_A, kFrac_A, kNtrkFcal_A
        };
        static const StageAlt src_stage[6] = {
            kNoCut_A, kNoCut_A, kC1Pass_A, kC2Pass_A, kC3Pass_A, kC4Pass_A
        };
        HTypeAlt  ht = cut_htype[cut_idx];
        StageAlt  ss = src_stage[cut_idx];
        TH2D* hd = (TH2D*)hh_[ss][ht]->Clone(Form("__alt_hsa_%d", cut_idx));
        hd->SetDirectory(nullptr);

        TCanvas* c = new TCanvas(Form("c_alt_cut%d", cut_idx), "", 800, 650);
        c->SetLogz(); c->SetRightMargin(0.14);
        hd->SetContour(99); hd->Draw("COLZ");

        TLatex tl; tl.SetNDC(); tl.SetTextSize(0.036);
        AddLabel(tl, 0.15, 0.88);

        if (cut_idx == 1 && g_cut1_) {
            TGraph* gd = (TGraph*)g_cut1_->Clone("__alt_gc1");
            gd->SetLineColor(kRed); gd->SetLineWidth(2);
            gd->SetMarkerColor(kRed); gd->SetMarkerStyle(20); gd->SetMarkerSize(0.5);
            gd->Draw("LP SAME");
            tl.DrawLatex(0.15, 0.82, "red: ZDC banana cut (#mu+5#sigma)");
        } else if (cut_idx == 2) {
            TLine l; l.SetLineColor(kRed); l.SetLineWidth(2); l.SetLineStyle(2);
            l.DrawLine(-CUT2_T_NS_ALT, -10., -CUT2_T_NS_ALT, 10.);
            l.DrawLine( CUT2_T_NS_ALT, -10.,  CUT2_T_NS_ALT, 10.);
            l.DrawLine(-10., -CUT2_T_NS_ALT, 10., -CUT2_T_NS_ALT);
            l.DrawLine(-10.,  CUT2_T_NS_ALT, 10.,  CUT2_T_NS_ALT);
            tl.DrawLatex(0.15, 0.82, Form("red: |t^{A}| < %.1f ns, |t^{C}| < %.1f ns",
                                          CUT2_T_NS_ALT, CUT2_T_NS_ALT));
        } else if (cut_idx == 4 && g_cut4_) {
            TGraph* gd = (TGraph*)g_cut4_->Clone("__alt_gc4");
            gd->SetLineColor(kRed); gd->SetLineWidth(2);
            gd->SetMarkerColor(kRed); gd->SetMarkerStyle(20); gd->SetMarkerSize(0.5);
            gd->Draw("LP SAME");
            tl.DrawLatex(0.15, 0.82, Form("red: frac lower cut (#mu-%.0f#sigma)", CUT4_N_SIGMA_ALT));
        } else if (cut_idx == 5 && g_cut5_lo_ && g_cut5_hi_) {
            TGraph* glo = (TGraph*)g_cut5_lo_->Clone("__alt_gc5lo");
            TGraph* ghi = (TGraph*)g_cut5_hi_->Clone("__alt_gc5hi");
            glo->SetLineColor(kRed); glo->SetLineWidth(2);
            glo->SetMarkerColor(kRed); glo->SetMarkerStyle(20); glo->SetMarkerSize(0.5);
            ghi->SetLineColor(kRed); ghi->SetLineWidth(2);
            ghi->SetMarkerColor(kRed); ghi->SetMarkerStyle(20); ghi->SetMarkerSize(0.5);
            glo->Draw("LP SAME"); ghi->Draw("LP SAME");
            tl.DrawLatex(0.15, 0.82, Form("red: linear fit (#mu#pm%.0f#sigma, FCal>%.1f TeV)", CUT5_N_SIGMA_ALT, FCAL_CUT5_LIN_MIN_ALT));
        }

        c->SaveAs(OutPath(Form("event_sel_cut%d_%s_standalone", cut_idx, kPbPbEvSelCutLabel(cut_idx))).c_str());
        delete c; delete hd;
    }

    // -------------------------------------------------------------------------
    void DrawSurvivalPlot() {
        static const StageAlt pass_stages[5] = {kC1Pass_A, kC2Pass_A, kC3Pass_A, kC4Pass_A, kC5Pass_A};
        double n_before[5];
        n_before[0] = (double)n_ctr80_[kNoCut_A];
        for (int k = 1; k < 5; ++k)
            n_before[k] = (double)n_ctr80_[pass_stages[k-1]];

        double xv[5], yv[5], exv[5], eyv[5];
        for (int k = 0; k < 5; ++k) {
            double nk = (double)n_ctr80_[pass_stages[k]];
            double nb = n_before[k];
            double p  = (nb > 0) ? nk / nb : 0.;
            xv[k]  = k + 1; yv[k] = p; exv[k] = 0.;
            eyv[k] = (nb > 0) ? std::sqrt(p * (1. - p) / nb) : 0.;
            std::cout << Form("  Cut %d (%s): %.4f survival (%lld / %lld ctr<80)  err=%.2e\n",
                              k+1, kPbPbEvSelCutLabel(k+1), p, (Long64_t)nk, (Long64_t)nb, eyv[k]);
        }

        TH1D* hframe = new TH1D("h_surv_frame_alt",
            ";;Survival fraction per cut (centrality 0#minus79)",
            5, 0.5, 5.5);
        hframe->SetDirectory(nullptr);
        for (int k = 0; k < 5; ++k)
            hframe->GetXaxis()->SetBinLabel(k + 1, Form("C%d  %s", k+1, kPbPbEvSelCutLabel(k+1)));
        hframe->SetFillStyle(0); hframe->SetLineColor(0);
        hframe->GetYaxis()->SetRangeUser(0.92, 1.002);
        hframe->GetXaxis()->SetLabelSize(0.038);

        TGraphErrors* gr = new TGraphErrors(5, xv, yv, exv, eyv);
        gr->SetMarkerStyle(20); gr->SetMarkerSize(1.2);
        gr->SetMarkerColor(kBlue+1); gr->SetLineColor(kBlue+1); gr->SetLineWidth(1);

        TCanvas* c = new TCanvas("c_alt_survival", "", 800, 600);
        c->SetLeftMargin(0.14); c->SetBottomMargin(0.18);
        hframe->Draw("AXIS");
        gr->Draw("PE SAME");

        TLatex tl; tl.SetNDC(); tl.SetTextSize(0.032);
        AddLabel(tl, 0.63, 0.80);
        tl.SetTextSize(0.027);
        for (int k = 0; k < 5; ++k)
            tl.DrawLatex(0.63, 0.74 - 0.055*k, Form("C%d: %.3f%%", k+1, yv[k]*100.));

        c->SaveAs(OutPath("event_sel_survival").c_str());
        delete c; delete hframe; delete gr;
    }

    // -------------------------------------------------------------------------
    void DrawCumulativeSurvivalPlot() {
        static const StageAlt pass_stages[5] = {kC1Pass_A, kC2Pass_A, kC3Pass_A, kC4Pass_A, kC5Pass_A};
        const double n_total = (double)n_ctr80_[kNoCut_A];

        double xv[5], yv[5], exv[5], eyv[5];
        for (int k = 0; k < 5; ++k) {
            double nk = (double)n_ctr80_[pass_stages[k]];
            double p  = (n_total > 0) ? nk / n_total : 0.;
            xv[k]  = k + 1; yv[k] = p; exv[k] = 0.;
            eyv[k] = (n_total > 0) ? std::sqrt(p * (1. - p) / n_total) : 0.;
            std::cout << Form("  Cut %d (%s): cumulative %.4f (%lld / %lld ctr<80)\n",
                              k+1, kPbPbEvSelCutLabel(k+1), p, (Long64_t)nk, (Long64_t)n_total);
        }

        double y_min = *std::min_element(yv, yv + 5);
        TH1D* hframe = new TH1D("h_surv_accum_frame_alt",
            ";;Cumulative survival fraction (centrality 0#minus79)",
            5, 0.5, 5.5);
        hframe->SetDirectory(nullptr);
        for (int k = 0; k < 5; ++k)
            hframe->GetXaxis()->SetBinLabel(k + 1, Form("C%d  %s", k+1, kPbPbEvSelCutLabel(k+1)));
        hframe->SetFillStyle(0); hframe->SetLineColor(0);
        hframe->GetYaxis()->SetRangeUser(std::max(0., y_min - 0.02), 1.002);
        hframe->GetXaxis()->SetLabelSize(0.038);

        TGraphErrors* gr = new TGraphErrors(5, xv, yv, exv, eyv);
        gr->SetMarkerStyle(20); gr->SetMarkerSize(1.2);
        gr->SetMarkerColor(kBlue+1); gr->SetLineColor(kBlue+1); gr->SetLineWidth(1);

        TCanvas* c = new TCanvas("c_alt_survival_accum", "", 800, 600);
        c->SetLeftMargin(0.14); c->SetBottomMargin(0.18);
        hframe->Draw("AXIS");
        gr->Draw("PE SAME");

        TLatex tl; tl.SetNDC(); tl.SetTextSize(0.032);
        AddLabel(tl, 0.63, 0.80);
        tl.SetTextSize(0.027);
        for (int k = 0; k < 5; ++k)
            tl.DrawLatex(0.63, 0.74 - 0.055*k, Form("C%d: %.3f%%", k+1, yv[k]*100.));

        c->SaveAs(OutPath("event_sel_survival_accum").c_str());
        delete c; delete hframe; delete gr;
    }

    // -------------------------------------------------------------------------
    void DrawZdcTimeCorrByCentr() {
        // 2x3 canvas: pads 1-5 = top 1% each (FCal ET-based), pad 6 = 5-10%
        static const char* ctr_labels[6] = {"0-1%","1-2%","2-3%","3-4%","4-5%","5-10%"};
        printf("ZDC time corr per centrality bin GetEntries (PbPb2023 Glauber thresholds):\n");
        for (int k = 0; k < 6; ++k)
            printf("  bin %d (%s): %.0f entries  FCal lower=%.5f TeV\n",
                   k, ctr_labels[k], h_zdctime_ctr_[k]->GetEntries(), (double)kCtrThresh[k]);
        TCanvas* c = new TCanvas("c_zdctime_ctr", "", 1500, 900);
        c->Divide(3, 2, 0.003, 0.003);
        for (int k = 0; k < 6; ++k) {
            c->cd(k + 1);
            gPad->SetLogz();
            gPad->SetRightMargin(0.16);
            gPad->SetLeftMargin(0.13);
            gPad->SetBottomMargin(0.13);
            h_zdctime_ctr_[k]->SetContour(99);
            h_zdctime_ctr_[k]->SetTitle(Form(";ZDC t^{A} [ns];ZDC t^{C} [ns]"));
            h_zdctime_ctr_[k]->Draw("COLZ");

            // Time-box cut lines
            TLine l; l.SetLineColor(kRed); l.SetLineWidth(2); l.SetLineStyle(2);
            l.DrawLine(-CUT2_T_NS_ALT, -10., -CUT2_T_NS_ALT, 10.);
            l.DrawLine( CUT2_T_NS_ALT, -10.,  CUT2_T_NS_ALT, 10.);
            l.DrawLine(-10., -CUT2_T_NS_ALT, 10., -CUT2_T_NS_ALT);
            l.DrawLine(-10.,  CUT2_T_NS_ALT, 10.,  CUT2_T_NS_ALT);

            TLatex tl; tl.SetNDC(); tl.SetTextSize(0.044);
            tl.DrawLatex(0.17, 0.89, Form("Pb+Pb 20%s  ctr %s (FCal)", yr_.c_str(), ctr_labels[k]));
            tl.SetTextSize(0.036);
            if (k == 0)
                tl.DrawLatex(0.17, 0.83, Form("FCal > %.3f TeV", (double)kCtrThresh[0]));
            else if (k < 5)
                tl.DrawLatex(0.17, 0.83, Form("FCal [%.3f, %.3f) TeV", (double)kCtrThresh[k], (double)kCtrThresh[k-1]));
            else
                tl.DrawLatex(0.17, 0.83, Form("FCal [%.3f, %.3f) TeV", (double)kCtrThresh[5], (double)kCtrThresh[4]));
        }
        c->SaveAs(OutPath("ZDC_time_AC_corr_top5_centr").c_str());
        delete c;
    }

    // -------------------------------------------------------------------------
    void DrawCentralityRatio() {
        TH1D* ratio = (TH1D*)h_centr_after_cuts_->Clone("h_centr_ratio");
        ratio->SetDirectory(nullptr);
        ratio->Divide(h_centr_before_cuts_);
        ratio->SetTitle(";Centrality percentile;(after cuts) / (before cuts)");
        ratio->SetMarkerStyle(20);
        ratio->SetMarkerSize(0.7);

        TCanvas* c = new TCanvas("c_centr_ratio", "", 800, 600);
        c->SetLeftMargin(0.13); c->SetBottomMargin(0.13);
        ratio->GetYaxis()->SetRangeUser(0., 1.05);
        ratio->Draw("E");
        TLine lone(-.5, 1., 79.5, 1.);
        lone.SetLineColor(kRed); lone.SetLineStyle(2); lone.SetLineWidth(1);
        lone.DrawClone();
        TLatex tl; tl.SetNDC(); tl.SetTextSize(0.035);
        AddLabel(tl, 0.17, 0.25);
        if (run_year_ == 25)
            tl.DrawLatex(0.17, 0.19, "Centrality from FCal E_{T} (pbpb2023 Glauber thresholds)");
        c->SaveAs(OutPath("centrality_ratio_after_before_cuts").c_str());
        delete c; delete ratio;
    }

    // -------------------------------------------------------------------------
    void DrawCentralityBeforeCuts() {
        TCanvas* c = new TCanvas("c_centr_before_cuts", "", 800, 600);
        c->SetLeftMargin(0.13); c->SetBottomMargin(0.13);
        h_centr_before_cuts_->SetMarkerStyle(20);
        h_centr_before_cuts_->SetMarkerSize(0.7);
        h_centr_before_cuts_->Draw("E");
        TLatex tl; tl.SetNDC(); tl.SetTextSize(0.035);
        AddLabel(tl, 0.17, 0.88);
        if (run_year_ == 25)
            tl.DrawLatex(0.17, 0.82, "Centrality from FCal E_{T} (pbpb2023 Glauber thresholds)");
        c->SaveAs(OutPath("centrality_before_cuts").c_str());
        delete c;
    }

    // -------------------------------------------------------------------------
    void DrawCentralityAfterCuts() {
        TCanvas* c = new TCanvas("c_centr_after_cuts", "", 800, 600);
        c->SetLeftMargin(0.13); c->SetBottomMargin(0.13);
        h_centr_after_cuts_->SetMarkerStyle(20);
        h_centr_after_cuts_->SetMarkerSize(0.7);
        h_centr_after_cuts_->Draw("E");
        TLatex tl; tl.SetNDC(); tl.SetTextSize(0.035);
        AddLabel(tl, 0.17, 0.88);
        if (run_year_ == 25)
            tl.DrawLatex(0.17, 0.82, "Centrality from FCal E_{T} (pbpb2023 Glauber thresholds)");
        c->SaveAs(OutPath("centrality_after_cuts").c_str());
        delete c;
    }

    // -------------------------------------------------------------------------
    // Part III: 3-panel preamp A vs C correlation for Cut-3 pass/fail groups.
    // For each year: events after Cuts 1+2 are split into:
    //   [0] Both A and C pass Cut 3   [1] Exactly one side fails   [2] Both fail
    void DrawPreampCorrelationGroups() {
        if (n_after_c12_ == 0) {
            std::cerr << "DrawPreampCorrelationGroups: no events after Cuts 1+2" << std::endl;
            return;
        }

        gStyle->SetOptStat(0);
        gStyle->SetPalette(kBird);

        static const char* grp_title[3] = {
            "Both A & C pass Cut 3",
            "Exactly one side fails Cut 3",
            "Both A & C fail Cut 3"
        };
        const bool is_perrun = (run_year_ == 25);
        const std::string cut_desc = is_perrun
            ? "per-run #mu+7#sigma"
            : Form("A=%.0f ADC, C=%.0f ADC", (double)cut3_preamp_A_, (double)cut3_preamp_C_);

        // Zoom range: show a bit beyond the hard cut; cap at histogram max
        const double zoom_lo = -800.;
        const double zoom_hi = std::min(3000., (double)std::max(cut3_preamp_A_, cut3_preamp_C_) * 2.5 + 300.);

        TCanvas* cv = new TCanvas(Form("c_pgrp_%d", run_year_), "", 2100, 700);
        cv->Divide(3, 1, 0.004, 0.004);

        for (int g = 0; g < 3; ++g) {
            cv->cd(g + 1);
            gPad->SetLogz();
            gPad->SetRightMargin(0.16);
            gPad->SetLeftMargin(0.13);
            gPad->SetBottomMargin(0.13);
            gPad->SetTopMargin(0.17);

            TH2D* hd = (TH2D*)h_preamp_grp_[g]->Clone(
                Form("hd_pgrp_%d_%d", run_year_, g));
            hd->GetXaxis()->SetRangeUser(zoom_lo, zoom_hi);
            hd->GetYaxis()->SetRangeUser(zoom_lo, zoom_hi);
            hd->SetContour(99);
            hd->SetMinimum(1.);
            hd->Draw("COLZ");

            // Cut lines for fixed cuts (23/24)
            if (!is_perrun) {
                TLine lA(cut3_preamp_A_, zoom_lo, cut3_preamp_A_, zoom_hi);
                TLine lC(zoom_lo, cut3_preamp_C_, zoom_hi, cut3_preamp_C_);
                for (TLine* l : {&lA, &lC}) {
                    l->SetLineColor(kRed); l->SetLineWidth(1); l->SetLineStyle(2);
                    l->DrawClone();
                }
            }

            double pct = 100. * n_preamp_grp_[g] / (double)n_after_c12_;

            // Text box at top of pad
            TLatex tl; tl.SetNDC(); tl.SetTextAlign(13);
            tl.SetTextSize(0.046);
            tl.DrawLatex(0.14, 0.98, Form("Pb+Pb 20%s  |  %s", yr_.c_str(), grp_title[g]));
            tl.SetTextSize(0.040);
            tl.DrawLatex(0.14, 0.92,
                Form("%.3f%%  (%lld / %lld events)", pct, n_preamp_grp_[g], n_after_c12_));
            tl.DrawLatex(0.14, 0.87, Form("[Cut 3: %s]", cut_desc.c_str()));
        }

        cv->SaveAs(OutPath("event_sel_cut3_preamp_AC_corr_groups").c_str());
        std::cout << "Saved: " << OutPath("event_sel_cut3_preamp_AC_corr_groups") << std::endl;
        delete cv;
        for (int g = 0; g < 3; ++g) {
            std::cout << Form("  Group %d (%s): %lld events (%.3f%%)\n",
                g, grp_title[g], n_preamp_grp_[g], 100.*n_preamp_grp_[g]/(double)n_after_c12_);
        }

        for (int g = 0; g < 3; ++g) { delete h_preamp_grp_[g]; h_preamp_grp_[g] = nullptr; }
    }

    // -------------------------------------------------------------------------
    void SavePlots() {
        gStyle->SetOptStat(0);
        gStyle->SetPalette(kBird);

        DrawFivePanelCanvas(kNoCut_A, OutPath("event_sel_nocuts_5panel"));

        static const StageAlt pass_stages[5] = {kC1Pass_A, kC2Pass_A, kC3Pass_A, kC4Pass_A, kC5Pass_A};
        static const StageAlt fail_stages[5] = {kC1Fail_A, kC2Fail_A, kC3Fail_A, kC4Fail_A, kC5Fail_A};
        for (int cut = 1; cut <= 5; ++cut) {
            DrawStandaloneCut(cut);
            DrawFivePanelCanvas(pass_stages[cut-1],
                OutPath(Form("event_sel_cut%d_%s_pass_5panel", cut, kPbPbEvSelCutLabel(cut))));
            DrawFivePanelCanvas(fail_stages[cut-1],
                OutPath(Form("event_sel_cut%d_%s_fail_5panel", cut, kPbPbEvSelCutLabel(cut))));
        }

        DrawFivePanelCanvas(kC3FailOne_A,
            OutPath("event_sel_cut3_ZDC_preamp_fail_one_side_5panel"));
        DrawFivePanelCanvas(kC3FailBoth_A,
            OutPath("event_sel_cut3_ZDC_preamp_fail_both_sides_5panel"));

        DrawSurvivalPlot();
        DrawCumulativeSurvivalPlot();
        DrawZdcTimeCorrByCentr();
        DrawCentralityBeforeCuts();
        DrawCentralityAfterCuts();
        DrawCentralityRatio();
        DrawPreampCorrelationGroups();

        std::cout << "All plots saved to: " << out_dir_ << std::endl;
    }
};

// =============================================================================
void plot_pbpb_event_sel_cuts(int run_year = 24) {
    PbPbEventSelCutsAlt plotter(run_year, /*is_alt=*/false);
    plotter.Run();
}

// Alt banana-cut variant: loads Cut 1 from *_alt.root, saves all derived cuts
// there (UPDATE), and writes plots to event_selection_alternative_banana/.
void plot_pbpb_event_sel_cuts_alt(int run_year = 24) {
    PbPbEventSelCutsAlt plotter(run_year, /*is_alt=*/true);
    plotter.Run();
}
