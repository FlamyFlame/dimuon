// plot_pbpb_event_sel_cuts.cxx
// Sequential 5-cut event selection pipeline for PbPb data.
//   Cut 1: ZDC vs FCal banana (loaded from event_sel_cuts_pbpb_YYYY.root)
//   Cut 2: ZDC time |tA| < 2 ns AND |tC| < 2 ns
//   Cut 3: ZDC preamp sum: A < +181 ADC AND C < +205 ADC
//   Cut 4: nTrk HItight fraction vs total nTrk — lower bound at mu-5sigma (per-slice Gauss)
//   Cut 5: nTrk HItight vs FCal ET — band [mu-5sigma, mu+5sigma] (per-slice Gauss)
// Derived cuts 4+5 saved to event_sel_cuts_pbpb_YYYY.root.
//
// Output per cut: standalone 2D-with-cut canvas, 5-panel pass canvas, 5-panel fail canvas.
// Additionally: 5-panel no-cuts canvas, per-cut survival rate plot.
//
// Events passing cuts 1-3 are cached in memory to allow correct per-event fills
// for cuts 4+5 (preserves correlations across histogram variables).
//
// Usage:
//   .L plot_pbpb_event_sel_cuts.cxx+
//   plot_pbpb_event_sel_cuts(24)

#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include "TChain.h"
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
#include "TLine.h"
#include "TPaveText.h"
#include "TBox.h"
#include "TParameter.h"
#include "TMath.h"
#include "../../NTupleProcessingCode/PbPbEventSelConfig.h"

// ---- Cut tuning -------------------------------------------------------------
static const double CUT2_T_NS        = 1.5;    // ZDC time window [ns]
// Cut 3 preamp bounds derived from EMG fit (mu+6sigma) — set in FitPreamp1D()
static const double CUT4_N_SIGMA     = 5.0;    // nTrk frac lower cut
static const double CUT5_N_SIGMA     = 5.0;    // nTrk-FCal band half-width

static const int    N_BINS_FRAC_SL   = 4;      // nTrk_total bins per slice (cut 4)
static const int    N_BINS_FCAL_SL   = 2;      // FCal bins per slice (cut 5)
static const int    MIN_ENTRIES_FIT  = 100;
static const double FCAL_CUT5_LIN_MIN = 0.2;  // FCal lower bound for cut-5 linear fit [TeV]

// ---- Input file map ---------------------------------------------------------
static std::map<int, std::vector<std::string>> BuildFileMap() {
    const std::string base = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/";
    return {
        {23, {
            base + "pbpb_2023/data_pbpb23_part1.root",
            base + "pbpb_2023/data_pbpb23_part2.root",
            base + "pbpb_2023/data_pbpb23_part3.root",
            base + "pbpb_2023/data_pbpb23_part4.root",
            base + "pbpb_2023/data_pbpb23_part5.root",
            base + "pbpb_2023/data_pbpb23_part6.root",
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
enum Stage {
    kNoCut  = 0,
    kC1Pass = 1, kC1Fail = 2,
    kC2Pass = 3, kC2Fail = 4,
    kC3Pass = 5, kC3Fail = 6,
    kC4Pass = 7, kC4Fail = 8,
    kC5Pass = 9, kC5Fail = 10,
    kNStages = 11
};
// Types 0-4 appear in the 5-panel canvases; type 5 (kPreampAC) is standalone-cut-3 only.
enum HType {
    kZdcFcal  = 0,  // ZDC E_tot vs FCal A+C
    kZdcTime  = 1,  // ZDC time A vs C
    kFCalAC   = 2,  // FCal ET A vs C
    kFrac     = 3,  // nTrk HItight fraction vs total nTrk
    kNtrkFcal = 4,  // nTrk HItight vs FCal A+C
    kPreampAC = 5,  // ZDC preamp A vs C  (standalone cut 3 canvas only)
    kNHTypes  = 6
};
static const int kNPanelTypes = 5;  // only types 0-4 go into 5-panel canvases

static const char* kStageLabel[kNStages] = {
    "No cuts",
    "Pass cut 1 (ZDC banana)", "Fail cut 1 (ZDC banana)",
    "Pass cut 2 (ZDC time)",   "Fail cut 2 (ZDC time)",
    "Pass cut 3 (ZDC preamp)", "Fail cut 3 (ZDC preamp)",
    "Pass cut 4 (nTrk frac)",  "Fail cut 4 (nTrk frac)",
    "Pass cut 5 (nTrk-FCal)",  "Fail cut 5 (nTrk-FCal)"
};

// ---- Helpers ----------------------------------------------------------------

struct SliceFit { double x, mu, sig, ex; };

static std::vector<SliceFit> FitSlices2D(TH2D* h2, int n_bins_per_slice) {
    std::vector<SliceFit> result;
    const int nx = h2->GetNbinsX();
    const int n_slices = nx / n_bins_per_slice;
    TF1 gf("__gf_slice", "gaus");
    for (int isl = 0; isl < n_slices; ++isl) {
        const int ix_lo = isl * n_bins_per_slice + 1;
        const int ix_hi = std::min((isl + 1) * n_bins_per_slice, nx);
        TH1D* hp = (TH1D*)h2->ProjectionY(Form("__hp_%s_%d", h2->GetName(), isl), ix_lo, ix_hi);
        hp->SetDirectory(nullptr);
        if (hp->GetEntries() < MIN_ENTRIES_FIT) { delete hp; continue; }
        const double mu0 = hp->GetBinCenter(hp->GetMaximumBin());
        const double rms = (hp->GetRMS() > 0.) ? hp->GetRMS() : 1.0;
        const double ymin = hp->GetXaxis()->GetXmin();
        const double ymax = hp->GetXaxis()->GetXmax();
        // Pass 1
        gf.SetRange(std::max(ymin, mu0 - 3.*rms), std::min(ymax, mu0 + 3.*rms));
        gf.SetParameters(hp->GetMaximum(), mu0, rms);
        if (hp->Fit(&gf, "RNQ") != 0 || gf.GetParameter(2) <= 0.) { delete hp; continue; }
        double mu1 = gf.GetParameter(1), sig1 = gf.GetParameter(2);
        // Pass 2
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

// Exponentially-Modified Gaussian: p[0]=A, p[1]=mu, p[2]=sigma, p[3]=lambda (>0, right tail)
static double EMGFunc(double* x, double* p) {
    double A   = p[0];
    double mu  = p[1];
    double sig = std::abs(p[2]);
    double lam = std::abs(p[3]);
    return A * TMath::Exp(lam*(mu + 0.5*lam*sig*sig - x[0]))
             * TMath::Erfc((mu + lam*sig*sig - x[0]) / (sig * TMath::Sqrt2()));
}

// Cut labels are provided by kPbPbEvSelCutLabel() from PbPbEventSelConfig.h.

// ---- Per-event data cache (events passing cuts 1-3) -------------------------
struct EvData {
    float fcal_A, fcal_C, fcal_AC, zdc_tot, tA, tC, preamp_A, preamp_C;
    int   ntrk_total, ntrk_tight;
    int   centrality;
};

// =============================================================================
class PbPbEventSelCuts {
public:
    explicit PbPbEventSelCuts(int run_year) : run_year_(run_year % 2000) {
        yr_      = std::to_string(run_year_);
        out_dir_ = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/"
                   "single_b_analysis/event_selection/pbpb_20" + yr_;
        cuts_dir_= "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_20" + yr_;
        auto fm = BuildFileMap();
        if (fm.find(run_year_) == fm.end())
            throw std::runtime_error("No files for year " + yr_);
        infiles_ = fm.at(run_year_);
    }

    void Run() {
        gStyle->SetOptStat(0);
        gStyle->SetPalette(kBird);
        gSystem->mkdir(out_dir_.c_str(), true);
        BookHists();
        LoadCut1();
        FillHists();
        FitPreamp1D();
        ApplyCut3();
        DeriveCut4();
        FillCuts45();
        SaveCuts();
        SavePlots();
    }

private:
    int run_year_;
    std::string yr_, out_dir_, cuts_dir_;
    std::vector<std::string> infiles_;

    TH2D* hh_[kNStages][kNHTypes] = {};
    TGraph* g_cut1_    = nullptr;
    TGraph* g_cut4_    = nullptr;
    TGraph* g_cut5_lo_ = nullptr;
    TGraph* g_cut5_hi_ = nullptr;

    TH1D* h1_preamp_A_ = nullptr;
    TH1D* h1_preamp_C_ = nullptr;
    TF1*  fit_preamp_A_ = nullptr;
    TF1*  fit_preamp_C_ = nullptr;

    Long64_t n_HLT_              = 0;
    Long64_t n_ctr80_[kNStages]  = {};  // per-stage count for centrality 0-79 only
    float cut3_preamp_A_         = 0.f; // derived from EMG fit: mu+6sigma, side A
    float cut3_preamp_C_         = 0.f; // derived from EMG fit: mu+6sigma, side C
    std::vector<EvData> stored_pre_c3_; // events passing cuts 1-2 (for cut 3 derivation)
    std::vector<EvData> stored_;        // events passing cuts 1-3

    // -------------------------------------------------------------------------
    void BookHists() {
        // Histogram specs per type
        struct Spec { int nx; double xlo, xhi; int ny; double ylo, yhi;
                      const char* xtitle; const char* ytitle; };
        Spec specs[kNHTypes] = {
            {120, -0.5, 5.5,   100,     0.,  400., "FCal E_{T}^{A+C} [TeV]",           "ZDC E_{total} [TeV]"},
            {100, -10., 10.,   100,   -10.,   10., "ZDC t^{A} [ns]",                   "ZDC t^{C} [ns]"},
            {110, -0.5, 3.0,   110,   -0.5,   3.0, "FCal E_{T}^{A} [TeV]",             "FCal E_{T}^{C} [TeV]"},
            {100,   0., 5000., 110,     0.,   1.1, "N_{trk}^{total} (p_{T}>400 MeV)",  "N_{trk}^{HItight} / N_{trk}^{total}"},
            {120, -0.5, 5.5,   100,     0., 3500., "FCal E_{T}^{A+C} [TeV]",           "N_{trk}^{HItight}"},
            {200,-1000., 3000., 200, -1000., 3000., "ZDC preamp sum side A [ADC]",      "ZDC preamp sum side C [ADC]"},
        };
        for (int s = 0; s < kNStages; ++s) {
            for (int t = 0; t < kNHTypes; ++t) {
                const Spec& sp = specs[t];
                hh_[s][t] = new TH2D(Form("hh_%d_%d", s, t),
                    Form(";%s;%s", sp.xtitle, sp.ytitle),
                    sp.nx, sp.xlo, sp.xhi, sp.ny, sp.ylo, sp.yhi);
                hh_[s][t]->SetDirectory(nullptr);
            }
        }
    }

    // -------------------------------------------------------------------------
    void LoadCut1() {
        const std::string path = cuts_dir_ + "/event_sel_cuts_pbpb_20" + yr_ + ".root";
        TFile* f = TFile::Open(path.c_str(), "READ");
        if (!f || f->IsZombie())
            throw std::runtime_error("Cannot open cuts file: " + path);
        TGraph* g = (TGraph*)f->Get(PbPbEvSelKey::kZDCFCalCut);
        if (!g) { f->Close(); throw std::runtime_error(std::string(PbPbEvSelKey::kZDCFCalCut) + " not found in " + path); }
        g_cut1_ = (TGraph*)g->Clone("g_cut1_local");
        f->Close();
        std::cout << "Loaded cut 1 TGraph: " << g_cut1_->GetN() << " points." << std::endl;
    }

    // -------------------------------------------------------------------------
    void FillOneStage(Stage s, const EvData& ev) {
        hh_[s][kZdcFcal] ->Fill(ev.fcal_AC, ev.zdc_tot);
        hh_[s][kZdcTime] ->Fill(ev.tA,      ev.tC);
        hh_[s][kFCalAC]  ->Fill(ev.fcal_A,  ev.fcal_C);
        hh_[s][kPreampAC]->Fill(ev.preamp_A, ev.preamp_C);
        if (ev.ntrk_total > 0) {
            double frac = (double)ev.ntrk_tight / ev.ntrk_total;
            hh_[s][kFrac]    ->Fill(ev.ntrk_total, frac);
            hh_[s][kNtrkFcal]->Fill(ev.fcal_AC, ev.ntrk_tight);
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

        Int_t   b_HLT = 0;
        Float_t FCal_Et_P = 0.f, FCal_Et_N = 0.f;
        Float_t zdc_E[2] = {}, zdc_t[2] = {};
        Float_t preamp[2][4] = {};
        Float_t centrality_f = 0.f;
        std::vector<int>* trk_numqual = nullptr;

        chain.SetBranchStatus("b_HLT_mu4_L1MU3V",          1);
        chain.SetBranchStatus("FCal_Et_P",                   1);
        chain.SetBranchStatus("FCal_Et_N",                   1);
        chain.SetBranchStatus("zdc_ZdcEnergy",               1);
        chain.SetBranchStatus("zdc_ZdcTime",                 1);
        chain.SetBranchStatus("zdc_ZdcModulePreSampleAmp",   1);
        chain.SetBranchStatus("trk_numqual",                 1);
        chain.SetBranchStatus("centrality",                  1);

        chain.SetBranchAddress("b_HLT_mu4_L1MU3V",          &b_HLT);
        chain.SetBranchAddress("FCal_Et_P",                  &FCal_Et_P);
        chain.SetBranchAddress("FCal_Et_N",                  &FCal_Et_N);
        chain.SetBranchAddress("zdc_ZdcEnergy",              zdc_E);
        chain.SetBranchAddress("zdc_ZdcTime",                zdc_t);
        chain.SetBranchAddress("zdc_ZdcModulePreSampleAmp",  preamp);
        chain.SetBranchAddress("trk_numqual",                &trk_numqual);
        chain.SetBranchAddress("centrality",                 &centrality_f);

        const Long64_t n_total = chain.GetEntries();
        std::cout << "FillHists: " << n_total << " total events (PbPb 20" << yr_ << ")" << std::endl;

        for (Long64_t i = 0; i < n_total; ++i) {
            chain.GetEntry(i);
            if (!b_HLT) continue;
            ++n_HLT_;

            EvData ev;
            ev.fcal_A   = FCal_Et_P * 1e-6f;
            ev.fcal_C   = FCal_Et_N * 1e-6f;
            ev.fcal_AC  = ev.fcal_A + ev.fcal_C;
            ev.zdc_tot  = (zdc_E[0] + zdc_E[1]) / 1000.f;
            ev.tA       = zdc_t[0];
            ev.tC       = zdc_t[1];
            ev.preamp_A = preamp[0][0] + preamp[0][1] + preamp[0][2] + preamp[0][3];
            ev.preamp_C = preamp[1][0] + preamp[1][1] + preamp[1][2] + preamp[1][3];
            ev.ntrk_total = 0; ev.ntrk_tight = 0;
            if (trk_numqual && (int)trk_numqual->size() >= 4) {
                ev.ntrk_total = (*trk_numqual)[0];
                ev.ntrk_tight = (*trk_numqual)[3];
            }
            ev.centrality = (int)centrality_f;
            const bool is_ctr80 = (ev.centrality >= 0 && ev.centrality < 80);

            FillOneStage(kNoCut, ev);
            if (is_ctr80) ++n_ctr80_[kNoCut];

            // Cut 1: ZDC banana
            const double cut1_val = PbPbEvSelEvalCut(g_cut1_, ev.fcal_AC);
            const bool   pass_c1  = (cut1_val <= 0. || ev.zdc_tot <= cut1_val);
            FillOneStage(pass_c1 ? kC1Pass : kC1Fail, ev);
            if (is_ctr80) ++n_ctr80_[pass_c1 ? kC1Pass : kC1Fail];
            if (!pass_c1) continue;

            // Cut 2: ZDC time box
            const bool pass_c2 = (std::abs(ev.tA) < CUT2_T_NS && std::abs(ev.tC) < CUT2_T_NS);
            FillOneStage(pass_c2 ? kC2Pass : kC2Fail, ev);
            if (is_ctr80) ++n_ctr80_[pass_c2 ? kC2Pass : kC2Fail];
            if (!pass_c2) continue;

            stored_pre_c3_.push_back(ev);  // cache for cut 3 derivation and beyond
        }
        std::cout << "HLT events: " << n_HLT_ << ", passing cuts 1-2: " << stored_pre_c3_.size() << std::endl;
    }

    // -------------------------------------------------------------------------
    // Project 1D preamp histograms from events entering cut 3 (passing cuts 1+2) and fit EMG.
    void FitPreamp1D() {
        h1_preamp_A_ = (TH1D*)hh_[kC2Pass][kPreampAC]->ProjectionX("h1_preamp_A");
        h1_preamp_C_ = (TH1D*)hh_[kC2Pass][kPreampAC]->ProjectionY("h1_preamp_C");
        h1_preamp_A_->SetDirectory(nullptr);
        h1_preamp_C_->SetDirectory(nullptr);
        h1_preamp_A_->GetYaxis()->SetTitle("Events");
        h1_preamp_C_->GetYaxis()->SetTitle("Events");

        auto doFit = [&](TH1D* h, const char* name) -> TF1* {
            // Pass 1: Gaussian seed
            TF1 g1(Form("__g1_%s", name), "gaus", -700., 350.);
            g1.SetParameters(h->GetMaximum(), -150., 200.);
            h->Fit(&g1, "RQN");
            double mu1  = g1.GetParameter(1);
            double sig1 = std::abs(g1.GetParameter(2));
            // Pass 2: EMG
            TF1* emg = new TF1(name, EMGFunc, mu1 - 3.*sig1, 400., 4);
            emg->SetParameters(g1.GetParameter(0), mu1, sig1, 0.004);
            emg->SetParLimits(2, 5., 500.);
            emg->SetParLimits(3, 1e-5, 0.05);
            h->Fit(emg, "RQN");
            emg->SetLineColor(kRed); emg->SetLineWidth(2);
            return emg;
        };

        fit_preamp_A_ = doFit(h1_preamp_A_, "emg_preamp_A_evsel");
        fit_preamp_C_ = doFit(h1_preamp_C_, "emg_preamp_C_evsel");
        cut3_preamp_A_ = (float)(fit_preamp_A_->GetParameter(1) + 6.*std::abs(fit_preamp_A_->GetParameter(2)));
        cut3_preamp_C_ = (float)(fit_preamp_C_->GetParameter(1) + 6.*std::abs(fit_preamp_C_->GetParameter(2)));
        std::cout << Form("Preamp fit A: mu=%.1f  sig=%.1f  lam=%.5f  cut(mu+6sig)=%.1f\n",
                          fit_preamp_A_->GetParameter(1), std::abs(fit_preamp_A_->GetParameter(2)),
                          std::abs(fit_preamp_A_->GetParameter(3)), (double)cut3_preamp_A_);
        std::cout << Form("Preamp fit C: mu=%.1f  sig=%.1f  lam=%.5f  cut(mu+6sig)=%.1f\n",
                          fit_preamp_C_->GetParameter(1), std::abs(fit_preamp_C_->GetParameter(2)),
                          std::abs(fit_preamp_C_->GetParameter(3)), (double)cut3_preamp_C_);
    }

    // -------------------------------------------------------------------------
    void ApplyCut3() {
        for (const auto& ev : stored_pre_c3_) {
            const bool pass_c3 = (ev.preamp_A < cut3_preamp_A_ && ev.preamp_C < cut3_preamp_C_);
            FillOneStage(pass_c3 ? kC3Pass : kC3Fail, ev);
            if (ev.centrality >= 0 && ev.centrality < 80)
                ++n_ctr80_[pass_c3 ? kC3Pass : kC3Fail];
            if (pass_c3) stored_.push_back(ev);
        }
        std::cout << "Cut 3 applied (preamp_A<" << cut3_preamp_A_
                  << ", preamp_C<" << cut3_preamp_C_ << "): pass=" << stored_.size()
                  << " fail=" << (stored_pre_c3_.size() - stored_.size()) << std::endl;
    }

    // -------------------------------------------------------------------------
    void DeriveCut4() {
        auto fits = FitSlices2D(hh_[kC3Pass][kFrac], N_BINS_FRAC_SL);
        if (fits.empty()) { std::cerr << "DeriveCut4: no converged slices!" << std::endl; return; }
        std::vector<double> gx, gy_lo;
        for (const auto& f : fits) {
            gx.push_back(f.x);
            gy_lo.push_back(f.mu - CUT4_N_SIGMA * f.sig);
            std::cout << Form("  cut4 slice x=%.0f  mu=%.4f  sig=%.4f  cut_lo=%.4f\n",
                              f.x, f.mu, f.sig, f.mu - CUT4_N_SIGMA * f.sig);
        }
        g_cut4_ = new TGraph((int)gx.size(), gx.data(), gy_lo.data());
        g_cut4_->SetTitle(Form("nTrk HItight frac lower cut (mu-%.0fsigma);N_{{trk}}^{{total}};frac cut", CUT4_N_SIGMA));
        std::cout << "Cut 4 derived: " << g_cut4_->GetN() << " slices." << std::endl;
    }

    // -------------------------------------------------------------------------
    void FillCuts45() {
        if (!g_cut4_) { std::cerr << "FillCuts45: cut 4 not derived!" << std::endl; return; }

        // First pass: fill stages kC4Pass/kC4Fail and accumulate nTrk-FCal for cut 5 derivation
        for (const auto& ev : stored_) {
            if (ev.ntrk_total <= 0) continue;
            const double frac  = (double)ev.ntrk_tight / ev.ntrk_total;
            const double cut4  = PbPbEvSelEvalCut(g_cut4_, ev.ntrk_total);
            const bool pass_c4 = (frac >= cut4);
            FillOneStage(pass_c4 ? kC4Pass : kC4Fail, ev);
            if (ev.centrality >= 0 && ev.centrality < 80)
                ++n_ctr80_[pass_c4 ? kC4Pass : kC4Fail];
        }
        std::cout << "Cut 4 filled: pass=" << (Long64_t)hh_[kC4Pass][kZdcFcal]->Integral()
                  << " fail=" << (Long64_t)hh_[kC4Fail][kZdcFcal]->Integral() << std::endl;

        // Derive cut 5 from kC4Pass nTrk-FCal distribution
        auto fits5 = FitSlices2D(hh_[kC4Pass][kNtrkFcal], N_BINS_FCAL_SL);
        if (fits5.empty()) { std::cerr << "DeriveCut5: no converged slices!" << std::endl; return; }
        // Collect all slice points; filter FCal > FCAL_CUT5_LIN_MIN for linear fit
        std::vector<double> gx, gy_lo, gy_hi, lin_x, lin_y_lo, lin_y_hi;
        for (const auto& f : fits5) {
            gx.push_back(f.x);
            gy_lo.push_back(f.mu - CUT5_N_SIGMA * f.sig);
            gy_hi.push_back(f.mu + CUT5_N_SIGMA * f.sig);
            std::cout << Form("  cut5 slice x=%.3f  mu=%.1f  sig=%.1f  lo=%.1f  hi=%.1f\n",
                              f.x, f.mu, f.sig,
                              f.mu - CUT5_N_SIGMA * f.sig, f.mu + CUT5_N_SIGMA * f.sig);
            if (f.x > FCAL_CUT5_LIN_MIN) {
                lin_x.push_back(f.x);
                lin_y_lo.push_back(f.mu - CUT5_N_SIGMA * f.sig);
                lin_y_hi.push_back(f.mu + CUT5_N_SIGMA * f.sig);
            }
        }
        if ((int)lin_x.size() < 2) { std::cerr << "DeriveCut5: too few points for linear fit!" << std::endl; return; }
        // Linear fit to filtered points; build dense TGraph over full FCal range
        TGraph g_lo_raw((int)lin_x.size(), lin_x.data(), lin_y_lo.data());
        TGraph g_hi_raw((int)lin_x.size(), lin_x.data(), lin_y_hi.data());
        TF1 f5_lo("__f5_lo", "pol1", lin_x.front(), lin_x.back());
        TF1 f5_hi("__f5_hi", "pol1", lin_x.front(), lin_x.back());
        g_lo_raw.Fit(&f5_lo, "RNQ");
        g_hi_raw.Fit(&f5_hi, "RNQ");
        std::cout << Form("  cut5 lin-fit lo: intercept=%.1f  slope=%.2f\n",
                          f5_lo.GetParameter(0), f5_lo.GetParameter(1));
        std::cout << Form("  cut5 lin-fit hi: intercept=%.1f  slope=%.2f\n",
                          f5_hi.GetParameter(0), f5_hi.GetParameter(1));
        const double fcal_lo_g = -0.5, fcal_hi_g = 5.5;
        const int n_lin_pts = 12;
        std::vector<double> lg5x(n_lin_pts), lg5y_lo(n_lin_pts), lg5y_hi(n_lin_pts);
        for (int k = 0; k < n_lin_pts; ++k) {
            double x = fcal_lo_g + k * (fcal_hi_g - fcal_lo_g) / (n_lin_pts - 1);
            lg5x[k]    = x;
            lg5y_lo[k] = f5_lo.Eval(x);
            lg5y_hi[k] = f5_hi.Eval(x);
        }
        g_cut5_lo_ = new TGraph(n_lin_pts, lg5x.data(), lg5y_lo.data());
        g_cut5_hi_ = new TGraph(n_lin_pts, lg5x.data(), lg5y_hi.data());
        g_cut5_lo_->SetTitle(Form("nTrk HItight lower cut linear fit (FCal>%.1f TeV);FCal E_{{T}} [TeV];nTrk cut", FCAL_CUT5_LIN_MIN));
        g_cut5_hi_->SetTitle(Form("nTrk HItight upper cut linear fit (FCal>%.1f TeV);FCal E_{{T}} [TeV];nTrk cut", FCAL_CUT5_LIN_MIN));
        std::cout << "Cut 5 linear fit done: " << (int)lin_x.size() << " points used." << std::endl;

        // Second pass: fill stages kC5Pass/kC5Fail
        for (const auto& ev : stored_) {
            if (ev.ntrk_total <= 0) continue;
            const double frac  = (double)ev.ntrk_tight / ev.ntrk_total;
            const double cut4  = PbPbEvSelEvalCut(g_cut4_, ev.ntrk_total);
            if (frac < cut4) continue;  // only events passing cut 4

            const double cut5_lo = PbPbEvSelEvalCut(g_cut5_lo_, ev.fcal_AC);
            const double cut5_hi = PbPbEvSelEvalCut(g_cut5_hi_, ev.fcal_AC);
            const bool pass_c5   = (ev.ntrk_tight >= cut5_lo && ev.ntrk_tight <= cut5_hi);
            FillOneStage(pass_c5 ? kC5Pass : kC5Fail, ev);
            if (ev.centrality >= 0 && ev.centrality < 80)
                ++n_ctr80_[pass_c5 ? kC5Pass : kC5Fail];
        }
        std::cout << "Cut 5 filled: pass=" << (Long64_t)hh_[kC5Pass][kZdcFcal]->Integral()
                  << " fail=" << (Long64_t)hh_[kC5Fail][kZdcFcal]->Integral() << std::endl;
    }

    // -------------------------------------------------------------------------
    void SaveCuts() {
        const std::string path = cuts_dir_ + "/event_sel_cuts_pbpb_20" + yr_ + ".root";
        TFile* f = TFile::Open(path.c_str(), "UPDATE");
        if (!f || f->IsZombie()) { std::cerr << "Cannot open " << path << std::endl; return; }
        if (g_cut4_) {
            TGraph* gc = (TGraph*)g_cut4_->Clone(PbPbEvSelKey::kNTrkFracCutLo);
            gc->Write(PbPbEvSelKey::kNTrkFracCutLo, TObject::kOverwrite);
        }
        if (g_cut5_lo_ && g_cut5_hi_) {
            TGraph* glo = (TGraph*)g_cut5_lo_->Clone(PbPbEvSelKey::kNTrkFCalCutLo);
            TGraph* ghi = (TGraph*)g_cut5_hi_->Clone(PbPbEvSelKey::kNTrkFCalCutHi);
            glo->Write(PbPbEvSelKey::kNTrkFCalCutLo, TObject::kOverwrite);
            ghi->Write(PbPbEvSelKey::kNTrkFCalCutHi, TObject::kOverwrite);
        }
        // Scalar cut values as TParameter<double> for retrieval by downstream code
        TParameter<double>(PbPbEvSelKey::kZDCTimeCutNs,  CUT2_T_NS)      .Write(PbPbEvSelKey::kZDCTimeCutNs,  TObject::kOverwrite);
        TParameter<double>(PbPbEvSelKey::kPreampACutADC, cut3_preamp_A_)  .Write(PbPbEvSelKey::kPreampACutADC, TObject::kOverwrite);
        TParameter<double>(PbPbEvSelKey::kPreampCCutADC, cut3_preamp_C_)  .Write(PbPbEvSelKey::kPreampCCutADC, TObject::kOverwrite);
        TParameter<double>("nTrk_frac_n_sigma",    CUT4_N_SIGMA) .Write("nTrk_frac_n_sigma",    TObject::kOverwrite);
        TParameter<double>("nTrk_FCal_band_n_sigma",CUT5_N_SIGMA).Write("nTrk_FCal_band_n_sigma",TObject::kOverwrite);
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

    // Draw a single pad for the 5-panel canvas
    void DrawPad(int pad_idx, Stage s) {
        gPad->SetLogz();
        gPad->SetRightMargin(0.16);
        gPad->SetLeftMargin(0.13);
        gPad->SetBottomMargin(0.13);
        TH2D* h = hh_[s][pad_idx];
        h->SetContour(99);
        h->Draw("COLZ");
    }

    // 5-panel canvas for a given stage (3 columns × 2 rows; 6th pad left empty)
    void DrawFivePanelCanvas(Stage s, const std::string& outpath) {
        TCanvas* c = new TCanvas(Form("c5_%d", s), "", 1500, 900);
        c->Divide(3, 2, 0.003, 0.003);
        for (int t = 0; t < kNPanelTypes; ++t) {
            c->cd(t + 1);
            DrawPad(t, s);
        }
        c->SaveAs(outpath.c_str());
        delete c;
    }

    // Debug: 1D preamp distributions after cut 2 (entering cut 3), no fitting, full range
    void DrawDebugPreamp() {
        TH2D* h2 = hh_[kC2Pass][kPreampAC];  // events entering cut 3
        TH1D* hA = (TH1D*)h2->ProjectionX("__dbg_preamp_A");
        TH1D* hC = (TH1D*)h2->ProjectionY("__dbg_preamp_C");
        hA->SetDirectory(nullptr); hC->SetDirectory(nullptr);

        TH1D*       hArr[2] = {hA,              hC};
        float       cutV[2] = {cut3_preamp_A_,  cut3_preamp_C_};
        const char* sides[2] = {"A",             "C"};

        TCanvas* c = new TCanvas("c_dbg_preamp", "", 1200, 600);
        c->Divide(2, 1, 0.005, 0.005);
        for (int i = 0; i < 2; ++i) {
            c->cd(i + 1);
            gPad->SetLogy();
            gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);

            hArr[i]->SetMarkerStyle(20); hArr[i]->SetMarkerSize(0.5);
            hArr[i]->GetXaxis()->SetTitle(Form("ZDC preamp sum side %s [ADC counts]", sides[i]));
            hArr[i]->GetYaxis()->SetTitle("Events");
            hArr[i]->Draw("E");

            double ymax = hArr[i]->GetMaximum();
            TLine lcut(cutV[i], hArr[i]->GetMinimum(0.1), cutV[i], ymax * 3.0);
            lcut.SetLineColor(kBlue+1); lcut.SetLineWidth(2); lcut.SetLineStyle(2);
            lcut.DrawClone();

            TLatex tl; tl.SetNDC(); tl.SetTextSize(0.038);
            tl.DrawLatex(0.17, 0.88, Form("Pb+Pb 20%s  side %s (after cut 2)", yr_.c_str(), sides[i]));
            tl.SetTextColor(kBlue+1); tl.SetTextSize(0.033);
            tl.DrawLatex(0.17, 0.82, Form("cut = %.0f ADC", (double)cutV[i]));
            tl.SetTextColor(kBlack);
        }
        c->SaveAs(OutPath("event_sel_cut3_ZDC_preamp_debug").c_str());
        delete c; delete hA; delete hC;
    }

    // Standalone cut 3: 1D preamp (A + C) with EMG fit + μ+6σ cut line — mirrors zoom fit plot
    void DrawStandaloneCut3() {
        if (!h1_preamp_A_ || !h1_preamp_C_ || !fit_preamp_A_ || !fit_preamp_C_) return;
        const int rebin = 3;
        TCanvas* c = new TCanvas("c_cut3_sa", "", 1200, 600);
        c->Divide(2, 1, 0.005, 0.005);

        TH1D*  hArr[2] = {h1_preamp_A_,   h1_preamp_C_};
        TF1*   fArr[2] = {fit_preamp_A_,  fit_preamp_C_};
        float  cutV[2] = {cut3_preamp_A_, cut3_preamp_C_};
        const char* sides[2] = {"A", "C"};
        TH1D* clones[2] = {nullptr, nullptr};

        for (int i = 0; i < 2; ++i) {
            c->cd(i + 1);
            gPad->SetLogy();
            gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);

            TH1D* hd = (TH1D*)hArr[i]->Clone(Form("__hd3_%d", i));
            hd->Rebin(rebin);
            hd->Scale(1.0 / rebin);
            clones[i] = hd;

            hd->GetXaxis()->SetRangeUser(-800., 800.);
            // Find peak max in displayed range for y-axis scaling
            int b1 = hd->FindBin(-800.), b2 = hd->FindBin(800.);
            double lmax = 0.;
            for (int b = b1; b <= b2; ++b) lmax = std::max(lmax, hd->GetBinContent(b));
            hd->GetYaxis()->SetRangeUser(0.5, lmax * 3.0);
            hd->GetXaxis()->SetTitle(Form("ZDC preamp sum side %s [ADC counts]", sides[i]));
            hd->GetYaxis()->SetTitle("Events / 60 ADC");
            hd->SetMarkerStyle(20); hd->SetMarkerSize(0.5);

            hd->Draw("E");
            fArr[i]->Draw("SAME");

            // Vertical cut line
            TLine lcut(cutV[i], 0.5, cutV[i], lmax * 3.0);
            lcut.SetLineColor(kBlue+1); lcut.SetLineWidth(2); lcut.SetLineStyle(2);
            lcut.DrawClone();

            TLatex tl; tl.SetNDC(); tl.SetTextSize(0.038);
            tl.DrawLatex(0.17, 0.88, Form("Pb+Pb 20%s  side %s", yr_.c_str(), sides[i]));
            tl.SetTextSize(0.033);
            tl.DrawLatex(0.17, 0.82, Form("#mu = %.1f ADC",   fArr[i]->GetParameter(1)));
            tl.DrawLatex(0.17, 0.76, Form("#sigma = %.1f ADC", std::abs(fArr[i]->GetParameter(2))));
            tl.DrawLatex(0.17, 0.70, Form("1/#lambda = %.0f ADC", 1./std::abs(fArr[i]->GetParameter(3))));
            tl.SetTextColor(kBlue+1);
            tl.DrawLatex(0.17, 0.64, Form("cut (#mu+6#sigma) = %.0f ADC", (double)cutV[i]));
            tl.SetTextColor(kBlack);

            // Equation box (top right) — matches ZDC_preamp_fit_zoom style
            double mu_f  = fArr[i]->GetParameter(1);
            double sig_f = std::abs(fArr[i]->GetParameter(2));
            double lam_f = std::abs(fArr[i]->GetParameter(3));
            double c1    = mu_f + 0.5*lam_f*sig_f*sig_f;
            double c2    = mu_f + lam_f*sig_f*sig_f;
            double c3    = sig_f * TMath::Sqrt2();
            TPaveText* pt = new TPaveText(0.55, 0.63, 0.97, 0.88, "NDC");
            pt->SetBorderSize(0); pt->SetFillStyle(0);
            pt->SetTextSize(0.031); pt->SetTextAlign(12);
            pt->AddText("f(x) = A e^{#lambda(c_{1}-x)} erfc((c_{2}-x)/c_{3})");
            pt->AddText(Form("c_{1} #equiv #mu+#lambda#sigma^{2}/2 = %.1f ADC", c1));
            pt->AddText(Form("c_{2} #equiv #mu+#lambda#sigma^{2}   = %.1f ADC", c2));
            pt->AddText(Form("c_{3} #equiv #sigma#sqrt{2}          = %.1f ADC", c3));
            pt->AddText(Form("#lambda = %.5f ADC^{-1}", lam_f));
            pt->Draw();
        }
        c->SaveAs(OutPath(Form("event_sel_cut3_%s_standalone", kPbPbEvSelCutLabel(3))).c_str());
        delete c;
        for (TH1D* cl : clones) delete cl;
    }

    // Standalone 2D cut canvas for cut index 1-5
    void DrawStandaloneCut(int cut_idx) {
        if (cut_idx == 3) { DrawStandaloneCut3(); return; }
        // Determine which histogram type and which source stage
        static const HType cut_htype[6] = {
            kZdcFcal,   // unused [0]
            kZdcFcal,   // cut 1
            kZdcTime,   // cut 2
            kPreampAC,  // cut 3 — uses the standalone-only preamp type (index 5)
            kFrac,      // cut 4
            kNtrkFcal   // cut 5
        };
        static const Stage src_stage[6] = {
            kNoCut,     // unused
            kNoCut,     // cut 1: before = no-cut
            kC1Pass,    // cut 2: before = pass cut 1
            kC2Pass,    // cut 3: before = pass cut 2
            kC3Pass,    // cut 4: before = pass cut 3
            kC4Pass     // cut 5: before = pass cut 4
        };
        HType ht = cut_htype[cut_idx];
        Stage ss = src_stage[cut_idx];
        TH2D* hd = (TH2D*)hh_[ss][ht]->Clone(Form("__hsa_%d", cut_idx));
        hd->SetDirectory(nullptr);

        TCanvas* c = new TCanvas(Form("c_cut%d", cut_idx), "", 800, 650);
        c->SetLogz();
        c->SetRightMargin(0.14);
        hd->SetContour(99);
        hd->Draw("COLZ");

        TLatex tl; tl.SetNDC(); tl.SetTextSize(0.036);
        AddLabel(tl, 0.15, 0.88);

        if (cut_idx == 1 && g_cut1_) {
            TGraph* gd = (TGraph*)g_cut1_->Clone("__gc1");
            gd->SetLineColor(kRed); gd->SetLineWidth(2);
            gd->SetMarkerColor(kRed); gd->SetMarkerStyle(20); gd->SetMarkerSize(0.5);
            gd->Draw("LP SAME");
            tl.DrawLatex(0.15, 0.82, "red: ZDC banana cut (#mu+5#sigma)");
        } else if (cut_idx == 2) {
            TLine l; l.SetLineColor(kRed); l.SetLineWidth(2); l.SetLineStyle(2);
            // Draw ±2 ns lines: vertical and horizontal
            l.DrawLine(-CUT2_T_NS, -10., -CUT2_T_NS, 10.);
            l.DrawLine( CUT2_T_NS, -10.,  CUT2_T_NS, 10.);
            l.DrawLine(-10., -CUT2_T_NS, 10., -CUT2_T_NS);
            l.DrawLine(-10.,  CUT2_T_NS, 10.,  CUT2_T_NS);
            tl.DrawLatex(0.15, 0.82, Form("red: |t^{A}| < %.1f ns, |t^{C}| < %.1f ns", CUT2_T_NS, CUT2_T_NS));
        } else if (cut_idx == 3) {
            TLine l; l.SetLineColor(kRed); l.SetLineWidth(2); l.SetLineStyle(2);
            l.DrawLine(cut3_preamp_A_, -1000., cut3_preamp_A_, 3000.);
            l.DrawLine(-1000., cut3_preamp_C_, 3000., cut3_preamp_C_);
            tl.DrawLatex(0.15, 0.82, Form("red: preamp A < %.0f, C < %.0f ADC",
                                          (double)cut3_preamp_A_, (double)cut3_preamp_C_));
        } else if (cut_idx == 4 && g_cut4_) {
            TGraph* gd = (TGraph*)g_cut4_->Clone("__gc4");
            gd->SetLineColor(kRed); gd->SetLineWidth(2);
            gd->SetMarkerColor(kRed); gd->SetMarkerStyle(20); gd->SetMarkerSize(0.5);
            gd->Draw("LP SAME");
            tl.DrawLatex(0.15, 0.82, Form("red: frac lower cut (#mu-%.0f#sigma)", CUT4_N_SIGMA));
        } else if (cut_idx == 5 && g_cut5_lo_ && g_cut5_hi_) {
            TGraph* glo = (TGraph*)g_cut5_lo_->Clone("__gc5lo");
            TGraph* ghi = (TGraph*)g_cut5_hi_->Clone("__gc5hi");
            glo->SetLineColor(kRed); glo->SetLineWidth(2);
            glo->SetMarkerColor(kRed); glo->SetMarkerStyle(20); glo->SetMarkerSize(0.5);
            ghi->SetLineColor(kRed); ghi->SetLineWidth(2);
            ghi->SetMarkerColor(kRed); ghi->SetMarkerStyle(20); ghi->SetMarkerSize(0.5);
            glo->Draw("LP SAME"); ghi->Draw("LP SAME");
            tl.DrawLatex(0.15, 0.82, Form("red: linear fit (#mu#pm%.0f#sigma, FCal>%.1f TeV)", CUT5_N_SIGMA, FCAL_CUT5_LIN_MIN));
        }

        c->SaveAs(OutPath(Form("event_sel_cut%d_%s_standalone", cut_idx, kPbPbEvSelCutLabel(cut_idx))).c_str());
        delete c; delete hd;
    }

    // -------------------------------------------------------------------------
    void DrawSurvivalPlot() {
        // Incremental survival using centrality 0-79 events only.
        // n_before[k] = events entering cut k (centrality 0-79); n_after[k] = events passing.
        static const Stage pass_stages[5] = {kC1Pass, kC2Pass, kC3Pass, kC4Pass, kC5Pass};
        double n_before[5];
        n_before[0] = (double)n_ctr80_[kNoCut];
        for (int k = 1; k < 5; ++k)
            n_before[k] = (double)n_ctr80_[pass_stages[k-1]];

        double xv[5], yv[5], exv[5], eyv[5];
        for (int k = 0; k < 5; ++k) {
            double nk = (double)n_ctr80_[pass_stages[k]];
            double nb = n_before[k];
            double p  = (nb > 0) ? nk / nb : 0.;
            xv[k]  = k + 1;
            yv[k]  = p;
            exv[k] = 0.;
            eyv[k] = (nb > 0) ? std::sqrt(p * (1. - p) / nb) : 0.;
            std::cout << Form("  Cut %d (%s): %.4f survival (%lld / %lld ctr<80)  err=%.2e\n",
                              k+1, kPbPbEvSelCutLabel(k+1), p, (Long64_t)nk, (Long64_t)nb, eyv[k]);
        }

        // Frame histogram: provides x-axis bin labels, no fill/line
        TH1D* hframe = new TH1D("h_surv_frame",
            ";;Survival fraction per cut (centrality 0#minus79)",
            5, 0.5, 5.5);
        hframe->SetDirectory(nullptr);
        for (int k = 0; k < 5; ++k)
            hframe->GetXaxis()->SetBinLabel(k + 1,
                Form("C%d  %s", k+1, kPbPbEvSelCutLabel(k+1)));
        hframe->SetFillStyle(0); hframe->SetLineColor(0);
        hframe->GetYaxis()->SetRangeUser(0.92, 1.002);
        hframe->GetXaxis()->SetLabelSize(0.038);

        TGraphErrors* gr = new TGraphErrors(5, xv, yv, exv, eyv);
        gr->SetMarkerStyle(20); gr->SetMarkerSize(1.2);
        gr->SetMarkerColor(kBlue+1); gr->SetLineColor(kBlue+1); gr->SetLineWidth(1);

        TCanvas* c = new TCanvas("c_survival", "", 800, 600);
        c->SetLeftMargin(0.14); c->SetBottomMargin(0.18);
        hframe->Draw("AXIS");
        gr->Draw("PE SAME");

        TLatex tl; tl.SetNDC(); tl.SetTextSize(0.032);
        AddLabel(tl, 0.63, 0.80);
        tl.SetTextSize(0.027);
        for (int k = 0; k < 5; ++k)
            tl.DrawLatex(0.63, 0.74 - 0.055*k,
                Form("C%d: %.3f%%", k+1, yv[k]*100.));

        c->SaveAs(OutPath("event_sel_survival").c_str());
        delete c; delete hframe; delete gr;
    }

    // -------------------------------------------------------------------------
    void SavePlots() {
        gStyle->SetOptStat(0);
        gStyle->SetPalette(kBird);

        // 5-panel: no cuts
        DrawFivePanelCanvas(kNoCut, OutPath("event_sel_nocuts_5panel"));

        // Per-cut: standalone + 5-panel pass + fail
        static const Stage pass_stages[5] = {kC1Pass, kC2Pass, kC3Pass, kC4Pass, kC5Pass};
        static const Stage fail_stages[5] = {kC1Fail, kC2Fail, kC3Fail, kC4Fail, kC5Fail};
        for (int cut = 1; cut <= 5; ++cut) {
            DrawStandaloneCut(cut);
            DrawFivePanelCanvas(pass_stages[cut-1],
                OutPath(Form("event_sel_cut%d_%s_pass_5panel", cut, kPbPbEvSelCutLabel(cut))));
            DrawFivePanelCanvas(fail_stages[cut-1],
                OutPath(Form("event_sel_cut%d_%s_fail_5panel", cut, kPbPbEvSelCutLabel(cut))));
        }

        DrawDebugPreamp();
        DrawSurvivalPlot();

        std::cout << "All plots saved to: " << out_dir_ << std::endl;
    }
};

// =============================================================================
void plot_pbpb_event_sel_cuts(int run_year = 24) {
    PbPbEventSelCuts plotter(run_year);
    plotter.Run();
}
