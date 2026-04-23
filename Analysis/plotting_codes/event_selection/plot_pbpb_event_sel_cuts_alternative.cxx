// plot_pbpb_event_sel_cuts_alternative.cxx
// Sequential 5-cut event selection pipeline for PbPb data — alternative cut 3.
//   Cut 1: ZDC vs FCal banana (loaded from event_sel_cuts_pbpb_YYYY.root)
//   Cut 2: ZDC time |tA| < 1.8 ns AND |tC| < 1.8 ns
//   Cut 3: ZDC preamp sum: A < 385 ADC AND C < 385 ADC  (fixed hard cut)
//   Cut 4: nTrk HItight fraction vs total nTrk — lower bound at mu-5sigma (per-slice Gauss)
//   Cut 5: nTrk HItight vs FCal ET — band [mu-5sigma, mu+5sigma] (per-slice Gauss)
// Derived cuts 4+5 saved to event_sel_cuts_pbpb_YYYY_alternative.root.
//
// Output per cut: standalone canvas, 5-panel pass canvas, 5-panel fail canvas.
// Additionally: 5-panel no-cuts canvas, per-cut survival rate plot (centrality 0-79).
//
// Usage:
//   .L plot_pbpb_event_sel_cuts_alternative.cxx+
//   plot_pbpb_event_sel_cuts_alternative(24)

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
#include "TBox.h"
#include "TParameter.h"
#include "TMath.h"
#include "../../NTupleProcessingCode/PbPbEventSelConfig.h"

// ---- Cut tuning -------------------------------------------------------------
static const double CUT2_T_NS_ALT    = 1.5;    // ZDC time window [ns]
static const float  CUT3_PREAMP_ALT  = 385.f;  // ZDC preamp upper bound [ADC], both sides
static const double CUT4_N_SIGMA_ALT = 5.0;    // nTrk frac lower cut
static const double CUT5_N_SIGMA_ALT = 5.0;    // nTrk-FCal band half-width

static const int    N_BINS_FRAC_SL_ALT  = 4;
static const int    N_BINS_FCAL_SL_ALT  = 2;
static const int    MIN_ENTRIES_FIT_ALT = 100;
static const double FCAL_CUT5_LIN_MIN_ALT = 0.2;  // FCal lower bound for cut-5 linear fit [TeV]

// ---- Input file map ---------------------------------------------------------
static std::map<int, std::vector<std::string>> BuildFileMapAlt() {
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
enum StageAlt {
    kNoCut_A  = 0,
    kC1Pass_A = 1, kC1Fail_A = 2,
    kC2Pass_A = 3, kC2Fail_A = 4,
    kC3Pass_A = 5, kC3Fail_A = 6,
    kC4Pass_A = 7, kC4Fail_A = 8,
    kC5Pass_A = 9, kC5Fail_A = 10,
    kNStages_A = 11
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
    explicit PbPbEventSelCutsAlt(int run_year) : run_year_(run_year % 2000) {
        yr_      = std::to_string(run_year_);
        out_dir_ = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/"
                   "single_b_analysis/event_selection_alternative_zdc_presample/pbpb_20" + yr_;
        cuts_dir_= "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_20" + yr_;
        auto fm = BuildFileMapAlt();
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
        DeriveCut4();
        FillCuts45();
        SaveCuts();
        SavePlots();
    }

private:
    int run_year_;
    std::string yr_, out_dir_, cuts_dir_;
    std::vector<std::string> infiles_;

    TH2D*   hh_[kNStages_A][kNHTypes_A] = {};
    TGraph* g_cut1_    = nullptr;
    TGraph* g_cut4_    = nullptr;
    TGraph* g_cut5_lo_ = nullptr;
    TGraph* g_cut5_hi_ = nullptr;

    Long64_t n_HLT_             = 0;
    Long64_t n_ctr80_[kNStages_A] = {};
    std::vector<EvDataAlt> stored_;  // events passing cuts 1-3

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
    }

    // -------------------------------------------------------------------------
    void LoadCut1() {
        const std::string path = cuts_dir_ + "/event_sel_cuts_pbpb_20" + yr_ + ".root";
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

        Int_t   b_HLT = 0;
        Float_t FCal_Et_P = 0.f, FCal_Et_N = 0.f;
        Float_t zdc_E[2] = {}, zdc_t[2] = {};
        Float_t preamp[2][4] = {};
        Float_t centrality_f = 0.f;
        std::vector<int>* trk_numqual = nullptr;

        chain.SetBranchStatus("b_HLT_mu4_L1MU3V",         1);
        chain.SetBranchStatus("FCal_Et_P",                  1);
        chain.SetBranchStatus("FCal_Et_N",                  1);
        chain.SetBranchStatus("zdc_ZdcEnergy",              1);
        chain.SetBranchStatus("zdc_ZdcTime",                1);
        chain.SetBranchStatus("zdc_ZdcModulePreSampleAmp",  1);
        chain.SetBranchStatus("trk_numqual",                1);
        chain.SetBranchStatus("centrality",                 1);

        chain.SetBranchAddress("b_HLT_mu4_L1MU3V",         &b_HLT);
        chain.SetBranchAddress("FCal_Et_P",                 &FCal_Et_P);
        chain.SetBranchAddress("FCal_Et_N",                 &FCal_Et_N);
        chain.SetBranchAddress("zdc_ZdcEnergy",             zdc_E);
        chain.SetBranchAddress("zdc_ZdcTime",               zdc_t);
        chain.SetBranchAddress("zdc_ZdcModulePreSampleAmp", preamp);
        chain.SetBranchAddress("trk_numqual",               &trk_numqual);
        chain.SetBranchAddress("centrality",                &centrality_f);

        const Long64_t n_total = chain.GetEntries();
        std::cout << "FillHists: " << n_total << " total events (PbPb 20" << yr_ << ")" << std::endl;

        for (Long64_t i = 0; i < n_total; ++i) {
            chain.GetEntry(i);
            if (!b_HLT) continue;
            ++n_HLT_;

            EvDataAlt ev;
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

            FillOneStage(kNoCut_A, ev);
            if (is_ctr80) ++n_ctr80_[kNoCut_A];

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

            // Cut 3: ZDC preamp — fixed hard cut at 385 ADC for both sides
            const bool pass_c3 = (ev.preamp_A < CUT3_PREAMP_ALT && ev.preamp_C < CUT3_PREAMP_ALT);
            FillOneStage(pass_c3 ? kC3Pass_A : kC3Fail_A, ev);
            if (is_ctr80) ++n_ctr80_[pass_c3 ? kC3Pass_A : kC3Fail_A];
            if (!pass_c3) continue;

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
            if (ev.centrality >= 0 && ev.centrality < 80)
                ++n_ctr80_[pass_c5 ? kC5Pass_A : kC5Fail_A];
        }
        std::cout << "Cut 5 filled: pass=" << (Long64_t)hh_[kC5Pass_A][kZdcFcal_A]->Integral()
                  << " fail=" << (Long64_t)hh_[kC5Fail_A][kZdcFcal_A]->Integral() << std::endl;
    }

    // -------------------------------------------------------------------------
    void SaveCuts() {
        const std::string path = cuts_dir_ + "/event_sel_cuts_pbpb_20" + yr_ + "_alternative.root";
        TFile* f = TFile::Open(path.c_str(), "RECREATE");
        if (!f || f->IsZombie()) { std::cerr << "Cannot open " << path << std::endl; return; }
        if (g_cut4_) {
            TGraph* gc = (TGraph*)g_cut4_->Clone(PbPbEvSelKey::kNTrkFracCutLo);
            gc->Write(PbPbEvSelKey::kNTrkFracCutLo);
        }
        if (g_cut5_lo_ && g_cut5_hi_) {
            TGraph* glo = (TGraph*)g_cut5_lo_->Clone(PbPbEvSelKey::kNTrkFCalCutLo);
            TGraph* ghi = (TGraph*)g_cut5_hi_->Clone(PbPbEvSelKey::kNTrkFCalCutHi);
            glo->Write(PbPbEvSelKey::kNTrkFCalCutLo);
            ghi->Write(PbPbEvSelKey::kNTrkFCalCutHi);
        }
        TParameter<double>(PbPbEvSelKey::kZDCTimeCutNs,  CUT2_T_NS_ALT)  .Write(PbPbEvSelKey::kZDCTimeCutNs);
        TParameter<double>(PbPbEvSelKey::kPreampACutADC, CUT3_PREAMP_ALT).Write(PbPbEvSelKey::kPreampACutADC);
        TParameter<double>(PbPbEvSelKey::kPreampCCutADC, CUT3_PREAMP_ALT).Write(PbPbEvSelKey::kPreampCCutADC);
        TParameter<double>("nTrk_frac_n_sigma",     CUT4_N_SIGMA_ALT).Write("nTrk_frac_n_sigma");
        TParameter<double>("nTrk_FCal_band_n_sigma", CUT5_N_SIGMA_ALT).Write("nTrk_FCal_band_n_sigma");
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

        TCanvas* c = new TCanvas("c_alt_cut3_sa", "", 1200, 600);
        c->Divide(2, 1, 0.005, 0.005);
        for (int i = 0; i < 2; ++i) {
            c->cd(i + 1);
            gPad->SetLogy();
            gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);

            hd[i] = (TH1D*)hArr[i]->Clone(Form("__alt_hd3_%d", i));
            hd[i]->Rebin(rebin);
            hd[i]->Scale(1.0 / rebin);
            hd[i]->GetXaxis()->SetRangeUser(-800., 800.);

            // Find peak max within displayed window for y-axis scaling
            double lmax = 0.;
            const int b_lo = hd[i]->FindBin(-800.);
            const int b_hi = hd[i]->FindBin( 800.);
            for (int b = b_lo; b <= b_hi; ++b)
                lmax = std::max(lmax, hd[i]->GetBinContent(b));
            hd[i]->GetYaxis()->SetRangeUser(0.5, lmax * 3.0);
            hd[i]->GetXaxis()->SetTitle(Form("ZDC preamp sum side %s [ADC counts]", sides[i]));
            hd[i]->GetYaxis()->SetTitle("Events / 60 ADC");
            hd[i]->SetMarkerStyle(20); hd[i]->SetMarkerSize(0.5);
            hd[i]->Draw("E");

            TLine lcut(CUT3_PREAMP_ALT, 0.5, CUT3_PREAMP_ALT, lmax * 3.0);
            lcut.SetLineColor(kBlue+1); lcut.SetLineWidth(2); lcut.SetLineStyle(2);
            lcut.DrawClone();

            TLatex tl; tl.SetNDC(); tl.SetTextSize(0.038);
            tl.DrawLatex(0.17, 0.88, Form("Pb+Pb 20%s  side %s", yr_.c_str(), sides[i]));
            tl.SetTextColor(kBlue+1); tl.SetTextSize(0.033);
            tl.DrawLatex(0.17, 0.82, Form("cut = %.0f ADC", (double)CUT3_PREAMP_ALT));
            tl.SetTextColor(kBlack);
        }
        c->SaveAs(OutPath(Form("event_sel_cut3_%s_standalone", kPbPbEvSelCutLabel(3))).c_str());
        delete c;
        for (int i = 0; i < 2; ++i) delete hd[i];
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

        DrawSurvivalPlot();
        DrawCumulativeSurvivalPlot();

        std::cout << "All plots saved to: " << out_dir_ << std::endl;
    }
};

// =============================================================================
void plot_pbpb_event_sel_cuts_alternative(int run_year = 24) {
    PbPbEventSelCutsAlt plotter(run_year);
    plotter.Run();
}
