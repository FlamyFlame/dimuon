// plot_pbpb_event_sel_event_level.cxx
// Event-level event selection diagnostics for PbPb data.
// Reads raw skim (HeavyIonD3PD tree), filters on HLT_mu4_L1MU3V event-level trigger.
// ZDC vs FCal cut: two-pass per-slice Gaussian fit with TGraph interpolation.
// Derived cut curves are saved to a ROOT file for use in the muon pair analysis.
//
// Usage:
//   .L plot_pbpb_event_sel_event_level.cxx+
//   plot_pbpb_event_sel_event_level(24)   // run_year = 23 or 24

#include <string>
#include <vector>
#include <map>
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
#include "TVectorD.h"

// ---- ZDC fit tuning ---------------------------------------------------------
// N_BINS_PER_SLICE=2 gives 0.10 TeV slices (120 bins × 0.05 TeV/bin).
// With ~6.8M trigger-passing events, 55+/60 slices converge (min ~110 entries),
// giving dense TGraph interpolation points vs the old 0.25 TeV step-function.
static const int    N_BINS_PER_SLICE    = 2;
static const int    MIN_ENTRIES         = 100;
static const double N_SIGMA_CUT         = 5.0;   // main band: cut at mu + N_SIGMA_CUT * sigma
static const double FCAL_BG_FIT_MIN_TEV = 2.0;   // FCal ET above which background band is also fitted
static const double BG_N_SIGMA_INNER    = 2.5;   // background band second-pass window half-width
static const double BG_N_SIGMA_CUT      = 3.0;   // background band: veto at mu - BG_N_SIGMA_CUT * sigma

// Per-year ZDC search window for the out-of-time pileup band (in h_zdc_fcal_ y-axis units).
static std::pair<double,double> GetBgSearchRange(int yr) {
    static const std::map<int, std::pair<double,double>> kBgMap = {
        {23, {220., 350.}},
        {24, {180., 300.}},
        {25, {170., 320.}},
    };
    auto it = kBgMap.find(yr);
    return (it != kBgMap.end()) ? it->second : std::make_pair(180., 300.);
}

// Per-year quadratic initial guess for background-band mu vs FCal ET [TeV].
// Coefficients a,b,c of  mu(x) = a*x^2 + b*x + c  fitted through 3 calibration points:
//   2023: (2.2,252), (2.6,306), (3.4,291)  -> a=-128.125,  b=750,       c=-777.875
//   2024: (2.0,272), (3.6,256), (5.4,225)  -> a=-325/153,  b=290/153,   c=272+720/153
//   2025: (2.2,270), (3.4,252), (4.6,218)  -> a=-50/9,     b=145/9,     c=270-77/9
static double GetBgMuGuess(int yr, double xcen) {
    struct Quad { double a, b, c; };
    static const std::map<int, Quad> kQuad = {
        {23, {-8./3.,          16./3.,          296.              }},  // pts: (3.5,282),(4.5,266),(5.0,256)
        {24, {-325./153.,     290./153.,       272. + 720./153. }},
        {25, {-50./9.,        145./9.,         270. -  77./9.   }},
    };
    auto it = kQuad.find(yr);
    if (it == kQuad.end()) return 250.;
    const auto& q = it->second;
    return q.a * xcen * xcen + q.b * xcen + q.c;
}

// ---- per-year input file lists ----------------------------------------------
static std::map<int, std::vector<std::string>> BuildFileMap() {
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

// =============================================================================
class PbPbEventSel {
public:
    explicit PbPbEventSel(int run_year) : run_year_(run_year % 2000) {
        yr_           = std::to_string(run_year_);
        out_dir_      = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/"
                        "single_b_analysis/event_selection/pbpb_20" + yr_;
        alt_out_dir_  = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/"
                        "single_b_analysis/event_selection_alternative_banana/pbpb_20" + yr_;
        cuts_dir_     = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_20" + yr_;
        auto file_map = BuildFileMap();
        if (file_map.find(run_year_) == file_map.end())
            throw std::runtime_error("No input files configured for run year " + yr_);
        infiles_ = file_map.at(run_year_);
    }

    void Run() {
        gStyle->SetOptStat(0);
        gStyle->SetPalette(kBird);
        gSystem->mkdir(out_dir_.c_str(), true);
        gSystem->mkdir(alt_out_dir_.c_str(), true);
        BookHists();
        FillHists();
        SavePlots();
    }

    // Run only the centrality ZDC time plot (loads banana cut from saved ROOT file).
    void RunCentralityOnly() {
        gStyle->SetOptStat(0);
        gStyle->SetPalette(kBird);
        gSystem->mkdir(out_dir_.c_str(), true);
        MakeZDCTimeCentralityPlot();
    }

private:
    int run_year_;
    std::string yr_, out_dir_, alt_out_dir_, cuts_dir_;
    std::vector<std::string> infiles_;

    TH2D *h_fcal_ac_    = nullptr;
    TH2D *h_zdc_fcal_   = nullptr;
    TH2D *h_ntrk_fcal_  = nullptr;
    TH2D *h_zdc_t_      = nullptr;
    TH2D *h_frac_       = nullptr;

    // -------------------------------------------------------------------------
    void BookHists() {
        h_fcal_ac_ = new TH2D("h_fcal_ac",
            ";FCal E_{T}^{A} [TeV];FCal E_{T}^{C} [TeV]",
            110, -0.5, 3.0, 110, -0.5, 3.0);

        h_zdc_fcal_ = new TH2D("h_zdc_fcal",
            ";FCal E_{T}^{A+C} [TeV];ZDC E_{total} [TeV]",
            120, -0.5, 5.5, 100, 0., 400.);

        h_ntrk_fcal_ = new TH2D("h_ntrk_fcal",
            ";FCal E_{T}^{A+C} [TeV];N_{trk}^{HItight}",
            120, -0.5, 5.5, 100, 0., 5000.);

        h_zdc_t_ = new TH2D("h_zdc_t",
            ";ZDC t^{A} [ns];ZDC t^{C} [ns]",
            100, -10., 10., 100, -10., 10.);

        h_frac_ = new TH2D("h_frac",
            ";N_{trk}^{total} (p_{T} > 400 MeV);N_{trk}^{HItight} / N_{trk}^{total}",
            100, 0., 5000., 110, 0., 1.1);

        for (TH2D* h : {h_fcal_ac_, h_zdc_fcal_, h_ntrk_fcal_, h_zdc_t_, h_frac_})
            h->SetDirectory(nullptr);
    }

    // -------------------------------------------------------------------------
    void FillHists() {
        TChain chain("HeavyIonD3PD", "HeavyIonD3PD");
        for (const auto& f : infiles_) {
            if (gSystem->AccessPathName(f.c_str())) {
                std::cerr << "File not found (skipping): " << f << std::endl;
                continue;
            }
            chain.Add(f.c_str());
        }
        if (chain.GetEntries() == 0)
            throw std::runtime_error("TChain is empty for run year 20" + yr_);

        chain.SetMakeClass(1);
        chain.SetBranchStatus("*", 0);

        Int_t   b_HLT_mu4 = 0;
        Float_t FCal_Et_P = 0.f, FCal_Et_N = 0.f;
        Float_t zdc_E[2]  = {}, zdc_t[2] = {};
        std::vector<int>* trk_numqual = nullptr;

        chain.SetBranchStatus("b_HLT_mu4_L1MU3V", 1);
        chain.SetBranchStatus("FCal_Et_P",         1);
        chain.SetBranchStatus("FCal_Et_N",         1);
        chain.SetBranchStatus("zdc_ZdcEnergy",     1);
        chain.SetBranchStatus("zdc_ZdcTime",       1);
        chain.SetBranchStatus("trk_numqual",       1);

        chain.SetBranchAddress("b_HLT_mu4_L1MU3V", &b_HLT_mu4);
        chain.SetBranchAddress("FCal_Et_P",         &FCal_Et_P);
        chain.SetBranchAddress("FCal_Et_N",         &FCal_Et_N);
        chain.SetBranchAddress("zdc_ZdcEnergy",     zdc_E);
        chain.SetBranchAddress("zdc_ZdcTime",       zdc_t);
        chain.SetBranchAddress("trk_numqual",       &trk_numqual);

        const Long64_t n = chain.GetEntries();
        std::cout << "Processing " << n << " events (Pb+Pb 20" << yr_ << ")" << std::endl;
        Long64_t n_pass = 0;

        for (Long64_t i = 0; i < n; ++i) {
            chain.GetEntry(i);
            if (!b_HLT_mu4) continue;
            ++n_pass;

            const float fcal_A  = FCal_Et_P * 1e-6f;       // MeV → TeV; [0]=A-side
            const float fcal_C  = FCal_Et_N * 1e-6f;
            const float fcal_AC = fcal_A + fcal_C;
            const float zdc_tot = (zdc_E[0] + zdc_E[1]) / 1000.f;  // GeV → TeV
            const float zdc_tA  = zdc_t[1];                // [1]=A-side
            const float zdc_tC  = zdc_t[0];                // [0]=C-side

            // trk_numqual[0] = all tracks pT > 400 MeV
            // trk_numqual[3] = HItight + pT > 400 MeV
            int ntrk_total = 0, ntrk_tight = 0;
            if (trk_numqual && (int)trk_numqual->size() >= 4) {
                ntrk_total = (*trk_numqual)[0];
                ntrk_tight = (*trk_numqual)[3];
            }

            h_fcal_ac_ ->Fill(fcal_A, fcal_C);
            h_zdc_fcal_->Fill(fcal_AC, zdc_tot);
            h_ntrk_fcal_->Fill(fcal_AC, ntrk_tight);
            h_zdc_t_   ->Fill(zdc_tA, zdc_tC);
            if (ntrk_total > 0)
                h_frac_->Fill(ntrk_total, (double)ntrk_tight / ntrk_total);
        }
        std::cout << "Trigger-passing events: " << n_pass << " / " << n << std::endl;
    }

    // -------------------------------------------------------------------------
    std::string OutPath(const std::string& stem) const {
        return out_dir_ + "/" + stem + "_pbpb_20" + yr_ + ".png";
    }
    void AddLabel(TLatex& tl, double x, double y) const {
        tl.DrawLatex(x, y, ("Pb+Pb 20" + yr_ + " data").c_str());
    }
    void DrawColz(TH2D* h) { h->SetContour(99); h->Draw("COLZ"); }

    // -------------------------------------------------------------------------
    void SavePlots() {
        TLatex tl; tl.SetNDC(); tl.SetTextSize(0.038);

        {
            TCanvas* c = new TCanvas("c_fcal_ac", "", 700, 600);
            c->SetLogz(); c->SetRightMargin(0.13);
            DrawColz(h_fcal_ac_);
            AddLabel(tl, 0.15, 0.87);
            c->SaveAs(OutPath("FCal_Et_A_vs_C").c_str());
            delete c;
        }
        {
            TH2D* hd = (TH2D*)h_zdc_fcal_->Clone("__h_zdc_plain");
            hd->SetDirectory(nullptr);
            TCanvas* c = new TCanvas("c_zdc_fcal", "", 700, 600);
            c->SetLogz(); c->SetRightMargin(0.13);
            DrawColz(hd);
            AddLabel(tl, 0.15, 0.87);
            c->SaveAs(OutPath("ZDC_E_tot_vs_FCal_Et_AC").c_str());
            delete c; delete hd;
        }
        {
            TCanvas* c = new TCanvas("c_ntrk_fcal", "", 700, 600);
            c->SetLogz(); c->SetRightMargin(0.13);
            DrawColz(h_ntrk_fcal_);
            AddLabel(tl, 0.15, 0.87);
            c->SaveAs(OutPath("ntrk_HItight_vs_FCal_Et_AC").c_str());
            delete c;
        }
        {
            TCanvas* c = new TCanvas("c_zdc_t", "", 700, 600);
            c->SetLogz(); c->SetRightMargin(0.13);
            DrawColz(h_zdc_t_);
            AddLabel(tl, 0.15, 0.87);
            c->SaveAs(OutPath("ZDC_t_A_vs_C").c_str());
            delete c;
        }
        {
            TCanvas* c = new TCanvas("c_frac", "", 700, 600);
            c->SetLogz(); c->SetRightMargin(0.13);
            DrawColz(h_frac_);
            AddLabel(tl, 0.15, 0.87);
            c->SaveAs(OutPath("ntrk_HItight_frac_vs_total").c_str());
            delete c;
        }
        FitAndPlotZDC();
        FitAndPlotZDCAlt();
        MakeZDCTimeCentralityPlot();
    }

    // -------------------------------------------------------------------------
    // Evaluate the cut TGraph with clamped constant extrapolation outside the
    // fitted FCal range (avoids unphysical linear extrapolation into tails).
    static double EvalCut(const TGraph* g, double x) {
        if (g->GetN() == 0) return -1.;
        if (x <= g->GetX()[0])           return g->GetY()[0];
        if (x >= g->GetX()[g->GetN()-1]) return g->GetY()[g->GetN()-1];
        return g->Eval(x);  // linear interpolation between fitted points
    }

    // -------------------------------------------------------------------------
    void FitAndPlotZDC() {
        const int nx       = h_zdc_fcal_->GetNbinsX();
        const int ny       = h_zdc_fcal_->GetNbinsY();
        const int n_slices = nx / N_BINS_PER_SLICE;

        const auto [bg_lo, bg_hi] = GetBgSearchRange(run_year_);

        // Main band fit results (all FCal slices that converge)
        std::vector<double> gr_x, gr_cut, gr_main_cut, gr_mu, gr_sigma, gr_ex;
        // Background band fit results (FCal > FCAL_BG_FIT_MIN_TEV slices only)
        std::vector<double> gr_bg_x, gr_bg_cut, gr_bg_mu, gr_bg_sig;
        // Initial guess used for each bg-region slice (for debug plot)
        std::vector<double> gr_guess_x, gr_guess_mu;
        // Slices in bg region where Gaussian didn't converge — extended by post-loop retry/extrapolation
        struct PendingBg { double xcen; int gr_idx; int isl; };
        std::vector<PendingBg> pending_bg;

        struct FirstSlice {
            bool   found = false;
            int    ix_lo = 0, ix_hi = 0;
            double xlo = 0., xhi = 0.;
            double amp2 = 0., mu2 = 0., sig2 = 0., lo2 = 0., hi2 = 0.;
        } first_sl;

        TF1 gfunc("gfit_evsel", "gaus");

        for (int isl = 0; isl < n_slices; ++isl) {
            const int    ix_lo = isl * N_BINS_PER_SLICE + 1;
            const int    ix_hi = std::min((isl + 1) * N_BINS_PER_SLICE, nx);
            const double xlo   = h_zdc_fcal_->GetXaxis()->GetBinLowEdge(ix_lo);
            const double xhi   = h_zdc_fcal_->GetXaxis()->GetBinUpEdge(ix_hi);
            const double xcen  = 0.5 * (xlo + xhi);

            TH1D* hpy = (TH1D*)h_zdc_fcal_->ProjectionY(
                Form("__hpy_ev_%d", isl), ix_lo, ix_hi);
            hpy->SetDirectory(nullptr);

            if (hpy->GetEntries() < MIN_ENTRIES) { delete hpy; continue; }

            const double mu0  = hpy->GetBinCenter(hpy->GetMaximumBin());
            const double rms  = (hpy->GetRMS() > 0.) ? hpy->GetRMS() : 1.0;
            const double ymin = hpy->GetXaxis()->GetXmin();
            const double ymax = hpy->GetXaxis()->GetXmax();

            // ---- Main banana band: two-pass Gaussian fit --------------------
            gfunc.SetRange(std::max(ymin, mu0 - 3.*rms), std::min(ymax, mu0 + 3.*rms));
            gfunc.SetParameters(hpy->GetMaximum(), mu0, rms);
            if (hpy->Fit(&gfunc, "RNQ") != 0 || gfunc.GetParameter(2) <= 0.) {
                delete hpy; continue;
            }
            const double mu1 = gfunc.GetParameter(1), sig1 = gfunc.GetParameter(2);
            const double lo2 = std::max(ymin, mu1 - 3.*sig1);
            const double hi2 = std::min(ymax, mu1 + 3.*sig1);
            gfunc.SetRange(lo2, hi2);
            gfunc.SetParameters(gfunc.GetParameter(0), mu1, sig1);
            if (hpy->Fit(&gfunc, "RNQ") != 0 || gfunc.GetParameter(2) <= 0.) {
                delete hpy; continue;
            }
            const double mu2       = gfunc.GetParameter(1);
            const double sig2      = gfunc.GetParameter(2);
            const double main_cut  = mu2 + N_SIGMA_CUT * sig2;

            if (!first_sl.found)
                first_sl = {true, ix_lo, ix_hi, xlo, xhi,
                            gfunc.GetParameter(0), mu2, sig2, lo2, hi2};

            // ---- Background (pileup) band: two-pass fit for FCal > 2 TeV ---
            double final_cut = main_cut;
            bool   bg_ok     = false;
            double bg_mu2 = 0., bg_sig2 = 0.;
            if (xcen > FCAL_BG_FIT_MIN_TEV) {
                // Quadratic initial guess for bg_mu (per-year polynomial through 3 calibration points)
                const double bg_mu0_raw = GetBgMuGuess(run_year_, xcen);
                const double bg_mu0    = std::max(bg_lo, std::min(bg_hi, bg_mu0_raw));
                gr_guess_x.push_back(xcen);
                gr_guess_mu.push_back(bg_mu0);
                const double sigma_init = (bg_hi - bg_lo) / 6.;
                // Pass-1 fit range: ±3σ_init centered on quadratic guess, clamped to [bg_lo,bg_hi].
                // This excludes the main-band tail at the lower edge of the search window,
                // which otherwise dominates and causes fit divergence at low FCal ET.
                const double p1_lo = std::max(bg_lo, bg_mu0 - 3.*sigma_init);
                const double p1_hi = std::min(bg_hi, bg_mu0 + 3.*sigma_init);
                const int bp1_lo = hpy->FindBin(p1_lo);
                const int bp1_hi = hpy->FindBin(p1_hi);
                double bg_amp0 = 1.;
                for (int b = bp1_lo; b <= bp1_hi; ++b)
                    bg_amp0 = std::max(bg_amp0, hpy->GetBinContent(b));
                const double bg_entries = hpy->Integral(bp1_lo, bp1_hi);
                if (bg_entries >= MIN_ENTRIES) {
                    // Pass 1: narrow window centered on quadratic guess
                    gfunc.SetRange(p1_lo, p1_hi);
                    gfunc.SetParameters(bg_amp0, bg_mu0, sigma_init);
                    if (hpy->Fit(&gfunc, "RNQ") == 0 && gfunc.GetParameter(2) > 0.) {
                        const double bg_mu1  = gfunc.GetParameter(1);
                        const double bg_sig1 = gfunc.GetParameter(2);
                        // Pass 2: mu1 ± BG_N_SIGMA_INNER * sigma1, clamped to [bg_lo, bg_hi]
                        gfunc.SetRange(std::max(bg_lo, bg_mu1 - BG_N_SIGMA_INNER*bg_sig1),
                                       std::min(bg_hi, bg_mu1 + BG_N_SIGMA_INNER*bg_sig1));
                        gfunc.SetParameters(gfunc.GetParameter(0), bg_mu1, bg_sig1);
                        if (hpy->Fit(&gfunc, "RNQ") == 0 && gfunc.GetParameter(2) > 0.) {
                            bg_mu2  = gfunc.GetParameter(1);
                            bg_sig2 = gfunc.GetParameter(2);
                            // Sanity: mu2 must remain within the search window
                            if (bg_mu2 >= bg_lo && bg_mu2 <= bg_hi) {
                                bg_ok   = true;
                                const double bg_cut = bg_mu2 - BG_N_SIGMA_CUT * bg_sig2;
                                gr_bg_x.push_back(xcen);
                                gr_bg_cut.push_back(bg_cut);
                                gr_bg_mu.push_back(bg_mu2);
                                gr_bg_sig.push_back(bg_sig2);
                                final_cut = std::max(main_cut, bg_cut);
                            }
                        }
                    }
                }
            }
            delete hpy;

            std::cout << Form("  sl %2d [%.2f-%.2f]  mu=%.2f sig=%.2f main=%.2f",
                              isl, xlo, xhi, mu2, sig2, main_cut);
            if (bg_ok)
                std::cout << Form("  bg_mu=%.2f bg_sig=%.2f bg_cut=%.2f  final=%.2f",
                                  bg_mu2, bg_sig2, bg_mu2 - BG_N_SIGMA_CUT*bg_sig2, final_cut);
            std::cout << "\n";

            if (xcen > FCAL_BG_FIT_MIN_TEV && !bg_ok)
                pending_bg.push_back({xcen, (int)gr_x.size(), isl});
            gr_x.push_back(xcen);
            gr_cut.push_back(final_cut);
            gr_main_cut.push_back(main_cut);
            gr_mu.push_back(mu2);
            gr_sigma.push_back(sig2);
            gr_ex.push_back(0.5*(xhi - xlo));
        }

        const int n_fit = (int)gr_x.size();
        const int n_bg  = (int)gr_bg_x.size();

        // ---- Post-loop: wider-slice retry then TGraph extrapolation for non-converged bg slices ----
        // Results extend the plotting TGraphs to FCAL_BG_FIT_MIN_TEV; gr_cut updated where tighter.
        std::vector<double> ext_bg_x, ext_bg_cut, ext_bg_mu, ext_bg_sig;
        int n_bg_interp = 0;

        if (!pending_bg.empty() && n_bg >= 2) {
            TF1 gfe("gfit_ext", "gaus");
            const double sigma_init = (bg_hi - bg_lo) / 6.;
            // Two-pass bg Gaussian fit on arbitrary projection h at xcen; returns true if converged.
            auto try_bg_fit = [&](TH1D* h, double xcen, double& mu_out, double& sig_out) -> bool {
                const double mu0   = std::max(bg_lo, std::min(bg_hi, GetBgMuGuess(run_year_, xcen)));
                const double p1_lo = std::max(bg_lo, mu0 - 3.*sigma_init);
                const double p1_hi = std::min(bg_hi, mu0 + 3.*sigma_init);
                const int bp1 = h->FindBin(p1_lo), bp2 = h->FindBin(p1_hi);
                if (h->Integral(bp1, bp2) < MIN_ENTRIES) return false;
                double amp = 1.;
                for (int b = bp1; b <= bp2; ++b) amp = std::max(amp, h->GetBinContent(b));
                gfe.SetRange(p1_lo, p1_hi);
                gfe.SetParameters(amp, mu0, sigma_init);
                if (h->Fit(&gfe, "RNQ") != 0 || gfe.GetParameter(2) <= 0.) return false;
                const double mu1 = gfe.GetParameter(1), sig1 = gfe.GetParameter(2);
                gfe.SetRange(std::max(bg_lo, mu1 - BG_N_SIGMA_INNER*sig1),
                             std::min(bg_hi, mu1 + BG_N_SIGMA_INNER*sig1));
                gfe.SetParameters(gfe.GetParameter(0), mu1, sig1);
                if (h->Fit(&gfe, "RNQ") != 0 || gfe.GetParameter(2) <= 0.) return false;
                mu_out  = gfe.GetParameter(1);
                sig_out = gfe.GetParameter(2);
                return (mu_out >= bg_lo && mu_out <= bg_hi);
            };

            // Build extrapolation TGraphs from converged points (clamped at boundaries)
            TGraph g_mu_e (n_bg, gr_bg_x.data(), gr_bg_mu.data());
            TGraph g_sig_e(n_bg, gr_bg_x.data(), gr_bg_sig.data());
            auto eval_e = [](const TGraph& g, double x) -> double {
                int n = g.GetN();
                if (x <= g.GetX()[0])   return g.GetY()[0];
                if (x >= g.GetX()[n-1]) return g.GetY()[n-1];
                return g.Eval(x);
            };

            for (auto& pb : pending_bg) {
                double mu_e = 0., sig_e = 0.;
                bool found = false;
                // Try 4× wider then 8× wider bins (0.2 and 0.4 TeV slices)
                for (int mult : {4, 8}) {
                    const int half    = mult / 2;
                    const int ix_w_lo = std::max(1,  (pb.isl - half) * N_BINS_PER_SLICE + 1);
                    const int ix_w_hi = std::min(nx, (pb.isl - half + mult) * N_BINS_PER_SLICE);
                    TH1D* hw = (TH1D*)h_zdc_fcal_->ProjectionY(
                        Form("__hpy_w%d_%d", mult, pb.isl), ix_w_lo, ix_w_hi);
                    hw->SetDirectory(nullptr);
                    found = try_bg_fit(hw, pb.xcen, mu_e, sig_e);
                    delete hw;
                    if (found) break;
                }
                // Fallback: extrapolate from converged TGraph
                if (!found) {
                    mu_e  = eval_e(g_mu_e,  pb.xcen);
                    sig_e = eval_e(g_sig_e, pb.xcen);
                }
                const double cut_e = mu_e - BG_N_SIGMA_CUT * sig_e;
                ext_bg_x.push_back(pb.xcen);
                ext_bg_mu.push_back(mu_e);
                ext_bg_sig.push_back(sig_e);
                ext_bg_cut.push_back(cut_e);

                const double new_cut = std::max(gr_main_cut[pb.gr_idx], cut_e);
                if (new_cut > gr_cut[pb.gr_idx]) {
                    gr_cut[pb.gr_idx] = new_cut;
                    ++n_bg_interp;
                }
            }
        }

        // Full bg vectors: extended (lower FCal) prepended to converged (higher FCal)
        std::vector<double> gr_bg_x_all = ext_bg_x,   gr_bg_cut_all = ext_bg_cut,
                            gr_bg_mu_all = ext_bg_mu, gr_bg_sig_all = ext_bg_sig;
        gr_bg_x_all.insert(gr_bg_x_all.end(),    gr_bg_x.begin(),   gr_bg_x.end());
        gr_bg_cut_all.insert(gr_bg_cut_all.end(), gr_bg_cut.begin(), gr_bg_cut.end());
        gr_bg_mu_all.insert(gr_bg_mu_all.end(),   gr_bg_mu.begin(),  gr_bg_mu.end());
        gr_bg_sig_all.insert(gr_bg_sig_all.end(), gr_bg_sig.begin(), gr_bg_sig.end());
        const int n_bg_all = (int)gr_bg_x_all.size();

        std::cout << n_fit << " main-band slices converged; "
                  << n_bg  << " bg-band converged; "
                  << (int)ext_bg_x.size() << " extended (wider-slice or extrapolated); "
                  << n_bg_interp << " tightened gr_cut." << std::endl;
        if (n_fit < 2) {
            std::cerr << "FitAndPlotZDC: too few converged slices!" << std::endl;
            return;
        }

        TGraph  g_cut     (n_fit, gr_x.data(),       gr_cut.data());
        TGraph  g_main    (n_fit, gr_x.data(),       gr_main_cut.data());
        TGraph  g_mu      (n_fit, gr_x.data(),       gr_mu.data());
        TGraphErrors g_mu_err(n_fit, gr_x.data(), gr_mu.data(),
                              gr_ex.data(), gr_sigma.data());
        // Background cut graph and mu±sigma graph (converged + extended)
        const bool have_bg = (n_bg_all > 0);
        TGraph       g_bg;
        TGraphErrors g_bg_mu_err;
        if (have_bg) {
            g_bg        = TGraph(n_bg_all, gr_bg_x_all.data(), gr_bg_cut_all.data());
            g_bg_mu_err = TGraphErrors(n_bg_all, gr_bg_x_all.data(), gr_bg_mu_all.data(),
                                       nullptr, gr_bg_sig_all.data());
        }

        // ---- save cut curves to ROOT file ------------------------------------
        const std::string cuts_path = cuts_dir_ + "/event_sel_cuts_pbpb_20" + yr_ + ".root";
        {
            TFile* fcuts = TFile::Open(cuts_path.c_str(), "UPDATE");
            if (!fcuts || fcuts->IsZombie()) {
                std::cerr << "Cannot open cuts file for writing: " << cuts_path << std::endl;
            } else {
                auto wg = [&](TGraph* g, const char* name, const char* title) {
                    g->SetTitle(title); g->Write(name, TObject::kOverwrite);
                };
                wg((TGraph*)g_cut.Clone ("g_ZDC_FCal_cut"),
                   "g_ZDC_FCal_cut",
                   Form("ZDC cut: max(main+%.0f#sigma, bg-%.0f#sigma);FCal E_{T} [TeV];ZDC cut",
                        N_SIGMA_CUT, BG_N_SIGMA_CUT));
                wg((TGraph*)g_main.Clone("g_ZDC_FCal_main_cut"),
                   "g_ZDC_FCal_main_cut",
                   Form("ZDC banana main band #mu+%.0f#sigma;FCal E_{T} [TeV];ZDC cut", N_SIGMA_CUT));
                wg((TGraph*)g_mu.Clone  ("g_ZDC_FCal_mu"),
                   "g_ZDC_FCal_mu",
                   "ZDC Gaussian mean per FCal slice;FCal E_{T} [TeV];ZDC mean");
                ((TGraphErrors*)g_mu_err.Clone("g_ZDC_FCal_mu_sigma"))
                    ->Write("g_ZDC_FCal_mu_sigma", TObject::kOverwrite);
                if (have_bg)
                    wg((TGraph*)g_bg.Clone("g_ZDC_FCal_bg_cut"),
                       "g_ZDC_FCal_bg_cut",
                       Form("ZDC bg band #mu-%.0f#sigma;FCal E_{T} [TeV];ZDC cut", BG_N_SIGMA_CUT));
                fcuts->Close();
                std::cout << "Cut TGraphs saved to: " << cuts_path << std::endl;
            }
        }

        const std::string lbl = "Pb+Pb 20" + yr_ + " data";
        TLatex tl; tl.SetNDC(); tl.SetTextSize(0.036);

        // ---- Canvas 1: 2D + final cut curve ---------------------------------
        {
            TH2D* hd = (TH2D*)h_zdc_fcal_->Clone("__h2_c1_ev");
            hd->SetDirectory(nullptr);
            TCanvas* c1 = new TCanvas("c_zdc_fcal_cut_ev", "", 800, 650);
            c1->SetLogz(); c1->SetRightMargin(0.13);
            hd->SetContour(99); hd->Draw("COLZ");
            {
                TGraph* gd = (TGraph*)g_cut.Clone("__gd_c1");
                gd->SetMarkerStyle(20); gd->SetMarkerSize(0.6);
                gd->SetMarkerColor(kRed); gd->SetLineColor(kRed); gd->SetLineWidth(2);
                gd->Draw("LP SAME");
            }
            tl.DrawLatex(0.15, 0.87, lbl.c_str());
            tl.DrawLatex(0.15, 0.81,
                Form("red: max(main #mu+%.0f#sigma, bg #mu-%.0f#sigma)", N_SIGMA_CUT, BG_N_SIGMA_CUT));
            c1->SaveAs(OutPath("ZDC_E_tot_vs_FCal_Et_AC_with_cut").c_str());
            delete c1; delete hd;
        }

        // ---- Canvas 2: pass | fail using TGraph interpolation ---------------
        {
            TH2D* h_pass = (TH2D*)h_zdc_fcal_->Clone("__h2_pass_ev");
            TH2D* h_fail = (TH2D*)h_zdc_fcal_->Clone("__h2_fail_ev");
            h_pass->Reset(); h_pass->SetDirectory(nullptr);
            h_fail->Reset(); h_fail->SetDirectory(nullptr);

            for (int ix = 1; ix <= nx; ++ix) {
                const double xc  = h_zdc_fcal_->GetXaxis()->GetBinCenter(ix);
                const double cut = EvalCut(&g_cut, xc);
                for (int iy = 1; iy <= ny; ++iy) {
                    const double w = h_zdc_fcal_->GetBinContent(ix, iy);
                    if (w == 0.) continue;
                    const double yc = h_zdc_fcal_->GetYaxis()->GetBinCenter(iy);
                    if (cut > 0. && yc > cut) h_fail->Fill(xc, yc, w);
                    else                       h_pass->Fill(xc, yc, w);
                }
            }

            TCanvas* c2 = new TCanvas("c_zdc_fcal_pf_ev", "", 1400, 620);
            c2->Divide(2, 1);
            auto draw_cut = [&](TGraph* src, const char* name) {
                TGraph* gd = (TGraph*)src->Clone(name);
                gd->SetMarkerStyle(20); gd->SetMarkerSize(0.6);
                gd->SetMarkerColor(kRed); gd->SetLineColor(kRed); gd->SetLineWidth(2);
                gd->Draw("LP SAME");
            };
            c2->cd(1); gPad->SetLogz(); gPad->SetRightMargin(0.15);
            h_pass->SetContour(99); h_pass->Draw("COLZ");
            draw_cut(&g_cut, "__gd_pass");
            { TLatex t2; t2.SetNDC(); t2.SetTextSize(0.036);
              t2.DrawLatex(0.15, 0.87, lbl.c_str());
              t2.DrawLatex(0.15, 0.81, "Passing cut"); }
            c2->cd(2); gPad->SetLogz(); gPad->SetRightMargin(0.15);
            h_fail->SetContour(99); h_fail->Draw("COLZ");
            draw_cut(&g_cut, "__gd_fail");
            { TLatex t2; t2.SetNDC(); t2.SetTextSize(0.036);
              t2.DrawLatex(0.15, 0.87, lbl.c_str());
              t2.DrawLatex(0.15, 0.81, "Failing cut"); }
            c2->SaveAs(OutPath("ZDC_E_tot_vs_FCal_Et_AC_pass_fail").c_str());
            delete c2; delete h_pass; delete h_fail;
        }

        // ---- Canvas 3 (debug): 2D + cut + mu±sigma --------------------------
        {
            TH2D* hd = (TH2D*)h_zdc_fcal_->Clone("__h2_c3_ev");
            hd->SetDirectory(nullptr);
            TCanvas* c3 = new TCanvas("c_zdc_fcal_dbg_ev", "", 800, 650);
            c3->SetLogz(); c3->SetRightMargin(0.13);
            hd->SetContour(99); hd->Draw("COLZ");
            {
                TGraph* gd = (TGraph*)g_cut.Clone("__gd_c3_cut");
                gd->SetMarkerStyle(20); gd->SetMarkerSize(0.6);
                gd->SetMarkerColor(kRed); gd->SetLineColor(kRed); gd->SetLineWidth(2);
                gd->Draw("LP SAME");
                TGraphErrors* gm = (TGraphErrors*)g_mu_err.Clone("__gd_c3_mu");
                gm->SetMarkerStyle(24); gm->SetMarkerSize(0.6);
                gm->SetMarkerColor(kRed); gm->SetLineColor(kRed); gm->SetLineWidth(1);
                gm->Draw("P SAME");
            }
            tl.DrawLatex(0.15, 0.87, lbl.c_str());
            tl.DrawLatex(0.15, 0.81, "#bullet final cut  #circ #mu_{2}#pm#sigma_{2}");
            c3->SaveAs(OutPath("ZDC_E_tot_vs_FCal_Et_AC_with_cut_debug").c_str());
            delete c3; delete hd;
        }

        // ---- Canvas 4 (debug): 1D ZDC projection of first converged slice ---
        if (first_sl.found) {
            TH1D* h1d = (TH1D*)h_zdc_fcal_->ProjectionY("__h1d_first_ev",
                first_sl.ix_lo, first_sl.ix_hi);
            h1d->SetDirectory(nullptr);
            h1d->SetTitle(";ZDC E_{total} [TeV];Counts");
            TF1 gfit2("gfit2_ev", "gaus");
            gfit2.SetRange(first_sl.lo2, first_sl.hi2);
            gfit2.SetParameters(first_sl.amp2, first_sl.mu2, first_sl.sig2);
            h1d->Fit(&gfit2, "RNQ");
            gfit2.SetLineColor(kRed); gfit2.SetLineWidth(2);
            TCanvas* c4 = new TCanvas("c_1d_first_ev", "", 750, 600);
            c4->SetLeftMargin(0.13); c4->SetBottomMargin(0.13);
            h1d->SetMarkerStyle(20); h1d->SetMarkerSize(0.7);
            h1d->Draw("E"); gfit2.Draw("SAME");
            TF1 gfit2_ext("gfit2_ext_ev", "gaus",
                first_sl.mu2 - 5.*first_sl.sig2,
                first_sl.mu2 + 5.*first_sl.sig2);
            gfit2_ext.SetParameters(gfit2.GetParameter(0),
                gfit2.GetParameter(1), gfit2.GetParameter(2));
            gfit2_ext.SetLineColor(kRed); gfit2_ext.SetLineWidth(2); gfit2_ext.SetLineStyle(2);
            gfit2_ext.Draw("SAME");
            TLatex tl4; tl4.SetNDC(); tl4.SetTextSize(0.038);
            tl4.DrawLatex(0.55, 0.82,
                Form("FCal [%.3f, %.3f] TeV", first_sl.xlo, first_sl.xhi));
            tl4.DrawLatex(0.55, 0.75, Form("#mu_{2} = %.2f", first_sl.mu2));
            tl4.DrawLatex(0.55, 0.68, Form("#sigma_{2} = %.2f", first_sl.sig2));
            tl4.DrawLatex(0.55, 0.61,
                Form("#mu_{2}+5#sigma_{2} = %.2f", first_sl.mu2 + N_SIGMA_CUT*first_sl.sig2));
            c4->SaveAs(OutPath("ZDC_1D_first_slice_fit_debug").c_str());
            delete c4; delete h1d;
        }

        // ---- Supporting plot: main-band and bg-band TGraphs separately ------
        {
            TH2D* hd = (TH2D*)h_zdc_fcal_->Clone("__h2_support");
            hd->SetDirectory(nullptr);
            TCanvas* cs = new TCanvas("c_zdc_support", "", 800, 650);
            cs->SetLogz(); cs->SetRightMargin(0.13);
            hd->SetContour(99); hd->Draw("COLZ");
            {
                TGraph* gm = (TGraph*)g_main.Clone("__gd_main_sup");
                gm->SetMarkerStyle(20); gm->SetMarkerSize(0.6);
                gm->SetMarkerColor(kRed); gm->SetLineColor(kRed); gm->SetLineWidth(2);
                gm->Draw("LP SAME");
            }
            if (have_bg) {
                TGraph* gb = (TGraph*)g_bg.Clone("__gd_bg_sup");
                gb->SetMarkerStyle(20); gb->SetMarkerSize(0.6);
                gb->SetMarkerColor(kBlue+1); gb->SetLineColor(kBlue+1); gb->SetLineWidth(2);
                gb->Draw("LP SAME");
                TGraphErrors* gmu = (TGraphErrors*)g_bg_mu_err.Clone("__gd_bg_mu_sup");
                gmu->SetMarkerStyle(21); gmu->SetMarkerSize(0.8);
                gmu->SetMarkerColor(kYellow+1); gmu->SetLineColor(kYellow+1); gmu->SetLineWidth(2);
                gmu->Draw("P SAME");
            }
            const int n_guess = (int)gr_guess_x.size();
            if (n_guess > 0) {
                TGraph* gg = new TGraph(n_guess, gr_guess_x.data(), gr_guess_mu.data());
                gg->SetMarkerStyle(22); gg->SetMarkerSize(0.7);
                gg->SetMarkerColor(kOrange+1); gg->SetLineColor(kOrange+1);
                gg->SetLineWidth(1); gg->SetLineStyle(2);
                gg->Draw("LP SAME");
            }
            TLatex ts; ts.SetNDC(); ts.SetTextSize(0.036);
            ts.DrawLatex(0.15, 0.87, lbl.c_str());
            ts.SetTextColor(kRed);
            ts.DrawLatex(0.15, 0.81, Form("red: main band #mu+%.0f#sigma", N_SIGMA_CUT));
            if (have_bg) {
                ts.SetTextColor(kBlue+1);
                ts.DrawLatex(0.15, 0.75,
                    Form("blue: bg band #mu-%.0f#sigma  (FCal > %.1f TeV)",
                         BG_N_SIGMA_CUT, FCAL_BG_FIT_MIN_TEV));
                ts.SetTextColor(kYellow+1);
                ts.DrawLatex(0.15, 0.69, "yellow: bg band #mu #pm #sigma");
            }
            if (n_guess > 0) {
                ts.SetTextColor(kOrange+1);
                ts.DrawLatex(0.15, 0.63, "orange: bg #mu initial guess");
            }
            cs->SaveAs(OutPath("cut1_ZDC_FCal_2graph_support").c_str());
            delete cs; delete hd;
        }
    }
    // -------------------------------------------------------------------------
    std::string OutPathAlt(const std::string& stem) const {
        return alt_out_dir_ + "/" + stem + "_pbpb_20" + yr_ + ".png";
    }

    // Hard-coded calibration points (x2,y2,x3,y3) for the alt-cut quadratic.
    static std::tuple<double,double,double,double> GetAltCutPts(int yr) {
        static const std::map<int, std::array<double,4>> kPts = {
            {23, {3.4, 225., 4.8, 165.}},
            {24, {3.4, 202., 4.8, 153.}},
            {25, {3.4, 200., 4.8, 152.}},
        };
        auto it = kPts.find(yr);
        if (it == kPts.end()) return {3.4, 202., 4.8, 153.};
        return {it->second[0], it->second[1], it->second[2], it->second[3]};
    }

    // Solve a*x^2 + b*x + c through three distinct points.
    static std::tuple<double,double,double> SolveQuad3Pts(
        double x1, double y1, double x2, double y2, double x3, double y3)
    {
        const double d21 = x2 - x1, d31 = x3 - x1, d32 = x3 - x2;
        const double a = ((y3 - y1)*d21 - (y2 - y1)*d31) / (d31 * d21 * d32);
        const double b = (y2 - y1)/d21 - (x2 + x1)*a;
        const double c = y1 - a*x1*x1 - b*x1;
        return {a, b, c};
    }

    // -------------------------------------------------------------------------
    // Alternative Cut 1: no background-band fitting.
    // For FCal < ALT_POLY_X1: use main-band TGraph (mu + N_SIGMA_CUT*sigma).
    // For FCal >= ALT_POLY_X1: quadratic through (ALT_POLY_X1, main_cut_at_x1)
    //   and two per-year hard-coded calibration points.
    void FitAndPlotZDCAlt() {
        static const double ALT_POLY_X1 = 2.0;  // TeV: cutoff between TGraph and quadratic

        const int nx       = h_zdc_fcal_->GetNbinsX();
        const int ny       = h_zdc_fcal_->GetNbinsY();
        const int n_slices = nx / N_BINS_PER_SLICE;

        std::vector<double> gr_x, gr_main_cut, gr_mu, gr_sigma, gr_ex;

        struct FirstSlice {
            bool found = false;
            int ix_lo = 0, ix_hi = 0;
            double xlo = 0., xhi = 0.;
            double amp2 = 0., mu2 = 0., sig2 = 0., lo2 = 0., hi2 = 0.;
        } first_sl;

        TF1 gfunc("gfit_alt", "gaus");

        for (int isl = 0; isl < n_slices; ++isl) {
            const int    ix_lo = isl * N_BINS_PER_SLICE + 1;
            const int    ix_hi = std::min((isl + 1) * N_BINS_PER_SLICE, nx);
            const double xlo   = h_zdc_fcal_->GetXaxis()->GetBinLowEdge(ix_lo);
            const double xhi   = h_zdc_fcal_->GetXaxis()->GetBinUpEdge(ix_hi);
            const double xcen  = 0.5 * (xlo + xhi);

            TH1D* hpy = (TH1D*)h_zdc_fcal_->ProjectionY(
                Form("__hpy_alt_%d", isl), ix_lo, ix_hi);
            hpy->SetDirectory(nullptr);

            if (hpy->GetEntries() < MIN_ENTRIES) { delete hpy; continue; }

            const double mu0  = hpy->GetBinCenter(hpy->GetMaximumBin());
            const double rms  = (hpy->GetRMS() > 0.) ? hpy->GetRMS() : 1.0;
            const double ymin = hpy->GetXaxis()->GetXmin();
            const double ymax = hpy->GetXaxis()->GetXmax();

            gfunc.SetRange(std::max(ymin, mu0 - 3.*rms), std::min(ymax, mu0 + 3.*rms));
            gfunc.SetParameters(hpy->GetMaximum(), mu0, rms);
            if (hpy->Fit(&gfunc, "RNQ") != 0 || gfunc.GetParameter(2) <= 0.) {
                delete hpy; continue;
            }
            const double mu1 = gfunc.GetParameter(1), sig1 = gfunc.GetParameter(2);
            const double lo2 = std::max(ymin, mu1 - 3.*sig1);
            const double hi2 = std::min(ymax, mu1 + 3.*sig1);
            gfunc.SetRange(lo2, hi2);
            gfunc.SetParameters(gfunc.GetParameter(0), mu1, sig1);
            if (hpy->Fit(&gfunc, "RNQ") != 0 || gfunc.GetParameter(2) <= 0.) {
                delete hpy; continue;
            }
            const double mu2      = gfunc.GetParameter(1);
            const double sig2     = gfunc.GetParameter(2);
            const double main_cut = mu2 + N_SIGMA_CUT * sig2;
            delete hpy;

            if (!first_sl.found)
                first_sl = {true, ix_lo, ix_hi, xlo, xhi,
                            gfunc.GetParameter(0), mu2, sig2, lo2, hi2};

            gr_x.push_back(xcen);
            gr_main_cut.push_back(main_cut);
            gr_mu.push_back(mu2);
            gr_sigma.push_back(sig2);
            gr_ex.push_back(0.5*(xhi - xlo));
        }

        const int n_fit = (int)gr_x.size();
        if (n_fit < 2) {
            std::cerr << "FitAndPlotZDCAlt: too few converged slices!" << std::endl;
            return;
        }

        TGraph       g_main  (n_fit, gr_x.data(), gr_main_cut.data());
        TGraph       g_mu    (n_fit, gr_x.data(), gr_mu.data());
        TGraphErrors g_mu_err(n_fit, gr_x.data(), gr_mu.data(),
                              gr_ex.data(), gr_sigma.data());

        // First quadratic point: main-band cut evaluated at ALT_POLY_X1
        const double y1 = EvalCut(&g_main, ALT_POLY_X1);
        const auto [x2, y2, x3, y3] = GetAltCutPts(run_year_);
        const auto [qa, qb, qc]      = SolveQuad3Pts(ALT_POLY_X1, y1, x2, y2, x3, y3);

        std::cout << Form("[Alt] poly through (%.1f,%.2f),(%.1f,%.2f),(%.1f,%.2f)"
                          "  =>  a=%.4f b=%.4f c=%.4f\n",
                          ALT_POLY_X1, y1, x2, y2, x3, y3, qa, qb, qc);

        // Build combined alt cut vector
        std::vector<double> gr_alt_cut;
        for (int i = 0; i < n_fit; ++i) {
            const double x = gr_x[i];
            gr_alt_cut.push_back(x < ALT_POLY_X1 ? gr_main_cut[i]
                                                  : qa*x*x + qb*x + qc);
        }
        TGraph g_cut(n_fit, gr_x.data(), gr_alt_cut.data());

        // ---- save to ROOT file -----------------------------------------------
        const std::string cuts_path = cuts_dir_ + "/event_sel_cuts_pbpb_20" + yr_ + "_alt.root";
        {
            TFile* fcuts = TFile::Open(cuts_path.c_str(), "RECREATE");
            if (!fcuts || fcuts->IsZombie()) {
                std::cerr << "Cannot open alt cuts file: " << cuts_path << std::endl;
            } else {
                auto wg = [&](TGraph* g, const char* name, const char* title) {
                    g->SetTitle(title); g->Write(name, TObject::kOverwrite);
                };
                wg((TGraph*)g_cut.Clone("g_ZDC_FCal_cut"),
                   "g_ZDC_FCal_cut",
                   Form("ZDC alt cut (main<%.1fTeV, quad>=%.1fTeV);FCal E_{T} [TeV];ZDC cut",
                        ALT_POLY_X1, ALT_POLY_X1));
                wg((TGraph*)g_main.Clone("g_ZDC_FCal_main_cut"),
                   "g_ZDC_FCal_main_cut",
                   Form("ZDC banana main band #mu+%.0f#sigma;FCal E_{T} [TeV];ZDC cut", N_SIGMA_CUT));
                wg((TGraph*)g_mu.Clone("g_ZDC_FCal_mu"),
                   "g_ZDC_FCal_mu",
                   "ZDC Gaussian mean per FCal slice;FCal E_{T} [TeV];ZDC mean");
                ((TGraphErrors*)g_mu_err.Clone("g_ZDC_FCal_mu_sigma"))
                    ->Write("g_ZDC_FCal_mu_sigma", TObject::kOverwrite);
                TVectorD poly(4);
                poly[0] = qa; poly[1] = qb; poly[2] = qc; poly[3] = ALT_POLY_X1;
                poly.Write("alt_poly_a_b_c_x1");
                fcuts->Close();
                std::cout << "Alt cut TGraphs saved to: " << cuts_path << std::endl;
            }
        }

        const std::string lbl = "Pb+Pb 20" + yr_ + " data (alt)";
        TLatex tl; tl.SetNDC(); tl.SetTextSize(0.036);

        // ---- Canvas 1: 2D + alt cut ------------------------------------------
        {
            TH2D* hd = (TH2D*)h_zdc_fcal_->Clone("__h2_c1_alt");
            hd->SetDirectory(nullptr);
            TCanvas* c1 = new TCanvas("c_zdc_fcal_cut_alt", "", 800, 650);
            c1->SetLogz(); c1->SetRightMargin(0.13);
            hd->SetContour(99); hd->Draw("COLZ");
            {
                TGraph* gd = (TGraph*)g_cut.Clone("__gd_c1_alt");
                gd->SetMarkerStyle(20); gd->SetMarkerSize(0.6);
                gd->SetMarkerColor(kRed); gd->SetLineColor(kRed); gd->SetLineWidth(2);
                gd->Draw("LP SAME");
            }
            tl.DrawLatex(0.15, 0.87, lbl.c_str());
            tl.DrawLatex(0.15, 0.81,
                Form("red: main #mu+%.0f#sigma (<%.0fTeV), quadratic (#geq%.0fTeV)",
                     N_SIGMA_CUT, ALT_POLY_X1, ALT_POLY_X1));
            c1->SaveAs(OutPathAlt("ZDC_E_tot_vs_FCal_Et_AC_with_cut").c_str());
            delete c1; delete hd;
        }

        // ---- Canvas 2: pass | fail -------------------------------------------
        {
            TH2D* h_pass = (TH2D*)h_zdc_fcal_->Clone("__h2_pass_alt");
            TH2D* h_fail = (TH2D*)h_zdc_fcal_->Clone("__h2_fail_alt");
            h_pass->Reset(); h_pass->SetDirectory(nullptr);
            h_fail->Reset(); h_fail->SetDirectory(nullptr);

            for (int ix = 1; ix <= nx; ++ix) {
                const double xc  = h_zdc_fcal_->GetXaxis()->GetBinCenter(ix);
                const double cut = EvalCut(&g_cut, xc);
                for (int iy = 1; iy <= ny; ++iy) {
                    const double w = h_zdc_fcal_->GetBinContent(ix, iy);
                    if (w == 0.) continue;
                    const double yc = h_zdc_fcal_->GetYaxis()->GetBinCenter(iy);
                    if (cut > 0. && yc > cut) h_fail->Fill(xc, yc, w);
                    else                       h_pass->Fill(xc, yc, w);
                }
            }

            TCanvas* c2 = new TCanvas("c_zdc_fcal_pf_alt", "", 1400, 620);
            c2->Divide(2, 1);
            auto draw_cut = [&](TGraph* src, const char* name) {
                TGraph* gd = (TGraph*)src->Clone(name);
                gd->SetMarkerStyle(20); gd->SetMarkerSize(0.6);
                gd->SetMarkerColor(kRed); gd->SetLineColor(kRed); gd->SetLineWidth(2);
                gd->Draw("LP SAME");
            };
            c2->cd(1); gPad->SetLogz(); gPad->SetRightMargin(0.15);
            h_pass->SetContour(99); h_pass->Draw("COLZ");
            draw_cut(&g_cut, "__gd_pass_alt");
            { TLatex t2; t2.SetNDC(); t2.SetTextSize(0.036);
              t2.DrawLatex(0.15, 0.87, lbl.c_str());
              t2.DrawLatex(0.15, 0.81, "Passing cut"); }
            c2->cd(2); gPad->SetLogz(); gPad->SetRightMargin(0.15);
            h_fail->SetContour(99); h_fail->Draw("COLZ");
            draw_cut(&g_cut, "__gd_fail_alt");
            { TLatex t2; t2.SetNDC(); t2.SetTextSize(0.036);
              t2.DrawLatex(0.15, 0.87, lbl.c_str());
              t2.DrawLatex(0.15, 0.81, "Failing cut"); }
            c2->SaveAs(OutPathAlt("ZDC_E_tot_vs_FCal_Et_AC_pass_fail").c_str());
            delete c2; delete h_pass; delete h_fail;
        }

    }

    // -------------------------------------------------------------------------
    // ZDC time A vs C for the top 5 individual centrality percent bins (0-4)
    // and the combined 5-10% bin (centrality 5-9), for events passing the
    // banana (ZDC vs FCal) cut.  Centrality is 0-based (0 = 0-1%, most central).
    void MakeZDCTimeCentralityPlot() {
        // Load banana cut TGraph from file
        const std::string cuts_path = cuts_dir_ + "/event_sel_cuts_pbpb_20" + yr_ + ".root";
        TFile* fcuts = TFile::Open(cuts_path.c_str(), "READ");
        if (!fcuts || fcuts->IsZombie()) {
            std::cerr << "Cannot open cuts file: " << cuts_path << std::endl;
            return;
        }
        TGraph* g_raw = (TGraph*)fcuts->Get("g_ZDC_FCal_cut");
        if (!g_raw) { fcuts->Close(); std::cerr << "g_ZDC_FCal_cut not found" << std::endl; return; }
        TGraph g_cut(*g_raw);   // clone before file close
        fcuts->Close();

        // Book: indices 0-4 = centrality 0-4 (1% each); index 5 = centrality 5-9 (5-10%)
        TH2D* hh[6] = {};
        const char* clabels[6] = {"0-1%","1-2%","2-3%","3-4%","4-5%","5-10%"};
        for (int ic = 0; ic < 6; ++ic) {
            hh[ic] = new TH2D(Form("h_zdc_t_cent%d", ic),
                Form(";ZDC t^{A} [ns];ZDC t^{C} [ns]"),
                100, -10., 10., 100, -10., 10.);
            hh[ic]->SetDirectory(nullptr);
        }

        // Event loop
        TChain chain("HeavyIonD3PD", "HeavyIonD3PD");
        for (const auto& f : infiles_) {
            if (!gSystem->AccessPathName(f.c_str())) chain.Add(f.c_str());
        }
        chain.SetMakeClass(1);
        chain.SetBranchStatus("*", 0);
        Int_t   b_HLT_mu4 = 0, centrality = -1;
        Float_t FCal_Et_P = 0.f, FCal_Et_N = 0.f;
        Float_t zdc_E[2] = {}, zdc_t[2] = {};
        chain.SetBranchStatus("b_HLT_mu4_L1MU3V", 1);
        chain.SetBranchStatus("centrality",        1);
        chain.SetBranchStatus("FCal_Et_P",         1);
        chain.SetBranchStatus("FCal_Et_N",         1);
        chain.SetBranchStatus("zdc_ZdcEnergy",     1);
        chain.SetBranchStatus("zdc_ZdcTime",       1);
        chain.SetBranchAddress("b_HLT_mu4_L1MU3V", &b_HLT_mu4);
        chain.SetBranchAddress("centrality",        &centrality);
        chain.SetBranchAddress("FCal_Et_P",         &FCal_Et_P);
        chain.SetBranchAddress("FCal_Et_N",         &FCal_Et_N);
        chain.SetBranchAddress("zdc_ZdcEnergy",     zdc_E);
        chain.SetBranchAddress("zdc_ZdcTime",       zdc_t);

        const Long64_t n = chain.GetEntries();
        std::cout << "ZDC time centrality pass: processing " << n << " events" << std::endl;

        for (Long64_t i = 0; i < n; ++i) {
            chain.GetEntry(i);
            if (!b_HLT_mu4) continue;
            if (centrality < 0 || centrality > 9) continue;  // keep 0-9 only

            const float fcal_AC = (FCal_Et_P + FCal_Et_N) * 1e-6f;
            const float zdc_tot = (zdc_E[0] + zdc_E[1]) / 1000.f;
            // Apply banana cut (interpolated)
            if (EvalCut(&g_cut, fcal_AC) > 0. && zdc_tot > EvalCut(&g_cut, fcal_AC)) continue;

            const float tA = zdc_t[1], tC = zdc_t[0];  // [1]=A, [0]=C
            int ic = (centrality <= 4) ? centrality : 5;  // 0-4 individual, 5-9 → bin 5
            hh[ic]->Fill(tA, tC);
        }
        std::cout << "ZDC time centrality histograms filled." << std::endl;

        // ---- 2×3 canvas ----
        gStyle->SetOptStat(0);
        gStyle->SetPalette(kBird);
        TCanvas* c = new TCanvas("c_zdc_t_centr", "", 1800, 1200);
        c->Divide(3, 2);
        for (int ic = 0; ic < 6; ++ic) {
            c->cd(ic + 1);
            gPad->SetLogz();
            gPad->SetRightMargin(0.14);
            hh[ic]->SetContour(99);
            hh[ic]->Draw("COLZ");
            TLatex tl; tl.SetNDC(); tl.SetTextSize(0.055);
            tl.DrawLatex(0.15, 0.87, Form("Pb+Pb 20%s, centr. %s", yr_.c_str(), clabels[ic]));
            tl.DrawLatex(0.15, 0.79, "Banana cut passed");
        }
        c->SaveAs(OutPath("ZDC_time_AC_corr_top_5_centr").c_str());
        std::cout << "Saved: " << OutPath("ZDC_time_AC_corr_top_5_centr") << std::endl;
        delete c;
        for (int ic = 0; ic < 6; ++ic) delete hh[ic];
    }
};

// =============================================================================
void plot_pbpb_event_sel_event_level(int run_year = 24) {
    PbPbEventSel plotter(run_year);
    plotter.Run();
}

// Run only the centrality ZDC time plot (requires cuts ROOT file from prior Run()).
void plot_pbpb_zdc_time_centr(int run_year = 24) {
    PbPbEventSel plotter(run_year);
    plotter.RunCentralityOnly();
}
