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

// ---- ZDC fit tuning ---------------------------------------------------------
// N_BINS_PER_SLICE=2 gives 0.10 TeV slices (120 bins × 0.05 TeV/bin).
// With ~6.8M trigger-passing events, 55+/60 slices converge (min ~110 entries),
// giving dense TGraph interpolation points vs the old 0.25 TeV step-function.
static const int    N_BINS_PER_SLICE = 2;
static const int    MIN_ENTRIES      = 100;
static const double N_SIGMA_CUT      = 5.0;
static const double FCAL_SPLIT_TEV   = 4.6;   // FCal ET above which ZDC bands split; extrapolate linearly beyond

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
        yr_       = std::to_string(run_year_);
        out_dir_  = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/"
                    "single_b_analysis/event_selection/pbpb_20" + yr_;
        cuts_dir_ = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_20" + yr_;
        auto file_map = BuildFileMap();
        if (file_map.find(run_year_) == file_map.end())
            throw std::runtime_error("No input files configured for run year " + yr_);
        infiles_ = file_map.at(run_year_);
    }

    void Run() {
        gStyle->SetOptStat(0);
        gStyle->SetPalette(kBird);
        gSystem->mkdir(out_dir_.c_str(), true);
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
    std::string yr_, out_dir_, cuts_dir_;
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
        const int nx      = h_zdc_fcal_->GetNbinsX();
        const int ny      = h_zdc_fcal_->GetNbinsY();
        const int n_slices = nx / N_BINS_PER_SLICE;

        std::vector<double> gr_x, gr_cut, gr_mu, gr_sigma, gr_ex;

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

            // Above FCAL_SPLIT_TEV the signal and pile-up bands split; skip fitting
            // entirely here — the cut is linearly extrapolated below.
            if (xcen > FCAL_SPLIT_TEV) continue;

            TH1D* hpy = (TH1D*)h_zdc_fcal_->ProjectionY(
                Form("__hpy_ev_%d", isl), ix_lo, ix_hi);
            hpy->SetDirectory(nullptr);

            if (hpy->GetEntries() < MIN_ENTRIES) { delete hpy; continue; }

            const double mu0  = hpy->GetBinCenter(hpy->GetMaximumBin());
            const double rms  = (hpy->GetRMS() > 0.) ? hpy->GetRMS() : 1.0;
            const double ymin = hpy->GetXaxis()->GetXmin();
            const double ymax = hpy->GetXaxis()->GetXmax();

            // Pass 1
            gfunc.SetRange(std::max(ymin, mu0 - 3.*rms), std::min(ymax, mu0 + 3.*rms));
            gfunc.SetParameters(hpy->GetMaximum(), mu0, rms);
            if (hpy->Fit(&gfunc, "RNQ") != 0 || gfunc.GetParameter(2) <= 0.) {
                delete hpy; continue;
            }
            const double mu1 = gfunc.GetParameter(1), sig1 = gfunc.GetParameter(2);

            // Pass 2
            const double lo2 = std::max(ymin, mu1 - 3.*sig1);
            const double hi2 = std::min(ymax, mu1 + 3.*sig1);
            gfunc.SetRange(lo2, hi2);
            gfunc.SetParameters(gfunc.GetParameter(0), mu1, sig1);
            if (hpy->Fit(&gfunc, "RNQ") != 0 || gfunc.GetParameter(2) <= 0.) {
                delete hpy; continue;
            }

            const double mu2  = gfunc.GetParameter(1);
            const double sig2 = gfunc.GetParameter(2);
            const double cut  = mu2 + N_SIGMA_CUT * sig2;

            if (!first_sl.found)
                first_sl = {true, ix_lo, ix_hi, xlo, xhi,
                            gfunc.GetParameter(0), mu2, sig2, lo2, hi2};
            delete hpy;

            gr_x.push_back(xcen);
            gr_cut.push_back(cut);
            gr_mu.push_back(mu2);
            gr_sigma.push_back(sig2);
            gr_ex.push_back(0.5*(xhi-xlo));

            std::cout << Form("  slice %2d [%.3f-%.3f TeV]  mu=%.3f  sig=%.3f  cut=%.3f TeV\n",
                              isl, xlo, xhi, mu2, sig2, cut);
        }
        const int n_fit = (int)gr_x.size();
        std::cout << "Fitting: " << n_fit << " slices converged (FCal <= "
                  << FCAL_SPLIT_TEV << " TeV)." << std::endl;

        // ---- linear extrapolation beyond FCAL_SPLIT_TEV using last 3 points ---
        if (n_fit < 3) {
            std::cerr << "FitAndPlotZDC: too few converged slices for extrapolation!" << std::endl;
            return;
        }
        {
            TGraph g_last3(3, gr_x.data() + n_fit - 3, gr_cut.data() + n_fit - 3);
            TF1 f_extrap("__f_extrap_c1", "pol1",
                         gr_x[n_fit - 3], gr_x[n_fit - 1]);
            g_last3.Fit(&f_extrap, "RNQ");
            const double slope     = f_extrap.GetParameter(1);
            const double intercept = f_extrap.GetParameter(0);
            std::cout << Form("Banana extrapolation (last 3 pts, linear): "
                              "intercept=%.3f  slope=%.4f\n", intercept, slope);
            std::cout << Form("  anchor pts: (%.3f, %.3f) (%.3f, %.3f) (%.3f, %.3f)\n",
                              gr_x[n_fit-3], gr_cut[n_fit-3],
                              gr_x[n_fit-2], gr_cut[n_fit-2],
                              gr_x[n_fit-1], gr_cut[n_fit-1]);
            // Append extrapolated points up to the histogram edge
            const double fcal_end = h_zdc_fcal_->GetXaxis()->GetXmax();
            const int n_ext = 8;
            for (int k = 0; k < n_ext; ++k) {
                double x = FCAL_SPLIT_TEV + (k + 1) * (fcal_end - FCAL_SPLIT_TEV) / n_ext;
                gr_x.push_back(x);
                gr_cut.push_back(f_extrap.Eval(x));
            }
        }

        // ---- build TGraph objects for interpolation and persistence ----------
        // g_cut includes the extrapolated tail; g_mu and g_mu_err use only fitted pts.
        TGraph  g_cut((int)gr_x.size(), gr_x.data(), gr_cut.data());
        TGraph  g_mu (n_fit, gr_x.data(), gr_mu.data());
        TGraphErrors g_mu_err(n_fit, gr_x.data(), gr_mu.data(),
                              gr_ex.data(), gr_sigma.data());

        // Sorted by x is guaranteed because slices are processed in order.

        // ---- save cut curves to ROOT file ------------------------------------
        const std::string cuts_path = cuts_dir_ + "/event_sel_cuts_pbpb_20" + yr_ + ".root";
        {
            TFile* fcuts = TFile::Open(cuts_path.c_str(), "UPDATE");
            if (!fcuts || fcuts->IsZombie()) {
                std::cerr << "Cannot open cuts file for writing: " << cuts_path << std::endl;
            } else {
                TGraph* g_c = (TGraph*)g_cut.Clone("g_ZDC_FCal_cut");
                TGraph* g_m = (TGraph*)g_mu.Clone("g_ZDC_FCal_mu");
                TGraphErrors* g_s = (TGraphErrors*)g_mu_err.Clone("g_ZDC_FCal_mu_sigma");
                g_c->SetTitle(Form("ZDC out-of-time pileup cut (mu+%.0fsigma);FCal E_{T}^{A+C} [TeV];ZDC E_{total} cut [TeV]", N_SIGMA_CUT));
                g_m->SetTitle("ZDC Gaussian mean per FCal slice;FCal E_{T}^{A+C} [TeV];ZDC mean [TeV]");
                g_s->SetTitle("ZDC Gaussian mean #pm sigma per FCal slice;FCal E_{T}^{A+C} [TeV];ZDC mean [TeV]");
                g_c->Write("g_ZDC_FCal_cut",  TObject::kOverwrite);
                g_m->Write("g_ZDC_FCal_mu",   TObject::kOverwrite);
                g_s->Write("g_ZDC_FCal_mu_sigma", TObject::kOverwrite);
                fcuts->Close();
                std::cout << "Cut TGraphs saved to: " << cuts_path << std::endl;
            }
        }

        const std::string lbl = "Pb+Pb 20" + yr_ + " data";
        TLatex tl; tl.SetNDC(); tl.SetTextSize(0.036);

        // ---- Canvas 1: 2D + interpolated cut curve --------------------------
        {
            TH2D* hd = (TH2D*)h_zdc_fcal_->Clone("__h2_c1_ev");
            hd->SetDirectory(nullptr);
            TCanvas* c1 = new TCanvas("c_zdc_fcal_cut_ev", "", 800, 650);
            c1->SetLogz(); c1->SetRightMargin(0.13);
            hd->SetContour(99); hd->Draw("COLZ");
            if (n_fit > 0) {
                TGraph* gd = (TGraph*)g_cut.Clone("__gd_c1");
                gd->SetMarkerStyle(20); gd->SetMarkerSize(0.6);
                gd->SetMarkerColor(kRed); gd->SetLineColor(kRed); gd->SetLineWidth(2);
                gd->Draw("LP SAME");  // line through interpolation points
            }
            tl.DrawLatex(0.15, 0.87, lbl.c_str());
            tl.DrawLatex(0.15, 0.81, Form("red: #mu_{2} + %.0f#sigma_{2} (TGraph interp.)", N_SIGMA_CUT));
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
                const double cut = (n_fit > 0) ? EvalCut(&g_cut, xc) : -1.;
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

            c2->cd(1); gPad->SetLogz(); gPad->SetRightMargin(0.15);
            h_pass->SetContour(99); h_pass->Draw("COLZ");
            if (n_fit > 0) {
                TGraph* gd = (TGraph*)g_cut.Clone("__gd_pass");
                gd->SetMarkerStyle(20); gd->SetMarkerSize(0.6);
                gd->SetMarkerColor(kRed); gd->SetLineColor(kRed); gd->SetLineWidth(2);
                gd->Draw("LP SAME");
            }
            {
                TLatex t2; t2.SetNDC(); t2.SetTextSize(0.036);
                t2.DrawLatex(0.15, 0.87, lbl.c_str());
                t2.DrawLatex(0.15, 0.81, Form("Passing: ZDC #leq interp. cut (#mu_{2}+%.0f#sigma_{2})", N_SIGMA_CUT));
            }

            c2->cd(2); gPad->SetLogz(); gPad->SetRightMargin(0.15);
            h_fail->SetContour(99); h_fail->Draw("COLZ");
            if (n_fit > 0) {
                TGraph* gd = (TGraph*)g_cut.Clone("__gd_fail");
                gd->SetMarkerStyle(20); gd->SetMarkerSize(0.6);
                gd->SetMarkerColor(kRed); gd->SetLineColor(kRed); gd->SetLineWidth(2);
                gd->Draw("LP SAME");
            }
            {
                TLatex t2; t2.SetNDC(); t2.SetTextSize(0.036);
                t2.DrawLatex(0.15, 0.87, lbl.c_str());
                t2.DrawLatex(0.15, 0.81, Form("Failing: ZDC > interp. cut (#mu_{2}+%.0f#sigma_{2})", N_SIGMA_CUT));
            }
            c2->SaveAs(OutPath("ZDC_E_tot_vs_FCal_Et_AC_pass_fail").c_str());
            delete c2; delete h_pass; delete h_fail;
        }

        // ---- Canvas 3 (debug): 2D + mu+5sigma line + mu±sigma error bars ----
        {
            TH2D* hd = (TH2D*)h_zdc_fcal_->Clone("__h2_c3_ev");
            hd->SetDirectory(nullptr);
            TCanvas* c3 = new TCanvas("c_zdc_fcal_dbg_ev", "", 800, 650);
            c3->SetLogz(); c3->SetRightMargin(0.13);
            hd->SetContour(99); hd->Draw("COLZ");
            if (n_fit > 0) {
                TGraph* gd_cut = (TGraph*)g_cut.Clone("__gd_c3_cut");
                gd_cut->SetMarkerStyle(20); gd_cut->SetMarkerSize(0.6);
                gd_cut->SetMarkerColor(kRed); gd_cut->SetLineColor(kRed); gd_cut->SetLineWidth(2);
                gd_cut->Draw("LP SAME");
                TGraphErrors* gd_mu = (TGraphErrors*)g_mu_err.Clone("__gd_c3_mu");
                gd_mu->SetMarkerStyle(24); gd_mu->SetMarkerSize(0.6);
                gd_mu->SetMarkerColor(kRed); gd_mu->SetLineColor(kRed); gd_mu->SetLineWidth(1);
                gd_mu->Draw("P SAME");
            }
            tl.DrawLatex(0.15, 0.87, lbl.c_str());
            tl.DrawLatex(0.15, 0.81, Form("#bullet  #mu_{2}+%.0f#sigma_{2}  (filled, interpolated)", N_SIGMA_CUT));
            tl.DrawLatex(0.15, 0.75, "#circ  #mu_{2} #pm #sigma_{2}  (open, x = slice half-width)");
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
                Form("FCal E_{T}^{A+C} #in [%.3f, %.3f] TeV", first_sl.xlo, first_sl.xhi));
            tl4.DrawLatex(0.55, 0.75, Form("#mu_{2} = %.2f TeV", first_sl.mu2));
            tl4.DrawLatex(0.55, 0.68, Form("#sigma_{2} = %.2f TeV", first_sl.sig2));
            tl4.DrawLatex(0.55, 0.61,
                Form("#mu_{2} + 5#sigma_{2} = %.2f TeV", first_sl.mu2 + N_SIGMA_CUT*first_sl.sig2));
            tl4.SetTextSize(0.032);
            tl4.DrawLatex(0.55, 0.53,
                Form("Fit range: [%.2f, %.2f] TeV", first_sl.lo2, first_sl.hi2));
            c4->SaveAs(OutPath("ZDC_1D_first_slice_fit_debug").c_str());
            delete c4; delete h1d;
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
