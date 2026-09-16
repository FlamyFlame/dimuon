// Single-muon reconstruction efficiency: R17618 (full HIJING truth -> Pythia/HIJING
// barcode collision) vs R17662 (signal-only truth -> no collision), 0-5% and 5-10%
// centrality, with an R17618/R17662 RATIO panel (the two are expected to agree, so the
// ratio is the quantity of interest -> it belongs on the plot).
//
// Both samples: pure prob>0.5 barcode matching (the ad-hoc dR fallback was deleted
// 2026-07-10) and the same pT-hat slice (pTH8_14 / kin0), so the muon populations match.
//
// NOTE on what the two r-tags are: they are the SAME generated events (100% eventNumber
// overlap, identical Pythia truth, identical overlaid HIJING event -> FCal_Et is
// bit-identical), but TWO DIFFERENT digitisation+reco passes: `digiSteeringConf:
// StandardSignalOnlyTruth` (the only AMI difference) shifts the tracking-detector RNG
// stream. Physics is unchanged -- detector response (momentum resolution), reco-muon
// multiplicity and efficiency all agree statistically (KS p = 0.75 / 0.94 / 1.00) --
// but the two are independent realisations, so the ratio errors are treated as
// uncorrelated and per-muon differences are expected.
//
// x-axis: 4-20 GeV in 6 LOG-spaced bins. The pTH8_14 sample has essentially no muons
// above ~20 GeV, so a wider axis would strand the data in a small fraction of the range.
// Efficiency = pass_{WP} / all fiducial truth muons vs truth pT (weighted, binomial
// errors), modelled on plot_single_muon_reco_effcy.cxx.
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TPad.h>
#include <TLine.h>
#include <cmath>
#include <ROOT/RDataFrame.hxx>
#include <string>
#include <vector>
#include <iostream>
#include "../../RDFBasedHistFilling/CommonEffcyConfig.h"

namespace {
    // The pTH8_14 sample populates muon pT only up to ~20 GeV -- above that the bins are
    // empty/statistics-starved, so the axis stops at 20 GeV (an axis running to 100 GeV
    // would leave the data in a small fraction of the range).
    // 6 bins, LOG-spaced (the axis is drawn log-x, so the binning must be log too):
    // edges = 4.5 * (20/4.5)^(i/6).
    const std::vector<double> pt_edges = [](){
        // low edge 4 -> 4.5 with the offline muon cut (2026-09-08)
        std::vector<double> e; const double lo = 4.5, hi = 20.0;
        const int nb = 6;
        for (int i = 0; i <= nb; ++i) e.push_back(lo * std::pow(hi / lo, double(i) / nb));
        return e;
    }();
    const int nPtBins = (int)pt_edges.size() - 1;

    // Both pure prob>0.5, pTH8_14 (kin0). R17618 uses only the kin0 slice so the
    // pT-hat mix matches R17662 exactly (fair comparison of the same muon population).
    const std::string R17618 =
        "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_hijing_overlay_test_sample_pbpb23/r17618_kin0_run/"
        "muon_pairs_pythia_fullsim_hijing_overlay_pbpb23_no_data_resonance_cuts_r17618kin0_single_muon.root";
    const std::string R17662 =
        "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_hijing_overlay_test_sample_pbpb23/r17662_run/"
        "muon_pairs_pythia_fullsim_hijing_overlay_pbpb23_no_data_resonance_cuts_r17662_single_muon.root";
    std::string OUTDIR =
        "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_hijing_overlay_test_sample_pbpb23/plots/r17618_vs_r17662_comparison/";

    // Muon working point (NOMINAL = Tight; set by the entry function). Only the quality bit
    // differs. make_eff() and the axis titles read these. Medium routes to a distinct outdir.
    std::string wp_filter = "pass_tight";
    std::string wp_label  = "tight";

    // Single-muon reco efficiency vs truth pT for one file, centrality [clo,chi),
    // and optional q*eta window (integrated=true ignores qeta). Same low-stat
    // suppression as plot_single_muon_reco_effcy.cxx.
    TH1D* make_eff(const std::string& file, int clo, int chi,
                   bool integrated, float qlo, float qhi, const char* name) {
        ROOT::RDataFrame df("muon_tree", file);
        auto d = df.Define("q_eta", "truth_charge * truth_eta")
                   .Filter([clo, chi](int c){ return c >= clo && c < chi; }, {"ev_centrality"});
        ROOT::RDF::RNode dq = d;
        if (!integrated)
            dq = d.Filter([qlo, qhi](float qe){ return qe >= qlo && qe < qhi; }, {"q_eta"});

        auto h_denom = dq.Histo1D({"h_denom", "", nPtBins, pt_edges.data()}, "truth_pt", "ev_weight");
        auto h_num   = dq.Filter(wp_filter).Histo1D({"h_num", "", nPtBins, pt_edges.data()}, "truth_pt", "ev_weight");

        TH1D* hd = (TH1D*)h_denom->Clone((std::string(name) + "_d").c_str());
        TH1D* hn = (TH1D*)h_num->Clone((std::string(name) + "_n").c_str());
        hd->SetDirectory(nullptr); hn->SetDirectory(nullptr);

        TH1D* eff = (TH1D*)hn->Clone(name);
        eff->SetDirectory(nullptr);
        eff->Divide(hn, hd, 1.0, 1.0, "B");
        for (int ib = 1; ib <= eff->GetNbinsX(); ++ib) {
            double den_val = hd->GetBinContent(ib);
            double den_err = hd->GetBinError(ib);
            if (den_err <= 0 || (den_val / den_err) * (den_val / den_err) < 10.0) {
                eff->SetBinContent(ib, 0);
                eff->SetBinError(ib, 0);
            }
        }
        delete hd; delete hn;
        return eff;
    }

    void style_eff(TH1D* h, int color, int mstyle) {
        h->SetLineColor(color); h->SetMarkerColor(color);
        h->SetMarkerStyle(mstyle); h->SetMarkerSize(1.0); h->SetLineWidth(2);
    }

    // R17618 / R17662 ratio, with errors propagated as independent (the two r-tags are
    // separate digitisation+reco realisations of the same generated events, so their
    // reconstruction fluctuations are independent). Bins where either input was
    // low-stat-suppressed (content 0) are dropped.
    TH1D* make_ratio(const TH1D* a, const TH1D* b, const char* name) {
        TH1D* r = (TH1D*)a->Clone(name);
        r->SetDirectory(nullptr);
        for (int i = 1; i <= r->GetNbinsX(); ++i) {
            double va = a->GetBinContent(i), vb = b->GetBinContent(i);
            double ea = a->GetBinError(i),   eb = b->GetBinError(i);
            if (va <= 0 || vb <= 0) { r->SetBinContent(i, 0); r->SetBinError(i, 0); continue; }
            double q = va / vb;
            r->SetBinContent(i, q);
            r->SetBinError(i, q * std::sqrt((ea/va)*(ea/va) + (eb/vb)*(eb/vb)));
        }
        r->SetLineColor(kBlack); r->SetMarkerColor(kBlack);
        r->SetMarkerStyle(20); r->SetMarkerSize(1.0); r->SetLineWidth(2);
        return r;
    }

    // Style a ratio-panel frame (bottom pad).
    // Auto-range from the ratio points (+errors) so none is ever clipped (R3).
    void set_ratio_range(TH1D* f, const TH1D* r) {
        double lo = 1.0, hi = 1.0;
        for (int i = 1; i <= r->GetNbinsX(); ++i) {
            double v = r->GetBinContent(i), e = r->GetBinError(i);
            if (v <= 0) continue;
            lo = std::min(lo, v - e); hi = std::max(hi, v + e);
        }
        double pad = 0.10 * std::max(hi - lo, 0.30);
        f->SetMinimum(std::max(0.0, std::min(lo - pad, 0.85)));
        f->SetMaximum(std::max(hi + pad, 1.15));
    }

    void style_ratio_frame(TH1D* f, double titsize, double labsize) {
        f->SetMinimum(0.70); f->SetMaximum(1.30);   // overridden by set_ratio_range()
        f->GetYaxis()->SetTitle("R17618 / R17662");
        f->GetYaxis()->SetNdivisions(505);
        f->GetXaxis()->SetTitle("Truth p_{T} [GeV]");
        f->GetXaxis()->SetTitleSize(titsize); f->GetYaxis()->SetTitleSize(titsize);
        f->GetXaxis()->SetLabelSize(labsize); f->GetYaxis()->SetLabelSize(labsize);
        f->GetXaxis()->SetMoreLogLabels(true); f->GetXaxis()->SetNoExponent(true);
    }

    TH1D* make_frame(const char* name, double titsize = 0.045, double labsize = 0.038) {
        TH1D* f = new TH1D(name, "", nPtBins, pt_edges.data());
        f->SetMinimum(0.0); f->SetMaximum(1.15);
        f->GetXaxis()->SetTitle("Truth p_{T} [GeV]");
        f->GetYaxis()->SetTitle(("Single-muon reco efficiency (" + wp_label + ")").c_str());
        f->GetXaxis()->SetTitleSize(titsize); f->GetYaxis()->SetTitleSize(titsize);
        f->GetXaxis()->SetLabelSize(labsize); f->GetYaxis()->SetLabelSize(labsize);
        f->GetXaxis()->SetMoreLogLabels(true); f->GetXaxis()->SetNoExponent(true);
        return f;
    }
}

void plot_single_muon_reco_effcy_r17618_vs_r17662(bool useTight = true) {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    // NOMINAL muon WP = Tight; pass false for Medium (routed to a distinct subdir).
    wp_filter = useTight ? "pass_tight" : "pass_medium";
    wp_label  = useTight ? "tight"      : "medium";
    if (!useTight) OUTDIR += "medium_wp/";
    gSystem->mkdir(OUTDIR.c_str(), kTRUE);

    static const CommonEffcyConfig ecfg;
    const QEtaBinning& qeta_bins = ecfg.q_eta_proj_ranges_coarse_incl_gap;
    const int nQEta = (int)qeta_bins.size();

    struct CtrBin { int lo, hi; std::string suffix, label; };
    std::vector<CtrBin> ctr_bins = { {0, 5, "ctr0_5", "0-5%"}, {5, 10, "ctr5_10", "5-10%"} };

    for (const auto& ctr : ctr_bins) {

        // ---- (A) q*eta-integrated overlay + R17618/R17662 ratio panel ----
        {
            TH1D* e618 = make_eff(R17618, ctr.lo, ctr.hi, true, 0, 0, ("eint618_" + ctr.suffix).c_str());
            TH1D* e662 = make_eff(R17662, ctr.lo, ctr.hi, true, 0, 0, ("eint662_" + ctr.suffix).c_str());
            style_eff(e618, kBlack,  20);
            style_eff(e662, kRed+1,  24);
            TH1D* ratio = make_ratio(e618, e662, ("rint_" + ctr.suffix).c_str());

            TCanvas c("c_int", "c_int", 900, 800);
            TPad* pTop = new TPad("pTop", "", 0.0, 0.32, 1.0, 1.0);
            TPad* pBot = new TPad("pBot", "", 0.0, 0.0,  1.0, 0.32);
            pTop->SetBottomMargin(0.02); pTop->SetTicks(1,1); pTop->SetLogx(true);
            pBot->SetTopMargin(0.03); pBot->SetBottomMargin(0.32); pBot->SetTicks(1,1); pBot->SetLogx(true);
            pTop->Draw(); pBot->Draw();

            pTop->cd();
            TH1D* frame = make_frame(("fint_" + ctr.suffix).c_str());
            frame->GetXaxis()->SetLabelSize(0);   // x labels live on the ratio pad
            frame->GetXaxis()->SetTitleSize(0);
            frame->Draw();
            e618->Draw("PE1 SAME");
            e662->Draw("PE1 SAME");

            TLegend leg(0.45, 0.13, 0.92, 0.30);
            leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.040);
            leg.AddEntry(e618, "R17618 (full HIJING truth)", "lpe");
            leg.AddEntry(e662, "R17662 (signal-only truth)", "lpe");
            leg.Draw("SAME");

            TLatex lat; lat.SetNDC();
            lat.SetTextSize(0.048); lat.DrawLatex(0.16, 0.93, "Pythia fullsim HIJING overlay");
            lat.SetTextSize(0.042); lat.DrawLatex(0.16, 0.87, "q#eta-integrated");
            lat.DrawLatex(0.16, 0.81, ("Centrality: " + ctr.label).c_str());

            pBot->cd();
            TH1D* rframe = new TH1D(("rf_" + ctr.suffix).c_str(), "", nPtBins, pt_edges.data());
            style_ratio_frame(rframe, 0.095, 0.080);
            set_ratio_range(rframe, ratio);
            rframe->GetYaxis()->SetTitleOffset(0.52);
            rframe->GetXaxis()->SetTitleOffset(1.25);
            rframe->Draw();
            TLine* l1 = new TLine(pt_edges.front(), 1.0, pt_edges.back(), 1.0);
            l1->SetLineStyle(2); l1->SetLineColor(kRed+1); l1->Draw("SAME");
            ratio->Draw("PE1 SAME");

            std::string out = OUTDIR + "single_muon_reco_effcy_vs_pt_integrated_" + ctr.suffix + "_r17618_vs_r17662.png";
            c.SaveAs(out.c_str());
            std::cout << "  Saved: " << out << "\n";
            delete frame; delete rframe; delete e618; delete e662; delete ratio;
        }

        // ---- (B) q*eta-binned: one bin per subplot (3x3), each with its own ratio panel ----
        {
            const int nrows = 3, ncols = 3;
            TCanvas c("c_sub", "c_sub", 1500, 1250);
            std::vector<TH1D*> frames, e618s, e662s, ratios;

            for (int iq = 0; iq < nQEta; ++iq) {
                float qlo = qeta_bins[iq].first, qhi = qeta_bins[iq].second;
                const std::string sfx = ctr.suffix + "_q" + std::to_string(iq);
                TH1D* e618 = make_eff(R17618, ctr.lo, ctr.hi, false, qlo, qhi, ("e618_" + sfx).c_str());
                TH1D* e662 = make_eff(R17662, ctr.lo, ctr.hi, false, qlo, qhi, ("e662_" + sfx).c_str());
                style_eff(e618, kBlack, 20); e618->SetMarkerSize(0.8);
                style_eff(e662, kRed+1, 24); e662->SetMarkerSize(0.8);
                TH1D* rat = make_ratio(e618, e662, ("rsub_" + sfx).c_str());
                rat->SetMarkerSize(0.8);
                e618s.push_back(e618); e662s.push_back(e662); ratios.push_back(rat);

                // cell geometry: each of the 9 cells is split into an efficiency pad (top)
                // and a ratio pad (bottom).
                const int col = iq % ncols, row = iq / ncols;
                const double x0 = col / double(ncols), x1 = (col + 1) / double(ncols);
                const double yTop = 1.0 - row / double(nrows);
                const double yBot = 1.0 - (row + 1) / double(nrows);
                const double ySplit = yBot + 0.32 * (yTop - yBot);

                c.cd(0);
                TPad* pe = new TPad(("pe" + sfx).c_str(), "", x0, ySplit, x1, yTop);
                TPad* pr = new TPad(("pr" + sfx).c_str(), "", x0, yBot,   x1, ySplit);
                pe->SetTicks(1,1); pe->SetLogx(true);
                pe->SetLeftMargin(0.17); pe->SetBottomMargin(0.02); pe->SetTopMargin(0.07); pe->SetRightMargin(0.04);
                pr->SetTicks(1,1); pr->SetLogx(true);
                pr->SetLeftMargin(0.17); pr->SetBottomMargin(0.34); pr->SetTopMargin(0.03); pr->SetRightMargin(0.04);
                pe->Draw(); pr->Draw();

                pe->cd();
                TH1D* hf = make_frame(("fsub_" + sfx).c_str(), 0.075, 0.062);
                hf->GetYaxis()->SetTitle(("#varepsilon_{reco} (" + wp_label + ")").c_str());
                hf->GetYaxis()->SetTitleOffset(1.05);
                hf->GetXaxis()->SetLabelSize(0); hf->GetXaxis()->SetTitleSize(0);
                hf->Draw();
                frames.push_back(hf);
                e618->Draw("PE1 SAME");
                e662->Draw("PE1 SAME");

                char buf[64];
                snprintf(buf, sizeof(buf), "%.1f < q#eta < %.1f", qlo, qhi);
                TLatex lab; lab.SetNDC(); lab.SetTextSize(0.075);
                lab.DrawLatex(0.22, 0.86, buf);

                if (iq == 0) {
                    TLegend* leg = new TLegend(0.24, 0.08, 0.97, 0.30);
                    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.068);
                    leg->AddEntry(e618, "R17618 (full HIJING truth)", "lpe");
                    leg->AddEntry(e662, "R17662 (signal-only truth)", "lpe");
                    leg->Draw("SAME");
                }

                pr->cd();
                TH1D* rf = new TH1D(("rfsub_" + sfx).c_str(), "", nPtBins, pt_edges.data());
                style_ratio_frame(rf, 0.135, 0.115);
                set_ratio_range(rf, rat);
                rf->GetYaxis()->SetTitleOffset(0.45);
                rf->GetXaxis()->SetTitleOffset(1.05);
                rf->Draw();
                frames.push_back(rf);
                TLine* ln = new TLine(pt_edges.front(), 1.0, pt_edges.back(), 1.0);
                ln->SetLineStyle(2); ln->SetLineColor(kRed+1); ln->Draw("SAME");
                rat->Draw("PE1 SAME");
            }

            c.cd(0);
            TLatex tit; tit.SetNDC();
            tit.SetTextSize(0.022);
            tit.DrawLatex(0.10, 0.985, "Pythia fullsim HIJING overlay  |  single-muon reco efficiency");
            tit.SetTextSize(0.020);
            tit.DrawLatex(0.62, 0.985, ("Centrality: " + ctr.label).c_str());

            std::string out = OUTDIR + "single_muon_reco_effcy_vs_pt_subplots_" + ctr.suffix + "_r17618_vs_r17662.png";
            c.SaveAs(out.c_str());
            std::cout << "  Saved: " << out << "\n";

            for (auto* h : frames) delete h;
            for (auto* h : e618s) delete h;
            for (auto* h : e662s) delete h;
            for (auto* h : ratios) delete h;
        }
    }
    std::cout << "Done.\n";
}
