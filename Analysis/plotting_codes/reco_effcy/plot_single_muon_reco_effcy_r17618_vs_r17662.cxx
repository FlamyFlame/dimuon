// Overlay R17662 (signal-truth-only) on R17618 (collision sample, +dR fix) for the
// single-muon reconstruction efficiency, for the 0-5% and 5-10% centrality bins.
// Two plot types per centrality: q-eta-integrated, and q-eta-binned (3x3 subplots).
// Efficiency = pass_medium / all fiducial truth muons vs truth pT (weighted, binomial
// errors), modelled on plot_single_muon_reco_effcy.cxx. Sample comparison: the
// single-muon ε is expected to agree within stats in the turn-on region (R17662 is
// pTH8_14 only, so it has no plateau stats above ~15 GeV).
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TSystem.h>
#include <TStyle.h>
#include <ROOT/RDataFrame.hxx>
#include <string>
#include <vector>
#include <iostream>
#include "../../RDFBasedHistFilling/CommonEffcyConfig.h"

namespace {
    const std::vector<double> pt_edges = {4, 5, 6, 7, 8, 10, 12, 15, 20, 30, 50, 80, 120};
    const int nPtBins = (int)pt_edges.size() - 1;

    const std::string R17618 =
        "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_hijing_overlay_test_sample/"
        "muon_pairs_pythia_fullsim_hijing_overlay_pp24_no_data_resonance_cuts_single_muon.root";
    const std::string R17662 =
        "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_hijing_overlay_test_sample/r17662_run/"
        "muon_pairs_pythia_fullsim_hijing_overlay_pp24_no_data_resonance_cuts_r17662_withdr_single_muon.root";
    const std::string OUTDIR =
        "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_hijing_overlay_test_sample/plots/r17618_vs_r17662_comparison/";

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
        auto h_num   = dq.Filter("pass_medium").Histo1D({"h_num", "", nPtBins, pt_edges.data()}, "truth_pt", "ev_weight");

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

    TH1D* make_frame(const char* name, double titsize = 0.045, double labsize = 0.038) {
        TH1D* f = new TH1D(name, "", nPtBins, pt_edges.data());
        f->SetMinimum(0.0); f->SetMaximum(1.15);
        f->GetXaxis()->SetTitle("Truth p_{T} [GeV]");
        f->GetYaxis()->SetTitle("Single-muon reco efficiency (medium)");
        f->GetXaxis()->SetTitleSize(titsize); f->GetYaxis()->SetTitleSize(titsize);
        f->GetXaxis()->SetLabelSize(labsize); f->GetYaxis()->SetLabelSize(labsize);
        f->GetXaxis()->SetMoreLogLabels(true); f->GetXaxis()->SetNoExponent(true);
        return f;
    }
}

void plot_single_muon_reco_effcy_r17618_vs_r17662() {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gSystem->mkdir(OUTDIR.c_str(), kTRUE);

    static const CommonEffcyConfig ecfg;
    const QEtaBinning& qeta_bins = ecfg.q_eta_proj_ranges_coarse_incl_gap;
    const int nQEta = (int)qeta_bins.size();

    struct CtrBin { int lo, hi; std::string suffix, label; };
    std::vector<CtrBin> ctr_bins = { {0, 5, "ctr0_5", "0-5%"}, {5, 10, "ctr5_10", "5-10%"} };

    for (const auto& ctr : ctr_bins) {

        // ---- (A) q*eta-integrated overlay ----
        {
            TH1D* e618 = make_eff(R17618, ctr.lo, ctr.hi, true, 0, 0, ("eint618_" + ctr.suffix).c_str());
            TH1D* e662 = make_eff(R17662, ctr.lo, ctr.hi, true, 0, 0, ("eint662_" + ctr.suffix).c_str());
            style_eff(e618, kBlack,  20);
            style_eff(e662, kRed+1,  24);

            TCanvas c("c_int", "c_int", 900, 700);
            c.SetTicks(1, 1); c.SetLogx(true);
            TH1D* frame = make_frame(("fint_" + ctr.suffix).c_str());
            frame->Draw();
            e618->Draw("PE1 SAME");
            e662->Draw("PE1 SAME");

            TLegend leg(0.45, 0.16, 0.90, 0.30);
            leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.032);
            leg.AddEntry(e618, "R17618 (+dR fix)", "lpe");
            leg.AddEntry(e662, "R17662 (signal-only truth)", "lpe");
            leg.Draw("SAME");

            TLatex lat; lat.SetNDC();
            lat.SetTextSize(0.040); lat.DrawLatex(0.16, 0.92, "Pythia fullsim HIJING overlay");
            lat.SetTextSize(0.036); lat.DrawLatex(0.16, 0.87, "q#eta-integrated");
            lat.DrawLatex(0.16, 0.82, ("Centrality: " + ctr.label).c_str());

            std::string out = OUTDIR + "single_muon_reco_effcy_vs_pt_integrated_" + ctr.suffix + "_r17618_vs_r17662.png";
            c.SaveAs(out.c_str());
            std::cout << "  Saved: " << out << "\n";
            delete frame; delete e618; delete e662;
        }

        // ---- (B) q*eta-binned: one bin per subplot (3x3) ----
        {
            const int nrows = 3, ncols = 3;
            TCanvas c("c_sub", "c_sub", 1500, 1000);
            c.Divide(ncols, nrows, 0.002, 0.002);
            std::vector<TH1D*> frames, e618s, e662s;

            for (int iq = 0; iq < nQEta; ++iq) {
                float qlo = qeta_bins[iq].first, qhi = qeta_bins[iq].second;
                TH1D* e618 = make_eff(R17618, ctr.lo, ctr.hi, false, qlo, qhi,
                                      ("e618_" + ctr.suffix + "_q" + std::to_string(iq)).c_str());
                TH1D* e662 = make_eff(R17662, ctr.lo, ctr.hi, false, qlo, qhi,
                                      ("e662_" + ctr.suffix + "_q" + std::to_string(iq)).c_str());
                style_eff(e618, kBlack, 20); e618->SetMarkerSize(0.9);
                style_eff(e662, kRed+1, 24); e662->SetMarkerSize(0.9);
                e618s.push_back(e618); e662s.push_back(e662);

                c.cd(iq + 1);
                gPad->SetTicks(1, 1); gPad->SetLogx(true);
                gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.13);
                gPad->SetTopMargin(0.06); gPad->SetRightMargin(0.04);

                TH1D* hf = make_frame(("fsub_" + ctr.suffix + "_q" + std::to_string(iq)).c_str(), 0.055, 0.048);
                hf->GetYaxis()->SetTitle("#varepsilon_{reco} (medium)");
                hf->GetXaxis()->SetTitleOffset(1.1); hf->GetYaxis()->SetTitleOffset(1.2);
                hf->Draw();
                frames.push_back(hf);

                e618->Draw("PE1 SAME");
                e662->Draw("PE1 SAME");

                char buf[64];
                snprintf(buf, sizeof(buf), "%.1f < q#eta < %.1f", qlo, qhi);
                TLatex lab; lab.SetNDC(); lab.SetTextSize(0.055);
                lab.DrawLatex(0.20, 0.88, buf);

                if (iq == 0) {
                    TLegend* leg = new TLegend(0.30, 0.13, 0.97, 0.34);
                    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.050);
                    leg->AddEntry(e618, "R17618 (+dR fix)", "lpe");
                    leg->AddEntry(e662, "R17662 (signal-only)", "lpe");
                    leg->Draw("SAME");
                }
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
        }
    }
    std::cout << "Done.\n";
}
