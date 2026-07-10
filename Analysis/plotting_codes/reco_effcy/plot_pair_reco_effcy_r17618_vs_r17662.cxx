// Overlay R17662 (signal-only truth, no barcode collision) on R17618 (full HIJING
// truth, has the Pythia/HIJING barcode collision) for the PAIR reconstruction
// efficiency, 0-5% and 5-10% centrality. BOTH samples: pure prob>0.5 barcode matching
// (dR fallback deleted 2026-07-10), pTH8_14 (R17618 = kin0 only), so the pT-hat mix
// matches. Any R17618<R17662 deficit is a genuine pair-reco-efficiency UNDERESTIMATE
// from the barcode collision.
//
// Uses the OFFICIAL RDF reco-efficiency histograms (produced by
// RDFBasedHistFillingPythiaFullsimOverlay): efficiency = h_<var>_<cat>_<ctr>_pass_<WP>
// / h_<var>_<cat>_<ctr>, for cat in {op, single_b}, var in {truth_pair_pt,
// truth_dr_zoomin}. Binomial errors.
#include <TROOT.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TSystem.h>
#include <TStyle.h>
#include <string>
#include <vector>
#include <iostream>

namespace {
    const std::string DIR =
        "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_hijing_overlay_test_sample/";
    const std::string H618 = DIR +
        "histograms_pythia_fullsim_hijing_overlay_pbpb23_no_data_resonance_cuts_r17618kin0_pair.root";
    const std::string H662 = DIR +
        "histograms_pythia_fullsim_hijing_overlay_pbpb23_no_data_resonance_cuts_r17662_pair.root";
    std::string OUTDIR = DIR + "plots/r17618_vs_r17662_comparison/pair/";

    std::string wp = "pass_tight";   // NOMINAL Tight; entry fn can switch to medium
    std::string wp_label = "tight";

    // Build a binomial-error efficiency TGraph from num/denom hists in one file.
    // Suppresses bins with < ~10 effective denom entries (matches single-muon macro).
    TGraphAsymmErrors* make_eff(TFile* f, const std::string& var, const std::string& cat,
                                const std::string& ctr, const char* gname, int color, int mstyle) {
        std::string dname = "h_" + var + "_" + cat + "_" + ctr;
        std::string nname = dname + "_" + wp;
        TH1D* hd = (TH1D*)f->Get(dname.c_str());
        TH1D* hn = (TH1D*)f->Get(nname.c_str());
        if (!hd || !hn) {
            std::cerr << "  MISSING: " << dname << " or " << nname << " in " << f->GetName() << "\n";
            return nullptr;
        }
        hd = (TH1D*)hd->Clone(Form("%s_d", gname)); hd->SetDirectory(nullptr);
        hn = (TH1D*)hn->Clone(Form("%s_n", gname)); hn->SetDirectory(nullptr);
        // zero out low-stat denom bins so TGraphAsymmErrors::Divide is well-defined
        for (int ib = 1; ib <= hd->GetNbinsX(); ++ib) {
            double v = hd->GetBinContent(ib), e = hd->GetBinError(ib);
            bool bad = (e <= 0) || ((v/e)*(v/e) < 10.0) || (hn->GetBinContent(ib) > v);
            if (bad) { hd->SetBinContent(ib,0); hd->SetBinError(ib,0);
                       hn->SetBinContent(ib,0); hn->SetBinError(ib,0); }
        }
        TGraphAsymmErrors* g = new TGraphAsymmErrors();
        g->Divide(hn, hd, "cl=0.683 b(1,1) mode");
        g->SetName(gname);
        g->SetLineColor(color); g->SetMarkerColor(color);
        g->SetMarkerStyle(mstyle); g->SetMarkerSize(1.1); g->SetLineWidth(2);
        delete hd; delete hn;
        return g;
    }

    void one_plot(TFile* f618, TFile* f662, const std::string& var, const std::string& xtitle,
                  const std::string& cat, const std::string& catlabel,
                  const std::string& ctr, const std::string& ctrlabel, bool logx) {
        TGraphAsymmErrors* g618 = make_eff(f618, var, cat, ctr, "g618", kBlack, 20);
        TGraphAsymmErrors* g662 = make_eff(f662, var, cat, ctr, "g662", kRed+1, 24);
        if (!g618 || !g662) return;

        TCanvas c("c", "c", 900, 700);
        c.SetTicks(1,1); if (logx) c.SetLogx(true);
        // explicit x-range from the graphs' points
        double xmin = 1e9, xmax = -1e9;
        for (auto* g : {g618, g662})
            for (int i = 0; i < g->GetN(); ++i) { double x=g->GetPointX(i); xmin=std::min(xmin,x); xmax=std::max(xmax,x); }
        TH1D* frame = new TH1D("frame","",1, xmin*0.9, xmax*1.1);
        frame->SetMinimum(0.0); frame->SetMaximum(1.15);
        frame->GetXaxis()->SetTitle(xtitle.c_str());
        frame->GetYaxis()->SetTitle(("Pair reco efficiency (" + wp_label + ")").c_str());
        frame->GetXaxis()->SetTitleSize(0.045); frame->GetYaxis()->SetTitleSize(0.045);
        frame->Draw();
        g618->Draw("PZ SAME");
        g662->Draw("PZ SAME");

        TLegend leg(0.42, 0.16, 0.90, 0.30);
        leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.032);
        leg.AddEntry(g618, "R17618 (full HIJING truth)", "lpe");
        leg.AddEntry(g662, "R17662 (signal-only truth)", "lpe");
        leg.Draw("SAME");

        TLatex lat; lat.SetNDC();
        lat.SetTextSize(0.038); lat.DrawLatex(0.16, 0.92, "Pythia fullsim HIJING overlay (pTH8_14)");
        lat.SetTextSize(0.034); lat.DrawLatex(0.16, 0.87, (catlabel + ",  Centrality: " + ctrlabel).c_str());

        std::string out = OUTDIR + "pair_reco_effcy_" + var + "_" + cat + "_" + ctr + "_r17618_vs_r17662.png";
        c.SaveAs(out.c_str());
        std::cout << "  Saved: " << out << "\n";
        delete frame; delete g618; delete g662;
    }
}

void plot_pair_reco_effcy_r17618_vs_r17662(bool useTight = true) {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    wp = useTight ? "pass_tight" : "pass_medium";
    wp_label = useTight ? "tight" : "medium";
    if (!useTight) OUTDIR += "medium_wp/";
    gSystem->mkdir(OUTDIR.c_str(), kTRUE);

    TFile* f618 = TFile::Open(H618.c_str());
    TFile* f662 = TFile::Open(H662.c_str());
    if (!f618 || f618->IsZombie() || !f662 || f662->IsZombie()) {
        std::cerr << "Cannot open input hist files\n"; return;
    }

    struct Ctr { std::string tag, label; };
    std::vector<Ctr> ctrs = { {"ctr0_5","0-5%"}, {"ctr5_10","5-10%"} };
    struct Cat { std::string tag, label; };
    std::vector<Cat> cats = { {"op","OS pairs"}, {"single_b","single-b pairs"} };

    for (const auto& ctr : ctrs)
        for (const auto& cat : cats) {
            one_plot(f618, f662, "truth_pair_pt", "Truth pair p_{T} [GeV]",
                     cat.tag, cat.label, ctr.tag, ctr.label, false);
            one_plot(f618, f662, "truth_dr_zoomin", "Truth #DeltaR(#mu#mu)",
                     cat.tag, cat.label, ctr.tag, ctr.label, false);
        }
    f618->Close(); f662->Close();
    std::cout << "Done.\n";
}
