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
#include <TPad.h>
#include <TLine.h>
#include <cmath>
#include <string>
#include <vector>
#include <iostream>

namespace {
    const std::string DIR =
        "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_hijing_overlay_test_sample_pbpb23/";
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

    // Binomial-error efficiency as a TH1D (used for the ratio panel; same low-stat
    // suppression as make_eff).
    TH1D* eff_hist(TFile* f, const std::string& var, const std::string& cat,
                   const std::string& ctr, const char* name) {
        std::string dname = "h_" + var + "_" + cat + "_" + ctr;
        TH1D* hd = (TH1D*)f->Get(dname.c_str());
        TH1D* hn = (TH1D*)f->Get((dname + "_" + wp).c_str());
        if (!hd || !hn) return nullptr;
        hd = (TH1D*)hd->Clone(Form("%s_d", name)); hd->SetDirectory(nullptr);
        hn = (TH1D*)hn->Clone(Form("%s_n", name)); hn->SetDirectory(nullptr);
        TH1D* e = (TH1D*)hn->Clone(name); e->SetDirectory(nullptr);
        e->Divide(hn, hd, 1.0, 1.0, "B");
        for (int i = 1; i <= e->GetNbinsX(); ++i) {
            double v = hd->GetBinContent(i), er = hd->GetBinError(i);
            if (er <= 0 || (v/er)*(v/er) < 10.0) { e->SetBinContent(i,0); e->SetBinError(i,0); }
        }
        delete hd; delete hn;
        return e;
    }

    // R17618/R17662 ratio (independent errors: the two r-tags are separate digi+reco passes).
    TH1D* make_ratio(const TH1D* a, const TH1D* b, const char* name) {
        TH1D* r = (TH1D*)a->Clone(name); r->SetDirectory(nullptr);
        for (int i = 1; i <= r->GetNbinsX(); ++i) {
            double va=a->GetBinContent(i), vb=b->GetBinContent(i);
            double ea=a->GetBinError(i),   eb=b->GetBinError(i);
            if (va<=0 || vb<=0) { r->SetBinContent(i,0); r->SetBinError(i,0); continue; }
            double q=va/vb; r->SetBinContent(i,q);
            r->SetBinError(i, q*std::sqrt((ea/va)*(ea/va)+(eb/vb)*(eb/vb)));
        }
        r->SetLineColor(kBlack); r->SetMarkerColor(kBlack);
        r->SetMarkerStyle(20); r->SetMarkerSize(1.0); r->SetLineWidth(2);
        return r;
    }

    // Nudge a graph slightly in x so two exactly-agreeing series do not hide one another
    // (a fully overlapping marker reads as a missing point).
    void x_offset(TGraphAsymmErrors* g, double frac) {
        for (int i = 0; i < g->GetN(); ++i) {
            const double w = g->GetErrorXlow(i) + g->GetErrorXhigh(i);   // bin width
            g->SetPointX(i, g->GetPointX(i) + frac * w);
        }
    }

    void one_plot(TFile* f618, TFile* f662, const std::string& var, const std::string& xtitle,
                  const std::string& cat, const std::string& catlabel,
                  const std::string& ctr, const std::string& ctrlabel, bool logx) {
        TGraphAsymmErrors* g618 = make_eff(f618, var, cat, ctr, "g618", kBlack, 20);
        TGraphAsymmErrors* g662 = make_eff(f662, var, cat, ctr, "g662", kRed+1, 24);
        TH1D* e618 = eff_hist(f618, var, cat, ctr, "e618r");
        TH1D* e662 = eff_hist(f662, var, cat, ctr, "e662r");
        if (!g618 || !g662 || !e618 || !e662) return;
        TH1D* ratio = make_ratio(e618, e662, "ratio_pair");
        x_offset(g662, 0.10);   // display-only nudge; the ratio is computed from the histograms

        TCanvas c("c", "c", 900, 800);
        TPad* pTop = new TPad("pTop","",0.0,0.32,1.0,1.0);
        TPad* pBot = new TPad("pBot","",0.0,0.0,1.0,0.32);
        pTop->SetBottomMargin(0.02); pTop->SetTicks(1,1); if (logx) pTop->SetLogx(true);
        pBot->SetTopMargin(0.03); pBot->SetBottomMargin(0.32); pBot->SetTicks(1,1); if (logx) pBot->SetLogx(true);
        pTop->Draw(); pBot->Draw();

        // x-range from the surviving bins' EDGES (R1: data fills the axis, error bars not clipped)
        double xmin = 1e30, xmax = -1e30;
        for (auto* g : {g618, g662})
            for (int i = 0; i < g->GetN(); ++i) {
                xmin = std::min(xmin, g->GetPointX(i) - g->GetErrorXlow(i));
                xmax = std::max(xmax, g->GetPointX(i) + g->GetErrorXhigh(i));
            }
        if (xmax <= xmin) {   // no surviving points -> do not write a broken canvas
            std::cerr << "  SKIP (no bins survive the low-stat cut): " << var << " " << cat << " " << ctr << "\n";
            delete g618; delete g662; delete e618; delete e662; delete ratio; return;
        }
        double lo, hi;
        if (logx) { const double r = xmax/xmin; lo = xmin/std::pow(r,0.06); hi = xmax*std::pow(r,0.06); }
        else      { lo = xmin - 0.04*(xmax-xmin); hi = xmax + 0.04*(xmax-xmin); }

        pTop->cd();
        TH1D* frame = new TH1D("frame","",1, lo, hi);
        frame->SetMinimum(0.0); frame->SetMaximum(1.15);
        frame->GetYaxis()->SetTitle(("Pair reco efficiency (" + wp_label + ")").c_str());
        frame->GetYaxis()->SetTitleSize(0.055); frame->GetYaxis()->SetLabelSize(0.048);
        frame->GetXaxis()->SetLabelSize(0); frame->GetXaxis()->SetTitleSize(0);
        if (logx) { frame->GetXaxis()->SetMoreLogLabels(true); frame->GetXaxis()->SetNoExponent(true); }
        frame->Draw();
        g618->Draw("PZ SAME");
        g662->Draw("PZ SAME");

        TLegend leg(0.52, 0.68, 0.94, 0.86);   // upper-right: the data sits low (eff ~0.1-0.6)
        leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.042);
        leg.AddEntry(g618, "R17618 (full HIJING truth)", "lpe");
        leg.AddEntry(g662, "R17662 (signal-only truth)", "lpe");
        leg.Draw("SAME");

        TLatex lat; lat.SetNDC();
        lat.SetTextSize(0.048); lat.DrawLatex(0.14, 0.93, "Pythia fullsim HIJING overlay");
        lat.SetTextSize(0.042); lat.DrawLatex(0.14, 0.87, (catlabel + ",  Centrality: " + ctrlabel).c_str());

        pBot->cd();
        TH1D* rframe = new TH1D("rframe","",1, lo, hi);
        // Auto-range so no ratio point is ever clipped (R3); keep a readable minimum span.
        double rlo = 1.0, rhi = 1.0;
        for (int i = 1; i <= ratio->GetNbinsX(); ++i) {
            double v = ratio->GetBinContent(i), e = ratio->GetBinError(i);
            if (v <= 0) continue;
            rlo = std::min(rlo, v - e); rhi = std::max(rhi, v + e);
        }
        double pad = 0.10 * std::max(rhi - rlo, 0.30);
        rlo = std::min(rlo - pad, 0.85); rhi = std::max(rhi + pad, 1.15);
        rframe->SetMinimum(std::max(0.0, rlo)); rframe->SetMaximum(rhi);
        rframe->GetYaxis()->SetTitle("R17618 / R17662");
        rframe->GetYaxis()->SetNdivisions(505);
        rframe->GetXaxis()->SetTitle(xtitle.c_str());
        rframe->GetXaxis()->SetTitleSize(0.11); rframe->GetYaxis()->SetTitleSize(0.095);
        rframe->GetXaxis()->SetLabelSize(0.09); rframe->GetYaxis()->SetLabelSize(0.080);
        rframe->GetYaxis()->SetTitleOffset(0.52); rframe->GetXaxis()->SetTitleOffset(1.25);
        if (logx) { rframe->GetXaxis()->SetMoreLogLabels(true); rframe->GetXaxis()->SetNoExponent(true); }
        rframe->Draw();
        TLine* l1 = new TLine(lo, 1.0, hi, 1.0);
        l1->SetLineStyle(2); l1->SetLineColor(kRed+1); l1->Draw("SAME");
        ratio->Draw("PE1 SAME");

        std::string out = OUTDIR + "pair_reco_effcy_" + var + "_" + cat + "_" + ctr + "_r17618_vs_r17662.png";
        c.SaveAs(out.c_str());
        std::cout << "  Saved: " << out << "\n";
        delete frame; delete rframe; delete g618; delete g662; delete e618; delete e662; delete ratio;
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
    std::vector<Cat> cats = { {"op","opposite-sign pairs"}, {"single_b","single-b pairs"} };

    for (const auto& ctr : ctrs)
        for (const auto& cat : cats) {
            one_plot(f618, f662, "truth_pair_pt", "Truth pair p_{T} [GeV]",
                     cat.tag, cat.label, ctr.tag, ctr.label, true);   // log bins -> log axis (R2)
            one_plot(f618, f662, "truth_dr_zoomin", "Truth #DeltaR(#mu#mu)",
                     cat.tag, cat.label, ctr.tag, ctr.label, false);
        }
    f618->Close(); f662->Close();
    std::cout << "Done.\n";
}
