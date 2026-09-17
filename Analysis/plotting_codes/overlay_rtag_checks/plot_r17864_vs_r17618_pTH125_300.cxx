// r-tag sanity comparison, HIJING overlay Pb+Pb 2024 conditions (r17864) vs Pb+Pb 2023
// conditions (r17618), SAME pT-hat slice 125-300 GeV (the only slice the r17864 test sample
// has), both processed by the SAME ntuple-processing code on the same day:
//   r17864 : pythia_fullsim_hijing_overlay_test_sample/            (OVERLAY_YEAR=24 run scripts)
//   r17618 : ..._test_sample_pbpb23/r17618_pTH125_300_run/          (single-slice override run,
//            NTupleProcessingCode/run_pythia_fullsim_overlay_pbpb23_pTH125_300_slice.sh)
// Tracking doc: docs/tracking/hijing_overlay_pbpb24_test_sample_skim.md, verification round.
//
// Every input is ntuple-processing OUTPUT (never a raw NTUP): the truth-seeded single-muon
// tree (`_single_muon`), its trigger-propagated twin (`_mc_trig_single_muon`) and the pair
// trees. Selections MIRROR the nominal chain:
//   recoeff  : plot_single_muon_reco_effcy_r17618_vs_r17662.cxx  -- eps = pass_<WP> / all truth
//              muons vs truth pT, per centrality bin (ev_centrality), weight ev_weight.
//   trigeff  : FillMCTrigEffHists.cxx Step 1 (mc_trigger_efficiency.md par. 3.1) -- denominator =
//              WP && pt > 4.5 && |eta| < 2.4 && truth pT > 4.5 && |truth eta| < 2.4 && detector
//              gap cut (ParamsSet::FiducialGapCutExpr) && 0-5 % centrality (doc D2); numerator =
//              passmu4 (the muon's own full mu4 chain); binning pT_bins_8 + pT_bins_60.
//   detresp  : (pT_reco - pT_truth)/pT_truth of WP muons, same single-muon selection as trigeff.
//   samb     : from_same_b fraction of opposite-sign reco pairs, pair_pass_<WP>, table only.
// CENTRALITY: the nominal 0-5 % (and 5-10 %) FCal bins are kept on purpose (user decision
// 2026-09-16) although the r17864 FCal is broken -- only 4.2 % of its events land in 0-5 %.
// Both samples are single-slice, so ev_weight is one constant per sample; it is applied anyway
// (it is what the nominal chain applies).
//
//   root -l -b -q 'plot_r17864_vs_r17618_pTH125_300.cxx+("recoeff")'
//   root -l -b -q 'plot_r17864_vs_r17618_pTH125_300.cxx+("trigeff")'
//   root -l -b -q 'plot_r17864_vs_r17618_pTH125_300.cxx+("detresp")'
//   root -l -b -q 'plot_r17864_vs_r17618_pTH125_300.cxx+("samb")'
// Second argument: useTight (default true = nominal Tight WP; false = Medium, own subdir).
#include <TFile.h>
#include <TH1D.h>
#include <TEfficiency.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TROOT.h>
#include <ROOT/RDataFrame.hxx>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>
#include "../../MuonObjectsParamsAndHelpers/ParamsSet.h"

namespace {
    const std::string D24 = "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_hijing_overlay_test_sample/";
    const std::string D23 = "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_hijing_overlay_test_sample_pbpb23/r17618_pTH125_300_run/";
    const std::string B24 = D24 + "muon_pairs_pythia_fullsim_hijing_overlay_pbpb24_no_data_resonance_cuts";
    const std::string B23 = D23 + "muon_pairs_pythia_fullsim_hijing_overlay_pbpb23_no_data_resonance_cuts";
    // ntuple-processing products (suffix order = PythiaAlgCoreT: ..._mc_trig + extra_output_suffix)
    const std::string SINGLE_24  = B24 + "_single_muon.root";
    const std::string SINGLE_23  = B23 + "_r17618pTH125_300_single_muon.root";
    const std::string TRIG_24    = B24 + "_mc_trig_single_muon.root";
    const std::string TRIG_23    = B23 + "_mc_trig_r17618pTH125_300_single_muon.root";
    const std::string PAIR_24    = B24 + ".root";
    const std::string PAIR_23    = B23 + "_r17618pTH125_300.root";
    std::string OUTDIR = D24 + "plots/r17864_rtag_sanity/";

    const std::string L24 = "Pb+Pb 2024 conditions (r17864)";
    const std::string L23 = "Pb+Pb 2023 conditions (r17618)";
    const int C24 = kRed + 1, C23 = kBlue + 1;

    std::string wp_col = "pass_tight", wp_label = "Tight";

    struct CtrBin { int lo, hi; std::string suffix, label; };
    // nominal FCal-centrality bins of the overlay chain (plot_single_muon_reco_effcy.cxx first two)
    const std::vector<CtrBin> ctr_bins = { {0, 5, "ctr0_5", "0-5%"}, {5, 10, "ctr5_10", "5-10%"} };

    // reco-efficiency pT axis = the canonical coarse single-muon binning (ParamsSet, never retyped):
    // the r17864 0-5 % bin holds only ~750 truth muons, so the coarse axis is the readable one.
    std::vector<double> coarse_pt_edges() { ParamsSet pms; return pms.single_mu_pt_coarse_bins; }
    // Step-1 trigger-efficiency binning, read from ParamsSet exactly as FillMCTrigEffHists does
    std::vector<double> trig_pt_edges() {
        ParamsSet pms;
        std::vector<double> e = pms.pT_bins_8;
        e.insert(e.end(), pms.pT_bins_60.begin(), pms.pT_bins_60.end());
        // the data convention duplicates the 8.0 edge (one zero-width bin); drop it for drawing
        std::vector<double> out;
        for (double x : e) if (out.empty() || x > out.back()) out.push_back(x);
        return out;
    }

    void style_pts(TH1D* h, int color, int mstyle) {
        h->SetLineColor(color); h->SetMarkerColor(color);
        h->SetMarkerStyle(mstyle); h->SetMarkerSize(1.0); h->SetLineWidth(2);
    }
    // Efficiency with Clopper-Pearson 68.3 % intervals (TEfficiency::ClopperPearson on the
    // EFFECTIVE entries N_eff = (sum w)^2 / sum w^2, k_eff = eps * N_eff), symmetrised as the larger
    // of the two half-widths. TH1::Divide(...,"B") gives a ZERO error for a bin at eps = 1 (21/21
    // muons at 11-12 GeV in the r17864 0-5 % sample), which made its ratio read as ~8 sigma; the
    // CP interval is the honest small-N error. Bins whose denominator has < 10 effective entries
    // are blanked (same floor as the reference macro).
    TH1D* make_eff(ROOT::RDF::RNode den, const std::string& num_sel, const std::string& var,
                   const std::vector<double>& edges, const char* name) {
        const int nb = (int)edges.size() - 1;
        auto hd = den.Histo1D({Form("%s_d", name), "", nb, edges.data()}, var, "ev_weight");
        auto hn = den.Filter(num_sel).Histo1D({Form("%s_n", name), "", nb, edges.data()}, var, "ev_weight");
        TH1D* eff = (TH1D*)hn->Clone(name); eff->SetDirectory(nullptr); eff->Reset();
        for (int ib = 1; ib <= nb; ++ib) {
            const double d = hd->GetBinContent(ib), de = hd->GetBinError(ib), n = hn->GetBinContent(ib);
            const double neff = (de > 0) ? (d * d) / (de * de) : 0.0;
            if (d <= 0 || neff < 10.0) { eff->SetBinContent(ib, 0); eff->SetBinError(ib, 0); continue; }
            const double e = n / d;
            const int N = (int)std::lround(neff), k = (int)std::lround(e * neff);
            const double up = TEfficiency::ClopperPearson(N, k, 0.683, true);
            const double lo = TEfficiency::ClopperPearson(N, k, 0.683, false);
            eff->SetBinContent(ib, e);
            eff->SetBinError(ib, std::max(up - e, e - lo));
        }
        return eff;
    }
    // ratio a/b with independent errors (two independent detector realisations); blanked bins dropped
    TH1D* make_ratio(const TH1D* a, const TH1D* b, const char* name) {
        TH1D* r = (TH1D*)a->Clone(name); r->SetDirectory(nullptr);
        for (int i = 1; i <= r->GetNbinsX(); ++i) {
            const double va = a->GetBinContent(i), vb = b->GetBinContent(i);
            const double ea = a->GetBinError(i),   eb = b->GetBinError(i);
            if (va <= 0 || vb <= 0) { r->SetBinContent(i, 0); r->SetBinError(i, 0); continue; }
            const double q = va / vb;
            r->SetBinContent(i, q); r->SetBinError(i, q * std::sqrt((ea/va)*(ea/va) + (eb/vb)*(eb/vb)));
        }
        style_pts(r, kBlack, 20);
        return r;
    }
    void set_ratio_range(TH1D* f, const TH1D* r) {
        double lo = 1.0, hi = 1.0;
        for (int i = 1; i <= r->GetNbinsX(); ++i) {
            const double v = r->GetBinContent(i), e = r->GetBinError(i);
            if (v <= 0) continue;
            lo = std::min(lo, v - e); hi = std::max(hi, v + e);
        }
        const double pad = 0.10 * std::max(hi - lo, 0.30);
        f->SetMinimum(std::max(0.0, std::min(lo - pad, 0.85)));
        f->SetMaximum(std::max(hi + pad, 1.15));
    }
    // efficiency-vs-pT canvas with ratio pad; returns after saving
    void draw_eff_canvas(TH1D* e24, TH1D* e23, const std::vector<double>& edges,
                         const std::string& ytitle, const std::string& xtitle,
                         const std::vector<std::string>& header_lines, const std::string& out) {
        TH1D* ratio = make_ratio(e24, e23, Form("%s_ratio", e24->GetName()));
        TCanvas c("c", "", 900, 800);
        TPad* pTop = new TPad("pTop", "", 0.0, 0.32, 1.0, 1.0);
        TPad* pBot = new TPad("pBot", "", 0.0, 0.0,  1.0, 0.32);
        pTop->SetBottomMargin(0.02); pTop->SetTicks(1,1); pTop->SetLogx(true);
        pBot->SetTopMargin(0.03); pBot->SetBottomMargin(0.32); pBot->SetTicks(1,1); pBot->SetLogx(true);
        pTop->Draw(); pBot->Draw();
        pTop->cd();
        TH1D* frame = new TH1D(Form("%s_frame", e24->GetName()), "", (int)edges.size() - 1, edges.data());
        frame->SetMinimum(0.0); frame->SetMaximum(1.25);
        frame->GetYaxis()->SetTitle(ytitle.c_str());
        frame->GetYaxis()->SetTitleSize(0.045); frame->GetYaxis()->SetLabelSize(0.038);
        frame->GetXaxis()->SetLabelSize(0); frame->GetXaxis()->SetTitleSize(0);
        frame->Draw();
        e24->Draw("PE1 SAME"); e23->Draw("PE1 SAME");
        TLegend leg(0.42, 0.14, 0.92, 0.30);
        leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.040);
        leg.AddEntry(e24, L24.c_str(), "lpe");
        leg.AddEntry(e23, L23.c_str(), "lpe");
        leg.Draw("SAME");
        TLatex lat; lat.SetNDC(); lat.SetTextSize(0.042);
        double y = 0.92;
        for (const auto& s : header_lines) { lat.DrawLatex(0.16, y, s.c_str()); y -= 0.06; }
        pBot->cd();
        TH1D* rframe = new TH1D(Form("%s_rframe", e24->GetName()), "", (int)edges.size() - 1, edges.data());
        rframe->GetYaxis()->SetTitle("2024 / 2023 conditions");
        rframe->GetYaxis()->SetNdivisions(505);
        rframe->GetXaxis()->SetTitle(xtitle.c_str());
        rframe->GetXaxis()->SetTitleSize(0.095); rframe->GetYaxis()->SetTitleSize(0.095);
        rframe->GetXaxis()->SetLabelSize(0.080); rframe->GetYaxis()->SetLabelSize(0.080);
        rframe->GetXaxis()->SetMoreLogLabels(true); rframe->GetXaxis()->SetNoExponent(true);
        rframe->GetYaxis()->SetTitleOffset(0.52); rframe->GetXaxis()->SetTitleOffset(1.25);
        set_ratio_range(rframe, ratio);
        rframe->Draw();
        TLine* l1 = new TLine(edges.front(), 1.0, edges.back(), 1.0);
        l1->SetLineStyle(2); l1->SetLineColor(kGray + 2); l1->Draw("SAME");
        ratio->Draw("PE1 SAME");
        c.SaveAs(out.c_str());
        std::printf("  Saved: %s\n", out.c_str());
        delete frame; delete rframe; delete ratio;
    }
    void style() { gROOT->SetBatch(kTRUE); gStyle->SetOptStat(0); gStyle->SetOptTitle(0); gSystem->mkdir(OUTDIR.c_str(), kTRUE); }
}

// (1) single-muon reconstruction efficiency vs truth pT, nominal centrality bins, with ratio
void plot_recoeff() {
    for (const auto& ctr : ctr_bins) {
        ROOT::RDataFrame d24("muon_tree", SINGLE_24), d23("muon_tree", SINGLE_23);
        auto sel = [&](ROOT::RDataFrame& d) {
            return d.Filter(Form("ev_centrality >= %d && ev_centrality < %d", ctr.lo, ctr.hi));
        };
        auto n24 = sel(d24), n23 = sel(d23);
        const std::vector<double> reco_pt_edges = coarse_pt_edges();
        TH1D* e24 = make_eff(n24, wp_col, "truth_pt", reco_pt_edges, ("reco24_" + ctr.suffix).c_str());
        TH1D* e23 = make_eff(n23, wp_col, "truth_pt", reco_pt_edges, ("reco23_" + ctr.suffix).c_str());
        style_pts(e24, C24, 20); style_pts(e23, C23, 24);
        const double N24 = *n24.Count(), N23 = *n23.Count();
        std::printf("recoeff %s: truth muons in bin: r17864 %.0f, r17618 %.0f\n", ctr.label.c_str(), N24, N23);
        draw_eff_canvas(e24, e23, reco_pt_edges,
                        "single-muon reconstruction efficiency (" + wp_label + ")", "truth p_{T} [GeV]",
                        {"HIJING overlay, #hat{p}_{T} 125-300 GeV", "FCal centrality " + ctr.label + ", q#eta-integrated"},
                        OUTDIR + "single_muon_reco_effcy_vs_pt_" + ctr.suffix + "_r17864_vs_r17618.png");
        delete e24; delete e23;
    }
}

// (2) MC single-muon mu4 efficiency vs pT, Step-1 mirror, 0-5 % centrality, with ratio
void plot_trigeff() {
    const std::vector<double> edges = trig_pt_edges();
    const std::string gap = ParamsSet::FiducialGapCutExpr("charge * eta");
    const std::string sel_single = wp_col + " && pt > 4.5 && fabs(eta) < 2.4"
                                 + " && truth_pt > 4.5 && fabs(truth_eta) < 2.4 && " + gap
                                 + " && ev_centrality >= 0 && ev_centrality < 5";
    ROOT::RDataFrame d24("muon_tree", TRIG_24), d23("muon_tree", TRIG_23);
    auto n24 = d24.Filter(sel_single), n23 = d23.Filter(sel_single);
    TH1D* e24 = make_eff(n24, "passmu4", "pt", edges, "trig24");
    TH1D* e23 = make_eff(n23, "passmu4", "pt", edges, "trig23");
    style_pts(e24, C24, 20); style_pts(e23, C23, 24);
    const double N24 = *n24.Count(), N23 = *n23.Count();
    const double P24 = *n24.Filter("passmu4").Count(), P23 = *n23.Filter("passmu4").Count();
    std::printf("trigeff 0-5%%: denominator muons r17864 %.0f (mu4 %.0f, %.3f), r17618 %.0f (mu4 %.0f, %.3f)\n",
                N24, P24, P24 / N24, N23, P23, P23 / N23);
    draw_eff_canvas(e24, e23, edges,
                    "single-muon mu4 efficiency (" + wp_label + ")", "p_{T} [GeV]",
                    {"HIJING overlay, #hat{p}_{T} 125-300 GeV", "FCal centrality 0-5%, q#eta-integrated"},
                    OUTDIR + "single_muon_mu4_effcy_vs_pt_ctr0_5_r17864_vs_r17618.png");
    delete e24; delete e23;
    // same efficiency on the canonical COARSE single-muon axis (ParamsSet::single_mu_pt_coarse_bins):
    // with ~420 r17864 muons in 0-5 % the Step-1 axis leaves most bins below the 10-entry floor.
    const std::vector<double> cedges = coarse_pt_edges();
    TH1D* c24 = make_eff(n24, "passmu4", "pt", cedges, "trig24_coarse");
    TH1D* c23 = make_eff(n23, "passmu4", "pt", cedges, "trig23_coarse");
    style_pts(c24, C24, 20); style_pts(c23, C23, 24);
    draw_eff_canvas(c24, c23, cedges,
                    "single-muon mu4 efficiency (" + wp_label + ")", "p_{T} [GeV]",
                    {"HIJING overlay, #hat{p}_{T} 125-300 GeV", "FCal centrality 0-5%, q#eta-integrated"},
                    OUTDIR + "single_muon_mu4_effcy_vs_pt_coarse_ctr0_5_r17864_vs_r17618.png");
    delete c24; delete c23;
}

// (3) detector response: pT residual of WP muons, same single-muon selection as (2)
void plot_detresp() {
    const std::string gap = ParamsSet::FiducialGapCutExpr("charge * eta");
    const std::string sel_single = wp_col + " && pt > 4.5 && fabs(eta) < 2.4"
                                 + " && truth_pt > 4.5 && fabs(truth_eta) < 2.4 && " + gap
                                 + " && ev_centrality >= 0 && ev_centrality < 5";
    ROOT::RDataFrame d24("muon_tree", TRIG_24), d23("muon_tree", TRIG_23);
    auto def = [&](ROOT::RDataFrame& d) {
        return d.Filter(sel_single).Define("dpt_rel", "(pt - truth_pt) / truth_pt")
                                   .Define("deta", "eta - truth_eta");
    };
    auto n24 = def(d24), n23 = def(d23);
    auto h24 = n24.Histo1D({"dpt24", "", 60, -0.15, 0.15}, "dpt_rel", "ev_weight");
    auto h23 = n23.Histo1D({"dpt23", "", 60, -0.15, 0.15}, "dpt_rel", "ev_weight");
    auto g24 = n24.Histo1D({"deta24", "", 60, -0.006, 0.006}, "deta", "ev_weight");
    auto g23 = n23.Histo1D({"deta23", "", 60, -0.006, 0.006}, "deta", "ev_weight");
    const double N24 = *n24.Count(), N23 = *n23.Count();
    std::printf("detresp 0-5%%: muons r17864 %.0f  dpt/pt mean %+.4f rms %.4f | r17618 %.0f  mean %+.4f rms %.4f\n",
                N24, h24->GetMean(), h24->GetRMS(), N23, h23->GetMean(), h23->GetRMS());
    TCanvas c("c", "", 1600, 700); c.Divide(2, 1);
    // `scale` = unit of the printed mean/RMS (1 for pT, 1e-3 for eta), `unit` its label
    auto panel = [&](int ipad, TH1* a, TH1* b, const char* xt, double bw_unit, double scale, const char* unit) {
        c.cd(ipad); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13); gPad->SetTicks(1, 1);
        a->Scale(1.0 / a->Integral()); b->Scale(1.0 / b->Integral());
        a->SetLineColor(C24); a->SetLineWidth(2); b->SetLineColor(C23); b->SetLineWidth(2);
        a->GetXaxis()->SetTitle(xt); a->GetYaxis()->SetTitle(Form("fraction of muons / %g", bw_unit));
        a->GetXaxis()->SetTitleSize(0.05); a->GetYaxis()->SetTitleSize(0.05);
        a->GetXaxis()->SetLabelSize(0.045); a->GetYaxis()->SetLabelSize(0.045);
        a->GetYaxis()->SetTitleOffset(1.3);
        a->SetMaximum(1.9 * std::max(a->GetBinContent(a->GetMaximumBin()), b->GetBinContent(b->GetMaximumBin())));
        a->SetMinimum(0);
        a->Draw("hist"); b->Draw("hist same");
        TLegend* l = new TLegend(0.16, 0.56, 0.88, 0.90); l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(0.031);
        l->SetHeader(("HIJING overlay, #hat{p}_{T} 125-300 GeV, FCal 0-5%, " + wp_label).c_str());
        l->AddEntry(a, L24.c_str(), "l");
        l->AddEntry((TObject*)nullptr, Form("mean %+.2f%s, RMS %.2f%s", a->GetMean() / scale, unit, a->GetRMS() / scale, unit), "");
        l->AddEntry(b, L23.c_str(), "l");
        l->AddEntry((TObject*)nullptr, Form("mean %+.2f%s, RMS %.2f%s", b->GetMean() / scale, unit, b->GetRMS() / scale, unit), "");
        l->Draw();
    };
    panel(1, h24.GetPtr(), h23.GetPtr(), "(p_{T}^{reco} - p_{T}^{truth}) / p_{T}^{truth}", 0.005, 1e-2, " %");
    panel(2, g24.GetPtr(), g23.GetPtr(), "#eta^{reco} - #eta^{truth}", 0.0002, 1e-3, "#times10^{-3}");
    const std::string out = OUTDIR + "single_muon_pt_eta_response_ctr0_5_r17864_vs_r17618.png";
    c.SaveAs(out.c_str());
    std::printf("  Saved: %s\n", out.c_str());
}

// (4) from_same_b fraction of opposite-sign reco pairs -- table + CSV, no plot
void print_samb() {
    const std::string pair_wp = (wp_col == "pass_tight") ? "pair_pass_tight" : "pair_pass_medium";
    struct Row { std::string sel_label, sel; };
    const std::vector<Row> rows = {
        {"all reco OS pairs",                          "true"},
        {"OS, " + wp_label + " pair, |y|<2.4, pT>4.5 both", pair_wp + " && m1.pt > 4.5 && m2.pt > 4.5 && fabs(m1.eta) < 2.4 && fabs(m2.eta) < 2.4"},
        {"OS, " + wp_label + " pair, 0-5% centrality",       pair_wp + " && m1.pt > 4.5 && m2.pt > 4.5 && fabs(m1.eta) < 2.4 && fabs(m2.eta) < 2.4 && avg_centrality >= 0 && avg_centrality < 5"},
    };
    std::ofstream csv(OUTDIR + "from_same_b_fraction_r17864_vs_r17618.csv");
    csv << "selection,sample,N_pairs,from_same_b_fraction,from_same_ancestors_fraction\n";
    std::printf("%-46s | %-8s | %7s | %12s | %16s\n", "selection (opposite-sign reco pairs)", "sample", "N", "from_same_b", "from_same_anc.");
    for (const auto& r : rows) {
        for (int s = 0; s < 2; ++s) {
            ROOT::RDataFrame d("muon_pair_tree_sign2", s == 0 ? PAIR_24 : PAIR_23);   // sign2 = opposite sign
            auto n = d.Filter(r.sel);
            const double N = *n.Count();
            const double fb = N > 0 ? *n.Filter("from_same_b").Count() / N : 0;
            const double fa = N > 0 ? *n.Filter("from_same_ancestors").Count() / N : 0;
            const double eb = N > 0 ? std::sqrt(fb * (1 - fb) / N) : 0;
            std::printf("%-46s | %-8s | %7.0f | %6.4f+-%.4f | %16.4f\n", r.sel_label.c_str(), s == 0 ? "r17864" : "r17618", N, fb, eb, fa);
            csv << r.sel_label << "," << (s == 0 ? "r17864" : "r17618") << "," << N << "," << fb << "," << fa << "\n";
        }
    }
    std::printf("  CSV: %sfrom_same_b_fraction_r17864_vs_r17618.csv\n", OUTDIR.c_str());
}

void plot_r17864_vs_r17618_pTH125_300(const char* what = "recoeff", bool useTight = true) {
    wp_col = useTight ? "pass_tight" : "pass_medium";
    wp_label = useTight ? "Tight" : "Medium";
    if (!useTight) OUTDIR += "medium_wp/";
    style();
    const std::string w = what;
    if      (w == "recoeff") plot_recoeff();
    else if (w == "trigeff") plot_trigeff();
    else if (w == "detresp") plot_detresp();
    else if (w == "samb")    print_samb();
    else std::printf("unknown mode '%s' (recoeff|trigeff|detresp|samb)\n", what);
}
