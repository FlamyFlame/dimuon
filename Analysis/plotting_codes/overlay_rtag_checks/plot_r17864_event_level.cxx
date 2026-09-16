// Event-level sanity plots for the HIJING-overlay Pb+Pb-2024-conditions test sample
// (r17864, pTH125_300 only, vtx z = -71.2 mm) against the Pb+Pb-2023-conditions test
// sample (r17618), SAME pT-hat slice so the two event populations match.
// Tracking doc: docs/tracking/hijing_overlay_pbpb24_test_sample_skim.md (R1).
//
// Reads the raw skim NTUPs (HeavyIonD3PD) directly: FCal_Et* and vtx_z / truth_vtx_z are
// per-event skim quantities with no selection procedure behind them, so nothing from the
// ntuple-processing chain is needed or bypassed.
//
//   root -l -b -q 'plot_r17864_event_level.cxx+("fcal")'   FCal sum-ET, both samples
//   root -l -b -q 'plot_r17864_event_level.cxx+("vtx")'    reco vs truth vertex z
//
// Every event is unweighted: both samples are ONE pT-hat slice (no slice mixing), and the
// quantities are underlying-event / beam-spot properties independent of the hard scatter.
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TLine.h>
#include <ROOT/RDataFrame.hxx>
#include <string>
#include <cstdio>

namespace {
    const std::string NEW_LABEL = "2024 conditions (r17864)";
    const std::string OLD_LABEL = "2023 conditions (r17618)";
    const std::string NEW_FILE =
        "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_hijing_overlay_test_sample/"
        "Pythia_5p36TeV_pp_hQCD_DiMu_pTH125_300.FullSimHIJINGOverlayPP24.NTUP.root";
    const std::string OLD_FILE =
        "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_hijing_overlay_test_sample_pbpb23/"
        "Pythia_5p36TeV_pp_hQCD_DiMu_pTH125_300.FullSimHIJINGOverlayPP24.NTUP.root";
    const std::string OUTDIR =
        "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_hijing_overlay_test_sample/plots/r17864_rtag_sanity/";
    const int kNew = kRed + 1, kOld = kBlue + 1;

    void style() {
        gStyle->SetOptStat(0);
        gStyle->SetOptTitle(0);
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);
        gSystem->mkdir(OUTDIR.c_str(), true);
    }
    void axes(TH1* h, const char* xt, const char* yt) {
        h->GetXaxis()->SetTitle(xt); h->GetYaxis()->SetTitle(yt);
        h->GetXaxis()->SetTitleSize(0.05); h->GetYaxis()->SetTitleSize(0.05);
        h->GetXaxis()->SetLabelSize(0.045); h->GetYaxis()->SetLabelSize(0.045);
        h->GetYaxis()->SetTitleOffset(1.3);
    }
    // FCal_Et is stored in MeV; TeV on the plots.
    ROOT::RDF::RNode fcal_node(ROOT::RDataFrame& df) {
        return df.Define("fcal_tev",   "FCal_Et   / 1e6f")
                 .Define("fcal_a_tev", "FCal_Et_P / 1e6f")   // A side = eta > 0
                 .Define("fcal_c_tev", "FCal_Et_N / 1e6f");  // C side = eta < 0
    }
}

void plot_fcal() {
    style();
    ROOT::RDataFrame dnew("HeavyIonD3PD", NEW_FILE), dold("HeavyIonD3PD", OLD_FILE);
    auto nn = fcal_node(dnew), no = fcal_node(dold);
    const ROOT::RDF::TH1DModel mtot("", "", 110, -3.0, 8.0), mside("", "", 90, -2.0, 4.0);
    auto h_new = nn.Histo1D(mtot, "fcal_tev"),    h_old = no.Histo1D(mtot, "fcal_tev");
    auto a_new = nn.Histo1D(mside, "fcal_a_tev"), c_new = nn.Histo1D(mside, "fcal_c_tev");
    auto a_old = no.Histo1D(mside, "fcal_a_tev"), c_old = no.Histo1D(mside, "fcal_c_tev");
    const double n_new = *nn.Count(), n_old = *no.Count();
    const double neg_new = *nn.Filter("fcal_tev < 0").Count(), neg_old = *no.Filter("fcal_tev < 0").Count();

    TCanvas c("c", "", 1600, 700);
    c.Divide(2, 1);
    // (a) total FCal sum-ET, per-event fraction so the two 10k samples overlay directly
    c.cd(1); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
    h_new->Scale(1.0 / n_new); h_old->Scale(1.0 / n_old);
    axes(h_new.GetPtr(), "FCal #Sigma E_{T} [TeV]", "fraction of events / 0.1 TeV");
    h_new->SetLineColor(kNew); h_new->SetLineWidth(2);
    h_old->SetLineColor(kOld); h_old->SetLineWidth(2);
    const double ymax = 1.75 * std::max(h_new->GetBinContent(h_new->GetMaximumBin()),
                                        h_old->GetBinContent(h_old->GetMaximumBin()));
    h_new->SetMaximum(ymax); h_new->SetMinimum(0);
    h_new->Draw("hist"); h_old->Draw("hist same");
    TLine zero(0, 0, 0, ymax); zero.SetLineStyle(2); zero.SetLineColor(kGray + 2); zero.Draw();
    TLegend l1(0.32, 0.62, 0.86, 0.89); l1.SetBorderSize(0); l1.SetFillStyle(0); l1.SetTextSize(0.031);
    l1.SetHeader("HIJING overlay, b < 5 fm, #hat{p}_{T} 125-300 GeV");
    l1.AddEntry(h_new.GetPtr(), Form("%s", NEW_LABEL.c_str()), "l");
    l1.AddEntry((TObject*)nullptr, Form("median %.2f TeV, %.0f%% negative", [&]{ double q, p = 0.5; h_new->GetQuantiles(1, &q, &p); return q; }(), 100.0 * neg_new / n_new), "");
    l1.AddEntry(h_old.GetPtr(), Form("%s", OLD_LABEL.c_str()), "l");
    l1.AddEntry((TObject*)nullptr, Form("median %.2f TeV, %.0f%% negative", [&]{ double q, p = 0.5; h_old->GetQuantiles(1, &q, &p); return q; }(), 100.0 * neg_old / n_old), "");
    l1.Draw();
    // (b) per side, A (eta>0) solid vs C (eta<0) dashed -- a timing/readout problem need not be symmetric
    c.cd(2); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
    for (auto* h : {a_new.GetPtr(), c_new.GetPtr()}) { h->Scale(1.0 / n_new); h->SetLineColor(kNew); h->SetLineWidth(2); }
    for (auto* h : {a_old.GetPtr(), c_old.GetPtr()}) { h->Scale(1.0 / n_old); h->SetLineColor(kOld); h->SetLineWidth(2); }
    c_new->SetLineStyle(2); c_old->SetLineStyle(2);
    axes(a_new.GetPtr(), "FCal #Sigma E_{T} per side [TeV]", "fraction of events / 0.067 TeV");
    double ymax2 = 0; for (auto* h : {a_new.GetPtr(), c_new.GetPtr(), a_old.GetPtr(), c_old.GetPtr()}) ymax2 = std::max(ymax2, h->GetBinContent(h->GetMaximumBin()));
    a_new->SetMaximum(1.7 * ymax2); a_new->SetMinimum(0);
    a_new->Draw("hist"); c_new->Draw("hist same"); a_old->Draw("hist same"); c_old->Draw("hist same");
    TLine zero2(0, 0, 0, 1.7 * ymax2); zero2.SetLineStyle(2); zero2.SetLineColor(kGray + 2); zero2.Draw();
    TLegend l2(0.30, 0.64, 0.89, 0.89); l2.SetBorderSize(0); l2.SetFillStyle(0); l2.SetTextSize(0.032);
    l2.SetHeader("side A: #eta > 0 (solid), side C: #eta < 0 (dashed)");
    l2.AddEntry(a_new.GetPtr(), "2024 conditions (r17864), side A", "l");
    l2.AddEntry(c_new.GetPtr(), "2024 conditions (r17864), side C", "l");
    l2.AddEntry(a_old.GetPtr(), "2023 conditions (r17618), side A", "l");
    l2.AddEntry(c_old.GetPtr(), "2023 conditions (r17618), side C", "l");
    l2.Draw();
    c.SaveAs((OUTDIR + "fcal_sum_et_r17864_vs_r17618_pTH125_300.png").c_str());
    std::printf("FCal sum-ET: new %.0f events (%.1f%% < 0), old %.0f events (%.1f%% < 0)\n",
                n_new, 100.0 * neg_new / n_new, n_old, 100.0 * neg_old / n_old);
}

void plot_vtx() {
    style();
    ROOT::RDataFrame dnew("HeavyIonD3PD", NEW_FILE), dold("HeavyIonD3PD", OLD_FILE);
    // vtx_z[0] = the reconstructed primary vertex; truth_vtx_z[0] = the generated interaction
    // point (every stored truth vertex sits at the signal PV, checked 2026-09-16). The generated
    // vertex is UNSMEARED: truth z is exactly -71.2000 mm in every r17864 event (the 0.01 mm in
    // r17864.txt is the reconstruction beam-spot width, not a truth smearing), so a 2D
    // reco-vs-truth view carries no information beyond the residual.
    auto def = [](ROOT::RDataFrame& d) {
        return d.Filter("vtx_z.size() > 0 && truth_vtx_z.size() > 0")
                .Define("z_reco",  "vtx_z[0]").Define("z_truth", "truth_vtx_z[0]")
                .Define("dz_um",   "(vtx_z[0] - truth_vtx_z[0]) * 1000.f");     // mm -> um
    };
    auto nn = def(dnew), no = def(dold);
    // edges offset by half a bin so the truth value -71.2000 mm is bin-CENTRED, not on an edge
    auto h_reco  = nn.Histo1D({"", "", 40, -71.2975, -71.0975}, "z_reco");   // data span -71.25..-71.16 mm
    auto h_truth = nn.Histo1D({"", "", 40, -71.2975, -71.0975}, "z_truth");
    auto r_new   = nn.Histo1D({"", "", 80, -80., 80.}, "dz_um");
    auto r_old   = no.Histo1D({"", "", 80, -80., 80.}, "dz_um");

    TCanvas c("c", "", 1600, 700);
    c.Divide(2, 1);
    // log y: the truth point is a delta (all 1e4 events in one 5 um bin) next to a reco
    // distribution falling over ~3 decades to single events; linear would hide the reco tails.
    c.cd(1); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13); gPad->SetLogy();
    axes(h_reco.GetPtr(), "vertex z [mm]", "events / 5 #mum");
    h_reco->GetXaxis()->SetNdivisions(505);
    h_reco->SetLineColor(kNew); h_reco->SetLineWidth(2);
    h_truth->SetLineColor(kBlack); h_truth->SetLineWidth(2); h_truth->SetLineStyle(2);
    h_reco->SetMaximum(30 * std::max(h_reco->GetBinContent(h_reco->GetMaximumBin()), h_truth->GetBinContent(h_truth->GetMaximumBin())));
    h_reco->SetMinimum(0.5);
    h_reco->Draw("hist"); h_truth->Draw("hist same");
    TLegend l1(0.16, 0.72, 0.92, 0.90); l1.SetBorderSize(0); l1.SetFillStyle(0); l1.SetTextSize(0.035);
    l1.SetHeader("Pb+Pb 2024 conditions (r17864), #hat{p}_{T} 125-300 GeV");
    l1.AddEntry(h_reco.GetPtr(),  "reconstructed primary vertex", "l");
    l1.AddEntry(h_truth.GetPtr(), "truth interaction point", "l");
    l1.Draw();
    c.cd(2); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
    axes(r_new.GetPtr(), "z_{reco} - z_{truth} [#mum]", "events / 2 #mum");
    r_new->SetLineColor(kNew); r_new->SetLineWidth(2);
    r_old->SetLineColor(kOld); r_old->SetLineWidth(2);
    r_new->SetMaximum(1.7 * std::max(r_new->GetBinContent(r_new->GetMaximumBin()), r_old->GetBinContent(r_old->GetMaximumBin())));
    r_new->SetMinimum(0);
    r_new->Draw("hist"); r_old->Draw("hist same");
    TLegend l3(0.16, 0.64, 0.92, 0.90); l3.SetBorderSize(0); l3.SetFillStyle(0); l3.SetTextSize(0.034);
    l3.AddEntry(r_new.GetPtr(), "2024 conditions (r17864), z_{vtx} = -71.2 mm", "l");
    l3.AddEntry((TObject*)nullptr, Form("mean %+.1f #mum, RMS %.1f #mum", r_new->GetMean(), r_new->GetRMS()), "");
    l3.AddEntry(r_old.GetPtr(), "2023 conditions (r17618), z_{vtx} = -3.3 mm", "l");
    l3.AddEntry((TObject*)nullptr, Form("mean %+.1f #mum, RMS %.1f #mum", r_old->GetMean(), r_old->GetRMS()), "");
    l3.Draw();
    c.SaveAs((OUTDIR + "vertex_z_reco_vs_truth_r17864.png").c_str());
    std::printf("vertex z: new reco mean %.4f mm, truth mean %.4f mm; residual mean %+.2f um RMS %.2f um (old: %+.2f / %.2f um)\n",
                h_reco->GetMean(), h_truth->GetMean(), r_new->GetMean(), r_new->GetRMS(), r_old->GetMean(), r_old->GetRMS());
}

void plot_r17864_event_level(const char* what = "fcal") {
    const std::string w = what;
    if      (w == "fcal") plot_fcal();
    else if (w == "vtx")  plot_vtx();
    else std::printf("unknown mode '%s' (fcal|vtx)\n", what);
}
