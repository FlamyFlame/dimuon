// Plot differential cross-section of pT / pT(last b-hadron) for single-B truth signal pairs.
// Selection: muon_pair_tree_sign2 with from_same_b == true.
// Three distributions: pair pT, lead muon pT, sublead muon pT, all normalised to the
// pT of the last b-hadron ancestor of the lead muon (m1_last_b_hadron_prt_pt_eta_phi_m[0]/1000 GeV).

#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TSystem.h>
#include <algorithm>

void plot_single_b_pt_ratio_to_last_b_hadron() {

    const std::string input_file =
        "/usatlas/u/yuhanguo/usatlasdata/pythia_truth_full_sample/pythia_5p36TeV/"
        "muon_pairs_pythia_5p36TeV_no_data_resonance_cuts.root";
    const std::string output_dir =
        "/usatlas/u/yuhanguo/usatlasdata/pythia_truth_full_sample/pythia_5p36TeV/plots/";
    gSystem->mkdir(output_dir.c_str(), true);

    ROOT::RDataFrame df("muon_pair_tree_sign2", input_file.c_str());

    // select single-b signal pairs; define b-hadron pT (GeV) from lead muon's last b ancestor
    auto df_ratios = df
        .Filter("from_same_b")
        .Define("b_pt",
            [](const std::vector<float>& v){ return v.empty() ? 0.f : v[0] / 1000.f; },
            {"m1_last_b_hadron_prt_pt_eta_phi_m"})
        .Filter("b_pt > 0")
        .Define("R_pair", "truth_pair_pt / b_pt")
        .Define("R_m1",   "m1.truth_pt    / b_pt")
        .Define("R_m2",   "m2.truth_pt    / b_pt");

    const int    nbins = 50;
    const double xmin  = 0., xmax = 1.;

    auto hp_pair = df_ratios.Histo1D({"h_R_pair", "", nbins, xmin, xmax}, "R_pair", "weight");
    auto hp_m1   = df_ratios.Histo1D({"h_R_m1",   "", nbins, xmin, xmax}, "R_m1",   "weight");
    auto hp_m2   = df_ratios.Histo1D({"h_R_m2",   "", nbins, xmin, xmax}, "R_m2",   "weight");

    TH1D* h_pair = (TH1D*) hp_pair->Clone("h_pair");
    TH1D* h_m1   = (TH1D*) hp_m1  ->Clone("h_m1");
    TH1D* h_m2   = (TH1D*) hp_m2  ->Clone("h_m2");

    // differential in mub: weight unit is mub; Scale("width") divides each bin by its width
    for (TH1D* h : {h_pair, h_m1, h_m2})
        h->Scale(1., "width");

    // style
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    auto style_hist = [](TH1D* h, int color, int marker) {
        h->SetLineColor(color);
        h->SetMarkerColor(color);
        h->SetMarkerStyle(marker);
        h->SetMarkerSize(0.8);
        h->SetLineWidth(1);
    };
    style_hist(h_pair, kRed,   20);
    style_hist(h_m1,   kBlack, 20);
    style_hist(h_m2,   kBlue,  20);

    double ymax = std::max({h_pair->GetMaximum(), h_m1->GetMaximum(), h_m2->GetMaximum()});

    TCanvas* c = new TCanvas("c", "c", 700, 600);
    c->SetLeftMargin(0.16);
    c->SetBottomMargin(0.13);
    c->SetRightMargin(0.05);

    h_pair->GetXaxis()->SetTitle("R = p_{T} / p_{T}^{last b hadron}");
    h_pair->GetYaxis()->SetTitle("d#sigma/dR [#mub]");
    h_pair->GetXaxis()->SetRangeUser(xmin, xmax);
    h_pair->GetYaxis()->SetRangeUser(0., ymax * 1.15);
    h_pair->Draw("E");
    h_m1->Draw("E same");
    h_m2->Draw("E same");

    TLegend* leg = new TLegend(0.52, 0.55, 0.93, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(43);
    leg->SetTextSize(16);
    leg->AddEntry(h_pair, "pair p_{T}",            "ep");
    leg->AddEntry(h_m1,   "p_{T}^{#mu, lead}",     "ep");
    leg->AddEntry(h_m2,   "p_{T}^{#mu, sublead}",  "ep");
    leg->Draw();

    TLatex txt;
    txt.SetNDC();
    txt.SetTextFont(43);
    txt.SetTextSize(15);
    txt.DrawLatex(0.52, 0.91, "Truth Single-B Signal");

    c->SaveAs((output_dir + "single_b_pt_ratio_to_last_b_hadron.png").c_str());
    delete c;
}
