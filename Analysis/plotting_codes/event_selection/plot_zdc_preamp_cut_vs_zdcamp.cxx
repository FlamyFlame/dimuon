// plot_zdc_preamp_cut_vs_zdcamp.cxx
//
// Ratio of the hard ZDC preamp cut (ADC) to the mean ZDC amplitude sum
// (zdc_ZdcAmp, same side) for events passing trigger + Cut1 + Cut2.
// Plotted for side A and side C across the three PbPb years.
//
// If the ratio is approximately constant, the inter-year cut variation tracks
// the amplitude scale (calibration shift).  A varying ratio signals a genuine
// change in the noise/pile-up threshold relative to the signal amplitude.
//
// Usage (run from Analysis/):
//   source /usatlas/u/yuhanguo/setup.sh
//   root -l -b -q 'plotting_codes/event_selection/plot_zdc_preamp_cut_vs_zdcamp.cxx+'

#include <string>
#include <vector>
#include <iostream>
#include "TChain.h"
#include "TFile.h"
#include "TParameter.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TSystem.h"
#include "TStyle.h"
#include "../../NTupleProcessingCode/PbPbEventSelConfig.h"

static const std::string kBase    = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/";
static const std::string kPlotDir = kBase + "plots/";

static std::vector<std::string> FilesForYear(int yr) {
    std::string base = kBase + "pbpb_20" + std::to_string(yr) + "/";
    if (yr == 23) return { base+"data_pbpb23_part1.root",
                           base+"data_pbpb23_part2.root",
                           base+"data_pbpb23_part3.root" };
    if (yr == 24) return { base+"data_pbpb24_part1.root",
                           base+"data_pbpb24_part2.root" };
    if (yr == 25) return { base+"data_pbpb25_part1.root",
                           base+"data_pbpb25_part2.root",
                           base+"data_pbpb25_part3.root",
                           base+"data_pbpb25_part4.root",
                           base+"data_pbpb25_part5.root",
                           base+"data_pbpb25_part6.root" };
    return {};
}

struct YearResult {
    double ratio_A{0.}, ratio_C{0.};  // cut_preamp / mean_zdcAmp
};

static YearResult ProcessYear(int yr) {
    const std::string cpath = PbPbEvSelCutsPath(2000 + yr);
    TFile* fc = TFile::Open(cpath.c_str(), "READ");
    if (!fc || fc->IsZombie())
        throw std::runtime_error("Cuts file not found: " + cpath);
    TGraph* g_cut1  = (TGraph*)((TGraph*)fc->Get(PbPbEvSelKey::kZDCFCalCut))->Clone();
    double cut2_ns  = ((TParameter<double>*)fc->Get(PbPbEvSelKey::kZDCTimeCutNs))->GetVal();
    double cut_A    = ((TParameter<double>*)fc->Get(PbPbEvSelKey::kPreampACutADC))->GetVal();
    double cut_C    = ((TParameter<double>*)fc->Get(PbPbEvSelKey::kPreampCCutADC))->GetVal();
    fc->Close();

    TChain chain("HeavyIonD3PD", "HeavyIonD3PD");
    for (const auto& f : FilesForYear(yr)) {
        if (gSystem->AccessPathName(f.c_str()))
            { std::cerr << "Skipping: " << f << std::endl; continue; }
        chain.Add(f.c_str());
    }
    chain.SetMakeClass(1);
    chain.SetBranchStatus("*", 0);
    Int_t   b_HLT  = 0;
    Float_t FCal_P = 0.f, FCal_N = 0.f;
    Float_t zdcE[2]{}, zdcT[2]{}, zdcAmp[2]{};
    chain.SetBranchStatus("b_HLT_mu4_L1MU3V", 1);
    chain.SetBranchStatus("FCal_Et_P",         1);
    chain.SetBranchStatus("FCal_Et_N",         1);
    chain.SetBranchStatus("zdc_ZdcEnergy",     1);
    chain.SetBranchStatus("zdc_ZdcTime",       1);
    chain.SetBranchStatus("zdc_ZdcAmp",        1);
    chain.SetBranchAddress("b_HLT_mu4_L1MU3V", &b_HLT);
    chain.SetBranchAddress("FCal_Et_P",         &FCal_P);
    chain.SetBranchAddress("FCal_Et_N",         &FCal_N);
    chain.SetBranchAddress("zdc_ZdcEnergy",     zdcE);
    chain.SetBranchAddress("zdc_ZdcTime",       zdcT);
    chain.SetBranchAddress("zdc_ZdcAmp",        zdcAmp);

    // apply trigger + Cut1 + Cut2 only (NOT Cut3 — that's what we study)
    double sumA = 0., sumC = 0.;
    long long n_sel = 0;
    const Long64_t n = chain.GetEntries();
    std::cout << "20" << yr << ": " << n << " events..." << std::flush;
    for (Long64_t i = 0; i < n; ++i) {
        chain.GetEntry(i);
        if (!b_HLT) continue;
        const float fcal_AC = (FCal_P + FCal_N) * 1e-6f;
        const float zdcTot  = (zdcE[0] + zdcE[1]) / 1000.f;
        if (zdcTot > (float)PbPbEvSelEvalCut(g_cut1, fcal_AC)) continue; // Cut1
        if (std::abs(zdcT[1]) >= (float)cut2_ns) continue;               // Cut2 A  [1]=A
        if (std::abs(zdcT[0]) >= (float)cut2_ns) continue;               // Cut2 C  [0]=C
        sumA += zdcAmp[1];   // [1] = side A
        sumC += zdcAmp[0];   // [0] = side C
        ++n_sel;
    }
    delete g_cut1;
    std::cout << "  " << n_sel << " pass trigger+Cut1+Cut2\n";

    const double meanA = (n_sel > 0) ? sumA / n_sel : 1.;
    const double meanC = (n_sel > 0) ? sumC / n_sel : 1.;
    YearResult res;
    res.ratio_A = cut_A / meanA;
    res.ratio_C = cut_C / meanC;
    std::cout << "  mean_zdcAmp_A=" << meanA << "  mean_zdcAmp_C=" << meanC
              << "  cut_A=" << cut_A << "  cut_C=" << cut_C
              << "  ratio_A=" << res.ratio_A << "  ratio_C=" << res.ratio_C << "\n";
    return res;
}

void plot_zdc_preamp_cut_vs_zdcamp() {
    gStyle->SetOptStat(0);
    gSystem->mkdir(kPlotDir.c_str(), true);

    const int years[3] = {23, 24, 25};
    YearResult R[3];
    for (int i = 0; i < 3; ++i) R[i] = ProcessYear(years[i]);

    const double xs[3]     = {0., 1., 2.};
    double ratioA[3], ratioC[3];
    for (int i = 0; i < 3; ++i) { ratioA[i] = R[i].ratio_A; ratioC[i] = R[i].ratio_C; }

    // y-axis range with 25% padding above/below
    double lo = *std::min_element(ratioA, ratioA+3);
    lo = std::min(lo, *std::min_element(ratioC, ratioC+3));
    double hi = *std::max_element(ratioA, ratioA+3);
    hi = std::max(hi, *std::max_element(ratioC, ratioC+3));
    lo = std::min(lo, 0.);   // always include zero
    const double pad = 0.25 * (hi - lo);
    const double ylo = lo - pad, yhi = hi + pad;

    TCanvas* c = new TCanvas("c", "ZDC preamp cut / ZdcAmp", 600, 550);
    c->SetLeftMargin(0.17);
    c->SetRightMargin(0.05);
    c->SetBottomMargin(0.18);
    c->SetTopMargin(0.07);

    TGraph* gA = new TGraph(3, xs, ratioA);
    TGraph* gC = new TGraph(3, xs, ratioC);
    gA->SetMarkerStyle(20); gA->SetMarkerSize(1.8);
    gA->SetMarkerColor(kBlack); gA->SetLineColor(kBlack); gA->SetLineWidth(2);
    gC->SetMarkerStyle(21); gC->SetMarkerSize(1.8);
    gC->SetMarkerColor(kBlue+1); gC->SetLineColor(kBlue+1); gC->SetLineWidth(2);

    // invisible frame to set axis range and labels
    TGraph* fr = new TGraph(3, xs, ratioA);
    fr->SetTitle(";PbPb year;Preamp cut / mean ZDC amplitude");
    fr->GetXaxis()->SetLimits(-0.5, 2.5);
    fr->GetYaxis()->SetRangeUser(ylo, yhi);
    fr->GetXaxis()->SetNdivisions(3);
    fr->GetXaxis()->SetLabelOffset(999);
    fr->GetYaxis()->SetTitleOffset(1.6);
    fr->GetYaxis()->SetTitleSize(0.050);
    fr->GetXaxis()->SetTitleSize(0.050);
    fr->SetMarkerStyle(1); fr->SetMarkerColor(0); fr->SetLineColor(0);
    fr->Draw("AP");
    fr->GetXaxis()->SetLimits(-0.5, 2.5);
    fr->GetYaxis()->SetRangeUser(ylo, yhi);

    // dashed line at y=0 if range crosses zero
    if (ylo < 0. && yhi > 0.) {
        TLine* zl = new TLine(-0.5, 0., 2.5, 0.);
        zl->SetLineStyle(2); zl->SetLineColor(kGray+1);
        zl->Draw();
    }

    gA->Draw("LP same");
    gC->Draw("LP same");

    // year labels on x-axis
    TLatex lab;
    lab.SetTextAlign(22); lab.SetTextSize(0.055);
    const double laby = ylo - 0.09*(yhi - ylo);
    lab.DrawLatex(0., laby, "2023");
    lab.DrawLatex(1., laby, "2024");
    lab.DrawLatex(2., laby, "2025");

    // value annotations
    TLatex val;
    val.SetTextSize(0.042); val.SetTextAlign(21);
    const double off = 0.04 * (yhi - ylo);
    for (int i = 0; i < 3; ++i) {
        val.SetTextColor(kBlack);
        val.DrawLatex(xs[i] - 0.09, ratioA[i] + off, Form("%.2f", ratioA[i]));
        val.SetTextColor(kBlue+1);
        val.DrawLatex(xs[i] + 0.09, ratioC[i] + off, Form("%.2f", ratioC[i]));
    }

    TLegend* leg = new TLegend(0.20, 0.75, 0.50, 0.91);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.050);
    leg->AddEntry(gA, "Side A", "lp");
    leg->AddEntry(gC, "Side C", "lp");
    leg->Draw();

    TLatex atl; atl.SetNDC(); atl.SetTextSize(0.038); atl.SetTextAlign(11);
    atl.DrawLatex(0.19, 0.09, "#bf{ATLAS Internal} | ZDC preamp cut vs ZDC amp");

    const std::string out = kPlotDir + "zdc_preamp_cut_vs_zdcamp.png";
    c->SaveAs(out.c_str());
    std::cout << "Saved: " << out << std::endl;
}
