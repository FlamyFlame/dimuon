// Per-run fraction of events kept by the NEW cut-2 window relative to the OLD |t^A| < 1.5 ns window
// (CSV from cut2_window_per_run.cxx; events passing cut 1, t^C inside its window). One marker per run, binomial errors.
// docs/tracking/pbpb_zdc_time_cut2_window.md
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "../../NTupleProcessingCode/PbPbEventSelConfig.h"

void plot_cut2_window_per_run(int yr, const std::string& csv, const std::string& outpng) {
    gStyle->SetOptStat(0);
    std::vector<int> run; std::vector<double> f, ef;
    double S_new = 0, S_old = 0;
    std::ifstream in(csv); std::string line; std::getline(in, line);
    while (std::getline(in, line)) {
        std::stringstream ss(line); std::string c; std::vector<std::string> v;
        while (std::getline(ss, c, ',')) v.push_back(c);
        const double n_old = std::stod(v[2]) + std::stod(v[3]), n_new = std::stod(v[3]);
        if (n_old <= 0) continue;
        const double p = n_new / n_old;
        run.push_back(std::stoi(v[0])); f.push_back(p); ef.push_back(std::sqrt(p * (1 - p) / n_old));
        S_new += n_new; S_old += n_old;
    }
    const int nr = run.size();
    const auto w = PbPbZdcTimeCut(yr);
    TH1D* fr = new TH1D("fr", "", nr, 0, nr);
    TGraphErrors* g = new TGraphErrors(nr);
    for (int i = 0; i < nr; ++i) {
        g->SetPoint(i, i + 0.5, f[i]); g->SetPointError(i, 0, ef[i]);
        fr->GetXaxis()->SetBinLabel(i + 1, std::to_string(run[i]).c_str());
    }
    TCanvas c("c", "", 1400, 750);
    c.SetBottomMargin(0.16); c.SetLeftMargin(0.09); c.SetRightMargin(0.03); c.SetTopMargin(0.07);
    fr->SetMinimum(0); fr->SetMaximum(1.15);
    fr->GetXaxis()->SetLabelSize(nr > 40 ? 0.022 : 0.03); fr->GetXaxis()->LabelsOption("v"); fr->GetXaxis()->SetTickLength(0);
    fr->GetXaxis()->SetTitle("Run number"); fr->GetXaxis()->SetTitleOffset(2.2); fr->GetXaxis()->SetTitleSize(0.04);
    fr->GetYaxis()->SetTitle(Form("N(%.1f < t^{A} < %.1f ns) / N(|t^{A}| < 1.5 ns)", w.a_lo, w.a_hi));
    fr->GetYaxis()->SetTitleSize(0.04); fr->GetYaxis()->SetTitleOffset(1.0);
    fr->Draw("AXIS");
    TLine* one = new TLine(0, 1, nr, 1); one->SetLineStyle(3); one->SetLineColor(kGray + 2); one->Draw();
    TLine* avg = new TLine(0, S_new / S_old, nr, S_new / S_old); avg->SetLineStyle(2); avg->SetLineColor(kRed + 1); avg->SetLineWidth(2); avg->Draw();
    g->SetMarkerStyle(20); g->SetMarkerSize(1.0); g->SetLineWidth(2); g->Draw("P SAME");
    TLatex t; t.SetNDC(); t.SetTextSize(0.04); t.SetTextFont(72); t.DrawLatex(0.10, 0.94, "ATLAS");
    t.SetTextFont(42);
    t.DrawLatex(0.17, 0.94, Form("Internal   Pb+Pb 20%d, HLT_mu4, after cut 1, %.1f < t^{C} < %.1f ns", yr, w.c_lo, w.c_hi));
    t.SetTextSize(0.035); t.SetTextColor(kRed + 1);
    t.DrawLatex(0.12, 0.25, Form("dataset average %.4f", S_new / S_old));
    c.SaveAs(outpng.c_str());
}
