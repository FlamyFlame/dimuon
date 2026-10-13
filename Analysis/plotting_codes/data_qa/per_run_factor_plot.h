// One marker per run, dataset-average line, same style as plot_signal_yield_per_lumi_per_run.
// Used by the factored per-run plots (run_norm_diagnostics.cxx for Pb+Pb, pp_norm_factors.cxx for pp24).
// y = num/den per run, stat error sqrt(num)/den; runs with NaN den are not drawn. The power of ten of the
// average is moved into the y title (no axis exponent label, which collides with the ATLAS label).
#pragma once
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"

inline void DrawPerRunFactor(const std::vector<int>& runs, const std::vector<double>& num, const std::vector<double>& den,
                             double mean, const std::string& ytitle, const std::string& unit, const std::string& label,
                             const std::string& outpath) {
    std::vector<int> idx;
    for (size_t i = 0; i < runs.size(); ++i) if (!std::isnan(den[i])) idx.push_back(i);
    const int nr = idx.size();
    const int pw = (int)std::floor(std::log10(mean));
    const double sc = std::pow(10., -pw);
    TH1D* fr = new TH1D(Form("frp_%s", outpath.c_str()), "", nr, 0, nr);
    TGraphErrors* g = new TGraphErrors(nr);
    double ymax = 0;
    for (int j = 0; j < nr; ++j) {
        const int i = idx[j];
        const double y = sc * num[i] / den[i], ey = sc * std::sqrt(std::max(num[i], 1.)) / den[i];
        g->SetPoint(j, j + 0.5, y); g->SetPointError(j, 0, ey); ymax = std::max(ymax, y + ey);
        fr->GetXaxis()->SetBinLabel(j + 1, std::to_string(runs[i]).c_str());
    }
    TCanvas* cc = new TCanvas(Form("cp_%s", outpath.c_str()), "", 1400, 750);
    cc->SetBottomMargin(0.16); cc->SetLeftMargin(0.09); cc->SetRightMargin(0.03); cc->SetTopMargin(0.07);
    fr->SetMinimum(0); fr->SetMaximum(1.35 * ymax);
    fr->GetXaxis()->SetLabelSize(nr > 40 ? 0.022 : 0.03); fr->GetXaxis()->LabelsOption("v"); fr->GetXaxis()->SetTickLength(0);
    fr->GetXaxis()->SetTitle("Run number"); fr->GetXaxis()->SetTitleOffset(2.2); fr->GetXaxis()->SetTitleSize(0.04);
    fr->GetYaxis()->SetTitle(Form("%s  [10^{%d}%s]", ytitle.c_str(), pw, unit.c_str()));
    fr->GetYaxis()->SetTitleSize(0.04); fr->GetYaxis()->SetTitleOffset(1.1);
    fr->Draw("AXIS");
    g->SetMarkerStyle(20); g->SetMarkerSize(1.0); g->SetLineWidth(2);
    TLine* ln = new TLine(0, sc * mean, nr, sc * mean); ln->SetLineStyle(2); ln->SetLineColor(kRed + 1); ln->SetLineWidth(2);
    ln->Draw(); g->Draw("P SAME");
    TLegend* lg = new TLegend(0.60, 0.78, 0.96, 0.90); lg->SetBorderSize(0); lg->SetFillStyle(0);
    lg->AddEntry(g, "per run (stat. error)", "pe");
    lg->AddEntry(ln, Form("dataset average %.4g%s", mean, unit.c_str()), "l"); lg->Draw();
    TLatex t; t.SetNDC(); t.SetTextSize(0.04); t.SetTextFont(72); t.DrawLatex(0.10, 0.94, "ATLAS");
    t.SetTextFont(42); t.DrawLatex(0.17, 0.94, ("Internal   " + label).c_str());
    cc->SaveAs(outpath.c_str());
}
