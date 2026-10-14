// plot_dr_correction_applied_demo.cxx
//
// Presentation figure: the Step-3 2mu4 dR correction of ONE pair-pT bin, drawn with the fit form
// the cross-section APPLIES in each cell (docs/tracking/pp24_trig_eff_hybrid_application.md §3b),
// for BOTH sign series -- i.e. what is applied once same-sign pairs get their own numbers (D2
// reverted). Per cell and sign the method comes from the application config, never retyped here:
//   opposite sign: DrCorrCrossxPolyMethod() in the cells of DrCorrCrossxPolyPtBins("os") x the
//                  LAST |eta^pair| group, DrCorrCrossxMethod() (expo) elsewhere;
//   same sign:     DrCorrCrossxPolyPtBins("ss") (empty) -> expo everywhere.
// Inputs = the persisted fits of fit_dr_corrections.cxx in the applied mode DrCorrCrossxMode()
// (no plateau correction, 3-group |eta^pair| fold): TF1 f_step3_pt<i>_eta<j> + the measured
// TGraphErrors g_step3_pt<i>_eta<j> it was fitted to; fit_ok from h_step3_fit_ok must be 1.
// Pure read-back + drawing, no physics computation.
//
// Usage: root -l -b -q 'plot_dr_correction_applied_demo.cxx+(3)'   // 1-based coarse pair-pT bin

#include "dr_correction_sample_cfg.h"
#include "dr_correction_cell_groups.h"
#include "../../../MuonObjectsParamsAndHelpers/ParamsSet.h"
#include "../../../RDFBasedHistFilling/CommonEffcyConfig.h"

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
struct Cell { TF1* f = nullptr; TGraphErrors* g = nullptr; std::string method; };

Cell LoadCell(const DrCorrSample& s, const std::string& method, const std::string& sign,
              int ipt, int ieta) {
    const std::string path = DrCorrFitFile(s, /*tight*/ true, 3, method, sign, DrCorrCrossxMode());
    TFile* f = TFile::Open(path.c_str());
    if (!f || f->IsZombie()) throw std::runtime_error("cannot open " + path);
    const std::string tag = Form("pt%d_eta%d", ipt, ieta);
    Cell c;
    c.method = method;
    c.f = dynamic_cast<TF1*>(f->Get(("f_step3_" + tag).c_str()));
    c.g = dynamic_cast<TGraphErrors*>(f->Get(("g_step3_" + tag).c_str()));
    auto* ok = dynamic_cast<TH2D*>(f->Get("h_step3_fit_ok"));
    if (!c.f || !c.g || !ok) throw std::runtime_error("missing f/g/fit_ok for " + tag + " in " + path);
    if (ok->GetBinContent(ipt, ieta) != 1)
        throw std::runtime_error(tag + " is not fit_ok in " + path + " (the analysis would fall back)");
    return c;
}
}  // namespace

void plot_dr_correction_applied_demo(int ipt = 3) {
    gROOT->SetBatch(true);
    gStyle->SetOptStat(0);
    const DrCorrSample s = GetDrCorrSample("pp_full", true);

    ParamsSet pms;
    const std::vector<double>& ptb = pms.pair_pt_coarse_bins;
    if (ipt < 1 || ipt > (int)ptb.size() - 1) throw std::runtime_error("pair-pT bin out of range");
    // |eta^pair| groups: 0, the fold's interior boundaries, the canonical outer edge
    CommonEffcyConfig cfg;
    std::vector<double> etab = {0.};
    for (double b : DrEtaAbsMergeInteriorBoundaries()) etab.push_back(b);
    etab.push_back(cfg.pair_eta_proj_ranges_coarse_incl_gap.back().second);
    const int neta = (int)etab.size() - 1;

    const auto& poly_os = DrCorrCrossxPolyPtBins("os");
    const auto& poly_ss = DrCorrCrossxPolyPtBins("ss");
    auto method_for = [&](const std::vector<int>& poly, int ieta) -> std::string {
        const bool is_poly = ieta == neta && std::find(poly.begin(), poly.end(), ipt) != poly.end();
        return is_poly ? DrCorrCrossxPolyMethod() : DrCorrCrossxMethod();
    };

    TCanvas c("c", "", 1560, 620);
    const double top = 0.70;
    const int colSS = kBlue, colOS = kRed;
    std::vector<TObject*> keep;
    TF1* leg_f[2] = {nullptr, nullptr};
    TGraphErrors* leg_g[2] = {nullptr, nullptr};
    for (int ie = 1; ie <= neta; ++ie) {
        c.cd();
        TPad* p = new TPad(Form("p%d", ie), "", (ie - 1) / double(neta), 0.0, ie / double(neta), top);
        p->SetLeftMargin(0.15); p->SetRightMargin(0.04); p->SetBottomMargin(0.16); p->SetTopMargin(0.40);
        p->SetTicks(1, 1); p->Draw(); p->cd();
        TH1D* fr = new TH1D(Form("fr%d", ie), "", 100, 0., 1.);
        fr->SetMinimum(0.); fr->SetMaximum(1.4);
        fr->GetXaxis()->SetTitle("#DeltaR(#mu_{1}, #mu_{2})");
        fr->GetYaxis()->SetTitle("#varepsilon_{#DeltaR}^{2mu4}");
        for (TAxis* a : {fr->GetXaxis(), fr->GetYaxis()}) { a->SetTitleSize(0.065); a->SetLabelSize(0.055); }
        fr->GetYaxis()->SetTitleOffset(1.0);
        fr->Draw("axis");
        TLine* one = new TLine(0., 1., 1., 1.); one->SetLineStyle(2); one->SetLineColor(kGray + 2); one->Draw();

        TLatex t; t.SetNDC(); t.SetTextFont(42);
        t.SetTextSize(0.060);
        t.DrawLatex(0.15, 0.92, Form("%.1f < |#eta^{pair}| < %.1f", etab[ie - 1], etab[ie]));

        for (int is = 0; is < 2; ++is) {   // 0 = same sign, 1 = opposite sign
            const std::string sign = is == 0 ? "ss" : "os";
            Cell cl = LoadCell(s, method_for(is == 0 ? poly_ss : poly_os, ie), sign, ipt, ie);
            const int col = is == 0 ? colSS : colOS;
            cl.g->SetMarkerStyle(is == 0 ? 22 : 20); cl.g->SetMarkerSize(1.0);
            cl.g->SetMarkerColor(col + 1); cl.g->SetLineColor(col + 1);
            cl.f->SetLineColor(col); cl.f->SetLineWidth(2); cl.f->SetNpx(500);
            cl.g->Draw("p same");
            cl.f->Draw("same");
            if (ie == 1) { leg_f[is] = cl.f; leg_g[is] = cl.g; }
            // parameters of the drawn curve, one column per sign
            t.SetTextSize(0.040); t.SetTextColor(col);
            const double x0 = is == 0 ? 0.42 : 0.70;
            double y = 0.855;
            t.DrawLatex(x0, y, cl.method == DrCorrCrossxPolyMethod() ? "polynomial" : "exponential");
            for (int ip = 0; ip < cl.f->GetNpar(); ++ip) {
                y -= 0.037;
                double lo, hi; cl.f->GetParLimits(ip, lo, hi);
                const double v = cl.f->GetParameter(ip);
                const bool fixed = (lo == hi && lo != 0.) || cl.f->GetParError(ip) == 0.;
                // at a limit: absolute 1e-4 test (a range-scaled tolerance is meaningless for wide limits)
                const bool at_limit = !fixed && lo < hi && (std::fabs(v - lo) < 1e-4 || std::fabs(v - hi) < 1e-4);
                t.DrawLatex(x0, y, fixed    ? Form("%s = %.3g (fixed)", cl.f->GetParName(ip), v)
                                 : at_limit ? Form("%s = %.2g (at limit)", cl.f->GetParName(ip), v)
                                            : Form("%s = %.3g #pm %.2g", cl.f->GetParName(ip), v, cl.f->GetParError(ip)));
            }
            y -= 0.037;
            t.DrawLatex(x0, y, Form("#chi^{2}/ndf = %.2f", cl.f->GetNDF() > 0 ? cl.f->GetChisquare() / cl.f->GetNDF() : -1.));
            t.SetTextColor(kBlack);
        }
        gPad->RedrawAxis();
    }

    // header: sample, cell, definition and the two fit equations
    c.cd();
    TLatex h; h.SetNDC(); h.SetTextFont(42);
    h.SetTextSize(0.034);
    h.DrawLatex(0.015, 0.94, Form("Pythia8 pp, #sqrt{s} = 5.36 TeV (2024 conditions), Tight muons,   %.1f < p_{T}^{pair} < %.1f GeV",
                                  ptb[ipt - 1], ptb[ipt]));
    h.SetTextSize(0.027);
    h.DrawLatex(0.015, 0.875, "#varepsilon_{#DeltaR}^{2mu4}(#DeltaR) #equiv P(pair passes 2mu4 | #DeltaR) / (#varepsilon_{1}#varepsilon_{2}),   applied correction = f(#DeltaR)/C for #DeltaR < 1");
    h.DrawLatex(0.015, 0.825, "exponential:  f(#DeltaR) = C + A #upoint exp[-(#DeltaR/#lambda)^{p}]");
    h.DrawLatex(0.015, 0.775, "polynomial:  f(#DeltaR) = C + u^{2}[A + a_{3}(u - 1) + a_{4}(u^{2} - 1)],   u #equiv max(0, 1 - #DeltaR/R_{p})");
    TLegend* l = new TLegend(0.55, 0.76, 0.99, 0.92);
    l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(0.030); l->SetNColumns(2);
    l->AddEntry(leg_g[0], "same sign", "p");
    l->AddEntry(leg_g[1], "opposite sign", "p");
    l->AddEntry(leg_f[0], "fit, same sign", "l");
    l->AddEntry(leg_f[1], "fit, opposite sign", "l");
    l->Draw();

    const std::string outdir = s.plot_root + "mc_based/step3_dr_fit/applied_demo/";
    gSystem->mkdir(outdir.c_str(), kTRUE);
    const std::string out = outdir + Form("step3_dr_fit_applied_pairpt_%.1f_%.1f.png", ptb[ipt - 1], ptb[ipt]);
    c.SaveAs(out.c_str());
}
