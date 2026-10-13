// Overlay of two event-selection acceptances (per-cut survival and cumulative), same definitions and style as
// plot_pbpb_event_sel_cuts.cxx DrawSurvivalPlot / DrawCumulativeSurvivalPlot (HLT_mu4, centrality 0-79).
// Inputs = the derivation's own printed counts ("Cut k (...): p survival (n_k / n_before ctr<80)" and
// "Cut k (...): cumulative p (n_k / n_total ctr<80)"), read from a log or a saved text file; the LAST block in the file is used.
// Used for the cut-2 window change (docs/tracking/pbpb_zdc_time_cut2_window.md): old vs new, one PNG each, per year.
#include <cmath>
#include <fstream>
#include <regex>
#include <stdexcept>
#include <string>
#include <vector>
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include "../../NTupleProcessingCode/PbPbEventSelConfig.h"

namespace evsov {
struct Counts { double surv_n[5], surv_d[5], cum_n[5], cum_d[5]; };
inline Counts Parse(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    const std::regex rs(R"(Cut (\d) \(\w+\): [0-9.]+ survival \((\d+) / (\d+) ctr<80\))");
    const std::regex rc(R"(Cut (\d) \(\w+\): cumulative [0-9.]+ \((\d+) / (\d+) ctr<80\))");
    Counts c{}; int ns = 0, nc = 0; std::string line; std::smatch m;
    while (std::getline(in, line)) {   // later blocks overwrite earlier ones -> last block wins
        if (std::regex_search(line, m, rs)) { int k = std::stoi(m[1]) - 1; c.surv_n[k] = std::stod(m[2]); c.surv_d[k] = std::stod(m[3]); ++ns; }
        if (std::regex_search(line, m, rc)) { int k = std::stoi(m[1]) - 1; c.cum_n[k] = std::stod(m[2]); c.cum_d[k] = std::stod(m[3]); ++nc; }
    }
    if (ns < 5 || nc < 5) throw std::runtime_error("incomplete survival block in " + path);
    return c;
}
inline TGraphErrors* Graph(const double* n, const double* d, double dx, int col, int mk) {
    auto* g = new TGraphErrors(5);
    for (int k = 0; k < 5; ++k) {
        const double p = n[k] / d[k];
        g->SetPoint(k, k + 1 + dx, p); g->SetPointError(k, 0, std::sqrt(p * (1 - p) / d[k]));
    }
    g->SetMarkerStyle(mk); g->SetMarkerSize(1.2); g->SetMarkerColor(col); g->SetLineColor(col);
    return g;
}
}

void plot_pbpb_event_sel_survival_old_vs_new(int yr, const std::string& old_counts, const std::string& new_counts,
                                             const std::string& old_label, const std::string& new_label,
                                             const std::string& outdir) {
    using namespace evsov;
    gStyle->SetOptStat(0);
    const Counts o = Parse(old_counts), n = Parse(new_counts);
    for (int cum = 0; cum < 2; ++cum) {
        const double* on = cum ? o.cum_n : o.surv_n; const double* od = cum ? o.cum_d : o.surv_d;
        const double* nn = cum ? n.cum_n : n.surv_n; const double* nd = cum ? n.cum_d : n.surv_d;
        double ymin = 1;
        for (int k = 0; k < 5; ++k) ymin = std::min({ymin, on[k] / od[k], nn[k] / nd[k]});
        TH1D fr(Form("fr_%d", cum), cum ? ";;Cumulative survival fraction (centrality 0#minus79)"
                                        : ";;Survival fraction per cut (centrality 0#minus79)", 5, 0.5, 5.5);
        for (int k = 0; k < 5; ++k) fr.GetXaxis()->SetBinLabel(k + 1, Form("C%d  %s", k + 1, kPbPbEvSelCutLabel(k + 1)));
        fr.GetYaxis()->SetRangeUser(std::max(0., ymin - 0.6 * (1.004 - ymin)), 1.004);   // room below the points for the legend/table
        fr.GetXaxis()->SetLabelSize(0.038);
        TCanvas c(Form("c_%d", cum), "", 800, 600); c.SetLeftMargin(0.14); c.SetBottomMargin(0.18);
        fr.Draw("AXIS");
        TGraphErrors* go = Graph(on, od, -0.08, kGray + 2, 24);
        TGraphErrors* gn = Graph(nn, nd, +0.08, kBlue + 1, 20);
        go->Draw("PE SAME"); gn->Draw("PE SAME");
        TLegend lg(0.16, 0.22, 0.60, 0.36); lg.SetBorderSize(0); lg.SetFillStyle(0); lg.SetTextSize(0.03);
        lg.AddEntry(go, old_label.c_str(), "pe"); lg.AddEntry(gn, new_label.c_str(), "pe"); lg.Draw();
        TLatex tl; tl.SetNDC(); tl.SetTextSize(0.032);
        const double y0 = cum ? 0.86 : 0.62;   // cumulative falls left->right: table top-right; per-cut: bottom-right
        tl.DrawLatex(0.60, y0, Form("Pb+Pb 20%d data", yr));
        tl.SetTextSize(0.026);
        tl.DrawLatex(0.60, y0 - 0.05, "         old        new");
        for (int k = 0; k < 5; ++k)
            tl.DrawLatex(0.60, y0 - 0.10 - 0.045 * k, Form("C%d:  %.3f%%   %.3f%%", k + 1, 100 * on[k] / od[k], 100 * nn[k] / nd[k]));
        c.SaveAs((outdir + Form("/event_sel_survival%s_old_vs_new_pbpb_20%d.png", cum ? "_accum" : "", yr)).c_str());
    }
}
