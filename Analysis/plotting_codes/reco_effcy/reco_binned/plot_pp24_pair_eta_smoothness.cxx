// =================================================================================================
// plot_pp24_pair_eta_smoothness.cxx
//
// eta^pair smoothness test of the reco-binned 3D bin-by-bin correction
// (docs/tracking/pair_reco_eff_bin_by_bin_reco_binned.md §3c, §4.1). Reads the outputs of
// fill_pp24_data_pair_eta_smoothness.C and fill_pp24_fullsim_pair_eta_smoothness_closure.C.
//
// Figures (dir SmoothnessPlotDir()/<wp>/), each 9 panels = all p_T^pair + the 8 coarse p_T^pair bins,
// coarse eta^pair cell edges as thin dashed lines:
//   data_corrected_<wp>WP.png   dsigma/deta^pair: trigger + reco-binned eps (nominal test), trigger only,
//                               trigger + the crossx-applied reco weight (eps^{reco,SF}, pre-D7 numerator);
//                               ratio pad: applied / reco-binned corrected
//   data_vs_mc_truth_<wp>WP.png data (trigger + reco-binned) vs Pythia truth T shape normalised to the
//                               data integral of the panel; ratio pad data / MC
//   mc_closure_<wp>WP.png       MC R (truth in T, D7) corrected by 1/eps^reco vs MC T; ratio pad corrected R / T
//   smoothness_values_<wp>WP.txt
// Smoothness numbers (per panel):
//   * boundary steps: at every interior coarse boundary, the step of r = data/MC-shape between the
//     last fine bin left and the first right, pull = step / sigma (data and MC stat); the same pulls
//     for adjacent fine bins NOT at a boundary are the reference distribution;
//   * within-cell chi2: in each coarse cell, data fine bins vs the MC truth shape normalised to the
//     data in THAT cell (ndf = nfine - 1) -- the within-cell shape only;
//   * MC closure: per coarse cell sum(corrected R)/sum(T) - 1 (must be ~0 by construction) and the
//     within-cell chi2 of corrected R vs T (indicative: R and T share events).
//
// Usage (from Analysis/plotting_codes/reco_effcy/reco_binned/):
//   root -l -b -q 'plot_pp24_pair_eta_smoothness.cxx+(true)'    // Tight
//   root -l -b -q 'plot_pp24_pair_eta_smoothness.cxx+(false)'   // Medium
//   root -l -b -q 'plot_pp24_pair_eta_smoothness.cxx+(true,true)'   // Tight, COARSE eta^pair view
// COARSE view (user 2026-10-11): the same figures with every eta^pair histogram summed into the
// canonical coarse cells (CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap, via the
// evaluator's EtaEdges(); exact because every coarse edge is a fine edge, checked in
// PairEtaFineSmoothnessBins()). Files *_coarse_eta_<wp>WP.png + coarse_eta_values_<wp>WP.txt
// (per-cell values; the fine-bin smoothness statistics do not apply).
// =================================================================================================
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TArrow.h>
#include <TLine.h>
#include <TMath.h>
#include <TNamed.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include "pair_eta_smoothness_common.h"

namespace {

TH1D* Get(TFile* f, const std::string& n)
{
    auto* h = dynamic_cast<TH1D*>(f->Get(n.c_str()));
    if (!h) throw std::runtime_error("missing " + n + " in " + f->GetName());
    auto* c = static_cast<TH1D*>(h->Clone((n + "_c").c_str()));
    c->SetDirectory(nullptr);
    return c;
}

void CheckFine(const TH1D* h, const std::vector<double>& fine)
{
    if (h->GetNbinsX() != (int)fine.size() - 1) throw std::runtime_error(std::string(h->GetName()) + ": not the fine binning");
    for (int i = 1; i <= h->GetNbinsX() + 1; ++i)
        if (std::fabs(h->GetXaxis()->GetBinLowEdge(i) - fine[i - 1]) > 1e-6)
            throw std::runtime_error(std::string(h->GetName()) + ": fine edge mismatch");
}

// coarse cell of each fine bin (by bin centre); cell -> [first, last] fine bin (1-based)
struct CellMap { std::vector<int> first, last; };
CellMap MakeCellMap(const std::vector<double>& coarse, const TH1D* h)
{
    CellMap m;
    for (size_t c = 0; c + 1 < coarse.size(); ++c) {
        int f = -1, l = -1;
        for (int b = 1; b <= h->GetNbinsX(); ++b) {
            const double x = h->GetBinCenter(b);
            if (x > coarse[c] && x < coarse[c + 1]) { if (f < 0) f = b; l = b; }
        }
        m.first.push_back(f); m.last.push_back(l);
    }
    return m;
}

// sum fine bins into the coarse cells (contents and sum of w^2 add)
TH1D* ToCoarse(const TH1D* h, const std::vector<double>& coarse)
{
    auto* r = static_cast<TH1D*>(const_cast<TH1D*>(h)->Rebin((int)coarse.size() - 1, (std::string(h->GetName()) + "_coarse").c_str(), coarse.data()));
    r->SetDirectory(nullptr);
    return r;
}

TH1D* Ratio(const TH1D* a, const TH1D* b, const std::string& n)   // uncorrelated errors
{
    auto* r = static_cast<TH1D*>(a->Clone(n.c_str()));
    r->SetDirectory(nullptr);
    r->Divide(b);
    return r;
}
TH1D* ScaledTo(const TH1D* t, double integral, const std::string& n)
{
    auto* s = static_cast<TH1D*>(t->Clone(n.c_str()));
    s->SetDirectory(nullptr);
    const double it = s->Integral();
    if (it > 0) s->Scale(integral / it);
    return s;
}

struct Series { TH1D* h; std::string label; Color_t col; Style_t mrk; bool hist = false; };
struct PanelData {
    std::string label;
    std::vector<Series> s;     // top pad
    std::vector<Series> ratios;   // bottom pad
};

void DrawFigure(const std::string& path, const std::string& header, const std::string& wp_text,
                const std::string& line2, const std::vector<PanelData>& P, const std::vector<double>& coarse,
                const std::string& ytitle, const std::string& rtitle, double rlo, double rhi)
{
    const int np = P.size();
    const int ncol = (np <= 3) ? np : (int)std::ceil(std::sqrt((double)np));
    const int nrow = (int)std::ceil((double)np / ncol);
    const int hstrip = 150;
    TCanvas c("c", "", 470 * ncol, 400 * nrow + hstrip);
    const double fh = (double)hstrip / (400.0 * nrow + hstrip);
    TPad* head = new TPad("head", "", 0, 1 - fh, 1, 1);
    TPad* body = new TPad("body", "", 0, 0, 1, 1 - fh);
    head->SetFillStyle(0); body->SetFillStyle(0); head->Draw(); body->Draw();
    head->cd();
    {
        const double sz = 0.24 * 110.0 / hstrip;
        TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(sz);
        double y = 0.80;
        t.DrawLatex(0.02, y, header.c_str());
        t.SetTextAlign(31); t.DrawLatex(0.98, y, wp_text.c_str()); t.SetTextAlign(11);
        y -= 1.35 * sz; t.DrawLatex(0.02, y, line2.c_str());
        TLegend* leg = new TLegend(0.02, 0.02, 0.98, y - 0.6 * sz);
        leg->SetNColumns((int)P[0].s.size()); leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(sz);
        for (const auto& s : P[0].s) leg->AddEntry(s.h, s.label.c_str(), s.hist ? "l" : "lp");
        leg->Draw();
    }
    body->cd();
    body->Divide(ncol, nrow);
    for (int p = 0; p < np; ++p) {
        body->cd(p + 1);
        TVirtualPad* cell = gPad;
        TPad* top = new TPad(Form("t%d", p), "", 0, 0.32, 1, 1);
        top->SetBottomMargin(0.02); top->SetTopMargin(0.12); top->SetLeftMargin(0.16); top->SetTicks(1, 1);
        top->Draw(); top->cd();
        double ymax = 0;
        for (const auto& s : P[p].s)
            for (int b = 1; b <= s.h->GetNbinsX(); ++b) ymax = std::max(ymax, s.h->GetBinContent(b) + s.h->GetBinError(b));
        bool first = true;
        for (const auto& s : P[p].s) {
            s.h->SetStats(0); s.h->SetTitle("");
            s.h->SetLineColor(s.col); s.h->SetMarkerColor(s.col); s.h->SetMarkerStyle(s.mrk);
            s.h->SetMarkerSize(0.8); s.h->SetLineWidth(s.hist ? 2 : 1);
            s.h->GetYaxis()->SetRangeUser(0, 1.25 * ymax);
            s.h->GetYaxis()->SetTitle(ytitle.c_str());
            s.h->GetYaxis()->SetTitleSize(0.060); s.h->GetYaxis()->SetLabelSize(0.050); s.h->GetYaxis()->SetTitleOffset(1.25);
            s.h->GetYaxis()->SetMaxDigits(3);
            s.h->GetXaxis()->SetLabelSize(0); s.h->GetXaxis()->SetTitleSize(0);
            s.h->Draw(first ? (s.hist ? "HIST" : "PE") : (s.hist ? "HIST SAME" : "PE SAME"));
            first = false;
        }
        for (size_t k = 1; k + 1 < coarse.size(); ++k) {
            TLine* l = new TLine(coarse[k], 0, coarse[k], 1.25 * ymax);
            l->SetLineStyle(2); l->SetLineColor(kGray + 1); l->SetLineWidth(1); l->Draw();
        }
        for (const auto& s : P[p].s) s.h->Draw(s.hist ? "HIST SAME" : "PE SAME");
        TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.062); t.SetTextAlign(31);
        t.DrawLatex(0.95, 0.91, P[p].label.c_str());
        cell->cd();
        TPad* bot = new TPad(Form("b%d", p), "", 0, 0, 1, 0.32);
        bot->SetTopMargin(0.03); bot->SetBottomMargin(0.34); bot->SetLeftMargin(0.16); bot->SetTicks(1, 1);
        bot->Draw(); bot->cd();
        // Ratio range: the nominal [rlo, rhi], widened per panel so that every point (with its error)
        // stays visible, capped at [0, kRatioCap]; points beyond the cap are marked by an arrow.
        const double kRatioCap = 3.0;
        double plo = rlo, phi = rhi;
        for (const auto& rs : P[p].ratios)
            for (int b = 1; b <= rs.h->GetNbinsX(); ++b) {
                const double v = rs.h->GetBinContent(b), e = rs.h->GetBinError(b);
                if (v <= 0) continue;
                plo = std::min(plo, v - e); phi = std::max(phi, v + e);
            }
        plo = std::max(0.0, std::floor(plo * 10) / 10);
        phi = std::min(kRatioCap, std::ceil(phi * 10) / 10);
        for (size_t q = 0; q < P[p].ratios.size(); ++q) {
            TH1D* r = P[p].ratios[q].h;
            r->SetStats(0); r->SetTitle("");
            r->SetLineColor(P[p].ratios[q].col); r->SetMarkerColor(P[p].ratios[q].col);
            r->SetMarkerStyle(P[p].ratios[q].mrk); r->SetMarkerSize(0.7);
            r->GetYaxis()->SetRangeUser(plo, phi); r->GetYaxis()->SetNdivisions(504);
            r->GetYaxis()->SetTitle(rtitle.c_str());
            r->GetYaxis()->SetTitleSize(0.12); r->GetYaxis()->SetLabelSize(0.11); r->GetYaxis()->SetTitleOffset(0.62);
            r->GetXaxis()->SetTitle("#eta^{pair}");
            r->GetXaxis()->SetTitleSize(0.14); r->GetXaxis()->SetLabelSize(0.12); r->GetXaxis()->SetTitleOffset(1.0);
            r->Draw(q == 0 ? "PE" : "PE SAME");
        }
        TLine* one = new TLine(coarse.front(), 1, coarse.back(), 1); one->SetLineStyle(2); one->Draw();
        for (size_t k = 1; k + 1 < coarse.size(); ++k) {
            TLine* l = new TLine(coarse[k], plo, coarse[k], phi);
            l->SetLineStyle(2); l->SetLineColor(kGray + 1); l->Draw();
        }
        for (const auto& rs : P[p].ratios) {
            rs.h->Draw("PE SAME");
            for (int b = 1; b <= rs.h->GetNbinsX(); ++b)
                if (rs.h->GetBinContent(b) > phi) {   // off-scale above the cap
                    const double x = rs.h->GetBinCenter(b), y0 = plo + 0.70 * (phi - plo);
                    TArrow* a = new TArrow(x, y0, x, phi - 0.02 * (phi - plo), 0.012, "|>");
                    a->SetLineColor(rs.h->GetLineColor()); a->SetFillColor(rs.h->GetLineColor()); a->Draw();
                }
        }
    }
    gSystem->mkdir(gSystem->DirName(path.c_str()), kTRUE);
    c.SaveAs(path.c_str());
}

// --- smoothness statistics ---
struct StepStats { std::vector<double> boundary_pulls; double ref_rms = 0; int ref_n = 0; };
StepStats Steps(const TH1D* r, const CellMap& m)
{
    StepStats s;
    std::vector<int> bnd_left;   // fine bin to the left of each interior boundary
    for (size_t c = 0; c + 1 < m.first.size(); ++c) bnd_left.push_back(m.last[c]);
    double ss = 0;
    for (int b = 1; b < r->GetNbinsX(); ++b) {
        const double e2 = std::pow(r->GetBinError(b), 2) + std::pow(r->GetBinError(b + 1), 2);
        if (e2 <= 0 || r->GetBinContent(b) <= 0 || r->GetBinContent(b + 1) <= 0) {
            if (std::find(bnd_left.begin(), bnd_left.end(), b) != bnd_left.end()) s.boundary_pulls.push_back(NAN);
            continue;
        }
        const double pull = (r->GetBinContent(b + 1) - r->GetBinContent(b)) / std::sqrt(e2);
        if (std::find(bnd_left.begin(), bnd_left.end(), b) != bnd_left.end()) s.boundary_pulls.push_back(pull);
        else { ss += pull * pull; ++s.ref_n; }
    }
    s.ref_rms = s.ref_n ? std::sqrt(ss / s.ref_n) : 0;
    return s;
}
// within-cell chi2 of d vs t normalised to d inside the cell; returns {chi2, ndf}
std::pair<double, int> CellChi2(const TH1D* d, const TH1D* t, int f, int l)
{
    double sd = 0, st = 0;
    for (int b = f; b <= l; ++b) { sd += d->GetBinContent(b); st += t->GetBinContent(b); }
    if (sd <= 0 || st <= 0) return {0, 0};
    const double k = sd / st;
    double chi2 = 0; int n = 0;
    for (int b = f; b <= l; ++b) {
        const double e2 = std::pow(d->GetBinError(b), 2) + k * k * std::pow(t->GetBinError(b), 2);
        if (e2 <= 0) continue;
        chi2 += std::pow(d->GetBinContent(b) - k * t->GetBinContent(b), 2) / e2;
        ++n;
    }
    return {chi2, n - 1};
}

}  // namespace

void plot_pp24_pair_eta_smoothness(bool tight_WP = true, bool coarse_eta = false)
{
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    const std::string wp = tight_WP ? "tight" : "medium";
    const std::string wp_text = tight_WP ? "Tight muons" : "Medium muons";
    const std::vector<double> fine = PairEtaFineSmoothnessBins();
    const std::vector<double> coarse = PairRecoEffRecoBinnedEvaluator::EtaEdges();
    const std::vector<double> ept = PairRecoEffRecoBinnedEvaluator::PtEdges();
    const int npt = ept.size() - 1;
    const std::string outdir = SmoothnessPlotDir() + wp + "/";

    std::unique_ptr<TFile> fd(TFile::Open(SmoothnessDataFile(tight_WP).c_str(), "READ"));
    std::unique_ptr<TFile> fm(TFile::Open(SmoothnessMCFile().c_str(), "READ"));
    if (!fd || fd->IsZombie() || !fm || fm->IsZombie()) throw std::runtime_error("cannot open inputs");
    {
        auto* w = dynamic_cast<TNamed*>(fd->Get("working_point"));
        if (!w || std::string(w->GetTitle()) != wp) throw std::runtime_error("data file WP mismatch");
    }

    std::vector<std::string> tags{""}, labels{Form("%.4g < p_{T}^{pair} < %.4g GeV", ept.front(), ept.back())};
    for (int k = 0; k < npt; ++k) {
        tags.push_back(Form("_pt%d", k));
        labels.push_back(Form("%.4g < p_{T}^{pair} < %.4g GeV", ept[k], ept[k + 1]));
    }

    gSystem->mkdir(outdir.c_str(), kTRUE);
    std::ofstream tab(coarse_eta ? std::string("/dev/null") : outdir + "smoothness_values_" + wp + "WP.txt");
    std::ofstream ctab(coarse_eta ? outdir + "coarse_eta_values_" + wp + "WP.txt" : std::string("/dev/null"));
    ctab << "Coarse eta^pair view, pp24 2mu4 data, " << wp_text << ". Data: " << SmoothnessDataFile(tight_WP)
         << "\nMC: " << SmoothnessMCFile() << "\nYields per coarse cell (integrated, not per unit eta).\n";
    tab << "eta^pair smoothness test, pp24 2mu4 data, " << wp_text << ". Data: " << SmoothnessDataFile(tight_WP)
        << "\nMC: " << SmoothnessMCFile() << "\nFine binning pair_eta_fine_smoothness (0.1); coarse cell edges:";
    for (double e : coarse) tab << " " << e;
    tab << "\nr = data(trig + reco-binned) / MC truth shape (normalised to the data integral of the panel).\n"
        << "Boundary pulls = step of r across each interior coarse boundary / sigma; reference = rms of the same pull for\n"
        << "adjacent fine bins NOT at a coarse boundary. Within-cell chi2: data vs MC truth shape normalised inside the cell.\n";

    std::vector<PanelData> P1, P2, P3;
    double max_closure_dev = 0, max_closure_dev_incl = 0;
    double chi2_rb_tot = 0, chi2_tb_tot = 0, chi2_mc_tot = 0, chi2_dr_tot = 0; int ndf_rb_tot = 0, ndf_tb_tot = 0, ndf_mc_tot = 0, ndf_dr_tot = 0;
    double max_closure_dev_unmeas = 0;
    PairRecoEffRecoBinnedEvaluator ev;
    ev.Load(RecoBinnedProductFile(), tight_WP);
    for (size_t p = 0; p < tags.size(); ++p) {
        const std::string& tg = tags[p];
        TH1D* d_rb  = Get(fd.get(), "h_pair_eta_fine_trig_reco_rb" + tg);
        TH1D* d_tr  = Get(fd.get(), "h_pair_eta_fine_trig" + tg);
        TH1D* d_tb  = Get(fd.get(), "h_pair_eta_fine_trig_reco_applied" + tg);
        TH1D* d_cnt = Get(fd.get(), "h_pair_eta_fine_counts" + tg);
        TH1D* mT    = Get(fm.get(), "h_mcT_fine" + tg);
        TH1D* mRc   = Get(fm.get(), "h_mcR_fine_corr_" + wp + tg);
        for (TH1D* h : {d_rb, d_tr, d_tb, d_cnt, mT, mRc}) CheckFine(h, fine);
        if (coarse_eta) {
            d_rb = ToCoarse(d_rb, coarse); d_tr = ToCoarse(d_tr, coarse); d_tb = ToCoarse(d_tb, coarse);
            d_cnt = ToCoarse(d_cnt, coarse); mT = ToCoarse(mT, coarse); mRc = ToCoarse(mRc, coarse);
            TH1D* d_raw = ToCoarse(Get(fd.get(), "h_pair_eta_fine_raw" + tg), coarse);
            TH1D* mT_n = ScaledTo(mT, d_rb->Integral(), "mTn" + tg);
            ctab << "\n[" << labels[p] << "]  cell | data pairs | raw [pb] | trig [pb] | trig+reco-binned [pb] | trig+applied eps^{reco,SF} [pb] | "
                    "applied/reco-binned | data(reco-binned)/MC truth (panel-normalised) | MC corrected R / T - 1\n";
            for (int c = 1; c <= d_rb->GetNbinsX(); ++c) {
                const double rr = mT_n->GetBinContent(c) > 0 ? d_rb->GetBinContent(c) / mT_n->GetBinContent(c) : 0;
                const double re = mT_n->GetBinContent(c) > 0 ? rr * std::hypot(d_rb->GetBinError(c) / std::max(d_rb->GetBinContent(c), 1e-300),
                                                                               mT_n->GetBinError(c) / mT_n->GetBinContent(c)) : 0;
                ctab << Form("    [%+.1f,%+.1f] %7.0f | %8.2f | %8.2f | %8.2f +- %6.2f | %8.2f +- %6.2f | %.4f | %.3f +- %.3f | %+.2e\n",
                             coarse[c - 1], coarse[c], d_cnt->GetBinContent(c), d_raw->GetBinContent(c), d_tr->GetBinContent(c),
                             d_rb->GetBinContent(c), d_rb->GetBinError(c), d_tb->GetBinContent(c), d_tb->GetBinError(c),
                             d_rb->GetBinContent(c) > 0 ? d_tb->GetBinContent(c) / d_rb->GetBinContent(c) : 0, rr, re,
                             mT->GetBinContent(c) > 0 ? mRc->GetBinContent(c) / mT->GetBinContent(c) - 1 : 0);
            }
        }
        const CellMap cm = MakeCellMap(coarse, d_rb);

        // ---- numbers (on yields, before the width scaling; widths are uniform anyway) ----
        TH1D* mT_d = ScaledTo(mT, d_rb->Integral(), "mT_d" + tg);
        TH1D* r_dm = Ratio(d_rb, mT_d, "r_dm" + tg);
        TH1D* mT_tb = ScaledTo(mT, d_tb->Integral(), "mT_tb" + tg);
        TH1D* r_tb = Ratio(d_tb, mT_tb, "r_tb" + tg);
        TH1D* mR_d = ScaledTo(mRc, d_rb->Integral(), "mR_d" + tg);
        TH1D* r_dr = Ratio(d_rb, mR_d, "r_dr" + tg);           // data / MC corrected reco (same sculpting)
        TH1D* r_mc = Ratio(mRc, mT, "r_mc" + tg);              // MC closure ratio
        if (!coarse_eta) {
        const StepStats st_dr = Steps(r_dr, cm), st_mc = Steps(r_mc, cm);
        const StepStats st_rb = Steps(r_dm, cm), st_tb = Steps(r_tb, cm);
        tab << "\n[" << labels[p] << "]  data pairs (unweighted) = " << d_cnt->Integral() << "\n";
        tab << "  boundary pulls (reco-binned eps):";
        for (size_t k = 0; k < st_rb.boundary_pulls.size(); ++k) tab << Form("  %+.4g:%+.2f", coarse[k + 1], st_rb.boundary_pulls[k]);
        tab << Form("   | reference rms %.2f (n=%d)\n", st_rb.ref_rms, st_rb.ref_n);
        tab << "  boundary pulls (crossx-applied eps^{reco,SF}, comparison):";
        for (size_t k = 0; k < st_tb.boundary_pulls.size(); ++k) tab << Form("  %+.4g:%+.2f", coarse[k + 1], st_tb.boundary_pulls[k]);
        tab << Form("   | reference rms %.2f (n=%d)\n", st_tb.ref_rms, st_tb.ref_n);
        tab << "  boundary pulls of data(reco-binned)/MC(corrected reco, normalised):";
        for (size_t k = 0; k < st_dr.boundary_pulls.size(); ++k) tab << Form("  %+.4g:%+.2f", coarse[k + 1], st_dr.boundary_pulls[k]);
        tab << Form("   | reference rms %.2f (n=%d)\n", st_dr.ref_rms, st_dr.ref_n);
        tab << "  boundary pulls of MC closure corrected-R/T (R,T share events: indicative):";
        for (size_t k = 0; k < st_mc.boundary_pulls.size(); ++k) tab << Form("  %+.4g:%+.2f", coarse[k + 1], st_mc.boundary_pulls[k]);
        tab << Form("   | reference rms %.2f (n=%d)\n", st_mc.ref_rms, st_mc.ref_n);
        tab << "  within-cell chi2/ndf  data(reco-binned) vs MC T | data(applied eps^{reco,SF}) vs MC T | MC corrected R vs MC T | data(reco-binned) vs MC corrected R | MC closure sum(corrR)/sum(T)-1 [* = cell contains a 3D cell without eps measurement]\n";
        double c2rb = 0, c2tb = 0, c2mc = 0, c2dr = 0; int nrb = 0, ntb = 0, nmc = 0, ndr = 0;
        for (size_t c = 0; c < cm.first.size(); ++c) {
            const auto a = CellChi2(d_rb, mT, cm.first[c], cm.last[c]);
            const auto b = CellChi2(d_tb, mT, cm.first[c], cm.last[c]);
            const auto m = CellChi2(mRc, mT, cm.first[c], cm.last[c]);
            const auto x = CellChi2(d_rb, mRc, cm.first[c], cm.last[c]);
            bool unmeas = false;   // any 3D (pT, eta, m) cell of this coarse eta cell without a measurement
            for (int ip = (p == 0 ? 0 : (int)p - 1); ip <= (p == 0 ? npt - 1 : (int)p - 1); ++ip)
                for (int im = 0; im < ev.h_eff->GetNbinsZ(); ++im)
                    if (ev.h_eff->GetBinContent(ip + 1, (int)c + 1, im + 1) <= 0) unmeas = true;
            double sr = 0, stt = 0;
            for (int bb = cm.first[c]; bb <= cm.last[c]; ++bb) { sr += mRc->GetBinContent(bb); stt += mT->GetBinContent(bb); }
            const double dev = stt > 0 ? sr / stt - 1 : 0;
            if (unmeas) max_closure_dev_unmeas = std::max(max_closure_dev_unmeas, std::fabs(dev));
            else if (p > 0) max_closure_dev = std::max(max_closure_dev, std::fabs(dev));
            else max_closure_dev_incl = std::max(max_closure_dev_incl, std::fabs(dev));
            c2rb += a.first; nrb += a.second; c2tb += b.first; ntb += b.second; c2mc += m.first; nmc += m.second;
            c2dr += x.first; ndr += x.second;
            tab << Form("    [%+.1f,%+.1f]  %6.1f/%d (p=%.3f) | %6.1f/%d (p=%.3f) | %6.1f/%d (p=%.3f) | %6.1f/%d (p=%.3f) | %+.2e%s\n",
                        coarse[c], coarse[c + 1], a.first, a.second, a.second > 0 ? TMath::Prob(a.first, a.second) : 1.,
                        b.first, b.second, b.second > 0 ? TMath::Prob(b.first, b.second) : 1.,
                        m.first, m.second, m.second > 0 ? TMath::Prob(m.first, m.second) : 1.,
                        x.first, x.second, x.second > 0 ? TMath::Prob(x.first, x.second) : 1., dev, unmeas ? " *" : "");
        }
        tab << Form("    panel total: %.1f/%d (p=%.3f) | %.1f/%d (p=%.3f) | %.1f/%d (p=%.3f) | %.1f/%d (p=%.3f)\n",
                    c2rb, nrb, nrb > 0 ? TMath::Prob(c2rb, nrb) : 1., c2tb, ntb, ntb > 0 ? TMath::Prob(c2tb, ntb) : 1.,
                    c2mc, nmc, nmc > 0 ? TMath::Prob(c2mc, nmc) : 1., c2dr, ndr, ndr > 0 ? TMath::Prob(c2dr, ndr) : 1.);
        if (p > 0) { chi2_rb_tot += c2rb; ndf_rb_tot += nrb; chi2_tb_tot += c2tb; ndf_tb_tot += ntb; chi2_mc_tot += c2mc; ndf_mc_tot += nmc;
                     chi2_dr_tot += c2dr; ndf_dr_tot += ndr; }
        }  // !coarse_eta

        // ---- figures: dsigma/deta (lumi factor already in the weights) ----
        for (TH1D* h : {d_rb, d_tr, d_tb, mT_d, mR_d}) h->Scale(1.0, "width");
        PanelData A; A.label = labels[p];
        A.s = {{d_rb, "trigger + reco-binned #varepsilon_{reco}", kBlack, 20},
               {d_tr, "trigger only", kBlue + 1, 21},   // filled: open square (25) renders randomly black in batch PNG (ROOT 6.34)
               {d_tb, "trigger + #varepsilon_{reco} #times SF (crossx)", kRed + 1, 24}};
        A.ratios = {{Ratio(d_tb, d_rb, "rA" + tg), "", kBlack, 20}};
        P1.push_back(A);
        PanelData B; B.label = labels[p];
        B.s = {{d_rb, "data, trigger + reco-binned #varepsilon_{reco}", kBlack, 20},
               {mT_d, "Pythia8 truth", kRed + 1, 1, true},
               {mR_d, "Pythia8 reco (truth in signal region) / reco-binned #varepsilon_{reco}", kBlue + 1, 1, true}};
        B.ratios = {{r_dm, "", kRed + 1, 20}, {r_dr, "", kBlue + 1, 24}};
        P2.push_back(B);
        TH1D* mTc = static_cast<TH1D*>(mT->Clone(("mTc" + tg).c_str())); mTc->SetDirectory(nullptr);
        TH1D* mRcc = static_cast<TH1D*>(mRc->Clone(("mRcc" + tg).c_str())); mRcc->SetDirectory(nullptr);
        if (coarse_eta) { mTc->Scale(1.0, "width"); mRcc->Scale(1.0, "width"); }   // unequal coarse widths
        PanelData C; C.label = labels[p];
        C.s = {{mRcc, "reco level (truth in signal region) / reco-binned #varepsilon_{reco}", kBlack, 20},
               {mTc, "truth level", kRed + 1, 1, true}};
        C.ratios = {{r_mc, "", kBlack, 20}};
        P3.push_back(C);
    }
    tab << Form("\nSUM over the 8 pT panels: within-cell chi2/ndf data(reco-binned) %.1f/%d (p=%.3g), data(applied eps^{reco,SF}) %.1f/%d (p=%.3g), "
                "MC corrected R vs T %.1f/%d (p=%.3g), data(reco-binned) vs MC corrected R %.1f/%d (p=%.3g)\n"
                "MC closure: max |sum(corrR)/sum(T) - 1| per coarse cell: pT panels %.3e, all pT %.3e; cells containing a 3D cell "
                "without eps measurement (marked *) %.3e\n",
                chi2_rb_tot, ndf_rb_tot, TMath::Prob(chi2_rb_tot, ndf_rb_tot), chi2_tb_tot, ndf_tb_tot, TMath::Prob(chi2_tb_tot, ndf_tb_tot),
                chi2_mc_tot, ndf_mc_tot, TMath::Prob(chi2_mc_tot, ndf_mc_tot), chi2_dr_tot, ndf_dr_tot, TMath::Prob(chi2_dr_tot, ndf_dr_tot),
                max_closure_dev, max_closure_dev_incl, max_closure_dev_unmeas);
    tab.close(); ctab.close();
    if (!coarse_eta) std::cout << Form("[closure] max |sum(corrR)/sum(T)-1| per coarse cell: pT panels %.3e, all pT %.3e, cells w/ unmeasured 3D cell %.3e\n",
                      max_closure_dev, max_closure_dev_incl, max_closure_dev_unmeas);
    if (!coarse_eta) std::cout << Form("[chi2] within-cell sum over pT panels: data rb vs T %.1f/%d, data applied vs T %.1f/%d, MC corrR vs T %.1f/%d, data rb vs MC corrR %.1f/%d\n",
                      chi2_rb_tot, ndf_rb_tot, chi2_tb_tot, ndf_tb_tot, chi2_mc_tot, ndf_mc_tot, chi2_dr_tot, ndf_dr_tot);

    const std::string hdr_data = "pp 2024 data, #sqrt{s} = 5.36 TeV, 2mu4";
    const std::string sel = "opposite-sign pairs, signal region";
    const std::string sfx = coarse_eta ? "_coarse_eta_" : "_";
    DrawFigure(outdir + "data_corrected" + sfx + wp + "WP.png", hdr_data, wp_text, sel, P1, coarse,
               "d#sigma/d#eta^{pair} [pb]", "crossx/reco-binned", 0.8, 1.2);
    DrawFigure(outdir + "data_vs_mc_truth" + sfx + wp + "WP.png", hdr_data, wp_text,
               sel + "; Pythia8 normalised to the data integral of each panel", P2, coarse,
               "d#sigma/d#eta^{pair} [pb]", "data / MC", coarse_eta ? 0.5 : 0.0, coarse_eta ? 1.5 : 2.0);
    DrawFigure(outdir + "mc_closure" + sfx + wp + "WP.png", "Pythia8 full simulation, pp 2024 conditions", wp_text,
               "all opposite-sign pairs, signal region", P3, coarse, coarse_eta ? "d(weighted pairs)/d#eta^{pair}" : "weighted pairs / 0.1", "reco / truth", 0.0, 2.0);
    std::cout << "[done] -> " << outdir << std::endl;
}
