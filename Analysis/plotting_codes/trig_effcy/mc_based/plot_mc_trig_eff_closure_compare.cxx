// =============================================================================
// plot_mc_trig_eff_closure_compare.cxx
//
// THE FOUR CELL-GROUPING APPROACHES, OVERLAID IN ONE FIGURE
// (docs/tracking/mc_trigeff_dr_binning_approaches.md §PP-4; inputs from FillMCTrigEffClosure.cxx,
//  one file per approach).
//
// The Step-3 dR correction eps_dR is fitted per (pair pT, pair eta) cell, and the four approaches
// differ ONLY in how that plane is partitioned before the fit:
//
//   A  8 pair-pT x 9 pair-eta   the un-merged reference                    kRed
//   B  7 x 9                    the last two pair-pT bins combined         kBlue
//   C  8 x 3                    pair eta combined into the 3 detector      kGreen+2
//                               regions (negative endcap, barrel,
//                               positive endcap)
//   D  7 x 3                    both                                       kMagenta
//
// Everything else is identical between them -- the same sample, the same selection, the same
// applied eps_MC, and the SAME delivered cascade (expo fit -> polyu fit -> interpolation ->
// measured values). So a difference between the five lines below is a difference between the
// GROUPINGS, which is the question the figure exists to answer: which of the four corrects the
// pp24 cross-section most accurately.
//
// SAME BINNING, ALWAYS. Every panel is drawn on the pp24 CROSS-SECTION's binning
// (ParamsSet::pT_bins_150, 15 log bins 8-150 GeV, x the 9
// CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap panels) -- that is what makes the
// overlay legitimate, and it is enforced here rather than assumed: the four denominators are
// compared bin by bin and a disagreement THROWS.
//
// THE NUMBERS ARE THE SAME OBJECT the per-approach figures show: this macro re-reads each
// approach's own closure file and re-forms the same ratio with the same conditional error, so a
// point here must equal the point in that approach's own subdirectory.
//
// TWO FIGURES PER SAMPLE VERSION:
//   closure_compare_pair_pt_*.png         the 9 pair-eta panels, 5 series each
//   closure_compare_nonclosure_chi2_*.png the same closure ratios compressed to ONE panel --
//                                         chi^2/ndof(pT) = [sum over pair eta of (C-1)^2/sigma_C^2]
//                                         / n_dof, one curve per approach. A PROPER chi^2 (each
//                                         term IS divided by its statistical variance, per the
//                                         user 2026-09-03 -- the earlier un-normalized sum of
//                                         (C-1)^2 let noisy, statistics-starved cells dominate the
//                                         total on equal footing with genuinely mis-corrected
//                                         ones); see the block that builds it.
//
// Usage (from Analysis/plotting_codes/trig_effcy/mc_based/):
//   root -l -b -q 'plot_mc_trig_eff_closure_compare.cxx+("pp_full", true)'
//   root -l -b -q 'plot_mc_trig_eff_closure_compare.cxx+("pp_full", false)'   // Medium
// =============================================================================

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TNamed.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include "dr_correction_sample_cfg.h"
#include "dr_correction_apply.h"
#include "dr_correction_ratio.h"
#include "../../../Utilities/CommonLogYRange.h"
#include "../../../Utilities/MCTrigEffPairPtBinning.h"
#include "../../../RDFBasedHistFilling/CommonEffcyConfig.h"
#include "../../../Utilities/PairEtaPanelBins.h"

namespace {

// One approach. `legend_shape` is the PHYSICS description of the grouping; the cell COUNTS in the
// legend are read from each file's own provenance stamp (below) rather than typed here, so a
// mislabelled figure is impossible.
struct Approach {
    std::string mode;            // plateau-mode token; never drawn
    std::string legend_shape;    // what was combined, in physics words
    Color_t     colour;
    Style_t     marker;
};
// Colours fixed by the user (2026-08-24): kRed / kBlue / kGreen+2 / kMagenta, in this order.
const std::vector<Approach> kApproaches = {
    {"nocorr",                  "no combination",                             kRed,      20},
    {"nocorr_ptmerge",          "top two p_{T}^{pair} cells combined",        kBlue,     22},
    {"nocorr_etamerge",         "#eta^{pair} combined into endcap / barrel / endcap",
                                                                              kGreen + 2, 23},
    {"nocorr_etamerge_ptmerge", "both combined",                              kMagenta,  33},
};
// The DELIVERED series key written by FillMCTrigEffClosure.cxx -- the same in every approach.
const char* kCascadeKey = "cascade";

const Color_t kUncorrColour = kBlack;
const Style_t kUncorrMarker = 21;

struct VersionCfg { std::string key, file, file_chi2, headline; };
const std::vector<VersionCfg> kVersions = {
    {"all_os", "closure_compare_pair_pt_all_opposite_sign",
               "closure_compare_nonclosure_chi2_all_opposite_sign",
     "all opposite-sign muon pairs"},
    {"signal", "closure_compare_pair_pt_single_b_signal_cuts",
               "closure_compare_nonclosure_chi2_single_b_signal_cuts",
     "opposite-sign pairs in the single-b signal region"},
};

std::string EtaLabel(const std::pair<float, float>& r)
{
    return Form("%.1f < #eta^{pair} < %.1f", r.first, r.second);
}

template <typename T>
T* Get(TFile* f, const std::string& n)
{
    T* o = dynamic_cast<T*>(f->Get(n.c_str()));
    if (!o) throw std::runtime_error("plot_mc_trig_eff_closure_compare: missing '" + n + "' in "
                                     + f->GetName());
    return o;
}

TH1D* Row(TH2D* h, int iz, const std::string& nm)
{
    TH1D* p = h->ProjectionX(nm.c_str(), iz, iz, "e");
    p->SetDirectory(nullptr);
    return p;
}

// "eps_dR CELLS: 8 pair pT x 9 pair eta" -> "8 #times 9 cells". Read from the file the points came
// from, so the legend cannot describe a grouping the histogram was not built with.
std::string CellCountFromProvenance(TFile* f)
{
    auto* prov = dynamic_cast<TNamed*>(f->Get("provenance"));
    if (!prov) return "";
    const TString t = prov->GetTitle();
    const Ssiz_t i = t.Index("eps_dR CELLS: ");
    if (i < 0) return "";
    const Ssiz_t j = t.Index(" (", i);
    if (j <= i) return "";
    TString cells = t(i + 14, j - i - 14);           // "8 pair pT x 9 pair eta"
    cells.ReplaceAll(" pair pT x ", " #times ");
    cells.ReplaceAll(" pair eta", " cells");
    return std::string(cells.Data());
}

}  // namespace

// =============================================================================
void plot_mc_trig_eff_closure_compare(const std::string& sample = "pp_full",
                                      bool use_tight_wp = true)
{
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);

    const DrCorrSample cfg = GetDrCorrSample(sample, use_tight_wp);
    const std::string wp_suf  = DrCorrWpSuffix(use_tight_wp);
    const std::string wp_text = use_tight_wp ? "Tight muons" : "Medium muons";

    static const CommonEffcyConfig ecfg{};
    const auto& eta_ranges = ecfg.pair_eta_proj_ranges_coarse_incl_gap;
    const int neta = static_cast<int>(eta_ranges.size());
    const int ncol = static_cast<int>(std::ceil(std::sqrt((double)neta)));
    const int nrow = static_cast<int>(std::ceil((double)neta / ncol));

    const std::string outdir = cfg.out_base + "closure/approach_comparison/";
    gSystem->mkdir(outdir.c_str(), kTRUE);

    // ---------------- open the four approach files ------------------------------------------
    std::vector<TFile*>      fin(kApproaches.size(), nullptr);
    std::vector<std::string> legend(kApproaches.size());
    for (size_t a = 0; a < kApproaches.size(); ++a) {
        const std::string p = cfg.mc_dir + "mc_trig_eff_closure_" + cfg.mc_label + wp_suf
                            + MCTrigEffPairPt::FileSuffix()
                            + DrCorrPlateauModeTag(kApproaches[a].mode) + ".root";
        fin[a] = TFile::Open(p.c_str(), "READ");
        if (!fin[a] || fin[a]->IsZombie())
            throw std::runtime_error("plot_mc_trig_eff_closure_compare: cannot open " + p
                                     + " -- run FillMCTrigEffClosure with plateau mode '"
                                     + kApproaches[a].mode + "' first");
        const std::string cells = CellCountFromProvenance(fin[a]);
        legend[a] = (cells.empty() ? std::string() : cells + ":  ") + kApproaches[a].legend_shape;
        std::cout << "  " << kApproaches[a].mode << " <- " << p << "  [" << legend[a] << "]"
                  << std::endl;
    }

    // ---------------- build everything first (the shared ratio range needs it all) -----------
    struct VersionData {
        TH2D* h_den = nullptr;
        std::vector<TH2D*> h_num;                       // [approach]
        std::vector<TH1D*> spec_unc;                    // [pair-eta bin]
        std::vector<std::vector<TH1D*>> spec_cor, ratio;  // [pair-eta bin][approach]
    };
    std::vector<VersionData> vdata(kVersions.size());
    std::vector<double> ratio_pts;

    for (size_t iv = 0; iv < kVersions.size(); ++iv) {
        const VersionCfg& V = kVersions[iv];
        VersionData&      D = vdata[iv];
        const std::string pre = "h_closure_" + V.key + "_";

        // THE DENOMINATOR IS THE SAME OBJECT IN ALL FOUR FILES -- same sample, same selection,
        // same crossx binning. Checked bin by bin rather than assumed: if it ever differs, the
        // four lines are not four corrections of ONE spectrum and the overlay is meaningless.
        D.h_den = Get<TH2D>(fin[0], pre + "den");
        for (size_t a = 1; a < kApproaches.size(); ++a) {
            TH2D* d = Get<TH2D>(fin[a], pre + "den");
            if (d->GetNbinsX() != D.h_den->GetNbinsX() || d->GetNbinsY() != D.h_den->GetNbinsY())
                throw std::runtime_error("plot_mc_trig_eff_closure_compare: approach '"
                    + kApproaches[a].mode + "' was filled on a DIFFERENT histogram binning than '"
                    + kApproaches[0].mode + "' -- the overlay would compare different cells");
            for (int bx = 1; bx <= d->GetNbinsX(); ++bx)
                for (int by = 1; by <= d->GetNbinsY(); ++by) {
                    const double u = D.h_den->GetBinContent(bx, by), w = d->GetBinContent(bx, by);
                    if (std::fabs(u - w) > 1e-6 * std::max(1.0, std::fabs(u)))
                        throw std::runtime_error(Form("plot_mc_trig_eff_closure_compare: the "
                            "no-trigger denominators of '%s' and '%s' differ in bin (%d,%d): "
                            "%g vs %g -- the four runs did not see the same sample",
                            kApproaches[0].mode.c_str(), kApproaches[a].mode.c_str(), bx, by, u, w));
                }
        }
        if (D.h_den->GetNbinsY() != neta)
            throw std::runtime_error("plot_mc_trig_eff_closure_compare: the histograms carry "
                + std::to_string(D.h_den->GetNbinsY()) + " pair-eta bins but CommonEffcyConfig has "
                + std::to_string(neta) + " -- stale input?");
        // BIN COUNT IS NOT ENOUGH. Panel `iz` is drawn with the label EtaLabel(eta_ranges[iz-1])
        // read LIVE from CommonEffcyConfig, while its content is Y-bin `iz` of a file written
        // earlier. A count check passes even when the two describe different ranges -- which is
        // exactly what happened on 2026-09-07, when the coarse pair-eta outer bins moved
        // +-(2.0,2.4) -> +-(2.0,2.2) to follow the pair-level gap cut: every closure file on disk
        // still had 9 bins over +-2.4 and would have been relabelled silently. Check the EDGES,
        // per panel, through the shared guard (muon_gap_cuts_acceptance.md F18).
        PairEtaPanels::CheckPanelsMatchFiducialCut();
        for (int iz = 1; iz <= neta; ++iz)
            PairEtaPanels::CheckAxisAligned(D.h_den->GetYaxis(), eta_ranges[iz - 1], iz, iz,
                                            "plot_mc_trig_eff_closure_compare");

        D.h_num.assign(kApproaches.size(), nullptr);
        std::vector<TH2D*> h_A(kApproaches.size(), nullptr), h_B(kApproaches.size(), nullptr);
        for (size_t a = 0; a < kApproaches.size(); ++a) {
            D.h_num[a] = Get<TH2D>(fin[a], pre + "num_epsmc_"  + kCascadeKey);
            h_A[a]     = Get<TH2D>(fin[a], pre + "numA_epsmc_" + kCascadeKey);
            h_B[a]     = Get<TH2D>(fin[a], pre + "numB_epsmc_" + kCascadeKey);
        }

        D.spec_unc.assign(neta, nullptr);
        D.spec_cor.assign(neta, std::vector<TH1D*>(kApproaches.size(), nullptr));
        D.ratio   .assign(neta, std::vector<TH1D*>(kApproaches.size(), nullptr));
        std::vector<TH1*> for_range;

        for (int iz = 1; iz <= neta; ++iz) {
            TH1D* den = Row(D.h_den, iz, Form("cmp_%s_den_eta%d", V.key.c_str(), iz));
            for (size_t a = 0; a < kApproaches.size(); ++a) {
                TH1D* num = Row(D.h_num[a], iz, Form("cmp_%s_num%zu_eta%d", V.key.c_str(), a, iz));
                TH1D* A   = Row(h_A[a],     iz, Form("cmp_%s_A%zu_eta%d",   V.key.c_str(), a, iz));
                TH1D* B   = Row(h_B[a],     iz, Form("cmp_%s_B%zu_eta%d",   V.key.c_str(), a, iz));
                // Formed exactly as in plot_mc_trig_eff_closure.cxx -- same order (ratio BEFORE
                // the width scaling, because A and B are sums of squared weights and do not scale
                // the same way) and the same single conditional-error implementation. That is what
                // makes a point here identical to the same point in the approach's own figure.
                TH1D* r = static_cast<TH1D*>(num->Clone(Form("cmp_%s_ratio%zu_eta%d",
                                                             V.key.c_str(), a, iz)));
                r->SetDirectory(nullptr);
                r->Divide(den);
                SetConditionalRatioErrors(r, den, A, B);
                delete A; delete B;
                D.ratio[iz - 1][a] = r;
                for (int b = 1; b <= r->GetNbinsX(); ++b)
                    if (r->GetBinContent(b) > 0.) ratio_pts.push_back(r->GetBinContent(b));
                num->Scale(1.0, "width");
                D.spec_cor[iz - 1][a] = num;
                for_range.push_back(num);
            }
            den->Scale(1.0, "width");
            D.spec_unc[iz - 1] = den;
            for_range.push_back(den);
        }
        ApplyCommonLogYRange(for_range);
    }

    // ---------------- the shared ratio range -------------------------------------------------
    double rmin = 1.0, rmax = 1.0;
    if (!ratio_pts.empty()) {
        rmin = *std::min_element(ratio_pts.begin(), ratio_pts.end());
        rmax = *std::max_element(ratio_pts.begin(), ratio_pts.end());
    }
    rmin = std::min(rmin, 1.0); rmax = std::max(rmax, 1.0);
    if (rmax - rmin < 0.10) { const double m = 0.5 * (rmin + rmax); rmin = m - 0.05; rmax = m + 0.05; }
    const double rpad = 0.08 * (rmax - rmin);
    rmin -= rpad; rmax += rpad;

    // ---------------- draw --------------------------------------------------------------------
    for (size_t iv = 0; iv < kVersions.size(); ++iv) {
        const VersionCfg& V = kVersions[iv];
        VersionData&      D = vdata[iv];

        int n_offscale = 0;
        for (int iz = 0; iz < neta; ++iz)
            for (size_t a = 0; a < kApproaches.size(); ++a) {
                TH1D* r = D.ratio[iz][a];
                for (int b = 1; b <= r->GetNbinsX(); ++b) {
                    const double y = r->GetBinContent(b);
                    if (y > 0. && (y < rmin || y > rmax)) ++n_offscale;
                }
            }

        // Five series need three legend rows, so the reserved header strip is taller than the
        // single-approach figure's.
        const double head = 0.21;
        TCanvas c(("c_cmp_" + V.key).c_str(), "", 1800, 1800);

        for (int iz = 1; iz <= neta; ++iz) {
            c.cd();   // a new TPad is adopted by gPad -- without this every panel nests in the last
            const int col = (iz - 1) % ncol, row = (iz - 1) / ncol;
            const double x1 = col / (double)ncol, x2 = (col + 1) / (double)ncol;
            const double ytop = 1.0 - head - row * (1.0 - head) / nrow;
            const double ybot = 1.0 - head - (row + 1) * (1.0 - head) / nrow;
            const double ysplit = ybot + 0.36 * (ytop - ybot);

            auto* pu = new TPad(Form("cmp_pu_%s_%d", V.key.c_str(), iz), "", x1, ysplit, x2, ytop);
            auto* pl = new TPad(Form("cmp_pl_%s_%d", V.key.c_str(), iz), "", x1, ybot,   x2, ysplit);
            for (TPad* p : {pu, pl}) {
                p->SetLeftMargin(0.235); p->SetRightMargin(0.02);
                p->SetLogx(1);
                p->SetTicks(1, 1);
            }
            pu->SetTopMargin(0.04); pu->SetBottomMargin(0.02);
            pl->SetTopMargin(0.02); pl->SetBottomMargin(0.32);
            pu->Draw(); pl->Draw();

            pu->cd();
            pu->SetLogy(1);
            TH1D* du = D.spec_unc[iz - 1];
            du->SetTitle("");
            du->GetYaxis()->SetTitle("d#sigma/dp_{T}^{pair} [nb/GeV]");
            du->GetYaxis()->SetTitleSize(0.070); du->GetYaxis()->SetLabelSize(0.060);
            du->GetYaxis()->SetTitleOffset(1.55);
            du->GetXaxis()->SetLabelSize(0.0);
            du->SetMarkerColor(kUncorrColour); du->SetLineColor(kUncorrColour);
            du->SetMarkerStyle(kUncorrMarker); du->SetMarkerSize(1.0);
            du->Draw("PE");
            for (size_t a = 0; a < kApproaches.size(); ++a) {
                TH1D* h = D.spec_cor[iz - 1][a];
                h->SetMarkerColor(kApproaches[a].colour); h->SetLineColor(kApproaches[a].colour);
                h->SetMarkerStyle(kApproaches[a].marker); h->SetMarkerSize(1.0);
                h->Draw("PE SAME");
            }
            auto* tl = new TLatex();
            tl->SetNDC(); tl->SetTextFont(42); tl->SetTextSize(0.068);
            tl->DrawLatex(0.28, 0.09, EtaLabel(eta_ranges[iz - 1]).c_str());

            pl->cd();
            for (size_t a = 0; a < kApproaches.size(); ++a) {
                TH1D* r = D.ratio[iz - 1][a];
                r->SetTitle("");
                r->SetMarkerColor(kApproaches[a].colour); r->SetLineColor(kApproaches[a].colour);
                r->SetMarkerStyle(kApproaches[a].marker); r->SetMarkerSize(1.0);
                r->GetYaxis()->SetRangeUser(rmin, rmax);
                r->GetYaxis()->SetTitle("corrected / no trigger");
                r->GetYaxis()->SetNdivisions(505);
                r->GetYaxis()->SetTitleSize(0.100); r->GetYaxis()->SetLabelSize(0.095);
                r->GetYaxis()->SetTitleOffset(1.08);
                r->GetXaxis()->SetTitle("p_{T}^{pair} [GeV]");
                r->GetXaxis()->SetTitleSize(0.125); r->GetXaxis()->SetLabelSize(0.105);
                r->GetXaxis()->SetTitleOffset(1.00);
                r->Draw(a == 0 ? "PE" : "PE SAME");
            }
            auto* ln = new TLine(D.spec_unc[iz - 1]->GetXaxis()->GetXmin(), 1.0,
                                 D.spec_unc[iz - 1]->GetXaxis()->GetXmax(), 1.0);
            ln->SetLineStyle(2); ln->SetLineColor(kGray + 2);
            ln->Draw();
        }

        // ---- header strip -------------------------------------------------------------------
        c.cd();
        auto* hd = new TLatex();
        hd->SetNDC(); hd->SetTextFont(42);
        hd->SetTextSize(0.0175);
        hd->DrawLatex(0.030, 0.984,
                      (cfg.sample_text + ",  " + wp_text + ",  " + V.headline).c_str());
        hd->SetTextSize(0.0145);
        hd->DrawLatex(0.030, 0.884,
                      "w = 1 / [#varepsilon(p_{T,1},q#eta_{1}) "
                      "#varepsilon(p_{T,2},q#eta_{2}) #varepsilon_{#DeltaR}(#DeltaR)],"
                      "   #varepsilon = single-muon mu4 efficiency from the MC turn-on");
        hd->DrawLatex(0.030, 0.864,
                      Form("#varepsilon_{#DeltaR}(#DeltaR) = f(#DeltaR)/C for #DeltaR < %g and 1 "
                           "above;  f = the exponential C + A e^{-(#DeltaR/#lambda)^{p}} where its "
                           "fit is accepted, the polynomial C + u^{2}(a_{2}+a_{3}u+a_{4}u^{2}) "
                           "where it is not,",
                           DrCorrectionEvaluator::kDrMax));
        hd->DrawLatex(0.030, 0.844,
                      "the linear interpolation of the measured points where neither is, and those "
                      "measured values themselves where none is #minus IDENTICALLY in all four "
                      "series, which differ only in the cells f is fitted in");
        if (n_offscale > 0)
            hd->DrawLatex(0.030, 0.824,
                          Form("%d ratio point(s) outside the lower-pad range", n_offscale));

        auto* leg = new TLegend(0.030, 0.900, 0.990, 0.976);
        leg->SetNColumns(2);
        leg->SetBorderSize(0); leg->SetFillStyle(0);
        leg->SetTextFont(42); leg->SetTextSize(0.0145);
        leg->AddEntry(D.spec_unc[0], "no trigger requirement", "PE");
        for (size_t a = 0; a < kApproaches.size(); ++a)
            leg->AddEntry(D.spec_cor[0][a],
                          ("2mu4, corrected  #minus  " + legend[a]).c_str(), "PE");
        leg->Draw();

        const std::string png = outdir + V.file + ".png";
        c.SaveAs(png.c_str());
        std::cout << "  wrote " << png << std::endl;

        // The inclusive number per approach -- it MUST equal the one printed by
        // plot_mc_trig_eff_closure.cxx for that approach's own figure.
        std::cout << "  inclusive closure (" << V.key << "):";
        const double d = D.h_den->Integral(1, D.h_den->GetNbinsX(), 1, neta);
        for (size_t a = 0; a < kApproaches.size(); ++a)
            std::cout << "  " << kApproaches[a].mode << " = "
                      << (d > 0 ? D.h_num[a]->Integral(1, D.h_num[a]->GetNbinsX(), 1, neta) / d
                                : -1.);
        std::cout << std::endl;

        // ===== the pair-eta-summed chi^2/ndof non-closure, chi^2(pT)/n_dof =====================
        //
        //     chi^2(p_T^pair)      = sum over pair-eta bins of (C - 1)^2 / sigma_C^2
        //     chi^2/ndof(p_T^pair) = chi^2(p_T^pair) / n_dof(p_T^pair)
        //
        // over exactly the closure ratios C the lower pads above draw, in the same crossx
        // pair-pT binning. It compresses the 9 pair-eta panels of one approach into one curve so
        // the four approaches can be read against each other bin by bin.
        //
        // THIS IS NOW A PROPER chi^2 (user, 2026-09-03 -- supersedes the un-normalized D^2 of
        // 2026-08-25, "sum (C-1)^2, not divided by errors"). Each term IS divided by its
        // statistical variance sigma_C^2, so a noisy, statistics-starved cell no longer
        // contributes on equal footing with a genuinely mis-corrected one: R2's own reviewer
        // check on the un-normalized version ("bin 14 alone is 6.998 +/- 10.8 -- a ~0.6 sigma
        // effect") is exactly the failure this fixes. chi^2/ndof ~ 1 is the closure reference;
        // its bar is the null-hypothesis width sqrt(2/ndof), NOT a propagated point error -- a
        // chi^2 statistic is compared to its own distribution, not to itself with a Gaussian bar.
        //
        // A pair-eta cell enters the sum only where the no-trigger DENOMINATOR is non-empty AND
        // sigma_C > 0: an empty cell has C = 0, sigma_C = 0 by construction (dr_correction_ratio.h
        // returns both as 0 for D <= 0) and is excluded rather than dividing by zero. A cell with
        // a filled denominator but nothing passing the trigger has C = 0 and sigma_C =
        // max(0,1)/n_eff > 0 (the k=0 boundary case in SetConditionalRatioErrors), so it correctly
        // contributes a large but finite term -- real total non-closure, not a NaN.
        //
        // ndof(p_T) does not depend on the approach: the four denominators were checked bin by
        // bin to be identical above, so the same cells are skipped in every approach and the four
        // raw chi^2 sums share one ndof per pair-pT bin.
        std::vector<TH1D*> chi2h(kApproaches.size(), nullptr);    // plotted: chi^2/ndof per pT bin
        std::vector<TH1D*> chi2raw(kApproaches.size(), nullptr); // chi^2 (unnormalized by ndof)
        for (size_t a = 0; a < kApproaches.size(); ++a) {
            chi2h[a] = static_cast<TH1D*>(D.ratio[0][a]->Clone(
                        Form("cmp_%s_chi2_%zu", V.key.c_str(), a)));
            chi2h[a]->SetDirectory(nullptr);
            chi2h[a]->Reset("ICES");
            chi2raw[a] = static_cast<TH1D*>(chi2h[a]->Clone(Form("cmp_%s_chi2raw_%zu",
                                                                  V.key.c_str(), a)));
            chi2raw[a]->SetDirectory(nullptr);
        }
        const int nptx = D.h_den->GetNbinsX();
        std::vector<int> ncell(nptx + 1, 0);
        double redmax = 0., redminpos = 1e300;
        for (int bx = 1; bx <= nptx; ++bx)
            for (size_t a = 0; a < kApproaches.size(); ++a) {
                double chi2 = 0.;
                int ncontrib = 0;
                for (int iz = 1; iz <= neta; ++iz) {
                    if (D.h_den->GetBinContent(bx, iz) <= 0.) continue;
                    const TH1D* r = D.ratio[iz - 1][a];
                    const double dev = r->GetBinContent(bx) - 1.0;
                    const double er  = r->GetBinError(bx);
                    if (er <= 0.) continue;   // no measured variance -> cannot enter a chi^2 term
                    chi2 += (dev * dev) / (er * er);
                    ++ncontrib;
                }
                chi2raw[a]->SetBinContent(bx, chi2);
                if (a == 0) ncell[bx] = ncontrib;   // the same in all four -- same denominators
                const double reduced = (ncontrib > 0) ? chi2 / ncontrib : 0.;
                const double err     = (ncontrib > 0) ? std::sqrt(2.0 / ncontrib) : 0.;
                chi2h[a]->SetBinContent(bx, reduced);
                chi2h[a]->SetBinError(bx, err);
                if (ncontrib > 0) {
                    if (reduced > redmax) redmax = reduced;
                    if (reduced > 0. && reduced < redminpos) redminpos = reduced;
                }
            }

        // Log y when the curves span more than ~2 decades, linear from 0 otherwise -- on a linear
        // axis the low-pT half of a many-decade curve is a flat line on the frame and unreadable.
        // Either way the range is widened, if needed, so the chi^2/ndof = 1 closure reference is
        // always on the frame.
        const bool chi2_logy = (redminpos < 1e299) && (redmax > 0.) && (redmax / redminpos > 100.);

        TCanvas c2(("c_chi2_" + V.key).c_str(), "", 1100, 850);
        c2.SetLeftMargin(0.145); c2.SetRightMargin(0.035);
        c2.SetTopMargin(0.105);  c2.SetBottomMargin(0.115);
        c2.SetLogx(1);
        c2.SetLogy(chi2_logy ? 1 : 0);
        c2.SetTicks(1, 1);

        for (size_t a = 0; a < kApproaches.size(); ++a) {
            TH1D* h = chi2h[a];
            h->SetTitle("");
            h->SetMarkerColor(kApproaches[a].colour); h->SetLineColor(kApproaches[a].colour);
            h->SetMarkerStyle(kApproaches[a].marker); h->SetMarkerSize(1.3);
            h->GetXaxis()->SetTitle("p_{T}^{pair} [GeV]");
            h->GetYaxis()->SetTitle("#chi^{2}_{#eta^{pair}} / n_{dof}");
            h->GetXaxis()->SetTitleSize(0.045); h->GetXaxis()->SetLabelSize(0.040);
            h->GetYaxis()->SetTitleSize(0.045); h->GetYaxis()->SetLabelSize(0.040);
            h->GetYaxis()->SetTitleOffset(1.50);
            if (chi2_logy) h->GetYaxis()->SetRangeUser(std::min(0.4 * redminpos, 0.5),
                                                        std::max(3.0 * redmax, 2.0));
            else           h->GetYaxis()->SetRangeUser(0.0, std::max(1.25 * redmax, 1.3));
            h->Draw(a == 0 ? "PE" : "PE SAME");
        }

        // Perfect-closure reference: chi^2/ndof = 1, the value the whole figure is measured
        // against (same dashed-red convention as the closure ratio pads above).
        TLine* lref = new TLine(chi2h[0]->GetXaxis()->GetXmin(), 1.0,
                                 chi2h[0]->GetXaxis()->GetXmax(), 1.0);
        lref->SetLineStyle(2); lref->SetLineColor(kRed + 1); lref->Draw("SAME");

        auto* h2 = new TLatex();
        h2->SetNDC(); h2->SetTextFont(42); h2->SetTextSize(0.0235);
        h2->DrawLatex(0.100, 0.957,
                      (cfg.sample_text + ",  " + wp_text + ",  " + V.headline).c_str());
        h2->SetTextSize(0.0215);
        h2->DrawLatex(0.100, 0.922,
                      // TLatex leaves a visible gap in "1/#varepsilon" -- the inverse power is
                      // both correct and renders cleanly.
                      "C(p_{T}^{pair}, #eta^{pair}) = "
                      "[2mu4 yield weighted by (#varepsilon_{trig}^{pair})^{-1}] / "
                      "[yield with no trigger requirement]");

        auto* leg2 = new TLegend(0.170, 0.630, 0.680, 0.880);
        leg2->SetBorderSize(0); leg2->SetFillStyle(0);
        leg2->SetTextFont(42); leg2->SetTextSize(0.028);
        for (size_t a = 0; a < kApproaches.size(); ++a)
            leg2->AddEntry(chi2h[a], legend[a].c_str(), "PE");
        leg2->Draw();

        const std::string png2 = outdir + V.file_chi2 + ".png";
        c2.SaveAs(png2.c_str());
        std::cout << "  wrote " << png2 << std::endl;

        // The pair-pT-summed total per approach -- the single reduced chi^2 the figure is a
        // differential view of -- and the number of pair-eta cells (= ndof) each pair-pT bin
        // actually summed over, which is NOT constant across the axis (the high-pT bins are empty
        // in some pair-eta panels) and is therefore needed to read the curve.
        //
        // AND THE SAME TOTAL SPLIT AT THE PAIR-pT MERGE EDGE. The merge fuses the top TWO of the
        // 8 eps_dR pair-pT cells, so below that edge a pair gets an identical correction whether
        // or not the merge is on: A and B coincide bin by bin, and so do C and D. Every difference
        // the pair-pT merge makes therefore lives in the few crossx bins above the edge. The split
        // bin is FOUND, not typed: it is the first bin where a merged-pT approach's raw chi^2
        // departs from its un-merged twin's.
        int bsep = nptx + 1;
        for (int bx = 1; bx <= nptx && bsep > nptx; ++bx)
            for (size_t a = 0; a + 1 < kApproaches.size(); a += 2)   // (A,B) then (C,D)
                if (std::fabs(chi2raw[a + 1]->GetBinContent(bx) - chi2raw[a]->GetBinContent(bx))
                        > 1e-9 * std::max(1.0, chi2raw[a]->GetBinContent(bx))) { bsep = bx; break; }

        double ndofAll = 0.;
        for (int bx = 1; bx <= nptx; ++bx) ndofAll += ncell[bx];
        std::cout << "  chi^2/ndof summed over pair pT (" << V.key << "), ndof=" << ndofAll << ":";
        for (size_t a = 0; a < kApproaches.size(); ++a) {
            const double tot = chi2raw[a]->Integral(1, nptx);   // bin counts, not a density
            std::cout << "  " << kApproaches[a].mode << " = "
                      << (ndofAll > 0 ? tot / ndofAll : 0.) << " (chi2=" << tot << ")";
        }
        if (bsep <= nptx) {
            double ndofLo = 0., ndofHi = 0.;
            for (int bx = 1; bx < bsep; ++bx) ndofLo += ncell[bx];
            for (int bx = bsep; bx <= nptx; ++bx) ndofHi += ncell[bx];
            std::cout << "\n    of which bins 1-" << bsep - 1 << " (p_T^pair < "
                      << chi2raw[0]->GetXaxis()->GetBinLowEdge(bsep)
                      << " GeV, below the pair-pT merge edge, ndof=" << ndofLo << "):";
            for (size_t a = 0; a < kApproaches.size(); ++a) {
                const double tot = chi2raw[a]->Integral(1, bsep - 1);
                std::cout << "  " << (ndofLo > 0 ? tot / ndofLo : 0.);
            }
            std::cout << "\n    and bins " << bsep << "-" << nptx << " (the rest, ndof=" << ndofHi
                      << "):";
            for (size_t a = 0; a < kApproaches.size(); ++a) {
                const double tot = chi2raw[a]->Integral(bsep, nptx);
                std::cout << "  " << (ndofHi > 0 ? tot / ndofHi : 0.);
            }
        }
        std::cout << "\n  pair-eta cells (= ndof) summed per pair-pT bin:";
        for (int bx = 1; bx <= nptx; ++bx) std::cout << " " << ncell[bx];
        std::cout << " (of " << neta << ")" << std::endl;
    }
    std::cout << "  shared ratio-pad range [" << rmin << ", " << rmax << "]" << std::endl;

    for (TFile* f : fin) f->Close();
    std::cout << "done." << std::endl;
}
