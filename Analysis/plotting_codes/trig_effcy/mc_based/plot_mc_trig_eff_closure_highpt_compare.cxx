// =================================================================================================
// plot_mc_trig_eff_closure_highpt_compare.cxx
//
// THE dR-CORRECTION PROCEDURE vs THE SINGLE-VALUE PAIR EFFICIENCY, ABOVE 52.2 GeV
// (docs/tracking/mc_trigeff_single_value_pair_eff.md PP-4 and R8; inputs from FillMCTrigEffClosure)
//
// TWO COMPARISONS, each into its OWN subdirectory of
// <plot base>/closure/single_value_highpt_comparison/ -- one PNG each per working point. Both keep
// the SAME two dR-procedure series as the reference against which the single-value procedure is
// judged, so the two can be read side by side:
//
//   mass_window_compr/   does the MASS WINDOW matter?
//       no trigger requirement                                                     black, open
//       dR, top two pair-pT cells merged        (`nocorr_ptmerge`)                 kBlue
//       dR, that plus the |eta^pair| fold       (`nocorr_etamerge_ptmerge`)        kMagenta
//       single value, 1.08 < m < 2.9 GeV (signal window)                           kRed
//       single value,    1 < m < 4   GeV (template-fit window)                     kGreen+2
//
//   pt_merge_compr/      does MERGING the top two pair-pT cells cost anything?
//       the same no-trigger and two dR series
//       single value, signal window, the canonical 8 pair-pT cells                 kRed
//       single value, signal window, top two combined into [74.24, 150) GeV        kGreen+2
//
// The dR series weight a firing pair by 1/[eps_MC(1) eps_MC(2) eps_dR(dR ; cell)]; the single-value
// ones replace that product by ONE measured number per (pair pT, |eta^pair|, sign) cell.
//
// THE RAW MEASURED NUMBER IS WHAT IS APPLIED (user, 2026-09-08). This is an MC-only test, so no
// data/MC correction belongs in it: the plan for the cross-section is to correct eps^pair by the
// product of the two single-muon data/MC scale factors, but such a factor cancels identically in an
// MC closure and would only obscure what the test measures. An earlier version of this macro also
// drew a "calibrated" series, w = 1/[eps(1) eps(2) K] with K = eps^pair / <eps_1 eps_2>; it is
// withdrawn, because writing the weight that way RE-INTRODUCES the very factorization the
// single-value procedure exists to avoid.
//
// WHY ONLY ABOVE 52.2 GeV. The single-value efficiency is delivered for the coarse pair-pT cells
// above 52.23 GeV, because that is where the dR-shape fit runs out of pairs. The lowest delivered
// cell is the CONTROL region, where the dR procedure still works and the two must agree.
//
// WHY ONLY THE SIGNAL SAMPLE VERSION. The mass window is part of the single-value efficiency's
// definition, so the all-opposite-sign sample is a different mass mixture and the number does not
// apply to it. FillMCTrigEffClosure books the single-value numerators for the `signal` version
// only, for the same reason.
//
// COVERAGE, NOT NON-CLOSURE. A presentation bin that straddles a cell edge, or that lies in a cell
// the delivery gate refuses, is only partly reached by the single-value correction. Such bins are
// identified from the closure file's own coverage denominator (`den_paireff_<window><mode suffix>`
// against `den`) and the single-value point is OMITTED there -- drawing it would show a coverage
// artefact as if the procedure did not close. The two cell modes have DIFFERENT coverage (that is
// the point of the merge), so the coverage key carries the mode suffix and the omitted-bin count
// is reported per series.
//
// SAME BINNING AS THE CROSS-SECTION, ALWAYS: ParamsSet::pT_bins_150 x the 9
// CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap panels (doc PP-1/D8), enforced here
// rather than assumed -- the two closure files' denominators are compared bin by bin and the panel
// edges are checked against the live ranges.
//
// Usage (from Analysis/plotting_codes/trig_effcy/mc_based/):
//   root -l -b -q 'plot_mc_trig_eff_closure_highpt_compare.cxx+("pp_full", true)'
//   root -l -b -q 'plot_mc_trig_eff_closure_highpt_compare.cxx+("pp_full", false)'   // Medium
// =================================================================================================

#include <algorithm>
#include <functional>
#include <map>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <TCanvas.h>
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

#include "dr_correction_sample_cfg.h"
#include "dr_correction_apply.h"
#include "dr_correction_ratio.h"
#include "../../../Utilities/CommonLogYRange.h"
#include "../../../Utilities/MCTrigEffPairPtBinning.h"
#include "../../../Utilities/PairTrigEffEvaluator.h"
#include "../../../RDFBasedHistFilling/CommonEffcyConfig.h"
#include "../../../Utilities/PairEtaPanelBins.h"

namespace {

// One drawn series. `num_key` is the histogram-name suffix FillMCTrigEffClosure wrote it under;
// `file_index` says which of the two closure files it comes from (0 = the pT-merged dR approach,
// which also carries the single-value numerators, 1 = the eta-folded one).
struct Series {
    std::string num_key;
    std::string legend;
    std::string cov;           // non-empty => a single-value series, masked by this coverage token
    std::string mode;          // cell mode of a single-value series ("" for the dR ones)
    int         file_index;
    Color_t     colour;
    Style_t     marker;
};

// TWO COMPARISONS, each answering a different question and each into its own subdirectory.
//   kMassWindow  does the mass window matter?  -> the two windows on the canonical cells
//   kPtMerge     does merging the top two pair-pT cells help? -> the signal window, both cell modes
// Both keep the SAME two dR-procedure series as the reference against which the single-value
// procedure is judged, so the two figures can be read side by side.
enum class Comparison { kMassWindow, kPtMerge };

const char* kCascadeKey = "cascade";
const std::vector<std::string> kModes = {"nocorr_ptmerge", "nocorr_etamerge_ptmerge"};

// The single-value numerators live in BOTH files (they do not depend on the dR cell grouping);
// they are read from file 0 and cross-checked against file 1, so a stale file cannot slip in.
std::vector<Series> BuildSeries(Comparison what)
{
    const std::string p = "paireff_";
    std::vector<Series> v = {
        {std::string("num_epsmc_") + kCascadeKey, "#varepsilon_{#DeltaR}(#DeltaR):  top two "
         "p_{T}^{pair} cells combined", "", "", 0, kBlue, 22},
        {std::string("num_epsmc_") + kCascadeKey, "#varepsilon_{#DeltaR}(#DeltaR):  that plus the "
         "|#eta^{pair}| fold", "", "", 1, kMagenta, 33},
    };
    if (what == Comparison::kMassWindow) {
        v.push_back({"num_" + p + "sig",
                     "single value:  1.08 < m_{#mu#mu} < 2.9 GeV (signal window)",
                     "sig", "nomerge", 0, kRed, 20});
        v.push_back({"num_" + p + "wide",
                     "single value:  1 < m_{#mu#mu} < 4 GeV (template-fit window)",
                     "wide", "nomerge", 0, kGreen + 2, 21});
    } else {
        v.push_back({"num_" + p + "sig",
                     "single value, signal window:  8 p_{T}^{pair} cells",
                     "sig", "nomerge", 0, kRed, 20});
        v.push_back({"num_" + p + "sig_ptmerge",
                     "single value, signal window:  top two p_{T}^{pair} cells combined",
                     "sig_ptmerge", "ptmerge", 0, kGreen + 2, 21});
    }
    return v;
}

struct ComparisonCfg { Comparison what; std::string subdir; std::string headline_extra; };
const std::vector<ComparisonCfg> kComparisons = {
    {Comparison::kMassWindow, "mass_window_compr", ""},
    {Comparison::kPtMerge,    "pt_merge_compr",    ""},
};

const Color_t kUncorrColour = kBlack;
const Style_t kUncorrMarker = 24;

// The sample version: the single-value numerators exist for this one only (header block).
const char* kVersion = "signal";

std::string EtaLabel(const std::pair<float, float>& r)
{
    return Form("%.1f < #eta^{pair} < %.1f", r.first, r.second);
}

template <typename T>
T* Get(TFile* f, const std::string& n)
{
    T* o = dynamic_cast<T*>(f->Get(n.c_str()));
    if (!o) throw std::runtime_error("plot_mc_trig_eff_closure_highpt_compare: missing '" + n
                                     + "' in " + f->GetName()
                                     + " -- rerun FillMCTrigEffClosure for that mode");
    return o;
}

TH1D* Row(TH2D* h, int iz, const std::string& nm)
{
    TH1D* p = h->ProjectionX(nm.c_str(), iz, iz, "e");
    p->SetDirectory(nullptr);
    return p;
}

}  // namespace

// =================================================================================================
void plot_mc_trig_eff_closure_highpt_compare(const std::string& sample = "pp_full",
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

    // The zoom edge is the CELL boundary the delivered range starts at, read from the canonical
    // coarse axis -- not a round 50 typed here (.claude/CLAUDE.md Binnings).
    const std::vector<double> coarse = PairTrigEff::PairPtEdges();
    const double pt_cell_lo = coarse[PairTrigEff::FirstDeliveredPtBin() - 1];
    const double pt_zoom_hi = coarse.back();
    // The DRAWN lower edge is a presentation-bin edge, not the cell edge: TAxis::SetRangeUser snaps
    // to whole bins, so a frame asked to start at the cell edge actually starts at the low edge of
    // the bin containing it. That edge -- read back from the axis, never assumed -- is what the
    // header states, so the label and the axis cannot disagree (.claude/CLAUDE.md Binnings item 5).
    double pt_zoom_lo = pt_cell_lo;   // resolved against the histogram axis once it is open

    // One SUBDIRECTORY per comparison, under the shared parent: the two figures answer different
    // questions on the same axes, so they belong beside each other rather than in one flat folder.
    const std::string base_outdir = cfg.out_base + "closure/single_value_highpt_comparison/";

    // ---------------- open the two dR-approach closure files ---------------------------------
    std::vector<TFile*> fin(kModes.size(), nullptr);
    for (size_t a = 0; a < kModes.size(); ++a) {
        const std::string p = cfg.mc_dir + "mc_trig_eff_closure_" + cfg.mc_label + wp_suf
                            + MCTrigEffPairPt::FileSuffix()
                            + DrCorrPlateauModeTag(kModes[a]) + ".root";
        fin[a] = TFile::Open(p.c_str(), "READ");
        if (!fin[a] || fin[a]->IsZombie())
            throw std::runtime_error("plot_mc_trig_eff_closure_highpt_compare: cannot open " + p
                                     + " -- run FillMCTrigEffClosure with plateau mode '"
                                     + kModes[a] + "' first");
        std::cout << "  " << kModes[a] << " <- " << p << std::endl;
    }

    const std::string pre = std::string("h_closure_") + kVersion + "_";

    // ---------------- the shared denominator, checked rather than assumed --------------------
    TH2D* h_den = Get<TH2D>(fin[0], pre + "den");
    for (size_t a = 1; a < kModes.size(); ++a) {
        TH2D* d = Get<TH2D>(fin[a], pre + "den");
        if (d->GetNbinsX() != h_den->GetNbinsX() || d->GetNbinsY() != h_den->GetNbinsY())
            throw std::runtime_error("plot_mc_trig_eff_closure_highpt_compare: the two closure "
                "files were filled on different histogram binnings -- the overlay would compare "
                "different cells");
        for (int bx = 1; bx <= d->GetNbinsX(); ++bx)
            for (int by = 1; by <= d->GetNbinsY(); ++by) {
                const double u = h_den->GetBinContent(bx, by), w = d->GetBinContent(bx, by);
                // Relative, never floored at 1: every bin is a weighted yield well below 1 nb, so
                // a max(1.0, |u|) floor would make this an absolute 1e-6 nb test that two quite
                // different files would pass (see the coverage mask below).
                if (std::fabs(u - w) > 1e-6 * std::fabs(u))
                    throw std::runtime_error(Form("plot_mc_trig_eff_closure_highpt_compare: the "
                        "no-trigger denominators of '%s' and '%s' differ in bin (%d,%d): %g vs %g "
                        "-- the two runs did not see the same sample",
                        kModes[0].c_str(), kModes[a].c_str(), bx, by, u, w));
            }
    }
    if (h_den->GetNbinsY() != neta)
        throw std::runtime_error("plot_mc_trig_eff_closure_highpt_compare: the histograms carry "
            + std::to_string(h_den->GetNbinsY()) + " pair-eta bins but CommonEffcyConfig has "
            + std::to_string(neta) + " -- stale input?");
    PairEtaPanels::CheckPanelsMatchFiducialCut();
    for (int iz = 1; iz <= neta; ++iz)
        PairEtaPanels::CheckAxisAligned(h_den->GetYaxis(), eta_ranges[iz - 1], iz, iz,
                                        "plot_mc_trig_eff_closure_highpt_compare");
    pt_zoom_lo = h_den->GetXaxis()->GetBinLowEdge(
        h_den->GetXaxis()->FindFixBin(pt_cell_lo * 1.0001));

    // The single-value numerators are written into every mode's file; they must be IDENTICAL, or
    // one of the two runs used a different pair-efficiency file.
    for (const auto& C : kComparisons)
        for (const auto& S : BuildSeries(C.what)) {
            if (S.cov.empty()) continue;
            TH2D* a0 = Get<TH2D>(fin[0], pre + S.num_key);
            TH2D* a1 = Get<TH2D>(fin[1], pre + S.num_key);
            for (int bx = 1; bx <= a0->GetNbinsX(); ++bx)
                for (int by = 1; by <= a0->GetNbinsY(); ++by) {
                    const double u = a0->GetBinContent(bx, by), w = a1->GetBinContent(bx, by);
                    if (std::fabs(u - w) > 1e-6 * std::fabs(u))
                        throw std::runtime_error(Form("plot_mc_trig_eff_closure_highpt_compare: "
                            "the single-value numerator '%s' differs between the two closure "
                            "files in bin (%d,%d): %g vs %g -- one of them was filled against a "
                            "different pair_trig_eff file", S.num_key.c_str(), bx, by, u, w));
                }
        }

    // ---------------- build and draw, once per applied form ----------------------------------
    for (const auto& C : kComparisons) {
      const std::string outdir = base_outdir + C.subdir + "/";
      gSystem->mkdir(outdir.c_str(), kTRUE);
      {
        const std::vector<Series> series = BuildSeries(C.what);

        // The cell MODES actually drawn. Computed here because the header needs one line per mode,
        // and the header's height sets where the panels start: with two modes the block is six
        // lines and a fixed 0.21 strip put the last one on top of the first panel row.
        std::vector<std::string> modes_drawn;
        for (const auto& S : series)
            if (!S.mode.empty()
                && std::find(modes_drawn.begin(), modes_drawn.end(), S.mode) == modes_drawn.end())
                modes_drawn.push_back(S.mode);

        // ---- per pair-eta panel: the spectra, the ratios, and the coverage mask -------------
        std::vector<TH1D*> spec_unc(neta, nullptr);
        std::vector<std::vector<TH1D*>> spec_cor(neta, std::vector<TH1D*>(series.size(), nullptr));
        std::vector<std::vector<TGraphErrors*>> gratio(neta,
                                                std::vector<TGraphErrors*>(series.size(), nullptr));
        std::vector<TH1D*> ratio_frame(neta, nullptr);
        std::vector<double> ratio_pts;
        std::vector<TH1*> for_range;
        const std::string tag = C.subdir;
        // The masked (panel, bin) pairs, as a UNION over the single-value windows: the header
        // states how many BINS lost their single-value points, so counting once per series would
        // report twice the number, and counting only the first window would undercount if the two
        // windows ever refused different cells.
        std::set<std::pair<int, int>> masked_bins;                 // the union, over all series
        std::map<std::string, std::set<std::pair<int, int>>> masked_by_series;

        for (int iz = 1; iz <= neta; ++iz) {
            TH1D* den = Row(h_den, iz, Form("hp_%s_den_eta%d", tag.c_str(), iz));

            for (size_t a = 0; a < series.size(); ++a) {
                const Series& S = series[a];
                TFile* f = fin[S.file_index];
                TH1D* num = Row(Get<TH2D>(f, pre + S.num_key),
                                iz, Form("hp_%s_num%zu_eta%d", tag.c_str(), a, iz));
                TH1D* A   = Row(Get<TH2D>(f, pre + "numA" + S.num_key.substr(3)),
                                iz, Form("hp_%s_A%zu_eta%d", tag.c_str(), a, iz));
                TH1D* B   = Row(Get<TH2D>(f, pre + "numB" + S.num_key.substr(3)),
                                iz, Form("hp_%s_B%zu_eta%d", tag.c_str(), a, iz));

                // The ratio, formed exactly as in the per-approach figures: divide first, then the
                // single conditional-error implementation (A and B are sums of squared weights and
                // do not scale like a density), so a point here equals the same point there.
                TH1D* r = static_cast<TH1D*>(num->Clone(Form("hp_%s_ratio%zu_eta%d",
                                                             tag.c_str(), a, iz)));
                r->SetDirectory(nullptr);
                r->Divide(den);
                SetConditionalRatioErrors(r, den, A, B);
                delete A; delete B;

                // COVERAGE MASK for the single-value series. A presentation bin counts only if the
                // window's own coverage denominator equals the full one there -- i.e. every pair
                // in the bin has a measured cell. This distinguishes "not covered" from a genuine
                // C = 0 (a bin where pairs exist and none fired), which a content-based test could
                // not.
                std::vector<bool> ok(r->GetNbinsX() + 2, true);
                if (!S.cov.empty()) {
                    TH1D* dcov = Row(Get<TH2D>(fin[0], pre + "den_paireff_" + S.cov),
                                     iz, Form("hp_%s_cov%zu_eta%d", tag.c_str(), a, iz));
                    for (int b = 1; b <= r->GetNbinsX(); ++b) {
                        const double dfull = den->GetBinContent(b), dc = dcov->GetBinContent(b);
                        // RELATIVE tolerance, and it must stay relative. Bin contents here are
                        // weighted yields in nb, all of them far below 1, so a `max(1.0, dfull)`
                        // floor turns this into an ABSOLUTE 1e-6 nb test that every bin passes --
                        // including a bin with ZERO coverage. That is what it did until
                        // 2026-09-08: the forward [105.53,150) x |eta| 2-2.2 cell is refused by the
                        // reader (4 raw pairs), so bin [101.47, 123.37) of the 2.0 < eta^pair < 2.2
                        // panel has den = 6.2e-07 and den_paireff = 0, yet was drawn -- at C = 0,
                        // exactly the coverage artefact PP-4 exists to keep off the figure.
                        const bool covered = dfull > 0.
                                          && std::fabs(dc - dfull) <= 1e-6 * dfull;
                        ok[b] = covered;
                        if (!covered && dfull > 0.
                            && den->GetXaxis()->GetBinUpEdge(b) > pt_zoom_lo)
                        { masked_bins.insert({iz, b}); masked_by_series[S.cov].insert({iz, b}); }
                    }
                    delete dcov;
                }

                // A SMALL x OFFSET PER SERIES. All four corrected series are defined at the same
                // bin centre, so drawn there they paint over one another and the later ones hide
                // the earlier -- measured: blue 0.8707 and green 0.8700 in one bin, the blue marker
                // completely covered. The offset is GEOMETRIC (the axis is log), a fixed fraction
                // of the bin's own width, so it is the same visual shift in every bin and cannot be
                // mistaken for a different p_T: the points still belong to the bin they sit in.
                const double shift = ((double)a - 0.5 * (series.size() - 1)) * 0.11;
                auto* g = new TGraphErrors();
                g->SetName(Form("hp_%s_g%zu_eta%d", tag.c_str(), a, iz));
                for (int b = 1; b <= r->GetNbinsX(); ++b) {
                    if (den->GetBinContent(b) <= 0. || !ok[b]) continue;
                    if (den->GetXaxis()->GetBinLowEdge(b) < pt_zoom_lo
                        && den->GetXaxis()->GetBinUpEdge(b) <= pt_zoom_lo) continue;
                    const double wlo = den->GetXaxis()->GetBinLowEdge(b);
                    const double whi = den->GetXaxis()->GetBinUpEdge(b);
                    const int n = g->GetN();
                    g->SetPoint(n, r->GetBinCenter(b) * std::pow(whi / wlo, shift),
                                r->GetBinContent(b));
                    g->SetPointError(n, 0.0, r->GetBinError(b));
                    ratio_pts.push_back(r->GetBinContent(b));
                }
                gratio[iz - 1][a] = g;
                if (a == 0) {
                    ratio_frame[iz - 1] = static_cast<TH1D*>(r->Clone(
                        Form("hp_%s_frame_eta%d", tag.c_str(), iz)));
                    ratio_frame[iz - 1]->SetDirectory(nullptr);
                    ratio_frame[iz - 1]->Reset();
                }
                delete r;

                // The spectrum pad is log-y, so a masked (zeroed) bin simply does not appear --
                // no artificial point at 0.
                for (int b = 1; b <= num->GetNbinsX(); ++b)
                    if (!ok[b]) { num->SetBinContent(b, 0.); num->SetBinError(b, 0.); }
                num->Scale(1.0, "width");
                spec_cor[iz - 1][a] = num;
                for_range.push_back(num);
            }
            den->Scale(1.0, "width");
            spec_unc[iz - 1] = den;
            for_range.push_back(den);
        }
        // The common log-y scale must be derived from the bins the figure actually DRAWS.
        // ApplyCommonLogYRange scans every bin of every histogram it is given, and SetRangeUser is
        // applied later, so scanning the FULL pair-pT spectrum put the frame ceiling ~2.8 decades
        // above the highest drawn point -- 37 % of every upper pad empty. Bins below the zoom are
        // zeroed here (they can never be displayed: the frame starts at the first drawn bin) so the
        // shared helper sees only the drawn range, and the one-scale-for-all-panels property is
        // kept.
        {
            const int b_first = spec_unc[0]->GetXaxis()->FindFixBin(pt_zoom_lo * 1.0001);
            for (TH1* h : for_range)
                for (int b = 1; b < b_first; ++b) { h->SetBinContent(b, 0.); h->SetBinError(b, 0.); }
        }
        ApplyCommonLogYRange(for_range);

        double rmin = 1.0, rmax = 1.0;
        if (!ratio_pts.empty()) {
            rmin = *std::min_element(ratio_pts.begin(), ratio_pts.end());
            rmax = *std::max_element(ratio_pts.begin(), ratio_pts.end());
        }
        rmin = std::min(rmin, 1.0); rmax = std::max(rmax, 1.0);
        if (rmax - rmin < 0.10) { const double m = 0.5 * (rmin + rmax); rmin = m - 0.05; rmax = m + 0.05; }
        const double rpad = 0.08 * (rmax - rmin);
        rmin -= rpad; rmax += rpad;

        // ---- draw ---------------------------------------------------------------------------
        // One extra 0.02 strip per additional cell-mode line, so the header can never reach down
        // into the panels.
        const double head = 0.21 + 0.02 * (double)(modes_drawn.size() > 1 ? modes_drawn.size() - 1
                                                                          : 0);
        TCanvas c(("c_hp_" + tag).c_str(), "", 1800, 1800);

        for (int iz = 1; iz <= neta; ++iz) {
            c.cd();
            const int col = (iz - 1) % ncol, row = (iz - 1) / ncol;
            const double x1 = col / (double)ncol, x2 = (col + 1) / (double)ncol;
            const double ytop = 1.0 - head - row * (1.0 - head) / nrow;
            const double ybot = 1.0 - head - (row + 1) * (1.0 - head) / nrow;
            const double ysplit = ybot + 0.36 * (ytop - ybot);

            auto* pu = new TPad(Form("hp_pu_%s_%d", tag.c_str(), iz), "", x1, ysplit, x2, ytop);
            auto* pl = new TPad(Form("hp_pl_%s_%d", tag.c_str(), iz), "", x1, ybot,   x2, ysplit);
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
            TH1D* du = spec_unc[iz - 1];
            du->SetTitle("");
            du->GetXaxis()->SetRangeUser(pt_zoom_lo, pt_zoom_hi);
            du->GetYaxis()->SetTitle("d#sigma/dp_{T}^{pair} [nb/GeV]");
            du->GetYaxis()->SetTitleSize(0.070); du->GetYaxis()->SetLabelSize(0.060);
            du->GetYaxis()->SetTitleOffset(1.55);
            du->GetXaxis()->SetLabelSize(0.0);
            du->SetMarkerColor(kUncorrColour); du->SetLineColor(kUncorrColour);
            du->SetMarkerStyle(kUncorrMarker); du->SetMarkerSize(1.1);
            du->Draw("PE");
            for (size_t a = 0; a < series.size(); ++a) {
                TH1D* h = spec_cor[iz - 1][a];
                h->GetXaxis()->SetRangeUser(pt_zoom_lo, pt_zoom_hi);
                h->SetMarkerColor(series[a].colour); h->SetLineColor(series[a].colour);
                h->SetMarkerStyle(series[a].marker); h->SetMarkerSize(1.1);
                h->Draw("PE SAME");
            }
            auto* tl = new TLatex();
            tl->SetNDC(); tl->SetTextFont(42); tl->SetTextSize(0.068);
            tl->DrawLatex(0.28, 0.09, EtaLabel(eta_ranges[iz - 1]).c_str());

            pl->cd();
            TH1D* fr = ratio_frame[iz - 1];
            fr->SetTitle("");
            fr->GetXaxis()->SetRangeUser(pt_zoom_lo, pt_zoom_hi);
            fr->GetYaxis()->SetRangeUser(rmin, rmax);
            fr->GetYaxis()->SetTitle("corrected / no trigger");
            fr->GetYaxis()->SetNdivisions(505);
            fr->GetYaxis()->SetTitleSize(0.100); fr->GetYaxis()->SetLabelSize(0.095);
            fr->GetYaxis()->SetTitleOffset(1.08);
            fr->GetXaxis()->SetTitle("p_{T}^{pair} [GeV]");
            fr->GetXaxis()->SetTitleSize(0.125); fr->GetXaxis()->SetLabelSize(0.105);
            fr->GetXaxis()->SetTitleOffset(1.00);
            // The drawn range is the top part of the pair-pT axis, inside which ROOT's default
            // log labelling prints only "10^{2}" -- so the reader cannot find the coarse cell
            // edges this figure is about. The edges are not retyped here (Binnings rule 1) and
            // pt_zoom_lo is derived from the axis, so both follow ParamsSet automatically.
            fr->GetXaxis()->SetMoreLogLabels();
            fr->GetXaxis()->SetNoExponent();
            fr->Draw("AXIS");
            auto* ln = new TLine(pt_zoom_lo, 1.0, pt_zoom_hi, 1.0);
            ln->SetLineStyle(2); ln->SetLineColor(kGray + 2);
            ln->Draw();
            for (size_t a = 0; a < series.size(); ++a) {
                TGraphErrors* g = gratio[iz - 1][a];
                g->SetMarkerColor(series[a].colour); g->SetLineColor(series[a].colour);
                g->SetMarkerStyle(series[a].marker); g->SetMarkerSize(1.1);
                if (g->GetN() > 0) g->Draw("P SAME");
            }
        }

        // ---- header strip --------------------------------------------------------------------
        c.cd();
        auto* hd = new TLatex();
        hd->SetNDC(); hd->SetTextFont(42);
        hd->SetTextSize(0.0175);
        hd->DrawLatex(0.030, 0.984,
                      (cfg.sample_text + ",  " + wp_text
                       + ",  opposite-sign pairs in the single-b signal region,  p_{T}^{pair} > "
                       + Form("%.1f", pt_zoom_lo) + " GeV").c_str());
        hd->SetTextSize(0.0145);
        hd->DrawLatex(0.030, 0.884,
                      "#DeltaR procedure:   w = 1 / [#varepsilon(p_{T,1},q#eta_{1}) "
                      "#varepsilon(p_{T,2},q#eta_{2}) #varepsilon_{#DeltaR}(#DeltaR)]");
        hd->DrawLatex(0.030, 0.864,
                      "single value:   w = 1 / #varepsilon_{2#mu4}^{pair},   "
                      "#varepsilon_{2#mu4}^{pair} = #Sigma_{2#mu4} w / #Sigma_{all} w "
                      "per (p_{T}^{pair}, |#eta^{pair}|, sign) cell");
        // EVERY edge in these lines is READ, never typed: the pair-pT cells from the canonical
        // coarse vector, the |eta^pair| groups from the live fold. The outer pair-eta edge already
        // moved 2.4 -> 2.2 on 2026-09-07; a typed label would now silently disagree with the axis
        // the points came from (.claude/CLAUDE.md Binnings items 1 and 5).
        // The cell list gets a LINE OF ITS OWN: appended to the epsilon line it ran off the right
        // edge of the canvas and was clipped mid-number, which is worse than no label at all.
        hd->DrawLatex(0.030, 0.844,
                      "#varepsilon = single-muon mu4 efficiency from the MC turn-on");
        // The cell list of every cell MODE actually drawn, each edge read from that mode's own
        // axis. In the pair-pT-merge comparison the two single-value series live on DIFFERENT
        // cells, and a figure that showed only one of the two lists would mislabel the other.
        std::string etacols;
        {
            const std::vector<double>& ge = PairTrigEff::AbsEtaGroups().edges;
            for (size_t g = 0; g + 1 < ge.size(); ++g)
                etacols += Form("%s %g#minus%g", g ? "," : "", ge[g], ge[g + 1]);
        }
        double ytext = 0.824;
        for (const std::string& m : modes_drawn) {
            const std::vector<double> me = PairTrigEff::PairPtEdges(m);
            const int lo = PairTrigEff::FirstDeliveredPtBin(m);
            std::string cells = "single-value cells";
            if (modes_drawn.size() > 1)
                cells += (m == "ptmerge" ? " (top two combined)" : " (8 p_{T}^{pair} cells)");
            cells += ": ";
            for (size_t i = lo - 1; i + 1 < me.size(); ++i)
                cells += Form("%s[%.2f, %.2f)", i == (size_t)lo - 1 ? "" : ", ", me[i], me[i + 1]);
            cells += Form(" GeV  #times  |#eta^{pair}|%s", etacols.c_str());
            hd->DrawLatex(0.030, ytext, cells.c_str());
            ytext -= 0.020;
        }
        // THE OMITTED-BIN COUNT IS NOT DRAWN (user, 2026-09-08). Which bins are absent is not
        // something the reader has to act on: "not fully inside a measured cell" is an artefact of
        // the presentation binning being finer than the correction cells, and the bins dropped for
        // want of statistics are already visible, cell by cell and with their raw pair counts, in
        // the statistics tables beside this figure. It is reported to the LOG instead, where it
        // still serves as a check that the mask did what it should.
        for (const auto& S : series) {
            if (S.cov.empty()) continue;
            std::cout << "  [" << tag << "] " << S.cov << ": "
                      << masked_by_series[S.cov].size()
                      << " presentation bins omitted (not fully inside a measured cell, or in a "
                         "cell with fewer than " << PairTrigEff::MinCellPairs() << " pairs)"
                      << std::endl;
        }

        auto* leg = new TLegend(0.030, 0.900, 0.990, 0.976);
        leg->SetNColumns(2);
        leg->SetBorderSize(0); leg->SetFillStyle(0);
        leg->SetTextFont(42); leg->SetTextSize(0.0145);
        leg->AddEntry(spec_unc[0], "no trigger requirement", "PE");
        for (size_t a = 0; a < series.size(); ++a)
            leg->AddEntry(spec_cor[0][a], series[a].legend.c_str(), "PE");
        leg->Draw();

        const std::string png = outdir + "closure_highpt_single_value_" + C.subdir + ".png";
        c.SaveAs(png.c_str());
        std::cout << "  wrote " << png << std::endl;

        // ---- the numbers behind the figure ---------------------------------------------------
        // ⚠ THE DENOMINATOR IS PER SERIES HERE, and deliberately so. The figure drops the
        // partially covered presentation bin from the single-value series (so the five lines keep
        // ONE common no-trigger reference); this table instead divides each single-value series by
        // its OWN coverage denominator `den_paireff_<window>`, which is the exactly matched
        // no-trigger yield of the pairs that series actually corrects. Dividing it by the FULL
        // denominator would mix a coverage fraction into the closure -- integrated above the cell
        // edge that mixture is large, because the straddling bin holds most of the yield -- and
        // would read as a 30 % non-closure that is nothing of the kind.
        {
            const int b_lo = h_den->GetXaxis()->FindFixBin(pt_zoom_lo * 1.0001);
            const int b_hi = h_den->GetXaxis()->GetNbins();
            std::cout << "\n===== closure above " << pt_zoom_lo << " GeV, " << C.subdir << ", "
                      << wp_text << " =====\n"
                      << "  (each series integrated over its own covered region of the panel; the "
                         "dR series cover every bin,\n   the single-value series only pair p_T >= "
                      << pt_cell_lo << " GeV, so they carry their own matched denominator)\n  "
                      << std::setw(24) << std::left << "eta^pair";
            for (const auto& S : series)
                std::cout << std::setw(22)
                          << (S.cov.empty()
                                  ? (S.file_index == 0 ? "dR: pT-merged" : "dR: +|eta| fold")
                                  : ("single value " + S.cov));
            std::cout << std::endl;

            auto row = [&](int iz_lo, int iz_hi, const char* label) {
                std::cout << "  " << std::setw(24) << std::left << label;
                for (const auto& S : series) {
                    TH2D* n = Get<TH2D>(fin[S.file_index], pre + S.num_key);
                    TH2D* d = S.cov.empty()
                                  ? h_den
                                  : Get<TH2D>(fin[0], pre + "den_paireff_" + S.cov);
                    const double D = d->Integral(b_lo, b_hi, iz_lo, iz_hi);
                    std::cout << std::setw(22)
                              << (D > 0 ? Form("%.4f", n->Integral(b_lo, b_hi, iz_lo, iz_hi) / D)
                                        : "-");
                }
                std::cout << std::endl;
            };
            for (int iz = 1; iz <= neta; ++iz)
                row(iz, iz, Form("[%.1f, %.1f)", eta_ranges[iz - 1].first,
                                 eta_ranges[iz - 1].second));
            row(1, neta, "ALL PANELS");
            // The coverage fraction itself, so the reader can see how much of the region above the
            // cell edge each single-value series reaches at all. Per DRAWN series, because merging
            // the top pair-pT cells changes the answer -- that is the point of the merge.
            for (const auto& S : series) {
                if (S.cov.empty()) continue;
                TH2D* d = Get<TH2D>(fin[0], pre + "den_paireff_" + S.cov);
                const double Dc = d->Integral(b_lo, b_hi, 1, neta);
                const double Df = h_den->Integral(b_lo, b_hi, 1, neta);
                std::cout << "  coverage (" << S.cov << "): the single-value cells reach "
                          << (Df > 0 ? 100.0 * Dc / Df : 0.0)
                          << "% of the no-trigger yield above " << pt_zoom_lo
                          << " GeV (the rest is the bin straddling the cell edge, plus any cell "
                             "the delivery gate refuses)" << std::endl;
            }
        }
      }
    }

    for (auto* f : fin) f->Close();
    std::cout << "done." << std::endl;
}
