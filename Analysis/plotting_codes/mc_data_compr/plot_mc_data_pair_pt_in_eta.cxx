// =================================================================================================
// plot_mc_data_pair_pt_in_eta.cxx  --  SIGNAL family
//
// pp24 data vs Pythia truth: the single-b dimuon pair-pT cross-section over 8-150 GeV, both
// pair-eta-integrated and split into the 9 canonical pair-eta bins. The MC counterpart of
// plots/single_b_analysis/pp24/pp24_crossx_pair_pt_in_eta_subplots.png (user, 2026-08-18).
// Both views project the SAME 2D histogram, so the integrated panel and the nine panels cannot
// describe different cells.
//
// THE MC PARTNER IS THE pp24-CONDITION FULLSIM FULL SAMPLE. It replaced the Pythia TRUTH full
// sample on 2026-08-25 for a physics reason, not a convenience one: the truth sample's weight
// sets the beam ratio to the Pb isospin mixture 4:6:6:9 unconditionally, so its absolute sigma is
// a Pb isospin-averaged NN cross-section -- not a pp cross-section, and therefore the wrong object
// to place beside pp24 data (16.6 % in the signal region). The fullsim pp24 sample is pp beam only
// with isospin weight 1, is AMI-weighted, and already implements the data signal region bit for
// bit on truth quantities.
//
//   data : dimuon_data/pp_2024/histograms_real_pairs_pp_2024_2mu4_nominal.root
//          h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts
//          = (1/L) * sum 1/(eps_trig * eps_reco), L = 400.412 pb^-1, single-b signal region
//            (m in (1.08,2.9), pair pT > ParamsSet::signal_pair_pt_min, BOTH muons passing the
//            fiducial gap cut and |eta^pair| < ParamsSet::pair_eta_fiducial_max), Tight WP.
//   Pythia: pythia_fullsim_full_sample/histograms_pythia_fullsim_pp24_no_data_resonance_cuts_full.root
//          h_truth_pair_eta_crossx_vs_truth_pair_pt_log_150_single_b_pass_signal_truth
//          = truth from_same_b OS pairs in the truth signal region, AMI-weighted.
//   POWHEG: powheg_full_sample/histograms_powheg_fullsim_pp17.root
//          h_truth_pair_eta_crossx_vs_truth_pair_pt_log_150_single_b_pass_signal_truth_gapcut
//          = the SAME truth selection (from_same_b + the data's mass window, pair-pT threshold
//            and fiducial gap windows on truth q*eta), per-sample cross-section weighted.
//
// POWHEG FULLSIM IS A **pp17-CONDITIONS** SAMPLE (sqrt(s_NN) = 5.02 TeV, Run-2 detector). No
// POWHEG FullSim pp24-conditions production exists, so this is the only POWHEG FullSim available
// and the difference in beam energy and detector conditions is REAL. The legend says so on every
// canvas; do not quote a POWHEG/data ratio from this figure without repeating it.
//
// NORMALIZATION: three samples, three different rules, none of them guessed.
//   data   : already d(sigma) in pb (the crossx histogram is 1/L-weighted at fill time).
//   Pythia : the AMI `crossSection` is in **nb** (PythiaAlgCoreT.c:549-557) -> exactly one
//            factor, nb -> pb = 1000.
//   POWHEG : `weight = EventWeights[0] * filter_effcy` is already a cross-section in **pb**
//            (PowhegAlgCoreT.c:209-212; read as nb it would be an 11 mb b-bbar cross-section,
//            55x the total sigma_bb at 5.02 TeV) -> factor **1**, and it must NOT be given
//            Pythia's x1000. The per-event division by the GENERATED statistics is done inside
//            the RDF stage, PER SAMPLE (RDFBasedHistFillingPowhegFullsim.cxx): the run reads the
//            bb and cc productions together and the class-wide `weight_norm` divides both by
//            N_bb + N_cc, which is an N-weighted average and not a sum -- a ~2x under-
//            normalization that is invisible in the reco-efficiency ratios the class was written
//            for but fatal here.
//
// KNOWN, DELIBERATE MISMATCHES -- state them wherever this plot is used, they are NOT bugs:
//   1. SIGNAL DEFINITION. MC is pure truth single-b; the data is raw OS with NO same-sign
//      subtraction and no template-fit background subtraction, so it still contains
//      gluon-splitting and combinatorial background.
//   2. The data is reco-level, corrected but NOT unfolded; the MC is truth-level.
// (The fiducial-region and isospin mismatches that the truth sample carried are GONE with it: the
// fullsim applies the same three gap windows to truth q*eta, and runs the pp beam alone.)
//
// BINNINGS: nothing here is retyped. The pair-pT and pair-eta axes come from the histograms
// themselves (both sides booked from `ParamsSet::pT_bins_150` and `ParamsSet::pair_eta_crossx_bins`
// by name), and the nine PANEL boundaries -- a different object -- are read from
// `CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap`.
// The fine pair-eta axis is 48 uniform bins on [-2.4, 2.4] (`ParamsSet::N_PAIR_ETA_CROSSX_BINS`),
// i.e. width exactly 0.1, so EVERY panel boundary -- the eight internal ones at +-0.5, 1.0, 1.5,
// 2.0 AND the outer ones at +-2.2 -- is a bin edge, and the projection below is exact by
// construction: no bin lands in two panels and the panel labels are literally the range drawn.
// (With the previous 44-bin axis none of the eight was an edge, eight bins were double-counted
// and the nine panels summed to +20.0 % of the true total.)
// SINCE 2026-09-07 the panels tile only |eta^pair| < ParamsSet::pair_eta_fiducial_max = 2.2, so
// the four outermost fine bins (|eta^pair| > 2.2) lie outside EVERY panel. They are empty here
// because all three inputs carry the pair-level gap cut. TWO different guards keep that honest:
// PairEtaPanels::Bins throws if a panel edge is not a bin edge of the histogram being projected
// (ALIGNMENT), and `assert_inside_fiducial` below throws if an input still has yield beyond +-2.2
// (EMPTINESS). The second is the one that matters here, because the pair-eta-INTEGRATED
// projection further down runs over the FULL axis while its label claims the fiducial region --
// the two views agree only through that emptiness.
//
// Output: /usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/mc_data_compr/signal/
// Usage:  root -l -b -q 'plot_mc_data_pair_pt_in_eta.cxx+()'           (Tight, nominal)
//         root -l -b -q 'plot_mc_data_pair_pt_in_eta.cxx+("medium")'   (Medium WP systematic)
// =================================================================================================

#include <algorithm>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>

#include "McDataComprColors.h"
#include "McDataComprConfig.h"
#include "McDataComprRatio.h"
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>

#include "../helper_functions.c"
#include "../../RDFBasedHistFilling/CommonEffcyConfig.h"
#include "../../Utilities/PairEtaPanelBins.h"

namespace {

// The NOMINAL data 2D: the UNSUFFIXED family, booked on ParamsSet::pT_bins_150 (16 log bins
// 9 -> 150 GeV) since 2026-09-09. It used to be the separate "pt_150" family; that family was
// demoted to the opt-in 9 -> 120 GeV alternative ("pt_120") because it covered only a handful of
// observables -- notably not the Pb+Pb R_AA global 3Ds -- so keeping it as nominal would have
// split the analysis across two pair-pT binnings.
const char* kDataHist = "h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts";

// The 2D naming convention is h_<y-variable>_vs_<x-variable><filter>
// (RDFBasedHistFillingBaseClass.cxx:561): x = pair pT, y = pair eta, same as the data histogram.
//
// The "_150" here is part of the var1D variable name `truth_pair_pt_log_150`
// (RDFBasedHistFilling/var1D_{pythia,powheg}_fullsim.json), which is bound to
// "binning": "pT_bins_150" -- i.e. the NOMINAL 9 -> 150 GeV axis, the same one the data histogram
// above now uses. There is no second truth pair-pT variable, so this side has no alternative view
// and nothing to select between. The token is a leftover of the old naming and says nothing about
// which axis is nominal; if the producers ever drop it, this name follows.
const char* kMcHist =
    "h_truth_pair_eta_crossx_vs_truth_pair_pt_log_150_single_b_pass_signal_truth";

// POWHEG FullSim pp17: the SAME 2D view on the SAME named axes, under the gap-cut truth signal
// filter added on 2026-09-03. (The class's pre-existing `_pass_signal_truth` was migrated off the
// retired one-sided q*eta < 2.2 on 2026-09-08 (D7), so the two now differ only by the extra
// truth_dr < 1.0 that `_single_b` carries.)
const char* kPowhegHist =
    "h_truth_pair_eta_crossx_vs_truth_pair_pt_log_150_single_b_pass_signal_truth_gapcut";

// AMI crossSection is in nb; the data is normalized by a luminosity in pb^-1.
const double kNbToPb = 1.0e3;
// POWHEG's weight is already a cross-section in pb. See the normalization block at the top.
const double kPowhegNorm = 1.0;

// Legend entries. POWHEG's names the sample's conditions because they are NOT the data's.
// The data entry names its beam energy too: a legend that says "5.02 TeV" on one curve and
// nothing on another leaves the contrast the POWHEG caveat exists for to be inferred.
const char* kDataLegend   = "pp data 2024 (5.36 TeV)";
const char* kPythiaLegend = "Pythia, single-b (pp24 cond.)";
// Two lines, not one: at nine-panel scale the one-line form ran past the pad edge and was
// clipped mid-word. The conditions line is what makes the curve honest, so it must always render.
const char* kPowhegLegend     = "POWHEG, single-b";
const char* kPowhegCondLegend = "  (pp17 cond., 5.02 TeV)";

const char* kOutDir = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/mc_data_compr/signal/";

TH2D* Fetch(const std::string& file, const char* hist)
{
    McDataComprConfig::AssertInputExists(file, hist);
    TFile* f = TFile::Open(file.c_str(), "READ");
    if (!f || f->IsZombie())
        throw std::runtime_error(std::string("cannot open ") + file);
    auto* h = dynamic_cast<TH2D*>(f->Get(hist));
    if (!h)
        throw std::runtime_error(std::string("no ") + hist + " in " + file);
    auto* c = static_cast<TH2D*>(h->Clone((std::string(hist) + "_clone").c_str()));
    c->SetDirectory(nullptr);
    f->Close();
    return c;
}

// The two axes MUST agree bin for bin: a panel cut at a different edge on each side would compare
// different eta ranges while every histogram still fills. Compared with a RELATIVE tolerance --
// one side is built as (nbins, min, max) and the other from a generated edge vector, and those
// agree only to double rounding.
void AssertSameAxes(const TH2D* a, const TH2D* b)
{
    auto cmp = [](const TAxis* x, const TAxis* y, const char* what) {
        if (x->GetNbins() != y->GetNbins())
            throw std::runtime_error(std::string("axis bin count differs on ") + what);
        for (int i = 1; i <= x->GetNbins() + 1; ++i) {
            const double ea = x->GetBinLowEdge(i);
            const double eb = y->GetBinLowEdge(i);
            if (std::fabs(ea - eb) > 1e-9 * std::max(1.0, std::fabs(ea)))
                throw std::runtime_error(std::string("axis edge ") + std::to_string(i)
                                         + " differs on " + what);
        }
    };
    cmp(a->GetXaxis(), b->GetXaxis(), "pair pT");
    cmp(a->GetYaxis(), b->GetYaxis(), "pair eta");
}

TH1D* Project(TH2D* h, double eta_lo, double eta_hi, const char* name, double scale,
              Color_t col, Style_t mstyle)
{
    const auto yb = PairEtaPanels::Bins(h->GetYaxis(),
                                       {static_cast<float>(eta_lo), static_cast<float>(eta_hi)},
                                       "plot_mc_data_pair_pt_in_eta::Project");
    const int y1 = yb.first, y2 = yb.second;
    TH1D* p = h->ProjectionX(name, y1, y2, "e");
    p->SetDirectory(nullptr);
    p->Scale(scale, "width");          // -> dsigma/dpT [pb/GeV]
    p->SetLineWidth(2);
    p->SetLineColor(col);
    p->SetMarkerColor(col);
    p->SetMarkerStyle(mstyle);
    p->SetMarkerSize(0.9);
    p->SetStats(0);
    p->GetXaxis()->SetTitle("p_{T}^{pair} [GeV]");
    // The observable is the PAIR pT, and the title has to say so: "d#sigma/dp_{T}" alone reads as
    // a single-muon spectrum, which this is not.
    p->GetYaxis()->SetTitle("d#sigma/dp_{T}^{pair} [pb GeV^{-1}]");
    return p;
}

void StyleFrame(TH1D* h, double lo, double hi)
{
    h->GetYaxis()->SetRangeUser(lo, hi);
    h->GetXaxis()->SetTitleSize(0.055);
    h->GetYaxis()->SetTitleSize(0.055);
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetTitleOffset(1.25);
    // On a log x axis the "10^{2}" label is tall; at the default offset the axis TITLE lands on
    // top of it. Push the title clear of the labels.
    h->GetXaxis()->SetTitleOffset(1.15);
}

// -------------------------------------------------------------------------------------------
// A common log-y range for the nine panels that is not mostly empty.
//
// The old rule was [global min x 0.3, global max x 3]. The global minimum is set by a couple of
// nearly-empty bins in the outermost |eta| panels, so the frame reached ~1e-5 while most panels
// stop near 1e-3 -- roughly two of the six drawn decades were blank in every panel. Taking a low
// PERCENTILE of the drawn values instead trims that sparse tail while keeping the bulk of every
// panel inside the frame. The number of points that fall below the frame is printed, so the
// trade-off is auditable rather than silent; it is deliberately NOT written on the canvas.
// -------------------------------------------------------------------------------------------
void CommonLogYRange(const std::vector<TH1D*>& hists, double pct, double& lo, double& hi)
{
    std::vector<double> v;
    for (const TH1D* h : hists)
        for (int b = 1; b <= h->GetNbinsX(); ++b) {
            const double y = h->GetBinContent(b);
            if (y > 0.) v.push_back(y);
        }
    if (v.empty()) { lo = 1e-6; hi = 1.; return; }
    std::sort(v.begin(), v.end());
    const size_t idx = static_cast<size_t>(pct * (v.size() - 1));
    lo = v[idx] * 0.5;
    hi = v.back() * 3.0;

    size_t below = 0;
    for (double y : v) if (y < lo) ++below;
    printf("[range] 9-panel common y = %.4g .. %.4g ; %zu of %zu positive points fall below it\n",
           lo, hi, below, v.size());
}

}  // namespace

// `wp` = "tight" (nominal) or "medium" (WP systematic). Repo rule: every plot set exposes a
// Medium/Tight config var and defaults to Tight. It switches the DATA input file only -- the MC
// here is a TRUTH histogram and has no reconstruction working point. See McDataComprConfig::MuonWP.
void plot_mc_data_pair_pt_in_eta(const char* wp = "tight")
{
    const McDataComprConfig::MuonWP muon_wp = McDataComprConfig::ParseWP(wp);
    std::cout << "plot_mc_data_pair_pt_in_eta: muon WP = " << McDataComprConfig::WPName(muon_wp)
              << std::endl;

    gStyle->SetOptStat(0);
    gSystem->mkdir(kOutDir, kTRUE);

    std::unique_ptr<TH2D> h_data(Fetch(McDataComprConfig::DataFile(muon_wp), kDataHist));
    std::unique_ptr<TH2D> h_mc  (Fetch(McDataComprConfig::PythiaFile(false), kMcHist));
    std::unique_ptr<TH2D> h_pow (Fetch(McDataComprConfig::PowhegFullsimFile(), kPowhegHist));
    // All three sides are booked from the SAME named binnings (`pT_bins_150`, `pair_eta_crossx`),
    // so a divergence is a producer bug, never something to silently plot.
    AssertSameAxes(h_data.get(), h_mc.get());
    AssertSameAxes(h_data.get(), h_pow.get());

    static const CommonEffcyConfig cfg{};
    const auto& eta_bins = cfg.pair_eta_proj_ranges_coarse_incl_gap;

    // THE FIDUCIAL CLAIM IS CHECKED, NOT ASSUMED. The pair-eta-integrated panel is labelled
    // |eta^pair| < ParamsSet::pair_eta_fiducial_max, and the 9 panels tile only that region, so
    // both are honest ONLY if the inputs actually carry the pair-level gap cut. An input filled
    // before that cut existed still has yield beyond +-2.2: it would be silently INCLUDED in the
    // integrated panel (which projects the whole axis) while the label said it was excluded, and
    // silently DROPPED from the 9-panel view. PairEtaPanels::Bins guards panel/axis ALIGNMENT;
    // only this guards EMPTINESS. Mirrors the residual check in plot_pp_counts_pair_pt_in_eta.cxx.
    auto assert_inside_fiducial = [&](TH2D* h, const char* which) {
        const TAxis* ay = h->GetYaxis();
        const auto lo = PairEtaPanels::Bins(ay, eta_bins.front(), "plot_mc_data_pair_pt_in_eta");
        const auto hi = PairEtaPanels::Bins(ay, eta_bins.back(),  "plot_mc_data_pair_pt_in_eta");
        const double all    = h->Integral();
        const double inside = h->Integral(1, h->GetNbinsX(), lo.first, hi.second);
        if (std::fabs(all - inside) > 1e-6 * std::max(1.0, std::fabs(all)))
            throw std::runtime_error(
                std::string(which) + ": " + std::to_string(100. * (all - inside) / all)
                + " % of the yield lies outside |eta^pair| < "
                + std::to_string(ParamsSet::pair_eta_fiducial_max)
                + ", so this file was filled BEFORE the pair-level gap cut was added (2026-09-07). "
                  "The integrated panel would include it while its label says it is excluded, and "
                  "the 9-panel view would drop it. Rerun the hist filling that produced this file.");
    };
    assert_inside_fiducial(h_data.get(), "pp24 data");
    assert_inside_fiducial(h_mc.get(),   "Pythia fullsim");
    assert_inside_fiducial(h_pow.get(),  "POWHEG fullsim");

    printf("data integral = %.6g pb ; Pythia = %.6g nb -> %.6g pb ; Pythia/data = %.4g\n",
           h_data->Integral(), h_mc->Integral(), h_mc->Integral() * kNbToPb,
           h_mc->Integral() * kNbToPb / h_data->Integral());
    printf("POWHEG FullSim pp17 integral = %.6g pb ; POWHEG/data = %.4g\n",
           h_pow->Integral() * kPowhegNorm,
           h_pow->Integral() * kPowhegNorm / h_data->Integral());

    // The full pair-eta range is read from the histogram axis, never retyped. It is the range
    // PROJECTED (the whole +-2.4 axis, which is correct -- nothing is left out); it is NOT the
    // range the spectrum COVERS, because every input carries the pair-level gap cut
    // |eta^pair| < ParamsSet::pair_eta_fiducial_max = 2.2. The panel is therefore LABELLED from
    // the physics, not from the axis -- otherwise a final-results figure would overstate its own
    // acceptance by 0.2 in |eta^pair| (.claude/CLAUDE.md Binnings rule 5).
    const double eta_full_lo = h_data->GetYaxis()->GetXmin();
    const double eta_full_hi = h_data->GetYaxis()->GetXmax();

    // ---------------------------------------------------------------- pair-eta integrated
    {
        TH1D* d = Project(h_data.get(), eta_full_lo, eta_full_hi, "hpt_int_data", 1.0,     McDataComprColors::kSignalData, 20);
        TH1D* m = Project(h_mc.get(),   eta_full_lo, eta_full_hi, "hpt_int_mc",   kNbToPb, McDataComprColors::kSignalMc, 21);
        TH1D* w = Project(h_pow.get(),  eta_full_lo, eta_full_hi, "hpt_int_pow",  kPowhegNorm, McDataComprColors::kSignalPowheg, 22);

        double lo = 1e300, hi = -1e300;
        for (TH1D* h : {d, m, w})
            for (int b = 1; b <= h->GetNbinsX(); ++b) {
                const double v = h->GetBinContent(b);
                if (v > 0.) { lo = std::min(lo, v); hi = std::max(hi, v); }
            }
        // Headroom above the curves so the legend and the eta label are clear of the markers.
        lo *= 0.3; hi *= 8.0;
        StyleFrame(d, lo, hi);
        McDataComprRatio::ApplyLogYLabelPolicy(d, lo, hi);
        McDataComprRatio::HideXAxis(d);

        TCanvas c("c_int", "pair pT, pair-eta integrated", 700, 800);
        TPad* pad_main = nullptr;
        TPad* pad_ratio = nullptr;
        McDataComprRatio::SplitPadForRatio(&c, pad_main, pad_ratio, 0.17, 0.04);

        pad_main->cd();
        gPad->SetLogx(); gPad->SetLogy();
        // DRAW ORDER = POWHEG, Pythia, DATA -- the data LAST so it is on top and can never be
        // hidden by an MC curve where the three cross (user instruction, 2026-09-03). The frame
        // is still defined by the data histogram, which is why it is the one that was styled and
        // range-set above; `d->Draw("E")` first paints the frame, then the MC curves go under
        // ROOT's painting order and the data is re-drawn last on top of them.
        d->Draw("AXIS");
        w->Draw("E,same");
        m->Draw("E,same");
        d->Draw("E,same");
        TLegend l(0.22, 0.10, 0.68, 0.33);
        l.SetBorderSize(0); l.SetFillStyle(0); l.SetTextFont(42); l.SetTextSize(0.040);
        // Legend order = draw order = POWHEG, Pythia, pp24.
        l.AddEntry(w, kPowhegLegend, "lp");
        l.AddEntry((TObject*)nullptr, kPowhegCondLegend, "");
        l.AddEntry(m, kPythiaLegend, "lp");
        l.AddEntry(d, kDataLegend,   "lp");
        l.Draw();
        TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.045);
        t.DrawLatex(0.22, 0.87, Form("|#eta^{pair}| < %.1f", ParamsSet::pair_eta_fiducial_max));

        pad_ratio->cd();
        gPad->SetLogx();
        TH1D* r  = McDataComprRatio::MakeRatio(m, d, "r_pt_int");
        TH1D* rw = McDataComprRatio::MakeRatio(w, d, "r_pt_int_pow");
        double rlo = 0., rhi = 2.; bool rlog = false;
        McDataComprRatio::AutoRatioRange(std::vector<TH1*>{r, rw}, rlo, rhi, rlog);
        printf("[ratio] signal pair_pt (eta integrated) : y-range %.4g .. %.4g%s\n",
               rlo, rhi, rlog ? " (log)" : "");
        gPad->SetLogy(rlog);
        // StyleFrame uses RELATIVE (precision-2) font sizes, so the ratio pad's text must be
        // scaled by the pad-height ratio to come out the same physical size as the main pad's.
        StyleFrame(r, rlo, rhi);
        McDataComprRatio::StyleRatioFrame(r, "p_{T}^{pair} [GeV]", rlo, rhi, rlog,
                                          McDataComprRatio::RelativeTextScale());
        r->Draw("E");
        McDataComprRatio::DrawUnityLine(r);
        rw->Draw("E,same");
        r->Draw("E,same");
        McDataComprRatio::DrawOutOfRangeMarkers(rw, rlo, rhi);
        McDataComprRatio::DrawOutOfRangeMarkers(r,  rlo, rhi);

        c.SaveAs((std::string(kOutDir) + "pair_pt_mc_data_compr.png").c_str());
    }

    // ---------------------------------------------------------------- 9 pair-eta panels
    {
        int nrow = 1, ncol = 1;
        DetermineSubplotGrid(static_cast<int>(eta_bins.size()), nrow, ncol);

        std::vector<TH1D*> ds, ms, ws, all;
        for (size_t i = 0; i < eta_bins.size(); ++i) {
            TH1D* d = Project(h_data.get(), eta_bins[i].first, eta_bins[i].second,
                              ("hd" + std::to_string(i)).c_str(), 1.0,     McDataComprColors::kSignalData, 20);
            TH1D* m = Project(h_mc.get(),   eta_bins[i].first, eta_bins[i].second,
                              ("hm" + std::to_string(i)).c_str(), kNbToPb, McDataComprColors::kSignalMc, 21);
            TH1D* w = Project(h_pow.get(),  eta_bins[i].first, eta_bins[i].second,
                              ("hw" + std::to_string(i)).c_str(), kPowhegNorm, McDataComprColors::kSignalPowheg, 22);
            ds.push_back(d); ms.push_back(m); ws.push_back(w);
            all.push_back(d); all.push_back(m); all.push_back(w);
        }

        // ONE common log-y range across all panels, so the nine are directly comparable by eye.
        double lo = 0., hi = 0.;
        CommonLogYRange(all, 0.03, lo, hi);

        // ONE common ratio range too, for the same reason as the common y range: nine panels each
        // on its own ratio scale cannot be compared by eye, which is the whole point of the view.
        std::vector<TH1*> rs, rws, r_all;
        for (size_t i = 0; i < eta_bins.size(); ++i) {
            rs.push_back(McDataComprRatio::MakeRatio(ms[i], ds[i],
                                                     "r_eta" + std::to_string(i)));
            rws.push_back(McDataComprRatio::MakeRatio(ws[i], ds[i],
                                                      "rw_eta" + std::to_string(i)));
            r_all.push_back(rs.back()); r_all.push_back(rws.back());
        }
        double rlo = 0., rhi = 2.; bool rlog = false;
        McDataComprRatio::AutoRatioRange(r_all, rlo, rhi, rlog);
        printf("[ratio] signal pair_pt_in_eta (9 panels, common) : y-range %.4g .. %.4g%s\n",
               rlo, rhi, rlog ? " (log)" : "");

        // Each cell is 450 x 450 px, not 450 x 350: the bottom 30 % is the ratio pad and the main
        // panel must not lose height to it.
        TCanvas c("c_eta", "pair pT in pair-eta bins", 450 * ncol, 450 * nrow);
        c.Divide(ncol, nrow);
        for (size_t i = 0; i < eta_bins.size(); ++i) {
            TPad* pad_main = nullptr;
            TPad* pad_ratio = nullptr;
            McDataComprRatio::SplitPadForRatio(c.cd(static_cast<int>(i) + 1),
                                               pad_main, pad_ratio, 0.20, 0.04);

            pad_main->cd();
            gPad->SetLogx(); gPad->SetLogy();
            StyleFrame(ds[i], lo, hi);
            McDataComprRatio::ApplyLogYLabelPolicy(ds[i], lo, hi);
            McDataComprRatio::HideXAxis(ds[i]);
            // POWHEG, Pythia, DATA -- data LAST so it sits on top (see the integrated panel).
            ds[i]->Draw("AXIS");
            ws[i]->Draw("E,same");
            ms[i]->Draw("E,same");
            ds[i]->Draw("E,same");

            TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.06);
            t.DrawLatex(0.24, 0.88, Form("#eta^{pair} #in [%.1f, %.1f]",
                                         eta_bins[i].first, eta_bins[i].second));
            if (i == 0) {
                // Upper RIGHT, not lower left: the pair-pT spectrum falls by ~5 decades, so the
                // bottom-left corner that used to hold the legend is exactly where the high-pT
                // tail of all three curves lives, and the POWHEG entry ran straight through its
                // own markers. The upper-right quadrant is empty in every panel by construction.
                TLegend* l = new TLegend(0.55, 0.60, 0.98, 0.88);
                l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextFont(42); l->SetTextSize(0.040);
                l->AddEntry(ws[i], kPowhegLegend, "lp");
                l->AddEntry((TObject*)nullptr, kPowhegCondLegend, "");
                l->AddEntry(ms[i], kPythiaLegend, "lp");
                l->AddEntry(ds[i], kDataLegend,   "lp");
                l->Draw();
            }

            pad_ratio->cd();
            gPad->SetLogx();
            gPad->SetLogy(rlog);
            TH1D* r  = static_cast<TH1D*>(rs[i]);
            TH1D* rw = static_cast<TH1D*>(rws[i]);
            StyleFrame(r, rlo, rhi);
            McDataComprRatio::StyleRatioFrame(r, "p_{T}^{pair} [GeV]", rlo, rhi, rlog,
                                              McDataComprRatio::RelativeTextScale());
            r->Draw("E");
            McDataComprRatio::DrawUnityLine(r);
            rw->Draw("E,same");
            r->Draw("E,same");
            McDataComprRatio::DrawOutOfRangeMarkers(rw, rlo, rhi);
            McDataComprRatio::DrawOutOfRangeMarkers(r,  rlo, rhi);
        }
        c.SaveAs((std::string(kOutDir) + "pair_pt_in_eta_subplots_mc_data_compr.png").c_str());
    }
}
