// plot_muon_q_eta_spectrum.cxx
//
// Single-muon q*eta spectrum after the nominal NTuple-processing selection, as the
// data input for choosing a detector-gap fiducial cut (+ acceptance efficiency).
// See Analysis/docs/tracking/muon_gap_cuts_acceptance.md.
//
// INPUT (NTuple-processing OUTPUT — never the raw NTUPs):
//   single_muon_trees_<sample>_part<N>_<trig>_mindR_0_02.root : TTree "muon_tree".
//   A muon is in that tree iff it belongs to a pair surviving the full nominal pair
//   chain (trigger match, Combined+Medium+IDCuts+MuonCuts, |eta|<2.4, pt>4.5 GeV,
//   one-sided dp/p, d0/z0, photoproduction veto (PbPb), resonance veto), de-duplicated
//   by muon index. It is therefore the trigger-biased single-muon projection of the
//   selected dimuon sample — which is exactly the population a gap cut would act on.
//   NO trigger-efficiency correction is applied here: the raw spectrum is to be read
//   together with the data-driven mu4 trigger efficiencies.
//
// BINNING: q*eta uses ParamsSet::makeEtaTrigEffcyBinning(1), the existing canonical
//   fine q*eta axis (0.01 for |q*eta|<0.2, 0.02 over [-1.3,1.4] and [2.2,2.4], 0.10
//   elsewhere) already used for the 2D trigger-efficiency maps. Not a new binning.
//   Bins are variable-width => everything is plotted as a DENSITY dN/d(q*eta)
//   via Scale(1., "width").
//   CAVEAT for the reader: over |q*eta| in [1.4,2.2] the canonical step is 0.10, coarser
//   than the real structure there, so a cut boundary cannot be read off this figure to
//   better than 0.1 in that region. The canonical binning is NOT changed to suit a plot.
// CENTRALITY: ParamsSet::ctrbins = {0,5,10,20,30,50,80}, with 0-5 and 5-10 merged
//   into 0-10 per the request.
//
// WORKING POINT: use_tight_wp (default TRUE) selects quality&16 (Tight); false selects
//   quality&8 (Medium). NOTE this is a PER-MUON working point, which is NOT identical to
//   the nominal pair-level requirement: the analysis requires the WP on BOTH muons of the
//   pair (m1.quality & m2.quality & 16, DimuonDataAlgCoreT.c:584). The single-muon tree is
//   de-duplicated and carries no partner information, so the per-muon form is the only one
//   available here; ~5% of these muons come from pairs the nominal Tight analysis drops.
//   That is the right population for choosing a SINGLE-MUON fiducial cut, but it is not
//   the nominal pair selection. use_tight_wp=false is effectively a no-op: the trees are
//   already written with pair-level Medium applied, so quality&8 holds for every entry.
//
// GAP-CUT SET (use_new_gap_cuts): which candidate fiducial definition is overlaid, and
//   which one the "fraction removed" numbers refer to. The spectra themselves are
//   IDENTICAL in both cases -- only the overlay and the quoted fractions change.
//   false -> OLD: the legacy PassSingleMuonGapCut (RDFBasedHistFillingData.cxx:141),
//            i.e. |eta| < eta_gap_cut1 at all pT, plus charge_eta_gap_cuts for pT < 6 GeV.
//            Output: muon_gap_cuts/old_gap_cuts/
//   true  -> NEW: the proposed pT-INDEPENDENT fiducial windows in
//            ParamsSet::single_mu_fiducial_gap_cuts. Output: muon_gap_cuts/new_gap_cuts/
//            The windows are READ from ParamsSet -- nothing is retyped or overridden
//            here -- and are printed at run time, so re-tuning the cut means editing
//            ParamsSet.h alone and re-running this macro.
//
// Usage:
//   root -l -b -q 'plot_muon_q_eta_spectrum.cxx+()'                // Tight, old cuts
//   root -l -b -q 'plot_muon_q_eta_spectrum.cxx+(false)'           // Medium, old cuts
//   root -l -b -q 'plot_muon_q_eta_spectrum.cxx+(true,  true)'     // Tight, new cuts
//   root -l -b -q 'plot_muon_q_eta_spectrum.cxx+(false, true)'     // Medium, new cuts

#include "../../MuonObjectsParamsAndHelpers/ParamsSet.h"

#include <TBox.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

namespace {

const char* kOutDirBase =
    "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/single_b_analysis/muon_gap_cuts";
const char* kDataBase = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/";

// Plot axis limits. NOT retyped: set once from ParamsSet::makeEtaTrigEffcyBinning(1)
// at the top of plot_muon_q_eta_spectrum().
double kAxisLo = 0., kAxisHi = 0.;

// ---------------------------------------------------------------------------
// Existing gap definitions in the code base, read from their sources.
// (1) ParamsSet::eta_gap_cut1 : reject |eta| < 0.135, at ALL pT. Since |q*eta| ==
//     |eta| this is the q*eta window (-0.135, +0.135) and can be drawn on the same axis.
// (2) ParamsSet::charge_eta_gap_cuts : reject q*eta inside any window, but ONLY for
//     muons with pT < 6 GeV (RDFBasedHistFillingData.cxx:144).
// Both are used by PassSingleMuonGapCut (RDFBasedHistFillingData.cxx:141), which today
// only feeds the diagnostic "_wgapcut" histograms, not the signal selection.
// ---------------------------------------------------------------------------
const double kGapCutPtThreshold = 6.0;  // RDFBasedHistFillingData.cxx:144

// Windows applied at all pT (the |eta| crack cut).
std::vector<std::pair<double, double>> GapWindowsAllPt() {
    return {{-ParamsSet::eta_gap_cut1, ParamsSet::eta_gap_cut1}};
}
// Windows applied only below kGapCutPtThreshold.
std::vector<std::pair<double, double>> GapWindowsLowPtOnly() {
    std::vector<std::pair<double, double>> w;
    for (const auto& g : ParamsSet::charge_eta_gap_cuts) w.emplace_back(g[0], g[1]);
    std::sort(w.begin(), w.end());
    return w;
}

// Faithful re-implementation of RDFBasedHistFillingData.cxx:141 PassSingleMuonGapCut,
// so the "fraction removed" number below is the real one, pT condition included.
bool PassSingleMuonGapCut(double eta, double pt, int charge) {
    if (std::fabs(eta) < ParamsSet::eta_gap_cut1) return false;
    if (pt < kGapCutPtThreshold)
        for (const auto& g : ParamsSet::charge_eta_gap_cuts)
            if (charge * eta > g[0] && charge * eta < g[1]) return false;
    return true;
}

// ---- the PROPOSED fiducial cut ----
// Windows are READ from ParamsSet::single_mu_fiducial_gap_cuts, the single source of truth;
// the live values are PRINTED at run time and deliberately not repeated here (a literal list in
// this comment went stale twice, 2026-08-12 and 2026-09-07). No value is retyped or overridden
// here, so re-tuning the cut means editing ParamsSet.h alone and re-running. Built once into
// g_new_windows because the test runs per muon over ~10 M entries.
bool kUseNewGapCuts = false;
std::vector<std::pair<double, double>> g_new_windows;

void BuildNewWindows() {
    g_new_windows.clear();
    for (const auto& g : ParamsSet::single_mu_fiducial_gap_cuts)
        g_new_windows.emplace_back(g.first, g.second);
    std::sort(g_new_windows.begin(), g_new_windows.end());
}

// pT-INDEPENDENT by design, so the acceptance factors as a pure function of q*eta.
bool PassNewFiducialCut(double eta, int charge) {
    const double x = charge * eta;
    for (const auto& w : g_new_windows)
        if (x > w.first && x < w.second) return false;
    return true;
}

// Dispatch used by the fill loop and by the drawing helpers.
bool PassSelectedGapCut(double eta, double pt, int charge) {
    return kUseNewGapCuts ? PassNewFiducialCut(eta, charge)
                          : PassSingleMuonGapCut(eta, pt, charge);
}


struct Sample {
    std::string dir;   // e.g. "pbpb_2023/"
    std::string tag;   // e.g. "pbpb_2023"
    std::string trig;  // e.g. "single_mu4"
    int nparts;
};

void AddParts(TChain& ch, const Sample& s) {
    for (int p = 1; p <= s.nparts; ++p)
        ch.Add((std::string(kDataBase) + s.dir + "single_muon_trees_" + s.tag + "_part" +
                std::to_string(p) + "_" + s.trig + "_mindR_0_02.root")
                   .c_str());
}

struct FillCounts {
    Long64_t n_wp = 0;      // muons in the plotted population (WP; PbPb also 0<=cent<80)
    Long64_t n_gapcut = 0;  // ... of which FAIL PassSingleMuonGapCut
    std::vector<Long64_t> n_wp_ctr, n_gapcut_ctr;  // same, per centrality bin
};

// Fill q*eta histograms from one sample.
// A PbPb muon with ev_centrality < 0 (the GetCentrality >85% sentinel) or >= 80 is dropped
// from the centrality-binned histograms, from the centrality-integrated histogram, AND from
// the counters — so every printed number describes exactly the plotted population.
FillCounts FillSample(const Sample& s, bool use_tight_wp, bool is_pbpb, TH1D* h_all,
                      TH1D* h_pos, TH1D* h_neg, TH1D* h_lopt, TH1D* h_hipt,
                      std::vector<TH1D*>& h_ctr, const std::vector<int>& ctr_lo,
                      const std::vector<int>& ctr_hi) {
    TChain ch("muon_tree");
    AddParts(ch, s);

    // SetMakeClass(1): the MuonObj branch is split, so leaves are addressable by name
    // without a compiled dictionary for MuonPbPb / MuonPP.
    ch.SetMakeClass(1);
    Float_t eta = 0.f, pt = 0.f;
    Int_t charge = 0, quality = 0, ev_centrality = -1;
    ch.SetBranchStatus("*", 0);
    for (const char* b : {"eta", "pt", "charge", "quality"}) ch.SetBranchStatus(b, 1);
    ch.SetBranchAddress("eta", &eta);
    ch.SetBranchAddress("pt", &pt);
    ch.SetBranchAddress("charge", &charge);
    ch.SetBranchAddress("quality", &quality);
    if (is_pbpb) {
        ch.SetBranchStatus("ev_centrality", 1);
        ch.SetBranchAddress("ev_centrality", &ev_centrality);
    }

    const int wp_bit = use_tight_wp ? 16 : 8;

    FillCounts fc;
    fc.n_wp_ctr.assign(ctr_lo.size(), 0);
    fc.n_gapcut_ctr.assign(ctr_lo.size(), 0);

    const Long64_t n_tree = ch.GetEntries();
    for (Long64_t i = 0; i < n_tree; ++i) {
        ch.GetEntry(i);
        if ((quality & wp_bit) == 0) continue;
        if (is_pbpb && (ev_centrality < 0 || ev_centrality >= ctr_hi.back())) continue;

        ++fc.n_wp;
        const bool fails_gap = !PassSelectedGapCut(eta, pt, charge);
        if (fails_gap) ++fc.n_gapcut;

        const double q_eta = charge * eta;
        if (is_pbpb)
            for (size_t k = 0; k < ctr_lo.size(); ++k)
                if (ev_centrality >= ctr_lo[k] && ev_centrality < ctr_hi[k]) {
                    h_ctr[k]->Fill(q_eta);
                    ++fc.n_wp_ctr[k];
                    if (fails_gap) ++fc.n_gapcut_ctr[k];
                }
        h_all->Fill(q_eta);
        (charge > 0 ? h_pos : h_neg)->Fill(q_eta);
        (pt < kGapCutPtThreshold ? h_lopt : h_hipt)->Fill(q_eta);
    }
    printf("  %-12s : %10lld muons in tree, %10lld in the plotted population (%s%s), "
           "%10lld (%.4f) fail the selected gap cut\n",
           s.tag.c_str(), n_tree, fc.n_wp, use_tight_wp ? "Tight" : "Medium",
           is_pbpb ? ", 0-80%" : "", fc.n_gapcut,
           fc.n_wp > 0 ? double(fc.n_gapcut) / fc.n_wp : 0.);
    return fc;
}

// ------------------------------ drawing helpers ------------------------------

void StyleSpectrum(TH1D* h, int color) {
    h->SetLineColor(color);
    h->SetLineWidth(2);
    h->SetStats(0);
    h->GetXaxis()->SetTitle("q #times #eta^{#mu}");
    h->GetYaxis()->SetTitle("dN_{#mu} / d(q#times#eta)");
    h->GetXaxis()->SetTitleSize(0.050);
    h->GetYaxis()->SetTitleSize(0.050);
    h->GetXaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetTitleOffset(1.70);
}

// Legends must be opaque: they sit over the shaded gap bands.
void StyleLegend(TLegend& l, double text_size) {
    l.SetBorderSize(0);
    l.SetFillColor(kWhite);
    l.SetFillStyle(1001);
    l.SetTextSize(text_size);
}

// Height a legend needs so its entries do not print on top of each other. ROOT divides the
// box evenly over the entries, so a box sized independently of the entry count silently
// overlaps its own text once the count grows (which is exactly what happened to the pp-vs-
// PbPb legend after the gap-window lines were added to it).
double LegendHeight(int n_entries, double text_size) {
    return n_entries * 1.8 * text_size + 0.02;
}

// Where to put a legend so it hides no curve.
//
// The legends here MUST be opaque (they sit over the shaded gap bands) and they are drawn
// LAST, so wherever one lands it erases whatever passes underneath. A fixed corner is
// therefore only safe for one particular spectrum: the same corner that is empty in pp is
// crossed by the Pb+Pb acceptance edge, which is how the Pb+Pb panels ended up with curve
// segments hidden. This scores the four in-frame corners by how much of the drawn curve
// falls inside the box and returns the emptiest, so the placement follows the data.
//
// MUST be called after the frame exists and gPad->Update() has run -- it reads the pad's
// user range. Scoring uses consecutive-bin SEGMENTS, not bin contents: a step histogram
// draws the vertical connector between two bins, so a steep plunge can cross a box without
// any single bin content landing inside it.
struct LegBox {
    double x1, y1, x2, y2;
};

LegBox AutoLegendBox(const std::vector<TH1D*>& hs, double w, double h) {
    const double lm = gPad->GetLeftMargin(), rm = gPad->GetRightMargin();
    const double bm = gPad->GetBottomMargin(), tm = gPad->GetTopMargin();
    const double inset = 0.02;
    const double ux0 = gPad->GetUxmin(), ux1 = gPad->GetUxmax();
    const double uy0 = gPad->GetUymin(), uy1 = gPad->GetUymax();  // log10 units if logy
    const bool logy = gPad->GetLogy();

    // Conventional order, so a tie keeps the familiar top-right placement.
    const std::vector<std::pair<double, double>> cands = {
        {1 - rm - inset - w, 1 - tm - inset - h},  // top right
        {lm + inset, 1 - tm - inset - h},          // top left
        {1 - rm - inset - w, bm + inset},          // bottom right
        {lm + inset, bm + inset}};                 // bottom left

    auto to_x = [&](double xn) { return ux0 + (xn - lm) / (1 - lm - rm) * (ux1 - ux0); };
    auto to_y = [&](double yn) {
        const double v = uy0 + (yn - bm) / (1 - bm - tm) * (uy1 - uy0);
        return logy ? std::pow(10., v) : v;
    };

    size_t best = 0;
    int best_score = -1;
    for (size_t i = 0; i < cands.size(); ++i) {
        const double bx0 = to_x(cands[i].first), bx1 = to_x(cands[i].first + w);
        const double by0 = to_y(cands[i].second), by1 = to_y(cands[i].second + h);
        int hits = 0;
        for (const TH1D* hh : hs)
            for (int b = 1; b < hh->GetNbinsX(); ++b) {
                const double xa = hh->GetBinCenter(b), xb = hh->GetBinCenter(b + 1);
                if (xb < bx0 || xa > bx1) continue;
                const double ca = hh->GetBinContent(b), cb = hh->GetBinContent(b + 1);
                if (std::max(ca, cb) >= by0 && std::min(ca, cb) <= by1) ++hits;
            }
        if (best_score < 0 || hits < best_score) {
            best_score = hits;
            best = i;
        }
    }
    return {cands[best].first, cands[best].second, cands[best].first + w,
            cands[best].second + h};
}

// Shaded bands for the gap-cut windows. Drawn under the curves, which are then redrawn on
// top. Darker band = rejected at all pT; lighter = rejected only for pT < 6 GeV.
// NOTE (2026-08-12): this used to also draw dashed edges bounding the q*eta regions with no
// fitted trigger efficiency. Those regions no longer exist -- the trigger-efficiency code
// now fits in the CONTIGUOUS CommonEffcyConfig::q_eta_proj_ranges_coarse_incl_gap binning,
// which includes the gaps and whose top edge tracks the fiducial cut, so every surviving
// muon has a fitted turn-on and there is nothing to mark.
void DrawGapMarkers(const std::vector<TH1D*>& redraw) {
    const double ymin = gPad->GetUymin(), ymax = gPad->GetUymax();
    const double y0 = gPad->GetLogy() ? std::pow(10., ymin) : ymin;
    const double y1 = gPad->GetLogy() ? std::pow(10., ymax) : ymax;

    if (kUseNewGapCuts) {
        for (const auto& w : g_new_windows) {
            TBox* b = new TBox(w.first, y0, w.second, y1);
            b->SetFillColorAlpha(kGreen - 7, 0.45);
            b->SetLineColor(kGreen - 7);
            b->Draw("same");
        }
    } else {
        for (const auto& w : GapWindowsAllPt()) {
            TBox* b = new TBox(w.first, y0, w.second, y1);
            b->SetFillColorAlpha(kRed - 7, 0.45);
            b->SetLineColor(kRed - 7);
            b->Draw("same");
        }
        for (const auto& w : GapWindowsLowPtOnly()) {
            TBox* b = new TBox(w.first, y0, w.second, y1);
            b->SetFillColorAlpha(kOrange - 9, 0.45);
            b->SetLineColor(kOrange - 9);
            b->Draw("same");
        }
    }
    for (TH1D* h : redraw) h->Draw("hist same");
    gPad->RedrawAxis();
}

// Common LINEAR-y range for a set of overlaid histograms. Floor pinned at 0 -- a linear
// axis must start at zero or the apparent depth of the gap structure this figure exists
// to show is visually distorted -- and ceiling at kHeadroom above the peak so the curve
// never touches the frame. Unlike a log floor, nothing is ever pushed off scale.
const double kHeadroom = 1.15;
void SetCommonLinRange(const std::vector<TH1D*>& hs) {
    double hi = 0.;
    for (TH1D* h : hs)
        for (int i = 1; i <= h->GetNbinsX(); ++i)
            hi = std::max(hi, h->GetBinContent(i));
    for (TH1D* h : hs) {
        h->SetMinimum(0.);
        h->SetMaximum(hi * kHeadroom);
    }
}

void Header(const std::string& text, const std::string& sub) {
    TLatex t;
    t.SetNDC();
    t.SetTextFont(42);
    // Both lines must clear the frame: baselines are placed above 1-topMargin. The
    // near-frame line (the only line when sub is empty; the sub line otherwise) is pushed
    // up to +0.055/+0.050 rather than the tighter +0.025/+0.020 used when the axis was log:
    // on a linear axis with large content, ROOT draws a "x10^n" multiplier box just above
    // the frame, and the tighter offsets used to sit right on top of it.
    const double frame_top = 1.0 - gPad->GetTopMargin();
    t.SetTextSize(0.042);
    t.DrawLatex(gPad->GetLeftMargin(), frame_top + (sub.empty() ? 0.055 : 0.095),
                text.c_str());
    if (!sub.empty()) {
        t.SetTextSize(0.035);
        t.DrawLatex(gPad->GetLeftMargin(), frame_top + 0.050, sub.c_str());
    }
}

// The key for the gap-band overlays. Every canvas that draws the overlays
// needs it (a reader of one figure alone cannot otherwise decode the shading). Values are
// FORMATTED from the constants, never retyped, so the text cannot drift from the bands.
// compact=true puts all proposed windows on ONE entry -- fine for the full-canvas key,
// where there is width; the per-pad legend (1/3 canvas wide) needs one line per window.
void AddGapKeyEntries(TLegend& l, bool compact = false) {
    if (kUseNewGapCuts) {
        // one line per window: the full list on a single line overflows a 1/3-width pad
        TBox* b = new TBox();
        b->SetFillColorAlpha(kGreen - 7, 0.45);
        if (compact) {
            std::string txt = "rejected at all p_{T}, q#times#eta in:";
            for (size_t i = 0; i < g_new_windows.size(); ++i)
                txt += Form("%s (%.2f, %.2f)", i ? "," : "", g_new_windows[i].first,
                            g_new_windows[i].second);
            l.AddEntry(b, txt.c_str(), "f");
        } else {
            l.AddEntry(b, "rejected (all p_{T}):", "f");
            for (const auto& w : g_new_windows)
                l.AddEntry((TObject*)nullptr,
                           Form("   %.2f < q#times#eta < %.2f", w.first, w.second), "");
        }
    } else {
        TBox* b1 = new TBox();
        b1->SetFillColorAlpha(kRed - 7, 0.45);
        l.AddEntry(b1,
                   Form("rejected: |#eta| < %.3f  (all p_{T})", ParamsSet::eta_gap_cut1),
                   "f");
        TBox* b2 = new TBox();
        b2->SetFillColorAlpha(kOrange - 9, 0.45);
        l.AddEntry(b2, Form("rejected: q#times#eta windows  (p_{T} < %.0f GeV only)",
                            kGapCutPtThreshold),
                   "f");
    }
}

// The gap key drawn ONCE across the top of a multi-panel canvas: in a half-width panel a
// per-pad legend clips these strings, and the overlay applies to every panel anyway.
void CanvasGapKey(TCanvas& c, double y1, double y2, double text_size = 0.011) {
    c.cd();
    TLegend* l = new TLegend(0.06, y1, 0.98, y2);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetTextSize(text_size);
    // Default margin reserves a QUARTER of the box width for the colour swatch; across a
    // full-canvas box that is a green slab wider than the text it labels.
    l->SetMargin(0.03);
    l->SetNColumns(kUseNewGapCuts ? 1 : 2);
    AddGapKeyEntries(*l, /*compact=*/true);
    l->Draw();
}

// Sample headline drawn ONCE across the top of the canvas. The per-panel pads are only
// 1/3 of the canvas wide, which clips a full "Pb+Pb sqrt(s_NN) = ..., HLT_mu4" string;
// the sample/working point is common to all panels anyway, so it belongs to the canvas.
void CanvasHeadline(TCanvas& c, const std::string& text, const std::string& sub) {
    c.cd();
    TLatex t;
    t.SetNDC();
    t.SetTextFont(42);
    t.SetTextSize(0.038);
    t.DrawLatex(0.065, 0.955, text.c_str());
    if (!sub.empty()) {
        t.SetTextSize(0.032);
        t.DrawLatex(0.065, 0.915, sub.c_str());
    }
}

// A main+ratio pad pair inside the rectangle [x1,x2] x [y1,y2] of the current canvas.
std::pair<TPad*, TPad*> MakeRatioPads(TCanvas& c, const std::string& name, double x1,
                                      double y1, double x2, double y2) {
    c.cd();  // pads belong to the CANVAS, not to whatever sub-pad is current
    const double ysplit = y1 + 0.34 * (y2 - y1);
    TPad* p_main = new TPad((name + "_main").c_str(), "", x1, ysplit, x2, y2);
    TPad* p_rat = new TPad((name + "_rat").c_str(), "", x1, y1, x2, ysplit);
    p_main->SetLeftMargin(0.20);
    p_main->SetRightMargin(0.04);
    p_main->SetTopMargin(0.14);
    p_main->SetBottomMargin(0.02);
    p_main->SetTickx(1);
    p_main->SetTicky(1);
    p_rat->SetLeftMargin(0.20);  // MUST match p_main's or the x-axes will not line up
    p_rat->SetRightMargin(0.04);
    p_rat->SetTopMargin(0.03);
    p_rat->SetBottomMargin(0.34);
    p_rat->SetTickx(1);
    p_rat->SetTicky(1);
    p_main->Draw();
    p_rat->Draw();
    return {p_main, p_rat};
}

// A single full-height pad in [x1,x2] x [y1,y2].
TPad* MakeSinglePad(TCanvas& c, const std::string& name, double x1, double y1, double x2,
                    double y2) {
    c.cd();  // pads belong to the CANVAS, not to whatever sub-pad is current
    TPad* p = new TPad(name.c_str(), "", x1, y1, x2, y2);
    p->SetLeftMargin(0.20);
    p->SetRightMargin(0.04);
    p->SetTopMargin(0.14);
    p->SetBottomMargin(0.14);
    p->SetTickx(1);
    p->SetTicky(1);
    p->Draw();
    return p;
}

// Build num/den. The two histograms are always DISJOINT samples here (mu+ vs mu-,
// low-pT vs high-pT, one centrality vs another, pp vs PbPb), so the default independent
// error propagation of TH1::Divide is correct — unlike a subset/superset ratio, which
// needs the conditional (binomial) form (cf. mc_trigger_efficiency.md R12).
TH1D* MakeRatio(const TH1D* num, const TH1D* den, const std::string& name, int color,
                const std::string& ytitle) {
    TH1D* r = (TH1D*)num->Clone(name.c_str());
    r->Divide(den);
    r->SetStats(0);
    r->SetLineColor(color);
    r->SetMarkerColor(color);
    r->SetLineWidth(2);
    r->GetYaxis()->SetTitle(ytitle.c_str());
    r->GetXaxis()->SetTitle("q #times #eta^{#mu}");
    // The ratio pad is ~1/3 the height of the main pad, so text/ticks need scaling up.
    r->GetXaxis()->SetTitleSize(0.115);
    r->GetXaxis()->SetLabelSize(0.100);
    r->GetXaxis()->SetTitleOffset(1.20);
    r->GetYaxis()->SetTitleSize(0.095);
    r->GetYaxis()->SetLabelSize(0.090);
    r->GetYaxis()->SetTitleOffset(0.80);
    r->GetYaxis()->SetNdivisions(505);
    return r;
}

// Y-range covering every plotted ratio point INCLUDING its error bar, so nothing is
// cropped (the failure mode fixed in raa_from_rdf_crossx.md).
void SetRatioRange(const std::vector<TH1D*>& rs) {
    double lo = 1e300, hi = -1e300;
    for (TH1D* r : rs)
        for (int i = 1; i <= r->GetNbinsX(); ++i) {
            if (r->GetBinContent(i) == 0 && r->GetBinError(i) == 0) continue;
            lo = std::min(lo, r->GetBinContent(i) - r->GetBinError(i));
            hi = std::max(hi, r->GetBinContent(i) + r->GetBinError(i));
        }
    if (lo > hi) return;
    // A ratio panel must show its unity reference; without this the dashed y=1 line is
    // clipped to the PAD (not the frame) and renders across the x-axis tick labels.
    lo = std::min(lo, 1.0);
    hi = std::max(hi, 1.0);
    const double pad = 0.08 * (hi - lo) + 1e-6;
    for (TH1D* r : rs) {
        r->SetMinimum(lo - pad);
        r->SetMaximum(hi + pad);
    }
}

void DrawRatioUnityLine() {
    TLine* l = new TLine(kAxisLo, 1., kAxisHi, 1.);
    l->SetLineStyle(2);
    l->SetLineColor(kGray + 2);
    l->Draw("same");
}

// Report the deepest bins inside/near the candidate gap regions, so the exact
// dip positions are readable without eyeballing the PNG.
void ReportDips(const TH1D* h, const char* label) {
    printf("\n--- %s : local minima of dN/d(q*eta) (bins below 60%% of the "
           "local 21-bin median) ---\n",
           label);
    const int n = h->GetNbinsX();
    for (int i = 1; i <= n; ++i) {
        std::vector<double> loc;
        for (int j = std::max(1, i - 10); j <= std::min(n, i + 10); ++j)
            loc.push_back(h->GetBinContent(j));
        std::sort(loc.begin(), loc.end());
        const double med = loc[loc.size() / 2];
        const double c = h->GetBinContent(i);
        if (med > 0 && c < 0.60 * med)
            printf("   q*eta in [%+.3f, %+.3f]  density = %10.1f   local median = "
                   "%10.1f   ratio = %.2f\n",
                   h->GetBinLowEdge(i), h->GetBinLowEdge(i + 1), c, med, c / med);
    }
}

}  // namespace

void plot_muon_q_eta_spectrum(bool use_tight_wp = true, bool use_new_gap_cuts = false) {
    gROOT->SetBatch(true);
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    const std::string wp_tag = use_tight_wp ? "" : "_medium_wp";
    const std::string wp_lbl = use_tight_wp ? "Tight" : "Medium";

    kUseNewGapCuts = use_new_gap_cuts;
    BuildNewWindows();
    std::string outdir = std::string(kOutDirBase) +
                         (use_new_gap_cuts ? "/new_gap_cuts" : "/old_gap_cuts") +
                         (use_tight_wp ? "" : "/medium");
    gSystem->mkdir(outdir.c_str(), kTRUE);

    printf("gap-cut set: %s\n", use_new_gap_cuts ? "NEW (proposed fiducial)"
                                                 : "OLD (legacy PassSingleMuonGapCut)");
    if (use_new_gap_cuts) {
        printf("  windows from ParamsSet::single_mu_fiducial_gap_cuts =");
        for (const auto& w : g_new_windows) printf(" (%.2f, %.2f)", w.first, w.second);
        printf("\n");
    }

    // ---- canonical binnings, read from their single sources of truth ----
    const std::vector<double> qeta_bins = ParamsSet::makeEtaTrigEffcyBinning(1);
    const int nb = static_cast<int>(qeta_bins.size()) - 1;
    kAxisLo = qeta_bins.front();  // never retyped
    kAxisHi = qeta_bins.back();
    printf("q*eta binning: %d bins over [%.2f, %.2f] "
           "(ParamsSet::makeEtaTrigEffcyBinning(1))\n",
           nb, qeta_bins.front(), qeta_bins.back());

    // ParamsSet::ctrbins = {0,5,10,20,30,50,80}; merge 0-5 and 5-10 into 0-10.
    std::vector<int> ctr_lo, ctr_hi;
    // Merge the FIRST TWO bins (by index, not by value) into one, per the request.
    for (size_t i = 1; i + 1 < ParamsSet::ctrbins.size(); ++i) {
        ctr_lo.push_back(i == 1 ? ParamsSet::ctrbins.front() : ParamsSet::ctrbins[i]);
        ctr_hi.push_back(ParamsSet::ctrbins[i + 1]);
    }
    const int nctr = static_cast<int>(ctr_lo.size());
    printf("centrality bins (ParamsSet::ctrbins, 0-5 + 5-10 merged):");
    for (int k = 0; k < nctr; ++k) printf(" %d-%d%%", ctr_lo[k], ctr_hi[k]);
    printf("\n");

    auto book = [&](const char* name) {
        TH1D* h = new TH1D(name, "", nb, qeta_bins.data());
        h->Sumw2();
        return h;
    };

    // ======================= pp24 =======================
    printf("\n=== pp24 ===\n");
    TH1D *pp_all = book("h_pp_all"), *pp_pos = book("h_pp_pos"), *pp_neg = book("h_pp_neg");
    TH1D *pp_lopt = book("h_pp_lopt"), *pp_hipt = book("h_pp_hipt");
    std::vector<TH1D*> pp_dummy;
    const FillCounts pp_fc =
        FillSample({"pp_2024/", "pp_2024", "2mu4", 12}, use_tight_wp, false, pp_all,
                   pp_pos, pp_neg, pp_lopt, pp_hipt, pp_dummy, ctr_lo, ctr_hi);

    // ================= PbPb 23+24+25 ====================
    printf("\n=== PbPb 2023 + 2024 + 2025 (combined) ===\n");
    TH1D *pb_all = book("h_pb_all"), *pb_pos = book("h_pb_pos"), *pb_neg = book("h_pb_neg");
    TH1D *pb_lopt = book("h_pb_lopt"), *pb_hipt = book("h_pb_hipt");
    std::vector<TH1D*> pb_ctr(nctr);
    for (int k = 0; k < nctr; ++k)
        pb_ctr[k] = book(Form("h_pb_ctr%d_%d", ctr_lo[k], ctr_hi[k]));

    const std::vector<Sample> pbpb = {{"pbpb_2023/", "pbpb_2023", "single_mu4", 4},
                                      {"pbpb_2024/", "pbpb_2024", "single_mu4", 2},
                                      {"pbpb_2025/", "pbpb_2025", "single_mu4", 6}};
    FillCounts pb_fc;
    pb_fc.n_wp_ctr.assign(nctr, 0);
    pb_fc.n_gapcut_ctr.assign(nctr, 0);
    for (const auto& s : pbpb) {
        const FillCounts fc = FillSample(s, use_tight_wp, true, pb_all, pb_pos, pb_neg,
                                         pb_lopt, pb_hipt, pb_ctr, ctr_lo, ctr_hi);
        pb_fc.n_wp += fc.n_wp;
        pb_fc.n_gapcut += fc.n_gapcut;
        for (int k = 0; k < nctr; ++k) {
            pb_fc.n_wp_ctr[k] += fc.n_wp_ctr[k];
            pb_fc.n_gapcut_ctr[k] += fc.n_gapcut_ctr[k];
        }
    }

    // ---- densities (variable-width bins) ----
    for (TH1D* h : {pp_all, pp_pos, pp_neg, pp_lopt, pp_hipt, pb_all, pb_pos, pb_neg,
                    pb_lopt, pb_hipt})
        h->Scale(1., "width");
    for (TH1D* h : pb_ctr) h->Scale(1., "width");

    // ========== Plots A (pp24) and B (PbPb combined, centrality-integrated) ==========
    // Three columns (N <= 3 => one row): the spectrum with the existing gap windows
    // overlaid; the same split by muon charge; and split at pT = 6 GeV, above which the
    // existing q*eta windows are NOT applied. Columns 2 and 3 carry a ratio pad, because
    // the ratio is exactly what those comparisons are about.
    struct PanelSet {
        TH1D *all, *pos, *neg, *lopt, *hipt;
        std::string header, sub, fname;
    };
    const std::vector<PanelSet> sets = {
        {pp_all, pp_pos, pp_neg, pp_lopt, pp_hipt,
         "pp #sqrt{s} = 5.36 TeV, 2024, HLT_2mu4",
         "single muons of selected pairs, " + wp_lbl + " WP (per muon)",
         "muon_q_eta_spectrum_pp24"},
        {pb_all, pb_pos, pb_neg, pb_lopt, pb_hipt,
         "Pb+Pb #sqrt{s_{NN}} = 5.36 TeV, 2023+2024+2025, HLT_mu4",
         "0-80% centrality, single muons of selected pairs, " + wp_lbl + " WP (per muon)",
         "muon_q_eta_spectrum_pbpb_combined"}};

    for (const auto& ps : sets) {
        TCanvas c(("c_" + ps.fname).c_str(), "", 2100, 840);

        // --- column 1: the spectrum + gap windows (single curve, no ratio) ---
        MakeSinglePad(c, ps.fname + "_p1", 0.000, 0.0, 0.334, 0.875)->cd();
        StyleSpectrum(ps.all, kBlack);
        SetCommonLinRange({ps.all});
        ps.all->Draw("hist");
        gPad->Update();
        DrawGapMarkers({ps.all});
        Header(kUseNewGapCuts ? "all muons; proposed gap cuts overlaid"
                              : "all muons; existing gap cuts overlaid",
               "");
        // NO in-frame legend here. It used to carry the whole gap-window key, which made it
        // the width of the frame and five entries tall -- with nowhere left to put it, it sat
        // on the Pb+Pb acceptance edge and on the q*eta = 0 dip. The key is now drawn once
        // across the canvas (below), and this panel has a single curve that the header names.

        // --- column 2: by muon charge, with mu+/mu- ratio ---
        {
            auto pads = MakeRatioPads(c, ps.fname + "_p2", 0.334, 0.0, 0.667, 0.875);
            pads.first->cd();
            StyleSpectrum(ps.pos, kRed + 1);
            StyleSpectrum(ps.neg, kAzure + 2);
            for (TH1D* h : {ps.pos, ps.neg}) { h->GetXaxis()->SetLabelSize(0.); h->GetXaxis()->SetTitleSize(0.); }
            SetCommonLinRange({ps.pos, ps.neg});
            ps.pos->Draw("hist");
            ps.neg->Draw("hist same");
            gPad->Update();
            DrawGapMarkers({ps.pos, ps.neg});
            Header("by muon charge", "");
            {
                const LegBox lb = AutoLegendBox({ps.pos, ps.neg}, 0.23,
                                                LegendHeight(2, 0.050));
                TLegend l(lb.x1, lb.y1, lb.x2, lb.y2);
                StyleLegend(l, 0.050);
                l.AddEntry(ps.pos, "#mu^{+}", "l");
                l.AddEntry(ps.neg, "#mu^{-}", "l");
                l.DrawClone();
            }
            pads.second->cd();
            TH1D* r = MakeRatio(ps.pos, ps.neg, ps.fname + "_r2", kBlack,
                                "#mu^{+} / #mu^{-}");
            SetRatioRange({r});
            r->Draw("hist");
            gPad->Update();
            DrawGapMarkers({r});
            DrawRatioUnityLine();
        }

        // --- column 3: split at pT = 6 GeV, with the ratio ---
        {
            auto pads = MakeRatioPads(c, ps.fname + "_p3", 0.667, 0.0, 1.000, 0.875);
            pads.first->cd();
            StyleSpectrum(ps.lopt, kMagenta + 1);
            StyleSpectrum(ps.hipt, kGreen + 2);
            for (TH1D* h : {ps.lopt, ps.hipt}) { h->GetXaxis()->SetLabelSize(0.); h->GetXaxis()->SetTitleSize(0.); }
            SetCommonLinRange({ps.lopt, ps.hipt});
            ps.hipt->Draw("hist");
            ps.lopt->Draw("hist same");
            gPad->Update();
            DrawGapMarkers({ps.hipt, ps.lopt});
            Header("split at p_{T} = 6 GeV", "");
            {
                const LegBox lb = AutoLegendBox({ps.lopt, ps.hipt}, 0.32,
                                                LegendHeight(2, 0.045));
                TLegend l(lb.x1, lb.y1, lb.x2, lb.y2);
                StyleLegend(l, 0.045);
                l.AddEntry(ps.lopt, "p_{T} < 6 GeV", "l");
                l.AddEntry(ps.hipt, "p_{T} #geq 6 GeV", "l");
                l.DrawClone();
            }
            pads.second->cd();
            TH1D* r = MakeRatio(ps.lopt, ps.hipt, ps.fname + "_r3", kBlack,
                                "#frac{p_{T} < 6 GeV}{p_{T} #geq 6 GeV}");
            SetRatioRange({r});
            r->Draw("hist");
            gPad->Update();
            DrawGapMarkers({r});
            DrawRatioUnityLine();
        }

        // One key for all three panels, in the band freed above them -- the overlay is the
        // same in every panel, and no in-frame box can be placed without covering data.
        CanvasGapKey(c, 0.878, 0.906, 0.020);
        CanvasHeadline(c, ps.header, ps.sub);
        c.SaveAs((outdir + "/" + ps.fname + wp_tag + ".png").c_str());
    }

    // ============ Plot C : PbPb combined, centrality-binned ============
    // 5 centrality panels + 1 shape-overlay panel; layout nrows >= ncols.
    // Each centrality panel gets its OWN y-range so the gap structure fills the frame;
    // cross-centrality comparison is the job of the unit-area overlay + its ratio pad.
    {
        TCanvas c("cC", "", 1400, 1760);
        const int ncol = 2, nrow = 3;
        const double kPanelTop = 0.965;  // leave a band for the canvas-level gap key
        for (int k = 0; k < nctr; ++k) {
            const int col = k % ncol, row = k / ncol;
            MakeSinglePad(c, Form("cC_p%d", k), col / double(ncol),
                          kPanelTop * (1.0 - (row + 1) / double(nrow)),
                          (col + 1) / double(ncol),
                          kPanelTop * (1.0 - row / double(nrow)))
                ->cd();
            StyleSpectrum(pb_ctr[k], kBlack);
            SetCommonLinRange({pb_ctr[k]});
            pb_ctr[k]->Draw("hist");
            gPad->Update();
            DrawGapMarkers({pb_ctr[k]});
            Header(Form("Pb+Pb #sqrt{s_{NN}} = 5.36 TeV, 23+24+25, %d-%d%%", ctr_lo[k],
                        ctr_hi[k]),
                   wp_lbl + " WP (per muon)");
        }

        // last panel: unit-normalised shape overlay + ratio to the 0-10% shape
        const int k = nctr, col = k % ncol, row = k / ncol;
        auto pads = MakeRatioPads(c, "cC_shape", col / double(ncol),
                                  kPanelTop * (1.0 - (row + 1) / double(nrow)),
                                  (col + 1) / double(ncol),
                                  kPanelTop * (1.0 - row / double(nrow)));
        const int cols[6] = {kBlack, kRed + 1, kOrange + 7, kGreen + 2, kAzure + 2,
                             kMagenta + 1};
        std::vector<TH1D*> shapes;
        for (int j = 0; j < nctr; ++j) {
            TH1D* h = (TH1D*)pb_ctr[j]->Clone(Form("h_shape_%d", j));
            const double integ = h->Integral("width");
            if (integ > 0) h->Scale(1. / integ);
            StyleSpectrum(h, cols[j % 6]);
            h->GetXaxis()->SetLabelSize(0.);
            h->GetXaxis()->SetTitleSize(0.);
            h->GetYaxis()->SetTitle("normalised dN_{#mu} / d(q#times#eta)");
            shapes.push_back(h);
        }
        pads.first->cd();
        SetCommonLinRange(shapes);
        for (size_t j = 0; j < shapes.size(); ++j)
            shapes[j]->Draw(j == 0 ? "hist" : "hist same");
        gPad->Update();
        DrawGapMarkers(shapes);
        Header("shape vs centrality (unit area)", "");
        {
            const LegBox lb = AutoLegendBox(shapes, 0.27, LegendHeight(nctr, 0.042));
            TLegend l(lb.x1, lb.y1, lb.x2, lb.y2);
            StyleLegend(l, 0.042);
            for (int j = 0; j < nctr; ++j)
                l.AddEntry(shapes[j], Form("%d-%d%%", ctr_lo[j], ctr_hi[j]), "l");
            l.DrawClone();
        }
        pads.second->cd();
        std::vector<TH1D*> srat;
        for (int j = 1; j < nctr; ++j)
            srat.push_back(MakeRatio(shapes[j], shapes[0], Form("h_shape_r%d", j),
                                     cols[j % 6],
                                     Form("shape / %d-%d%%", ctr_lo[0], ctr_hi[0])));
        SetRatioRange(srat);
        for (size_t j = 0; j < srat.size(); ++j)
            srat[j]->Draw(j == 0 ? "hist" : "hist same");
        gPad->Update();
        DrawGapMarkers(srat);
        DrawRatioUnityLine();

        CanvasGapKey(c, kPanelTop + 0.002, 0.999);
        c.SaveAs((outdir + "/muon_q_eta_spectrum_pbpb_combined_ctr_binned" + wp_tag +
                  ".png")
                     .c_str());
    }

    // ============ Plot D : pp24 vs PbPb shape comparison, with ratio ============
    {
        TCanvas c("cD", "", 950, 850);
        auto pads = MakeRatioPads(c, "cD", 0.0, 0.0, 1.0, 0.945);
        TH1D* pp_s = (TH1D*)pp_all->Clone("h_pp_shape");
        TH1D* pb_s = (TH1D*)pb_all->Clone("h_pb_shape");
        for (TH1D* h : {pp_s, pb_s}) {
            const double integ = h->Integral("width");
            if (integ > 0) h->Scale(1. / integ);
        }
        StyleSpectrum(pp_s, kBlack);
        StyleSpectrum(pb_s, kRed + 1);
        // AFTER StyleSpectrum, which would otherwise restore the x labels onto the
        // main pad (they belong to the ratio pad below).
        for (TH1D* h : {pp_s, pb_s}) {
            h->GetYaxis()->SetTitle("normalised dN_{#mu} / d(q#times#eta)");
            h->GetXaxis()->SetLabelSize(0.);
            h->GetXaxis()->SetTitleSize(0.);
        }
        pads.first->cd();
        SetCommonLinRange({pp_s, pb_s});
        pp_s->Draw("hist");
        pb_s->Draw("hist same");
        gPad->Update();
        DrawGapMarkers({pp_s, pb_s});
        Header("q#times#eta shape: pp vs Pb+Pb", "unit area, " + wp_lbl + " WP (per muon)");
        {
            // Two entries only: the gap-window key moved to the canvas band. Kept in this
            // box it needed five lines in a box sized for two, so the window text printed
            // over the sample text AND the box sat across the q*eta > 2.3 acceptance edge.
            const LegBox lb = AutoLegendBox({pp_s, pb_s}, 0.49, LegendHeight(2, 0.038));
            TLegend l(lb.x1, lb.y1, lb.x2, lb.y2);
            StyleLegend(l, 0.038);
            l.AddEntry(pp_s, "pp 2024 (HLT_2mu4)", "l");
            l.AddEntry(pb_s, "Pb+Pb 23+24+25, 0-80% (HLT_mu4)", "l");
            l.DrawClone();
        }
        pads.second->cd();
        TH1D* r = MakeRatio(pp_s, pb_s, "h_pppb_ratio", kBlack, "pp / Pb+Pb");
        SetRatioRange({r});
        r->Draw("hist");
        gPad->Update();
        DrawGapMarkers({r});
        DrawRatioUnityLine();
        CanvasGapKey(c, 0.950, 0.990, 0.020);
        c.SaveAs(
            (outdir + "/muon_q_eta_spectrum_pp_vs_pbpb_shape" + wp_tag + ".png").c_str());
    }

    // ---------------- numeric summary ----------------
    printf("\n================ gap definitions ================\n");
    if (use_new_gap_cuts) {
        printf("SELECTED (plotted): proposed pT-INDEPENDENT fiducial windows in q*eta\n");
        for (const auto& w : g_new_windows)
            printf("      [%+.4f, %+.4f]\n", w.first, w.second);
        printf("   (read from ParamsSet::single_mu_fiducial_gap_cuts)\n");
        printf("NOT selected, shown for reference --\n");
    } else {
        printf("SELECTED (plotted): legacy PassSingleMuonGapCut "
               "(RDFBasedHistFillingData.cxx:141)\n");
    }
    printf("PassSingleMuonGapCut (RDFBasedHistFillingData.cxx:141):\n");
    printf("  applied at ALL pT   -- reject |eta| < %.4f (ParamsSet::eta_gap_cut1)\n",
           ParamsSet::eta_gap_cut1);
    printf("  applied ONLY for pT < %.1f GeV -- reject q*eta in "
           "(ParamsSet::charge_eta_gap_cuts):\n",
           kGapCutPtThreshold);
    for (const auto& w : GapWindowsLowPtOnly())
        printf("      [%+.4f, %+.4f]\n", w.first, w.second);

    ReportDips(pp_all, "pp24");
    ReportDips(pb_all, "PbPb 23+24+25, 0-80%");

    // Fraction of single muons the EXISTING PassSingleMuonGapCut would remove.
    // Counted muon-by-muon in the fill loop over exactly the plotted population,
    // so the pT < 6 GeV condition and the 0-80% centrality gate are both exact.
    printf("\nFraction of single muons FAILING the SELECTED gap cut (%s):\n",
           use_new_gap_cuts ? "proposed fiducial" : "legacy PassSingleMuonGapCut");
    printf("   pp24                : %.4f  (%lld / %lld)\n",
           pp_fc.n_wp > 0 ? double(pp_fc.n_gapcut) / pp_fc.n_wp : 0., pp_fc.n_gapcut,
           pp_fc.n_wp);
    printf("   PbPb 23+24+25 0-80%% : %.4f  (%lld / %lld)\n",
           pb_fc.n_wp > 0 ? double(pb_fc.n_gapcut) / pb_fc.n_wp : 0., pb_fc.n_gapcut,
           pb_fc.n_wp);
    for (int k = 0; k < nctr; ++k)
        printf("   PbPb %2d-%2d%%         : %.4f  (%lld / %lld)\n", ctr_lo[k], ctr_hi[k],
               pb_fc.n_wp_ctr[k] > 0
                   ? double(pb_fc.n_gapcut_ctr[k]) / pb_fc.n_wp_ctr[k]
                   : 0.,
               pb_fc.n_gapcut_ctr[k], pb_fc.n_wp_ctr[k]);

    printf("\nPNGs written to %s\n", outdir.c_str());
}
