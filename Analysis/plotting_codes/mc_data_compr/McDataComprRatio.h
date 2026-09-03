#pragma once

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include <TAxis.h>
#include <TH1.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMarker.h>
#include <TPad.h>
#include <TVirtualPad.h>

// =================================================================================================
// McDataComprRatio -- the MC/data ratio pad, and the two axis-cosmetics rules that come with it.
//
// FREE inline functions at namespace scope, deliberately: all THREE producers of this plot set
// need them, and `plot_mc_data_pair_pt_in_eta.cxx` does not derive from PlotMCDataComprBaseClass
// (see McDataComprColors.h for why that matters).
//
// ERROR PROPAGATION. The ratio is built with TH1::Divide(num, den, 1, 1) and NO "B" option, i.e.
// the binomial treatment is deliberately NOT used: data and MC are INDEPENDENT samples (a pp24
// data yield and a Pythia/POWHEG prediction), not a subset and its parent, so the correct
// combination is the quadrature sum of the two relative uncertainties. Using "B" here would
// understate the error badly wherever the two are of comparable size.
// =================================================================================================
namespace McDataComprRatio {

// Fraction of the cell height given to the ratio pad. One constant, because the text-size rescale
// of a RELATIVE-font pad (see RelativeTextScale) has to agree with it exactly.
constexpr double kRatioFrac = 0.30;

// Legend boxes must not sit on the top frame line. ROOT's auto-placement fills whatever gap it
// finds, including the topmost strip, so the guard is applied AFTER the pad has been painted.
constexpr double kLegendTopGuard = 0.88;

// -------------------------------------------------------------------------------------------
// Split one pad (typically a cell of TCanvas::Divide) into a main pad and a ratio pad that share
// the x axis. Returns through the two out-parameters; the caller cd()s into them.
// -------------------------------------------------------------------------------------------
inline void SplitPadForRatio(TVirtualPad* cell, TPad*& main, TPad*& ratio,
                             double left_margin = 0.23, double right_margin = 0.04)
{
    cell->cd();
    const std::string base = std::string(cell->GetName());

    main = new TPad((base + "_main").c_str(), "", 0., kRatioFrac, 1., 1.);
    main->SetLeftMargin(left_margin);
    main->SetRightMargin(right_margin);
    main->SetTopMargin(0.06);
    main->SetBottomMargin(0.02);   // the x axis is drawn by the ratio pad, not here
    main->Draw();

    ratio = new TPad((base + "_ratio").c_str(), "", 0., 0., 1., kRatioFrac);
    ratio->SetLeftMargin(left_margin);
    ratio->SetRightMargin(right_margin);
    ratio->SetTopMargin(0.03);
    ratio->SetBottomMargin(0.38);
    // NO SetGridy(): on a log ratio axis with sub-decade labels the grid degenerates into a dense
    // hatch that competes with the markers. The dashed unity line is the reference this pad needs.
    ratio->Draw();
}

// The factor by which a RELATIVE-size (font precision 2) text must be scaled in the ratio pad to
// come out the same physical size as in the main pad. Not needed for font 43 (precision 3, pixel
// sizes), which is what plotting_codes/helper_functions.c::hist_helper sets.
inline double RelativeTextScale(){ return (1. - kRatioFrac) / kRatioFrac; }

// -------------------------------------------------------------------------------------------
// The main pad must not repeat the x axis: no numeric labels, no title. (The frame line stays --
// the two pads share it.)
// -------------------------------------------------------------------------------------------
inline void HideXAxis(TH1* h)
{
    if (!h) return;
    h->GetXaxis()->SetLabelSize(0.);
    h->GetXaxis()->SetTitleSize(0.);
    h->GetXaxis()->SetTitle("");
}

// -------------------------------------------------------------------------------------------
// MC / data, with independent-sample errors. `num` and `den` must already carry whatever
// normalization and width scaling they are drawn with, so that the ratio is the ratio of the two
// curves the reader sees.
// -------------------------------------------------------------------------------------------
inline TH1D* MakeRatio(const TH1* num, const TH1* den, const std::string& name)
{
    if (!num || !den) return nullptr;
    TH1D* r = static_cast<TH1D*>(num->Clone(name.c_str()));
    r->SetDirectory(nullptr);
    r->Divide(num, den, 1., 1.);          // no "B": independent samples, see the file header
    r->SetTitle("");
    r->SetStats(0);
    r->SetLineColor(num->GetLineColor());
    r->SetMarkerColor(num->GetMarkerColor());
    r->SetMarkerStyle(num->GetMarkerStyle());
    r->SetMarkerSize(num->GetMarkerSize());
    return r;
}

// -------------------------------------------------------------------------------------------
// A y range that shows the ratio rather than a few outliers. Bins with zero content on either
// side come back as exactly 0 from TH1::Divide and are NOT data points -- they are excluded.
//
// `use_log` is set when the drawn ratios span more than ~1.5 decades, which happens in the
// generic family (the ΔR -> 0 truth divergence reaches MC/data ~ 1e3): on a linear axis that
// single point would flatten every other bin onto the unity line.
// -------------------------------------------------------------------------------------------
inline void AutoRatioRange(const std::vector<TH1*>& ratios, double& lo, double& hi, bool& use_log)
{
    // Two passes. The FIRST ignores bins whose ratio has a relative uncertainty above 50 %: a
    // point that is barely distinguishable from zero carries no information and must not set the
    // scale for the bins that do (in the 9-panel pair-pT view a single 54 %-error bin at 150 GeV
    // stretched the common range to 0.034 .. 6.56 and squeezed the whole 1.3-1.5 bulk into a
    // sliver). Those points are still DRAWN; they may simply fall outside the frame. The second
    // pass, over every positive bin, is the fallback when too few precise bins remain to define a
    // range at all.
    auto scan = [&ratios](double max_rel_err, double& rmin, double& rmax, int& n){
        rmin = 1e300; rmax = -1e300; n = 0;
        for (const TH1* r : ratios){
            if (!r) continue;
            for (int b = 1; b <= r->GetNbinsX(); ++b){
                const double v = r->GetBinContent(b);
                if (v <= 0. || !std::isfinite(v)) continue;
                if (max_rel_err > 0. && r->GetBinError(b) > max_rel_err * v) continue;
                rmin = std::min(rmin, v);
                rmax = std::max(rmax, v);
                ++n;
            }
        }
    };

    double rmin = 0., rmax = 0.; int n = 0;
    scan(0.5, rmin, rmax, n);
    if (n < 3) scan(-1., rmin, rmax, n);
    if (n == 0){ lo = 0.; hi = 2.; use_log = false; return; }

    use_log = (rmax / rmin > 30.);
    if (use_log){
        lo = rmin / 2.5;
        hi = rmax * 2.5;
    } else {
        // Always keep unity inside the frame: the dashed line at 1 is the reference the pad
        // exists for, and a range that excludes it hides the sign of the discrepancy.
        const double span = std::max(rmax, 1.) - std::min(rmin, 1.);
        lo = std::max(0., std::min(rmin, 1.) - 0.20 * span);
        hi = std::max(rmax, 1.) + 0.20 * span;
        if (hi - lo < 1e-6){ lo = 0.5; hi = 1.5; }
    }
}

// -------------------------------------------------------------------------------------------
// Log-y numeric labels. A ROOT log axis labels only the powers of ten by default, so a drawn
// range narrower than one decade can contain NO power of ten and end up with ZERO numeric labels
// (measured on generic/Deta_zoomin SS: 1303 .. 8569 = 0.82 decades, zero labels; Dphi_zoomin SS:
// 0.98 decades, one label). SetMoreLogLabels() adds the 2,3,...,9 sub-labels; SetNoExponent()
// prints them as plain numbers, which is only readable while the values themselves are of
// moderate magnitude. Applied to EVERY log-y frame in this plot set, not just the two panels
// where the defect was noticed.
// -------------------------------------------------------------------------------------------
inline void ApplyLogYLabelPolicy(TH1* frame, double lo, double hi)
{
    if (!frame || lo <= 0. || hi <= 0.) return;
    const double decades = std::log10(hi / lo);
    if (decades < 3.0) frame->GetYaxis()->SetMoreLogLabels();
    if (decades < 1.5 && hi < 1e5 && lo > 1e-3) frame->GetYaxis()->SetNoExponent();
}

// -------------------------------------------------------------------------------------------
// Ratio-pad frame. `rel_scale` = RelativeTextScale() for a relative-font histogram, or 1 for a
// pixel-font (precision 3) one, where the size is already pad-independent.
// -------------------------------------------------------------------------------------------
inline void StyleRatioFrame(TH1* r, const std::string& xtitle, double lo, double hi,
                            bool use_log, double rel_scale = 1.,
                            const std::string& ytitle = "MC / data")
{
    if (!r) return;
    TAxis* xa = r->GetXaxis();
    TAxis* ya = r->GetYaxis();
    xa->SetTitle(xtitle.c_str());
    ya->SetTitle(ytitle.c_str());
    ya->SetRangeUser(lo, hi);
    ya->SetNdivisions(505);
    if (rel_scale != 1.){
        xa->SetLabelSize(xa->GetLabelSize() * rel_scale);
        xa->SetTitleSize(xa->GetTitleSize() * rel_scale);
        ya->SetLabelSize(ya->GetLabelSize() * rel_scale);
        ya->SetTitleSize(ya->GetTitleSize() * rel_scale);
        ya->SetTitleOffset(ya->GetTitleOffset() / rel_scale);
        // Right-aligned by default, so at offset 1 the title lands on the last tick label
        // ("10^{2}" on the log pair-pT axis). Push it clear.
        xa->SetTitleOffset(1.45);
    }
    if (use_log){
        // The RATIO pad is only ~30 % of the cell height, so its label budget is a third of the
        // main pad's: MoreLogLabels over a multi-decade range packs the sub-decade labels tight
        // enough to run into the "MC / data" title (seen on generic/minv_zoomin SS, 2.5 decades).
        // Add them only where the range is narrow enough that the plain powers of ten would leave
        // the axis under-labelled -- the same defect ApplyLogYLabelPolicy exists to fix.
        const double decades = std::log10(hi / lo);
        if (decades < 1.5){
            ya->SetMoreLogLabels();
            if (hi < 1e4 && lo > 1e-3) ya->SetNoExponent();
        }
    }
}

// The dashed reference at MC/data = 1. Returned so the caller can keep it alive.
inline TLine* DrawUnityLine(const TH1* frame)
{
    if (!frame) return nullptr;
    const double x1 = frame->GetXaxis()->GetXmin();
    const double x2 = frame->GetXaxis()->GetXmax();
    TLine* l = new TLine(x1, 1., x2, 1.);
    l->SetLineStyle(2);
    l->SetLineColor(kGray + 2);
    l->SetLineWidth(1);
    l->Draw("same");
    return l;
}

// -------------------------------------------------------------------------------------------
// OUT-OF-RANGE ratio points.
//
// `AutoRatioRange` deliberately ignores points with a >50 % relative error when it picks the
// frame, so a genuine but badly-measured ratio can land outside it. ROOT then draws the error
// bar stub and clips the marker, and the reader cannot tell an off-scale point from a missing
// one -- found in review 2026-09-03 on the nine-panel signal figure, where two green POWHEG
// points at ratio ~4.1 rendered as bare stubs.
//
// This paints a hollow triangle AT the frame edge, pointing the way the point went. It is a
// position marker, not a value: it says "this point is off-scale in this direction", which is
// exactly the information the clipped drawing destroys. Call it AFTER the histogram is drawn.
// -------------------------------------------------------------------------------------------
inline std::vector<TMarker*> DrawOutOfRangeMarkers(const TH1* r, double lo, double hi)
{
    std::vector<TMarker*> marks;
    if (!r) return marks;
    const double span = hi - lo;
    if (!(span > 0.)) return marks;
    for (int b = 1; b <= r->GetNbinsX(); ++b) {
        const double y = r->GetBinContent(b);
        if (y == 0. && r->GetBinError(b) == 0.) continue;   // empty bin, not an off-scale one
        if (y <= hi && y >= lo) continue;
        const bool above = (y > hi);
        TMarker* m = new TMarker(r->GetXaxis()->GetBinCenter(b),
                                 above ? hi - 0.03 * span : lo + 0.03 * span,
                                 above ? 22 : 23);          // triangle up / down
        m->SetMarkerColor(r->GetMarkerColor());
        m->SetMarkerSize(1.0);
        m->Draw();
        marks.push_back(m);
    }
    return marks;
}

// -------------------------------------------------------------------------------------------
// Keep an auto-placed legend off the top frame line. ROOT resolves a default-constructed
// TLegend's position only when the pad is PAINTED, so the pad must be updated first; after that
// the coordinates are fixed and a shift sticks.
// -------------------------------------------------------------------------------------------
inline void FixLegendTopOverlap(TLegend* l, TVirtualPad* pad, double guard = kLegendTopGuard)
{
    if (!l || !pad) return;
    pad->Modified();
    pad->Update();
    const double y2 = l->GetY2NDC();
    if (y2 <= guard) return;
    const double shift = y2 - guard;
    const double y1 = l->GetY1NDC() - shift;
    if (y1 < 0.02) return;             // no room to move it; leave it where ROOT put it
    l->SetY1NDC(y1);
    l->SetY2NDC(y2 - shift);
    pad->Modified();
    pad->Update();
}

}  // namespace McDataComprRatio
