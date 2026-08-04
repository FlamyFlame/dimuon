#pragma once

#include <algorithm>
#include <limits>
#include <vector>

#include "TH1.h"

// =============================================================================
// Common log-y range for every subplot of one multi-panel PNG.
//
// SINGLE SOURCE for the y-range convention of the crossx pair-pT-in-eta subplot
// figures (pp and Pb+Pb). Two properties, both required (user, 2026-08-04):
//
//  1. NO DATA POINT MAY FALL OFF THE BOTTOM OF THE FRAME. The floor is derived
//     from the data itself — the smallest positive bin content over ALL panels —
//     never from a fixed ratio to the maximum. A fixed floor such as
//     `SetMinimum(ymax * 1e-5)` silently deletes the highest-pair-pT bins: the
//     pp24 crossx spectra span max/min ≈ 2.3e5 within a single pair-eta panel,
//     so the last one or two points (the ones that show how far in pair pT the
//     data actually reach) fall below such a floor and are never drawn.
//
//  2. ALL PANELS OF ONE PNG SHARE ONE SCALE, so panels can be compared by eye.
//     ROOT's per-pad autoscale gives each panel its own range, which makes two
//     adjacent panels with the same y-position mean different values.
//
// The floor is set from the smallest positive bin CONTENT, not from the lower
// end of the error bar. In the highest-pT bins the content comes from one or two
// pairs, so the Poisson error reaches down to ~0 — a value no log axis can ever
// contain. Requiring the error bar to fit would drive the axis to an arbitrary
// depth and compress the whole spectrum; the marker being visible (with its bar
// clipped at the frame) is the correct rendering of "we have a point here, with
// a bad uncertainty".
//
// The padding factors reproduce the look of ROOT's own log autoscale (~1/3 of a
// decade below the lowest point, x1.45 above the highest) so the lowest marker
// sits clear of the frame edge instead of on it.
// =============================================================================

static constexpr double kCommonLogYFloorPad = 1.0 / 3.0;
static constexpr double kCommonLogYCeilPad  = 1.45;

// Scan every histogram, then give them all the same [min, max]. `hists` must be
// EVERY histogram drawn in the PNG, across all panels and all overlaid curves —
// a histogram left out of the scan can still be drawn off-frame.
// `ceil_pad` raises the top of the frame where a multi-entry legend needs more
// clearance than the default headroom (it only ADDS empty space; it can never
// hide a point).
// No-op if nothing positive is present (leaves ROOT's default autoscale).
inline void ApplyCommonLogYRange(const std::vector<TH1*>& hists,
                                 double ceil_pad = kCommonLogYCeilPad)
{
    double gmin = std::numeric_limits<double>::max();
    double gmax = 0.0;
    for (const TH1* h : hists) {
        if (!h) continue;
        for (int b = 1; b <= h->GetNbinsX(); ++b) {
            const double v = h->GetBinContent(b);
            if (v <= 0.0) continue;   // empty bins carry no information on a log axis
            gmin = std::min(gmin, v);
            gmax = std::max(gmax, v);
        }
    }
    if (gmax <= 0.0) return;

    const double ymin = gmin * kCommonLogYFloorPad;
    const double ymax = gmax * ceil_pad;
    for (TH1* h : hists) {
        if (!h) continue;
        h->SetMinimum(ymin);
        h->SetMaximum(ymax);
    }
}
