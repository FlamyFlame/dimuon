// dr_correction_plateau.h
//
// THE ONE implementation of the per-cell large-dR PLATEAU of a dR-correction ratio: the weighted
// mean of eps_dR (or eps_single) over the plateau window, its statistical error, its weighted RMS
// scatter, and the plateau-WINDOW systematic |p[nominal] - p[retired [1,4]]|.
//
// It lives in a header for the same reason dr_correction_ratio.h does: this number NORMALIZES the
// nominal fit, so a second, drifted copy of the weighted mean would silently change every
// normalized curve and every chi2 built on it. Round 9 (mc_trigger_efficiency.md R17) is the
// standing example of what two coexisting definitions of one quantity cost.
//
// Producers/consumers:
//   plot_mc_trig_eff.cxx    measures it per (pair pT, pair eta) cell and writes the plateau ROOT
//                           file the fit stage reads (the .txt tables beside the plots are for
//                           humans; no code may parse them).
//   fit_dr_corrections.cxx  re-measures it on the MERGED cells when the pair-pT grouping merges
//                           the last two bins (plateau mode "nocorr_ptmerge"), because the plateau
//                           file on disk describes the UN-merged grid. Same code, same window,
//                           same estimator -- only the cell it is measured over differs.
//
// The window edges come from MCTrigEffPlateauWindow.h and are never retyped.

#ifndef DR_CORRECTION_PLATEAU_H
#define DR_CORRECTION_PLATEAU_H

#include <cmath>
#include <utility>
#include <vector>

#include <TH1D.h>

#include "../../../Utilities/MCTrigEffPlateauWindow.h"

// `syst` is < 0 when the alternative window has no usable dR bin, i.e. "not evaluable" -- the
// consumer must not silently read that as zero uncertainty.
struct PlateauCell { double mean, err, rms; int nb; double syst = -1.; };

// The per-cell large-dR plateau of a dR-correction ratio histogram, measured in the nominal window
// and in the retired [1,4] one (same histogram, so their difference is purely the window).
inline PlateauCell PlateauFromRatio(const TH1D* r)
{
    auto window = [&](double lo, double hi) -> PlateauCell {
        double sw = 0, swv = 0;
        std::vector<std::pair<double,double>> vw;  // (value, weight)
        for (int i = 1; i <= r->GetNbinsX(); ++i) {
            const double xc = r->GetBinCenter(i);
            if (xc < lo || xc > hi) continue;
            const double v = r->GetBinContent(i), e = r->GetBinError(i);
            if (e <= 0. || v == 0.) continue;
            const double w = 1. / (e * e); sw += w; swv += w * v; vw.push_back({v, w});
        }
        if (sw <= 0.) return PlateauCell{ -1, -1, -1, 0 };
        const double mean = swv / sw, err = std::sqrt(1. / sw);
        double swd = 0;
        for (auto& q : vw) swd += q.second * (q.first - mean) * (q.first - mean);
        const double rms = std::sqrt(swd / sw);  // weighted RMS scatter about the mean
        return PlateauCell{ mean, err, rms, (int)vw.size() };
    };
    PlateauCell nom = window(MCTrigEffPlateau::kLo, MCTrigEffPlateau::kHi);
    const PlateauCell alt = window(MCTrigEffPlateau::kSystLo, MCTrigEffPlateau::kSystHi);
    if (nom.nb > 0 && alt.nb > 0) nom.syst = std::fabs(nom.mean - alt.mean);
    return nom;
}

#endif  // DR_CORRECTION_PLATEAU_H
