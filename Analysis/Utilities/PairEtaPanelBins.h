#ifndef PAIR_ETA_PANEL_BINS_H
#define PAIR_ETA_PANEL_BINS_H

#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

#include "TAxis.h"

#include "../MuonObjectsParamsAndHelpers/ParamsSet.h"
#include "../RDFBasedHistFilling/CommonEffcyConfig.h"

// =================================================================================================
// PANEL/AXIS ALIGNMENT GUARD for the coarse pair-eta panels
// (CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap).
//
// WHY IT EXISTS. Every coarse pair-eta panel is produced by projecting a bin RANGE of a
// histogram's eta axis. That is only the panel it claims to be if BOTH panel edges are bin edges
// of THAT histogram's axis. When they are not, the projection silently integrates a different
// range than the label says -- and at the outermost panel it silently DROPS every pair beyond the
// last included bin, with no warning and no visible artefact.
//
// This became a live hazard on 2026-09-07, when the coarse panels moved their outer edges to
// +-2.2 to follow the new pair-level fiducial cut |eta^pair| < ParamsSet::pair_eta_fiducial_max:
//   * pp24 books its pair-eta axis on ParamsSet::N_PAIR_ETA_CROSSX_BINS = 48 over [-2.4, 2.4],
//     width exactly 0.1, on which 2.2 IS a bin edge -> unaffected.
//   * Pb+Pb (RDFBasedHistFillingPbPb.cxx) and the TRUTH signal-acceptance producers
//     (RDFBasedHistFillingPythiaTruth.cxx, RDFBasedHistFillingPowhegTruth.cxx) USED to book a
//     RETYPED `44, -2.4, 2.4` axis, width 0.10909..., on which 2.2 is NOT a bin edge: the outer
//     panels quietly stopped at +-2.29091, discarding pairs in |eta^pair| in [2.29091, 2.4]
//     that no cut in those samples removed, while still being labelled (2.0, 2.2).
//     ALL THREE were migrated onto N_PAIR_ETA_CROSSX_BINS on 2026-09-08
//     (mu_pt45_gap125_pairpt9_adoption.md D1/D7), so the CODE no longer produces such an axis.
//     This guard stays, and is now aimed at STALE FILES: any histogram written before that date
//     still carries the 44-bin axis, and reading one must fail rather than mislabel.
// So: fail loudly. The fix for those samples is to book the axis on N_PAIR_ETA_CROSSX_BINS AND to
// adopt the fiducial + pair-level cuts, rerunning their hist filling in the SAME step
// (docs/signal_selection_change_impact.md; docs/tracking/muon_gap_cuts_acceptance.md F18).
//
// This header exists because the guard has MANY consumers, spread across
// plotting_codes/{single_b_analysis,mc_data_compr,reco_effcy,trig_effcy} and the reco-eff writer,
// and a guard that only some of them call is not a guard. Deliberately no count is maintained
// here: the authoritative list is whatever `grep -rn PairEtaPanels::` returns. The hazard is not
// only `FindBin` -- code that equates histogram BIN i with PANEL i (the MC-trig-eff closure
// plotters) is the same association spelled without FindBin, and needs the same check.
// =================================================================================================
namespace PairEtaPanels {

// The panels must tile exactly the region the pair-level fiducial cut leaves. This is the
// pair-eta analogue of the q*eta coupling that RDFBasedHistFillingData::SetIOPaths enforces with a
// startup throw (coarse q*eta top edge == forward gap window). Checked once per process.
inline void CheckPanelsMatchFiducialCut()
{
    static bool checked = false;
    if (checked) return;
    checked = true;

    const CommonEffcyConfig cfg;
    const auto& panels = cfg.pair_eta_proj_ranges_coarse_incl_gap;
    if (panels.empty()) throw std::runtime_error("PairEtaPanels: coarse pair-eta panels are empty");

    const float lo  = panels.front().first;
    const float hi  = panels.back().second;
    const float max = ParamsSet::pair_eta_fiducial_max;
    if (std::fabs(hi - max) > 1e-4f || std::fabs(lo + max) > 1e-4f)
        throw std::runtime_error(
            "PairEtaPanels: CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap spans ["
            + std::to_string(lo) + ", " + std::to_string(hi) + "] but the pair-level fiducial cut "
            "keeps |eta^pair| < ParamsSet::pair_eta_fiducial_max = " + std::to_string(max)
            + ". The panels must tile exactly the surviving region: if they are WIDER the outer "
              "panels are diluted by an unpopulated slice while their labels claim it, and if they "
              "are NARROWER they silently drop surviving pairs. Move one to match the other.");
}

inline void CheckAxisAligned(const TAxis* ax, const std::pair<float, float>& panel,
                             int b1, int b2, const std::string& context)
{
    const double tol = 1e-4;
    const double lo_edge = ax->GetBinLowEdge(b1);
    const double hi_edge = ax->GetBinUpEdge(b2);
    if (std::fabs(lo_edge - panel.first) > tol || std::fabs(hi_edge - panel.second) > tol)
        throw std::runtime_error(
            context + ": coarse pair-eta panel (" + std::to_string(panel.first) + ", "
            + std::to_string(panel.second) + ") is NOT bin-aligned with this histogram's eta axis "
            "-- the projection would actually cover (" + std::to_string(lo_edge) + ", "
            + std::to_string(hi_edge) + "), so the panel would be mislabelled and, at the outermost "
            "panel, would silently lose yield. This histogram's eta axis has "
            + std::to_string(ax->GetNbins()) + " bins over [" + std::to_string(ax->GetXmin())
            + ", " + std::to_string(ax->GetXmax()) + "]. Fix depends on which axis this is: "
            "if it is a FINE crossx pair-eta axis, book it on ParamsSet::N_PAIR_ETA_CROSSX_BINS "
            "(48 over [-2.4, 2.4], width 0.1, on which every coarse panel edge is a bin edge) and "
            "rerun the hist filling that produced this file, in the same step as adopting the "
            "fiducial + pair-level gap cuts; if it IS the 9-bin panel axis itself (booked by "
            "RangesToEdges from these same ranges, e.g. the MC-trig-eff closure files), the file "
            "simply predates the panel edges that are live now -- rerun its producer. See "
            "docs/tracking/muon_gap_cuts_acceptance.md F18.");
}

// Panel bin range + both checks in one call, so no call site can forget them.
inline std::pair<int, int> Bins(const TAxis* ax, const std::pair<float, float>& panel,
                                const std::string& context)
{
    CheckPanelsMatchFiducialCut();
    const int b1 = ax->FindBin(panel.first  + 1e-6);
    const int b2 = ax->FindBin(panel.second - 1e-6);
    CheckAxisAligned(ax, panel, b1, b2, context);
    return {b1, b2};
}

}  // namespace PairEtaPanels

#endif  // PAIR_ETA_PANEL_BINS_H
