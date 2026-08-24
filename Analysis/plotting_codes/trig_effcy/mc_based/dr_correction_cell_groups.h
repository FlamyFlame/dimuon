// dr_correction_cell_groups.h   (was dr_correction_pt_groups.h until 2026-08-24)
//
// THE CELL GROUPING of the DeltaR-correction fit cells -- the ONE place that knows a fit cell may
// span more than one FILLED bin, on EITHER axis. Shared by fit_dr_corrections.cxx (which fits the
// cells), plot_dr_correction_fits.cxx (which re-projects the measured points into them) and
// dr_correction_apply.h (whose raw-bin fallback projects them again).
//
// WHY IT EXISTS. Two plateau-mode families group cells rather than fitting the filled binning:
//
//   * PAIR-pT MERGE ("nocorr_ptmerge", user 2026-08-17) -- the LAST TWO pair-pT bins become one
//     cell. On the 8-bin log axis the top two cells, p_T^pair in [72.1,104) and [104,150) GeV, run
//     past where pp Pythia has yield: they hold the plateau-guard failures
//     (mc_trigger_efficiency.md R24/R27) and the closure collapse (mc_trig_eff_closure.md R5), and
//     their dR fits are noise-dominated.
//   * PAIR-ETA MERGE ("nocorr_etamerge", user 2026-08-24) -- the 9 pair-eta bins become THREE
//     PHYSICAL DETECTOR REGIONS: negative-eta endcap (-2.4,-1.0), barrel (-1.0,1.0), positive-eta
//     endcap (1.0,2.4). Same trade-off on the second axis: the 9-bin pair-eta grid is the
//     CROSS-SECTION's presentation binning, never chosen for eps_dR's statistics, and grouping it
//     triples the pairs per fit while keeping the one distinction that is physically motivated
//     (endcap L1 geometry differs from the barrel's). The two ENDCAPS are deliberately NOT merged
//     with each other: the r16578 forward anomaly is NEGATIVE-eta only (parent R8/R10/R14).
//   * "nocorr_etamerge_ptmerge" applies both.
//
// NEITHER IS A NEW BINNING (.claude/CLAUDE.md 'Binnings'). The FILLED histograms keep
// `ParamsSet::pair_pt_coarse_bins` and `CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap`
// untouched; a merged cell is its source bins PROJECTED TOGETHER, which is numerically identical
// to having filled a coarser axis (the projection sums num / denom / errA / errB before any ratio
// is formed, exactly as the fill would have). Two consequences that make this the safe
// construction:
//   * every group edge is READ from the source histogram's own axis, never retyped, so the
//     grouping cannot drift from the binning it groups. For pair eta, where the grouping cannot be
//     expressed as an index rule, the requested INTERIOR boundaries are LOOKED UP in the axis and
//     a boundary that is not an existing edge THROWS -- it never rebins;
//   * every variant is opt-in and suffixed (`_nocorr_ptmerge`, `_nocorr_etamerge`,
//     `_nocorr_etamerge_ptmerge`), so it can neither overwrite nor be mistaken for another tree.
//
// SCOPE GUARD (pair-pT merge only): defined for the 8-bin NOMINAL axis. The `_pt4bin` comparison
// variant already merges the top cells by construction, and merging its last two would leave 3
// cells over 8-150 GeV -- a third pair-pT view nobody asked for. Asking for it THROWS.
// The pair-eta merge has no such restriction: it is orthogonal to the pair-pT axis.

#ifndef DR_CORRECTION_CELL_GROUPS_H
#define DR_CORRECTION_CELL_GROUPS_H

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include <TAxis.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TString.h>

#include "dr_correction_ratio.h"
#include "dr_correction_sample_cfg.h"
#include "../../../Utilities/MCTrigEffPairPtBinning.h"

// One entry per FIT CELL along one axis: the source-histogram bin range it covers, plus the group
// edges. `n == source bins` and every range is a single bin unless that axis is grouped.
struct DrAxisGroups {
    int                 n = 0;        // number of groups = number of fit cells on this axis
    std::vector<int>    lo, hi;       // 1-based source-bin range of each group, size n
    std::vector<double> edges;        // size n+1, read from the source axis
    bool                merged = false;
};

// The historical name, kept so the pair-pT call sites read the same as before.
using DrPtGroups = DrAxisGroups;

// The IDENTITY grouping of an axis: one group per filled bin. Every mode uses it on the axes it
// does not group, so the call sites never branch on "is this axis grouped".
inline DrAxisGroups MakeDrIdentityGroups(const TAxis* ax)
{
    if (!ax) throw std::runtime_error("MakeDrIdentityGroups: null axis");
    DrAxisGroups G;
    for (int i = 1; i <= ax->GetNbins(); ++i) { G.lo.push_back(i); G.hi.push_back(i); }
    G.n = (int)G.lo.size();
    G.edges.resize(G.n + 1);
    for (int g = 0; g < G.n; ++g) G.edges[g] = ax->GetBinLowEdge(G.lo[g]);
    G.edges[G.n] = ax->GetBinUpEdge(G.hi[G.n - 1]);
    return G;
}

// Build the pair-pT grouping from the axis of the histogram that is actually being fitted.
// `merge_last_two` comes from DrCorrModeMergeLastTwoPt(plateau_mode) -- never from a literal.
inline DrAxisGroups MakeDrPtGroups(const TAxis* pt_axis, bool merge_last_two)
{
    if (!pt_axis) throw std::runtime_error("MakeDrPtGroups: null pair-pT axis");
    const int nsrc = pt_axis->GetNbins();

    if (!merge_last_two) return MakeDrIdentityGroups(pt_axis);

    // The merge is a property of the 8-bin nominal axis (see the SCOPE GUARD above). Both
    // conditions are checked: the switch, and the axis itself -- the switch says what the job
    // asked for, the axis says what it actually got.
    if (MCTrigEffPairPt::UseFourBin())
        throw std::runtime_error(
            "MakeDrPtGroups: the merged-last-two-pair-pT-bins mode is defined for the 8-bin "
            "nominal axis ONLY, but MCTRIGEFF_PAIRPT_4BIN is set. The _pt4bin variant already "
            "merges the top cells by construction.");
    if (nsrc != ParamsSet::N_COARSE_PAIR_PT_BINS)
        throw std::runtime_error(
            "MakeDrPtGroups: the merged mode expects the "
            + std::to_string(ParamsSet::N_COARSE_PAIR_PT_BINS)
            + "-bin nominal pair-pT axis, got " + std::to_string(nsrc) + " bins");

    DrAxisGroups G;
    G.merged = true;
    for (int i = 1; i <= nsrc; ++i) {
        if (i == nsrc) { G.hi.back() = i; continue; }   // fold the last bin into the previous group
        G.lo.push_back(i);
        G.hi.push_back(i);
    }
    G.n = (int)G.lo.size();
    G.edges.resize(G.n + 1);
    for (int g = 0; g < G.n; ++g) G.edges[g] = pt_axis->GetBinLowEdge(G.lo[g]);
    G.edges[G.n] = pt_axis->GetBinUpEdge(G.hi[G.n - 1]);
    return G;
}

// THE PAIR-ETA GROUP BOUNDARIES (user, 2026-08-24). INTERIOR boundaries only -- the outer ones are
// whatever the source axis's outer edges are, so nothing about the acceptance is retyped here.
// The three regions they produce are the physical ones: negative-eta endcap | barrel | positive-eta
// endcap. They are NOT a binning: each must already be an edge of
// CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap, and MakeDrEtaGroups throws if it is not.
inline const std::vector<double>& DrEtaMergeInteriorBoundaries()
{
    static const std::vector<double> b = {-1.0, 1.0};
    return b;
}

// Build the pair-eta grouping. `merge_eta` comes from DrCorrModeMergeEta(plateau_mode).
inline DrAxisGroups MakeDrEtaGroups(const TAxis* eta_axis, bool merge_eta)
{
    if (!eta_axis) throw std::runtime_error("MakeDrEtaGroups: null pair-eta axis");
    if (!merge_eta) return MakeDrIdentityGroups(eta_axis);

    const int nsrc = eta_axis->GetNbins();
    const auto& bnd = DrEtaMergeInteriorBoundaries();

    // LOOK THE BOUNDARIES UP in the axis; never rebin to them. `FindBin` alone would happily
    // return the bin CONTAINING a boundary that is not an edge, which is exactly how a grouping
    // silently stops describing the binning it groups.
    std::vector<int> cut_bins;              // first source bin of each group after the first
    for (double b : bnd) {
        int found = -1;
        for (int i = 1; i <= nsrc; ++i)
            if (std::fabs(eta_axis->GetBinLowEdge(i) - b) < 1e-6) { found = i; break; }
        if (found < 0)
            throw std::runtime_error(
                "MakeDrEtaGroups: the requested pair-eta group boundary " + std::to_string(b)
                + " is not an edge of the filled pair-eta axis. The merge groups EXISTING bins; "
                  "it must never invent an edge (.claude/CLAUDE.md 'Binnings').");
        if (found == 1)
            throw std::runtime_error(
                "MakeDrEtaGroups: pair-eta group boundary " + std::to_string(b)
                + " coincides with the axis's lower edge, which would produce an empty group");
        cut_bins.push_back(found);
    }
    for (size_t i = 1; i < cut_bins.size(); ++i)
        if (cut_bins[i] <= cut_bins[i - 1])
            throw std::runtime_error("MakeDrEtaGroups: the pair-eta group boundaries are not "
                                     "strictly increasing along the axis");

    DrAxisGroups G;
    G.merged = true;
    size_t next_cut = 0;
    for (int i = 1; i <= nsrc; ++i) {
        const bool starts_group = (i == 1) ||
            (next_cut < cut_bins.size() && i == cut_bins[next_cut]);
        if (starts_group) {
            if (i != 1) ++next_cut;
            G.lo.push_back(i);
            G.hi.push_back(i);
        } else {
            G.hi.back() = i;
        }
    }
    G.n = (int)G.lo.size();
    if (G.n != (int)bnd.size() + 1)
        throw std::runtime_error("MakeDrEtaGroups: internal inconsistency -- "
                                 + std::to_string(bnd.size()) + " boundaries produced "
                                 + std::to_string(G.n) + " groups");
    G.edges.resize(G.n + 1);
    for (int g = 0; g < G.n; ++g) G.edges[g] = eta_axis->GetBinLowEdge(G.lo[g]);
    G.edges[G.n] = eta_axis->GetBinUpEdge(G.hi[G.n - 1]);
    return G;
}

// One line for a log / report, so a run always states what it grouped.
inline std::string DrGroupsDescribe(const DrAxisGroups& G, const char* axis_name, const char* unit)
{
    if (!G.merged) return std::string(axis_name) + ": " + std::to_string(G.n)
                        + " cells (filled binning, no grouping)";
    std::string s = std::string(axis_name) + ": MERGED into " + std::to_string(G.n) + " cells --";
    for (int g = 0; g < G.n; ++g)
        s += Form("%s [%.4g, %.4g)%s", g ? "," : "", G.edges[g], G.edges[g + 1], unit);
    return s;
}

// The (pair pT, pair eta) cell ratio for pair-pT GROUP `iy` and pair-eta GROUP `iz` (both 1-based;
// iy = 0 integrates over pair pT and iz = 0 over pair eta, i.e. the inclusive cell). Identical to
// DrCellRatio for identity groupings -- it IS DrCellRatio with the groups' source-bin ranges, so a
// merged cell is built by summing its source bins' num/denom/errA/errB before the ratio, never by
// averaging ratios.
inline TH1D* DrGroupCellRatio(TH3D* hn, TH3D* hd, TH3D* ha, TH3D* hb, TH3D* hp, TH3D* hq,
                              const DrAxisGroups& Gpt, const DrAxisGroups& Geta,
                              int iy, int iz, const char* nm)
{
    if (iy < 0 || iy > Gpt.n)
        throw std::runtime_error("DrGroupCellRatio: pair-pT group " + std::to_string(iy)
                                 + " out of range [0," + std::to_string(Gpt.n) + "]");
    if (iz < 0 || iz > Geta.n)
        throw std::runtime_error("DrGroupCellRatio: pair-eta group " + std::to_string(iz)
                                 + " out of range [0," + std::to_string(Geta.n) + "]");
    const int ylo = (iy == 0) ? 0 : Gpt.lo[iy - 1];    // 0 -> DrCellRatioRange integrates the axis
    const int yhi = (iy == 0) ? 0 : Gpt.hi[iy - 1];
    const int zlo = (iz == 0) ? 0 : Geta.lo[iz - 1];
    const int zhi = (iz == 0) ? 0 : Geta.hi[iz - 1];
    return DrCellRatioRange(hn, hd, ha, hb, hp, hq, ylo, yhi, zlo, zhi, nm);
}

// A (pair pT, pair eta) map on the GROUPED axes. Used for every per-cell map the fit stage writes,
// so the maps and the fitted cells cannot describe different grids.
inline TH2D* BookDrGroupMap(const DrAxisGroups& Gpt, const DrAxisGroups& Geta,
                            const std::string& name, const std::string& ztitle)
{
    auto* h = new TH2D(name.c_str(),
                       (";p_{T}^{pair} [GeV];#eta^{pair};" + ztitle).c_str(),
                       Gpt.n,  const_cast<double*>(Gpt.edges.data()),
                       Geta.n, const_cast<double*>(Geta.edges.data()));
    h->SetDirectory(nullptr);
    return h;
}

#endif  // DR_CORRECTION_CELL_GROUPS_H
