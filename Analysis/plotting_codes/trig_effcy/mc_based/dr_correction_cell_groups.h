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
//   * PAIR-ETA MERGE ("nocorr_etamerge", user 2026-09-03; SUPERSEDES the 2026-08-24 signed
//     negative-endcap/barrel/positive-endcap grouping) -- the 9 pair-eta bins become THREE
//     SIGN-INDEPENDENT |eta^pair| BINS: |eta| < 1.0 (barrel), 1.0 <= |eta| < 2.0, and 2.0 up to
//     the source axis's top edge -- which is 2.2 since 2026-09-07, because that axis's outer
//     edges track ParamsSet::pair_eta_fiducial_max. The group boundaries {1.0, 2.0} are looked
//     up in the source axis rather than retyped, so the fold follows the axis automatically.
//     The dR correlation was found to barely depend on the SIGN of pair eta,
//     while the >2.0-vs-<2.0 split inside the endcap is a much bigger effect than any
//     negative/positive asymmetry -- so the fold that buys back statistics is now in |eta|, not in
//     the detector-region sign. Same trade-off as the pair-pT merge: the 9-bin pair-eta grid is
//     the CROSS-SECTION's presentation binning, never chosen for eps_dR's statistics, and folding
//     it triples the pairs per fit while keeping the one distinction that is physically motivated.
//     The barrel group is ONE contiguous run of source bins; the two forward groups FOLD their
//     negative- and positive-eta source bins TOGETHER -- the one grouping on this axis a single
//     contiguous range cannot describe (see DrAxisGroups below).
//   * "nocorr_etamerge_ptmerge" applies both.
//
// NEITHER IS A NEW BINNING (.claude/CLAUDE.md 'Binnings'). The FILLED histograms keep
// `ParamsSet::pair_pt_coarse_bins` and `CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap`
// untouched; a merged cell is its source bins PROJECTED TOGETHER, which is numerically identical
// to having filled a coarser/folded axis (the projection sums num / denom / errA / errB over every
// source sub-range before any ratio is formed, exactly as the fill would have). Two consequences
// that make this the safe construction:
//   * every group edge is READ from the source histogram's own axis, never retyped, so the
//     grouping cannot drift from the binning it groups. For pair eta, where the grouping cannot be
//     expressed as an index rule, the requested |eta| boundaries are LOOKED UP as a SYMMETRIC PAIR
//     of existing edges (one on each side of 0), and a boundary that fails that THROWS -- it never
//     rebins or invents an edge;
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
#include <utility>
#include <vector>

#include <TAxis.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TString.h>

#include "dr_correction_ratio.h"
#include "dr_correction_sample_cfg.h"
#include "../../../Utilities/MCTrigEffPairPtBinning.h"

// One entry per FIT CELL along one axis: the source-histogram sub-range(s) it covers, plus the
// group edges. `n == source bins` (identity) or the number of merged cells. A group is USUALLY one
// contiguous 1-based inclusive [lo,hi] sub-range of the source axis; the sign-independent |eta|
// fold is the one grouping that needs MORE than one -- a forward |eta| group is the UNION of a
// negative-eta and a positive-eta source sub-range, which a single [lo,hi] cannot express.
struct DrAxisGroups {
    int                 n = 0;        // number of groups = number of fit cells on this axis
    std::vector<std::vector<std::pair<int,int>>> ranges;  // size n; each group's source sub-range(s)
    std::vector<double> edges;        // size n+1, the OUTPUT axis (read from / derived off the source)
    bool                merged = false;
    bool                folded = false;   // true if any group folds >1 disjoint source sub-range
};

// The historical name, kept so the pair-pT call sites read the same as before.
using DrPtGroups = DrAxisGroups;

// The IDENTITY grouping of an axis: one group per filled bin. Every mode uses it on the axes it
// does not group, so the call sites never branch on "is this axis grouped".
inline DrAxisGroups MakeDrIdentityGroups(const TAxis* ax)
{
    if (!ax) throw std::runtime_error("MakeDrIdentityGroups: null axis");
    DrAxisGroups G;
    const int nb = ax->GetNbins();
    G.ranges.resize(nb);
    for (int i = 1; i <= nb; ++i) G.ranges[i - 1] = {{i, i}};
    G.n = nb;
    G.edges.resize(G.n + 1);
    for (int g = 0; g < G.n; ++g) G.edges[g] = ax->GetBinLowEdge(g + 1);
    G.edges[G.n] = ax->GetBinUpEdge(G.n);
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
        if (i == nsrc) { G.ranges.back()[0].second = i; continue; }   // fold last bin into previous
        G.ranges.push_back({{i, i}});
    }
    G.n = (int)G.ranges.size();
    G.edges.resize(G.n + 1);
    for (int g = 0; g < G.n; ++g) G.edges[g] = pt_axis->GetBinLowEdge(G.ranges[g][0].first);
    G.edges[G.n] = pt_axis->GetBinUpEdge(G.ranges[G.n - 1][0].second);
    return G;
}

// THE PAIR-ETA GROUP BOUNDARIES (user, 2026-09-03; SUPERSEDES the 2026-08-24 signed
// negative-endcap/barrel/positive-endcap grouping -- see the header comment and
// docs/tracking/mc_trigeff_dr_binning_approaches.md). Sign-independent |eta^pair| INTERIOR
// boundaries; the outer one is read from the source axis. Each boundary MUST be an existing edge
// on BOTH sides of 0 (the source axis must be symmetric about 0 for a sign-independent fold to be
// well defined); MakeDrEtaGroups throws otherwise, never inventing an edge.
inline const std::vector<double>& DrEtaAbsMergeInteriorBoundaries()
{
    static const std::vector<double> b = {1.0, 2.0};
    return b;
}

// Build the pair-eta grouping. `merge_eta` comes from DrCorrModeMergeEta(plateau_mode). The
// barrel group (|eta| < the first boundary) straddles zero and is therefore ONE contiguous source
// range; every other group is the UNION of a negative-eta and a positive-eta source sub-range
// (e.g. |eta| in [1,2) = source bins covering [-2,-1) UNION source bins covering [1,2)).
inline DrAxisGroups MakeDrEtaGroups(const TAxis* eta_axis, bool merge_eta)
{
    if (!eta_axis) throw std::runtime_error("MakeDrEtaGroups: null pair-eta axis");
    if (!merge_eta) return MakeDrIdentityGroups(eta_axis);

    const int nsrc = eta_axis->GetNbins();

    // LOOK EDGES UP in the axis; never rebin to them. `FindBin` alone would happily return the
    // bin CONTAINING a boundary that is not an edge, which is exactly how a grouping silently
    // stops describing the binning it groups.
    auto find_bin_with_low_edge = [&](double val) -> int {
        for (int i = 1; i <= nsrc; ++i)
            if (std::fabs(eta_axis->GetBinLowEdge(i) - val) < 1e-6) return i;
        return -1;
    };
    auto find_bin_with_up_edge = [&](double val) -> int {
        for (int i = 1; i <= nsrc; ++i)
            if (std::fabs(eta_axis->GetBinUpEdge(i) - val) < 1e-6) return i;
        return -1;
    };

    const double eta_lo = eta_axis->GetBinLowEdge(1), eta_hi = eta_axis->GetBinUpEdge(nsrc);
    if (std::fabs(eta_lo + eta_hi) > 1e-6)
        throw std::runtime_error(
            "MakeDrEtaGroups: the source pair-eta axis [" + std::to_string(eta_lo) + ", "
            + std::to_string(eta_hi) + "] is not symmetric about 0 -- a sign-independent |eta| "
              "fold is not well defined for it");
    const double eta_max = eta_hi;

    std::vector<double> ab = {0.0};                                  // |eta| boundaries
    for (double b : DrEtaAbsMergeInteriorBoundaries()) ab.push_back(b);
    ab.push_back(eta_max);
    for (size_t i = 1; i < ab.size(); ++i)
        if (ab[i] <= ab[i - 1])
            throw std::runtime_error("MakeDrEtaGroups: |eta| group boundaries are not strictly "
                                     "increasing");

    DrAxisGroups G;
    G.merged = true;
    const int ngroup = (int)ab.size() - 1;
    G.ranges.resize(ngroup);
    G.edges.resize(ngroup + 1);
    for (int g = 0; g < ngroup; ++g) {
        const double lo = ab[g], hi = ab[g + 1];
        G.edges[g] = lo;
        if (g == 0) {
            // straddles zero -> one contiguous source range [-hi, +hi)
            const int blo = find_bin_with_low_edge(-hi);
            const int bhi = find_bin_with_up_edge(hi);
            if (blo < 0 || bhi < 0)
                throw std::runtime_error(
                    "MakeDrEtaGroups: the |eta| boundary " + std::to_string(hi) + " is not a "
                    "symmetric pair of existing edges on the source pair-eta axis. The merge "
                    "groups EXISTING bins; it must never invent an edge (.claude/CLAUDE.md "
                    "'Binnings').");
            G.ranges[g].push_back({blo, bhi});
        } else {
            const int neg_lo = find_bin_with_low_edge(-hi);
            const int neg_hi = find_bin_with_up_edge(-lo);
            const int pos_lo = find_bin_with_low_edge(lo);
            const int pos_hi = find_bin_with_up_edge(hi);
            if (neg_lo < 0 || neg_hi < 0 || pos_lo < 0 || pos_hi < 0)
                throw std::runtime_error(
                    "MakeDrEtaGroups: the |eta| boundary [" + std::to_string(lo) + ", "
                    + std::to_string(hi) + ") is not a symmetric pair of existing edges on the "
                      "source pair-eta axis. The merge groups EXISTING bins; it must never invent "
                      "an edge (.claude/CLAUDE.md 'Binnings').");
            G.ranges[g].push_back({neg_lo, neg_hi});
            G.ranges[g].push_back({pos_lo, pos_hi});
            G.folded = true;
        }
    }
    G.edges[ngroup] = ab.back();
    G.n = ngroup;

    // Every source bin must be covered EXACTLY ONCE -- the fold must be a partition of the source
    // axis, not an approximation of one.
    std::vector<int> cover(nsrc + 1, 0);
    for (const auto& rs : G.ranges)
        for (const auto& r : rs)
            for (int i = r.first; i <= r.second; ++i) ++cover[i];
    for (int i = 1; i <= nsrc; ++i)
        if (cover[i] != 1)
            throw std::runtime_error("MakeDrEtaGroups: internal inconsistency -- source bin "
                                     + std::to_string(i) + " covered " + std::to_string(cover[i])
                                     + " times (expected 1)");
    return G;
}

// One line for a log / report, so a run always states what it grouped. The OUTPUT axis (`edges`)
// is what is printed; for the |eta| fold that axis runs 0 -> eta_max, so the unit string should
// read "|eta^{pair}|" at call sites that pass a folded grouping (fit_dr_corrections.cxx).
inline std::string DrGroupsDescribe(const DrAxisGroups& G, const char* axis_name, const char* unit)
{
    if (!G.merged) return std::string(axis_name) + ": " + std::to_string(G.n)
                        + " cells (filled binning, no grouping)";
    std::string s = std::string(axis_name) + ": MERGED into " + std::to_string(G.n) + " cells";
    if (G.folded) s += " (sign-independent fold)";
    s += " --";
    for (int g = 0; g < G.n; ++g)
        s += Form("%s [%.4g, %.4g)%s", g ? "," : "", G.edges[g], G.edges[g + 1], unit);
    return s;
}

// The (pair pT, pair eta) cell ratio for pair-pT GROUP `iy` and pair-eta GROUP `iz` (both 1-based;
// iy = 0 integrates over pair pT and iz = 0 over pair eta, i.e. the inclusive cell). Identical to
// DrCellRatio for identity groupings -- it sums the group's source sub-range(s)' num/denom/errA/
// errB before the ratio, never averages ratios, and handles a group folded from more than one
// disjoint source sub-range (the sign-independent |eta| groups) via DrCellRatioMultiRange.
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
    const std::vector<std::pair<int,int>> yranges =
        (iy == 0) ? std::vector<std::pair<int,int>>{} : Gpt.ranges[iy - 1];
    const std::vector<std::pair<int,int>> zranges =
        (iz == 0) ? std::vector<std::pair<int,int>>{} : Geta.ranges[iz - 1];
    return DrCellRatioMultiRange(hn, hd, ha, hb, hp, hq, yranges, zranges, nm);
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
