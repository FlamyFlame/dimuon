// dr_correction_pt_groups.h
//
// THE PAIR-pT GROUPING of the DeltaR-correction cells -- the ONE place that knows a fit cell may
// span more than one filled pair-pT bin. Shared by fit_dr_corrections.cxx (which fits the cells),
// plot_dr_correction_fits.cxx (which re-projects the measured points into them) and
// dr_correction_apply.h (whose raw-bin fallback projects them again).
//
// WHY IT EXISTS. The plateau mode "nocorr_ptmerge" (user, 2026-08-17) fits the SAME histograms as
// "nocorr" but with the LAST TWO pair-pT bins merged into a single cell: on the 8-bin log axis the
// top two cells, p_T^pair in [72.1,104) and [104,150) GeV, run past where pp Pythia has yield --
// they hold the plateau-guard failures (mc_trigger_efficiency.md R24/R27) and the closure collapse
// (mc_trig_eff_closure.md R1b), and their dR fits are noise-dominated.
//
// IT IS NOT A NEW BINNING (.claude/CLAUDE.md 'Binnings'). The FILLED histograms keep
// `ParamsSet::pair_pt_coarse_bins` untouched; a merged cell is the two source bins PROJECTED
// TOGETHER, which is numerically identical to having filled a 7-bin axis (the projection sums
// num / denom / errA / errB before any ratio is formed, exactly as the fill would have). Two
// consequences that make this the safe construction:
//   * every group edge is READ from the source histogram's own axis, never retyped, so the
//     grouping cannot drift from the binning it groups;
//   * the variant is opt-in and suffixed (`_nocorr_ptmerge`), so it can neither overwrite nor be
//     mistaken for the un-merged trees.
//
// SCOPE GUARD: the merge is defined for the 8-bin NOMINAL axis only. The `_pt4bin` comparison
// variant already merges the top cells by construction, and merging its last two would leave 3
// cells over 8-150 GeV -- a third pair-pT view nobody asked for. Asking for it THROWS.

#ifndef DR_CORRECTION_PT_GROUPS_H
#define DR_CORRECTION_PT_GROUPS_H

#include <stdexcept>
#include <string>
#include <vector>

#include <TAxis.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>

#include "dr_correction_ratio.h"
#include "dr_correction_sample_cfg.h"
#include "../../../Utilities/MCTrigEffPairPtBinning.h"

// One entry per FIT CELL along pair pT: the source-histogram bin range it covers, plus the group
// edges. `n == source bins` and every range is a single bin unless the last two were merged.
struct DrPtGroups {
    int                 n = 0;        // number of groups = number of fit cells in pair pT
    std::vector<int>    lo, hi;       // 1-based source-bin range of each group, size n
    std::vector<double> edges;        // size n+1, read from the source axis
    bool                merged = false;
};

// Build the grouping from the axis of the histogram that is actually being fitted.
// `merge_last_two` comes from DrCorrModeMergeLastTwoPt(plateau_mode) -- never from a literal.
inline DrPtGroups MakeDrPtGroups(const TAxis* pt_axis, bool merge_last_two)
{
    if (!pt_axis) throw std::runtime_error("MakeDrPtGroups: null pair-pT axis");
    const int nsrc = pt_axis->GetNbins();

    if (merge_last_two) {
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
    }

    DrPtGroups G;
    G.merged = merge_last_two;
    for (int i = 1; i <= nsrc; ++i) {
        if (merge_last_two && i == nsrc) {            // fold the last bin into the previous group
            G.hi.back() = i;
            continue;
        }
        G.lo.push_back(i);
        G.hi.push_back(i);
    }
    G.n = (int)G.lo.size();
    G.edges.resize(G.n + 1);
    for (int g = 0; g < G.n; ++g) G.edges[g] = pt_axis->GetBinLowEdge(G.lo[g]);
    G.edges[G.n] = pt_axis->GetBinUpEdge(G.hi[G.n - 1]);
    return G;
}

// The (pair pT, pair eta) cell ratio for GROUP iy (1-based; 0 = integrate over pair pT, i.e. the
// inclusive cell). Identical to DrCellRatio for an un-merged grouping -- it IS DrCellRatio with the
// group's source-bin range, so the merged cell is built by summing the two bins' num/denom/errA/
// errB before the ratio, never by averaging two ratios.
inline TH1D* DrGroupCellRatio(TH3D* hn, TH3D* hd, TH3D* ha, TH3D* hb, TH3D* hp, TH3D* hq,
                              const DrPtGroups& G, int iy, int iz, const char* nm)
{
    if (iy < 0 || iy > G.n)
        throw std::runtime_error("DrGroupCellRatio: pair-pT group " + std::to_string(iy)
                                 + " out of range [0," + std::to_string(G.n) + "]");
    const int ylo = (iy == 0) ? 0 : G.lo[iy - 1];     // 0 -> DrCellRatioRange integrates the axis
    const int yhi = (iy == 0) ? 0 : G.hi[iy - 1];
    return DrCellRatioRange(hn, hd, ha, hb, hp, hq, ylo, yhi, iz, iz, nm);
}

// A (pair pT, pair eta) map on the GROUPED pair-pT axis, with the pair-eta axis copied from `like`.
// Used for every per-cell map the fit stage writes, so the maps and the fitted cells cannot
// describe different grids.
inline TH2D* BookDrGroupMap(const TH2D* like, const DrPtGroups& G, const std::string& name,
                            const std::string& ztitle)
{
    const TAxis* ay = like->GetYaxis();               // pair eta -- unchanged by the grouping
    const int ny = ay->GetNbins();
    std::vector<double> ye(ny + 1);
    for (int i = 0; i < ny; ++i) ye[i] = ay->GetBinLowEdge(i + 1);
    ye[ny] = ay->GetBinUpEdge(ny);
    auto* h = new TH2D(name.c_str(),
                       (";p_{T}^{pair} [GeV];#eta^{pair};" + ztitle).c_str(),
                       G.n, G.edges.data(), ny, ye.data());
    h->SetDirectory(nullptr);
    return h;
}

#endif  // DR_CORRECTION_PT_GROUPS_H
