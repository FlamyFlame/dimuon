#ifndef DR_CORRECTION_CASCADE_EVALUATOR_H
#define DR_CORRECTION_CASCADE_EVALUATOR_H

#include <atomic>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <TAxis.h>
#include <TString.h>

#include "MCTrigEffPairPtBinning.h"
#include "../MuonObjectsParamsAndHelpers/ParamsSet.h"
#include "../RDFBasedHistFilling/CommonEffcyConfig.h"
#include "../plotting_codes/trig_effcy/mc_based/dr_correction_apply.h"
#include "../plotting_codes/trig_effcy/mc_based/dr_correction_cell_groups.h"

// =================================================================================================
// DrCorrectionCascadeEvaluator -- the pair-trigger dR correlation correction eps_dR(dR) delivered
// by a CROSS-METHOD CASCADE, for ANY sign series and ANY no-plateau-correction plateau mode.
//
// WHAT IT ADDS OVER DrCorrectionEvaluator. That struct delivers ONE fit method per instance, with
// a per-cell chain  curve f(dR)/C  ->  raw measured bins / C_raw  ->  no correction. This one puts
// SEVERAL methods in front of that chain, per (pair pT, pair eta) cell:
//
//     1. the EXPONENTIAL fit        f = C + A exp[-(dR/lambda)^p]      if it is accepted
//     2. else the POLYNOMIAL fit    f = C + u^2 [A + a3(u-1) + a4(u^2-1)]  if IT is accepted
//                                   (the same quartic as C + a2 u^2 + a3 u^3 + a4 u^4, with
//                                   a2 = A - a3 - a4; written in A = f(0)-C since 2026-09-08 so
//                                   the Step-3 requirement f(0) <= C is one fit limit, A <= 0)
//     3. else the INTERPOLATION     linear through the measured points below R_p, flat at the last
//                                   measured knot above it            if IT is accepted
//     4. else the raw measured bins of the FIRST tier, normalised by their own dR in [0.5,1] mean
//        (DrCorrectionEvaluator's own TEMPORARY placeholder), and finally eps_dR = 1.
//
// Every tier delivers the SAME form  eps_dR(dR) = f(dR)/C  for dR < 1 and exactly 1 above, with C
// the tier's own baseline (the fitted free asymptote for 1-2, the flat-branch value for 3). That
// is why the interpolation can sit in the chain at all rather than being a different kind of
// object -- see the header of dr_correction_apply.h.
//
// "Accepted" is the producer's own verdict `h_step3_fit_ok == 1` AND DrCorrectionEvaluator's extra
// screen on the baseline C (`DrCorrPlateauUsable`), because C is what the cell is divided by. Both
// screens already live inside DrCorrectionEvaluator: a cell it could not use has `HasCurve()`
// false, which is exactly the signal this cascade branches on. Nothing is re-implemented here.
//
// WHY IT EXISTS (user, 2026-08-24; docs/tracking/mc_trigeff_dr_binning_approaches.md §PP-3 / D9).
// The four cell-grouping approaches being compared -- no merge / pair-pT merge / pair-eta merge /
// both -- must all be corrected by the SAME delivered cascade, on their own mode, or the
// comparison measures the cascade rather than the grouping. `interp` replaces the raw-bin
// placeholder as tier 3: it is one of the three resolutions mc_trigger_efficiency.md R26 names.
//
// ⚠ MIRROR NOTICE -- Utilities/DrCorrectionCrossxEvaluator.h IMPLEMENTS THE SAME IDEA.
// That class is what the pp24 CROSS-SECTION currently applies: the crossx defaults
// (DrCorrCrossxMethod/Sign/Mode) with the tier list `expo -> polyu_fixedRp -> raw bins`, i.e. THIS
// class with `{"expo","polyu_fixedRp"}` and no interpolation tier. It is deliberately left
// UNTOUCHED (user, 2026-08-24) so this comparison cannot move the cross-section before an approach
// is chosen. The two must not coexist beyond that decision: Remaining Work 1 of
// mc_trigeff_dr_binning_approaches.md is to REPLACE DrCorrectionCrossxEvaluator with this class
// configured for the winning mode and tier list -- not to edit it into agreement.
//
// BINNING GUARD. The correction cells must be the canonical binning, GROUPED exactly as the mode
// says, or a pair is corrected by another cell's curve while every histogram still fills
// (.claude/CLAUDE.md §Binnings). CheckCanonicalBinning below compares every fit file's own axes,
// edge by edge, against ParamsSet::pair_pt_coarse_bins and
// CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap put through the SAME grouping helpers
// the fit stage used, and THROWS on any mismatch -- it never silently rebins.
// =================================================================================================
// The default tier list. Named here rather than in dr_correction_sample_cfg.h because it is a
// property of THIS cascade, not of the fit production: the producer makes all three methods
// regardless of which of them any consumer chains.
inline const std::vector<std::string>& DrCorrCascadeMethods()
{
    static const std::vector<std::string> m = {"expo", "polyu_fixedRp", "interp"};
    return m;
}

struct DrCorrectionCascadeEvaluator {

    static constexpr int kRawTier = -1;          // route value: no tier has a curve

    std::vector<std::unique_ptr<DrCorrectionEvaluator>> tiers;
    std::vector<std::string> methods;
    std::string mode, sign, label;

    // Per-cell routing decided ONCE at load time: the tier index, or kRawTier.
    std::vector<std::vector<int>> route;

    // ATOMIC: Eval() runs inside RDF lambdas and RDataFrame is used with ImplicitMT enabled
    // (RDFBasedHistFillingBaseClass.cxx). A plain ++ here is a data race -- harmless to the
    // returned value (every lookup is a const read of a histogram, a TF1 or a TGraph) but it would
    // corrupt exactly the census that tells us how much of the correction is a placeholder.
    std::atomic<long long> n_eval{0}, n_raw{0}, n_outside{0};
    std::vector<std::atomic<long long>> n_tier;
    // Per-CELL census, decided once at load time. Distinct from the per-EVALUATION counters above:
    // a consumer that only holds a produced file must be able to state how much of the MAP is the
    // TEMPORARY raw-bin placeholder, which is a property of the cells, not of the pairs.
    std::vector<int> n_cells_tier;
    int n_cells_raw = 0;

    // The cell edges this correction is defined on, AFTER the mode's grouping. Exposed so the
    // caller can state them in its log rather than re-deriving them.
    std::vector<double> pt_cell_edges, eta_cell_edges;

    void Load(const DrCorrSample& cfg, bool use_tight_wp, const std::string& sign_in,
              const std::string& plateau_mode,
              const std::vector<std::string>& method_list = DrCorrCascadeMethods())
    {
        if (method_list.empty())
            throw std::runtime_error("DrCorrectionCascadeEvaluator: empty tier list");
        mode    = plateau_mode;
        sign    = sign_in;
        label   = cfg.mc_label;
        methods = method_list;

        for (const auto& m : methods) {
            auto e = std::make_unique<DrCorrectionEvaluator>();
            e->Load(cfg, use_tight_wp, m, sign, mode);
            tiers.push_back(std::move(e));
        }
        CheckSameGrid();
        CheckCanonicalBinning();

        const int npt  = Grid()->GetNbinsX();
        const int neta = Grid()->GetNbinsY();
        route.assign(npt, std::vector<int>(neta, kRawTier));
        n_cells_tier.assign(tiers.size(), 0);
        n_cells_raw = 0;
        for (int iy = 0; iy < npt; ++iy)
            for (int iz = 0; iz < neta; ++iz) {
                int r = kRawTier;
                for (size_t t = 0; t < tiers.size(); ++t)
                    if (tiers[t]->cells[iy][iz].HasCurve()) { r = (int)t; break; }
                route[iy][iz] = r;
                if (r == kRawTier) ++n_cells_raw; else ++n_cells_tier[r];
            }
        // std::atomic is neither copyable nor movable, so the counter vector is SIZED here rather
        // than assigned -- a resize on a live vector would try to move its elements.
        std::vector<std::atomic<long long>> fresh(tiers.size());
        for (auto& a : fresh) a.store(0);
        n_tier.swap(fresh);

        std::cout << "DrCorrectionCascadeEvaluator [" << label << " / "
                  << (sign.empty() ? "sign-integrated" : DrCorrSignText(sign)) << " / " << mode
                  << "]: " << npt << "x" << neta << " cells routed --";
        for (size_t t = 0; t < tiers.size(); ++t)
            std::cout << (t ? "," : "") << " " << n_cells_tier[t] << " on " << methods[t];
        std::cout << ", " << n_cells_raw << " on the RAW measured bins (every tier rejected)"
                  << std::endl;
        if (n_cells_raw > 0)
            std::cout << "  ** the raw-bin branch is a TEMPORARY placeholder "
                         "(mc_trigger_efficiency.md R26 is OPEN)" << std::endl;
        ReportDeliveredExtrema(npt, neta);
    }

    // The cell grid (axes + acceptance flags) of the FIRST tier -- CheckSameGrid has already
    // established that every tier carries the same one.
    const TH2D* Grid() const { return tiers.front()->h_fit_ok.get(); }

    // eps_dR for one pair. Arguments in the SAME order as DrCorrectionEvaluator::Eval.
    double Eval(double dr, double pair_pt, double pair_eta)
    {
        ++n_eval;
        const TH2D* g = Grid();
        const int iy = g->GetXaxis()->FindBin(pair_pt);
        // FOLD-AWARE (fixed 2026-09-08). A folded |eta| mode's Y axis runs 0 -> eta_max, so the
        // signed lookup this line used to do put EVERY negative-eta pair in the underflow bin,
        // where it was counted "outside the cell grid" and returned eps_dR = 1: half the sample
        // silently uncorrected in every `*_etamerge*` closure (50.67 % outside, against 0 % for an
        // un-folded mode). Same bug as mc_trigeff_dr_binning_approaches.md D11, in a class that had
        // reimplemented the lookup; the value now comes from the one shared helper.
        const int iz = g->GetYaxis()->FindBin(
            DrCorrectionEvaluator::CellLookupEta(pair_eta, tiers.front()->eta_folded));
        // OUTSIDE THE CELL GRID is its own outcome and must be counted as such. It is NOT the
        // raw-bin branch: the pair gets eps_dR = 1 (no correction) because no cell covers it --
        // typically pair pT below the coarse grid's bottom edge (ParamsSet::pair_pt_coarse_bins
        // .front(), 9 GeV since 2026-09-08), which the signal region excludes but the
        // GENERIC sample does not. Lumping the two together made the census read "53 % placeholder"
        // for a map with one placeholder cell.
        if (iy < 1 || iy > g->GetNbinsX() || iz < 1 || iz > g->GetNbinsY()) {
            ++n_outside;
            return 1.0;
        }
        const int r = route[iy - 1][iz - 1];
        if (r == kRawTier) {
            // Every tier rejected -> the raw branch, which lives inside each evaluator (they read
            // the same measured histograms and normalise them the same way). Take the first.
            ++n_raw;
            return tiers.front()->Eval(dr, pair_pt, pair_eta);
        }
        ++n_tier[r];
        return tiers[r]->Eval(dr, pair_pt, pair_eta);
    }

    // One line naming the cascade for a provenance stamp or a log.
    std::string Describe() const
    {
        std::string s = "cascade ";
        for (size_t t = 0; t < methods.size(); ++t) s += (t ? " -> " : "") + methods[t];
        s += " -> raw bins; cells routed:";
        for (size_t t = 0; t < methods.size(); ++t)
            s += Form("%s %d on %s", t ? "," : "", n_cells_tier[t], methods[t].c_str());
        s += Form(", %d on the RAW-BIN PLACEHOLDER", n_cells_raw);
        return s;
    }

    void PrintStats()
    {
        const long long N = n_eval;
        auto pct = [&](long long n) { return N ? 100.0 * n / N : 0.0; };
        std::cout << "DrCorrectionCascadeEvaluator: " << N << " evaluations --";
        for (size_t t = 0; t < tiers.size(); ++t)
            std::cout << (t ? "," : "") << " " << n_tier[t].load() << " routed to " << methods[t]
                      << " (" << pct(n_tier[t].load()) << "%)";
        std::cout << ", " << n_raw << " to the RAW-BIN PLACEHOLDER (" << pct(n_raw)
                  << "%  <-- the cascade's per-PAIR placeholder share; the per-method lines below "
                     "count each sub-evaluator's OWN calls, which is a different question), "
                  << n_outside << " OUTSIDE the cell grid (" << pct(n_outside)
                  << "%, eps_dR = 1 -- pairs the correction does not cover, e.g. pair pT < "
                  << ParamsSet{}.pair_pt_coarse_bins.front() << " GeV)"
                  << std::endl;
        for (size_t t = 0; t < tiers.size(); ++t)
            if (t == 0 || n_tier[t].load() > 0) tiers[t]->PrintStats();
    }

    // THE EXTREMUM OF WHAT THIS CASCADE ACTUALLY DELIVERS, over EVERY cell and every route.
    // Each DrCorrectionEvaluator prints its own span, but only over the cells IT accepted -- the
    // raw-bin branch bypasses that scan entirely, and the cascade can route a cell to a different
    // tier than the one that reported it. So no component's printout bounds what is applied.
    //
    // WHY THE BOUND MATTERS PHYSICALLY: eps_dR^2mu4 is the probability that a CLOSE-BY pair fires
    // 2mu4 relative to two independent legs. Close-by muons can only LOSE efficiency (shared RoIs,
    // the L1 two-muon separation), so eps_dR <= 1 is a physics requirement. A delivered value above
    // 1 does not make the pair more likely to fire; it means that cell's correction is an artefact
    // -- and since w_trig = 1/(e1 e2 eps_dR), it SUPPRESSES those pairs. Nothing is clipped here
    // (that is the user's call, coupled to the OPEN R26) but it can never be invisible.
    void ReportDeliveredExtrema(int npt, int neta)
    {
        // The scan below is a DIAGNOSTIC over cell CENTRES, not over pairs. Its Eval calls would
        // otherwise land in the sub-evaluators' counters and inflate the very census that says how
        // much of the delivered map is the TEMPORARY placeholder. Snapshot, scan, restore.
        std::vector<DrCorrectionEvaluator::Counters> snap;
        for (auto& t : tiers) snap.push_back(t->SnapshotCounters());
        double lo = 1e300, hi = -1e300;
        int lo_y = 0, lo_z = 0, hi_y = 0, hi_z = 0; double hi_dr = 0.;
        const TAxis* ax = Grid()->GetXaxis();
        const TAxis* ay = Grid()->GetYaxis();
        for (int iy = 1; iy <= npt; ++iy)
            for (int iz = 1; iz <= neta; ++iz) {
                const double pt  = ax->GetBinCenter(iy);
                const double eta = ay->GetBinCenter(iz);
                const int r = route[iy - 1][iz - 1];
                DrCorrectionEvaluator& ev = (r == kRawTier) ? *tiers.front() : *tiers[r];
                for (int k = 0; k < 200; ++k) {
                    const double dr = DrCorrectionEvaluator::kDrMax * k / 200.0;
                    const double v  = ev.Eval(dr, pt, eta);
                    if (v < lo) { lo = v; lo_y = iy; lo_z = iz; }
                    if (v > hi) { hi = v; hi_y = iy; hi_z = iz; hi_dr = dr; }
                }
            }
        std::cout << "     DELIVERED eps_dR over ALL " << npt * neta << " cells and all routes "
                     "spans [" << lo << " (pair pT bin " << lo_y << ", pair eta bin " << lo_z
                  << "), " << hi << " (pair pT bin " << hi_y << ", pair eta bin " << hi_z
                  << ", dR = " << hi_dr << ")]" << std::endl;
        for (size_t t = 0; t < tiers.size(); ++t) tiers[t]->RestoreCounters(snap[t]);
        if (hi > 1.0)
            std::cout << "  ** WARNING: a delivered eps_dR ABOVE 1. A 2mu4 close-by correction is "
                         "a LOSS, so eps_dR <= 1 is a physics requirement -- that cell is a fit / "
                         "raw-bin ARTEFACT, and w_trig = 1/(e1 e2 eps_dR) SUPPRESSES its pairs. "
                         "Nothing is clipped here; this is OPEN (mc_trigger_efficiency.md R26)."
                      << std::endl;
    }

private:
    void CheckSameGrid() const
    {
        for (size_t t = 1; t < tiers.size(); ++t)
            if (tiers[t]->h_fit_ok->GetNbinsX() != tiers[0]->h_fit_ok->GetNbinsX() ||
                tiers[t]->h_fit_ok->GetNbinsY() != tiers[0]->h_fit_ok->GetNbinsY())
                throw std::runtime_error("DrCorrectionCascadeEvaluator: the " + methods[0]
                    + " and " + methods[t] + " fit files have different cell grids -- one is stale");
    }

    // The cell grid the fits were produced on MUST be the canonical binning put through the mode's
    // own grouping. Read from ParamsSet / CommonEffcyConfig and grouped by the SAME helpers the fit
    // stage used (dr_correction_cell_groups.h); never retyped.
    void CheckCanonicalBinning()
    {
        static const ParamsSet pms{};
        static const CommonEffcyConfig ecfg{};

        // ParamsSet::pair_pt_coarse_bins DIRECTLY when the merge is on, because the merge is
        // DEFINED on the 8-bin nominal axis; otherwise MCTrigEffPairPt::Edges, so the `_pt4bin`
        // comparison variant validates against its own axis rather than against the nominal one.
        std::vector<double> pt = DrCorrModeMergeLastTwoPt(mode) ? pms.pair_pt_coarse_bins
                                                                : MCTrigEffPairPt::Edges(pms);
        if (DrCorrModeMergeLastTwoPt(mode)) {
            if (pt.size() < 3)
                throw std::runtime_error("DrCorrectionCascadeEvaluator: cannot merge the last two "
                                         "pair-pT bins of a binning with fewer than 2 bins");
            pt.erase(pt.end() - 2);          // drop the edge BETWEEN the last two bins
        }
        std::vector<double> eta;
        eta.push_back(ecfg.pair_eta_proj_ranges_coarse_incl_gap.front().first);
        for (const auto& r : ecfg.pair_eta_proj_ranges_coarse_incl_gap) eta.push_back(r.second);
        if (DrCorrModeMergeEta(mode)) {
            // Grouped by the SAME helper the fit stage used, so the guard cannot describe a
            // different grouping from the one that was fitted. A temporary axis carries the
            // canonical edges into MakeDrEtaGroups, which then LOOKS UP its boundaries in them.
            TAxis tmp((int)eta.size() - 1, eta.data());
            eta = MakeDrEtaGroups(&tmp, true).edges;
        }

        auto check = [&](const TAxis* ax, const std::vector<double>& edges, const char* what) {
            if (ax->GetNbins() != (int)edges.size() - 1)
                throw std::runtime_error(std::string("DrCorrectionCascadeEvaluator: the fit file "
                    "has ") + std::to_string(ax->GetNbins()) + " " + what + " cells but plateau "
                    "mode '" + mode + "' expects " + std::to_string(edges.size() - 1)
                    + " -- stale fit file, or the wrong plateau mode");
            for (size_t i = 0; i < edges.size(); ++i) {
                const double got = (i + 1 <= (size_t)ax->GetNbins())
                                 ? ax->GetBinLowEdge(i + 1) : ax->GetBinUpEdge(ax->GetNbins());
                if (std::fabs(got - edges[i]) > 1e-6)
                    throw std::runtime_error(std::string("DrCorrectionCascadeEvaluator: ") + what
                        + " edge " + std::to_string(i) + " is " + std::to_string(got)
                        + " in the fit file but " + std::to_string(edges[i])
                        + " canonically -- stale fit file");
            }
        };
        // EVERY tier, not just the first: CheckSameGrid compares only bin COUNTS, so a stale tier
        // with the same number of cells but different edges would correct a pair with another
        // cell's curve -- silently, since every histogram still fills.
        for (size_t t = 0; t < tiers.size(); ++t) {
            check(tiers[t]->h_fit_ok->GetXaxis(), pt,  "pair-pT");
            check(tiers[t]->h_fit_ok->GetYaxis(), eta, "pair-eta");
        }
        pt_cell_edges  = pt;
        eta_cell_edges = eta;
    }
};

#endif  // DR_CORRECTION_CASCADE_EVALUATOR_H
