#ifndef DR_CORRECTION_CROSSX_EVALUATOR_H
#define DR_CORRECTION_CROSSX_EVALUATOR_H

#include <atomic>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <TAxis.h>

#include "MCTrigEffPairPtBinning.h"
#include "../MuonObjectsParamsAndHelpers/ParamsSet.h"
#include "../RDFBasedHistFilling/CommonEffcyConfig.h"
#include "../plotting_codes/trig_effcy/mc_based/dr_correction_apply.h"

// =================================================================================================
// DrCorrectionCrossxEvaluator -- the pair-trigger dR correlation correction eps_dR(dR) AS THE
// CROSS-SECTION APPLIES IT.
//
// WHAT IT ADDS OVER DrCorrectionEvaluator. That struct delivers ONE fit method per instance, with
// a per-cell fallback chain  fitted f(dR)/C  ->  raw measured bins / C_raw  ->  no correction.
// The pp24 cross-section is required (user, 2026-08-17) to use a CROSS-METHOD cascade instead:
//
//     per (pair pT, pair eta) cell:
//        1. the EXPONENTIAL fit           f = C + A exp[-(dR/lambda)^p]        if it is accepted
//        2. else the POLYNOMIAL fit       f = C + u^2 (a2 + a3 u + a4 u^2)     if IT is accepted
//        3. else the RAW (absolute) measured bins, normalised by their own dR in [0.5,1] mean
//
// "Accepted" is the producer's own verdict `h_step3_fit_ok == 1` (fit converged, enough
// informative points, f > 0 over the fit domain) AND DrCorrectionEvaluator's extra screen on the
// fitted baseline C (`DrCorrPlateauUsable`), because C is what the cell is divided by. Both
// screens already live inside DrCorrectionEvaluator: a cell it could not fit has `fit == nullptr`,
// which is exactly the signal this cascade branches on. Nothing is re-implemented here.
//
// Step 3 of the cascade is inherited, not added: DrCorrectionEvaluator's own raw-bin fallback IS
// the "use absolute values" branch, and it is flagged TEMPORARY there (mc_trigger_efficiency.md
// R26 is OPEN on whether a chi2 screen belongs in `fit_ok`). When R26 is decided, this cascade
// changes with it.
//
// CONFIGURATION is not a free parameter here: it is read from the named crossx defaults in
// dr_correction_sample_cfg.h -- DrCorrCrossxMethod() = "expo" (the primary),
// DrCorrCrossxSign() = "os" (opposite sign, matching the OS cross-section), DrCorrCrossxMode() =
// "nocorr_ptmerge" (no plateau correction; the canonical 8 pair-pT bins with the last two merged
// into one cell, [72.1, 150) GeV). The BACKUP method is the repo's other parametric form,
// "polyu_fixedRp". Nothing is retyped.
//
// BINNING GUARD. The correction cells and the cross-section's pair kinematics must be the same
// binning, or a pair is corrected by another cell's curve while every histogram still fills
// (.claude/CLAUDE.md §Binnings). The guard below compares the fit file's own axes, edge by edge,
// against ParamsSet::pair_pt_coarse_bins (merged where the mode merges) and
// CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap, and THROWS on any mismatch -- it never
// silently rebins.
//
// The correction is defined for dR < 1 only (the fit domain) and is exactly 1 above; that, the
// dR clamp into the stored TF1 range, and the floor/cap counters all come from
// DrCorrectionEvaluator and are deliberately not duplicated.
// =================================================================================================
// The BACKUP parametric form, used where the primary (DrCorrCrossxMethod()) fit is rejected.
// Named here rather than in dr_correction_sample_cfg.h because it is a property of THIS cascade,
// not of the fit production: the producer makes all three methods regardless.
inline const char* DrCorrCrossxBackupMethod() { return "polyu_fixedRp"; }

struct DrCorrectionCrossxEvaluator {

    DrCorrectionEvaluator primary;   // expo
    DrCorrectionEvaluator backup;    // polyu_fixedRp
    std::string mode, sign;

    // Per-cell routing decided ONCE at load time: 0 = primary fit, 1 = backup fit, 2 = raw/none.
    std::vector<std::vector<int>> route;
    // ATOMIC: Eval() runs inside RDF lambdas and RDataFrame is used with ImplicitMT enabled
    // (RDFBasedHistFillingBaseClass.cxx). A plain ++ here is a data race -- harmless to the
    // returned value (every lookup is a const read of a histogram or TF1) but it would corrupt
    // exactly the census that tells us how much of the correction is a placeholder.
    std::atomic<long long> n_eval{0}, n_primary{0}, n_backup{0}, n_raw{0}, n_outside{0};
    // Per-CELL census, decided once at load time. Distinct from the per-EVALUATION counters above:
    // a consumer that only holds a produced file must be able to state how much of the MAP is the
    // TEMPORARY raw-bin placeholder, which is a property of the cells, not of the pairs.
    int n_cells_primary = 0, n_cells_backup = 0, n_cells_raw = 0;

    // The pair-pT cell edges this correction is defined on, AFTER the mode's merging. Exposed so
    // the caller can state them in its log rather than re-deriving them.
    std::vector<double> pt_cell_edges, eta_cell_edges;

    void Load(const DrCorrSample& cfg, bool use_tight_wp)
    {
        mode = DrCorrCrossxMode();
        sign = DrCorrCrossxSign();

        primary.Load(cfg, use_tight_wp, DrCorrCrossxMethod(), sign, mode);
        backup .Load(cfg, use_tight_wp, DrCorrCrossxBackupMethod(), sign, mode);

        CheckSameGrid();
        CheckCanonicalBinning();

        const int npt  = primary.h_fit_ok->GetNbinsX();
        const int neta = primary.h_fit_ok->GetNbinsY();
        route.assign(npt, std::vector<int>(neta, 2));
        int n0 = 0, n1 = 0, n2 = 0;
        for (int iy = 0; iy < npt; ++iy)
            for (int iz = 0; iz < neta; ++iz) {
                if      (primary.cells[iy][iz].fit) { route[iy][iz] = 0; ++n0; }
                else if (backup .cells[iy][iz].fit) { route[iy][iz] = 1; ++n1; }
                else                                { route[iy][iz] = 2; ++n2; }
            }
        n_cells_primary = n0; n_cells_backup = n1; n_cells_raw = n2;
        std::cout << "DrCorrectionCrossxEvaluator [" << cfg.mc_label << " / "
                  << DrCorrSignText(sign) << " / " << mode << "]: " << npt << "x" << neta
                  << " cells routed -- " << n0 << " on the " << DrCorrCrossxMethod()
                  << " fit, " << n1 << " on the " << DrCorrCrossxBackupMethod()
                  << " fit (exponential rejected), " << n2
                  << " on the RAW measured bins (both fits rejected)" << std::endl;
        if (n2 > 0)
            std::cout << "  ** the raw-bin branch is a TEMPORARY placeholder "
                         "(mc_trigger_efficiency.md R26 is OPEN)" << std::endl;
        ReportDeliveredExtrema(npt, neta);
    }

    // eps_dR for one pair. Arguments in the SAME order as DrCorrectionEvaluator::Eval.
    double Eval(double dr, double pair_pt, double pair_eta)
    {
        ++n_eval;
        const int iy = primary.h_fit_ok->GetXaxis()->FindBin(pair_pt);
        const int iz = primary.h_fit_ok->GetYaxis()->FindBin(pair_eta);
        // OUTSIDE THE CELL GRID is its own outcome and must be counted as such. It is NOT the
        // raw-bin branch: the pair gets eps_dR = 1 (no correction) because no cell covers it --
        // typically pair pT below the 8 GeV bottom edge, which the signal region excludes but the
        // GENERIC sample does not. Lumping it in with "both fits rejected" made the census read
        // 53 % placeholder when the map itself has only ONE such cell out of 63.
        if (iy < 1 || iy > primary.h_fit_ok->GetNbinsX() ||
            iz < 1 || iz > primary.h_fit_ok->GetNbinsY()) {
            ++n_outside;
            return 1.0;
        }
        const int r = route[iy - 1][iz - 1];
        if (r == 1) { ++n_backup;  return backup .Eval(dr, pair_pt, pair_eta); }
        if (r == 0) { ++n_primary; return primary.Eval(dr, pair_pt, pair_eta); }
        ++n_raw;
        // Both fits rejected -> the raw branch, which lives inside either evaluator (they read the
        // same measured histograms and normalise the same way). Take it from the primary.
        return primary.Eval(dr, pair_pt, pair_eta);
    }

    void PrintStats()
    {
        const long long N = n_eval;
        auto pct = [&](long long n) { return N ? 100.0 * n / N : 0.0; };
        std::cout << "DrCorrectionCrossxEvaluator: " << N << " evaluations -- "
                  << n_primary << " routed to " << DrCorrCrossxMethod() << " (" << pct(n_primary)
                  << "%), " << n_backup << " to " << DrCorrCrossxBackupMethod() << " ("
                  << pct(n_backup) << "%), " << n_raw
                  << " to the RAW-BIN PLACEHOLDER (" << pct(n_raw)
                  << "%  <-- the cascade's per-PAIR placeholder share; the per-method lines below "
                     "count each sub-evaluator's OWN calls, which is a different question), "
                  << n_outside << " OUTSIDE the cell grid (" << pct(n_outside)
                  << "%, eps_dR = 1 -- pairs the correction does not cover, e.g. pair pT < 8 GeV)"
                  << std::endl;
        primary.PrintStats();
        if (n_backup > 0) backup.PrintStats();
    }

    // THE EXTREMUM OF WHAT THIS CASCADE ACTUALLY DELIVERS, over EVERY cell and every route.
    // DrCorrectionEvaluator prints its own span, but only over the cells IT fitted -- the raw-bin
    // branch bypasses that scan entirely, and the cascade can route a cell to a different
    // evaluator than the one that reported it. So neither component's printout bounds what the
    // cross-section is actually divided by.
    //
    // WHY THE BOUND MATTERS PHYSICALLY: eps_dR^2mu4 is the probability that a CLOSE-BY pair fires
    // 2mu4 relative to two independent legs. Close-by muons can only LOSE efficiency (shared RoIs,
    // the L1 two-muon separation), so eps_dR <= 1 is a physics requirement. A delivered value
    // above 1 does not make the pair more likely to fire; it means that cell's correction is an
    // artefact -- and since w_trig = 1/(e1 e2 eps_dR), it SUPPRESSES those pairs in the
    // cross-section. Nothing is clipped here (that is the user's call, coupled to the OPEN R26
    // in mc_trigger_efficiency.md) but it can never again be invisible.
    void ReportDeliveredExtrema(int npt, int neta)
    {
        // The scan below is a DIAGNOSTIC over cell CENTRES, not over pairs. Its Eval calls would
        // otherwise land in the sub-evaluators' counters and inflate the very census that says how
        // much of the delivered map is the TEMPORARY placeholder (measured: the scan added 99
        // raw-bin + 101 empty-bin evaluations to a real 71, i.e. the quoted placeholder share came
        // out 2.2x too high). Snapshot, scan, restore.
        const auto snap_p = primary.SnapshotCounters();
        const auto snap_b = backup.SnapshotCounters();
        double lo = 1e300, hi = -1e300;
        int lo_y = 0, lo_z = 0, hi_y = 0, hi_z = 0; double hi_dr = 0.;
        const TAxis* ax = primary.h_fit_ok->GetXaxis();
        const TAxis* ay = primary.h_fit_ok->GetYaxis();
        for (int iy = 1; iy <= npt; ++iy)
            for (int iz = 1; iz <= neta; ++iz) {
                const double pt  = ax->GetBinCenter(iy);
                const double eta = ay->GetBinCenter(iz);
                for (int k = 0; k < 200; ++k) {
                    const double dr = DrCorrectionEvaluator::kDrMax * k / 200.0;
                    const double v  = EvalNoCount(dr, pt, eta);
                    if (v < lo) { lo = v; lo_y = iy; lo_z = iz; }
                    if (v > hi) { hi = v; hi_y = iy; hi_z = iz; hi_dr = dr; }
                }
            }
        std::cout << "     DELIVERED eps_dR over ALL " << npt * neta << " cells and all routes "
                     "spans [" << lo << " (pair pT bin " << lo_y << ", pair eta bin " << lo_z
                  << "), " << hi << " (pair pT bin " << hi_y << ", pair eta bin " << hi_z
                  << ", dR = " << hi_dr << ")]" << std::endl;
        primary.RestoreCounters(snap_p);
        backup .RestoreCounters(snap_b);
        if (hi > 1.0)
            std::cout << "  ** WARNING: a delivered eps_dR ABOVE 1. A 2mu4 close-by correction is "
                         "a LOSS, so eps_dR <= 1 is a physics requirement -- that cell is a fit / "
                         "raw-bin ARTEFACT, and w_trig = 1/(e1 e2 eps_dR) SUPPRESSES its pairs. "
                         "Nothing is clipped here; this is OPEN (mc_trigger_efficiency.md R26) and "
                         "must be stated wherever the cross-section is quoted." << std::endl;
    }

private:
    // Eval without touching the counters, for the load-time scan.
    double EvalNoCount(double dr, double pair_pt, double pair_eta)
    {
        const int iy = primary.h_fit_ok->GetXaxis()->FindBin(pair_pt);
        const int iz = primary.h_fit_ok->GetYaxis()->FindBin(pair_eta);
        if (iy < 1 || iy > primary.h_fit_ok->GetNbinsX() ||
            iz < 1 || iz > primary.h_fit_ok->GetNbinsY()) return 1.0;
        return (route[iy - 1][iz - 1] == 1) ? backup.Eval(dr, pair_pt, pair_eta)
                                            : primary.Eval(dr, pair_pt, pair_eta);
    }

    void CheckSameGrid() const
    {
        if (primary.h_fit_ok->GetNbinsX() != backup.h_fit_ok->GetNbinsX() ||
            primary.h_fit_ok->GetNbinsY() != backup.h_fit_ok->GetNbinsY())
            throw std::runtime_error(std::string("DrCorrectionCrossxEvaluator: the ")
                + DrCorrCrossxMethod() + " and " + DrCorrCrossxBackupMethod()
                + " fit files have different cell grids -- one of them is stale");
    }

    // The cell grid the fits were produced on MUST be the canonical binning, with the last two
    // pair-pT bins merged iff the mode says so. Read from ParamsSet / CommonEffcyConfig; never
    // retyped.
    void CheckCanonicalBinning()
    {
        static const ParamsSet pms{};
        static const CommonEffcyConfig cfg{};

        // ParamsSet::pair_pt_coarse_bins DIRECTLY, not MCTrigEffPairPt::Edges(pms): the latter
        // returns the 4-bin comparison variant when MCTRIGEFF_PAIRPT_4BIN is set in the
        // environment, and a trigger-efficiency STUDY's environment variable must never be able to
        // move the cross-section's binning guard.
        std::vector<double> pt = pms.pair_pt_coarse_bins;
        if (DrCorrModeMergeLastTwoPt(mode)) {
            if (pt.size() < 3)
                throw std::runtime_error("DrCorrectionCrossxEvaluator: cannot merge the last two "
                                         "pair-pT bins of a binning with fewer than 2 bins");
            pt.erase(pt.end() - 2);          // drop the edge BETWEEN the last two bins
        }
        std::vector<double> eta;
        eta.push_back(cfg.pair_eta_proj_ranges_coarse_incl_gap.front().first);
        for (const auto& r : cfg.pair_eta_proj_ranges_coarse_incl_gap) eta.push_back(r.second);

        auto check = [](const TAxis* ax, const std::vector<double>& edges, const char* what) {
            if (ax->GetNbins() != (int)edges.size() - 1)
                throw std::runtime_error(std::string("DrCorrectionCrossxEvaluator: the fit file "
                    "has ") + std::to_string(ax->GetNbins()) + " " + what + " cells but the "
                    "canonical binning has " + std::to_string(edges.size() - 1)
                    + " -- stale fit file, or the wrong plateau mode");
            for (size_t i = 0; i < edges.size(); ++i) {
                const double got = (i + 1 <= (size_t)ax->GetNbins())
                                 ? ax->GetBinLowEdge(i + 1) : ax->GetBinUpEdge(ax->GetNbins());
                if (std::fabs(got - edges[i]) > 1e-6)
                    throw std::runtime_error(std::string("DrCorrectionCrossxEvaluator: ") + what
                        + " edge " + std::to_string(i) + " is " + std::to_string(got)
                        + " in the fit file but " + std::to_string(edges[i])
                        + " canonically -- stale fit file");
            }
        };
        check(primary.h_fit_ok->GetXaxis(), pt,  "pair-pT");
        check(primary.h_fit_ok->GetYaxis(), eta, "pair-eta");
        // The BACKUP too: CheckSameGrid compares only bin COUNTS, so a stale backup file with the
        // same number of cells but different edges would otherwise correct a pair with another
        // cell's curve -- silently, since every histogram still fills (.claude/CLAUDE.md §Binnings).
        check(backup.h_fit_ok->GetXaxis(), pt,  "pair-pT (backup)");
        check(backup.h_fit_ok->GetYaxis(), eta, "pair-eta (backup)");
        pt_cell_edges  = pt;
        eta_cell_edges = eta;
    }
};

#endif  // DR_CORRECTION_CROSSX_EVALUATOR_H
