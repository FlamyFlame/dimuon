#ifndef DR_CORRECTION_CROSSX_EVALUATOR_H
#define DR_CORRECTION_CROSSX_EVALUATOR_H

#include <algorithm>
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
// pp24 CROSS-SECTION APPLIES IT, in REGION A of the hybrid weight
// (docs/tracking/pp24_trig_eff_hybrid_application.md Physics Procedure §2, §3a-b):
//
//     eps_trig^pair = eps^nc_data(1) * eps^nc_data(2) * eps_dR(dR; cell)     pair pT in bins 1..N-2
//
// of the canonical coarse pair-pT axis (ParamsSet::pair_pt_coarse_bins, N = 8 -> [9, 74.24) GeV).
// The top two bins are REGION B, served by the single-value pair efficiency instead
// (Utilities/PairTrigEffCrossxEvaluator.h), and are NOT delivered by this struct: asking it for a
// pair there THROWS.
//
// WHAT IT ADDS OVER DrCorrectionEvaluator. That struct delivers ONE fit method per instance. The
// cross-section is required (user, 2026-09-17) to use a per-cell PRIMARY FORM plus ONE fallback:
//
//     per region-A (pair pT, |eta^pair|) cell:
//        primary  = the EXPONENTIAL   f = C + A exp[-(dR/lambda)^p]              (DrCorrCrossxMethod)
//                   EXCEPT the cells DrCorrCrossxPolyPtBins() x the LAST |eta| group, whose
//                   primary is the constrained POLYNOMIAL
//                   f = C + u^2 [A + a3(u-1) + a4(u^2-1)]                          (DrCorrCrossxPolyMethod)
//        fallback = the linear INTERPOLATION through the measured points          (DrCorrCrossxFallbackMethod)
//                   if the primary is rejected
//        nothing else: a cell with neither accepted THROWS at load time. The raw-bin placeholder
//                   of DrCorrectionEvaluator is UNREACHABLE from here.
//
// "Accepted" is the producer's own verdict `h_step3_fit_ok == 1` AND DrCorrectionEvaluator's two
// baseline screens (`DrCorrPlateauUsable`, `DrCorrBaselineConsistent`), because C is what the cell
// is divided by. All of that lives inside DrCorrectionEvaluator: a cell it could not deliver has
// `HasCurve() == false`, which is exactly the signal this cascade branches on.
//
// CONFIGURATION is not a free parameter here: it is read from the named crossx choices in
// dr_correction_sample_cfg.h -- DrCorrCrossxMode() = "nocorr_etamerge" (no plateau correction,
// the canonical 8 pair-pT bins x the 3 sign-independent |eta^pair| groups), DrCorrCrossxSign() =
// "os" (opposite sign, matching the OS cross-section; same-sign pairs are TEMPORARILY weighted
// with these same numbers -- doc D2), and the three method tokens. Nothing is retyped.
//
// BINNING GUARD. The correction cells and the cross-section's pair kinematics must be the same
// binning, or a pair is corrected by another cell's curve while every histogram still fills
// (.claude/CLAUDE.md §Binnings). The guard below compares each fit file's own axes, edge by edge,
// against ParamsSet::pair_pt_coarse_bins and the live MakeDrEtaGroups fold of
// CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap, and THROWS on any mismatch -- it never
// silently rebins.
//
// The correction is defined for dR < 1 only (the fit domain) and is exactly 1 above; that, the dR
// clamp into the stored TF1 range, and the floor/cap counters all come from DrCorrectionEvaluator
// and are deliberately not duplicated.
// =================================================================================================
struct DrCorrectionCrossxEvaluator {

    DrCorrectionEvaluator expo;      // DrCorrCrossxMethod()
    DrCorrectionEvaluator poly;      // DrCorrCrossxPolyMethod()
    DrCorrectionEvaluator interp;    // DrCorrCrossxFallbackMethod()
    std::string mode, sign;

    enum Route { kExpo = 0, kPoly = 1, kInterp = 2, kRegionB = -1 };
    // Per-cell routing decided ONCE at load time; kRegionB marks the top two pair-pT bins, which
    // this struct does not serve.
    std::vector<std::vector<int>> route;
    int n_pt_region_a = 0;           // = N - 2: the pair-pT bins this struct delivers

    // ATOMIC: Eval() runs inside RDF lambdas and RDataFrame is used with ImplicitMT enabled
    // (RDFBasedHistFillingBaseClass.cxx). A plain ++ here is a data race -- harmless to the
    // returned value (every lookup is a const read of a histogram or TF1) but it would corrupt
    // exactly the census that tells us which curve each pair was corrected by.
    std::atomic<long long> n_eval{0}, n_expo{0}, n_poly{0}, n_interp{0}, n_outside{0};
    // Per-CELL census, decided once at load time.
    int n_cells_expo = 0, n_cells_poly = 0, n_cells_interp = 0;   // as ROUTED
    int n_cells_poly_primary = 0;                                   // how many cells ASK for poly

    // The cell edges this correction is defined on, exposed so the caller can state them.
    std::vector<double> pt_cell_edges, eta_cell_edges;

    void Load(const DrCorrSample& cfg, bool use_tight_wp)
    {
        mode = DrCorrCrossxMode();
        sign = DrCorrCrossxSign();

        // Each loader prints its OWN census over EVERY cell of its fit file, including the two
        // region-B pair-pT bins this struct never serves (there the etamerge fit files do carry
        // rejected cells and DrCorrectionEvaluator builds its raw-bin placeholder for them). The
        // routing census printed AFTER them is the one that describes what the cross-section uses.
        expo  .Load(cfg, use_tight_wp, DrCorrCrossxMethod(),         sign, mode);
        poly  .Load(cfg, use_tight_wp, DrCorrCrossxPolyMethod(),     sign, mode);
        interp.Load(cfg, use_tight_wp, DrCorrCrossxFallbackMethod(), sign, mode);

        CheckSameGrid();
        CheckCanonicalBinning();

        const int npt  = expo.h_fit_ok->GetNbinsX();
        const int neta = expo.h_fit_ok->GetNbinsY();
        if (npt < 3)
            throw std::runtime_error("DrCorrectionCrossxEvaluator: the coarse pair-pT axis has "
                                     + std::to_string(npt) + " bins; the hybrid needs at least 3 "
                                     "(region A = all but the last two)");
        n_pt_region_a = npt - 2;
        const TAxis* ax = expo.h_fit_ok->GetXaxis();
        const TAxis* ay = expo.h_fit_ok->GetYaxis();

        // The poly-primary cells: named pair-pT bins x the LAST |eta| group. Validated against the
        // axis they index (a bin index past region A would silently name nothing).
        for (int b : DrCorrCrossxPolyPtBins())
            if (b < 1 || b > n_pt_region_a)
                throw std::runtime_error("DrCorrectionCrossxEvaluator: DrCorrCrossxPolyPtBins names "
                    "pair-pT bin " + std::to_string(b) + ", outside region A (1.."
                    + std::to_string(n_pt_region_a) + ")");
        auto is_poly_primary = [&](int iy, int iz) {
            const auto& pb = DrCorrCrossxPolyPtBins();
            return iz == neta && std::find(pb.begin(), pb.end(), iy) != pb.end();
        };

        route.assign(npt, std::vector<int>(neta, kRegionB));
        std::cout << "DrCorrectionCrossxEvaluator [" << cfg.mc_label << " / " << DrCorrSignText(sign)
                  << " / " << mode << "]: region A = pair-pT bins 1.." << n_pt_region_a << " = ["
                  << ax->GetBinLowEdge(1) << ", " << ax->GetBinUpEdge(n_pt_region_a)
                  << ") GeV x " << neta << " |eta^pair| groups; bins " << n_pt_region_a + 1 << ".."
                  << npt << " = [" << ax->GetBinLowEdge(n_pt_region_a + 1) << ", "
                  << ax->GetBinUpEdge(npt) << ") GeV are REGION B (single-value pair efficiency, "
                     "not served here)" << std::endl;
        for (int iy = 1; iy <= n_pt_region_a; ++iy)
            for (int iz = 1; iz <= neta; ++iz) {
                const bool want_poly = is_poly_primary(iy, iz);
                if (want_poly) ++n_cells_poly_primary;
                DrCorrectionEvaluator& primary = want_poly ? poly : expo;
                const char* primary_name = want_poly ? DrCorrCrossxPolyMethod() : DrCorrCrossxMethod();
                const std::string cell = Form("pair pT [%g,%g) x |eta^pair| [%g,%g)",
                                              ax->GetBinLowEdge(iy), ax->GetBinUpEdge(iy),
                                              ay->GetBinLowEdge(iz), ay->GetBinUpEdge(iz));
                if (primary.cells[iy - 1][iz - 1].HasCurve()) {
                    route[iy - 1][iz - 1] = want_poly ? kPoly : kExpo;
                    if (want_poly) { ++n_cells_poly; std::cout << "  cell " << cell << ": POLYNOMIAL primary (user-named cell)" << std::endl; }
                    else             ++n_cells_expo;
                } else if (interp.cells[iy - 1][iz - 1].HasCurve()) {
                    route[iy - 1][iz - 1] = kInterp;
                    ++n_cells_interp;
                    std::cout << "  ** cell " << cell << ": primary " << primary_name
                              << " REJECTED -> FALLBACK to " << DrCorrCrossxFallbackMethod() << std::endl;
                } else {
                    throw std::runtime_error("DrCorrectionCrossxEvaluator: region-A cell " + cell
                        + " has neither its primary (" + primary_name + ") nor the fallback ("
                        + DrCorrCrossxFallbackMethod() + ") accepted. There is deliberately NO "
                        "raw-bin tier (user, 2026-09-17): refit / inspect the cell instead of "
                        "delivering an unmeasured correction.");
                }
            }
        std::cout << "DrCorrectionCrossxEvaluator: " << n_pt_region_a * neta << " region-A cells routed -- "
                  << n_cells_expo << " on " << DrCorrCrossxMethod() << ", " << n_cells_poly << " on "
                  << DrCorrCrossxPolyMethod() << " (" << n_cells_poly_primary << " asked for it), "
                  << n_cells_interp << " on the " << DrCorrCrossxFallbackMethod() << " fallback"
                  << std::endl;
        ReportDeliveredExtrema(neta);
    }

    // eps_dR for one REGION-A pair. Arguments in the SAME order as DrCorrectionEvaluator::Eval.
    double Eval(double dr, double pair_pt, double pair_eta)
    {
        ++n_eval;
        const int iy = expo.h_fit_ok->GetXaxis()->FindBin(pair_pt);
        // FOLD-AWARE: the mode is a `*_etamerge*` one, so the Y axis is |eta^pair| (0 -> 2.2) and a
        // SIGNED lookup would send every negative-eta pair to the underflow bin -- half the sample
        // silently uncorrected (mc_trigeff_single_value_pair_eff.md R5). Shared helper, one definition.
        const int iz = expo.h_fit_ok->GetYaxis()->FindBin(
            DrCorrectionEvaluator::CellLookupEta(pair_eta, expo.eta_folded));
        if (iy < 1 || iy > expo.h_fit_ok->GetNbinsX() ||
            iz < 1 || iz > expo.h_fit_ok->GetNbinsY()) {
            // Outside the cell grid there is no correction; the caller (the hybrid) is expected
            // to have restricted region A to the grid, so this is COUNTED and must come out 0.
            ++n_outside;
            return 1.0;
        }
        const int r = route[iy - 1][iz - 1];
        if (r == kRegionB)
            throw std::runtime_error(Form("DrCorrectionCrossxEvaluator::Eval: pair pT = %.2f GeV is in "
                "region B (the single-value pair efficiency), which this struct does not serve. "
                "Route the pair through PairTrigEffCrossxEvaluator.", pair_pt));
        if (r == kPoly)   { ++n_poly;   return poly  .Eval(dr, pair_pt, pair_eta); }
        if (r == kInterp) { ++n_interp; return interp.Eval(dr, pair_pt, pair_eta); }
        ++n_expo;
        return expo.Eval(dr, pair_pt, pair_eta);
    }

    void PrintStats()
    {
        const long long N = n_eval;
        auto pct = [&](long long n) { return N ? 100.0 * n / N : 0.0; };
        std::cout << "DrCorrectionCrossxEvaluator (region A): " << N << " evaluations -- "
                  << n_expo << " on " << DrCorrCrossxMethod() << " (" << pct(n_expo) << "%), "
                  << n_poly << " on " << DrCorrCrossxPolyMethod() << " (" << pct(n_poly) << "%), "
                  << n_interp << " on the " << DrCorrCrossxFallbackMethod() << " fallback ("
                  << pct(n_interp) << "%), " << n_outside << " OUTSIDE the cell grid ("
                  << pct(n_outside) << "%, eps_dR = 1; must be 0 -- the hybrid restricts region A "
                     "to the grid)" << std::endl;
        expo.PrintStats();
        if (n_poly   > 0) poly  .PrintStats();
        if (n_interp > 0) interp.PrintStats();
    }

    // THE EXTREMUM OF WHAT THIS CASCADE ACTUALLY DELIVERS over every region-A cell and route.
    // Each DrCorrectionEvaluator prints its own span, but over ALL the cells it fitted -- including
    // region-B cells of the fit file that this struct never uses, and cells routed elsewhere. So
    // neither component's printout bounds what the cross-section is actually divided by.
    //
    // WHY THE BOUND MATTERS PHYSICALLY: eps_dR^2mu4 is the probability that a CLOSE-BY pair fires
    // 2mu4 relative to two independent legs. Close-by muons can only LOSE efficiency, so
    // eps_dR <= 1 is a physics requirement; a delivered value above 1 is an artefact that
    // SUPPRESSES those pairs in the cross-section. Nothing is clipped here (the bound on the
    // delivered f/C is OPEN, mc_trigger_efficiency.md R12/R35) but it can never be invisible.
    void ReportDeliveredExtrema(int neta)
    {
        // The scan is a DIAGNOSTIC over cell centres, not over pairs; keep it out of the census.
        const auto se = expo.SnapshotCounters();
        const auto sp = poly.SnapshotCounters();
        const auto si = interp.SnapshotCounters();
        double lo = 1e300, hi = -1e300;
        int lo_y = 0, lo_z = 0, hi_y = 0, hi_z = 0; double hi_dr = 0.;
        const TAxis* ax = expo.h_fit_ok->GetXaxis();
        const TAxis* ay = expo.h_fit_ok->GetYaxis();
        for (int iy = 1; iy <= n_pt_region_a; ++iy)
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
        std::cout << "     DELIVERED eps_dR over the " << n_pt_region_a * neta << " region-A cells and "
                     "all routes spans [" << lo << " (pair pT bin " << lo_y << ", |eta| group " << lo_z
                  << "), " << hi << " (pair pT bin " << hi_y << ", |eta| group " << hi_z
                  << ", dR = " << hi_dr << ")]" << std::endl;
        expo.RestoreCounters(se);
        poly.RestoreCounters(sp);
        interp.RestoreCounters(si);
        if (hi > 1.0)
            std::cout << "  ** WARNING: a delivered eps_dR ABOVE 1. A 2mu4 close-by correction is "
                         "a LOSS, so eps_dR <= 1 is a physics requirement -- that cell is a fit "
                         "ARTEFACT, and w_trig = 1/(e1 e2 eps_dR) SUPPRESSES its pairs. Nothing is "
                         "clipped here; this is OPEN (mc_trigger_efficiency.md R26/R35) and must be "
                         "stated wherever the cross-section is quoted." << std::endl;
    }

private:
    // Eval without touching the counters, for the load-time scan. Region-A cells only.
    double EvalNoCount(double dr, double pair_pt, double pair_eta)
    {
        const int iy = expo.h_fit_ok->GetXaxis()->FindBin(pair_pt);
        const int iz = expo.h_fit_ok->GetYaxis()->FindBin(
            DrCorrectionEvaluator::CellLookupEta(pair_eta, expo.eta_folded));
        if (iy < 1 || iy > n_pt_region_a || iz < 1 || iz > expo.h_fit_ok->GetNbinsY()) return 1.0;
        const int r = route[iy - 1][iz - 1];
        if (r == kPoly)   return poly  .Eval(dr, pair_pt, pair_eta);
        if (r == kInterp) return interp.Eval(dr, pair_pt, pair_eta);
        return expo.Eval(dr, pair_pt, pair_eta);
    }

    void CheckSameGrid() const
    {
        auto same = [&](const DrCorrectionEvaluator& a, const DrCorrectionEvaluator& b) {
            return a.h_fit_ok->GetNbinsX() == b.h_fit_ok->GetNbinsX() &&
                   a.h_fit_ok->GetNbinsY() == b.h_fit_ok->GetNbinsY();
        };
        if (!same(expo, poly) || !same(expo, interp))
            throw std::runtime_error(std::string("DrCorrectionCrossxEvaluator: the ")
                + DrCorrCrossxMethod() + ", " + DrCorrCrossxPolyMethod() + " and "
                + DrCorrCrossxFallbackMethod()
                + " fit files have different cell grids -- one of them is stale");
    }

    // The cell grid the fits were produced on MUST be the canonical binning. Read from ParamsSet /
    // CommonEffcyConfig; never retyped. Checked on ALL THREE files: CheckSameGrid compares only
    // bin COUNTS, so a stale file with the same number of cells but different edges would
    // otherwise correct a pair with another cell's curve -- silently, since every histogram still
    // fills (.claude/CLAUDE.md §Binnings).
    void CheckCanonicalBinning()
    {
        static const ParamsSet pms{};
        static const CommonEffcyConfig cfg{};

        // ParamsSet::pair_pt_coarse_bins DIRECTLY, not MCTrigEffPairPt::Edges(pms): the latter
        // returns the 4-bin comparison variant when MCTRIGEFF_PAIRPT_4BIN is set in the
        // environment, and a trigger-efficiency STUDY's environment variable must never be able to
        // move the cross-section's binning guard.
        std::vector<double> pt = pms.pair_pt_coarse_bins;
        if (DrCorrModeMergeLastTwoPt(mode))
            throw std::runtime_error("DrCorrectionCrossxEvaluator: the hybrid weight serves the last "
                                     "two pair-pT bins from the single-value pair efficiency, so "
                                     "DrCorrCrossxMode() must NOT be a pT-merged mode (got '" + mode + "')");
        std::vector<double> eta;
        eta.push_back(cfg.pair_eta_proj_ranges_coarse_incl_gap.front().first);
        for (const auto& r : cfg.pair_eta_proj_ranges_coarse_incl_gap) eta.push_back(r.second);
        // A FOLDED mode's fit file carries the |eta| group axis (0 -> eta_max), not the signed
        // 9-bin one. Built with the SAME MakeDrEtaGroups fold the fit stage used, so it cannot
        // describe a different grouping.
        if (DrCorrModeMergeEta(mode)) {
            TAxis src((int)eta.size() - 1, eta.data());
            eta = MakeDrEtaGroups(&src, true).edges;
        }

        auto check = [](const TAxis* ax, const std::vector<double>& edges, const std::string& what) {
            if (ax->GetNbins() != (int)edges.size() - 1)
                throw std::runtime_error("DrCorrectionCrossxEvaluator: the fit file has "
                    + std::to_string(ax->GetNbins()) + " " + what + " cells but the canonical "
                    "binning has " + std::to_string(edges.size() - 1)
                    + " -- stale fit file, or the wrong plateau mode");
            for (size_t i = 0; i < edges.size(); ++i) {
                const double got = (i + 1 <= (size_t)ax->GetNbins())
                                 ? ax->GetBinLowEdge(i + 1) : ax->GetBinUpEdge(ax->GetNbins());
                if (std::fabs(got - edges[i]) > 1e-6)
                    throw std::runtime_error("DrCorrectionCrossxEvaluator: " + what + " edge "
                        + std::to_string(i) + " is " + std::to_string(got)
                        + " in the fit file but " + std::to_string(edges[i])
                        + " canonically -- stale fit file");
            }
        };
        for (const DrCorrectionEvaluator* e : {&expo, &poly, &interp}) {
            check(e->h_fit_ok->GetXaxis(), pt,  "pair-pT (" + e->method + ")");
            check(e->h_fit_ok->GetYaxis(), eta, "pair-eta (" + e->method + ")");
        }
        pt_cell_edges  = pt;
        eta_cell_edges = eta;
    }
};

#endif  // DR_CORRECTION_CROSSX_EVALUATOR_H
