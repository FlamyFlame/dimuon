#ifndef PAIR_TRIG_EFF_CROSSX_EVALUATOR_H
#define PAIR_TRIG_EFF_CROSSX_EVALUATOR_H

#include <atomic>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <TAxis.h>
#include <TString.h>

#include "../MuonObjectsParamsAndHelpers/ParamsSet.h"
#include "../plotting_codes/trig_effcy/mc_based/dr_correction_sample_cfg.h"
#include "DrCorrectionCrossxEvaluator.h"
#include "PairTrigEffEvaluator.h"
#include "SingleMuEffEvaluator.h"

// =================================================================================================
// PairTrigEffCrossxEvaluator -- THE pp24 PER-PAIR 2mu4 TRIGGER EFFICIENCY, as the cross-section
// applies it (docs/tracking/pp24_trig_eff_hybrid_application.md, Physics Procedure §2):
//
//   pair pT in [p_lo, p_split)  REGION A:  eps_trig = eps^nc_data(1) eps^nc_data(2) eps_dR(dR; cell)
//   pair pT in [p_split, p_hi)  REGION B:  eps_trig = eps^pair_MC(cell) * SF(1) * SF(2)
//                                          SF(i) = eps^nc_data(pT_i, q_i eta_i) / eps_MC(pT_i, q_i eta_i)
//   w_trig = 1 / eps_trig
//
// with p_lo, p_split, p_hi = ParamsSet::pair_pt_coarse_bins edges 0, N-2 and N (9, 74.24, 150 GeV
// on the 8-bin axis): region B is the LAST TWO canonical coarse bins. Never retyped.
//
// REGION A is DrCorrectionCrossxEvaluator (expo primary / polynomial primary in the user-named
// forward cells / interpolation fallback / no raw tier; opposite sign; 3-group |eta^pair| fold;
// no a-priori plateau). The single-muon efficiencies are the DATA tag-and-probe turn-ons, passed
// in by the caller (they are RDF columns already).
//
// REGION B is the SINGLE-VALUE pair 2mu4 efficiency measured on the pp24 fullsim in the SIGNAL
// mass window (Utilities/PairTrigEffEvaluator.h, PURE form S1/S0) -- a PAIR-level efficiency, NOT
// a dR correction. It therefore replaces the whole eps1*eps2*eps_dR product and must never be
// multiplied by the single-muon efficiencies. Being pure MC it carries the MC/data pair
// over-efficiency (<r1 r2> = 1.27, mc_trig_eff_closure.md R1); that is removed by the product of
// the two single-muon DATA/MC SCALE FACTORS, evaluated per leg at the exact (pT, q*eta) from the
// SAME two turn-on sets the analysis already has (data T&P; MC direct = FitMCSinglesEffcy, the
// eps_MC the Step-3 inverse weighting divides by). Per pair, not cell-averaged: the user's
// decision (mc_trigeff_single_value_pair_eff.md D4) is explicitly NOT the calibrated K form.
//
// FOUR STATED LIMITATIONS, all user decisions of 2026-09-17 (doc D1-D4):
//   D1 TEMPORARY: one region-B cell ([105.53,150) x |eta| [2.0,2.2), 3 raw MC pairs) is refused by
//      the delivery gate and is served from the pT-MERGED cell of its |eta| group
//      (PairTrigEffCrossxRefusedCellFallbackMode). Revert when the requested high-pT MC arrives.
//   D2 TEMPORARY: SAME-SIGN pairs get these OPPOSITE-sign numbers (the SS cells above 74 GeV are
//      unmeasurable in the present MC). Revert to the SS series when the SS-filtered MC arrives.
//   D3: the region-B number is the SIGNAL-WINDOW one for EVERY pair, including pairs outside
//      1.08-2.9 GeV in the template-fit and generic passes (the 1-4 GeV window differs by 10-21 %).
//   D4: no raw-bin tier anywhere; a region-A cell with nothing usable throws at load time.
//
// GUARDS: the region split is checked to be a cell edge of BOTH inputs; the delivered eps_trig is
// capped at 1 (an efficiency) with the cap COUNTED; pairs outside [p_lo, p_hi) get NO pair-level
// correction (eps1*eps2 only, as before this change) and are COUNTED -- the signal region has no
// upper pair-pT cut, so a handful of pairs above 150 GeV exist and land in histogram overflow.
// =================================================================================================
struct PairTrigEffCrossxEvaluator {

    DrCorrectionCrossxEvaluator dr;                 // region A
    PairTrigEffEvaluator        sv;                 // region B, the nominal (un-merged) cells
    PairTrigEffEvaluator        sv_fallback;        // region B, D1 -- loaded only if configured
    bool                        have_fallback = false;
    SingleMuEffEvaluator        eps_mc;             // for the SF

    double p_lo = 0., p_split = 0., p_hi = 0.;
    int    n_pt_total = 0;                          // N
    // Region-B routing, decided ONCE at load time, indexed [pt bin - (N-2) - 1][eta group - 1]:
    // 0 = nominal cell, 1 = the D1 merged fallback.
    std::vector<std::vector<int>> route_b;
    // D1 is a decision about ONE named cell (doc §3c: "Any other refusal THROWS"). More than one
    // refused cell means the single-value file changed underneath this configuration, and the
    // user must re-decide rather than have the merged cell quietly serve a second one.
    static constexpr int kMaxFallbackCells = 1;
    int n_cells_fallback = 0;

    // ATOMIC: Eval() runs inside RDF lambdas under ImplicitMT.
    std::atomic<long long> n_eval{0}, n_a{0}, n_b{0}, n_b_fallback{0}, n_outside{0}, n_cap{0};

    void Load(const DrCorrSample& cfg, bool use_tight_wp)
    {
        static const ParamsSet pms{};
        const std::vector<double>& e = pms.pair_pt_coarse_bins;
        n_pt_total = (int)e.size() - 1;
        if (n_pt_total < 3)
            throw std::runtime_error("PairTrigEffCrossxEvaluator: pair_pt_coarse_bins has fewer than 3 bins");
        p_lo = e.front(); p_split = e[n_pt_total - 2]; p_hi = e.back();

        // ---- region A
        dr.Load(cfg, use_tight_wp);
        if (dr.n_pt_region_a != n_pt_total - 2 ||
            std::fabs(dr.pt_cell_edges[dr.n_pt_region_a] - p_split) > 1e-6)
            throw std::runtime_error("PairTrigEffCrossxEvaluator: the dR cascade's region A does not end "
                                     "at the canonical split edge -- stale fit file");

        // ---- region B
        const std::string sv_file = PairTrigEff::FileName(cfg.sample_dir, cfg.mc_label,
                                                          DrCorrWpSuffix(use_tight_wp));
        sv.Load(sv_file, PairTrigEffCrossxSign(), PairTrigEffCrossxWindow(),
                PairTrigEffEvaluator::ApplyForm::kPure, PairTrigEffCrossxMode());
        const std::string fb_mode = PairTrigEffCrossxRefusedCellFallbackMode();
        have_fallback = !fb_mode.empty();
        if (have_fallback)
            sv_fallback.Load(sv_file, PairTrigEffCrossxSign(), PairTrigEffCrossxWindow(),
                             PairTrigEffEvaluator::ApplyForm::kPure, fb_mode);
        // The single-value file's pair-pT axis is ParamsSet::pair_pt_coarse_bins (checked inside
        // its loader), so the split edge is one of its edges by construction; the first DELIVERED
        // bin of that file must not lie above region B's first bin, or region B would be refused
        // wholesale.
        {
            const TAxis* ax = sv.Hist()->GetXaxis();
            const int first_b = ax->FindFixBin(p_split);      // bins are [lo, hi): the bin STARTING at the split
            if (first_b != n_pt_total - 1 || std::fabs(ax->GetBinLowEdge(first_b) - p_split) > 1e-6)
                throw std::runtime_error("PairTrigEffCrossxEvaluator: the region split " + std::to_string(p_split)
                    + " is not a cell edge of the single-value file");
            if (PairTrigEff::FirstDeliveredPtBin(PairTrigEffCrossxMode()) > first_b)
                throw std::runtime_error("PairTrigEffCrossxEvaluator: the single-value file delivers from bin "
                    + std::to_string(PairTrigEff::FirstDeliveredPtBin(PairTrigEffCrossxMode()))
                    + " but region B starts at bin " + std::to_string(first_b));
        }

        // ---- the MC single-muon efficiency for the SF
        eps_mc.Load(SingleMuEffEvaluator::Src::kMCDirect, DrCorrSinglesFitFile(cfg, use_tight_wp));

        // ---- region-B census: every cell either delivered, served by the D1 fallback, or a throw.
        const TAxis* ax = sv.Hist()->GetXaxis();
        const TAxis* ay = sv.Hist()->GetYaxis();
        const int neta = ay->GetNbins();
        route_b.assign(2, std::vector<int>(neta, 0));
        std::cout << "PairTrigEffCrossxEvaluator [" << cfg.mc_label << "]: region A = [" << p_lo << ", "
                  << p_split << ") GeV (dR procedure), region B = [" << p_split << ", " << p_hi
                  << ") GeV (single-value pair efficiency x SF1 x SF2, " << sv.Describe() << ")" << std::endl;
        for (int k = 0; k < 2; ++k) {
            const int ix = n_pt_total - 1 + k;             // 1-based bin of region B
            const double pt = ax->GetBinCenter(ix);
            for (int iz = 1; iz <= neta; ++iz) {
                const double eta = ay->GetBinCenter(iz);
                const std::string cell = Form("pair pT [%g,%g) x |eta^pair| [%g,%g)",
                                              ax->GetBinLowEdge(ix), ax->GetBinUpEdge(ix),
                                              ay->GetBinLowEdge(iz), ay->GetBinUpEdge(iz));
                const auto st = sv.Status(pt, eta);
                if (st == PairTrigEffEvaluator::Reject::kDelivered) {
                    std::cout << "  cell " << cell << ": eps^pair = " << sv.Eval(pt, eta) << " +- "
                              << sv.Error(pt, eta) << " (" << sv.NRawPass(pt, eta) << " / "
                              << sv.NRaw(pt, eta) << " raw pairs)" << std::endl;
                    continue;
                }
                if (have_fallback && sv_fallback.Covered(pt, eta)) {
                    route_b[k][iz - 1] = 1;
                    if (++n_cells_fallback > kMaxFallbackCells)
                        throw std::runtime_error("PairTrigEffCrossxEvaluator: a SECOND region-B cell ("
                            + cell + ") is refused by the delivery gate. The D1 fallback is a "
                            "user decision about exactly one cell -- the single-value file has "
                            "changed; re-decide (docs/tracking/pp24_trig_eff_hybrid_application.md D1).");
                    std::cout << "  ** cell " << cell << ": REFUSED by the delivery gate ("
                              << PairTrigEffEvaluator::StatusText(st) << ", " << sv.NRaw(pt, eta)
                              << " raw pairs) -> TEMPORARY D1 fallback to the '" << fb_mode
                              << "' cell of this |eta| group: eps^pair = " << sv_fallback.Eval(pt, eta)
                              << " +- " << sv_fallback.Error(pt, eta) << " (" << sv_fallback.NRawPass(pt, eta)
                              << " / " << sv_fallback.NRaw(pt, eta) << " raw pairs). Revert when the "
                                 "requested high-pT MC statistics arrive." << std::endl;
                    continue;
                }
                throw std::runtime_error("PairTrigEffCrossxEvaluator: region-B cell " + cell
                    + " is not delivered (" + PairTrigEffEvaluator::StatusText(st) + ")"
                    + (have_fallback ? " and the '" + fb_mode + "' fallback cell is not delivered either"
                                     : " and no fallback is configured")
                    + " -- the cross-section cannot be corrected there.");
            }
        }
    }

    // eps_trig^pair for one pair. `e1d`/`e2d` are the DATA single-muon turn-ons the caller already
    // evaluated (RDF columns effcy1/effcy2); the leg kinematics are needed for the MC twin.
    double Eval(double dr_, double pair_pt, double pair_eta, double e1d, double e2d,
                float pt1, float eta1, int q1, float pt2, float eta2, int q2)
    {
        ++n_eval;
        if (!(pair_pt >= p_lo && pair_pt < p_hi)) {
            // No cell covers the pair: no PAIR-level correction, singles product only (the
            // behaviour every pair outside the coarse grid had before this change). Counted.
            ++n_outside;
            return e1d * e2d;
        }
        double v;
        if (pair_pt < p_split) {
            ++n_a;
            v = e1d * e2d * dr.Eval(dr_, pair_pt, pair_eta);
        } else {
            ++n_b;
            const int ix = sv.Hist()->GetXaxis()->FindFixBin(pair_pt);
            const int iz = sv.Hist()->GetYaxis()->FindFixBin(std::fabs(pair_eta));
            const int k  = ix - (n_pt_total - 1);
            // The |eta^pair| < 2.2 pair-level fiducial cut puts every selected pair inside the
            // group axis; a pair outside it has no cell and no correction was ever defined for it.
            if (k < 0 || k > 1 || iz < 1 || iz > (int)route_b[0].size()) { ++n_outside; return e1d * e2d; }
            const bool fb = route_b[k][iz - 1] == 1;
            if (fb) ++n_b_fallback;
            const double eps_pair = fb ? sv_fallback.Eval(pair_pt, pair_eta) : sv.Eval(pair_pt, pair_eta);
            // Numerator and denominator are evaluated by their own conventions: the data turn-on
            // (caller) is the TFormula at the exact pT, uncapped in pT, cap 1, floor 0.01; the MC
            // twin clamps pT into its fit range [4.5, 60], cap 1, floor 0.02. Measured effect for
            // legs above 60 GeV: <= 1 % per leg (doc R5); floors never active (legs >= 4.5 GeV).
            const double sf1 = e1d / eps_mc.Eval(pt1, eta1, q1);
            const double sf2 = e2d / eps_mc.Eval(pt2, eta2, q2);
            v = eps_pair * sf1 * sf2;
        }
        if (v > 1.0) { v = 1.0; ++n_cap; }
        return v;
    }

    void PrintStats()
    {
        const long long N = n_eval;
        auto pct = [&](long long n) { return N ? 100.0 * n / N : 0.0; };
        std::cout << "PairTrigEffCrossxEvaluator: " << N << " evaluations -- region A " << n_a << " ("
                  << pct(n_a) << "%), region B " << n_b << " (" << pct(n_b) << "%; of which "
                  << n_b_fallback << " on the TEMPORARY D1 pT-merged fallback cell), "
                  << n_outside << " outside [" << p_lo << ", " << p_hi << ") GeV or the |eta| groups ("
                  << pct(n_outside) << "%, singles product only), " << n_cap
                  << " capped at eps_trig = 1 (" << pct(n_cap) << "%)" << std::endl;
        dr.PrintStats();
        eps_mc.PrintStats("region-B SF denominator");
    }
};

#endif  // PAIR_TRIG_EFF_CROSSX_EVALUATOR_H
