#ifndef DR_CORRECTION_APPLY_H
#define DR_CORRECTION_APPLY_H

#include <algorithm>
#include <atomic>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TString.h>

#include "dr_correction_sample_cfg.h"
#include "dr_correction_ratio.h"
#include "dr_correction_cell_groups.h"

// =================================================================================================
// APPLYING eps_dR: turn the per-(pair pT, pair eta) Step-3 fits into a number for ONE pair.
//
// This is the CONSUMER side of fit_dr_corrections.cxx. The fit stage measures and fits; this reads
// the persisted result back and answers "what is eps_dR for a pair at this dR, in this cell?".
//
// THE FORM APPLIED (user, 2026-08-13; docs/tracking/mc_trig_eff_closure.md §3.2):
//
//     eps_dR(dR) = f_cell(dR) / C_cell     for dR <  1
//     eps_dR(dR) = 1                       for dR >= 1
//
// taken from the NO-PLATEAU-CORRECTION fit, whose shape f(dR) = C + A exp[-(dR/lambda)^p] carries
// its own free asymptote C determined by the dR < 1 data alone. Dividing by C afterwards is what
// puts the plateau at 1, which is what the factorised per-pair weight eps1*eps2*eps_dR assumes.
// Because f(dR) -> C well before dR = 1, f(1)/C ~ 1 and there is NO STEP at the dR = 1 boundary
// where the correction hands over to unity.
//
// THE THREE CURVE FLAVOURS THIS STRUCT CAN DELIVER (one per `method`), all in the SAME form
// f(dR)/C, so a cascade can put them in one chain without special-casing any of them:
//   `expo`           f = C + A exp[-(dR/lambda)^p]          C = the fitted free asymptote
//   `polyu_fixedRp`  f = C + u^2[A + a3(u-1) + a4(u^2-1)]   C = the fitted free asymptote
//                    (the same quartic as C + a2 u^2 + a3 u^3 + a4 u^4, with
//                    a2 = A - a3 - a4; carried in A = f(0)-C since 2026-09-08 so the
//                    Step-3 restriction f(0) <= C is a limit on one parameter)
//   `interp`         linear interpolation through the measured points below R_p, FLAT at the last
//                    measured knot above it (fit_dr_corrections.cxx's `flat_val`)
//                                                          C = that flat-branch value
// The interpolation was added 2026-08-24 (user; docs/tracking/mc_trigeff_dr_binning_approaches.md
// §PP-3) because it is the third tier of the delivered cascade, REPLACING the raw-bin placeholder
// below. It is one of the three resolutions mc_trigger_efficiency.md R26 names ("per-cell
// interpolation"), and it enters on exactly the same footing as the two fits: a curve divided by
// its own baseline. It carries no parameter errors, so its baseline screen is
// DrCorrPlateauUsable(C, 0) -- i.e. the "is C a sane normalization?" half of the test, which is
// the only half that is defined without an error.
//
// WHICH VARIANT pp24 CROSSX APPLIES (user, 2026-09-17; named in dr_correction_sample_cfg.h so no
// consumer retypes it): plateau mode `nocorr_etamerge` = DrCorrCrossxMode() (the
// no-plateau-correction fit on the 3-group |eta^pair| fold, un-merged pair pT), series `os` =
// DrCorrCrossxSign(), primary `expo` = DrCorrCrossxMethod() with `polyu_fixedRp` primary in three
// named forward cells and `interp` as the ONLY fallback -- and ONLY in coarse pair-pT bins 1..N-2;
// the last two bins are served by the single-value pair efficiency instead
// (Utilities/PairTrigEffCrossxEvaluator.h, docs/tracking/pp24_trig_eff_hybrid_application.md).
// The raw-bin placeholder below is therefore UNREACHABLE from the cross-section; it remains for
// the MC closure. Load() still DEFAULTS to the un-merged "nocorr", so the closure thread keeps
// consuming exactly what it consumed before.
//
// WHY THE NO-PLATEAU-CORRECTION VARIANT and not the nominal plateau-normalized one: the measured
// eps_dR carries structure out to large dR, worst in the pair-eta bins that enclose the detector
// gap, so the far-away [2, 3.5] plateau is not necessarily the right baseline for the small-dR
// region the correction is about (mc_trigger_efficiency.md R27: one cell has C = 0.849 +- 0.013
// against a plateau of 0.997 -- a 15% difference in exactly that region).
//
// ---------------------------------------------------------------------------------------------
// ⚠ TEMPORARY PLACEHOLDER -- THE RAW-BIN FALLBACK. A cell whose fit was REJECTED
// (h_step3_fit_ok == 0) falls back to the RAW MEASURED eps_dR of that cell, normalized by its own
// mean over dR in [0.5, 1.0]:
//
//     eps_dR(dR) = eps_dR^raw(bin containing dR) / C_raw ,  dR < 1
//
// C_raw plays exactly the role the fitted C plays -- without it the raw values keep the cell's own
// normalization offset and the correction would step at dR = 1.
//
// THIS IS A STOP-GAP, NOT THE FINAL SOLUTION. It exists only so the closure can be produced at all
// while mc_trigger_efficiency.md R26 is OPEN and awaiting a user decision (both parametric fit
// forms fail badly in a substantial minority of cells, and `usable` carries no chi2 term, so a
// badly-fitted cell is published as fit_ok = 1 while some perfectly measurable cells are rejected).
// Once R26 is decided -- a chi2 screen in `usable`, per-cell interpolation, or a fit form that can
// describe the measured shape -- this branch must be REPLACED, not extended.
// ---------------------------------------------------------------------------------------------
//
// Guards: dR is clamped into the stored TF1 range before evaluation (a TF1 read back from a file
// returns 0 outside it -- pp_trig_eff_highpt_jump.md), and the returned correction is floored at
// kMinCorr and capped at kMaxCorr with every firing COUNTED and printed. A correction <= 0 would
// flip the sign of a per-pair weight; a silent clamp would hide it.
//
// THE BASELINE SCREEN. `fit_ok == 1` is necessary but NOT sufficient here: it says the fit
// converged and is positive, and it carries no chi2 term (parent doc R26, OPEN). Since C is what
// the whole cell is DIVIDED BY, it gets the repo's existing unusable-normalization screen
// `DrCorrPlateauUsable(C, err_C)` on top -- the same test, for the same reason, that
// dr_correction_sample_cfg.h applies to a measured plateau. A cell that fails it is NOT silently
// used: it is announced and falls through to the raw-bin placeholder.
// =================================================================================================
struct DrCorrectionEvaluator {

    // Beyond this dR no correction is measured and none is assumed. MUST equal the fit domain's
    // upper edge (kFitHi in fit_dr_corrections.cxx) -- the fits contain no information above it.
    static constexpr double kDrMax   = 1.0;
    static constexpr double kMinCorr = 0.02;   // a trigger correction can never be <= 0
    static constexpr double kMaxCorr = 5.0;    // pathology guard, not a physics choice

    struct Cell {
        TF1*    fit   = nullptr;   // parametric methods; owned by the open TFile
        TGraph* knots = nullptr;   // `interp`; owned by the open TFile
        double  C     = 0.;        // the baseline of whichever of the two is set (see the header)
        TH1D*   raw   = nullptr;   // TEMPORARY PLACEHOLDER: raw measured eps_dR of this cell
        double  C_raw = 0.;        // its dR in [0.5, 1.0] mean
        // "does this cell deliver a measured CURVE?" -- the question every cascade branches on.
        bool HasCurve() const { return fit != nullptr || knots != nullptr; }
    };

    std::vector<std::vector<Cell>> cells;      // [iy-1][iz-1]
    std::unique_ptr<TH2D> h_fit_ok;            // also supplies the cell axes
    std::unique_ptr<TH2D> h_plateau;           // the cell's OWN measured plateau; screens C
    std::string method, sign, label, mode;
    // Set from DrCorrModeMergeEta(mode) at Load() time. A folded |eta| grouping's Y axis (of
    // h_fit_ok) runs 0 -> eta_max, NOT -eta_max -> eta_max, so Eval() below must look the cell up
    // by |pair_eta|, never by the signed value -- a signed lookup on that axis would send every
    // negative-eta pair to the underflow bin and silently deliver eps_dR = 1 (no correction) to
    // half the sample.
    bool eta_folded = false;

    // ATOMIC (2026-09-17): Eval() runs inside RDF lambdas under ImplicitMT (the crossx and the
    // closure); plain ++ lost ~0.4 % of the counts (outer atomic census 1250918 vs inner 1245904)
    // and made the "0 floored / 0 capped / 0 outside (must be 0)" guards untrustworthy. The
    // delivered values were never affected (const TF1 / histogram reads).
    std::atomic<long long> n_eval{0}, n_flat{0}, n_fit{0}, n_raw{0}, n_rawempty{0}, n_none{0};
    std::atomic<long long> n_floor{0}, n_cap{0}, n_outside{0};
    // Per-CELL census (as opposed to per-evaluation), so a consumer of the output file can be told
    // how much of the map is a placeholder rather than a fit.
    int n_cells_fitted = 0, n_cells_raw = 0, n_cells_dead = 0;
    int n_cells_bad_C  = 0;   // fit_ok = 1 but the fitted baseline C failed the screen
    // Extrema of the DELIVERED correction over the accepted cells (reported, not screened).
    double max_corr = -1e300, min_corr = 1e300;
    int max_corr_iy = 0, max_corr_iz = 0, min_corr_iy = 0, min_corr_iz = 0;

    // ------------------------------------------------------------------ construction
    // `sign` is the series token ("os" for opposite sign). `plateau_mode` must be one of the
    // NO-PLATEAU-CORRECTION modes -- the applied form above is defined in terms of the FREE
    // baseline C, which only those fits have:
    //   "nocorr"          the un-merged cells (the default; the reference of the 4-approach
    //                     comparison in docs/tracking/mc_trigeff_dr_binning_approaches.md)
    //   "nocorr_etamerge" the same fit with the 9 pair-eta bins merged into 3 sign-independent
    //                     |eta^pair| bins: |eta| < 1.0 (barrel), 1.0 <= |eta| < 2.0 and
    //                     2.0 <= |eta| < the source axis's top edge (2.2 since 2026-09-07,
    //                     tracking ParamsSet::pair_eta_fiducial_max; see
    //                     dr_correction_cell_groups.h). Eval() below looks
    //                     the cell up by |pair_eta| whenever this mode is loaded.
    //   "nocorr_etamerge_ptmerge"  both merges at once
    //   "nocorr_ptmerge"  the same fit with the LAST TWO pair-pT bins merged into one cell. It
    //                     WAS the pp24 crossx variant from 2026-08-17 to 2026-09-17 (its top cell
    //                     covered p_T^pair in [72.1, 150) GeV, where the un-merged pair of cells
    //                     runs past the sample's yield -- mc_trigger_efficiency.md R24/R27, and
    //                     the closure collapse in mc_trig_eff_closure.md R1b); the crossx now
    //                     uses `nocorr_etamerge` below 74 GeV and the single-value pair
    //                     efficiency above (DrCorrCrossxMode(), PairTrigEffCrossxEvaluator.h).
    // The DEFAULT is deliberately the un-merged mode: changing it would silently move every
    // existing consumer (the MC closure) onto a different correction.
    void Load(const DrCorrSample& cfg, bool use_tight_wp, const std::string& method_in,
              const std::string& sign_in, const std::string& plateau_mode = "nocorr")
    {
        method = method_in;
        sign   = sign_in;
        mode   = plateau_mode;
        label  = cfg.mc_label;
        eta_folded = DrCorrModeMergeEta(plateau_mode);

        if (!DrCorrModeNoPlateau(plateau_mode))
            throw std::runtime_error("DrCorrectionEvaluator: plateau mode '" + plateau_mode +
                                     "' has no free baseline C, so the f(dR)/C form is undefined "
                                     "for it. Use 'nocorr' or 'nocorr_ptmerge'.");

        const std::string fit_path =
            DrCorrFitFile(cfg, use_tight_wp, 3, method, sign, plateau_mode);
        TFile* ff = TFile::Open(fit_path.c_str(), "READ");
        if (!ff || ff->IsZombie())
            throw std::runtime_error("DrCorrectionEvaluator: cannot open " + fit_path +
                                     " -- run fit_dr_corrections with plateau mode '"
                                     + plateau_mode + "' and sign '" + sign + "' first");

        auto* ok = dynamic_cast<TH2D*>(ff->Get("h_step3_fit_ok"));
        if (!ok)
            throw std::runtime_error("DrCorrectionEvaluator: no h_step3_fit_ok in " + fit_path);
        h_fit_ok.reset(static_cast<TH2D*>(ok->Clone("h_fit_ok_clone")));
        h_fit_ok->SetDirectory(nullptr);

        // The cell's OWN MEASURED plateau, written by the fit stage beside fit_ok. It is the
        // reference the fitted baseline C is screened against below (user decision 2026-09-10).
        // Absent only in a fit file older than that stage, in which case the ratio screen is
        // skipped rather than guessed -- see the comment at the screen itself.
        auto* plat_h = dynamic_cast<TH2D*>(ff->Get("h_step3_plateau"));
        if (plat_h) {
            h_plateau.reset(static_cast<TH2D*>(plat_h->Clone("h_plateau_clone")));
            h_plateau->SetDirectory(nullptr);
        } else {
            std::cout << "  ** WARNING: no h_step3_plateau in " << fit_path
                      << " -- the C-vs-measured-plateau screen is INACTIVE for this file."
                      << std::endl;
        }

        const int npt = h_fit_ok->GetNbinsX(), neta = h_fit_ok->GetNbinsY();
        cells.assign(npt, std::vector<Cell>(neta));
        // The cell's own measured plateau, or 0 when the fit file predates it (screen inapplicable).
        auto measured_plateau = [&](int ix, int iz_) -> double {
            return h_plateau ? h_plateau->GetBinContent(ix, iz_) : 0.;
        };

        // The RAW measured curves come from the SAME histogram file the fit stage read, projected
        // with the SAME per-cell projection + conditional errors (dr_correction_ratio.h). Opened
        // lazily -- a run in which every cell is fitted never needs it.
        TFile* fh = nullptr;
        TH3D *h3n = nullptr, *h3d = nullptr, *h3a = nullptr, *h3b = nullptr;
        // The pair-pT GROUPING of the fit cells, needed by the raw-bin fallback below: with the
        // last two bins merged, cell `iy` at the top spans TWO filled bins, and projecting bin
        // `iy` alone would silently deliver half of it. Built from the histograms' own axis, so it
        // is filled in together with them (dr_correction_cell_groups.h).
        // ... and the pair-ETA grouping, for the same reason: with the 9 filled pair-eta bins
        // folded into 3 sign-independent |eta| bins, a forward cell `iz` spans TWO disjoint
        // filled-bin sub-ranges (a negative-eta one and a positive-eta one).
        DrAxisGroups Gpt, Geta;
        const std::string hist_path = DrCorrHistFile(cfg, use_tight_wp,
                                                     MCTrigEffPairPt::FileSuffix() + "_step3");
        const std::string hp = "h_mc_dr_" + (sign.empty() ? std::string() : sign + "_")
                             + "zoom_vs_pt_eta_";

        int& n_fitted  = n_cells_fitted;
        int& n_fallback = n_cells_raw;
        int& n_dead    = n_cells_dead;
        for (int iy = 1; iy <= npt; ++iy) {
            for (int iz = 1; iz <= neta; ++iz) {
                Cell& c = cells[iy - 1][iz - 1];
                if (h_fit_ok->GetBinContent(iy, iz) > 0.5 && method == "interp") {
                    // The interpolation is persisted as a TGraph of knots, not a TF1. Its baseline
                    // is the FLAT BRANCH the producer pinned to the last measured knot below R_p,
                    // which is what the fitted C is for the parametric forms -- so `f(dR)/C` means
                    // the same thing here, and there is no step at dR = 1.
                    auto* gk = dynamic_cast<TGraph*>(ff->Get(Form("gknots_step3_pt%d_eta%d",
                                                                  iy, iz)));
                    if (gk && gk->GetN() > 1) {
                        const double C = gk->Eval(kDrMax);   // the flat branch, by construction
                        // No parameter errors exist for an interpolation, so the screen is the
                        // error-free half of DrCorrPlateauUsable: "is C a sane normalization?".
                        // Passing 0 makes its `err >= plateau` clause vacuous, which is the honest
                        // reading -- not a silently weaker test smuggled in as the same one.
                        if (DrCorrPlateauUsable(C, 0.) && DrCorrBaselineConsistent(C, measured_plateau(iy, iz))) {
                            c.knots = gk; c.C = C; ++n_cells_fitted;
                            for (int k = 0; k <= 200; ++k) {
                                const double x = kDrMax * k / 200.0;
                                const double v = gk->Eval(x) / C;
                                if (v > max_corr) { max_corr = v; max_corr_iy = iy; max_corr_iz = iz; }
                                if (v < min_corr) { min_corr = v; min_corr_iy = iy; min_corr_iz = iz; }
                            }
                            continue;
                        }
                        ++n_cells_bad_C;
                        std::cout << "  ** cell (pair pT bin " << iy << ", pair eta bin " << iz
                                  << ") has fit_ok = 1 but an UNUSABLE interpolation baseline C = "
                                  << C << " -> rejected, falling back to the raw bins" << std::endl;
                    }
                } else if (h_fit_ok->GetBinContent(iy, iz) > 0.5) {
                    auto* f = dynamic_cast<TF1*>(ff->Get(Form("f_step3_pt%d_eta%d", iy, iz)));
                    // The free baseline is the LAST parameter in both nocorr formulas
                    // (fit_dr_corrections.cxx MakeMethodCfg: expo -> [3], polyu_fixedRp -> [4]).
                    if (f && f->GetNpar() > 0) {
                        const int    ip = f->GetNpar() - 1;
                        const double C  = f->GetParameter(ip);
                        const double eC = f->GetParError(ip);
                        // C IS THE NORMALIZATION, so it needs the SAME sanity screen the repo
                        // already applies to a measured plateau -- `DrCorrPlateauUsable`, whose
                        // own comment states the reason exactly: "eps_dR is a ratio that must tend
                        // to 1 at large dR, so a plateau far below 1 is not a normalization offset
                        // but an empty cell; dividing by 0.01 inflates the curve 100x rather than
                        // normalizing it."  `C > 0` alone is NOT enough: a `polyu_fixedRp` cell was
                        // accepted with C = 0.0534, which inflated f/C to 21.1 over dR < 1, was
                        // truncated by kMaxCorr, and drove one delivered closure point to
                        // 0.193 +- 0.031 -- 26 sigma from 1, an artefact of the fit, not a
                        // measurement (mc_trig_eff_closure.md R2; parent doc R26: `usable` carries
                        // no chi2 term, so such a fit is published with fit_ok = 1).
                        // A cell rejected here falls through to the raw-bin placeholder below,
                        // which is the same route any other rejected fit takes.
                        // ...and C must also be CONSISTENT WITH THE PLATEAU THIS CELL MEASURED.
                        // DrCorrPlateauUsable bounds C only from below; a runaway C passes it and
                        // silently inflates 1/eps_dR. See DrCorrBaselineConsistent for the case
                        // that forced this (C = 3.98 against a measured plateau of 0.82).
                        if (DrCorrPlateauUsable(C, eC) && DrCorrBaselineConsistent(C, measured_plateau(iy, iz))) {
                            c.fit = f; c.C = C; ++n_fitted;
                            // The screen bounds the DENOMINATOR C, not the delivered correction
                            // f/C. Those are different questions, and a C that only just clears
                            // the screen (e.g. 0.553 +- 0.339, a 61 % relative error) can still
                            // put f/C above 2. Whether a hard bound on the correction itself
                            // should be imposed is a PHYSICS decision entangled with the parent
                            // doc's open R26, so nothing is rejected on it here -- but the
                            // extremum is measured and PRINTED, so it can never be invisible.
                            for (int k = 0; k <= 200; ++k) {
                                const double x = kDrMax * k / 200.0;
                                const double v = f->Eval(x) / C;
                                if (v > max_corr) { max_corr = v; max_corr_iy = iy; max_corr_iz = iz; }
                                if (v < min_corr) { min_corr = v; min_corr_iy = iy; min_corr_iz = iz; }
                            }
                            continue;
                        }
                        ++n_cells_bad_C;
                        const double mp = measured_plateau(iy, iz);
                        const bool   runaway = !DrCorrBaselineConsistent(C, mp);
                        std::cout << "  ** cell (pair pT bin " << iy << ", pair eta bin " << iz
                                  << ") has fit_ok = 1 but an UNUSABLE fitted baseline C = "
                                  << C << " +- " << eC;
                        if (runaway)
                            std::cout << " -- INCONSISTENT with the plateau this cell MEASURED ("
                                      << mp << "; ratio " << (mp > 0. ? C / mp : -1.)
                                      << ", allowed " << 1.0 / kDrCorrBaselineMaxRatio << " - "
                                      << kDrCorrBaselineMaxRatio << ")";
                        std::cout << " -> rejected, falling back to the next tier" << std::endl;
                    }
                }
                // ---- TEMPORARY PLACEHOLDER (see the header comment) ----
                if (!fh) {
                    fh = TFile::Open(hist_path.c_str(), "READ");
                    if (!fh || fh->IsZombie())
                        throw std::runtime_error("DrCorrectionEvaluator: a cell has no usable fit "
                                                 "and the raw fallback needs " + hist_path);
                    h3n = dynamic_cast<TH3D*>(fh->Get((hp + "num").c_str()));
                    h3d = dynamic_cast<TH3D*>(fh->Get((hp + "denom").c_str()));
                    h3a = dynamic_cast<TH3D*>(fh->Get((hp + "errA").c_str()));
                    h3b = dynamic_cast<TH3D*>(fh->Get((hp + "errB").c_str()));
                    if (!h3n || !h3d || !h3a || !h3b)
                        throw std::runtime_error("DrCorrectionEvaluator: missing " + hp +
                                                 "{num,denom,errA,errB} in " + hist_path);
                    Gpt  = MakeDrPtGroups (h3n->GetYaxis(),
                                           DrCorrModeMergeLastTwoPt(plateau_mode));
                    Geta = MakeDrEtaGroups(h3n->GetZaxis(),
                                           DrCorrModeMergeEta(plateau_mode));
                    if (Gpt.n != npt || Geta.n != neta)
                        throw std::runtime_error("DrCorrectionEvaluator: the fit file has "
                            + std::to_string(npt) + "x" + std::to_string(neta)
                            + " cells but plateau mode '" + plateau_mode + "' groups " + hist_path
                            + " into " + std::to_string(Gpt.n) + "x" + std::to_string(Geta.n)
                            + " -- fit file and mode disagree");
                }
                TH1D* r = DrGroupCellRatio(h3n, h3d, h3a, h3b, nullptr, nullptr, Gpt, Geta, iy, iz,
                                           Form("rawcorr_pt%d_eta%d", iy, iz));
                double sum = 0.; int nb = 0;
                for (int i = 1; i <= r->GetNbinsX(); ++i) {
                    const double x = r->GetBinCenter(i);
                    if (x < 0.5 || x > kDrMax) continue;
                    if (r->GetBinError(i) <= 0.) continue;   // no information in this bin
                    sum += r->GetBinContent(i); ++nb;
                }
                if (nb > 0 && sum > 0.) { c.raw = r; c.C_raw = sum / nb; ++n_fallback; }
                else                    { delete r; ++n_dead; }
            }
        }
        std::cout << "DrCorrectionEvaluator [" << label << " / " << method << " / "
                  << (sign.empty() ? "sign-integrated" : DrCorrSignText(sign)) << " / "
                  << plateau_mode
                  << "]: " << npt << "x" << neta << " cells -- " << n_fitted
                  << (method == "interp" ? " interpolated, " : " fitted, ")
                  << n_fallback << " on the RAW-BIN PLACEHOLDER, " << n_dead
                  << " with no correction at all (eps_dR = 1); of the rejected, "
                  << n_cells_bad_C << " had fit_ok = 1 but an unusable baseline C" << std::endl;
        if (n_fitted > 0)
            std::cout << "     delivered eps_dR over the ACCEPTED cells spans ["
                      << min_corr << " (pair pT bin " << min_corr_iy << ", pair eta bin "
                      << min_corr_iz << "), " << max_corr << " (pair pT bin " << max_corr_iy
                      << ", pair eta bin " << max_corr_iz << ")]"
                      << (max_corr > 1.5 ? "  <-- ABOVE 1.5: a 2mu4 close-by correction is a LOSS,"
                                           " so this is a fit artefact, not a measurement"
                                         : "")
                      << std::endl;
        if (n_fallback > 0)
            std::cout << "  ** " << n_fallback << " cell(s) use the TEMPORARY raw-bin fallback "
                         "(mc_trig_eff_closure.md §3.3) -- not a final solution." << std::endl;
    }

    // ------------------------------------------------------------------ evaluation
    // THE ONE PLACE THE pair-eta CELL LOOKUP VALUE IS FORMED, for every consumer of a Step-3 fit
    // map. A folded |eta| grouping's Y axis runs 0 -> eta_max, so a SIGNED lookup sends every
    // negative-eta pair into the underflow bin, where it is counted "outside the cell grid" and
    // silently gets eps_dR = 1 -- half the sample left uncorrected, with no crash and no warning.
    // That bug was fixed here on 2026-09-03 (mc_trigeff_dr_binning_approaches.md D11) but NOT in
    // DrCorrectionCascadeEvaluator or DrCorrectionCrossxEvaluator, which reimplemented the same
    // two lines; it was found again on 2026-09-08 in the cascade class, where it had corrupted
    // every `*_etamerge*` closure. Hence a shared helper: the fold can no longer be forgotten by
    // the next class that needs a cell index.
    static double CellLookupEta(double pair_eta, bool folded)
    {
        return folded ? std::fabs(pair_eta) : pair_eta;
    }

    double Eval(double dr, double pair_pt, double pair_eta)
    {
        ++n_eval;
        if (dr >= kDrMax) { ++n_flat; return 1.0; }

        const int iy = h_fit_ok->GetXaxis()->FindBin(pair_pt);
        // A folded |eta| grouping's Y axis runs 0 -> eta_max (see `eta_folded` above); the
        // un-merged / pair-pT-only-merged modes keep the signed -eta_max -> eta_max axis.
        const int iz = h_fit_ok->GetYaxis()->FindBin(CellLookupEta(pair_eta, eta_folded));
        if (iy < 1 || iy > h_fit_ok->GetNbinsX() || iz < 1 || iz > h_fit_ok->GetNbinsY()) {
            // Outside the measured cells there is no correction. The caller is expected to have
            // restricted the sample to the cell grid, so this is COUNTED and must come out 0.
            ++n_outside;
            return 1.0;
        }

        const Cell& c = cells[iy - 1][iz - 1];
        double v = 1.0;
        if (c.fit) {
            const double x = std::min(std::max(dr, c.fit->GetXmin()), c.fit->GetXmax());
            v = c.fit->Eval(x) / c.C;
            ++n_fit;
        } else if (c.knots) {
            // A TGraph EXTRAPOLATES LINEARLY outside its knots, so dR is clamped into [0, kDrMax]
            // exactly as it is clamped into a TF1's range above. The producer already pinned
            // dR = 0 to the first measured value for the same reason.
            const double x = std::min(std::max(dr, 0.0), kDrMax);
            v = c.knots->Eval(x) / c.C;
            ++n_fit;
        } else if (c.raw) {
            // A cell on the placeholder, but the individual raw BIN this pair lands in may still
            // carry no information (zero denominator) -> no correction. Counted separately, so the
            // printed placeholder fraction is not inflated by pairs that got eps_dR = 1 anyway.
            const int b = c.raw->FindBin(dr);
            const double y = c.raw->GetBinContent(b);
            if (b >= 1 && b <= c.raw->GetNbinsX() && c.raw->GetBinError(b) > 0. && y > 0.) {
                v = y / c.C_raw;
                ++n_raw;
            } else {
                ++n_rawempty;
            }
        } else {
            ++n_none;
        }

        if (v < kMinCorr) { v = kMinCorr; ++n_floor; }
        if (v > kMaxCorr) { v = kMaxCorr; ++n_cap; }
        return v;
    }

    // The load-time scans (DrCorrectionCrossxEvaluator::ReportDeliveredExtrema) call Eval over a
    // grid of cell centres to bound the DELIVERED correction. Those calls are diagnostics, not
    // pairs: left in the counters they inflate exactly the census that says how much of the map is
    // a placeholder (measured: 99 raw-bin + 101 empty-bin evaluations added to a real 71, i.e. the
    // placeholder share came out 2.2x too high). A caller that scans therefore snapshots the
    // counters and restores them afterwards.
    struct Counters {
        long long n_eval, n_flat, n_fit, n_raw, n_rawempty, n_none, n_floor, n_cap, n_outside;
    };
    Counters SnapshotCounters() const
    {
        return {n_eval, n_flat, n_fit, n_raw, n_rawempty, n_none, n_floor, n_cap, n_outside};
    }
    void RestoreCounters(const Counters& c)
    {
        n_eval = c.n_eval; n_flat = c.n_flat; n_fit = c.n_fit; n_raw = c.n_raw;
        n_rawempty = c.n_rawempty; n_none = c.n_none; n_floor = c.n_floor; n_cap = c.n_cap;
        n_outside = c.n_outside;
    }

    void PrintStats() const
    {
        const Counters c = SnapshotCounters();
        auto pct = [&](long long n) { return c.n_eval ? 100.0 * n / c.n_eval : 0.0; };
        std::cout << "DrCorrectionEvaluator [" << label << " / " << method << "]: "
                  << c.n_eval << " evaluations -- "
                  << c.n_flat << " at dR >= " << kDrMax << " (" << pct(c.n_flat) << "%, no correction), "
                  << c.n_fit  << " from a fit (" << pct(c.n_fit) << "%), "
                  << c.n_raw  << " from the RAW-BIN PLACEHOLDER (" << pct(c.n_raw) << "%), "
                  << c.n_rawempty << " in an empty bin of a placeholder cell (" << pct(c.n_rawempty)
                  << "%, no correction), "
                  << c.n_none << " with no correction available (" << pct(c.n_none) << "%); "
                  << c.n_floor << " floored at " << kMinCorr << ", " << c.n_cap << " capped at "
                  << kMaxCorr << ", " << c.n_outside << " outside the cell grid (must be 0)"
                  << std::endl;
    }
};

#endif  // DR_CORRECTION_APPLY_H
