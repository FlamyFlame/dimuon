#ifndef DR_CORRECTION_APPLY_H
#define DR_CORRECTION_APPLY_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TString.h>

#include "dr_correction_sample_cfg.h"
#include "dr_correction_ratio.h"
#include "dr_correction_pt_groups.h"

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
// WHICH VARIANT pp24 CROSSX APPLIES (user, 2026-08-17; TEMPORARY, named in
// dr_correction_sample_cfg.h so no consumer retypes it): method `expo` = DrCorrCrossxMethod(),
// series `os` (opposite sign) = DrCorrCrossxSign(), plateau mode `nocorr_ptmerge` =
// DrCorrCrossxMode() -- the no-plateau-correction fit with the LAST TWO pair-pT bins merged into a
// single cell covering p_T^pair in [72.1, 150) GeV. Load() still DEFAULTS to the un-merged
// "nocorr", so the MC closure thread keeps consuming exactly what it consumed before; the merged
// variant is requested explicitly.
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
        TF1*   fit = nullptr;      // owned by the open TFile
        double C   = 0.;           // the fit's own free asymptote (last parameter)
        TH1D*  raw = nullptr;      // TEMPORARY PLACEHOLDER: raw measured eps_dR of this cell
        double C_raw = 0.;         // its dR in [0.5, 1.0] mean
    };

    std::vector<std::vector<Cell>> cells;      // [iy-1][iz-1]
    std::unique_ptr<TH2D> h_fit_ok;            // also supplies the cell axes
    std::string method, sign, label, mode;

    long long n_eval = 0, n_flat = 0, n_fit = 0, n_raw = 0, n_rawempty = 0, n_none = 0;
    long long n_floor = 0, n_cap = 0, n_outside = 0;
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
    //   "nocorr"          the un-merged pair-pT cells (the default; what the MC closure consumes)
    //   "nocorr_ptmerge"  the same fit with the LAST TWO pair-pT bins merged into one cell. This
    //                     is the variant the pp24 crossx application uses (DrCorrCrossxMode(),
    //                     with DrCorrCrossxMethod() / DrCorrCrossxSign()); its top cell covers
    //                     p_T^pair in [72.1, 150) GeV, where the un-merged pair of cells runs past
    //                     the sample's yield (mc_trigger_efficiency.md R24/R27, and the closure
    //                     collapse in mc_trig_eff_closure.md R1b).
    // The DEFAULT is deliberately the un-merged mode: changing it would silently move every
    // existing consumer (the MC closure) onto a different correction.
    void Load(const DrCorrSample& cfg, bool use_tight_wp, const std::string& method_in,
              const std::string& sign_in, const std::string& plateau_mode = "nocorr")
    {
        method = method_in;
        sign   = sign_in;
        mode   = plateau_mode;
        label  = cfg.mc_label;

        if (!DrCorrModeNoPlateau(plateau_mode))
            throw std::runtime_error("DrCorrectionEvaluator: plateau mode '" + plateau_mode +
                                     "' has no free baseline C, so the f(dR)/C form is undefined "
                                     "for it. Use 'nocorr' or 'nocorr_ptmerge'.");

        if (method == "interp")
            throw std::runtime_error("DrCorrectionEvaluator: method 'interp' has no fitted "
                                     "baseline C, so the f(dR)/C form is undefined for it. Use a "
                                     "parametric method (expo | polyu_fixedRp).");

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

        const int npt = h_fit_ok->GetNbinsX(), neta = h_fit_ok->GetNbinsY();
        cells.assign(npt, std::vector<Cell>(neta));

        // The RAW measured curves come from the SAME histogram file the fit stage read, projected
        // with the SAME per-cell projection + conditional errors (dr_correction_ratio.h). Opened
        // lazily -- a run in which every cell is fitted never needs it.
        TFile* fh = nullptr;
        TH3D *h3n = nullptr, *h3d = nullptr, *h3a = nullptr, *h3b = nullptr;
        // The pair-pT GROUPING of the fit cells, needed by the raw-bin fallback below: with the
        // last two bins merged, cell `iy` at the top spans TWO filled bins, and projecting bin
        // `iy` alone would silently deliver half of it. Built from the histograms' own axis, so it
        // is filled in together with them (dr_correction_pt_groups.h).
        DrPtGroups Gpt;
        const std::string hist_path = cfg.mc_dir + "mc_trig_eff_hists_" + cfg.mc_label
                                    + DrCorrWpSuffix(use_tight_wp)
                                    + MCTrigEffPairPt::FileSuffix() + "_step3.root";
        const std::string hp = "h_mc_dr_" + (sign.empty() ? std::string() : sign + "_")
                             + "zoom_vs_pt_eta_";

        int& n_fitted  = n_cells_fitted;
        int& n_fallback = n_cells_raw;
        int& n_dead    = n_cells_dead;
        for (int iy = 1; iy <= npt; ++iy) {
            for (int iz = 1; iz <= neta; ++iz) {
                Cell& c = cells[iy - 1][iz - 1];
                if (h_fit_ok->GetBinContent(iy, iz) > 0.5) {
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
                        if (DrCorrPlateauUsable(C, eC)) {
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
                        std::cout << "  ** cell (pair pT bin " << iy << ", pair eta bin " << iz
                                  << ") has fit_ok = 1 but an UNUSABLE fitted baseline C = "
                                  << C << " +- " << eC
                                  << " -> rejected, falling back to the raw bins" << std::endl;
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
                    Gpt = MakeDrPtGroups(h3n->GetYaxis(),
                                         DrCorrModeMergeLastTwoPt(plateau_mode));
                    if (Gpt.n != npt)
                        throw std::runtime_error("DrCorrectionEvaluator: the fit file has "
                            + std::to_string(npt) + " pair-pT cells but plateau mode '"
                            + plateau_mode + "' groups " + hist_path + " into "
                            + std::to_string(Gpt.n) + " -- fit file and mode disagree");
                }
                TH1D* r = DrGroupCellRatio(h3n, h3d, h3a, h3b, nullptr, nullptr, Gpt, iy, iz,
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
                  << "]: " << npt << "x" << neta << " cells -- " << n_fitted << " fitted, "
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
    double Eval(double dr, double pair_pt, double pair_eta)
    {
        ++n_eval;
        if (dr >= kDrMax) { ++n_flat; return 1.0; }

        const int iy = h_fit_ok->GetXaxis()->FindBin(pair_pt);
        const int iz = h_fit_ok->GetYaxis()->FindBin(pair_eta);
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

    void PrintStats() const
    {
        auto pct = [&](long long n) { return n_eval ? 100.0 * n / n_eval : 0.0; };
        std::cout << "DrCorrectionEvaluator [" << label << " / " << method << "]: "
                  << n_eval << " evaluations -- "
                  << n_flat << " at dR >= " << kDrMax << " (" << pct(n_flat) << "%, no correction), "
                  << n_fit  << " from a fit (" << pct(n_fit) << "%), "
                  << n_raw  << " from the RAW-BIN PLACEHOLDER (" << pct(n_raw) << "%), "
                  << n_rawempty << " in an empty bin of a placeholder cell (" << pct(n_rawempty)
                  << "%, no correction), "
                  << n_none << " with no correction available (" << pct(n_none) << "%); "
                  << n_floor << " floored at " << kMinCorr << ", " << n_cap << " capped at "
                  << kMaxCorr << ", " << n_outside << " outside the cell grid (must be 0)"
                  << std::endl;
    }
};

#endif  // DR_CORRECTION_APPLY_H
