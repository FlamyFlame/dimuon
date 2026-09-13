#pragma once

#include <stdexcept>
#include <string>
#include <iostream>

// =============================================================================
// PbPb HLT_mu4 prescale-corrected SAMPLED LUMINOSITY per year [nb^-1].
//
// SINGLE SOURCE for the per-year luminosities used to luminosity-weight the
// combined-year Pb+Pb results. These MUST match the L_int factors baked into
// PbPbBaseClass::make_crossx_factors_pbpb_20YY (the crossx normalization).
//
// Year combination procedure (HF-muon R_AA internal note HION-2019-58, §4.1
// Eq.3): the combined Pb+Pb central value is the LUMINOSITY-WEIGHTED AVERAGE of
// the per-year yields,  combined = Sum_y( L_y * h_y ) / Sum_y( L_y ).
// Because each per-year histogram h_y is already weighted by crossx_factor_y
// (proportional to 1/L_y), this equals the correct combined normalization
// Sum_y(N_y) / (f * sigma * T_AA * Sum_y L_y) — i.e. total counts over total
// luminosity, NOT a naive sum of the per-year normalized histograms.
//
// Source of values: PbPbBaseClass.h crossx-factor functions / IntNotes
// analysis_metadata.md (2023 excludes the two b-hadron runs 461674 + 462964).
// 2024/2025/2026 T_AA are 2023 placeholders, but the luminosities here are the
// real per-year values.
// =============================================================================
inline double PbPbMu4SampledLumiNb(int run_year){
    switch (run_year % 2000){
        // 2023 corrected 2026-06-19: 1.02426 -> 1.17576 nb^-1. The old value used a
        // wrong/old GRL. Correct-GRL Prescale-Corrected Total = 1183.650457 ub^-1
        // (1.18365 nb^-1); the R_AA luminosity subtracts the two b-hadron runs
        // 461674 (2.623332) + 462964 (5.262407) ub^-1 => 1175.764718 ub^-1 =
        // 1.17576 nb^-1. (Old value only excluded 462964; 461674 was below the old
        // table range.) MUST match PbPbBaseClass.h make_crossx_factors_pbpb_2023.
        case 23: return 1.17576;
        // 2024 corrected 2026-06-19: 1.59663 -> 0.85112 nb^-1. Old value used a
        // stale GRL that did not exclude bad runs <489703; corrected GRL
        // (physics_HI2024_50ns.xml, runs >=489703). MUST match PbPbBaseClass.h
        // make_crossx_factors_pbpb_2024.
        case 24: return 0.85112;
        case 25: return 2.59933;
        // 2026 (added 2026-09-10): Prescale-Corrected Total = 2623.16 ub^-1 =
        // 2.62316 nb^-1 from lumitable_pbpb_26_HLT_mu4.csv (35 runs, 522041-523437),
        // whose run list is byte-identical to the skim GRL
        // physics_HI2026_50ns_noIBL.xml. No 2026 run is excluded at event level
        // (PbPbBadRuns has no 26 entry), so the R_AA luminosity is the GRL total.
        // MUST match PbPbBaseClass.h make_crossx_factors_pbpb_2026.
        case 26: {
            // ⚠ PROVISIONAL while the 2026 skim is still recovering (as of 2026-09-13).
            //
            // 2.62316 nb^-1 is the GRL total, i.e. the luminosity of ALL 35 runs.  It is the
            // right denominator ONLY once the skim covers all of them.  A skim bug
            // (TrigRates::ProcessZdc reading a ZDC aux item that only exists under
            // StoreZdc & 2) killed jobs on lumiblocks whose ZDC reco produced no RPD data, so
            // several runs are currently only PARTIALLY skimmed -- e.g. run 522200 had
            // 578/2192 files, run 522721 0/890.  The bug is fixed and recovery tasks are
            // running.
            //
            // Until every run reads 100 %, the numerator covers fewer events than this
            // denominator describes, so a 2026 cross-section computed now is biased LOW
            // (equivalently: the luminosity is too high for the yield it is divided into).
            // This cannot be detected downstream -- every histogram fills and every plot
            // renders -- so warn once, loudly, rather than returning the number silently.
            static bool warned = false;
            if (!warned) {
                warned = true;
                std::cerr << "\n*** WARNING: PbPb 2026 luminosity 2.62316 nb^-1 is the GRL total, "
                             "but the 2026 skim is INCOMPLETE (recovery in progress).\n"
                             "    Any 2026 luminosity-normalised result is PROVISIONAL and biased "
                             "low until every run is fully skimmed.\n"
                             "    See Analysis/docs/tracking/pbpb2026_analysis_support.md "
                             "(registry P14) before using it.***\n" << std::endl;
            }
            return 2.62316;
        }
        // Run 2 (2015 + 2018) is still a live code path: RDFBasedHistFillingPbPb keeps a
        // full `run_year == 15 || run_year == 18` branch and PbPbBaseClass registers
        // {18,"default"} crossx factors, so FillHistogramsCrossx can reach this function
        // with 15/18 (trigger_mode 2/3; mode 1 short-circuits via trigger_effcy_calc).
        // The value is the SAME combined 2015+2018 luminosity baked into
        // make_crossx_factors_pbpb_run2(): 0.436882 + 1.3803 = 1.817182 nb^-1.
        // Both years return it because the Run-2 crossx factors are the COMBINED ones
        // (see the "2015 + 2018 are combined" comment in PbPbBaseClass::BuildPbPbMaps).
        case 15:
        case 18: return 1.817182;
        // Throw, never return 0.  A zero here was silently catastrophic in two ways:
        // RDFBasedHistFillingPbPb computes 1.0/L, so 0 gave an INFINITE differential
        // cross-section weight; and in the luminosity-weighted year combination
        // Sum_y(L_y h_y)/Sum_y(L_y) a new year would be dropped at weight zero while
        // every plot still rendered with the year in its legend.
        default:
            throw std::runtime_error(
                "PbPbMu4SampledLumiNb: no sampled luminosity for Pb+Pb run year 20" +
                std::to_string(run_year % 2000) +
                " (known: 2015, 2018, 2023, 2024, 2025, 2026). Add it here AND in "
                "PbPbBaseClass.h::make_crossx_factors_pbpb_<yr>().");
    }
}
