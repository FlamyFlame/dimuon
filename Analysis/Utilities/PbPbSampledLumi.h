#pragma once

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
// analysis_metadata.md (2023 excludes bad run 462964). 2024/2025 T_AA are 2023
// placeholders, but the luminosities here are the real per-year values.
// =============================================================================
inline double PbPbMu4SampledLumiNb(int run_year){
    switch (run_year % 2000){
        case 23: return 1.02426;
        case 24: return 1.59663;
        case 25: return 2.59933;
        default: return 0.;   // unknown year -> 0 weight (caller guards Sum L > 0)
    }
}
