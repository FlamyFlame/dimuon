#pragma once

#include <string>
#include <utility>
#include <vector>

using QEtaBinning = std::vector<std::pair<float, float>>;

std::string pairToSuffix(const std::pair<float, float>& p);

struct CommonEffcyConfig {
    // NOTE: all q_eta_proj_ranges_* below are for SINGLE-MUON TRIGGER EFFICIENCY
    // fitting using q×eta bins. Do NOT use these for pair eta binning.

    QEtaBinning q_eta_proj_ranges_fine_excl_gap = { // Run 3; -2.4 <= q*eta < 2.2 (one-sided), gap excluded
        // (-2.4,-2.0) split into two (round-5 change #2): the forward-endcap anomaly
        // (mc_trigger_efficiency.md R3/R8/R10) has q*eta sub-structure inside this wide bin.
        // The fine 2D trig-eff axis (makeEtaTrigEffcyBinning) already has an edge at -2.2.
        {-2.4f, -2.2f},
        {-2.2f, -2.0f},
        {-2.0f, -1.6f},
        {-1.6f, -1.3f},
        {-0.9f, -0.5f},
        {-0.5f, -0.1f},
        {0.1f, 0.5f},
        {0.5f, 1.0f},
        {1.3f, 1.6f},
        {1.6f, 2.0f},
        {2.0f, 2.2f}
    };

    QEtaBinning q_eta_proj_ranges_fine_excl_gap_run2 = { // Run 2; gap excluded
        {-2.4f, -2.0f},
        {-2.0f, -1.6f},
        {-1.6f, -1.3f},
        {-0.9f, -0.5f},
        {-0.5f, -0.1f},
        {0.1f, 0.5f},
        {0.5f, 1.0f},
        {1.3f, 1.6f},
        {1.6f, 2.0f},
        {2.0f, 2.4f}
    };

    // RUN-3 COARSE q*eta binning -- the NOMINAL binning for the single-muon mu4 turn-on fits
    // since 2026-08-04 (user, advisor feedback). 10 contiguous bins spanning -2.4 <= q*eta < 2.3
    // (one-sided; the forward slice q*eta > 2.3 is removed by the fiducial gap cut,
    // ParamsSet::single_mu_fiducial_gap_cuts).
    //
    // WHY COARSE, AND WHY IT MATTERS: the fine `q_eta_proj_ranges_fine_excl_gap` binning EXCLUDES
    // the gap regions, so muons landing in a gap had no fitted turn-on and fell back to an
    // unfitted 2D (pT, q*eta) histogram ratio -- and to a -1.0 "no efficiency" sentinel when that
    // 2D bin was empty, which silently dropped the pair from every corrected histogram
    // (pp_trig_eff_highpt_jump.md). The hypothesis behind this change is that the same
    // binning/fluctuation effects in those gap regions are also what pushed the dR-correction
    // plateaus away from 1. A CONTIGUOUS binning has no holes, so there is no fallback and no
    // sentinel: every surviving muon has a fitted turn-on.
    //
    // The only difference from the pre-2026-08-04 version is that (-0.5, 0.5) is SPLIT at 0.
    // The split matters because q*eta folds the two charges through the toroid bending
    // direction, and the eta ~ 0 crack is not symmetric in q*eta (muon_gap_cuts_acceptance.md F7).
    // TOP EDGE TRACKS THE GAP CUT: it must equal the lower edge of the forward window in
    // ParamsSet::single_mu_fiducial_gap_cuts (2.30 since 2026-08-04). If they disagree, muons
    // between the two survive the cut with no fitted turn-on and the evaluator throws.
    QEtaBinning q_eta_proj_ranges_coarse_incl_gap = {
        {-2.4f, -2.0f},
        {-2.0f, -1.5f},
        {-1.5f, -1.0f},
        {-1.0f, -0.5f},
        {-0.5f, 0.0f},
        {0.0f, 0.5f},
        {0.5f, 1.0f},
        {1.0f, 1.5f},
        {1.5f, 2.0f},
        {2.0f, 2.3f}
    };

    QEtaBinning q_eta_proj_ranges_coarse_incl_gap_run2 = { // Run 2; gap included
        {-2.4f, -2.0f},
        {-2.0f, -1.5f},
        {-1.5f, -1.0f},
        {-1.0f, -0.5f},
        {-0.5f, 0.5f},
        {0.5f, 1.0f},
        {1.0f, 1.5f},
        {1.5f, 2.0f},
        {2.0f, 2.4f}
    };

    // For PAIR ETA binning (reco efficiency, cross-section, signal acceptance).
    // Last bin extends to 2.4: signal cuts are per-muon q*eta < 2.2 (one-sided,
    // -2.4 <= q*eta < 2.2), so pair eta reaches 2.4. No run-year split needed.
    QEtaBinning pair_eta_proj_ranges_coarse_incl_gap = {
        {-2.4f, -2.0f},
        {-2.0f, -1.5f},
        {-1.5f, -1.0f},
        {-1.0f, -0.5f},
        {-0.5f, 0.5f},
        {0.5f, 1.0f},
        {1.0f, 1.5f},
        {1.5f, 2.0f},
        {2.0f, 2.4f}
    };

};

inline void SetQEtaProjRanges(
    int run_year,
    QEtaBinning& q_eta_proj_ranges,
    std::vector<std::string>& q_eta_ranges_str,
    bool useCoarseQEtaBin = false)
{
    static const CommonEffcyConfig cfg{};

    q_eta_proj_ranges = (run_year > 20)
        ? (useCoarseQEtaBin ? cfg.q_eta_proj_ranges_coarse_incl_gap      : cfg.q_eta_proj_ranges_fine_excl_gap)
        : (useCoarseQEtaBin ? cfg.q_eta_proj_ranges_coarse_incl_gap_run2 : cfg.q_eta_proj_ranges_fine_excl_gap_run2);

    q_eta_ranges_str.clear();
    for (const auto& pair : q_eta_proj_ranges){
        q_eta_ranges_str.push_back("_q_eta_" + pairToSuffix(pair));
    }
}
