#pragma once

#include <string>
#include <utility>
#include <vector>

using QEtaBinning = std::vector<std::pair<float, float>>;

std::string pairToSuffix(const std::pair<float, float>& p);

struct CommonEffcyConfig {
    // NOTE: all q_eta_proj_ranges_* below are for SINGLE-MUON TRIGGER EFFICIENCY
    // fitting using q×eta bins. Do NOT use these for pair eta binning.

    QEtaBinning q_eta_proj_ranges_fine_excl_gap = { // Run 3; |q×eta| < 2.2, gap excluded
        {-2.4f, -2.0f},
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

    QEtaBinning q_eta_proj_ranges_coarse_incl_gap = { // Run 3; |q×eta| < 2.2, gap included
        {-2.4f, -2.0f},
        {-2.0f, -1.5f},
        {-1.5f, -1.0f},
        {-1.0f, -0.5f},
        {-0.5f, 0.5f},
        {0.5f, 1.0f},
        {1.0f, 1.5f},
        {1.5f, 2.0f},
        {2.0f, 2.2f}
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
    // Last bin extends to 2.4: signal cuts are per-muon |q×eta| < 2.2, so pair
    // eta reaches 2.4. No run-year split needed.
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
