#pragma once

#include <string>
#include <utility>
#include <vector>

using QEtaBinning = std::vector<std::pair<float, float>>;

std::string pairToSuffix(const std::pair<float, float>& p);

struct CommonEffcyConfig {
    // maps of q_eta bins to pT projection ranges serving single-muon effciency fitting
    QEtaBinning q_eta_proj_ranges_fine_excl_gap = {
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

    QEtaBinning q_eta_proj_ranges_fine_excl_gap_run2 = {
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

    QEtaBinning q_eta_proj_ranges_coarse_incl_gap = { // coarse bins including gaps
        {-2.4f, -2.0f}, 
        {-2.0f, -1.5}, 
        {-1.5, -1.0f}, 
        {-1.0f, -0.5f}, 
        {-0.5f, 0.5f}, 
        {0.5f, 1.0f}, 
        {1.0f, 1.5f}, 
        {1.5f, 2.0f}, 
        {2.0f, 2.2f}
    };

    QEtaBinning q_eta_proj_ranges_coarse_incl_gap_run2 = { // including [2.2, 2.4]
        {-2.4f, -2.0f}, 
        {-2.0f, -1.5}, 
        {-1.5, -1.0f}, 
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
