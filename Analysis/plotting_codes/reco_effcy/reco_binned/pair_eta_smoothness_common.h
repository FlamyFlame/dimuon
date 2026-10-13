#ifndef PAIR_ETA_SMOOTHNESS_COMMON_H
#define PAIR_ETA_SMOOTHNESS_COMMON_H
// =================================================================================================
// Shared configuration of the eta^pair smoothness test (docs/tracking/pair_reco_eff_bin_by_bin_reco_binned.md
// §3c, §4.1): the FINE plot-only eta^pair binning and the output locations.
//
// FINE BINNING `pair_eta_fine_smoothness` = uniform 0.1 in [-2.2, 2.2], 44 bins (user 2026-10-11).
// It is NOT retyped: it IS ParamsSet::pair_eta_measured_bins (the existing opt-in measured-range
// axis, whose edges are exactly the crossx fine-axis edges inside |eta^pair| < 2.2), checked here
// against ParamsSet::pair_eta_fiducial_max and the approved width. Every coarse
// CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap edge is a fine edge (checked).
// The fine binning is ONLY for the eta^pair histograms; eps^reco is always evaluated on the
// canonical COARSE axis (§4.1).
// =================================================================================================
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../../MuonObjectsParamsAndHelpers/ParamsSet.h"
#include "../../../MuonObjectsParamsAndHelpers/FullSimSampleType.h"
#include "../../../Utilities/PairRecoEffRecoBinnedEvaluator.h"

inline std::vector<double> PairEtaFineSmoothnessBins()
{
    static const ParamsSet pms{};
    const std::vector<double>& e = pms.pair_eta_measured_bins;
    const double kApprovedWidth = 0.1;   // user-approved fine width (2026-10-11)
    if (std::fabs(e.front() + ParamsSet::pair_eta_fiducial_max) > 1e-6
        || std::fabs(e.back() - ParamsSet::pair_eta_fiducial_max) > 1e-6)
        throw std::runtime_error("pair_eta_fine_smoothness: endpoints differ from +-pair_eta_fiducial_max");
    for (size_t i = 1; i < e.size(); ++i)
        if (std::fabs(e[i] - e[i - 1] - kApprovedWidth) > 1e-9)
            throw std::runtime_error("pair_eta_fine_smoothness: bin width is not the approved 0.1");
    for (double c : PairRecoEffRecoBinnedEvaluator::EtaEdges()) {
        bool found = false;
        for (double f : e) if (std::fabs(f - c) < 1e-5) found = true;
        if (!found) throw std::runtime_error("pair_eta_fine_smoothness: coarse edge " + std::to_string(c)
                                             + " is not a fine edge");
    }
    return e;
}

inline std::string SmoothnessDataDir() { return "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/pair_eta_smoothness/"; }
inline std::string SmoothnessPlotDir() { return "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/single_b_analysis/pp24/pair_eta_smoothness/"; }
inline std::string WpSuffix(bool tight) { return tight ? "" : "_medium_wp"; }   // crossx convention
inline std::string SmoothnessDataFile(bool tight) { return SmoothnessDataDir() + "pair_eta_smoothness_pp24_2mu4" + WpSuffix(tight) + ".root"; }
inline std::string SmoothnessMCFile() { return SmoothnessDataDir() + "pair_eta_smoothness_mc_closure_pp24_fullsim_full.root"; }
inline std::string RecoBinnedProductFile()
{
    return FullSimRecoEffDir(FullSimSampleInputDir(FullSimSampleType::pp, /*is_test_sample=*/false))
         + "pair_reco_eff_reco_binned_pp24_full.root";
}

#endif
