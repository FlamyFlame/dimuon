// The PROPOSED generator-level loose dimuon mass filter of the additional-statistics MC request
// (docs/tracking/pthat_slice_mass_stats_sample_request.md). ONE definition, shared by the
// per-pT-hat-slice mass/statistics macro (mc_pthat_slice_mass_statistics.cxx) and the fullsim
// per-slice statistics plots (plot_pythia_fullsim_kn_pt_crossx.cxx, same-sign variant), so the
// two cannot drift apart on the value being argued about.
#pragma once
namespace MCRequest {
constexpr double kLooseMassMax = 10.0;   // GeV
}
