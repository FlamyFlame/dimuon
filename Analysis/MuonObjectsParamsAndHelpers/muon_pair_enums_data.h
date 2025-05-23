#ifndef __muon_pair_enums_data_h__
#define __muon_pair_enums_data_h__

#include <array>
#include <string>

enum class CutsCommon : int {
  nocut = 0,
  pass_trigger_match,
  pass_muon_quality,
  pass_muon_eta,
  pass_muon_pt,
  pass_muon_dP_overP,
  pass_muon_d0_z0,
  pass_muon_trk_charge,
  NumCommonCuts
};

enum class CutsPP : int {
  pass_resonance = static_cast<int>(CutsCommon::NumCommonCuts),
  nCuts_pp_data
};

enum class CutsPbPb : int {
  pass_photoprod = static_cast<int>(CutsCommon::NumCommonCuts),
  pass_resonance,
  nCuts_PbPb_data
};

inline const std::array<std::string, static_cast<int>(CutsPP::nCuts_pp_data)> cutLabelsPP = {
  "no cut", "trigger match", "muon quality", "muon eta", "muon pT",
  "dPoverP", "d0,z0", "muon+track charge", "resonance"
};

inline const std::array<std::string, static_cast<int>(CutsPbPb::nCuts_PbPb_data)> cutLabelsPbPb = {
  "no cut", "trigger match", "muon quality", "muon eta", "muon pT",
  "dPoverP", "d0,z0", "muon+track charge", "photoproduction", "resonance"
};

#endif
