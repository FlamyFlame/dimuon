#ifndef __muon_pair_enums_data_h__
#define __muon_pair_enums_data_h__

enum cuts_data{
  nocut, 
  // pass_trigger, // whether pass dimuon trigger for current event
  pass_trigger_match, // trigger match for current pair
  pass_muon_quality, // single-muon quality cut
  pass_muon_eta, 
  pass_muon_pt, 
  pass_muon_dP_overP, // single-muon dP over P cut to exclude backgrounds from s-and-light-hadron decay
  pass_photoprod,
  pass_resonance 
};
std::vector<std::string> cutLabels_data = {"no cut", "trigger match", "muon quality", "muon eta", "muon pT", "muon dPoverP", "photoproduction", "resonance"};
const int numCuts_data = cutLabels_data.size();

#endif