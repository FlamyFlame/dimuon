#pragma once

enum cuts_MC{
  nocut, 
  pass_muon_eta, 
  pass_muon_pt, 
  // pass_photoprod,
  pass_resonance 
};

// enum class CutsCommon : int {
//   nocut = 0,
//   pass_trigger_match,
//   pass_muon_quality,
//   pass_muon_eta,
//   pass_muon_pt,
//   pass_muon_dP_overP,
//   pass_muon_d0_z0,
//   pass_muon_trk_charge,
//   NumCommonCuts
// };


std::vector<std::string> cutLabels_MC = {"no cut", "muon eta", "muon pT", "resonance"};

enum single_muon_parent_group{ // parent groups for a single muon; used for pythia + powheg
  resonance_decay,
  direct_b, // excluding J/Psi
  b_to_c,
  direct_c,
  s_light,
  single_photon,
  prt_drell_yan
};

enum hard_scatter_categories{
  gg_gg,
  qqbar_gg,
  gq_gq,
  flavor_creat, // 0 to 2
  flavor_excit, // 1 to 1
  hq_scatt, // 2 to 2
  gg_qqbar,
  qqprime_qqprime,
  hard_scatt_drell_yan // never to be filled for HF muon pairs
};

enum pair_flavor_index{
  from_resonance = 0,
  resonance_contaminated = 1,
  from_single_b = 2,
  bb = 3,
  cc = 4,
  one_b_one_c = 5,
  photon_to_dimuon_splitting = 6, // photon to dimuon splitting, NOT photon to QQbar splitting + semi-leptonic decays
  drell_yan = 7,
  other_flavor = 8, // should just be combinatoric pairs
  nFlavors = 9
};

enum muon_pair_both_from_open_HF_origin_catgr{ // will only be applied to [both from b] or [both fromc c] cases
  fc = 0,
  same_gs_fsr = 1, // including photon splitting
  same_phs_fsr = 2, // including photon splitting
  same_gs_isr_zero_hard_scatt = 3, // both semi-hard pT
  same_gs_isr_one_hard_scatt = 4, // one FE, one semi-hard pT
  diff_gs_same_hard_scatt = 5,
  others = 6, // should just be combinatorid + be uncorrelated - plot Dphi & print out to check
  nOrigins = 7, // put this in front of non-HF --> do NOT fill hists for non-HF
  not_both_from_open_and_different_HF_hadrons = 8 // pairs from single-b / resonances / resonance-contamined also goes here
};

enum powheg_ancestor_categories{
  gg,
  gq,
  single_gluon,
  qqbar,
  incoming
};

enum powheg_origin_categories{
  gg_gQQbar,
  qg_qQQbar,
  qqbar_gQQbar,
  not_in_powheg
};
