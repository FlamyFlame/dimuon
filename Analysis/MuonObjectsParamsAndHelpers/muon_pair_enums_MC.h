#ifndef __muon_pair_enums_MC_h__
#define __muon_pair_enums_MC_h__


enum cuts_MC{
  nocut, 
  pass_muon_eta, 
  pass_muon_pt, 
  // pass_photoprod,
  pass_resonance 
};

std::vector<std::string> cutLabels_MC = {"no cut", "muon eta", "muon pT", "resonance"};
const int numCuts_MC = cutLabels_MC.size();

enum parent_groups{ // parent groups for a single muon; used for pythia + powheg
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
  from_resonance,
  from_single_b,
  bb,
  cc,
  other_flavor
};

enum muon_pair_both_from_open_b_or_c_origin_categories{ // will only be applied to [both from b] or [both fromc c] cases
  fc,
  same_gs_phs_fsr, // including photon splitting
  same_gs_isr_zero_hard_scatt, // both semi-hard pT
  same_gs_isr_one_hard_scatt, // one FE, one semi-hard pT
  same_gs_isr_both_hard_scatt, // different hard scatt, should be uncorrelated
  diff_gs_same_hard_scatt,
  // diff_gs_same_hard_scatt_XX_to_gg, // both from FSR gs
  // diff_gs_same_hard_scatt_FE, // one from FSR gs, the other from ISR gs / is incoming
  // diff_gs_same_hard_scatt_QQprime, // both are either from ISR gs or incoming
  others, // should be uncorrelated - plot Dphi & print out to check
  drell_yan // is a muon pair origin category, but never filled for HF muon pairs
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

#endif