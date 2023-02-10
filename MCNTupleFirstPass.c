#ifndef __MCNTupleFirstPass_C__
#define __MCNTupleFirstPass_C__

#include "MCNTupleFirstPass.h"
#include "Riostream.h"
#include "TLorentzVector.h"
#include "time.h"
#include<bits/stdc++.h>
#include <assert.h>

// Benefits of using enum: like labeling equations 
// and referring to them using the labels rather than numbers
// we do not need to remember which number corresponds to which cut
// also, if we add new cuts, we don't need to change the part of the code
// corresponding to the existing cuts
enum cuts_MC{
  nocut, 
  pass_muon_eta, 
  pass_muon_pt, 
  pass_resonance, 
  pass_photoprod
};

bool MCNTupleFirstPass::PassCuts(){
  //Apply ALL CUTS but for resonances

  // NO quality cuts on the MC truth muons
  // if((mpair->m1.quality&mpair->m2.quality&1  )==0) return false;//compair->m2ined muon
  // if((mpair->m1.quality&mpair->m2.quality&8  )==0) return false;//Medium muon
  // if((mpair->m1.quality&mpair->m2.quality&32 )==0) return false;//IDCuts
  //if((mpair->m1.quality&mpair->m2.quality&256)==0) return false;//MuonCuts

  if (fabs(mpair->m1.eta) > 2.4 || fabs(mpair->m2.eta) > 2.4) return false;
  
  // passed muon eta cut
  // if (mpair->same_sign) h_cutAcceptance[0]->Fill(pass_muon_eta + 0.5);
  // else h_cutAcceptance[1]->Fill(pass_muon_eta + 0.5);
  h_cutAcceptance[!mpair->same_sign]->Fill(pass_muon_eta + 0.5); // if same sign: fill the h_cutAcceptance[0] histogram; if opposite sign, fill the h_cutAcceptance[1] histogram

  if (mpair->m1.pt < 4 || mpair->m2.pt < 4) return false;
  
  // passed muon pt cut
  // if (mpair->same_sign) h_cutAcceptance[0]->Fill(pass_muon_pt + 0.5);
  // else h_cutAcceptance[1]->Fill(pass_muon_pt + 0.5);
  h_cutAcceptance[!mpair->same_sign]->Fill(pass_muon_pt + 0.5);

  // if( fabs(mpair->m1.dP_overP) > pms.deltaP_overP_thrsh || fabs(mpair->m2.dP_overP) > pms.deltaP_overP_thrsh ) pass = false;
 
  return true;
}

bool MCNTupleFirstPass::IsResonance(){
  // assuming we have already filled up the muon pair mpair
  if (mpair->same_sign) h_cutAcceptance[0]->Fill(pass_resonance + 0.5);
  else{ // opposite sign
    if (mpair->minv > pms.minv_upper) return true; // upper cut at 60 GeV

    // bool isresonance = false;
    
    for (array<float,2> ires : pms.minv_cuts){
      if (mpair->minv > ires[0] && mpair->minv < ires[1]){
        return true;
      }
    }
    
    h_cutAcceptance[1]->Fill(pass_resonance + 0.5);
  }
  
  return false;
}

bool MCNTupleFirstPass::IsPhotoProduction(){
  // A := |Delta PT| / (sum pT) < 0.05 && alpha := (pi-Dphi)/pi < 0.01
  // float A = (mpair->m1.pt - mpair->m2.pt) / (mpair->m1.pt + mpair->m2.pt);
  // assert (A >= 0);
  // float alpha = (pms.PI - fabs(mpair->dphi)) / pms.PI;
  // return (!(mpair->same_sign) && A < 0.05 && alpha < 0.01);
  // return (!(mpair->same_sign) && mpair->asym < 0.05 && mpair->acop < 0.01);
  
  if (mpair->same_sign) h_cutAcceptance[0]->Fill(pass_photoprod + 0.5);
  else if (mpair->asym < 0.05 && mpair->acop < 0.01) return true; // is photoproduction 
  else h_cutAcceptance[1]->Fill(pass_photoprod + 0.5);

  return false;
}

void MCNTupleFirstPass::FillSingleMuonTree(){
  muonOutTree->Fill();
}

int MCNTupleFirstPass::ParentGrouping(int parent_id, bool c_tag){

  assert(parent_id > 0);

  if ((parent_id >= 500 && parent_id < 600) || (parent_id >= 5000 && parent_id < 6000)){
    if (! c_tag) return 1; // direct b
    else         return 2; // b -> c
  }
  if (c_tag) return 3; // c not from b
  if ((parent_id >= 100 && parent_id < 400) || (parent_id >= 1000 && parent_id < 4000)) // light and s-hadrons
    return 4; // light and s-flavored hadrons
  if (parent_id == 22) // photons
    return 5;

  return -1; // others
}

int MCNTupleFirstPass::UpdateCurParents(bool isMuon1, std::vector<int>& cur_prt_bars, std::vector<int>& cur_prt_ids, bool record_history = true, int hf_quark_index = -1){
  
  // // print out if there are more than one hadron-level parents
  // check before updating so that the earliest hadron-level parents don't get printed out
  bool cur_is_quark_gluon = (abs(cur_prt_ids[0]) <= 5 || cur_prt_ids[0] == 21);
  // if (cur_prt_bars.size() > 1 && !(cur_is_quark_gluon)){
  if (cur_prt_bars.size() > 1 && hf_quark_index < 0){
    // LATER REMEMBER TO ADD [EXCLUDE 21 21 / 21 2]
    m_unspecified_parent_file << "Event#: " << mpair->m1.ev_num << std::endl;
    m_unspecified_parent_file << "More than one parent at hadronic level:" << std::endl;
    for (int parent_id : cur_prt_ids) m_unspecified_parent_file << parent_id << " ";
    m_unspecified_parent_file << std::endl << std::endl;
  }

  int prev_first_prt_id = cur_prt_ids[0];
  cur_prt_ids.clear();

  std::vector<int>::iterator itbar;
  if (hf_quark_index > 0)
    itbar = std::find(truth_barcode->begin(),truth_barcode->end(),cur_prt_bars[hf_quark_index]);
  else
    itbar = std::find(truth_barcode->begin(),truth_barcode->end(),cur_prt_bars[0]);
  if (itbar == truth_barcode->end()){ // found the barcode of muon1 among all truth particles
    std::cout << "Error:: Barcode of the current parent not found among all truth particles, quitting" << std::endl;
    throw std::exception();
  }
  cur_prt_bars = truth_parents->at(itbar - truth_barcode->begin());

  for (int parent_bar : cur_prt_bars){
    itbar = std::find(truth_barcode->begin(),truth_barcode->end(),parent_bar);
    if (itbar == truth_barcode->end()){
      std::cout << "Error:: Barcode of the current parent not found among all truth particles, quitting" << std::endl;
      throw std::exception();
    }
    // cur_prt_ids.push_back(abs(truth_id->at(itbar - truth_barcode->begin())) % 10000);
    cur_prt_ids.push_back(truth_id->at(itbar - truth_barcode->begin()) % 10000);
  }

  // save the vector of the updated parent ids into the muon history chain
  if (record_history){
    if (isMuon1) m1_history_before_hadrons.push_back(cur_prt_ids);
    else         m2_history_before_hadrons.push_back(cur_prt_ids);
  }

  // record if there is oscillation
  if (cur_prt_ids[0] == (-1) * prev_first_prt_id){
    if (isMuon1) mpair->m1_osc = true;
    else         mpair->m2_osc = true;
  }

  if (cur_prt_ids[0] == 441 || cur_prt_ids[0] == 443 || cur_prt_ids[0] == 445){
    if (isMuon1) m1_from_J_psi = true;
    else         m2_from_J_psi = true;
  }

  if (cur_prt_ids[0] == 551 || cur_prt_ids[0] == 553 || cur_prt_ids[0] == 555){
    if (isMuon1) m1_from_Upsilon = true;
    else         m2_from_Upsilon = true;
  }


  //c-tag
  if ((abs(cur_prt_ids[0]) >= 400 && abs(cur_prt_ids[0]) < 500) || (abs(cur_prt_ids[0]) >= 4000 && abs(cur_prt_ids[0]) < 5000)){ // c-hadrons
    if (isMuon1) mpair->m1_c_tag = true;
    else         mpair->m2_c_tag = true;
  }

  // if (isFirstTime) h_numParents->Fill(cur_prt_bars.size() - 0.5);

  bool prev_is_hf_hadron = ((abs(prev_first_prt_id) >= 400 && abs(prev_first_prt_id) < 600) || (abs(prev_first_prt_id) >= 4000 && abs(prev_first_prt_id) < 6000));
  if (prev_is_hf_hadron && ((abs(cur_prt_ids[0]) < 4000 && abs(cur_prt_ids[0]) % 100 <= 5) || cur_prt_ids[0] == 21)) // including
    return prev_first_prt_id;
  return 0;
}


int MCNTupleFirstPass::FindHeavyQuarks(std::vector<int>& cur_prt_ids, int quark_type, bool isMuon1, int hadron_child_id = 0){
  if (quark_type != 4 && quark_type != 5){
    std::cout << "Error:: the parameter quark_type must take value of 4 (c) or 5 (b), quitting" << std::endl;
    throw std::exception();
  }

  if (cur_prt_ids.size() == 0){
    std::cout << "Error:: parent list is empty, quitting" << std::endl;
    throw std::exception();
  }

  if (cur_prt_ids.size() == 1 && cur_prt_ids[0] == 2212)
    return -2; // case II

  std::vector<int>::iterator it_q = std::find(cur_prt_ids.begin(),cur_prt_ids.end(),quark_type);
  std::vector<int>::iterator it_qbar = std::find(cur_prt_ids.begin(),cur_prt_ids.end(), (-1) * quark_type);
  if (it_q != cur_prt_ids.end()){ // found q
    if (it_qbar != cur_prt_ids.end()){ // found qbar
      int sign_correction = (quark_type == 4)? +1 : -1; // b quark to hadron changes sign
      if (hadron_child_id * sign_correction > 0)
        return it_q - cur_prt_ids.begin();
      if (hadron_child_id * sign_correction < 0)
        return it_qbar - cur_prt_ids.begin();
      // not the immediate parent of a hadron (hadron_child_id == 0)
      if (isMuon1) m1_multi_hf_quark_ids = cur_prt_ids;
      else         m2_multi_hf_quark_ids = cur_prt_ids;
    }
    return it_q - cur_prt_ids.begin();
  }
  if (it_qbar != cur_prt_ids.end()){
    return it_qbar - cur_prt_ids.begin();
  }

  // if (it_q == cur_prt_ids.end() && it_qbar == cur_prt_ids.end()){ // not find either
  //   return -1; // case I: non-empty parents of the heavy quarks
  // }

  return -1;
}


void MCNTupleFirstPass::FindSingleMuonParents(bool isMuon1){
  std::vector<int> parent_bars;
  std::vector<int> parent_ids;
  int first_hadron_id = 0;

  int ind = (isMuon1)? mpair->m1.ind : mpair->m2.ind;
  int charge = (isMuon1)? mpair->m1.charge : mpair->m2.charge;
  parent_bars.push_back(ind);
  parent_ids.push_back(-13 * charge);
  if (isMuon1) m1_history_before_hadrons.push_back(parent_ids);
  else         m2_history_before_hadrons.push_back(parent_ids);

  while(abs(parent_ids[0]) == 13 || abs(parent_ids[0]) == 15 || (abs(parent_ids[0]) >= 400 && abs(parent_ids[0]) < 500) || (abs(parent_ids[0]) >= 4000 && abs(parent_ids[0]) < 5000)){
    first_hadron_id = UpdateCurParents(isMuon1,parent_bars,parent_ids);
  }

  if (isMuon1){
    mpair->m1_earliest_parent_id = parent_ids[0];
    cur_m1_earliest_parent_barcode = parent_bars[0];
  } 
  else{
    mpair->m2_earliest_parent_id = parent_ids[0];
    cur_m2_earliest_parent_barcode = parent_bars[0];
  }        

  if (mc_mode == "mc_truth_bb"){
    if (!((abs(parent_ids[0]) >= 500 && abs(parent_ids[0]) < 600) || (abs(parent_ids[0]) >= 5000 && abs(parent_ids[0]) < 6000)))
      return;
    while((abs(parent_ids[0]) >= 500 && abs(parent_ids[0]) < 600) || (abs(parent_ids[0]) >= 5000 && abs(parent_ids[0]) < 6000)){
      first_hadron_id = UpdateCurParents(isMuon1,parent_bars,parent_ids);
    }
  }else{ // cc
    bool first_hadron_is_c_hadron = ((abs(first_hadron_id) >= 400 && abs(first_hadron_id) < 500) || (abs(first_hadron_id) >= 4000 && abs(first_hadron_id) < 5000));
    if (!first_hadron_is_c_hadron)
      return;
  }

  int quark = (mc_mode == "mc_truth_bb")? 5:4;
  // if (mpair->m1.ev_num < 80) std::cout << first_hadron_id << std::endl;
  int quark_index = FindHeavyQuarks(parent_ids, quark, isMuon1, first_hadron_id);

  // bool alreadyFinished = false;
  // if ((isMuon1 && m1_from_J_psi) || (!isMuon1 && m2_from_J_psi)){
  //   if (parent_ids.size() == 2 && parent_ids[0] == 21 && parent_ids[1] == 21)
  //     alreadyFinished = true; // skip the following
  //   if (parent_ids.size() == 2 && ((parent_ids[0] == 21 && abs(parent_ids[1]) <= 3) || (parent_ids[1] == 21 && abs(parent_ids[0]) <= 3)))
  //     alreadyFinished = true;
  // }

  if ((isMuon1 && m1_from_J_psi) || (!isMuon1 && m2_from_J_psi)){
    while(parent_ids.size() == 1 && (abs(parent_ids[0]) == 3 || abs(parent_ids[0]) == 1003 || abs(parent_ids[0]) == 2003)){ // 9940003
      UpdateCurParents(isMuon1,parent_bars,parent_ids); // need to update at least one step up anyway
    }
  }

  if ((isMuon1 && m1_from_Upsilon) || (!isMuon1 && m2_from_Upsilon)){
    while(parent_ids.size() == 1 && (abs(parent_ids[0]) == 1103 || abs(parent_ids[0]) == 203)){ // 9940003
      UpdateCurParents(isMuon1,parent_bars,parent_ids); // need to update at least one step up anyway
    }
  }

  while(quark_index >= 0){
    UpdateCurParents(isMuon1,parent_bars,parent_ids,quark_index); // need to update at least one step up anyway
    quark_index = FindHeavyQuarks(parent_ids, quark, isMuon1);
  }

  if (isMuon1){
    m1_ancestor_is_incoming = (quark_index == -2); // if Case II: then is incoming
    cur_m1_ancestor_ids = parent_ids;
    cur_m1_ancestor_bars = parent_bars;
  }
  else{
    m2_ancestor_is_incoming = (quark_index == -2);
    cur_m2_ancestor_ids = parent_ids;
    cur_m2_ancestor_bars = parent_bars;
  }
}


int MCNTupleFirstPass::AncestorGrouping(std::vector<int>& ancestor_ids, bool sameprts){
  // incoming
  if (ancestor_ids.size() == 1 && ancestor_ids[0] == 2212)
    return 0;

  // gg
  if (ancestor_ids.size() == 2 && ancestor_ids[0] == 21 && ancestor_ids[1] == 21)
    return 1;
  
  // gq
  if(ancestor_ids.size() == 2 && 
    ((ancestor_ids[0] == 21 && (abs(ancestor_ids[1]) == 1 || abs(ancestor_ids[1]) == 2 || abs(ancestor_ids[1]) == 3)) || 
     (ancestor_ids[1] == 21 && (abs(ancestor_ids[0]) == 1 || abs(ancestor_ids[0]) == 2 || abs(ancestor_ids[0]) == 3))))
    return 2;
  if(mc_mode == "mc_truth_bb" && ancestor_ids.size() == 2 && 
    ((ancestor_ids[0] == 21 && abs(ancestor_ids[1]) == 4) ||
     (ancestor_ids[1] == 21 && abs(ancestor_ids[0]) == 4)))
    return 2;

  // single gluon
  if (ancestor_ids.size() == 1 && ancestor_ids[0] == 21)
    return 3;
  
  // q qbar
  if (ancestor_ids.size() == 2 && (ancestor_ids[0] == (-1) * ancestor_ids[1])){
    if (abs(ancestor_ids[0]) == 1 || abs(ancestor_ids[0]) == 2 || abs(ancestor_ids[0]) == 3)
      return 4;
    if (mc_mode == "mc_truth_bb" && abs(ancestor_ids[0]) == 4)
      return 4;
  }

  // others
  m_unspecified_parent_file << "Event#: " << mpair->m1.ev_num << std::endl;
  if (sameprts)
    m_unspecified_parent_file << "Same parents. Ancestor not in any group." << std::endl;
  else
    m_unspecified_parent_file << "Different parents. One ancestor not in any group." << std::endl;
  for (auto v : ancestor_ids) m_unspecified_parent_file << v << " ";
  m_unspecified_parent_file << std::endl << std::endl;
  return -1;
}


void MCNTupleFirstPass::HardScatteringAnalysis(std::vector<int>& ancestor_bars, std::vector<int>& ancestor_ids, bool sameprts = true){
  std::vector<int>::iterator itbar = std::find(truth_barcode->begin(),truth_barcode->end(),ancestor_bars[0]);
  if (itbar == truth_barcode->end()){
    std::cout << "Error:: Barcode of an ancestor not found among all truth particles, quitting" << std::endl;
    throw std::exception();
  }
  std::vector<int> hard_scatt_out_bars = truth_children->at(itbar - truth_barcode->begin());
  std::vector<int> hard_scatt_out_ids = {};
  for (int child_bar : hard_scatt_out_bars){
    itbar = std::find(truth_barcode->begin(),truth_barcode->end(),child_bar);
    if (itbar == truth_barcode->end()){
      std::cout << "Error:: Barcode of an outgoing particle of the hard scattering not found among all truth particles, quitting" << std::endl;
      throw std::exception();
    }
    hard_scatt_out_ids.push_back(truth_id->at(itbar - truth_barcode->begin()) % 10000);
  }

  h_bb_ss_near_num_hard_scatt_out->Fill(hard_scatt_out_ids.size() + 0.5);
}

void MCNTupleFirstPass::PrintHistory(std::ofstream& f, bool same_ancestors, bool state_same_ancestors = false){
  f << "Event #: " << mpair->m1.ev_num << std::endl;
  if (same_ancestors){
    if (state_same_ancestors) f << "Both b's have the same ancestors: " << std::endl;

    for (int vv : m1_history_before_hadrons[0]) f << vv << " ";
    for (std::vector<std::vector<int>>::iterator it = m1_history_before_hadrons.begin() + 1; it < m1_history_before_hadrons.end(); it++){
      f << "<--- ";
      for (int vv : *it){
        f << vv << " ";
      }
    }
    f << std::endl;

    for (int vv : m2_history_before_hadrons[0]) f << vv << " ";
    for (std::vector<std::vector<int>>::iterator it = m2_history_before_hadrons.begin() + 1; it < m2_history_before_hadrons.end(); it++){
      f << "<--- ";
      for (int vv : *it){
        f << vv << " ";
      }
    }
    f << std::endl << std::endl;
  }else{ // different ancestors

    f << "Leading muon: " << std::endl;
    for (int vv : m1_history_before_hadrons[0]) f << vv << " ";
    for (std::vector<std::vector<int>>::iterator it = m1_history_before_hadrons.begin() + 1; it < m1_history_before_hadrons.end(); it++){
      f << "<--- ";
      for (int vv : *it){
        f << vv << " ";
      }
    }
    f << std::endl;

    f << "Subleading muon: " << std::endl;
    for (int vv : m2_history_before_hadrons[0]) f << vv << " ";
    for (std::vector<std::vector<int>>::iterator it = m2_history_before_hadrons.begin() + 1; it < m2_history_before_hadrons.end(); it++){
      f << "<--- ";
      for (int vv : *it){
        f << vv << " ";
      }
    }
    f << std::endl << std::endl;
  }
}


void MCNTupleFirstPass::FillMuonPairParents(){
  bool not_near = !(abs(mpair->dphi) < pms.PI / 2.);
  mpair->m1_c_tag = false;
  mpair->m2_c_tag = false;
  mpair->m1_osc = false;
  mpair->m2_osc = false;
  m1_from_J_psi = false;
  m2_from_J_psi = false;
  m1_from_Upsilon = false;
  m2_from_Upsilon = false;

  for (std::vector<int> v : m1_history_before_hadrons) v.clear();
  for (std::vector<int> v : m2_history_before_hadrons) v.clear();
  m1_history_before_hadrons.clear();
  m2_history_before_hadrons.clear();

  m1_multi_hf_quark_ids.clear();
  m2_multi_hf_quark_ids.clear();

  // m1_multi_hadronic_parents_ids.clear();
  // m2_multi_hadronic_parents_ids.clear();

  FindSingleMuonParents(true);
  FindSingleMuonParents(false);
  
  // parent groups: {direct b, c from b, c not from b, strange & light hadrons, photons};
  mpair->m1_parent_group = ParentGrouping(abs(mpair->m1_earliest_parent_id), mpair->m1_c_tag);
  mpair->m2_parent_group = ParentGrouping(abs(mpair->m2_earliest_parent_id), mpair->m2_c_tag);

  if (mpair->m1_parent_group < 0) m_unspecified_parent_file << "*** Ungrouped parent: " << mpair->m1_earliest_parent_id << " ***" << std::endl << std::endl;
  if (mpair->m2_parent_group < 0) m_unspecified_parent_file << "*** Ungrouped parent: " << mpair->m2_earliest_parent_id << " ***" << std::endl << std::endl;

  h_unweighted_parent_groups[!mpair->same_sign][not_near]->Fill(mpair->m1_parent_group - 0.5, mpair->m2_parent_group - 0.5);
  h_parent_groups[!mpair->same_sign][not_near]->Fill(mpair->m1_parent_group - 0.5, mpair->m2_parent_group - 0.5, mpair->weight);


  // the case of opposite sign muon pair, where one muon is directly from b, the other is from b to c
  // check if both are initially from the same b
  bool one_b_one_btoc = ((mpair->m1_parent_group == 1 && mpair->m2_parent_group == 2) || (mpair->m1_parent_group == 2 && mpair->m2_parent_group == 1));
  bool bb_op_sign_one_b_one_btoc = (mc_mode == "mc_truth_bb" && (!mpair->same_sign) && one_b_one_btoc);
  
  if (bb_op_sign_one_b_one_btoc){
    if (cur_m1_earliest_parent_barcode == cur_m2_earliest_parent_barcode){
      // case I: from same b
      h_unweighted_bb_op_one_b_one_btoc[not_near]->Fill(0.5);
      h_bb_op_one_b_one_btoc[not_near]->Fill(0.5,mpair->weight);
    }else if (mpair->m1_osc || mpair->m2_osc){
      // case II: involve oscillation
      h_unweighted_bb_op_one_b_one_btoc[not_near]->Fill(1.5);
      h_bb_op_one_b_one_btoc[not_near]->Fill(1.5,mpair->weight);
    }else if (cur_m1_ancestor_bars == cur_m2_ancestor_bars){
      // case III: from different 
      h_unweighted_bb_op_one_b_one_btoc[not_near]->Fill(2.5);
      h_bb_op_one_b_one_btoc[not_near]->Fill(2.5,mpair->weight);
    }else{
      h_unweighted_bb_op_one_b_one_btoc[not_near]->Fill(3.5);
      h_bb_op_one_b_one_btoc[not_near]->Fill(3.5,mpair->weight);
    }

    m_unspecified_parent_file << "Event#: " << mpair->m1.ev_num << std::endl;
    m_unspecified_parent_file << "bb op one b one b-to-c:" << std::endl;
    m_unspecified_parent_file << "Leading muon b from: ";
    for (int pid : cur_m1_ancestor_ids) m_unspecified_parent_file << pid << " ";
    m_unspecified_parent_file << std::endl << "Subleading muon b from: ";
    for (int pid : cur_m2_ancestor_ids) m_unspecified_parent_file << pid << " ";
    m_unspecified_parent_file << std::endl << std::endl;
  }

  // define useful booleans
  bool both_from_b = ((mpair->m1_parent_group == 1 || mpair->m1_parent_group == 2) && (mpair->m2_parent_group == 1 || mpair->m2_parent_group == 2));
  bool both_from_c = (mpair->m1_parent_group == 3 && mpair->m2_parent_group == 3);
  bool same_ancestors = (cur_m1_ancestor_bars == cur_m2_ancestor_bars);
  bool near_side = (abs(mpair->dphi) < pms.PI / 2.);
  bool small_dphi = (abs(mpair->dphi) < 0.4);

  bool bb_op_sign_one_b_one_btoc_from_same_b = bb_op_sign_one_b_one_btoc && (cur_m1_earliest_parent_barcode == cur_m2_earliest_parent_barcode);
  mpair->from_same_b = bb_op_sign_one_b_one_btoc_from_same_b;


  if (mc_mode == "mc_truth_bb"){

    // near-side dphi plots
    if (near_side){
      if (mpair->same_sign){  // near-side & same-sign
        h_dphi_bb_ss_near->Fill(mpair->dphi, mpair->weight); // fill in dphi regardless of parents & ancestors
      }else{ // opposite sign
        h_dphi_bb_op_near->Fill(mpair->dphi, mpair->weight);
        if (mpair->m1_parent_group == 1 && mpair->m2_parent_group == 1) // both from direct b
          h_dphi_bb_op_near_both_from_b->Fill(mpair->dphi, mpair->weight);
        else if (bb_op_sign_one_b_one_btoc){
          h_dphi_bb_op_near_one_b_one_btoc->Fill(mpair->dphi, mpair->weight);
          if (cur_m1_earliest_parent_barcode == cur_m2_earliest_parent_barcode)
            h_dphi_bb_op_near_from_same_b->Fill(mpair->dphi, mpair->weight);
        }
      }
    }

    // analysis of muon pairs both from b & NOT from the same b
    if (both_from_b && (!bb_op_sign_one_b_one_btoc_from_same_b)){
      
      // print out unspecified cases where an ancestor vector contains both b and -b
      if (!same_ancestors && (m1_multi_hf_quark_ids.size() != 0 || m2_multi_hf_quark_ids.size() != 0)){
        std::vector<int> multi_hf_quark_ids = (m1_multi_hf_quark_ids.size() != 0)? m1_multi_hf_quark_ids : m2_multi_hf_quark_ids;
        m_unspecified_parent_file << "Event#: " << mpair->m1.ev_num << std::endl;
        m_unspecified_parent_file << "Unexpected. Both b and -b found." << std::endl;
        for (auto v : multi_hf_quark_ids) m_unspecified_parent_file << v << " ";
        m_unspecified_parent_file << std::endl << std::endl;
      }

      // fill in histograms
      h_unweighted_bb_both_from_b_same_ancestors[!mpair->same_sign][not_near]->Fill(!same_ancestors + 0.5); // 0.5 for same parents; 1.5 for different parents
      h_bb_both_from_b_same_ancestors[!mpair->same_sign][not_near]->Fill(!same_ancestors + 0.5, mpair->weight); // 0.5 for same parents; 1.5 for different parents
      if (same_ancestors){
        h_unweighted_bb_both_from_b_ancestor_sp[!mpair->same_sign][not_near]->Fill(AncestorGrouping(cur_m1_ancestor_ids, false) + 0.5);
        h_bb_both_from_b_ancestor_sp[!mpair->same_sign][not_near]->Fill(AncestorGrouping(cur_m1_ancestor_ids, false) + 0.5, mpair->weight);
      }else{
        h_unweighted_bb_both_from_b_ancestor_dp[!mpair->same_sign][not_near]->Fill(AncestorGrouping(cur_m1_ancestor_ids, false) + 0.5, AncestorGrouping(cur_m2_ancestor_ids, false) + 0.5);
        h_bb_both_from_b_ancestor_dp[!mpair->same_sign][not_near]->Fill(AncestorGrouping(cur_m1_ancestor_ids, false) + 0.5, AncestorGrouping(cur_m2_ancestor_ids, false) + 0.5, mpair->weight);
      }

      // print history
      if (print_prt_history && mpair->m1.ev_num < 1000){
      // if (print_prt_history){
        if (same_ancestors)  PrintHistory(m_b_parent_file[!mpair->same_sign][not_near], true, true);
        else            PrintHistory(m_b_parent_file[!mpair->same_sign][not_near], false, true);
      }

      // analysis of the [same sign, near side, same ancestors] case
      if (mpair->same_sign && near_side && same_ancestors){
        if (mpair->m1_osc && mpair->m2_osc){
          if (one_b_one_btoc)
            h_bb_ss_near_involv_osc->Fill(4.5,mpair->weight);
          else{
            h_bb_ss_near_involv_osc->Fill(5.5,mpair->weight);
            PrintHistory(m_bb_ss_near_file, true);
          }
        }else if ((mpair->m1_osc && (!mpair->m2_osc)) || (mpair->m2_osc && (!mpair->m1_osc))){
          if (one_b_one_btoc){
            h_bb_ss_near_involv_osc->Fill(2.5,mpair->weight);
            PrintHistory(m_bb_ss_near_file, true);
          }
          else
            h_bb_ss_near_involv_osc->Fill(3.5,mpair->weight);
        }else{
          if (one_b_one_btoc)
            h_bb_ss_near_involv_osc->Fill(0.5,mpair->weight);
          else{
            h_bb_ss_near_involv_osc->Fill(1.5,mpair->weight);
            PrintHistory(m_bb_ss_near_file, true);
          }
        }

        // find the outgoing particles from the hard scattering
        // record the number of hard-scattering products
        // and their kinematics (& correlation with the muon kinematics)
        HardScatteringAnalysis(cur_m1_ancestor_bars, cur_m1_ancestor_ids);
      }
    }
  }


  // case of c-cbar same-sign muon pairs both from c (directly)
  if (mc_mode == "mc_truth_cc" && both_from_c){
    
    // print out unspecified cases where an ancestor vector contains both c and -c
    if (!same_ancestors && (m1_multi_hf_quark_ids.size() != 0 || m2_multi_hf_quark_ids.size() != 0)){
      std::vector<int> multi_hf_quark_ids = (m1_multi_hf_quark_ids.size() != 0)? m1_multi_hf_quark_ids : m2_multi_hf_quark_ids;
      m_unspecified_parent_file << "Event#: " << mpair->m1.ev_num << std::endl;
      m_unspecified_parent_file << "Unexpected. Both c and -c found." << std::endl;
      for (auto v : multi_hf_quark_ids) m_unspecified_parent_file << v << " ";
      m_unspecified_parent_file << std::endl << std::endl;
    }

    // fill histograms
    h_unweighted_cc_both_from_c_same_ancestors[!mpair->same_sign][not_near]->Fill(!same_ancestors + 0.5); // 0.5 for same parents; 1.5 for different parents
    h_cc_both_from_c_same_ancestors[!mpair->same_sign][not_near]->Fill(!same_ancestors + 0.5, mpair->weight); // 0.5 for same parents; 1.5 for different parents
    if (same_ancestors){
      h_unweighted_cc_both_from_c_ancestor_sp[!mpair->same_sign][not_near]->Fill(AncestorGrouping(cur_m1_ancestor_ids, true) + 0.5);
      h_cc_both_from_c_ancestor_sp[!mpair->same_sign][not_near]->Fill(AncestorGrouping(cur_m1_ancestor_ids, true) + 0.5, mpair->weight);
    }else{
      h_unweighted_cc_both_from_c_ancestor_dp[!mpair->same_sign][not_near]->Fill(AncestorGrouping(cur_m1_ancestor_ids, false) + 0.5, AncestorGrouping(cur_m2_ancestor_ids, false) + 0.5);
      h_cc_both_from_c_ancestor_dp[!mpair->same_sign][not_near]->Fill(AncestorGrouping(cur_m1_ancestor_ids, false) + 0.5, AncestorGrouping(cur_m2_ancestor_ids, false) + 0.5, mpair->weight);
    }

    // print history
    if (print_prt_history && mpair->m1.ev_num < 1000){
    // if (print_prt_history){
      if (same_ancestors)  PrintHistory(m_c_parent_file[!mpair->same_sign][not_near], true, true);
      else            PrintHistory(m_b_parent_file[!mpair->same_sign][not_near], false, true);
    }
  }

  if (mc_mode == "mc_truth_cc" && mpair->same_sign){
    if (small_dphi){
      h_cc_ss_small_dphi_prt_gps->Fill(mpair->m1_parent_group - 0.5, mpair->m2_parent_group - 0.5,mpair->weight);
      h_cc_ss_small_dphi_same_ancestors->Fill(!same_ancestors + 0.5,mpair->weight);
      h_cc_ss_small_dphi_sp->Fill(AncestorGrouping(cur_m1_ancestor_ids, true) + 0.5,mpair->weight);
      h_cc_ss_small_dphi_dp->Fill(AncestorGrouping(cur_m1_ancestor_ids, false) + 0.5, AncestorGrouping(cur_m2_ancestor_ids, false) + 0.5,mpair->weight);
      m_cc_ss_small_dphi_file << mpair->m1.ev_num << "\t" << mpair->m1_parent_group << "\t" << mpair->m2_parent_group << std::endl;
    }else if(abs(mpair->dphi) > 0.7 && abs(mpair->dphi) < 1.2){
      h_cc_ss_plateau_prt_gps->Fill(mpair->m1_parent_group - 0.5, mpair->m2_parent_group - 0.5,mpair->weight);
      h_cc_ss_plateau_same_ancestors->Fill(!same_ancestors + 0.5,mpair->weight);
      h_cc_ss_plateau_sp->Fill(AncestorGrouping(cur_m1_ancestor_ids, true) + 0.5,mpair->weight);
      h_cc_ss_plateau_dp->Fill(AncestorGrouping(cur_m1_ancestor_ids, false) + 0.5, AncestorGrouping(cur_m2_ancestor_ids, false) + 0.5,mpair->weight);
    }
  }
}


void MCNTupleFirstPass::FillMuonPairTree(){
  int nsign = (mpair->same_sign)? 0:1;
  muonPairOutTree[nsign]->Fill();
  for (unsigned int idr = 0; idr < pms.ndRselcs; idr++){
    if (mpair->dr < pms.deltaR_thrsh[idr]){
      muonPairOutTreeBinned[idr][nsign]->Fill();
    }
  }
}


void MCNTupleFirstPass::ProcessData(){

  mpair = new MuonPairMC();
  // Long64_t nentries = 100000;
  Long64_t nentries = fChain->GetEntries();//number of events
  for (Long64_t jentry=0; jentry<nentries; jentry++) {//loop over the events
  // for (Long64_t jentry=0; jentry<5;jentry++) {//loop over the events

    if(jentry%10000==0) cout<<"Processing "<<jentry<<" event out of "<<nentries<<" events"<<std::endl;

    int num_bytes = fChain->GetEntry(jentry);//read in an event
    if(num_bytes==0){
      std::cout<<"Error:: Read in event has size of zero bytes,  quitting"<<std::endl;
      throw std::exception();
    }
 

    //trigger requirement for event
    // if(!b_HLT_2mu4) continue;

    // std::vector<int> muon_index_list = {};
    // std::vector<int>::iterator it;

    int NPairs = muon_pair_muon1_pt->size();//number of muon pairs in the event
    // h_numMuonPairsRaw->Fill(NPairs - 0.5);
    int NPairsAfter = 0;

    for(int i=0;i<NPairs;i++){//loop over all muon-pairs in the event

      // // use ind to record barcode instead
      mpair->m1.ind     = static_cast<int>(muon_pair_muon1_bar->at(i));
      mpair->m2.ind     = static_cast<int>(muon_pair_muon2_bar->at(i));

      mpair->m1.pt    = fabs(muon_pair_muon1_pt->at(i))/1000.0;//pt of the first muon in the pair
      mpair->m2.pt    = fabs(muon_pair_muon2_pt->at(i))/1000.0;//pt of the second muon in the pair
      mpair->m1.eta   = muon_pair_muon1_eta->at(i);
      mpair->m2.eta   = muon_pair_muon2_eta->at(i);
      mpair->m1.phi   = muon_pair_muon1_phi->at(i);
      mpair->m2.phi   = muon_pair_muon2_phi->at(i);
      mpair->m1.charge  = static_cast<int>(muon_pair_muon1_ch->at(i));//sign of pt stores charge
      mpair->m2.charge  = static_cast<int>(muon_pair_muon2_ch->at(i));//sign of pt stores charge
      // mpair->m1.quality = muon_pair_muon1_quality->at(i);
      // mpair->m2.quality = muon_pair_muon2_quality->at(i);
      // mpair->m1.dP_overP = muon_deltaP_overP->at(mpair->m1.ind);
      // mpair->m2.dP_overP = muon_deltaP_overP->at(mpair->m2.ind);
      mpair->m1.ev_num = jentry;
      mpair->m2.ev_num = jentry;

      mpair->weight     = EventWeights->at(0) * filter_effcy / nentries;
      mpair->minv       = truth_mupair_m->at(i)/1000.;
      mpair->pair_pt    = truth_mupair_pt->at(i)/1000.;
      mpair->pair_y     = truth_mupair_y->at(i);
      mpair->asym       = abs(truth_mupair_asym->at(i));
      mpair->acop       = abs(truth_mupair_acop->at(i));
      assert (mpair->asym >= 0 && mpair->acop >= 0);

      // if charge unequal gives 1 (opposite sign); otherwise gives 0 (same sign)
      h_cutAcceptance[mpair->m1.charge != mpair->m2.charge]->Fill(nocut + 0.5);
      // ------------------------------------------------------------

      // Apply cuts

      // Trigger match for muon pair
      // if(dimuon_b_HLT_2mu4->at(i)==false) continue;

      if (!PassCuts())continue;
    
      //------------------------------------------------------------

      //Two things at this step: 
      //1) sort pt, eta, phi by pt
      //2) update the muon-pair values
      // mpair->Update();
      mpair->UpdateShort();

      // resonance cut
      if (IsResonance()) continue;

      // photo-production cut
      if (IsPhotoProduction()) continue;
      
      //------------------------------------------------------------

      // fill in muon parents (now that we have finished applying all the cuts)
      FillMuonPairParents();

      //------------------------------------------------------------

      NPairsAfter++;
      FillMuonPairTree(); // mode = 2 for MC truth
      
      h_crossx->Fill(EventWeights->at(0));
      if(mpair->pt_lead < 8) // 4-8
        h_crossx_pt_binned[0]->Fill(EventWeights->at(0));
      else if(mpair->pt_lead < 15) // 8-15
        h_crossx_pt_binned[1]->Fill(EventWeights->at(0));
      else // > 15
        h_crossx_pt_binned[2]->Fill(EventWeights->at(0));
    }

    h_numMuonPairs->Fill(NPairsAfter - 0.5);
  }//loop over events

  delete mpair;
}


void MCNTupleFirstPass::Run(){
  clock_t start, end;
  double cpu_time_used;
  start = clock();

  if (mc_mode == "mc_truth_bb") filter_effcy = filter_effcy_bb;
  else if (mc_mode == "mc_truth_cc") filter_effcy = filter_effcy_cc;
  else{
      filter_effcy = -1000;
      std::cout << "MC mode has to be 'mc_truth_bb' or 'mc_truth_cc.' Else program will terminate." << std::endl;
  }
  std::cout << "The MC mode is " << mc_mode << ". Filter efficiency is " << filter_effcy << std::endl;

  // std::cout << "Mode = " << mode << ". Output file is " << m_outfile << std::endl;
  InitInput();
  InitOutput();
  ProcessData();
  m_unspecified_parent_file.close();

  if (mc_mode == "mc_truth_bb") m_bb_ss_near_file.close();

  if (print_prt_history){
    if (mc_mode == "mc_truth_bb"){
      for (int isign = 0; isign < ParamsSet::nSigns; isign++){
        m_b_parent_file[isign][0].close();
        m_b_parent_file[isign][1].close();
      }      
    }else{
      for (int isign = 0; isign < ParamsSet::nSigns; isign++){
        m_c_parent_file[isign][0].close();
        m_c_parent_file[isign][1].close();
        m_cc_ss_small_dphi_file.close();
      }
    }
  }

  // for loop over ibin: h->GetXaxis()->SetAxisLabel(ibin,label[ibin]);
  for (int ibin = 0; ibin < numCuts; ibin++){
    h_cutAcceptance[0]->GetXaxis()->SetBinLabel(ibin+1,cutLabels[ibin].c_str());
    h_cutAcceptance[1]->GetXaxis()->SetBinLabel(ibin+1,cutLabels[ibin].c_str());
  }

  for (int ibin = 0; ibin < nParentGroups; ibin++){
    for (int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
      for (int lphi = 0; lphi < 2; lphi++){
        h_unweighted_parent_groups[ksign][lphi]->GetXaxis()->SetBinLabel(ibin+1,parentGroupLabels[ibin].c_str());
        h_unweighted_parent_groups[ksign][lphi]->GetYaxis()->SetBinLabel(ibin+1,parentGroupLabels[ibin].c_str());
        h_parent_groups[ksign][lphi]->GetXaxis()->SetBinLabel(ibin+1,parentGroupLabels[ibin].c_str());
        h_parent_groups[ksign][lphi]->GetYaxis()->SetBinLabel(ibin+1,parentGroupLabels[ibin].c_str());
      }
    }
  }

  if (mc_mode == "mc_truth_bb"){
    for (int ibin = 0; ibin < 6; ibin++){
      h_bb_ss_near_involv_osc->GetXaxis()->SetBinLabel(ibin+1,osc_labels[ibin].c_str());
    }
    for (int ibin = 0; ibin < 3; ibin++){
      h_bb_ss_near_num_hard_scatt_out->GetXaxis()->SetBinLabel(ibin+1,num_hard_scatt_out_labels[ibin].c_str());
    }
    for (int lphi = 0; lphi < 2; lphi++){
      for (int ibin = 0; ibin < 4; ibin++){
        h_unweighted_bb_op_one_b_one_btoc[lphi]->GetXaxis()->SetBinLabel(ibin+1,bb_op_one_b_one_btoc_labels[ibin].c_str());
        h_bb_op_one_b_one_btoc[lphi]->GetXaxis()->SetBinLabel(ibin+1,bb_op_one_b_one_btoc_labels[ibin].c_str());
      }
      for (int isign = 0; isign < ParamsSet::nSigns; isign++){
        for (int ibin = 0; ibin < 2; ibin++){
          h_unweighted_bb_both_from_b_same_ancestors[isign][lphi]->GetXaxis()->SetBinLabel(ibin+1,samePrtsLabels[ibin].c_str());
          h_bb_both_from_b_same_ancestors[isign][lphi]->GetXaxis()->SetBinLabel(ibin+1,samePrtsLabels[ibin].c_str());
        }
        for (int ibin = 0; ibin < nAncestorGroups; ibin++){
          h_unweighted_bb_both_from_b_ancestor_sp[isign][lphi]->GetXaxis()->SetBinLabel(ibin+1,ancestor_labels[ibin].c_str());
          h_unweighted_bb_both_from_b_ancestor_dp[isign][lphi]->GetXaxis()->SetBinLabel(ibin+1,ancestor_labels[ibin].c_str());
          h_unweighted_bb_both_from_b_ancestor_dp[isign][lphi]->GetYaxis()->SetBinLabel(ibin+1,ancestor_labels[ibin].c_str());
          h_bb_both_from_b_ancestor_sp[isign][lphi]->GetXaxis()->SetBinLabel(ibin+1,ancestor_labels[ibin].c_str());
          h_bb_both_from_b_ancestor_dp[isign][lphi]->GetXaxis()->SetBinLabel(ibin+1,ancestor_labels[ibin].c_str());
          h_bb_both_from_b_ancestor_dp[isign][lphi]->GetYaxis()->SetBinLabel(ibin+1,ancestor_labels[ibin].c_str());
        }
      }
    }
  }
  else{
    for (int ibin = 0; ibin < nParentGroups; ibin++){
      h_cc_ss_small_dphi_prt_gps->GetXaxis()->SetBinLabel(ibin+1,parentGroupLabels[ibin].c_str());
      h_cc_ss_small_dphi_prt_gps->GetYaxis()->SetBinLabel(ibin+1,parentGroupLabels[ibin].c_str());
      h_cc_ss_plateau_prt_gps->GetXaxis()->SetBinLabel(ibin+1,parentGroupLabels[ibin].c_str());
      h_cc_ss_plateau_prt_gps->GetYaxis()->SetBinLabel(ibin+1,parentGroupLabels[ibin].c_str());
    }
    for (int ibin = 0; ibin < 2; ibin++){
      h_cc_ss_small_dphi_same_ancestors->GetXaxis()->SetBinLabel(ibin+1,samePrtsLabels[ibin].c_str());
      h_cc_ss_plateau_same_ancestors->GetXaxis()->SetBinLabel(ibin+1,samePrtsLabels[ibin].c_str());
    }
    for (int ibin = 0; ibin < nAncestorGroups; ibin++){
      h_cc_ss_small_dphi_sp->GetXaxis()->SetBinLabel(ibin+1,ancestor_labels[ibin].c_str());
      h_cc_ss_small_dphi_dp->GetXaxis()->SetBinLabel(ibin+1,ancestor_labels[ibin].c_str());
      h_cc_ss_small_dphi_dp->GetYaxis()->SetBinLabel(ibin+1,ancestor_labels[ibin].c_str());
      h_cc_ss_plateau_sp->GetXaxis()->SetBinLabel(ibin+1,ancestor_labels[ibin].c_str());
      h_cc_ss_plateau_dp->GetXaxis()->SetBinLabel(ibin+1,ancestor_labels[ibin].c_str());
      h_cc_ss_plateau_dp->GetYaxis()->SetBinLabel(ibin+1,ancestor_labels[ibin].c_str());
    }

    for (int lphi = 0; lphi < 2; lphi++){
      for (int isign = 0; isign < ParamsSet::nSigns; isign++){
        for (int ibin = 0; ibin < 2; ibin++){
          h_unweighted_cc_both_from_c_same_ancestors[isign][lphi]->GetXaxis()->SetBinLabel(ibin+1,samePrtsLabels[ibin].c_str());
          h_cc_both_from_c_same_ancestors[isign][lphi]->GetXaxis()->SetBinLabel(ibin+1,samePrtsLabels[ibin].c_str());
        }
        for (int ibin = 0; ibin < nAncestorGroups; ibin++){
          h_unweighted_cc_both_from_c_ancestor_sp[isign][lphi]->GetXaxis()->SetBinLabel(ibin+1,ancestor_labels[ibin].c_str());
          h_unweighted_cc_both_from_c_ancestor_dp[isign][lphi]->GetXaxis()->SetBinLabel(ibin+1,ancestor_labels[ibin].c_str());
          h_unweighted_cc_both_from_c_ancestor_dp[isign][lphi]->GetYaxis()->SetBinLabel(ibin+1,ancestor_labels[ibin].c_str());
          h_cc_both_from_c_ancestor_sp[isign][lphi]->GetXaxis()->SetBinLabel(ibin+1,ancestor_labels[ibin].c_str());
          h_cc_both_from_c_ancestor_dp[isign][lphi]->GetXaxis()->SetBinLabel(ibin+1,ancestor_labels[ibin].c_str());
          h_cc_both_from_c_ancestor_dp[isign][lphi]->GetYaxis()->SetBinLabel(ibin+1,ancestor_labels[ibin].c_str());
        }
      }
    }

    h_cc_ss_small_dphi_prt_gps->Scale(1./h_cc_ss_small_dphi_prt_gps->Integral());
    h_cc_ss_small_dphi_same_ancestors->Scale(1./h_cc_ss_small_dphi_same_ancestors->Integral());
    h_cc_ss_small_dphi_sp->Scale(1./h_cc_ss_small_dphi_sp->Integral());
    h_cc_ss_small_dphi_dp->Scale(1./h_cc_ss_small_dphi_dp->Integral());
    h_cc_ss_plateau_prt_gps->Scale(1./h_cc_ss_plateau_prt_gps->Integral());
    h_cc_ss_plateau_same_ancestors->Scale(1./h_cc_ss_plateau_same_ancestors->Integral());
    h_cc_ss_plateau_sp->Scale(1./h_cc_ss_plateau_sp->Integral());
    h_cc_ss_plateau_dp->Scale(1./h_cc_ss_plateau_dp->Integral());
  }

  h_cutAcceptance[0]->Scale(1./h_cutAcceptance[0]->GetBinContent(1));
  h_cutAcceptance[1]->Scale(1./h_cutAcceptance[1]->GetBinContent(1));

  m_outfile->Write();


  // h_numParents->Write();
  // h_numMuonPairs->Write();
  // if (mc_mode == "mc_truth_bb"){
  //   h_crossx->Write();
  //   for (int ipt = 0; ipt < npTbins; ipt++){
  //     h_crossx_pt_binned[ipt]->Write();
  //   }
  // }
  // for (unsigned int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
  //   muonPairOutTree[ksign]->Write();
  //   for (unsigned int idr = 0; idr < ParamsSet::ndRselcs; idr++){
  //     muonPairOutTreeBinned[idr][ksign]->Write();
  //   }
  // }

  // if (mc_mode == "mc_truth_bb"){
  //   h_unweighted_bb_op_one_b_one_btoc[0]->Write();
  //   h_unweighted_bb_op_one_b_one_btoc[1]->Write();
  //   h_bb_op_one_b_one_btoc[0]->Write();
  //   h_bb_op_one_b_one_btoc[1]->Write();
  //   for (int isign = 0; isign < ParamsSet::nSigns; isign++){
  //     h_unweighted_bb_both_from_b_same_ancestors[isign][0]->Write();
  //     h_unweighted_bb_both_from_b_same_ancestors[isign][1]->Write();
  //     h_bb_both_from_b_same_ancestors[isign][0]->Write();
  //     h_bb_both_from_b_same_ancestors[isign][1]->Write();
  //     h_unweighted_bb_both_from_b_ancestor_sp[isign][0]->Write();
  //     h_unweighted_bb_both_from_b_ancestor_sp[isign][1]->Write();
  //     h_unweighted_bb_both_from_b_ancestor_dp[isign][0]->Write();
  //     h_unweighted_bb_both_from_b_ancestor_dp[isign][1]->Write();
  //     h_bb_both_from_b_ancestor_sp[isign][0]->Write();
  //     h_bb_both_from_b_ancestor_sp[isign][1]->Write();
  //     h_bb_both_from_b_ancestor_dp[isign][0]->Write();
  //     h_bb_both_from_b_ancestor_dp[isign][1]->Write();
  //   }
  // }










  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  std::cout << "CPU time used (in seconds) for mode " << mode << " is " << cpu_time_used << std::endl;

}

#endif
