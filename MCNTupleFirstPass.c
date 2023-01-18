#ifndef __MCNTupleFirstPass_C__
#define __MCNTupleFirstPass_C__

#include "MCNTupleFirstPass.h"
#include "Riostream.h"
#include "TLorentzVector.h"
#include "time.h"
#include<bits/stdc++.h>

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

void MCNTupleFirstPass::UpdateCurParents(bool isFirstTime, bool isMuon1, std::vector<int>& cur_prt_bars, std::vector<int>& cur_prt_ids, int hf_quark_index = -1){
  
  // // print out if there are more than one hadron-level parents
  // bool prev_id_c_hadron = ((abs(prev_first_prt_id) >= 400 && abs(prev_first_prt_id) < 500) || (abs(prev_first_prt_id) >= 4000 && abs(prev_first_prt_id) < 5000));
  // bool prev_id_b_hadron = ((abs(prev_first_prt_id) >= 500 && abs(prev_first_prt_id) < 600) || (abs(prev_first_prt_id) >= 5000 && abs(prev_first_prt_id) < 6000));
  // bool cur_is_quark_gluon = (abs(cur_prt_ids[0]) <= 5 || cur_prt_ids[0] == 21);
  // check before updating so that the earliest hadron-level parents don't get printed out
  if (cur_prt_bars.size() > 1 && hf_quark_index < 0){
    if (isMuon1) m1_multi_hadronic_parents_ids = cur_prt_ids;
    else         m2_multi_hadronic_parents_ids = cur_prt_ids;
    // m_unspecified_parent_file << "Event#: " << mpair->m1.ev_num << std::endl;
    // m_unspecified_parent_file << "More than more parent at hadronic level:" << std::endl;
    // m_unspecified_parent_file << prev_first_prt_id << " <- ";
    // for (int parent_id : cur_prt_ids) m_unspecified_parent_file << parent_id << " ";
    // m_unspecified_parent_file << std::endl << std::endl;
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

  // record if there is oscillation
  if (cur_prt_ids[0] == (-1) * prev_first_prt_id){
    if (isMuon1) mpair->m1_osc = true;
    else         mpair->m2_osc = true;
  }

  //c-tag
  if ((abs(cur_prt_ids[0]) >= 400 && abs(cur_prt_ids[0]) < 500) || (abs(cur_prt_ids[0]) >= 4000 && abs(cur_prt_ids[0]) < 5000)){ // c-hadrons
    if (isMuon1) mpair->m1_c_tag = true;
    else         mpair->m2_c_tag = true;
  }

  // store the first parent id if isFirstTime
  if (isFirstTime){
    if (isMuon1)
      mpair->m1_parent_id = cur_prt_ids[0];
    else
      mpair->m2_parent_id = cur_prt_ids[0];

    h_numParents->Fill(cur_prt_bars.size() - 0.5);
  }
}


int MCNTupleFirstPass::FindHeavyQuarks(std::vector<int>& cur_prt_ids, int quark_type){
  if (quark_type != 4 && quark_type != 5){
    std::cout << "Error:: the parameter quark_type must take value of 4 (c) or 5 (b), quitting" << std::endl;
    throw std::exception();
  }

  if (cur_prt_ids.size() == 0) return -2; // case II

  std::vector<int>::iterator it_q = std::find(cur_prt_ids.begin(),cur_prt_ids.end(),quark_type);
  std::vector<int>::iterator it_qbar = std::find(cur_prt_ids.begin(),cur_prt_ids.end(), (-1) * quark_type);
  if (it_q != cur_prt_ids.end()){
    if (it_qbar != cur_prt_ids.end()){
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

  int ind = (isMuon1)? mpair->m1.ind : mpair->m2.ind;
  parent_bars.push_back(ind);
  parent_ids.push_back(13);
  UpdateCurParents(true,isMuon1,parent_bars,parent_ids);

  while(abs(parent_ids[0]) == 13 || abs(parent_ids[0]) == 15 || (abs(parent_ids[0]) >= 400 && abs(parent_ids[0]) < 500) || (abs(parent_ids[0]) >= 4000 && abs(parent_ids[0]) < 5000)){
    UpdateCurParents(false,isMuon1,parent_bars,parent_ids);
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
    while((abs(parent_ids[0]) >= 500 && abs(parent_ids[0]) < 600) || (abs(parent_ids[0]) >= 5000 && abs(parent_ids[0]) < 6000)){
      UpdateCurParents(false,isMuon1,parent_bars,parent_ids);
    }
  }

  int quark = (mc_mode == "mc_truth_bb")? 5:4;
  int quark_index = FindHeavyQuarks(parent_ids, quark);
  if (isMuon1) m1_history_before_hadrons.push_back(parent_ids);
  else         m2_history_before_hadrons.push_back(parent_ids);

  do{
    UpdateCurParents(false,isMuon1,parent_bars,parent_ids,quark_index); // need to update at least one step up anyway
    quark_index = FindHeavyQuarks(parent_ids, quark);
    if (isMuon1) m1_history_before_hadrons.push_back(parent_ids);
    else         m2_history_before_hadrons.push_back(parent_ids);
  }while (quark_index >= 0);

  // Remark: this does not always give us the correct result
  // But for the events we care about (which we record and plot)
  // (where both muons come from b or c for b-bbar or c-cbar sample, respectively)
  // this gives the correct result
  // the others (say if a muon comes from s/light hadron, 
  // or if a muon in the c-cbar sample comes from b initially)
  // won't give us any error: FindHeavyQuarks will just return 1 for not finding any HF quark
  // so the m*_ancestor_is_incoming, cur_m*_ancestor_ids, cur_m*_ancestor_bars are not informative
  // but we won't use them anyway

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


void MCNTupleFirstPass::FillMuonPairParents(){
  bool not_near = !(mpair->dphi < pms.PI / 2.);
  mpair->m1_c_tag = false;
  mpair->m2_c_tag = false;
  mpair->m1_osc = false;
  mpair->m2_osc = false;

  for (std::vector<int> v : m1_history_before_hadrons) v.clear();
  for (std::vector<int> v : m2_history_before_hadrons) v.clear();
  m1_history_before_hadrons.clear();
  m2_history_before_hadrons.clear();

  m1_multi_hf_quark_ids.clear();
  m2_multi_hf_quark_ids.clear();

  m1_multi_hadronic_parents_ids.clear();
  m2_multi_hadronic_parents_ids.clear();

  FindSingleMuonParents(true);
  FindSingleMuonParents(false);
  
  // parent groups: {direct b, c from b, c not from b, strange & light hadrons, photons};
  mpair->m1_parent_group = ParentGrouping(abs(mpair->m1_earliest_parent_id), mpair->m1_c_tag);
  mpair->m2_parent_group = ParentGrouping(abs(mpair->m2_earliest_parent_id), mpair->m2_c_tag);

  if (mpair->m1_parent_group < 0) m_unspecified_parent_file << "*** Ungrouped parent: " << mpair->m1_earliest_parent_id << " ***" << std::endl << std::endl;
  if (mpair->m2_parent_group < 0) m_unspecified_parent_file << "*** Ungrouped parent: " << mpair->m2_earliest_parent_id << " ***" << std::endl << std::endl;

  h_MuonPairParentGroups[!mpair->same_sign][not_near]->Fill(mpair->m1_parent_group - 0.5, mpair->m2_parent_group - 0.5);
  h_MuonPairParentGroups_weighted[!mpair->same_sign][not_near]->Fill(mpair->m1_parent_group - 0.5, mpair->m2_parent_group - 0.5, mpair->weight);


  // outprint the case where either leading-muon or subleading-muon have more than one hadronic-level parents
  // if ((m1_multi_hadronic_parents_ids.size() != 0 || m2_multi_hadronic_parents_ids.size() != 0) && (m1_multi_hadronic_parents_ids != m2_multi_hadronic_parents_ids)){
  if (m1_multi_hadronic_parents_ids.size() != 0 || m2_multi_hadronic_parents_ids.size() != 0){
    std::vector<int> multi_hadronic_parents_ids = (m1_multi_hadronic_parents_ids.size() != 0)? m1_multi_hadronic_parents_ids : m2_multi_hadronic_parents_ids;
    m_unspecified_parent_file << "Event#: " << mpair->m1.ev_num << std::endl;
    m_unspecified_parent_file << "More than more parents at hadronic level:" << std::endl;
    // m_unspecified_parent_file << prev_first_prt_id << " <- ";
    for (int parent_id : multi_hadronic_parents_ids) m_unspecified_parent_file << parent_id << " ";
    m_unspecified_parent_file << std::endl << std::endl;
  }


  // the case of opposite sign muon pair, where one muon is directly from b, the other is from b to c
  // check if both are initially from the same b
  bool bb_op_sign_one_b_one_btoc = (mc_mode == "mc_truth_bb" && !mpair->same_sign) && ((mpair->m1_parent_group == 1 && mpair->m2_parent_group == 2) || (mpair->m1_parent_group == 2 && mpair->m2_parent_group == 1));
  
  if (bb_op_sign_one_b_one_btoc){
    if (cur_m1_earliest_parent_barcode == cur_m2_earliest_parent_barcode){
      // case I: from same b
      h_bb_op_one_b_one_btoc[not_near]->Fill(0.5);
      h_weighted_bb_op_one_b_one_btoc[not_near]->Fill(0.5,mpair->weight);
    }else if (mpair->m1_c_tag || mpair->m2_c_tag){
      // case II: involve oscillation
      h_bb_op_one_b_one_btoc[not_near]->Fill(1.5);
      h_weighted_bb_op_one_b_one_btoc[not_near]->Fill(1.5,mpair->weight);
    }else{
      // case III: others
      h_bb_op_one_b_one_btoc[not_near]->Fill(2.5);
      h_weighted_bb_op_one_b_one_btoc[not_near]->Fill(2.5,mpair->weight);

      m_unspecified_parent_file << "Event#: " << mpair->m1.ev_num << std::endl;
      m_unspecified_parent_file << "bb op one b one b-to-c:" << std::endl;
      m_unspecified_parent_file << "Leading muon b from: ";
      for (int pid : cur_m1_ancestor_ids) m_unspecified_parent_file << pid << " ";
      m_unspecified_parent_file << std::endl << "Subleading muon b from: ";
      for (int pid : cur_m2_ancestor_ids) m_unspecified_parent_file << pid << " ";
      m_unspecified_parent_file << std::endl << std::endl;
    }
  }

  bool bb_op_sign_one_b_one_btoc_from_same_b = bb_op_sign_one_b_one_btoc && (cur_m1_earliest_parent_barcode == cur_m2_earliest_parent_barcode);
  // case of b-bbar sample & both muons are from b initially
  if (mc_mode == "mc_truth_bb" && (mpair->m1_parent_group == 1 || mpair->m1_parent_group == 2) && (mpair->m2_parent_group == 1 || mpair->m2_parent_group == 2) && (!bb_op_sign_one_b_one_btoc_from_same_b)){
    // first check if same parents 
    bool same_prts = (!m1_ancestor_is_incoming && !m2_ancestor_is_incoming && cur_m1_ancestor_bars == cur_m2_ancestor_bars);
    if (!same_prts && (m1_multi_hf_quark_ids.size() != 0 || m2_multi_hf_quark_ids.size() != 0)){
      std::vector<int> multi_hf_quark_ids = (m1_multi_hf_quark_ids.size() != 0)? m1_multi_hf_quark_ids : m2_multi_hf_quark_ids;
      m_unspecified_parent_file << "Event#: " << mpair->m1.ev_num << std::endl;
      m_unspecified_parent_file << "Unexpected. Both " << quark_type << " and -" << quark_type << " found." << std::endl;
      for (auto v : multi_hf_quark_ids) m_unspecified_parent_file << v << " ";
      m_unspecified_parent_file << std::endl << std::endl;
    }

    h_bb_both_from_b_same_prts[not_near]->Fill(!same_prts + 0.5); // 0.5 for same parents; 1.5 for different parents
    h_bb_both_from_b_ancestor[not_near]->Fill(m1_ancestor_is_incoming + 0.5, m2_ancestor_is_incoming + 0.5);
    // print out the parents for near/away side
    if (mpair->m1.ev_num < 200){ // only print out for the first 100 events
      m_b_parent_file[not_near] << "Event#: " << mpair->m1.ev_num << std::endl;
      if (same_prts){
        // m_b_parent_file[not_near] << "Both b's have the same ancestors: " << std::endl;
        // m_b_parent_file[not_near] << mpair->m1_earliest_parent_id << ", " << mpair->m2_earliest_parent_id << " <------ ";
        // for (int pid : cur_m1_ancestor_ids) m_b_parent_file[not_near] << pid << " ";
        // m_b_parent_file[not_near] << std::endl << std::endl;
        m_b_parent_file[not_near] << "Both b's have the same ancestors: " << std::endl;
        m_b_parent_file[not_near] << mpair->m1_earliest_parent_id << ", " << mpair->m2_earliest_parent_id;
        for (std::vector<int> v : m1_history_before_hadrons){
          m_b_parent_file[not_near] << " <--- ";
          for (int vv : v) m_b_parent_file[not_near] << vv << " ";
        }
        m_b_parent_file[not_near] << std::endl << std::endl;

      }
      else{
        if (m1_ancestor_is_incoming) m_b_parent_file[not_near] << "Leading muon ancestor is incoming." << std::endl;
        else{
          // m_b_parent_file[not_near] << "Leading muon b ancestor: " << mpair->m1_earliest_parent_id << " <------ ";
          // for (int pid : cur_m1_ancestor_ids) m_b_parent_file[not_near] << pid << " ";
          // m_b_parent_file[not_near] << std::endl;
          m_b_parent_file[not_near] << "Leading muon b ancestor: " << mpair->m1_earliest_parent_id;
          for (std::vector<int> v : m1_history_before_hadrons){
            m_b_parent_file[not_near] << " <--- ";
            for (int vv : v) m_b_parent_file[not_near] << vv << " ";
          }
          m_b_parent_file[not_near] << std::endl;
        }
        if (m2_ancestor_is_incoming) m_b_parent_file[not_near] << "Subleading muon ancestor is incoming." << std::endl;
        else{
          // m_b_parent_file[not_near] << "Subleading muon b ancestor: " << mpair->m2_earliest_parent_id << " <------ ";
          // for (int pid : cur_m2_ancestor_ids) m_b_parent_file[not_near] << pid << " ";
          // m_b_parent_file[not_near] << std::endl << std::endl;
          m_b_parent_file[not_near] << "Subleading muon b ancestor: " << mpair->m2_earliest_parent_id;
          for (std::vector<int> v : m2_history_before_hadrons){
            m_b_parent_file[not_near] << " <--- ";
            for (int vv : v) m_b_parent_file[not_near] << vv << " ";
          }
          m_b_parent_file[not_near] << std::endl << std::endl;
        }
      }
    }
  }


  // case of c-cbar same-sign muon pairs both from c (directly)
  if (mc_mode == "mc_truth_cc" && mpair->m1_parent_group == 3 && mpair->m2_parent_group == 3){
    bool same_prts = (!m1_ancestor_is_incoming && !m2_ancestor_is_incoming && cur_m1_ancestor_bars == cur_m2_ancestor_bars);
    h_cc_both_from_c_same_prts[not_near]->Fill(!same_prts + 0.5); // 0.5 for same parents; 1.5 for different parents
    h_cc_both_from_c_ancestor[!mpair->same_sign][not_near]->Fill(m1_ancestor_is_incoming + 0.5, m2_ancestor_is_incoming + 0.5);

    if (mpair->m1.ev_num < 200){ // only print out for the first 100 events
      m_c_parent_file[!mpair->same_sign][not_near] << "Event#: " << mpair->m1.ev_num << std::endl;
      if (same_prts){
        m_c_parent_file[!mpair->same_sign][not_near] << "Both c's have the same ancestors: " << std::endl << "hadrons";
        for (std::vector<int> v : m1_history_before_hadrons){
          m_c_parent_file[!mpair->same_sign][not_near] << " <--- ";
          for (int vv : v) m_c_parent_file[!mpair->same_sign][not_near] << vv << " ";
        }
        m_c_parent_file[!mpair->same_sign][not_near] << std::endl << std::endl;
      }
      else{
        if (m1_ancestor_is_incoming) m_c_parent_file[!mpair->same_sign][not_near] << "Leading muon ancestor is incoming." << std::endl;
        else{
          m_c_parent_file[!mpair->same_sign][not_near] << "Leading muon c ancestor: " << std::endl << "hadrons";
          for (std::vector<int> v : m1_history_before_hadrons){
            m_c_parent_file[!mpair->same_sign][not_near] << " <--- ";
            for (int vv : v) m_c_parent_file[!mpair->same_sign][not_near] << vv << " ";
          }
          m_c_parent_file[!mpair->same_sign][not_near] << std::endl;
        }
        if (m2_ancestor_is_incoming) m_c_parent_file[!mpair->same_sign][not_near] << "Subleading muon ancestor is incoming." << std::endl;
        else{
          m_c_parent_file[!mpair->same_sign][not_near] << "Subleading muon c ancestor: " << std::endl << "hadrons";
          for (std::vector<int> v : m2_history_before_hadrons){
            m_c_parent_file[!mpair->same_sign][not_near] << " <--- ";
            for (int vv : v) m_c_parent_file[!mpair->same_sign][not_near] << vv << " ";
          }
          m_c_parent_file[!mpair->same_sign][not_near] << std::endl << std::endl;
        }
      }
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

    if(jentry%100000==0) cout<<"Processing "<<jentry<<" event out of "<<nentries<<" events"<<std::endl;

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

      mpair->weight     = fabs(EventWeights->at(0)) * filter_effcy / nentries;
      mpair->minv       = truth_mupair_m->at(i)/1000.;
      mpair->pair_pt    = truth_mupair_pt->at(i)/1000.;
      mpair->pair_y     = truth_mupair_y->at(i);
      mpair->asym       = truth_mupair_asym->at(i);
      mpair->acop       = truth_mupair_acop->at(i);

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

      // NPairsAfter++;
      FillMuonPairTree(); // mode = 2 for MC truth
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
  if (mc_mode == "mc_truth_bb"){    
    m_b_parent_file[0].close();
    m_b_parent_file[1].close();
  }else{
    for (int isign = 0; isign < ParamsSet::nSigns; isign++){
      m_c_parent_file[isign][0].close();
      m_c_parent_file[isign][1].close();
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
        h_MuonPairParentGroups[ksign][lphi]->GetXaxis()->SetBinLabel(ibin+1,parentGroupLabels[ibin].c_str());
        h_MuonPairParentGroups[ksign][lphi]->GetYaxis()->SetBinLabel(ibin+1,parentGroupLabels[ibin].c_str());
      }
    }
  }

  if (mc_mode == "mc_truth_bb"){
    for (int ibin = 0; ibin < 3; ibin++){
      for (int lphi = 0; lphi < 2; lphi++){
        h_bb_op_one_b_one_btoc[lphi]->GetXaxis()->SetBinLabel(ibin+1,bb_op_one_b_one_btoc_labels[ibin].c_str());
        h_weighted_bb_op_one_b_one_btoc[lphi]->GetXaxis()->SetBinLabel(ibin+1,bb_op_one_b_one_btoc_labels[ibin].c_str());
      }
    }
    for (int ibin = 0; ibin < 2; ibin++){
      for (int lphi = 0; lphi < 2; lphi++){
        h_bb_both_from_b_same_prts[lphi]->GetXaxis()->SetBinLabel(ibin+1,samePrtsLabels[ibin].c_str());
        h_bb_both_from_b_ancestor[lphi]->GetXaxis()->SetBinLabel(ibin+1,bb_both_from_b_ancestor_labels[ibin].c_str());
        h_bb_both_from_b_ancestor[lphi]->GetYaxis()->SetBinLabel(ibin+1,bb_both_from_b_ancestor_labels[ibin].c_str());
      }
    }
  }
  else{
    for (int ibin = 0; ibin < 2; ibin++){
      for (int lphi = 0; lphi < 2; lphi++){
        h_cc_both_from_c_same_prts[lphi]->GetXaxis()->SetBinLabel(ibin+1,samePrtsLabels[ibin].c_str());
        for (int isign = 0; isign < ParamsSet::nSigns; isign++){
          h_cc_both_from_c_ancestor[isign][lphi]->GetXaxis()->SetBinLabel(ibin+1,cc_both_from_c_ancestor_labels[ibin].c_str());
          h_cc_both_from_c_ancestor[isign][lphi]->GetYaxis()->SetBinLabel(ibin+1,cc_both_from_c_ancestor_labels[ibin].c_str());
        }
      }
    }
  }

  h_cutAcceptance[0]->Scale(1./h_cutAcceptance[0]->GetBinContent(1));
  h_cutAcceptance[1]->Scale(1./h_cutAcceptance[1]->GetBinContent(1));

  m_outfile->Write();

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  std::cout << "CPU time used (in seconds) for mode " << mode << " is " << cpu_time_used << std::endl;

}

#endif
