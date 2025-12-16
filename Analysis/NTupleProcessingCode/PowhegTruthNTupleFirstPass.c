#include "PowhegNTupleFirstPass.c"

void PowhegTruthNTupleFirstPass::ProcessData(){
    int quark = (mc_mode == "bb")? 5:4;
    qqpair = new TruthQQPair(quark);

    PowhegNTupleFirstPass::ProcessData();
    
    delete qqpair;

    std::cout << "The total integral of events where the hardest scattering is 2-to-2 is: " << crossx_2_to_2 << std::endl;
    std::cout << "The total integral of events where the hardest scattering is 2-to-3 is: " << crossx_2_to_3 << std::endl;
    std::cout << "Total percentage in crossx of events where the relevant hard scattering for a same-ancestor muon pair is not the hardest scattering is: " << crossx_relevant_hard_isnt_hardest / (crossx_2_to_2 + crossx_2_to_3) << std::endl;
    std::cout << "Total percentage in crossx of skipped HF events (not filled in HF-muon-pair-ancestor-category trees): " << skipped_event_crossx / (crossx_2_to_2 + crossx_2_to_3) << std::endl;
}

void PowhegTruthNTupleFirstPass::InitTempVariables(){
    single_gluon_history = new std::vector<std::vector<int>>();
    m1_history = new std::vector<std::vector<int>>();
    m2_history = new std::vector<std::vector<int>>();
    m1_history_particle = new std::vector<std::vector<Particle>>();
    m2_history_particle = new std::vector<std::vector<Particle>>();
    // m1_last_hf_hadron_prt_pt_eta_phi_m = new std::vector<float>();
    // m2_last_hf_hadron_prt_pt_eta_phi_m = new std::vector<float>();
    // m1_last_b_hadron_prt_pt_eta_phi_m = new std::vector<float>();
    // m2_last_b_hadron_prt_pt_eta_phi_m = new std::vector<float>();
    m1_first_hf_hadron_prt_pt_eta_phi_m = new std::vector<float>();
    m2_first_hf_hadron_prt_pt_eta_phi_m = new std::vector<float>();
    // m1_hq_ancestor_pt_eta_phi_m = new std::vector<float>();
    // m2_hq_ancestor_pt_eta_phi_m = new std::vector<float>();
}

void PowhegTruthNTupleFirstPass::InitOutput(){
    InitOutStreamFiles();
    PowhegNTupleFirstPass::InitOutput();
}


void PowhegTruthNTupleFirstPass::InitOutStreamFiles(){
    // ---------------------------------------------------------------------------------------------------------------------------

    if (print_specific_prt_history){
        if (mc_mode == "cc"){

            m_cc_ss_small_dphi_file = new std::ofstream(Form("%scc_ss_small_dphi.txt", mcdir.c_str()));
            *m_cc_ss_small_dphi_file << "Event#\tm1-grp\tm2-grp" << std::endl;
        }else{
            m_bb_ss_near_file = new std::ofstream(Form("%sbb_ss_near.txt", mcdir.c_str()));
            m_bb_ss_away_file = new std::ofstream(Form("%sbb_ss_away.txt", mcdir.c_str()));
            m_bb_op_near_one_b_one_btoc_others_file = new std::ofstream(Form("%sbb_op_near_one_b_one_btoc_others.txt", mcdir.c_str()));
        }
    }

    // ---------------------------------------------------------------------------------------------------------------------------

    if (mc_mode == "bb"){
        m_unspecified_parent_file = new std::ofstream(mcdir + "unspecified_parents_bb.txt");
        if (print_prt_history){
            for (int isign = 0; isign < ParamsSet::nSigns; isign++){
                for (int jdphi = 0; jdphi < 2; jdphi++){
                    m_b_parent_file[isign][jdphi] = new std::ofstream(Form("%sb_parents_%s%s.txt", mcdir.c_str(), sign_labels[isign].c_str(), dphis[jdphi].c_str()));
                }
            }
        }
    }else{
        m_unspecified_parent_file = new std::ofstream(mcdir + "unspecified_parents_cc.txt");
        if (print_prt_history){
            for (int isign = 0; isign < ParamsSet::nSigns; isign++){
                for (int jdphi = 0; jdphi < 2; jdphi++){
                    m_c_parent_file[isign][jdphi] = new std::ofstream(Form("%sc_parents_%s%s.txt", mcdir.c_str(), sign_labels[isign].c_str(), dphis[jdphi].c_str()));
                }
            }
        }
    } 
}

void PowhegTruthNTupleFirstPass::InitOutputTrees(){

    PowhegNTupleFirstPass::InitOutputTrees();
    
    if (output_QQpair_tree){
        for (unsigned int isign = 0; isign < ParamsSet::nSigns; isign++){
            for (int jdphi = 0; jdphi < 2; jdphi++){
                for (int kgrp = 0; kgrp < nAncestorGroups; kgrp++){
                    QQPairOutTree[isign][jdphi][kgrp] = new TTree(Form("QQ_pair_tree%s%s%s",signs[isign].c_str(),dphis[jdphi].c_str(),ancestor_grps[kgrp].c_str()),Form("QQ pairs, muon %s%s%s",signs[isign].c_str(),dphis[jdphi].c_str(),ancestor_grps[kgrp].c_str()));
                    QQPairOutTree[isign][jdphi][kgrp]->Branch("QQPairObj",&qqpair);
                }
            }
        }
    }
}

void PowhegTruthNTupleFirstPass::InitOutputHists(){
    PowhegNTupleFirstPass::InitOutputHists();

    if (output_truth_hists){

        for (int isign = 0; isign < ParamsSet::nSigns; isign++){
            for (int jdphi = 0; jdphi < 2; jdphi++){
                h_parent_groups[isign][jdphi] = new TH2D(Form("h_parent_groups_sign%d%s",isign+1, dphis[jdphi].c_str()),Form("h_parent_groups_sign%d%s",isign+1, dphis[jdphi].c_str()),nParentGroups,0,nParentGroups,nParentGroups,0,nParentGroups);
                h_num_hard_scatt_out[isign][jdphi] = new TH1D(Form("h_num_hard_scatt_out_sign%d%s",isign+1, dphis[jdphi].c_str()),Form("h_num_hard_scatt_out_sign%d%s",isign+1, dphis[jdphi].c_str()),3,2,5);
                h_pt_muon_pt_closest_hadr_ratio[isign][jdphi] = new TH1D(Form("h_pt_muon_pt_closest_hadr_ratio_sign%d%s",isign+1, dphis[jdphi].c_str()),Form("h_pt_muon_pt_closest_hadr_ratio_sign%d%s",isign+1, dphis[jdphi].c_str()),40,0,1.);
                h_pt_closest_hadr_pt_furthest_hadr_ratio[isign][jdphi] = new TH1D(Form("h_pt_closest_hadr_pt_furthest_hadr_ratio_sign%d%s",isign+1, dphis[jdphi].c_str()),Form("h_pt_closest_hadr_pt_furthest_hadr_ratio_sign%d%s",isign+1, dphis[jdphi].c_str()),40,0,1.);
                h_pt_hadr_hq_ratio[isign][jdphi] = new TH1D(Form("h_pt_hadr_hq_ratio_sign%d%s",isign+1, dphis[jdphi].c_str()),Form("h_pt_hadr_hq_ratio_sign%d%s",isign+1, dphis[jdphi].c_str()),40,0,1.);
                h_dphi_muon_closest_hadr[isign][jdphi] = new TH1D(Form("h_dphi_muon_closest_hadr_sign%d%s",isign+1, dphis[jdphi].c_str()),Form("h_dphi_muon_closest_hadr_sign%d%s",isign+1, dphis[jdphi].c_str()),32,-pms.PI,pms.PI);
                
                for (int kgrp = 0; kgrp < nAncestorGroups; kgrp++){
                    h_QQ_DR[isign][jdphi][kgrp] = new TH1D(Form("h_QQ_DR_sign%d%s%s",isign+1, dphis[jdphi].c_str(), ancestor_grps[kgrp].c_str()),";#Delta R;d#sigma/d#Delta R", 50,0,5.75);
                    h_QQ_Dphi[isign][jdphi][kgrp] = new TH1D(Form("h_QQ_Dphi_sign%d%s%s",isign+1, dphis[jdphi].c_str(), ancestor_grps[kgrp].c_str()),";#Delta#phi;d#sigma/d#Delta#phi", 32,-pms.PI,pms.PI);
                    h_QQ_minv[isign][jdphi][kgrp] = new TH1D(Form("h_QQ_minv_sign%d%s%s",isign+1, dphis[jdphi].c_str(), ancestor_grps[kgrp].c_str()),";m_{QQ} [GeV]; d#sigma/dm_{QQ}",pms.n_hq_minv_bins,pms.hq_minvBins);
                    h_QQ_pair_pt_ptlead_ratio[isign][jdphi][kgrp] = new TH1D(Form("h_QQ_pair_pt_ptlead_ratio_sign%d%s%s",isign+1, dphis[jdphi].c_str(), ancestor_grps[kgrp].c_str()),";#frac{p_{T}^{pair}}{p_{T}^{lead}};d#sigma/d#frac{p_{T}^{pair}}{p_{T}^{lead}}", 25,0,2.);
                    h_QQ_pt_avg[isign][jdphi][kgrp] = new TH1D(Form("h_QQ_pt_avg_sign%d%s%s",isign+1, dphis[jdphi].c_str(), ancestor_grps[kgrp].c_str()),";p_{T}^{avg};d#sigma/dp_{T}^{avg}", pms.n_hq_pt_bins,pms.hq_pTBins);
                    h_QQ_asym[isign][jdphi][kgrp] = new TH1D(Form("h_QQ_asym_sign%d%s%s",isign+1, dphis[jdphi].c_str(), ancestor_grps[kgrp].c_str()),";A = (pT1 - pT2)/(pT1 + pT2);d#sigma/dA", 25,0,1);
                    h_QQ_minv_mHard_ratio[isign][jdphi][kgrp] = new TH1D(Form("h_QQ_minv_mHard_ratio_sign%d%s%s",isign+1, dphis[jdphi].c_str(), ancestor_grps[kgrp].c_str()),";#frac{m_{QQ}}{#hat{s}};d#sigma/d#frac{m_{QQ}}{#hat{s}}", 25,0,1);
                    
                    // h_QQ_ptlead_pair_pt[isign][jdphi] = new TH2D(Form("h_QQ_ptlead_pair_pt_sign%d%s",isign+1, dphis[jdphi].c_str()),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",pms.npairPT_bins,pms.pairPTBins[isign][2],pms.n_hq_pt_bins,pms.hq_pTBins);
                    h_QQ_ptlead_pair_pt[isign][jdphi][kgrp] = new TH2D(Form("h_QQ_ptlead_pair_pt_sign%d%s%s",isign+1, dphis[jdphi].c_str(), ancestor_grps[kgrp].c_str()),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",pms.n_hq_pt_bins,pms.hq_pTBins,pms.n_hq_pt_bins,pms.hq_pTBins);
                    h_QQ_pt1_pt2[isign][jdphi][kgrp] = new TH2D(Form("h_QQ_pt1_pt2_sign%d%s%s",isign+1, dphis[jdphi].c_str(), ancestor_grps[kgrp].c_str()),";p_{T}^{sublead} [GeV];p_{T}^{lead} [GeV]",pms.n_hq_pt_bins,pms.hq_pTBins,pms.n_hq_pt_bins,pms.hq_pTBins);
                    h_QQ_Deta_Dphi[isign][jdphi][kgrp] = new TH2D(Form("h_QQ_Deta_Dphi_sign%d%s%s",isign+1, dphis[jdphi].c_str(), ancestor_grps[kgrp].c_str()),";#Delta#phi;#Delta#eta", 32,-pms.PI,pms.PI,40,-4.8,4.8);
                    h_QQ_eta1_eta2[isign][jdphi][kgrp] = new TH2D(Form("h_QQ_eta1_eta2_sign%d%s%s",isign+1, dphis[jdphi].c_str(), ancestor_grps[kgrp].c_str()),";#eta_{sublead};#eta_{lead}",40,-2.4,2.4, 40,-2.4,2.4);
                    h_QQ_minv_pair_pt[isign][jdphi][kgrp] = new TH2D(Form("h_QQ_minv_pair_pt_sign%d%s%s",isign+1, dphis[jdphi].c_str(), ancestor_grps[kgrp].c_str()),";p_{T}^{pair} [GeV];m_{QQ} [GeV]",pms.n_hq_pt_bins,pms.hq_pTBins,pms.n_hq_minv_bins,pms.hq_minvBins);
                    h_QQ_minv_Dphi[isign][jdphi][kgrp] = new TH2D(Form("h_QQ_minv_Dphi_sign%d%s%s",isign+1, dphis[jdphi].c_str(), ancestor_grps[kgrp].c_str()),";#Delta#phi;m_{QQ} [GeV]",pms.npt_bins,pms.pTBins,pms.n_hq_minv_bins,pms.hq_minvBins);
                }
            }
        }

        for (int isign = 0; isign < ParamsSet::nSigns; isign++){
            h_QQ_both_from_Q_same_ancestors[isign][0] = new TH1D(Form("h_QQ_both_from_Q_same_ancestors_sign%d_near",isign+1),Form("h_QQ_both_from_Q_same_ancestors_sign%d_near",isign+1),2,0,2);
            h_QQ_both_from_Q_same_ancestors[isign][1] = new TH1D(Form("h_QQ_both_from_Q_same_ancestors_sign%d_away",isign+1),Form("h_QQ_both_from_Q_same_ancestors_sign%d_away",isign+1),2,0,2);
            h_QQ_both_from_Q_ancestor_sp[isign][0] = new TH1D(Form("h_QQ_both_from_Q_ancestor_sp_sign%d_near",isign+1),Form("h_QQ_both_from_Q_ancestor_sp_sign%d_near",isign+1),nAncestorGroups,0,nAncestorGroups);
            h_QQ_both_from_Q_ancestor_sp[isign][1] = new TH1D(Form("h_QQ_both_from_Q_ancestor_sp_sign%d_away",isign+1),Form("h_QQ_both_from_Q_ancestor_sp_sign%d_away",isign+1),nAncestorGroups,0,nAncestorGroups);
            h_QQ_both_from_Q_ancestor_dp[isign][0] = new TH2D(Form("h_QQ_both_from_Q_ancestor_dp_sign%d_near",isign+1),Form("h_QQ_both_from_Q_ancestor_dp_sign%d_near",isign+1),nAncestorGroups,0,nAncestorGroups,nAncestorGroups,0,nAncestorGroups);
            h_QQ_both_from_Q_ancestor_dp[isign][1] = new TH2D(Form("h_QQ_both_from_Q_ancestor_dp_sign%d_away",isign+1),Form("h_QQ_both_from_Q_ancestor_dp_sign%d_away",isign+1),nAncestorGroups,0,nAncestorGroups,nAncestorGroups,0,nAncestorGroups);
        }
    }

    // crossx histograms
    static const int ncrossx_bins = 40;
    float crossx_logpow_bb = 0.0711;
    float crossx_logpow_cc = 0.1288;
    float crossx_max_bb = 1.52 *pow(10,9);
    float crossx_max_cc = 2.9 * pow(10,10);
    double crossx_bins_bb[ncrossx_bins+1];
    double crossx_bins_cc[ncrossx_bins+1];

    for(int icrossx = 0; icrossx <= ncrossx_bins; icrossx++){
        crossx_bins_bb[icrossx] = crossx_max_bb * pow(10.0, ((float)(icrossx - ncrossx_bins))*crossx_logpow_bb);
        crossx_bins_cc[icrossx] = crossx_max_cc * pow(10.0, ((float)(icrossx - ncrossx_bins))*crossx_logpow_cc);
    }

    if (mc_mode == "bb"){
        h_crossx = new TH1D("h_crossx","h_crossx",ncrossx_bins,crossx_bins_bb);
        for (int ipt = 0; ipt < npTbins; ipt++){
            // h_crossx_pt_binned[ipt] = new TH1D(Form("h_crossx_pt%d",ipt+1),Form("h_crossx_pt%d",ipt+1),ncrossx_bins,crossx_bins_bb);
            h_crossx_pt_binned[ipt] = new TH1D(Form("h_crossx_pt%d",ipt+1),Form("h_crossx_pt%d",ipt+1),ncrossx_bins,crossx_bins_cc);
        }
    }else{
        h_crossx = new TH1D("h_crossx","h_crossx",ncrossx_bins,crossx_bins_cc);
        for (int ipt = 0; ipt < npTbins; ipt++){
            h_crossx_pt_binned[ipt] = new TH1D(Form("h_crossx_pt%d",ipt+1),Form("h_crossx_pt%d",ipt+1),ncrossx_bins,crossx_bins_cc);
        }
    }
}

void PowhegTruthNTupleFirstPass::PerformTruthPairAnalysis(){
    if (truth_children->at(2).size() == 2){
      crossx_2_to_2 += mpair->weight;
    }else if (truth_children->at(2).size() == 3){
      crossx_2_to_3 += mpair->weight;
    }else{
        std::cout << "The current event (Event# " << mpair->m1.ev_num << ")'s hardest scattering is neither 2-to-2 nor 2-to-3." << std::endl;
        std::cout << "Number of outgoing particles from the hard scattering is:" << truth_children->at(2).size() << std::endl;
        for (auto child : truth_children->at(2)){
            cout << child << " ";
        }
        cout << endl;
    }

    // fill in muon parents (now that we have finished applying all the cuts)
    MuonPairAncestorTracing();
    if (m1_ancestor_is_incoming && m2_ancestor_is_incoming){
        std::cout << "Both muons are from quarks that are incoming partons." << std::endl;
        PrintHistory(&std::cout, false, mpair->from_same_ancestors);
    }

    //------------------------------------------------------------

    // cout << "mQQ " << mpair->mQQ << ", mHard_relevant " << mpair->mHard_relevant << endl;
    FillMuonPairTree(); // mode = 2 for MC truth
    
    h_crossx->Fill(EventWeights->at(0));
    if(mpair->pt_lead < 8) // 4-8
        h_crossx_pt_binned[0]->Fill(EventWeights->at(0));
    else if(mpair->pt_lead < 15) // 8-15
        h_crossx_pt_binned[1]->Fill(EventWeights->at(0));
    else // > 15
        h_crossx_pt_binned[2]->Fill(EventWeights->at(0));

    int isign = !(mpair->same_sign);
    int jdphi = (abs(mpair->dphi) >= pms.PI / 2.);
}



int PowhegTruthNTupleFirstPass::ParentGrouping(std::vector<int>& parent_ids, bool c_tag, bool prev_is_lepton){

  if (parent_ids.size() == 2 && abs(parent_ids[0]) <= 5 && parent_ids[1] == (-1) * parent_ids[0] && prev_is_lepton){
    // std::cout << "For study purpose - Drell-Yan:" << std::endl;
    // PrintHistory(&std::cout, false, true);
    return prt_drell_yan;
  }

  int parent_id = abs(parent_ids[0]);

  if ((parent_id >= 500 && parent_id < 600) || (parent_id >= 5000 && parent_id < 6000)){
    if (! c_tag) return direct_b; // direct b
    else         return b_to_c; // b -> c
  }
  if (c_tag) return direct_c; // c not from b
  if ((parent_id >= 100 && parent_id < 400) || (parent_id >= 1000 && parent_id < 4000)) // light and s-hadrons
    return s_light; // light and s-flavored hadrons

  if (parent_id == 22) // photons
    return single_photon;

  std::cout << "Not in the given set of parent groups. Parent ids: " << std::endl;
  for (auto id : parent_ids) std::cout << id << " ";
  std::cout << std::endl;
  return -1; // others
}

void PowhegTruthNTupleFirstPass::GetPtEtaPhiMFromBarcode(int barcode, std::vector<float>* pt_eta_phi_m){
  std::vector<int>::iterator itbar = std::find(truth_barcode->begin(), truth_barcode->end(), barcode);

  if (itbar == truth_barcode->end()){ // found the barcode of muon1 among all truth particles
    std::cout << "Error:: Barcode of the current parent not found among all truth particles, quitting" << std::endl;
    throw std::exception();
  }
  
  pt_eta_phi_m->clear();
  pt_eta_phi_m->push_back(abs(truth_pt->at(itbar - truth_barcode->begin()))/1000.);
  pt_eta_phi_m->push_back(truth_eta->at(itbar - truth_barcode->begin()));
  pt_eta_phi_m->push_back(truth_phi->at(itbar - truth_barcode->begin()));
}

int PowhegTruthNTupleFirstPass::UpdateCurParents(bool isMuon1, std::vector<int>& cur_prt_bars, std::vector<int>& cur_prt_ids, int hf_quark_index = -1){
  // if (mpair->m1.ev_num == 1 && !isMuon1){
  //   PrintHistory(&std::cout, true, isMuon1);
  // }

  // // print out if there are more than one hadron-level parents
  // check before updating so that the earliest hadron-level parents don't get printed out
  bool cur_is_quark_gluon = (abs(cur_prt_ids[0]) <= 5 || cur_prt_ids[0] == 21);
  // if (cur_prt_bars.size() > 1 && !(cur_is_quark_gluon)){
  if (cur_prt_bars.size() > 1 && hf_quark_index < 0){
    // LATER REMEMBER TO ADD [EXCLUDE 21 21 / 21 2]
    *m_unspecified_parent_file << "Event#: " << mpair->m1.ev_num << std::endl;
    *m_unspecified_parent_file << "More than one parent at hadronic level:" << std::endl;
    for (int parent_id : cur_prt_ids) *m_unspecified_parent_file << parent_id << " ";
    *m_unspecified_parent_file << std::endl << std::endl;
  }

  int prev_first_prt_id = cur_prt_ids[0];
  int prev_first_prt_bar = cur_prt_bars[0];
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

  // these two are difference: e.g, the barcode can be larger than 10000
  // but the index can only go to #particles - 1
  cur_prt_bars = truth_parents->at(itbar - truth_barcode->begin());
  std::vector<int> cur_prt_indices;

  for (int parent_bar : cur_prt_bars){
    itbar = std::find(truth_barcode->begin(),truth_barcode->end(),parent_bar);
    if (itbar == truth_barcode->end()){
      std::cout << "Error:: Barcode of the current parent not found among all truth particles, quitting" << std::endl;
      throw std::exception();
    }
    // cur_prt_ids.push_back(abs(truth_id->at(itbar - truth_barcode->begin())) % 10000);
    int index = itbar - truth_barcode->begin();
    cur_prt_indices.push_back(index);
    cur_prt_ids.push_back(truth_id->at(index) % 10000);
  }

  // std::cout << "Update cur parents: find before saving the Particles" << std::endl;

  std::vector<Particle> cur_prt_profiles {};
  for (int index : cur_prt_indices){
    float pt = abs(truth_pt->at(index)) / 1000.;
    float eta = truth_eta->at(index);
    float phi = truth_phi->at(index);
    int bar = truth_barcode->at(index);
    int id = truth_id->at(index);
    int status = -10;
    Particle p {pt, eta, phi, bar, id, status};
    cur_prt_profiles.push_back(p);
  }

  // std::cout << "Update cur parents: find before saving the history" << std::endl;

  // save the vector of the updated parent ids into the muon history chain
  if (isMuon1){
    m1_history->push_back(cur_prt_ids);
    m1_history_particle->push_back(cur_prt_profiles);
  }
  else{
    m2_history->push_back(cur_prt_ids);
    m2_history_particle->push_back(cur_prt_profiles);
  }

  // record if there is oscillation
  if (cur_prt_ids[0] == (-1) * prev_first_prt_id && abs(prev_first_prt_id) != 4 && abs(prev_first_prt_id) != 5){
    if (isMuon1) m1_osc = true;
    else         m2_osc = true;
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
    if (isMuon1) m1_c_tag = true;
    else         m2_c_tag = true;
  }

  // std::cout << "Update cur parents: find beofre getting kinematics" << std::endl;

  // from muon to hadron stage
  if (abs(prev_first_prt_id) == 13 || abs(prev_first_prt_id) == 15){
    if ((abs(cur_prt_ids[0]) >= 400 && abs(cur_prt_ids[0]) < 600) || (abs(cur_prt_ids[0]) >= 4000 && abs(cur_prt_ids[0]) < 6000)){ // c-hadrons
      if (isMuon1)  GetPtEtaPhiMFromBarcode(cur_prt_bars[0], &mpair->m1_last_hf_hadron_prt_pt_eta_phi_m);
      else          GetPtEtaPhiMFromBarcode(cur_prt_bars[0], &mpair->m2_last_hf_hadron_prt_pt_eta_phi_m);
    }
  }

  // last b-flavored hadron
  bool prev_is_b_hadron = ((abs(prev_first_prt_id) >= 500 && abs(prev_first_prt_id) < 600) || (abs(prev_first_prt_id) >= 5000 && abs(prev_first_prt_id) < 6000));
  bool current_is_b_hadron = ((abs(cur_prt_ids[0]) >= 500 && abs(cur_prt_ids[0]) < 600) || (abs(cur_prt_ids[0]) >= 5000 && abs(cur_prt_ids[0]) < 6000));
  if ((!prev_is_b_hadron) && current_is_b_hadron){ // the last b hadron
    if (isMuon1)  GetPtEtaPhiMFromBarcode(cur_prt_bars[0], &mpair->m1_last_b_hadron_prt_pt_eta_phi_m);
    else          GetPtEtaPhiMFromBarcode(cur_prt_bars[0], &mpair->m2_last_b_hadron_prt_pt_eta_phi_m);
  }

  // from hadron to quark/gluon stage
  bool prev_is_hf_hadron = ((abs(prev_first_prt_id) >= 400 && abs(prev_first_prt_id) < 600) || (abs(prev_first_prt_id) >= 4000 && abs(prev_first_prt_id) < 6000));
  if (prev_is_hf_hadron && ((abs(cur_prt_ids[0]) < 4000 && abs(cur_prt_ids[0]) % 100 <= 5) || cur_prt_ids[0] == 21)){ // including light diquarks
    if (isMuon1)  GetPtEtaPhiMFromBarcode(prev_first_prt_bar, m1_first_hf_hadron_prt_pt_eta_phi_m);
    else          GetPtEtaPhiMFromBarcode(prev_first_prt_bar, m2_first_hf_hadron_prt_pt_eta_phi_m);
    return prev_first_prt_id;
  }
  return 0;
}

int PowhegTruthNTupleFirstPass::FindHeavyQuarks(std::vector<int>& cur_prt_ids, int quark_type, bool isMuon1, int hadron_child_id = 0){
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
      int sign_correction = (abs(hadron_child_id) >= 500 && abs(hadron_child_id) < 600)? -1 : +1; // b quark to B meson changes sign; c to D meson or b/c(bar) to hadron does not change sign
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

void PowhegTruthNTupleFirstPass::SingleMuonAncestorTracing(bool isMuon1){
  // cout << "starting on m1/m2 ancestor tracing" << endl;

  std::vector<int> parent_bars;
  std::vector<int> parent_ids;
  int first_hadron_id = 0;
  int prev_first_prt_id = -1;

  float pt = (isMuon1)? mpair->m1.pt : mpair->m2.pt;
  float eta = (isMuon1)? mpair->m1.eta : mpair->m2.eta;
  float phi = (isMuon1)? mpair->m1.phi : mpair->m2.phi;
  int ind = (isMuon1)? mpair->m1.ind : mpair->m2.ind;
  int charge = (isMuon1)? mpair->m1.charge : mpair->m2.charge;
  int id = -13 * charge;
  int status = -10;
  Particle p {pt, eta, phi, ind, id, status};

  parent_bars.push_back(ind);
  parent_ids.push_back(id);
  std::vector<Particle> cur_muon_profile {};
  cur_muon_profile.push_back(p);

  if (isMuon1){
    m1_history->push_back(parent_ids);
    m1_history_particle->push_back(cur_muon_profile);
  }
  else{
    m2_history->push_back(parent_ids);
    m2_history_particle->push_back(cur_muon_profile);
  }

  // cout << "before lepton/c-tracing step" << endl;

  while(abs(parent_ids[0]) == 13 || abs(parent_ids[0]) == 15 || (abs(parent_ids[0]) >= 400 && abs(parent_ids[0]) < 500) || (abs(parent_ids[0]) >= 4000 && abs(parent_ids[0]) < 5000)){
    prev_first_prt_id = parent_ids[0];
    first_hadron_id = UpdateCurParents(isMuon1,parent_bars,parent_ids);
  }

  // cout << "finished lepton/c-tracing step" << endl;

  bool prev_is_lepton = (abs(prev_first_prt_id) == 13 || abs(prev_first_prt_id) == 15);

  int cur_parent_group; // necessary to determine at this stage for the heavy-quark finding step
  if (isMuon1){
    m1_earliest_parent_id = parent_ids[0];
    cur_m1_earliest_parent_barcode = parent_bars[0];
    cur_parent_group = ParentGrouping(parent_ids, m1_c_tag, prev_is_lepton);
    mpair->m1_parent_group = cur_parent_group;
  } 
  else{
    m2_earliest_parent_id = parent_ids[0];
    cur_m2_earliest_parent_barcode = parent_bars[0];
    cur_parent_group = ParentGrouping(parent_ids, m2_c_tag, prev_is_lepton);
    mpair->m2_parent_group = cur_parent_group;
  }

  // cout << "finished hadronic parent tagging" << endl;

  if (mc_mode == "bb"){
    if (cur_parent_group != direct_b && cur_parent_group != b_to_c){
      mpair->both_from_b = false;
      return;
    }
    while((abs(parent_ids[0]) >= 500 && abs(parent_ids[0]) < 600) || (abs(parent_ids[0]) >= 5000 && abs(parent_ids[0]) < 6000)){
      first_hadron_id = UpdateCurParents(isMuon1,parent_bars,parent_ids);
    }
  }else{ // cc
    // bool first_hadron_is_c_hadron = ((abs(first_hadron_id) >= 400 && abs(first_hadron_id) < 500) || (abs(first_hadron_id) >= 4000 && abs(first_hadron_id) < 5000));
    if (cur_parent_group != direct_c){
      mpair->both_from_c = false;
      return;
    }
  }

  // std::cout << mpair->m1.ev_num << "\t" << first_hadron_id << ",\t";

  int quark = (mc_mode == "bb")? 5:4;
  // if (mpair->m1.ev_num < 80) std::cout << first_hadron_id << std::endl;
  int quark_index = FindHeavyQuarks(parent_ids, quark, isMuon1, first_hadron_id);

  // cout << "finished HQ index initializaiton" << endl;

  // for (int pid : parent_ids) cout << pid << " ";
  // std::cout << "(" << quark_index << "),\t";

  if ((isMuon1 && m1_from_J_psi) || (!isMuon1 && m2_from_J_psi)){
    while(parent_ids.size() == 1 && (abs(parent_ids[0]) == 3 || abs(parent_ids[0]) == 1003 || abs(parent_ids[0]) == 2003)){ // 9940003
      UpdateCurParents(isMuon1,parent_bars,parent_ids);
    }
  }

  if ((isMuon1 && m1_from_Upsilon) || (!isMuon1 && m2_from_Upsilon)){
    while(parent_ids.size() == 1 && (abs(parent_ids[0]) == 1103 || abs(parent_ids[0]) == 203)){ // 9940003
      UpdateCurParents(isMuon1,parent_bars,parent_ids);
    }
  }

  // cout << "finished resonance additional tracing" << endl;

  int prev_hq_bar = -1;
  while(quark_index >= 0){
    prev_hq_bar = parent_bars[quark_index];
    UpdateCurParents(isMuon1,parent_bars,parent_ids,quark_index);
    quark_index = FindHeavyQuarks(parent_ids, quark, isMuon1);
  }

  // cout << "finished all HQ-level tracing" << endl;

  if (prev_hq_bar < 0){
    if (!((isMuon1 && m1_from_J_psi) || (!isMuon1 && m2_from_J_psi) || (isMuon1 && m1_from_Upsilon) || (!isMuon1 && m2_from_Upsilon))){
      std::cout << "Previous HQ barcode is -1: heavy quark never found." << std::endl;
      std::cout << "Event #: " << mpair->m1.ev_num << std::endl;
      PrintHistory(&std::cout, true, isMuon1);
      skip_event = true;

      // if (isMuon1){
      //   for (int vv : (*m1_history)[0]) std::cout << vv << " ";
      //   for (std::vector<std::vector<int>>::iterator it = m1_history->begin() + 1; it < m1_history->end(); it++){
      //     std::cout << "<--- ";
      //     for (int vv : *it){
      //       std::cout << vv << " ";
      //     }
      //   }
      //   std::cout << std::endl;
      // }else{
      //   for (int vv : (*m2_history)[0]) std::cout << vv << " ";
      //   for (std::vector<std::vector<int>>::iterator it = m2_history->begin() + 1; it < m2_history->end(); it++){
      //     std::cout << "<--- ";
      //     for (int vv : *it){
      //       std::cout << vv << " ";
      //     }
      //   }
      //   std::cout << std::endl;
      // }
    }

    return;
  }

  if (isMuon1){
    m1_ancestor_is_incoming = (quark_index == -2); // if Case II: then is incoming
    cur_m1_ancestor_ids = parent_ids;
    cur_m1_ancestor_bars = parent_bars;
    GetPtEtaPhiMFromBarcode(prev_hq_bar, &mpair->m1_first_hq_ancestor_pt_eta_phi_m);
  }
  else{
    m2_ancestor_is_incoming = (quark_index == -2);
    cur_m2_ancestor_ids = parent_ids;
    cur_m2_ancestor_bars = parent_bars;
    GetPtEtaPhiMFromBarcode(prev_hq_bar, &mpair->m2_first_hq_ancestor_pt_eta_phi_m);
  }

  // cout << "finished pt_eta_phi_m_recording" << endl;

  // if (quark_index == -2){ // the current muon is incoming
  //   *m_unspecified_parent_file << "Event#: " << mpair->m1.ev_num << std::endl;
  //   std::string which_muon = (isMuon1)? "1" : "2";
  //   *m_unspecified_parent_file << "Muon " << which_muon << " is incoming:" << std::endl;
  //   *m_unspecified_parent_file << "Sea HQ barcode: " << prev_hq_bar;
  //   *m_unspecified_parent_file << std::endl << "History: ";
  //   std::vector<std::vector<int>>* m_history = (isMuon1)? m1_history : m2_history;
  //   for (int vv : (*m_history)[0]) *m_unspecified_parent_file << vv << " ";
  //   for (std::vector<std::vector<int>>::iterator it = m_history->begin() + 1; it < m_history->end(); it++){
  //     *m_unspecified_parent_file << "<--- ";
  //     for (int vv : *it){
  //       *m_unspecified_parent_file << vv << " ";
  //     }
  //   }
  //   *m_unspecified_parent_file << std::endl << std::endl;
  // }
}

int PowhegTruthNTupleFirstPass::AncestorGrouping(std::vector<int>& ancestor_ids){
  // incoming
  if (ancestor_ids.size() == 1 && ancestor_ids[0] == 2212)
    return incoming;

  // gg
  if (ancestor_ids.size() == 2 && ancestor_ids[0] == 21 && ancestor_ids[1] == 21)
    return gg;
  
  // gq
  if(ancestor_ids.size() == 2 && 
    ((ancestor_ids[0] == 21 && (abs(ancestor_ids[1]) == 1 || abs(ancestor_ids[1]) == 2 || abs(ancestor_ids[1]) == 3)) || 
     (ancestor_ids[1] == 21 && (abs(ancestor_ids[0]) == 1 || abs(ancestor_ids[0]) == 2 || abs(ancestor_ids[0]) == 3))))
    return gq;
  if(mc_mode == "bb" && ancestor_ids.size() == 2 && 
    ((ancestor_ids[0] == 21 && abs(ancestor_ids[1]) == 4) ||
     (ancestor_ids[1] == 21 && abs(ancestor_ids[0]) == 4)))
    return gq;

  // single gluon
  if (ancestor_ids.size() == 1 && (ancestor_ids[0] == 21 || ancestor_ids[0] == 22))
    return single_gluon;
  
  // q qbar
  if (ancestor_ids.size() == 2 && (ancestor_ids[0] == (-1) * ancestor_ids[1])){
    if (abs(ancestor_ids[0]) == 1 || abs(ancestor_ids[0]) == 2 || abs(ancestor_ids[0]) == 3)
      return qqbar;
    if (mc_mode == "bb" && abs(ancestor_ids[0]) == 4)
      return qqbar;
  }

  // others
  std::cout << "Event#: " << mpair->m1.ev_num << std::endl;
  // if (sameprts)
  std::cout << "Ancestor not in any group. Printing out the history of both for better understanding:" << std::endl;
  PrintHistory(&std::cout, false, mpair->from_same_ancestors);
  // for (auto v : ancestor_ids) std::cout << v << " ";
  // std::cout << std::endl << std::endl;
  return -1;
}

void PowhegTruthNTupleFirstPass::HardScatteringAnalysis(std::vector<int>& ancestor_bars, std::vector<int>& ancestor_ids, int sign_dphi_mode, int ancestor_grp){
  // find the outgoing particles from the hard scattering
  // record the number, and kinematics, of the hard-scattering outproducts
  // the latter in the different muon-pair-sign/dphi regions

  int quark = (mc_mode == "bb")? 5:4;

  std::vector<int>::iterator itbar = std::find(truth_barcode->begin(),truth_barcode->end(),ancestor_bars[0]);
  if (itbar == truth_barcode->end()){
    std::cout << "Error:: Barcode of an ancestor not found among all truth particles, quitting" << std::endl;
    throw std::exception();
  }
  std::vector<int> hard_out_bars = truth_children->at(itbar - truth_barcode->begin());
  std::vector<int> hard_out_ids = {}; // this is to store the ids of all outgoing particles, not just the heavy quarks
  std::vector<float> hard_out_qq_pts = {};
  std::vector<float> hard_out_qq_etas = {};
  std::vector<float> hard_out_qq_phis = {};
  std::vector<TLorentzVector> hard_out_lorentz_vecs = {};

  for (int child_bar : hard_out_bars){
    itbar = std::find(truth_barcode->begin(),truth_barcode->end(),child_bar);
    if (itbar == truth_barcode->end()){
      std::cout << "Error:: Barcode of an outgoing particle of the hard scattering not found among all truth particles, quitting" << std::endl;
      throw std::exception();
    }

    int cur_ind = itbar - truth_barcode->begin();
    int cur_id = truth_id->at(cur_ind) % 10000;
    float cur_pt = abs(truth_pt->at(cur_ind))/1000.0;
    float cur_m = abs(truth_m->at(cur_ind))/1000.0;
    float cur_eta = truth_eta->at(cur_ind);
    float cur_phi = truth_phi->at(cur_ind);

    TLorentzVector M;
    M.SetPtEtaPhiM(cur_pt,cur_eta,cur_phi,cur_m);
    hard_out_lorentz_vecs.push_back(M);

    hard_out_ids.push_back(cur_id);
    
    if (abs(cur_id) == quark){
      hard_out_qq_pts.push_back(cur_pt);
      hard_out_qq_etas.push_back(cur_eta);
      hard_out_qq_phis.push_back(cur_phi);
    }
  }

  assert(hard_out_lorentz_vecs.size() != 0);
  TLorentzVector Mtotal = hard_out_lorentz_vecs[0];
  for (auto ii = hard_out_lorentz_vecs.begin()+1; ii < hard_out_lorentz_vecs.end(); ii++){
    Mtotal += *ii;
  }

  float s_cm = Mtotal.M();

  int isign = (sign_dphi_mode >= 2);
  int idphi = sign_dphi_mode % 2;
  if (output_truth_hists) h_num_hard_scatt_out[isign][idphi]->Fill(hard_out_ids.size() + 0.5, mpair->weight);

  if (hard_out_qq_pts.size() != 2){
    std::cout << "Error:: Not exactly 2 heavy quarks found from the hard scattering." << std::endl;
    std::cout << "Event number: " << mpair->m1.ev_num << std::endl;
    std::cout << "Ancestor1 barcode: " << ancestor_bars[0] << std::endl;
    // throw std::exception();
    return;
  }

  qqpair->q1.pt = hard_out_qq_pts[0];
  qqpair->q1.eta = hard_out_qq_etas[0];
  qqpair->q1.phi = hard_out_qq_phis[0];
  qqpair->q2.pt = hard_out_qq_pts[1];
  qqpair->q2.eta = hard_out_qq_etas[1];
  qqpair->q2.phi = hard_out_qq_phis[1];
  qqpair->ev_num = mpair->m1.ev_num;
  qqpair->weight = mpair->weight;
  qqpair->weight = mpair->weight;
  qqpair->ancestor_group = ancestor_grp;
  qqpair->Update();

  if (output_QQpair_tree) QQPairOutTree[isign][idphi][ancestor_grp]->Fill();

  mpair->mQQ = qqpair->minv;
  mpair->mHard_relevant = s_cm;
  // cout << "mQQ " << mpair->mQQ << ", mHard_relevant " << mpair->mHard_relevant << endl;
  
  if (output_truth_hists){    
      h_QQ_Dphi[isign][idphi][ancestor_grp]->Fill(qqpair->dphi, qqpair->weight);
      h_QQ_DR[isign][idphi][ancestor_grp]->Fill(qqpair->dr, qqpair->weight);
      h_QQ_minv[isign][idphi][ancestor_grp]->Fill(qqpair->minv, qqpair->weight);
      h_QQ_pair_pt_ptlead_ratio[isign][idphi][ancestor_grp]->Fill(qqpair->pair_pt / qqpair->pt_lead, qqpair->weight);
      h_QQ_pt_avg[isign][idphi][ancestor_grp]->Fill((qqpair->q1.pt + qqpair->q2.pt) / 2, qqpair->weight);
      h_QQ_asym[isign][idphi][ancestor_grp]->Fill(qqpair->asym, qqpair->weight);
      if (ancestor_grp != 2) // not fill for single gluon case --> treat separately later
        h_QQ_minv_mHard_ratio[isign][idphi][ancestor_grp]->Fill(qqpair->minv / s_cm, qqpair->weight);

      h_QQ_ptlead_pair_pt[isign][idphi][ancestor_grp]->Fill(qqpair->pair_pt, qqpair->pt_lead, qqpair->weight);
      h_QQ_pt1_pt2[isign][idphi][ancestor_grp]->Fill(qqpair->q2.pt, qqpair->pt_lead, qqpair->weight);
      h_QQ_Deta_Dphi[isign][idphi][ancestor_grp]->Fill(qqpair->dphi, qqpair->deta, qqpair->weight);
      h_QQ_eta1_eta2[isign][idphi][ancestor_grp]->Fill(qqpair->q2.eta, qqpair->q1.eta, qqpair->weight);
      h_QQ_minv_pair_pt[isign][idphi][ancestor_grp]->Fill(qqpair->pair_pt, qqpair->minv, qqpair->weight);
      h_QQ_minv_Dphi[isign][idphi][ancestor_grp]->Fill(qqpair->dphi, qqpair->minv, qqpair->weight);
  }

  if (ancestor_grp == 2){
    assert (truth_eta->at(2) == 1000. && truth_eta->at(3) == 1000. && abs(truth_eta)->at(4) < 20);
    int num_hard_out = truth_children->at(2).size();

    std::vector<TLorentzVector> hard_out_lorentz_vecs = {};

    for (int cur_ind = 4; cur_ind < 4 + num_hard_out; cur_ind++){

      float cur_id = truth_id->at(cur_ind);
      float cur_pt = abs(truth_pt->at(cur_ind))/1000.0;
      float cur_m = abs(truth_m->at(cur_ind))/1000.0;
      float cur_eta = truth_eta->at(cur_ind);
      float cur_phi = truth_phi->at(cur_ind);

      TLorentzVector M;
      M.SetPtEtaPhiM(cur_pt,cur_eta,cur_phi,cur_m);
      hard_out_lorentz_vecs.push_back(M);      
    }

    assert(hard_out_lorentz_vecs.size() != 0);
    TLorentzVector Mtotal = hard_out_lorentz_vecs[0];
    for (auto ii = hard_out_lorentz_vecs.begin()+1; ii < hard_out_lorentz_vecs.end(); ii++){
      Mtotal += *ii;
    }

    float s_cm = Mtotal.M();

    mpair->mQQ = qqpair->minv;
    mpair->mHard_relevant = s_cm;

    if (output_truth_hists) h_QQ_minv_mHard_ratio[isign][idphi][ancestor_grp]->Fill(qqpair->minv / s_cm, qqpair->weight);
  }
}


void PowhegTruthNTupleFirstPass::PrintHistory(std::ostream* f, bool print_single, bool muon1_sameancestor){
  // print_single: if True then print single muon history; else then print both muons' history
  // muon1_sameancestor: meaning depends on print_single
  // if (print_single): True = muon1, False = muon2
  // if (!print_single): True = same partonic ancestors, False = different partonic ancestors
  // if True then print single muon history; else then print both muons' history
  *f << "Event #: " << mpair->m1.ev_num << std::endl;

  if (print_single){
    std::vector<std::vector<Particle>>* m_history_particle = (muon1_sameancestor)? m1_history_particle : m2_history_particle;
    // *f << "id (barcode, status, pt, eta, phi)" << std::endl;
    *f << "id (barcode, pt)" << std::endl;
    for (size_t i = 0; i < (*m_history_particle).size(); i++){
      for (size_t j = 0; j < (*m_history_particle)[i].size(); j++) {
        *f << (*m_history_particle)[i][j].id << " (" 
        << (*m_history_particle)[i][j].barcode  << ", "
        << (*m_history_particle)[i][j].pt   << ") ";
        // << (*m_history_particle)[i][j].pt       << ", "
        // << (*m_history_particle)[i][j].eta      << ", "
        // << (*m_history_particle)[i][j].phi      << ") ";
      }
      *f << "<--- ";
    }
    *f << std::endl << std::endl;
  }
  else{ // print out both muons' history

    // if (print_if_same_ancestor){
    //   if (muon1_sameancestor) // same ancestor
    //     *f << "Same ancestors." << std::endl;
    //   else
    //     *f << "Different ancestors." << std::endl;
    // }

    // *f << "id (barcode, status, pt, eta, phi)" << std::endl;
    *f << "id (barcode, pt)" << std::endl;
    for (size_t i = 0; i < (*m1_history_particle).size(); i++){
      for (size_t j = 0; j < (*m1_history_particle)[i].size(); j++) {
        *f << (*m1_history_particle)[i][j].id << " (" 
        << (*m1_history_particle)[i][j].barcode  << ", "
        << (*m1_history_particle)[i][j].pt   << ") ";
        // << (*m1_history_particle)[i][j].pt       << ", "
        // << (*m1_history_particle)[i][j].eta      << ", "
        // << (*m1_history_particle)[i][j].phi      << ") ";
      }
      *f << "<--- ";
    }
    *f << std::endl << std::endl;

    // *f << "id (barcode, status, pt, eta, phi)" << std::endl;
    // *f << "id (barcode, pt)" << std::endl;
    for (size_t i = 0; i < (*m2_history_particle).size(); i++){
      for (size_t j = 0; j < (*m2_history_particle)[i].size(); j++) {
        *f << (*m2_history_particle)[i][j].id << " (" 
        << (*m2_history_particle)[i][j].barcode  << ", "
        << (*m2_history_particle)[i][j].pt   << ") ";
        // << (*m2_history_particle)[i][j].pt       << ", "
        // << (*m2_history_particle)[i][j].eta      << ", "
        // << (*m2_history_particle)[i][j].phi      << ") ";
      }
      *f << "<--- ";
    }
    *f << std::endl << std::endl;
  }
}

void PowhegTruthNTupleFirstPass::MuonPairTagsReinit(){
  m1_c_tag = false;
  m2_c_tag = false;
  m1_osc = false;
  m2_osc = false;
  m1_from_J_psi = false;
  m2_from_J_psi = false;
  m1_from_Upsilon = false;
  m2_from_Upsilon = false;

  skip_event = false;

  mpair->from_same_b = false;
  // mpair->from_same_ancestors = false;
  mpair->both_from_b = true;
  // mpair->one_from_b_one_from_c = false;
  mpair->both_from_c = true;
  mpair->m1_ancestor_category = -10;
  mpair->m2_ancestor_category = -10;
  mpair->mQQ                  = -10.;
  mpair->mHard_relevant       = -10.;

  for (std::vector<int> v : *m1_history) v.clear();
  for (std::vector<int> v : *m2_history) v.clear();
  m1_history->clear();
  m2_history->clear();

  for (std::vector<Particle> v : *m1_history_particle) v.clear();
  for (std::vector<Particle> v : *m2_history_particle) v.clear();
  m1_history_particle->clear();
  m2_history_particle->clear();

  m1_multi_hf_quark_ids.clear();
  m2_multi_hf_quark_ids.clear();

  m1_first_hf_hadron_prt_pt_eta_phi_m->clear();
  m2_first_hf_hadron_prt_pt_eta_phi_m->clear();

}


void PowhegTruthNTupleFirstPass::CheckIfFromSameB(){
  // check if from same b
  // not necessarily [op sign + one from direct b, one from b to c]
  // can also be from some 1-to-n hadronic weak decay (so perhaps both are b to c)
  if (cur_m1_earliest_parent_barcode == cur_m2_earliest_parent_barcode){
    std::vector<int>::iterator it = std::find(truth_barcode->begin(),truth_barcode->end(),cur_m1_earliest_parent_barcode);
    if (it == truth_barcode->end()){
      std::cout << "Error:: Barcode of the first non-c hadronic barcode not found among all truth particles, quitting" << std::endl;
      throw std::exception();
    }

    int prt_ind = it - truth_barcode->begin();

    int prt_id = abs(truth_id->at(prt_ind));
    // necessary to make sure the "hadronic parent" is actually a b-flavored hadron
    // since it's possible that the first not-c-hadron parent vectors is quark level
    // and contains both 4 and -4, where the 4's barcode is recorded as cur_m1_earliest_parent_barcode
    if ((prt_id >= 500 && prt_id < 600) || (prt_id >= 5000 && prt_id < 6000)){
      mpair->from_same_b = true;

      // if (mpair->same_sign){
      //   std::cout << "The same b-flavored hadron give a SAME-SIGN muon pair. How does this happen?" << std::endl;
      //   std::cout << "muon 1 & muon2 hardonic parent barcode: " << cur_m1_earliest_parent_barcode << std::endl;
      //   PrintHistory(&std::cout, false, true);
      // }
    }
  }
}


void PowhegTruthNTupleFirstPass::MuonPairAncestorTracing(){

  MuonPairTagsReinit();
  
  SingleMuonAncestorTracing(true);
  SingleMuonAncestorTracing(false);
  
  bool not_near = !(abs(mpair->dphi) < pms.PI / 2.);

  if (output_truth_hists) h_parent_groups[!mpair->same_sign][not_near]->Fill(mpair->m1_parent_group + 0.5, mpair->m2_parent_group + 0.5, mpair->weight);

  if (skip_event){ // if skip event: return without recording any ancestor-categorizing tags, since these would be inaccurate
    skipped_event_crossx += mpair->weight;
    return; // return without filling in the ancestor-category histograms
  }

  // return if for the bb/cc sample not both from b/c
  // for these muon pairs, at least one's ancestor tracing stops at the hadronic level
  // (which has to happend since FindHeavyQuark() is b/c mode-dependent)
  // hence, cur_m1_ancestor_bars and cur_m1_ancestor_ids, etc. are NOT accurate
  if (mc_mode == "bb" && !mpair->both_from_b) return;
  if (mc_mode == "cc" && !mpair->both_from_c) return;

  bool same_ancestors = (cur_m1_ancestor_bars == cur_m2_ancestor_bars);
  mpair->from_same_ancestors = same_ancestors;

  mpair->m1_ancestor_category = AncestorGrouping(cur_m1_ancestor_ids);
  mpair->m2_ancestor_category = AncestorGrouping(cur_m2_ancestor_ids);

  if (same_ancestors && mpair->m1_ancestor_category != single_gluon && mpair->m1_ancestor_category != -1 && cur_m1_ancestor_bars[0] != 3){
    if (mpair->m1.ev_num < 10000){
      std::cout << "The relevant hard scattering is not the hardest scattering. Printing out history to make sure:" << std::endl;
      PrintHistory(&std::cout, false, true);
    }
    crossx_relevant_hard_isnt_hardest += mpair->weight;
  }
  
  CheckIfFromSameB();
  if (mpair->from_same_b) return;

  bool near_side = (abs(mpair->dphi) < pms.PI / 2.);
  bool small_dphi = (abs(mpair->dphi) < 0.4);

      
  // print out unspecified cases where an ancestor vector contains both b and -b
  if (!same_ancestors && (m1_multi_hf_quark_ids.size() != 0 || m2_multi_hf_quark_ids.size() != 0)){
    std::vector<int> multi_hf_quark_ids = (m1_multi_hf_quark_ids.size() != 0)? m1_multi_hf_quark_ids : m2_multi_hf_quark_ids;
    *m_unspecified_parent_file << "Event#: " << mpair->m1.ev_num << std::endl;
    *m_unspecified_parent_file << "Unexpected. Both Q and -Q found." << std::endl;
    for (auto v : multi_hf_quark_ids) *m_unspecified_parent_file << v << " ";
    *m_unspecified_parent_file << std::endl << std::endl;
  }

  // fill in histograms
  if (output_truth_hists){    
      h_QQ_both_from_Q_same_ancestors[!mpair->same_sign][not_near]->Fill(!same_ancestors + 0.5, mpair->weight); // 0.5 for same parents; 1.5 for different parents
      if (same_ancestors){
        h_QQ_both_from_Q_ancestor_sp[!mpair->same_sign][not_near]->Fill(mpair->m1_ancestor_category + 0.5, mpair->weight);
      }else{
        h_QQ_both_from_Q_ancestor_dp[!mpair->same_sign][not_near]->Fill(mpair->m1_ancestor_category + 0.5, mpair->m2_ancestor_category + 0.5, mpair->weight);
      }
  }

  // print history
  if (print_prt_history && mpair->m1.ev_num < 1000){
    std::ofstream* m_prt_file = (mc_mode == "bb")? m_b_parent_file[!mpair->same_sign][not_near] : m_c_parent_file[!mpair->same_sign][not_near];
    if (same_ancestors)  PrintHistory(m_prt_file, false, true);
    else            PrintHistory(m_prt_file, false, false);
  }

  // analysis of the [QQ sample, both both Q, same ancestors] case --> both muons from hard scattering
  if (same_ancestors){
    // // same sign - categorization
    // if (mpair->same_sign) SameSignSameAncestorsAnalysis(near_side, one_b_one_btoc);

    // analysis of hard scattering outproducts
    int isign = !(mpair->same_sign);
    int jdphi = !(near_side);
    HardScatteringAnalysis(cur_m1_ancestor_bars, cur_m1_ancestor_ids, 2 * isign + jdphi, mpair->m1_ancestor_category);
    if (output_truth_hists) KinematicCorrPlots(isign, jdphi);
  }
}


void PowhegTruthNTupleFirstPass::KinematicCorrPlots(int isign, int idphi){
  h_pt_muon_pt_closest_hadr_ratio[isign][idphi]->Fill(mpair->m1.pt / (mpair->m1_last_hf_hadron_prt_pt_eta_phi_m)[0],mpair->weight);
  h_pt_muon_pt_closest_hadr_ratio[isign][idphi]->Fill(mpair->m2.pt / (mpair->m2_last_hf_hadron_prt_pt_eta_phi_m)[0],mpair->weight);
  h_pt_closest_hadr_pt_furthest_hadr_ratio[isign][idphi]->Fill((mpair->m1_last_hf_hadron_prt_pt_eta_phi_m)[0] / (*m1_first_hf_hadron_prt_pt_eta_phi_m)[0], mpair->weight);
  h_pt_closest_hadr_pt_furthest_hadr_ratio[isign][idphi]->Fill((mpair->m2_last_hf_hadron_prt_pt_eta_phi_m)[0] / (*m2_first_hf_hadron_prt_pt_eta_phi_m)[0], mpair->weight);
  h_pt_hadr_hq_ratio[isign][idphi]->Fill((mpair->m1_last_hf_hadron_prt_pt_eta_phi_m)[0] / (mpair->m1_first_hq_ancestor_pt_eta_phi_m)[0], mpair->weight);
  h_pt_hadr_hq_ratio[isign][idphi]->Fill((mpair->m2_last_hf_hadron_prt_pt_eta_phi_m)[0] / (mpair->m2_first_hq_ancestor_pt_eta_phi_m)[0], mpair->weight);
  float muon_hadr_dphi = mpair->m1.phi - (mpair->m1_last_hf_hadron_prt_pt_eta_phi_m)[2];
  muon_hadr_dphi = atan2(sin(muon_hadr_dphi),cos(muon_hadr_dphi));//fold muon_hadr_dphi to [-pi,pi]
  h_dphi_muon_closest_hadr[isign][idphi]->Fill(abs(muon_hadr_dphi), mpair->weight);
  muon_hadr_dphi = mpair->m2.phi - (mpair->m2_last_hf_hadron_prt_pt_eta_phi_m)[2];
  muon_hadr_dphi = atan2(sin(muon_hadr_dphi),cos(muon_hadr_dphi));//fold muon_hadr_dphi to [-pi,pi]
  h_dphi_muon_closest_hadr[isign][idphi]->Fill(abs(muon_hadr_dphi), mpair->weight);
}



void PowhegTruthNTupleFirstPass::HistAdjust(){
  PowhegNTupleFirstPass::HistAdjust();

  for (int ibin = 0; ibin < 3; ibin++){
    for (int isign = 0; isign < ParamsSet::nSigns; isign++){
      for (int jdphi = 0; jdphi < 2; jdphi++){
        h_num_hard_scatt_out[isign][jdphi]->GetXaxis()->SetBinLabel(ibin+1,num_hard_scatt_out_labels[ibin].c_str());
      }
    }
  }

  for (int ibin = 0; ibin < nParentGroups; ibin++){
    for (int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
      for (int lphi = 0; lphi < 2; lphi++){
        h_parent_groups[ksign][lphi]->GetXaxis()->SetBinLabel(ibin+1,parentGroupLabels[ibin].c_str());
        h_parent_groups[ksign][lphi]->GetYaxis()->SetBinLabel(ibin+1,parentGroupLabels[ibin].c_str());
      }
    }
  }


  for (int lphi = 0; lphi < 2; lphi++){
    for (int isign = 0; isign < ParamsSet::nSigns; isign++){
      for (int ibin = 0; ibin < 2; ibin++){
        h_QQ_both_from_Q_same_ancestors[isign][lphi]->GetXaxis()->SetBinLabel(ibin+1,samePrtsLabels[ibin].c_str());
      }
      for (int ibin = 0; ibin < nAncestorGroups; ibin++){
        h_QQ_both_from_Q_ancestor_sp[isign][lphi]->GetXaxis()->SetBinLabel(ibin+1,ancestor_labels[ibin].c_str());
        h_QQ_both_from_Q_ancestor_dp[isign][lphi]->GetXaxis()->SetBinLabel(ibin+1,ancestor_labels[ibin].c_str());
        h_QQ_both_from_Q_ancestor_dp[isign][lphi]->GetYaxis()->SetBinLabel(ibin+1,ancestor_labels[ibin].c_str());
      }
    }
  }
}
void PowhegTruthNTupleFirstPass::Finalize(){

    if (output_truth_hists){
        HistAdjust();
    }
    
    m_outfile->Write();
    m_outHistFile->Write();

    delete m1_history;
    delete m2_history;
    delete m1_history_particle;
    delete m2_history_particle;
    delete m1_first_hf_hadron_prt_pt_eta_phi_m;
    delete m2_first_hf_hadron_prt_pt_eta_phi_m;

    m_unspecified_parent_file->close();
    delete m_unspecified_parent_file;

    if (print_prt_history){
        if (mc_mode == "bb"){
            for (int isign = 0; isign < ParamsSet::nSigns; isign++){
                m_b_parent_file[isign][0]->close();
                m_b_parent_file[isign][1]->close();
                delete m_b_parent_file[isign][0];
                delete m_b_parent_file[isign][1];
            } 
        }else{
            for (int isign = 0; isign < ParamsSet::nSigns; isign++){
                m_c_parent_file[isign][0]->close();
                m_c_parent_file[isign][1]->close();
                delete m_c_parent_file[isign][0];
                delete m_c_parent_file[isign][1];
            }
        }
    }

    if (print_specific_prt_history){
        if (mc_mode == "bb"){
            m_bb_ss_near_file->close();
            m_bb_ss_away_file->close();
            m_bb_op_near_one_b_one_btoc_others_file->close();
            delete m_bb_ss_near_file;
            delete m_bb_ss_away_file;
            delete m_bb_op_near_one_b_one_btoc_others_file;
        }else{
            m_cc_ss_small_dphi_file->close();
            delete m_cc_ss_small_dphi_file;
        }
    }

    for (int jdphi = 0; jdphi < 2; jdphi++){
        for (int isign = 0; isign < ParamsSet::nSigns; isign++){
            delete h_parent_groups[isign][jdphi];
            // delete h_n_to_2_ancestors[isign][jdphi];

            for (int kgrp = 0; kgrp < nAncestorGroups; kgrp++){

                delete h_QQ_DR[isign][jdphi][kgrp];
                delete h_QQ_Dphi[isign][jdphi][kgrp];
                delete h_QQ_minv[isign][jdphi][kgrp];

                delete h_QQ_ptlead_pair_pt[isign][jdphi][kgrp];
                delete h_QQ_pt1_pt2[isign][jdphi][kgrp];
                delete h_QQ_Deta_Dphi[isign][jdphi][kgrp];
                delete h_QQ_eta1_eta2[isign][jdphi][kgrp];
                delete h_QQ_minv_pair_pt[isign][jdphi][kgrp];
                delete h_QQ_minv_Dphi[isign][jdphi][kgrp];
            }
        }
    }

    for (int jdphi = 0; jdphi < 2; jdphi++){
        for (int isign = 0; isign < ParamsSet::nSigns; isign++){
            delete h_QQ_both_from_Q_same_ancestors[isign][jdphi];
            delete h_QQ_both_from_Q_ancestor_sp[isign][jdphi];
            delete h_QQ_both_from_Q_ancestor_dp[isign][jdphi];
        }
    }
}
