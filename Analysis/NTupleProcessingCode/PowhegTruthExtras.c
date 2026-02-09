#include "PowhegTruthExtras.h"
#include "../Utilities/tchain_helpers.h"

template <class PairT, class Derived>
void PowhegTruthExtras<PairT, Derived>::InitializeExtra(){
    int quark = (self().mcModeRef() == "bb")? 5:4;
    if (output_QQpair_tree){
        qqpair = new TruthQQPair(quark);
    }
}

template <class PairT, class Derived>
void PowhegTruthExtras<PairT, Derived>::InitTempVariablesExtra(){
    single_gluon_history = new std::vector<std::vector<int>>();
    m1_history = new std::vector<std::vector<int>>();
    m2_history = new std::vector<std::vector<int>>();
    m1_history_particle = new std::vector<std::vector<Particle>>();
    m2_history_particle = new std::vector<std::vector<Particle>>();
    m1_first_hf_hadron_prt_pt_eta_phi_m = new std::vector<float>();
    m2_first_hf_hadron_prt_pt_eta_phi_m = new std::vector<float>();
}


template <class PairT, class Derived>
void PowhegTruthExtras<PairT, Derived>::InitInputExtra(){

    if (!self().getIsFullsim()){ // not saved in fullsim
        enable_and_bind(self().fChainRef(), "Q"                        , &Q);
        enable_and_bind(self().fChainRef(), "truth_m"                  , &truth_m);
        enable_and_bind(self().fChainRef(), "truth_children"           , &truth_children);        
    }

    enable_and_bind(self().fChainRef(), "truth_id"                   , &truth_id);
    enable_and_bind(self().fChainRef(), "truth_barcode"              , &truth_barcode);
    enable_and_bind(self().fChainRef(), "truth_qual"                 , &truth_qual);
    enable_and_bind(self().fChainRef(), "truth_pt"                 , &truth_pt);
    enable_and_bind(self().fChainRef(), "truth_eta"                 , &truth_eta);
    enable_and_bind(self().fChainRef(), "truth_phi"                 , &truth_phi);
    enable_and_bind(self().fChainRef(), "truth_parents"              , &truth_parents);
}

template <class PairT, class Derived>
void PowhegTruthExtras<PairT, Derived>::CheckBranchPtrsExtra(){
    if (self().debug_mode) std::cout << "Calling PowhegTruthExtras::CheckBranchPtrsExtra" << std::endl;
    
    auto require = [&](auto* p, const char* name){
        if(!p) throw std::runtime_error(std::string("Null branch pointer: ") + name);
    };

    require(truth_id, "truth_id");
    require(truth_barcode, "truth_barcode");
    require(truth_qual, "truth_qual");
    require(truth_pt, "truth_pt");
    require(truth_eta, "truth_eta");
    require(truth_phi, "truth_phi");
    require(truth_parents, "truth_parents");

    if (!self().getIsFullsim()){ // Q not saved in fullsim
        require(truth_m, "truth_m");
        require(truth_children, "truth_children");
    }

    if (self().debug_mode) std::cout << "truth_id size " << truth_id->size() << ", " << truth_id->at(0) << std::endl;
    if (self().debug_mode) std::cout << "truth_pt size " << truth_pt->size() << ", " << truth_pt->at(0) << std::endl;

}


template <class PairT, class Derived>
void PowhegTruthExtras<PairT, Derived>::InitOutputExtra(){
    InitOutStreamFiles();
}

template <class PairT, class Derived>
void PowhegTruthExtras<PairT, Derived>::InitOutStreamFiles(){
    // ---------------------------------------------------------------------------------------------------------------------------

    if (print_specific_prt_history){
        if (self().mcModeRef() == "cc"){

            m_cc_ss_small_dphi_file = new std::ofstream(Form("%scc_ss_small_dphi.txt", self().mcdirRef().c_str()));
            *m_cc_ss_small_dphi_file << "Event#\tm1-grp\tm2-grp" << std::endl;
        }else{
            m_bb_ss_near_file = new std::ofstream(Form("%sbb_ss_near.txt", self().mcdirRef().c_str()));
            m_bb_ss_away_file = new std::ofstream(Form("%sbb_ss_away.txt", self().mcdirRef().c_str()));
            m_bb_op_near_one_b_one_btoc_others_file = new std::ofstream(Form("%sbb_op_near_one_b_one_btoc_others.txt", self().mcdirRef().c_str()));
        }
    }

    // ---------------------------------------------------------------------------------------------------------------------------

    if (self().mcModeRef() == "bb"){
        m_unspecified_parent_file = new std::ofstream(self().mcdirRef() + "unspecified_parents_bb.txt");
        if (print_prt_history){
            for (int isign = 0; isign < ParamsSet::nSigns; isign++){
                for (int jdphi = 0; jdphi < 2; jdphi++){
                    m_b_parent_file[isign][jdphi] = new std::ofstream(Form("%sb_parents_%s%s.txt", self().mcdirRef().c_str(), sign_labels[isign].c_str(), dphis[jdphi].c_str()));
                }
            }
        }
    }else{
        m_unspecified_parent_file = new std::ofstream(self().mcdirRef() + "unspecified_parents_cc.txt");
        if (print_prt_history){
            for (int isign = 0; isign < ParamsSet::nSigns; isign++){
                for (int jdphi = 0; jdphi < 2; jdphi++){
                    m_c_parent_file[isign][jdphi] = new std::ofstream(Form("%sc_parents_%s%s.txt", self().mcdirRef().c_str(), sign_labels[isign].c_str(), dphis[jdphi].c_str()));
                }
            }
        }
    } 
}

template <class PairT, class Derived>
void PowhegTruthExtras<PairT, Derived>::InitOutputTreesExtra(){    
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

template <class PairT, class Derived>
void PowhegTruthExtras<PairT, Derived>::InitOutputHistsExtra(){
    if (output_truth_hists){

        for (int isign = 0; isign < ParamsSet::nSigns; isign++){
            for (int jdphi = 0; jdphi < 2; jdphi++){
                h_parent_groups[isign][jdphi] = new TH2D(Form("h_parent_groups_sign%d%s",isign+1, dphis[jdphi].c_str()),Form("h_parent_groups_sign%d%s",isign+1, dphis[jdphi].c_str()),nParentGroups,0,nParentGroups,nParentGroups,0,nParentGroups);
                h_num_hard_scatt_out[isign][jdphi] = new TH1D(Form("h_num_hard_scatt_out_sign%d%s",isign+1, dphis[jdphi].c_str()),Form("h_num_hard_scatt_out_sign%d%s",isign+1, dphis[jdphi].c_str()),3,2,5);
                h_pt_muon_pt_closest_hadr_ratio[isign][jdphi] = new TH1D(Form("h_pt_muon_pt_closest_hadr_ratio_sign%d%s",isign+1, dphis[jdphi].c_str()),Form("h_pt_muon_pt_closest_hadr_ratio_sign%d%s",isign+1, dphis[jdphi].c_str()),40,0,1.);
                h_pt_closest_hadr_pt_furthest_hadr_ratio[isign][jdphi] = new TH1D(Form("h_pt_closest_hadr_pt_furthest_hadr_ratio_sign%d%s",isign+1, dphis[jdphi].c_str()),Form("h_pt_closest_hadr_pt_furthest_hadr_ratio_sign%d%s",isign+1, dphis[jdphi].c_str()),40,0,1.);
                h_pt_hadr_hq_ratio[isign][jdphi] = new TH1D(Form("h_pt_hadr_hq_ratio_sign%d%s",isign+1, dphis[jdphi].c_str()),Form("h_pt_hadr_hq_ratio_sign%d%s",isign+1, dphis[jdphi].c_str()),40,0,1.);
                h_dphi_muon_closest_hadr[isign][jdphi] = new TH1D(Form("h_dphi_muon_closest_hadr_sign%d%s",isign+1, dphis[jdphi].c_str()),Form("h_dphi_muon_closest_hadr_sign%d%s",isign+1, dphis[jdphi].c_str()),32,-self().pmsRef().PI,self().pmsRef().PI);
                
                for (int kgrp = 0; kgrp < nAncestorGroups; kgrp++){
                    h_QQ_DR[isign][jdphi][kgrp] = new TH1D(Form("h_QQ_DR_sign%d%s%s",isign+1, dphis[jdphi].c_str(), ancestor_grps[kgrp].c_str()),";#Delta R;d#sigma/d#Delta R", 50,0,5.75);
                    h_QQ_Dphi[isign][jdphi][kgrp] = new TH1D(Form("h_QQ_Dphi_sign%d%s%s",isign+1, dphis[jdphi].c_str(), ancestor_grps[kgrp].c_str()),";#Delta#phi;d#sigma/d#Delta#phi", 32,-self().pmsRef().PI,self().pmsRef().PI);
                    h_QQ_minv[isign][jdphi][kgrp] = new TH1D(Form("h_QQ_minv_sign%d%s%s",isign+1, dphis[jdphi].c_str(), ancestor_grps[kgrp].c_str()),";m_{QQ} [GeV]; d#sigma/dm_{QQ}",self().pmsRef().n_hq_minv_bins,self().pmsRef().hq_minvBins);
                    h_QQ_pair_pt_ptlead_ratio[isign][jdphi][kgrp] = new TH1D(Form("h_QQ_pair_pt_ptlead_ratio_sign%d%s%s",isign+1, dphis[jdphi].c_str(), ancestor_grps[kgrp].c_str()),";#frac{p_{T}^{pair}}{p_{T}^{lead}};d#sigma/d#frac{p_{T}^{pair}}{p_{T}^{lead}}", 25,0,2.);
                    h_QQ_pt_avg[isign][jdphi][kgrp] = new TH1D(Form("h_QQ_pt_avg_sign%d%s%s",isign+1, dphis[jdphi].c_str(), ancestor_grps[kgrp].c_str()),";p_{T}^{avg};d#sigma/dp_{T}^{avg}", self().pmsRef().n_hq_pt_bins,self().pmsRef().hq_pTBins);
                    h_QQ_asym[isign][jdphi][kgrp] = new TH1D(Form("h_QQ_asym_sign%d%s%s",isign+1, dphis[jdphi].c_str(), ancestor_grps[kgrp].c_str()),";A = (pT1 - pT2)/(pT1 + pT2);d#sigma/dA", 25,0,1);
                    h_QQ_minv_mHard_ratio[isign][jdphi][kgrp] = new TH1D(Form("h_QQ_minv_mHard_ratio_sign%d%s%s",isign+1, dphis[jdphi].c_str(), ancestor_grps[kgrp].c_str()),";#frac{m_{QQ}}{#hat{s}};d#sigma/d#frac{m_{QQ}}{#hat{s}}", 25,0,1);
                    
                    // h_QQ_ptlead_pair_pt[isign][jdphi] = new TH2D(Form("h_QQ_ptlead_pair_pt_sign%d%s",isign+1, dphis[jdphi].c_str()),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",self().pmsRef().npairPT_bins,self().pmsRef().pairPTBins[isign][2],self().pmsRef().n_hq_pt_bins,self().pmsRef().hq_pTBins);
                    h_QQ_ptlead_pair_pt[isign][jdphi][kgrp] = new TH2D(Form("h_QQ_ptlead_pair_pt_sign%d%s%s",isign+1, dphis[jdphi].c_str(), ancestor_grps[kgrp].c_str()),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",self().pmsRef().n_hq_pt_bins,self().pmsRef().hq_pTBins,self().pmsRef().n_hq_pt_bins,self().pmsRef().hq_pTBins);
                    h_QQ_pt1_pt2[isign][jdphi][kgrp] = new TH2D(Form("h_QQ_pt1_pt2_sign%d%s%s",isign+1, dphis[jdphi].c_str(), ancestor_grps[kgrp].c_str()),";p_{T}^{sublead} [GeV];p_{T}^{lead} [GeV]",self().pmsRef().n_hq_pt_bins,self().pmsRef().hq_pTBins,self().pmsRef().n_hq_pt_bins,self().pmsRef().hq_pTBins);
                    h_QQ_Deta_Dphi[isign][jdphi][kgrp] = new TH2D(Form("h_QQ_Deta_Dphi_sign%d%s%s",isign+1, dphis[jdphi].c_str(), ancestor_grps[kgrp].c_str()),";#Delta#phi;#Delta#eta", 32,-self().pmsRef().PI,self().pmsRef().PI,40,-4.8,4.8);
                    h_QQ_eta1_eta2[isign][jdphi][kgrp] = new TH2D(Form("h_QQ_eta1_eta2_sign%d%s%s",isign+1, dphis[jdphi].c_str(), ancestor_grps[kgrp].c_str()),";#eta_{sublead};#eta_{lead}",40,-2.4,2.4, 40,-2.4,2.4);
                    h_QQ_minv_pair_pt[isign][jdphi][kgrp] = new TH2D(Form("h_QQ_minv_pair_pt_sign%d%s%s",isign+1, dphis[jdphi].c_str(), ancestor_grps[kgrp].c_str()),";p_{T}^{pair} [GeV];m_{QQ} [GeV]",self().pmsRef().n_hq_pt_bins,self().pmsRef().hq_pTBins,self().pmsRef().n_hq_minv_bins,self().pmsRef().hq_minvBins);
                    h_QQ_minv_Dphi[isign][jdphi][kgrp] = new TH2D(Form("h_QQ_minv_Dphi_sign%d%s%s",isign+1, dphis[jdphi].c_str(), ancestor_grps[kgrp].c_str()),";#Delta#phi;m_{QQ} [GeV]",self().pmsRef().npt_bins,self().pmsRef().pTBins,self().pmsRef().n_hq_minv_bins,self().pmsRef().hq_minvBins);
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

    if (self().mcModeRef() == "bb"){
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

template <class PairT, class Derived>
void PowhegTruthExtras<PairT, Derived>::FillMuonPairExtra(int pair_ind){
    self().mpairRef()->Q          = (!self().getIsFullsim())? Q : -1000.;
}

template <class PairT, class Derived>
void PowhegTruthExtras<PairT, Derived>::PerformTruthPairAnalysis(){

    if (truth_children && truth_children->size() > 2){ // if truth_children branch is skimmed        
        if (truth_children->at(2).size() == 2){
          crossx_2_to_2 += self().mpairRef()->weight;
        }else if (truth_children->at(2).size() == 3){
          crossx_2_to_3 += self().mpairRef()->weight;
        }else{
            std::cout << "The current event (Event# " << self().mpairRef()->m1.ev_num << ")'s hardest scattering is neither 2-to-2 nor 2-to-3." << std::endl;
            std::cout << "Number of outgoing particles from the hard scattering is:" << truth_children->at(2).size() << std::endl;
            for (auto child : truth_children->at(2)){
                cout << child << " ";
            }
            cout << endl;
        }
    }

    // fill in muon parents (now that we have finished applying all the cuts)
    MuonPairAncestorTracing();
    if (m1_ancestor_is_incoming && m2_ancestor_is_incoming){
        if (debug_mode_cout_warnings){            
            std::cout << "Both muons are from quarks that are incoming partons." << std::endl;
            PrintHistory(&std::cout, false, self().mpairRef()->from_same_ancestors);
        }
    }

    //------------------------------------------------------------
    
    if (self().EventWeightsRef()->size() > 0){        
        h_crossx->Fill(self().EventWeightsRef()->at(0));
        if(self().mpairRef()->truth_pt_lead < 8) // 4-8
            h_crossx_pt_binned[0]->Fill(self().EventWeightsRef()->at(0));
        else if(self().mpairRef()->truth_pt_lead < 15) // 8-15
            h_crossx_pt_binned[1]->Fill(self().EventWeightsRef()->at(0));
        else // > 15
            h_crossx_pt_binned[2]->Fill(self().EventWeightsRef()->at(0));
    }

    // int isign = !(self().mpairRef()->truth_same_sign);
    // int jdphi = (abs(self().mpairRef()->truth_dphi) >= self().pmsRef().PI / 2.);
}

template <class PairT, class Derived>
void PowhegTruthExtras<PairT, Derived>::MuonPairAncestorTracing(){

  MuonPairTagsReinit();
  
  SingleMuonAncestorTracing(true);
  SingleMuonAncestorTracing(false);
  
  bool not_near = !(abs(self().mpairRef()->truth_dphi) < self().pmsRef().PI / 2.);

  if (output_truth_hists) h_parent_groups[!self().mpairRef()->truth_same_sign][not_near]->Fill(self().mpairRef()->m1_parent_group + 0.5, self().mpairRef()->m2_parent_group + 0.5, self().mpairRef()->weight);

  if (skip_event){ // if skip event: return without recording any ancestor-categorizing tags, since these would be inaccurate
    skipped_event_crossx += self().mpairRef()->weight;
    return; // return without filling in the ancestor-category histograms
  }

  // return if for the bb/cc sample not both from b/c
  // for these muon pairs, at least one's ancestor tracing stops at the hadronic level
  // (which has to happend since FindHeavyQuark() is b/c mode-dependent)
  // hence, cur_m1_ancestor_bars and cur_m1_ancestor_ids, etc. are NOT accurate
  if (self().mcModeRef() == "bb" && !self().mpairRef()->both_from_b) return;
  if (self().mcModeRef() == "cc" && !self().mpairRef()->both_from_c) return;

  bool same_ancestors = (cur_m1_ancestor_bars == cur_m2_ancestor_bars);
  self().mpairRef()->from_same_ancestors = same_ancestors;

  if    (m1_ancestor_is_incoming) self().mpairRef()->m1_ancestor_category = incoming;
  else  self().mpairRef()->m1_ancestor_category = AncestorGrouping(cur_m1_ancestor_ids);
  
  if    (m2_ancestor_is_incoming) self().mpairRef()->m2_ancestor_category = incoming;
  else  self().mpairRef()->m2_ancestor_category = AncestorGrouping(cur_m2_ancestor_ids);

  if (same_ancestors && self().mpairRef()->m1_ancestor_category != single_gluon && self().mpairRef()->m1_ancestor_category != -1 && cur_m1_ancestor_bars.size() > 0 && cur_m1_ancestor_bars.at(0) != 3){
    if (self().mpairRef()->m1.ev_num < 10000 && debug_mode_cout_warnings){
      std::cout << "The relevant hard scattering is not the hardest scattering. Printing out history to make sure:" << std::endl;
      PrintHistory(&std::cout, false, true);
    }
    crossx_relevant_hard_isnt_hardest += self().mpairRef()->weight;
  }
  
  CheckIfFromSameB();
  if (self().mpairRef()->from_same_b) return;

  bool near_side = (abs(self().mpairRef()->truth_dphi) < self().pmsRef().PI / 2.);
  bool small_dphi = (abs(self().mpairRef()->truth_dphi) < 0.4);

      
  // print out unspecified cases where an ancestor vector contains both b and -b
  if (!same_ancestors && (m1_multi_hf_quark_ids.size() != 0 || m2_multi_hf_quark_ids.size() != 0) && m_unspecified_parent_file) {
    std::vector<int> multi_hf_quark_ids = (m1_multi_hf_quark_ids.size() != 0)? m1_multi_hf_quark_ids : m2_multi_hf_quark_ids;
    *m_unspecified_parent_file << "Event#: " << self().mpairRef()->m1.ev_num << std::endl;
    *m_unspecified_parent_file << "Unexpected. Both Q and -Q found." << std::endl;
    for (auto v : multi_hf_quark_ids) *m_unspecified_parent_file << v << " ";
    *m_unspecified_parent_file << std::endl << std::endl;
  }

  // fill in histograms
  if (output_truth_hists){    
      h_QQ_both_from_Q_same_ancestors[!self().mpairRef()->truth_same_sign][not_near]->Fill(!same_ancestors + 0.5, self().mpairRef()->weight); // 0.5 for same parents; 1.5 for different parents
      if (same_ancestors){
        h_QQ_both_from_Q_ancestor_sp[!self().mpairRef()->truth_same_sign][not_near]->Fill(self().mpairRef()->m1_ancestor_category + 0.5, self().mpairRef()->weight);
      }else{
        h_QQ_both_from_Q_ancestor_dp[!self().mpairRef()->truth_same_sign][not_near]->Fill(self().mpairRef()->m1_ancestor_category + 0.5, self().mpairRef()->m2_ancestor_category + 0.5, self().mpairRef()->weight);
      }
  }

  // print history
  if (print_prt_history && self().mpairRef()->m1.ev_num < 1000){
    std::ofstream* m_prt_file = (self().mcModeRef() == "bb")? m_b_parent_file[!self().mpairRef()->truth_same_sign][not_near] : m_c_parent_file[!self().mpairRef()->truth_same_sign][not_near];
    if (same_ancestors)  PrintHistory(m_prt_file, false, true);
    else            PrintHistory(m_prt_file, false, false);
  }

  // analysis of the [QQ sample, both both Q, same ancestors] case --> both muons from hard scattering
  if (same_ancestors){
    // analysis of hard scattering outproducts
    int isign = !(self().mpairRef()->truth_same_sign);
    int jdphi = !(near_side);
    HardScatteringAnalysis(cur_m1_ancestor_bars, cur_m1_ancestor_ids, 2 * isign + jdphi, self().mpairRef()->m1_ancestor_category);
    if (output_truth_hists) KinematicCorrPlots(isign, jdphi);
  }
}

template <class PairT, class Derived>
void PowhegTruthExtras<PairT, Derived>::SingleMuonAncestorTracing(bool isMuon1){
  if(self().debug_mode) std::cout << "starting on m1/m2 ancestor tracing" << std::endl;

  std::vector<int> parent_bars;
  std::vector<int> parent_ids;
  int first_hadron_id = 0;
  int prev_first_prt_id = -1;

  float pt = (isMuon1)? self().mpairRef()->m1.truth_pt : self().mpairRef()->m2.truth_pt;
  float eta = (isMuon1)? self().mpairRef()->m1.truth_eta : self().mpairRef()->m2.truth_eta;
  float phi = (isMuon1)? self().mpairRef()->m1.truth_phi : self().mpairRef()->m2.truth_phi;
  int ind = (isMuon1)? self().mpairRef()->m1.truth_bar : self().mpairRef()->m2.truth_bar;
  int charge = (isMuon1)? self().mpairRef()->m1.truth_charge : self().mpairRef()->m2.truth_charge;
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

  if(self().debug_mode) std::cout << "before lepton/c-tracing step" << std::endl;

  while(!parent_ids.empty() && (abs(parent_ids.at(0)) == 13 || abs(parent_ids.at(0)) == 15 || (abs(parent_ids.at(0)) >= 400 && abs(parent_ids.at(0)) < 500) || (abs(parent_ids.at(0)) >= 4000 && abs(parent_ids.at(0)) < 5000))) {
    prev_first_prt_id = parent_ids.at(0);
    first_hadron_id = UpdateCurParents(isMuon1,parent_bars,parent_ids);
  }

  if(self().debug_mode) std::cout << "finished lepton/c-tracing step" << std::endl;

  bool prev_is_lepton = (abs(prev_first_prt_id) == 13 || abs(prev_first_prt_id) == 15);

  if (parent_ids.empty()){
    if (self().debug_mode){
      std::cout << "Muon history ends at the lepton/c-tracing step. Print out history for sanity check!" << std::endl;
      PrintHistory(&std::cout, true, isMuon1);
    }
    return;
  }

  int cur_parent_group; // necessary to determine at this stage for the heavy-quark finding step
  if (isMuon1){
    m1_earliest_parent_id = parent_ids.at(0);
    cur_m1_earliest_parent_barcode = parent_bars.at(0);
    cur_parent_group = ParentGrouping(parent_ids, m1_c_tag, prev_is_lepton);
    self().mpairRef()->m1_parent_group = cur_parent_group;
  } 
  else{
    m2_earliest_parent_id = parent_ids.at(0);
    cur_m2_earliest_parent_barcode = parent_bars.at(0);
    cur_parent_group = ParentGrouping(parent_ids, m2_c_tag, prev_is_lepton);
    self().mpairRef()->m2_parent_group = cur_parent_group;
  }

  if(self().debug_mode) std::cout << "finished hadronic parent tagging" << std::endl;

  if (self().mcModeRef() == "bb"){
    if (cur_parent_group != direct_b && cur_parent_group != b_to_c){
      self().mpairRef()->both_from_b = false;
      return;
    }
    while(!parent_ids.empty() && ((abs(parent_ids.at(0)) >= 500 && abs(parent_ids.at(0)) < 600) || (abs(parent_ids.at(0)) >= 5000 && abs(parent_ids.at(0)) < 6000))){
      first_hadron_id = UpdateCurParents(isMuon1,parent_bars,parent_ids);
    }
  }else{ // cc
    // bool first_hadron_is_c_hadron = ((abs(first_hadron_id) >= 400 && abs(first_hadron_id) < 500) || (abs(first_hadron_id) >= 4000 && abs(first_hadron_id) < 5000));
    if (cur_parent_group != direct_c){
      self().mpairRef()->both_from_c = false;
      return;
    }
  }

  if (parent_ids.empty()){
    if (self().debug_mode){
      std::cout << "Muon history ends at the b-tracing step. Print out history for sanity check!" << std::endl;
      PrintHistory(&std::cout, true, isMuon1);
    }
    return;
  }

  int quark = (self().mcModeRef() == "bb")? 5:4;
  int quark_index = FindHeavyQuarks(parent_ids, quark, isMuon1, first_hadron_id);

  if(self().debug_mode) std::cout << "finished HQ index initializaiton" << std::endl;


  if ((isMuon1 && m1_from_J_psi) || (!isMuon1 && m2_from_J_psi)){
    while(parent_ids.size() == 1 && (abs(parent_ids.at(0)) == 3 || abs(parent_ids.at(0)) == 1003 || abs(parent_ids.at(0)) == 2003)){ // 9940003
      UpdateCurParents(isMuon1,parent_bars,parent_ids);
    }
  }

  if ((isMuon1 && m1_from_Upsilon) || (!isMuon1 && m2_from_Upsilon)){
    while(parent_ids.size() == 1 && (abs(parent_ids.at(0)) == 1103 || abs(parent_ids.at(0)) == 203)){ // 9940003
      UpdateCurParents(isMuon1,parent_bars,parent_ids);
    }
  }

  if(self().debug_mode) std::cout << "finished resonance additional tracing" << std::endl;

  int prev_hq_bar = -1;
  while(quark_index >= 0){
    if (parent_bars.size() > quark_index){
        prev_hq_bar = parent_bars.at(quark_index);    
    } else{
        std::cout << "WARNING: during heavy quark loop, parent_bars size is NOT large enough for quark_index!" << std::endl;
        std::cout << "parent_bars size: " << parent_bars.size() << ", quark_index: " << quark_index << std::endl;
        break;
    }
    
    UpdateCurParents(isMuon1,parent_bars,parent_ids,quark_index);
    if (parent_bars.empty()){
        quark_index = -2;
        break;
    }
    quark_index = FindHeavyQuarks(parent_ids, quark, isMuon1);
  }

  if(self().debug_mode) std::cout << "finished all HQ-level tracing" << std::endl;

  if (prev_hq_bar < 0){
    if (!((isMuon1 && m1_from_J_psi) || (!isMuon1 && m2_from_J_psi) || (isMuon1 && m1_from_Upsilon) || (!isMuon1 && m2_from_Upsilon))){
      std::cout << "Previous HQ barcode is -1: heavy quark never found." << std::endl;
      std::cout << "Event #: " << self().mpairRef()->m1.ev_num << std::endl;
      PrintHistory(&std::cout, true, isMuon1);
      skip_event = true;
    }

    return;
  }

  if (isMuon1){
    m1_ancestor_is_incoming = (quark_index == -2); // if Case II: then is incoming
    cur_m1_ancestor_ids = parent_ids;
    cur_m1_ancestor_bars = parent_bars;
    GetPtEtaPhiMFromBarcode(prev_hq_bar, &self().mpairRef()->m1_first_hq_ancestor_pt_eta_phi_m);
  }
  else{
    m2_ancestor_is_incoming = (quark_index == -2);
    cur_m2_ancestor_ids = parent_ids;
    cur_m2_ancestor_bars = parent_bars;
    GetPtEtaPhiMFromBarcode(prev_hq_bar, &self().mpairRef()->m2_first_hq_ancestor_pt_eta_phi_m);
  }

  if(self().debug_mode) std::cout << "finished pt_eta_phi_m_recording" << std::endl;
}


template <class PairT, class Derived>
int PowhegTruthExtras<PairT, Derived>::ParentGrouping(std::vector<int>& parent_ids, bool c_tag, bool prev_is_lepton){

  if (parent_ids.size() == 2 && abs(parent_ids.at(0)) <= 5 && parent_ids.at(1) == (-1) * parent_ids.at(0) && prev_is_lepton){
    // std::cout << "For study purpose - Drell-Yan:" << std::endl;
    // PrintHistory(&std::cout, false, true);
    return prt_drell_yan;
  }

  int parent_id = abs(parent_ids.at(0));

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

template <class PairT, class Derived>
void PowhegTruthExtras<PairT, Derived>::GetPtEtaPhiMFromBarcode(int barcode, std::vector<float>* pt_eta_phi_m){
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

template <class PairT, class Derived>
int PowhegTruthExtras<PairT, Derived>::UpdateCurParents(bool isMuon1, std::vector<int>& cur_prt_bars, std::vector<int>& cur_prt_ids, int hf_quark_index){
    
    if (cur_prt_bars.empty() || cur_prt_ids.empty()){
        std::cout << "UpdateCurParents:: cur_prt_bars or cur_prt_ids empty at start of calling UpdateCurParents!" << std::endl;
        PrintHistory(&std::cout, true, isMuon1);        return 0;
    }
    
    // // print out if there are more than one hadron-level parents
    // check before updating so that the earliest hadron-level parents don't get printed out
    bool cur_is_quark_gluon = (!cur_prt_ids.empty() && (abs(cur_prt_ids.at(0)) <= 5 || cur_prt_ids.at(0) == 21));
    // if (cur_prt_bars.size() > 1 && !(cur_is_quark_gluon)){
    if (cur_prt_bars.size() > 1 && hf_quark_index < 0 && m_unspecified_parent_file){
      // LATER REMEMBER TO ADD [EXCLUDE 21 21 / 21 2]
      *m_unspecified_parent_file << "Event#: " << self().mpairRef()->m1.ev_num << std::endl;
      *m_unspecified_parent_file << "More than one parent at hadronic level:" << std::endl;
      for (int parent_id : cur_prt_ids) *m_unspecified_parent_file << parent_id << " ";
      *m_unspecified_parent_file << std::endl << std::endl;
    }
    
    int prev_first_prt_id = cur_prt_ids.at(0);
    int prev_first_prt_bar = cur_prt_bars.at(0);
    
    cur_prt_ids.clear();
    
    std::vector<int>::iterator itbar;
    if (hf_quark_index > 0){
        try{
            itbar = std::find(truth_barcode->begin(),truth_barcode->end(),cur_prt_bars.at(hf_quark_index));
        } catch(const std::out_of_range& e){
            std::cerr << "UpdateCurParents:: ERROR: out_of_range exception caught when accessing cur_prt_bars.at(hf_quark_index): " << e.what() << std::endl;
            std::cerr << "Throwing exception for hf_quark_index: " << hf_quark_index << " and cur_prt_bars size: " << cur_prt_bars.size() << std::endl;
            std::cerr << "History:" << std::endl;
            PrintHistory(&std::cout, true, isMuon1);
            throw std::exception();
        }
  
    }
    else
      itbar = std::find(truth_barcode->begin(),truth_barcode->end(),cur_prt_bars.at(0));
    
    if (itbar == truth_barcode->end()){ // found the barcode of muon1 among all truth particles
      std::cout << "Error:: Barcode of the current parent not found among all truth particles, quitting" << std::endl;
      throw std::exception();
    }

    // these two are difference: e.g, the barcode can be larger than 10000
    // but the index can only go to #particles - 1
    cur_prt_bars = truth_parents->at(itbar - truth_barcode->begin());
    std::vector<int> cur_prt_indices;
    
    if (cur_prt_bars.empty()){      return 0;  
    } 

    for (int parent_bar : cur_prt_bars){
      itbar = std::find(truth_barcode->begin(),truth_barcode->end(),parent_bar);
      if (itbar == truth_barcode->end()){
        std::cout << "Error:: Barcode of the current parent not found among all truth particles, quitting" << std::endl;
        throw std::exception();
      }
      // cur_prt_ids.push_back(abs(truth_id->at(itbar - truth_barcode->begin())) % 10000);
      int index = itbar - truth_barcode->begin();
      cur_prt_indices.push_back(index);
      try{
        cur_prt_ids.push_back(truth_id->at(index) % 10000);
      } catch(const std::out_of_range& e){
        std::cout << "Error:: Truth_id->at(index) not valid in UpdateCurParents(), quitting" << std::endl;
        throw std::exception();
      }
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
    if (cur_prt_ids.at(0) == (-1) * prev_first_prt_id && abs(prev_first_prt_id) != 4 && abs(prev_first_prt_id) != 5){
      if (isMuon1) m1_osc = true;
      else         m2_osc = true;
    }

    if (cur_prt_ids.at(0) == 441 || cur_prt_ids.at(0) == 443 || cur_prt_ids.at(0) == 445){
      if (isMuon1) m1_from_J_psi = true;
      else         m2_from_J_psi = true;
    }

    if (cur_prt_ids.at(0) == 551 || cur_prt_ids.at(0) == 553 || cur_prt_ids.at(0) == 555){
      if (isMuon1) m1_from_Upsilon = true;
      else         m2_from_Upsilon = true;
    }

    //c-tag
    if ((abs(cur_prt_ids.at(0)) >= 400 && abs(cur_prt_ids.at(0)) < 500) || (abs(cur_prt_ids.at(0)) >= 4000 && abs(cur_prt_ids.at(0)) < 5000)){ // c-hadrons
      if (isMuon1) m1_c_tag = true;
      else         m2_c_tag = true;
    }

    // std::cout << "Update cur parents: find beofre getting kinematics" << std::endl;

    // from muon to hadron stage
    if (abs(prev_first_prt_id) == 13 || abs(prev_first_prt_id) == 15){
      if ((abs(cur_prt_ids.at(0)) >= 400 && abs(cur_prt_ids.at(0)) < 600) || (abs(cur_prt_ids.at(0)) >= 4000 && abs(cur_prt_ids.at(0)) < 6000)){ // c-hadrons
        if (isMuon1)  GetPtEtaPhiMFromBarcode(cur_prt_bars.at(0), &self().mpairRef()->m1_last_hf_hadron_prt_pt_eta_phi_m);
        else          GetPtEtaPhiMFromBarcode(cur_prt_bars.at(0), &self().mpairRef()->m2_last_hf_hadron_prt_pt_eta_phi_m);
      }
    }

    // last b-flavored hadron
    bool prev_is_b_hadron = ((abs(prev_first_prt_id) >= 500 && abs(prev_first_prt_id) < 600) || (abs(prev_first_prt_id) >= 5000 && abs(prev_first_prt_id) < 6000));
    bool current_is_b_hadron = ((abs(cur_prt_ids.at(0)) >= 500 && abs(cur_prt_ids.at(0)) < 600) || (abs(cur_prt_ids.at(0)) >= 5000 && abs(cur_prt_ids.at(0)) < 6000));
    if ((!prev_is_b_hadron) && current_is_b_hadron){ // the last b hadron
      if (isMuon1)  GetPtEtaPhiMFromBarcode(cur_prt_bars.at(0), &self().mpairRef()->m1_last_b_hadron_prt_pt_eta_phi_m);
      else          GetPtEtaPhiMFromBarcode(cur_prt_bars.at(0), &self().mpairRef()->m2_last_b_hadron_prt_pt_eta_phi_m);
    }

    // from hadron to quark/gluon stage
    bool prev_is_hf_hadron = ((abs(prev_first_prt_id) >= 400 && abs(prev_first_prt_id) < 600) || (abs(prev_first_prt_id) >= 4000 && abs(prev_first_prt_id) < 6000));
    if (prev_is_hf_hadron && ((abs(cur_prt_ids.at(0)) < 4000 && abs(cur_prt_ids.at(0)) % 100 <= 5) || cur_prt_ids.at(0) == 21)){ // including light diquarks
      if (isMuon1)  GetPtEtaPhiMFromBarcode(prev_first_prt_bar, m1_first_hf_hadron_prt_pt_eta_phi_m);
      else          GetPtEtaPhiMFromBarcode(prev_first_prt_bar, m2_first_hf_hadron_prt_pt_eta_phi_m);      return prev_first_prt_id;
    }

    return 0;
}

template <class PairT, class Derived>
int PowhegTruthExtras<PairT, Derived>::FindHeavyQuarks(std::vector<int>& cur_prt_ids, int quark_type, bool isMuon1, int hadron_child_id){
    
    if (quark_type != 4 && quark_type != 5){
      std::cout << "Error:: the parameter quark_type must take value of 4 (c) or 5 (b), quitting" << std::endl;
      throw std::exception();
    }

    if (cur_prt_ids.size() == 0){
        if (self().getIsFullsim()){
            if (debug_mode_cout_warnings){
                std::cout << "FindHeavyQuarks: WARNING:: parent list is empty, assume incoming" << std::endl;
                PrintHistory(&std::cout, true, isMuon1);
            }
            return -2;
        }else{
            std::cout << "FindHeavyQuarks: Error:: parent list is empty, quitting" << std::endl;
            PrintHistory(&std::cout, true, isMuon1);
            throw std::exception();
        }
    }

    if (cur_prt_ids.size() == 1 && cur_prt_ids.at(0) == 2212)
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

    return -1;
}

template <class PairT, class Derived>
int PowhegTruthExtras<PairT, Derived>::AncestorGrouping(std::vector<int>& ancestor_ids){
    // incoming
    if (ancestor_ids.size() == 1 && ancestor_ids.at(0) == 2212)
        return incoming;

    // gg
    if (ancestor_ids.size() == 2 && ancestor_ids.at(0) == 21 && ancestor_ids.at(1) == 21)
        return gg;
    
    // gq
    if(ancestor_ids.size() == 2 && 
      ((ancestor_ids.at(0) == 21 && (abs(ancestor_ids.at(1)) == 1 || abs(ancestor_ids.at(1)) == 2 || abs(ancestor_ids.at(1)) == 3)) || 
       (ancestor_ids.at(1) == 21 && (abs(ancestor_ids.at(0)) == 1 || abs(ancestor_ids.at(0)) == 2 || abs(ancestor_ids.at(0)) == 3))))
        return gq;
    if(self().mcModeRef() == "bb" && ancestor_ids.size() == 2 && 
      ((ancestor_ids.at(0) == 21 && abs(ancestor_ids.at(1)) == 4) ||
       (ancestor_ids.at(1) == 21 && abs(ancestor_ids.at(0)) == 4)))
        return gq;

    // single gluon
    if (ancestor_ids.size() == 1 && (ancestor_ids.at(0) == 21 || ancestor_ids.at(0) == 22))
        return single_gluon;
    
    // q qbar
    if (ancestor_ids.size() == 2 && (ancestor_ids.at(0) == (-1) * ancestor_ids.at(1))){
        if (abs(ancestor_ids.at(0)) == 1 || abs(ancestor_ids.at(0)) == 2 || abs(ancestor_ids.at(0)) == 3)
            return qqbar;
        if (self().mcModeRef() == "bb" && abs(ancestor_ids.at(0)) == 4)
            return qqbar;
    }

    // others
    // if (sameprts)
    if (debug_mode_cout_warnings){
        std::cout << "Event#: " << self().mpairRef()->m1.ev_num << std::endl;
        std::cout << "Ancestor not in any group. Printing out the history of both for better understanding:" << std::endl;
        PrintHistory(&std::cout, false, self().mpairRef()->from_same_ancestors);    
    }

    return -1;
}

template <class PairT, class Derived>
void PowhegTruthExtras<PairT, Derived>::HardScatteringAnalysis(std::vector<int>& ancestor_bars, std::vector<int>& ancestor_ids, int sign_dphi_mode, int ancestor_grp){
    if (self().debug_mode) std::cout << "Calling HardScatteringAnalysis" << std::endl;

    if (ancestor_grp < 0 || ancestor_grp >= powheg_ancestor_categories::NANCESTOERS){
        std::cerr << "HardScatteringAnalysis:: WARNING: ancestor group passed in (" << ancestor_grp << ") is INVALID! Return without hard scattering analysis!" << std::endl;
        std::cerr << "Print out history for checks:" << std::endl;
        PrintHistory(&std::cerr, false, true);
        return;
    }

    if (ancestor_grp == powheg_ancestor_categories::incoming){
        std::cout << "HardScatteringAnalysis:: INFO: both muons' ancestors are incoming. Return without hard scattering analysis!" << std::endl;
        std::cout << "Print out pair history for checks:" << std::endl;
        PrintHistory(&std::cout, false, true);
        return;
    }


    // find the outgoing particles from the hard scattering
    // record the number, and kinematics, of the hard-scattering outproducts
    // the latter in the different muon-pair-sign/dphi regions

    // if truth_children branch is not available, cannot perform hard scattering analysis --> return
    if (!truth_children || !truth_m) return;

    if (ancestor_bars.empty()){
      std::cout << "HardScatteringAnalysis:: WARNING: ancestor_bars is empty: return without performing analysis!" << std::endl;
      return;
    }

    int quark = (self().mcModeRef() == "bb")? 5:4;

    std::vector<int>::iterator itbar = std::find(truth_barcode->begin(),truth_barcode->end(),ancestor_bars.at(0));
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
    TLorentzVector Mtotal = hard_out_lorentz_vecs.at(0);
    for (auto ii = hard_out_lorentz_vecs.begin()+1; ii < hard_out_lorentz_vecs.end(); ii++){
      Mtotal += *ii;
    }

    float s_cm = Mtotal.M();

    int isign = (sign_dphi_mode >= 2);
    int idphi = sign_dphi_mode % 2;
    if (output_truth_hists) h_num_hard_scatt_out[isign][idphi]->Fill(hard_out_ids.size() + 0.5, self().mpairRef()->weight);

    if (hard_out_qq_pts.size() != 2){
      std::cout << "Error:: Not exactly 2 heavy quarks found from the hard scattering." << std::endl;
      std::cout << "Event number: " << self().mpairRef()->m1.ev_num << std::endl;
      std::cout << "Ancestor1 barcode: " << ancestor_bars.at(0) << std::endl;
      // throw std::exception();
      return;
    }

    if (qqpair){
        qqpair->Clear();
        qqpair->q1.pt = hard_out_qq_pts.at(0);
        qqpair->q1.eta = hard_out_qq_etas.at(0);
        qqpair->q1.phi = hard_out_qq_phis.at(0);
        qqpair->q2.pt = hard_out_qq_pts.at(1);
        qqpair->q2.eta = hard_out_qq_etas.at(1);
        qqpair->q2.phi = hard_out_qq_phis.at(1);
        qqpair->ev_num = self().mpairRef()->m1.ev_num;
        qqpair->weight = self().mpairRef()->weight;
        qqpair->weight = self().mpairRef()->weight;
        qqpair->ancestor_group = ancestor_grp;
        qqpair->Update();

        if (output_QQpair_tree){
            QQPairOutTree[isign][idphi][ancestor_grp]->Fill();
        }

        self().mpairRef()->mQQ = qqpair->minv;
    }

    self().mpairRef()->mHard_relevant = s_cm;
    // cout << "mQQ " << self().mpairRef()->mQQ << ", mHard_relevant " << self().mpairRef()->mHard_relevant << endl;
    
    if (qqpair && output_truth_hists){    
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

      self().mpairRef()->mHard_relevant = s_cm;

      if (qqpair){
          self().mpairRef()->mQQ = qqpair->minv;
          if (output_truth_hists) h_QQ_minv_mHard_ratio[isign][idphi][ancestor_grp]->Fill(qqpair->minv / s_cm, qqpair->weight);        
      }
    }

    if (self().debug_mode) std::cout << "Finished calling HardScatteringAnalysis" << std::endl;
}

template <class PairT, class Derived>
void PowhegTruthExtras<PairT, Derived>::PrintHistory(std::ostream* f, bool print_single, bool muon1_sameancestor){
  // print_single: if True then print single muon history; else then print both muons' history
  // muon1_sameancestor: meaning depends on print_single
  // if (print_single): True = muon1, False = muon2
  // if (!print_single): True = same partonic ancestors, False = different partonic ancestors
  // if True then print single muon history; else then print both muons' history
  *f << "Event #: " << self().mpairRef()->m1.ev_num << std::endl;

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

template <class PairT, class Derived>
void PowhegTruthExtras<PairT, Derived>::MuonPairTagsReinit(){
    m1_c_tag = false;
    m2_c_tag = false;
    m1_osc = false;
    m2_osc = false;
    m1_from_J_psi = false;
    m2_from_J_psi = false;
    m1_from_Upsilon = false;
    m2_from_Upsilon = false;

    skip_event = false;

    self().mpairRef()->from_same_b = false;
    // self().mpairRef()->from_same_ancestors = false;
    self().mpairRef()->both_from_b = true;
    // self().mpairRef()->one_from_b_one_from_c = false;
    self().mpairRef()->both_from_c = true;
    self().mpairRef()->m1_ancestor_category = -10;
    self().mpairRef()->m2_ancestor_category = -10;
    self().mpairRef()->mQQ                  = -10.;
    self().mpairRef()->mHard_relevant       = -10.;

    cur_m1_earliest_parent_barcode = -10;
    cur_m2_earliest_parent_barcode = -10;
    m1_earliest_parent_id = -10;
    m2_earliest_parent_id = -10;
    
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

    cur_m1_ancestor_ids.clear();
    cur_m1_ancestor_bars.clear();
    cur_m2_ancestor_ids.clear();
    cur_m2_ancestor_bars.clear();
}

template <class PairT, class Derived>
void PowhegTruthExtras<PairT, Derived>::CheckIfFromSameB(){
  // check if from same b
  // not necessarily [op sign + one from direct b, one from b to c]
  // can also be from some 1-to-n hadronic weak decay (so perhaps both are b to c)
  if (cur_m1_earliest_parent_barcode == cur_m2_earliest_parent_barcode && cur_m1_earliest_parent_barcode != -10){
    std::vector<int>::iterator it = std::find(truth_barcode->begin(),truth_barcode->end(),cur_m1_earliest_parent_barcode);
    if (it == truth_barcode->end()){
      std::cerr << "Error:: Barcode of the first non-c hadronic barcode not found among all truth particles, quitting" << std::endl;
      std::cerr << "cur_m1_earliest_parent_barcode: " << cur_m1_earliest_parent_barcode << std::endl;
      std::cerr << "cur_m2_earliest_parent_barcode: " << cur_m2_earliest_parent_barcode << std::endl;
      std::cerr << "m1_earliest_parent_id: " << m1_earliest_parent_id << std::endl;
      std::cerr << "m2_earliest_parent_id: " << m2_earliest_parent_id << std::endl;
      PrintHistory(&std::cerr, true, true);
      throw std::exception();
    }

    int prt_ind = it - truth_barcode->begin();

    int prt_id = abs(truth_id->at(prt_ind));
    // necessary to make sure the "hadronic parent" is actually a b-flavored hadron
    // since it's possible that the first not-c-hadron parent vectors is quark level
    // and contains both 4 and -4, where the 4's barcode is recorded as cur_m1_earliest_parent_barcode
    if ((prt_id >= 500 && prt_id < 600) || (prt_id >= 5000 && prt_id < 6000)){
      self().mpairRef()->from_same_b = true;

      // if (self().mpairRef()->truth_same_sign){
      //   std::cout << "The same b-flavored hadron give a SAME-SIGN muon pair. How does this happen?" << std::endl;
      //   std::cout << "muon 1 & muon2 hardonic parent barcode: " << cur_m1_earliest_parent_barcode << std::endl;
      //   PrintHistory(&std::cout, false, true);
      // }
    }
  }
}

template <class PairT, class Derived>
void PowhegTruthExtras<PairT, Derived>::KinematicCorrPlots(int isign, int idphi){
    try{
        h_pt_muon_pt_closest_hadr_ratio[isign][idphi]->Fill(self().mpairRef()->m1.truth_pt / (self().mpairRef()->m1_last_hf_hadron_prt_pt_eta_phi_m).at(0),self().mpairRef()->weight);
        h_pt_muon_pt_closest_hadr_ratio[isign][idphi]->Fill(self().mpairRef()->m2.truth_pt / (self().mpairRef()->m2_last_hf_hadron_prt_pt_eta_phi_m).at(0),self().mpairRef()->weight);
        h_pt_closest_hadr_pt_furthest_hadr_ratio[isign][idphi]->Fill((self().mpairRef()->m1_last_hf_hadron_prt_pt_eta_phi_m).at(0) / (*m1_first_hf_hadron_prt_pt_eta_phi_m).at(0), self().mpairRef()->weight);
        h_pt_closest_hadr_pt_furthest_hadr_ratio[isign][idphi]->Fill((self().mpairRef()->m2_last_hf_hadron_prt_pt_eta_phi_m).at(0) / (*m2_first_hf_hadron_prt_pt_eta_phi_m).at(0), self().mpairRef()->weight);
        h_pt_hadr_hq_ratio[isign][idphi]->Fill((self().mpairRef()->m1_last_hf_hadron_prt_pt_eta_phi_m).at(0) / (self().mpairRef()->m1_first_hq_ancestor_pt_eta_phi_m).at(0), self().mpairRef()->weight);
        h_pt_hadr_hq_ratio[isign][idphi]->Fill((self().mpairRef()->m2_last_hf_hadron_prt_pt_eta_phi_m).at(0) / (self().mpairRef()->m2_first_hq_ancestor_pt_eta_phi_m).at(0), self().mpairRef()->weight);
        float muon_hadr_dphi = self().mpairRef()->m1.truth_phi - (self().mpairRef()->m1_last_hf_hadron_prt_pt_eta_phi_m).at(2);
        muon_hadr_dphi = atan2(sin(muon_hadr_dphi),cos(muon_hadr_dphi));//fold muon_hadr_dphi to [-pi,pi]
        h_dphi_muon_closest_hadr[isign][idphi]->Fill(abs(muon_hadr_dphi), self().mpairRef()->weight);
        muon_hadr_dphi = self().mpairRef()->m2.truth_phi - (self().mpairRef()->m2_last_hf_hadron_prt_pt_eta_phi_m).at(2);
        muon_hadr_dphi = atan2(sin(muon_hadr_dphi),cos(muon_hadr_dphi));//fold muon_hadr_dphi to [-pi,pi]
        h_dphi_muon_closest_hadr[isign][idphi]->Fill(abs(muon_hadr_dphi), self().mpairRef()->weight);
    } catch(const std::out_of_range& e){
        if (self().debug_mode){
            std::cerr << "KinematicCorrPlots:: out_of_range exception caught: " << e.what() << std::endl;
            std::cerr << "Return without filling KinematicCorr plots" << e.what() << std::endl;
        }
        return;
    }
}

template <class PairT, class Derived>
void PowhegTruthExtras<PairT, Derived>::HistAdjustExtra(){
    if (!output_truth_hists) return;

    try{
        for (int ibin = 0; ibin < 3; ibin++){
            for (int isign = 0; isign < ParamsSet::nSigns; isign++){
                for (int jdphi = 0; jdphi < 2; jdphi++){
                    h_num_hard_scatt_out[isign][jdphi]->GetXaxis()->SetBinLabel(ibin+1,num_hard_scatt_out_labels.at(ibin).c_str());
                }
            }
        }

        for (int ibin = 0; ibin < nParentGroups; ibin++){
            for (int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
                for (int lphi = 0; lphi < 2; lphi++){
                    h_parent_groups[ksign][lphi]->GetXaxis()->SetBinLabel(ibin+1,parentGroupLabels.at(ibin).c_str());
                    h_parent_groups[ksign][lphi]->GetYaxis()->SetBinLabel(ibin+1,parentGroupLabels.at(ibin).c_str());
                }
            }
        }


        for (int lphi = 0; lphi < 2; lphi++){
            for (int isign = 0; isign < ParamsSet::nSigns; isign++){
                for (int ibin = 0; ibin < 2; ibin++){
                    h_QQ_both_from_Q_same_ancestors[isign][lphi]->GetXaxis()->SetBinLabel(ibin+1,samePrtsLabels.at(ibin).c_str());
                }
                for (int ibin = 0; ibin < nAncestorGroups; ibin++){
                    h_QQ_both_from_Q_ancestor_sp[isign][lphi]->GetXaxis()->SetBinLabel(ibin+1,ancestor_labels.at(ibin).c_str());
                    h_QQ_both_from_Q_ancestor_dp[isign][lphi]->GetXaxis()->SetBinLabel(ibin+1,ancestor_labels.at(ibin).c_str());
                    h_QQ_both_from_Q_ancestor_dp[isign][lphi]->GetYaxis()->SetBinLabel(ibin+1,ancestor_labels.at(ibin).c_str());
                }
            }
        }
    } catch(const std::out_of_range& e){
        std::cout << "PowhegTruthExtras::HistAdjustExtra: out_of_range error caught when setting bin labels!" << std::endl;
        std::cout << "Continue program without setting labels correctly!" << std::endl;
    }
}

template <class PairT, class Derived>
void PowhegTruthExtras<PairT, Derived>::FinalizeExtra(){    
    delete qqpair;

    std::cout << "The total integral of events where the hardest scattering is 2-to-2 is: " << crossx_2_to_2 << std::endl;
    std::cout << "The total integral of events where the hardest scattering is 2-to-3 is: " << crossx_2_to_3 << std::endl;
    std::cout << "Total percentage in crossx of events where the relevant hard scattering for a same-ancestor muon pair is not the hardest scattering is: " << crossx_relevant_hard_isnt_hardest / (crossx_2_to_2 + crossx_2_to_3) << std::endl;
    std::cout << "Total percentage in crossx of skipped HF events (not filled in HF-muon-pair-ancestor-category trees): " << skipped_event_crossx / (crossx_2_to_2 + crossx_2_to_3) << std::endl;

    delete m1_history;
    delete m2_history;
    delete m1_history_particle;
    delete m2_history_particle;
    delete m1_first_hf_hadron_prt_pt_eta_phi_m;
    delete m2_first_hf_hadron_prt_pt_eta_phi_m;

    m_unspecified_parent_file->close();
    delete m_unspecified_parent_file;

    if (print_prt_history){
        if (self().mcModeRef() == "bb"){
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
        if (self().mcModeRef() == "bb"){
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
}
