#ifndef __MCNTupleFirstPass_h__
#define __MCNTupleFirstPass_h__

#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <string.h>
#include <fstream>
// #include <cmath>  
#include "MuonPairPowheg.h"
#include "TruthQQPair.h"
#include "ParamsSet.h"
#include "/usatlas/u/yuhanguo/workarea/pythia/struct_particle.h"

// #include "vector"
#include "TH1D.h"
#include "TH2D.h"

class MCNTupleFirstPass{
    // Read through the N-tuple, apply appropriate cuts
    // Then fill in histograms and/or output trees
    // mode = 1: output single-muon information into a TTree
    // mode = 2: output muon-pair information into a TTree
    // NOW ONLY HAVE MODE = 1, 2


private:
// --------------------- general settings ---------------------------

    double crossx_cut;
    double filter_effcy;
    double filter_effcy_bb = 0.003;
    double filter_effcy_cc = 0.001108;
    // static const int nMCmodes = 2;

    ParamsSet pms;

    std::string mcdir;

// --------------------- input files & trees & data for setting branches---------------------------
  
    TChain          *fChain;   //!pointer to the analyzed TTree or TChain
  
    // Declaration of leaf types
    // UInt_t          RunNumber;
    float                 Q;
    std::vector<float>   *EventWeights           =nullptr;

    std::vector<int>* truth_id = nullptr;
    std::vector<int>* truth_barcode = nullptr;
    std::vector<int>* truth_qual = nullptr;
    std::vector<float>* truth_pt = nullptr;
    std::vector<float>* truth_eta = nullptr;
    std::vector<float>* truth_phi = nullptr;
    std::vector<float>* truth_m = nullptr;
    std::vector<std::vector<int>>* truth_parents = nullptr;
    std::vector<std::vector<int>>* truth_children = nullptr;
    // std::vector<std::vector<int>>* truth_muon_parent_ids = nullptr;
    // std::vector<std::vector<int>>* truth_muon_parent_bars = nullptr;

    std::vector<float>   *muon_pair_muon1_pt           =nullptr;
    std::vector<float>   *muon_pair_muon1_eta       =nullptr;
    std::vector<float>   *muon_pair_muon1_phi          =nullptr;
    std::vector<int>     *muon_pair_muon1_ch         =nullptr;
    std::vector<int>     *muon_pair_muon1_bar         =nullptr;

    std::vector<float>   *muon_pair_muon2_pt       =nullptr;
    std::vector<float>   *muon_pair_muon2_eta          =nullptr;
    std::vector<float>   *muon_pair_muon2_phi         =nullptr;
    std::vector<int>     *muon_pair_muon2_ch         =nullptr;
    std::vector<int>     *muon_pair_muon2_bar         =nullptr;

    std::vector<float>   *truth_mupair_asym         =nullptr;
    std::vector<float>   *truth_mupair_acop         =nullptr;
    std::vector<float>   *truth_mupair_pt         =nullptr;
    std::vector<float>   *truth_mupair_y         =nullptr;
    // std::vector<float>   *truth_mupair_phi         =nullptr;
    std::vector<float>   *truth_mupair_m         =nullptr;
    
// --------------------- temporary variables (muon, muonpair objects, vectors, etc.) ---------------------------
  
    Long64_t nentries;
    std::string sub_dir;

    Muon* tempmuon = nullptr;
    MuonPairPowheg* mpair = nullptr;
    TruthQQPair* qqpair = nullptr;
    std::vector<int> resonance_tagged_muon_index_list {};

    int cur_m1_earliest_parent_barcode;
    int cur_m2_earliest_parent_barcode;

    // bool same_ancestors;
    bool m1_ancestor_is_incoming;
    bool m2_ancestor_is_incoming;

    bool m1_from_J_psi;
    bool m2_from_J_psi;

    bool m1_from_Upsilon;
    bool m2_from_Upsilon;

    float crossx_2_to_2 = 0.;
    float crossx_2_to_3 = 0.;
    float crossx_relevant_hard_isnt_hardest = 0.;

    int nmuonpairs;

    double skipped_event_crossx = 0;
    bool skip_event = false;
    // boolean to tag whether 
    bool pre_return = false;
    
    std::vector<std::vector<int>>* single_gluon_history;
    std::vector<std::vector<int>>* m1_history;
    std::vector<std::vector<int>>* m2_history;
    std::vector<std::vector<Particle>>* m1_history_particle;
    std::vector<std::vector<Particle>>* m2_history_particle;


    std::vector<int> cur_m1_ancestor_ids;
    std::vector<int> cur_m2_ancestor_ids;
    std::vector<int> cur_m1_ancestor_bars;
    std::vector<int> cur_m2_ancestor_bars;

    std::vector<int> m1_multi_hf_quark_ids;
    std::vector<int> m2_multi_hf_quark_ids;

    bool m1_c_tag;
    bool m2_c_tag;
    bool m1_osc;
    bool m2_osc;
    int m1_earliest_parent_id;
    int m2_earliest_parent_id;

    // std::vector<float>* m1_last_b_hadron_prt_pt_eta_phi_m;
    // std::vector<float>* m2_last_b_hadron_prt_pt_eta_phi_m;
    // std::vector<float>* m1_last_hf_hadron_prt_pt_eta_phi_m;
    // std::vector<float>* m2_last_hf_hadron_prt_pt_eta_phi_m;
    std::vector<float>* m1_first_hf_hadron_prt_pt_eta_phi_m;
    std::vector<float>* m2_first_hf_hadron_prt_pt_eta_phi_m;
    // std::vector<float>* m1_hq_ancestor_pt_eta_phi_m;
    // std::vector<float>* m2_hq_ancestor_pt_eta_phi_m;

  
// --------------------- output file, histograms & trees ---------------------------
  
    TFile *m_outfile = nullptr;
    std::ofstream* m_unspecified_parent_file = nullptr;
    std::ofstream* m_b_parent_file[ParamsSet::nSigns][2];
    std::ofstream* m_c_parent_file[ParamsSet::nSigns][2];
    std::ofstream* m_cc_ss_small_dphi_file = nullptr;
    std::ofstream* m_bb_ss_near_file = nullptr;
    std::ofstream* m_bb_ss_away_file = nullptr;
    std::ofstream* m_bb_op_near_one_b_one_btoc_others_file = nullptr;

    TTree* meta_tree;
    TTree* muonPairOutTree[ParamsSet::nSigns];
    static const int nAncestorGroups = 4;
    TTree* QQPairOutTree[ParamsSet::nSigns][2][nAncestorGroups];
    TTree* muonPairOutTreeBinned[ParamsSet::ndRselcs][ParamsSet::nSigns];

    TH1D* h_cutAcceptance[ParamsSet::nSigns];
    // TH1D* h_numParents;
    TH1D* h_numMuonPairs;
    TH2D* h_ptlead_pair_pt[ParamsSet::nSigns][2];
    TH2D* h_parent_groups[ParamsSet::nSigns][2];

    TH1D* h_QQ_DR[ParamsSet::nSigns][2][4];
    TH1D* h_QQ_Dphi[ParamsSet::nSigns][2][4];
    TH1D* h_QQ_minv[ParamsSet::nSigns][2][4];
    TH1D* h_QQ_pair_pt_ptlead_ratio[ParamsSet::nSigns][2][4];
    TH1D* h_QQ_pt_avg[ParamsSet::nSigns][2][4];
    TH1D* h_QQ_asym[ParamsSet::nSigns][2][4];
    TH1D* h_QQ_minv_mHard_ratio[ParamsSet::nSigns][2][4];

    TH2D* h_QQ_ptlead_pair_pt[ParamsSet::nSigns][2][4];
    TH2D* h_QQ_pt1_pt2[ParamsSet::nSigns][2][4];
    TH2D* h_QQ_Deta_Dphi[ParamsSet::nSigns][2][4];
    TH2D* h_QQ_eta1_eta2[ParamsSet::nSigns][2][4];
    TH2D* h_QQ_minv_pair_pt[ParamsSet::nSigns][2][4];
    TH2D* h_QQ_minv_Dphi[ParamsSet::nSigns][2][4];

    static const int npTbins = 3;
    TH1D* h_crossx;
    TH1D* h_crossx_pt_binned[npTbins];

    // TH1D* h_bb_op_one_b_one_btoc[2];

    TH1D* h_QQ_both_from_Q_same_ancestors[ParamsSet::nSigns][2];
    TH1D* h_QQ_both_from_Q_ancestor_sp[ParamsSet::nSigns][2];
    TH2D* h_QQ_both_from_Q_ancestor_dp[ParamsSet::nSigns][2];

    TH1D* h_num_hard_scatt_out[ParamsSet::nSigns][2];
    TH1D* h_pt_muon_pt_closest_hadr_ratio[ParamsSet::nSigns][2];
    TH1D* h_pt_closest_hadr_pt_furthest_hadr_ratio[ParamsSet::nSigns][2];
    TH1D* h_pt_hadr_hq_ratio[ParamsSet::nSigns][2];
    TH1D* h_dphi_muon_closest_hadr[ParamsSet::nSigns][2];

    // std::vector<std::string> parentGroupLabels = {"direct b","b to c","direct c","s/light","single photon", "Drell-Yan"};
    std::vector<std::string> parentGroupLabels = {"direct b","b to c","direct c","s/light","single photon"};
    int nParentGroups = parentGroupLabels.size();    
    
    std::vector<std::string> ancestor_labels = {"gg", "gq", "single g", "q qbar"};
    
    // static const int numCuts = 5;
    static const int numCuts = 6;
    std::string cutLabels[numCuts] = {"no cut", "prev resonance tag", "muon eta", "muon pT", "resonance", "photoproduction"};

    std::vector<std::string> samePrtsLabels = {"Same Parents", "Different Parents"};
    std::vector<std::string> bb_op_one_b_one_btoc_labels = {"Same b", "Involve osc(s)", "From different ancestors", "Others"};

    // std::vector<std::string> osc_labels = {"0 osc, one b one b-to-c", "0 osc, others", "1 osc, one b one b-to-c", "1 osc, regular", "2 oscs, one b one b-to-c", "2 oscs, others"};
    std::vector<std::string> osc_labels = {"0 osc, one c-tag", "0 osc, others(*)", "1 osc, one c-tag(*)", "1 osc, regular", "2 oscs, one c-tag", "2 oscs, others(*)"};
    std::vector<std::string> num_hard_scatt_out_labels = {"2","3","more"};

// --------------------- class methods ---------------------------
  
    void InitInput();
    void InitTempVariables();
    void InitOutput();
    void ProcessData();
    bool PassCuts();
    bool IsResonance();
    bool IsPhotoProduction();
    int  ParentGrouping(std::vector<int>& parent_ids, bool c_tag, bool prev_is_lepton);
    void GetPtEtaPhiMFromBarcode(int barcode, std::vector<float>* pt_eta_phi_m);
    int  UpdateCurParents(bool isMuon1, std::vector<int>& cur_prt_bars, std::vector<int>& cur_prt_ids, int hf_quark_index = -1);
    int  FindHeavyQuarks(std::vector<int>& cur_prt_ids, int quark_type, bool isMuon1, int hadron_child_id = 0);
    void SingleMuonAncestorTracing(bool isMuon1);
    int  AncestorGrouping(std::vector<int>& ancestor_ids);
    void HardScatteringAnalysis(std::vector<int>& ancestor_bars, std::vector<int>& ancestor_ids, int sign_dphi_mode, int ancestor_grp);
    void PrintHistory(std::ostream* f, bool print_single, bool muon1_sameancestor);
    // void SameSignSameAncestorsAnalysis(bool near_side, bool one_b_one_btoc, bool print_history = false);
    void KinematicCorrPlots(int isign, int jdphi);
    void MuonPairTagsReinit();
    void CheckIfFromSameB();
    void MuonPairAncestorTracing();
    void FillMuonPairTree();
    void HistAdjust();
    void Finalize();
  

public :
    // int mode = 2;
    bool is_full_sample = true;
    int full_sample_batch_num;
    std::string mc_mode = "mc_truth_cc";
    bool print_prt_history = false;
    bool print_specific_prt_history = false;
    MCNTupleFirstPass(){
        crossx_cut = 5 * pow(10,8);
    }
    ~MCNTupleFirstPass(){}
    void Run();
};


//initialize the TChain
void MCNTupleFirstPass::InitInput(){

    mcdir = (is_full_sample)? "/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample/" : "/usatlas/u/yuhanguo/usatlasdata/athena/runMCV2/";
    fChain = new TChain("HeavyIonD3PD","HeavyIonD3PD");
    fChain->SetMakeClass(1);
    if (is_full_sample){
        sub_dir = (mc_mode == "mc_truth_bb")? "bb_full_sample/" : "cc_full_sample/";
        for (int ifile = 5 * (full_sample_batch_num - 1); ifile < 5 * full_sample_batch_num; ifile++){
            char filename[100];
            std::sprintf(filename, "%s%s%s_%02d.root", mcdir.c_str(), sub_dir.c_str(), mc_mode.c_str(), ifile+1);
            std::cout << filename << std::endl;
            std::ifstream in_file(filename);
            
            if (!in_file.good()){
                std::cout << "Warning: File " << filename << " not found. Skip.\n";
                continue; // skip this file
            }
            fChain->Add(filename);
        }
    }else{
        fChain->Add(Form("%s%s.root", mcdir.c_str(), mc_mode.c_str()));
    }
    
    cout << "nentries: " << fChain->GetEntries() << endl;
    // MC muons have no quality, d0 or z0 recorded
    fChain->SetBranchAddress("truth_id"                   , &truth_id);
    fChain->SetBranchAddress("truth_barcode"              , &truth_barcode);
    fChain->SetBranchAddress("truth_qual"                 , &truth_qual);
    fChain->SetBranchAddress("truth_pt"                 , &truth_pt);
    fChain->SetBranchAddress("truth_eta"                 , &truth_eta);
    fChain->SetBranchAddress("truth_phi"                 , &truth_phi);
    fChain->SetBranchAddress("truth_m"                 , &truth_m);
    fChain->SetBranchAddress("truth_parents"              , &truth_parents);
    fChain->SetBranchAddress("truth_children"              , &truth_children);
    fChain->SetBranchAddress("EventWeights"               , &EventWeights);
    fChain->SetBranchAddress("Q"               , &Q);

    fChain->SetBranchAddress("truth_mupair_pt1"           , &muon_pair_muon1_pt);
    fChain->SetBranchAddress("truth_mupair_eta1"          , &muon_pair_muon1_eta);
    fChain->SetBranchAddress("truth_mupair_phi1"          , &muon_pair_muon1_phi);
    fChain->SetBranchAddress("truth_mupair_ch1"           , &muon_pair_muon1_ch);
    fChain->SetBranchAddress("truth_mupair_bar1"          , &muon_pair_muon1_bar);

    fChain->SetBranchAddress("truth_mupair_pt2"           , &muon_pair_muon2_pt);
    fChain->SetBranchAddress("truth_mupair_eta2"          , &muon_pair_muon2_eta);
    fChain->SetBranchAddress("truth_mupair_phi2"          , &muon_pair_muon2_phi);
    fChain->SetBranchAddress("truth_mupair_ch2"           , &muon_pair_muon2_ch);
    fChain->SetBranchAddress("truth_mupair_bar2"          , &muon_pair_muon2_bar);

    fChain->SetBranchAddress("truth_mupair_asym"          , &truth_mupair_asym);
    fChain->SetBranchAddress("truth_mupair_acop"          , &truth_mupair_acop);
    fChain->SetBranchAddress("truth_mupair_pt"            , &truth_mupair_pt);
    fChain->SetBranchAddress("truth_mupair_y"             , &truth_mupair_y);
    // fChain->SetBranchAddress("truth_mupair_phi"           , &truth_mupair_phi);
    fChain->SetBranchAddress("truth_mupair_m"             , &truth_mupair_m);


    //SetBranch Status
    fChain->SetBranchStatus("*"                               ,0);//switch off all branches, then enable just the ones that we need
  
    fChain->SetBranchStatus("truth_id"           ,1);
    fChain->SetBranchStatus("truth_barcode"              ,1);
    fChain->SetBranchStatus("truth_qual"             ,1);
    fChain->SetBranchStatus("truth_pt"             ,1);
    fChain->SetBranchStatus("truth_eta"             ,1);
    fChain->SetBranchStatus("truth_phi"             ,1);
    fChain->SetBranchStatus("truth_m"             ,1);
    fChain->SetBranchStatus("truth_parents"             ,1);
    fChain->SetBranchStatus("truth_children"             ,1);
    fChain->SetBranchStatus("EventWeights"             ,1);
    fChain->SetBranchStatus("Q"             ,1);

    fChain->SetBranchStatus("truth_mupair_pt1"           ,1);
    fChain->SetBranchStatus("truth_mupair_eta1"              ,1);
    fChain->SetBranchStatus("truth_mupair_phi1"             ,1);
    fChain->SetBranchStatus("truth_mupair_ch1"             ,1);
    fChain->SetBranchStatus("truth_mupair_bar1"             ,1);

    fChain->SetBranchStatus("truth_mupair_pt2"             ,1);
    fChain->SetBranchStatus("truth_mupair_eta2"         ,1);
    fChain->SetBranchStatus("truth_mupair_phi2"              ,1);
    fChain->SetBranchStatus("truth_mupair_ch2"              ,1);
    fChain->SetBranchStatus("truth_mupair_bar2"              ,1);

    fChain->SetBranchStatus("truth_mupair_asym"              ,1);
    fChain->SetBranchStatus("truth_mupair_acop"           ,1);
    fChain->SetBranchStatus("truth_mupair_pt"             ,1);
    fChain->SetBranchStatus("truth_mupair_y"             ,1);
    // fChain->SetBranchStatus("truth_mupair_phi"         ,1);
    fChain->SetBranchStatus("truth_mupair_m"              ,1);   
}

void MCNTupleFirstPass::InitTempVariables(){
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


void MCNTupleFirstPass::InitOutput(){

    if (print_specific_prt_history){
        if (mc_mode == "mc_truth_cc"){

            m_cc_ss_small_dphi_file = new std::ofstream(Form("%scc_ss_small_dphi.txt", mcdir.c_str()));
            *m_cc_ss_small_dphi_file << "Event#\tm1-grp\tm2-grp" << std::endl;
        }else{
            m_bb_ss_near_file = new std::ofstream(Form("%sbb_ss_near.txt", mcdir.c_str()));
            m_bb_ss_away_file = new std::ofstream(Form("%sbb_ss_away.txt", mcdir.c_str()));
            m_bb_op_near_one_b_one_btoc_others_file = new std::ofstream(Form("%sbb_op_near_one_b_one_btoc_others.txt", mcdir.c_str()));
        }
    }

    std::string sign_labels[ParamsSet::nSigns] = {"same_sign", "op_sign"};
    std::string signs[ParamsSet::nSigns] = {"_ss", "_op"};
    std::string dphis[2] = {"_near", "_away"};
    std::string ancestor_grps[nAncestorGroups] = {"_gg", "_qg","_single_g","_qq"};

    if (mc_mode == "mc_truth_bb"){
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

    // m_outfile=new TFile(Form("%smuon_pairs_%s.root", mcdir.c_str(), mc_mode.c_str()),"recreate");
    if (is_full_sample){
        std::string sub_dir = (mc_mode == "mc_truth_bb")? "bb_full_sample/" : "cc_full_sample/";
        m_outfile=new TFile(Form("%s%smuon_pairs_%s_%d-%d.root", mcdir.c_str(), sub_dir.c_str(), mc_mode.c_str(), 5 * full_sample_batch_num - 4, 5 * full_sample_batch_num),"recreate");
    }else{
        m_outfile=new TFile(Form("%smuon_pairs_%s.root", mcdir.c_str(), mc_mode.c_str()),"recreate");
    }
    
    meta_tree = new TTree("meta_tree","meta_tree");
    meta_tree->Branch("nentries_before_cuts", &nentries, "nentries_before_cuts/L");

    for (unsigned int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
        muonPairOutTree[ksign] = new TTree(Form("muon_pair_tree_sign%u",ksign+1),Form("all muon pairs, sign%u",ksign+1));
        muonPairOutTree[ksign]->Branch("MuonPairObj",&mpair);
        for (unsigned int idr = 0; idr < ParamsSet::ndRselcs; idr++){
            muonPairOutTreeBinned[idr][ksign] = new TTree(Form("muon_pair_tree_dr%u_sign%u",idr+1,ksign+1),Form("all muon pairs, dr%u, sign%u",idr+1,ksign+1));
            muonPairOutTreeBinned[idr][ksign]->Branch("MuonPairObj",&mpair);
        }
    }

    for (unsigned int isign = 0; isign < ParamsSet::nSigns; isign++){
        for (int jdphi = 0; jdphi < 2; jdphi++){
            for (int kgrp = 0; kgrp < nAncestorGroups; kgrp++){
                QQPairOutTree[isign][jdphi][kgrp] = new TTree(Form("QQ_pair_tree%s%s%s",signs[isign].c_str(),dphis[jdphi].c_str(),ancestor_grps[kgrp].c_str()),Form("QQ pairs, muon %s%s%s",signs[isign].c_str(),dphis[jdphi].c_str(),ancestor_grps[kgrp].c_str()));
                QQPairOutTree[isign][jdphi][kgrp]->Branch("QQPairObj",&qqpair);
            }
        }
    }

    // h_numParents = new TH1D("h_numParents","h_numParents",3,0,3);
    h_numMuonPairs = new TH1D("h_numMuonPairs","h_numMuonPairs",6,0,6);

    for (int isign = 0; isign < ParamsSet::nSigns; isign++){
        h_cutAcceptance[isign] = new TH1D(Form("h_cutAcceptance_sign%d",isign+1),Form("h_cutAcceptance_sign%d",isign+1),numCuts,0,numCuts);
        for (int jdphi = 0; jdphi < 2; jdphi++){
            h_parent_groups[isign][jdphi] = new TH2D(Form("h_parent_groups_sign%d%s",isign+1, dphis[jdphi].c_str()),Form("h_parent_groups_sign%d%s",isign+1, dphis[jdphi].c_str()),nParentGroups,0,nParentGroups,nParentGroups,0,nParentGroups);
            h_ptlead_pair_pt[isign][jdphi] = new TH2D(Form("h_ptlead_pair_pt_sign%d%s",isign+1, dphis[jdphi].c_str()),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",pms.npairPT_bins,pms.pairPTBins[isign][2],pms.npt_bins,pms.pTBins);
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

    // static const int nweight_bins = 40;
    // float weight_logpow_bb = 0.0711;
    // float weight_logpow_cc = 0.1254;
    // float weight_max_bb = 5.74;
    // float weight_max_cc = 1.26;
    // double weight_bins_bb[nweight_bins+1];
    // double weight_bins_cc[nweight_bins+1];

    // for(int iweight = 0; iweight <= nweight_bins; iweight++){
    //     weight_bins_bb[iweight] = weight_max_bb * pow(10.0, ((float)(iweight - nweight_bins))*weight_logpow_bb);
    //     weight_bins_cc[iweight] = weight_max_cc * pow(10.0, ((float)(iweight - nweight_bins))*weight_logpow_cc);
    // }

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


    if (mc_mode == "mc_truth_bb"){
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

    for (int isign = 0; isign < ParamsSet::nSigns; isign++){
        h_QQ_both_from_Q_same_ancestors[isign][0] = new TH1D(Form("h_QQ_both_from_Q_same_ancestors_sign%d_near",isign+1),Form("h_QQ_both_from_Q_same_ancestors_sign%d_near",isign+1),2,0,2);
        h_QQ_both_from_Q_same_ancestors[isign][1] = new TH1D(Form("h_QQ_both_from_Q_same_ancestors_sign%d_away",isign+1),Form("h_QQ_both_from_Q_same_ancestors_sign%d_away",isign+1),2,0,2);
        h_QQ_both_from_Q_ancestor_sp[isign][0] = new TH1D(Form("h_QQ_both_from_Q_ancestor_sp_sign%d_near",isign+1),Form("h_QQ_both_from_Q_ancestor_sp_sign%d_near",isign+1),nAncestorGroups,0,nAncestorGroups);
        h_QQ_both_from_Q_ancestor_sp[isign][1] = new TH1D(Form("h_QQ_both_from_Q_ancestor_sp_sign%d_away",isign+1),Form("h_QQ_both_from_Q_ancestor_sp_sign%d_away",isign+1),nAncestorGroups,0,nAncestorGroups);
        h_QQ_both_from_Q_ancestor_dp[isign][0] = new TH2D(Form("h_QQ_both_from_Q_ancestor_dp_sign%d_near",isign+1),Form("h_QQ_both_from_Q_ancestor_dp_sign%d_near",isign+1),nAncestorGroups,0,nAncestorGroups,nAncestorGroups,0,nAncestorGroups);
        h_QQ_both_from_Q_ancestor_dp[isign][1] = new TH2D(Form("h_QQ_both_from_Q_ancestor_dp_sign%d_away",isign+1),Form("h_QQ_both_from_Q_ancestor_dp_sign%d_away",isign+1),nAncestorGroups,0,nAncestorGroups,nAncestorGroups,0,nAncestorGroups);
    }
}

void MCNTupleFirstPass::Finalize(){
    delete m1_history;
    delete m2_history;
    delete m1_history_particle;
    delete m2_history_particle;
    // delete m1_last_b_hadron_prt_pt_eta_phi_m;
    // delete m2_last_b_hadron_prt_pt_eta_phi_m;
    // delete m1_last_hf_hadron_prt_pt_eta_phi_m;
    // delete m2_last_hf_hadron_prt_pt_eta_phi_m;
    delete m1_first_hf_hadron_prt_pt_eta_phi_m;
    delete m2_first_hf_hadron_prt_pt_eta_phi_m;
    // delete m1_hq_ancestor_pt_eta_phi_m;
    // delete m2_hq_ancestor_pt_eta_phi_m;

    m_unspecified_parent_file->close();
    delete m_unspecified_parent_file;

    if (print_prt_history){
        if (mc_mode == "mc_truth_bb"){
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
        if (mc_mode == "mc_truth_bb"){
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
            delete h_ptlead_pair_pt[isign][jdphi];
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

#endif


