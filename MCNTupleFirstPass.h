#ifndef __MCNTupleFirstPass_h__
#define __MCNTupleFirstPass_h__

#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <string.h>
#include <fstream>
// #include <cmath>  
#include "MuonPairMC.h"
#include "ParamsSet.h"
#include "vector"
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

    float filter_effcy_bb = 0.0003774031;
    float filter_effcy_cc = 0.000005964574;
    // static const int nMCmodes = 2;

    ParamsSet pms;

    std::string mcdir = "/usatlas/u/yuhanguo/usatlasdata/athena/runMCV2/";

// --------------------- input files & trees & data for setting branches---------------------------
  
    TChain          *fChain;   //!pointer to the analyzed TTree or TChain
  
    // Declaration of leaf types
    UInt_t          RunNumber;

    std::vector<int>* truth_id = nullptr;
    std::vector<int>* truth_barcode = nullptr;
    std::vector<int>* truth_qual = nullptr;
    std::vector<std::vector<int>>* truth_parents = nullptr;
    std::vector<std::vector<int>>* truth_children = nullptr;
    // std::vector<std::vector<int>>* truth_muon_parent_ids = nullptr;
    // std::vector<std::vector<int>>* truth_muon_parent_bars = nullptr;
    std::vector<float>   *EventWeights           =nullptr;

    std::vector<float>   *muon_pair_muon1_pt           =nullptr;
    std::vector<float>   *muon_pair_muon1_eta       =nullptr;
    std::vector<float>   *muon_pair_muon1_phi          =nullptr;
    std::vector<float>     *muon_pair_muon1_ch         =nullptr;
    std::vector<float>     *muon_pair_muon1_bar         =nullptr;

    std::vector<float>   *muon_pair_muon2_pt       =nullptr;
    std::vector<float>   *muon_pair_muon2_eta          =nullptr;
    std::vector<float>   *muon_pair_muon2_phi         =nullptr;
    std::vector<float>     *muon_pair_muon2_ch         =nullptr;
    std::vector<float>     *muon_pair_muon2_bar         =nullptr;

    std::vector<float>   *truth_mupair_asym         =nullptr;
    std::vector<float>   *truth_mupair_acop         =nullptr;
    std::vector<float>   *truth_mupair_pt         =nullptr;
    std::vector<float>   *truth_mupair_y         =nullptr;
    // std::vector<float>   *truth_mupair_phi         =nullptr;
    std::vector<float>   *truth_mupair_m         =nullptr;
    
// --------------------- temporary variables (muon, muonpair objects, vectors, etc.) ---------------------------
  
    Muon* tempmuon = nullptr;
    MuonPairMC* mpair = nullptr;

    int cur_m1_earliest_parent_barcode;
    int cur_m2_earliest_parent_barcode;

    // bool same_ancestors;
    bool m1_ancestor_is_incoming;
    bool m2_ancestor_is_incoming;

    bool m1_from_J_psi;
    bool m2_from_J_psi;

    bool m1_from_Upsilon;
    bool m2_from_Upsilon;

    // int cur_m1_beyond_earliest_parent_barcode;
    // int cur_m2_beyond_earliest_parent_barcode;
    // std::vector<int> cur_m1_beyond_earliest_parent_ids;
    // std::vector<int> cur_m2_beyond_earliest_parent_ids;

    int nmuonpairs;
    
    std::vector<std::vector<int>> m1_history_before_hadrons;
    std::vector<std::vector<int>> m2_history_before_hadrons;

    std::vector<int> cur_m1_ancestor_ids;
    std::vector<int> cur_m2_ancestor_ids;
    std::vector<int> cur_m1_ancestor_bars;
    std::vector<int> cur_m2_ancestor_bars;

    std::vector<int> m1_multi_hf_quark_ids;
    std::vector<int> m2_multi_hf_quark_ids;

    // std::vector<int> m1_multi_hadronic_parents_ids;
    // std::vector<int> m2_multi_hadronic_parents_ids;
  
// --------------------- output file, histograms & trees ---------------------------
  
    TFile *m_outfile = nullptr;
    ofstream m_unspecified_parent_file;
    ofstream m_b_parent_file[ParamsSet::nSigns][2];
    ofstream m_c_parent_file[ParamsSet::nSigns][2];
    ofstream m_cc_ss_small_dphi_file;
    ofstream m_bb_ss_near_file;

    TTree* muonOutTree;
    TTree* muonPairOutTree[ParamsSet::nSigns];
    TTree* muonPairOutTreeBinned[ParamsSet::ndRselcs][ParamsSet::nSigns];

    TH1D* h_cutAcceptance[ParamsSet::nSigns];
    // TH1D* h_numParents;
    TH1D* h_numMuonPairs;
    TH2D* h_unweighted_parent_groups[ParamsSet::nSigns][2];
    TH2D* h_parent_groups[ParamsSet::nSigns][2];

    static const int npTbins = 3;
    TH1D* h_crossx;
    TH1D* h_crossx_pt_binned[npTbins];

    TH1D* h_unweighted_bb_op_one_b_one_btoc[2];
    TH1D* h_bb_op_one_b_one_btoc[2];

    TH1D* h_unweighted_bb_both_from_b_same_ancestors[ParamsSet::nSigns][2];
    TH1D* h_unweighted_cc_both_from_c_same_ancestors[ParamsSet::nSigns][2];
    TH1D* h_bb_both_from_b_same_ancestors[ParamsSet::nSigns][2];
    TH1D* h_cc_both_from_c_same_ancestors[ParamsSet::nSigns][2];

    TH1D* h_unweighted_bb_both_from_b_ancestor_sp[ParamsSet::nSigns][2];
    TH2D* h_unweighted_bb_both_from_b_ancestor_dp[ParamsSet::nSigns][2];
    TH1D* h_bb_both_from_b_ancestor_sp[ParamsSet::nSigns][2];
    TH2D* h_bb_both_from_b_ancestor_dp[ParamsSet::nSigns][2];

    TH1D* h_unweighted_cc_both_from_c_ancestor_sp[ParamsSet::nSigns][2];
    TH2D* h_unweighted_cc_both_from_c_ancestor_dp[ParamsSet::nSigns][2];
    TH1D* h_cc_both_from_c_ancestor_sp[ParamsSet::nSigns][2];
    TH2D* h_cc_both_from_c_ancestor_dp[ParamsSet::nSigns][2];

    TH1D* h_dphi_bb_op_near_both_from_b;
    TH1D* h_dphi_bb_op_near_one_b_one_btoc;
    TH1D* h_dphi_bb_ss_near;
    TH1D* h_dphi_bb_op_near;
    TH1D* h_dphi_bb_op_near_from_same_b;

    TH1D* h_bb_ss_near_involv_osc;
    TH1D* h_bb_ss_near_num_hard_scatt_out;

    TH2D* h_cc_ss_small_dphi_prt_gps;
    TH1D* h_cc_ss_small_dphi_same_ancestors;
    TH1D* h_cc_ss_small_dphi_sp;
    TH2D* h_cc_ss_small_dphi_dp;
    TH2D* h_cc_ss_plateau_prt_gps;
    TH1D* h_cc_ss_plateau_same_ancestors;
    TH1D* h_cc_ss_plateau_sp;
    TH2D* h_cc_ss_plateau_dp;

    const int nParentGroups = 5;
    std::vector<std::string> parentGroupLabels = {"direct b","b to c","direct c","s/light","photon"};
    
    const int nAncestorGroups = 5;
    std::vector<std::string> ancestor_labels = {"incoming", "gg", "gq", "single g", "q qbar"};
    
    static const int numCuts = 5;
    std::string cutLabels[numCuts] = {"no cut", "muon pT", "muon eta", "resonance", "photoproduction"};

    std::vector<std::string> samePrtsLabels = {"Same Parents", "Different Parents"};
    std::vector<std::string> bb_op_one_b_one_btoc_labels = {"Same b", "Involve osc(s)", "From different ancestors", "Others"};

    // std::vector<std::string> osc_labels = {"0 osc, one b one b-to-c", "0 osc, others", "1 osc, one b one b-to-c", "1 osc, regular", "2 oscs, one b one b-to-c", "2 oscs, others"};
    std::vector<std::string> osc_labels = {"0 osc, one c-tag", "0 osc, others", "1 osc, one c-tag", "1 osc, regular", "2 oscs, one c-tag", "2 oscs, others"};
    std::vector<std::string> num_hard_scatt_out_labels = {"2","3","more"};

// --------------------- class methods ---------------------------
  
    void InitInput();
    void InitOutput();
    void ProcessData();
    bool PassCuts();
    bool IsResonance();
    bool IsPhotoProduction();
    int  ParentGrouping(int parent_id, bool c_tag);
    int  UpdateCurParents(bool isMuon1, std::vector<int>& cur_prt_bars, std::vector<int>& cur_prt_ids, bool record_history = true, int hf_quark_index = -1);
    int  FindHeavyQuarks(std::vector<int>& cur_prt_ids, int quark_type, bool isMuon1, int hadron_child_id = 0);
    void FindSingleMuonParents(bool isMuon1);
    int  AncestorGrouping(std::vector<int>& ancestor_ids, bool sameprts = true);
    void HardScatteringAnalysis(std::vector<int>& ancestor_bars, std::vector<int>& ancestor_ids, bool sameprts = true);
    void PrintHistory(std::ofstream& f, bool same_ancestors, bool state_same_ancestors = false);
    void FillMuonPairParents();
    void FillSingleMuonTree();
    void FillMuonPairTree();
  

public :
    int mode = 2;
    std::string mc_mode = "mc_truth_cc";
    bool print_prt_history = true;
    MCNTupleFirstPass();
    ~MCNTupleFirstPass(){}
    void Run();
    float filter_effcy;
};


MCNTupleFirstPass::MCNTupleFirstPass(){
    // if(mode != 1 && mode != 2){
    //     std::cout<<"Error:: Mode has to equal 1 or 2; code is used for outputting muon / muon-pair trees only."<<std::endl;
    //     throw std::exception();
    // }
    if(mode != 2){
        std::cout << "Error:: Mode has to equal 2; MC truth code is used for outputting muon-pair trees only." << std::endl;
        throw std::exception();
    }
}

//initialize the TChain
void MCNTupleFirstPass::InitInput(){

    fChain = new TChain("HeavyIonD3PD","HeavyIonD3PD");
    fChain->SetMakeClass(1);
    fChain->Add(Form("%s%s.root", mcdir.c_str(), mc_mode.c_str()));
    
    // MC muons have no quality, d0 or z0 recorded
    fChain->SetBranchAddress("truth_id"                   , &truth_id);
    fChain->SetBranchAddress("truth_barcode"              , &truth_barcode);
    fChain->SetBranchAddress("truth_qual"                 , &truth_qual);
    fChain->SetBranchAddress("truth_parents"              , &truth_parents);
    fChain->SetBranchAddress("truth_children"              , &truth_children);
    fChain->SetBranchAddress("EventWeights"               , &EventWeights);

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
    fChain->SetBranchStatus("truth_parents"             ,1);
    fChain->SetBranchStatus("truth_children"             ,1);
    fChain->SetBranchStatus("EventWeights"             ,1);

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


void MCNTupleFirstPass::InitOutput(){

    if (mc_mode == "mc_truth_cc"){
        m_cc_ss_small_dphi_file.open(Form("%scc_ss_small_dphi.txt", mcdir.c_str()));
        m_cc_ss_small_dphi_file << "Event#\tm1-grp\tm2-grp" << std::endl;
    }else{
        m_bb_ss_near_file.open(Form("%sbb_ss_near.txt", mcdir.c_str()));
    }


    std::string sign_labels[ParamsSet::nSigns] = {"same_sign", "op_sign"};
    if (mc_mode == "mc_truth_bb"){
        m_unspecified_parent_file.open(mcdir + "unspecified_parents_bb.txt");
        if (print_prt_history){
            for (int isign = 0; isign < ParamsSet::nSigns; isign++){
                m_b_parent_file[isign][0].open(Form("%sb_parents_%s_near.txt", mcdir.c_str(), sign_labels[isign].c_str()));
                m_b_parent_file[isign][1].open(Form("%sb_parents_%s_away.txt", mcdir.c_str(), sign_labels[isign].c_str()));
            }
        }
    }else{
        m_unspecified_parent_file.open(mcdir + "unspecified_parents_cc.txt");
        if (print_prt_history){
            for (int isign = 0; isign < ParamsSet::nSigns; isign++){
                m_c_parent_file[isign][0].open(Form("%sc_parents_%s_near.txt", mcdir.c_str(), sign_labels[isign].c_str()));
                m_c_parent_file[isign][1].open(Form("%sc_parents_%s_away.txt", mcdir.c_str(), sign_labels[isign].c_str()));
            }
        }
    } 

    if (mode == 1){
        m_outfile=new TFile(Form("%ssingle_muon_trees_%s.root", mcdir.c_str(), mc_mode.c_str()),"recreate");
        muonOutTree = new TTree("muon_tree","all single muons");
        muonOutTree->Branch("MuonObj",&tempmuon);
    }
    else{ //output muon pair trees
        m_outfile=new TFile(Form("%smuon_pairs_%s.root", mcdir.c_str(), mc_mode.c_str()),"recreate");

        for (unsigned int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
            muonPairOutTree[ksign] = new TTree(Form("muon_pair_tree_sign%u",ksign+1),Form("all muon pairs, sign%u",ksign+1));
            muonPairOutTree[ksign]->Branch("MuonPairObj",&mpair);
            for (unsigned int idr = 0; idr < ParamsSet::ndRselcs; idr++){
                muonPairOutTreeBinned[idr][ksign] = new TTree(Form("muon_pair_tree_dr%u_sign%u",idr+1,ksign+1),Form("all muon pairs, dr%u, sign%u",idr+1,ksign+1));
                muonPairOutTreeBinned[idr][ksign]->Branch("MuonPairObj",&mpair);
            }
        }
    }

    // h_numParents = new TH1D("h_numParents","h_numParents",3,0,3);
    h_numMuonPairs = new TH1D("h_numMuonPairs","h_numMuonPairs",6,0,6);

    static const int nweight_bins = 40;
    float weight_logpow_bb = 0.0711;
    float weight_logpow_cc = 0.1254;
    float weight_max_bb = 5.74;
    float weight_max_cc = 1.26;
    // std::vector<double> weight_bins_bb;
    // std::vector<double> weight_bins_cc;
    // weight_bins_bb.reserve(nweight_bins+1);
    // weight_bins_cc.reserve(nweight_bins+1);
    double weight_bins_bb[nweight_bins+1];
    double weight_bins_cc[nweight_bins+1];

    for(int iweight = 0; iweight <= nweight_bins; iweight++){
        weight_bins_bb[iweight] = weight_max_bb * pow(10.0, ((float)(iweight - nweight_bins))*weight_logpow_bb);
        weight_bins_cc[iweight] = weight_max_cc * pow(10.0, ((float)(iweight - nweight_bins))*weight_logpow_cc);
    }

    // std::vector<double> weight_bins = (mc_mode == "mc_truth_bb")? weight_bins_bb : weight_bins_cc;    

    if (mc_mode == "mc_truth_bb"){
        h_crossx = new TH1D("h_crossx","h_crossx",nweight_bins,weight_bins_bb);
        for (int ipt = 0; ipt < npTbins; ipt++){
            h_crossx_pt_binned[ipt] = new TH1D(Form("h_crossx_pt%d",ipt+1),Form("h_crossx_pt%d",ipt+1),nweight_bins,weight_bins_bb);
        }
    }else{
        h_crossx = new TH1D("h_crossx","h_crossx",nweight_bins,weight_bins_cc);
        for (int ipt = 0; ipt < npTbins; ipt++){
            h_crossx_pt_binned[ipt] = new TH1D(Form("h_crossx_pt%d",ipt+1),Form("h_crossx_pt%d",ipt+1),nweight_bins,weight_bins_cc);
        }
    }

    for (int isign = 0; isign < ParamsSet::nSigns; isign++){
        h_cutAcceptance[isign] = new TH1D(Form("h_cutAcceptance_sign%d",isign+1),Form("h_cutAcceptance_sign%d",isign+1),numCuts,0,numCuts);
        h_unweighted_parent_groups[isign][0] = new TH2D(Form("h_unweighted_parent_groups_sign%d_near",isign+1),Form("h_unweighted_parent_groups_sign%d_near",isign+1),nParentGroups,0,nParentGroups,nParentGroups,0,nParentGroups);
        h_unweighted_parent_groups[isign][1] = new TH2D(Form("h_unweighted_parent_groups_sign%d_away",isign+1),Form("h_unweighted_parent_groups_sign%d_away",isign+1),nParentGroups,0,nParentGroups,nParentGroups,0,nParentGroups);
        h_parent_groups[isign][0] = new TH2D(Form("h_parent_groups_sign%d_near",isign+1),Form("h_parent_groups_sign%d_near",isign+1),nParentGroups,0,nParentGroups,nParentGroups,0,nParentGroups);
        h_parent_groups[isign][1] = new TH2D(Form("h_parent_groups_sign%d_away",isign+1),Form("h_parent_groups_sign%d_away",isign+1),nParentGroups,0,nParentGroups,nParentGroups,0,nParentGroups);
    }

    if (mc_mode == "mc_truth_bb"){

        h_unweighted_bb_op_one_b_one_btoc[0] = new TH1D("h_unweighted_bb_op_one_b_one_btoc_near","h_unweighted_bb_op_one_b_one_btoc_near",3,0,3);
        h_unweighted_bb_op_one_b_one_btoc[1] = new TH1D("h_unweighted_bb_op_one_b_one_btoc_away","h_unweighted_bb_op_one_b_one_btoc_away",3,0,3);
        h_bb_op_one_b_one_btoc[0] = new TH1D("h_bb_op_one_b_one_btoc_near","h_bb_op_one_b_one_btoc_near",4,0,4);
        h_bb_op_one_b_one_btoc[1] = new TH1D("h_bb_op_one_b_one_btoc_away","h_bb_op_one_b_one_btoc_away",4,0,4);
    
        h_bb_ss_near_involv_osc = new TH1D("h_bb_ss_near_involv_osc","h_bb_ss_near_involv_osc",6,0,6);
        h_bb_ss_near_num_hard_scatt_out = new TH1D("h_bb_ss_near_num_hard_scatt_out","h_bb_ss_near_num_hard_scatt_out",3,2,5);

        h_dphi_bb_ss_near = new TH1D("h_dphi_bb_ss_near",";#Delta#phi;d#sigma/d#Delta#phi", 32,-pms.PI/2.,pms.PI/2.);
        h_dphi_bb_op_near = new TH1D("h_dphi_bb_op_near",";#Delta#phi;d#sigma/d#Delta#phi", 32,-pms.PI/2.,pms.PI/2.);
        h_dphi_bb_op_near_from_same_b = new TH1D("h_dphi_bb_op_near_from_same_b",";#Delta#phi;d#sigma/d#Delta#phi", 32,-pms.PI/2.,pms.PI/2.);
        
        h_dphi_bb_op_near_both_from_b = new TH1D("h_dphi_bb_op_near_both_from_b",";#Delta#phi;d#sigma/d#Delta#phi", 32,-pms.PI/2.,pms.PI/2.);
        h_dphi_bb_op_near_one_b_one_btoc = new TH1D("h_dphi_bb_op_near_one_b_one_btoc",";#Delta#phi;d#sigma/d#Delta#phi", 32,-pms.PI/2.,pms.PI/2.);

        for (int isign = 0; isign < ParamsSet::nSigns; isign++){
            h_unweighted_bb_both_from_b_same_ancestors[isign][0] = new TH1D(Form("h_unweighted_bb_both_from_b_same_ancestors_sign%d_near",isign+1),Form("h_unweighted_bb_both_from_b_same_ancestors_sign%d_near",isign+1),2,0,2);
            h_unweighted_bb_both_from_b_same_ancestors[isign][1] = new TH1D(Form("h_unweighted_bb_both_from_b_same_ancestors_sign%d_away",isign+1),Form("h_unweighted_bb_both_from_b_same_ancestors_sign%d_away",isign+1),2,0,2);
            h_bb_both_from_b_same_ancestors[isign][0] = new TH1D(Form("h_bb_both_from_b_same_ancestors_sign%d_near",isign+1),Form("h_bb_both_from_b_same_ancestors_sign%d_near",isign+1),2,0,2);
            h_bb_both_from_b_same_ancestors[isign][1] = new TH1D(Form("h_bb_both_from_b_same_ancestors_sign%d_away",isign+1),Form("h_bb_both_from_b_same_ancestors_sign%d_away",isign+1),2,0,2);
            h_unweighted_bb_both_from_b_ancestor_sp[isign][0] = new TH1D(Form("h_unweighted_bb_both_from_b_ancestor_sp_sign%d_near",isign+1),Form("h_unweighted_bb_both_from_b_ancestor_sp_sign%d_near",isign+1),nAncestorGroups,0,nAncestorGroups);
            h_unweighted_bb_both_from_b_ancestor_sp[isign][1] = new TH1D(Form("h_unweighted_bb_both_from_b_ancestor_sp_sign%d_away",isign+1),Form("h_unweighted_bb_both_from_b_ancestor_sp_sign%d_away",isign+1),nAncestorGroups,0,nAncestorGroups);
            h_unweighted_bb_both_from_b_ancestor_dp[isign][0] = new TH2D(Form("h_unweighted_bb_both_from_b_ancestor_dp_sign%d_near",isign+1),Form("h_unweighted_bb_both_from_b_ancestor_dp_sign%d_near",isign+1),nAncestorGroups,0,nAncestorGroups,nAncestorGroups,0,nAncestorGroups);
            h_unweighted_bb_both_from_b_ancestor_dp[isign][1] = new TH2D(Form("h_unweighted_bb_both_from_b_ancestor_dp_sign%d_away",isign+1),Form("h_unweighted_bb_both_from_b_ancestor_dp_sign%d_away",isign+1),nAncestorGroups,0,nAncestorGroups,nAncestorGroups,0,nAncestorGroups);
            h_bb_both_from_b_ancestor_sp[isign][0] = new TH1D(Form("h_bb_both_from_b_ancestor_sp_sign%d_near",isign+1),Form("h_bb_both_from_b_ancestor_sp_sign%d_near",isign+1),nAncestorGroups,0,nAncestorGroups);
            h_bb_both_from_b_ancestor_sp[isign][1] = new TH1D(Form("h_bb_both_from_b_ancestor_sp_sign%d_away",isign+1),Form("h_bb_both_from_b_ancestor_sp_sign%d_away",isign+1),nAncestorGroups,0,nAncestorGroups);
            h_bb_both_from_b_ancestor_dp[isign][0] = new TH2D(Form("h_bb_both_from_b_ancestor_dp_sign%d_near",isign+1),Form("h_bb_both_from_b_ancestor_dp_sign%d_near",isign+1),nAncestorGroups,0,nAncestorGroups,nAncestorGroups,0,nAncestorGroups);
            h_bb_both_from_b_ancestor_dp[isign][1] = new TH2D(Form("h_bb_both_from_b_ancestor_dp_sign%d_away",isign+1),Form("h_bb_both_from_b_ancestor_dp_sign%d_away",isign+1),nAncestorGroups,0,nAncestorGroups,nAncestorGroups,0,nAncestorGroups);
        }
    }else if (mc_mode == "mc_truth_cc"){
        h_cc_ss_small_dphi_prt_gps = new TH2D("h_cc_ss_small_dphi_prt_gps","h_cc_ss_small_dphi_prt_gps",nParentGroups,0,nParentGroups,nParentGroups,0,nParentGroups);
        h_cc_ss_small_dphi_same_ancestors = new TH1D("h_cc_ss_small_dphi_same_ancestors","h_cc_ss_small_dphi_same_ancestors",2,0,2);
        h_cc_ss_small_dphi_sp = new TH1D("h_cc_ss_small_dphi_sp","h_cc_ss_small_dphi_sp",nAncestorGroups,0,nAncestorGroups);
        h_cc_ss_small_dphi_dp = new TH2D("h_cc_ss_small_dphi_dp","h_cc_ss_small_dphi_dp",nAncestorGroups,0,nAncestorGroups,nAncestorGroups,0,nAncestorGroups);
        h_cc_ss_plateau_prt_gps = new TH2D("h_cc_ss_plateau_prt_gps","h_cc_ss_plateau_prt_gps",nParentGroups,0,nParentGroups,nParentGroups,0,nParentGroups);
        h_cc_ss_plateau_same_ancestors = new TH1D("h_cc_ss_plateau_same_ancestors","h_cc_ss_plateau_same_ancestors",2,0,2);
        h_cc_ss_plateau_sp = new TH1D("h_cc_ss_plateau_sp","h_cc_ss_plateau_sp",nAncestorGroups,0,nAncestorGroups);
        h_cc_ss_plateau_dp = new TH2D("h_cc_ss_plateau_dp","h_cc_ss_plateau_dp",nAncestorGroups,0,nAncestorGroups,nAncestorGroups,0,nAncestorGroups);


        for (int isign = 0; isign < ParamsSet::nSigns; isign++){
            h_unweighted_cc_both_from_c_same_ancestors[isign][0] = new TH1D(Form("h_unweighted_cc_both_from_c_same_ancestors_sign%d_near",isign+1),Form("h_unweighted_cc_both_from_c_same_ancestors_sign%d_near",isign+1),2,0,2);
            h_unweighted_cc_both_from_c_same_ancestors[isign][1] = new TH1D(Form("h_unweighted_cc_both_from_c_same_ancestors_sign%d_away",isign+1),Form("h_unweighted_cc_both_from_c_same_ancestors_sign%d_away",isign+1),2,0,2);
            h_cc_both_from_c_same_ancestors[isign][0] = new TH1D(Form("h_cc_both_from_c_same_ancestors_sign%d_near",isign+1),Form("h_cc_both_from_c_same_ancestors_sign%d_near",isign+1),2,0,2);
            h_cc_both_from_c_same_ancestors[isign][1] = new TH1D(Form("h_cc_both_from_c_same_ancestors_sign%d_away",isign+1),Form("h_cc_both_from_c_same_ancestors_sign%d_away",isign+1),2,0,2);
            h_unweighted_cc_both_from_c_ancestor_sp[isign][0] = new TH1D(Form("h_unweighted_cc_both_from_c_ancestor_sp_sign%d_near",isign+1),Form("h_unweighted_cc_both_from_c_ancestor_sp_sign%d_near",isign+1),nAncestorGroups,0,nAncestorGroups);
            h_unweighted_cc_both_from_c_ancestor_sp[isign][1] = new TH1D(Form("h_unweighted_cc_both_from_c_ancestor_sp_sign%d_away",isign+1),Form("h_unweighted_cc_both_from_c_ancestor_sp_sign%d_away",isign+1),nAncestorGroups,0,nAncestorGroups);
            h_unweighted_cc_both_from_c_ancestor_dp[isign][0] = new TH2D(Form("h_unweighted_cc_both_from_c_ancestor_dp_sign%d_near",isign+1),Form("h_unweighted_cc_both_from_c_ancestor_dp_sign%d_near",isign+1),nAncestorGroups,0,nAncestorGroups,nAncestorGroups,0,nAncestorGroups);
            h_unweighted_cc_both_from_c_ancestor_dp[isign][1] = new TH2D(Form("h_unweighted_cc_both_from_c_ancestor_dp_sign%d_away",isign+1),Form("h_unweighted_cc_both_from_c_ancestor_dp_sign%d_away",isign+1),nAncestorGroups,0,nAncestorGroups,nAncestorGroups,0,nAncestorGroups);
            h_cc_both_from_c_ancestor_sp[isign][0] = new TH1D(Form("h_cc_both_from_c_ancestor_sp_sign%d_near",isign+1),Form("h_cc_both_from_c_ancestor_sp_sign%d_near",isign+1),nAncestorGroups,0,nAncestorGroups);
            h_cc_both_from_c_ancestor_sp[isign][1] = new TH1D(Form("h_cc_both_from_c_ancestor_sp_sign%d_away",isign+1),Form("h_cc_both_from_c_ancestor_sp_sign%d_away",isign+1),nAncestorGroups,0,nAncestorGroups);
            h_cc_both_from_c_ancestor_dp[isign][0] = new TH2D(Form("h_cc_both_from_c_ancestor_dp_sign%d_near",isign+1),Form("h_cc_both_from_c_ancestor_dp_sign%d_near",isign+1),nAncestorGroups,0,nAncestorGroups,nAncestorGroups,0,nAncestorGroups);
            h_cc_both_from_c_ancestor_dp[isign][1] = new TH2D(Form("h_cc_both_from_c_ancestor_dp_sign%d_away",isign+1),Form("h_cc_both_from_c_ancestor_dp_sign%d_away",isign+1),nAncestorGroups,0,nAncestorGroups,nAncestorGroups,0,nAncestorGroups);
        }
    }
}

#endif


