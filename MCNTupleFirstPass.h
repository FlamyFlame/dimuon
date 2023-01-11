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

    const int nParentGroups = 5;
    std::vector<std::string> parentGroupLabels = {"direct b","b to c","direct c","s/light","photon"};

    static const int numCuts = 5;
    std::string cutLabels[numCuts] = {"no cut", "muon pT", "muon eta", "resonance", "photoproduction"};
    
// --------------------- input files & trees & data for setting branches---------------------------
  
    TChain          *fChain;   //!pointer to the analyzed TTree or TChain
  
    // Declaration of leaf types
    UInt_t          RunNumber;

    std::vector<int>* truth_id = nullptr;
    std::vector<int>* truth_barcode = nullptr;
    std::vector<int>* truth_qual = nullptr;
    std::vector<std::vector<int>>* truth_parents = nullptr;
    // std::vector<std::vector<int>>* truth_muon_parent_ids = nullptr;
    // std::vector<std::vector<int>>* truth_muon_parent_bars = nullptr;
    std::vector<float>   *EventWeights           =nullptr;

    std::vector<float>   *muon_pair_muon1_pt           =nullptr;
    std::vector<float>   *muon_pair_muon1_eta       =nullptr;
    std::vector<float>   *muon_pair_muon1_phi          =nullptr;
    // std::vector<int>     *muon_pair_muon1_ch         =nullptr;
    std::vector<float>     *muon_pair_muon1_ch         =nullptr;
    // std::vector<int>     *muon_pair_muon1_bar         =nullptr;
    std::vector<float>     *muon_pair_muon1_bar         =nullptr;

    std::vector<float>   *muon_pair_muon2_pt       =nullptr;
    std::vector<float>   *muon_pair_muon2_eta          =nullptr;
    std::vector<float>   *muon_pair_muon2_phi         =nullptr;
    // std::vector<int>     *muon_pair_muon2_ch         =nullptr;
    std::vector<float>     *muon_pair_muon2_ch         =nullptr;
    // std::vector<int>     *muon_pair_muon2_bar         =nullptr;
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

    // bool same_parents;
    bool m1_ancestor_is_incoming;
    bool m2_ancestor_is_incoming;

    // int cur_m1_beyond_earliest_parent_barcode;
    // int cur_m2_beyond_earliest_parent_barcode;
    // std::vector<int> cur_m1_beyond_earliest_parent_ids;
    // std::vector<int> cur_m2_beyond_earliest_parent_ids;
    
    std::vector<int> cur_m1_ancestor_ids;
    std::vector<int> cur_m2_ancestor_ids;
    std::vector<int> cur_m1_ancestor_bars;
    std::vector<int> cur_m2_ancestor_bars;
  
// --------------------- output file, histograms & trees ---------------------------
  
    TFile *m_outfile = nullptr;
    ofstream m_unspecified_parent_file;
    ofstream m_b_parent_file[2];
    // ofstream m_away_b_parent_file;
    ofstream m_c_parent_file[ParamsSet::nSigns][2];

    TTree* muonOutTree;
    TTree* muonPairOutTree[ParamsSet::nSigns];
    TTree* muonPairOutTreeBinned[ParamsSet::ndRselcs][ParamsSet::nSigns];

    TH1D* h_cutAcceptance[ParamsSet::nSigns];
    TH1D* h_numParents;
    TH1D* h_numMuonPairs;
    TH2D* h_MuonPairParentGroups[ParamsSet::nSigns][2];
    TH2D* h_MuonPairParentGroups_weighted[ParamsSet::nSigns][2];

    TH1D* h_bb_op_one_b_one_btoc[2];
    TH1D* h_weighted_bb_op_one_b_one_btoc[2];

    TH1D* h_cc_ss_both_from_c_same_prts[2];
    TH1D* h_cc_ss_both_from_c[2];
    // TH1D* h_weighted_cc_ss_both_from_c[2];

    TH1D* h_bb_both_from_b_same_prts[2];
    TH1D* h_bb_both_from_b[2];
    // TH1D* h_weighted_bb_both_from_b[2];

    TH2D* h_bb_both_from_b_ancestor[2];
    // TH2D* h_bb_both_from_b_ancestor_weighted[2];
    TH2D* h_cc_both_from_c_ancestor[ParamsSet::nSigns][2];
    // TH2D* h_cc_both_from_c_ancestor_weighted[ParamsSet::nSigns][2];

    std::vector<std::string> samePrtsLabels = {"Same Parents", "Different Parents"};
    std::vector<std::string> bb_op_one_b_one_btoc_labels = {"Same b", "Involve osc(s)", "Others"};
    std::vector<std::string> bb_both_from_b_ancestor_labels = {"Not incoming", "Incoming"};
    std::vector<std::string> cc_both_from_c_ancestor_labels = {"Not incoming", "Incoming"};

// --------------------- class methods ---------------------------
  
    void InitInput();
    void InitOutput();
    void ProcessData();
    bool PassCuts();
    bool IsResonance();
    bool IsPhotoProduction();
    int  ParentGrouping(int parent_id, bool c_tag);
    void UpdateCurParents(bool isFirstTime, bool isMuon1, std::vector<int>& cur_prt_bars, std::vector<int>& cur_prt_ids);
    int  CheckIfContainHeavyQuarks(std::vector<int>& cur_prt_ids, int quark_type);
    void FindSingleMuonParents(bool isMuon1);
    void FillMuonPairParents();
    void FillSingleMuonTree();
    void FillMuonPairTree();
  

public :
    int mode = 2;
    bool isTruth = true;
    // std::string mc_mode = "mc_truth_bb";
    std::string mc_mode = "mc_truth_cc";
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

    std::string sign_labels[ParamsSet::nSigns] = {"same_sign", "op_sign"};
    m_unspecified_parent_file.open(mcdir + "unspecified_parents.txt");
    if (mc_mode == "mc_truth_bb"){
        m_b_parent_file[0].open(mcdir + "b_parents_near.txt");
        m_b_parent_file[1].open(mcdir + "b_parents_away.txt");
    }else{
        for (int isign = 0; isign < ParamsSet::nSigns; isign++){
            m_c_parent_file[isign][0].open(Form("%sc_parents_%s_near.txt", mcdir.c_str(), sign_labels[isign].c_str()));
            m_c_parent_file[isign][1].open(Form("%sc_parents_%s_near.txt", mcdir.c_str(), sign_labels[isign].c_str()));
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

    h_numParents = new TH1D("h_numParents","h_numParents",3,0,3);
    h_numMuonPairs = new TH1D("h_numMuonPairs","h_numMuonPairs",6,0,6);
    for (int isign = 0; isign < ParamsSet::nSigns; isign++){
        h_cutAcceptance[isign] = new TH1D(Form("h_cutAcceptance_sign%d",isign+1),Form("h_cutAcceptance_sign%d",isign+1),numCuts,0,numCuts);
        h_MuonPairParentGroups[isign][0] = new TH2D(Form("h_MuonPairParentGroups_sign%d_near",isign+1),Form("h_MuonPairParentGroups_sign%d_near",isign+1),nParentGroups,0,nParentGroups,nParentGroups,0,nParentGroups);
        h_MuonPairParentGroups[isign][1] = new TH2D(Form("h_MuonPairParentGroups_sign%d_away",isign+1),Form("h_MuonPairParentGroups_sign%d_away",isign+1),nParentGroups,0,nParentGroups,nParentGroups,0,nParentGroups);
        h_MuonPairParentGroups_weighted[isign][0] = new TH2D(Form("h_MuonPairParentGroups_weighted_sign%d_near",isign+1),Form("h_MuonPairParentGroups_weighted_sign%d_near",isign+1),nParentGroups,0,nParentGroups,nParentGroups,0,nParentGroups);
        h_MuonPairParentGroups_weighted[isign][1] = new TH2D(Form("h_MuonPairParentGroups_weighted_sign%d_away",isign+1),Form("h_MuonPairParentGroups_weighted_sign%d_away",isign+1),nParentGroups,0,nParentGroups,nParentGroups,0,nParentGroups);
    }

    if (mc_mode == "mc_truth_bb"){
        h_bb_op_one_b_one_btoc[0] = new TH1D("h_bb_op_one_b_one_btoc_near","h_bb_op_one_b_one_btoc_near",3,0,3);
        h_bb_op_one_b_one_btoc[1] = new TH1D("h_bb_op_one_b_one_btoc_away","h_bb_op_one_b_one_btoc_away",3,0,3);
        h_weighted_bb_op_one_b_one_btoc[0] = new TH1D("h_weighted_bb_op_one_b_one_btoc_near","h_weighted_bb_op_one_b_one_btoc_near",3,0,3);
        h_weighted_bb_op_one_b_one_btoc[1] = new TH1D("h_weighted_bb_op_one_b_one_btoc_away","h_weighted_bb_op_one_b_one_btoc_away",3,0,3);
    
        h_bb_both_from_b_same_prts[0] = new TH1D("h_bb_both_from_b_same_prts_near","h_bb_both_from_b_same_prts_near",2,0,2);
        h_bb_both_from_b_same_prts[1] = new TH1D("h_bb_both_from_b_same_prts_away","h_bb_both_from_b_same_prts_away",2,0,2);
        // h_bb_both_from_b[0] = new TH1D("h_bb_both_from_b_near","h_bb_both_from_b_near",3,0,3);
        // h_bb_both_from_b[1] = new TH1D("h_bb_both_from_b_away","h_bb_both_from_b_away",3,0,3);
        // h_weighted_bb_both_from_b[0] = new TH1D("h_weighted_bb_both_from_b_near","h_weighted_bb_both_from_b_near",3,0,3);
        // h_weighted_bb_both_from_b[1] = new TH1D("h_weighted_bb_both_from_b_away","h_weighted_bb_both_from_b_away",3,0,3);
        h_bb_both_from_b_ancestor[0] = new TH2D("h_bb_both_from_b_ancestor_near","h_bb_both_from_b_ancestor_near",2,0,2,2,0,2);
        h_bb_both_from_b_ancestor[1] = new TH2D("h_bb_both_from_b_ancestor_away","h_bb_both_from_b_ancestor_away",2,0,2,2,0,2);
    }else if (mc_mode == "mc_truth_cc"){
        h_cc_ss_both_from_c_same_prts[0] = new TH1D("h_cc_ss_both_from_c_same_prts_near","h_cc_ss_both_from_c_same_prts_near",2,0,2);
        h_cc_ss_both_from_c_same_prts[1] = new TH1D("h_cc_ss_both_from_c_same_prts_away","h_cc_ss_both_from_c_same_prts_away",2,0,2);
        // h_cc_ss_both_from_c[0] = new TH1D("h_cc_ss_both_from_c_near","h_cc_ss_both_from_c_near",3,0,3);
        // h_cc_ss_both_from_c[1] = new TH1D("h_cc_ss_both_from_c_away","h_cc_ss_both_from_c_away",3,0,3);
        // h_weighted_cc_ss_both_from_c[0] = new TH1D("h_weighted_cc_ss_both_from_c_near","h_weighted_cc_ss_both_from_c_near",3,0,3);
        // h_weighted_cc_ss_both_from_c[1] = new TH1D("h_weighted_cc_ss_both_from_c_away","h_weighted_cc_ss_both_from_c_away",3,0,3);
        for (int isign = 0; isign < ParamsSet::nSigns; isign++){
            h_cc_both_from_c_ancestor[isign][0] = new TH2D(Form("h_cc_both_from_c_ancestor_sign%d_near",isign+1),Form("h_cc_both_from_c_ancestor_sign%d_near",isign+1),2,0,2,2,0,2);
            h_cc_both_from_c_ancestor[isign][1] = new TH2D(Form("h_cc_both_from_c_ancestor_sign%d_away",isign+1),Form("h_cc_both_from_c_ancestor_sign%d_away",isign+1),2,0,2,2,0,2);
        }
    }
}

#endif


