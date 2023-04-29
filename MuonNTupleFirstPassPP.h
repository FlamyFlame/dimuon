#ifndef __MuonNTupleFirstPassPP_h__
#define __MuonNTupleFirstPassPP_h__

#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <string.h>
#include <fstream>
// #include <cmath>  
#include "MuonPairData.h"
#include "ParamsSet.h"
#include "vector"
#include "TH1D.h"
#include "TH2D.h"

class MuonNTupleFirstPassPP{
    // Read through the N-tuple, apply appropriate cuts
    // Then fill in histograms and/or output trees
    // mode = 3: kill resonances & output single-muon information into a TTree
    // mode = 4: kill resonances & output muon-pair information into a TTree
    // NOW ONLY HAVE MODE = 3, 4


private:
// --------------------- general settings ---------------------------

    int mode = 4;
    ParamsSet pms;
    
// --------------------- input files & trees & data for setting branches---------------------------
  
    TChain          *fChain;   //!pointer to the analyzed TTree or TChain
  
    // Declaration of leaf types
    UInt_t          RunNumber;
    std::vector<float>   *muon_deltaP_overP           =nullptr;
  
    std::vector<int>     *muon_pair_muon1_index       =nullptr;
    std::vector<float>   *muon_pair_muon1_pt          =nullptr;
    std::vector<float>   *muon_pair_muon1_eta         =nullptr;
    std::vector<float>   *muon_pair_muon1_phi         =nullptr;
    std::vector<float>   *muon_pair_muon1_trk_pt          =nullptr;
    std::vector<float>   *muon_pair_muon1_trk_eta         =nullptr;
    std::vector<float>   *muon_pair_muon1_trk_phi         =nullptr;
    std::vector<int>     *muon_pair_muon1_quality     =nullptr;
    std::vector<float>   *muon_pair_muon1_d0          =nullptr;
    std::vector<float>   *muon_pair_muon1_z0          =nullptr;
  
    std::vector<int>     *muon_pair_muon2_index       =nullptr;
    std::vector<float>   *muon_pair_muon2_pt          =nullptr;
    std::vector<float>   *muon_pair_muon2_eta         =nullptr;
    std::vector<float>   *muon_pair_muon2_phi         =nullptr;
    std::vector<float>   *muon_pair_muon2_trk_pt          =nullptr;
    std::vector<float>   *muon_pair_muon2_trk_eta         =nullptr;
    std::vector<float>   *muon_pair_muon2_trk_phi         =nullptr;
    std::vector<int>     *muon_pair_muon2_quality     =nullptr;
    std::vector<float>   *muon_pair_muon2_d0          =nullptr;
    std::vector<float>   *muon_pair_muon2_z0          =nullptr;
  
    Bool_t               b_HLT_2mu4;
    std::vector<bool>    *dimuon_b_HLT_2mu4    =nullptr;
  
// --------------------- temporary muon and muonpair objects ---------------------------
  
    Muon* tempmuon = nullptr;
    MuonPairData* mpair = nullptr;
    std::vector<int> resonance_tagged_muon_index_list {};

  
// --------------------- class methods ---------------------------
  
    void InitInput();
    void InitOutput();
    void ProcessData();
    bool PassCuts();
    bool IsResonance();
    bool IsPhotoProduction();
    void FillSingleMuonTree();
    void FillMuonPairTree();
  
  
// --------------------- output file, histograms & trees ---------------------------
  
    TFile *m_outfile = nullptr;
    
    TTree* muonOutTree;
    // TTree* muonPairOutTree[ParamsSet::nSigns];
    TTree* muonPairOutTreeBinned[ParamsSet::ndRselcs][ParamsSet::nSigns];
    // TH1D* h_dphi_failing_photoprod;
    // TH2D* h_asym_acop_failing_photoprod;

public :
    bool isTight;
    MuonNTupleFirstPassPP();
    ~MuonNTupleFirstPassPP(){}
    void Run();
};


MuonNTupleFirstPassPP::MuonNTupleFirstPassPP(){
    isTight = false;
    if(mode != 3 && mode != 4){
        std::cout<<"Error:: Mode has to equal 3 or 4; code is used for outputting muon / muon-pair trees only."<<std::endl;
        throw std::exception();
    }
}

//initialize the TChain
void MuonNTupleFirstPassPP::InitInput(){

    fChain = new TChain("HeavyIonD3PD","HeavyIonD3PD");
    fChain->SetMakeClass(1);
    fChain->Add("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/New_All_DataPP2017_5TeV_Dec2021.root");
  
    fChain->SetBranchAddress("RunNumber"                   , &RunNumber);
    fChain->SetBranchAddress("muon_deltaP_overP"           , &muon_deltaP_overP);
  
    fChain->SetBranchAddress("muon_pair_muon1_index"       , &muon_pair_muon1_index);
    fChain->SetBranchAddress("muon_pair_muon1_pt"          , &muon_pair_muon1_pt);
    fChain->SetBranchAddress("muon_pair_muon1_eta"         , &muon_pair_muon1_eta);
    fChain->SetBranchAddress("muon_pair_muon1_phi"         , &muon_pair_muon1_phi);
    fChain->SetBranchAddress("muon_pair_muon1_quality"     , &muon_pair_muon1_quality);
    fChain->SetBranchAddress("muon_pair_muon1_d0"          , &muon_pair_muon1_d0);
    fChain->SetBranchAddress("muon_pair_muon1_z0"          , &muon_pair_muon1_z0);
    fChain->SetBranchAddress("muon_pair_muon1_trk_pt"      , &muon_pair_muon1_trk_pt);
    fChain->SetBranchAddress("muon_pair_muon1_trk_eta"     , &muon_pair_muon1_trk_eta);
    fChain->SetBranchAddress("muon_pair_muon1_trk_phi"     , &muon_pair_muon1_trk_phi);
  
    fChain->SetBranchAddress("muon_pair_muon2_index"       , &muon_pair_muon2_index);
    fChain->SetBranchAddress("muon_pair_muon2_pt"          , &muon_pair_muon2_pt);
    fChain->SetBranchAddress("muon_pair_muon2_eta"         , &muon_pair_muon2_eta);
    fChain->SetBranchAddress("muon_pair_muon2_phi"         , &muon_pair_muon2_phi);
    fChain->SetBranchAddress("muon_pair_muon2_quality"     , &muon_pair_muon2_quality);
    fChain->SetBranchAddress("muon_pair_muon2_d0"          , &muon_pair_muon2_d0);
    fChain->SetBranchAddress("muon_pair_muon2_z0"          , &muon_pair_muon2_z0);
    fChain->SetBranchAddress("muon_pair_muon2_trk_pt"      , &muon_pair_muon2_trk_pt);
    fChain->SetBranchAddress("muon_pair_muon2_trk_eta"     , &muon_pair_muon2_trk_eta);
    fChain->SetBranchAddress("muon_pair_muon2_trk_phi"     , &muon_pair_muon2_trk_phi);

    fChain->SetBranchAddress("b_HLT_2mu4"                  , &b_HLT_2mu4);
    fChain->SetBranchAddress("dimuon_b_HLT_2mu4"           , &dimuon_b_HLT_2mu4);


    //SetBranch Status
    fChain->SetBranchStatus("*"                               ,0);//switch off all branches, then enable just the ones that we need
    fChain->SetBranchStatus("RunNumber"                       ,1);
    fChain->SetBranchStatus("muon_deltaP_overP"               ,1);
  
    fChain->SetBranchStatus("muon_pair_muon1_index"           ,1);
    fChain->SetBranchStatus("muon_pair_muon1_pt"              ,1);
    fChain->SetBranchStatus("muon_pair_muon1_eta"             ,1);
    fChain->SetBranchStatus("muon_pair_muon1_phi"             ,1);
    fChain->SetBranchStatus("muon_pair_muon1_trk_pt"              ,1);
    fChain->SetBranchStatus("muon_pair_muon1_trk_eta"             ,1);
    fChain->SetBranchStatus("muon_pair_muon1_trk_phi"             ,1);
    fChain->SetBranchStatus("muon_pair_muon1_quality"         ,1);
    fChain->SetBranchStatus("muon_pair_muon1_d0"              ,1);
    fChain->SetBranchStatus("muon_pair_muon1_z0"              ,1);
  
    fChain->SetBranchStatus("muon_pair_muon2_index"           ,1);
    fChain->SetBranchStatus("muon_pair_muon2_pt"              ,1);
    fChain->SetBranchStatus("muon_pair_muon2_eta"             ,1);
    fChain->SetBranchStatus("muon_pair_muon2_phi"             ,1);
    fChain->SetBranchStatus("muon_pair_muon2_trk_pt"              ,1);
    fChain->SetBranchStatus("muon_pair_muon2_trk_eta"             ,1);
    fChain->SetBranchStatus("muon_pair_muon2_trk_phi"             ,1);
    fChain->SetBranchStatus("muon_pair_muon2_quality"         ,1);
    fChain->SetBranchStatus("muon_pair_muon2_d0"              ,1);
    fChain->SetBranchStatus("muon_pair_muon2_z0"              ,1);
  
    fChain->SetBranchStatus("b_HLT_2mu4"               ,1);
    fChain->SetBranchStatus("dimuon_b_HLT_2mu4"        ,1);
}


void MuonNTupleFirstPassPP::InitOutput(){

    std::string tight_postfix = (isTight)? "_tight" : "";
    if (mode == 3){
        m_outfile=new TFile(("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/single_muon_trees_pp" + tight_postfix + ".root").c_str(),"recreate");
        muonOutTree = new TTree("muon_tree","all single muons");
        muonOutTree->Branch("MuonObj",&tempmuon);
    }
    else{ //output muon pair trees
        m_outfile=new TFile(("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/muon_pairs_pp" + tight_postfix + ".root").c_str(),"recreate");

        for (unsigned int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
            // muonPairOutTree[ksign] = new TTree(Form("muon_pair_tree_sign%u",ksign+1),Form("all muon pairs, sign%u",ksign+1));
            // muonPairOutTree[ksign]->Branch("MuonPairObj",&mpair);
            for (unsigned int idr = 0; idr < ParamsSet::ndRselcs; idr++){
                muonPairOutTreeBinned[idr][ksign] = new TTree(Form("muon_pair_tree_dr%u_sign%u",idr+1,ksign+1),Form("all muon pairs, dr%u, sign%u",idr+1,ksign+1));
                muonPairOutTreeBinned[idr][ksign]->Branch("MuonPairObj",&mpair);
            }
        }
    }

    // h_dphi_failing_photoprod = new TH1D("h_dphi_failing_photoprod",";#Delta#phi;1/N_{evt} dN/d#Delta#phi", 128, -pms.PI, pms.PI);
    // h_asym_acop_failing_photoprod = new TH2D("h_asym_acop_failing_photoprod", ";#frac{#pi-|#Delta#phi|}{#pi};#frac{p_{T}^{lead} - p_{T}^{sublead}}{p_{T}^{lead} + p_{T}^{sublead}}",20,0,0.01,20,0,0.05);
}

#endif


