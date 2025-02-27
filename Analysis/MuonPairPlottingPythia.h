// Assumes that we exclude resonances

#ifndef __MuonPairPlottingPythia_h__
#define __MuonPairPlottingPythia_h__

#include <TROOT.h>
// #include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <string.h>
#include  <stdlib.h>
#include <fstream>
#include "ParamsSet.h"
// #include "MuonPairPythia.h"
#include "muon_pair_enums_MC.h"

#include "vector"
#include "TH1D.h"
#include "TH2D.h"

class MuonPairPlottingPythia{
  	//updates histograms - centrality binned
    //if need histograms with all centrality bins added together
    //can add at the histograms level


private:
// --------------------- general settings ---------------------------

    // int mode = 1;
    // bool isScram = true;
    // bool isScram = false;

  	ParamsSet pms;
    // static const int nFiles = 3;
    static const int nFiles = 2;
    static const int nKinRanges = 5;
    static const int nDphis = 2;
    // static const int nFlavors = ; // both from b, both from c, one b one c
    std::vector<float> kinRanges = {5., 10., 25., 60., 120., 3200.};


// --------------------- output file & histograms ---------------------------

  	TFile *outFile = nullptr;
    
    std::string dphi_regions[nDphis] = {"near", "away"};

    TH1D* h_pair_dP_overP[nDphis][ParamsSet::nSigns];
    TH1D* h_Dphi[nDphis][ParamsSet::nSigns];
    TH1D* h_DR[nDphis][ParamsSet::nSigns];
    TH1D* h_DR_zoomin[nDphis][ParamsSet::nSigns];
    TH1D* h_DR_jacobian_corrected[nDphis][ParamsSet::nSigns];
    TH1D* h_DR_zoomin_jacobian_corrected[nDphis][ParamsSet::nSigns];
    TH1D* h_pair_y[nDphis][ParamsSet::nSigns];
    TH1D* h_pt_asym[nDphis][ParamsSet::nSigns];
    TH1D* h_psrapidity_ordered_pt_asym[nDphis][ParamsSet::nSigns];
    TH1D* h_pair_pt_ptlead_ratio[nDphis][ParamsSet::nSigns];
    TH2D* h_eta_avg_Dphi[nDphis][ParamsSet::nSigns];
    TH2D* h_Deta_Dphi[nDphis][ParamsSet::nSigns];
    TH2D* h_eta1_eta2[nDphis][ParamsSet::nSigns];
    TH2D* h_eta_avg_Deta[nDphis][ParamsSet::nSigns];
    TH2D* h_pt1_pt2[nDphis][ParamsSet::nSigns];
    TH2D* h_ptlead_pair_pt[nDphis][ParamsSet::nSigns];
    TH2D* h_ptlead_pair_pt_zoomin[nDphis][ParamsSet::nSigns];
    TH2D* h_ptlead_pair_pt_log[nDphis][ParamsSet::nSigns];
    TH2D* h_minv_pair_pt[nDphis][ParamsSet::nSigns];
    TH2D* h_minv_pair_pt_zoomin[nDphis][ParamsSet::nSigns];
    TH2D* h_minv_pair_pt_jacobian_corrected[nDphis][ParamsSet::nSigns];
    TH2D* h_minv_pair_pt_zoomin_jacobian_corrected[nDphis][ParamsSet::nSigns];
    TH2D* h_minv_pair_pt_log[nDphis][ParamsSet::nSigns];

    TH2D* h_eta1_eta2_dphicut[ParamsSet::ndphiselcs][ParamsSet::nSigns];

    static const int nAncestorGroups = 5;
    std::string ancestor_grp_labels[nAncestorGroups + 2] = {"_GS_ISR_no_HS", "_gs_ISR_one_hard_scatt", "_diff_GS_same_HS", "_FC", "_gs_FSR", "_from_same_b", "_others"};
    int ancestor_grps[nAncestorGroups] = {same_gs_isr_zero_hard_scatt, same_gs_isr_one_hard_scatt, diff_gs_same_hard_scatt, fc, same_gs_phs_fsr};

    static const int nFlavors = 4;
    std::string flavor_grps[nFlavors] = {"_single_b", "_bb", "_cc", "_other_flavors"};

    TH1D* h_DR_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 2];
    TH1D* h_DR_zoomin_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 2];
    TH1D* h_DR_jacobian_corrected_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 2];
    TH1D* h_DR_zoomin_jacobian_corrected_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 2];
    TH1D* h_psrapidity_ordered_pt_asym_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 2];
    TH1D* h_pt_asym_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 2];
    TH1D* h_pair_pt_ptlead_ratio_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 2];
    
    TH1D* h_Qsplit_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 2];
    TH1D* h_Qsplit_pTHat_ratio_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 2];
    TH1D* h_Qsplit_mHat_ratio_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 2];
    
    TH2D* h_Deta_Dphi_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 2];
    TH2D* h_minv_pair_pt_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 2];
    TH2D* h_minv_pair_pt_zoomin_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 2];
    TH2D* h_minv_pair_pt_jacobian_corrected_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 2];
    TH2D* h_minv_pair_pt_zoomin_jacobian_corrected_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 2];
    TH2D* h_minv_pair_pt_log_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 2];
    TH2D* h_ptlead_pair_pt_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 2];
    TH2D* h_ptlead_pair_pt_zoomin_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 2];
    TH2D* h_ptlead_pair_pt_log_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 2];

    TH1D* h_DR_flavor_binned[ParamsSet::nSigns][nFlavors];
    TH1D* h_DR_zoomin_flavor_binned[ParamsSet::nSigns][nFlavors];
    TH1D* h_DR_jacobian_corrected_flavor_binned[ParamsSet::nSigns][nFlavors];
    TH1D* h_DR_zoomin_jacobian_corrected_flavor_binned[ParamsSet::nSigns][nFlavors];
    TH1D* h_psrapidity_ordered_pt_asym_flavor_binned[ParamsSet::nSigns][nFlavors];
    TH1D* h_pt_asym_flavor_binned[ParamsSet::nSigns][nFlavors];
    TH1D* h_pair_pt_ptlead_ratio_flavor_binned[ParamsSet::nSigns][nFlavors];
    
    TH2D* h_Deta_Dphi_flavor_binned[ParamsSet::nSigns][nFlavors];
    TH2D* h_minv_pair_pt_flavor_binned[ParamsSet::nSigns][nFlavors];
    TH2D* h_minv_pair_pt_zoomin_flavor_binned[ParamsSet::nSigns][nFlavors];
    TH2D* h_minv_pair_pt_jacobian_corrected_flavor_binned[ParamsSet::nSigns][nFlavors];
    TH2D* h_minv_pair_pt_zoomin_jacobian_corrected_flavor_binned[ParamsSet::nSigns][nFlavors];
    TH2D* h_minv_pair_pt_log_flavor_binned[ParamsSet::nSigns][nFlavors];
    TH2D* h_ptlead_pair_pt_flavor_binned[ParamsSet::nSigns][nFlavors];
    TH2D* h_ptlead_pair_pt_zoomin_flavor_binned[ParamsSet::nSigns][nFlavors];
    TH2D* h_ptlead_pair_pt_log_flavor_binned[ParamsSet::nSigns][nFlavors];


    TH1D* h_kinbin_pair_dP_overP[nKinRanges][nDphis][ParamsSet::nSigns];
    TH1D* h_kinbin_Dphi[nKinRanges][nDphis][ParamsSet::nSigns];
    TH1D* h_kinbin_DR[nKinRanges][nDphis][ParamsSet::nSigns];
    TH1D* h_kinbin_pair_y[nKinRanges][nDphis][ParamsSet::nSigns];
    TH1D* h_kinbin_pt_asym[nKinRanges][nDphis][ParamsSet::nSigns];
    TH1D* h_kinbin_pair_pt_ptlead_ratio[nKinRanges][nDphis][ParamsSet::nSigns];
    TH2D* h_kinbin_eta_avg_Dphi[nKinRanges][nDphis][ParamsSet::nSigns];
    TH2D* h_kinbin_Deta_Dphi[nKinRanges][nDphis][ParamsSet::nSigns];
    TH2D* h_kinbin_eta1_eta2[nKinRanges][nDphis][ParamsSet::nSigns];
    TH2D* h_kinbin_eta_avg_Deta[nKinRanges][nDphis][ParamsSet::nSigns];
    TH2D* h_kinbin_pt1_pt2[nKinRanges][nDphis][ParamsSet::nSigns];
    TH2D* h_kinbin_ptlead_pair_pt[nKinRanges][nDphis][ParamsSet::nSigns];
    TH2D* h_kinbin_minv_pair_pt[nKinRanges][nDphis][ParamsSet::nSigns];
    TH2D* h_kinbin_minv_pair_pt_log[nKinRanges][nDphis][ParamsSet::nSigns];


// --------------------- input files & trees & data for setting branches---------------------------

   	TFile* inFile[nFiles];
    // NEED TO CHECK IF NEW AFTER0322 CONTAINS 0429 + IF THE nevents_before_cuts ARE CORRECT
    // std::vector<std::string> file_names = {"/usatlas/u/yuhanguo/usatlasdata/pythia/muon_pairs_pythia_allto0318.root", "/usatlas/u/yuhanguo/usatlasdata/pythia/muon_pairs_pythia_after0322.root", "/usatlas/u/yuhanguo/usatlasdata/pythia/muon_pairs_pythia_0429.root"};
    std::vector<std::string> file_names = {"/usatlas/u/yuhanguo/usatlasdata/pythia/muon_pairs_pythia_allto0318_with_data_resonance_cuts.root", "/usatlas/u/yuhanguo/usatlasdata/pythia/muon_pairs_pythia_after0322_with_data_resonance_cuts.root"};
    
    // we need to multiply the previous weight for each of the five different kinematic ranges (k_i)
    // in each of the two files (file A, B) by the corresponding reweight factor
    // N^i_A / (N^i_A + N^i_B) or N^i_B / (N^i_A + N^i_B), respectively
    // std::vector<std::vector<int>> nevents_before_cuts = {{20000, 50000, 1000000, 1000000, 1000000}, {79990, 300000, 1000000, 1000000, 1000000}, {19920, 100000, 500000, 500000, 500000}};
    std::vector<std::vector<int>> nevents_before_cuts = {{20000, 50000, 1000000, 1000000, 1000000}, {99910, 400000, 1500000, 1500000, 1500000}};


    std::vector<int> nevents_before_cuts_combined {};

    // TTree *inTree_ctrbin[nKinRanges][ParamsSet::nSigns];
    // TTree *inTree[ParamsSet::nSigns];
    TTree *inTree[nFiles][nKinRanges][ParamsSet::nSigns];

    double weight[nFiles][nKinRanges][ParamsSet::nSigns];
    double pTHat[nFiles][nKinRanges][ParamsSet::nSigns];
    double Qsplit[nFiles][nKinRanges][ParamsSet::nSigns];
    double mHard_relevant[nFiles][nKinRanges][ParamsSet::nSigns];

    bool from_same_b[nFiles][nKinRanges][ParamsSet::nSigns];
    bool both_from_b[nFiles][nKinRanges][ParamsSet::nSigns];
    bool both_from_c[nFiles][nKinRanges][ParamsSet::nSigns];
    int muon_pair_origin_category[nFiles][nKinRanges][ParamsSet::nSigns];

    float pair_dPoverP[nFiles][nKinRanges][ParamsSet::nSigns];
    float pt_lead[nFiles][nKinRanges][ParamsSet::nSigns];
    float pair_pt[nFiles][nKinRanges][ParamsSet::nSigns];
    // float pair_eta[nFiles][nKinRanges][ParamsSet::nSigns];
    float pair_y[nFiles][nKinRanges][ParamsSet::nSigns];
    float asym[nFiles][nKinRanges][ParamsSet::nSigns];
    float dpt[nFiles][nKinRanges][ParamsSet::nSigns];
    float deta[nFiles][nKinRanges][ParamsSet::nSigns];
    float etaavg[nFiles][nKinRanges][ParamsSet::nSigns];
    float phiavg[nFiles][nKinRanges][ParamsSet::nSigns];
    float dphi[nFiles][nKinRanges][ParamsSet::nSigns];
    float dr[nFiles][nKinRanges][ParamsSet::nSigns];
    float minv[nFiles][nKinRanges][ParamsSet::nSigns];
    float m1pt[nFiles][nKinRanges][ParamsSet::nSigns];
    float m2pt[nFiles][nKinRanges][ParamsSet::nSigns];
    float m1eta[nFiles][nKinRanges][ParamsSet::nSigns];
    float m2eta[nFiles][nKinRanges][ParamsSet::nSigns];
    float m1phi[nFiles][nKinRanges][ParamsSet::nSigns];
    float m2phi[nFiles][nKinRanges][ParamsSet::nSigns];
    // int m1charge[nFiles][nKinRanges][ParamsSet::nSigns];
    // int m2charge[nFiles][nKinRanges][ParamsSet::nSigns];

// --------------------- class methods ---------------------------

   	void InitInput();
    void InitOutput();
   	void InitHists();
   	void ProcessData();
    void WriteOutput();
    bool PassSingleMuonGapCut(float meta, float mpt, int mcharge);
   	void FillHistograms(int nfile, int nkin, int nsign);
   	void FillPtBinnedHistograms(int nkin, int npt, int nsign);

public:
    int mode = 1;
    // std::string file_suffix = "_allto0318";
    // std::string file_suffix = "_after0322";
  	MuonPairPlottingPythia();
  	~MuonPairPlottingPythia(){}
  	void Run();

};

MuonPairPlottingPythia::MuonPairPlottingPythia(){
    // nFiles = file_names.size();

    if (nevents_before_cuts.size() != nFiles || file_names.size() != nFiles){
        std::cout << "The vectors nevents_before_cuts and file_names must have the same size: ";
        std::cout << "equals number of the batches/files to be added." << std::endl;
        throw std::exception();
    }else{
        for (auto nevents_cur_file : nevents_before_cuts){
            if (nevents_cur_file.size() != nKinRanges){
                std::cout << "All vectors in nevents_before_cuts must have the same size: ";
                std::cout << "equals number of kinematic ranges." << std::endl;
                throw std::exception();
            }
        }
    }

    for (int ikin = 0; ikin < nKinRanges; ikin++){
        int nevents_before_cuts_combined_kn = 0;
        for (int jfile = 0; jfile < nFiles; jfile++){
            nevents_before_cuts_combined_kn += nevents_before_cuts[jfile][ikin];
        }
        cout << nevents_before_cuts[0][ikin] << " " << nevents_before_cuts[1][ikin] << " " << nevents_before_cuts_combined_kn << endl;
        nevents_before_cuts_combined.push_back(nevents_before_cuts_combined_kn);
    }
}

void MuonPairPlottingPythia::InitInput(){

    for (int ifile = 0; ifile < nFiles; ifile++){
        inFile[ifile] = new TFile(file_names[ifile].c_str(),"read");
        for (int jkin = 0; jkin < nKinRanges; jkin++){
            for (int ksign = 0; ksign < ParamsSet::nSigns; ksign++){

                inTree[ifile][jkin][ksign] = (TTree*) inFile[ifile]->Get(Form("muon_pair_tree_kin%d_sign%d",jkin,ksign+1));  
                inTree[ifile][jkin][ksign]->SetBranchAddress("weight"          , &weight[ifile][jkin][ksign]);
                inTree[ifile][jkin][ksign]->SetBranchAddress("pTHat"          , &pTHat[ifile][jkin][ksign]);
                inTree[ifile][jkin][ksign]->SetBranchAddress("Qsplit"          , &Qsplit[ifile][jkin][ksign]);
                inTree[ifile][jkin][ksign]->SetBranchAddress("mHard_relevant"          , &mHard_relevant[ifile][jkin][ksign]);
                inTree[ifile][jkin][ksign]->SetBranchAddress("from_same_b"          , &from_same_b[ifile][jkin][ksign]);
                inTree[ifile][jkin][ksign]->SetBranchAddress("both_from_b"          , &both_from_b[ifile][jkin][ksign]);
                inTree[ifile][jkin][ksign]->SetBranchAddress("both_from_c"          , &both_from_c[ifile][jkin][ksign]);
                inTree[ifile][jkin][ksign]->SetBranchAddress("muon_pair_origin_category"          , &muon_pair_origin_category[ifile][jkin][ksign]);
            	inTree[ifile][jkin][ksign]->SetBranchAddress("pair_dPoverP"    , &pair_dPoverP[ifile][jkin][ksign]);
            	inTree[ifile][jkin][ksign]->SetBranchAddress("pt_lead"          , &pt_lead[ifile][jkin][ksign]);
            	inTree[ifile][jkin][ksign]->SetBranchAddress("pair_pt"          , &pair_pt[ifile][jkin][ksign]);
            	// inTree[ifile][jkin][ksign]->SetBranchAddress("pair_eta"     , &pair_eta[ifile][jkin][ksign]);
            	inTree[ifile][jkin][ksign]->SetBranchAddress("pair_y"           , &pair_y[ifile][jkin][ksign]);
                inTree[ifile][jkin][ksign]->SetBranchAddress("asym"             , &asym[ifile][jkin][ksign]);
            	inTree[ifile][jkin][ksign]->SetBranchAddress("dpt"           , &dpt[ifile][jkin][ksign]);
            	inTree[ifile][jkin][ksign]->SetBranchAddress("deta"       , &deta[ifile][jkin][ksign]);
            	inTree[ifile][jkin][ksign]->SetBranchAddress("etaavg"      , &etaavg[ifile][jkin][ksign]);
            	inTree[ifile][jkin][ksign]->SetBranchAddress("phiavg"            , &phiavg[ifile][jkin][ksign]);
            	inTree[ifile][jkin][ksign]->SetBranchAddress("dphi"     , &dphi[ifile][jkin][ksign]);
            	inTree[ifile][jkin][ksign]->SetBranchAddress("dr"        , &dr[ifile][jkin][ksign]);
            	inTree[ifile][jkin][ksign]->SetBranchAddress("minv"        , &minv[ifile][jkin][ksign]);
            	inTree[ifile][jkin][ksign]->SetBranchAddress("m1.pt"           , &m1pt[ifile][jkin][ksign]);
            	inTree[ifile][jkin][ksign]->SetBranchAddress("m2.pt"           , &m2pt[ifile][jkin][ksign]);
            	inTree[ifile][jkin][ksign]->SetBranchAddress("m1.eta"       , &m1eta[ifile][jkin][ksign]);
            	inTree[ifile][jkin][ksign]->SetBranchAddress("m2.eta"       , &m2eta[ifile][jkin][ksign]);
            	inTree[ifile][jkin][ksign]->SetBranchAddress("m1.phi"     	, &m1phi[ifile][jkin][ksign]);
            	inTree[ifile][jkin][ksign]->SetBranchAddress("m2.phi"     	, &m2phi[ifile][jkin][ksign]);
                // inTree[ifile][jkin][ksign]->SetBranchAddress("m1.charge"           , &m1charge[ifile][jkin][ksign]);
                // inTree[ifile][jkin][ksign]->SetBranchAddress("m2.charge"           , &m2charge[ifile][jkin][ksign]);

                inTree[ifile][jkin][ksign]->SetBranchStatus("*"                  ,0);//switch off all branches, then enable just the ones that we need
                inTree[ifile][jkin][ksign]->SetBranchStatus("weight"           ,1);
                inTree[ifile][jkin][ksign]->SetBranchStatus("pTHat"           ,1);
                inTree[ifile][jkin][ksign]->SetBranchStatus("Qsplit"           ,1);
                inTree[ifile][jkin][ksign]->SetBranchStatus("mHard_relevant"           ,1);
                inTree[ifile][jkin][ksign]->SetBranchStatus("from_same_b"           ,1);
                inTree[ifile][jkin][ksign]->SetBranchStatus("both_from_b"           ,1);
                inTree[ifile][jkin][ksign]->SetBranchStatus("both_from_c"           ,1);
                inTree[ifile][jkin][ksign]->SetBranchStatus("muon_pair_origin_category"           ,1);
                inTree[ifile][jkin][ksign]->SetBranchStatus("pair_dPoverP"           ,1);
                inTree[ifile][jkin][ksign]->SetBranchStatus("pt_lead"           ,1);
                inTree[ifile][jkin][ksign]->SetBranchStatus("pair_pt"           ,1);
                inTree[ifile][jkin][ksign]->SetBranchStatus("pair_y"           ,1);
                inTree[ifile][jkin][ksign]->SetBranchStatus("asym"           ,1);
                inTree[ifile][jkin][ksign]->SetBranchStatus("dpt"           ,1);
                inTree[ifile][jkin][ksign]->SetBranchStatus("deta"           ,1);
                inTree[ifile][jkin][ksign]->SetBranchStatus("etaavg"           ,1);
                inTree[ifile][jkin][ksign]->SetBranchStatus("phiavg"           ,1);
                inTree[ifile][jkin][ksign]->SetBranchStatus("dphi"           ,1);
                inTree[ifile][jkin][ksign]->SetBranchStatus("dr"           ,1);
                inTree[ifile][jkin][ksign]->SetBranchStatus("minv"           ,1);
                inTree[ifile][jkin][ksign]->SetBranchStatus("m1.pt"           ,1);
                inTree[ifile][jkin][ksign]->SetBranchStatus("m2.pt"           ,1);
                inTree[ifile][jkin][ksign]->SetBranchStatus("m1.eta"           ,1);
                inTree[ifile][jkin][ksign]->SetBranchStatus("m2.eta"           ,1);
                inTree[ifile][jkin][ksign]->SetBranchStatus("m1.phi"           ,1);
                inTree[ifile][jkin][ksign]->SetBranchStatus("m2.phi"           ,1);
                // inTree[ifile][jkin][ksign]->SetBranchStatus("m1.charge"           ,1);
                // inTree[ifile][jkin][ksign]->SetBranchStatus("m2.charge"           ,1);
            }
        }
    }
    
}

void MuonPairPlottingPythia::InitOutput(){
    outFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/pythia/histograms_pythia_combined.root","recreate");
}

void MuonPairPlottingPythia::InitHists(){
    // int nDR_bins = (isMCTruthBB || isMCTruthCC)? 80 : 200;
    // int nDphi_bins = (isMCTruthBB || isMCTruthCC)? 64 : 128;
    // int neta_bins = (isMCTruthBB || isMCTruthCC)? 50 : 100;
    // int nDeta_bins = (isMCTruthBB || isMCTruthCC)? 100 : 200;
    // int npair_y_bins = (isMCTruthBB || isMCTruthCC)? 45 : 90;
    // int npt_asym_bins = (isMCTruthBB || isMCTruthCC)? 50 : 100;
    // int npair_pt_ptlead_ratio_bins = (isMCTruthBB || isMCTruthCC)? 50 : 100;
    int nDR_bins =                          100;
    int nDR_zoomin_bins =                   100;
    int nDphi_bins =                        64;
    int neta_bins =                         50;
    int nDeta_bins =                        100;
    int npair_y_bins =                      45;
    int npt_asym_bins =                     50;
    int npair_pt_ptlead_ratio_bins =        50;
    int nminv_bins_linear = 50;
    int npair_pT_bins_linear = 50;
    int npT_lead_bins_linear = 50;

    
    static const int nminv_bins_log = 40;
    float minv_logpow[ParamsSet::nSigns];
    minv_logpow[0] = 0.062;
    minv_logpow[1] = 0.045;
    float minv_max = 60;
    double minv_bins_log[ParamsSet::nSigns][nminv_bins_log+1];

    for (unsigned int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
        for(int iminv = 0; iminv <= nminv_bins_log; iminv++){
            minv_bins_log[ksign][iminv] = minv_max * pow(10.0, (static_cast<float>(iminv - nminv_bins_log)) * minv_logpow[ksign]);
        }
    }

    // h_pt_asym[jflavor][ksign] = new TH1D(Form("h_pt_asym_%s_sign%d",flavor_labels[jflavor].c_str(),ksign+1),";A = (pT1 - pT2)/(pT1 + pT2);d#sigma/dA", npt_asym_bins,0,1.);
    // h_psrapidity_ordered_pt_asym[jflavor][ksign] = new TH1D(Form("h_psrapidity_ordered_pt_asym_%s_sign%d",flavor_labels[jflavor].c_str(),ksign+1),";A' = (pT_{large |#eta|} - pT_{small |#eta|})/(pT_{large |#eta|} + pT_{small |#eta|});d#sigma/dA'", npt_asym_bins*2,-1.,1.);
    // h_pair_pt_ptlead_ratio[jflavor][ksign] = new TH1D(Form("h_pair_pt_ptlead_ratio_%s_sign%d",flavor_labels[jflavor].c_str(),ksign+1),";#frac{p_{T}^{pair}}{p_{T}^{lead}};d#sigma/d#frac{p_{T}^{pair}}{p_{T}^{lead}}", npair_pt_ptlead_ratio_bins,0,2.);
   	// h_Deta_Dphi[jflavor][ksign] = new TH2D(Form("h_Deta_Dphi_%s_sign%d",flavor_labels[jflavor].c_str(),ksign+1),";#Delta#phi;#Delta#eta", nDphi_bins,-pms.PI,pms.PI,nDeta_bins,-4.8,4.8);
    // h_ptlead_pair_pt[jflavor][ksign] = new TH2D(Form("h_ptlead_pair_pt_%s_sign%d",flavor_labels[jflavor].c_str(),ksign+1),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",pms.npairPT_bins,pms.pairPTBins[ksign][2],pms.npt_bins,pms.pTBins);
    // h_minv_pair_pt[jflavor][ksign] = new TH2D(Form("h_minv_pair_pt_%s_sign%d",flavor_labels[jflavor].c_str(),ksign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",pms.npairPT_bins,pms.pairPTBins[ksign][2],nminv_bins,minv_bins[ksign]);
      
    if (mode == 1){
        for (unsigned int ksign = 0; ksign < ParamsSet::nSigns; ksign++){

            for (int kflav = 0; kflav < nFlavors; kflav++){
                h_DR_flavor_binned[ksign][kflav] = new TH1D(Form("h_DR_sign%d%s",ksign+1,flavor_grps[kflav].c_str()),";#Delta R;1/N_{evt} dN/d#Delta R", nDR_bins,0,pms.deltaR_thrsh[2]);
                h_DR_zoomin_flavor_binned[ksign][kflav] = new TH1D(Form("h_DR_sign%d%s",ksign+1,flavor_grps[kflav].c_str()),";#Delta R;1/N_{evt} dN/d#Delta R", nDR_zoomin_bins,0,pms.deltaR_thrsh_zoomin);
                h_DR_jacobian_corrected_flavor_binned[ksign][kflav] = new TH1D(Form("h_DR_sign%d%s",ksign+1,flavor_grps[kflav].c_str()),";#Delta R;#frac{1}{N_{evt}} #frac{1}{#Delta R} #frac{dN}{d#Delta R}", nDR_bins,0,pms.deltaR_thrsh[2]);
                h_DR_zoomin_jacobian_corrected_flavor_binned[ksign][kflav] = new TH1D(Form("h_DR_sign%d%s",ksign+1,flavor_grps[kflav].c_str()),";#Delta R;#frac{1}{N_{evt}} #frac{1}{#Delta R} #frac{dN}{d#Delta R}", nDR_zoomin_bins,0,pms.deltaR_thrsh_zoomin);
                h_pt_asym_flavor_binned[ksign][kflav] = new TH1D(Form("h_pt_asym_sign%d%s",ksign+1,flavor_grps[kflav].c_str()),";A = (pT1 - pT2)/(pT1 + pT2);d#sigma/dA", npt_asym_bins,0,1.);
                h_psrapidity_ordered_pt_asym_flavor_binned[ksign][kflav] = new TH1D(Form("h_psrapidity_ordered_pt_asym_sign%d%s",ksign+1,flavor_grps[kflav].c_str()),";A' = (pT_{large |#eta|} - pT_{small |#eta|})/(pT_{large |#eta|} + pT_{small |#eta|});d#sigma/dA'", npt_asym_bins*2,-1.,1.);
                h_pair_pt_ptlead_ratio_flavor_binned[ksign][kflav] = new TH1D(Form("h_pair_pt_ptlead_ratio_sign%d%s",ksign+1,flavor_grps[kflav].c_str()),";#frac{p_{T}^{pair}}{p_{T}^{lead}};d#sigma/d#frac{p_{T}^{pair}}{p_{T}^{lead}}", npair_pt_ptlead_ratio_bins,0,2.);
                
                h_Deta_Dphi_flavor_binned[ksign][kflav] = new TH2D(Form("h_Deta_Dphi_sign%d%s",ksign+1,flavor_grps[kflav].c_str()),";#Delta#phi;#Delta#eta", nDphi_bins,-pms.PI,pms.PI,nDeta_bins,-4.8,4.8);
                h_minv_pair_pt_flavor_binned[ksign][kflav] = new TH2D(Form("h_minv_pair_pt_sign%d%s",ksign+1,flavor_grps[kflav].c_str()),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",npair_pT_bins_linear,0,30,nminv_bins_linear,0,30);
                h_minv_pair_pt_zoomin_flavor_binned[ksign][kflav] = new TH2D(Form("h_minv_pair_pt_zoomin_sign%d%s",ksign+1,flavor_grps[kflav].c_str()),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",npair_pT_bins_linear,0,20,nminv_bins_linear,0,4);
                h_minv_pair_pt_jacobian_corrected_flavor_binned[ksign][kflav] = new TH2D(Form("h_minv_pair_pt_jacobian_corrected_sign%d%s",ksign+1,flavor_grps[kflav].c_str()),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",npair_pT_bins_linear,0,30,nminv_bins_linear,0,30);
                h_minv_pair_pt_zoomin_jacobian_corrected_flavor_binned[ksign][kflav] = new TH2D(Form("h_minv_pair_pt_zoomin_jacobian_corrected_sign%d%s",ksign+1,flavor_grps[kflav].c_str()),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",npair_pT_bins_linear,0,20,nminv_bins_linear,0,4);
                h_minv_pair_pt_log_flavor_binned[ksign][kflav] = new TH2D(Form("h_minv_pair_pt_log_sign%d%s",ksign+1,flavor_grps[kflav].c_str()),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",pms.npairPT_bins,pms.pairPTBins[0][2],nminv_bins_log,minv_bins_log[0]);
                h_ptlead_pair_pt_flavor_binned[ksign][kflav] = new TH2D(Form("h_ptlead_pair_pt_sign%d%s",ksign+1,flavor_grps[kflav].c_str()),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",npair_pT_bins_linear,0,30,npT_lead_bins_linear,0,30);
                h_ptlead_pair_pt_zoomin_flavor_binned[ksign][kflav] = new TH2D(Form("h_ptlead_pair_pt_zoomin_sign%d%s",ksign+1,flavor_grps[kflav].c_str()),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",npair_pT_bins_linear,0,20,npT_lead_bins_linear,0,15);
                h_ptlead_pair_pt_log_flavor_binned[ksign][kflav] = new TH2D(Form("h_ptlead_pair_pt_log_sign%d%s",ksign+1,flavor_grps[kflav].c_str()),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",pms.npairPT_bins,pms.pairPTBins[0][2],pms.npt_bins,pms.pTBins);

                // h_Qsplit_flavor_binned[ksign][kflav] = new TH1D(Form("h_Qsplit_sign%d%s",ksign+1,flavor_grps[kflav].c_str()),";Q_{split}", 160, 0, 80.);
                // h_Qsplit_pTHat_ratio_flavor_binned[ksign][kflav] = new TH1D(Form("h_Qsplit_pTHat_ratio_sign%d%s",ksign+1,flavor_grps[kflav].c_str()),";#frac{Q_{split}}{#hat{p}_{T}}};d#sigma/d#frac{Q_{split}}{#hat{p}_{T}}}", 100,-2.,2.);
                // h_Qsplit_mHat_ratio_flavor_binned[ksign][kflav] = new TH1D(Form("h_Qsplit_mHat_ratio_sign%d%s",ksign+1,flavor_grps[kflav].c_str()),";#frac{Q_{split}}{#sqrt{#hat{s}}};d#sigma/d#frac{Q_{split}}{#sqrt{#hat{s}}}", 140,-0.7,0.7);
                
                h_DR_flavor_binned[ksign][kflav]->Sumw2();
                h_DR_zoomin_flavor_binned[ksign][kflav]->Sumw2();
                h_DR_jacobian_corrected_flavor_binned[ksign][kflav]->Sumw2();
                h_DR_zoomin_jacobian_corrected_flavor_binned[ksign][kflav]->Sumw2();
                h_pt_asym_flavor_binned[ksign][kflav]->Sumw2();
                h_psrapidity_ordered_pt_asym_flavor_binned[ksign][kflav]->Sumw2();
                h_pair_pt_ptlead_ratio_flavor_binned[ksign][kflav]->Sumw2();
                h_Deta_Dphi_flavor_binned[ksign][kflav]->Sumw2();
                h_minv_pair_pt_flavor_binned[ksign][kflav]->Sumw2();
                h_minv_pair_pt_zoomin_flavor_binned[ksign][kflav]->Sumw2();
                h_minv_pair_pt_jacobian_corrected_flavor_binned[ksign][kflav]->Sumw2();
                h_minv_pair_pt_zoomin_jacobian_corrected_flavor_binned[ksign][kflav]->Sumw2();
                h_minv_pair_pt_log_flavor_binned[ksign][kflav]->Sumw2();
                h_ptlead_pair_pt_flavor_binned[ksign][kflav]->Sumw2();
                h_ptlead_pair_pt_zoomin_flavor_binned[ksign][kflav]->Sumw2();
                h_ptlead_pair_pt_log_flavor_binned[ksign][kflav]->Sumw2();
                // h_Qsplit_flavor_binned[ksign][kflav]->Sumw2();
                // h_Qsplit_pTHat_ratio_flavor_binned[ksign][kflav]->Sumw2();
                // h_Qsplit_mHat_ratio_flavor_binned[ksign][kflav]->Sumw2();
            }


            for (int kgrp = 0; kgrp < nAncestorGroups + 2; kgrp++){
                h_DR_ancestor_binned[ksign][kgrp] = new TH1D(Form("h_DR_sign%d%s",ksign+1,ancestor_grp_labels[kgrp].c_str()),";#Delta R;1/N_{evt} dN/d#Delta R", nDR_bins,0,pms.deltaR_thrsh[2]);
                h_DR_zoomin_ancestor_binned[ksign][kgrp] = new TH1D(Form("h_DR_zoomin_sign%d%s",ksign+1,ancestor_grp_labels[kgrp].c_str()),";#Delta R;1/N_{evt} dN/d#Delta R", nDR_zoomin_bins,0,pms.deltaR_thrsh_zoomin);
                h_DR_jacobian_corrected_ancestor_binned[ksign][kgrp] = new TH1D(Form("h_DR_sign%d%s",ksign+1,ancestor_grp_labels[kgrp].c_str()),";#Delta R;#frac{1}{N_{evt}} #frac{1}{#Delta R} #frac{dN}{d#Delta R}", nDR_bins,0,pms.deltaR_thrsh[2]);
                h_DR_zoomin_jacobian_corrected_ancestor_binned[ksign][kgrp] = new TH1D(Form("h_DR_zoomin_sign%d%s",ksign+1,ancestor_grp_labels[kgrp].c_str()),";#Delta R;#frac{1}{N_{evt}} #frac{1}{#Delta R} #frac{dN}{d#Delta R}", nDR_zoomin_bins,0,pms.deltaR_thrsh_zoomin);
                h_pt_asym_ancestor_binned[ksign][kgrp] = new TH1D(Form("h_pt_asym_sign%d%s",ksign+1,ancestor_grp_labels[kgrp].c_str()),";A = (pT1 - pT2)/(pT1 + pT2);d#sigma/dA", npt_asym_bins,0,1.);
                h_psrapidity_ordered_pt_asym_ancestor_binned[ksign][kgrp] = new TH1D(Form("h_psrapidity_ordered_pt_asym_sign%d%s",ksign+1,ancestor_grp_labels[kgrp].c_str()),";A' = (pT_{large |#eta|} - pT_{small |#eta|})/(pT_{large |#eta|} + pT_{small |#eta|});d#sigma/dA'", npt_asym_bins*2,-1.,1.);
                h_pair_pt_ptlead_ratio_ancestor_binned[ksign][kgrp] = new TH1D(Form("h_pair_pt_ptlead_ratio_sign%d%s",ksign+1,ancestor_grp_labels[kgrp].c_str()),";#frac{p_{T}^{pair}}{p_{T}^{lead}};d#sigma/d#frac{p_{T}^{pair}}{p_{T}^{lead}}", npair_pt_ptlead_ratio_bins,0,2.);
                
                h_Deta_Dphi_ancestor_binned[ksign][kgrp] = new TH2D(Form("h_Deta_Dphi_sign%d%s",ksign+1,ancestor_grp_labels[kgrp].c_str()),";#Delta#phi;#Delta#eta", nDphi_bins,-pms.PI,pms.PI,nDeta_bins,-4.8,4.8);
                h_minv_pair_pt_ancestor_binned[ksign][kgrp] = new TH2D(Form("h_minv_pair_pt_sign%d%s",ksign+1,ancestor_grp_labels[kgrp].c_str()),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",npair_pT_bins_linear,0,30,nminv_bins_linear,0,30);
                h_minv_pair_pt_zoomin_ancestor_binned[ksign][kgrp] = new TH2D(Form("h_minv_pair_pt_zoomin_sign%d%s",ksign+1,ancestor_grp_labels[kgrp].c_str()),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",npair_pT_bins_linear,0,20,nminv_bins_linear,0,4);
                h_minv_pair_pt_jacobian_corrected_ancestor_binned[ksign][kgrp] = new TH2D(Form("h_minv_pair_pt_jacobian_corrected_sign%d%s",ksign+1,ancestor_grp_labels[kgrp].c_str()),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",npair_pT_bins_linear,0,30,nminv_bins_linear,0,30);
                h_minv_pair_pt_zoomin_jacobian_corrected_ancestor_binned[ksign][kgrp] = new TH2D(Form("h_minv_pair_pt_zoomin_jacobian_corrected_sign%d%s",ksign+1,ancestor_grp_labels[kgrp].c_str()),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",npair_pT_bins_linear,0,20,nminv_bins_linear,0,4);
                h_minv_pair_pt_log_ancestor_binned[ksign][kgrp] = new TH2D(Form("h_minv_pair_pt_log_sign%d%s",ksign+1,ancestor_grp_labels[kgrp].c_str()),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",pms.npairPT_bins,pms.pairPTBins[0][2],nminv_bins_log,minv_bins_log[0]);
                h_ptlead_pair_pt_ancestor_binned[ksign][kgrp] = new TH2D(Form("h_ptlead_pair_pt_sign%d%s",ksign+1,ancestor_grp_labels[kgrp].c_str()),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",npair_pT_bins_linear,0,30,npT_lead_bins_linear,0,30);
                h_ptlead_pair_pt_zoomin_ancestor_binned[ksign][kgrp] = new TH2D(Form("h_ptlead_pair_pt_zoomin_sign%d%s",ksign+1,ancestor_grp_labels[kgrp].c_str()),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",npair_pT_bins_linear,0,20,npT_lead_bins_linear,0,15);
                h_ptlead_pair_pt_log_ancestor_binned[ksign][kgrp] = new TH2D(Form("h_ptlead_pair_pt_log_sign%d%s",ksign+1,ancestor_grp_labels[kgrp].c_str()),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",pms.npairPT_bins,pms.pairPTBins[0][2],pms.npt_bins,pms.pTBins);

                h_Qsplit_ancestor_binned[ksign][kgrp] = new TH1D(Form("h_Qsplit_sign%d%s",ksign+1,ancestor_grp_labels[kgrp].c_str()),";Q_{split}", 160, 0, 80.);
                h_Qsplit_pTHat_ratio_ancestor_binned[ksign][kgrp] = new TH1D(Form("h_Qsplit_pTHat_ratio_sign%d%s",ksign+1,ancestor_grp_labels[kgrp].c_str()),";#frac{Q_{split}}{#hat{p}_{T}}};d#sigma/d#frac{Q_{split}}{#hat{p}_{T}}}", 100,-2.,2.);
                h_Qsplit_mHat_ratio_ancestor_binned[ksign][kgrp] = new TH1D(Form("h_Qsplit_mHat_ratio_sign%d%s",ksign+1,ancestor_grp_labels[kgrp].c_str()),";#frac{Q_{split}}{#sqrt{#hat{s}}};d#sigma/d#frac{Q_{split}}{#sqrt{#hat{s}}}", 140,-0.7,0.7);
                
                h_DR_ancestor_binned[ksign][kgrp]->Sumw2();
                h_DR_zoomin_ancestor_binned[ksign][kgrp]->Sumw2();
                h_DR_jacobian_corrected_ancestor_binned[ksign][kgrp]->Sumw2();
                h_DR_zoomin_jacobian_corrected_ancestor_binned[ksign][kgrp]->Sumw2();
                h_pt_asym_ancestor_binned[ksign][kgrp]->Sumw2();
                h_psrapidity_ordered_pt_asym_ancestor_binned[ksign][kgrp]->Sumw2();
                h_pair_pt_ptlead_ratio_ancestor_binned[ksign][kgrp]->Sumw2();
                h_Deta_Dphi_ancestor_binned[ksign][kgrp]->Sumw2();
                h_minv_pair_pt_ancestor_binned[ksign][kgrp]->Sumw2();
                h_minv_pair_pt_zoomin_ancestor_binned[ksign][kgrp]->Sumw2();
                h_minv_pair_pt_jacobian_corrected_ancestor_binned[ksign][kgrp]->Sumw2();
                h_minv_pair_pt_zoomin_jacobian_corrected_ancestor_binned[ksign][kgrp]->Sumw2();
                h_minv_pair_pt_log_ancestor_binned[ksign][kgrp]->Sumw2();
                h_ptlead_pair_pt_ancestor_binned[ksign][kgrp]->Sumw2();
                h_ptlead_pair_pt_zoomin_ancestor_binned[ksign][kgrp]->Sumw2();
                h_ptlead_pair_pt_log_ancestor_binned[ksign][kgrp]->Sumw2();
                h_Qsplit_ancestor_binned[ksign][kgrp]->Sumw2();
                h_Qsplit_pTHat_ratio_ancestor_binned[ksign][kgrp]->Sumw2();
                h_Qsplit_mHat_ratio_ancestor_binned[ksign][kgrp]->Sumw2();
            }

            
            for (unsigned int idphi= 0; idphi < ParamsSet::ndphiselcs; idphi++){
                h_eta1_eta2_dphicut[idphi][ksign] = new TH2D(Form("h_eta1_eta2_DPHI%d_sign%d",idphi+1,ksign+1),";#eta_{sublead};#eta_{lead}",50,-2.4,2.4, 50,-2.4,2.4);
                h_eta1_eta2_dphicut[idphi][ksign]->Sumw2();
            }



//"_%s_%s_%s_%s_%s_%s",pms.pt_labels[ipt], pms.ctr_labels[ictr], pms.signs[isign], pms.dphi_regions[idphi], pms.deta_regions[ideta], pms.gapcut_labels[igap], pms.photocut_labels[iphoto])
            
            for (unsigned int jdphi = 0; jdphi < nDphis; jdphi++){

                h_pair_dP_overP[jdphi][ksign] = new TH1D(Form("h_pair_dP_overP_%s_sign%d",dphi_regions[jdphi].c_str(),ksign+1),";(#Delta p / p)_{pair};1/N_{evt} dN/d(#Delta p / p)_{pair}" ,pms.deltaP_overP_nbins,0,pms.deltaP_overP_max);
                h_DR[jdphi][ksign] = new TH1D(Form("h_DR_%s_sign%d",dphi_regions[jdphi].c_str(),ksign+1),";#Delta R;1/N_{evt} dN/d#Delta R", nDR_bins,0,pms.deltaR_thrsh[2]);
                h_DR_zoomin[jdphi][ksign] = new TH1D(Form("h_DR_zoomin_%s_sign%d",dphi_regions[jdphi].c_str(),ksign+1),";#Delta R;1/N_{evt} dN/d#Delta R", nDR_zoomin_bins,0,pms.deltaR_thrsh_zoomin);
                h_DR_jacobian_corrected[jdphi][ksign] = new TH1D(Form("h_DR_jacobian_corrected_%s_sign%d",dphi_regions[jdphi].c_str(),ksign+1),";#Delta R;#frac{1}{N_{evt}} #frac{1}{#Delta R} #frac{dN}{d#Delta R}", nDR_bins,0,pms.deltaR_thrsh[2]);
                h_DR_zoomin_jacobian_corrected[jdphi][ksign] = new TH1D(Form("h_DR_zoomin_jacobian_corrected_%s_sign%d",dphi_regions[jdphi].c_str(),ksign+1),";#Delta R;#frac{1}{N_{evt}} #frac{1}{#Delta R} #frac{dN}{d#Delta R}", nDR_zoomin_bins,0,pms.deltaR_thrsh_zoomin);
                h_Dphi[jdphi][ksign] = new TH1D(Form("h_Dphi_%s_sign%d",dphi_regions[jdphi].c_str(),ksign+1),";#Delta#phi;1/N_{evt} dN/d#Delta#phi", nDphi_bins,-pms.PI,pms.PI);
                h_pair_y[jdphi][ksign] = new TH1D(Form("h_pair_y_%s_sign%d",dphi_regions[jdphi].c_str(),ksign+1),";y_{pair};1/N_{evt} dN/dy_{pair}" ,npair_y_bins,-3,3);
                h_pt_asym[jdphi][ksign] = new TH1D(Form("h_pt_asym_%s_sign%d",dphi_regions[jdphi].c_str(),ksign+1),";A = (pT1 - pT2)/(pT1 + pT2);d#sigma/dA", npt_asym_bins,0,1.);
                h_psrapidity_ordered_pt_asym[jdphi][ksign] = new TH1D(Form("h_psrapidity_ordered_pt_asym_%s_sign%d",dphi_regions[jdphi].c_str(),ksign+1),";A' = (pT_{large |#eta|} - pT_{small |#eta|})/(pT_{large |#eta|} + pT_{small |#eta|});d#sigma/dA'", npt_asym_bins*2,-1.,1.);
                h_pair_pt_ptlead_ratio[jdphi][ksign] = new TH1D(Form("h_pair_pt_ptlead_ratio_%s_sign%d",dphi_regions[jdphi].c_str(),ksign+1),";#frac{p_{T}^{pair}}{p_{T}^{lead}};d#sigma/d#frac{p_{T}^{pair}}{p_{T}^{lead}}", npair_pt_ptlead_ratio_bins,0,2.);

                h_eta_avg_Dphi[jdphi][ksign] = new TH2D(Form("h_eta_avg_Dphi_%s_sign%d",dphi_regions[jdphi].c_str(),ksign+1),";#Delta#phi;#bar{#eta}", nDphi_bins,-pms.PI,pms.PI,neta_bins,-2.4,2.4);
                h_Deta_Dphi[jdphi][ksign] = new TH2D(Form("h_Deta_Dphi_%s_sign%d",dphi_regions[jdphi].c_str(),ksign+1),";#Delta#phi;#Delta#eta", nDphi_bins,-pms.PI,pms.PI,nDeta_bins,-4.8,4.8);
                h_eta1_eta2[jdphi][ksign] = new TH2D(Form("h_eta1_eta2_%s_sign%d",dphi_regions[jdphi].c_str(),ksign+1),";#eta_{sublead};#eta_{lead}",neta_bins,-2.4,2.4, neta_bins,-2.4,2.4);
                h_eta_avg_Deta[jdphi][ksign] = new TH2D(Form("h_eta_avg_Deta_%s_sign%d",dphi_regions[jdphi].c_str(),ksign+1),";#Delta#eta;#bar{#eta}",nDeta_bins,-4.8,4.8,neta_bins,-2.4,2.4);
                h_pt1_pt2[jdphi][ksign] = new TH2D(Form("h_pt1_pt2_%s_sign%d",dphi_regions[jdphi].c_str(),ksign+1),";p_{T}^{sublead} [GeV];p_{T}^{lead} [GeV]",pms.npt_bins,pms.pTBins,pms.npt_bins,pms.pTBins);
                // h_ptlead_pair_pt[jdphi][ksign] = new TH2D(Form("h_ptlead_pair_pt_%s_sign%d",dphi_regions[jdphi].c_str(),ksign+1),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",pms.npairPT_bins,pms.pairPTBins[ksign][2],pms.npt_bins,pms.pTBins);
                // h_minv_pair_pt[jdphi][ksign] = new TH2D(Form("h_minv_pair_pt_%s_sign%d",dphi_regions[jdphi].c_str(),ksign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",pms.npairPT_bins,pms.pairPTBins[ksign][2],nminv_bins,minv_bins[ksign]);
                h_ptlead_pair_pt[jdphi][ksign] = new TH2D(Form("h_ptlead_pair_pt_%s_sign%d",dphi_regions[jdphi].c_str(),ksign+1),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",npair_pT_bins_linear,0,30,npT_lead_bins_linear,0,30);
                h_ptlead_pair_pt_zoomin[jdphi][ksign] = new TH2D(Form("h_ptlead_pair_pt_zoomin_%s_sign%d",dphi_regions[jdphi].c_str(),ksign+1),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",npair_pT_bins_linear,0,20,npT_lead_bins_linear,0,15);
                h_ptlead_pair_pt_log[jdphi][ksign] = new TH2D(Form("h_ptlead_pair_pt_log_%s_sign%d",dphi_regions[jdphi].c_str(),ksign+1),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",pms.npairPT_bins,pms.pairPTBins[ksign][2],pms.npt_bins,pms.pTBins);
                h_minv_pair_pt[jdphi][ksign] = new TH2D(Form("h_minv_pair_pt_%s_sign%d",dphi_regions[jdphi].c_str(),ksign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",npair_pT_bins_linear,0,30,nminv_bins_linear,0,30);
                h_minv_pair_pt_zoomin[jdphi][ksign] = new TH2D(Form("h_minv_pair_pt_zoomin_%s_sign%d",dphi_regions[jdphi].c_str(),ksign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",npair_pT_bins_linear,0,20,nminv_bins_linear,0,4);
                h_minv_pair_pt_jacobian_corrected[jdphi][ksign] = new TH2D(Form("h_minv_pair_pt_jacobian_corrected_%s_sign%d",dphi_regions[jdphi].c_str(),ksign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",npair_pT_bins_linear,0,30,nminv_bins_linear,0,30);
                h_minv_pair_pt_zoomin_jacobian_corrected[jdphi][ksign] = new TH2D(Form("h_minv_pair_pt_zoomin_jacobian_corrected_%s_sign%d",dphi_regions[jdphi].c_str(),ksign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",npair_pT_bins_linear,0,20,nminv_bins_linear,0,4);
                h_minv_pair_pt_log[jdphi][ksign] = new TH2D(Form("h_minv_pair_pt_log_%s_sign%d",dphi_regions[jdphi].c_str(),ksign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",pms.npairPT_bins,pms.pairPTBins[ksign][2],nminv_bins_log,minv_bins_log[ksign]);


                h_pair_dP_overP[jdphi][ksign]->Sumw2();
                h_DR[jdphi][ksign]->Sumw2();
                h_DR_zoomin[jdphi][ksign]->Sumw2();
                h_DR_jacobian_corrected[jdphi][ksign]->Sumw2();
                h_DR_zoomin_jacobian_corrected[jdphi][ksign]->Sumw2();
                h_Dphi[jdphi][ksign]->Sumw2();
                h_pair_y[jdphi][ksign]->Sumw2();
                h_pt_asym[jdphi][ksign]->Sumw2();
                h_pair_pt_ptlead_ratio[jdphi][ksign]->Sumw2();
                h_eta_avg_Dphi[jdphi][ksign]->Sumw2();
                h_Deta_Dphi[jdphi][ksign]->Sumw2();
                h_eta1_eta2[jdphi][ksign]->Sumw2();
                h_eta_avg_Deta[jdphi][ksign]->Sumw2();
                h_pt1_pt2[jdphi][ksign]->Sumw2();
                h_ptlead_pair_pt[jdphi][ksign]->Sumw2();
                h_ptlead_pair_pt_zoomin[jdphi][ksign]->Sumw2();
                h_ptlead_pair_pt_log[jdphi][ksign]->Sumw2();
                h_minv_pair_pt[jdphi][ksign]->Sumw2();
                h_minv_pair_pt_zoomin[jdphi][ksign]->Sumw2();
                h_minv_pair_pt_jacobian_corrected[jdphi][ksign]->Sumw2();
                h_minv_pair_pt_zoomin_jacobian_corrected[jdphi][ksign]->Sumw2();
                h_minv_pair_pt_log[jdphi][ksign]->Sumw2();

                for (unsigned int ikin = 0; ikin < nKinRanges; ikin++){

                    h_kinbin_pair_dP_overP[ikin][jdphi][ksign] = new TH1D(Form("h_pair_dP_overP_k%d_%s_sign%d",ikin,dphi_regions[jdphi].c_str(),ksign+1),";(#Delta p / p)_{pair};1/N_{evt} dN/d(#Delta p / p)_{pair}" ,pms.deltaP_overP_nbins,0,pms.deltaP_overP_max);
                    h_kinbin_DR[ikin][jdphi][ksign] = new TH1D(Form("h_DR_k%d_%s_sign%d",ikin,dphi_regions[jdphi].c_str(),ksign+1),";#Delta R;1/N_{evt} dN/d#Delta R", nDR_bins,0,pms.deltaR_thrsh[2]);
                    h_kinbin_Dphi[ikin][jdphi][ksign] = new TH1D(Form("h_Dphi_k%d_%s_sign%d",ikin,dphi_regions[jdphi].c_str(),ksign+1),";#Delta#phi;1/N_{evt} dN/d#Delta#phi", nDphi_bins,-pms.PI,pms.PI);
                    h_kinbin_pair_y[ikin][jdphi][ksign] = new TH1D(Form("h_pair_y_k%d_%s_sign%d",ikin,dphi_regions[jdphi].c_str(),ksign+1),";y_{pair};1/N_{evt} dN/dy_{pair}" ,npair_y_bins,-3,3);
                    h_kinbin_pt_asym[ikin][jdphi][ksign] = new TH1D(Form("h_pt_asym_k%d_%s_sign%d",ikin,dphi_regions[jdphi].c_str(),ksign+1),";A = (pT1 - pT2)/(pT1 + pT2);d#sigma/dA", npt_asym_bins,0,1.);
                    h_kinbin_pair_pt_ptlead_ratio[ikin][jdphi][ksign] = new TH1D(Form("h_pair_pt_ptlead_ratio_k%d_%s_sign%d",ikin,dphi_regions[jdphi].c_str(),ksign+1),";#frac{p_{T}^{pair}}{p_{T}^{lead}};d#sigma/d#frac{p_{T}^{pair}}{p_{T}^{lead}}", npair_pt_ptlead_ratio_bins,0,2.);
                    
                    h_kinbin_eta_avg_Dphi[ikin][jdphi][ksign] = new TH2D(Form("h_eta_avg_Dphi_k%d_%s_sign%d",ikin,dphi_regions[jdphi].c_str(),ksign+1),";#Delta#phi;#bar{#eta}", nDphi_bins,-pms.PI,pms.PI,neta_bins,-2.4,2.4);
                    h_kinbin_Deta_Dphi[ikin][jdphi][ksign] = new TH2D(Form("h_Deta_Dphi_k%d_%s_sign%d",ikin,dphi_regions[jdphi].c_str(),ksign+1),";#Delta#phi;#Delta#eta", nDphi_bins,-pms.PI,pms.PI,nDeta_bins,-4.8,4.8);
                    h_kinbin_eta1_eta2[ikin][jdphi][ksign] = new TH2D(Form("h_eta1_eta2_k%d_%s_sign%d",ikin,dphi_regions[jdphi].c_str(),ksign+1),";#eta_{sublead};#eta_{lead}",neta_bins,-2.4,2.4, neta_bins,-2.4,2.4);
                    h_kinbin_eta_avg_Deta[ikin][jdphi][ksign] = new TH2D(Form("h_eta_avg_Deta_k%d_%s_sign%d",ikin,dphi_regions[jdphi].c_str(),ksign+1),";#Delta#eta;#bar{#eta}",nDeta_bins,-4.8,4.8,neta_bins,-2.4,2.4);
                    h_kinbin_pt1_pt2[ikin][jdphi][ksign] = new TH2D(Form("h_pt1_pt2_k%d_%s_sign%d",ikin,dphi_regions[jdphi].c_str(),ksign+1),";p_{T}^{sublead} [GeV];p_{T}^{lead} [GeV]",pms.npt_bins,pms.pTBins,pms.npt_bins,pms.pTBins);
                    h_kinbin_ptlead_pair_pt[ikin][jdphi][ksign] = new TH2D(Form("h_ptlead_pair_pt_k%d_%s_sign%d",ikin,dphi_regions[jdphi].c_str(),ksign+1),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",pms.npairPT_bins,pms.pairPTBins[ksign][2],pms.npt_bins,pms.pTBins);
                    h_kinbin_minv_pair_pt[ikin][jdphi][ksign] = new TH2D(Form("h_minv_pair_pt_k%d_%s_sign%d",ikin,dphi_regions[jdphi].c_str(),ksign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",pms.npairPT_bins,pms.pairPTBins[ksign][2],nminv_bins_linear,0,30);
                    h_kinbin_minv_pair_pt_log[ikin][jdphi][ksign] = new TH2D(Form("h_minv_pair_pt_log_k%d_%s_sign%d",ikin,dphi_regions[jdphi].c_str(),ksign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",pms.npairPT_bins,pms.pairPTBins[ksign][2],nminv_bins_log,minv_bins_log[ksign]);
                    
                    h_kinbin_pair_dP_overP[ikin][jdphi][ksign]->Sumw2();
                    h_kinbin_DR[ikin][jdphi][ksign]->Sumw2();
                    h_kinbin_Dphi[ikin][jdphi][ksign]->Sumw2();
                    h_kinbin_pair_y[ikin][jdphi][ksign]->Sumw2();
                    h_kinbin_pt_asym[ikin][jdphi][ksign]->Sumw2();
                    h_kinbin_pair_pt_ptlead_ratio[ikin][jdphi][ksign]->Sumw2();
                    h_kinbin_eta_avg_Dphi[ikin][jdphi][ksign]->Sumw2();
                    h_kinbin_Deta_Dphi[ikin][jdphi][ksign]->Sumw2();
                    h_kinbin_eta1_eta2[ikin][jdphi][ksign]->Sumw2();
                    h_kinbin_eta_avg_Deta[ikin][jdphi][ksign]->Sumw2();
                    h_kinbin_pt1_pt2[ikin][jdphi][ksign]->Sumw2();
                    h_kinbin_ptlead_pair_pt[ikin][jdphi][ksign]->Sumw2();
                    h_kinbin_minv_pair_pt[ikin][jdphi][ksign]->Sumw2();
                    h_kinbin_minv_pair_pt_log[ikin][jdphi][ksign]->Sumw2();
                }
            }
        }
    }
}


#endif