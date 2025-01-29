// Assumes that we exclude resonances

#ifndef __MuonPairPlottingPP_h__
#define __MuonPairPlottingPP_h__

#include <TROOT.h>
// #include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <string.h>
#include  <stdlib.h>
#include <fstream>
#include "ParamsSet.h"
#include "MuonPair.h"
#include "vector"
#include "TH1D.h"
#include "TH2D.h"

class MuonPairPlottingPP{
  	//updates histograms - centrality binned
    //if need histograms with all centrality bins added together
    //can add at the histograms level


private:
  	// --------------------- general settings ---------------------------

    // int mode = 1;
    // bool isScram = true;
    // bool isScram = false;

  	ParamsSet pms;

  	// --------------------- output file & histograms ---------------------------

  	TFile *outFile = nullptr;

    
    std::string dphi_regions[2] = {"near", "away"};

    // TH1D* h_pair_dP_overP[ParamsSet::nPtBins][ParamsSet::nSigns][ParamsSet::ndphiRegions][ParamsSet::ndetaRegions][ParamsSet::nPhotoProdCuts][ParamsSet::nGapCuts];
    // TH1D* h_Dphi[ParamsSet::nPtBins][ParamsSet::nSigns][ParamsSet::ndphiRegions][ParamsSet::ndetaRegions][ParamsSet::nPhotoProdCuts][ParamsSet::nGapCuts];
    // TH1D* h_DR[ParamsSet::nPtBins][ParamsSet::nSigns][ParamsSet::ndphiRegions][ParamsSet::ndetaRegions][ParamsSet::nPhotoProdCuts][ParamsSet::nGapCuts];
    // TH1D* h_pair_y[ParamsSet::nPtBins][ParamsSet::nSigns][ParamsSet::ndphiRegions][ParamsSet::ndetaRegions][ParamsSet::nPhotoProdCuts][ParamsSet::nGapCuts];
    // TH1D* h_pt_asym[ParamsSet::nPtBins][ParamsSet::nSigns][ParamsSet::ndphiRegions][ParamsSet::ndetaRegions][ParamsSet::nPhotoProdCuts][ParamsSet::nGapCuts];
    // TH1D* h_pair_pt_ptlead_ratio[ParamsSet::nPtBins][ParamsSet::nSigns][ParamsSet::ndphiRegions][ParamsSet::ndetaRegions][ParamsSet::nPhotoProdCuts][ParamsSet::nGapCuts];
    // TH2D* h_eta_avg_Dphi[ParamsSet::nPtBins][ParamsSet::nSigns][ParamsSet::ndphiRegions][ParamsSet::ndetaRegions][ParamsSet::nPhotoProdCuts][ParamsSet::nGapCuts];
    // TH2D* h_Deta_Dphi[ParamsSet::nPtBins][ParamsSet::nSigns][ParamsSet::ndphiRegions][ParamsSet::ndetaRegions][ParamsSet::nPhotoProdCuts][ParamsSet::nGapCuts];
    // TH2D* h_eta1_eta2[ParamsSet::nPtBins][ParamsSet::nSigns][ParamsSet::ndphiRegions][ParamsSet::ndetaRegions][ParamsSet::nPhotoProdCuts][ParamsSet::nGapCuts];
    // TH2D* h_eta_avg_Deta[ParamsSet::nPtBins][ParamsSet::nSigns][ParamsSet::ndphiRegions][ParamsSet::ndetaRegions][ParamsSet::nPhotoProdCuts][ParamsSet::nGapCuts];
    // TH2D* h_pt1_pt2[ParamsSet::nPtBins][ParamsSet::nSigns][ParamsSet::ndphiRegions][ParamsSet::ndetaRegions][ParamsSet::nPhotoProdCuts][ParamsSet::nGapCuts];
    // TH2D* h_ptlead_pair_pt[ParamsSet::nPtBins][ParamsSet::nSigns][ParamsSet::ndphiRegions][ParamsSet::ndetaRegions][ParamsSet::nPhotoProdCuts][ParamsSet::nGapCuts];
    // TH2D* h_ptlead_pair_pt_zoomin[ParamsSet::nPtBins][ParamsSet::nSigns][ParamsSet::ndphiRegions][ParamsSet::ndetaRegions][ParamsSet::nPhotoProdCuts][ParamsSet::nGapCuts];
    // TH2D* h_ptlead_pair_pt_log[ParamsSet::nPtBins][ParamsSet::nSigns][ParamsSet::ndphiRegions][ParamsSet::ndetaRegions][ParamsSet::nPhotoProdCuts][ParamsSet::nGapCuts];
    // TH2D* h_minv_pair_pt[ParamsSet::nPtBins][ParamsSet::nSigns][ParamsSet::ndphiRegions][ParamsSet::ndetaRegions][ParamsSet::nPhotoProdCuts][ParamsSet::nGapCuts];
    // TH2D* h_minv_pair_pt_zoomin[ParamsSet::nPtBins][ParamsSet::nSigns][ParamsSet::ndphiRegions][ParamsSet::ndetaRegions][ParamsSet::nPhotoProdCuts][ParamsSet::nGapCuts];
    // TH2D* h_minv_pair_pt_log[ParamsSet::nPtBins][ParamsSet::nSigns][ParamsSet::ndphiRegions][ParamsSet::ndetaRegions][ParamsSet::nPhotoProdCuts][ParamsSet::nGapCuts];

    // TH1D* h_pair_dP_overP[ParamsSet::ndphiRegions][ParamsSet::nSigns][ParamsSet::nGapCuts];
    // TH1D* h_Dphi[ParamsSet::ndphiRegions][ParamsSet::nSigns][ParamsSet::nGapCuts];
    // TH1D* h_DR[ParamsSet::ndphiRegions][ParamsSet::nSigns][ParamsSet::nGapCuts];
    // TH1D* h_pair_y[ParamsSet::ndphiRegions][ParamsSet::nSigns][ParamsSet::nGapCuts];
    // TH1D* h_pt_asym[ParamsSet::ndphiRegions][ParamsSet::nSigns][ParamsSet::nGapCuts];
    // TH1D* h_pair_pt_ptlead_ratio[ParamsSet::ndphiRegions][ParamsSet::nSigns][ParamsSet::nGapCuts];
    // TH2D* h_eta_avg_Dphi[ParamsSet::ndphiRegions][ParamsSet::nSigns][ParamsSet::nGapCuts];
    // TH2D* h_Deta_Dphi[ParamsSet::ndphiRegions][ParamsSet::nSigns][ParamsSet::nGapCuts];
    // TH2D* h_eta1_eta2[ParamsSet::ndphiRegions][ParamsSet::nSigns][ParamsSet::nGapCuts];
    // TH2D* h_eta_avg_Deta[ParamsSet::ndphiRegions][ParamsSet::nSigns][ParamsSet::nGapCuts];
    // TH2D* h_pt1_pt2[ParamsSet::ndphiRegions][ParamsSet::nSigns][ParamsSet::nGapCuts];
    // TH2D* h_ptlead_pair_pt[ParamsSet::ndphiRegions][ParamsSet::nSigns][ParamsSet::nGapCuts];
    // TH2D* h_ptlead_pair_pt_zoomin[ParamsSet::ndphiRegions][ParamsSet::nSigns][ParamsSet::nGapCuts];
    // TH2D* h_ptlead_pair_pt_log[ParamsSet::ndphiRegions][ParamsSet::nSigns][ParamsSet::nGapCuts];
    // TH2D* h_minv_pair_pt[ParamsSet::ndphiRegions][ParamsSet::nSigns][ParamsSet::nGapCuts];
    // TH2D* h_minv_pair_pt_zoomin[ParamsSet::ndphiRegions][ParamsSet::nSigns][ParamsSet::nGapCuts];
    // TH2D* h_minv_pair_pt_log[ParamsSet::ndphiRegions][ParamsSet::nSigns][ParamsSet::nGapCuts];

    TH1D* h_pair_dP_overP[2][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH1D* h_Dphi[2][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH1D* h_DR[2][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH1D* h_pair_y[2][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH1D* h_pt_asym[2][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH1D* h_pair_pt_ptlead_ratio[2][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_eta_avg_Dphi[2][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_Deta_Dphi[2][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_eta1_eta2[2][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_eta_avg_Deta[2][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_pt1_pt2[2][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_ptlead_pair_pt[2][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_ptlead_pair_pt_zoomin[2][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_ptlead_pair_pt_log[2][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_minv_pair_pt[2][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_minv_pair_pt_zoomin[2][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_minv_pair_pt_log[2][ParamsSet::nSigns][ParamsSet::nGapCuts];

    static const int nAncestorGroups = 4;

    // TH1D* h_DR_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 1];
    // TH1D* h_pt_asym_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 1];
    // // TH1D* h_psrapidity_ordered_pt_asym_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 1];
    // TH1D* h_pair_pt_ptlead_ratio_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 1];
    // TH2D* h_ptlead_pair_pt_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 1];
    // TH2D* h_Deta_Dphi_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 1];
    // TH2D* h_minv_pair_pt_ancestor_binned[ParamsSet::nSigns][nAncestorGroups + 1];

    // TH2D* h_unweighted_eta1_eta2_dphicut[ParamsSet::ndphiselcs][ParamsSet::nSigns][ParamsSet::nGapCuts];

    // --------------------- input files & trees & data for setting branches---------------------------

   	TFile *inFile;
    // TTree *inTree_ctrbin[ParamsSet::nSigns];
    // TTree *inTree[ParamsSet::nSigns];
    TTree *inTree[ParamsSet::nSigns];

    double weight[ParamsSet::nSigns];
  	float pair_dPoverP[ParamsSet::nSigns];
  	float pt_lead[ParamsSet::nSigns];
  	float pair_pt[ParamsSet::nSigns];
  	// float pair_eta[ParamsSet::nSigns];
  	float pair_y[ParamsSet::nSigns];
    float asym[ParamsSet::nSigns];
  	float dpt[ParamsSet::nSigns];
  	float deta[ParamsSet::nSigns];
  	float etaavg[ParamsSet::nSigns];
  	float phiavg[ParamsSet::nSigns];
  	float dphi[ParamsSet::nSigns];
  	float dr[ParamsSet::nSigns];
  	float minv[ParamsSet::nSigns];
  	float m1pt[ParamsSet::nSigns];
  	float m2pt[ParamsSet::nSigns];
  	float m1eta[ParamsSet::nSigns];
  	float m2eta[ParamsSet::nSigns];
  	float m1phi[ParamsSet::nSigns];
  	float m2phi[ParamsSet::nSigns];
    int m1charge[ParamsSet::nSigns];
    int m2charge[ParamsSet::nSigns];


    // --------------------- class methods ---------------------------

   	void InitInput();
    void InitOutput();
   	void InitHists();
   	void ProcessData();
    void WriteOutput();
    bool PassSingleMuonGapCut(float meta, float mpt, int mcharge);
   	void FillHistograms(int nsign);

public:
    int mode = 1;
    bool isScram;
    bool isTight;
    bool isRun3;
  	MuonPairPlottingPP();
  	~MuonPairPlottingPP(){}
  	void Run();

};

MuonPairPlottingPP::MuonPairPlottingPP(){
    // if (mode != 1 && mode != 3){
    //     std::cout<<"Error:: Mode has to be 1 (no binning) or 3 (binning by pT),  quitting"<<std::endl;
    //     throw std::exception();
    // }
    isScram = false;
    isTight = false;
    
    isRun3 = true;
}

void MuonPairPlottingPP::InitInput(){


    if (isScram){
        inFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/scrambled_muon_pairs_pp.root","read");
    }else if (isTight){
        inFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/muon_pairs_pp_tight.root","read");
    }else{
        if (isRun3){
            // inFile = new TFile("/eos/user/y/yuhang/data/pp_24/muon_pairs_pp_2024_part1.root","read");            
            inFile = new TFile("/eos/user/y/yuhang/data/pp_24/muon_pairs_pp_2024_part2.root","read");            
            // inFile = new TFile("/eos/user/y/yuhang/data/pp_24/muon_pairs_pp_2024_part3.root","read");            
        }else{
            inFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/muon_pairs_pp.root","read");            
        }
    }

    for (int isign = 0; isign < ParamsSet::nSigns; isign++){
        if (isScram){
            inTree[isign] = (TTree*) inFile->Get(Form("scramb_muon_pair_tree_sign%d",isign+1));  
        }else{
    	    inTree[isign] = (TTree*) inFile->Get(Form("muon_pair_tree_sign%d",isign+1));  
        }
        inTree[isign]->SetBranchAddress("weight"          , &weight[isign]);
    	inTree[isign]->SetBranchAddress("pair_dPoverP"    , &pair_dPoverP[isign]);
    	inTree[isign]->SetBranchAddress("pt_lead"          , &pt_lead[isign]);
    	inTree[isign]->SetBranchAddress("pair_pt"          , &pair_pt[isign]);
    	// inTree[isign]->SetBranchAddress("pair_eta"     , &pair_eta[isign]);
    	inTree[isign]->SetBranchAddress("pair_y"           , &pair_y[isign]);
        inTree[isign]->SetBranchAddress("asym"             , &asym[isign]);
    	inTree[isign]->SetBranchAddress("dpt"           , &dpt[isign]);
    	inTree[isign]->SetBranchAddress("deta"       , &deta[isign]);
    	inTree[isign]->SetBranchAddress("etaavg"      , &etaavg[isign]);
    	inTree[isign]->SetBranchAddress("phiavg"            , &phiavg[isign]);
    	inTree[isign]->SetBranchAddress("dphi"     , &dphi[isign]);
    	inTree[isign]->SetBranchAddress("dr"        , &dr[isign]);
    	inTree[isign]->SetBranchAddress("minv"        , &minv[isign]);
    	inTree[isign]->SetBranchAddress("m1.pt"           , &m1pt[isign]);
    	inTree[isign]->SetBranchAddress("m2.pt"           , &m2pt[isign]);
    	inTree[isign]->SetBranchAddress("m1.eta"       , &m1eta[isign]);
    	inTree[isign]->SetBranchAddress("m2.eta"       , &m2eta[isign]);
    	inTree[isign]->SetBranchAddress("m1.phi"     	, &m1phi[isign]);
    	inTree[isign]->SetBranchAddress("m2.phi"     	, &m2phi[isign]);
        inTree[isign]->SetBranchAddress("m1.charge"           , &m1charge[isign]);
        inTree[isign]->SetBranchAddress("m2.charge"           , &m2charge[isign]);

        inTree[isign]->SetBranchStatus("*"                  ,0);//switch off all branches, then enable just the ones that we need
        inTree[isign]->SetBranchStatus("weight"           ,1);
        inTree[isign]->SetBranchStatus("pair_dPoverP"           ,1);
        inTree[isign]->SetBranchStatus("pt_lead"           ,1);
        inTree[isign]->SetBranchStatus("pair_pt"           ,1);
        inTree[isign]->SetBranchStatus("pair_y"           ,1);
        inTree[isign]->SetBranchStatus("asym"           ,1);
        inTree[isign]->SetBranchStatus("dpt"           ,1);
        inTree[isign]->SetBranchStatus("deta"           ,1);
        inTree[isign]->SetBranchStatus("etaavg"           ,1);
        inTree[isign]->SetBranchStatus("phiavg"           ,1);
        inTree[isign]->SetBranchStatus("dphi"           ,1);
        inTree[isign]->SetBranchStatus("dr"           ,1);
        inTree[isign]->SetBranchStatus("minv"           ,1);
        inTree[isign]->SetBranchStatus("m1.pt"           ,1);
        inTree[isign]->SetBranchStatus("m2.pt"           ,1);
        inTree[isign]->SetBranchStatus("m1.eta"           ,1);
        inTree[isign]->SetBranchStatus("m2.eta"           ,1);
        inTree[isign]->SetBranchStatus("m1.phi"           ,1);
        inTree[isign]->SetBranchStatus("m2.phi"           ,1);
        inTree[isign]->SetBranchStatus("m1.charge"           ,1);
        inTree[isign]->SetBranchStatus("m2.charge"           ,1);

    }
}

void MuonPairPlottingPP::InitOutput(){
    // function to define the output file
    // this is needed since the histograms, once defined, belong to the last file before them
    // need to create the output files after creating TFile objects for the input file and between defining TH1D, TH2D objects for the histograms
    if (isScram){
        // outFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/histograms_scrambled_pairs_pp.root","update");
        outFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/histograms_scrambled_pairs_pp.root","recreate");
    }else if (isTight){
        outFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/histograms_real_pairs_pp_tight.root","recreate");
    }else{
        if (isRun3){
            // outFile = new TFile("/eos/user/y/yuhang/data/pp_24/histograms_real_pairs_pp_2024_part1.root","recreate");
            outFile = new TFile("/eos/user/y/yuhang/data/pp_24/histograms_real_pairs_pp_2024_part2.root","recreate");
            // outFile = new TFile("/eos/user/y/yuhang/data/pp_24/histograms_real_pairs_pp_2024_part3.root","recreate");
        }else{
            outFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/histograms_real_pairs_pp.root","recreate");
        }
    }
}

void MuonPairPlottingPP::InitHists(){
    int nDR_bins = 200;
    int nDphi_bins = 128;
    int neta_bins = 100;
    int nDeta_bins =  200;
    int npair_y_bins = 90;
    int npt_asym_bins = 100;
    int npair_pt_ptlead_ratio_bins = 100;
    int nminv_bins_linear = 100;
    int npair_pT_bins_linear = 100;
    int npT_lead_bins_linear = 100;

    static const int nminv_bins_log = 40;
    float minv_logpow[ParamsSet::nSigns];
    minv_logpow[0] = 0.062;
    minv_logpow[1] = 0.045;
    float minv_max = 60;
    double minv_bins_log[ParamsSet::nSigns][nminv_bins_log+1];

    std::string dphi_regions[2] = {"near", "away"};

    for (unsigned int isign = 0; isign < ParamsSet::nSigns; isign++){
        for(int iminv = 0; iminv <= nminv_bins_log; iminv++){
            minv_bins_log[isign][iminv] = minv_max * pow(10.0, (static_cast<float>(iminv - nminv_bins_log)) * minv_logpow[isign]);
        }
    }

    std::string ancestor_grps[nAncestorGroups + 1] = {"_from_same_b", "_gg", "_qg","_single_g","_qq"};

   	if (mode == 1){
        // for (unsigned int isign = 0; isign < ParamsSet::nSigns; isign++){
                
        // for (unsigned int ipt = 0; ipt < ParamsSet::nPtBins; ipt++){
        // for (unsigned int ictr = 0; ictr < ParamsSet::nCtrBins; ictr++){
        // for (unsigned int idphi = 0; idphi < ParamsSet::ndphiRegions; idphi++){
        // for (unsigned int ideta = 0; ideta < ParamsSet::ndetaRegions; ideta++){
        // for (unsigned int igap = 0; igap < ParamsSet::nGapCuts; igap++){
        // for (unsigned int iphoto = 0; iphoto < ParamsSet::nPhotoProdCuts; iphoto++){
        for (unsigned int idphi = 0; idphi < 2; idphi++){
            for (unsigned int isign = 0; isign < ParamsSet::nSigns; isign++){
                for (unsigned int igap = 0; igap < ParamsSet::nGapCuts; igap++){
                    h_pair_dP_overP[idphi][isign][igap] = new TH1D(Form("h_pair_dP_overP_%s_sign%d_gapcut%d",dphi_regions[idphi].c_str(),isign+1,igap+1),";(#Delta p / p)_{pair};1/N_{evt} dN/d(#Delta p / p)_{pair}" ,pms.deltaP_overP_nbins,0,pms.deltaP_overP_max);
                    h_DR[idphi][isign][igap] = new TH1D(Form("h_DR_%s_sign%d_gapcut%d",dphi_regions[idphi].c_str(),isign+1,igap+1),";#Delta R;1/N_{evt} dN/d#Delta R", nDR_bins,0,pms.deltaR_thrsh[2]);
                    h_Dphi[idphi][isign][igap] = new TH1D(Form("h_Dphi_%s_sign%d_gapcut%d",dphi_regions[idphi].c_str(),isign+1,igap+1),";#Delta#phi;1/N_{evt} dN/d#Delta#phi", nDphi_bins,-pms.PI,pms.PI);
                    h_pair_y[idphi][isign][igap] = new TH1D(Form("h_pair_y_%s_sign%d_gapcut%d",dphi_regions[idphi].c_str(),isign+1,igap+1),";y_{pair};1/N_{evt} dN/dy_{pair}" ,npair_y_bins,-3,3);
                    h_pt_asym[idphi][isign][igap] = new TH1D(Form("h_pt_asym_%s_sign%d_gapcut%d",dphi_regions[idphi].c_str(),isign+1,igap+1),";A = (pT1 - pT2)/(pT1 + pT2);d#sigma/dA", npt_asym_bins,0,1.);
                    h_pair_pt_ptlead_ratio[idphi][isign][igap] = new TH1D(Form("h_pair_pt_ptlead_ratio_%s_sign%d_gapcut%d",dphi_regions[idphi].c_str(),isign+1,igap+1),";#frac{p_{T}^{pair}}{p_{T}^{lead}};d#sigma/d#frac{p_{T}^{pair}}{p_{T}^{lead}}", npair_pt_ptlead_ratio_bins,0,2.);

                    h_eta_avg_Dphi[idphi][isign][igap] = new TH2D(Form("h_eta_avg_Dphi_%s_sign%d_gapcut%d",dphi_regions[idphi].c_str(),isign+1,igap+1),";#Delta#phi;#bar{#eta}", nDphi_bins,-pms.PI,pms.PI,neta_bins,-2.4,2.4);
                    h_Deta_Dphi[idphi][isign][igap] = new TH2D(Form("h_Deta_Dphi_%s_sign%d_gapcut%d",dphi_regions[idphi].c_str(),isign+1,igap+1),";#Delta#phi;#Delta#eta", nDphi_bins,-pms.PI,pms.PI,nDeta_bins,-4.8,4.8);
                    h_eta1_eta2[idphi][isign][igap] = new TH2D(Form("h_eta1_eta2_%s_sign%d_gapcut%d",dphi_regions[idphi].c_str(),isign+1,igap+1),";#eta_{sublead};#eta_{lead}",neta_bins,-2.4,2.4, neta_bins,-2.4,2.4);
                    h_eta_avg_Deta[idphi][isign][igap] = new TH2D(Form("h_eta_avg_Deta_%s_sign%d_gapcut%d",dphi_regions[idphi].c_str(),isign+1,igap+1),";#Delta#eta;#bar{#eta}",nDeta_bins,-4.8,4.8,neta_bins,-2.4,2.4);
                    h_pt1_pt2[idphi][isign][igap] = new TH2D(Form("h_pt1_pt2_%s_sign%d_gapcut%d",dphi_regions[idphi].c_str(),isign+1,igap+1),";p_{T}^{sublead} [GeV];p_{T}^{lead} [GeV]",pms.npt_bins,pms.pTBins,pms.npt_bins,pms.pTBins);
                    h_ptlead_pair_pt[idphi][isign][igap] = new TH2D(Form("h_ptlead_pair_pt_%s_sign%d_gapcut%d",dphi_regions[idphi].c_str(),isign+1,igap+1),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",npair_pT_bins_linear,0,30,npT_lead_bins_linear,0,30);
                    h_ptlead_pair_pt_zoomin[idphi][isign][igap] = new TH2D(Form("h_ptlead_pair_pt_zoomin_%s_sign%d_gapcut%d",dphi_regions[idphi].c_str(),isign+1,igap+1),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",npair_pT_bins_linear,0,20,npT_lead_bins_linear,0,15);
                    h_ptlead_pair_pt_log[idphi][isign][igap] = new TH2D(Form("h_ptlead_pair_pt_log_%s_sign%d_gapcut%d",dphi_regions[idphi].c_str(),isign+1,igap+1),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",pms.npairPT_bins,pms.pairPTBins[isign][2],pms.npt_bins,pms.pTBins);
                    h_minv_pair_pt[idphi][isign][igap] = new TH2D(Form("h_minv_pair_pt_%s_sign%d_gapcut%d",dphi_regions[idphi].c_str(),isign+1,igap+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",npair_pT_bins_linear,0,30,nminv_bins_linear,0,30);
                    h_minv_pair_pt_zoomin[idphi][isign][igap] = new TH2D(Form("h_minv_pair_pt_zoomin_%s_sign%d_gapcut%d",dphi_regions[idphi].c_str(),isign+1,igap+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",npair_pT_bins_linear,0,20,nminv_bins_linear,5,15);
                    h_minv_pair_pt_log[idphi][isign][igap] = new TH2D(Form("h_minv_pair_pt_log_%s_sign%d_gapcut%d",dphi_regions[idphi].c_str(),isign+1,igap+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",pms.npairPT_bins,pms.pairPTBins[isign][2],nminv_bins_log,minv_bins_log[isign]);
                }

                // for (unsigned int idphi= 0; idphi < ParamsSet::ndphiselcs; idphi++){
                //     // h_eta1_eta2_dphicut[idphi][isign][igap] = new TH2D(Form("h_eta1_eta2_dphi%d_sign%d_gapcut%d",idphi+1,isign+1,igap+1),";#eta_{sublead};#eta_{lead}",100,-2.4,2.4, 100,-2.4,2.4);
                //     h_eta1_eta2_dphicut[idphi][isign][igap] = new TH2D(Form("h_eta1_eta2_dphi%d_sign%d_gapcut%d",idphi+1,isign+1,igap+1),";#eta_{sublead};#eta_{lead}",50,-2.4,2.4, 50,-2.4,2.4);
                //     // if (isMCTruthBB || isMCTruthCC) h_unweighted_eta1_eta2_dphicut[idphi][isign][igap] = new TH2D(Form("h_unweighted_eta1_eta2_dphi%d_sign%d_gapcut%d",idphi+1,isign+1,igap+1),";#eta_{sublead};#eta_{lead}",50,-2.4,2.4, 50,-2.4,2.4);
                // }
            }
        }
    }
}


#endif