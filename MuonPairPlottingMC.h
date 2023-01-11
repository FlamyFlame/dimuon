#ifndef __MuonPairPlottingMC_h__
#define __MuonPairPlottingMC_h__

#include "MuonPairPlottingPP.h"

class MuonPairPlottingMC : public MuonPairPlottingPP{
private:
  	ParamsSet pms;

  	// --------------------- NEW output file & histograms ---------------------------

  	TFile *outFile = nullptr;

    TH1D* h_pair_dP_overP[ParamsSet::ndRselcs][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH1D* h_Dphi[ParamsSet::ndRselcs][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH1D* h_DR[ParamsSet::ndRselcs][ParamsSet::nSigns][ParamsSet::nGapCuts];
    // TH1D* h_Minv[ParamsSet::ndRselcs][ParamsSet::nSigns][ParamsSet::nGapCuts];
    // TH1D* h_pt_lead[ParamsSet::ndRselcs][ParamsSet::nSigns][ParamsSet::nGapCuts];
    // TH1D* h_eta_avg[ParamsSet::ndRselcs][ParamsSet::nSigns][ParamsSet::nGapCuts];
    // TH1D* h_pair_pt[ParamsSet::ndRselcs][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH1D* h_pair_y[ParamsSet::ndRselcs][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_eta_avg_Dphi[ParamsSet::ndRselcs][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_Deta_Dphi[ParamsSet::ndRselcs][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_eta1_eta2[ParamsSet::ndRselcs][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_eta_avg_Deta[ParamsSet::ndRselcs][ParamsSet::nSigns][ParamsSet::nGapCuts];
    // TH2D* h_eta_avg_pair_eta[ParamsSet::ndRselcs][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_pt1_pt2[ParamsSet::ndRselcs][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_ptlead_pair_pt[ParamsSet::ndRselcs][ParamsSet::nSigns][ParamsSet::nGapCuts];
    TH2D* h_minv_pair_pt[ParamsSet::ndRselcs][ParamsSet::nSigns][ParamsSet::nGapCuts];

    TH2D* h_eta1_eta2_dphicut[ParamsSet::ndphiselcs][ParamsSet::nSigns][ParamsSet::nGapCuts];

    // --------------------- NEW data variables for setting branches---------------------------

   	// TFile *inFile;
    // TTree *inTree[ParamsSet::ndRselcs][ParamsSet::nSigns];

   //  float weight[ParamsSet::ndRselcs][ParamsSet::nSigns];
  	// float pair_dPoverP[ParamsSet::ndRselcs][ParamsSet::nSigns];
  	// float pt_lead[ParamsSet::ndRselcs][ParamsSet::nSigns];
  	// float pair_pt[ParamsSet::ndRselcs][ParamsSet::nSigns];
  	// float pair_y[ParamsSet::ndRselcs][ParamsSet::nSigns];
  	// float dpt[ParamsSet::ndRselcs][ParamsSet::nSigns];
  	// float deta[ParamsSet::ndRselcs][ParamsSet::nSigns];
  	// float etaavg[ParamsSet::ndRselcs][ParamsSet::nSigns];
  	// float phiavg[ParamsSet::ndRselcs][ParamsSet::nSigns];
  	// float dphi[ParamsSet::ndRselcs][ParamsSet::nSigns];
  	// float dr[ParamsSet::ndRselcs][ParamsSet::nSigns];
  	// float minv[ParamsSet::ndRselcs][ParamsSet::nSigns];
  	// float m1pt[ParamsSet::ndRselcs][ParamsSet::nSigns];
  	// float m2pt[ParamsSet::ndRselcs][ParamsSet::nSigns];
  	// float m1eta[ParamsSet::ndRselcs][ParamsSet::nSigns];
  	// float m2eta[ParamsSet::ndRselcs][ParamsSet::nSigns];
  	// float m1phi[ParamsSet::ndRselcs][ParamsSet::nSigns];
  	// float m2phi[ParamsSet::ndRselcs][ParamsSet::nSigns];
   //  int m1charge[ParamsSet::ndRselcs][ParamsSet::nSigns];
    // int m2charge[ParamsSet::ndRselcs][ParamsSet::nSigns];
    int m1_parent_group[ParamsSet::ndRselcs][ParamsSet::nSigns];
    int m2_parent_group[ParamsSet::ndRselcs][ParamsSet::nSigns];


    // --------------------- NEW/OVERRIDEN class methods ---------------------------

   	void InitInput();
   	void InitHists();
   	void ProcessData();
    // void WriteOutput();
    bool PassSingleMuonGapCut(float meta, float mpt, int mcharge);
   	void FillParentGroupedHistograms(int ndr, int nsign);
   	void FillPtBinnedHistograms(int ndr, int npt, int nsign);

public:
    int mode = 1;
    bool isScram = false;
    bool isMCTruthBB = false;
    bool isMCTruthCC = true;
    // assert(isScram & isMCTruthBB & isMCTruthCC == 1 || isScram & isMCTruthBB & isMCTruthCC == 1); // at most one can be true (if all false: real data)
  	MuonPairPlottingMC();
  	~MuonPairPlottingMC(){}
  	void Run();

};

MuonPairPlottingMC::MuonPairPlottingMC(){
    // if (mode != 1 && mode != 3){
    //     std::cout<<"Error:: Mode has to be 1 (no binning) or 3 (binning by pT),  quitting"<<std::endl;
    //     throw std::exception();
    // }
}

void MuonPairPlottingMC::InitInput(){

    if (isScram){
        inFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/scrambled_muon_pairs_pp.root","read");
    }else if(isMCTruthBB){
        inFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/athena/runMCV2/muon_pairs_mc_truth_bb.root","read");
    }else if(isMCTruthCC){
        inFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/athena/runMCV2/muon_pairs_mc_truth_cc.root","read");
    }else{
   	    inFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/muon_pairs_pp.root","read");
    }

    for (int idr = 0; idr < ParamsSet::ndRselcs; idr++){
        for (int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
            if (isScram){
                inTree[idr][ksign] = (TTree*) inFile->Get(Form("scramb_muon_pair_tree_dr%d_sign%d",idr+1,ksign+1));  
            }else{
        	    inTree[idr][ksign] = (TTree*) inFile->Get(Form("muon_pair_tree_dr%d_sign%d",idr+1,ksign+1));  
            }
            inTree[idr][ksign]->SetBranchAddress("weight"          , &weight[idr][ksign]);
        	inTree[idr][ksign]->SetBranchAddress("pair_dPoverP"           , &pair_dPoverP[idr][ksign]);
        	inTree[idr][ksign]->SetBranchAddress("pt_lead"          , &pt_lead[idr][ksign]);
        	inTree[idr][ksign]->SetBranchAddress("pair_pt"          , &pair_pt[idr][ksign]);
        	// inTree[idr][ksign]->SetBranchAddress("pair_eta"     , &pair_eta[idr][ksign]);
        	inTree[idr][ksign]->SetBranchAddress("pair_y"           , &pair_y[idr][ksign]);
        	inTree[idr][ksign]->SetBranchAddress("dpt"           , &dpt[idr][ksign]);
        	inTree[idr][ksign]->SetBranchAddress("deta"       , &deta[idr][ksign]);
        	inTree[idr][ksign]->SetBranchAddress("etaavg"      , &etaavg[idr][ksign]);
        	inTree[idr][ksign]->SetBranchAddress("phiavg"            , &phiavg[idr][ksign]);
        	inTree[idr][ksign]->SetBranchAddress("dphi"     , &dphi[idr][ksign]);
        	inTree[idr][ksign]->SetBranchAddress("dr"        , &dr[idr][ksign]);
        	inTree[idr][ksign]->SetBranchAddress("minv"        , &minv[idr][ksign]);
        	inTree[idr][ksign]->SetBranchAddress("m1.pt"           , &m1pt[idr][ksign]);
        	inTree[idr][ksign]->SetBranchAddress("m2.pt"           , &m2pt[idr][ksign]);
        	inTree[idr][ksign]->SetBranchAddress("m1.eta"       , &m1eta[idr][ksign]);
        	inTree[idr][ksign]->SetBranchAddress("m2.eta"       , &m2eta[idr][ksign]);
        	inTree[idr][ksign]->SetBranchAddress("m1.phi"     	, &m1phi[idr][ksign]);
        	inTree[idr][ksign]->SetBranchAddress("m2.phi"     	, &m2phi[idr][ksign]);
            inTree[idr][ksign]->SetBranchAddress("m1.charge"           , &m1charge[idr][ksign]);
            inTree[idr][ksign]->SetBranchAddress("m2.charge"           , &m2charge[idr][ksign]);

        }
    }
}

void MuonPairPlottingMC::InitHists(){

   	if (mode == 1){
        for (unsigned int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
            for (unsigned int lgapcut = 0; lgapcut < ParamsSet::nGapCuts; lgapcut++){
   	            for (unsigned int idr = 0; idr < ParamsSet::ndRselcs; idr++){
                    h_pair_dP_overP[idr][ksign][lgapcut] = new TH1D(Form("h_pair_dP_overP_dr%d_sign%d_gapcut%d",idr+1,ksign+1,lgapcut+1),";(#Delta p / p)_{pair};1/N_{evt} dN/d(#Delta p / p)_{pair}" ,pms.deltaP_overP_nbins,0,pms.deltaP_overP_max);
                    h_DR[idr][ksign][lgapcut] = new TH1D(Form("h_DR_dr%d_sign%d_gapcut%d",idr+1,ksign+1,lgapcut+1),";#Delta R;1/N_{evt} dN/d#Delta R", pms.deltaR_nbins[idr],0,pms.deltaR_thrsh[idr]);
                    h_Dphi[idr][ksign][lgapcut] = new TH1D(Form("h_Dphi_dr%d_sign%d_gapcut%d",idr+1,ksign+1,lgapcut+1),";#Delta#phi;1/N_{evt} dN/d#Delta#phi", 128,-pms.PI,pms.PI);
                    h_pair_y[idr][ksign][lgapcut] = new TH1D(Form("h_pair_y_dr%d_sign%d_gapcut%d",idr+1,ksign+1,lgapcut+1),";y_{pair};1/N_{evt} dN/dy_{pair}" ,90,-3,3);
                    h_eta_avg_Dphi[idr][ksign][lgapcut] = new TH2D(Form("h_eta_avg_Dphi_dr%d_sign%d_gapcut%d",idr+1,ksign+1,lgapcut+1),";#Delta#phi;#bar{#eta}", 128,-pms.PI,pms.PI,100,-2.4,2.4);
                    h_Deta_Dphi[idr][ksign][lgapcut] = new TH2D(Form("h_Deta_Dphi_dr%d_sign%d_gapcut%d",idr+1,ksign+1,lgapcut+1),";#Delta#phi;#Delta#eta", 128,-pms.PI,pms.PI,200,-4.8,4.8);
                    h_eta1_eta2[idr][ksign][lgapcut] = new TH2D(Form("h_eta1_eta2_dr%d_sign%d_gapcut%d",idr+1,ksign+1,lgapcut+1),";#eta_{sublead};#eta_{lead}",100,-2.4,2.4, 100,-2.4,2.4);
                    h_eta_avg_Deta[idr][ksign][lgapcut] = new TH2D(Form("h_eta_avg_Deta_dr%d_sign%d_gapcut%d",idr+1,ksign+1,lgapcut+1),";#Delta#eta;#bar{#eta}",200,-4.8,4.8,100,-2.4,2.4);
                    // h_eta_avg_pair_eta[idr][ksign][lgapcut] = new TH2D(Form("h_eta_avg_pair_eta_dr%d_sign%d_gapcut%d",idr+1,ksign+1,lgapcut+1),";#eta_{pair};#bar{#eta}",100,-2.4,2.4, 100,-2.4,2.4);
                    h_pt1_pt2[idr][ksign][lgapcut] = new TH2D(Form("h_pt1_pt2_dr%d_sign%d_gapcut%d",idr+1,ksign+1,lgapcut+1),";p_{T}^{sublead} [GeV];p_{T}^{lead} [GeV]",pms.npt_bins,pms.pTBins,pms.npt_bins,pms.pTBins);
                    h_ptlead_pair_pt[idr][ksign][lgapcut] = new TH2D(Form("h_ptlead_pair_pt_dr%d_sign%d_gapcut%d",idr+1,ksign+1,lgapcut+1),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",pms.npairPT_bins,pms.pairPTBins[ksign][idr],pms.npt_bins,pms.pTBins);
                    h_minv_pair_pt[idr][ksign][lgapcut] = new TH2D(Form("h_minv_pair_pt_dr%d_sign%d_gapcut%d",idr+1,ksign+1,lgapcut+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",pms.npairPT_bins,pms.pairPTBins[ksign][idr],pms.minv_nbins[idr],0,pms.minv_max[idr]);
                    // double minv_arr[pms.minv_nbins[idr]+1];
                    // std::copy(pms.minv_bins[ksign][idr].begin(),pms.minv_bins[ksign][idr].end(),minv_arr);
                    // h_minv_pair_pt[idr][ksign][lgapcut] = new TH2D(Form("h_minv_pair_pt_dr%d_sign%d_gapcut%d",idr+1,ksign+1,lgapcut+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",pms.npairPT_bins,pms.pairPTBins[ksign][idr],pms.minv_nbins[idr],minv_arr);
                }

                for (unsigned int idphi= 0; idphi < ParamsSet::ndphiselcs; idphi++){
                    h_eta1_eta2_dphicut[idphi][ksign][lgapcut] = new TH2D(Form("h_eta1_eta2_dphi%d_sign%d_gapcut%d",idphi+1,ksign+1,lgapcut+1),";#eta_{sublead};#eta_{lead}",100,-2.4,2.4, 100,-2.4,2.4);
                }
            }
        }
    }
}



#endif