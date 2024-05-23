#include "MuonPairPlottingOLDPP.h"
#include "time.h"

void MuonPairPlottingOLDPP::ProcessData(){
	
  	for (int idr = 0; idr < ParamsSet::ndRselcs; idr++){
        for (int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
            
    	    Long64_t nentries = inTree[idr][ksign]->GetEntries(); //#muon pairs
            for (Long64_t lentry=0; lentry<nentries; lentry++) {//loop over the events
            // for (Long64_t lentry=0; lentry<100; lentry++) {//loop over the events
            if(lentry%100000==0) cout<<"Processing "<<lentry<<" event out of "<<nentries<<" events"<<std::endl;
          		int num_bytes = inTree[idr][ksign]->GetEntry(lentry);//read in an event
                if(num_bytes==0){
                	std::cout<<"Error:: Read in event has size of zero bytes,  quitting"<<std::endl;
                	throw std::exception();
                }
                if(mode == 1){
                	FillHistograms(idr, ksign);
                }else{ // mode has to be 1 or 3
                	// FillPtBinnedHistograms(idr, jpt, ksign);
                }
            }
        }
    }
}


bool MuonPairPlottingOLDPP::PassSingleMuonGapCut(float meta, float mpt, int mcharge){
    // parameters: eta and pT of the same muon
    if (fabs(meta) < pms.eta_gap_cut1) return false;
    if (mpt < 6){
        for (array<float,2> charge_eta_gap_cut : pms.charge_eta_gap_cuts){
            if (mcharge * meta > charge_eta_gap_cut[0] && mcharge * meta < charge_eta_gap_cut[1]) return false;
        }
    }
    // if (mpt < 6 && fabs(meta) > pms.eta_gap_cut2[0] && fabs(meta) < pms.eta_gap_cut2[1]) return false;
    return true;
}


void MuonPairPlottingOLDPP::FillHistograms(int ndr, int nsign){

    // ngapcut = 0: all; = 1: only those that pass

    //  bool pass_gapcut = fabs(m1eta[ndr][nsign]) > pms.eta_gap_cut && fabs(m2eta[ndr][nsign]) > pms.eta_gap_cut;
    bool pass_gapcut = PassSingleMuonGapCut(m1eta[ndr][nsign], m1pt[ndr][nsign], m1charge[ndr][nsign]) && PassSingleMuonGapCut(m2eta[ndr][nsign], m2pt[ndr][nsign], m2charge[ndr][nsign]);

    if (pass_gapcut){
        h_pair_dP_overP[ndr][nsign][1]->Fill(pair_dPoverP[ndr][nsign],weight[ndr][nsign]);
        h_pair_y[ndr][nsign][1]->Fill(pair_y[ndr][nsign],weight[ndr][nsign]);
        h_DR[ndr][nsign][1]->Fill(dr[ndr][nsign],weight[ndr][nsign]);
        h_Dphi[ndr][nsign][1]->Fill(dphi[ndr][nsign],weight[ndr][nsign]);
        h_eta_avg_Dphi[ndr][nsign][1]->Fill(dphi[ndr][nsign],etaavg[ndr][nsign],weight[ndr][nsign]);
        h_eta1_eta2[ndr][nsign][1]->Fill(m2eta[ndr][nsign],m1eta[ndr][nsign],weight[ndr][nsign]);
        h_pt1_pt2[ndr][nsign][1]->Fill(m2pt[ndr][nsign],m1pt[ndr][nsign],weight[ndr][nsign]);
        h_eta_avg_Deta[ndr][nsign][1]->Fill(deta[ndr][nsign],etaavg[ndr][nsign],weight[ndr][nsign]);
        h_ptlead_pair_pt[ndr][nsign][1]->Fill(pair_pt[ndr][nsign],m1pt[ndr][nsign],weight[ndr][nsign]);
        h_minv_pair_pt[ndr][nsign][1]->Fill(pair_pt[ndr][nsign],minv[ndr][nsign],weight[ndr][nsign]);
        h_Deta_Dphi[ndr][nsign][1]->Fill(dphi[ndr][nsign],deta[ndr][nsign],weight[ndr][nsign]);
            
        if (isMCTruthBB || isMCTruthCC){
            // h_unweighted_pair_dP_overP[ndr][nsign][1]->Fill(pair_dPoverP[ndr][nsign]);
            // h_unweighted_pair_y[ndr][nsign][1]->Fill(pair_y[ndr][nsign]);
            // h_unweighted_DR[ndr][nsign][1]->Fill(dr[ndr][nsign]);
            // h_unweighted_Dphi[ndr][nsign][1]->Fill(dphi[ndr][nsign]);
            // h_unweighted_eta_avg_Dphi[ndr][nsign][1]->Fill(dphi[ndr][nsign],etaavg[ndr][nsign]);
            // h_unweighted_eta1_eta2[ndr][nsign][1]->Fill(m2eta[ndr][nsign],m1eta[ndr][nsign]);
            // h_unweighted_pt1_pt2[ndr][nsign][1]->Fill(m2pt[ndr][nsign],m1pt[ndr][nsign]);
            // h_unweighted_eta_avg_Deta[ndr][nsign][1]->Fill(deta[ndr][nsign],etaavg[ndr][nsign]);
            // h_unweighted_ptlead_pair_pt[ndr][nsign][1]->Fill(pair_pt[ndr][nsign],m1pt[ndr][nsign]);
            // h_unweighted_minv_pair_pt[ndr][nsign][1]->Fill(pair_pt[ndr][nsign],minv[ndr][nsign]);
            h_unweighted_Deta_Dphi[ndr][nsign][1]->Fill(dphi[ndr][nsign],deta[ndr][nsign]);
        }
        
        if (ndr == 2){ //no delta R cut
            h_eta1_eta2_dphicut[2][nsign][1]->Fill(m2eta[ndr][nsign],m1eta[ndr][nsign],weight[ndr][nsign]);
            if (abs(dphi[ndr][nsign]) < 1){
                h_eta1_eta2_dphicut[0][nsign][1]->Fill(m2eta[ndr][nsign],m1eta[ndr][nsign],weight[ndr][nsign]);
            }else if(abs(dphi[ndr][nsign]) > pms.PI-1){
                h_eta1_eta2_dphicut[1][nsign][1]->Fill(m2eta[ndr][nsign],m1eta[ndr][nsign],weight[ndr][nsign]);
            }
        }
    }

    //fill [0] for everyone

    h_pair_dP_overP[ndr][nsign][0]->Fill(pair_dPoverP[ndr][nsign],weight[ndr][nsign]);
    h_pair_y[ndr][nsign][0]->Fill(pair_y[ndr][nsign],weight[ndr][nsign]);
    h_DR[ndr][nsign][0]->Fill(dr[ndr][nsign],weight[ndr][nsign]);
    h_Dphi[ndr][nsign][0]->Fill(dphi[ndr][nsign],weight[ndr][nsign]);
    h_eta_avg_Dphi[ndr][nsign][0]->Fill(dphi[ndr][nsign],etaavg[ndr][nsign],weight[ndr][nsign]);
    h_Deta_Dphi[ndr][nsign][0]->Fill(dphi[ndr][nsign],deta[ndr][nsign],weight[ndr][nsign]);
    h_eta1_eta2[ndr][nsign][0]->Fill(m2eta[ndr][nsign],m1eta[ndr][nsign],weight[ndr][nsign]);
    h_pt1_pt2[ndr][nsign][0]->Fill(m2pt[ndr][nsign],m1pt[ndr][nsign],weight[ndr][nsign]);
    h_eta_avg_Deta[ndr][nsign][0]->Fill(deta[ndr][nsign],etaavg[ndr][nsign],weight[ndr][nsign]);
    h_ptlead_pair_pt[ndr][nsign][0]->Fill(pair_pt[ndr][nsign],m1pt[ndr][nsign],weight[ndr][nsign]);
    h_minv_pair_pt[ndr][nsign][0]->Fill(pair_pt[ndr][nsign],minv[ndr][nsign],weight[ndr][nsign]);

    if (isMCTruthBB || isMCTruthCC){
        // h_unweighted_pair_dP_overP[ndr][nsign][0]->Fill(pair_dPoverP[ndr][nsign]);
        // h_unweighted_pair_y[ndr][nsign][0]->Fill(pair_y[ndr][nsign]);
        // h_unweighted_DR[ndr][nsign][0]->Fill(dr[ndr][nsign]);
        // h_unweighted_Dphi[ndr][nsign][0]->Fill(dphi[ndr][nsign]);
        // h_unweighted_eta_avg_Dphi[ndr][nsign][0]->Fill(dphi[ndr][nsign],etaavg[ndr][nsign]);
        h_unweighted_Deta_Dphi[ndr][nsign][0]->Fill(dphi[ndr][nsign],deta[ndr][nsign]);
        // h_unweighted_eta1_eta2[ndr][nsign][0]->Fill(m2eta[ndr][nsign],m1eta[ndr][nsign]);
        // h_unweighted_pt1_pt2[ndr][nsign][0]->Fill(m2pt[ndr][nsign],m1pt[ndr][nsign]);
        // h_unweighted_eta_avg_Deta[ndr][nsign][0]->Fill(deta[ndr][nsign],etaavg[ndr][nsign]);
        // h_unweighted_ptlead_pair_pt[ndr][nsign][0]->Fill(pair_pt[ndr][nsign],m1pt[ndr][nsign]);
        // h_unweighted_minv_pair_pt[ndr][nsign][0]->Fill(pair_pt[ndr][nsign],minv[ndr][nsign]);
    }


    if (ndr == 2){ //no delta R cut
        h_eta1_eta2_dphicut[2][nsign][0]->Fill(m2eta[ndr][nsign],m1eta[ndr][nsign],weight[ndr][nsign]);
        if (abs(dphi[ndr][nsign]) < 1){
            h_eta1_eta2_dphicut[0][nsign][0]->Fill(m2eta[ndr][nsign],m1eta[ndr][nsign],weight[ndr][nsign]);
        }else if(abs(dphi[ndr][nsign]) > pms.PI-1){
            h_eta1_eta2_dphicut[1][nsign][0]->Fill(m2eta[ndr][nsign],m1eta[ndr][nsign],weight[ndr][nsign]);
        }
    }
}

void MuonPairPlottingOLDPP::FillPtBinnedHistograms(int ndr, int npt, int nsign){}


void MuonPairPlottingOLDPP::WriteOutput(){
    if (isScram){
        // outFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/histograms_scrambled_pairs_pp_old.root","update");
        outFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/histograms_scrambled_pairs_pp_old.root","recreate");
    }else if (isTight){
        outFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/histograms_real_pairs_pp_tight_old.root","recreate");
    }else if(isMCTruthBB){
        outFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/athena/runMCV2/histograms_mc_truth_bb_old.root","recreate");
    }else if(isMCTruthCC){
        outFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/athena/runMCV2/histograms_mc_truth_cc_old.root","recreate");
    }else{
        // outFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/histograms_real_pairs_pp_no_resn_cut_old.root","recreate");
        // outFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/histograms_real_pairs_pp_no_resn_no_photo_old.root","recreate");
        outFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/histograms_real_pairs_pp_old.root","recreate");
        // outFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/histograms_real_pairs_pp_1s_only_old.root","recreate");
    }


    if (mode == 1){
        // outFile->Write();
        // outFile->cd();
        // gDirectory->Delete("h_*");
    
        for (unsigned int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
            for (unsigned int lgapcut = 0; lgapcut < ParamsSet::nGapCuts; lgapcut++){
                for (unsigned int idr = 0; idr < ParamsSet::ndRselcs; idr++){
                    // std::cout << "dR " << idr << ", sign " << ksign << "gapcut " << lgapcut << ", #entries " << h_pair_y[idr][ksign][lgapcut]->GetEntries() << std::endl;
                    h_pair_dP_overP[idr][ksign][lgapcut]->Write();
                    h_pair_y[idr][ksign][lgapcut]->Write();
                    h_DR[idr][ksign][lgapcut]->Write();
                    h_Dphi[idr][ksign][lgapcut]->Write();
                    h_eta_avg_Dphi[idr][ksign][lgapcut]->Write();
                    h_Deta_Dphi[idr][ksign][lgapcut]->Write();
                    h_eta1_eta2[idr][ksign][lgapcut]->Write();
                    h_eta_avg_Deta[idr][ksign][lgapcut]->Write();
                    h_pt1_pt2[idr][ksign][lgapcut]->Write();
                    h_ptlead_pair_pt[idr][ksign][lgapcut]->Write();
                    h_minv_pair_pt[idr][ksign][lgapcut]->Write();


                    if (isMCTruthBB || isMCTruthCC){
                        // h_unweighted_pair_dP_overP[idr][ksign][lgapcut]->Write();
                        // h_unweighted_pair_y[idr][ksign][lgapcut]->Write();
                        // h_unweighted_DR[idr][ksign][lgapcut]->Write();
                        // h_unweighted_Dphi[idr][ksign][lgapcut]->Write();
                        // h_unweighted_eta_avg_Dphi[idr][ksign][lgapcut]->Write();
                        h_unweighted_Deta_Dphi[idr][ksign][lgapcut]->Write();
                        // h_unweighted_eta1_eta2[idr][ksign][lgapcut]->Write();
                        // h_unweighted_eta_avg_Deta[idr][ksign][lgapcut]->Write();
                        // h_unweighted_pt1_pt2[idr][ksign][lgapcut]->Write();
                        // h_unweighted_ptlead_pair_pt[idr][ksign][lgapcut]->Write();
                        // h_unweighted_minv_pair_pt[idr][ksign][lgapcut]->Write();
                    }
                }
                for (unsigned int idphi= 0; idphi < ParamsSet::ndphiselcs; idphi++){
                    h_eta1_eta2_dphicut[idphi][ksign][lgapcut]->Write();
                }
            }
        }
    }else{ // mode = 3

    }
}

void MuonPairPlottingOLDPP::Run(){
    clock_t start, end;
    double cpu_time_used;
    start = clock();

    if (isScram) std::cout << "Data being processed: Scrambled" << std::endl;
    else if (isTight) std::cout << "Data being processed: Tight" << std::endl;
    else if (isMCTruthBB) std::cout << "Data being processed: MC Truth bb" << std::endl;
    else if (isMCTruthCC) std::cout << "Data being processed: MC Truth cc" << std::endl;
    else std::cout << "Data being processed: Real" << std::endl;

  	InitInput();
  	InitHists();
  	ProcessData();
  	WriteOutput();

    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;;
    std::cout << "CPU time used for running mode " << mode << " is " << cpu_time_used << " seconds" << std::endl;

}