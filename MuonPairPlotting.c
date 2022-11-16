#include "MuonPairPlotting.h"
#include "time.h"

void MuonPairPlotting::ProcessData(){
	
    int cur_ctr_intvl = 0;
    for (int jctr = 0; jctr < nCtrBins; jctr++){
        for (unsigned int jj = 0; jj < ctrBins[jctr].size(); jj++){
  	        for (int idr = 0; idr < ParamsSet::ndRselcs; idr++){
                for (int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
                    
        	        Long64_t nentries = inTree[idr][cur_ctr_intvl][ksign]->GetEntries(); //#muon pairs
                    for (Long64_t lentry=0; lentry<nentries; lentry++) {//loop over the events
                    // for (Long64_t lentry=0; lentry<100; lentry++) {//loop over the events
                    if(lentry%100000==0) cout<<"Processing "<<lentry<<" event out of "<<nentries<<" events"<<std::endl;

                  		int num_bytes = inTree[idr][cur_ctr_intvl][ksign]->GetEntry(lentry);//read in an event
                        if(num_bytes==0){
                        	std::cout<<"Error:: Read in event has size of zero bytes,  quitting"<<std::endl;
                        	throw std::exception();
                        }
                        if(mode == 1){
                        	FillUnbinnedHistograms(idr, cur_ctr_intvl, ksign);
                        }else if(mode == 2){
                        	FillCtrBinnedHistograms(idr, jctr, cur_ctr_intvl, ksign);
                        }else if(mode == 3){
                        	FillPtBinnedHistograms(idr, cur_ctr_intvl, ksign);
                        }
                    }
                }
            }
            cur_ctr_intvl++;
        }
    }
    assert(cur_ctr_intvl == 16);
}

void MuonPairPlotting::FillUnbinnedHistograms(int ndr, int nctr_intvl, int nsign){
    h_pair_dP_overP[ndr][nsign]->Fill(pair_dPoverP[ndr][nctr_intvl][nsign]);
    // h_Minv[ndr][nsign]   ->Fill(minv[ndr][nctr_intvl][nsign]);
    // h_Dphi[ndr][nsign]   ->Fill(dphi[ndr][nctr_intvl][nsign]);
    // h_Deta[ndr][nsign]   ->Fill(deta[ndr][nctr_intvl][nsign]);
    // h_DR[ndr][nsign]     ->Fill(dr[ndr][nctr_intvl][nsign]);
    // h_pt_lead[ndr][nsign] ->Fill(m1pt[ndr][nctr_intvl][nsign]);
    // h_eta_avg[ndr][nsign] ->Fill(etaavg[ndr][nctr_intvl][nsign]);
    // h_pair_pt[ndr][nsign]->Fill(pair_pt[ndr][nctr_intvl][nsign]);
    // h_pair_eta[ndr][nsign]->Fill(pair_eta[ndr][nctr_intvl][nsign]);
    h_pair_y[ndr][nsign]->Fill(pair_y[ndr][nctr_intvl][nsign]);

    h_eta_avg_Dphi[ndr][nsign]->Fill(phiavg[ndr][nctr_intvl][nsign],etaavg[ndr][nctr_intvl][nsign]);
    h_eta1_eta2[ndr][nsign]->Fill(m2eta[ndr][nctr_intvl][nsign],m1eta[ndr][nctr_intvl][nsign]);
    h_pt1_pt2[ndr][nsign]->Fill(m2pt[ndr][nctr_intvl][nsign],m1pt[ndr][nctr_intvl][nsign]);
    h_eta_avg_Deta[ndr][nsign]->Fill(deta[ndr][nctr_intvl][nsign],etaavg[ndr][nctr_intvl][nsign]);  
    h_eta_avg_pair_eta[ndr][nsign]->Fill(pair_eta[ndr][nctr_intvl][nsign],etaavg[ndr][nctr_intvl][nsign]);
    h_ptlead_pair_pt[ndr][nsign]->Fill(pair_pt[ndr][nctr_intvl][nsign],m1pt[ndr][nctr_intvl][nsign]);
    h_minv_pair_pt[ndr][nsign]->Fill(pair_pt[ndr][nctr_intvl][nsign],minv[ndr][nctr_intvl][nsign]);
}


bool MuonPairPlotting::PassSingleMuonGapCut(float meta, float mpt, int mcharge){
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

void MuonPairPlotting::FillCtrBinnedHistograms(int ndr, int nctr_bin, int nctr_intvl, int nsign){
    // ngapcut = 0: all; = 1: only those that pass

    //if pass eta gap cut
    bool pass_gapcut = PassSingleMuonGapCut(m1eta[ndr][nctr_intvl][nsign], m1pt[ndr][nctr_intvl][nsign], m1charge[ndr][nctr_intvl][nsign]) && PassSingleMuonGapCut(m2eta[ndr][nctr_intvl][nsign], m2pt[ndr][nctr_intvl][nsign], m2charge[ndr][nctr_intvl][nsign]);
    if (pass_gapcut){
        h_ctrbin_pair_dP_overP[ndr][nctr_bin][nsign][1]->Fill(pair_dPoverP[ndr][nctr_intvl][nsign]);
        h_ctrbin_pair_y[ndr][nctr_bin][nsign][1]->Fill(pair_y[ndr][nctr_intvl][nsign]);
        h_ctrbin_DR[ndr][nctr_bin][nsign][1]->Fill(dr[ndr][nctr_intvl][nsign]);
        h_ctrbin_eta_avg_Dphi[ndr][nctr_bin][nsign][1]->Fill(dphi[ndr][nctr_intvl][nsign],etaavg[ndr][nctr_intvl][nsign]);
        h_ctrbin_Deta_Dphi[ndr][nctr_bin][nsign][1]->Fill(dphi[ndr][nctr_intvl][nsign],deta[ndr][nctr_intvl][nsign]);
        h_ctrbin_eta1_eta2[ndr][nctr_bin][nsign][1]->Fill(m2eta[ndr][nctr_intvl][nsign],m1eta[ndr][nctr_intvl][nsign]);
        h_ctrbin_pt1_pt2[ndr][nctr_bin][nsign][1]->Fill(m2pt[ndr][nctr_intvl][nsign],m1pt[ndr][nctr_intvl][nsign]);
        h_ctrbin_eta_avg_Deta[ndr][nctr_bin][nsign][1]->Fill(deta[ndr][nctr_intvl][nsign],etaavg[ndr][nctr_intvl][nsign]);
        h_ctrbin_ptlead_pair_pt[ndr][nctr_bin][nsign][1]->Fill(pair_pt[ndr][nctr_intvl][nsign],m1pt[ndr][nctr_intvl][nsign]);
        h_ctrbin_minv_pair_pt[ndr][nctr_bin][nsign][1]->Fill(pair_pt[ndr][nctr_intvl][nsign],minv[ndr][nctr_intvl][nsign]);
        h_ctrbin_eta_avg_pair_eta[ndr][nctr_bin][nsign][1]->Fill(pair_eta[ndr][nctr_intvl][nsign],etaavg[ndr][nctr_intvl][nsign]);

        if (ndr == 2){ //no delta R cut
            h_ctrbin_eta1_eta2_dphicut[2][nctr_bin][nsign][1]->Fill(m2eta[ndr][nctr_intvl][nsign],m1eta[ndr][nctr_intvl][nsign]);
            if (dphi[ndr][nctr_intvl][nsign] < 1){
                h_ctrbin_eta1_eta2_dphicut[0][nctr_bin][nsign][1]->Fill(m2eta[ndr][nctr_intvl][nsign],m1eta[ndr][nctr_intvl][nsign]);
            }else if(dphi[ndr][nctr_intvl][nsign] > pms.PI-1){
                h_ctrbin_eta1_eta2_dphicut[1][nctr_bin][nsign][1]->Fill(m2eta[ndr][nctr_intvl][nsign],m1eta[ndr][nctr_intvl][nsign]);
            }
        }
    }else{
        // h_ctrbin_pair_eta_Dphi_failgapcut[ndr][nctr_bin][nsign]->Fill(dphi[ndr][nctr_intvl][nsign],pair_eta[ndr][nctr_intvl][nsign]);
        // h_ctrbin_Dphi_failgapcut[ndr][nctr_bin][nsign]->Fill(dphi[ndr][nctr_intvl][nsign]);
    }

    //fill [1] for everyone
    h_ctrbin_pair_dP_overP[ndr][nctr_bin][nsign][0]->Fill(pair_dPoverP[ndr][nctr_intvl][nsign]);
    h_ctrbin_pair_y[ndr][nctr_bin][nsign][0]->Fill(pair_y[ndr][nctr_intvl][nsign]);
    h_ctrbin_DR[ndr][nctr_bin][nsign][0]->Fill(dr[ndr][nctr_intvl][nsign]);
    h_ctrbin_Deta_Dphi[ndr][nctr_bin][nsign][0]->Fill(dphi[ndr][nctr_intvl][nsign],deta[ndr][nctr_intvl][nsign]);
    h_ctrbin_eta_avg_Dphi[ndr][nctr_bin][nsign][0]->Fill(dphi[ndr][nctr_intvl][nsign],etaavg[ndr][nctr_intvl][nsign]);
    h_ctrbin_eta1_eta2[ndr][nctr_bin][nsign][0]->Fill(m2eta[ndr][nctr_intvl][nsign],m1eta[ndr][nctr_intvl][nsign]);
    h_ctrbin_pt1_pt2[ndr][nctr_bin][nsign][0]->Fill(m2pt[ndr][nctr_intvl][nsign],m1pt[ndr][nctr_intvl][nsign]);
    h_ctrbin_eta_avg_Deta[ndr][nctr_bin][nsign][0]->Fill(deta[ndr][nctr_intvl][nsign],etaavg[ndr][nctr_intvl][nsign]);
    h_ctrbin_eta_avg_pair_eta[ndr][nctr_bin][nsign][0]->Fill(pair_eta[ndr][nctr_intvl][nsign],etaavg[ndr][nctr_intvl][nsign]);
    h_ctrbin_ptlead_pair_pt[ndr][nctr_bin][nsign][0]->Fill(pair_pt[ndr][nctr_intvl][nsign],m1pt[ndr][nctr_intvl][nsign]);
    h_ctrbin_minv_pair_pt[ndr][nctr_bin][nsign][0]->Fill(pair_pt[ndr][nctr_intvl][nsign],minv[ndr][nctr_intvl][nsign]);

    if (ndr == 2){ //no delta R cut
        h_ctrbin_eta1_eta2_dphicut[2][nctr_bin][nsign][0]->Fill(m2eta[ndr][nctr_intvl][nsign],m1eta[ndr][nctr_intvl][nsign]);
        if (dphi[ndr][nctr_intvl][nsign] < 1){
            h_ctrbin_eta1_eta2_dphicut[0][nctr_bin][nsign][0]->Fill(m2eta[ndr][nctr_intvl][nsign],m1eta[ndr][nctr_intvl][nsign]);
        }else if(dphi[ndr][nctr_intvl][nsign] > pms.PI-1){
            h_ctrbin_eta1_eta2_dphicut[1][nctr_bin][nsign][0]->Fill(m2eta[ndr][nctr_intvl][nsign],m1eta[ndr][nctr_intvl][nsign]);
        }
    }
}

void MuonPairPlotting::FillPtBinnedHistograms(int ndr, int nctr_intvl, int nsign){}


void MuonPairPlotting::WriteOutput(){
    if (isScram){
        outFile = new TFile(("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/histograms_scrambled_pairs" + postfix + ".root").c_str(),"update");
    }else{
        outFile = new TFile(("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/histograms_real_pairs" + postfix + ".root").c_str(),"update");
    }


    if (mode == 2){
        outFile->cd();
        // if (outFile->GetDirectory("ctr-binned") == nullptr){}
        if (gDirectory->GetDirectory("ctr-binned") == nullptr){
          gDirectory->mkdir("ctr-binned");
        }
        gDirectory->cd("ctr-binned");
        gDirectory->Delete("h_*");
    
        for (unsigned int jctr = 0; jctr < nCtrBins; jctr++){
            for (unsigned int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
                for (unsigned int idr = 0; idr < ParamsSet::ndRselcs; idr++){
                    // h_ctrbin_pair_eta_Dphi_failgapcut[idr][jctr][ksign]->Write();
                    // h_ctrbin_Dphi_failgapcut[idr][jctr][ksign]->Write();
                    for (unsigned int lgapcut = 0; lgapcut < ParamsSet::nGapCuts; lgapcut++){
                        h_ctrbin_pair_dP_overP[idr][jctr][ksign][lgapcut]->Write();
                        h_ctrbin_pair_y[idr][jctr][ksign][lgapcut]->Write();
                        h_ctrbin_DR[idr][jctr][ksign][lgapcut]->Write();
                        h_ctrbin_Deta_Dphi[idr][jctr][ksign][lgapcut]->Write();
                        h_ctrbin_eta_avg_Dphi[idr][jctr][ksign][lgapcut]->Write();
                        h_ctrbin_eta1_eta2[idr][jctr][ksign][lgapcut]->Write();
                        h_ctrbin_eta_avg_Deta[idr][jctr][ksign][lgapcut]->Write();
                        h_ctrbin_eta_avg_pair_eta[idr][jctr][ksign][lgapcut]->Write();
                        h_ctrbin_pt1_pt2[idr][jctr][ksign][lgapcut]->Write();
                        h_ctrbin_ptlead_pair_pt[idr][jctr][ksign][lgapcut]->Write();
                        h_ctrbin_minv_pair_pt[idr][jctr][ksign][lgapcut]->Write();
                    }
                }
                
                for (unsigned int idphi= 0; idphi < ParamsSet::ndphiselcs; idphi++){
                    for (unsigned int lgapcut = 0; lgapcut < ParamsSet::nGapCuts; lgapcut++){
                        h_ctrbin_eta1_eta2_dphicut[idphi][jctr][ksign][lgapcut]->Write();
                    }
                }
            }            
        }
    }else if (mode == 1){
        outFile->cd();
        if (gDirectory->Get("h_*") != nullptr) gDirectory->Delete("h_*");
    
        for (unsigned int idr = 0; idr < pms.ndRselcs; idr++){
            for (unsigned int ksign = 0; ksign < pms.nSigns; ksign++){
                h_pair_dP_overP[idr][ksign]->Write();
                h_DR[idr][ksign]->Write();
                h_pair_y[idr][ksign]->Write();
                h_eta_avg_Dphi[idr][ksign]->Write();
                h_Deta_Dphi[idr][ksign]->Write();
                h_eta1_eta2[idr][ksign]->Write();
                h_pt1_pt2[idr][ksign]->Write();
                h_eta_avg_Deta[idr][ksign]->Write();
                h_eta_avg_pair_eta[idr][ksign]->Write();
                h_ptlead_pair_pt[idr][ksign]->Write();
                h_minv_pair_pt[idr][ksign]->Write();
            }
        }
    }else if (mode == 3){

    }else{
        std::cout << "Mode = " << mode << "; mode should take a value in 1, 2, 3." << std::endl;
    }
}

void MuonPairPlotting::Run(){
    clock_t start, end;
    double cpu_time_used;
    start = clock();

  	InitInput();
  	InitHists();
  	ProcessData();
  	WriteOutput();

    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;;
    std::cout << "CPU time used for running mode " << mode << " is " << cpu_time_used << " seconds" << std::endl;

}