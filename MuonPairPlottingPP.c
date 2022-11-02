#include "MuonPairPlottingPP.h"
#include "time.h"

void MuonPairPlottingPP::ProcessData(){
	
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


void MuonPairPlottingPP::FillHistograms(int ndr, int nsign){

    // ngapcut = 0: all; = 1: only those that pass

    bool pass_gapcut = fabs(m1eta[ndr][nsign]) > pms.eta_gap_cut && fabs(m2eta[ndr][nsign]) > pms.eta_gap_cut;
    if (pass_gapcut){
        h_pair_dP_overP[ndr][nsign][1]->Fill(pair_dPoverP[ndr][nsign]);
        h_pair_y[ndr][nsign][1]->Fill(pair_y[ndr][nsign]);
        h_DR[ndr][nsign][1]->Fill(dr[ndr][nsign]);
        h_eta_avg_Dphi[ndr][nsign][1]->Fill(dphi[ndr][nsign],etaavg[ndr][nsign]);
        h_eta1_eta2[ndr][nsign][1]->Fill(m2eta[ndr][nsign],m1eta[ndr][nsign]);
        h_pt1_pt2[ndr][nsign][1]->Fill(m2pt[ndr][nsign],m1pt[ndr][nsign]);
        h_eta_avg_Deta[ndr][nsign][1]->Fill(deta[ndr][nsign],etaavg[ndr][nsign]);
        h_ptlead_pair_pt[ndr][nsign][1]->Fill(pair_pt[ndr][nsign],m1pt[ndr][nsign]);
        h_minv_pair_pt[ndr][nsign][1]->Fill(pair_pt[ndr][nsign],minv[ndr][nsign]);
        h_Deta_Dphi[ndr][nsign][1]->Fill(dphi[ndr][nsign],deta[ndr][nsign]);
        h_eta_avg_pair_eta[ndr][nsign][1]->Fill(pair_eta[ndr][nsign],etaavg[ndr][nsign]);


        if (ndr == 2){ //no delta R cut
            h_eta1_eta2_dphicut[2][nsign][1]->Fill(m2eta[ndr][nsign],m1eta[ndr][nsign]);
            if (dphi[ndr][nsign] < 1){
                h_eta1_eta2_dphicut[0][nsign][1]->Fill(m2eta[ndr][nsign],m1eta[ndr][nsign]);
            }else if(dphi[ndr][nsign] > pms.PI-1){
                h_eta1_eta2_dphicut[1][nsign][1]->Fill(m2eta[ndr][nsign],m1eta[ndr][nsign]);
            }
        }
    }

    //fill [0] for everyone
    h_pair_dP_overP[ndr][nsign][0]->Fill(pair_dPoverP[ndr][nsign]);
    h_pair_y[ndr][nsign][0]->Fill(pair_y[ndr][nsign]);
    h_DR[ndr][nsign][0]->Fill(dr[ndr][nsign]);
    h_eta_avg_Dphi[ndr][nsign][0]->Fill(dphi[ndr][nsign],etaavg[ndr][nsign]);
    h_Deta_Dphi[ndr][nsign][0]->Fill(dphi[ndr][nsign],deta[ndr][nsign]);
    h_eta1_eta2[ndr][nsign][0]->Fill(m2eta[ndr][nsign],m1eta[ndr][nsign]);
    h_pt1_pt2[ndr][nsign][0]->Fill(m2pt[ndr][nsign],m1pt[ndr][nsign]);
    h_eta_avg_Deta[ndr][nsign][0]->Fill(deta[ndr][nsign],etaavg[ndr][nsign]);
    h_eta_avg_pair_eta[ndr][nsign][0]->Fill(pair_eta[ndr][nsign],etaavg[ndr][nsign]);
    h_ptlead_pair_pt[ndr][nsign][0]->Fill(pair_pt[ndr][nsign],m1pt[ndr][nsign]);
    h_minv_pair_pt[ndr][nsign][0]->Fill(pair_pt[ndr][nsign],minv[ndr][nsign]);

    if (ndr == 2){ //no delta R cut
        h_eta1_eta2_dphicut[2][nsign][0]->Fill(m2eta[ndr][nsign],m1eta[ndr][nsign]);
        if (dphi[ndr][nsign] < 1){
            h_eta1_eta2_dphicut[0][nsign][0]->Fill(m2eta[ndr][nsign],m1eta[ndr][nsign]);
        }else if(dphi[ndr][nsign] > pms.PI-1){
            h_eta1_eta2_dphicut[1][nsign][0]->Fill(m2eta[ndr][nsign],m1eta[ndr][nsign]);
        }
    }
}

void MuonPairPlottingPP::FillPtBinnedHistograms(int ndr, int npt, int nsign){}


void MuonPairPlottingPP::WriteOutput(){
    if (isScram){
        // outFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/histograms_scrambled_pairs_pp.root","update");
        outFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/histograms_scrambled_pairs_pp.root","recreate");
    }else{
        // outFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/histograms_real_pairs_pp.root","update");
        outFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/histograms_real_pairs_pp.root","recreate");
    }


    if (mode == 1){
        // outFile->Write();
        // outFile->cd();
        // gDirectory->Delete("h_*");
    
        for (unsigned int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
            for (unsigned int lgapcut = 0; lgapcut < ParamsSet::nGapCuts; lgapcut++){
                for (unsigned int idr = 0; idr < ParamsSet::ndRselcs; idr++){
                    std::cout << "dR " << idr << ", sign " << ksign << "gapcut " << lgapcut << ", #entries " << h_pair_y[idr][ksign][lgapcut]->GetEntries() << std::endl;
                    h_pair_dP_overP[idr][ksign][lgapcut]->Write();
                    h_pair_y[idr][ksign][lgapcut]->Write();
                    h_DR[idr][ksign][lgapcut]->Write();
                    h_eta_avg_Dphi[idr][ksign][lgapcut]->Write();
                    h_Deta_Dphi[idr][ksign][lgapcut]->Write();
                    h_eta1_eta2[idr][ksign][lgapcut]->Write();
                    h_eta_avg_Deta[idr][ksign][lgapcut]->Write();
                    h_eta_avg_pair_eta[idr][ksign][lgapcut]->Write();
                    h_pt1_pt2[idr][ksign][lgapcut]->Write();
                    h_ptlead_pair_pt[idr][ksign][lgapcut]->Write();
                    h_minv_pair_pt[idr][ksign][lgapcut]->Write();
                }
                for (unsigned int idphi= 0; idphi < ParamsSet::ndphiselcs; idphi++){
                    h_eta1_eta2_dphicut[idphi][ksign][lgapcut]->Write();
                }
            }
        }
    }else{ // mode = 3

    }
}

void MuonPairPlottingPP::Run(){
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