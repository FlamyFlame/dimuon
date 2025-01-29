#include "MuonPairPlottingPP.h"
#include "time.h"

void MuonPairPlottingPP::ProcessData(){
	
    for (int isign = 0; isign < ParamsSet::nSigns; isign++){
        
        Long64_t nentries = inTree[isign]->GetEntries(); //#muon pairs
        for (Long64_t lentry=0; lentry<nentries; lentry++) {//loop over the events
        // for (Long64_t lentry=0; lentry<100; lentry++) {//loop over the events
        if(lentry%100000==0) cout<<"Processing "<<lentry<<" event out of "<<nentries<<" events"<<std::endl;
      		int num_bytes = inTree[isign]->GetEntry(lentry);//read in an event
            if(num_bytes==0){
            	std::cout<<"Error:: Read in event has size of zero bytes,  quitting"<<std::endl;
            	throw std::exception();
            }
            if(mode == 1){
            	FillHistograms(isign);
            }else{ // mode has to be 1 or 3
            	// FillPtBinnedHistograms(idr, jpt, isign);
            }
        }
    }
}


bool MuonPairPlottingPP::PassSingleMuonGapCut(float meta, float mpt, int mcharge){
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


void MuonPairPlottingPP::FillHistograms(int nsign){

    // ngapcut = 0: all; = 1: only those that pass

    //  bool pass_gapcut = fabs(m1eta[nsign]) > pms.eta_gap_cut && fabs(m2eta[nsign]) > pms.eta_gap_cut;
    bool pass_gapcut = PassSingleMuonGapCut(m1eta[nsign], m1pt[nsign], m1charge[nsign]) && PassSingleMuonGapCut(m2eta[nsign], m2pt[nsign], m2charge[nsign]);
    bool away_side = (abs(dphi[nsign]) >= pms.PI / 2.);

    if (pass_gapcut){
        h_pair_dP_overP[away_side][nsign][1]->Fill(pair_dPoverP[nsign],weight[nsign]);
        h_pair_y[away_side][nsign][1]->Fill(pair_y[nsign],weight[nsign]);
        h_DR[away_side][nsign][1]->Fill(dr[nsign],weight[nsign]);
        h_Dphi[away_side][nsign][1]->Fill(dphi[nsign],weight[nsign]);
        h_pt_asym[away_side][nsign][1]->Fill(asym[nsign],weight[nsign]);
        h_pair_pt_ptlead_ratio[away_side][nsign][1]->Fill(pair_pt[nsign]/pt_lead[nsign],weight[nsign]);

        h_eta_avg_Dphi[away_side][nsign][1]->Fill(dphi[nsign],etaavg[nsign],weight[nsign]);
        h_eta1_eta2[away_side][nsign][1]->Fill(m2eta[nsign],m1eta[nsign],weight[nsign]);
        h_pt1_pt2[away_side][nsign][1]->Fill(m2pt[nsign],m1pt[nsign],weight[nsign]);
        h_eta_avg_Deta[away_side][nsign][1]->Fill(deta[nsign],etaavg[nsign],weight[nsign]);
        h_Deta_Dphi[away_side][nsign][1]->Fill(dphi[nsign],deta[nsign],weight[nsign]);
        h_ptlead_pair_pt[away_side][nsign][1]->Fill(pair_pt[nsign],m1pt[nsign],weight[nsign]);
        h_ptlead_pair_pt_zoomin[away_side][nsign][1]->Fill(pair_pt[nsign],m1pt[nsign],weight[nsign]);
        h_ptlead_pair_pt_log[away_side][nsign][1]->Fill(pair_pt[nsign],m1pt[nsign],weight[nsign]);
        h_minv_pair_pt[away_side][nsign][1]->Fill(pair_pt[nsign],minv[nsign],weight[nsign]);
        h_minv_pair_pt_zoomin[away_side][nsign][1]->Fill(pair_pt[nsign],minv[nsign],weight[nsign]);
        h_minv_pair_pt_log[away_side][nsign][1]->Fill(pair_pt[nsign],minv[nsign],weight[nsign]);

        // h_eta1_eta2_dphicut[2][nsign][1]->Fill(m2eta[nsign],m1eta[nsign],weight[nsign]);
        // if (abs(dphi[nsign]) < 1){
        //     h_eta1_eta2_dphicut[0][nsign][1]->Fill(m2eta[nsign],m1eta[nsign],weight[nsign]);
        // }else if(abs(dphi[nsign]) > pms.PI-1){
        //     h_eta1_eta2_dphicut[1][nsign][1]->Fill(m2eta[nsign],m1eta[nsign],weight[nsign]);
        // }
    }

    h_pair_dP_overP[away_side][nsign][0]->Fill(pair_dPoverP[nsign],weight[nsign]);
    h_pair_y[away_side][nsign][0]->Fill(pair_y[nsign],weight[nsign]);
    h_DR[away_side][nsign][0]->Fill(dr[nsign],weight[nsign]);
    h_Dphi[away_side][nsign][0]->Fill(dphi[nsign],weight[nsign]);
    h_pt_asym[away_side][nsign][0]->Fill(asym[nsign],weight[nsign]);
    h_pair_pt_ptlead_ratio[away_side][nsign][0]->Fill(pair_pt[nsign]/pt_lead[nsign],weight[nsign]);

    h_eta_avg_Dphi[away_side][nsign][0]->Fill(dphi[nsign],etaavg[nsign],weight[nsign]);
    h_Deta_Dphi[away_side][nsign][0]->Fill(dphi[nsign],deta[nsign],weight[nsign]);
    h_eta1_eta2[away_side][nsign][0]->Fill(m2eta[nsign],m1eta[nsign],weight[nsign]);
    h_pt1_pt2[away_side][nsign][0]->Fill(m2pt[nsign],m1pt[nsign],weight[nsign]);
    h_eta_avg_Deta[away_side][nsign][0]->Fill(deta[nsign],etaavg[nsign],weight[nsign]);
    h_ptlead_pair_pt[away_side][nsign][0]->Fill(pair_pt[nsign],m1pt[nsign],weight[nsign]);
    h_ptlead_pair_pt_zoomin[away_side][nsign][0]->Fill(pair_pt[nsign],m1pt[nsign],weight[nsign]);
    h_ptlead_pair_pt_log[away_side][nsign][0]->Fill(pair_pt[nsign],m1pt[nsign],weight[nsign]);
    h_minv_pair_pt[away_side][nsign][0]->Fill(pair_pt[nsign],minv[nsign],weight[nsign]);
    h_minv_pair_pt_zoomin[away_side][nsign][0]->Fill(pair_pt[nsign],minv[nsign],weight[nsign]);
    h_minv_pair_pt_log[away_side][nsign][0]->Fill(pair_pt[nsign],minv[nsign],weight[nsign]);

    // h_eta1_eta2_dphicut[2][nsign][0]->Fill(m2eta[nsign],m1eta[nsign],weight[nsign]);
    // if (abs(dphi[nsign]) < 1){
    //     h_eta1_eta2_dphicut[0][nsign][0]->Fill(m2eta[nsign],m1eta[nsign],weight[nsign]);
    // }else if(abs(dphi[nsign]) > pms.PI-1){
    //     h_eta1_eta2_dphicut[1][nsign][0]->Fill(m2eta[nsign],m1eta[nsign],weight[nsign]);
    // }
}



void MuonPairPlottingPP::WriteOutput(){

    if (mode == 1){
        outFile->cd();
        outFile->Write();
        outFile->Close();

        // outFile->cd();
        // gDirectory->Delete("h_*");

        // for (unsigned int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
        //     for (unsigned int lgapcut = 0; lgapcut < ParamsSet::nGapCuts; lgapcut++){
        //         for (unsigned int jdphi = 0; jdphi < 2; jdphi++){
        //             h_pair_dP_overP[jdphi][ksign][lgapcut]->Write();
        //             h_pair_y[jdphi][ksign][lgapcut]->Write();
        //             h_DR[jdphi][ksign][lgapcut]->Write();
        //             h_Dphi[jdphi][ksign][lgapcut]->Write();
        //             h_pt_asym[jdphi][ksign][lgapcut]->Write();
        //             h_pair_pt_ptlead_ratio[jdphi][ksign][lgapcut]->Write();

        //             h_eta_avg_Dphi[jdphi][ksign][lgapcut]->Write();
        //             h_Deta_Dphi[jdphi][ksign][lgapcut]->Write();
        //             h_eta1_eta2[jdphi][ksign][lgapcut]->Write();
        //             h_eta_avg_Deta[jdphi][ksign][lgapcut]->Write();
        //             h_pt1_pt2[jdphi][ksign][lgapcut]->Write();
        //             h_ptlead_pair_pt[jdphi][ksign][lgapcut]->Write();
        //             h_ptlead_pair_pt_zoomin[jdphi][ksign][lgapcut]->Write();
        //             h_ptlead_pair_pt_log[jdphi][ksign][lgapcut]->Write();
        //             h_minv_pair_pt[jdphi][ksign][lgapcut]->Write();
        //             h_minv_pair_pt_zoomin[jdphi][ksign][lgapcut]->Write();
        //             h_minv_pair_pt_log[jdphi][ksign][lgapcut]->Write();
        //         }

        //         // for (unsigned int idphi= 0; idphi < ParamsSet::ndphiselcs; idphi++){
        //         //     h_eta1_eta2_dphicut[idphi][ksign][lgapcut]->Write();
        //         // }
        //     }
        // }
        // outFile->Close();

    }else{ // mode = 3

    }
}

void MuonPairPlottingPP::Run(){
    clock_t start, end;
    double cpu_time_used;
    start = clock();

    if (isScram) std::cout << "Data being processed: Scrambled" << std::endl;
    else if (isTight) std::cout << "Data being processed: Tight" << std::endl;
    else std::cout << "Data being processed: Real" << std::endl;

  	InitInput();
    InitOutput();
  	InitHists();
  	ProcessData();
  	WriteOutput();

    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;;
    std::cout << "CPU time used for running mode " << mode << " is " << cpu_time_used << " seconds" << std::endl;

}