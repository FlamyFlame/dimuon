#include "MuonPairPlotting.h"
#include "time.h"
#include "EfficiencyCorrs/TrigAndRecoEff.C"

void MuonPairPlotting::ProcessData(){
	
    int cur_ctr_intvl = 0;
    for (int jctr = 0; jctr < nCtrBins; jctr++){
        for (unsigned int jj = 0; jj < ctrBins[jctr].size(); jj++){
            for (int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
                
        	    Long64_t nentries = inTree[cur_ctr_intvl][ksign]->GetEntries(); //#muon pairs
                for (Long64_t lentry=0; lentry<nentries; lentry++) {//loop over the events
                // for (Long64_t lentry=0; lentry<100; lentry++) {//loop over the events
                if(lentry%100000==0) cout<<"Processing "<<lentry<<" event out of "<<nentries<<" events"<<std::endl;

              		int num_bytes = inTree[cur_ctr_intvl][ksign]->GetEntry(lentry);//read in an event
                    if(num_bytes==0){
                    	std::cout<<"Error:: Read in event has size of zero bytes,  quitting"<<std::endl;
                    	throw std::exception();
                    }
                    FillHistograms(jctr, cur_ctr_intvl, ksign);
                }
            }
            cur_ctr_intvl++;
        }
    }
    assert(cur_ctr_intvl == 16);
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


void MuonPairPlotting::FillHistograms(int nctr_bin, int nctr_intvl, int nsign){
    // 

    // ngapcut = 0: all; = 1: only those that pass

    //if pass eta gap cut
    bool pass_gapcut = PassSingleMuonGapCut(m1eta[nctr_intvl][nsign], m1pt[nctr_intvl][nsign], m1charge[nctr_intvl][nsign]) && PassSingleMuonGapCut(m2eta[nctr_intvl][nsign], m2pt[nctr_intvl][nsign], m2charge[nctr_intvl][nsign]);
    if (pass_gapcut){
        h_pair_dP_overP[nctr_bin][nsign][1]->Fill(pair_dPoverP[nctr_intvl][nsign]);
        h_pair_y[nctr_bin][nsign][1]->Fill(pair_y[nctr_intvl][nsign]);
        h_DR[nctr_bin][nsign][1]->Fill(dr[nctr_intvl][nsign]);
        h_pt_asym[nctr_bin][nsign][1]->Fill(asym[nctr_intvl][nsign]);
        h_pair_pt_ptlead_ratio[nctr_bin][nsign][1]->Fill(pair_pt[nctr_intvl][nsign]/m1pt[nctr_intvl][nsign]);
        h_eta_avg_Dphi[nctr_bin][nsign][1]->Fill(dphi[nctr_intvl][nsign],etaavg[nctr_intvl][nsign]);
        h_Deta_Dphi[nctr_bin][nsign][1]->Fill(dphi[nctr_intvl][nsign],deta[nctr_intvl][nsign]);
        h_eta1_eta2[nctr_bin][nsign][1]->Fill(m2eta[nctr_intvl][nsign],m1eta[nctr_intvl][nsign]);
        h_pt1_pt2[nctr_bin][nsign][1]->Fill(m2pt[nctr_intvl][nsign],m1pt[nctr_intvl][nsign]);
        h_eta_avg_Deta[nctr_bin][nsign][1]->Fill(deta[nctr_intvl][nsign],etaavg[nctr_intvl][nsign]);
        h_ptlead_pair_pt[nctr_bin][nsign][1]->Fill(pair_pt[nctr_intvl][nsign],m1pt[nctr_intvl][nsign]);
        h_minv_pair_pt[nctr_bin][nsign][1]->Fill(pair_pt[nctr_intvl][nsign],minv[nctr_intvl][nsign]);
        h_minv_pair_pt_zoomin[nctr_bin][nsign][1]->Fill(pair_pt[nctr_intvl][nsign],minv[nctr_intvl][nsign]);
        h_minv_pair_pt_log[nctr_bin][nsign][1]->Fill(pair_pt[nctr_intvl][nsign],minv[nctr_intvl][nsign]);
        h_eta_avg_pair_eta[nctr_bin][nsign][1]->Fill(pair_eta[nctr_intvl][nsign],etaavg[nctr_intvl][nsign]);

        // if (ndr == 2){ //no delta R cut
        //     h_eta1_eta2_dphicut[2][nctr_bin][nsign][1]->Fill(m2eta[nctr_intvl][nsign],m1eta[nctr_intvl][nsign]);
        //     if (dphi[nctr_intvl][nsign] < 1){
        //         h_eta1_eta2_dphicut[0][nctr_bin][nsign][1]->Fill(m2eta[nctr_intvl][nsign],m1eta[nctr_intvl][nsign]);
        //     }else if(dphi[nctr_intvl][nsign] > pms.PI-1){
        //         h_eta1_eta2_dphicut[1][nctr_bin][nsign][1]->Fill(m2eta[nctr_intvl][nsign],m1eta[nctr_intvl][nsign]);
        //     }
        // }
    }
    // else{
        // h_pair_eta_Dphi_failgapcut[nctr_bin][nsign]->Fill(dphi[nctr_intvl][nsign],pair_eta[nctr_intvl][nsign]);
        // h_Dphi_failgapcut[nctr_bin][nsign]->Fill(dphi[nctr_intvl][nsign]);
    // }

    //fill [1] for everyone
    h_pair_dP_overP[nctr_bin][nsign][0]->Fill(pair_dPoverP[nctr_intvl][nsign]);
    h_pair_y[nctr_bin][nsign][0]->Fill(pair_y[nctr_intvl][nsign]);
    h_DR[nctr_bin][nsign][0]->Fill(dr[nctr_intvl][nsign]);
    h_pt_asym[nctr_bin][nsign][0]->Fill(asym[nctr_intvl][nsign]);
    h_pair_pt_ptlead_ratio[nctr_bin][nsign][0]->Fill(pair_pt[nctr_intvl][nsign]/m1pt[nctr_intvl][nsign]);

    h_Deta_Dphi[nctr_bin][nsign][0]->Fill(dphi[nctr_intvl][nsign],deta[nctr_intvl][nsign]);
    h_eta_avg_Dphi[nctr_bin][nsign][0]->Fill(dphi[nctr_intvl][nsign],etaavg[nctr_intvl][nsign]);
    h_eta1_eta2[nctr_bin][nsign][0]->Fill(m2eta[nctr_intvl][nsign],m1eta[nctr_intvl][nsign]);
    h_pt1_pt2[nctr_bin][nsign][0]->Fill(m2pt[nctr_intvl][nsign],m1pt[nctr_intvl][nsign]);
    h_eta_avg_Deta[nctr_bin][nsign][0]->Fill(deta[nctr_intvl][nsign],etaavg[nctr_intvl][nsign]);
    h_eta_avg_pair_eta[nctr_bin][nsign][0]->Fill(pair_eta[nctr_intvl][nsign],etaavg[nctr_intvl][nsign]);
    h_ptlead_pair_pt[nctr_bin][nsign][0]->Fill(pair_pt[nctr_intvl][nsign],m1pt[nctr_intvl][nsign]);
    h_minv_pair_pt[nctr_bin][nsign][0]->Fill(pair_pt[nctr_intvl][nsign],minv[nctr_intvl][nsign]);

    // if (ndr == 2){ //no delta R cut
    //     h_eta1_eta2_dphicut[2][nctr_bin][nsign][0]->Fill(m2eta[nctr_intvl][nsign],m1eta[nctr_intvl][nsign]);
    //     if (dphi[nctr_intvl][nsign] < 1){
    //         h_eta1_eta2_dphicut[0][nctr_bin][nsign][0]->Fill(m2eta[nctr_intvl][nsign],m1eta[nctr_intvl][nsign]);
    //     }else if(dphi[nctr_intvl][nsign] > pms.PI-1){
    //         h_eta1_eta2_dphicut[1][nctr_bin][nsign][0]->Fill(m2eta[nctr_intvl][nsign],m1eta[nctr_intvl][nsign]);
    //     }
    // }
}


void MuonPairPlotting::WriteOutput(){
    if (isScram){
        outFile = new TFile(("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/histograms_scrambled_pairs" + postfix + ".root").c_str(),"recreate");
    }else{
        outFile = new TFile(("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/histograms_real_pairs" + postfix + ".root").c_str(),"recreate");
    }

    outFile->Write();

    // outFile->cd();
    // if (gDirectory->Get("h_*") != nullptr) gDirectory->Delete("h_*");

    // if (outFile->GetDirectory("ctr-binned") == nullptr){}
    // if (gDirectory->GetDirectory("ctr-binned") == nullptr){
    //   gDirectory->mkdir("ctr-binned");
    // }
    // gDirectory->cd("ctr-binned");
    // gDirectory->Delete("h_*");
    
    // for (unsigned int ipt = 0; ipt < ParamsSet::nPtBins; ipt++){
    //     for (unsigned int ictr = 0; ictr < ParamsSet::nCtrBins; ictr++){
    //         for (unsigned int isign = 0; isign < ParamsSet::nSigns; isign++){
    //             for (unsigned int idphi = 0; idphi < ParamsSet::ndphiRegions; idphi++){
    //                 for (unsigned int ideta = 0; ideta < ParamsSet::ndetaRegions; ideta++){
    //                     for (unsigned int igap = 0; igap < ParamsSet::nGapCuts; igap++){
    //                         h_pair_dP_overP[ipt][ictr][isign][idphi][ideta][igap]->Write();
    //                         h_pair_y[ipt][ictr][isign][idphi][ideta][igap]->Write();
    //                         h_DR[ipt][ictr][isign][idphi][ideta][igap]->Write();
    //                         h_pt_asym[ipt][ictr][isign][idphi][ideta][igap]->Write();
    //                         h_pair_pt_ptlead_ratio[ipt][ictr][isign][idphi][ideta][igap]->Write();
    //                         h_Deta_Dphi[ipt][ictr][isign][idphi][ideta][igap]->Write();
    //                         h_eta_avg_Dphi[ipt][ictr][isign][idphi][ideta][igap]->Write();
    //                         h_eta1_eta2[ipt][ictr][isign][idphi][ideta][igap]->Write();
    //                         h_eta_avg_Deta[ipt][ictr][isign][idphi][ideta][igap]->Write();
    //                         h_eta_avg_pair_eta[ipt][ictr][isign][idphi][ideta][igap]->Write();
    //                         h_pt1_pt2[ipt][ictr][isign][idphi][ideta][igap]->Write();
    //                         h_ptlead_pair_pt[ipt][ictr][isign][idphi][ideta][igap]->Write();
    //                         h_minv_pair_pt[ipt][ictr][isign][idphi][ideta][igap]->Write();
    //                     }
    //                 }
    //             }
    //         }
    //     }            
    // }

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