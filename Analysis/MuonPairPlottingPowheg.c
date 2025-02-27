#include "MuonPairPlottingPowheg.h"
#include "time.h"

void MuonPairPlottingPowheg::ProcessData(){
	
    for (int ibatch = 0; ibatch < nBatches; ibatch++){
        int num_bytes_meta_tree = meta_tree[ibatch]->GetEntry(0);//read in an event
        if (num_bytes_meta_tree == 0){
            std::cout<<"Error:: Read in event has size of zero bytes,  quitting"<<std::endl;
            throw std::exception();
        }
        std::cout << "nevents in batch#" << ibatch << ": " << nevents_before_cuts[ibatch] << std::endl;
        nevents_before_cuts_total += nevents_before_cuts[ibatch];
    }
    
    std::cout << "total number of events (before cuts): " << nevents_before_cuts_total << std::endl;

    for (int ibatch = 0; ibatch < nBatches; ibatch++){
        for (int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
            
    	    Long64_t nentries = inTree[ibatch][ksign]->GetEntries(); //#muon pairs
            for (Long64_t lentry=0; lentry<nentries; lentry++) {//loop over the events
            // for (Long64_t lentry=0; lentry<100; lentry++) {//loop over the events
            if(lentry%100000==0) cout<<"Processing "<<lentry<<" event out of "<<nentries<<" events"<<std::endl;
          		int num_bytes = inTree[ibatch][ksign]->GetEntry(lentry);//read in an event
                if(num_bytes==0){
                	std::cout<<"Error:: Read in event has size of zero bytes,  quitting"<<std::endl;
                	throw std::exception();
                }
                if(mode == 1){
                	FillHistograms(ibatch, ksign);
                }else{ // mode has to be 1 or 3
                	// FillPtBinnedHistograms(ibatch, jpt, ksign);
                }
            }
        }
    }
}


bool MuonPairPlottingPowheg::PassSingleMuonGapCut(float meta, float mpt, int mcharge){
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


void MuonPairPlottingPowheg::FillHistograms(int nbatch, int nsign){

    if (abs(weight[nbatch][nsign]) > crossx_cut * filter_effcy) return; // return without filling in the histograms
    double ev_weight = weight[nbatch][nsign] / nevents_before_cuts_total;
    // ngapcut = 0: all; = 1: only those that pass

    bool pass_gapcut = PassSingleMuonGapCut(m1eta[nbatch][nsign], m1pt[nbatch][nsign], m1charge[nbatch][nsign]) && PassSingleMuonGapCut(m2eta[nbatch][nsign], m2pt[nbatch][nsign], m2charge[nbatch][nsign]);
    bool ndphi = (abs(dphi[nbatch][nsign]) >= pms.PI / 2.);

    // if (pass_gapcut){
        
    //     if (nbatch == 2){ //no delta R cut
    //         h_pair_dP_overP[ndphi][nsign][1]->Fill(pair_dPoverP[nbatch][nsign],ev_weight);
    //         h_pair_y[ndphi][nsign][1]->Fill(pair_y[nbatch][nsign],ev_weight);
    //         h_DR[ndphi][nsign][1]->Fill(dr[nbatch][nsign],ev_weight);
    //         h_Dphi[ndphi][nsign][1]->Fill(dphi[nbatch][nsign],ev_weight);
    //         h_pt_asym[ndphi][nsign][1]->Fill(asym[nbatch][nsign],ev_weight);
    //         h_pair_pt_ptlead_ratio[ndphi][nsign][1]->Fill(pair_pt[nbatch][nsign]/pt_lead[nbatch][nsign],ev_weight);

    //         h_eta_avg_Dphi[ndphi][nsign][1]->Fill(dphi[nbatch][nsign],etaavg[nbatch][nsign],ev_weight);
    //         h_eta1_eta2[ndphi][nsign][1]->Fill(m2eta[nbatch][nsign],m1eta[nbatch][nsign],ev_weight);
    //         h_pt1_pt2[ndphi][nsign][1]->Fill(m2pt[nbatch][nsign],m1pt[nbatch][nsign],ev_weight);
    //         h_eta_avg_Deta[ndphi][nsign][1]->Fill(deta[nbatch][nsign],etaavg[nbatch][nsign],ev_weight);
    //         h_ptlead_pair_pt[ndphi][nsign][1]->Fill(pair_pt[nbatch][nsign],m1pt[nbatch][nsign],ev_weight);
    //         h_minv_pair_pt[ndphi][nsign][1]->Fill(pair_pt[nbatch][nsign],minv[nbatch][nsign],ev_weight);
    //         h_Deta_Dphi[ndphi][nsign][1]->Fill(dphi[nbatch][nsign],deta[nbatch][nsign],ev_weight);

    //         h_eta1_eta2_dphicut[2][nsign][1]->Fill(m2eta[nbatch][nsign],m1eta[nbatch][nsign],ev_weight);
    //         if (abs(dphi[nbatch][nsign]) < 1){
    //             h_eta1_eta2_dphicut[0][nsign][1]->Fill(m2eta[nbatch][nsign],m1eta[nbatch][nsign],ev_weight);
    //         }else if(abs(dphi[nbatch][nsign]) > pms.PI-1){
    //             h_eta1_eta2_dphicut[1][nsign][1]->Fill(m2eta[nbatch][nsign],m1eta[nbatch][nsign],ev_weight);
    //         }
    //     }
    // }

    bool is_filled = false;
    // if ((isMCTruthBB && both_from_b[nbatch][nsign]) || (isMCTruthCC && both_from_c[nbatch][nsign])){
    if (from_same_b[nbatch][nsign]){
        h_DR_ancestor_binned[nsign][nAncestorGroups]->Fill(dr[nbatch][nsign], ev_weight);
        h_pt_asym_ancestor_binned[nsign][nAncestorGroups]->Fill(asym[nbatch][nsign], ev_weight);
        h_pair_pt_ptlead_ratio_ancestor_binned[nsign][nAncestorGroups]->Fill(pair_pt[nbatch][nsign]/pt_lead[nbatch][nsign], ev_weight);
        h_ptlead_pair_pt_ancestor_binned[nsign][nAncestorGroups]->Fill(pair_pt[nbatch][nsign], pt_lead[nbatch][nsign], ev_weight);
        h_Deta_Dphi_ancestor_binned[nsign][nAncestorGroups]->Fill(dphi[nbatch][nsign], deta[nbatch][nsign], ev_weight);
        h_minv_pair_pt_ancestor_binned[nsign][nAncestorGroups]->Fill(pair_pt[nbatch][nsign], minv[nbatch][nsign], ev_weight);
        is_filled = true;

        h_DR_flavor_binned[nsign][0]->Fill(dr[nbatch][nsign], ev_weight);
        h_pt_asym_flavor_binned[nsign][0]->Fill(asym[nbatch][nsign], ev_weight);
        h_pair_pt_ptlead_ratio_flavor_binned[nsign][0]->Fill(pair_pt[nbatch][nsign]/pt_lead[nbatch][nsign], ev_weight);
        h_ptlead_pair_pt_flavor_binned[nsign][0]->Fill(pair_pt[nbatch][nsign], pt_lead[nbatch][nsign], ev_weight);
        h_Deta_Dphi_flavor_binned[nsign][0]->Fill(dphi[nbatch][nsign], deta[nbatch][nsign], ev_weight);
        h_minv_pair_pt_flavor_binned[nsign][0]->Fill(pair_pt[nbatch][nsign], minv[nbatch][nsign], ev_weight);
    }else{

        if (both_from_b[nbatch][nsign]){
            h_DR_flavor_binned[nsign][1]->Fill(dr[nbatch][nsign], ev_weight);
            h_pt_asym_flavor_binned[nsign][1]->Fill(asym[nbatch][nsign], ev_weight);
            h_pair_pt_ptlead_ratio_flavor_binned[nsign][1]->Fill(pair_pt[nbatch][nsign]/pt_lead[nbatch][nsign], ev_weight);
            h_ptlead_pair_pt_flavor_binned[nsign][1]->Fill(pair_pt[nbatch][nsign], pt_lead[nbatch][nsign], ev_weight);
            h_Deta_Dphi_flavor_binned[nsign][1]->Fill(dphi[nbatch][nsign], deta[nbatch][nsign], ev_weight);
            h_minv_pair_pt_flavor_binned[nsign][1]->Fill(pair_pt[nbatch][nsign], minv[nbatch][nsign], ev_weight);
        }else if (both_from_c[nbatch][nsign]){
            h_DR_flavor_binned[nsign][2]->Fill(dr[nbatch][nsign], ev_weight);
            h_pt_asym_flavor_binned[nsign][2]->Fill(asym[nbatch][nsign], ev_weight);
            h_pair_pt_ptlead_ratio_flavor_binned[nsign][2]->Fill(pair_pt[nbatch][nsign]/pt_lead[nbatch][nsign], ev_weight);
            h_ptlead_pair_pt_flavor_binned[nsign][2]->Fill(pair_pt[nbatch][nsign], pt_lead[nbatch][nsign], ev_weight);
            h_Deta_Dphi_flavor_binned[nsign][2]->Fill(dphi[nbatch][nsign], deta[nbatch][nsign], ev_weight);
            h_minv_pair_pt_flavor_binned[nsign][2]->Fill(pair_pt[nbatch][nsign], minv[nbatch][nsign], ev_weight);
        }else{
            h_DR_flavor_binned[nsign][3]->Fill(dr[nbatch][nsign], ev_weight);
            h_pt_asym_flavor_binned[nsign][3]->Fill(asym[nbatch][nsign], ev_weight);
            h_pair_pt_ptlead_ratio_flavor_binned[nsign][3]->Fill(pair_pt[nbatch][nsign]/pt_lead[nbatch][nsign], ev_weight);
            h_ptlead_pair_pt_flavor_binned[nsign][3]->Fill(pair_pt[nbatch][nsign], pt_lead[nbatch][nsign], ev_weight);
            h_Deta_Dphi_flavor_binned[nsign][3]->Fill(dphi[nbatch][nsign], deta[nbatch][nsign], ev_weight);
            h_minv_pair_pt_flavor_binned[nsign][3]->Fill(pair_pt[nbatch][nsign], minv[nbatch][nsign], ev_weight);
        }

        if (from_same_ancestors[nbatch][nsign]){ // from the same ancestors and not the same b
            for (int kgrp = 0; kgrp < nAncestorGroups; kgrp++){
                if (m1_ancestor_category[nbatch][nsign] == kgrp){
                    h_DR_ancestor_binned[nsign][kgrp]->Fill(dr[nbatch][nsign], ev_weight);
                    h_pt_asym_ancestor_binned[nsign][kgrp]->Fill(asym[nbatch][nsign], ev_weight);
                    h_pair_pt_ptlead_ratio_ancestor_binned[nsign][kgrp]->Fill(pair_pt[nbatch][nsign]/pt_lead[nbatch][nsign], ev_weight);
                    h_ptlead_pair_pt_ancestor_binned[nsign][kgrp]->Fill(pair_pt[nbatch][nsign], pt_lead[nbatch][nsign], ev_weight);
                    h_Deta_Dphi_ancestor_binned[nsign][kgrp]->Fill(dphi[nbatch][nsign], deta[nbatch][nsign], ev_weight);
                    h_minv_pair_pt_ancestor_binned[nsign][kgrp]->Fill(pair_pt[nbatch][nsign], minv[nbatch][nsign], ev_weight);
                    h_mQQ_ancestor_binned[nsign][kgrp]->Fill(mQQ[nbatch][nsign],ev_weight);
                    h_mQQ_Q_ratio_ancestor_binned[nsign][kgrp]->Fill(mQQ[nbatch][nsign] / Q[nbatch][nsign], ev_weight);
                    h_mQQ_mHard_ratio_ancestor_binned[nsign][kgrp]->Fill(mQQ[nbatch][nsign] / mHard_relevant[nbatch][nsign], ev_weight);
                    is_filled = true;
                    break;
                }
            }
        }
    }
    
    if (!is_filled){
        h_DR_ancestor_binned[nsign][nAncestorGroups+1]->Fill(dr[nbatch][nsign], ev_weight);
        h_pt_asym_ancestor_binned[nsign][nAncestorGroups+1]->Fill(asym[nbatch][nsign], ev_weight);
        h_pair_pt_ptlead_ratio_ancestor_binned[nsign][nAncestorGroups+1]->Fill(pair_pt[nbatch][nsign]/pt_lead[nbatch][nsign], ev_weight);
        h_ptlead_pair_pt_ancestor_binned[nsign][nAncestorGroups+1]->Fill(pair_pt[nbatch][nsign], pt_lead[nbatch][nsign], ev_weight);
        h_Deta_Dphi_ancestor_binned[nsign][nAncestorGroups+1]->Fill(dphi[nbatch][nsign], deta[nbatch][nsign], ev_weight);
        h_minv_pair_pt_ancestor_binned[nsign][nAncestorGroups+1]->Fill(pair_pt[nbatch][nsign], minv[nbatch][nsign], ev_weight);
    }

    h_pair_dP_overP[ndphi][nsign]->Fill(pair_dPoverP[nbatch][nsign],ev_weight);
    h_pair_y[ndphi][nsign]->Fill(pair_y[nbatch][nsign],ev_weight);
    h_DR[ndphi][nsign]->Fill(dr[nbatch][nsign],ev_weight);
    h_DR_zoomin[ndphi][nsign]->Fill(dr[nbatch][nsign],ev_weight);
    h_Dphi[ndphi][nsign]->Fill(dphi[nbatch][nsign],ev_weight);
    h_pt_asym[ndphi][nsign]->Fill(asym[nbatch][nsign],ev_weight);
    h_pair_pt_ptlead_ratio[ndphi][nsign]->Fill(pair_pt[nbatch][nsign]/pt_lead[nbatch][nsign],ev_weight);
    if (mQQ[nbatch][nsign] != -10 && mHard_relevant[nbatch][nsign] != -10){
        h_mQQ[ndphi][nsign]->Fill(mQQ[nbatch][nsign],ev_weight);
        h_mQQ_Q_ratio[ndphi][nsign]->Fill(mQQ[nbatch][nsign] / Q[nbatch][nsign], ev_weight);
        h_mQQ_mHard_ratio[ndphi][nsign]->Fill(mQQ[nbatch][nsign] / mHard_relevant[nbatch][nsign], ev_weight);
    }

    h_eta_avg_Dphi[ndphi][nsign]->Fill(dphi[nbatch][nsign],etaavg[nbatch][nsign],ev_weight);
    h_Deta_Dphi[ndphi][nsign]->Fill(dphi[nbatch][nsign],deta[nbatch][nsign],ev_weight);
    h_eta1_eta2[ndphi][nsign]->Fill(m2eta[nbatch][nsign],m1eta[nbatch][nsign],ev_weight);
    h_pt1_pt2[ndphi][nsign]->Fill(m2pt[nbatch][nsign],m1pt[nbatch][nsign],ev_weight);
    h_eta_avg_Deta[ndphi][nsign]->Fill(deta[nbatch][nsign],etaavg[nbatch][nsign],ev_weight);
    h_ptlead_pair_pt[ndphi][nsign]->Fill(pair_pt[nbatch][nsign],m1pt[nbatch][nsign],ev_weight);
    h_ptlead_pair_pt_zoomin[ndphi][nsign]->Fill(pair_pt[nbatch][nsign],m1pt[nbatch][nsign],ev_weight);
    h_ptlead_pair_pt_log[ndphi][nsign]->Fill(pair_pt[nbatch][nsign],m1pt[nbatch][nsign],ev_weight);
    h_minv_pair_pt[ndphi][nsign]->Fill(pair_pt[nbatch][nsign],minv[nbatch][nsign],ev_weight);
    h_minv_pair_pt_zoomin[ndphi][nsign]->Fill(pair_pt[nbatch][nsign],minv[nbatch][nsign],ev_weight);
    h_minv_pair_pt_log[ndphi][nsign]->Fill(pair_pt[nbatch][nsign],minv[nbatch][nsign],ev_weight);

    // h_eta1_eta2_dphicut[2][nsign]->Fill(m2eta[nbatch][nsign],m1eta[nbatch][nsign],ev_weight);
    // if (abs(dphi[nbatch][nsign]) < 1){
    //     h_eta1_eta2_dphicut[0][nsign]->Fill(m2eta[nbatch][nsign],m1eta[nbatch][nsign],ev_weight);
    // }else if(abs(dphi[nbatch][nsign]) > pms.PI-1){
    //     h_eta1_eta2_dphicut[1][nsign]->Fill(m2eta[nbatch][nsign],m1eta[nbatch][nsign],ev_weight);
    // }
}

void MuonPairPlottingPowheg::FillPtBinnedHistograms(int nbatch, int npt, int nsign){}


void MuonPairPlottingPowheg::WriteOutput(){
    // sub_dir = (isMCTruthBB)? "bb_full_sample/" : "cc_full_sample/";
    // outFile = new TFile(Form("%s%shistograms_%s_combined.root", mcdir.c_str(), sub_dir.c_str(), mc_mode.c_str()),"recreate");

    if (mode == 1){
        // outFile->Write();
        // outFile->cd();
        // gDirectory->Delete("h_*");
    
        for (unsigned int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
            // for (unsigned int lgapcut = 0; lgapcut < ParamsSet::nGapCuts; lgapcut++){

            for (int kgrp = 0; kgrp < nAncestorGroups + 2; kgrp++){
                h_DR_ancestor_binned[ksign][kgrp]->Write();
                h_pt_asym_ancestor_binned[ksign][kgrp]->Write();
                // h_psrapidity_ordered_pt_asym_ancestor_binned
                h_pair_pt_ptlead_ratio_ancestor_binned[ksign][kgrp]->Write();
                h_ptlead_pair_pt_ancestor_binned[ksign][kgrp]->Write();
                h_Deta_Dphi_ancestor_binned[ksign][kgrp]->Write();
                h_minv_pair_pt_ancestor_binned[ksign][kgrp]->Write();
                h_mQQ_ancestor_binned[ksign][kgrp]->Write();
                h_mQQ_Q_ratio_ancestor_binned[ksign][kgrp]->Write();
                h_mQQ_mHard_ratio_ancestor_binned[ksign][kgrp]->Write();
            }
            for (int kflav = 0; kflav < nFlavors; kflav++){
                h_DR_flavor_binned[ksign][kflav]->Write();
                h_pt_asym_flavor_binned[ksign][kflav]->Write();
                h_pair_pt_ptlead_ratio_flavor_binned[ksign][kflav]->Write();
                h_ptlead_pair_pt_flavor_binned[ksign][kflav]->Write();
                h_Deta_Dphi_flavor_binned[ksign][kflav]->Write();
                h_minv_pair_pt_flavor_binned[ksign][kflav]->Write();
                // h_mQQ_flavor_binned[ksign][kflav]->Write();
                // h_mQQ_Q_ratio_flavor_binned[ksign][kflav]->Write();
                // h_mQQ_mHard_ratio_flavor_binned[ksign][kflav]->Write();
            }

            for (unsigned int jdphi = 0; jdphi < 2; jdphi++){
                h_pair_dP_overP[jdphi][ksign]->Write();
                h_pair_y[jdphi][ksign]->Write();
                h_DR[jdphi][ksign]->Write();
                h_DR_zoomin[jdphi][ksign]->Write();
                h_Dphi[jdphi][ksign]->Write();
                h_pt_asym[jdphi][ksign]->Write();
                h_pair_pt_ptlead_ratio[jdphi][ksign]->Write();
                h_mQQ[jdphi][ksign]->Write();
                h_mQQ_Q_ratio[jdphi][ksign]->Write();
                h_mQQ_mHard_ratio[jdphi][ksign]->Write();

                h_eta_avg_Dphi[jdphi][ksign]->Write();
                h_Deta_Dphi[jdphi][ksign]->Write();
                h_eta1_eta2[jdphi][ksign]->Write();
                h_eta_avg_Deta[jdphi][ksign]->Write();
                h_pt1_pt2[jdphi][ksign]->Write();
                h_ptlead_pair_pt[jdphi][ksign]->Write();
                h_ptlead_pair_pt_zoomin[jdphi][ksign]->Write();
                h_ptlead_pair_pt_log[jdphi][ksign]->Write();
                h_minv_pair_pt[jdphi][ksign]->Write();
                h_minv_pair_pt_zoomin[jdphi][ksign]->Write();
                h_minv_pair_pt_log[jdphi][ksign]->Write();
            }

            // for (unsigned int idphi= 0; idphi < ParamsSet::ndphiselcs; idphi++){
            //     h_eta1_eta2_dphicut[idphi][ksign]->Write();
            // }
            // }
        }
    }else{ // mode = 3

    }
}

void MuonPairPlottingPowheg::Run(){
    clock_t start, end;
    double cpu_time_used;
    start = clock();

    if ((isMCTruthBB && isMCTruthCC) || (!isMCTruthBB && !isMCTruthCC)){
        std::cout << "one and only one of the booleans bb or cc needs to be true." << std::endl;
        throw std::exception();
    }
    if (isMCTruthBB) std::cout << "Data being processed: MC Truth bb" << std::endl;
    else std::cout << "Data being processed: MC Truth cc" << std::endl;

    filter_effcy = (isMCTruthBB)? filter_effcy_bb : filter_effcy_cc;

    sub_dir = (isMCTruthBB)? "bb_full_sample/" : "cc_full_sample/";
    mc_mode = (isMCTruthBB)? "mc_truth_bb" : "mc_truth_cc";

  	InitInput();
    InitOutput();
  	InitHists();
  	ProcessData();
  	WriteOutput();

    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;;
    std::cout << "CPU time used for running mode " << mode << " is " << cpu_time_used << " seconds" << std::endl;

}