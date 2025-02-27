#include "MuonPairPlottingPythia.h"
#include "time.h"

void MuonPairPlottingPythia::ProcessData(){
    for (int ifile = 0; ifile < nFiles; ifile++){
      	for (int jkin = 0; jkin < nKinRanges; jkin++){
            for (int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
                
        	    Long64_t nentries = inTree[ifile][jkin][ksign]->GetEntries(); //#muon pairs
                for (Long64_t lentry=0; lentry<nentries; lentry++) {//loop over the events
                // for (Long64_t lentry=0; lentry<100; lentry++) {//loop over the events
                if(lentry%100000==0) cout<<"Processing "<<lentry<<" event out of "<<nentries<<" events"<<std::endl;
              		int num_bytes = inTree[ifile][jkin][ksign]->GetEntry(lentry);//read in an event
                    if(num_bytes==0){
                    	std::cout<<"Error:: Read in event has size of zero bytes,  quitting"<<std::endl;
                    	throw std::exception();
                    }
                    if(mode == 1){
                    	FillHistograms(ifile, jkin, ksign);
                    }else{ // mode has to be 1 or 3
                    	// FillPtBinnedHistograms(jkin, jpt, ksign);
                    }
                }
            }
        }
    }
}


bool MuonPairPlottingPythia::PassSingleMuonGapCut(float meta, float mpt, int mcharge){
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


void MuonPairPlottingPythia::FillHistograms(int nfile, int nkin, int nsign){

    double ev_weight = weight[nfile][nkin][nsign] * nevents_before_cuts[nfile][nkin] / nevents_before_cuts_combined[nkin];

    // ngapcut = 0: all; = 1: only those that pass

    // bool pass_gapcut = PassSingleMuonGapCut(m1eta[nfile][nkin][nsign], m1pt[nfile][nkin][nsign], m1charge[nfile][nkin][nsign]) && PassSingleMuonGapCut(m2eta[nfile][nkin][nsign], m2pt[nfile][nkin][nsign], m2charge[nfile][nkin][nsign]);
    bool away_side = (abs(dphi[nfile][nkin][nsign]) >= pms.PI / 2.);

    float pT_large_eta = (abs(m1eta[nfile][nkin][nsign]) > abs(m2eta[nfile][nkin][nsign]))? m1pt[nfile][nkin][nsign] : m2pt[nfile][nkin][nsign];
    float pT_small_eta = (abs(m1eta[nfile][nkin][nsign]) > abs(m2eta[nfile][nkin][nsign]))? m2pt[nfile][nkin][nsign] : m1pt[nfile][nkin][nsign];
    float psrapidity_ordered_pt_asym = (pT_large_eta - pT_small_eta)/(pT_large_eta + pT_small_eta); // expect to be peaked towards negative value if there is asymmetry

    bool is_filled = false;
    if (from_same_b[nfile][nkin][nsign]){
        h_DR_ancestor_binned[nsign][nAncestorGroups]->Fill(dr[nfile][nkin][nsign],ev_weight);
        h_DR_zoomin_ancestor_binned[nsign][nAncestorGroups]->Fill(dr[nfile][nkin][nsign],ev_weight);
        h_DR_jacobian_corrected_ancestor_binned[nsign][nAncestorGroups]->Fill(dr[nfile][nkin][nsign],ev_weight * 1. / dr[nfile][nkin][nsign]);
        h_DR_zoomin_jacobian_corrected_ancestor_binned[nsign][nAncestorGroups]->Fill(dr[nfile][nkin][nsign],ev_weight * 1. / dr[nfile][nkin][nsign]);
        h_pt_asym_ancestor_binned[nsign][nAncestorGroups]->Fill(asym[nfile][nkin][nsign],ev_weight);
        h_psrapidity_ordered_pt_asym_ancestor_binned[nsign][nAncestorGroups]->Fill(psrapidity_ordered_pt_asym,ev_weight); 
        h_pair_pt_ptlead_ratio_ancestor_binned[nsign][nAncestorGroups]->Fill(pair_pt[nfile][nkin][nsign]/pt_lead[nfile][nkin][nsign],ev_weight);
        h_Deta_Dphi_ancestor_binned[nsign][nAncestorGroups]->Fill(dphi[nfile][nkin][nsign],deta[nfile][nkin][nsign],ev_weight);
        h_minv_pair_pt_ancestor_binned[nsign][nAncestorGroups]->Fill(pair_pt[nfile][nkin][nsign],minv[nfile][nkin][nsign],ev_weight);
        h_minv_pair_pt_zoomin_ancestor_binned[nsign][nAncestorGroups]->Fill(pair_pt[nfile][nkin][nsign],minv[nfile][nkin][nsign],ev_weight);
        h_minv_pair_pt_jacobian_corrected_ancestor_binned[nsign][nAncestorGroups]->Fill(pair_pt[nfile][nkin][nsign],minv[nfile][nkin][nsign],ev_weight * 1. / minv[nfile][nkin][nsign]);
        h_minv_pair_pt_zoomin_jacobian_corrected_ancestor_binned[nsign][nAncestorGroups]->Fill(pair_pt[nfile][nkin][nsign],minv[nfile][nkin][nsign],ev_weight * 1. / minv[nfile][nkin][nsign]);
        h_minv_pair_pt_log_ancestor_binned[nsign][nAncestorGroups]->Fill(pair_pt[nfile][nkin][nsign],minv[nfile][nkin][nsign],ev_weight);
        h_ptlead_pair_pt_ancestor_binned[nsign][nAncestorGroups]->Fill(pair_pt[nfile][nkin][nsign],pt_lead[nfile][nkin][nsign],ev_weight); 
        h_ptlead_pair_pt_zoomin_ancestor_binned[nsign][nAncestorGroups]->Fill(pair_pt[nfile][nkin][nsign],pt_lead[nfile][nkin][nsign],ev_weight); 
        h_ptlead_pair_pt_log_ancestor_binned[nsign][nAncestorGroups]->Fill(pair_pt[nfile][nkin][nsign],pt_lead[nfile][nkin][nsign],ev_weight); 
        
        h_DR_flavor_binned[nsign][0]->Fill(dr[nfile][nkin][nsign],ev_weight);
        h_DR_zoomin_flavor_binned[nsign][0]->Fill(dr[nfile][nkin][nsign],ev_weight);
        h_DR_jacobian_corrected_flavor_binned[nsign][0]->Fill(dr[nfile][nkin][nsign],ev_weight * 1. / dr[nfile][nkin][nsign]);
        h_DR_zoomin_jacobian_corrected_flavor_binned[nsign][0]->Fill(dr[nfile][nkin][nsign],ev_weight * 1. / dr[nfile][nkin][nsign]);
        h_pt_asym_flavor_binned[nsign][0]->Fill(asym[nfile][nkin][nsign],ev_weight);
        h_psrapidity_ordered_pt_asym_flavor_binned[nsign][0]->Fill(psrapidity_ordered_pt_asym,ev_weight); 
        h_pair_pt_ptlead_ratio_flavor_binned[nsign][0]->Fill(pair_pt[nfile][nkin][nsign]/pt_lead[nfile][nkin][nsign],ev_weight);
        h_Deta_Dphi_flavor_binned[nsign][0]->Fill(dphi[nfile][nkin][nsign],deta[nfile][nkin][nsign],ev_weight);
        h_minv_pair_pt_flavor_binned[nsign][0]->Fill(pair_pt[nfile][nkin][nsign],minv[nfile][nkin][nsign],ev_weight);
        h_minv_pair_pt_zoomin_flavor_binned[nsign][0]->Fill(pair_pt[nfile][nkin][nsign],minv[nfile][nkin][nsign],ev_weight);
        h_minv_pair_pt_jacobian_corrected_flavor_binned[nsign][0]->Fill(pair_pt[nfile][nkin][nsign],minv[nfile][nkin][nsign],ev_weight * 1. / minv[nfile][nkin][nsign]);
        h_minv_pair_pt_zoomin_jacobian_corrected_flavor_binned[nsign][0]->Fill(pair_pt[nfile][nkin][nsign],minv[nfile][nkin][nsign],ev_weight * 1. / minv[nfile][nkin][nsign]);
        h_minv_pair_pt_log_flavor_binned[nsign][0]->Fill(pair_pt[nfile][nkin][nsign],minv[nfile][nkin][nsign],ev_weight);
        h_ptlead_pair_pt_flavor_binned[nsign][0]->Fill(pair_pt[nfile][nkin][nsign],pt_lead[nfile][nkin][nsign],ev_weight); 
        h_ptlead_pair_pt_zoomin_flavor_binned[nsign][0]->Fill(pair_pt[nfile][nkin][nsign],pt_lead[nfile][nkin][nsign],ev_weight); 
        h_ptlead_pair_pt_log_flavor_binned[nsign][0]->Fill(pair_pt[nfile][nkin][nsign],pt_lead[nfile][nkin][nsign],ev_weight); 

        is_filled = true;
    }else{

        int flavor_ind = -1;
        if (both_from_b[nfile][nkin][nsign]){
            flavor_ind = 1;
        }else if (both_from_c[nfile][nkin][nsign]){
            flavor_ind = 2;
        }else{
            flavor_ind = 3;
        }

        if (flavor_ind > 0){
            h_DR_flavor_binned[nsign][flavor_ind]->Fill(dr[nfile][nkin][nsign],ev_weight);
            h_DR_zoomin_flavor_binned[nsign][flavor_ind]->Fill(dr[nfile][nkin][nsign],ev_weight);
            h_DR_jacobian_corrected_flavor_binned[nsign][flavor_ind]->Fill(dr[nfile][nkin][nsign],ev_weight * 1. / dr[nfile][nkin][nsign]);
            h_DR_zoomin_jacobian_corrected_flavor_binned[nsign][flavor_ind]->Fill(dr[nfile][nkin][nsign],ev_weight * 1. / dr[nfile][nkin][nsign]);
            h_pt_asym_flavor_binned[nsign][flavor_ind]->Fill(asym[nfile][nkin][nsign],ev_weight);
            h_psrapidity_ordered_pt_asym_flavor_binned[nsign][flavor_ind]->Fill(psrapidity_ordered_pt_asym,ev_weight); 
            h_pair_pt_ptlead_ratio_flavor_binned[nsign][flavor_ind]->Fill(pair_pt[nfile][nkin][nsign]/pt_lead[nfile][nkin][nsign],ev_weight);
            h_Deta_Dphi_flavor_binned[nsign][flavor_ind]->Fill(dphi[nfile][nkin][nsign],deta[nfile][nkin][nsign],ev_weight);
            h_minv_pair_pt_flavor_binned[nsign][flavor_ind]->Fill(pair_pt[nfile][nkin][nsign],minv[nfile][nkin][nsign],ev_weight);
            h_minv_pair_pt_zoomin_flavor_binned[nsign][flavor_ind]->Fill(pair_pt[nfile][nkin][nsign],minv[nfile][nkin][nsign],ev_weight);
            h_minv_pair_pt_jacobian_corrected_flavor_binned[nsign][flavor_ind]->Fill(pair_pt[nfile][nkin][nsign],minv[nfile][nkin][nsign],ev_weight * 1. / minv[nfile][nkin][nsign]);
            h_minv_pair_pt_zoomin_jacobian_corrected_flavor_binned[nsign][flavor_ind]->Fill(pair_pt[nfile][nkin][nsign],minv[nfile][nkin][nsign],ev_weight * 1. / minv[nfile][nkin][nsign]);
            h_minv_pair_pt_log_flavor_binned[nsign][flavor_ind]->Fill(pair_pt[nfile][nkin][nsign],minv[nfile][nkin][nsign],ev_weight);
            h_ptlead_pair_pt_flavor_binned[nsign][flavor_ind]->Fill(pair_pt[nfile][nkin][nsign],pt_lead[nfile][nkin][nsign],ev_weight); 
            h_ptlead_pair_pt_zoomin_flavor_binned[nsign][flavor_ind]->Fill(pair_pt[nfile][nkin][nsign],pt_lead[nfile][nkin][nsign],ev_weight); 
            h_ptlead_pair_pt_log_flavor_binned[nsign][flavor_ind]->Fill(pair_pt[nfile][nkin][nsign],pt_lead[nfile][nkin][nsign],ev_weight);             
        }

        for (int igrp = 0; igrp < nAncestorGroups; igrp++){
            if (muon_pair_origin_category[nfile][nkin][nsign] == ancestor_grps[igrp]){
                h_DR_ancestor_binned[nsign][igrp]->Fill(dr[nfile][nkin][nsign],ev_weight);
                h_DR_zoomin_ancestor_binned[nsign][igrp]->Fill(dr[nfile][nkin][nsign],ev_weight);
                h_DR_jacobian_corrected_ancestor_binned[nsign][igrp]->Fill(dr[nfile][nkin][nsign],ev_weight * 1. / dr[nfile][nkin][nsign]);
                h_DR_zoomin_jacobian_corrected_ancestor_binned[nsign][igrp]->Fill(dr[nfile][nkin][nsign],ev_weight * 1. / dr[nfile][nkin][nsign]);
                h_pt_asym_ancestor_binned[nsign][igrp]->Fill(asym[nfile][nkin][nsign],ev_weight);
                h_psrapidity_ordered_pt_asym_ancestor_binned[nsign][igrp]->Fill(psrapidity_ordered_pt_asym,ev_weight); 
                h_pair_pt_ptlead_ratio_ancestor_binned[nsign][igrp]->Fill(pair_pt[nfile][nkin][nsign]/pt_lead[nfile][nkin][nsign],ev_weight);
                h_Deta_Dphi_ancestor_binned[nsign][igrp]->Fill(dphi[nfile][nkin][nsign],deta[nfile][nkin][nsign],ev_weight);
                h_minv_pair_pt_ancestor_binned[nsign][igrp]->Fill(pair_pt[nfile][nkin][nsign],minv[nfile][nkin][nsign],ev_weight);
                h_minv_pair_pt_zoomin_ancestor_binned[nsign][igrp]->Fill(pair_pt[nfile][nkin][nsign],minv[nfile][nkin][nsign],ev_weight);
                h_minv_pair_pt_jacobian_corrected_ancestor_binned[nsign][igrp]->Fill(pair_pt[nfile][nkin][nsign],minv[nfile][nkin][nsign],ev_weight * 1. / minv[nfile][nkin][nsign]);
                h_minv_pair_pt_zoomin_jacobian_corrected_ancestor_binned[nsign][igrp]->Fill(pair_pt[nfile][nkin][nsign],minv[nfile][nkin][nsign],ev_weight * 1. / minv[nfile][nkin][nsign]);
                h_minv_pair_pt_log_ancestor_binned[nsign][igrp]->Fill(pair_pt[nfile][nkin][nsign],minv[nfile][nkin][nsign],ev_weight);
                h_ptlead_pair_pt_ancestor_binned[nsign][igrp]->Fill(pair_pt[nfile][nkin][nsign],pt_lead[nfile][nkin][nsign],ev_weight); 
                h_ptlead_pair_pt_zoomin_ancestor_binned[nsign][igrp]->Fill(pair_pt[nfile][nkin][nsign],pt_lead[nfile][nkin][nsign],ev_weight); 
                h_ptlead_pair_pt_log_ancestor_binned[nsign][igrp]->Fill(pair_pt[nfile][nkin][nsign],pt_lead[nfile][nkin][nsign],ev_weight); 
                
                if (Qsplit[nfile][nkin][nsign] != -10.){
                    h_Qsplit_ancestor_binned[nsign][igrp]->Fill(Qsplit[nfile][nkin][nsign],ev_weight);
                    h_Qsplit_pTHat_ratio_ancestor_binned[nsign][igrp]->Fill(Qsplit[nfile][nkin][nsign] / pTHat[nfile][nkin][nsign], ev_weight);
                    if (mHard_relevant[nfile][nkin][nsign] != -10.){
                        h_Qsplit_mHat_ratio_ancestor_binned[nsign][igrp]->Fill(Qsplit[nfile][nkin][nsign] / mHard_relevant[nfile][nkin][nsign], ev_weight);
                    }
                }
                is_filled = true;
            }
        }
    }

    if (!is_filled){
        h_DR_ancestor_binned[nsign][nAncestorGroups+1]->Fill(dr[nfile][nkin][nsign],ev_weight);
        h_DR_zoomin_ancestor_binned[nsign][nAncestorGroups+1]->Fill(dr[nfile][nkin][nsign],ev_weight);
        h_DR_jacobian_corrected_ancestor_binned[nsign][nAncestorGroups+1]->Fill(dr[nfile][nkin][nsign],ev_weight * 1. / dr[nfile][nkin][nsign]);
        h_DR_zoomin_jacobian_corrected_ancestor_binned[nsign][nAncestorGroups+1]->Fill(dr[nfile][nkin][nsign],ev_weight * 1. / dr[nfile][nkin][nsign]);
        h_pt_asym_ancestor_binned[nsign][nAncestorGroups+1]->Fill(asym[nfile][nkin][nsign],ev_weight);
        h_psrapidity_ordered_pt_asym_ancestor_binned[nsign][nAncestorGroups+1]->Fill(psrapidity_ordered_pt_asym,ev_weight); 
        h_pair_pt_ptlead_ratio_ancestor_binned[nsign][nAncestorGroups+1]->Fill(pair_pt[nfile][nkin][nsign]/pt_lead[nfile][nkin][nsign],ev_weight);
        h_Deta_Dphi_ancestor_binned[nsign][nAncestorGroups+1]->Fill(dphi[nfile][nkin][nsign],deta[nfile][nkin][nsign],ev_weight);
        h_minv_pair_pt_ancestor_binned[nsign][nAncestorGroups+1]->Fill(pair_pt[nfile][nkin][nsign],minv[nfile][nkin][nsign],ev_weight);
        h_minv_pair_pt_zoomin_ancestor_binned[nsign][nAncestorGroups+1]->Fill(pair_pt[nfile][nkin][nsign],minv[nfile][nkin][nsign],ev_weight);
        h_minv_pair_pt_jacobian_corrected_ancestor_binned[nsign][nAncestorGroups+1]->Fill(pair_pt[nfile][nkin][nsign],minv[nfile][nkin][nsign],ev_weight * 1. / minv[nfile][nkin][nsign]);
        h_minv_pair_pt_zoomin_jacobian_corrected_ancestor_binned[nsign][nAncestorGroups+1]->Fill(pair_pt[nfile][nkin][nsign],minv[nfile][nkin][nsign],ev_weight * 1. / minv[nfile][nkin][nsign]);
        h_minv_pair_pt_log_ancestor_binned[nsign][nAncestorGroups+1]->Fill(pair_pt[nfile][nkin][nsign],minv[nfile][nkin][nsign],ev_weight);
        h_ptlead_pair_pt_ancestor_binned[nsign][nAncestorGroups+1]->Fill(pair_pt[nfile][nkin][nsign],pt_lead[nfile][nkin][nsign],ev_weight); 
        h_ptlead_pair_pt_zoomin_ancestor_binned[nsign][nAncestorGroups+1]->Fill(pair_pt[nfile][nkin][nsign],pt_lead[nfile][nkin][nsign],ev_weight); 
        h_ptlead_pair_pt_log_ancestor_binned[nsign][nAncestorGroups+1]->Fill(pair_pt[nfile][nkin][nsign],pt_lead[nfile][nkin][nsign],ev_weight); 
    }

    //fill [0] for everyone
    h_pair_dP_overP[away_side][nsign]->Fill(pair_dPoverP[nfile][nkin][nsign],ev_weight);
    h_pair_y[away_side][nsign]->Fill(pair_y[nfile][nkin][nsign],ev_weight);
    h_DR[away_side][nsign]->Fill(dr[nfile][nkin][nsign],ev_weight);
    h_DR_zoomin[away_side][nsign]->Fill(dr[nfile][nkin][nsign],ev_weight);
    h_DR_jacobian_corrected[away_side][nsign]->Fill(dr[nfile][nkin][nsign],ev_weight * 1. / dr[nfile][nkin][nsign]);
    h_DR_zoomin_jacobian_corrected[away_side][nsign]->Fill(dr[nfile][nkin][nsign],ev_weight * 1. / dr[nfile][nkin][nsign]);
    h_Dphi[away_side][nsign]->Fill(dphi[nfile][nkin][nsign],ev_weight);
    h_pt_asym[away_side][nsign]->Fill(asym[nfile][nkin][nsign],ev_weight);
    h_pair_pt_ptlead_ratio[away_side][nsign]->Fill(pair_pt[nfile][nkin][nsign]/pt_lead[nfile][nkin][nsign],ev_weight);
    h_eta_avg_Dphi[away_side][nsign]->Fill(dphi[nfile][nkin][nsign],etaavg[nfile][nkin][nsign],ev_weight);
    h_Deta_Dphi[away_side][nsign]->Fill(dphi[nfile][nkin][nsign],deta[nfile][nkin][nsign],ev_weight);
    h_eta1_eta2[away_side][nsign]->Fill(m2eta[nfile][nkin][nsign],m1eta[nfile][nkin][nsign],ev_weight);
    h_pt1_pt2[away_side][nsign]->Fill(m2pt[nfile][nkin][nsign],m1pt[nfile][nkin][nsign],ev_weight);
    h_eta_avg_Deta[away_side][nsign]->Fill(deta[nfile][nkin][nsign],etaavg[nfile][nkin][nsign],ev_weight);
    h_ptlead_pair_pt[away_side][nsign]->Fill(pair_pt[nfile][nkin][nsign],m1pt[nfile][nkin][nsign],ev_weight);
    h_ptlead_pair_pt_zoomin[away_side][nsign]->Fill(pair_pt[nfile][nkin][nsign],m1pt[nfile][nkin][nsign],ev_weight);
    h_ptlead_pair_pt_log[away_side][nsign]->Fill(pair_pt[nfile][nkin][nsign],m1pt[nfile][nkin][nsign],ev_weight);
    h_minv_pair_pt[away_side][nsign]->Fill(pair_pt[nfile][nkin][nsign],minv[nfile][nkin][nsign],ev_weight);
    h_minv_pair_pt_zoomin[away_side][nsign]->Fill(pair_pt[nfile][nkin][nsign],minv[nfile][nkin][nsign],ev_weight);
    h_minv_pair_pt_jacobian_corrected[away_side][nsign]->Fill(pair_pt[nfile][nkin][nsign],minv[nfile][nkin][nsign],ev_weight * 1. / minv[nfile][nkin][nsign]);
    h_minv_pair_pt_zoomin_jacobian_corrected[away_side][nsign]->Fill(pair_pt[nfile][nkin][nsign],minv[nfile][nkin][nsign],ev_weight * 1. / minv[nfile][nkin][nsign]);
    h_minv_pair_pt_log[away_side][nsign]->Fill(pair_pt[nfile][nkin][nsign],minv[nfile][nkin][nsign],ev_weight);

    h_kinbin_pair_dP_overP[nkin][away_side][nsign]->Fill(pair_dPoverP[nfile][nkin][nsign],ev_weight);
    h_kinbin_pair_y[nkin][away_side][nsign]->Fill(pair_y[nfile][nkin][nsign],ev_weight);
    h_kinbin_DR[nkin][away_side][nsign]->Fill(dr[nfile][nkin][nsign],ev_weight);
    h_kinbin_Dphi[nkin][away_side][nsign]->Fill(dphi[nfile][nkin][nsign],ev_weight);
    h_kinbin_pt_asym[nkin][away_side][nsign]->Fill(asym[nfile][nkin][nsign],ev_weight);
    h_kinbin_pair_pt_ptlead_ratio[nkin][away_side][nsign]->Fill(pair_pt[nfile][nkin][nsign]/pt_lead[nfile][nkin][nsign],ev_weight);
    h_kinbin_eta_avg_Dphi[nkin][away_side][nsign]->Fill(dphi[nfile][nkin][nsign],etaavg[nfile][nkin][nsign],ev_weight);
    h_kinbin_Deta_Dphi[nkin][away_side][nsign]->Fill(dphi[nfile][nkin][nsign],deta[nfile][nkin][nsign],ev_weight);
    h_kinbin_eta1_eta2[nkin][away_side][nsign]->Fill(m2eta[nfile][nkin][nsign],m1eta[nfile][nkin][nsign],ev_weight);
    h_kinbin_pt1_pt2[nkin][away_side][nsign]->Fill(m2pt[nfile][nkin][nsign],m1pt[nfile][nkin][nsign],ev_weight);
    h_kinbin_eta_avg_Deta[nkin][away_side][nsign]->Fill(deta[nfile][nkin][nsign],etaavg[nfile][nkin][nsign],ev_weight);
    h_kinbin_ptlead_pair_pt[nkin][away_side][nsign]->Fill(pair_pt[nfile][nkin][nsign],m1pt[nfile][nkin][nsign],ev_weight);
    h_kinbin_minv_pair_pt[nkin][away_side][nsign]->Fill(pair_pt[nfile][nkin][nsign],minv[nfile][nkin][nsign],ev_weight);
    
    if (abs(dphi[nfile][nkin][nsign]) < 1){
        h_eta1_eta2_dphicut[0][nsign]->Fill(m2eta[nfile][nkin][nsign],m1eta[nfile][nkin][nsign],ev_weight);
    }else if(abs(dphi[nfile][nkin][nsign]) > pms.PI-1){
        h_eta1_eta2_dphicut[1][nsign]->Fill(m2eta[nfile][nkin][nsign],m1eta[nfile][nkin][nsign],ev_weight);
    }
}

void MuonPairPlottingPythia::FillPtBinnedHistograms(int nkin, int npt, int nsign){}


void MuonPairPlottingPythia::WriteOutput(){
    // outFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/pythia/histograms_pythia_combined.root","recreate");

    if (mode == 1){
        // outFile->Write();
        // outFile->cd();
        // gDirectory->Delete("h_*");
    

        for (unsigned int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
            for (int jgrp = 0; jgrp < nAncestorGroups + 2; jgrp++){
                h_DR_ancestor_binned[ksign][jgrp]->Write();
                h_DR_zoomin_ancestor_binned[ksign][jgrp]->Write();
                h_DR_jacobian_corrected_ancestor_binned[ksign][jgrp]->Write();
                h_DR_zoomin_jacobian_corrected_ancestor_binned[ksign][jgrp]->Write();
                h_pt_asym_ancestor_binned[ksign][jgrp]->Write();
                h_psrapidity_ordered_pt_asym_ancestor_binned[ksign][jgrp]->Write();
                h_pair_pt_ptlead_ratio_ancestor_binned[ksign][jgrp]->Write();
                h_Deta_Dphi_ancestor_binned[ksign][jgrp]->Write();
                h_minv_pair_pt_ancestor_binned[ksign][jgrp]->Write();
                h_minv_pair_pt_zoomin_ancestor_binned[ksign][jgrp]->Write();
                h_minv_pair_pt_jacobian_corrected_ancestor_binned[ksign][jgrp]->Write();
                h_minv_pair_pt_zoomin_jacobian_corrected_ancestor_binned[ksign][jgrp]->Write();
                h_minv_pair_pt_log_ancestor_binned[ksign][jgrp]->Write();
                h_ptlead_pair_pt_ancestor_binned[ksign][jgrp]->Write();
                h_ptlead_pair_pt_zoomin_ancestor_binned[ksign][jgrp]->Write();
                h_ptlead_pair_pt_log_ancestor_binned[ksign][jgrp]->Write();
                h_Qsplit_ancestor_binned[ksign][jgrp]->Write();
                h_Qsplit_pTHat_ratio_ancestor_binned[ksign][jgrp]->Write();
                h_Qsplit_mHat_ratio_ancestor_binned[ksign][jgrp]->Write();
            }

            for (int jflav = 0; jflav < nFlavors; jflav++){
                h_DR_flavor_binned[ksign][jflav]->Write();
                h_DR_zoomin_flavor_binned[ksign][jflav]->Write();
                h_DR_jacobian_corrected_flavor_binned[ksign][jflav]->Write();
                h_DR_zoomin_jacobian_corrected_flavor_binned[ksign][jflav]->Write();
                h_pt_asym_flavor_binned[ksign][jflav]->Write();
                h_psrapidity_ordered_pt_asym_flavor_binned[ksign][jflav]->Write();
                h_pair_pt_ptlead_ratio_flavor_binned[ksign][jflav]->Write();
                h_Deta_Dphi_flavor_binned[ksign][jflav]->Write();
                h_minv_pair_pt_flavor_binned[ksign][jflav]->Write();
                h_minv_pair_pt_zoomin_flavor_binned[ksign][jflav]->Write();
                h_minv_pair_pt_jacobian_corrected_flavor_binned[ksign][jflav]->Write();
                h_minv_pair_pt_zoomin_jacobian_corrected_flavor_binned[ksign][jflav]->Write();
                h_minv_pair_pt_log_flavor_binned[ksign][jflav]->Write();
                h_ptlead_pair_pt_flavor_binned[ksign][jflav]->Write();
                h_ptlead_pair_pt_zoomin_flavor_binned[ksign][jflav]->Write();
                h_ptlead_pair_pt_log_flavor_binned[ksign][jflav]->Write();
                // h_Qsplit_flavor_binned[ksign][jflav]->Write();
                // h_Qsplit_pTHat_ratio_flavor_binned[ksign][jflav]->Write();
                // h_Qsplit_mHat_ratio_flavor_binned[ksign][jflav]->Write();
            }

            for (unsigned int jdphi = 0; jdphi < 2; jdphi++){
                h_pair_dP_overP[jdphi][ksign]->Write();
                h_pair_y[jdphi][ksign]->Write();
                h_DR[jdphi][ksign]->Write();
                h_DR_zoomin[jdphi][ksign]->Write();
                h_DR_jacobian_corrected[jdphi][ksign]->Write();
                h_DR_zoomin_jacobian_corrected[jdphi][ksign]->Write();
                h_Dphi[jdphi][ksign]->Write();
                h_pt_asym[jdphi][ksign]->Write();
                h_pair_pt_ptlead_ratio[jdphi][ksign]->Write();
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
                h_minv_pair_pt_jacobian_corrected[jdphi][ksign]->Write();
                h_minv_pair_pt_zoomin_jacobian_corrected[jdphi][ksign]->Write();
                h_minv_pair_pt_log[jdphi][ksign]->Write();

                for (int ikin = 0; ikin < nKinRanges; ikin++){
                    h_kinbin_pair_dP_overP[ikin][jdphi][ksign]->Write();
                    h_kinbin_pair_y[ikin][jdphi][ksign]->Write();
                    h_kinbin_DR[ikin][jdphi][ksign]->Write();
                    h_kinbin_Dphi[ikin][jdphi][ksign]->Write();
                    h_kinbin_pt_asym[ikin][jdphi][ksign]->Write();
                    h_kinbin_pair_pt_ptlead_ratio[ikin][jdphi][ksign]->Write();
                    h_kinbin_eta_avg_Dphi[ikin][jdphi][ksign]->Write();
                    h_kinbin_Deta_Dphi[ikin][jdphi][ksign]->Write();
                    h_kinbin_eta1_eta2[ikin][jdphi][ksign]->Write();
                    h_kinbin_eta_avg_Deta[ikin][jdphi][ksign]->Write();
                    h_kinbin_pt1_pt2[ikin][jdphi][ksign]->Write();
                    h_kinbin_ptlead_pair_pt[ikin][jdphi][ksign]->Write();
                    h_kinbin_minv_pair_pt[ikin][jdphi][ksign]->Write();
                }
            }
            for (unsigned int idphi= 0; idphi < ParamsSet::ndphiselcs; idphi++){
                h_eta1_eta2_dphicut[idphi][ksign]->Write();
            }
        }
    }
}

void MuonPairPlottingPythia::Run(){
    clock_t start, end;
    double cpu_time_used;
    start = clock();

  	InitInput();
    InitOutput();
  	InitHists();
  	ProcessData();
  	WriteOutput();

    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;;
    std::cout << "CPU time used for running mode " << mode << " is " << cpu_time_used << " seconds" << std::endl;

}