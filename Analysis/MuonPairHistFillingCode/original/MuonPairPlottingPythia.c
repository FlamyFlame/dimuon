#include "MuonPairPlottingPythia.h"
#include "time.h"
#include "TLorentzVector.h"

MuonPairPlottingPythia::MuonPairPlottingPythia(){
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

    with_data_resonance_cuts_suffix = (with_data_resonance_cuts)? "_with_data_resonance_cuts" : "_no_data_resonance_cuts";
    for (std::string& file_name : file_names){
      file_name += with_data_resonance_cuts_suffix + ".root";
    } 

    for (int ifile = 0; ifile < nFiles; ifile++){
        inFile[ifile] = new TFile(file_names[ifile].c_str(),"read");
        if (!inFile[ifile]){
            std::cout << "File with the name " << file_names[ifile] << "does not exist or cannot be opened";
            throw std::exception();
        }


        for (int jkin = 0; jkin < nKinRanges; jkin++){
            for (int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
                std::string tree_name = "muon_pair_tree_kin" + std::to_string(jkin) + "_sign" + std::to_string(ksign+1);
                inTree[ifile][jkin][ksign] = (TTree*) inFile[ifile]->Get(Form("muon_pair_tree_kin%d_sign%d",jkin,ksign+1));
                if (!inTree[ifile][jkin][ksign]){
                    std::cout << "File with the name " << file_names[ifile] << "does not have a tree named" << tree_name << std::endl;
                    throw std::exception();
                }
                inTree[ifile][jkin][ksign]->SetBranchAddress("weight"          , &weight[ifile][jkin][ksign]);
                inTree[ifile][jkin][ksign]->SetBranchAddress("pTHat"          , &pTHat[ifile][jkin][ksign]);
                inTree[ifile][jkin][ksign]->SetBranchAddress("Qsplit"          , &Qsplit[ifile][jkin][ksign]);
                inTree[ifile][jkin][ksign]->SetBranchAddress("mHard_relevant"          , &mHard_relevant[ifile][jkin][ksign]);
                inTree[ifile][jkin][ksign]->SetBranchAddress("from_same_resonance"          , &from_same_resonance[ifile][jkin][ksign]);
                inTree[ifile][jkin][ksign]->SetBranchAddress("resonance_contaminated"          , &resonance_contaminated[ifile][jkin][ksign]);
                inTree[ifile][jkin][ksign]->SetBranchAddress("from_same_b"          , &from_same_b[ifile][jkin][ksign]);
                inTree[ifile][jkin][ksign]->SetBranchAddress("both_from_b"          , &both_from_b[ifile][jkin][ksign]);
                inTree[ifile][jkin][ksign]->SetBranchAddress("both_from_c"          , &both_from_c[ifile][jkin][ksign]);
                inTree[ifile][jkin][ksign]->SetBranchAddress("data_resonance_or_reso_contam_tagged_old"          , &data_resonance_or_reso_contam_tagged_old[ifile][jkin][ksign]);
                inTree[ifile][jkin][ksign]->SetBranchAddress("data_resonance_or_reso_contam_tagged_new"          , &data_resonance_or_reso_contam_tagged_new[ifile][jkin][ksign]);
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
                inTree[ifile][jkin][ksign]->SetBranchAddress("m1.phi"       , &m1phi[ifile][jkin][ksign]);
                inTree[ifile][jkin][ksign]->SetBranchAddress("m2.phi"       , &m2phi[ifile][jkin][ksign]);
                // inTree[ifile][jkin][ksign]->SetBranchAddress("m1.charge"           , &m1charge[ifile][jkin][ksign]);
                // inTree[ifile][jkin][ksign]->SetBranchAddress("m2.charge"           , &m2charge[ifile][jkin][ksign]);
        }
    }
}

void MuonPairPlottingPythia::InitOutput(){
    output_file_path = "/usatlas/u/yuhanguo/usatlasdata/pythia/histograms_pythia_combined" + with_data_resonance_cuts_suffix + ".root";
    outFile = new TFile(output_file_path.c_str(),"recreate");
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

    flavor_grp_map[pair_flavor_index::from_resonance] = "_resonance";
    flavor_grp_map[pair_flavor_index::resonance_contaminated] = "_resonance_contaminated";
    flavor_grp_map[pair_flavor_index::from_single_b] = "_single_b";
    flavor_grp_map[pair_flavor_index::bb] = "_bb";
    flavor_grp_map[pair_flavor_index::cc] = "_cc";
    flavor_grp_map[pair_flavor_index::other_flavor] = "_other_flavors";

    // h_pt_asym[jflavor][ksign] = new TH1D(Form("h_pt_asym_%s_sign%d",flavor_labels[jflavor].c_str(),ksign+1),";A = (pT1 - pT2)/(pT1 + pT2);d#sigma/dA", npt_asym_bins,0,1.);
    // h_psrapidity_ordered_pt_asym[jflavor][ksign] = new TH1D(Form("h_psrapidity_ordered_pt_asym_%s_sign%d",flavor_labels[jflavor].c_str(),ksign+1),";A' = (pT_{large |#eta|} - pT_{small |#eta|})/(pT_{large |#eta|} + pT_{small |#eta|});d#sigma/dA'", npt_asym_bins*2,-1.,1.);
    // h_pair_pt_ptlead_ratio[jflavor][ksign] = new TH1D(Form("h_pair_pt_ptlead_ratio_%s_sign%d",flavor_labels[jflavor].c_str(),ksign+1),";#frac{p_{T}^{pair}}{p_{T}^{lead}};d#sigma/d#frac{p_{T}^{pair}}{p_{T}^{lead}}", npair_pt_ptlead_ratio_bins,0,2.);
    // h_Deta_Dphi[jflavor][ksign] = new TH2D(Form("h_Deta_Dphi_%s_sign%d",flavor_labels[jflavor].c_str(),ksign+1),";#Delta#phi;#Delta#eta", nDphi_bins,-pms.PI,pms.PI,nDeta_bins,-4.8,4.8);
    // h_ptlead_pair_pt[jflavor][ksign] = new TH2D(Form("h_ptlead_pair_pt_%s_sign%d",flavor_labels[jflavor].c_str(),ksign+1),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",pms.npairPT_bins,pms.pairPTBins[ksign][2],pms.npt_bins,pms.pTBins);
    // h_minv_pair_pt[jflavor][ksign] = new TH2D(Form("h_minv_pair_pt_%s_sign%d",flavor_labels[jflavor].c_str(),ksign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",pms.npairPT_bins,pms.pairPTBins[ksign][2],nminv_bins,minv_bins[ksign]);
      

    // ------------------------------------------ cross sections ------------------------------------------
    h_crossx_truth_from_single_b_vs_pair_pt_pair_eta = new TH2D("h_crossx_truth_from_single_b_vs_pair_pt_pair_eta",";#eta^{pair};p_{T}^{pair} [GeV];#sigma^{truth}",npair_eta_bins_coarse,pair_eta_min,pair_eta_max,npair_pT_bins_coarse,pair_pt_min,pair_pt_max);
    h_crossx_truth_from_single_b_vs_pair_pt_pair_eta->Sumw2();

    // ------------------------------------------ kinematic distributions ------------------------------------------
    if (mode == 1){
        for (unsigned int ksign = 0; ksign < ParamsSet::nSigns; ksign++){

            for (int kflav = 0; kflav < pair_flavor_index::nFlavors; kflav++){
                h_DR_flavor_binned[ksign][kflav] = new TH1D(Form("h_DR_sign%d%s",ksign+1,flavor_grp_map[kflav].c_str()),";#Delta R;1/N_{evt} dN/d#Delta R", nDR_bins,0,pms.deltaR_thrsh[2]);
                h_DR_zoomin_flavor_binned[ksign][kflav] = new TH1D(Form("h_DR_zoomin_sign%d%s",ksign+1,flavor_grp_map[kflav].c_str()),";#Delta R;1/N_{evt} dN/d#Delta R", nDR_zoomin_bins,0,pms.deltaR_thrsh_zoomin);
                h_DR_jacobian_corrected_flavor_binned[ksign][kflav] = new TH1D(Form("h_DR_jacobian_corrected_sign%d%s",ksign+1,flavor_grp_map[kflav].c_str()),";#Delta R;#frac{1}{N_{evt}} #frac{1}{#Delta R} #frac{dN}{d#Delta R}", nDR_bins,0,pms.deltaR_thrsh[2]);
                h_DR_zoomin_jacobian_corrected_flavor_binned[ksign][kflav] = new TH1D(Form("h_DR_zoomin_jacobian_corrected_sign%d%s",ksign+1,flavor_grp_map[kflav].c_str()),";#Delta R;#frac{1}{N_{evt}} #frac{1}{#Delta R} #frac{dN}{d#Delta R}", nDR_zoomin_bins,0,pms.deltaR_thrsh_zoomin);
                h_pt_asym_flavor_binned[ksign][kflav] = new TH1D(Form("h_pt_asym_sign%d%s",ksign+1,flavor_grp_map[kflav].c_str()),";A = (pT1 - pT2)/(pT1 + pT2);d#sigma/dA", npt_asym_bins,0,1.);
                h_psrapidity_ordered_pt_asym_flavor_binned[ksign][kflav] = new TH1D(Form("h_psrapidity_ordered_pt_asym_sign%d%s",ksign+1,flavor_grp_map[kflav].c_str()),";A' = (pT_{large |#eta|} - pT_{small |#eta|})/(pT_{large |#eta|} + pT_{small |#eta|});d#sigma/dA'", npt_asym_bins*2,-1.,1.);
                h_pair_pt_ptlead_ratio_flavor_binned[ksign][kflav] = new TH1D(Form("h_pair_pt_ptlead_ratio_sign%d%s",ksign+1,flavor_grp_map[kflav].c_str()),";#frac{p_{T}^{pair}}{p_{T}^{lead}};d#sigma/d#frac{p_{T}^{pair}}{p_{T}^{lead}}", npair_pt_ptlead_ratio_bins,0,2.);
                h_pair_pt_jacobian_corrected_flavor_binned[ksign][kflav] = new TH1D(Form("h_pair_pt_jacobian_corrected_sign%d%s",ksign+1,flavor_grp_map[kflav].c_str()),";p_{T}^{pair} [GeV];#frac{1}{N_{evt}} #frac{1}{p_{T}^{pair}} #frac{dN}{dp_{T}^{pair}}",npair_pT_bins_linear,0,30);
                
                h_Deta_Dphi_flavor_binned[ksign][kflav] = new TH2D(Form("h_Deta_Dphi_sign%d%s",ksign+1,flavor_grp_map[kflav].c_str()),";#Delta#phi;#Delta#eta", nDphi_bins,-pms.PI,pms.PI,nDeta_bins,-4.8,4.8);
                h_minv_pair_pt_flavor_binned[ksign][kflav] = new TH2D(Form("h_minv_pair_pt_sign%d%s",ksign+1,flavor_grp_map[kflav].c_str()),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",npair_pT_bins_linear,0,30,nminv_bins_linear,0,30);
                h_minv_pair_pt_zoomin_flavor_binned[ksign][kflav] = new TH2D(Form("h_minv_pair_pt_zoomin_sign%d%s",ksign+1,flavor_grp_map[kflav].c_str()),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",npair_pT_bins_linear,0,20,nminv_bins_linear,0,4);
                h_minv_pair_pt_jacobian_corrected_flavor_binned[ksign][kflav] = new TH2D(Form("h_minv_pair_pt_jacobian_corrected_sign%d%s",ksign+1,flavor_grp_map[kflav].c_str()),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",npair_pT_bins_linear,0,30,nminv_bins_linear,0,30);
                h_minv_pair_pt_zoomin_jacobian_corrected_flavor_binned[ksign][kflav] = new TH2D(Form("h_minv_pair_pt_zoomin_jacobian_corrected_sign%d%s",ksign+1,flavor_grp_map[kflav].c_str()),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",npair_pT_bins_linear,0,20,nminv_bins_linear,0,4);
                h_minv_pair_pt_log_flavor_binned[ksign][kflav] = new TH2D(Form("h_minv_pair_pt_log_sign%d%s",ksign+1,flavor_grp_map[kflav].c_str()),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",pms.npairPT_bins,pms.pairPTBins[0][2],nminv_bins_log,minv_bins_log[0]);
                h_ptlead_pair_pt_flavor_binned[ksign][kflav] = new TH2D(Form("h_ptlead_pair_pt_sign%d%s",ksign+1,flavor_grp_map[kflav].c_str()),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",npair_pT_bins_linear,0,30,npT_lead_bins_linear,0,30);
                h_ptlead_pair_pt_zoomin_flavor_binned[ksign][kflav] = new TH2D(Form("h_ptlead_pair_pt_zoomin_sign%d%s",ksign+1,flavor_grp_map[kflav].c_str()),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",npair_pT_bins_linear,0,20,npT_lead_bins_linear,0,15);
                h_ptlead_pair_pt_log_flavor_binned[ksign][kflav] = new TH2D(Form("h_ptlead_pair_pt_log_sign%d%s",ksign+1,flavor_grp_map[kflav].c_str()),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",pms.npairPT_bins,pms.pairPTBins[0][2],pms.npt_bins,pms.pTBins);
                
                h_DR_flavor_binned[ksign][kflav]->Sumw2();
                h_DR_zoomin_flavor_binned[ksign][kflav]->Sumw2();
                h_DR_jacobian_corrected_flavor_binned[ksign][kflav]->Sumw2();
                h_DR_zoomin_jacobian_corrected_flavor_binned[ksign][kflav]->Sumw2();
                h_pt_asym_flavor_binned[ksign][kflav]->Sumw2();
                h_psrapidity_ordered_pt_asym_flavor_binned[ksign][kflav]->Sumw2();
                h_pair_pt_ptlead_ratio_flavor_binned[ksign][kflav]->Sumw2();
                h_pair_pt_jacobian_corrected_flavor_binned[ksign][kflav]->Sumw2();
                h_Deta_Dphi_flavor_binned[ksign][kflav]->Sumw2();
                h_minv_pair_pt_flavor_binned[ksign][kflav]->Sumw2();
                h_minv_pair_pt_zoomin_flavor_binned[ksign][kflav]->Sumw2();
                h_minv_pair_pt_jacobian_corrected_flavor_binned[ksign][kflav]->Sumw2();
                h_pair_pt_jacobian_corrected_flavor_binned[ksign][kflav]->Sumw2();
                h_minv_pair_pt_zoomin_jacobian_corrected_flavor_binned[ksign][kflav]->Sumw2();
                h_minv_pair_pt_log_flavor_binned[ksign][kflav]->Sumw2();
                h_ptlead_pair_pt_flavor_binned[ksign][kflav]->Sumw2();
                h_ptlead_pair_pt_zoomin_flavor_binned[ksign][kflav]->Sumw2();
                h_ptlead_pair_pt_log_flavor_binned[ksign][kflav]->Sumw2();
            }


            for (int kgrp = 0; kgrp < nAncestorGroupsTotal; kgrp++){
                h_DR_ancestor_binned[ksign][kgrp] = new TH1D(Form("h_DR_sign%d%s",ksign+1,ancestor_grp_labels[kgrp].c_str()),";#Delta R;1/N_{evt} dN/d#Delta R", nDR_bins,0,pms.deltaR_thrsh[2]);
                h_DR_zoomin_ancestor_binned[ksign][kgrp] = new TH1D(Form("h_DR_zoomin_sign%d%s",ksign+1,ancestor_grp_labels[kgrp].c_str()),";#Delta R;1/N_{evt} dN/d#Delta R", nDR_zoomin_bins,0,pms.deltaR_thrsh_zoomin);
                h_DR_jacobian_corrected_ancestor_binned[ksign][kgrp] = new TH1D(Form("h_DR_jacobian_corrected_sign%d%s",ksign+1,ancestor_grp_labels[kgrp].c_str()),";#Delta R;#frac{1}{N_{evt}} #frac{1}{#Delta R} #frac{dN}{d#Delta R}", nDR_bins,0,pms.deltaR_thrsh[2]);
                h_DR_zoomin_jacobian_corrected_ancestor_binned[ksign][kgrp] = new TH1D(Form("h_DR_zoomin_jacobian_corrected_sign%d%s",ksign+1,ancestor_grp_labels[kgrp].c_str()),";#Delta R;#frac{1}{N_{evt}} #frac{1}{#Delta R} #frac{dN}{d#Delta R}", nDR_zoomin_bins,0,pms.deltaR_thrsh_zoomin);
                h_pt_asym_ancestor_binned[ksign][kgrp] = new TH1D(Form("h_pt_asym_sign%d%s",ksign+1,ancestor_grp_labels[kgrp].c_str()),";A = (pT1 - pT2)/(pT1 + pT2);d#sigma/dA", npt_asym_bins,0,1.);
                h_psrapidity_ordered_pt_asym_ancestor_binned[ksign][kgrp] = new TH1D(Form("h_psrapidity_ordered_pt_asym_sign%d%s",ksign+1,ancestor_grp_labels[kgrp].c_str()),";A' = (pT_{large |#eta|} - pT_{small |#eta|})/(pT_{large |#eta|} + pT_{small |#eta|});d#sigma/dA'", npt_asym_bins*2,-1.,1.);
                h_pair_pt_ptlead_ratio_ancestor_binned[ksign][kgrp] = new TH1D(Form("h_pair_pt_ptlead_ratio_sign%d%s",ksign+1,ancestor_grp_labels[kgrp].c_str()),";#frac{p_{T}^{pair}}{p_{T}^{lead}};d#sigma/d#frac{p_{T}^{pair}}{p_{T}^{lead}}", npair_pt_ptlead_ratio_bins,0,2.);
                h_pair_pt_jacobian_corrected_ancestor_binned[ksign][kgrp] = new TH1D(Form("h_pair_pt_jacobian_corrected_sign%d%s",ksign+1,ancestor_grp_labels[kgrp].c_str()),";p_{T}^{pair} [GeV];#frac{1}{N_{evt}} #frac{1}{p_{T}^{pair}} #frac{dN}{dp_{T}^{pair}}",npair_pT_bins_linear,0,30);
                
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
                h_pair_pt_jacobian_corrected_ancestor_binned[ksign][kgrp]->Sumw2();
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
                h_pair_pt_jacobian_corrected[jdphi][ksign] = new TH1D(Form("h_pair_pt_jacobian_corrected_%s_sign%d",dphi_regions[jdphi].c_str(),ksign+1),";p_{T}^{pair} [GeV];#frac{1}{N_{evt}} #frac{1}{p_{T}^{pair}} #frac{dN}{dp_{T}^{pair}}",npair_pT_bins_linear,0,30);

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
                h_pair_pt_jacobian_corrected[jdphi][ksign]->Sumw2();
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

                if (plot_kin_binned_histograms){
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
}

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

        // -------------------- TEMPORARY!!! NEED TO WRITE PAIR ETA IN NTUPLE-PROCESSING CODE!!! --------------------
        TLorentzVector M1, M2, M3;
        M1.SetPtEtaPhiM(m1pt[nfile][nkin][nsign],m1eta[nfile][nkin][nsign],m1phi[nfile][nkin][nsign], 0.105658);
        M2.SetPtEtaPhiM(m2pt[nfile][nkin][nsign],m2eta[nfile][nkin][nsign],m2phi[nfile][nkin][nsign], 0.105658);
        M3=M1+M2;
        float pair_eta = M3.Eta();
        // ------------------------------------------------------------------------------------------------------------------------

        h_crossx_truth_from_single_b_vs_pair_pt_pair_eta->Fill(pair_eta,pair_pt[nfile][nkin][nsign],ev_weight);

    }        

    // ------------------------------------------------------------------------------------------------------------------------
    // flavor-categorized histogram filling
    // ------------------------------------------------------------------------------------------------------------------------
    int flavor_ind = pair_flavor_index::other_flavor;

    if (from_same_resonance[nfile][nkin][nsign]){
        flavor_ind = pair_flavor_index::from_resonance;
    }else if (resonance_contaminated[nfile][nkin][nsign]){
        flavor_ind = pair_flavor_index::resonance_contaminated;
    }else if (from_same_b[nfile][nkin][nsign]){
        flavor_ind = pair_flavor_index::from_single_b;
    }else if (both_from_b[nfile][nkin][nsign]){
        flavor_ind = pair_flavor_index::bb;
    }else if (both_from_c[nfile][nkin][nsign]){
        flavor_ind = pair_flavor_index::cc;
    }

    h_DR_flavor_binned[nsign][flavor_ind]->Fill(dr[nfile][nkin][nsign],ev_weight);
    h_DR_zoomin_flavor_binned[nsign][flavor_ind]->Fill(dr[nfile][nkin][nsign],ev_weight);
    h_DR_jacobian_corrected_flavor_binned[nsign][flavor_ind]->Fill(dr[nfile][nkin][nsign],ev_weight * 1. / dr[nfile][nkin][nsign]);
    h_DR_zoomin_jacobian_corrected_flavor_binned[nsign][flavor_ind]->Fill(dr[nfile][nkin][nsign],ev_weight * 1. / dr[nfile][nkin][nsign]);
    h_pt_asym_flavor_binned[nsign][flavor_ind]->Fill(asym[nfile][nkin][nsign],ev_weight);
    h_psrapidity_ordered_pt_asym_flavor_binned[nsign][flavor_ind]->Fill(psrapidity_ordered_pt_asym,ev_weight); 
    h_pair_pt_ptlead_ratio_flavor_binned[nsign][flavor_ind]->Fill(pair_pt[nfile][nkin][nsign]/pt_lead[nfile][nkin][nsign],ev_weight);
    h_pair_pt_jacobian_corrected_flavor_binned[nsign][flavor_ind]->Fill(pair_pt[nfile][nkin][nsign], ev_weight * 1. / pair_pt[nfile][nkin][nsign]);
    h_Deta_Dphi_flavor_binned[nsign][flavor_ind]->Fill(dphi[nfile][nkin][nsign],deta[nfile][nkin][nsign],ev_weight);
    h_minv_pair_pt_flavor_binned[nsign][flavor_ind]->Fill(pair_pt[nfile][nkin][nsign],minv[nfile][nkin][nsign],ev_weight);
    h_minv_pair_pt_zoomin_flavor_binned[nsign][flavor_ind]->Fill(pair_pt[nfile][nkin][nsign],minv[nfile][nkin][nsign],ev_weight);
    h_minv_pair_pt_jacobian_corrected_flavor_binned[nsign][flavor_ind]->Fill(pair_pt[nfile][nkin][nsign],minv[nfile][nkin][nsign],ev_weight * 1. / dr[nfile][nkin][nsign]);
    h_pair_pt_jacobian_corrected_flavor_binned[nsign][flavor_ind]->Fill(pair_pt[nfile][nkin][nsign],ev_weight * 1. / pair_pt[nfile][nkin][nsign]);
    h_minv_pair_pt_zoomin_jacobian_corrected_flavor_binned[nsign][flavor_ind]->Fill(pair_pt[nfile][nkin][nsign],minv[nfile][nkin][nsign],ev_weight * 1. / dr[nfile][nkin][nsign]);
    h_minv_pair_pt_log_flavor_binned[nsign][flavor_ind]->Fill(pair_pt[nfile][nkin][nsign],minv[nfile][nkin][nsign],ev_weight);
    h_ptlead_pair_pt_flavor_binned[nsign][flavor_ind]->Fill(pair_pt[nfile][nkin][nsign],pt_lead[nfile][nkin][nsign],ev_weight); 
    h_ptlead_pair_pt_zoomin_flavor_binned[nsign][flavor_ind]->Fill(pair_pt[nfile][nkin][nsign],pt_lead[nfile][nkin][nsign],ev_weight); 
    h_ptlead_pair_pt_log_flavor_binned[nsign][flavor_ind]->Fill(pair_pt[nfile][nkin][nsign],pt_lead[nfile][nkin][nsign],ev_weight);             

    // ------------------------------------------------------------------------------------------------------------------------
    // ancestor/origin-categorized histogram filling
    // ------------------------------------------------------------------------------------------------------------------------
    int ancestor_ind = nAncestorGroups; // other origin categories

    if (from_same_resonance[nfile][nkin][nsign] || resonance_contaminated[nfile][nkin][nsign] || from_same_b[nfile][nkin][nsign]){
        ancestor_ind = -1;  // do not fill: "single-b", "resonance" and "resonance_contaminated" are accounted for by flavor categories; we do not want them to show up in "others"
    } else{
     // the "proper ancestor groups" are filled only for pairs both from (different) open HF's - no overlap with resonance or single-b
        for (int igrp = 0; igrp < nAncestorGroups; igrp++){
            if (muon_pair_origin_category[nfile][nkin][nsign] == ancestor_grps[igrp]){
                ancestor_ind = igrp;

                if (Qsplit[nfile][nkin][nsign] != -10.){
                    h_Qsplit_ancestor_binned[nsign][igrp]->Fill(Qsplit[nfile][nkin][nsign],ev_weight);
                    h_Qsplit_pTHat_ratio_ancestor_binned[nsign][igrp]->Fill(Qsplit[nfile][nkin][nsign] / pTHat[nfile][nkin][nsign], ev_weight);
                    if (mHard_relevant[nfile][nkin][nsign] != -10.){
                        h_Qsplit_mHat_ratio_ancestor_binned[nsign][igrp]->Fill(Qsplit[nfile][nkin][nsign] / mHard_relevant[nfile][nkin][nsign], ev_weight);
                    }
                }

                break;
            }
        }
    }

    if (ancestor_ind > 0){ // not from resonance or single-b
        h_DR_ancestor_binned[nsign][ancestor_ind]->Fill(dr[nfile][nkin][nsign],ev_weight);
        h_DR_zoomin_ancestor_binned[nsign][ancestor_ind]->Fill(dr[nfile][nkin][nsign],ev_weight);
        h_DR_jacobian_corrected_ancestor_binned[nsign][ancestor_ind]->Fill(dr[nfile][nkin][nsign],ev_weight * 1. / dr[nfile][nkin][nsign]);
        h_DR_zoomin_jacobian_corrected_ancestor_binned[nsign][ancestor_ind]->Fill(dr[nfile][nkin][nsign],ev_weight * 1. / dr[nfile][nkin][nsign]);
        h_pt_asym_ancestor_binned[nsign][ancestor_ind]->Fill(asym[nfile][nkin][nsign],ev_weight);
        h_psrapidity_ordered_pt_asym_ancestor_binned[nsign][ancestor_ind]->Fill(psrapidity_ordered_pt_asym,ev_weight); 
        h_pair_pt_ptlead_ratio_ancestor_binned[nsign][ancestor_ind]->Fill(pair_pt[nfile][nkin][nsign]/pt_lead[nfile][nkin][nsign],ev_weight);
        h_pair_pt_jacobian_corrected_ancestor_binned[nsign][ancestor_ind]->Fill(pair_pt[nfile][nkin][nsign], ev_weight * 1. / pair_pt[nfile][nkin][nsign]);
        h_Deta_Dphi_ancestor_binned[nsign][ancestor_ind]->Fill(dphi[nfile][nkin][nsign],deta[nfile][nkin][nsign],ev_weight);
        h_minv_pair_pt_ancestor_binned[nsign][ancestor_ind]->Fill(pair_pt[nfile][nkin][nsign],minv[nfile][nkin][nsign],ev_weight);
        h_minv_pair_pt_zoomin_ancestor_binned[nsign][ancestor_ind]->Fill(pair_pt[nfile][nkin][nsign],minv[nfile][nkin][nsign],ev_weight);
        h_minv_pair_pt_jacobian_corrected_ancestor_binned[nsign][ancestor_ind]->Fill(pair_pt[nfile][nkin][nsign],minv[nfile][nkin][nsign],ev_weight * 1. / dr[nfile][nkin][nsign]);
        h_minv_pair_pt_zoomin_jacobian_corrected_ancestor_binned[nsign][ancestor_ind]->Fill(pair_pt[nfile][nkin][nsign],minv[nfile][nkin][nsign],ev_weight * 1. / dr[nfile][nkin][nsign]);
        h_minv_pair_pt_log_ancestor_binned[nsign][ancestor_ind]->Fill(pair_pt[nfile][nkin][nsign],minv[nfile][nkin][nsign],ev_weight);
        h_ptlead_pair_pt_ancestor_binned[nsign][ancestor_ind]->Fill(pair_pt[nfile][nkin][nsign],pt_lead[nfile][nkin][nsign],ev_weight); 
        h_ptlead_pair_pt_zoomin_ancestor_binned[nsign][ancestor_ind]->Fill(pair_pt[nfile][nkin][nsign],pt_lead[nfile][nkin][nsign],ev_weight); 
        h_ptlead_pair_pt_log_ancestor_binned[nsign][ancestor_ind]->Fill(pair_pt[nfile][nkin][nsign],pt_lead[nfile][nkin][nsign],ev_weight);         
    }

    // ------------------------------------------------------------------------------------------------------------------------
    //un-flavor-or-origin-categorized histogram filling
    // ------------------------------------------------------------------------------------------------------------------------
    h_pair_dP_overP[away_side][nsign]->Fill(pair_dPoverP[nfile][nkin][nsign],ev_weight);
    h_pair_y[away_side][nsign]->Fill(pair_y[nfile][nkin][nsign],ev_weight);
    h_DR[away_side][nsign]->Fill(dr[nfile][nkin][nsign],ev_weight);
    h_DR_zoomin[away_side][nsign]->Fill(dr[nfile][nkin][nsign],ev_weight);
    h_DR_jacobian_corrected[away_side][nsign]->Fill(dr[nfile][nkin][nsign],ev_weight * 1. / dr[nfile][nkin][nsign]);
    h_DR_zoomin_jacobian_corrected[away_side][nsign]->Fill(dr[nfile][nkin][nsign],ev_weight * 1. / dr[nfile][nkin][nsign]);
    h_Dphi[away_side][nsign]->Fill(dphi[nfile][nkin][nsign],ev_weight);
    h_pt_asym[away_side][nsign]->Fill(asym[nfile][nkin][nsign],ev_weight);
    h_pair_pt_ptlead_ratio[away_side][nsign]->Fill(pair_pt[nfile][nkin][nsign]/pt_lead[nfile][nkin][nsign],ev_weight);
    h_pair_pt_jacobian_corrected[away_side][nsign]->Fill(pair_pt[nfile][nkin][nsign], ev_weight * 1. / pair_pt[nfile][nkin][nsign]);
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
    h_minv_pair_pt_jacobian_corrected[away_side][nsign]->Fill(pair_pt[nfile][nkin][nsign],minv[nfile][nkin][nsign],ev_weight * 1. / dr[nfile][nkin][nsign]);
    h_minv_pair_pt_zoomin_jacobian_corrected[away_side][nsign]->Fill(pair_pt[nfile][nkin][nsign],minv[nfile][nkin][nsign],ev_weight * 1. / dr[nfile][nkin][nsign]);
    h_minv_pair_pt_log[away_side][nsign]->Fill(pair_pt[nfile][nkin][nsign],minv[nfile][nkin][nsign],ev_weight);

    if (plot_kin_binned_histograms){
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
    }
    
    if (abs(dphi[nfile][nkin][nsign]) < 1){
        h_eta1_eta2_dphicut[0][nsign]->Fill(m2eta[nfile][nkin][nsign],m1eta[nfile][nkin][nsign],ev_weight);
    }else if(abs(dphi[nfile][nkin][nsign]) > pms.PI-1){
        h_eta1_eta2_dphicut[1][nsign]->Fill(m2eta[nfile][nkin][nsign],m1eta[nfile][nkin][nsign],ev_weight);
    }
}

void MuonPairPlottingPythia::FillPtBinnedHistograms(int nkin, int npt, int nsign){}


void MuonPairPlottingPythia::WriteOutput(){
    // outFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/pythia/histograms_pythia_combined.root","recreate");

    h_crossx_truth_from_single_b_vs_pair_pt_pair_eta->Write();
    
    if (mode == 1){
        // outFile->Write();
        // outFile->cd();
        // gDirectory->Delete("h_*");
    

        for (unsigned int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
            for (int jgrp = 0; jgrp < nAncestorGroupsTotal; jgrp++){
                h_DR_ancestor_binned[ksign][jgrp]->Write();
                h_DR_zoomin_ancestor_binned[ksign][jgrp]->Write();
                h_DR_jacobian_corrected_ancestor_binned[ksign][jgrp]->Write();
                h_DR_zoomin_jacobian_corrected_ancestor_binned[ksign][jgrp]->Write();
                h_pt_asym_ancestor_binned[ksign][jgrp]->Write();
                h_psrapidity_ordered_pt_asym_ancestor_binned[ksign][jgrp]->Write();
                h_pair_pt_ptlead_ratio_ancestor_binned[ksign][jgrp]->Write();
                h_pair_pt_jacobian_corrected_ancestor_binned[ksign][jgrp]->Write();
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

            for (int jflav = 0; jflav < pair_flavor_index::nFlavors; jflav++){
                h_DR_flavor_binned[ksign][jflav]->Write();
                h_DR_zoomin_flavor_binned[ksign][jflav]->Write();
                h_DR_jacobian_corrected_flavor_binned[ksign][jflav]->Write();
                h_DR_zoomin_jacobian_corrected_flavor_binned[ksign][jflav]->Write();
                h_pt_asym_flavor_binned[ksign][jflav]->Write();
                h_psrapidity_ordered_pt_asym_flavor_binned[ksign][jflav]->Write();
                h_pair_pt_ptlead_ratio_flavor_binned[ksign][jflav]->Write();
                h_pair_pt_jacobian_corrected_flavor_binned[ksign][jflav]->Write();
                h_Deta_Dphi_flavor_binned[ksign][jflav]->Write();
                h_minv_pair_pt_flavor_binned[ksign][jflav]->Write();
                h_minv_pair_pt_zoomin_flavor_binned[ksign][jflav]->Write();
                h_minv_pair_pt_jacobian_corrected_flavor_binned[ksign][jflav]->Write();
                h_minv_pair_pt_zoomin_jacobian_corrected_flavor_binned[ksign][jflav]->Write();
                h_minv_pair_pt_log_flavor_binned[ksign][jflav]->Write();
                h_ptlead_pair_pt_flavor_binned[ksign][jflav]->Write();
                h_ptlead_pair_pt_zoomin_flavor_binned[ksign][jflav]->Write();
                h_ptlead_pair_pt_log_flavor_binned[ksign][jflav]->Write();
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
                h_pair_pt_jacobian_corrected[jdphi][ksign]->Write();
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

                if (plot_kin_binned_histograms){
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