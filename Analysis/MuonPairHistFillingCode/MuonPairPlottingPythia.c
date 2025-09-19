#include "MuonPairPlottingPythia.h"
#include "time.h"
#include "TLorentzVector.h"

MuonPairPlottingPythia::MuonPairPlottingPythia(){}

void MuonPairPlottingPythia::InitInput(){

    with_data_resonance_cuts_suffix = (with_data_resonance_cuts)? "_with_data_resonance_cuts" : "_no_data_resonance_cuts";
    file_name += with_data_resonance_cuts_suffix + ".root";

    inFile = new TFile(file_name.c_str(),"read");
    if (!inFile){
        std::cout << "File with the name " << file_name << "does not exist or cannot be opened";
        throw std::exception();
    }


    for (int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
        std::string tree_name = "muon_pair_tree_sign" + std::to_string(ksign+1);
        inTree[ksign] = (TTree*) inFile->Get(tree_name.c_str());
        if (!inTree[ksign]){
            std::cout << "File with the name " << file_name << "does not have a tree named" << tree_name << std::endl;
            throw std::exception();
        }
        inTree[ksign]->SetBranchAddress("weight"          , &weight[ksign]);
        inTree[ksign]->SetBranchAddress("pTHat"          , &pTHat[ksign]);
        inTree[ksign]->SetBranchAddress("Qsplit"          , &Qsplit[ksign]);
        inTree[ksign]->SetBranchAddress("mHard_relevant"          , &mHard_relevant[ksign]);
        inTree[ksign]->SetBranchAddress("from_same_resonance"          , &from_same_resonance[ksign]);
        inTree[ksign]->SetBranchAddress("resonance_contaminated"          , &resonance_contaminated[ksign]);
        inTree[ksign]->SetBranchAddress("from_same_b"          , &from_same_b[ksign]);
        inTree[ksign]->SetBranchAddress("muon_pair_flavor_category"          , &muon_pair_flavor_category[ksign]);
        inTree[ksign]->SetBranchAddress("muon_pair_origin_category"          , &muon_pair_origin_category[ksign]);
        inTree[ksign]->SetBranchAddress("data_resonance_or_reso_contam_tagged_old"          , &data_resonance_or_reso_contam_tagged_old[ksign]);
        inTree[ksign]->SetBranchAddress("data_resonance_or_reso_contam_tagged_new"          , &data_resonance_or_reso_contam_tagged_new[ksign]);
        inTree[ksign]->SetBranchAddress("muon_pair_origin_category"          , &muon_pair_origin_category[ksign]);
        inTree[ksign]->SetBranchAddress("pair_dPoverP"    , &pair_dPoverP[ksign]);
        inTree[ksign]->SetBranchAddress("pt_lead"          , &pt_lead[ksign]);
        inTree[ksign]->SetBranchAddress("pair_pt"          , &pair_pt[ksign]);
        inTree[ksign]->SetBranchAddress("pair_eta"     , &pair_eta[ksign]);
        inTree[ksign]->SetBranchAddress("pair_y"           , &pair_y[ksign]);
        inTree[ksign]->SetBranchAddress("asym"             , &asym[ksign]);
        inTree[ksign]->SetBranchAddress("dpt"           , &dpt[ksign]);
        inTree[ksign]->SetBranchAddress("deta"       , &deta[ksign]);
        inTree[ksign]->SetBranchAddress("etaavg"      , &etaavg[ksign]);
        inTree[ksign]->SetBranchAddress("phiavg"            , &phiavg[ksign]);
        inTree[ksign]->SetBranchAddress("dphi"     , &dphi[ksign]);
        inTree[ksign]->SetBranchAddress("dr"        , &dr[ksign]);
        inTree[ksign]->SetBranchAddress("minv"        , &minv[ksign]);
        inTree[ksign]->SetBranchAddress("m1.pt"           , &m1pt[ksign]);
        inTree[ksign]->SetBranchAddress("m2.pt"           , &m2pt[ksign]);
        inTree[ksign]->SetBranchAddress("m1.eta"       , &m1eta[ksign]);
        inTree[ksign]->SetBranchAddress("m2.eta"       , &m2eta[ksign]);
        inTree[ksign]->SetBranchAddress("m1.phi"       , &m1phi[ksign]);
        inTree[ksign]->SetBranchAddress("m2.phi"       , &m2phi[ksign]);
    }
}


void MuonPairPlottingPythia::InitOutput(){
    output_file_path = "/usatlas/u/yuhanguo/usatlasdata/pythia_private_sample/histograms_pythia_combined" + with_data_resonance_cuts_suffix + ".root";
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

    // ------------------------------------------ flavor + origin hist suffices ------------------------------------------

    origin_grp_map_build(origin_grp_map);
    flavor_grp_map_build(flavor_grp_map);      

    // ------------------------------------------ cross sections ------------------------------------------
    h_crossx_truth_from_single_b_vs_pair_pt_pair_eta = new TH2D("h_crossx_truth_from_single_b_vs_pair_pt_pair_eta",";#eta^{pair};p_{T}^{pair} [GeV];#sigma^{truth}",npair_eta_bins_coarse,pair_eta_min,pair_eta_max,npair_pT_bins_coarse,pair_pt_min,pair_pt_max);
    h_crossx_truth_from_single_b_vs_pair_pt_pair_eta->Sumw2();

    // ------------------------------------------ resonance-cut study ------------------------------------------
    int minv_sub_GeV_nbins = 120;
    double minv_sub_GeV_max = 1.2;
    int minv_single_b_region_nbins = 80;
    double minv_single_b_region_max = 3.2;

    h_minv_sub_GeV_signal_no_res_cut = new TH1D("h_minv_sub_GeV_signal_no_res_cut",";m_{#mu#mu} [GeV];1/N_{evt} dN/dm_{#mu#mu}",minv_sub_GeV_nbins,0,minv_sub_GeV_max);
    h_minv_sub_GeV_signal_old_res_cut = new TH1D("h_minv_sub_GeV_signal_old_res_cut",";m_{#mu#mu} [GeV];1/N_{evt} dN/dm_{#mu#mu}",minv_sub_GeV_nbins,0,minv_sub_GeV_max);
    h_minv_sub_GeV_signal_new_res_cut = new TH1D("h_minv_sub_GeV_signal_new_res_cut",";m_{#mu#mu} [GeV];1/N_{evt} dN/dm_{#mu#mu}",minv_sub_GeV_nbins,0,minv_sub_GeV_max);
    h_minv_sub_GeV_jacobian_corrected_signal_no_res_cut = new TH1D("h_minv_sub_GeV_jacobian_corrected_signal_no_res_cut",";m_{#mu#mu};#frac{1}{N_{evt}} #frac{1}{#delta R} #frac{dN}{dm_{#mu#mu}}",minv_sub_GeV_nbins,0,minv_sub_GeV_max);
    h_minv_sub_GeV_jacobian_corrected_signal_old_res_cut = new TH1D("h_minv_sub_GeV_jacobian_corrected_signal_old_res_cut",";m_{#mu#mu};#frac{1}{N_{evt}} #frac{1}{#delta R} #frac{dN}{dm_{#mu#mu}}",minv_sub_GeV_nbins,0,minv_sub_GeV_max);
    h_minv_sub_GeV_jacobian_corrected_signal_new_res_cut = new TH1D("h_minv_sub_GeV_jacobian_corrected_signal_new_res_cut",";m_{#mu#mu};#frac{1}{N_{evt}} #frac{1}{#delta R} #frac{dN}{dm_{#mu#mu}}",minv_sub_GeV_nbins,0,minv_sub_GeV_max);

    h_minv_sub_GeV_resonance_and_res_contam_bkg_no_res_cut = new TH1D("h_minv_sub_GeV_resonance_and_res_contam_bkg_no_res_cut",";m_{#mu#mu} [GeV];1/N_{evt} dN/dm_{#mu#mu}",minv_sub_GeV_nbins,0,minv_sub_GeV_max);
    h_minv_sub_GeV_resonance_and_res_contam_bkg_old_res_cut = new TH1D("h_minv_sub_GeV_resonance_and_res_contam_bkg_old_res_cut",";m_{#mu#mu} [GeV];1/N_{evt} dN/dm_{#mu#mu}",minv_sub_GeV_nbins,0,minv_sub_GeV_max);
    h_minv_sub_GeV_resonance_and_res_contam_bkg_new_res_cut = new TH1D("h_minv_sub_GeV_resonance_and_res_contam_bkg_new_res_cut",";m_{#mu#mu} [GeV];1/N_{evt} dN/dm_{#mu#mu}",minv_sub_GeV_nbins,0,minv_sub_GeV_max);
    h_minv_sub_GeV_jacobian_corrected_resonance_and_res_contam_bkg_no_res_cut = new TH1D("h_minv_sub_GeV_jacobian_corrected_resonance_and_res_contam_bkg_no_res_cut",";m_{#mu#mu};#frac{1}{N_{evt}} #frac{1}{#delta R} #frac{dN}{dm_{#mu#mu}}",minv_sub_GeV_nbins,0,minv_sub_GeV_max);
    h_minv_sub_GeV_jacobian_corrected_resonance_and_res_contam_bkg_old_res_cut = new TH1D("h_minv_sub_GeV_jacobian_corrected_resonance_and_res_contam_bkg_old_res_cut",";m_{#mu#mu};#frac{1}{N_{evt}} #frac{1}{#delta R} #frac{dN}{dm_{#mu#mu}}",minv_sub_GeV_nbins,0,minv_sub_GeV_max);
    h_minv_sub_GeV_jacobian_corrected_resonance_and_res_contam_bkg_new_res_cut = new TH1D("h_minv_sub_GeV_jacobian_corrected_resonance_and_res_contam_bkg_new_res_cut",";m_{#mu#mu};#frac{1}{N_{evt}} #frac{1}{#delta R} #frac{dN}{dm_{#mu#mu}}",minv_sub_GeV_nbins,0,minv_sub_GeV_max);

    h_minv_single_b_region_signal_no_res_cut = new TH1D("h_minv_single_b_region_signal_no_res_cut",";m_{#mu#mu} [GeV];1/N_{evt} dN/dm_{#mu#mu}",minv_single_b_region_nbins,0,minv_single_b_region_max);
    h_minv_single_b_region_signal_old_res_cut = new TH1D("h_minv_single_b_region_signal_old_res_cut",";m_{#mu#mu} [GeV];1/N_{evt} dN/dm_{#mu#mu}",minv_single_b_region_nbins,0,minv_single_b_region_max);
    h_minv_single_b_region_signal_new_res_cut = new TH1D("h_minv_single_b_region_signal_new_res_cut",";m_{#mu#mu} [GeV];1/N_{evt} dN/dm_{#mu#mu}",minv_single_b_region_nbins,0,minv_single_b_region_max);
    h_minv_single_b_region_jacobian_corrected_signal_no_res_cut = new TH1D("h_minv_single_b_region_jacobian_corrected_signal_no_res_cut",";m_{#mu#mu};#frac{1}{N_{evt}} #frac{1}{#delta R} #frac{dN}{dm_{#mu#mu}}",minv_single_b_region_nbins,0,minv_single_b_region_max);
    h_minv_single_b_region_jacobian_corrected_signal_old_res_cut = new TH1D("h_minv_single_b_region_jacobian_corrected_signal_old_res_cut",";m_{#mu#mu};#frac{1}{N_{evt}} #frac{1}{#delta R} #frac{dN}{dm_{#mu#mu}}",minv_single_b_region_nbins,0,minv_single_b_region_max);
    h_minv_single_b_region_jacobian_corrected_signal_new_res_cut = new TH1D("h_minv_single_b_region_jacobian_corrected_signal_new_res_cut",";m_{#mu#mu};#frac{1}{N_{evt}} #frac{1}{#delta R} #frac{dN}{dm_{#mu#mu}}",minv_single_b_region_nbins,0,minv_single_b_region_max);

    h_minv_single_b_region_resonance_and_res_contam_bkg_no_res_cut = new TH1D("h_minv_single_b_region_resonance_and_res_contam_bkg_no_res_cut",";m_{#mu#mu} [GeV];1/N_{evt} dN/dm_{#mu#mu}",minv_single_b_region_nbins,0,minv_single_b_region_max);
    h_minv_single_b_region_resonance_and_res_contam_bkg_old_res_cut = new TH1D("h_minv_single_b_region_resonance_and_res_contam_bkg_old_res_cut",";m_{#mu#mu} [GeV];1/N_{evt} dN/dm_{#mu#mu}",minv_single_b_region_nbins,0,minv_single_b_region_max);
    h_minv_single_b_region_resonance_and_res_contam_bkg_new_res_cut = new TH1D("h_minv_single_b_region_resonance_and_res_contam_bkg_new_res_cut",";m_{#mu#mu} [GeV];1/N_{evt} dN/dm_{#mu#mu}",minv_single_b_region_nbins,0,minv_single_b_region_max);
    h_minv_single_b_region_jacobian_corrected_resonance_and_res_contam_bkg_no_res_cut = new TH1D("h_minv_single_b_region_jacobian_corrected_resonance_and_res_contam_bkg_no_res_cut",";m_{#mu#mu};#frac{1}{N_{evt}} #frac{1}{#delta R} #frac{dN}{dm_{#mu#mu}}",minv_single_b_region_nbins,0,minv_single_b_region_max);
    h_minv_single_b_region_jacobian_corrected_resonance_and_res_contam_bkg_old_res_cut = new TH1D("h_minv_single_b_region_jacobian_corrected_resonance_and_res_contam_bkg_old_res_cut",";m_{#mu#mu};#frac{1}{N_{evt}} #frac{1}{#delta R} #frac{dN}{dm_{#mu#mu}}",minv_single_b_region_nbins,0,minv_single_b_region_max);
    h_minv_single_b_region_jacobian_corrected_resonance_and_res_contam_bkg_new_res_cut = new TH1D("h_minv_single_b_region_jacobian_corrected_resonance_and_res_contam_bkg_new_res_cut",";m_{#mu#mu};#frac{1}{N_{evt}} #frac{1}{#delta R} #frac{dN}{dm_{#mu#mu}}",minv_single_b_region_nbins,0,minv_single_b_region_max);

    // ------------------------------------------ kinematic distributions ------------------------------------------
    if (mode == 1){
        for (unsigned int ksign = 0; ksign < ParamsSet::nSigns; ksign++){

            for (int kflav = 0; kflav < pair_flavor_index::nFlavors; kflav++){
                h_DR_flavor_binned[ksign][kflav] = new TH1D(Form("h_DR_sign%d%s",ksign+1,flavor_grp_map[kflav].c_str()),";#Delta R;1/N_{evt} dN/d#Delta R", nDR_bins,0,pms.deltaR_thrsh[2]);
                h_DR_zoomin_flavor_binned[ksign][kflav] = new TH1D(Form("h_DR_zoomin_sign%d%s",ksign+1,flavor_grp_map[kflav].c_str()),";#Delta R;1/N_{evt} dN/d#Delta R", nDR_zoomin_bins,0,pms.deltaR_thrsh_zoomin);
                h_Deta_zoomin_flavor_binned[ksign][kflav] = new TH1D(Form("h_Deta_zoomin_sign%d%s",ksign+1,flavor_grp_map[kflav].c_str()),";#Delta #eta;1/N_{evt} dN/d#Delta #eta", nDR_zoomin_bins,-1.,1.);
                h_Dphi_zoomin_flavor_binned[ksign][kflav] = new TH1D(Form("h_Dphi_zoomin_sign%d%s",ksign+1,flavor_grp_map[kflav].c_str()),";#Delta #phi;1/N_{evt} dN/d#Delta #phi", nDR_zoomin_bins,-1.,1.);
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
                h_Deta_zoomin_flavor_binned[ksign][kflav]->Sumw2();
                h_Dphi_zoomin_flavor_binned[ksign][kflav]->Sumw2();
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


            for (int korigin = 0; korigin < muon_pair_both_from_open_HF_origin_catgr::nOrigins; korigin++){
                h_DR_origin_binned[ksign][korigin] = new TH1D(Form("h_DR_sign%d%s",ksign+1,origin_grp_map[korigin].c_str()),";#Delta R;1/N_{evt} dN/d#Delta R", nDR_bins,0,pms.deltaR_thrsh[2]);
                h_DR_zoomin_origin_binned[ksign][korigin] = new TH1D(Form("h_DR_zoomin_sign%d%s",ksign+1,origin_grp_map[korigin].c_str()),";#Delta R;1/N_{evt} dN/d#Delta R", nDR_zoomin_bins,0,pms.deltaR_thrsh_zoomin);
                h_Deta_zoomin_origin_binned[ksign][korigin] = new TH1D(Form("h_Deta_zoomin_sign%d%s",ksign+1,origin_grp_map[korigin].c_str()),";#Delta #eta;1/N_{evt} dN/d#Delta #eta", nDR_zoomin_bins,-1,1);
                h_Dphi_zoomin_origin_binned[ksign][korigin] = new TH1D(Form("h_Dphi_zoomin_sign%d%s",ksign+1,origin_grp_map[korigin].c_str()),";#Delta #phi;1/N_{evt} dN/d#Delta #phi", nDR_zoomin_bins,-1,1);
                h_DR_jacobian_corrected_origin_binned[ksign][korigin] = new TH1D(Form("h_DR_jacobian_corrected_sign%d%s",ksign+1,origin_grp_map[korigin].c_str()),";#Delta R;#frac{1}{N_{evt}} #frac{1}{#Delta R} #frac{dN}{d#Delta R}", nDR_bins,0,pms.deltaR_thrsh[2]);
                h_DR_zoomin_jacobian_corrected_origin_binned[ksign][korigin] = new TH1D(Form("h_DR_zoomin_jacobian_corrected_sign%d%s",ksign+1,origin_grp_map[korigin].c_str()),";#Delta R;#frac{1}{N_{evt}} #frac{1}{#Delta R} #frac{dN}{d#Delta R}", nDR_zoomin_bins,0,pms.deltaR_thrsh_zoomin);
                h_pt_asym_origin_binned[ksign][korigin] = new TH1D(Form("h_pt_asym_sign%d%s",ksign+1,origin_grp_map[korigin].c_str()),";A = (pT1 - pT2)/(pT1 + pT2);d#sigma/dA", npt_asym_bins,0,1.);
                h_psrapidity_ordered_pt_asym_origin_binned[ksign][korigin] = new TH1D(Form("h_psrapidity_ordered_pt_asym_sign%d%s",ksign+1,origin_grp_map[korigin].c_str()),";A' = (pT_{large |#eta|} - pT_{small |#eta|})/(pT_{large |#eta|} + pT_{small |#eta|});d#sigma/dA'", npt_asym_bins*2,-1.,1.);
                h_pair_pt_ptlead_ratio_origin_binned[ksign][korigin] = new TH1D(Form("h_pair_pt_ptlead_ratio_sign%d%s",ksign+1,origin_grp_map[korigin].c_str()),";#frac{p_{T}^{pair}}{p_{T}^{lead}};d#sigma/d#frac{p_{T}^{pair}}{p_{T}^{lead}}", npair_pt_ptlead_ratio_bins,0,2.);
                h_pair_pt_jacobian_corrected_origin_binned[ksign][korigin] = new TH1D(Form("h_pair_pt_jacobian_corrected_sign%d%s",ksign+1,origin_grp_map[korigin].c_str()),";p_{T}^{pair} [GeV];#frac{1}{N_{evt}} #frac{1}{p_{T}^{pair}} #frac{dN}{dp_{T}^{pair}}",npair_pT_bins_linear,0,30);
                
                h_Deta_Dphi_origin_binned[ksign][korigin] = new TH2D(Form("h_Deta_Dphi_sign%d%s",ksign+1,origin_grp_map[korigin].c_str()),";#Delta#phi;#Delta#eta", nDphi_bins,-pms.PI,pms.PI,nDeta_bins,-4.8,4.8);
                h_minv_pair_pt_origin_binned[ksign][korigin] = new TH2D(Form("h_minv_pair_pt_sign%d%s",ksign+1,origin_grp_map[korigin].c_str()),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",npair_pT_bins_linear,0,30,nminv_bins_linear,0,30);
                h_minv_pair_pt_zoomin_origin_binned[ksign][korigin] = new TH2D(Form("h_minv_pair_pt_zoomin_sign%d%s",ksign+1,origin_grp_map[korigin].c_str()),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",npair_pT_bins_linear,0,20,nminv_bins_linear,0,4);
                h_minv_pair_pt_jacobian_corrected_origin_binned[ksign][korigin] = new TH2D(Form("h_minv_pair_pt_jacobian_corrected_sign%d%s",ksign+1,origin_grp_map[korigin].c_str()),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",npair_pT_bins_linear,0,30,nminv_bins_linear,0,30);
                h_minv_pair_pt_zoomin_jacobian_corrected_origin_binned[ksign][korigin] = new TH2D(Form("h_minv_pair_pt_zoomin_jacobian_corrected_sign%d%s",ksign+1,origin_grp_map[korigin].c_str()),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",npair_pT_bins_linear,0,20,nminv_bins_linear,0,4);
                h_minv_pair_pt_log_origin_binned[ksign][korigin] = new TH2D(Form("h_minv_pair_pt_log_sign%d%s",ksign+1,origin_grp_map[korigin].c_str()),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",pms.npairPT_bins,pms.pairPTBins[0][2],nminv_bins_log,minv_bins_log[0]);
                h_ptlead_pair_pt_origin_binned[ksign][korigin] = new TH2D(Form("h_ptlead_pair_pt_sign%d%s",ksign+1,origin_grp_map[korigin].c_str()),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",npair_pT_bins_linear,0,30,npT_lead_bins_linear,0,30);
                h_ptlead_pair_pt_zoomin_origin_binned[ksign][korigin] = new TH2D(Form("h_ptlead_pair_pt_zoomin_sign%d%s",ksign+1,origin_grp_map[korigin].c_str()),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",npair_pT_bins_linear,0,20,npT_lead_bins_linear,0,15);
                h_ptlead_pair_pt_log_origin_binned[ksign][korigin] = new TH2D(Form("h_ptlead_pair_pt_log_sign%d%s",ksign+1,origin_grp_map[korigin].c_str()),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",pms.npairPT_bins,pms.pairPTBins[0][2],pms.npt_bins,pms.pTBins);

                h_Qsplit_origin_binned[ksign][korigin] = new TH1D(Form("h_Qsplit_sign%d%s",ksign+1,origin_grp_map[korigin].c_str()),";Q_{split}", 160, 0, 80.);
                h_Qsplit_pTHat_ratio_origin_binned[ksign][korigin] = new TH1D(Form("h_Qsplit_pTHat_ratio_sign%d%s",ksign+1,origin_grp_map[korigin].c_str()),";#frac{Q_{split}}{#hat{p}_{T}}};d#sigma/d#frac{Q_{split}}{#hat{p}_{T}}}", 100,-2.,2.);
                h_Qsplit_mHat_ratio_origin_binned[ksign][korigin] = new TH1D(Form("h_Qsplit_mHat_ratio_sign%d%s",ksign+1,origin_grp_map[korigin].c_str()),";#frac{Q_{split}}{#sqrt{#hat{s}}};d#sigma/d#frac{Q_{split}}{#sqrt{#hat{s}}}", 140,-0.7,0.7);
                
                h_DR_origin_binned[ksign][korigin]->Sumw2();
                h_DR_zoomin_origin_binned[ksign][korigin]->Sumw2();
                h_Deta_zoomin_origin_binned[ksign][korigin]->Sumw2();
                h_Dphi_zoomin_origin_binned[ksign][korigin]->Sumw2();
                h_DR_jacobian_corrected_origin_binned[ksign][korigin]->Sumw2();
                h_DR_zoomin_jacobian_corrected_origin_binned[ksign][korigin]->Sumw2();
                h_pt_asym_origin_binned[ksign][korigin]->Sumw2();
                h_psrapidity_ordered_pt_asym_origin_binned[ksign][korigin]->Sumw2();
                h_pair_pt_ptlead_ratio_origin_binned[ksign][korigin]->Sumw2();
                h_pair_pt_jacobian_corrected_origin_binned[ksign][korigin]->Sumw2();
                h_Deta_Dphi_origin_binned[ksign][korigin]->Sumw2();
                h_minv_pair_pt_origin_binned[ksign][korigin]->Sumw2();
                h_minv_pair_pt_zoomin_origin_binned[ksign][korigin]->Sumw2();
                h_minv_pair_pt_jacobian_corrected_origin_binned[ksign][korigin]->Sumw2();
                h_minv_pair_pt_zoomin_jacobian_corrected_origin_binned[ksign][korigin]->Sumw2();
                h_minv_pair_pt_log_origin_binned[ksign][korigin]->Sumw2();
                h_ptlead_pair_pt_origin_binned[ksign][korigin]->Sumw2();
                h_ptlead_pair_pt_zoomin_origin_binned[ksign][korigin]->Sumw2();
                h_ptlead_pair_pt_log_origin_binned[ksign][korigin]->Sumw2();
                h_Qsplit_origin_binned[ksign][korigin]->Sumw2();
                h_Qsplit_pTHat_ratio_origin_binned[ksign][korigin]->Sumw2();
                h_Qsplit_mHat_ratio_origin_binned[ksign][korigin]->Sumw2();
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
                h_Deta_zoomin[jdphi][ksign] = new TH1D(Form("h_Deta_zoomin_%s_sign%d",dphi_regions[jdphi].c_str(),ksign+1),";#Delta #eta;1/N_{evt} dN/d#Delta #eta", nDR_zoomin_bins,-1,1);
                h_Dphi_zoomin[jdphi][ksign] = new TH1D(Form("h_Dphi_zoomin_%s_sign%d",dphi_regions[jdphi].c_str(),ksign+1),";#Delta #phi;1/N_{evt} dN/d#Delta #phi", nDR_zoomin_bins,-1,1);
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
                h_Deta_zoomin[jdphi][ksign]->Sumw2();
                h_Dphi_zoomin[jdphi][ksign]->Sumw2();
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
            }
        }
    }
}

void MuonPairPlottingPythia::ProcessData(){
    for (int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
        
        Long64_t nentries = inTree[ksign]->GetEntries(); //#muon pairs
        for (Long64_t lentry=0; lentry<nentries; lentry++) {//loop over the events
        // for (Long64_t lentry=0; lentry<100; lentry++) {//loop over the events
        if(lentry%100000==0) cout<<"Processing "<<lentry<<" event out of "<<nentries<<" events"<<std::endl;
      		int num_bytes = inTree[ksign]->GetEntry(lentry);//read in an event
            if(num_bytes==0){
            	std::cout<<"Error:: Read in event has size of zero bytes,  quitting"<<std::endl;
            	throw std::exception();
            }
            if(mode == 1){
            	FillHistograms(ksign);
            }else{ // mode has to be 1 or 3
            	// FillPtBinnedHistograms(jpt, ksign);
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


void MuonPairPlottingPythia::FillHistograms(int nsign){

    double ev_weight = weight[nsign];

    // ngapcut = 0: all; = 1: only those that pass

    // bool pass_gapcut = PassSingleMuonGapCut(m1eta[nsign], m1pt[nsign], m1charge[nsign]) && PassSingleMuonGapCut(m2eta[nsign], m2pt[nsign], m2charge[nsign]);
    bool away_side = (abs(dphi[nsign]) >= pms.PI / 2.);

    float pT_large_eta = (abs(m1eta[nsign]) > abs(m2eta[nsign]))? m1pt[nsign] : m2pt[nsign];
    float pT_small_eta = (abs(m1eta[nsign]) > abs(m2eta[nsign]))? m2pt[nsign] : m1pt[nsign];
    float psrapidity_ordered_pt_asym = (pT_large_eta - pT_small_eta)/(pT_large_eta + pT_small_eta); // expect to be peaked towards negative value if there is asymmetry

    bool is_filled = false;
    if (from_same_b[nsign]){

        h_crossx_truth_from_single_b_vs_pair_pt_pair_eta->Fill(pair_eta[nsign],pair_pt[nsign],ev_weight);

        h_minv_sub_GeV_signal_no_res_cut->Fill(minv[nsign],ev_weight);
        h_minv_sub_GeV_jacobian_corrected_signal_no_res_cut->Fill(minv[nsign],ev_weight * 1. / dr[nsign]);
        h_minv_single_b_region_signal_no_res_cut->Fill(minv[nsign],ev_weight);
        h_minv_single_b_region_jacobian_corrected_signal_no_res_cut->Fill(minv[nsign],ev_weight * 1. / dr[nsign]);

        if (!data_resonance_or_reso_contam_tagged_old[nsign]){
            h_minv_sub_GeV_signal_old_res_cut->Fill(minv[nsign],ev_weight);
            h_minv_sub_GeV_jacobian_corrected_signal_old_res_cut->Fill(minv[nsign],ev_weight * 1. / dr[nsign]);
            h_minv_single_b_region_signal_old_res_cut->Fill(minv[nsign],ev_weight);
            h_minv_single_b_region_jacobian_corrected_signal_old_res_cut->Fill(minv[nsign],ev_weight * 1. / dr[nsign]);
        }
        if (!data_resonance_or_reso_contam_tagged_new[nsign]){
            h_minv_sub_GeV_signal_new_res_cut->Fill(minv[nsign],ev_weight);
            h_minv_sub_GeV_jacobian_corrected_signal_new_res_cut->Fill(minv[nsign],ev_weight * 1. / dr[nsign]);
            h_minv_single_b_region_signal_new_res_cut->Fill(minv[nsign],ev_weight);
            h_minv_single_b_region_jacobian_corrected_signal_new_res_cut->Fill(minv[nsign],ev_weight * 1. / dr[nsign]);
        }
    }

    if (from_same_resonance[nsign] || resonance_contaminated[nsign]){
        h_minv_sub_GeV_resonance_and_res_contam_bkg_no_res_cut->Fill(minv[nsign],ev_weight);
        h_minv_sub_GeV_jacobian_corrected_resonance_and_res_contam_bkg_no_res_cut->Fill(minv[nsign],ev_weight * 1. / dr[nsign]);
        h_minv_single_b_region_resonance_and_res_contam_bkg_no_res_cut->Fill(minv[nsign],ev_weight);
        h_minv_single_b_region_jacobian_corrected_resonance_and_res_contam_bkg_no_res_cut->Fill(minv[nsign],ev_weight * 1. / dr[nsign]);

        if (!data_resonance_or_reso_contam_tagged_old[nsign]){
            h_minv_sub_GeV_resonance_and_res_contam_bkg_old_res_cut->Fill(minv[nsign],ev_weight);
            h_minv_sub_GeV_jacobian_corrected_resonance_and_res_contam_bkg_old_res_cut->Fill(minv[nsign],ev_weight * 1. / dr[nsign]);
            h_minv_single_b_region_resonance_and_res_contam_bkg_old_res_cut->Fill(minv[nsign],ev_weight);
            h_minv_single_b_region_jacobian_corrected_resonance_and_res_contam_bkg_old_res_cut->Fill(minv[nsign],ev_weight * 1. / dr[nsign]);
        }
        if (!data_resonance_or_reso_contam_tagged_new[nsign]){
            h_minv_sub_GeV_resonance_and_res_contam_bkg_new_res_cut->Fill(minv[nsign],ev_weight);
            h_minv_sub_GeV_jacobian_corrected_resonance_and_res_contam_bkg_new_res_cut->Fill(minv[nsign],ev_weight * 1. / dr[nsign]);
            h_minv_single_b_region_resonance_and_res_contam_bkg_new_res_cut->Fill(minv[nsign],ev_weight);
            h_minv_single_b_region_jacobian_corrected_resonance_and_res_contam_bkg_new_res_cut->Fill(minv[nsign],ev_weight * 1. / dr[nsign]);
        }
    }

    // ------------------------------------------------------------------------------------------------------------------------
    // flavor-categorized histogram filling
    // ------------------------------------------------------------------------------------------------------------------------
    
    int iflavor = muon_pair_flavor_category[nsign];

    if (iflavor < pair_flavor_index::nFlavors){ // safety guard
        h_DR_flavor_binned[nsign][iflavor]->Fill(dr[nsign],ev_weight);
        h_DR_zoomin_flavor_binned[nsign][iflavor]->Fill(dr[nsign],ev_weight);
        h_Deta_zoomin_flavor_binned[nsign][iflavor]->Fill(deta[nsign],ev_weight);
        h_Dphi_zoomin_flavor_binned[nsign][iflavor]->Fill(dphi[nsign],ev_weight);
        h_DR_jacobian_corrected_flavor_binned[nsign][iflavor]->Fill(dr[nsign],ev_weight * 1. / dr[nsign]);
        h_DR_zoomin_jacobian_corrected_flavor_binned[nsign][iflavor]->Fill(dr[nsign],ev_weight * 1. / dr[nsign]);
        h_pt_asym_flavor_binned[nsign][iflavor]->Fill(asym[nsign],ev_weight);
        h_psrapidity_ordered_pt_asym_flavor_binned[nsign][iflavor]->Fill(psrapidity_ordered_pt_asym,ev_weight); 
        h_pair_pt_ptlead_ratio_flavor_binned[nsign][iflavor]->Fill(pair_pt[nsign]/pt_lead[nsign],ev_weight);
        h_pair_pt_jacobian_corrected_flavor_binned[nsign][iflavor]->Fill(pair_pt[nsign], ev_weight * 1. / pair_pt[nsign]);
        h_Deta_Dphi_flavor_binned[nsign][iflavor]->Fill(dphi[nsign],deta[nsign],ev_weight);
        h_minv_pair_pt_flavor_binned[nsign][iflavor]->Fill(pair_pt[nsign],minv[nsign],ev_weight);
        h_minv_pair_pt_zoomin_flavor_binned[nsign][iflavor]->Fill(pair_pt[nsign],minv[nsign],ev_weight);
        h_minv_pair_pt_jacobian_corrected_flavor_binned[nsign][iflavor]->Fill(pair_pt[nsign],minv[nsign],ev_weight * 1. / dr[nsign]);
        h_pair_pt_jacobian_corrected_flavor_binned[nsign][iflavor]->Fill(pair_pt[nsign],ev_weight * 1. / pair_pt[nsign]);
        h_minv_pair_pt_zoomin_jacobian_corrected_flavor_binned[nsign][iflavor]->Fill(pair_pt[nsign],minv[nsign],ev_weight * 1. / dr[nsign]);
        h_minv_pair_pt_log_flavor_binned[nsign][iflavor]->Fill(pair_pt[nsign],minv[nsign],ev_weight);
        h_ptlead_pair_pt_flavor_binned[nsign][iflavor]->Fill(pair_pt[nsign],pt_lead[nsign],ev_weight);
        h_ptlead_pair_pt_zoomin_flavor_binned[nsign][iflavor]->Fill(pair_pt[nsign],pt_lead[nsign],ev_weight);
        h_ptlead_pair_pt_log_flavor_binned[nsign][iflavor]->Fill(pair_pt[nsign],pt_lead[nsign],ev_weight);
    }

    // ------------------------------------------------------------------------------------------------------------------------
    // ancestor/origin-categorized histogram filling
    // ------------------------------------------------------------------------------------------------------------------------

    int iorigin = muon_pair_origin_category[nsign];

    if (iorigin < muon_pair_both_from_open_HF_origin_catgr::nOrigins){
        h_DR_origin_binned[nsign][iorigin]->Fill(dr[nsign],ev_weight);
        h_DR_zoomin_origin_binned[nsign][iorigin]->Fill(dr[nsign],ev_weight);
        h_Deta_zoomin_origin_binned[nsign][iorigin]->Fill(deta[nsign],ev_weight);
        h_Dphi_zoomin_origin_binned[nsign][iorigin]->Fill(dphi[nsign],ev_weight);
        h_DR_jacobian_corrected_origin_binned[nsign][iorigin]->Fill(dr[nsign],ev_weight * 1. / dr[nsign]);
        h_DR_zoomin_jacobian_corrected_origin_binned[nsign][iorigin]->Fill(dr[nsign],ev_weight * 1. / dr[nsign]);
        h_pt_asym_origin_binned[nsign][iorigin]->Fill(asym[nsign],ev_weight);
        h_psrapidity_ordered_pt_asym_origin_binned[nsign][iorigin]->Fill(psrapidity_ordered_pt_asym,ev_weight); 
        h_pair_pt_ptlead_ratio_origin_binned[nsign][iorigin]->Fill(pair_pt[nsign]/pt_lead[nsign],ev_weight);
        h_pair_pt_jacobian_corrected_origin_binned[nsign][iorigin]->Fill(pair_pt[nsign], ev_weight * 1. / pair_pt[nsign]);
        h_Deta_Dphi_origin_binned[nsign][iorigin]->Fill(dphi[nsign],deta[nsign],ev_weight);
        h_minv_pair_pt_origin_binned[nsign][iorigin]->Fill(pair_pt[nsign],minv[nsign],ev_weight);
        h_minv_pair_pt_zoomin_origin_binned[nsign][iorigin]->Fill(pair_pt[nsign],minv[nsign],ev_weight);
        h_minv_pair_pt_jacobian_corrected_origin_binned[nsign][iorigin]->Fill(pair_pt[nsign],minv[nsign],ev_weight * 1. / dr[nsign]);
        h_minv_pair_pt_zoomin_jacobian_corrected_origin_binned[nsign][iorigin]->Fill(pair_pt[nsign],minv[nsign],ev_weight * 1. / dr[nsign]);
        h_minv_pair_pt_log_origin_binned[nsign][iorigin]->Fill(pair_pt[nsign],minv[nsign],ev_weight);
        h_ptlead_pair_pt_origin_binned[nsign][iorigin]->Fill(pair_pt[nsign],pt_lead[nsign],ev_weight); 
        h_ptlead_pair_pt_zoomin_origin_binned[nsign][iorigin]->Fill(pair_pt[nsign],pt_lead[nsign],ev_weight); 
        h_ptlead_pair_pt_log_origin_binned[nsign][iorigin]->Fill(pair_pt[nsign],pt_lead[nsign],ev_weight);         

        if (Qsplit[nsign] != -10.){
            h_Qsplit_origin_binned[nsign][iorigin]->Fill(Qsplit[nsign],ev_weight);
            h_Qsplit_pTHat_ratio_origin_binned[nsign][iorigin]->Fill(Qsplit[nsign] / pTHat[nsign], ev_weight);
            if (mHard_relevant[nsign] != -10.){
                h_Qsplit_mHat_ratio_origin_binned[nsign][iorigin]->Fill(Qsplit[nsign] / mHard_relevant[nsign], ev_weight);
            }
        }
    }

    // ------------------------------------------------------------------------------------------------------------------------
    //un-flavor-or-origin-categorized histogram filling
    // ------------------------------------------------------------------------------------------------------------------------
    h_pair_dP_overP[away_side][nsign]->Fill(pair_dPoverP[nsign],ev_weight);
    h_pair_y[away_side][nsign]->Fill(pair_y[nsign],ev_weight);
    h_DR[away_side][nsign]->Fill(dr[nsign],ev_weight);
    h_DR_zoomin[away_side][nsign]->Fill(dr[nsign],ev_weight);
    h_Deta_zoomin[away_side][nsign]->Fill(deta[nsign],ev_weight);
    h_Dphi_zoomin[away_side][nsign]->Fill(dphi[nsign],ev_weight);
    h_DR_jacobian_corrected[away_side][nsign]->Fill(dr[nsign],ev_weight * 1. / dr[nsign]);
    h_DR_zoomin_jacobian_corrected[away_side][nsign]->Fill(dr[nsign],ev_weight * 1. / dr[nsign]);
    h_Dphi[away_side][nsign]->Fill(dphi[nsign],ev_weight);
    h_pt_asym[away_side][nsign]->Fill(asym[nsign],ev_weight);
    h_pair_pt_ptlead_ratio[away_side][nsign]->Fill(pair_pt[nsign]/pt_lead[nsign],ev_weight);
    h_pair_pt_jacobian_corrected[away_side][nsign]->Fill(pair_pt[nsign], ev_weight * 1. / pair_pt[nsign]);
    h_eta_avg_Dphi[away_side][nsign]->Fill(dphi[nsign],etaavg[nsign],ev_weight);
    h_Deta_Dphi[away_side][nsign]->Fill(dphi[nsign],deta[nsign],ev_weight);
    h_eta1_eta2[away_side][nsign]->Fill(m2eta[nsign],m1eta[nsign],ev_weight);
    h_pt1_pt2[away_side][nsign]->Fill(m2pt[nsign],m1pt[nsign],ev_weight);
    h_eta_avg_Deta[away_side][nsign]->Fill(deta[nsign],etaavg[nsign],ev_weight);
    h_ptlead_pair_pt[away_side][nsign]->Fill(pair_pt[nsign],m1pt[nsign],ev_weight);
    h_ptlead_pair_pt_zoomin[away_side][nsign]->Fill(pair_pt[nsign],m1pt[nsign],ev_weight);
    h_ptlead_pair_pt_log[away_side][nsign]->Fill(pair_pt[nsign],m1pt[nsign],ev_weight);
    h_minv_pair_pt[away_side][nsign]->Fill(pair_pt[nsign],minv[nsign],ev_weight);
    h_minv_pair_pt_zoomin[away_side][nsign]->Fill(pair_pt[nsign],minv[nsign],ev_weight);
    h_minv_pair_pt_jacobian_corrected[away_side][nsign]->Fill(pair_pt[nsign],minv[nsign],ev_weight * 1. / dr[nsign]);
    h_minv_pair_pt_zoomin_jacobian_corrected[away_side][nsign]->Fill(pair_pt[nsign],minv[nsign],ev_weight * 1. / dr[nsign]);
    h_minv_pair_pt_log[away_side][nsign]->Fill(pair_pt[nsign],minv[nsign],ev_weight);

    if (abs(dphi[nsign]) < 1){
        h_eta1_eta2_dphicut[0][nsign]->Fill(m2eta[nsign],m1eta[nsign],ev_weight);
    }else if(abs(dphi[nsign]) > pms.PI-1){
        h_eta1_eta2_dphicut[1][nsign]->Fill(m2eta[nsign],m1eta[nsign],ev_weight);
    }
}

void MuonPairPlottingPythia::FillPtBinnedHistograms(int npt, int nsign){}


void MuonPairPlottingPythia::WriteOutput(){
    // outFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/pythia_private_sample/histograms_pythia_combined.root","recreate");

    h_crossx_truth_from_single_b_vs_pair_pt_pair_eta->Write();
    
    if (mode == 1){
        outFile->Write();
        outFile->Close();
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