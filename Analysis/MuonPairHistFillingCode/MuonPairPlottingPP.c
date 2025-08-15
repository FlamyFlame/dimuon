#include "MuonPairPlottingPP.h"
#include "time.h"

MuonPairPlottingPP::MuonPairPlottingPP(){
    isScram = false;
    isTight = false;
    
    isRun3 = true;
    // isMu4mu4noL1 = true;
}

void MuonPairPlottingPP::InitInput(){


    if (isScram){
        inFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/scrambled_muon_pairs_pp.root","read");
    }else if (isTight){
        inFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/muon_pairs_pp_tight.root","read");
    }else{
        if (isRun3){
            inFile = new TFile(Form("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/muon_pairs_pp_2024_combined%s.root", trig_suffix.c_str()),"read");            
        }else{
            inFile = new TFile(Form("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_run2/muon_pairs_pp_2017_combined%s.root", trig_suffix.c_str()),"read");
        }
    }

    if (!inFile || inFile->IsZombie()) {
        std::cerr << "Error opening input file!" << std::endl;
        throw std::exception();
    }

    for (int isign = 0; isign < ParamsSet::nSigns; isign++){
        if (isScram){
            inTree[isign] = (TTree*) inFile->Get(Form("scramb_muon_pair_tree_sign%d",isign+1));  
        }else{
            inTree[isign] = (TTree*) inFile->Get(Form("muon_pair_tree_sign%d",isign+1));  
        }
        inTree[isign]->SetBranchAddress("weight"          , &weight[isign]);
        inTree[isign]->SetBranchAddress("pair_dPoverP"    , &pair_dPoverP[isign]);
        inTree[isign]->SetBranchAddress("pt_lead"          , &pt_lead[isign]);
        inTree[isign]->SetBranchAddress("pair_pt"          , &pair_pt[isign]);
        inTree[isign]->SetBranchAddress("pair_eta"     , &pair_eta[isign]);
        inTree[isign]->SetBranchAddress("pair_y"           , &pair_y[isign]);
        inTree[isign]->SetBranchAddress("asym"             , &asym[isign]);
        inTree[isign]->SetBranchAddress("dpt"           , &dpt[isign]);
        inTree[isign]->SetBranchAddress("deta"       , &deta[isign]);
        inTree[isign]->SetBranchAddress("etaavg"      , &etaavg[isign]);
        inTree[isign]->SetBranchAddress("phiavg"            , &phiavg[isign]);
        inTree[isign]->SetBranchAddress("dphi"     , &dphi[isign]);
        inTree[isign]->SetBranchAddress("dr"        , &dr[isign]);
        inTree[isign]->SetBranchAddress("minv"        , &minv[isign]);
        inTree[isign]->SetBranchAddress("m1.pt"           , &m1pt[isign]);
        inTree[isign]->SetBranchAddress("m2.pt"           , &m2pt[isign]);
        inTree[isign]->SetBranchAddress("m1.eta"       , &m1eta[isign]);
        inTree[isign]->SetBranchAddress("m2.eta"       , &m2eta[isign]);
        inTree[isign]->SetBranchAddress("m1.phi"        , &m1phi[isign]);
        inTree[isign]->SetBranchAddress("m2.phi"        , &m2phi[isign]);
        inTree[isign]->SetBranchAddress("m1.charge"           , &m1charge[isign]);
        inTree[isign]->SetBranchAddress("m2.charge"           , &m2charge[isign]);
        inTree[isign]->SetBranchAddress("passmu4mu4noL1"           , &passmu4mu4noL1[isign]);
        inTree[isign]->SetBranchAddress("pass2mu4"           , &pass2mu4[isign]);
        inTree[isign]->SetBranchAddress("mu1PassSingle"           , &mu1PassSingle[isign]);
        inTree[isign]->SetBranchAddress("mu2PassSingle"           , &mu2PassSingle[isign]);
        inTree[isign]->SetBranchAddress("passSeparated"           , &passSeparated[isign]);
        inTree[isign]->SetBranchAddress("passSeparatedDeta"           , &passSeparatedDeta[isign]);

        inTree[isign]->SetBranchStatus("*"                  ,0);//switch off all branches, then enable just the ones that we need
        inTree[isign]->SetBranchStatus("weight"           ,1);
        inTree[isign]->SetBranchStatus("pair_dPoverP"           ,1);
        inTree[isign]->SetBranchStatus("pt_lead"           ,1);
        inTree[isign]->SetBranchStatus("pair_pt"           ,1);
        inTree[isign]->SetBranchStatus("pair_y"           ,1);
        inTree[isign]->SetBranchStatus("asym"           ,1);
        inTree[isign]->SetBranchStatus("dpt"           ,1);
        inTree[isign]->SetBranchStatus("deta"           ,1);
        inTree[isign]->SetBranchStatus("etaavg"           ,1);
        inTree[isign]->SetBranchStatus("phiavg"           ,1);
        inTree[isign]->SetBranchStatus("dphi"           ,1);
        inTree[isign]->SetBranchStatus("dr"           ,1);
        inTree[isign]->SetBranchStatus("minv"           ,1);
        inTree[isign]->SetBranchStatus("m1.pt"           ,1);
        inTree[isign]->SetBranchStatus("m2.pt"           ,1);
        inTree[isign]->SetBranchStatus("m1.eta"           ,1);
        inTree[isign]->SetBranchStatus("m2.eta"           ,1);
        inTree[isign]->SetBranchStatus("m1.phi"           ,1);
        inTree[isign]->SetBranchStatus("m2.phi"           ,1);
        inTree[isign]->SetBranchStatus("m1.charge"           ,1);
        inTree[isign]->SetBranchStatus("m2.charge"           ,1);
        inTree[isign]->SetBranchStatus("passmu4mu4noL1"           ,1);
        inTree[isign]->SetBranchStatus("pass2mu4"           ,1);
        inTree[isign]->SetBranchStatus("mu1PassSingle"           ,1);
        inTree[isign]->SetBranchStatus("mu2PassSingle"           ,1);
        inTree[isign]->SetBranchStatus("passSeparated"           ,1);
        inTree[isign]->SetBranchStatus("passSeparatedDeta"           ,1);
    }
}

void MuonPairPlottingPP::InitOutput(){
    // function to define the output file
    // this is needed since the histograms, once defined, belong to the last file before them
    // need to create the output files after creating TFile objects for the input file and between defining TH1D, TH2D objects for the histograms
    if (isScram){
        outFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/histograms_scrambled_pairs_pp.root","recreate");
    }else if (isTight){
        outFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/histograms_real_pairs_pp_tight.root","recreate");
    }else{
        if (isRun3){
            outFile = new TFile(Form("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/histograms_real_pairs_pp_2024%s.root",trig_suffix.c_str()),"recreate");
        }else{
            outFile = new TFile(Form("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_run2/histograms_real_pairs_pp_2017%s.root",trig_suffix.c_str()),"recreate");
        }
    }
}

void MuonPairPlottingPP::InitHists(){
    int nDR_bins = 200;
    int nDR_zoomin_bins = 100;
    int nDphi_bins = 128;
    int neta_bins = 100;
    int nDeta_bins =  200;
    int npair_y_bins = 90;
    int npair_eta_bins = 100;
    int npt_asym_bins = 100;
    int npair_pt_ptlead_ratio_bins = 100;
    int nminv_bins_linear = 100;
    int npair_pT_bins_linear = 100;
    int npT_lead_bins_linear = 100;

    int neta_bins_trig_effcy = 48;
    int nDphi_bins_trig_effcy = 32;
    int nDeta_bins_trig_effcy = 40;
    int nDR_bins_trig_effcy = 40;
    int nminv_bins_trig_effcy = 40;
    int nDR_deta_dphi_zoomin_bins_trig_effcy = 20;

    std::vector<double> pT_bins_single_muon (pms.pT_bins_60); // make a copy of a suitable set of single-muon pT bins (adjustable) --> use the copy for histogram settings

    static const int nminv_bins_log = 40;
    float minv_logpow[ParamsSet::nSigns];
    minv_logpow[0] = 0.062;
    minv_logpow[1] = 0.045;
    float minv_max = 60;
    double minv_bins_log[ParamsSet::nSigns][nminv_bins_log+1];

    std::string dphi_regions[2] = {"near", "away"};

    for (unsigned int isign = 0; isign < ParamsSet::nSigns; isign++){
        for(int iminv = 0; iminv <= nminv_bins_log; iminv++){
            minv_bins_log[isign][iminv] = minv_max * pow(10.0, (static_cast<float>(iminv - nminv_bins_log)) * minv_logpow[isign]);
        }
    }


    h_crossx_dR_cut_vs_pair_pt_pair_eta = new TH2D("h_crossx_dR_cut_vs_pair_pt_pair_eta",";#eta^{pair};p_{T}^{pair} [GeV];#sigma^{truth}",npair_eta_bins_coarse,pair_eta_min,pair_eta_max,npair_pT_bins_coarse,pair_pt_min,pair_pt_max);
    h_crossx_dR_cut_vs_pair_pt_pair_eta->Sumw2();

    // ------------------- trigger efficiency histograms with single-b signal selection -------------------

    h_Deta_mu4_w_single_b_sig_sel = new TH1D("h_Deta_mu4_w_single_b_sig_sel",";#Delta#eta;",nDeta_bins_trig_effcy,-4.8,4.8);
    h_Deta_zoomin_mu4_w_single_b_sig_sel = new TH1D("h_Deta_zoomin_mu4_w_single_b_sig_sel",";#Delta#eta;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
    h_Dphi_mu4_w_single_b_sig_sel = new TH1D("h_Dphi_mu4_w_single_b_sig_sel",";#Delta#phi;",nDphi_bins_trig_effcy,-pms.PI,pms.PI);
    h_Dphi_zoomin_mu4_w_single_b_sig_sel = new TH1D("h_Dphi_zoomin_mu4_w_single_b_sig_sel",";#Delta#phi;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
    h_DR_mu4_w_single_b_sig_sel = new TH1D("h_DR_mu4_w_single_b_sig_sel",";#DeltaR;",nDR_bins_trig_effcy,0,pms.deltaR_thrsh[2]);
    h_DR_zoomin_mu4_w_single_b_sig_sel = new TH1D("h_DR_zoomin_mu4_w_single_b_sig_sel",";#DeltaR;",nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
    h_DR_zoomin_vs_pt2nd_mu4_w_single_b_sig_sel = new TH2D("h_DR_zoomin_vs_pt2nd_mu4_w_single_b_sig_sel",";p_{T,2nd} [GeV];#DeltaR",int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data(), nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
    h_DR_0_2_mu4_w_single_b_sig_sel = new TH1D("h_DR_0_2_mu4_w_single_b_sig_sel",";#DeltaR;",20,0,0.2);
    h_DR_0_2_vs_pt2nd_mu4_w_single_b_sig_sel = new TH2D("h_DR_0_2_vs_pt2nd_mu4_w_single_b_sig_sel",";p_{T,2nd} [GeV];#DeltaR",int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data(), 20,0,0.2);

    h_minv_zoomin_mu4_w_single_b_sig_sel = new TH1D("h_minv_zoomin_mu4_w_single_b_sig_sel",";m_{#mu#mu} [GeV];",nminv_bins_trig_effcy,0,3.);
    h_pair_pt_log_mu4_w_single_b_sig_sel = new TH1D("h_pair_pt_log_mu4_w_single_b_sig_sel",";p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
    h_pt2nd_mu4_w_single_b_sig_sel = new TH1D("h_pt2nd_mu4_w_single_b_sig_sel",";p_{T,2nd} [GeV];",int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
    h_pt2nd_vs_q_eta_2nd_mu4_w_single_b_sig_sel = new TH2D("h_pt2nd_vs_q_eta_2nd_mu4_w_single_b_sig_sel",";q*#eta;p_{T} [GeV]",neta_bins_trig_effcy,-2.4,2.4, int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
    h_pair_eta_vs_pair_pT_mu4_w_single_b_sig_sel = new TH2D("h_pair_eta_vs_pair_pT_mu4_w_single_b_sig_sel",";p_{T}^{pair} [GeV];#eta^{pair}",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),neta_bins_trig_effcy,-2.4,2.4);
    h_Deta_Dphi_mu4_w_single_b_sig_sel = new TH2D("h_Deta_Dphi_mu4_w_single_b_sig_sel",";#Delta#phi;#Delta#eta", nDphi_bins_trig_effcy,-pms.PI,pms.PI,nDeta_bins_trig_effcy,-4.8,4.8);
    h_eta1_eta2_mu4_w_single_b_sig_sel = new TH2D("h_eta1_eta2_mu4_w_single_b_sig_sel",";#eta_{sublead};#eta_{lead}",neta_bins_trig_effcy,-2.4,2.4, neta_bins_trig_effcy,-2.4,2.4);
    h_eta_avg_Deta_mu4_w_single_b_sig_sel = new TH2D("h_eta_avg_Deta_mu4_w_single_b_sig_sel",";#bar{#eta};#Delta#eta", neta_bins_trig_effcy,-2.4,2.4,nDeta_bins_trig_effcy,-4.8,4.8);
    h_eta_avg_Dphi_mu4_w_single_b_sig_sel = new TH2D("h_eta_avg_Dphi_mu4_w_single_b_sig_sel",";#Delta#phi;#bar{#eta}", nDphi_bins_trig_effcy,-pms.PI,pms.PI,neta_bins_trig_effcy,-2.4,2.4);
    h_minv_pair_pt_log_mu4_w_single_b_sig_sel = new TH2D("h_minv_pair_pt_log_mu4_w_single_b_sig_sel",";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),nminv_bins_log,minv_bins_log[1]);

    h_DR_zoomin_vs_pair_pT_mu4_w_single_b_sig_sel = new TH2D("h_DR_zoomin_vs_pair_pT_mu4_w_single_b_sig_sel",";p_{T}^{pair} [GeV];#DeltaR;", int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(), nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);

    h_Deta_mu4_mu4noL1_w_single_b_sig_sel = new TH1D("h_Deta_mu4_mu4noL1_w_single_b_sig_sel",";#Delta#eta;",nDeta_bins_trig_effcy,-4.8,4.8);
    h_Deta_zoomin_mu4_mu4noL1_w_single_b_sig_sel = new TH1D("h_Deta_zoomin_mu4_mu4noL1_w_single_b_sig_sel",";#Delta#eta;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
    h_Dphi_mu4_mu4noL1_w_single_b_sig_sel = new TH1D("h_Dphi_mu4_mu4noL1_w_single_b_sig_sel",";#Delta#phi;",nDphi_bins_trig_effcy,-pms.PI,pms.PI);
    h_Dphi_zoomin_mu4_mu4noL1_w_single_b_sig_sel = new TH1D("h_Dphi_zoomin_mu4_mu4noL1_w_single_b_sig_sel",";#Delta#phi;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
    h_DR_mu4_mu4noL1_w_single_b_sig_sel = new TH1D("h_DR_mu4_mu4noL1_w_single_b_sig_sel",";#DeltaR;",nDR_bins_trig_effcy,0,pms.deltaR_thrsh[2]);
    h_DR_zoomin_mu4_mu4noL1_w_single_b_sig_sel = new TH1D("h_DR_zoomin_mu4_mu4noL1_w_single_b_sig_sel",";#DeltaR;",nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
    h_DR_zoomin_vs_pt2nd_mu4_mu4noL1_w_single_b_sig_sel = new TH2D("h_DR_zoomin_vs_pt2nd_mu4_mu4noL1_w_single_b_sig_sel",";p_{T,2nd} [GeV];#DeltaR",int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data(), nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
    h_DR_0_2_mu4_mu4noL1_w_single_b_sig_sel = new TH1D("h_DR_0_2_mu4_mu4noL1_w_single_b_sig_sel",";#DeltaR;",20,0,0.2);
    h_DR_0_2_vs_pt2nd_mu4_mu4noL1_w_single_b_sig_sel = new TH2D("h_DR_0_2_vs_pt2nd_mu4_mu4noL1_w_single_b_sig_sel",";p_{T,2nd} [GeV];#DeltaR",int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data(), 20,0,0.2);

    h_minv_zoomin_mu4_mu4noL1_w_single_b_sig_sel = new TH1D("h_minv_zoomin_mu4_mu4noL1_w_single_b_sig_sel",";m_{#mu#mu} [GeV];",nminv_bins_trig_effcy,0,3.);
    h_pair_pt_log_mu4_mu4noL1_w_single_b_sig_sel = new TH1D("h_pair_pt_log_mu4_mu4noL1_w_single_b_sig_sel",";p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
    h_pt2nd_mu4_mu4noL1_w_single_b_sig_sel = new TH1D("h_pt2nd_mu4_mu4noL1_w_single_b_sig_sel",";p_{T,2nd} [GeV];",int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
    h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_w_single_b_sig_sel = new TH2D("h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_w_single_b_sig_sel",";q*#eta;p_{T} [GeV]",neta_bins_trig_effcy,-2.4,2.4, int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
    h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sepr_divided_w_single_b_sig_sel = new TH2D("h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sepr_divided_w_single_b_sig_sel",";q*#eta;p_{T} [GeV]",neta_bins_trig_effcy,-2.4,2.4, int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
    h_pair_eta_vs_pair_pT_mu4_mu4noL1_w_single_b_sig_sel = new TH2D("h_pair_eta_vs_pair_pT_mu4_mu4noL1_w_single_b_sig_sel",";p_{T}^{pair} [GeV];#eta^{pair}",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),neta_bins_trig_effcy,-2.4,2.4);
    h_Deta_Dphi_mu4_mu4noL1_w_single_b_sig_sel = new TH2D("h_Deta_Dphi_mu4_mu4noL1_w_single_b_sig_sel",";#Delta#phi;#Delta#eta", nDphi_bins_trig_effcy,-pms.PI,pms.PI,nDeta_bins_trig_effcy,-4.8,4.8);
    h_eta1_eta2_mu4_mu4noL1_w_single_b_sig_sel = new TH2D("h_eta1_eta2_mu4_mu4noL1_w_single_b_sig_sel",";#eta_{sublead};#eta_{lead}",neta_bins_trig_effcy,-2.4,2.4, neta_bins_trig_effcy,-2.4,2.4);
    h_eta_avg_Deta_mu4_mu4noL1_w_single_b_sig_sel = new TH2D("h_eta_avg_Deta_mu4_mu4noL1_w_single_b_sig_sel",";#bar{#eta};#Delta#eta", neta_bins_trig_effcy,-2.4,2.4,nDeta_bins_trig_effcy,-4.8,4.8);
    h_eta_avg_Dphi_mu4_mu4noL1_w_single_b_sig_sel = new TH2D("h_eta_avg_Dphi_mu4_mu4noL1_w_single_b_sig_sel",";#Delta#phi;#bar{#eta}", nDphi_bins_trig_effcy,-pms.PI,pms.PI,neta_bins_trig_effcy,-2.4,2.4);
    h_minv_pair_pt_log_mu4_mu4noL1_w_single_b_sig_sel = new TH2D("h_minv_pair_pt_log_mu4_mu4noL1_w_single_b_sig_sel",";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),nminv_bins_log,minv_bins_log[1]);

    h_Deta_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_w_single_b_sig_sel = new TH1D("h_Deta_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_w_single_b_sig_sel",";#Delta#eta;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
    h_Dphi_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_w_single_b_sig_sel = new TH1D("h_Dphi_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_w_single_b_sig_sel",";#Delta#phi;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
    h_DR_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_w_single_b_sig_sel = new TH1D("h_DR_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_w_single_b_sig_sel",";#DeltaR;",nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
    h_DR_0_2_mu4_mu4noL1_inv_w_by_single_mu_effcy_w_single_b_sig_sel = new TH1D("h_DR_0_2_mu4_mu4noL1_inv_w_by_single_mu_effcy_w_single_b_sig_sel",";#DeltaR;",20,0,0.2);
    h_minv_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_w_single_b_sig_sel = new TH1D("h_minv_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_w_single_b_sig_sel",";m_{#mu#mu} [GeV];",nminv_bins_trig_effcy,0,3.);
    h_pair_pt_log_mu4_mu4noL1_inv_w_by_single_mu_effcy_w_single_b_sig_sel = new TH1D("h_pair_pt_log_mu4_mu4noL1_inv_w_by_single_mu_effcy_w_single_b_sig_sel",";p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
    h_pt2nd_mu4_mu4noL1_inv_w_by_single_mu_effcy_w_single_b_sig_sel = new TH1D("h_pt2nd_mu4_mu4noL1_inv_w_by_single_mu_effcy_w_single_b_sig_sel",";p_{T,2nd} [GeV];",int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());

    h_Deta_mu4_mu4noL1_excl_w_single_b_sig_sel = new TH1D("h_Deta_mu4_mu4noL1_excl_w_single_b_sig_sel",";#Delta#eta;",nDeta_bins_trig_effcy,-4.8,4.8);
    h_Deta_zoomin_mu4_mu4noL1_excl_w_single_b_sig_sel = new TH1D("h_Deta_zoomin_mu4_mu4noL1_excl_w_single_b_sig_sel",";#Delta#eta;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
    h_Dphi_mu4_mu4noL1_excl_w_single_b_sig_sel = new TH1D("h_Dphi_mu4_mu4noL1_excl_w_single_b_sig_sel",";#Delta#phi;",nDphi_bins_trig_effcy,-pms.PI,pms.PI);
    h_Dphi_zoomin_mu4_mu4noL1_excl_w_single_b_sig_sel = new TH1D("h_Dphi_zoomin_mu4_mu4noL1_excl_w_single_b_sig_sel",";#Delta#phi;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
    h_DR_mu4_mu4noL1_excl_w_single_b_sig_sel = new TH1D("h_DR_mu4_mu4noL1_excl_w_single_b_sig_sel",";#DeltaR;",nDR_bins_trig_effcy,0,pms.deltaR_thrsh[2]);
    h_DR_zoomin_mu4_mu4noL1_excl_w_single_b_sig_sel = new TH1D("h_DR_zoomin_mu4_mu4noL1_excl_w_single_b_sig_sel",";#DeltaR;",nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
    h_minv_zoomin_mu4_mu4noL1_excl_w_single_b_sig_sel = new TH1D("h_minv_zoomin_mu4_mu4noL1_excl_w_single_b_sig_sel",";m_{#mu#mu} [GeV];",nminv_bins_trig_effcy,0,3.);
    h_pair_pt_log_mu4_mu4noL1_excl_w_single_b_sig_sel = new TH1D("h_pair_pt_log_mu4_mu4noL1_excl_w_single_b_sig_sel",";p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
    h_pt2nd_mu4_mu4noL1_excl_w_single_b_sig_sel = new TH1D("h_pt2nd_mu4_mu4noL1_excl_w_single_b_sig_sel",";p_{T,2nd} [GeV];",int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
    h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_excl_w_single_b_sig_sel = new TH2D("h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_excl_w_single_b_sig_sel",";q*#eta;p_{T} [GeV]",neta_bins_trig_effcy,-2.4,2.4, int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
    h_pair_eta_vs_pair_pT_mu4_mu4noL1_excl_w_single_b_sig_sel = new TH2D("h_pair_eta_vs_pair_pT_mu4_mu4noL1_excl_w_single_b_sig_sel",";p_{T}^{pair} [GeV];#eta^{pair}",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),neta_bins_trig_effcy,-2.4,2.4);
    h_Deta_Dphi_mu4_mu4noL1_excl_w_single_b_sig_sel = new TH2D("h_Deta_Dphi_mu4_mu4noL1_excl_w_single_b_sig_sel",";#Delta#phi;#Delta#eta", nDphi_bins_trig_effcy,-pms.PI,pms.PI,nDeta_bins_trig_effcy,-4.8,4.8);
    h_eta1_eta2_mu4_mu4noL1_excl_w_single_b_sig_sel = new TH2D("h_eta1_eta2_mu4_mu4noL1_excl_w_single_b_sig_sel",";#eta_{sublead};#eta_{lead}",neta_bins_trig_effcy,-2.4,2.4, neta_bins_trig_effcy,-2.4,2.4);
    h_eta_avg_Deta_mu4_mu4noL1_excl_w_single_b_sig_sel = new TH2D("h_eta_avg_Deta_mu4_mu4noL1_excl_w_single_b_sig_sel",";#bar{#eta};#Delta#eta", neta_bins_trig_effcy,-2.4,2.4,nDeta_bins_trig_effcy,-4.8,4.8);
    h_eta_avg_Dphi_mu4_mu4noL1_excl_w_single_b_sig_sel = new TH2D("h_eta_avg_Dphi_mu4_mu4noL1_excl_w_single_b_sig_sel",";#Delta#phi;#bar{#eta}", nDphi_bins_trig_effcy,-pms.PI,pms.PI,neta_bins_trig_effcy,-2.4,2.4);
    h_minv_pair_pt_log_mu4_mu4noL1_excl_w_single_b_sig_sel = new TH2D("h_minv_pair_pt_log_mu4_mu4noL1_excl_w_single_b_sig_sel",";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),nminv_bins_log,minv_bins_log[1]);

    h_Deta_2mu4_w_single_b_sig_sel = new TH1D("h_Deta_2mu4_w_single_b_sig_sel",";#Delta#eta;",nDeta_bins_trig_effcy,-4.8,4.8);
    h_Deta_zoomin_2mu4_w_single_b_sig_sel = new TH1D("h_Deta_zoomin_2mu4_w_single_b_sig_sel",";#Delta#eta;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
    h_Dphi_2mu4_w_single_b_sig_sel = new TH1D("h_Dphi_2mu4_w_single_b_sig_sel",";#Delta#phi;",nDphi_bins_trig_effcy,-pms.PI,pms.PI);
    h_Dphi_zoomin_2mu4_w_single_b_sig_sel = new TH1D("h_Dphi_zoomin_2mu4_w_single_b_sig_sel",";#Delta#phi;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
    h_DR_2mu4_w_single_b_sig_sel = new TH1D("h_DR_2mu4_w_single_b_sig_sel",";#DeltaR;",nDR_bins_trig_effcy,0,pms.deltaR_thrsh[2]);
    h_DR_zoomin_2mu4_w_single_b_sig_sel = new TH1D("h_DR_zoomin_2mu4_w_single_b_sig_sel",";#DeltaR;",nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
    h_DR_zoomin_vs_pt2nd_2mu4_w_single_b_sig_sel = new TH2D("h_DR_zoomin_vs_pt2nd_2mu4_w_single_b_sig_sel",";p_{T,2nd} [GeV];#DeltaR",int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data(), nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
    h_DR_0_2_2mu4_w_single_b_sig_sel = new TH1D("h_DR_0_2_2mu4_w_single_b_sig_sel",";#DeltaR;",20,0,0.2);
    h_DR_0_2_vs_pt2nd_2mu4_w_single_b_sig_sel = new TH2D("h_DR_0_2_vs_pt2nd_2mu4_w_single_b_sig_sel",";p_{T,2nd} [GeV];#DeltaR",int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data(), 20,0,0.2);

    h_minv_zoomin_2mu4_w_single_b_sig_sel = new TH1D("h_minv_zoomin_2mu4_w_single_b_sig_sel",";m_{#mu#mu} [GeV];",nminv_bins_trig_effcy,0,3.);
    h_pair_pt_log_2mu4_w_single_b_sig_sel = new TH1D("h_pair_pt_log_2mu4_w_single_b_sig_sel",";p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
    h_pt2nd_2mu4_w_single_b_sig_sel = new TH1D("h_pt2nd_2mu4_w_single_b_sig_sel",";p_{T,2nd} [GeV];",int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
    h_pt2nd_vs_q_eta_2nd_2mu4_w_single_b_sig_sel = new TH2D("h_pt2nd_vs_q_eta_2nd_2mu4_w_single_b_sig_sel",";q*#eta;p_{T} [GeV]",neta_bins_trig_effcy,-2.4,2.4, int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
    h_pt2nd_vs_q_eta_2nd_2mu4_sepr_divided_w_single_b_sig_sel = new TH2D("h_pt2nd_vs_q_eta_2nd_2mu4_sepr_divided_w_single_b_sig_sel",";q*#eta;p_{T} [GeV]",neta_bins_trig_effcy,-2.4,2.4, int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
    h_pair_eta_vs_pair_pT_2mu4_w_single_b_sig_sel = new TH2D("h_pair_eta_vs_pair_pT_2mu4_w_single_b_sig_sel",";p_{T}^{pair} [GeV];#eta^{pair}",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),neta_bins_trig_effcy,-2.4,2.4);
    h_Deta_Dphi_2mu4_w_single_b_sig_sel = new TH2D("h_Deta_Dphi_2mu4_w_single_b_sig_sel",";#Delta#phi;#Delta#eta", nDphi_bins_trig_effcy,-pms.PI,pms.PI,nDeta_bins_trig_effcy,-4.8,4.8);
    h_eta1_eta2_2mu4_w_single_b_sig_sel = new TH2D("h_eta1_eta2_2mu4_w_single_b_sig_sel",";#eta_{sublead};#eta_{lead}",neta_bins_trig_effcy,-2.4,2.4, neta_bins_trig_effcy,-2.4,2.4);
    h_eta_avg_Deta_2mu4_w_single_b_sig_sel = new TH2D("h_eta_avg_Deta_2mu4_w_single_b_sig_sel",";#bar{#eta};#Delta#eta", neta_bins_trig_effcy,-2.4,2.4,nDeta_bins_trig_effcy,-4.8,4.8);
    h_eta_avg_Dphi_2mu4_w_single_b_sig_sel = new TH2D("h_eta_avg_Dphi_2mu4_w_single_b_sig_sel",";#Delta#phi;#bar{#eta}", nDphi_bins_trig_effcy,-pms.PI,pms.PI,neta_bins_trig_effcy,-2.4,2.4);
    h_minv_pair_pt_log_2mu4_w_single_b_sig_sel = new TH2D("h_minv_pair_pt_log_2mu4_w_single_b_sig_sel",";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),nminv_bins_log,minv_bins_log[1]);
    
    h_Deta_zoomin_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel = new TH1D("h_Deta_zoomin_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel",";#Delta#eta;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
    h_Dphi_zoomin_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel = new TH1D("h_Dphi_zoomin_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel",";#Delta#phi;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
    h_DR_zoomin_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel = new TH1D("h_DR_zoomin_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel",";#DeltaR;",nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
    h_DR_0_2_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel = new TH1D("h_DR_0_2_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel",";#DeltaR;",20,0,0.2);
    h_minv_zoomin_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel = new TH1D("h_minv_zoomin_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel",";m_{#mu#mu} [GeV];",nminv_bins_trig_effcy,0,3.);
    h_pair_pt_log_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel = new TH1D("h_pair_pt_log_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel",";p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
    h_pt2nd_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel = new TH1D("h_pt2nd_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel",";p_{T,2nd} [GeV];",int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());

    for (unsigned int isign = 0; isign < ParamsSet::nSigns; isign++){
        
        // ------------------- trigger efficiency histograms full range (without signal selection) -------------------
        
        // pair kinematics requiring only single mu4 trigger
        h_Deta_mu4[isign] = new TH1D(Form("h_Deta_mu4_sign%d",isign+1),";#Delta#eta;",nDeta_bins_trig_effcy,-4.8,4.8);
        h_Deta_zoomin_mu4[isign] = new TH1D(Form("h_Deta_zoomin_mu4_sign%d",isign+1),";#Delta#eta;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
        h_Dphi_mu4[isign] = new TH1D(Form("h_Dphi_mu4_sign%d",isign+1),";#Delta#phi;",nDphi_bins_trig_effcy,-pms.PI,pms.PI);
        h_Dphi_zoomin_mu4[isign] = new TH1D(Form("h_Dphi_zoomin_mu4_sign%d",isign+1),";#Delta#phi;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
        h_DR_mu4[isign] = new TH1D(Form("h_DR_mu4_sign%d",isign+1),";#DeltaR;",nDR_bins_trig_effcy,0,pms.deltaR_thrsh[2]);
        h_DR_zoomin_mu4[isign] = new TH1D(Form("h_DR_zoomin_mu4_sign%d",isign+1),";#DeltaR;",nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
        h_DR_zoomin_vs_pt2nd_mu4[isign] = new TH2D(Form("h_DR_zoomin_vs_pt2nd_mu4_sign%d",isign+1),";p_{T,2nd} [GeV];#DeltaR",int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data(), nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
        h_DR_0_2_mu4[isign] = new TH1D(Form("h_DR_0_2_mu4_sign%d",isign+1),";#DeltaR;",20,0,0.2);
        h_DR_0_2_vs_pt2nd_mu4[isign] = new TH2D(Form("h_DR_0_2_vs_pt2nd_mu4_sign%d",isign+1),";p_{T,2nd} [GeV];#DeltaR",int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data(), 20,0,0.2);

        h_minv_zoomin_mu4[isign] = new TH1D(Form("h_minv_zoomin_mu4_sign%d",isign+1),";m_{#mu#mu} [GeV];",nminv_bins_trig_effcy,0,3.);
        h_pair_pt_log_mu4[isign] = new TH1D(Form("h_pair_pt_log_mu4_sign%d",isign+1),";p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
        h_pt2nd_mu4[isign] = new TH1D(Form("h_pt2nd_mu4_sign%d",isign+1),";p_{T,2nd} [GeV];",int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());

        h_pt2nd_vs_q_eta_2nd_mu4[isign] = new TH2D(Form("h_pt2nd_vs_q_eta_2nd_mu4_sign%d",isign+1),";q*#eta;p_{T} [GeV]",neta_bins_trig_effcy,-2.4,2.4, int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
        h_pt2nd_vs_phi2nd_mu4[isign] = new TH2D(Form("h_pt2nd_vs_phi2nd_mu4_sign%d",isign+1),";#phi;p_{T} [GeV]",nDphi_bins_trig_effcy,-pms.PI,pms.PI,int(pT_bins_single_muon.size())-1,pT_bins_single_muon.data());
        h_phi2nd_vs_q_eta_2nd_mu4[isign] =new TH2D(Form("h_phi2nd_vs_q_eta_2nd_mu4_sign%d",isign+1),";q*#eta;#phi",neta_bins_trig_effcy,-2.4,2.4,nDphi_bins_trig_effcy,-pms.PI,pms.PI);
        h_pair_eta_vs_pair_pT_mu4[isign] = new TH2D(Form("h_pair_eta_vs_pair_pT_mu4_sign%d",isign+1),";p_{T}^{pair} [GeV];#eta^{pair}",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),neta_bins_trig_effcy,-2.4,2.4);
        h_Deta_Dphi_mu4[isign] = new TH2D(Form("h_Deta_Dphi_mu4_sign%d",isign+1),";#Delta#phi;#Delta#eta", nDphi_bins_trig_effcy,-pms.PI,pms.PI,nDeta_bins_trig_effcy,-4.8,4.8);
        h_eta1_eta2_mu4[isign] = new TH2D(Form("h_eta1_eta2_mu4_sign%d",isign+1),";#eta_{sublead};#eta_{lead}",neta_bins_trig_effcy,-2.4,2.4, neta_bins_trig_effcy,-2.4,2.4);
        h_eta_avg_Deta_mu4[isign] = new TH2D(Form("h_eta_avg_Deta_mu4_sign%d",isign+1),";#bar{#eta};#Delta#eta", neta_bins_trig_effcy,-2.4,2.4,nDeta_bins_trig_effcy,-4.8,4.8);
        h_eta_avg_Dphi_mu4[isign] = new TH2D(Form("h_eta_avg_Dphi_mu4_sign%d",isign+1),";#Delta#phi;#bar{#eta}", nDphi_bins_trig_effcy,-pms.PI,pms.PI,neta_bins_trig_effcy,-2.4,2.4);
        h_minv_pair_pt_log_mu4[isign] = new TH2D(Form("h_minv_pair_pt_log_mu4_sign%d",isign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),nminv_bins_log,minv_bins_log[isign]);

        // pair kinematics requiring mu4_mu4noL1 trigger
        h_Deta_mu4_mu4noL1[isign] = new TH1D(Form("h_Deta_mu4_mu4noL1_sign%d",isign+1),";#Delta#eta;",nDeta_bins_trig_effcy,-4.8,4.8);
        h_Deta_zoomin_mu4_mu4noL1[isign] = new TH1D(Form("h_Deta_zoomin_mu4_mu4noL1_sign%d",isign+1),";#Delta#eta;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
        h_Dphi_mu4_mu4noL1[isign] = new TH1D(Form("h_Dphi_mu4_mu4noL1_sign%d",isign+1),";#Delta#phi;",nDphi_bins_trig_effcy,-pms.PI,pms.PI);
        h_Dphi_zoomin_mu4_mu4noL1[isign] = new TH1D(Form("h_Dphi_zoomin_mu4_mu4noL1_sign%d",isign+1),";#Delta#phi;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
        h_DR_mu4_mu4noL1[isign] = new TH1D(Form("h_DR_mu4_mu4noL1_sign%d",isign+1),";#DeltaR;",nDR_bins_trig_effcy,0,pms.deltaR_thrsh[2]);
        h_DR_zoomin_mu4_mu4noL1[isign] = new TH1D(Form("h_DR_zoomin_mu4_mu4noL1_sign%d",isign+1),";#DeltaR;",nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
        h_DR_zoomin_vs_pt2nd_mu4_mu4noL1[isign] = new TH2D(Form("h_DR_zoomin_vs_pt2nd_mu4_mu4noL1_sign%d",isign+1),";p_{T,2nd} [GeV];#DeltaR",int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data(), nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
        h_DR_0_2_mu4_mu4noL1[isign] = new TH1D(Form("h_DR_0_2_mu4_mu4noL1_sign%d",isign+1),";#DeltaR;",20,0,0.2);
        h_DR_0_2_vs_pt2nd_mu4_mu4noL1[isign] = new TH2D(Form("h_DR_0_2_vs_pt2nd_mu4_mu4noL1_sign%d",isign+1),";p_{T,2nd} [GeV];#DeltaR",int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data(), 20,0,0.2);

        h_minv_zoomin_mu4_mu4noL1[isign] = new TH1D(Form("h_minv_zoomin_mu4_mu4noL1_sign%d",isign+1),";m_{#mu#mu} [GeV];",nminv_bins_trig_effcy,0,3.);
        h_pair_pt_log_mu4_mu4noL1[isign] = new TH1D(Form("h_pair_pt_log_mu4_mu4noL1_sign%d",isign+1),";p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
        h_pt2nd_mu4_mu4noL1[isign] = new TH1D(Form("h_pt2nd_mu4_mu4noL1_sign%d",isign+1),";p_{T,2nd} [GeV];",int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());

        h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1[isign] = new TH2D(Form("h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sign%d",isign+1),";q*#eta;p_{T} [GeV]",neta_bins_trig_effcy,-2.4,2.4, int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
        h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sepr_divided[isign] = new TH2D(Form("h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sepr_divided_sign%d",isign+1),";q*#eta;p_{T} [GeV]",neta_bins_trig_effcy,-2.4,2.4, int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
        h_pt2nd_vs_phi2nd_mu4_mu4noL1[isign] = new TH2D(Form("h_pt2nd_vs_phi2nd_mu4_mu4noL1_sign%d",isign+1),";#phi;p_{T} [GeV]",nDphi_bins_trig_effcy,-pms.PI,pms.PI,int(pT_bins_single_muon.size())-1,pT_bins_single_muon.data());
        h_phi2nd_vs_q_eta_2nd_mu4_mu4noL1[isign] =new TH2D(Form("h_phi2nd_vs_q_eta_2nd_mu4_mu4noL1_sign%d",isign+1),";q*#eta;#phi",neta_bins_trig_effcy,-2.4,2.4,nDphi_bins_trig_effcy,-pms.PI,pms.PI);
        h_pair_eta_vs_pair_pT_mu4_mu4noL1[isign] = new TH2D(Form("h_pair_eta_vs_pair_pT_mu4_mu4noL1_sign%d",isign+1),";p_{T}^{pair} [GeV];#eta^{pair}",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),neta_bins_trig_effcy,-2.4,2.4);
        h_Deta_Dphi_mu4_mu4noL1[isign] = new TH2D(Form("h_Deta_Dphi_mu4_mu4noL1_sign%d",isign+1),";#Delta#phi;#Delta#eta", nDphi_bins_trig_effcy,-pms.PI,pms.PI,nDeta_bins_trig_effcy,-4.8,4.8);
        h_eta1_eta2_mu4_mu4noL1[isign] = new TH2D(Form("h_eta1_eta2_mu4_mu4noL1_sign%d",isign+1),";#eta_{sublead};#eta_{lead}",neta_bins_trig_effcy,-2.4,2.4, neta_bins_trig_effcy,-2.4,2.4);
        h_eta_avg_Deta_mu4_mu4noL1[isign] = new TH2D(Form("h_eta_avg_Deta_mu4_mu4noL1_sign%d",isign+1),";#bar{#eta};#Delta#eta", neta_bins_trig_effcy,-2.4,2.4,nDeta_bins_trig_effcy,-4.8,4.8);
        h_eta_avg_Dphi_mu4_mu4noL1[isign] = new TH2D(Form("h_eta_avg_Dphi_mu4_mu4noL1_sign%d",isign+1),";#Delta#phi;#bar{#eta}", nDphi_bins_trig_effcy,-pms.PI,pms.PI,neta_bins_trig_effcy,-2.4,2.4);
        h_minv_pair_pt_log_mu4_mu4noL1[isign] = new TH2D(Form("h_minv_pair_pt_log_mu4_mu4noL1_sign%d",isign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),nminv_bins_log,minv_bins_log[isign]);

        h_Deta_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy[isign] = new TH1D(Form("h_Deta_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_sign%d",isign+1),";#Delta#eta;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
        h_Dphi_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy[isign] = new TH1D(Form("h_Dphi_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_sign%d",isign+1),";#Delta#phi;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
        h_DR_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy[isign] = new TH1D(Form("h_DR_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_sign%d",isign+1),";#DeltaR;",nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
        h_DR_0_2_mu4_mu4noL1_inv_w_by_single_mu_effcy[isign] = new TH1D(Form("h_DR_0_2_mu4_mu4noL1_inv_w_by_single_mu_effcy_sign%d",isign+1),";#DeltaR;",20,0,0.2);
        h_minv_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy[isign] = new TH1D(Form("h_minv_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_sign%d",isign+1),";m_{#mu#mu} [GeV];",nminv_bins_trig_effcy,0,3.);
        h_pair_pt_log_mu4_mu4noL1_inv_w_by_single_mu_effcy[isign] = new TH1D(Form("h_pair_pt_log_mu4_mu4noL1_inv_w_by_single_mu_effcy_sign%d",isign+1),";p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
        h_pt2nd_mu4_mu4noL1_inv_w_by_single_mu_effcy[isign] = new TH1D(Form("h_pt2nd_mu4_mu4noL1_inv_w_by_single_mu_effcy_sign%d",isign+1),";p_{T,2nd} [GeV];",int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());

        h_Deta_mu4_mu4noL1_excl[isign] = new TH1D(Form("h_Deta_mu4_mu4noL1_excl_sign%d",isign+1),";#Delta#eta;",nDeta_bins_trig_effcy,-4.8,4.8);
        h_Deta_zoomin_mu4_mu4noL1_excl[isign] = new TH1D(Form("h_Deta_zoomin_mu4_mu4noL1_excl_sign%d",isign+1),";#Delta#eta;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
        h_Dphi_mu4_mu4noL1_excl[isign] = new TH1D(Form("h_Dphi_mu4_mu4noL1_excl_sign%d",isign+1),";#Delta#phi;",nDphi_bins_trig_effcy,-pms.PI,pms.PI);
        h_Dphi_zoomin_mu4_mu4noL1_excl[isign] = new TH1D(Form("h_Dphi_zoomin_mu4_mu4noL1_excl_sign%d",isign+1),";#Delta#phi;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
        h_DR_mu4_mu4noL1_excl[isign] = new TH1D(Form("h_DR_mu4_mu4noL1_excl_sign%d",isign+1),";#DeltaR;",nDR_bins_trig_effcy,0,pms.deltaR_thrsh[2]);
        h_DR_zoomin_mu4_mu4noL1_excl[isign] = new TH1D(Form("h_DR_zoomin_mu4_mu4noL1_excl_sign%d",isign+1),";#DeltaR;",nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
        h_minv_zoomin_mu4_mu4noL1_excl[isign] = new TH1D(Form("h_minv_zoomin_mu4_mu4noL1_excl_sign%d",isign+1),";m_{#mu#mu} [GeV];",nminv_bins_trig_effcy,0,3.);
        h_pair_pt_log_mu4_mu4noL1_excl[isign] = new TH1D(Form("h_pair_pt_log_mu4_mu4noL1_excl_sign%d",isign+1),";p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
        h_pt2nd_mu4_mu4noL1_excl[isign] = new TH1D(Form("h_pt2nd_mu4_mu4noL1_excl_sign%d",isign+1),";p_{T,2nd} [GeV];",int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());

        h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_excl[isign] = new TH2D(Form("h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_excl_sign%d",isign+1),";q*#eta;p_{T} [GeV]",neta_bins_trig_effcy,-2.4,2.4, int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
        h_pair_eta_vs_pair_pT_mu4_mu4noL1_excl[isign] = new TH2D(Form("h_pair_eta_vs_pair_pT_mu4_mu4noL1_excl_sign%d",isign+1),";p_{T}^{pair} [GeV];#eta^{pair}",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),neta_bins_trig_effcy,-2.4,2.4);
        h_Deta_Dphi_mu4_mu4noL1_excl[isign] = new TH2D(Form("h_Deta_Dphi_mu4_mu4noL1_excl_sign%d",isign+1),";#Delta#phi;#Delta#eta", nDphi_bins_trig_effcy,-pms.PI,pms.PI,nDeta_bins_trig_effcy,-4.8,4.8);
        h_eta1_eta2_mu4_mu4noL1_excl[isign] = new TH2D(Form("h_eta1_eta2_mu4_mu4noL1_excl_sign%d",isign+1),";#eta_{sublead};#eta_{lead}",neta_bins_trig_effcy,-2.4,2.4, neta_bins_trig_effcy,-2.4,2.4);
        h_eta_avg_Deta_mu4_mu4noL1_excl[isign] = new TH2D(Form("h_eta_avg_Deta_mu4_mu4noL1_excl_sign%d",isign+1),";#bar{#eta};#Delta#eta", neta_bins_trig_effcy,-2.4,2.4,nDeta_bins_trig_effcy,-4.8,4.8);
        h_eta_avg_Dphi_mu4_mu4noL1_excl[isign] = new TH2D(Form("h_eta_avg_Dphi_mu4_mu4noL1_excl_sign%d",isign+1),";#Delta#phi;#bar{#eta}", nDphi_bins_trig_effcy,-pms.PI,pms.PI,neta_bins_trig_effcy,-2.4,2.4);
        h_minv_pair_pt_log_mu4_mu4noL1_excl[isign] = new TH2D(Form("h_minv_pair_pt_log_mu4_mu4noL1_excl_sign%d",isign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),nminv_bins_log,minv_bins_log[isign]);

        // pair kinematics requiring 2mu4 trigger
        h_Deta_2mu4[isign] = new TH1D(Form("h_Deta_2mu4_sign%d",isign+1),";#Delta#eta;",nDeta_bins_trig_effcy,-4.8,4.8);
        h_Deta_zoomin_2mu4[isign] = new TH1D(Form("h_Deta_zoomin_2mu4_sign%d",isign+1),";#Delta#eta;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
        h_Dphi_2mu4[isign] = new TH1D(Form("h_Dphi_2mu4_sign%d",isign+1),";#Delta#phi;",nDphi_bins_trig_effcy,-pms.PI,pms.PI);
        h_Dphi_zoomin_2mu4[isign] = new TH1D(Form("h_Dphi_zoomin_2mu4_sign%d",isign+1),";#Delta#phi;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
        h_DR_2mu4[isign] = new TH1D(Form("h_DR_2mu4_sign%d",isign+1),";#DeltaR;",nDR_bins_trig_effcy,0,pms.deltaR_thrsh[2]);
        h_DR_zoomin_2mu4[isign] = new TH1D(Form("h_DR_zoomin_2mu4_sign%d",isign+1),";#DeltaR;",nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
        h_DR_zoomin_vs_pt2nd_2mu4[isign] = new TH2D(Form("h_DR_zoomin_vs_pt2nd_2mu4_sign%d",isign+1),";p_{T,2nd} [GeV];#DeltaR",int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data(), nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
        h_DR_0_2_2mu4[isign] = new TH1D(Form("h_DR_0_2_2mu4_sign%d",isign+1),";#DeltaR;",20,0,0.2);
        h_DR_0_2_vs_pt2nd_2mu4[isign] = new TH2D(Form("h_DR_0_2_vs_pt2nd_2mu4_sign%d",isign+1),";p_{T,2nd} [GeV];#DeltaR",int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data(), 20,0,0.2);

        h_minv_zoomin_2mu4[isign] = new TH1D(Form("h_minv_zoomin_2mu4_sign%d",isign+1),";m_{#mu#mu} [GeV];",nminv_bins_trig_effcy,0,3.);
        h_pair_pt_log_2mu4[isign] = new TH1D(Form("h_pair_pt_log_2mu4_sign%d",isign+1),";p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
        h_pt2nd_2mu4[isign] = new TH1D(Form("h_pt2nd_2mu4_sign%d",isign+1),";p_{T,2nd} [GeV];",int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());

        h_pt2nd_vs_q_eta_2nd_2mu4[isign] = new TH2D(Form("h_pt2nd_vs_q_eta_2nd_2mu4_sign%d",isign+1),";q*#eta;p_{T} [GeV]",neta_bins_trig_effcy,-2.4,2.4, int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
        h_pt2nd_vs_q_eta_2nd_2mu4_sepr_divided[isign] = new TH2D(Form("h_pt2nd_vs_q_eta_2nd_2mu4_sepr_divided_sign%d",isign+1),";q*#eta;p_{T} [GeV]",neta_bins_trig_effcy,-2.4,2.4, int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
        h_pt2nd_vs_phi2nd_2mu4[isign] = new TH2D(Form("h_pt2nd_vs_phi2nd_2mu4_sign%d",isign+1),";#phi;p_{T} [GeV]",nDphi_bins_trig_effcy,-pms.PI,pms.PI,int(pT_bins_single_muon.size())-1,pT_bins_single_muon.data());
        h_phi2nd_vs_q_eta_2nd_2mu4[isign] =new TH2D(Form("h_phi2nd_vs_q_eta_2nd_2mu4_sign%d",isign+1),";q*#eta;#phi",neta_bins_trig_effcy,-2.4,2.4,nDphi_bins_trig_effcy,-pms.PI,pms.PI);
        h_pair_eta_vs_pair_pT_2mu4[isign] = new TH2D(Form("h_pair_eta_vs_pair_pT_2mu4_sign%d",isign+1),";p_{T}^{pair} [GeV];#eta^{pair}",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),neta_bins_trig_effcy,-2.4,2.4);
        h_Deta_Dphi_2mu4[isign] = new TH2D(Form("h_Deta_Dphi_2mu4_sign%d",isign+1),";#Delta#phi;#Delta#eta", nDphi_bins_trig_effcy,-pms.PI,pms.PI,nDeta_bins_trig_effcy,-4.8,4.8);
        h_eta1_eta2_2mu4[isign] = new TH2D(Form("h_eta1_eta2_2mu4_sign%d",isign+1),";#eta_{sublead};#eta_{lead}",neta_bins_trig_effcy,-2.4,2.4, neta_bins_trig_effcy,-2.4,2.4);
        h_eta_avg_Deta_2mu4[isign] = new TH2D(Form("h_eta_avg_Deta_2mu4_sign%d",isign+1),";#bar{#eta};#Delta#eta", neta_bins_trig_effcy,-2.4,2.4,nDeta_bins_trig_effcy,-4.8,4.8);
        h_eta_avg_Dphi_2mu4[isign] = new TH2D(Form("h_eta_avg_Dphi_2mu4_sign%d",isign+1),";#Delta#phi;#bar{#eta}", nDphi_bins_trig_effcy,-pms.PI,pms.PI,neta_bins_trig_effcy,-2.4,2.4);
        h_minv_pair_pt_log_2mu4[isign] = new TH2D(Form("h_minv_pair_pt_log_2mu4_sign%d",isign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),nminv_bins_log,minv_bins_log[isign]);

        h_Deta_zoomin_2mu4_inv_w_by_single_mu_effcy[isign] = new TH1D(Form("h_Deta_zoomin_2mu4_inv_w_by_single_mu_effcy_sign%d",isign+1),";#Delta#eta;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
        h_Dphi_zoomin_2mu4_inv_w_by_single_mu_effcy[isign] = new TH1D(Form("h_Dphi_zoomin_2mu4_inv_w_by_single_mu_effcy_sign%d",isign+1),";#Delta#phi;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
        h_DR_zoomin_2mu4_inv_w_by_single_mu_effcy[isign] = new TH1D(Form("h_DR_zoomin_2mu4_inv_w_by_single_mu_effcy_sign%d",isign+1),";#DeltaR;",nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
        h_DR_0_2_2mu4_inv_w_by_single_mu_effcy[isign] = new TH1D(Form("h_DR_0_2_2mu4_inv_w_by_single_mu_effcy_sign%d",isign+1),";#DeltaR;",20,0,0.2);
        h_minv_zoomin_2mu4_inv_w_by_single_mu_effcy[isign] = new TH1D(Form("h_minv_zoomin_2mu4_inv_w_by_single_mu_effcy_sign%d",isign+1),";m_{#mu#mu} [GeV];",nminv_bins_trig_effcy,0,3.);
        h_pair_pt_log_2mu4_inv_w_by_single_mu_effcy[isign] = new TH1D(Form("h_pair_pt_log_2mu4_inv_w_by_single_mu_effcy_sign%d",isign+1),";p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
        h_pt2nd_2mu4_inv_w_by_single_mu_effcy[isign] = new TH1D(Form("h_pt2nd_2mu4_inv_w_by_single_mu_effcy_sign%d",isign+1),";p_{T,2nd} [GeV];",int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
        // ------------------- trigger efficiency histograms in regions with good single-muon acceptance (without signal selection) -------------------
        
        // pair kinematics requiring only single mu4 trigger
        h_Deta_mu4_good_accept[isign] = new TH1D(Form("h_Deta_mu4_good_accept_sign%d",isign+1),";#Delta#eta;",nDeta_bins_trig_effcy,-4.8,4.8);
        h_Deta_zoomin_mu4_good_accept[isign] = new TH1D(Form("h_Deta_zoomin_mu4_good_accept_sign%d",isign+1),";#Delta#eta;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
        h_Dphi_mu4_good_accept[isign] = new TH1D(Form("h_Dphi_mu4_good_accept_sign%d",isign+1),";#Delta#phi;",nDphi_bins_trig_effcy,-pms.PI,pms.PI);
        h_Dphi_zoomin_mu4_good_accept[isign] = new TH1D(Form("h_Dphi_zoomin_mu4_good_accept_sign%d",isign+1),";#Delta#phi;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
        h_DR_mu4_good_accept[isign] = new TH1D(Form("h_DR_mu4_good_accept_sign%d",isign+1),";#DeltaR;",nDR_bins_trig_effcy,0,pms.deltaR_thrsh[2]);
        h_DR_zoomin_mu4_good_accept[isign] = new TH1D(Form("h_DR_zoomin_mu4_good_accept_sign%d",isign+1),";#DeltaR;",nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
        h_minv_zoomin_mu4_good_accept[isign] = new TH1D(Form("h_minv_zoomin_mu4_good_accept_sign%d",isign+1),";m_{#mu#mu} [GeV];",nminv_bins_trig_effcy,0,3.);
        h_pair_pt_log_mu4_good_accept[isign] = new TH1D(Form("h_pair_pt_log_mu4_good_accept_sign%d",isign+1),";p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());

        // pair kinematics requiring mu4_mu4noL1 trigger
        h_Deta_mu4_mu4noL1_good_accept[isign] = new TH1D(Form("h_Deta_mu4_mu4noL1_good_accept_sign%d",isign+1),";#Delta#eta;",nDeta_bins_trig_effcy,-4.8,4.8);
        h_Deta_zoomin_mu4_mu4noL1_good_accept[isign] = new TH1D(Form("h_Deta_zoomin_mu4_mu4noL1_good_accept_sign%d",isign+1),";#Delta#eta;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
        h_Dphi_mu4_mu4noL1_good_accept[isign] = new TH1D(Form("h_Dphi_mu4_mu4noL1_good_accept_sign%d",isign+1),";#Delta#phi;",nDphi_bins_trig_effcy,-pms.PI,pms.PI);
        h_Dphi_zoomin_mu4_mu4noL1_good_accept[isign] = new TH1D(Form("h_Dphi_zoomin_mu4_mu4noL1_good_accept_sign%d",isign+1),";#Delta#phi;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
        h_DR_mu4_mu4noL1_good_accept[isign] = new TH1D(Form("h_DR_mu4_mu4noL1_good_accept_sign%d",isign+1),";#DeltaR;",nDR_bins_trig_effcy,0,pms.deltaR_thrsh[2]);
        h_DR_zoomin_mu4_mu4noL1_good_accept[isign] = new TH1D(Form("h_DR_zoomin_mu4_mu4noL1_good_accept_sign%d",isign+1),";#DeltaR;",nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
        h_minv_zoomin_mu4_mu4noL1_good_accept[isign] = new TH1D(Form("h_minv_zoomin_mu4_mu4noL1_good_accept_sign%d",isign+1),";m_{#mu#mu} [GeV];",nminv_bins_trig_effcy,0,3.);
        h_pair_pt_log_mu4_mu4noL1_good_accept[isign] = new TH1D(Form("h_pair_pt_log_mu4_mu4noL1_good_accept_sign%d",isign+1),";p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());

        // pair kinematics requiring 2mu4 trigger
        h_Deta_2mu4_good_accept[isign] = new TH1D(Form("h_Deta_2mu4_good_accept_sign%d",isign+1),";#Delta#eta;",nDeta_bins_trig_effcy,-4.8,4.8);
        h_Deta_zoomin_2mu4_good_accept[isign] = new TH1D(Form("h_Deta_zoomin_2mu4_good_accept_sign%d",isign+1),";#Delta#eta;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
        h_Dphi_2mu4_good_accept[isign] = new TH1D(Form("h_Dphi_2mu4_good_accept_sign%d",isign+1),";#Delta#phi;",nDphi_bins_trig_effcy,-pms.PI,pms.PI);
        h_Dphi_zoomin_2mu4_good_accept[isign] = new TH1D(Form("h_Dphi_zoomin_2mu4_good_accept_sign%d",isign+1),";#Delta#phi;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
        h_DR_2mu4_good_accept[isign] = new TH1D(Form("h_DR_2mu4_good_accept_sign%d",isign+1),";#DeltaR;",nDR_bins_trig_effcy,0,pms.deltaR_thrsh[2]);
        h_DR_zoomin_2mu4_good_accept[isign] = new TH1D(Form("h_DR_zoomin_2mu4_good_accept_sign%d",isign+1),";#DeltaR;",nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
        h_minv_zoomin_2mu4_good_accept[isign] = new TH1D(Form("h_minv_zoomin_2mu4_good_accept_sign%d",isign+1),";m_{#mu#mu} [GeV];",nminv_bins_trig_effcy,0,3.);
        h_pair_pt_log_2mu4_good_accept[isign] = new TH1D(Form("h_pair_pt_log_2mu4_good_accept_sign%d",isign+1),";p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());

        // ------------------- trigger efficiency vs pair pT or pT1st histograms (full range) -------------------
        
        // pair kinematics requiring only single mu4 trigger
        h_Deta_vs_pT_1st_mu4[isign] = new TH2D(Form("h_Deta_vs_pT_1st_mu4_sign%d",isign+1),";p_{T,1st} [GeV];#Delta#eta;",  int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data(),  nDeta_bins_trig_effcy,-4.8,4.8);
        h_Deta_zoomin_vs_pT_1st_mu4[isign] = new TH2D(Form("h_Deta_zoomin_vs_pT_1st_mu4_sign%d",isign+1),";p_{T,1st} [GeV];#Delta#eta;",  int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data(),  nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
        h_Dphi_vs_pT_1st_mu4[isign] = new TH2D(Form("h_Dphi_vs_pT_1st_mu4_sign%d",isign+1),";p_{T,1st} [GeV];#Delta#phi;",  int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data(),  nDphi_bins_trig_effcy,-pms.PI,pms.PI);
        h_Dphi_zoomin_vs_pT_1st_mu4[isign] = new TH2D(Form("h_Dphi_zoomin_vs_pT_1st_mu4_sign%d",isign+1),";p_{T,1st} [GeV];#Delta#phi;",  int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data(),  nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
        h_DR_vs_pT_1st_mu4[isign] = new TH2D(Form("h_DR_vs_pT_1st_mu4_sign%d",isign+1),";p_{T,1st} [GeV];#DeltaR;",  int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data(),  nDR_bins_trig_effcy,0,pms.deltaR_thrsh[2]);
        h_DR_zoomin_vs_pT_1st_mu4[isign] = new TH2D(Form("h_DR_zoomin_vs_pT_1st_mu4_sign%d",isign+1),";p_{T,1st} [GeV];#DeltaR;",  int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data(),  nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
        h_minv_zoomin_vs_pT_1st_mu4[isign] = new TH2D(Form("h_minv_zoomin_vs_pT_1st_mu4_sign%d",isign+1),";p_{T,1st} [GeV];m_{#mu#mu} [GeV];",  int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data(),  nminv_bins_trig_effcy,0,3.);
        h_pair_pt_log_vs_pT_1st_mu4[isign] = new TH2D(Form("h_pair_pt_log_vs_pT_1st_mu4_sign%d",isign+1),";p_{T,1st} [GeV];p_{T}^{pair} [GeV];",  int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data(),  int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
        h_pt2nd_vs_pT_1st_mu4[isign] = new TH2D(Form("h_pt2nd_vs_pT_1st_mu4_sign%d",isign+1),";p_{T,1st} [GeV];p_{T,2nd} [GeV];",  int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data(),  int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());

        h_Deta_vs_pair_pT_mu4[isign] = new TH2D(Form("h_Deta_vs_pair_pT_mu4_sign%d",isign+1),";p_{T}^{pair} [GeV];#Delta#eta;", int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(), nDeta_bins_trig_effcy,-4.8,4.8);
        h_Deta_zoomin_vs_pair_pT_mu4[isign] = new TH2D(Form("h_Deta_zoomin_vs_pair_pT_mu4_sign%d",isign+1),";p_{T}^{pair} [GeV];#Delta#eta;", int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(), nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
        h_Dphi_vs_pair_pT_mu4[isign] = new TH2D(Form("h_Dphi_vs_pair_pT_mu4_sign%d",isign+1),";p_{T}^{pair} [GeV];#Delta#phi;", int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(), nDphi_bins_trig_effcy,-pms.PI,pms.PI);
        h_Dphi_zoomin_vs_pair_pT_mu4[isign] = new TH2D(Form("h_Dphi_zoomin_vs_pair_pT_mu4_sign%d",isign+1),";p_{T}^{pair} [GeV];#Delta#phi;", int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(), nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
        h_DR_vs_pair_pT_mu4[isign] = new TH2D(Form("h_DR_vs_pair_pT_mu4_sign%d",isign+1),";p_{T}^{pair} [GeV];#DeltaR;", int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(), nDR_bins_trig_effcy,0,pms.deltaR_thrsh[2]);
        h_DR_zoomin_vs_pair_pT_mu4[isign] = new TH2D(Form("h_DR_zoomin_vs_pair_pT_mu4_sign%d",isign+1),";p_{T}^{pair} [GeV];#DeltaR;", int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(), nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
        h_minv_zoomin_vs_pair_pT_mu4[isign] = new TH2D(Form("h_minv_zoomin_vs_pair_pT_mu4_sign%d",isign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV];", int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(), nminv_bins_trig_effcy,0,3.);
        h_pair_pt_log_vs_pair_pT_mu4[isign] = new TH2D(Form("h_pair_pt_log_vs_pair_pT_mu4_sign%d",isign+1),";p_{T}^{pair} [GeV];p_{T}^{pair} [GeV];", int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(), int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
        h_pt2nd_vs_pair_pT_mu4[isign] = new TH2D(Form("h_pt2nd_vs_pair_pT_mu4_sign%d",isign+1),";p_{T}^{pair} [GeV];p_{T,2nd} [GeV];", int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(), int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());

        // pair kinematics requiring mu4_mu4noL1 trigger
        h_Deta_vs_pT_1st_mu4_mu4noL1[isign] = new TH2D(Form("h_Deta_vs_pT_1st_mu4_mu4noL1_sign%d",isign+1),";p_{T,1st} [GeV];#Delta#eta;",nDeta_bins_trig_effcy,-4.8,4.8, int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
        h_Deta_zoomin_vs_pT_1st_mu4_mu4noL1[isign] = new TH2D(Form("h_Deta_zoomin_vs_pT_1st_mu4_mu4noL1_sign%d",isign+1),";p_{T,1st} [GeV];#Delta#eta;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8, int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
        h_Dphi_vs_pT_1st_mu4_mu4noL1[isign] = new TH2D(Form("h_Dphi_vs_pT_1st_mu4_mu4noL1_sign%d",isign+1),";p_{T,1st} [GeV];#Delta#phi;",nDphi_bins_trig_effcy,-pms.PI,pms.PI, int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
        h_Dphi_zoomin_vs_pT_1st_mu4_mu4noL1[isign] = new TH2D(Form("h_Dphi_zoomin_vs_pT_1st_mu4_mu4noL1_sign%d",isign+1),";p_{T,1st} [GeV];#Delta#phi;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8, int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
        h_DR_vs_pT_1st_mu4_mu4noL1[isign] = new TH2D(Form("h_DR_vs_pT_1st_mu4_mu4noL1_sign%d",isign+1),";p_{T,1st} [GeV];#DeltaR;",nDR_bins_trig_effcy,0,pms.deltaR_thrsh[2], int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
        h_DR_zoomin_vs_pT_1st_mu4_mu4noL1[isign] = new TH2D(Form("h_DR_zoomin_vs_pT_1st_mu4_mu4noL1_sign%d",isign+1),";p_{T,1st} [GeV];#DeltaR;",nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8, int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
        h_minv_zoomin_vs_pT_1st_mu4_mu4noL1[isign] = new TH2D(Form("h_minv_zoomin_vs_pT_1st_mu4_mu4noL1_sign%d",isign+1),";p_{T,1st} [GeV];m_{#mu#mu} [GeV];",nminv_bins_trig_effcy,0,3., int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
        h_pair_pt_log_vs_pT_1st_mu4_mu4noL1[isign] = new TH2D(Form("h_pair_pt_log_vs_pT_1st_mu4_mu4noL1_sign%d",isign+1),";p_{T,1st} [GeV];p_{T}^{pair} [GeV];",int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data(), int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
        h_pt2nd_vs_pT_1st_mu4_mu4noL1[isign] = new TH2D(Form("h_pt2nd_vs_pT_1st_mu4_mu4noL1_sign%d",isign+1),";p_{T,1st} [GeV];p_{T,2nd} [GeV];",int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data(), int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());

        h_Deta_vs_pair_pT_mu4_mu4noL1[isign] = new TH2D(Form("h_Deta_vs_pair_pT_mu4_mu4noL1_sign%d",isign+1),";p_{T}^{pair} [GeV];#Delta#eta;",nDeta_bins_trig_effcy,-4.8,4.8, int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
        h_Deta_zoomin_vs_pair_pT_mu4_mu4noL1[isign] = new TH2D(Form("h_Deta_zoomin_vs_pair_pT_mu4_mu4noL1_sign%d",isign+1),";p_{T}^{pair} [GeV];#Delta#eta;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8, int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
        h_Dphi_vs_pair_pT_mu4_mu4noL1[isign] = new TH2D(Form("h_Dphi_vs_pair_pT_mu4_mu4noL1_sign%d",isign+1),";p_{T}^{pair} [GeV];#Delta#phi;",nDphi_bins_trig_effcy,-pms.PI,pms.PI, int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
        h_Dphi_zoomin_vs_pair_pT_mu4_mu4noL1[isign] = new TH2D(Form("h_Dphi_zoomin_vs_pair_pT_mu4_mu4noL1_sign%d",isign+1),";p_{T}^{pair} [GeV];#Delta#phi;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8, int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
        h_DR_vs_pair_pT_mu4_mu4noL1[isign] = new TH2D(Form("h_DR_vs_pair_pT_mu4_mu4noL1_sign%d",isign+1),";p_{T}^{pair} [GeV];#DeltaR;",nDR_bins_trig_effcy,0,pms.deltaR_thrsh[2], int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
        h_DR_zoomin_vs_pair_pT_mu4_mu4noL1[isign] = new TH2D(Form("h_DR_zoomin_vs_pair_pT_mu4_mu4noL1_sign%d",isign+1),";p_{T}^{pair} [GeV];#DeltaR;",nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8, int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
        h_minv_zoomin_vs_pair_pT_mu4_mu4noL1[isign] = new TH2D(Form("h_minv_zoomin_vs_pair_pT_mu4_mu4noL1_sign%d",isign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV];",nminv_bins_trig_effcy,0,3., int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
        h_pair_pt_log_vs_pair_pT_mu4_mu4noL1[isign] = new TH2D(Form("h_pair_pt_log_vs_pair_pT_mu4_mu4noL1_sign%d",isign+1),";p_{T}^{pair} [GeV];p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(), int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
        h_pt2nd_vs_pair_pT_mu4_mu4noL1[isign] = new TH2D(Form("h_pt2nd_vs_pair_pT_mu4_mu4noL1_sign%d",isign+1),";p_{T}^{pair} [GeV];p_{T,2nd} [GeV];",int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data(), int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());

        h_Deta_vs_pT_1st_mu4_mu4noL1_excl[isign] = new TH2D(Form("h_Deta_vs_pT_1st_mu4_mu4noL1_excl_sign%d",isign+1),";p_{T,1st} [GeV];#Delta#eta;",nDeta_bins_trig_effcy,-4.8,4.8,  int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
        h_Deta_zoomin_vs_pT_1st_mu4_mu4noL1_excl[isign] = new TH2D(Form("h_Deta_zoomin_vs_pT_1st_mu4_mu4noL1_excl_sign%d",isign+1),";p_{T,1st} [GeV];#Delta#eta;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8,  int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
        h_Dphi_vs_pT_1st_mu4_mu4noL1_excl[isign] = new TH2D(Form("h_Dphi_vs_pT_1st_mu4_mu4noL1_excl_sign%d",isign+1),";p_{T,1st} [GeV];#Delta#phi;",nDphi_bins_trig_effcy,-pms.PI,pms.PI,  int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
        h_Dphi_zoomin_vs_pT_1st_mu4_mu4noL1_excl[isign] = new TH2D(Form("h_Dphi_zoomin_vs_pT_1st_mu4_mu4noL1_excl_sign%d",isign+1),";p_{T,1st} [GeV];#Delta#phi;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8,  int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
        h_DR_vs_pT_1st_mu4_mu4noL1_excl[isign] = new TH2D(Form("h_DR_vs_pT_1st_mu4_mu4noL1_excl_sign%d",isign+1),";p_{T,1st} [GeV];#DeltaR;",nDR_bins_trig_effcy,0,pms.deltaR_thrsh[2],  int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
        h_DR_zoomin_vs_pT_1st_mu4_mu4noL1_excl[isign] = new TH2D(Form("h_DR_zoomin_vs_pT_1st_mu4_mu4noL1_excl_sign%d",isign+1),";p_{T,1st} [GeV];#DeltaR;",nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8,  int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
        h_minv_zoomin_vs_pT_1st_mu4_mu4noL1_excl[isign] = new TH2D(Form("h_minv_zoomin_vs_pT_1st_mu4_mu4noL1_excl_sign%d",isign+1),";p_{T,1st} [GeV];m_{#mu#mu} [GeV];",nminv_bins_trig_effcy,0,3.,  int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
        h_pair_pt_log_vs_pT_1st_mu4_mu4noL1_excl[isign] = new TH2D(Form("h_pair_pt_log_vs_pT_1st_mu4_mu4noL1_excl_sign%d",isign+1),";p_{T,1st} [GeV];p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),  int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
        h_pt2nd_vs_pT_1st_mu4_mu4noL1_excl[isign] = new TH2D(Form("h_pt2nd_vs_pT_1st_mu4_mu4noL1_excl_sign%d",isign+1),";p_{T,1st} [GeV];p_{T,2nd} [GeV];",int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data(),  int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());

        h_Deta_vs_pair_pT_mu4_mu4noL1_excl[isign] = new TH2D(Form("h_Deta_vs_pair_pT_mu4_mu4noL1_excl_sign%d",isign+1),";p_{T}^{pair} [GeV];#Delta#eta;",nDeta_bins_trig_effcy,-4.8,4.8,  int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
        h_Deta_zoomin_vs_pair_pT_mu4_mu4noL1_excl[isign] = new TH2D(Form("h_Deta_zoomin_vs_pair_pT_mu4_mu4noL1_excl_sign%d",isign+1),";p_{T}^{pair} [GeV];#Delta#eta;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8,  int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
        h_Dphi_vs_pair_pT_mu4_mu4noL1_excl[isign] = new TH2D(Form("h_Dphi_vs_pair_pT_mu4_mu4noL1_excl_sign%d",isign+1),";p_{T}^{pair} [GeV];#Delta#phi;",nDphi_bins_trig_effcy,-pms.PI,pms.PI,  int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
        h_Dphi_zoomin_vs_pair_pT_mu4_mu4noL1_excl[isign] = new TH2D(Form("h_Dphi_zoomin_vs_pair_pT_mu4_mu4noL1_excl_sign%d",isign+1),";p_{T}^{pair} [GeV];#Delta#phi;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8,  int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
        h_DR_vs_pair_pT_mu4_mu4noL1_excl[isign] = new TH2D(Form("h_DR_vs_pair_pT_mu4_mu4noL1_excl_sign%d",isign+1),";p_{T}^{pair} [GeV];#DeltaR;",nDR_bins_trig_effcy,0,pms.deltaR_thrsh[2],  int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
        h_DR_zoomin_vs_pair_pT_mu4_mu4noL1_excl[isign] = new TH2D(Form("h_DR_zoomin_vs_pair_pT_mu4_mu4noL1_excl_sign%d",isign+1),";p_{T}^{pair} [GeV];#DeltaR;",nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8,  int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
        h_minv_zoomin_vs_pair_pT_mu4_mu4noL1_excl[isign] = new TH2D(Form("h_minv_zoomin_vs_pair_pT_mu4_mu4noL1_excl_sign%d",isign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV];",nminv_bins_trig_effcy,0,3.,  int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
        h_pair_pt_log_vs_pair_pT_mu4_mu4noL1_excl[isign] = new TH2D(Form("h_pair_pt_log_vs_pair_pT_mu4_mu4noL1_excl_sign%d",isign+1),";p_{T}^{pair} [GeV];p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),  int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
        h_pt2nd_vs_pair_pT_mu4_mu4noL1_excl[isign] = new TH2D(Form("h_pt2nd_vs_pair_pT_mu4_mu4noL1_excl_sign%d",isign+1),";p_{T}^{pair} [GeV];p_{T,2nd} [GeV];",int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data(),  int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());

        // pair kinematics requiring 2mu4 trigger
        h_Deta_vs_pT_1st_2mu4[isign] = new TH2D(Form("h_Deta_vs_pT_1st_2mu4_sign%d",isign+1),";p_{T,1st} [GeV];#Delta#eta;",nDeta_bins_trig_effcy,-4.8,4.8,  int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
        h_Deta_zoomin_vs_pT_1st_2mu4[isign] = new TH2D(Form("h_Deta_zoomin_vs_pT_1st_2mu4_sign%d",isign+1),";p_{T,1st} [GeV];#Delta#eta;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8,  int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
        h_Dphi_vs_pT_1st_2mu4[isign] = new TH2D(Form("h_Dphi_vs_pT_1st_2mu4_sign%d",isign+1),";p_{T,1st} [GeV];#Delta#phi;",nDphi_bins_trig_effcy,-pms.PI,pms.PI,  int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
        h_Dphi_zoomin_vs_pT_1st_2mu4[isign] = new TH2D(Form("h_Dphi_zoomin_vs_pT_1st_2mu4_sign%d",isign+1),";p_{T,1st} [GeV];#Delta#phi;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8,  int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
        h_DR_vs_pT_1st_2mu4[isign] = new TH2D(Form("h_DR_vs_pT_1st_2mu4_sign%d",isign+1),";p_{T,1st} [GeV];#DeltaR;",nDR_bins_trig_effcy,0,pms.deltaR_thrsh[2],  int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
        h_DR_zoomin_vs_pT_1st_2mu4[isign] = new TH2D(Form("h_DR_zoomin_vs_pT_1st_2mu4_sign%d",isign+1),";p_{T,1st} [GeV];#DeltaR;",nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8,  int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
        h_minv_zoomin_vs_pT_1st_2mu4[isign] = new TH2D(Form("h_minv_zoomin_vs_pT_1st_2mu4_sign%d",isign+1),";p_{T,1st} [GeV];m_{#mu#mu} [GeV];",nminv_bins_trig_effcy,0,3.,  int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
        h_pair_pt_log_vs_pT_1st_2mu4[isign] = new TH2D(Form("h_pair_pt_log_vs_pT_1st_2mu4_sign%d",isign+1),";p_{T,1st} [GeV];p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),  int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
        h_pt2nd_vs_pT_1st_2mu4[isign] = new TH2D(Form("h_pt2nd_vs_pT_1st_2mu4_sign%d",isign+1),";p_{T,1st} [GeV];p_{T,2nd} [GeV];",int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data(),  int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());

        h_Deta_vs_pair_pT_2mu4[isign] = new TH2D(Form("h_Deta_vs_pair_pT_2mu4_sign%d",isign+1),";p_{T}^{pair} [GeV];#Delta#eta;",nDeta_bins_trig_effcy,-4.8,4.8, int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
        h_Deta_zoomin_vs_pair_pT_2mu4[isign] = new TH2D(Form("h_Deta_zoomin_vs_pair_pT_2mu4_sign%d",isign+1),";p_{T}^{pair} [GeV];#Delta#eta;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8, int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
        h_Dphi_vs_pair_pT_2mu4[isign] = new TH2D(Form("h_Dphi_vs_pair_pT_2mu4_sign%d",isign+1),";p_{T}^{pair} [GeV];#Delta#phi;",nDphi_bins_trig_effcy,-pms.PI,pms.PI, int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
        h_Dphi_zoomin_vs_pair_pT_2mu4[isign] = new TH2D(Form("h_Dphi_zoomin_vs_pair_pT_2mu4_sign%d",isign+1),";p_{T}^{pair} [GeV];#Delta#phi;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8, int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
        h_DR_vs_pair_pT_2mu4[isign] = new TH2D(Form("h_DR_vs_pair_pT_2mu4_sign%d",isign+1),";p_{T}^{pair} [GeV];#DeltaR;",nDR_bins_trig_effcy,0,pms.deltaR_thrsh[2], int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
        h_DR_zoomin_vs_pair_pT_2mu4[isign] = new TH2D(Form("h_DR_zoomin_vs_pair_pT_2mu4_sign%d",isign+1),";p_{T}^{pair} [GeV];#DeltaR;",nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8, int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
        h_minv_zoomin_vs_pair_pT_2mu4[isign] = new TH2D(Form("h_minv_zoomin_vs_pair_pT_2mu4_sign%d",isign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV];",nminv_bins_trig_effcy,0,3., int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
        h_pair_pt_log_vs_pair_pT_2mu4[isign] = new TH2D(Form("h_pair_pt_log_vs_pair_pT_2mu4_sign%d",isign+1),";p_{T}^{pair} [GeV];p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(), int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
        h_pt2nd_vs_pair_pT_2mu4[isign] = new TH2D(Form("h_pt2nd_vs_pair_pT_2mu4_sign%d",isign+1),";p_{T}^{pair} [GeV];p_{T,2nd} [GeV];",int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data(), int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());



        // ----------------- separated (dR > 0.8) -----------------

        // pair kinematics requiring only single mu4 trigger
        h_Deta_mu4_sepr[isign] = new TH1D(Form("h_Deta_mu4_sepr_sign%d",isign+1),";#Delta#eta;",nDeta_bins_trig_effcy,-4.8,4.8);
        h_Deta_zoomin_mu4_sepr[isign] = new TH1D(Form("h_Deta_zoomin_mu4_sepr_sign%d",isign+1),";#Delta#eta;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
        h_Dphi_mu4_sepr[isign] = new TH1D(Form("h_Dphi_mu4_sepr_sign%d",isign+1),";#Delta#phi;",nDphi_bins_trig_effcy,-pms.PI,pms.PI);
        h_Dphi_zoomin_mu4_sepr[isign] = new TH1D(Form("h_Dphi_zoomin_mu4_sepr_sign%d",isign+1),";#Delta#phi;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
        h_DR_mu4_sepr[isign] = new TH1D(Form("h_DR_mu4_sepr_sign%d",isign+1),";#DeltaR;",nDR_bins_trig_effcy,0,pms.deltaR_thrsh[2]);
        h_DR_zoomin_mu4_sepr[isign] = new TH1D(Form("h_DR_zoomin_mu4_sepr_sign%d",isign+1),";#DeltaR;",nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
        h_minv_zoomin_mu4_sepr[isign] = new TH1D(Form("h_minv_zoomin_mu4_sepr_sign%d",isign+1),";m_{#mu#mu} [GeV];",nminv_bins_trig_effcy,0,3.);
        h_pair_pt_log_mu4_sepr[isign] = new TH1D(Form("h_pair_pt_log_mu4_sepr_sign%d",isign+1),";p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
        h_pt2nd_mu4_sepr[isign] = new TH1D(Form("h_pt2nd_mu4_sepr_sign%d",isign+1),";p_{T,2nd} [GeV];",int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());

        h_pt2nd_vs_q_eta_2nd_mu4_sepr[isign] = new TH2D(Form("h_pt2nd_vs_q_eta_2nd_mu4_sepr_sign%d",isign+1),";q*#eta;p_{T} [GeV]",neta_bins_trig_effcy,-2.4,2.4, int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
        h_pair_eta_vs_pair_pT_mu4_sepr[isign] = new TH2D(Form("h_pair_eta_vs_pair_pT_mu4_sepr_sign%d",isign+1),";p_{T}^{pair} [GeV];#eta^{pair}",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),neta_bins_trig_effcy,-2.4,2.4);
        h_Deta_Dphi_mu4_sepr[isign] = new TH2D(Form("h_Deta_Dphi_mu4_sepr_sign%d",isign+1),";#Delta#phi;#Delta#eta", nDphi_bins_trig_effcy,-pms.PI,pms.PI,nDeta_bins_trig_effcy,-4.8,4.8);
        h_eta1_eta2_mu4_sepr[isign] = new TH2D(Form("h_eta1_eta2_mu4_sepr_sign%d",isign+1),";#eta_{sublead};#eta_{lead}",neta_bins_trig_effcy,-2.4,2.4, neta_bins_trig_effcy,-2.4,2.4);
        h_eta_avg_Deta_mu4_sepr[isign] = new TH2D(Form("h_eta_avg_Deta_mu4_sepr_sign%d",isign+1),";#bar{#eta};#Delta#eta", neta_bins_trig_effcy,-2.4,2.4,nDeta_bins_trig_effcy,-4.8,4.8);
        h_eta_avg_Dphi_mu4_sepr[isign] = new TH2D(Form("h_eta_avg_Dphi_mu4_sepr_sign%d",isign+1),";#Delta#phi;#bar{#eta}", nDphi_bins_trig_effcy,-pms.PI,pms.PI,neta_bins_trig_effcy,-2.4,2.4);
        h_minv_pair_pt_log_mu4_sepr[isign] = new TH2D(Form("h_minv_pair_pt_log_mu4_sepr_sign%d",isign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),nminv_bins_log,minv_bins_log[isign]);

        // pair kinematics requiring mu4_mu4noL1 trigger
        h_Deta_mu4_mu4noL1_sepr[isign] = new TH1D(Form("h_Deta_mu4_mu4noL1_sepr_sign%d",isign+1),";#Delta#eta;",nDeta_bins_trig_effcy,-4.8,4.8);
        h_Deta_zoomin_mu4_mu4noL1_sepr[isign] = new TH1D(Form("h_Deta_zoomin_mu4_mu4noL1_sepr_sign%d",isign+1),";#Delta#eta;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
        h_Dphi_mu4_mu4noL1_sepr[isign] = new TH1D(Form("h_Dphi_mu4_mu4noL1_sepr_sign%d",isign+1),";#Delta#phi;",nDphi_bins_trig_effcy,-pms.PI,pms.PI);
        h_Dphi_zoomin_mu4_mu4noL1_sepr[isign] = new TH1D(Form("h_Dphi_zoomin_mu4_mu4noL1_sepr_sign%d",isign+1),";#Delta#phi;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
        h_DR_mu4_mu4noL1_sepr[isign] = new TH1D(Form("h_DR_mu4_mu4noL1_sepr_sign%d",isign+1),";#DeltaR;",nDR_bins_trig_effcy,0,pms.deltaR_thrsh[2]);
        h_DR_zoomin_mu4_mu4noL1_sepr[isign] = new TH1D(Form("h_DR_zoomin_mu4_mu4noL1_sepr_sign%d",isign+1),";#DeltaR;",nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
        h_minv_zoomin_mu4_mu4noL1_sepr[isign] = new TH1D(Form("h_minv_zoomin_mu4_mu4noL1_sepr_sign%d",isign+1),";m_{#mu#mu} [GeV];",nminv_bins_trig_effcy,0,3.);
        h_pair_pt_log_mu4_mu4noL1_sepr[isign] = new TH1D(Form("h_pair_pt_log_mu4_mu4noL1_sepr_sign%d",isign+1),";p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
        h_pt2nd_mu4_mu4noL1_sepr[isign] = new TH1D(Form("h_pt2nd_mu4_mu4noL1_sepr_sign%d",isign+1),";p_{T,2nd} [GeV];",int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());

        h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sepr[isign] = new TH2D(Form("h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sepr_sign%d",isign+1),";q*#eta;p_{T} [GeV]",neta_bins_trig_effcy,-2.4,2.4, int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
        h_pair_eta_vs_pair_pT_mu4_mu4noL1_sepr[isign] = new TH2D(Form("h_pair_eta_vs_pair_pT_mu4_mu4noL1_sepr_sign%d",isign+1),";p_{T}^{pair} [GeV];#eta^{pair}",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),neta_bins_trig_effcy,-2.4,2.4);
        h_Deta_Dphi_mu4_mu4noL1_sepr[isign] = new TH2D(Form("h_Deta_Dphi_mu4_mu4noL1_sepr_sign%d",isign+1),";#Delta#phi;#Delta#eta", nDphi_bins_trig_effcy,-pms.PI,pms.PI,nDeta_bins_trig_effcy,-4.8,4.8);
        h_eta1_eta2_mu4_mu4noL1_sepr[isign] = new TH2D(Form("h_eta1_eta2_mu4_mu4noL1_sepr_sign%d",isign+1),";#eta_{sublead};#eta_{lead}",neta_bins_trig_effcy,-2.4,2.4, neta_bins_trig_effcy,-2.4,2.4);
        h_eta_avg_Deta_mu4_mu4noL1_sepr[isign] = new TH2D(Form("h_eta_avg_Deta_mu4_mu4noL1_sepr_sign%d",isign+1),";#bar{#eta};#Delta#eta", neta_bins_trig_effcy,-2.4,2.4,nDeta_bins_trig_effcy,-4.8,4.8);
        h_eta_avg_Dphi_mu4_mu4noL1_sepr[isign] = new TH2D(Form("h_eta_avg_Dphi_mu4_mu4noL1_sepr_sign%d",isign+1),";#Delta#phi;#bar{#eta}", nDphi_bins_trig_effcy,-pms.PI,pms.PI,neta_bins_trig_effcy,-2.4,2.4);
        h_minv_pair_pt_log_mu4_mu4noL1_sepr[isign] = new TH2D(Form("h_minv_pair_pt_log_mu4_mu4noL1_sepr_sign%d",isign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),nminv_bins_log,minv_bins_log[isign]);

        h_Deta_mu4_mu4noL1_excl_sepr[isign] = new TH1D(Form("h_Deta_mu4_mu4noL1_excl_sepr_sign%d",isign+1),";#Delta#eta;",nDeta_bins_trig_effcy,-4.8,4.8);
        h_Deta_zoomin_mu4_mu4noL1_excl_sepr[isign] = new TH1D(Form("h_Deta_zoomin_mu4_mu4noL1_excl_sepr_sign%d",isign+1),";#Delta#eta;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
        h_Dphi_mu4_mu4noL1_excl_sepr[isign] = new TH1D(Form("h_Dphi_mu4_mu4noL1_excl_sepr_sign%d",isign+1),";#Delta#phi;",nDphi_bins_trig_effcy,-pms.PI,pms.PI);
        h_Dphi_zoomin_mu4_mu4noL1_excl_sepr[isign] = new TH1D(Form("h_Dphi_zoomin_mu4_mu4noL1_excl_sepr_sign%d",isign+1),";#Delta#phi;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
        h_DR_mu4_mu4noL1_excl_sepr[isign] = new TH1D(Form("h_DR_mu4_mu4noL1_excl_sepr_sign%d",isign+1),";#DeltaR;",nDR_bins_trig_effcy,0,pms.deltaR_thrsh[2]);
        h_DR_zoomin_mu4_mu4noL1_excl_sepr[isign] = new TH1D(Form("h_DR_zoomin_mu4_mu4noL1_excl_sepr_sign%d",isign+1),";#DeltaR;",nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
        h_minv_zoomin_mu4_mu4noL1_excl_sepr[isign] = new TH1D(Form("h_minv_zoomin_mu4_mu4noL1_excl_sepr_sign%d",isign+1),";m_{#mu#mu} [GeV];",nminv_bins_trig_effcy,0,3.);
        h_pair_pt_log_mu4_mu4noL1_excl_sepr[isign] = new TH1D(Form("h_pair_pt_log_mu4_mu4noL1_excl_sepr_sign%d",isign+1),";p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
        h_pt2nd_mu4_mu4noL1_excl_sepr[isign] = new TH1D(Form("h_pt2nd_mu4_mu4noL1_excl_sepr_sign%d",isign+1),";p_{T,2nd} [GeV];",int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());

        h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_excl_sepr[isign] = new TH2D(Form("h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_excl_sepr_sign%d",isign+1),";q*#eta;p_{T} [GeV]",neta_bins_trig_effcy,-2.4,2.4, int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
        h_pair_eta_vs_pair_pT_mu4_mu4noL1_excl_sepr[isign] = new TH2D(Form("h_pair_eta_vs_pair_pT_mu4_mu4noL1_excl_sepr_sign%d",isign+1),";p_{T}^{pair} [GeV];#eta^{pair}",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),neta_bins_trig_effcy,-2.4,2.4);
        h_Deta_Dphi_mu4_mu4noL1_excl_sepr[isign] = new TH2D(Form("h_Deta_Dphi_mu4_mu4noL1_excl_sepr_sign%d",isign+1),";#Delta#phi;#Delta#eta", nDphi_bins_trig_effcy,-pms.PI,pms.PI,nDeta_bins_trig_effcy,-4.8,4.8);
        h_eta1_eta2_mu4_mu4noL1_excl_sepr[isign] = new TH2D(Form("h_eta1_eta2_mu4_mu4noL1_excl_sepr_sign%d",isign+1),";#eta_{sublead};#eta_{lead}",neta_bins_trig_effcy,-2.4,2.4, neta_bins_trig_effcy,-2.4,2.4);
        h_eta_avg_Deta_mu4_mu4noL1_excl_sepr[isign] = new TH2D(Form("h_eta_avg_Deta_mu4_mu4noL1_excl_sepr_sign%d",isign+1),";#bar{#eta};#Delta#eta", neta_bins_trig_effcy,-2.4,2.4,nDeta_bins_trig_effcy,-4.8,4.8);
        h_eta_avg_Dphi_mu4_mu4noL1_excl_sepr[isign] = new TH2D(Form("h_eta_avg_Dphi_mu4_mu4noL1_excl_sepr_sign%d",isign+1),";#Delta#phi;#bar{#eta}", nDphi_bins_trig_effcy,-pms.PI,pms.PI,neta_bins_trig_effcy,-2.4,2.4);
        h_minv_pair_pt_log_mu4_mu4noL1_excl_sepr[isign] = new TH2D(Form("h_minv_pair_pt_log_mu4_mu4noL1_excl_sepr_sign%d",isign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),nminv_bins_log,minv_bins_log[isign]);

        // pair kinematics requiring 2mu4 trigger
        h_Deta_2mu4_sepr[isign] = new TH1D(Form("h_Deta_2mu4_sepr_sign%d",isign+1),";#Delta#eta;",nDeta_bins_trig_effcy,-4.8,4.8);
        h_Deta_zoomin_2mu4_sepr[isign] = new TH1D(Form("h_Deta_zoomin_2mu4_sepr_sign%d",isign+1),";#Delta#eta;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
        h_Dphi_2mu4_sepr[isign] = new TH1D(Form("h_Dphi_2mu4_sepr_sign%d",isign+1),";#Delta#phi;",nDphi_bins_trig_effcy,-pms.PI,pms.PI);
        h_Dphi_zoomin_2mu4_sepr[isign] = new TH1D(Form("h_Dphi_zoomin_2mu4_sepr_sign%d",isign+1),";#Delta#phi;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
        h_DR_2mu4_sepr[isign] = new TH1D(Form("h_DR_2mu4_sepr_sign%d",isign+1),";#DeltaR;",nDR_bins_trig_effcy,0,pms.deltaR_thrsh[2]);
        h_DR_zoomin_2mu4_sepr[isign] = new TH1D(Form("h_DR_zoomin_2mu4_sepr_sign%d",isign+1),";#DeltaR;",nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
        h_minv_zoomin_2mu4_sepr[isign] = new TH1D(Form("h_minv_zoomin_2mu4_sepr_sign%d",isign+1),";m_{#mu#mu} [GeV];",nminv_bins_trig_effcy,0,3.);
        h_pair_pt_log_2mu4_sepr[isign] = new TH1D(Form("h_pair_pt_log_2mu4_sepr_sign%d",isign+1),";p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
        h_pt2nd_2mu4_sepr[isign] = new TH1D(Form("h_pt2nd_2mu4_sepr_sign%d",isign+1),";p_{T,2nd} [GeV];",int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());

        h_pt2nd_vs_q_eta_2nd_2mu4_sepr[isign] = new TH2D(Form("h_pt2nd_vs_q_eta_2nd_2mu4_sepr_sign%d",isign+1),";q*#eta;p_{T} [GeV]",neta_bins_trig_effcy,-2.4,2.4, int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
        h_pair_eta_vs_pair_pT_2mu4_sepr[isign] = new TH2D(Form("h_pair_eta_vs_pair_pT_2mu4_sepr_sign%d",isign+1),";p_{T}^{pair} [GeV];#eta^{pair}",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),neta_bins_trig_effcy,-2.4,2.4);
        h_Deta_Dphi_2mu4_sepr[isign] = new TH2D(Form("h_Deta_Dphi_2mu4_sepr_sign%d",isign+1),";#Delta#phi;#Delta#eta", nDphi_bins_trig_effcy,-pms.PI,pms.PI,nDeta_bins_trig_effcy,-4.8,4.8);
        h_eta1_eta2_2mu4_sepr[isign] = new TH2D(Form("h_eta1_eta2_2mu4_sepr_sign%d",isign+1),";#eta_{sublead};#eta_{lead}",neta_bins_trig_effcy,-2.4,2.4, neta_bins_trig_effcy,-2.4,2.4);
        h_eta_avg_Deta_2mu4_sepr[isign] = new TH2D(Form("h_eta_avg_Deta_2mu4_sepr_sign%d",isign+1),";#bar{#eta};#Delta#eta", neta_bins_trig_effcy,-2.4,2.4,nDeta_bins_trig_effcy,-4.8,4.8);
        h_eta_avg_Dphi_2mu4_sepr[isign] = new TH2D(Form("h_eta_avg_Dphi_2mu4_sepr_sign%d",isign+1),";#Delta#phi;#bar{#eta}", nDphi_bins_trig_effcy,-pms.PI,pms.PI,neta_bins_trig_effcy,-2.4,2.4);
        h_minv_pair_pt_log_2mu4_sepr[isign] = new TH2D(Form("h_minv_pair_pt_log_2mu4_sepr_sign%d",isign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),nminv_bins_log,minv_bins_log[isign]);
        
        if (output_non_trig_effcy_hists){
            for (unsigned int idphi = 0; idphi < 2; idphi++){
                for (unsigned int igap = 0; igap < ParamsSet::nGapCuts; igap++){
                    h_pair_dP_overP[idphi][isign][igap] = new TH1D(Form("h_pair_dP_overP_%s_sign%d_gapcut%d",dphi_regions[idphi].c_str(),isign+1,igap+1),";(#Delta p / p)_{pair};1/N_{evt} dN/d(#Delta p / p)_{pair}" ,pms.deltaP_overP_nbins,0,pms.deltaP_overP_max);
                    h_DR[idphi][isign][igap] = new TH1D(Form("h_DR_%s_sign%d_gapcut%d",dphi_regions[idphi].c_str(),isign+1,igap+1),";#Delta R;1/N_{evt} dN/d#Delta R", nDR_bins,0,pms.deltaR_thrsh[2]);
                    h_DR_zoomin[idphi][isign][igap] = new TH1D(Form("h_DR_zoomin_%s_sign%d_gapcut%d",dphi_regions[idphi].c_str(),isign+1,igap+1),";#Delta R;1/N_{evt} dN/d#Delta R", nDR_zoomin_bins,0,pms.deltaR_thrsh_zoomin);
                    h_DR_jacobian_corrected[idphi][isign][igap] = new TH1D(Form("h_DR_jacobian_corrected_%s_sign%d_gapcut%d",dphi_regions[idphi].c_str(),isign+1,igap+1),";#Delta R;#frac{1}{N_{evt}} #frac{1}{#Delta R} #frac{dN}{d#Delta R}", nDR_bins,0,pms.deltaR_thrsh[2]);
                    h_DR_zoomin_jacobian_corrected[idphi][isign][igap] = new TH1D(Form("h_DR_zoomin_jacobian_corrected_%s_sign%d_gapcut%d",dphi_regions[idphi].c_str(),isign+1,igap+1),";#Delta R;#frac{1}{N_{evt}} #frac{1}{#Delta R} #frac{dN}{d#Delta R}", nDR_zoomin_bins,0,pms.deltaR_thrsh_zoomin);
                    h_Dphi[idphi][isign][igap] = new TH1D(Form("h_Dphi_%s_sign%d_gapcut%d",dphi_regions[idphi].c_str(),isign+1,igap+1),";#Delta#phi;1/N_{evt} dN/d#Delta#phi", nDphi_bins,-pms.PI,pms.PI);
                    h_pair_y[idphi][isign][igap] = new TH1D(Form("h_pair_y_%s_sign%d_gapcut%d",dphi_regions[idphi].c_str(),isign+1,igap+1),";y_{pair};1/N_{evt} dN/dy_{pair}" ,npair_y_bins,-3,3);
                    h_pt_asym[idphi][isign][igap] = new TH1D(Form("h_pt_asym_%s_sign%d_gapcut%d",dphi_regions[idphi].c_str(),isign+1,igap+1),";A = (pT1 - pT2)/(pT1 + pT2);d#sigma/dA", npt_asym_bins,0,1.);
                    h_pair_pt_ptlead_ratio[idphi][isign][igap] = new TH1D(Form("h_pair_pt_ptlead_ratio_%s_sign%d_gapcut%d",dphi_regions[idphi].c_str(),isign+1,igap+1),";#frac{p_{T}^{pair}}{p_{T}^{lead}};d#sigma/d#frac{p_{T}^{pair}}{p_{T}^{lead}}", npair_pt_ptlead_ratio_bins,0,2.);

                    h_eta_avg_Dphi[idphi][isign][igap] = new TH2D(Form("h_eta_avg_Dphi_%s_sign%d_gapcut%d",dphi_regions[idphi].c_str(),isign+1,igap+1),";#Delta#phi;#bar{#eta}", nDphi_bins,-pms.PI,pms.PI,neta_bins,-2.4,2.4);
                    h_Deta_Dphi[idphi][isign][igap] = new TH2D(Form("h_Deta_Dphi_%s_sign%d_gapcut%d",dphi_regions[idphi].c_str(),isign+1,igap+1),";#Delta#phi;#Delta#eta", nDphi_bins,-pms.PI,pms.PI,nDeta_bins,-4.8,4.8);
                    h_eta1_eta2[idphi][isign][igap] = new TH2D(Form("h_eta1_eta2_%s_sign%d_gapcut%d",dphi_regions[idphi].c_str(),isign+1,igap+1),";#eta_{sublead};#eta_{lead}",neta_bins,-2.4,2.4, neta_bins,-2.4,2.4);
                    h_eta_avg_Deta[idphi][isign][igap] = new TH2D(Form("h_eta_avg_Deta_%s_sign%d_gapcut%d",dphi_regions[idphi].c_str(),isign+1,igap+1),";#Delta#eta;#bar{#eta}",nDeta_bins,-4.8,4.8,neta_bins,-2.4,2.4);
                    h_pt1_pt2[idphi][isign][igap] = new TH2D(Form("h_pt1_pt2_%s_sign%d_gapcut%d",dphi_regions[idphi].c_str(),isign+1,igap+1),";p_{T}^{sublead} [GeV];p_{T}^{lead} [GeV]",pms.npt_bins,pms.pTBins,pms.npt_bins,pms.pTBins);
                    h_ptlead_pair_pt[idphi][isign][igap] = new TH2D(Form("h_ptlead_pair_pt_%s_sign%d_gapcut%d",dphi_regions[idphi].c_str(),isign+1,igap+1),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",npair_pT_bins_linear,0,30,npT_lead_bins_linear,0,30);
                    h_ptlead_pair_pt_zoomin[idphi][isign][igap] = new TH2D(Form("h_ptlead_pair_pt_zoomin_%s_sign%d_gapcut%d",dphi_regions[idphi].c_str(),isign+1,igap+1),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",npair_pT_bins_linear,0,20,npT_lead_bins_linear,0,15);
                    h_ptlead_pair_pt_log[idphi][isign][igap] = new TH2D(Form("h_ptlead_pair_pt_log_%s_sign%d_gapcut%d",dphi_regions[idphi].c_str(),isign+1,igap+1),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",pms.npairPT_bins,pms.pairPTBins[isign][2],pms.npt_bins,pms.pTBins);
                    h_minv_pair_pt[idphi][isign][igap] = new TH2D(Form("h_minv_pair_pt_%s_sign%d_gapcut%d",dphi_regions[idphi].c_str(),isign+1,igap+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",npair_pT_bins_linear,0,30,nminv_bins_linear,0,30);
                    h_minv_pair_pt_zoomin[idphi][isign][igap] = new TH2D(Form("h_minv_pair_pt_zoomin_%s_sign%d_gapcut%d",dphi_regions[idphi].c_str(),isign+1,igap+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",npair_pT_bins_linear,0,20,nminv_bins_linear,0,4);
                    h_minv_pair_pt_jacobian_corrected[idphi][isign][igap] = new TH2D(Form("h_minv_pair_pt_jacobian_corrected_%s_sign%d_gapcut%d",dphi_regions[idphi].c_str(),isign+1,igap+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",npair_pT_bins_linear,0,30,nminv_bins_linear,0,30);
                    h_minv_pair_pt_zoomin_jacobian_corrected[idphi][isign][igap] = new TH2D(Form("h_minv_pair_pt_zoomin_jacobian_corrected_%s_sign%d_gapcut%d",dphi_regions[idphi].c_str(),isign+1,igap+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",npair_pT_bins_linear,0,20,nminv_bins_linear,0,4);
                    h_minv_pair_pt_log[idphi][isign][igap] = new TH2D(Form("h_minv_pair_pt_log_%s_sign%d_gapcut%d",dphi_regions[idphi].c_str(),isign+1,igap+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",pms.npairPT_bins,pms.pairPTBins[isign][2],nminv_bins_log,minv_bins_log[isign]);
                }
            }
        }
    }
}


void MuonPairPlottingPP::ProcessData(){
	
    cout << "1st loop over muon pairs" << endl;

    for (int isign = 0; isign < ParamsSet::nSigns; isign++){ // 1st loop: fill all nominal histograms
        
        Long64_t nentries = inTree[isign]->GetEntries(); //#muon pairs
        for (Long64_t lentry=0; lentry<nentries; lentry++) {//loop over the events
        // for (Long64_t lentry=0; lentry<100; lentry++) {//loop over the events
        if(lentry%100000==0) cout<<"Processing "<<lentry<<" event out of "<<nentries<<" events"<<std::endl;
      		int num_bytes = inTree[isign]->GetEntry(lentry);//read in an event
            if(num_bytes==0){
            	std::cout<<"Error:: Read in event has size of zero bytes,  quitting"<<std::endl;
            	throw std::exception();
            }
            
            FillHistograms(isign);
        }
    }

    CalculateSingleMuonTrigEffcyRatios();

    cout << "2nd loop over muon pairs" << endl;

    for (int isign = 0; isign < ParamsSet::nSigns; isign++){ // 2nd loop: fill pair efficiency histograms inversely weighted by single-muon trigger efficiencies
        
        Long64_t nentries = inTree[isign]->GetEntries(); //#muon pairs
        for (Long64_t lentry=0; lentry<nentries; lentry++) {//loop over the events
        // for (Long64_t lentry=0; lentry<100; lentry++) {//loop over the events
        if(lentry%100000==0) cout<<"Processing "<<lentry<<" event out of "<<nentries<<" events"<<std::endl;
            int num_bytes = inTree[isign]->GetEntry(lentry);//read in an event
            if(num_bytes==0){
                std::cout<<"Error:: Read in event has size of zero bytes,  quitting"<<std::endl;
                throw std::exception();
            }
            
            FillTrigEffcyHistsInvWeightedbySingleMuonEffcies(isign);
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

    bool pass_single_b_signal_selection = (nsign == 1 && minv[nsign] > 1.08 && minv[nsign] < 2.9 && pair_pt[nsign] > 8);
    bool pass_single_muon_good_acceptance_selection[ParamsSet::nSigns];
    // 1st muon good accept
    // pass_single_muon_good_acceptance_selection[0] = (mu1PassSingle[nsign])? (m1pt[nsign] >= 6 && ((m1eta[nsign] > 1.1 && m1eta[nsign] < 2.3) || (m1eta[nsign] > -2.3 && m1eta[nsign] < -1.2))) : false;
    // pass_single_muon_good_acceptance_selection[1] = (mu2PassSingle[nsign])? (m2pt[nsign] >= 6 && ((m2eta[nsign] > 1.1 && m2eta[nsign] < 2.3) || (m2eta[nsign] > -2.3 && m2eta[nsign] < -1.2))) : false;
    
    // 2nd muon good accept
    pass_single_muon_good_acceptance_selection[0] = (mu2PassSingle[nsign])? (m1pt[nsign] >= 6 && ((m1eta[nsign] > 1.1 && m1eta[nsign] < 2.3) || (m1eta[nsign] > -2.3 && m1eta[nsign] < -1.2))) : false;
    pass_single_muon_good_acceptance_selection[1] = (mu1PassSingle[nsign])? (m2pt[nsign] >= 6 && ((m2eta[nsign] > 1.1 && m2eta[nsign] < 2.3) || (m2eta[nsign] > -2.3 && m2eta[nsign] < -1.2))) : false;

    // ngapcut = 0: all; = 1: only those that pass

    //  bool pass_gapcut = fabs(m1eta[nsign]) > pms.eta_gap_cut && fabs(m2eta[nsign]) > pms.eta_gap_cut;
    bool pass_gapcut = PassSingleMuonGapCut(m1eta[nsign], m1pt[nsign], m1charge[nsign]) && PassSingleMuonGapCut(m2eta[nsign], m2pt[nsign], m2charge[nsign]);
    bool away_side = (abs(dphi[nsign]) >= pms.PI / 2.);
    if (pass_single_b_signal_selection){
        h_crossx_dR_cut_vs_pair_pt_pair_eta->Fill(pair_eta[nsign],pair_pt[nsign],weight[nsign]);
    }

    if (output_non_trig_effcy_hists){
        if (pass_gapcut){

            h_pair_dP_overP[away_side][nsign][1]->Fill(pair_dPoverP[nsign],weight[nsign]);
            h_pair_y[away_side][nsign][1]->Fill(pair_y[nsign],weight[nsign]);
            h_DR[away_side][nsign][1]->Fill(dr[nsign],weight[nsign]);
            h_DR_zoomin[away_side][nsign][1]->Fill(dr[nsign],weight[nsign]);
            h_DR_jacobian_corrected[away_side][nsign][1]->Fill(dr[nsign],weight[nsign] * 1. / dr[nsign]);
            h_DR_zoomin_jacobian_corrected[away_side][nsign][1]->Fill(dr[nsign],weight[nsign] * 1. / dr[nsign]);
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
            h_minv_pair_pt_jacobian_corrected[away_side][nsign][1]->Fill(pair_pt[nsign],minv[nsign],weight[nsign] * 1. / dr[nsign]);
            h_minv_pair_pt_zoomin_jacobian_corrected[away_side][nsign][1]->Fill(pair_pt[nsign],minv[nsign],weight[nsign] * 1. / dr[nsign]);
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
        h_DR_zoomin[away_side][nsign][0]->Fill(dr[nsign],weight[nsign]);
        h_DR_jacobian_corrected[away_side][nsign][0]->Fill(dr[nsign],weight[nsign] * 1. / dr[nsign]);
        h_DR_zoomin_jacobian_corrected[away_side][nsign][0]->Fill(dr[nsign],weight[nsign] * 1. / dr[nsign]);
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
        h_minv_pair_pt_jacobian_corrected[away_side][nsign][0]->Fill(pair_pt[nsign],minv[nsign],weight[nsign] * 1. / dr[nsign]);
        h_minv_pair_pt_zoomin_jacobian_corrected[away_side][nsign][0]->Fill(pair_pt[nsign],minv[nsign],weight[nsign] * 1. / dr[nsign]);
        h_minv_pair_pt_log[away_side][nsign][0]->Fill(pair_pt[nsign],minv[nsign],weight[nsign]);
    }

    h_Deta_mu4[nsign]->Fill(deta[nsign],weight[nsign]);
    h_Deta_zoomin_mu4[nsign]->Fill(deta[nsign],weight[nsign]);
    h_Dphi_mu4[nsign]->Fill(dphi[nsign],weight[nsign]);
    h_Dphi_zoomin_mu4[nsign]->Fill(dphi[nsign],weight[nsign]);
    h_DR_mu4[nsign]->Fill(dr[nsign],weight[nsign]);
    h_DR_zoomin_mu4[nsign]->Fill(dr[nsign],weight[nsign]);
    h_DR_0_2_mu4[nsign]->Fill(dr[nsign],weight[nsign]);

    h_Deta_vs_pair_pT_mu4[nsign]->Fill(pair_pt[nsign], deta[nsign],weight[nsign]);
    h_Deta_zoomin_vs_pair_pT_mu4[nsign]->Fill(pair_pt[nsign], deta[nsign],weight[nsign]);
    h_Dphi_vs_pair_pT_mu4[nsign]->Fill(pair_pt[nsign], dphi[nsign],weight[nsign]);
    h_Dphi_zoomin_vs_pair_pT_mu4[nsign]->Fill(pair_pt[nsign], dphi[nsign],weight[nsign]);
    h_DR_vs_pair_pT_mu4[nsign]->Fill(pair_pt[nsign], dr[nsign],weight[nsign]);
    h_DR_zoomin_vs_pair_pT_mu4[nsign]->Fill(pair_pt[nsign], dr[nsign],weight[nsign]);

    if (mu1PassSingle[nsign]){
        h_pt2nd_mu4[nsign]->Fill(m2pt[nsign],weight[nsign]);
        h_pt2nd_vs_q_eta_2nd_mu4[nsign]->Fill(m2charge[nsign] * m2eta[nsign], m2pt[nsign],weight[nsign]);
        h_pt2nd_vs_phi2nd_mu4[nsign]->Fill(m2phi[nsign], m2pt[nsign], weight[nsign]);
        h_phi2nd_vs_q_eta_2nd_mu4[nsign]->Fill(m2charge[nsign]*m2eta[nsign], m2phi[nsign], weight[nsign]);
        h_DR_zoomin_vs_pt2nd_mu4[nsign]->Fill(m2pt[nsign],dr[nsign],weight[nsign]);
        h_DR_0_2_vs_pt2nd_mu4[nsign]->Fill(m2pt[nsign],dr[nsign],weight[nsign]);

        h_Deta_vs_pT_1st_mu4[nsign]->Fill(m2pt[nsign], deta[nsign],weight[nsign]);
        h_Deta_zoomin_vs_pT_1st_mu4[nsign]->Fill(m2pt[nsign], deta[nsign],weight[nsign]);
        h_Dphi_vs_pT_1st_mu4[nsign]->Fill(m2pt[nsign], dphi[nsign],weight[nsign]);
        h_Dphi_zoomin_vs_pT_1st_mu4[nsign]->Fill(m2pt[nsign], dphi[nsign],weight[nsign]);
        h_DR_vs_pT_1st_mu4[nsign]->Fill(m2pt[nsign], dr[nsign],weight[nsign]);
        h_DR_zoomin_vs_pT_1st_mu4[nsign]->Fill(m2pt[nsign], dr[nsign],weight[nsign]);
    }
    if (mu2PassSingle[nsign]){
        h_pt2nd_mu4[nsign]->Fill(m1pt[nsign],weight[nsign]);
        h_pt2nd_vs_q_eta_2nd_mu4[nsign]->Fill(m1charge[nsign] * m1eta[nsign], m1pt[nsign],weight[nsign]);
        h_pt2nd_vs_phi2nd_mu4[nsign]->Fill(m1phi[nsign], m1pt[nsign], weight[nsign]);
        h_phi2nd_vs_q_eta_2nd_mu4[nsign]->Fill(m1charge[nsign]*m1eta[nsign], m1phi[nsign], weight[nsign]);
        h_DR_zoomin_vs_pt2nd_mu4[nsign]->Fill(m1pt[nsign],dr[nsign],weight[nsign]);
        h_DR_0_2_vs_pt2nd_mu4[nsign]->Fill(m1pt[nsign],dr[nsign],weight[nsign]);

        h_Deta_vs_pT_1st_mu4[nsign]->Fill(m1pt[nsign], deta[nsign],weight[nsign]);
        h_Deta_zoomin_vs_pT_1st_mu4[nsign]->Fill(m1pt[nsign], deta[nsign],weight[nsign]);
        h_Dphi_vs_pT_1st_mu4[nsign]->Fill(m1pt[nsign], dphi[nsign],weight[nsign]);
        h_Dphi_zoomin_vs_pT_1st_mu4[nsign]->Fill(m1pt[nsign], dphi[nsign],weight[nsign]);
        h_DR_vs_pT_1st_mu4[nsign]->Fill(m1pt[nsign], dr[nsign],weight[nsign]);
        h_DR_zoomin_vs_pT_1st_mu4[nsign]->Fill(m1pt[nsign], dr[nsign],weight[nsign]);
    }

    h_pair_eta_vs_pair_pT_mu4[nsign]->Fill(pair_pt[nsign],pair_eta[nsign],weight[nsign]);
    h_Deta_Dphi_mu4[nsign]->Fill(dphi[nsign],deta[nsign],weight[nsign]);
    h_eta1_eta2_mu4[nsign]->Fill(m2eta[nsign],m1eta[nsign],weight[nsign]);
    h_eta_avg_Deta_mu4[nsign]->Fill(deta[nsign],etaavg[nsign],weight[nsign]);
    h_eta_avg_Dphi_mu4[nsign]->Fill(dphi[nsign],etaavg[nsign],weight[nsign]);
    h_minv_zoomin_mu4[nsign]->Fill(minv[nsign],weight[nsign]);
    h_pair_pt_log_mu4[nsign]->Fill(pair_pt[nsign],weight[nsign]);
    h_minv_pair_pt_log_mu4[nsign]->Fill(pair_pt[nsign],minv[nsign],weight[nsign]);

    for (int i = 0; i < ParamsSet::nSigns; i++){
        if (pass_single_muon_good_acceptance_selection[i]){ // fill for both
            h_Deta_mu4_good_accept[nsign]->Fill(deta[nsign],weight[nsign]);
            h_Deta_zoomin_mu4_good_accept[nsign]->Fill(deta[nsign],weight[nsign]);
            h_Dphi_mu4_good_accept[nsign]->Fill(dphi[nsign],weight[nsign]);
            h_Dphi_zoomin_mu4_good_accept[nsign]->Fill(dphi[nsign],weight[nsign]);
            h_DR_mu4_good_accept[nsign]->Fill(dr[nsign],weight[nsign]);
            h_DR_zoomin_mu4_good_accept[nsign]->Fill(dr[nsign],weight[nsign]);
            h_minv_zoomin_mu4_good_accept[nsign]->Fill(minv[nsign],weight[nsign]);
            h_pair_pt_log_mu4_good_accept[nsign]->Fill(pair_pt[nsign],weight[nsign]);    
        }
    }

    if (passSeparated[nsign]){
        h_Deta_mu4_sepr[nsign]->Fill(deta[nsign],weight[nsign]);
        h_Deta_zoomin_mu4_sepr[nsign]->Fill(deta[nsign],weight[nsign]);
        h_Deta_Dphi_mu4_sepr[nsign]->Fill(dphi[nsign],deta[nsign],weight[nsign]);
        h_eta_avg_Deta_mu4_sepr[nsign]->Fill(deta[nsign],etaavg[nsign],weight[nsign]);
    }

    if (passSeparated[nsign]){
        h_Dphi_mu4_sepr[nsign]->Fill(dphi[nsign],weight[nsign]);
        h_Dphi_zoomin_mu4_sepr[nsign]->Fill(dphi[nsign],weight[nsign]);
        h_DR_mu4_sepr[nsign]->Fill(dr[nsign],weight[nsign]);
        h_DR_zoomin_mu4_sepr[nsign]->Fill(dr[nsign],weight[nsign]);
        if (mu1PassSingle[nsign]){
            h_pt2nd_mu4_sepr[nsign]->Fill(m2pt[nsign],weight[nsign]);
            h_pt2nd_vs_q_eta_2nd_mu4_sepr[nsign]->Fill(m2charge[nsign] * m2eta[nsign], m2pt[nsign],weight[nsign]);
        }
        if (mu2PassSingle[nsign]){
            h_pt2nd_mu4_sepr[nsign]->Fill(m1pt[nsign],weight[nsign]);
            h_pt2nd_vs_q_eta_2nd_mu4_sepr[nsign]->Fill(m1charge[nsign] * m1eta[nsign], m1pt[nsign],weight[nsign]);
        }

        h_pair_eta_vs_pair_pT_mu4_sepr[nsign]->Fill(pair_pt[nsign],pair_eta[nsign],weight[nsign]);
        h_eta1_eta2_mu4_sepr[nsign]->Fill(m2eta[nsign],m1eta[nsign],weight[nsign]);
        h_eta_avg_Dphi_mu4_sepr[nsign]->Fill(dphi[nsign],etaavg[nsign],weight[nsign]);
        h_minv_zoomin_mu4_sepr[nsign]->Fill(minv[nsign],weight[nsign]);
        h_pair_pt_log_mu4_sepr[nsign]->Fill(pair_pt[nsign],weight[nsign]);
        h_minv_pair_pt_log_mu4_sepr[nsign]->Fill(pair_pt[nsign],minv[nsign],weight[nsign]);
    }

    if (pass_single_b_signal_selection){
        h_DR_zoomin_vs_pair_pT_mu4_w_single_b_sig_sel->Fill(pair_pt[nsign], dr[nsign],weight[nsign]);

        h_Deta_mu4_w_single_b_sig_sel->Fill(deta[nsign],weight[nsign]);
        h_Deta_zoomin_mu4_w_single_b_sig_sel->Fill(deta[nsign],weight[nsign]);
        h_Dphi_mu4_w_single_b_sig_sel->Fill(dphi[nsign],weight[nsign]);
        h_Dphi_zoomin_mu4_w_single_b_sig_sel->Fill(dphi[nsign],weight[nsign]);
        h_DR_mu4_w_single_b_sig_sel->Fill(dr[nsign],weight[nsign]);
        h_DR_zoomin_mu4_w_single_b_sig_sel->Fill(dr[nsign],weight[nsign]);
        h_DR_0_2_mu4_w_single_b_sig_sel->Fill(dr[nsign],weight[nsign]);

        if (mu1PassSingle[nsign]){
            h_pt2nd_mu4_w_single_b_sig_sel->Fill(m2pt[nsign],weight[nsign]);
            h_pt2nd_vs_q_eta_2nd_mu4_w_single_b_sig_sel->Fill(m2charge[nsign] * m2eta[nsign], m2pt[nsign],weight[nsign]);
            h_DR_zoomin_vs_pt2nd_mu4_w_single_b_sig_sel->Fill(m2pt[nsign],dr[nsign],weight[nsign]);
            h_DR_0_2_vs_pt2nd_mu4_w_single_b_sig_sel->Fill(m2pt[nsign],dr[nsign],weight[nsign]);
        }
        if (mu2PassSingle[nsign]){
            h_pt2nd_mu4_w_single_b_sig_sel->Fill(m1pt[nsign],weight[nsign]);
            h_pt2nd_vs_q_eta_2nd_mu4_w_single_b_sig_sel->Fill(m1charge[nsign] * m1eta[nsign], m1pt[nsign],weight[nsign]);
            h_DR_zoomin_vs_pt2nd_mu4_w_single_b_sig_sel->Fill(m1pt[nsign],dr[nsign],weight[nsign]);
            h_DR_0_2_vs_pt2nd_mu4_w_single_b_sig_sel->Fill(m1pt[nsign],dr[nsign],weight[nsign]);
        }

        h_pair_eta_vs_pair_pT_mu4_w_single_b_sig_sel->Fill(pair_pt[nsign],pair_eta[nsign],weight[nsign]);
        h_Deta_Dphi_mu4_w_single_b_sig_sel->Fill(dphi[nsign],deta[nsign],weight[nsign]);
        h_eta1_eta2_mu4_w_single_b_sig_sel->Fill(m2eta[nsign],m1eta[nsign],weight[nsign]);
        h_eta_avg_Deta_mu4_w_single_b_sig_sel->Fill(deta[nsign],etaavg[nsign],weight[nsign]);
        h_eta_avg_Dphi_mu4_w_single_b_sig_sel->Fill(dphi[nsign],etaavg[nsign],weight[nsign]);
        h_minv_zoomin_mu4_w_single_b_sig_sel->Fill(minv[nsign],weight[nsign]);
        h_pair_pt_log_mu4_w_single_b_sig_sel->Fill(pair_pt[nsign],weight[nsign]);
        h_minv_pair_pt_log_mu4_w_single_b_sig_sel->Fill(pair_pt[nsign],minv[nsign],weight[nsign]);
    }

    if (passmu4mu4noL1[nsign]){
        h_Deta_mu4_mu4noL1[nsign]->Fill(deta[nsign],weight[nsign]);
        h_Deta_zoomin_mu4_mu4noL1[nsign]->Fill(deta[nsign],weight[nsign]);
        h_Dphi_mu4_mu4noL1[nsign]->Fill(dphi[nsign],weight[nsign]);
        h_Dphi_zoomin_mu4_mu4noL1[nsign]->Fill(dphi[nsign],weight[nsign]);
        h_DR_mu4_mu4noL1[nsign]->Fill(dr[nsign],weight[nsign]);
        h_DR_zoomin_mu4_mu4noL1[nsign]->Fill(dr[nsign],weight[nsign]);
        h_DR_0_2_mu4_mu4noL1[nsign]->Fill(dr[nsign],weight[nsign]);

        h_Deta_vs_pair_pT_mu4_mu4noL1[nsign]->Fill(pair_pt[nsign], deta[nsign],weight[nsign]);
        h_Deta_zoomin_vs_pair_pT_mu4_mu4noL1[nsign]->Fill(pair_pt[nsign], deta[nsign],weight[nsign]);
        h_Dphi_vs_pair_pT_mu4_mu4noL1[nsign]->Fill(pair_pt[nsign], dphi[nsign],weight[nsign]);
        h_Dphi_zoomin_vs_pair_pT_mu4_mu4noL1[nsign]->Fill(pair_pt[nsign], dphi[nsign],weight[nsign]);
        h_DR_vs_pair_pT_mu4_mu4noL1[nsign]->Fill(pair_pt[nsign], dr[nsign],weight[nsign]);
        h_DR_zoomin_vs_pair_pT_mu4_mu4noL1[nsign]->Fill(pair_pt[nsign], dr[nsign],weight[nsign]);

        if (mu1PassSingle[nsign]){
            h_pt2nd_mu4_mu4noL1[nsign]->Fill(m2pt[nsign],weight[nsign]);
            h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1[nsign]->Fill(m2charge[nsign] * m2eta[nsign], m2pt[nsign],weight[nsign]);
            h_pt2nd_vs_phi2nd_mu4_mu4noL1[nsign]->Fill(m2phi[nsign], m2pt[nsign], weight[nsign]);
            h_phi2nd_vs_q_eta_2nd_mu4_mu4noL1[nsign]->Fill(m2charge[nsign]*m2eta[nsign], m2phi[nsign], weight[nsign]);
            h_DR_zoomin_vs_pt2nd_mu4_mu4noL1[nsign]->Fill(m2pt[nsign],dr[nsign],weight[nsign]);
            h_DR_0_2_vs_pt2nd_mu4_mu4noL1[nsign]->Fill(m2pt[nsign],dr[nsign],weight[nsign]);

            h_Deta_vs_pT_1st_mu4_mu4noL1[nsign]->Fill(m2pt[nsign], deta[nsign],weight[nsign]);
            h_Deta_zoomin_vs_pT_1st_mu4_mu4noL1[nsign]->Fill(m2pt[nsign], deta[nsign],weight[nsign]);
            h_Dphi_vs_pT_1st_mu4_mu4noL1[nsign]->Fill(m2pt[nsign], dphi[nsign],weight[nsign]);
            h_Dphi_zoomin_vs_pT_1st_mu4_mu4noL1[nsign]->Fill(m2pt[nsign], dphi[nsign],weight[nsign]);
            h_DR_vs_pT_1st_mu4_mu4noL1[nsign]->Fill(m2pt[nsign], dr[nsign],weight[nsign]);
            h_DR_zoomin_vs_pT_1st_mu4_mu4noL1[nsign]->Fill(m2pt[nsign], dr[nsign],weight[nsign]);

        }
        if (mu2PassSingle[nsign]){
            h_pt2nd_mu4_mu4noL1[nsign]->Fill(m1pt[nsign],weight[nsign]);
            h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1[nsign]->Fill(m1charge[nsign] * m1eta[nsign], m1pt[nsign],weight[nsign]);
            h_pt2nd_vs_phi2nd_mu4_mu4noL1[nsign]->Fill(m1phi[nsign], m1pt[nsign], weight[nsign]);
            h_phi2nd_vs_q_eta_2nd_mu4_mu4noL1[nsign]->Fill(m1charge[nsign]*m1eta[nsign], m1phi[nsign], weight[nsign]);
            h_DR_zoomin_vs_pt2nd_mu4_mu4noL1[nsign]->Fill(m1pt[nsign],dr[nsign],weight[nsign]);
            h_DR_0_2_vs_pt2nd_mu4_mu4noL1[nsign]->Fill(m1pt[nsign],dr[nsign],weight[nsign]);

            h_Deta_vs_pT_1st_mu4_mu4noL1[nsign]->Fill(m1pt[nsign], deta[nsign],weight[nsign]);
            h_Deta_zoomin_vs_pT_1st_mu4_mu4noL1[nsign]->Fill(m1pt[nsign], deta[nsign],weight[nsign]);
            h_Dphi_vs_pT_1st_mu4_mu4noL1[nsign]->Fill(m1pt[nsign], dphi[nsign],weight[nsign]);
            h_Dphi_zoomin_vs_pT_1st_mu4_mu4noL1[nsign]->Fill(m1pt[nsign], dphi[nsign],weight[nsign]);
            h_DR_vs_pT_1st_mu4_mu4noL1[nsign]->Fill(m1pt[nsign], dr[nsign],weight[nsign]);
            h_DR_zoomin_vs_pT_1st_mu4_mu4noL1[nsign]->Fill(m1pt[nsign], dr[nsign],weight[nsign]);
        }

        h_pair_eta_vs_pair_pT_mu4_mu4noL1[nsign]->Fill(pair_pt[nsign],pair_eta[nsign],weight[nsign]);
        h_Deta_Dphi_mu4_mu4noL1[nsign]->Fill(dphi[nsign],deta[nsign],weight[nsign]);
        h_eta1_eta2_mu4_mu4noL1[nsign]->Fill(m2eta[nsign],m1eta[nsign],weight[nsign]);
        h_eta_avg_Deta_mu4_mu4noL1[nsign]->Fill(deta[nsign],etaavg[nsign],weight[nsign]);
        h_eta_avg_Dphi_mu4_mu4noL1[nsign]->Fill(dphi[nsign],etaavg[nsign],weight[nsign]);
        h_minv_zoomin_mu4_mu4noL1[nsign]->Fill(minv[nsign],weight[nsign]);
        h_pair_pt_log_mu4_mu4noL1[nsign]->Fill(pair_pt[nsign],weight[nsign]);
        h_minv_pair_pt_log_mu4_mu4noL1[nsign]->Fill(pair_pt[nsign],minv[nsign],weight[nsign]);

        for (int i = 0; i < ParamsSet::nSigns; i++){
            if (pass_single_muon_good_acceptance_selection[i]){ // fill for both
                h_Deta_mu4_mu4noL1_good_accept[nsign]->Fill(deta[nsign],weight[nsign]);
                h_Deta_zoomin_mu4_mu4noL1_good_accept[nsign]->Fill(deta[nsign],weight[nsign]);
                h_Dphi_mu4_mu4noL1_good_accept[nsign]->Fill(dphi[nsign],weight[nsign]);
                h_Dphi_zoomin_mu4_mu4noL1_good_accept[nsign]->Fill(dphi[nsign],weight[nsign]);
                h_DR_mu4_mu4noL1_good_accept[nsign]->Fill(dr[nsign],weight[nsign]);
                h_DR_zoomin_mu4_mu4noL1_good_accept[nsign]->Fill(dr[nsign],weight[nsign]);
                h_minv_zoomin_mu4_mu4noL1_good_accept[nsign]->Fill(minv[nsign],weight[nsign]);
                h_pair_pt_log_mu4_mu4noL1_good_accept[nsign]->Fill(pair_pt[nsign],weight[nsign]);    
            }
        }

        if ((!mu1PassSingle[nsign] || !mu2PassSingle[nsign])){ // exclusive: only one muon passes mu4
            h_Deta_mu4_mu4noL1_excl[nsign]->Fill(deta[nsign],weight[nsign]);
            h_Deta_zoomin_mu4_mu4noL1_excl[nsign]->Fill(deta[nsign],weight[nsign]);
            h_Dphi_mu4_mu4noL1_excl[nsign]->Fill(dphi[nsign],weight[nsign]);
            h_Dphi_zoomin_mu4_mu4noL1_excl[nsign]->Fill(dphi[nsign],weight[nsign]);
            h_DR_mu4_mu4noL1_excl[nsign]->Fill(dr[nsign],weight[nsign]);
            h_DR_zoomin_mu4_mu4noL1_excl[nsign]->Fill(dr[nsign],weight[nsign]);
            if (mu1PassSingle[nsign]){
                h_pt2nd_mu4_mu4noL1_excl[nsign]->Fill(m2pt[nsign],weight[nsign]);
                h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_excl[nsign]->Fill(m2charge[nsign] * m2eta[nsign], m2pt[nsign],weight[nsign]);
            }
            if (mu2PassSingle[nsign]){
                h_pt2nd_mu4_mu4noL1_excl[nsign]->Fill(m1pt[nsign],weight[nsign]);
                h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_excl[nsign]->Fill(m1charge[nsign] * m1eta[nsign], m1pt[nsign],weight[nsign]);
            }

            h_pair_eta_vs_pair_pT_mu4_mu4noL1_excl[nsign]->Fill(pair_pt[nsign],pair_eta[nsign],weight[nsign]);
            h_Deta_Dphi_mu4_mu4noL1_excl[nsign]->Fill(dphi[nsign],deta[nsign],weight[nsign]);
            h_eta1_eta2_mu4_mu4noL1_excl[nsign]->Fill(m2eta[nsign],m1eta[nsign],weight[nsign]);
            h_eta_avg_Deta_mu4_mu4noL1_excl[nsign]->Fill(deta[nsign],etaavg[nsign],weight[nsign]);
            h_eta_avg_Dphi_mu4_mu4noL1_excl[nsign]->Fill(dphi[nsign],etaavg[nsign],weight[nsign]);
            h_minv_zoomin_mu4_mu4noL1_excl[nsign]->Fill(minv[nsign],weight[nsign]);
            h_pair_pt_log_mu4_mu4noL1_excl[nsign]->Fill(pair_pt[nsign],weight[nsign]);
            h_minv_pair_pt_log_mu4_mu4noL1_excl[nsign]->Fill(pair_pt[nsign],minv[nsign],weight[nsign]);
        }

        if (passSeparated[nsign]){
            h_Deta_mu4_mu4noL1_sepr[nsign]->Fill(deta[nsign],weight[nsign]);
            h_Deta_zoomin_mu4_mu4noL1_sepr[nsign]->Fill(deta[nsign],weight[nsign]);
            h_Deta_Dphi_mu4_mu4noL1_sepr[nsign]->Fill(dphi[nsign],deta[nsign],weight[nsign]);
            h_eta_avg_Deta_mu4_mu4noL1_sepr[nsign]->Fill(deta[nsign],etaavg[nsign],weight[nsign]);
            if ((!mu1PassSingle[nsign] || !mu2PassSingle[nsign])){ // exclusive: only one muon passes mu4
                h_Deta_mu4_mu4noL1_excl_sepr[nsign]->Fill(deta[nsign],weight[nsign]);
                h_Deta_zoomin_mu4_mu4noL1_excl_sepr[nsign]->Fill(deta[nsign],weight[nsign]);
                h_Deta_Dphi_mu4_mu4noL1_excl_sepr[nsign]->Fill(dphi[nsign],deta[nsign],weight[nsign]);
                h_eta_avg_Deta_mu4_mu4noL1_excl_sepr[nsign]->Fill(deta[nsign],etaavg[nsign],weight[nsign]);
            }
        }

        if (passSeparated[nsign]){
            h_Dphi_mu4_mu4noL1_sepr[nsign]->Fill(dphi[nsign],weight[nsign]);
            h_Dphi_zoomin_mu4_mu4noL1_sepr[nsign]->Fill(dphi[nsign],weight[nsign]);
            h_DR_mu4_mu4noL1_sepr[nsign]->Fill(dr[nsign],weight[nsign]);
            h_DR_zoomin_mu4_mu4noL1_sepr[nsign]->Fill(dr[nsign],weight[nsign]);
            if (mu1PassSingle[nsign]){
                h_pt2nd_mu4_mu4noL1_sepr[nsign]->Fill(m2pt[nsign],weight[nsign]);
                h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sepr[nsign]->Fill(m2charge[nsign] * m2eta[nsign], m2pt[nsign],weight[nsign]);
            }
            if (mu2PassSingle[nsign]){
                h_pt2nd_mu4_mu4noL1_sepr[nsign]->Fill(m1pt[nsign],weight[nsign]);
                h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sepr[nsign]->Fill(m1charge[nsign] * m1eta[nsign], m1pt[nsign],weight[nsign]);
            }

            h_pair_eta_vs_pair_pT_mu4_mu4noL1_sepr[nsign]->Fill(pair_pt[nsign],pair_eta[nsign],weight[nsign]);
            h_eta1_eta2_mu4_mu4noL1_sepr[nsign]->Fill(m2eta[nsign],m1eta[nsign],weight[nsign]);
            h_eta_avg_Dphi_mu4_mu4noL1_sepr[nsign]->Fill(dphi[nsign],etaavg[nsign],weight[nsign]);
            h_minv_zoomin_mu4_mu4noL1_sepr[nsign]->Fill(minv[nsign],weight[nsign]);
            h_pair_pt_log_mu4_mu4noL1_sepr[nsign]->Fill(pair_pt[nsign],weight[nsign]);
            h_minv_pair_pt_log_mu4_mu4noL1_sepr[nsign]->Fill(pair_pt[nsign],minv[nsign],weight[nsign]);

            if ((!mu1PassSingle[nsign] || !mu2PassSingle[nsign])){ // exclusive: only one muon passes mu4
                h_Dphi_mu4_mu4noL1_excl_sepr[nsign]->Fill(dphi[nsign],weight[nsign]);
                h_Dphi_zoomin_mu4_mu4noL1_excl_sepr[nsign]->Fill(dphi[nsign],weight[nsign]);
                h_DR_mu4_mu4noL1_excl_sepr[nsign]->Fill(dr[nsign],weight[nsign]);
                h_DR_zoomin_mu4_mu4noL1_excl_sepr[nsign]->Fill(dr[nsign],weight[nsign]);
                if (mu1PassSingle[nsign]){
                    h_pt2nd_mu4_mu4noL1_excl_sepr[nsign]->Fill(m2pt[nsign],weight[nsign]);
                    h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_excl_sepr[nsign]->Fill(m2charge[nsign] * m2eta[nsign], m2pt[nsign],weight[nsign]);
                }
                if (mu2PassSingle[nsign]){
                    h_pt2nd_mu4_mu4noL1_excl_sepr[nsign]->Fill(m1pt[nsign],weight[nsign]);
                    h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_excl_sepr[nsign]->Fill(m1charge[nsign] * m1eta[nsign], m1pt[nsign],weight[nsign]);
                }

                h_pair_eta_vs_pair_pT_mu4_mu4noL1_excl_sepr[nsign]->Fill(pair_pt[nsign],pair_eta[nsign],weight[nsign]);
                h_eta1_eta2_mu4_mu4noL1_excl_sepr[nsign]->Fill(m2eta[nsign],m1eta[nsign],weight[nsign]);
                h_eta_avg_Dphi_mu4_mu4noL1_excl_sepr[nsign]->Fill(dphi[nsign],etaavg[nsign],weight[nsign]);
                h_minv_zoomin_mu4_mu4noL1_excl_sepr[nsign]->Fill(minv[nsign],weight[nsign]);
                h_pair_pt_log_mu4_mu4noL1_excl_sepr[nsign]->Fill(pair_pt[nsign],weight[nsign]);
                h_minv_pair_pt_log_mu4_mu4noL1_excl_sepr[nsign]->Fill(pair_pt[nsign],minv[nsign],weight[nsign]);
            }
        }

        if (pass_single_b_signal_selection){
            h_Deta_mu4_mu4noL1_w_single_b_sig_sel->Fill(deta[nsign],weight[nsign]);
            h_Deta_zoomin_mu4_mu4noL1_w_single_b_sig_sel->Fill(deta[nsign],weight[nsign]);
            h_Dphi_mu4_mu4noL1_w_single_b_sig_sel->Fill(dphi[nsign],weight[nsign]);
            h_Dphi_zoomin_mu4_mu4noL1_w_single_b_sig_sel->Fill(dphi[nsign],weight[nsign]);
            h_DR_mu4_mu4noL1_w_single_b_sig_sel->Fill(dr[nsign],weight[nsign]);
            h_DR_zoomin_mu4_mu4noL1_w_single_b_sig_sel->Fill(dr[nsign],weight[nsign]);
            h_DR_0_2_mu4_mu4noL1_w_single_b_sig_sel->Fill(dr[nsign],weight[nsign]);

            if (mu1PassSingle[nsign]){
                h_pt2nd_mu4_mu4noL1_w_single_b_sig_sel->Fill(m2pt[nsign],weight[nsign]);
                h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_w_single_b_sig_sel->Fill(m2charge[nsign] * m2eta[nsign], m2pt[nsign],weight[nsign]);
                h_DR_zoomin_vs_pt2nd_mu4_mu4noL1_w_single_b_sig_sel->Fill(m2pt[nsign],dr[nsign],weight[nsign]);
                h_DR_0_2_vs_pt2nd_mu4_mu4noL1_w_single_b_sig_sel->Fill(m2pt[nsign],dr[nsign],weight[nsign]);

            }
            if (mu2PassSingle[nsign]){
                h_pt2nd_mu4_mu4noL1_w_single_b_sig_sel->Fill(m1pt[nsign],weight[nsign]);
                h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_w_single_b_sig_sel->Fill(m1charge[nsign] * m1eta[nsign], m1pt[nsign],weight[nsign]);
                h_DR_zoomin_vs_pt2nd_mu4_mu4noL1_w_single_b_sig_sel->Fill(m1pt[nsign],dr[nsign],weight[nsign]);
                h_DR_0_2_vs_pt2nd_mu4_mu4noL1_w_single_b_sig_sel->Fill(m1pt[nsign],dr[nsign],weight[nsign]);
            }

            h_pair_eta_vs_pair_pT_mu4_mu4noL1_w_single_b_sig_sel->Fill(pair_pt[nsign],pair_eta[nsign],weight[nsign]);
            h_Deta_Dphi_mu4_mu4noL1_w_single_b_sig_sel->Fill(dphi[nsign],deta[nsign],weight[nsign]);
            h_eta1_eta2_mu4_mu4noL1_w_single_b_sig_sel->Fill(m2eta[nsign],m1eta[nsign],weight[nsign]);
            h_eta_avg_Deta_mu4_mu4noL1_w_single_b_sig_sel->Fill(deta[nsign],etaavg[nsign],weight[nsign]);
            h_eta_avg_Dphi_mu4_mu4noL1_w_single_b_sig_sel->Fill(dphi[nsign],etaavg[nsign],weight[nsign]);
            h_minv_zoomin_mu4_mu4noL1_w_single_b_sig_sel->Fill(minv[nsign],weight[nsign]);
            h_pair_pt_log_mu4_mu4noL1_w_single_b_sig_sel->Fill(pair_pt[nsign],weight[nsign]);
            h_minv_pair_pt_log_mu4_mu4noL1_w_single_b_sig_sel->Fill(pair_pt[nsign],minv[nsign],weight[nsign]);

            if ((!mu1PassSingle[nsign] || !mu2PassSingle[nsign])){ // exclusive: only one muon passes mu4
                h_Deta_mu4_mu4noL1_excl_w_single_b_sig_sel->Fill(deta[nsign],weight[nsign]);
                h_Deta_zoomin_mu4_mu4noL1_excl_w_single_b_sig_sel->Fill(deta[nsign],weight[nsign]);
                h_Dphi_mu4_mu4noL1_excl_w_single_b_sig_sel->Fill(dphi[nsign],weight[nsign]);
                h_Dphi_zoomin_mu4_mu4noL1_excl_w_single_b_sig_sel->Fill(dphi[nsign],weight[nsign]);
                h_DR_mu4_mu4noL1_excl_w_single_b_sig_sel->Fill(dr[nsign],weight[nsign]);
                h_DR_zoomin_mu4_mu4noL1_excl_w_single_b_sig_sel->Fill(dr[nsign],weight[nsign]);
                if (mu1PassSingle[nsign]){
                    h_pt2nd_mu4_mu4noL1_excl_w_single_b_sig_sel->Fill(m2pt[nsign],weight[nsign]);
                    h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_excl_w_single_b_sig_sel->Fill(m2charge[nsign] * m2eta[nsign], m2pt[nsign],weight[nsign]);
                }
                if (mu2PassSingle[nsign]){
                    h_pt2nd_mu4_mu4noL1_excl_w_single_b_sig_sel->Fill(m1pt[nsign],weight[nsign]);
                    h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_excl_w_single_b_sig_sel->Fill(m1charge[nsign] * m1eta[nsign], m1pt[nsign],weight[nsign]);
                }

                h_pair_eta_vs_pair_pT_mu4_mu4noL1_excl_w_single_b_sig_sel->Fill(pair_pt[nsign],pair_eta[nsign],weight[nsign]);
                h_Deta_Dphi_mu4_mu4noL1_excl_w_single_b_sig_sel->Fill(dphi[nsign],deta[nsign],weight[nsign]);
                h_eta1_eta2_mu4_mu4noL1_excl_w_single_b_sig_sel->Fill(m2eta[nsign],m1eta[nsign],weight[nsign]);
                h_eta_avg_Deta_mu4_mu4noL1_excl_w_single_b_sig_sel->Fill(deta[nsign],etaavg[nsign],weight[nsign]);
                h_eta_avg_Dphi_mu4_mu4noL1_excl_w_single_b_sig_sel->Fill(dphi[nsign],etaavg[nsign],weight[nsign]);
                h_minv_zoomin_mu4_mu4noL1_excl_w_single_b_sig_sel->Fill(minv[nsign],weight[nsign]);
                h_pair_pt_log_mu4_mu4noL1_excl_w_single_b_sig_sel->Fill(pair_pt[nsign],weight[nsign]);
                h_minv_pair_pt_log_mu4_mu4noL1_excl_w_single_b_sig_sel->Fill(pair_pt[nsign],minv[nsign],weight[nsign]);
            }
        }
    }
    if (pass2mu4[nsign]){
        h_Deta_2mu4[nsign]->Fill(deta[nsign],weight[nsign]);
        h_Deta_zoomin_2mu4[nsign]->Fill(deta[nsign],weight[nsign]);
        h_Dphi_2mu4[nsign]->Fill(dphi[nsign],weight[nsign]);
        h_Dphi_zoomin_2mu4[nsign]->Fill(dphi[nsign],weight[nsign]);
        h_DR_2mu4[nsign]->Fill(dr[nsign],weight[nsign]);
        h_DR_zoomin_2mu4[nsign]->Fill(dr[nsign],weight[nsign]);
        h_DR_0_2_2mu4[nsign]->Fill(dr[nsign],weight[nsign]);

        h_Deta_vs_pair_pT_2mu4[nsign]->Fill(pair_pt[nsign], deta[nsign],weight[nsign]);
        h_Deta_zoomin_vs_pair_pT_2mu4[nsign]->Fill(pair_pt[nsign], deta[nsign],weight[nsign]);
        h_Dphi_vs_pair_pT_2mu4[nsign]->Fill(pair_pt[nsign], dphi[nsign],weight[nsign]);
        h_Dphi_zoomin_vs_pair_pT_2mu4[nsign]->Fill(pair_pt[nsign], dphi[nsign],weight[nsign]);
        h_DR_vs_pair_pT_2mu4[nsign]->Fill(pair_pt[nsign], dr[nsign],weight[nsign]);
        h_DR_zoomin_vs_pair_pT_2mu4[nsign]->Fill(pair_pt[nsign], dr[nsign],weight[nsign]);

        if (mu1PassSingle[nsign]){
            h_pt2nd_2mu4[nsign]->Fill(m2pt[nsign],weight[nsign]);
            h_pt2nd_vs_q_eta_2nd_2mu4[nsign]->Fill(m2charge[nsign] * m2eta[nsign], m2pt[nsign],weight[nsign]);
            h_pt2nd_vs_phi2nd_2mu4[nsign]->Fill(m2phi[nsign], m2pt[nsign], weight[nsign]);
            h_phi2nd_vs_q_eta_2nd_2mu4[nsign]->Fill(m2charge[nsign]*m2eta[nsign], m2phi[nsign], weight[nsign]);
            h_DR_zoomin_vs_pt2nd_2mu4[nsign]->Fill(m2pt[nsign],dr[nsign],weight[nsign]);
            h_DR_0_2_vs_pt2nd_2mu4[nsign]->Fill(m2pt[nsign],dr[nsign],weight[nsign]);

            h_Deta_vs_pT_1st_2mu4[nsign]->Fill(m2pt[nsign], deta[nsign],weight[nsign]);
            h_Deta_zoomin_vs_pT_1st_2mu4[nsign]->Fill(m2pt[nsign], deta[nsign],weight[nsign]);
            h_Dphi_vs_pT_1st_2mu4[nsign]->Fill(m2pt[nsign], dphi[nsign],weight[nsign]);
            h_Dphi_zoomin_vs_pT_1st_2mu4[nsign]->Fill(m2pt[nsign], dphi[nsign],weight[nsign]);
            h_DR_vs_pT_1st_2mu4[nsign]->Fill(m2pt[nsign], dr[nsign],weight[nsign]);
            h_DR_zoomin_vs_pT_1st_2mu4[nsign]->Fill(m2pt[nsign], dr[nsign],weight[nsign]);
        }
        if (mu2PassSingle[nsign]){
            h_pt2nd_2mu4[nsign]->Fill(m1pt[nsign],weight[nsign]);
            h_pt2nd_vs_q_eta_2nd_2mu4[nsign]->Fill(m1charge[nsign] * m1eta[nsign], m1pt[nsign],weight[nsign]);
            h_pt2nd_vs_phi2nd_2mu4[nsign]->Fill(m1phi[nsign], m1pt[nsign], weight[nsign]);
            h_phi2nd_vs_q_eta_2nd_2mu4[nsign]->Fill(m1charge[nsign]*m1eta[nsign], m1phi[nsign], weight[nsign]);
            h_DR_zoomin_vs_pt2nd_2mu4[nsign]->Fill(m1pt[nsign],dr[nsign],weight[nsign]);
            h_DR_0_2_vs_pt2nd_2mu4[nsign]->Fill(m1pt[nsign],dr[nsign],weight[nsign]);

            h_Deta_vs_pT_1st_2mu4[nsign]->Fill(m1pt[nsign], deta[nsign],weight[nsign]);
            h_Deta_zoomin_vs_pT_1st_2mu4[nsign]->Fill(m1pt[nsign], deta[nsign],weight[nsign]);
            h_Dphi_vs_pT_1st_2mu4[nsign]->Fill(m1pt[nsign], dphi[nsign],weight[nsign]);
            h_Dphi_zoomin_vs_pT_1st_2mu4[nsign]->Fill(m1pt[nsign], dphi[nsign],weight[nsign]);
            h_DR_vs_pT_1st_2mu4[nsign]->Fill(m1pt[nsign], dr[nsign],weight[nsign]);
            h_DR_zoomin_vs_pT_1st_2mu4[nsign]->Fill(m1pt[nsign], dr[nsign],weight[nsign]);

        }

        h_pair_eta_vs_pair_pT_2mu4[nsign]->Fill(pair_pt[nsign],pair_eta[nsign],weight[nsign]);
        h_Deta_Dphi_2mu4[nsign]->Fill(dphi[nsign],deta[nsign],weight[nsign]);
        h_eta1_eta2_2mu4[nsign]->Fill(m2eta[nsign],m1eta[nsign],weight[nsign]);
        h_eta_avg_Deta_2mu4[nsign]->Fill(deta[nsign],etaavg[nsign],weight[nsign]);
        h_eta_avg_Dphi_2mu4[nsign]->Fill(dphi[nsign],etaavg[nsign],weight[nsign]);
        h_minv_zoomin_2mu4[nsign]->Fill(minv[nsign],weight[nsign]);
        h_pair_pt_log_2mu4[nsign]->Fill(pair_pt[nsign],weight[nsign]);
        h_minv_pair_pt_log_2mu4[nsign]->Fill(pair_pt[nsign],minv[nsign],weight[nsign]);

        for (int i = 0; i < ParamsSet::nSigns; i++){
            if (pass_single_muon_good_acceptance_selection[i]){ // fill for both
                h_Deta_2mu4_good_accept[nsign]->Fill(deta[nsign],weight[nsign]);
                h_Deta_zoomin_2mu4_good_accept[nsign]->Fill(deta[nsign],weight[nsign]);
                h_Dphi_2mu4_good_accept[nsign]->Fill(dphi[nsign],weight[nsign]);
                h_Dphi_zoomin_2mu4_good_accept[nsign]->Fill(dphi[nsign],weight[nsign]);
                h_DR_2mu4_good_accept[nsign]->Fill(dr[nsign],weight[nsign]);
                h_DR_zoomin_2mu4_good_accept[nsign]->Fill(dr[nsign],weight[nsign]);
                h_minv_zoomin_2mu4_good_accept[nsign]->Fill(minv[nsign],weight[nsign]);
                h_pair_pt_log_2mu4_good_accept[nsign]->Fill(pair_pt[nsign],weight[nsign]);    
            }
        }

        if (passSeparated[nsign]){
            h_Deta_2mu4_sepr[nsign]->Fill(deta[nsign],weight[nsign]);
            h_Deta_zoomin_2mu4_sepr[nsign]->Fill(deta[nsign],weight[nsign]);
            h_Deta_Dphi_2mu4_sepr[nsign]->Fill(dphi[nsign],deta[nsign],weight[nsign]);
            h_eta_avg_Deta_2mu4_sepr[nsign]->Fill(deta[nsign],etaavg[nsign],weight[nsign]);
        }

        if (passSeparated[nsign]){
            h_Dphi_2mu4_sepr[nsign]->Fill(dphi[nsign],weight[nsign]);
            h_Dphi_zoomin_2mu4_sepr[nsign]->Fill(dphi[nsign],weight[nsign]);
            h_DR_2mu4_sepr[nsign]->Fill(dr[nsign],weight[nsign]);
            h_DR_zoomin_2mu4_sepr[nsign]->Fill(dr[nsign],weight[nsign]);
            if (mu1PassSingle[nsign]){
                h_pt2nd_2mu4_sepr[nsign]->Fill(m2pt[nsign],weight[nsign]);
                h_pt2nd_vs_q_eta_2nd_2mu4_sepr[nsign]->Fill(m2charge[nsign] * m2eta[nsign], m2pt[nsign],weight[nsign]);
            }
            if (mu2PassSingle[nsign]){
                h_pt2nd_2mu4_sepr[nsign]->Fill(m1pt[nsign],weight[nsign]);
                h_pt2nd_vs_q_eta_2nd_2mu4_sepr[nsign]->Fill(m1charge[nsign] * m1eta[nsign], m1pt[nsign],weight[nsign]);
            }
    
            h_pair_eta_vs_pair_pT_2mu4_sepr[nsign]->Fill(pair_pt[nsign],pair_eta[nsign],weight[nsign]);
            h_eta1_eta2_2mu4_sepr[nsign]->Fill(m2eta[nsign],m1eta[nsign],weight[nsign]);
            h_eta_avg_Dphi_2mu4_sepr[nsign]->Fill(dphi[nsign],etaavg[nsign],weight[nsign]);
            h_minv_zoomin_2mu4_sepr[nsign]->Fill(minv[nsign],weight[nsign]);
            h_pair_pt_log_2mu4_sepr[nsign]->Fill(pair_pt[nsign],weight[nsign]);
            h_minv_pair_pt_log_2mu4_sepr[nsign]->Fill(pair_pt[nsign],minv[nsign],weight[nsign]);
        }

        if (pass_single_b_signal_selection){
            h_Deta_2mu4_w_single_b_sig_sel->Fill(deta[nsign],weight[nsign]);
            h_Deta_zoomin_2mu4_w_single_b_sig_sel->Fill(deta[nsign],weight[nsign]);
            h_Dphi_2mu4_w_single_b_sig_sel->Fill(dphi[nsign],weight[nsign]);
            h_Dphi_zoomin_2mu4_w_single_b_sig_sel->Fill(dphi[nsign],weight[nsign]);
            h_DR_2mu4_w_single_b_sig_sel->Fill(dr[nsign],weight[nsign]);
            h_DR_zoomin_2mu4_w_single_b_sig_sel->Fill(dr[nsign],weight[nsign]);
            h_DR_0_2_2mu4_w_single_b_sig_sel->Fill(dr[nsign],weight[nsign]);

            if (mu1PassSingle[nsign]){
                h_pt2nd_2mu4_w_single_b_sig_sel->Fill(m2pt[nsign],weight[nsign]);
                h_pt2nd_vs_q_eta_2nd_2mu4_w_single_b_sig_sel->Fill(m2charge[nsign] * m2eta[nsign], m2pt[nsign],weight[nsign]);
                h_DR_zoomin_vs_pt2nd_2mu4_w_single_b_sig_sel->Fill(m2pt[nsign],dr[nsign],weight[nsign]);
                h_DR_0_2_vs_pt2nd_2mu4_w_single_b_sig_sel->Fill(m2pt[nsign],dr[nsign],weight[nsign]);
            }
            if (mu2PassSingle[nsign]){
                h_pt2nd_2mu4_w_single_b_sig_sel->Fill(m1pt[nsign],weight[nsign]);
                h_pt2nd_vs_q_eta_2nd_2mu4_w_single_b_sig_sel->Fill(m1charge[nsign] * m1eta[nsign], m1pt[nsign],weight[nsign]);
                h_DR_zoomin_vs_pt2nd_2mu4_w_single_b_sig_sel->Fill(m1pt[nsign],dr[nsign],weight[nsign]);
                h_DR_0_2_vs_pt2nd_2mu4_w_single_b_sig_sel->Fill(m1pt[nsign],dr[nsign],weight[nsign]);
            }

            h_pair_eta_vs_pair_pT_2mu4_w_single_b_sig_sel->Fill(pair_pt[nsign],pair_eta[nsign],weight[nsign]);
            h_Deta_Dphi_2mu4_w_single_b_sig_sel->Fill(dphi[nsign],deta[nsign],weight[nsign]);
            h_eta1_eta2_2mu4_w_single_b_sig_sel->Fill(m2eta[nsign],m1eta[nsign],weight[nsign]);
            h_eta_avg_Deta_2mu4_w_single_b_sig_sel->Fill(deta[nsign],etaavg[nsign],weight[nsign]);
            h_eta_avg_Dphi_2mu4_w_single_b_sig_sel->Fill(dphi[nsign],etaavg[nsign],weight[nsign]);
            h_minv_zoomin_2mu4_w_single_b_sig_sel->Fill(minv[nsign],weight[nsign]);
            h_pair_pt_log_2mu4_w_single_b_sig_sel->Fill(pair_pt[nsign],weight[nsign]);
            h_minv_pair_pt_log_2mu4_w_single_b_sig_sel->Fill(pair_pt[nsign],minv[nsign],weight[nsign]);
        }
    }

}


void MuonPairPlottingPP::FillTrigEffcyHistsInvWeightedbySingleMuonEffcies(int nsign){
    bool pass_single_b_signal_selection = (nsign == 1 && minv[nsign] > 1.08 && minv[nsign] < 2.9 && pair_pt[nsign] > 8);

    double single_muon_pt_max = h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sepr_divided[nsign]->GetYaxis()->GetXmax(); // max value single-muon pT distribution reaches: pairs with 2nd muon pT above this value will not appear in the weighted trigger efficiency histogramms

    double pt_2nd = mu1PassSingle[nsign]? m2pt[nsign] : m1pt[nsign];
    if (pt_2nd >= single_muon_pt_max) return; // return without filling the inversed-weighted-by-single-muon-efficiency histograms

    double q_eta_2nd = mu1PassSingle[nsign]? m2charge[nsign] * m2eta[nsign] : m1charge[nsign] * m1eta[nsign];

    int bin_num = h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sepr_divided[nsign]->FindBin(q_eta_2nd, pt_2nd);

    if (passmu4mu4noL1[nsign]){
        h_Deta_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy[nsign]->Fill(deta[nsign],weight[nsign] * 1. / h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sepr_divided[nsign]->GetBinContent(bin_num));
        h_Dphi_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy[nsign]->Fill(dphi[nsign],weight[nsign] * 1. / h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sepr_divided[nsign]->GetBinContent(bin_num));
        h_DR_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy[nsign]->Fill(dr[nsign],weight[nsign] * 1. / h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sepr_divided[nsign]->GetBinContent(bin_num));
        h_DR_0_2_mu4_mu4noL1_inv_w_by_single_mu_effcy[nsign]->Fill(dr[nsign],weight[nsign] * 1. / h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sepr_divided[nsign]->GetBinContent(bin_num));
        h_minv_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy[nsign]->Fill(minv[nsign],weight[nsign] * 1. / h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sepr_divided[nsign]->GetBinContent(bin_num));
        h_pair_pt_log_mu4_mu4noL1_inv_w_by_single_mu_effcy[nsign]->Fill(pair_pt[nsign],weight[nsign] * 1. / h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sepr_divided[nsign]->GetBinContent(bin_num));
        h_pt2nd_mu4_mu4noL1_inv_w_by_single_mu_effcy[nsign]->Fill(pt_2nd,weight[nsign] * 1. / h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sepr_divided[nsign]->GetBinContent(bin_num));

        if (pass_single_b_signal_selection){
            h_Deta_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_w_single_b_sig_sel->Fill(deta[nsign],weight[nsign] * 1. / h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sepr_divided_w_single_b_sig_sel->GetBinContent(bin_num));
            h_Dphi_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_w_single_b_sig_sel->Fill(dphi[nsign],weight[nsign] * 1. / h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sepr_divided_w_single_b_sig_sel->GetBinContent(bin_num));
            h_DR_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_w_single_b_sig_sel->Fill(dr[nsign],weight[nsign] * 1. / h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sepr_divided_w_single_b_sig_sel->GetBinContent(bin_num));
            h_DR_0_2_mu4_mu4noL1_inv_w_by_single_mu_effcy_w_single_b_sig_sel->Fill(dr[nsign],weight[nsign] * 1. / h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sepr_divided_w_single_b_sig_sel->GetBinContent(bin_num));
            h_minv_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_w_single_b_sig_sel->Fill(minv[nsign],weight[nsign] * 1. / h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sepr_divided_w_single_b_sig_sel->GetBinContent(bin_num));
            h_pair_pt_log_mu4_mu4noL1_inv_w_by_single_mu_effcy_w_single_b_sig_sel->Fill(pair_pt[nsign],weight[nsign] * 1. / h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sepr_divided_w_single_b_sig_sel->GetBinContent(bin_num));
            h_pt2nd_mu4_mu4noL1_inv_w_by_single_mu_effcy_w_single_b_sig_sel->Fill(pt_2nd,weight[nsign] * 1. / h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sepr_divided_w_single_b_sig_sel->GetBinContent(bin_num));
        }
    }

    if (pass2mu4[nsign]){
        h_Deta_zoomin_2mu4_inv_w_by_single_mu_effcy[nsign]->Fill(deta[nsign],weight[nsign] * 1. / h_pt2nd_vs_q_eta_2nd_2mu4_sepr_divided[nsign]->GetBinContent(bin_num));
        h_Dphi_zoomin_2mu4_inv_w_by_single_mu_effcy[nsign]->Fill(dphi[nsign],weight[nsign] * 1. / h_pt2nd_vs_q_eta_2nd_2mu4_sepr_divided[nsign]->GetBinContent(bin_num));
        h_DR_zoomin_2mu4_inv_w_by_single_mu_effcy[nsign]->Fill(dr[nsign],weight[nsign] * 1. / h_pt2nd_vs_q_eta_2nd_2mu4_sepr_divided[nsign]->GetBinContent(bin_num));
        h_DR_0_2_2mu4_inv_w_by_single_mu_effcy[nsign]->Fill(dr[nsign],weight[nsign] * 1. / h_pt2nd_vs_q_eta_2nd_2mu4_sepr_divided[nsign]->GetBinContent(bin_num));
        h_minv_zoomin_2mu4_inv_w_by_single_mu_effcy[nsign]->Fill(minv[nsign],weight[nsign] * 1. / h_pt2nd_vs_q_eta_2nd_2mu4_sepr_divided[nsign]->GetBinContent(bin_num));
        h_pair_pt_log_2mu4_inv_w_by_single_mu_effcy[nsign]->Fill(pair_pt[nsign],weight[nsign] * 1. / h_pt2nd_vs_q_eta_2nd_2mu4_sepr_divided[nsign]->GetBinContent(bin_num));
        h_pt2nd_2mu4_inv_w_by_single_mu_effcy[nsign]->Fill(pt_2nd,weight[nsign] * 1. / h_pt2nd_vs_q_eta_2nd_2mu4_sepr_divided[nsign]->GetBinContent(bin_num));

        if (pass_single_b_signal_selection){
            h_Deta_zoomin_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel->Fill(deta[nsign],weight[nsign] * 1. / h_pt2nd_vs_q_eta_2nd_2mu4_sepr_divided_w_single_b_sig_sel->GetBinContent(bin_num));
            h_Dphi_zoomin_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel->Fill(dphi[nsign],weight[nsign] * 1. / h_pt2nd_vs_q_eta_2nd_2mu4_sepr_divided_w_single_b_sig_sel->GetBinContent(bin_num));
            h_DR_zoomin_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel->Fill(dr[nsign],weight[nsign] * 1. / h_pt2nd_vs_q_eta_2nd_2mu4_sepr_divided_w_single_b_sig_sel->GetBinContent(bin_num));
            h_DR_0_2_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel->Fill(dr[nsign],weight[nsign] * 1. / h_pt2nd_vs_q_eta_2nd_2mu4_sepr_divided_w_single_b_sig_sel->GetBinContent(bin_num));
            h_minv_zoomin_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel->Fill(minv[nsign],weight[nsign] * 1. / h_pt2nd_vs_q_eta_2nd_2mu4_sepr_divided_w_single_b_sig_sel->GetBinContent(bin_num));
            h_pair_pt_log_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel->Fill(pair_pt[nsign],weight[nsign] * 1. / h_pt2nd_vs_q_eta_2nd_2mu4_sepr_divided_w_single_b_sig_sel->GetBinContent(bin_num));
            h_pt2nd_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel->Fill(pt_2nd,weight[nsign] * 1. / h_pt2nd_vs_q_eta_2nd_2mu4_sepr_divided_w_single_b_sig_sel->GetBinContent(bin_num));
        }
    }
}


void MuonPairPlottingPP::CalculateSingleMuonTrigEffcyRatios(){
    for (int isign = 0; isign < ParamsSet::nSigns; isign++){
        h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sepr_divided[isign] = (TH2D*)h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sepr[isign]->Clone(Form("%s_clone", h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sepr[isign]->GetName()));
        h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sepr_divided[isign]->Divide(h_pt2nd_vs_q_eta_2nd_mu4_sepr[isign]);
        h_pt2nd_vs_q_eta_2nd_2mu4_sepr_divided[isign] = (TH2D*)h_pt2nd_vs_q_eta_2nd_2mu4_sepr[isign]->Clone(Form("%s_clone", h_pt2nd_vs_q_eta_2nd_2mu4_sepr[isign]->GetName()));
        h_pt2nd_vs_q_eta_2nd_2mu4_sepr_divided[isign]->Divide(h_pt2nd_vs_q_eta_2nd_mu4_sepr[isign]);
    }

    h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sepr_divided_w_single_b_sig_sel = (TH2D*)h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_w_single_b_sig_sel->Clone(Form("%s_clone", h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_w_single_b_sig_sel->GetName()));
    h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sepr_divided_w_single_b_sig_sel->Divide(h_pt2nd_vs_q_eta_2nd_mu4_w_single_b_sig_sel);
    h_pt2nd_vs_q_eta_2nd_2mu4_sepr_divided_w_single_b_sig_sel = (TH2D*)h_pt2nd_vs_q_eta_2nd_2mu4_w_single_b_sig_sel->Clone(Form("%s_clone", h_pt2nd_vs_q_eta_2nd_2mu4_w_single_b_sig_sel->GetName()));
    h_pt2nd_vs_q_eta_2nd_2mu4_sepr_divided_w_single_b_sig_sel->Divide(h_pt2nd_vs_q_eta_2nd_mu4_w_single_b_sig_sel);
}

void MuonPairPlottingPP::WriteOutput(){

    outFile->cd();
    outFile->Write();
    outFile->Close();
}

void MuonPairPlottingPP::Run(){
    clock_t start, end;
    double cpu_time_used;
    start = clock();

    if (trigger_mode != 0 && trigger_mode != 1 && doTrigEffcy){
        std::cerr << "If doTrigEffcy, trigger_mode must be 1 (single-muon trigger selection) or 0 (minbias trigger selection)!!" << std::endl;
        throw std::exception();
    }

    if (isScram) std::cout << "Data being processed: Scrambled" << std::endl;
    else if (isTight) std::cout << "Data being processed: Tight" << std::endl;
    else std::cout << "Data being processed: Real" << std::endl;

    switch(trigger_mode){
    case 0:
        trig_suffix = "_min_bias";
        break;
    case 1:
        trig_suffix = "_single_mu4";
        break;
    case 2:
        trig_suffix = "_mu4_mu4noL1";
        break;
    case 3:
        trig_suffix = "_2mu4";
        break;        
    default:
        std::cerr << "Trigger mode INVALID: must be 0 / 1 / 2 / 3!" << std::endl;
        throw std::exception();
    }

    // if (doTrigEffcy && !filter_out_photo_resn_for_trig_effcy) trig_suffix += "_no_photo_resn_cuts"; // if not filter out photoprod/resn pairs for trigger efficiency study
    if (!filter_out_photo_resn_for_trig_effcy) trig_suffix += "_no_photo_resn_cuts"; // if not filter out photoprod/resn pairs for trigger efficiency study

    output_non_trig_effcy_hists = !(trigger_mode == 0 || trigger_mode == 1);

  	InitInput();
    InitOutput();
  	InitHists();
  	ProcessData();
    // CalculateTrigEffcyRatio();
  	WriteOutput();

    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;;
    std::cout << "CPU time used is " << cpu_time_used << " seconds" << std::endl;

}