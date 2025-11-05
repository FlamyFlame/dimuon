#include "MuonPairPlottingPP.h"
#include "time.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"

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
        // inTree[isign]->SetBranchAddress("m1.passmu4"           , &mu1PassSingle[isign]);
        // inTree[isign]->SetBranchAddress("m2.passmu4"           , &mu2PassSingle[isign]);
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
        // inTree[isign]->SetBranchStatus("m1.passmu4"           ,1);
        // inTree[isign]->SetBranchStatus("m2.passmu4"           ,1);
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

    std::vector<double> eta_bins_trig_effcy = ParamsSet::makeEtaTrigEffcyBinning();
    int neta_bins_trig_effcy = eta_bins_trig_effcy.size() - 1;
    int nphi_bins_trig_effcy = 128; // phi 2nd muon
    int nDphi_bins_trig_effcy = 64;
    int nDeta_bins_trig_effcy = 40;
    int nDR_bins_trig_effcy = 40;
    int nminv_bins_trig_effcy = 40;
    int nDR_deta_dphi_zoomin_bins_trig_effcy = 20;

    // Build uniform phi edges so we can use the (xbins, ybins, zbins) TH3D ctor
    std::vector<double> phi2nd_bins(nphi_bins_trig_effcy + 1);
    for (int i = 0; i <= nphi_bins_trig_effcy; ++i) {
        phi2nd_bins[i] = -pms.PI + (2.0 * pms.PI) * (static_cast<double>(i) / nphi_bins_trig_effcy);
    }

    std::vector<double> pT_bins_single_muon (pms.pT_bins_8); // make a copy of a suitable set of single-muon pT bins (adjustable) --> use the copy for histogram settings

    pT_bins_single_muon.insert(pT_bins_single_muon.end(), pms.pT_bins_60.begin(), pms.pT_bins_60.end());

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

    // turn on Sumw2 for all histograms
    TH1::SetDefaultSumw2(kTRUE);

    h_crossx_dR_cut_vs_pair_pt_pair_eta = new TH2D("h_crossx_dR_cut_vs_pair_pt_pair_eta",";#eta^{pair};p_{T}^{pair} [GeV];#sigma^{truth}",npair_eta_bins_coarse,pair_eta_min,pair_eta_max,npair_pT_bins_coarse,pair_pt_min,pair_pt_max);

    if (trigger_mode == 1){
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
        h_pair_eta_vs_pair_pT_mu4_w_single_b_sig_sel = new TH2D("h_pair_eta_vs_pair_pT_mu4_w_single_b_sig_sel",";p_{T}^{pair} [GeV];#eta^{pair}",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),neta_bins_trig_effcy,eta_bins_trig_effcy.data());
        h_Deta_Dphi_mu4_w_single_b_sig_sel = new TH2D("h_Deta_Dphi_mu4_w_single_b_sig_sel",";#Delta#phi;#Delta#eta", nDphi_bins_trig_effcy,-pms.PI,pms.PI,nDeta_bins_trig_effcy,-4.8,4.8);
        h_eta1_eta2_mu4_w_single_b_sig_sel = new TH2D("h_eta1_eta2_mu4_w_single_b_sig_sel",";#eta_{sublead};#eta_{lead}",neta_bins_trig_effcy,eta_bins_trig_effcy.data(), neta_bins_trig_effcy,eta_bins_trig_effcy.data());
        h_eta_avg_Deta_mu4_w_single_b_sig_sel = new TH2D("h_eta_avg_Deta_mu4_w_single_b_sig_sel",";#bar{#eta};#Delta#eta", neta_bins_trig_effcy,eta_bins_trig_effcy.data(),nDeta_bins_trig_effcy,-4.8,4.8);
        h_eta_avg_Dphi_mu4_w_single_b_sig_sel = new TH2D("h_eta_avg_Dphi_mu4_w_single_b_sig_sel",";#Delta#phi;#bar{#eta}", nDphi_bins_trig_effcy,-pms.PI,pms.PI,neta_bins_trig_effcy,eta_bins_trig_effcy.data());
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
        h_pair_eta_vs_pair_pT_mu4_mu4noL1_w_single_b_sig_sel = new TH2D("h_pair_eta_vs_pair_pT_mu4_mu4noL1_w_single_b_sig_sel",";p_{T}^{pair} [GeV];#eta^{pair}",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),neta_bins_trig_effcy,eta_bins_trig_effcy.data());
        h_Deta_Dphi_mu4_mu4noL1_w_single_b_sig_sel = new TH2D("h_Deta_Dphi_mu4_mu4noL1_w_single_b_sig_sel",";#Delta#phi;#Delta#eta", nDphi_bins_trig_effcy,-pms.PI,pms.PI,nDeta_bins_trig_effcy,-4.8,4.8);
        h_eta1_eta2_mu4_mu4noL1_w_single_b_sig_sel = new TH2D("h_eta1_eta2_mu4_mu4noL1_w_single_b_sig_sel",";#eta_{sublead};#eta_{lead}",neta_bins_trig_effcy,eta_bins_trig_effcy.data(), neta_bins_trig_effcy,eta_bins_trig_effcy.data());
        h_eta_avg_Deta_mu4_mu4noL1_w_single_b_sig_sel = new TH2D("h_eta_avg_Deta_mu4_mu4noL1_w_single_b_sig_sel",";#bar{#eta};#Delta#eta", neta_bins_trig_effcy,eta_bins_trig_effcy.data(),nDeta_bins_trig_effcy,-4.8,4.8);
        h_eta_avg_Dphi_mu4_mu4noL1_w_single_b_sig_sel = new TH2D("h_eta_avg_Dphi_mu4_mu4noL1_w_single_b_sig_sel",";#Delta#phi;#bar{#eta}", nDphi_bins_trig_effcy,-pms.PI,pms.PI,neta_bins_trig_effcy,eta_bins_trig_effcy.data());
        h_minv_pair_pt_log_mu4_mu4noL1_w_single_b_sig_sel = new TH2D("h_minv_pair_pt_log_mu4_mu4noL1_w_single_b_sig_sel",";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),nminv_bins_log,minv_bins_log[1]);

        h_Deta_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_w_single_b_sig_sel = new TH1D("h_Deta_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_w_single_b_sig_sel",";#Delta#eta;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
        h_Dphi_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_w_single_b_sig_sel = new TH1D("h_Dphi_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_w_single_b_sig_sel",";#Delta#phi;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
        h_DR_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_w_single_b_sig_sel = new TH1D("h_DR_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_w_single_b_sig_sel",";#DeltaR;",nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
        h_DR_0_2_mu4_mu4noL1_inv_w_by_single_mu_effcy_w_single_b_sig_sel = new TH1D("h_DR_0_2_mu4_mu4noL1_inv_w_by_single_mu_effcy_w_single_b_sig_sel",";#DeltaR;",20,0,0.2);
        h_minv_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_w_single_b_sig_sel = new TH1D("h_minv_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_w_single_b_sig_sel",";m_{#mu#mu} [GeV];",nminv_bins_trig_effcy,0,3.);
        h_pair_pt_log_mu4_mu4noL1_inv_w_by_single_mu_effcy_w_single_b_sig_sel = new TH1D("h_pair_pt_log_mu4_mu4noL1_inv_w_by_single_mu_effcy_w_single_b_sig_sel",";p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
        h_DR_zoomin_vs_pair_pt_log_mu4_mu4noL1_inv_w_by_single_mu_effcy_w_single_b_sig_sel = new TH2D("h_DR_zoomin_vs_pair_pt_log_mu4_mu4noL1_inv_w_by_single_mu_effcy_w_single_b_sig_sel",";p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);

        h_Deta_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom_w_single_b_sig_sel = new TH1D("h_Deta_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom_w_single_b_sig_sel",";#Delta#eta;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
        h_Dphi_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom_w_single_b_sig_sel = new TH1D("h_Dphi_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom_w_single_b_sig_sel",";#Delta#phi;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
        h_DR_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom_w_single_b_sig_sel = new TH1D("h_DR_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom_w_single_b_sig_sel",";#DeltaR;",nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
        h_DR_0_2_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom_w_single_b_sig_sel = new TH1D("h_DR_0_2_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom_w_single_b_sig_sel",";#DeltaR;",20,0,0.2);
        h_minv_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom_w_single_b_sig_sel = new TH1D("h_minv_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom_w_single_b_sig_sel",";m_{#mu#mu} [GeV];",nminv_bins_trig_effcy,0,3.);
        h_pair_pt_log_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom_w_single_b_sig_sel = new TH1D("h_pair_pt_log_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom_w_single_b_sig_sel",";p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
        h_DR_zoomin_vs_pair_pt_log_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom_w_single_b_sig_sel = new TH2D("h_DR_zoomin_vs_pair_pt_log_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom_w_single_b_sig_sel",";p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);




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
        h_pair_eta_vs_pair_pT_2mu4_w_single_b_sig_sel = new TH2D("h_pair_eta_vs_pair_pT_2mu4_w_single_b_sig_sel",";p_{T}^{pair} [GeV];#eta^{pair}",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),neta_bins_trig_effcy,eta_bins_trig_effcy.data());
        h_Deta_Dphi_2mu4_w_single_b_sig_sel = new TH2D("h_Deta_Dphi_2mu4_w_single_b_sig_sel",";#Delta#phi;#Delta#eta", nDphi_bins_trig_effcy,-pms.PI,pms.PI,nDeta_bins_trig_effcy,-4.8,4.8);
        h_eta1_eta2_2mu4_w_single_b_sig_sel = new TH2D("h_eta1_eta2_2mu4_w_single_b_sig_sel",";#eta_{sublead};#eta_{lead}",neta_bins_trig_effcy,eta_bins_trig_effcy.data(), neta_bins_trig_effcy,eta_bins_trig_effcy.data());
        h_eta_avg_Deta_2mu4_w_single_b_sig_sel = new TH2D("h_eta_avg_Deta_2mu4_w_single_b_sig_sel",";#bar{#eta};#Delta#eta", neta_bins_trig_effcy,eta_bins_trig_effcy.data(),nDeta_bins_trig_effcy,-4.8,4.8);
        h_eta_avg_Dphi_2mu4_w_single_b_sig_sel = new TH2D("h_eta_avg_Dphi_2mu4_w_single_b_sig_sel",";#Delta#phi;#bar{#eta}", nDphi_bins_trig_effcy,-pms.PI,pms.PI,neta_bins_trig_effcy,eta_bins_trig_effcy.data());
        h_minv_pair_pt_log_2mu4_w_single_b_sig_sel = new TH2D("h_minv_pair_pt_log_2mu4_w_single_b_sig_sel",";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),nminv_bins_log,minv_bins_log[1]);
        
        h_Deta_zoomin_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel = new TH1D("h_Deta_zoomin_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel",";#Delta#eta;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
        h_Dphi_zoomin_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel = new TH1D("h_Dphi_zoomin_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel",";#Delta#phi;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
        h_DR_zoomin_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel = new TH1D("h_DR_zoomin_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel",";#DeltaR;",nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
        h_DR_0_2_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel = new TH1D("h_DR_0_2_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel",";#DeltaR;",20,0,0.2);
        h_minv_zoomin_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel = new TH1D("h_minv_zoomin_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel",";m_{#mu#mu} [GeV];",nminv_bins_trig_effcy,0,3.);
        h_pair_pt_log_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel = new TH1D("h_pair_pt_log_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel",";p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
        h_DR_zoomin_vs_pair_pt_log_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel = new TH2D("h_DR_zoomin_vs_pair_pt_log_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel",";p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);

        h_Deta_zoomin_2mu4_inv_w_by_single_mu_effcy_denom_w_single_b_sig_sel = new TH1D("h_Deta_zoomin_2mu4_inv_w_by_single_mu_effcy_denom_w_single_b_sig_sel",";#Delta#eta;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
        h_Dphi_zoomin_2mu4_inv_w_by_single_mu_effcy_denom_w_single_b_sig_sel = new TH1D("h_Dphi_zoomin_2mu4_inv_w_by_single_mu_effcy_denom_w_single_b_sig_sel",";#Delta#phi;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
        h_DR_zoomin_2mu4_inv_w_by_single_mu_effcy_denom_w_single_b_sig_sel = new TH1D("h_DR_zoomin_2mu4_inv_w_by_single_mu_effcy_denom_w_single_b_sig_sel",";#DeltaR;",nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
        h_DR_0_2_2mu4_inv_w_by_single_mu_effcy_denom_w_single_b_sig_sel = new TH1D("h_DR_0_2_2mu4_inv_w_by_single_mu_effcy_denom_w_single_b_sig_sel",";#DeltaR;",20,0,0.2);
        h_minv_zoomin_2mu4_inv_w_by_single_mu_effcy_denom_w_single_b_sig_sel = new TH1D("h_minv_zoomin_2mu4_inv_w_by_single_mu_effcy_denom_w_single_b_sig_sel",";m_{#mu#mu} [GeV];",nminv_bins_trig_effcy,0,3.);
        h_pair_pt_log_2mu4_inv_w_by_single_mu_effcy_denom_w_single_b_sig_sel = new TH1D("h_pair_pt_log_2mu4_inv_w_by_single_mu_effcy_denom_w_single_b_sig_sel",";p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
        h_DR_zoomin_vs_pair_pt_log_2mu4_inv_w_by_single_mu_effcy_denom_w_single_b_sig_sel = new TH2D("h_DR_zoomin_vs_pair_pt_log_2mu4_inv_w_by_single_mu_effcy_denom_w_single_b_sig_sel",";p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);


        // --------------------- signed histograms below ---------------------------

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

            h_pair_eta_vs_pair_pT_mu4[isign] = new TH2D(Form("h_pair_eta_vs_pair_pT_mu4_sign%d",isign+1),";p_{T}^{pair} [GeV];#eta^{pair}",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),neta_bins_trig_effcy,eta_bins_trig_effcy.data());
            h_Deta_Dphi_mu4[isign] = new TH2D(Form("h_Deta_Dphi_mu4_sign%d",isign+1),";#Delta#phi;#Delta#eta", nDphi_bins_trig_effcy,-pms.PI,pms.PI,nDeta_bins_trig_effcy,-4.8,4.8);
            h_eta1_eta2_mu4[isign] = new TH2D(Form("h_eta1_eta2_mu4_sign%d",isign+1),";#eta_{sublead};#eta_{lead}",neta_bins_trig_effcy,eta_bins_trig_effcy.data(), neta_bins_trig_effcy,eta_bins_trig_effcy.data());
            h_eta_avg_Deta_mu4[isign] = new TH2D(Form("h_eta_avg_Deta_mu4_sign%d",isign+1),";#bar{#eta};#Delta#eta", neta_bins_trig_effcy,eta_bins_trig_effcy.data(),nDeta_bins_trig_effcy,-4.8,4.8);
            h_eta_avg_Dphi_mu4[isign] = new TH2D(Form("h_eta_avg_Dphi_mu4_sign%d",isign+1),";#Delta#phi;#bar{#eta}", nDphi_bins_trig_effcy,-pms.PI,pms.PI,neta_bins_trig_effcy,eta_bins_trig_effcy.data());
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

            h_pair_eta_vs_pair_pT_mu4_mu4noL1[isign] = new TH2D(Form("h_pair_eta_vs_pair_pT_mu4_mu4noL1_sign%d",isign+1),";p_{T}^{pair} [GeV];#eta^{pair}",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),neta_bins_trig_effcy,eta_bins_trig_effcy.data());
            h_Deta_Dphi_mu4_mu4noL1[isign] = new TH2D(Form("h_Deta_Dphi_mu4_mu4noL1_sign%d",isign+1),";#Delta#phi;#Delta#eta", nDphi_bins_trig_effcy,-pms.PI,pms.PI,nDeta_bins_trig_effcy,-4.8,4.8);
            h_eta1_eta2_mu4_mu4noL1[isign] = new TH2D(Form("h_eta1_eta2_mu4_mu4noL1_sign%d",isign+1),";#eta_{sublead};#eta_{lead}",neta_bins_trig_effcy,eta_bins_trig_effcy.data(), neta_bins_trig_effcy,eta_bins_trig_effcy.data());
            h_eta_avg_Deta_mu4_mu4noL1[isign] = new TH2D(Form("h_eta_avg_Deta_mu4_mu4noL1_sign%d",isign+1),";#bar{#eta};#Delta#eta", neta_bins_trig_effcy,eta_bins_trig_effcy.data(),nDeta_bins_trig_effcy,-4.8,4.8);
            h_eta_avg_Dphi_mu4_mu4noL1[isign] = new TH2D(Form("h_eta_avg_Dphi_mu4_mu4noL1_sign%d",isign+1),";#Delta#phi;#bar{#eta}", nDphi_bins_trig_effcy,-pms.PI,pms.PI,neta_bins_trig_effcy,eta_bins_trig_effcy.data());
            h_minv_pair_pt_log_mu4_mu4noL1[isign] = new TH2D(Form("h_minv_pair_pt_log_mu4_mu4noL1_sign%d",isign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),nminv_bins_log,minv_bins_log[isign]);

            h_Deta_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy[isign] = new TH1D(Form("h_Deta_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_sign%d",isign+1),";#Delta#eta;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
            h_Dphi_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy[isign] = new TH1D(Form("h_Dphi_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_sign%d",isign+1),";#Delta#phi;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
            h_DR_mu4_mu4noL1_inv_w_by_single_mu_effcy[isign] = new TH1D(Form("h_DR_mu4_mu4noL1_inv_w_by_single_mu_effcy_sign%d",isign+1),";#DeltaR;",nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
            h_DR_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy[isign] = new TH1D(Form("h_DR_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_sign%d",isign+1),";#DeltaR;",nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
            h_DR_0_2_mu4_mu4noL1_inv_w_by_single_mu_effcy[isign] = new TH1D(Form("h_DR_0_2_mu4_mu4noL1_inv_w_by_single_mu_effcy_sign%d",isign+1),";#DeltaR;",20,0,0.2);
            h_minv_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy[isign] = new TH1D(Form("h_minv_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_sign%d",isign+1),";m_{#mu#mu} [GeV];",nminv_bins_trig_effcy,0,3.);
            h_pair_pt_log_mu4_mu4noL1_inv_w_by_single_mu_effcy[isign] = new TH1D(Form("h_pair_pt_log_mu4_mu4noL1_inv_w_by_single_mu_effcy_sign%d",isign+1),";p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
            h_DR_zoomin_vs_pair_pt_log_mu4_mu4noL1_inv_w_by_single_mu_effcy[isign] = new TH2D(Form("h_DR_zoomin_vs_pair_pt_log_mu4_mu4noL1_inv_w_by_single_mu_effcy_sign%d",isign+1),";p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);

            h_Deta_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom[isign] = new TH1D(Form("h_Deta_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom_sign%d",isign+1),";#Delta#eta;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
            h_Dphi_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom[isign] = new TH1D(Form("h_Dphi_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom_sign%d",isign+1),";#Delta#phi;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
            h_DR_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom[isign] = new TH1D(Form("h_DR_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom_sign%d",isign+1),";#DeltaR;",nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
            h_DR_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom[isign] = new TH1D(Form("h_DR_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom_sign%d",isign+1),";#DeltaR;",nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
            h_DR_0_2_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom[isign] = new TH1D(Form("h_DR_0_2_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom_sign%d",isign+1),";#DeltaR;",20,0,0.2);
            h_minv_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom[isign] = new TH1D(Form("h_minv_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom_sign%d",isign+1),";m_{#mu#mu} [GeV];",nminv_bins_trig_effcy,0,3.);
            h_pair_pt_log_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom[isign] = new TH1D(Form("h_pair_pt_log_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom_sign%d",isign+1),";p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
            h_DR_zoomin_vs_pair_pt_log_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom[isign] = new TH2D(Form("h_DR_zoomin_vs_pair_pt_log_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom_sign%d",isign+1),";p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);




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

            h_pair_eta_vs_pair_pT_2mu4[isign] = new TH2D(Form("h_pair_eta_vs_pair_pT_2mu4_sign%d",isign+1),";p_{T}^{pair} [GeV];#eta^{pair}",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),neta_bins_trig_effcy,eta_bins_trig_effcy.data());
            h_Deta_Dphi_2mu4[isign] = new TH2D(Form("h_Deta_Dphi_2mu4_sign%d",isign+1),";#Delta#phi;#Delta#eta", nDphi_bins_trig_effcy,-pms.PI,pms.PI,nDeta_bins_trig_effcy,-4.8,4.8);
            h_eta1_eta2_2mu4[isign] = new TH2D(Form("h_eta1_eta2_2mu4_sign%d",isign+1),";#eta_{sublead};#eta_{lead}",neta_bins_trig_effcy,eta_bins_trig_effcy.data(), neta_bins_trig_effcy,eta_bins_trig_effcy.data());
            h_eta_avg_Deta_2mu4[isign] = new TH2D(Form("h_eta_avg_Deta_2mu4_sign%d",isign+1),";#bar{#eta};#Delta#eta", neta_bins_trig_effcy,eta_bins_trig_effcy.data(),nDeta_bins_trig_effcy,-4.8,4.8);
            h_eta_avg_Dphi_2mu4[isign] = new TH2D(Form("h_eta_avg_Dphi_2mu4_sign%d",isign+1),";#Delta#phi;#bar{#eta}", nDphi_bins_trig_effcy,-pms.PI,pms.PI,neta_bins_trig_effcy,eta_bins_trig_effcy.data());
            h_minv_pair_pt_log_2mu4[isign] = new TH2D(Form("h_minv_pair_pt_log_2mu4_sign%d",isign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),nminv_bins_log,minv_bins_log[isign]);

            h_Deta_zoomin_2mu4_inv_w_by_single_mu_effcy[isign] = new TH1D(Form("h_Deta_zoomin_2mu4_inv_w_by_single_mu_effcy_sign%d",isign+1),";#Delta#eta;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
            h_Dphi_zoomin_2mu4_inv_w_by_single_mu_effcy[isign] = new TH1D(Form("h_Dphi_zoomin_2mu4_inv_w_by_single_mu_effcy_sign%d",isign+1),";#Delta#phi;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
            h_DR_2mu4_inv_w_by_single_mu_effcy[isign] = new TH1D(Form("h_DR_2mu4_inv_w_by_single_mu_effcy_sign%d",isign+1),";#DeltaR;",nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
            h_DR_zoomin_2mu4_inv_w_by_single_mu_effcy[isign] = new TH1D(Form("h_DR_zoomin_2mu4_inv_w_by_single_mu_effcy_sign%d",isign+1),";#DeltaR;",nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
            h_DR_0_2_2mu4_inv_w_by_single_mu_effcy[isign] = new TH1D(Form("h_DR_0_2_2mu4_inv_w_by_single_mu_effcy_sign%d",isign+1),";#DeltaR;",20,0,0.2);
            h_minv_zoomin_2mu4_inv_w_by_single_mu_effcy[isign] = new TH1D(Form("h_minv_zoomin_2mu4_inv_w_by_single_mu_effcy_sign%d",isign+1),";m_{#mu#mu} [GeV];",nminv_bins_trig_effcy,0,3.);
            h_pair_pt_log_2mu4_inv_w_by_single_mu_effcy[isign] = new TH1D(Form("h_pair_pt_log_2mu4_inv_w_by_single_mu_effcy_sign%d",isign+1),";p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
            h_DR_zoomin_vs_pair_pt_log_2mu4_inv_w_by_single_mu_effcy[isign] = new TH2D(Form("h_DR_zoomin_vs_pair_pt_log_2mu4_inv_w_by_single_mu_effcy_sign%d",isign+1),";p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);

            h_Deta_zoomin_2mu4_inv_w_by_single_mu_effcy_denom[isign] = new TH1D(Form("h_Deta_zoomin_2mu4_inv_w_by_single_mu_effcy_denom_sign%d",isign+1),";#Delta#eta;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
            h_Dphi_zoomin_2mu4_inv_w_by_single_mu_effcy_denom[isign] = new TH1D(Form("h_Dphi_zoomin_2mu4_inv_w_by_single_mu_effcy_denom_sign%d",isign+1),";#Delta#phi;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
            h_DR_2mu4_inv_w_by_single_mu_effcy_denom[isign] = new TH1D(Form("h_DR_2mu4_inv_w_by_single_mu_effcy_denom_sign%d",isign+1),";#DeltaR;",nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
            h_DR_zoomin_2mu4_inv_w_by_single_mu_effcy_denom[isign] = new TH1D(Form("h_DR_zoomin_2mu4_inv_w_by_single_mu_effcy_denom_sign%d",isign+1),";#DeltaR;",nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
            h_DR_0_2_2mu4_inv_w_by_single_mu_effcy_denom[isign] = new TH1D(Form("h_DR_0_2_2mu4_inv_w_by_single_mu_effcy_denom_sign%d",isign+1),";#DeltaR;",20,0,0.2);
            h_minv_zoomin_2mu4_inv_w_by_single_mu_effcy_denom[isign] = new TH1D(Form("h_minv_zoomin_2mu4_inv_w_by_single_mu_effcy_denom_sign%d",isign+1),";m_{#mu#mu} [GeV];",nminv_bins_trig_effcy,0,3.);
            h_pair_pt_log_2mu4_inv_w_by_single_mu_effcy_denom[isign] = new TH1D(Form("h_pair_pt_log_2mu4_inv_w_by_single_mu_effcy_denom_sign%d",isign+1),";p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
            h_DR_zoomin_vs_pair_pt_log_2mu4_inv_w_by_single_mu_effcy_denom[isign] = new TH2D(Form("h_DR_zoomin_vs_pair_pt_log_2mu4_inv_w_by_single_mu_effcy_denom_sign%d",isign+1),";p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
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


            // --------------------- 2D 2nd-muon histograms for trigger efficiency ---------------------------
            // mu4
            h_pt2nd_vs_q_eta_2nd_mu4[isign] = new TH2D(Form("h_pt2nd_vs_q_eta_2nd_mu4_sign%d",isign+1),";q*#eta;p_{T} [GeV]",neta_bins_trig_effcy,eta_bins_trig_effcy.data(), int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
            h_pt2nd_vs_phi2nd_mu4[isign] = new TH2D(Form("h_pt2nd_vs_phi2nd_mu4_sign%d",isign+1),";#phi;p_{T} [GeV]",nphi_bins_trig_effcy,-pms.PI,pms.PI,int(pT_bins_single_muon.size())-1,pT_bins_single_muon.data());
            h_phi2nd_vs_q_eta_2nd_mu4[isign] =new TH2D(Form("h_phi2nd_vs_q_eta_2nd_mu4_sign%d",isign+1),";q*#eta;#phi",neta_bins_trig_effcy,eta_bins_trig_effcy.data(),nphi_bins_trig_effcy,-pms.PI,pms.PI);

            // mu4_mu4noL1
            h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1[isign] = new TH2D(Form("h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sign%d",isign+1),";q*#eta;p_{T} [GeV]",neta_bins_trig_effcy,eta_bins_trig_effcy.data(), int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
            h_pt2nd_vs_phi2nd_mu4_mu4noL1[isign] = new TH2D(Form("h_pt2nd_vs_phi2nd_mu4_mu4noL1_sign%d",isign+1),";#phi;p_{T} [GeV]",nphi_bins_trig_effcy,-pms.PI,pms.PI,int(pT_bins_single_muon.size())-1,pT_bins_single_muon.data());
            h_phi2nd_vs_q_eta_2nd_mu4_mu4noL1[isign] =new TH2D(Form("h_phi2nd_vs_q_eta_2nd_mu4_mu4noL1_sign%d",isign+1),";q*#eta;#phi",neta_bins_trig_effcy,eta_bins_trig_effcy.data(),nphi_bins_trig_effcy,-pms.PI,pms.PI);

            // 2mu4
            h_pt2nd_vs_q_eta_2nd_2mu4[isign] = new TH2D(Form("h_pt2nd_vs_q_eta_2nd_2mu4_sign%d",isign+1),";q*#eta;p_{T} [GeV]",neta_bins_trig_effcy,eta_bins_trig_effcy.data(), int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
            h_pt2nd_vs_phi2nd_2mu4[isign] = new TH2D(Form("h_pt2nd_vs_phi2nd_2mu4_sign%d",isign+1),";#phi;p_{T} [GeV]",nphi_bins_trig_effcy,-pms.PI,pms.PI,int(pT_bins_single_muon.size())-1,pT_bins_single_muon.data());
            h_phi2nd_vs_q_eta_2nd_2mu4[isign] =new TH2D(Form("h_phi2nd_vs_q_eta_2nd_2mu4_sign%d",isign+1),";q*#eta;#phi",neta_bins_trig_effcy,eta_bins_trig_effcy.data(),nphi_bins_trig_effcy,-pms.PI,pms.PI);

            // separated (dR > 0.8)
            h_pt2nd_vs_q_eta_2nd_mu4_sepr[isign] = new TH2D(Form("h_pt2nd_vs_q_eta_2nd_mu4_sepr_sign%d",isign+1),";q*#eta;p_{T} [GeV]",neta_bins_trig_effcy,eta_bins_trig_effcy.data(), int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
            h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sepr[isign] = new TH2D(Form("h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sepr_sign%d",isign+1),";q*#eta;p_{T} [GeV]",neta_bins_trig_effcy,eta_bins_trig_effcy.data(), int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
            h_pt2nd_vs_q_eta_2nd_2mu4_sepr[isign] = new TH2D(Form("h_pt2nd_vs_q_eta_2nd_2mu4_sepr_sign%d",isign+1),";q*#eta;p_{T} [GeV]",neta_bins_trig_effcy,eta_bins_trig_effcy.data(), int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());

            h_pt2nd_vs_phi2nd_mu4_sepr[isign] = new TH2D(Form("h_pt2nd_vs_phi2nd_mu4_sepr_sign%d",isign+1),";#phi;p_{T} [GeV]",nphi_bins_trig_effcy,-pms.PI,pms.PI,int(pT_bins_single_muon.size())-1,pT_bins_single_muon.data());
            h_pt2nd_vs_phi2nd_mu4_mu4noL1_sepr[isign] = new TH2D(Form("h_pt2nd_vs_phi2nd_mu4_mu4noL1_sepr_sign%d",isign+1),";#phi;p_{T} [GeV]",nphi_bins_trig_effcy,-pms.PI,pms.PI,int(pT_bins_single_muon.size())-1,pT_bins_single_muon.data());
            h_pt2nd_vs_phi2nd_2mu4_sepr[isign] = new TH2D(Form("h_pt2nd_vs_phi2nd_2mu4_sepr_sign%d",isign+1),";#phi;p_{T} [GeV]",nphi_bins_trig_effcy,-pms.PI,pms.PI,int(pT_bins_single_muon.size())-1,pT_bins_single_muon.data());
            
            // separated divided
            h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sepr_divided[isign] = new TH2D(Form("h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sepr_divided_sign%d",isign+1),";q*#eta;p_{T} [GeV]",neta_bins_trig_effcy,eta_bins_trig_effcy.data(), int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
            h_pt2nd_vs_q_eta_2nd_2mu4_sepr_divided[isign] = new TH2D(Form("h_pt2nd_vs_q_eta_2nd_2mu4_sepr_divided_sign%d",isign+1),";q*#eta;p_{T} [GeV]",neta_bins_trig_effcy,eta_bins_trig_effcy.data(), int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());

            // signal
            h_pt2nd_vs_q_eta_2nd_mu4_w_single_b_sig_sel[isign] = new TH2D(Form("h_pt2nd_vs_q_eta_2nd_mu4_w_single_b_sig_sel_sign%d",isign+1),";q*#eta;p_{T} [GeV]",neta_bins_trig_effcy,eta_bins_trig_effcy.data(), int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
            h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_w_single_b_sig_sel[isign] = new TH2D(Form("h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_w_single_b_sig_sel_sign%d",isign+1),";q*#eta;p_{T} [GeV]",neta_bins_trig_effcy,eta_bins_trig_effcy.data(), int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
            h_pt2nd_vs_q_eta_2nd_2mu4_w_single_b_sig_sel[isign] = new TH2D(Form("h_pt2nd_vs_q_eta_2nd_2mu4_w_single_b_sig_sel_sign%d",isign+1),";q*#eta;p_{T} [GeV]",neta_bins_trig_effcy,eta_bins_trig_effcy.data(), int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
            
            h_pt2nd_vs_phi2nd_mu4_w_single_b_sig_sel[isign] = new TH2D(Form("h_pt2nd_vs_phi2nd_mu4_w_single_b_sig_sel_sign%d",isign+1),";#phi;p_{T} [GeV]",nphi_bins_trig_effcy,-pms.PI,pms.PI,int(pT_bins_single_muon.size())-1,pT_bins_single_muon.data());
            h_pt2nd_vs_phi2nd_mu4_mu4noL1_w_single_b_sig_sel[isign] = new TH2D(Form("h_pt2nd_vs_phi2nd_mu4_mu4noL1_w_single_b_sig_sel_sign%d",isign+1),";#phi;p_{T} [GeV]",nphi_bins_trig_effcy,-pms.PI,pms.PI,int(pT_bins_single_muon.size())-1,pT_bins_single_muon.data());
            h_pt2nd_vs_phi2nd_2mu4_w_single_b_sig_sel[isign] = new TH2D(Form("h_pt2nd_vs_phi2nd_2mu4_w_single_b_sig_sel_sign%d",isign+1),";#phi;p_{T} [GeV]",nphi_bins_trig_effcy,-pms.PI,pms.PI,int(pT_bins_single_muon.size())-1,pT_bins_single_muon.data());

            // signal divided
            h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_divided_w_single_b_sig_sel[isign] = new TH2D(Form("h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_divided_w_single_b_sig_sel_sign%d",isign+1),";q*#eta;p_{T} [GeV]",neta_bins_trig_effcy,eta_bins_trig_effcy.data(), int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
            h_pt2nd_vs_q_eta_2nd_2mu4_divided_w_single_b_sig_sel[isign] = new TH2D(Form("h_pt2nd_vs_q_eta_2nd_2mu4_divided_w_single_b_sig_sel_sign%d",isign+1),";q*#eta;p_{T} [GeV]",neta_bins_trig_effcy,eta_bins_trig_effcy.data(), int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());


            // --------------------- 3D 2nd-muon histograms for trigger efficiency ---------------------------

            // mu4 sepr
            h_pt2nd_vs_q_eta_2nd_vs_phi2nd_mu4_sepr[isign] =
                new TH3D(Form("h_pt2nd_vs_q_eta_2nd_vs_phi2nd_mu4_sepr_sign%d",isign+1),
                         ";#phi_{2nd};q#eta_{2nd};p_{T,2nd} [GeV]",
                         nphi_bins_trig_effcy, phi2nd_bins.data(),
                         neta_bins_trig_effcy, eta_bins_trig_effcy.data(),
                         int(pT_bins_single_muon.size())-1, pT_bins_single_muon.data());

            // mu4_mu4noL1 sepr
            h_pt2nd_vs_q_eta_2nd_vs_phi2nd_mu4_mu4noL1_sepr[isign] =
                new TH3D(Form("h_pt2nd_vs_q_eta_2nd_vs_phi2nd_mu4_mu4noL1_sepr_sign%d",isign+1),
                         ";#phi_{2nd};q#eta_{2nd};p_{T,2nd} [GeV]",
                         nphi_bins_trig_effcy, phi2nd_bins.data(),
                         neta_bins_trig_effcy, eta_bins_trig_effcy.data(),
                         int(pT_bins_single_muon.size())-1, pT_bins_single_muon.data());

            // 2mu4 sepr
            h_pt2nd_vs_q_eta_2nd_vs_phi2nd_2mu4_sepr[isign] =
                new TH3D(Form("h_pt2nd_vs_q_eta_2nd_vs_phi2nd_2mu4_sepr_sign%d",isign+1),
                         ";#phi_{2nd};q#eta_{2nd};p_{T,2nd} [GeV]",
                         nphi_bins_trig_effcy, phi2nd_bins.data(),
                         neta_bins_trig_effcy, eta_bins_trig_effcy.data(),
                         int(pT_bins_single_muon.size())-1, pT_bins_single_muon.data());

            // mu4_mu4noL1 sepr divided
            h_pt2nd_vs_q_eta_2nd_vs_phi2nd_mu4_mu4noL1_sepr_divided[isign] =
                new TH3D(Form("h_pt2nd_vs_q_eta_2nd_vs_phi2nd_mu4_mu4noL1_sepr_divided_sign%d",isign+1),
                         ";#phi_{2nd};q#eta_{2nd};p_{T,2nd} [GeV]",
                         nphi_bins_trig_effcy, phi2nd_bins.data(),
                         neta_bins_trig_effcy, eta_bins_trig_effcy.data(),
                         int(pT_bins_single_muon.size())-1, pT_bins_single_muon.data());

            // 2mu4 sepr divided
            h_pt2nd_vs_q_eta_2nd_vs_phi2nd_2mu4_sepr_divided[isign] =
                new TH3D(Form("h_pt2nd_vs_q_eta_2nd_vs_phi2nd_2mu4_sepr_divided_sign%d",isign+1),
                         ";#phi_{2nd};q#eta_{2nd};p_{T,2nd} [GeV]",
                         nphi_bins_trig_effcy, phi2nd_bins.data(),
                         neta_bins_trig_effcy, eta_bins_trig_effcy.data(),
                         int(pT_bins_single_muon.size())-1, pT_bins_single_muon.data());


            // mu4 signal
            h_pt2nd_vs_q_eta_2nd_vs_phi2nd_mu4_w_single_b_sig_sel[isign] =
                new TH3D(Form("h_pt2nd_vs_q_eta_2nd_vs_phi2nd_mu4_w_single_b_sig_sel_sign%d",isign+1),
                         ";#phi_{2nd};q#eta_{2nd};p_{T,2nd} [GeV]",
                         nphi_bins_trig_effcy, phi2nd_bins.data(),
                         neta_bins_trig_effcy, eta_bins_trig_effcy.data(),
                         int(pT_bins_single_muon.size())-1, pT_bins_single_muon.data());

            // mu4_mu4noL1 signal
            h_pt2nd_vs_q_eta_2nd_vs_phi2nd_mu4_mu4noL1_w_single_b_sig_sel[isign] =
                new TH3D(Form("h_pt2nd_vs_q_eta_2nd_vs_phi2nd_mu4_mu4noL1_w_single_b_sig_sel_sign%d",isign+1),
                         ";#phi_{2nd};q#eta_{2nd};p_{T,2nd} [GeV]",
                         nphi_bins_trig_effcy, phi2nd_bins.data(),
                         neta_bins_trig_effcy, eta_bins_trig_effcy.data(),
                         int(pT_bins_single_muon.size())-1, pT_bins_single_muon.data());

            // 2mu4 signal
            h_pt2nd_vs_q_eta_2nd_vs_phi2nd_2mu4_w_single_b_sig_sel[isign] =
                new TH3D(Form("h_pt2nd_vs_q_eta_2nd_vs_phi2nd_2mu4_w_single_b_sig_sel_sign%d",isign+1),
                         ";#phi_{2nd};q#eta_{2nd};p_{T,2nd} [GeV]",
                         nphi_bins_trig_effcy, phi2nd_bins.data(),
                         neta_bins_trig_effcy, eta_bins_trig_effcy.data(),
                         int(pT_bins_single_muon.size())-1, pT_bins_single_muon.data());

            // mu4_mu4noL1 signal divided
            h_pt2nd_vs_q_eta_2nd_vs_phi2nd_mu4_mu4noL1_divided_w_single_b_sig_sel[isign] =
                new TH3D(Form("h_pt2nd_vs_q_eta_2nd_vs_phi2nd_mu4_mu4noL1_divided_w_single_b_sig_sel_sign%d",isign+1),
                         ";#phi_{2nd};q#eta_{2nd};p_{T,2nd} [GeV]",
                         nphi_bins_trig_effcy, phi2nd_bins.data(),
                         neta_bins_trig_effcy, eta_bins_trig_effcy.data(),
                         int(pT_bins_single_muon.size())-1, pT_bins_single_muon.data());

            // 2mu4 signal divided
            h_pt2nd_vs_q_eta_2nd_vs_phi2nd_2mu4_divided_w_single_b_sig_sel[isign] =
                new TH3D(Form("h_pt2nd_vs_q_eta_2nd_vs_phi2nd_2mu4_divided_w_single_b_sig_sel_sign%d",isign+1),
                         ";#phi_{2nd};q#eta_{2nd};p_{T,2nd} [GeV]",
                         nphi_bins_trig_effcy, phi2nd_bins.data(),
                         neta_bins_trig_effcy, eta_bins_trig_effcy.data(),
                         int(pT_bins_single_muon.size())-1, pT_bins_single_muon.data());
            
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

            h_pair_eta_vs_pair_pT_mu4_sepr[isign] = new TH2D(Form("h_pair_eta_vs_pair_pT_mu4_sepr_sign%d",isign+1),";p_{T}^{pair} [GeV];#eta^{pair}",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),neta_bins_trig_effcy,eta_bins_trig_effcy.data());
            h_Deta_Dphi_mu4_sepr[isign] = new TH2D(Form("h_Deta_Dphi_mu4_sepr_sign%d",isign+1),";#Delta#phi;#Delta#eta", nDphi_bins_trig_effcy,-pms.PI,pms.PI,nDeta_bins_trig_effcy,-4.8,4.8);
            h_eta1_eta2_mu4_sepr[isign] = new TH2D(Form("h_eta1_eta2_mu4_sepr_sign%d",isign+1),";#eta_{sublead};#eta_{lead}",neta_bins_trig_effcy,eta_bins_trig_effcy.data(), neta_bins_trig_effcy,eta_bins_trig_effcy.data());
            h_eta_avg_Deta_mu4_sepr[isign] = new TH2D(Form("h_eta_avg_Deta_mu4_sepr_sign%d",isign+1),";#bar{#eta};#Delta#eta", neta_bins_trig_effcy,eta_bins_trig_effcy.data(),nDeta_bins_trig_effcy,-4.8,4.8);
            h_eta_avg_Dphi_mu4_sepr[isign] = new TH2D(Form("h_eta_avg_Dphi_mu4_sepr_sign%d",isign+1),";#Delta#phi;#bar{#eta}", nDphi_bins_trig_effcy,-pms.PI,pms.PI,neta_bins_trig_effcy,eta_bins_trig_effcy.data());
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

            h_pair_eta_vs_pair_pT_mu4_mu4noL1_sepr[isign] = new TH2D(Form("h_pair_eta_vs_pair_pT_mu4_mu4noL1_sepr_sign%d",isign+1),";p_{T}^{pair} [GeV];#eta^{pair}",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),neta_bins_trig_effcy,eta_bins_trig_effcy.data());
            h_Deta_Dphi_mu4_mu4noL1_sepr[isign] = new TH2D(Form("h_Deta_Dphi_mu4_mu4noL1_sepr_sign%d",isign+1),";#Delta#phi;#Delta#eta", nDphi_bins_trig_effcy,-pms.PI,pms.PI,nDeta_bins_trig_effcy,-4.8,4.8);
            h_eta1_eta2_mu4_mu4noL1_sepr[isign] = new TH2D(Form("h_eta1_eta2_mu4_mu4noL1_sepr_sign%d",isign+1),";#eta_{sublead};#eta_{lead}",neta_bins_trig_effcy,eta_bins_trig_effcy.data(), neta_bins_trig_effcy,eta_bins_trig_effcy.data());
            h_eta_avg_Deta_mu4_mu4noL1_sepr[isign] = new TH2D(Form("h_eta_avg_Deta_mu4_mu4noL1_sepr_sign%d",isign+1),";#bar{#eta};#Delta#eta", neta_bins_trig_effcy,eta_bins_trig_effcy.data(),nDeta_bins_trig_effcy,-4.8,4.8);
            h_eta_avg_Dphi_mu4_mu4noL1_sepr[isign] = new TH2D(Form("h_eta_avg_Dphi_mu4_mu4noL1_sepr_sign%d",isign+1),";#Delta#phi;#bar{#eta}", nDphi_bins_trig_effcy,-pms.PI,pms.PI,neta_bins_trig_effcy,eta_bins_trig_effcy.data());
            h_minv_pair_pt_log_mu4_mu4noL1_sepr[isign] = new TH2D(Form("h_minv_pair_pt_log_mu4_mu4noL1_sepr_sign%d",isign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),nminv_bins_log,minv_bins_log[isign]);


            

            // pair kinematics requiring 2mu4 trigger
            h_Deta_2mu4_sepr[isign] = new TH1D(Form("h_Deta_2mu4_sepr_sign%d",isign+1),";#Delta#eta;",nDeta_bins_trig_effcy,-4.8,4.8);
            h_Deta_zoomin_2mu4_sepr[isign] = new TH1D(Form("h_Deta_zoomin_2mu4_sepr_sign%d",isign+1),";#Delta#eta;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
            h_Dphi_2mu4_sepr[isign] = new TH1D(Form("h_Dphi_2mu4_sepr_sign%d",isign+1),";#Delta#phi;",nDphi_bins_trig_effcy,-pms.PI,pms.PI);
            h_Dphi_zoomin_2mu4_sepr[isign] = new TH1D(Form("h_Dphi_zoomin_2mu4_sepr_sign%d",isign+1),";#Delta#phi;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
            h_DR_2mu4_sepr[isign] = new TH1D(Form("h_DR_2mu4_sepr_sign%d",isign+1),";#DeltaR;",nDR_bins_trig_effcy,0,pms.deltaR_thrsh[2]);
            h_DR_zoomin_2mu4_sepr[isign] = new TH1D(Form("h_DR_zoomin_2mu4_sepr_sign%d",isign+1),";#DeltaR;",nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
            h_minv_zoomin_2mu4_sepr[isign] = new TH1D(Form("h_minv_zoomin_2mu4_sepr_sign%d",isign+1),";m_{#mu#mu} [GeV];",nminv_bins_trig_effcy,0,3.);
            h_pair_pt_log_2mu4_sepr[isign] = new TH1D(Form("h_pair_pt_log_2mu4_sepr_sign%d",isign+1),";p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());

            h_pair_eta_vs_pair_pT_2mu4_sepr[isign] = new TH2D(Form("h_pair_eta_vs_pair_pT_2mu4_sepr_sign%d",isign+1),";p_{T}^{pair} [GeV];#eta^{pair}",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),neta_bins_trig_effcy,eta_bins_trig_effcy.data());
            h_Deta_Dphi_2mu4_sepr[isign] = new TH2D(Form("h_Deta_Dphi_2mu4_sepr_sign%d",isign+1),";#Delta#phi;#Delta#eta", nDphi_bins_trig_effcy,-pms.PI,pms.PI,nDeta_bins_trig_effcy,-4.8,4.8);
            h_eta1_eta2_2mu4_sepr[isign] = new TH2D(Form("h_eta1_eta2_2mu4_sepr_sign%d",isign+1),";#eta_{sublead};#eta_{lead}",neta_bins_trig_effcy,eta_bins_trig_effcy.data(), neta_bins_trig_effcy,eta_bins_trig_effcy.data());
            h_eta_avg_Deta_2mu4_sepr[isign] = new TH2D(Form("h_eta_avg_Deta_2mu4_sepr_sign%d",isign+1),";#bar{#eta};#Delta#eta", neta_bins_trig_effcy,eta_bins_trig_effcy.data(),nDeta_bins_trig_effcy,-4.8,4.8);
            h_eta_avg_Dphi_2mu4_sepr[isign] = new TH2D(Form("h_eta_avg_Dphi_2mu4_sepr_sign%d",isign+1),";#Delta#phi;#bar{#eta}", nDphi_bins_trig_effcy,-pms.PI,pms.PI,neta_bins_trig_effcy,eta_bins_trig_effcy.data());
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
    } // end if statement - trigger mode = 1 (MUON-PAIR trigger | mu4)
    else if (trigger_mode == 0){
        for (unsigned int isign = 0; isign < ParamsSet::nSigns; isign++){
            // MB
            h_pt2nd_vs_q_eta_2nd_MB[isign] = new TH2D(Form("h_pt2nd_vs_q_eta_2nd_MB_sign%d",isign+1),";q*#eta;p_{T} [GeV]",neta_bins_trig_effcy,eta_bins_trig_effcy.data(), int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
            h_pt2nd_vs_phi2nd_MB[isign] = new TH2D(Form("h_pt2nd_vs_phi2nd_MB_sign%d",isign+1),";#phi;p_{T} [GeV]",nphi_bins_trig_effcy,-pms.PI,pms.PI,int(pT_bins_single_muon.size())-1,pT_bins_single_muon.data());
            h_phi2nd_vs_q_eta_2nd_MB[isign] =new TH2D(Form("h_phi2nd_vs_q_eta_2nd_MB_sign%d",isign+1),";q*#eta;#phi",neta_bins_trig_effcy,eta_bins_trig_effcy.data(),nphi_bins_trig_effcy,-pms.PI,pms.PI);

            // mu4
            h_pt2nd_vs_q_eta_2nd_mu4[isign] = new TH2D(Form("h_pt2nd_vs_q_eta_2nd_mu4_sign%d",isign+1),";q*#eta;p_{T} [GeV]",neta_bins_trig_effcy,eta_bins_trig_effcy.data(), int(pT_bins_single_muon.size() - 1), pT_bins_single_muon.data());
            h_pt2nd_vs_phi2nd_mu4[isign] = new TH2D(Form("h_pt2nd_vs_phi2nd_mu4_sign%d",isign+1),";#phi;p_{T} [GeV]",nphi_bins_trig_effcy,-pms.PI,pms.PI,int(pT_bins_single_muon.size())-1,pT_bins_single_muon.data());
            h_phi2nd_vs_q_eta_2nd_mu4[isign] =new TH2D(Form("h_phi2nd_vs_q_eta_2nd_mu4_sign%d",isign+1),";q*#eta;#phi",neta_bins_trig_effcy,eta_bins_trig_effcy.data(),nphi_bins_trig_effcy,-pms.PI,pms.PI);
        }
    } // end if statement - trigger mode = 0 (mu4 | MB)
}

void MuonPairPlottingPP::OpenEffcyPtFitFile(){
    // Open the ROOT file
    file_effcy_pTfits = TFile::Open(file_name_effcy_pTfits.c_str(), "READ");
    if (!file_effcy_pTfits || file_effcy_pTfits->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }
}

void MuonPairPlottingPP::ProcessData(){
	
    cout << "1st loop over muon pairs" << endl;

    for (int isign = 0; isign < ParamsSet::nSigns; isign++){ // 1st loop: fill all nominal histograms
        
        Long64_t nentries = inTree[isign]->GetEntries(); //#muon pairs
        for (Long64_t lentry=0; lentry<nentries; lentry++) {//loop over the events
        // for (Long64_t lentry=0; lentry<1000; lentry++) {//loop over the events
        if(lentry%100000==0) cout<<"Processing "<<lentry<<" event out of "<<nentries<<" events"<<std::endl;
      		int num_bytes = inTree[isign]->GetEntry(lentry);//read in an event
            if(num_bytes==0){
            	std::cout<<"Error:: Read in event has size of zero bytes,  quitting"<<std::endl;
            	throw std::exception();
            }
            
            FillHistograms(isign);
        }
    }

    if (trigger_mode != 1) return; // only continue if calculate muon-pair trigger efficiency selecting on mu4 trigger

    MakeAndWriteSingleMuonPtTrigEffGraphs();

    CalculateSingleMuonTrigEffcyRatios();

    cout << "2nd loop over muon pairs" << endl;

    for (int isign = 0; isign < ParamsSet::nSigns; isign++){ // 2nd loop: fill pair efficiency histograms inversely weighted by single-muon trigger efficiencies
        
        Long64_t nentries = inTree[isign]->GetEntries(); //#muon pairs
        for (Long64_t lentry=0; lentry<nentries; lentry++) {//loop over the events
        // for (Long64_t lentry=0; lentry<1000; lentry++) {//loop over the events
        if(lentry%100000==0) cout<<"Processing "<<lentry<<" event out of "<<nentries<<" events"<<std::endl;
            int num_bytes = inTree[isign]->GetEntry(lentry);//read in an event
            if(num_bytes==0){
                std::cout<<"Error:: Read in event has size of zero bytes,  quitting"<<std::endl;
                throw std::exception();
            }
            
            FillTrigEffcyHistsInvWeightedbySingleMuonEffcies(isign);
        }
    }

    MakeAndWriteDRTrigEffGraphs();

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

// helper function to find a bin (range) containing a number (q*eta value) and return the string given by the pairToSuffix function
std::string MuonPairPlottingPP::FindBinReturnStr(float number, const std::vector<std::pair<float, float>>& ranges) {
    for (const auto& r : ranges) {
        if (number >= r.first && number <= r.second) {
            return pairToSuffix(r); // found the enclosing range
        }
    }
    // q*eta value is not in any bin used for pT fitting --> return empty string & use binned value (not fitted) for trigger efficiencies
    return "";
}

void MuonPairPlottingPP::CalculateSingleMuonTrigEffcyRatios(){
    for (int isign = 0; isign < 2; isign++){ // looping over the 2nd-muon charge sign
        // mu4_mu4noL1 2D
        h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sepr_divided[isign] = (TH2D*)h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sepr[isign]->Clone(Form("%s_clone", h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sepr[isign]->GetName()));
        h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sepr_divided[isign]->Divide(h_pt2nd_vs_q_eta_2nd_mu4_sepr[isign]);
        
        // 2mu4 2D
        h_pt2nd_vs_q_eta_2nd_2mu4_sepr_divided[isign] = (TH2D*)h_pt2nd_vs_q_eta_2nd_2mu4_sepr[isign]->Clone(Form("%s_clone", h_pt2nd_vs_q_eta_2nd_2mu4_sepr[isign]->GetName()));
        h_pt2nd_vs_q_eta_2nd_2mu4_sepr_divided[isign]->Divide(h_pt2nd_vs_q_eta_2nd_mu4_sepr[isign]);

        // mu4_mu4noL1 signal 2D
        h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_divided_w_single_b_sig_sel[isign] = (TH2D*)h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_w_single_b_sig_sel[isign]->Clone(Form("%s_clone", h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_w_single_b_sig_sel[isign]->GetName()));
        h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_divided_w_single_b_sig_sel[isign]->Divide(h_pt2nd_vs_q_eta_2nd_mu4_w_single_b_sig_sel[isign]);
        
        // 2mu4 signal 2D
        h_pt2nd_vs_q_eta_2nd_2mu4_divided_w_single_b_sig_sel[isign] = (TH2D*)h_pt2nd_vs_q_eta_2nd_2mu4_w_single_b_sig_sel[isign]->Clone(Form("%s_clone", h_pt2nd_vs_q_eta_2nd_2mu4_w_single_b_sig_sel[isign]->GetName()));
        h_pt2nd_vs_q_eta_2nd_2mu4_divided_w_single_b_sig_sel[isign]->Divide(h_pt2nd_vs_q_eta_2nd_mu4_w_single_b_sig_sel[isign]);

        // mu4_mu4noL1 3D
        h_pt2nd_vs_q_eta_2nd_vs_phi2nd_mu4_mu4noL1_sepr_divided[isign] = (TH3D*) h_pt2nd_vs_q_eta_2nd_vs_phi2nd_mu4_mu4noL1_sepr[isign]->Clone(Form("%s_clone", h_pt2nd_vs_q_eta_2nd_vs_phi2nd_mu4_mu4noL1_sepr[isign]->GetName()));
        h_pt2nd_vs_q_eta_2nd_vs_phi2nd_mu4_mu4noL1_sepr_divided[isign]->Divide(h_pt2nd_vs_q_eta_2nd_vs_phi2nd_mu4_sepr[isign]);

        // 2mu4 3D
        h_pt2nd_vs_q_eta_2nd_vs_phi2nd_2mu4_sepr_divided[isign] = (TH3D*) h_pt2nd_vs_q_eta_2nd_vs_phi2nd_2mu4_sepr[isign]->Clone(Form("%s_clone", h_pt2nd_vs_q_eta_2nd_vs_phi2nd_2mu4_sepr[isign]->GetName()));
        h_pt2nd_vs_q_eta_2nd_vs_phi2nd_2mu4_sepr_divided[isign]->Divide(h_pt2nd_vs_q_eta_2nd_vs_phi2nd_mu4_sepr[isign]);

        // mu4_mu4noL1 signal 3D
        h_pt2nd_vs_q_eta_2nd_vs_phi2nd_mu4_mu4noL1_divided_w_single_b_sig_sel[isign] = (TH3D*) h_pt2nd_vs_q_eta_2nd_vs_phi2nd_mu4_mu4noL1_w_single_b_sig_sel[isign]->Clone(Form("%s_clone", h_pt2nd_vs_q_eta_2nd_vs_phi2nd_mu4_mu4noL1_w_single_b_sig_sel[isign]->GetName()));
        h_pt2nd_vs_q_eta_2nd_vs_phi2nd_mu4_mu4noL1_divided_w_single_b_sig_sel[isign]->Divide(h_pt2nd_vs_q_eta_2nd_vs_phi2nd_mu4_w_single_b_sig_sel[isign]);

        // 2mu4 signal 3D
        h_pt2nd_vs_q_eta_2nd_vs_phi2nd_2mu4_divided_w_single_b_sig_sel[isign] = (TH3D*) h_pt2nd_vs_q_eta_2nd_vs_phi2nd_2mu4_w_single_b_sig_sel[isign]->Clone(Form("%s_clone", h_pt2nd_vs_q_eta_2nd_vs_phi2nd_2mu4_w_single_b_sig_sel[isign]->GetName()));
        h_pt2nd_vs_q_eta_2nd_vs_phi2nd_2mu4_divided_w_single_b_sig_sel[isign]->Divide(h_pt2nd_vs_q_eta_2nd_vs_phi2nd_mu4_w_single_b_sig_sel[isign]);
    }
}

// helper function to evaluate single muon efficiency given the trigger of interest, the 2nd muon sign, pT, q*eta & phi 
// return -1 if using 3D (currently no fitting done) or if q*eta does not lie in range (crack area with strong q * eta dependence)
float MuonPairPlottingPP::EvaluateSingleMuonEffcyPtFitted(bool charge_sign, std::string trg, float pt_2nd, float q_eta_2nd, float phi_2nd){
    if (use_3D_2nd_muon) return -1;
        
    auto q_eta_bin_str = FindBinReturnStr(q_eta_2nd, q_eta_proj_ranges_for_single_muon_effcy_pT_fitting);

    std::string musign = charge_sign? "sign2" : "sign1";

    if (q_eta_bin_str != "") { // q*eta value found in a bin (range) --> evaluate pT-fitted functions at pT_2nd
        std::string fname = "f_pt2nd_vs_q_eta_2nd_";
        
        fname += trg + "_sepr_" + musign + "_py_" + q_eta_bin_str + "_divided";

        TF1 *func_pTfit = dynamic_cast<TF1*>(file_effcy_pTfits->Get(fname.c_str()));  // replace with actual name
        if (!func_pTfit) {
            std::cerr << "Warning: TF1 " << fname << " not found in file." << std::endl;
            return -1;
        }

        return func_pTfit->Eval(pt_2nd);
    } else return -1;   
}

float MuonPairPlottingPP::EvaluateSingleMuonEffcy(bool charge_sign, std::string trg, float pt_2nd, float q_eta_2nd, float phi_2nd){
    if (trg != "mu4_mu4noL1" && trg != "2mu4"){
        std::cerr << "EvaluateSingleMuonEffcy: The string trg passed as argument MUST equal mu4_mu4noL1 or 2mu4!" << std::endl;
        return -1;
    }

    float singleMuonEffcyPtFitted = EvaluateSingleMuonEffcyPtFitted(charge_sign, trg, pt_2nd, q_eta_2nd, phi_2nd);
    if (singleMuonEffcyPtFitted != -1) return singleMuonEffcyPtFitted;

    float effcy = -1;
    if (use_3D_2nd_muon){ // 3D binned value
        TH3D* h_effcy = (trg == "mu4_mu4noL1")? h_pt2nd_vs_q_eta_2nd_vs_phi2nd_mu4_mu4noL1_sepr_divided[charge_sign] : h_pt2nd_vs_q_eta_2nd_vs_phi2nd_2mu4_sepr_divided[charge_sign];
        int bin_num = h_effcy->FindBin(phi_2nd, q_eta_2nd, pt_2nd);
        effcy = h_effcy->GetBinContent(bin_num);
    } else{ // 2D binned value
        TH2D* h_effcy = (trg == "mu4_mu4noL1")? h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sepr_divided[charge_sign] : h_pt2nd_vs_q_eta_2nd_2mu4_sepr_divided[charge_sign];
        int bin_num = h_effcy->FindBin(q_eta_2nd, pt_2nd);
        effcy = h_effcy->GetBinContent(bin_num);
    }

    if (effcy < 0){
        std::cout << "EvaluateSingleMuonEffcy: WARNING:: Efficiency for charge sign " << charge_sign << ", trigger " << trg 
                  << ", pT " << pt_2nd << ", q*eta " << q_eta_2nd << ", phi " << phi_2nd << ", use_3D_2nd_muon? " << use_3D_2nd_muon
                  << "is negative!" << std::endl;
    }

    // if (singleMuonEffcyPtFitted != -1){
    //     cout << "pT: " << pt_2nd << ", pT fitted: " << singleMuonEffcyPtFitted << ", binned: " << effcy << endl;
    //     return singleMuonEffcyPtFitted;
    // } 
    return effcy;
}

void MuonPairPlottingPP::FillHistograms(int nsign){

    bool pass_single_mu4 [2] = {mu1PassSingle[nsign], mu2PassSingle[nsign]};
    double pt_2nd [2] = {m2pt[nsign], m1pt[nsign]};
    double pt_1st [2] = {m1pt[nsign], m2pt[nsign]};
    double q_eta_2nd [2] = {m2charge[nsign] * m2eta[nsign], m1charge[nsign] * m1eta[nsign]};
    int charge_2nd [2] = {m2charge[nsign], m1charge[nsign]};
    double eta_2nd [2] = {m2eta[nsign], m1eta[nsign]};
    double phi_2nd [2] = {m2phi[nsign], m1phi[nsign]};

    bool pass_single_b_signal_selection = (nsign == 1 && minv[nsign] > 1.08 && minv[nsign] < 2.9 && pair_pt[nsign] > 8);
    bool pass_single_muon_good_acceptance_selection[2];

    pass_single_muon_good_acceptance_selection[0] = (mu1PassSingle[nsign])? (m2pt[nsign] >= 6 && ((m2eta[nsign] > 1.1 && m2eta[nsign] < 2.3) || (m2eta[nsign] > -2.3 && m2eta[nsign] < -1.2))) : false;
    pass_single_muon_good_acceptance_selection[1] = (mu2PassSingle[nsign])? (m1pt[nsign] >= 6 && ((m1eta[nsign] > 1.1 && m1eta[nsign] < 2.3) || (m1eta[nsign] > -2.3 && m1eta[nsign] < -1.2))) : false;

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

    if (trigger_mode == 1){ // evaluate MUON-PAIR trigger efficiencies, selecting on mu4
        for (int muon_ind = 0; muon_ind < 2; muon_ind++){ // loop over the two muons + check if either passes the single-muon mu4
            if (!pass_single_mu4[muon_ind]) continue;
            bool charge_sign_2nd_muon = (charge_2nd[muon_ind] > 0)? 0: 1;

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

            // pt-2nd-dependent histograms
            h_pt2nd_vs_q_eta_2nd_mu4[charge_sign_2nd_muon]->Fill(q_eta_2nd[muon_ind], pt_2nd[muon_ind],weight[nsign]);
            h_pt2nd_vs_phi2nd_mu4[charge_sign_2nd_muon]->Fill(phi_2nd[muon_ind], pt_2nd[muon_ind], weight[nsign]);
            h_phi2nd_vs_q_eta_2nd_mu4[charge_sign_2nd_muon]->Fill(q_eta_2nd[muon_ind], phi_2nd[muon_ind], weight[nsign]);
            h_DR_zoomin_vs_pt2nd_mu4[nsign]->Fill(pt_2nd[muon_ind],dr[nsign],weight[nsign]);
            h_DR_0_2_vs_pt2nd_mu4[nsign]->Fill(pt_2nd[muon_ind],dr[nsign],weight[nsign]);

            // pt-1st-dependent histograms
            h_Deta_vs_pT_1st_mu4[nsign]->Fill(pt_1st[muon_ind], deta[nsign],weight[nsign]);
            h_Deta_zoomin_vs_pT_1st_mu4[nsign]->Fill(pt_1st[muon_ind], deta[nsign],weight[nsign]);
            h_Dphi_vs_pT_1st_mu4[nsign]->Fill(pt_1st[muon_ind], dphi[nsign],weight[nsign]);
            h_Dphi_zoomin_vs_pT_1st_mu4[nsign]->Fill(pt_1st[muon_ind], dphi[nsign],weight[nsign]);
            h_DR_vs_pT_1st_mu4[nsign]->Fill(pt_1st[muon_ind], dr[nsign],weight[nsign]);
            h_DR_zoomin_vs_pT_1st_mu4[nsign]->Fill(pt_1st[muon_ind], dr[nsign],weight[nsign]);

            // pair-observable 2D histograms
            h_pair_eta_vs_pair_pT_mu4[nsign]->Fill(pair_pt[nsign],pair_eta[nsign],weight[nsign]);
            h_Deta_Dphi_mu4[nsign]->Fill(dphi[nsign],deta[nsign],weight[nsign]);
            h_eta1_eta2_mu4[nsign]->Fill(m2eta[nsign],m1eta[nsign],weight[nsign]);
            h_eta_avg_Deta_mu4[nsign]->Fill(deta[nsign],etaavg[nsign],weight[nsign]);
            h_eta_avg_Dphi_mu4[nsign]->Fill(dphi[nsign],etaavg[nsign],weight[nsign]);
            h_minv_zoomin_mu4[nsign]->Fill(minv[nsign],weight[nsign]);
            h_pair_pt_log_mu4[nsign]->Fill(pair_pt[nsign],weight[nsign]);
            h_minv_pair_pt_log_mu4[nsign]->Fill(pair_pt[nsign],minv[nsign],weight[nsign]);
           
            // good-acceptance histograms
            if (pass_single_muon_good_acceptance_selection[muon_ind]){ // fill for both
                h_Deta_mu4_good_accept[nsign]->Fill(deta[nsign],weight[nsign]);
                h_Deta_zoomin_mu4_good_accept[nsign]->Fill(deta[nsign],weight[nsign]);
                h_Dphi_mu4_good_accept[nsign]->Fill(dphi[nsign],weight[nsign]);
                h_Dphi_zoomin_mu4_good_accept[nsign]->Fill(dphi[nsign],weight[nsign]);
                h_DR_mu4_good_accept[nsign]->Fill(dr[nsign],weight[nsign]);
                h_DR_zoomin_mu4_good_accept[nsign]->Fill(dr[nsign],weight[nsign]);
                h_minv_zoomin_mu4_good_accept[nsign]->Fill(minv[nsign],weight[nsign]);
                h_pair_pt_log_mu4_good_accept[nsign]->Fill(pair_pt[nsign],weight[nsign]);    
            }

            // well-separated histograms
            if (passSeparated[nsign]){
                h_Deta_mu4_sepr[nsign]->Fill(deta[nsign],weight[nsign]);
                h_Deta_zoomin_mu4_sepr[nsign]->Fill(deta[nsign],weight[nsign]);
                h_Deta_Dphi_mu4_sepr[nsign]->Fill(dphi[nsign],deta[nsign],weight[nsign]);
                h_eta_avg_Deta_mu4_sepr[nsign]->Fill(deta[nsign],etaavg[nsign],weight[nsign]);

                h_Dphi_mu4_sepr[nsign]->Fill(dphi[nsign],weight[nsign]);
                h_Dphi_zoomin_mu4_sepr[nsign]->Fill(dphi[nsign],weight[nsign]);
                h_DR_mu4_sepr[nsign]->Fill(dr[nsign],weight[nsign]);
                h_DR_zoomin_mu4_sepr[nsign]->Fill(dr[nsign],weight[nsign]);
                
                h_pt2nd_vs_q_eta_2nd_mu4_sepr[charge_sign_2nd_muon]->Fill(q_eta_2nd[muon_ind], pt_2nd[muon_ind],weight[nsign]);
                h_pt2nd_vs_phi2nd_mu4_sepr[charge_sign_2nd_muon]->Fill(phi_2nd[muon_ind], pt_2nd[muon_ind],weight[nsign]);
                h_pt2nd_vs_q_eta_2nd_vs_phi2nd_mu4_sepr[charge_sign_2nd_muon]->Fill(phi_2nd[muon_ind], q_eta_2nd[muon_ind], pt_2nd[muon_ind], weight[nsign]);

                h_pair_eta_vs_pair_pT_mu4_sepr[nsign]->Fill(pair_pt[nsign],pair_eta[nsign],weight[nsign]);
                h_eta1_eta2_mu4_sepr[nsign]->Fill(m2eta[nsign],m1eta[nsign],weight[nsign]);
                h_eta_avg_Dphi_mu4_sepr[nsign]->Fill(dphi[nsign],etaavg[nsign],weight[nsign]);
                h_minv_zoomin_mu4_sepr[nsign]->Fill(minv[nsign],weight[nsign]);
                h_pair_pt_log_mu4_sepr[nsign]->Fill(pair_pt[nsign],weight[nsign]);
                h_minv_pair_pt_log_mu4_sepr[nsign]->Fill(pair_pt[nsign],minv[nsign],weight[nsign]);
                    
            } // end if statement for well-separated pairs

            // signal-selected-pair histograms
            if (pass_single_b_signal_selection){
                h_DR_zoomin_vs_pair_pT_mu4_w_single_b_sig_sel->Fill(pair_pt[nsign], dr[nsign],weight[nsign]);

                h_Deta_mu4_w_single_b_sig_sel->Fill(deta[nsign],weight[nsign]);
                h_Deta_zoomin_mu4_w_single_b_sig_sel->Fill(deta[nsign],weight[nsign]);
                h_Dphi_mu4_w_single_b_sig_sel->Fill(dphi[nsign],weight[nsign]);
                h_Dphi_zoomin_mu4_w_single_b_sig_sel->Fill(dphi[nsign],weight[nsign]);
                h_DR_mu4_w_single_b_sig_sel->Fill(dr[nsign],weight[nsign]);
                h_DR_zoomin_mu4_w_single_b_sig_sel->Fill(dr[nsign],weight[nsign]);
                h_DR_0_2_mu4_w_single_b_sig_sel->Fill(dr[nsign],weight[nsign]);

                h_pt2nd_vs_q_eta_2nd_mu4_w_single_b_sig_sel[charge_sign_2nd_muon]->Fill(q_eta_2nd[muon_ind], pt_2nd[muon_ind],weight[nsign]);
                h_pt2nd_vs_phi2nd_mu4_w_single_b_sig_sel[charge_sign_2nd_muon]->Fill(phi_2nd[muon_ind], pt_2nd[muon_ind],weight[nsign]);
                h_pt2nd_vs_q_eta_2nd_vs_phi2nd_mu4_w_single_b_sig_sel[charge_sign_2nd_muon]->Fill(phi_2nd[muon_ind], q_eta_2nd[muon_ind], pt_2nd[muon_ind], weight[nsign]);
                h_DR_zoomin_vs_pt2nd_mu4_w_single_b_sig_sel->Fill(pt_2nd[muon_ind],dr[nsign],weight[nsign]);
                h_DR_0_2_vs_pt2nd_mu4_w_single_b_sig_sel->Fill(pt_2nd[muon_ind],dr[nsign],weight[nsign]);

                h_pair_eta_vs_pair_pT_mu4_w_single_b_sig_sel->Fill(pair_pt[nsign],pair_eta[nsign],weight[nsign]);
                h_Deta_Dphi_mu4_w_single_b_sig_sel->Fill(dphi[nsign],deta[nsign],weight[nsign]);
                h_eta1_eta2_mu4_w_single_b_sig_sel->Fill(m2eta[nsign],m1eta[nsign],weight[nsign]);
                h_eta_avg_Deta_mu4_w_single_b_sig_sel->Fill(deta[nsign],etaavg[nsign],weight[nsign]);
                h_eta_avg_Dphi_mu4_w_single_b_sig_sel->Fill(dphi[nsign],etaavg[nsign],weight[nsign]);
                h_minv_zoomin_mu4_w_single_b_sig_sel->Fill(minv[nsign],weight[nsign]);
                h_pair_pt_log_mu4_w_single_b_sig_sel->Fill(pair_pt[nsign],weight[nsign]);
                h_minv_pair_pt_log_mu4_w_single_b_sig_sel->Fill(pair_pt[nsign],minv[nsign],weight[nsign]);
            } // end if statement for single-b signal pairs
        } // end loop over muon index

        if (passmu4mu4noL1[nsign]){
            for (int muon_ind = 0; muon_ind < 2; muon_ind++){ // loop over the two muons + check if either passes the single-muon mu4
                if (!pass_single_mu4[muon_ind]) continue;
                bool charge_sign_2nd_muon = (charge_2nd[muon_ind] > 0)? 0: 1;

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

                // pt-2nd-dependent histograms
                h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1[charge_sign_2nd_muon]->Fill(q_eta_2nd[muon_ind], pt_2nd[muon_ind],weight[nsign]);
                h_pt2nd_vs_phi2nd_mu4_mu4noL1[charge_sign_2nd_muon]->Fill(phi_2nd[muon_ind], pt_2nd[muon_ind], weight[nsign]);
                h_phi2nd_vs_q_eta_2nd_mu4_mu4noL1[charge_sign_2nd_muon]->Fill(q_eta_2nd[muon_ind], phi_2nd[muon_ind], weight[nsign]);
                h_DR_zoomin_vs_pt2nd_mu4_mu4noL1[nsign]->Fill(pt_2nd[muon_ind],dr[nsign],weight[nsign]);
                h_DR_0_2_vs_pt2nd_mu4_mu4noL1[nsign]->Fill(pt_2nd[muon_ind],dr[nsign],weight[nsign]);

                // pt-1st-dependent histograms
                h_Deta_vs_pT_1st_mu4_mu4noL1[nsign]->Fill(pt_1st[muon_ind], deta[nsign],weight[nsign]);
                h_Deta_zoomin_vs_pT_1st_mu4_mu4noL1[nsign]->Fill(pt_1st[muon_ind], deta[nsign],weight[nsign]);
                h_Dphi_vs_pT_1st_mu4_mu4noL1[nsign]->Fill(pt_1st[muon_ind], dphi[nsign],weight[nsign]);
                h_Dphi_zoomin_vs_pT_1st_mu4_mu4noL1[nsign]->Fill(pt_1st[muon_ind], dphi[nsign],weight[nsign]);
                h_DR_vs_pT_1st_mu4_mu4noL1[nsign]->Fill(pt_1st[muon_ind], dr[nsign],weight[nsign]);
                h_DR_zoomin_vs_pT_1st_mu4_mu4noL1[nsign]->Fill(pt_1st[muon_ind], dr[nsign],weight[nsign]);
                
                // pair-observable 2D histograms
                h_pair_eta_vs_pair_pT_mu4_mu4noL1[nsign]->Fill(pair_pt[nsign],pair_eta[nsign],weight[nsign]);
                h_Deta_Dphi_mu4_mu4noL1[nsign]->Fill(dphi[nsign],deta[nsign],weight[nsign]);
                h_eta1_eta2_mu4_mu4noL1[nsign]->Fill(m2eta[nsign],m1eta[nsign],weight[nsign]);
                h_eta_avg_Deta_mu4_mu4noL1[nsign]->Fill(deta[nsign],etaavg[nsign],weight[nsign]);
                h_eta_avg_Dphi_mu4_mu4noL1[nsign]->Fill(dphi[nsign],etaavg[nsign],weight[nsign]);
                h_minv_zoomin_mu4_mu4noL1[nsign]->Fill(minv[nsign],weight[nsign]);
                h_pair_pt_log_mu4_mu4noL1[nsign]->Fill(pair_pt[nsign],weight[nsign]);
                h_minv_pair_pt_log_mu4_mu4noL1[nsign]->Fill(pair_pt[nsign],minv[nsign],weight[nsign]);

                // good-acceptance histograms
                if (pass_single_muon_good_acceptance_selection[muon_ind]){ // fill for both
                    h_Deta_mu4_mu4noL1_good_accept[nsign]->Fill(deta[nsign],weight[nsign]);
                    h_Deta_zoomin_mu4_mu4noL1_good_accept[nsign]->Fill(deta[nsign],weight[nsign]);
                    h_Dphi_mu4_mu4noL1_good_accept[nsign]->Fill(dphi[nsign],weight[nsign]);
                    h_Dphi_zoomin_mu4_mu4noL1_good_accept[nsign]->Fill(dphi[nsign],weight[nsign]);
                    h_DR_mu4_mu4noL1_good_accept[nsign]->Fill(dr[nsign],weight[nsign]);
                    h_DR_zoomin_mu4_mu4noL1_good_accept[nsign]->Fill(dr[nsign],weight[nsign]);
                    h_minv_zoomin_mu4_mu4noL1_good_accept[nsign]->Fill(minv[nsign],weight[nsign]);
                    h_pair_pt_log_mu4_mu4noL1_good_accept[nsign]->Fill(pair_pt[nsign],weight[nsign]);    
                }

                // well-separated histograms
                if (passSeparated[nsign]){
                    h_Deta_mu4_mu4noL1_sepr[nsign]->Fill(deta[nsign],weight[nsign]);
                    h_Deta_zoomin_mu4_mu4noL1_sepr[nsign]->Fill(deta[nsign],weight[nsign]);
                    h_Deta_Dphi_mu4_mu4noL1_sepr[nsign]->Fill(dphi[nsign],deta[nsign],weight[nsign]);
                    h_eta_avg_Deta_mu4_mu4noL1_sepr[nsign]->Fill(deta[nsign],etaavg[nsign],weight[nsign]);

                    h_Dphi_mu4_mu4noL1_sepr[nsign]->Fill(dphi[nsign],weight[nsign]);
                    h_Dphi_zoomin_mu4_mu4noL1_sepr[nsign]->Fill(dphi[nsign],weight[nsign]);
                    h_DR_mu4_mu4noL1_sepr[nsign]->Fill(dr[nsign],weight[nsign]);
                    h_DR_zoomin_mu4_mu4noL1_sepr[nsign]->Fill(dr[nsign],weight[nsign]);
                    
                    h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sepr[charge_sign_2nd_muon]->Fill(q_eta_2nd[muon_ind], pt_2nd[muon_ind],weight[nsign]);
                    h_pt2nd_vs_phi2nd_mu4_mu4noL1_sepr[charge_sign_2nd_muon]->Fill(phi_2nd[muon_ind], pt_2nd[muon_ind],weight[nsign]);
                    h_pt2nd_vs_q_eta_2nd_vs_phi2nd_mu4_mu4noL1_sepr[charge_sign_2nd_muon]->Fill(phi_2nd[muon_ind], q_eta_2nd[muon_ind], pt_2nd[muon_ind], weight[nsign]);

                    h_pair_eta_vs_pair_pT_mu4_mu4noL1_sepr[nsign]->Fill(pair_pt[nsign],pair_eta[nsign],weight[nsign]);
                    h_eta1_eta2_mu4_mu4noL1_sepr[nsign]->Fill(m2eta[nsign],m1eta[nsign],weight[nsign]);
                    h_eta_avg_Dphi_mu4_mu4noL1_sepr[nsign]->Fill(dphi[nsign],etaavg[nsign],weight[nsign]);
                    h_minv_zoomin_mu4_mu4noL1_sepr[nsign]->Fill(minv[nsign],weight[nsign]);
                    h_pair_pt_log_mu4_mu4noL1_sepr[nsign]->Fill(pair_pt[nsign],weight[nsign]);
                    h_minv_pair_pt_log_mu4_mu4noL1_sepr[nsign]->Fill(pair_pt[nsign],minv[nsign],weight[nsign]);

                } // end if statement for well-separated pairs

                // signal-selected-pair histograms
                if (pass_single_b_signal_selection){
                    h_Deta_mu4_mu4noL1_w_single_b_sig_sel->Fill(deta[nsign],weight[nsign]);
                    h_Deta_zoomin_mu4_mu4noL1_w_single_b_sig_sel->Fill(deta[nsign],weight[nsign]);
                    h_Dphi_mu4_mu4noL1_w_single_b_sig_sel->Fill(dphi[nsign],weight[nsign]);
                    h_Dphi_zoomin_mu4_mu4noL1_w_single_b_sig_sel->Fill(dphi[nsign],weight[nsign]);
                    h_DR_mu4_mu4noL1_w_single_b_sig_sel->Fill(dr[nsign],weight[nsign]);
                    h_DR_zoomin_mu4_mu4noL1_w_single_b_sig_sel->Fill(dr[nsign],weight[nsign]);
                    h_DR_0_2_mu4_mu4noL1_w_single_b_sig_sel->Fill(dr[nsign],weight[nsign]);

                    h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_w_single_b_sig_sel[charge_sign_2nd_muon]->Fill(q_eta_2nd[muon_ind], pt_2nd[muon_ind],weight[nsign]);
                    h_pt2nd_vs_phi2nd_mu4_mu4noL1_w_single_b_sig_sel[charge_sign_2nd_muon]->Fill(phi_2nd[muon_ind], pt_2nd[muon_ind],weight[nsign]);
                    h_pt2nd_vs_q_eta_2nd_vs_phi2nd_mu4_mu4noL1_w_single_b_sig_sel[charge_sign_2nd_muon]->Fill(phi_2nd[muon_ind], q_eta_2nd[muon_ind], pt_2nd[muon_ind], weight[nsign]);
                    h_DR_zoomin_vs_pt2nd_mu4_mu4noL1_w_single_b_sig_sel->Fill(pt_2nd[muon_ind],dr[nsign],weight[nsign]);
                    h_DR_0_2_vs_pt2nd_mu4_mu4noL1_w_single_b_sig_sel->Fill(pt_2nd[muon_ind],dr[nsign],weight[nsign]);

                    h_pair_eta_vs_pair_pT_mu4_mu4noL1_w_single_b_sig_sel->Fill(pair_pt[nsign],pair_eta[nsign],weight[nsign]);
                    h_Deta_Dphi_mu4_mu4noL1_w_single_b_sig_sel->Fill(dphi[nsign],deta[nsign],weight[nsign]);
                    h_eta1_eta2_mu4_mu4noL1_w_single_b_sig_sel->Fill(m2eta[nsign],m1eta[nsign],weight[nsign]);
                    h_eta_avg_Deta_mu4_mu4noL1_w_single_b_sig_sel->Fill(deta[nsign],etaavg[nsign],weight[nsign]);
                    h_eta_avg_Dphi_mu4_mu4noL1_w_single_b_sig_sel->Fill(dphi[nsign],etaavg[nsign],weight[nsign]);
                    h_minv_zoomin_mu4_mu4noL1_w_single_b_sig_sel->Fill(minv[nsign],weight[nsign]);
                    h_pair_pt_log_mu4_mu4noL1_w_single_b_sig_sel->Fill(pair_pt[nsign],weight[nsign]);
                    h_minv_pair_pt_log_mu4_mu4noL1_w_single_b_sig_sel->Fill(pair_pt[nsign],minv[nsign],weight[nsign]);
                } // end if statement for single-b signal pairs
            } // end loop over muon index
        } // end if statement for mu4_mu4noL1 trigger

        if (pass2mu4[nsign]){

            for (int muon_ind = 0; muon_ind < 2; muon_ind++){ // loop over the two muons + check if either passes the single-muon mu4
                if (!pass_single_mu4[muon_ind]) continue;
                bool charge_sign_2nd_muon = (charge_2nd[muon_ind] > 0)? 0: 1;
                
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

                // pt-2nd-dependent histograms
                h_pt2nd_vs_q_eta_2nd_2mu4[charge_sign_2nd_muon]->Fill(q_eta_2nd[muon_ind], pt_2nd[muon_ind],weight[nsign]);
                h_pt2nd_vs_phi2nd_2mu4[charge_sign_2nd_muon]->Fill(phi_2nd[muon_ind], pt_2nd[muon_ind], weight[nsign]);
                h_phi2nd_vs_q_eta_2nd_2mu4[charge_sign_2nd_muon]->Fill(q_eta_2nd[muon_ind], phi_2nd[muon_ind], weight[nsign]);
                h_DR_zoomin_vs_pt2nd_2mu4[nsign]->Fill(pt_2nd[muon_ind],dr[nsign],weight[nsign]);
                h_DR_0_2_vs_pt2nd_2mu4[nsign]->Fill(pt_2nd[muon_ind],dr[nsign],weight[nsign]);

                // pt-1st-dependent histograms
                h_Deta_vs_pT_1st_2mu4[nsign]->Fill(pt_1st[muon_ind], deta[nsign],weight[nsign]);
                h_Deta_zoomin_vs_pT_1st_2mu4[nsign]->Fill(pt_1st[muon_ind], deta[nsign],weight[nsign]);
                h_Dphi_vs_pT_1st_2mu4[nsign]->Fill(pt_1st[muon_ind], dphi[nsign],weight[nsign]);
                h_Dphi_zoomin_vs_pT_1st_2mu4[nsign]->Fill(pt_1st[muon_ind], dphi[nsign],weight[nsign]);
                h_DR_vs_pT_1st_2mu4[nsign]->Fill(pt_1st[muon_ind], dr[nsign],weight[nsign]);
                h_DR_zoomin_vs_pT_1st_2mu4[nsign]->Fill(pt_1st[muon_ind], dr[nsign],weight[nsign]);

                // pair-observable 2D histograms
                h_pair_eta_vs_pair_pT_2mu4[nsign]->Fill(pair_pt[nsign],pair_eta[nsign],weight[nsign]);
                h_Deta_Dphi_2mu4[nsign]->Fill(dphi[nsign],deta[nsign],weight[nsign]);
                h_eta1_eta2_2mu4[nsign]->Fill(m2eta[nsign],m1eta[nsign],weight[nsign]);
                h_eta_avg_Deta_2mu4[nsign]->Fill(deta[nsign],etaavg[nsign],weight[nsign]);
                h_eta_avg_Dphi_2mu4[nsign]->Fill(dphi[nsign],etaavg[nsign],weight[nsign]);
                h_minv_zoomin_2mu4[nsign]->Fill(minv[nsign],weight[nsign]);
                h_pair_pt_log_2mu4[nsign]->Fill(pair_pt[nsign],weight[nsign]);
                h_minv_pair_pt_log_2mu4[nsign]->Fill(pair_pt[nsign],minv[nsign],weight[nsign]);

                // good-acceptance histograms
                if (pass_single_muon_good_acceptance_selection[muon_ind]){ // fill for both
                    h_Deta_2mu4_good_accept[nsign]->Fill(deta[nsign],weight[nsign]);
                    h_Deta_zoomin_2mu4_good_accept[nsign]->Fill(deta[nsign],weight[nsign]);
                    h_Dphi_2mu4_good_accept[nsign]->Fill(dphi[nsign],weight[nsign]);
                    h_Dphi_zoomin_2mu4_good_accept[nsign]->Fill(dphi[nsign],weight[nsign]);
                    h_DR_2mu4_good_accept[nsign]->Fill(dr[nsign],weight[nsign]);
                    h_DR_zoomin_2mu4_good_accept[nsign]->Fill(dr[nsign],weight[nsign]);
                    h_minv_zoomin_2mu4_good_accept[nsign]->Fill(minv[nsign],weight[nsign]);
                    h_pair_pt_log_2mu4_good_accept[nsign]->Fill(pair_pt[nsign],weight[nsign]);    
                }

                // well-separated histograms
                if (passSeparated[nsign]){
                    h_Deta_2mu4_sepr[nsign]->Fill(deta[nsign],weight[nsign]);
                    h_Deta_zoomin_2mu4_sepr[nsign]->Fill(deta[nsign],weight[nsign]);
                    h_Deta_Dphi_2mu4_sepr[nsign]->Fill(dphi[nsign],deta[nsign],weight[nsign]);
                    h_eta_avg_Deta_2mu4_sepr[nsign]->Fill(deta[nsign],etaavg[nsign],weight[nsign]);
                
                    h_Dphi_2mu4_sepr[nsign]->Fill(dphi[nsign],weight[nsign]);
                    h_Dphi_zoomin_2mu4_sepr[nsign]->Fill(dphi[nsign],weight[nsign]);
                    h_DR_2mu4_sepr[nsign]->Fill(dr[nsign],weight[nsign]);
                    h_DR_zoomin_2mu4_sepr[nsign]->Fill(dr[nsign],weight[nsign]);
                    
                    h_pt2nd_vs_q_eta_2nd_2mu4_sepr[charge_sign_2nd_muon]->Fill(q_eta_2nd[muon_ind], pt_2nd[muon_ind],weight[nsign]);
                    h_pt2nd_vs_phi2nd_2mu4_sepr[charge_sign_2nd_muon]->Fill(phi_2nd[muon_ind], pt_2nd[muon_ind],weight[nsign]);
                    h_pt2nd_vs_q_eta_2nd_vs_phi2nd_2mu4_sepr[charge_sign_2nd_muon]->Fill(phi_2nd[muon_ind], q_eta_2nd[muon_ind], pt_2nd[muon_ind], weight[nsign]);
            
                    h_pair_eta_vs_pair_pT_2mu4_sepr[nsign]->Fill(pair_pt[nsign],pair_eta[nsign],weight[nsign]);
                    h_eta1_eta2_2mu4_sepr[nsign]->Fill(m2eta[nsign],m1eta[nsign],weight[nsign]);
                    h_eta_avg_Dphi_2mu4_sepr[nsign]->Fill(dphi[nsign],etaavg[nsign],weight[nsign]);
                    h_minv_zoomin_2mu4_sepr[nsign]->Fill(minv[nsign],weight[nsign]);
                    h_pair_pt_log_2mu4_sepr[nsign]->Fill(pair_pt[nsign],weight[nsign]);
                    h_minv_pair_pt_log_2mu4_sepr[nsign]->Fill(pair_pt[nsign],minv[nsign],weight[nsign]);
                } // end if statement for well-separated pairs

                // signal-selected-pair histograms
                if (pass_single_b_signal_selection){
                    h_Deta_2mu4_w_single_b_sig_sel->Fill(deta[nsign],weight[nsign]);
                    h_Deta_zoomin_2mu4_w_single_b_sig_sel->Fill(deta[nsign],weight[nsign]);
                    h_Dphi_2mu4_w_single_b_sig_sel->Fill(dphi[nsign],weight[nsign]);
                    h_Dphi_zoomin_2mu4_w_single_b_sig_sel->Fill(dphi[nsign],weight[nsign]);
                    h_DR_2mu4_w_single_b_sig_sel->Fill(dr[nsign],weight[nsign]);
                    h_DR_zoomin_2mu4_w_single_b_sig_sel->Fill(dr[nsign],weight[nsign]);
                    h_DR_0_2_2mu4_w_single_b_sig_sel->Fill(dr[nsign],weight[nsign]);

                    h_pt2nd_vs_q_eta_2nd_2mu4_w_single_b_sig_sel[charge_sign_2nd_muon]->Fill(q_eta_2nd[muon_ind], pt_2nd[muon_ind],weight[nsign]);
                    h_pt2nd_vs_phi2nd_2mu4_w_single_b_sig_sel[charge_sign_2nd_muon]->Fill(phi_2nd[muon_ind], pt_2nd[muon_ind],weight[nsign]);
                    h_pt2nd_vs_q_eta_2nd_vs_phi2nd_2mu4_w_single_b_sig_sel[charge_sign_2nd_muon]->Fill(phi_2nd[muon_ind], q_eta_2nd[muon_ind], pt_2nd[muon_ind], weight[nsign]);
                    h_DR_zoomin_vs_pt2nd_2mu4_w_single_b_sig_sel->Fill(pt_2nd[muon_ind],dr[nsign],weight[nsign]);
                    h_DR_0_2_vs_pt2nd_2mu4_w_single_b_sig_sel->Fill(pt_2nd[muon_ind],dr[nsign],weight[nsign]);

                    h_pair_eta_vs_pair_pT_2mu4_w_single_b_sig_sel->Fill(pair_pt[nsign],pair_eta[nsign],weight[nsign]);
                    h_Deta_Dphi_2mu4_w_single_b_sig_sel->Fill(dphi[nsign],deta[nsign],weight[nsign]);
                    h_eta1_eta2_2mu4_w_single_b_sig_sel->Fill(m2eta[nsign],m1eta[nsign],weight[nsign]);
                    h_eta_avg_Deta_2mu4_w_single_b_sig_sel->Fill(deta[nsign],etaavg[nsign],weight[nsign]);
                    h_eta_avg_Dphi_2mu4_w_single_b_sig_sel->Fill(dphi[nsign],etaavg[nsign],weight[nsign]);
                    h_minv_zoomin_2mu4_w_single_b_sig_sel->Fill(minv[nsign],weight[nsign]);
                    h_pair_pt_log_2mu4_w_single_b_sig_sel->Fill(pair_pt[nsign],weight[nsign]);
                    h_minv_pair_pt_log_2mu4_w_single_b_sig_sel->Fill(pair_pt[nsign],minv[nsign],weight[nsign]);
                } // end if statement for single-b signal pairs
            } // end loop over muon index
        } // end if statement for 2mu4 trigger
    } // end if statement - trigger mode = 1 (MUON-PAIR trigger | mu4)
    else if (trigger_mode == 0){

    } // end if statement - trigger mode = 0 (mu4 | MB)
}

void MuonPairPlottingPP::FillTrigEffcyHistsInvWeightedbySingleMuonEffcies(int nsign){
    bool pass_single_b_signal_selection = (nsign == 1 && minv[nsign] > 1.08 && minv[nsign] < 2.9 && pair_pt[nsign] > 8);

    bool pass_single_mu4 [2] = {mu1PassSingle[nsign], mu2PassSingle[nsign]};
    double pt_2nd [2] = {m2pt[nsign], m1pt[nsign]};
    double q_eta_2nd [2] = {m2charge[nsign] * m2eta[nsign], m1charge[nsign] * m1eta[nsign]};
    double phi_2nd [2] = {m2phi[nsign], m1phi[nsign]};
    int charge_2nd [2] = {m2charge[nsign], m1charge[nsign]};

    // ---------- loop over the muon indices ----------
    for (int muon_ind = 0; muon_ind < 2; muon_ind++){ // loop over the two muons + check if either passes the single-muon mu4
        bool charge_sign_2nd_muon = (charge_2nd[muon_ind] > 0)? 0: 1;

        // ---------- require passing single mu4 & 2nd muon pT < pT_max (60GeV) ----------
        double single_muon_pt_max =
            use_3D_2nd_muon
            ? h_pt2nd_vs_q_eta_2nd_vs_phi2nd_mu4_mu4noL1_sepr_divided[charge_sign_2nd_muon]->GetZaxis()->GetXmax()
            : h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sepr_divided[charge_sign_2nd_muon]->GetYaxis()->GetXmax();

        if (!pass_single_mu4[muon_ind]) continue;
        if (pt_2nd[muon_ind] >= single_muon_pt_max) continue; // skip to next muon without filling the inversed-weighted-by-single-muon-efficiency histograms

        float eff_2nd_muon_mu4_mu4noL1 = EvaluateSingleMuonEffcy(charge_sign_2nd_muon, "mu4_mu4noL1", pt_2nd[muon_ind], q_eta_2nd[muon_ind], phi_2nd[muon_ind]);
        float eff_2nd_muon_2mu4 = EvaluateSingleMuonEffcy(charge_sign_2nd_muon, "2mu4", pt_2nd[muon_ind], q_eta_2nd[muon_ind], phi_2nd[muon_ind]);
        // ---------- fill histograms ----------
        if (eff_2nd_muon_mu4_mu4noL1 > 0){ // datapoints exist in the 2nd-muon kinematics bin to evaluate single-muon trigger efficiency
            
            h_Deta_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom[nsign]->Fill(deta[nsign],weight[nsign]);
            h_Dphi_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom[nsign]->Fill(dphi[nsign],weight[nsign]);
            h_DR_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom[nsign]->Fill(dr[nsign],weight[nsign]);
            h_DR_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom[nsign]->Fill(dr[nsign],weight[nsign]);
            h_DR_0_2_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom[nsign]->Fill(dr[nsign],weight[nsign]);
            h_minv_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom[nsign]->Fill(minv[nsign],weight[nsign]);
            h_pair_pt_log_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom[nsign]->Fill(pair_pt[nsign],weight[nsign]);
            h_DR_zoomin_vs_pair_pt_log_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom[nsign]->Fill(pair_pt[nsign],dr[nsign],weight[nsign] * 1. / eff_2nd_muon_mu4_mu4noL1);
            

            if (passmu4mu4noL1[nsign]){
                h_Deta_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy[nsign]->Fill(deta[nsign],weight[nsign] * 1. / eff_2nd_muon_mu4_mu4noL1);
                h_Dphi_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy[nsign]->Fill(dphi[nsign],weight[nsign] * 1. / eff_2nd_muon_mu4_mu4noL1);
                h_DR_mu4_mu4noL1_inv_w_by_single_mu_effcy[nsign]->Fill(dr[nsign],weight[nsign] * 1. / eff_2nd_muon_mu4_mu4noL1);
                h_DR_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy[nsign]->Fill(dr[nsign],weight[nsign] * 1. / eff_2nd_muon_mu4_mu4noL1);
                h_DR_0_2_mu4_mu4noL1_inv_w_by_single_mu_effcy[nsign]->Fill(dr[nsign],weight[nsign] * 1. / eff_2nd_muon_mu4_mu4noL1);
                h_minv_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy[nsign]->Fill(minv[nsign],weight[nsign] * 1. / eff_2nd_muon_mu4_mu4noL1);
                h_pair_pt_log_mu4_mu4noL1_inv_w_by_single_mu_effcy[nsign]->Fill(pair_pt[nsign],weight[nsign] * 1. / eff_2nd_muon_mu4_mu4noL1);
                h_DR_zoomin_vs_pair_pt_log_mu4_mu4noL1_inv_w_by_single_mu_effcy[nsign]->Fill(pair_pt[nsign],dr[nsign],weight[nsign] * 1. / eff_2nd_muon_mu4_mu4noL1);
            }
        }

        if (pass_single_b_signal_selection && eff_2nd_muon_mu4_mu4noL1 > 0){
            
            h_Deta_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom_w_single_b_sig_sel->Fill(deta[nsign],weight[nsign]);
            h_Dphi_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom_w_single_b_sig_sel->Fill(dphi[nsign],weight[nsign]);
            h_DR_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom_w_single_b_sig_sel->Fill(dr[nsign],weight[nsign]);
            h_DR_0_2_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom_w_single_b_sig_sel->Fill(dr[nsign],weight[nsign]);
            h_minv_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom_w_single_b_sig_sel->Fill(minv[nsign],weight[nsign]);
            h_pair_pt_log_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom_w_single_b_sig_sel->Fill(pair_pt[nsign],weight[nsign]);
            h_DR_zoomin_vs_pair_pt_log_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom_w_single_b_sig_sel->Fill(pair_pt[nsign],dr[nsign],weight[nsign] * 1. / eff_2nd_muon_mu4_mu4noL1);

            if (passmu4mu4noL1[nsign]){
                h_Deta_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_w_single_b_sig_sel->Fill(deta[nsign],weight[nsign] * 1. / eff_2nd_muon_mu4_mu4noL1);
                h_Dphi_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_w_single_b_sig_sel->Fill(dphi[nsign],weight[nsign] * 1. / eff_2nd_muon_mu4_mu4noL1);
                h_DR_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_w_single_b_sig_sel->Fill(dr[nsign],weight[nsign] * 1. / eff_2nd_muon_mu4_mu4noL1);
                h_DR_0_2_mu4_mu4noL1_inv_w_by_single_mu_effcy_w_single_b_sig_sel->Fill(dr[nsign],weight[nsign] * 1. / eff_2nd_muon_mu4_mu4noL1);
                h_minv_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_w_single_b_sig_sel->Fill(minv[nsign],weight[nsign] * 1. / eff_2nd_muon_mu4_mu4noL1);
                h_pair_pt_log_mu4_mu4noL1_inv_w_by_single_mu_effcy_w_single_b_sig_sel->Fill(pair_pt[nsign],weight[nsign] * 1. / eff_2nd_muon_mu4_mu4noL1);
                h_DR_zoomin_vs_pair_pt_log_mu4_mu4noL1_inv_w_by_single_mu_effcy_w_single_b_sig_sel->Fill(pair_pt[nsign],dr[nsign],weight[nsign] * 1. / eff_2nd_muon_mu4_mu4noL1);
            }
        }


        if (eff_2nd_muon_2mu4 > 0){

            h_Deta_zoomin_2mu4_inv_w_by_single_mu_effcy_denom[nsign]->Fill(deta[nsign],weight[nsign]);
            h_Dphi_zoomin_2mu4_inv_w_by_single_mu_effcy_denom[nsign]->Fill(dphi[nsign],weight[nsign]);
            h_DR_2mu4_inv_w_by_single_mu_effcy_denom[nsign]->Fill(dr[nsign],weight[nsign]);
            h_DR_zoomin_2mu4_inv_w_by_single_mu_effcy_denom[nsign]->Fill(dr[nsign],weight[nsign]);
            h_DR_0_2_2mu4_inv_w_by_single_mu_effcy_denom[nsign]->Fill(dr[nsign],weight[nsign]);
            h_minv_zoomin_2mu4_inv_w_by_single_mu_effcy_denom[nsign]->Fill(minv[nsign],weight[nsign]);
            h_pair_pt_log_2mu4_inv_w_by_single_mu_effcy_denom[nsign]->Fill(pair_pt[nsign],weight[nsign]);
            h_DR_zoomin_vs_pair_pt_log_2mu4_inv_w_by_single_mu_effcy_denom[nsign]->Fill(pair_pt[nsign],dr[nsign],weight[nsign] * 1. / eff_2nd_muon_mu4_mu4noL1);

            if (pass2mu4[nsign]){
                h_Deta_zoomin_2mu4_inv_w_by_single_mu_effcy[nsign]->Fill(deta[nsign],weight[nsign] * 1. / eff_2nd_muon_2mu4);
                h_Dphi_zoomin_2mu4_inv_w_by_single_mu_effcy[nsign]->Fill(dphi[nsign],weight[nsign] * 1. / eff_2nd_muon_2mu4);
                h_DR_2mu4_inv_w_by_single_mu_effcy[nsign]->Fill(dr[nsign],weight[nsign] * 1. / eff_2nd_muon_2mu4);
                h_DR_zoomin_2mu4_inv_w_by_single_mu_effcy[nsign]->Fill(dr[nsign],weight[nsign] * 1. / eff_2nd_muon_2mu4);
                h_DR_0_2_2mu4_inv_w_by_single_mu_effcy[nsign]->Fill(dr[nsign],weight[nsign] * 1. / eff_2nd_muon_2mu4);
                h_minv_zoomin_2mu4_inv_w_by_single_mu_effcy[nsign]->Fill(minv[nsign],weight[nsign] * 1. / eff_2nd_muon_2mu4);
                h_pair_pt_log_2mu4_inv_w_by_single_mu_effcy[nsign]->Fill(pair_pt[nsign],weight[nsign] * 1. / eff_2nd_muon_2mu4);
                h_DR_zoomin_vs_pair_pt_log_2mu4_inv_w_by_single_mu_effcy[nsign]->Fill(pair_pt[nsign],dr[nsign],weight[nsign] * 1. / eff_2nd_muon_mu4_mu4noL1);
            }
        }

        if (pass_single_b_signal_selection && eff_2nd_muon_2mu4 > 0){

            h_Deta_zoomin_2mu4_inv_w_by_single_mu_effcy_denom_w_single_b_sig_sel->Fill(deta[nsign],weight[nsign]);
            h_Dphi_zoomin_2mu4_inv_w_by_single_mu_effcy_denom_w_single_b_sig_sel->Fill(dphi[nsign],weight[nsign]);
            h_DR_zoomin_2mu4_inv_w_by_single_mu_effcy_denom_w_single_b_sig_sel->Fill(dr[nsign],weight[nsign]);
            h_DR_0_2_2mu4_inv_w_by_single_mu_effcy_denom_w_single_b_sig_sel->Fill(dr[nsign],weight[nsign]);
            h_minv_zoomin_2mu4_inv_w_by_single_mu_effcy_denom_w_single_b_sig_sel->Fill(minv[nsign],weight[nsign]);
            h_pair_pt_log_2mu4_inv_w_by_single_mu_effcy_denom_w_single_b_sig_sel->Fill(pair_pt[nsign],weight[nsign]);
            h_DR_zoomin_vs_pair_pt_log_2mu4_inv_w_by_single_mu_effcy_denom_w_single_b_sig_sel->Fill(pair_pt[nsign],dr[nsign],weight[nsign] * 1. / eff_2nd_muon_mu4_mu4noL1);

            if (pass2mu4[nsign]){
                h_Deta_zoomin_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel->Fill(deta[nsign],weight[nsign] * 1. / eff_2nd_muon_2mu4);
                h_Dphi_zoomin_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel->Fill(dphi[nsign],weight[nsign] * 1. / eff_2nd_muon_2mu4);
                h_DR_zoomin_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel->Fill(dr[nsign],weight[nsign] * 1. / eff_2nd_muon_2mu4);
                h_DR_0_2_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel->Fill(dr[nsign],weight[nsign] * 1. / eff_2nd_muon_2mu4);
                h_minv_zoomin_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel->Fill(minv[nsign],weight[nsign] * 1. / eff_2nd_muon_2mu4);
                h_pair_pt_log_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel->Fill(pair_pt[nsign],weight[nsign] * 1. / eff_2nd_muon_2mu4);
                h_DR_zoomin_vs_pair_pt_log_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel->Fill(pair_pt[nsign],dr[nsign],weight[nsign] * 1. / eff_2nd_muon_mu4_mu4noL1);
            }
        }
    } // end loop over the two muons
}


// Build 12 graphs for DR variables and write them out.
// Uses a local analog of clipNumeratorIfInvW and calls it before BayesDivide.
void MuonPairPlottingPP::MakeAndWriteDRTrigEffGraphs()
{

    auto make_and_write = [&](TH1* num_src, TH1* den_src){
        if (!num_src || !den_src) return (TH1*)nullptr;

        // Work on clones so originals remain untouched
        TH1* num((TH1*)num_src->Clone(Form("%s_divided", num_src->GetName())));
        TH1* den((TH1*)den_src->Clone(Form("%s_divided", den_src->GetName())));

        num->Divide(den);

        num->Write();
        return num;
    };

    // helper for projection, making & writing of TEfficiency graphs
    auto proj_make_and_write = [&](TH2* hNum2D, TH2* hDen2D, bool projy = true, int firstbin = 1, int lastbin = -1, std::string proj_range_str = ""){
        if (!hNum2D || !hDen2D){
            cout << "Either the num or the den TH2D* is nullptr!" << endl;
            return (TH1*)nullptr;
        }

        cout << "names: num: " << hNum2D->GetName() << "; den: " << hDen2D->GetName() << endl;
        cout << "2d numerator #events: " << hNum2D->GetEntries() << "; den #events: " << hDen2D->GetEntries() << endl;
        cout << "h_DR_zoomin_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel #events: " << h_DR_zoomin_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel->GetEntries() << endl;

        // suffix that captures projection axis & range
        std::string proj_suffix = projy? "_py" : "_px";
        if (proj_range_str != "") proj_suffix += "_" + proj_range_str;

        std::unique_ptr<TH1> hNum1D(
            projy ? hNum2D->ProjectionY(Form("%s%s", hNum2D->GetName(), proj_suffix.c_str()), firstbin, lastbin, "e")
                  : hNum2D->ProjectionX(Form("%s%s", hNum2D->GetName(), proj_suffix.c_str()), firstbin, lastbin, "e")
        );

        std::unique_ptr<TH1> hDen1D(
            projy ? hDen2D->ProjectionY(Form("%s%s", hDen2D->GetName(), proj_suffix.c_str()), firstbin, lastbin, "e")
                  : hDen2D->ProjectionX(Form("%s%s", hDen2D->GetName(), proj_suffix.c_str()), firstbin, lastbin, "e")
        );

        // Work on clones so originals remain untouched
        TH1* num((TH1*)hNum1D->Clone(Form("%s_divided", hNum1D->GetName())));
        TH1* den((TH1*)hDen1D->Clone(Form("%s_divided", hDen1D->GetName())));

        cout << "numerator #events: " << num->GetEntries() << "; den #events: " << den->GetEntries() << endl;
        cout << "numerator bin 3 content: " << num->GetBinContent(3) << "; den bin 3 content: " << den->GetBinContent(3) << endl;
        num->Divide(den);
        cout << "numerator bin 3 content AFTER DIVIDING: " << num->GetBinContent(3) << endl;

        num->Write();
        return num;
    };

    // -------- per-sign --------
    for (int sign = 0; sign < ParamsSet::nSigns; ++sign) {
        make_and_write(h_DR_zoomin_2mu4_inv_w_by_single_mu_effcy[sign],        h_DR_zoomin_mu4[sign]);
        make_and_write(h_DR_0_2_2mu4_inv_w_by_single_mu_effcy[sign],           h_DR_0_2_mu4[sign]);
        make_and_write(h_DR_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy[sign], h_DR_zoomin_mu4[sign]);
        make_and_write(h_DR_0_2_mu4_mu4noL1_inv_w_by_single_mu_effcy[sign],    h_DR_0_2_mu4[sign]);
    }

    // -------- signal pairs --------
    make_and_write(h_DR_zoomin_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel,        h_DR_zoomin_mu4_w_single_b_sig_sel);
    make_and_write(h_DR_0_2_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel,           h_DR_0_2_mu4_w_single_b_sig_sel);
    make_and_write(h_DR_zoomin_mu4_mu4noL1_inv_w_by_single_mu_effcy_w_single_b_sig_sel, h_DR_zoomin_mu4_w_single_b_sig_sel);
    make_and_write(h_DR_0_2_mu4_mu4noL1_inv_w_by_single_mu_effcy_w_single_b_sig_sel,    h_DR_0_2_mu4_w_single_b_sig_sel);

    // -------- signed & signal pairs in pair pT ranges --------
    for (auto range : pair_pT_ranges_for_weighted_effcy_dR_fitting){
        int bin_first = bin_number(range.first, pms.pT_bins_80) + 1; // + 1 pushes into next bin (bin lower end agree with range edge if range edge founded)
        int bin_last = bin_number(range.second, pms.pT_bins_80); // bin higher end agree with range edge if range edge founded
        std::string proj_suffix = pairToSuffix(range);

        std::cout << "pair pT Bin range: " << bin_first << ", " << bin_last << std::endl;

        for (int sign = 0; sign < ParamsSet::nSigns; ++sign) {
            proj_make_and_write(h_DR_zoomin_vs_pair_pt_log_mu4_mu4noL1_inv_w_by_single_mu_effcy[sign],     h_DR_zoomin_vs_pair_pt_log_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom[sign], true, bin_first, bin_last, proj_suffix);
            proj_make_and_write(h_DR_zoomin_vs_pair_pt_log_2mu4_inv_w_by_single_mu_effcy[sign],     h_DR_zoomin_vs_pair_pt_log_2mu4_inv_w_by_single_mu_effcy_denom[sign], true, bin_first, bin_last, proj_suffix);
        }

        proj_make_and_write(h_DR_zoomin_vs_pair_pt_log_mu4_mu4noL1_inv_w_by_single_mu_effcy_w_single_b_sig_sel,     h_DR_zoomin_vs_pair_pt_log_mu4_mu4noL1_inv_w_by_single_mu_effcy_denom_w_single_b_sig_sel, true, bin_first, bin_last, proj_suffix);
        proj_make_and_write(h_DR_zoomin_vs_pair_pt_log_2mu4_inv_w_by_single_mu_effcy_w_single_b_sig_sel,     h_DR_zoomin_vs_pair_pt_log_2mu4_inv_w_by_single_mu_effcy_denom_w_single_b_sig_sel, true, bin_first, bin_last, proj_suffix);
    }

}

void MuonPairPlottingPP::MakeAndWriteSingleMuonPtTrigEffGraphs()
{
    // helper for projection, making & writing of TEfficiency graphs
    auto proj_make_and_write = [&](TH2* hNum2D, TH2* hDen2D, bool projy = true, int firstbin = 1, int lastbin = -1, std::string proj_range_str = ""){
        if (!hNum2D || !hDen2D) return (TGraphAsymmErrors*)nullptr;

        // suffix that captures projection axis & range
        std::string proj_suffix = projy? "_py" : "_px";
        if (proj_range_str != "") proj_suffix += "_" + proj_range_str;

        std::unique_ptr<TH1> hNum1D(
            projy ? hNum2D->ProjectionY(Form("%s%s", hNum2D->GetName(), proj_suffix.c_str()), firstbin, lastbin, "e")
                  : hNum2D->ProjectionX(Form("%s%s", hNum2D->GetName(), proj_suffix.c_str()), firstbin, lastbin, "e")
        );

        std::unique_ptr<TH1> hDen1D(
            projy ? hDen2D->ProjectionY(Form("%s%s", hDen2D->GetName(), proj_suffix.c_str()), firstbin, lastbin, "e")
                  : hDen2D->ProjectionX(Form("%s%s", hDen2D->GetName(), proj_suffix.c_str()), firstbin, lastbin, "e")
        );

        auto g = new TGraphAsymmErrors();
        g->BayesDivide(hNum1D.get(), hDen1D.get());

        // Name: "h_<...>" → "g_<...>" from the *numerator* histo name
        std::string n = hNum1D->GetName() ? hNum1D->GetName() : "graph";
        
        n += "_divided";
        if (n.rfind("h_", 0) == 0) n.replace(0, 2, "g_"); else n = "g_" + n;
        g->SetName(n.c_str());

        g->Write(g->GetName(), TObject::kOverwrite);
        return g;
    };

    // q-eta bins
    std::vector<double> eta_bins_trig_effcy = ParamsSet::makeEtaTrigEffcyBinning();

    for (int sign = 0; sign < 2; ++sign) {
        proj_make_and_write(h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sepr[sign],     h_pt2nd_vs_q_eta_2nd_mu4_sepr[sign]);
        proj_make_and_write(h_pt2nd_vs_q_eta_2nd_2mu4_sepr[sign],            h_pt2nd_vs_q_eta_2nd_mu4_sepr[sign]);

        for (auto range : q_eta_proj_ranges_for_single_muon_effcy_pT_fitting){
            int bin_first = bin_number(range.first, eta_bins_trig_effcy) + 1; // + 1 pushes into next bin (bin lower end agree with range edge if range edge founded)
            int bin_last = bin_number(range.second, eta_bins_trig_effcy); // bin higher end agree with range edge if range edge founded
            std::string proj_suffix = pairToSuffix(range);

            proj_make_and_write(h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_sepr[sign],     h_pt2nd_vs_q_eta_2nd_mu4_sepr[sign], true, bin_first, bin_last, proj_suffix);
            proj_make_and_write(h_pt2nd_vs_q_eta_2nd_2mu4_sepr[sign],            h_pt2nd_vs_q_eta_2nd_mu4_sepr[sign], true, bin_first, bin_last, proj_suffix);
        }
    }
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
    if (trigger_mode == 1) OpenEffcyPtFitFile();
    InitOutput();
  	InitHists();
  	ProcessData();
  	WriteOutput();

    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;;
    std::cout << "CPU time used is " << cpu_time_used << " seconds" << std::endl;

}