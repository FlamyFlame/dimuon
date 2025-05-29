#include "MuonPairPlottingPP.h"
#include "time.h"

MuonPairPlottingPP::MuonPairPlottingPP(){
    // if (mode != 1 && mode != 3){
    //     std::cout<<"Error:: Mode has to be 1 (no binning) or 3 (binning by pT),  quitting"<<std::endl;
    //     throw std::exception();
    // }
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
            outFile = new TFile("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_run2/histograms_real_pairs_pp_2017.root","recreate");
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

    int neta_bins_trig_effcy = 40;
    int nDphi_bins_trig_effcy = 32;
    int nDeta_bins_trig_effcy = 40;
    int nDR_bins_trig_effcy = 40;
    int nminv_bins_trig_effcy = 40;
    int nDR_deta_dphi_zoomin_bins_trig_effcy = 20;


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

    if (mode == 1){
        // for (unsigned int isign = 0; isign < ParamsSet::nSigns; isign++){
                
        // for (unsigned int ipt = 0; ipt < ParamsSet::nPtBins; ipt++){
        // for (unsigned int ictr = 0; ictr < ParamsSet::nCtrBins; ictr++){
        // for (unsigned int idphi = 0; idphi < ParamsSet::ndphiRegions; idphi++){
        // for (unsigned int ideta = 0; ideta < ParamsSet::ndetaRegions; ideta++){
        // for (unsigned int igap = 0; igap < ParamsSet::nGapCuts; igap++){
        // for (unsigned int iphoto = 0; iphoto < ParamsSet::nPhotoProdCuts; iphoto++){
        for (unsigned int isign = 0; isign < ParamsSet::nSigns; isign++){
            h_Deta_mu4[isign] = new TH1D(Form("h_Deta_mu4_sign%d",isign+1),";#Delta#eta;",nDeta_bins_trig_effcy,-4.8,4.8);
            h_Deta_zoomin_mu4[isign] = new TH1D(Form("h_Deta_zoomin_mu4_sign%d",isign+1),";#Delta#eta;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
            h_Dphi_mu4[isign] = new TH1D(Form("h_Dphi_mu4_sign%d",isign+1),";#Delta#phi;",nDphi_bins_trig_effcy,-pms.PI,pms.PI);
            h_Dphi_zoomin_mu4[isign] = new TH1D(Form("h_Dphi_zoomin_mu4_sign%d",isign+1),";#Delta#phi;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
            h_DR_mu4[isign] = new TH1D(Form("h_DR_mu4_sign%d",isign+1),";#DeltaR;",nDR_bins_trig_effcy,0,pms.deltaR_thrsh[2]);
            h_DR_zoomin_mu4[isign] = new TH1D(Form("h_DR_zoomin_mu4_sign%d",isign+1),";#DeltaR;",nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
            h_minv_zoomin_mu4[isign] = new TH1D(Form("h_minv_zoomin_mu4_sign%d",isign+1),";m_{#mu#mu} [GeV];",nminv_bins_trig_effcy,0,3.);
            h_pair_pt_log_mu4[isign] = new TH1D(Form("h_pair_pt_log_mu4_sign%d",isign+1),";p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
            h_pt2nd_mu4[isign] = new TH1D(Form("h_pt2nd_mu4_sign%d",isign+1),";p_{T,2nd} [GeV];",int(pms.pT_bins_40.size() - 1), pms.pT_bins_40.data());

            h_pt2nd_vs_q_eta_2nd_mu4[isign] = new TH2D(Form("h_pt2nd_vs_q_eta_2nd_mu4_sign%d",isign+1),";q*#eta;p_{T} [GeV]",neta_bins_trig_effcy,-2.4,2.4, int(pms.pT_bins_40.size() - 1), pms.pT_bins_40.data());
            h_pair_eta_vs_pair_pT_mu4[isign] = new TH2D(Form("h_pair_eta_vs_pair_pT_mu4_sign%d",isign+1),";p_{T}^{pair} [GeV];#eta^{pair}",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),neta_bins_trig_effcy,-2.4,2.4);
            h_Deta_Dphi_mu4[isign] = new TH2D(Form("h_Deta_Dphi_mu4_sign%d",isign+1),";#Delta#phi;#Delta#eta", nDphi_bins_trig_effcy,-pms.PI,pms.PI,nDeta_bins_trig_effcy,-4.8,4.8);
            h_eta1_eta2_mu4[isign] = new TH2D(Form("h_eta1_eta2_mu4_sign%d",isign+1),";#eta_{sublead};#eta_{lead}",neta_bins_trig_effcy,-2.4,2.4, neta_bins_trig_effcy,-2.4,2.4);
            h_eta_avg_Deta_mu4[isign] = new TH2D(Form("h_eta_avg_Deta_mu4_sign%d",isign+1),";#bar{#eta};#Delta#eta", neta_bins_trig_effcy,-2.4,2.4,nDeta_bins_trig_effcy,-4.8,4.8);
            h_eta_avg_Dphi_mu4[isign] = new TH2D(Form("h_eta_avg_Dphi_mu4_sign%d",isign+1),";#Delta#phi;#bar{#eta}", nDphi_bins_trig_effcy,-pms.PI,pms.PI,neta_bins_trig_effcy,-2.4,2.4);
            h_minv_pair_pt_log_mu4[isign] = new TH2D(Form("h_minv_pair_pt_log_mu4_sign%d",isign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),nminv_bins_log,minv_bins_log[isign]);
            // h_minv_pair_pt_log_mu4[isign] = new TH2D(Form("h_minv_pair_pt_log_mu4_sign%d",isign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",pms.npairPT_bins,pms.pairPTBins[isign][2],nminv_bins_log,minv_bins_log[isign]);

            h_Deta_mu4_w_sig_sel[isign] = new TH1D(Form("h_Deta_mu4_w_sig_sel_sign%d",isign+1),";#Delta#eta;",nDeta_bins_trig_effcy,-4.8,4.8);
            h_Deta_zoomin_mu4_w_sig_sel[isign] = new TH1D(Form("h_Deta_zoomin_mu4_w_sig_sel_sign%d",isign+1),";#Delta#eta;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
            h_Dphi_mu4_w_sig_sel[isign] = new TH1D(Form("h_Dphi_mu4_w_sig_sel_sign%d",isign+1),";#Delta#phi;",nDphi_bins_trig_effcy,-pms.PI,pms.PI);
            h_Dphi_zoomin_mu4_w_sig_sel[isign] = new TH1D(Form("h_Dphi_zoomin_mu4_w_sig_sel_sign%d",isign+1),";#Delta#phi;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
            h_DR_mu4_w_sig_sel[isign] = new TH1D(Form("h_DR_mu4_w_sig_sel_sign%d",isign+1),";#DeltaR;",nDR_bins_trig_effcy,0,pms.deltaR_thrsh[2]);
            h_DR_zoomin_mu4_w_sig_sel[isign] = new TH1D(Form("h_DR_zoomin_mu4_w_sig_sel_sign%d",isign+1),";#DeltaR;",nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
            h_minv_zoomin_mu4_w_sig_sel[isign] = new TH1D(Form("h_minv_zoomin_mu4_w_sig_sel_sign%d",isign+1),";m_{#mu#mu} [GeV];",nminv_bins_trig_effcy,0,3.);
            h_pair_pt_log_mu4_w_sig_sel[isign] = new TH1D(Form("h_pair_pt_log_mu4_w_sig_sel_sign%d",isign+1),";p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
            h_pt2nd_mu4_w_sig_sel[isign] = new TH1D(Form("h_pt2nd_mu4_w_sig_sel_sign%d",isign+1),";p_{T,2nd} [GeV];",int(pms.pT_bins_40.size() - 1), pms.pT_bins_40.data());

            h_pt2nd_vs_q_eta_2nd_mu4_w_sig_sel[isign] = new TH2D(Form("h_pt2nd_vs_q_eta_2nd_mu4_w_sig_sel_sign%d",isign+1),";q*#eta;p_{T} [GeV]",neta_bins_trig_effcy,-2.4,2.4, int(pms.pT_bins_40.size() - 1), pms.pT_bins_40.data());
            h_pair_eta_vs_pair_pT_mu4_w_sig_sel[isign] = new TH2D(Form("h_pair_eta_vs_pair_pT_mu4_w_sig_sel_sign%d",isign+1),";p_{T}^{pair} [GeV];#eta^{pair}",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),neta_bins_trig_effcy,-2.4,2.4);
            h_Deta_Dphi_mu4_w_sig_sel[isign] = new TH2D(Form("h_Deta_Dphi_mu4_w_sig_sel_sign%d",isign+1),";#Delta#phi;#Delta#eta", nDphi_bins_trig_effcy,-pms.PI,pms.PI,nDeta_bins_trig_effcy,-4.8,4.8);
            h_eta1_eta2_mu4_w_sig_sel[isign] = new TH2D(Form("h_eta1_eta2_mu4_w_sig_sel_sign%d",isign+1),";#eta_{sublead};#eta_{lead}",neta_bins_trig_effcy,-2.4,2.4, neta_bins_trig_effcy,-2.4,2.4);
            h_eta_avg_Deta_mu4_w_sig_sel[isign] = new TH2D(Form("h_eta_avg_Deta_mu4_w_sig_sel_sign%d",isign+1),";#bar{#eta};#Delta#eta", neta_bins_trig_effcy,-2.4,2.4,nDeta_bins_trig_effcy,-4.8,4.8);
            h_eta_avg_Dphi_mu4_w_sig_sel[isign] = new TH2D(Form("h_eta_avg_Dphi_mu4_w_sig_sel_sign%d",isign+1),";#Delta#phi;#bar{#eta}", nDphi_bins_trig_effcy,-pms.PI,pms.PI,neta_bins_trig_effcy,-2.4,2.4);
            h_minv_pair_pt_log_mu4_w_sig_sel[isign] = new TH2D(Form("h_minv_pair_pt_log_mu4_w_sig_sel_sign%d",isign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),nminv_bins_log,minv_bins_log[isign]);
            // h_minv_pair_pt_log_mu4_w_sig_sel[isign] = new TH2D(Form("h_minv_pair_pt_log_mu4_w_sig_sel_sign%d",isign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",pms.npairPT_bins,pms.pairPTBins[isign][2],nminv_bins_log,minv_bins_log[isign]);

            h_Deta_mu4_mu4noL1_given_mu4[isign] = new TH1D(Form("h_Deta_mu4_mu4noL1_given_mu4_sign%d",isign+1),";#Delta#eta;",nDeta_bins_trig_effcy,-4.8,4.8);
            h_Deta_zoomin_mu4_mu4noL1_given_mu4[isign] = new TH1D(Form("h_Deta_zoomin_mu4_mu4noL1_given_mu4_sign%d",isign+1),";#Delta#eta;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
            h_Dphi_mu4_mu4noL1_given_mu4[isign] = new TH1D(Form("h_Dphi_mu4_mu4noL1_given_mu4_sign%d",isign+1),";#Delta#phi;",nDphi_bins_trig_effcy,-pms.PI,pms.PI);
            h_Dphi_zoomin_mu4_mu4noL1_given_mu4[isign] = new TH1D(Form("h_Dphi_zoomin_mu4_mu4noL1_given_mu4_sign%d",isign+1),";#Delta#phi;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
            h_DR_mu4_mu4noL1_given_mu4[isign] = new TH1D(Form("h_DR_mu4_mu4noL1_given_mu4_sign%d",isign+1),";#DeltaR;",nDR_bins_trig_effcy,0,pms.deltaR_thrsh[2]);
            h_DR_zoomin_mu4_mu4noL1_given_mu4[isign] = new TH1D(Form("h_DR_zoomin_mu4_mu4noL1_given_mu4_sign%d",isign+1),";#DeltaR;",nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
            h_minv_zoomin_mu4_mu4noL1_given_mu4[isign] = new TH1D(Form("h_minv_zoomin_mu4_mu4noL1_given_mu4_sign%d",isign+1),";m_{#mu#mu} [GeV];",nminv_bins_trig_effcy,0,3.);
            h_pair_pt_log_mu4_mu4noL1_given_mu4[isign] = new TH1D(Form("h_pair_pt_log_mu4_mu4noL1_given_mu4_sign%d",isign+1),";p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
            h_pt2nd_mu4_mu4noL1_given_mu4[isign] = new TH1D(Form("h_pt2nd_mu4_mu4noL1_given_mu4_sign%d",isign+1),";p_{T,2nd} [GeV];",int(pms.pT_bins_40.size() - 1), pms.pT_bins_40.data());

            h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_given_mu4[isign] = new TH2D(Form("h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_given_mu4_sign%d",isign+1),";q*#eta;p_{T} [GeV]",neta_bins_trig_effcy,-2.4,2.4, int(pms.pT_bins_40.size() - 1), pms.pT_bins_40.data());
            h_pair_eta_vs_pair_pT_mu4_mu4noL1_given_mu4[isign] = new TH2D(Form("h_pair_eta_vs_pair_pT_mu4_mu4noL1_given_mu4_sign%d",isign+1),";p_{T}^{pair} [GeV];#eta^{pair}",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),neta_bins_trig_effcy,-2.4,2.4);
            h_Deta_Dphi_mu4_mu4noL1_given_mu4[isign] = new TH2D(Form("h_Deta_Dphi_mu4_mu4noL1_given_mu4_sign%d",isign+1),";#Delta#phi;#Delta#eta", nDphi_bins_trig_effcy,-pms.PI,pms.PI,nDeta_bins_trig_effcy,-4.8,4.8);
            h_eta1_eta2_mu4_mu4noL1_given_mu4[isign] = new TH2D(Form("h_eta1_eta2_mu4_mu4noL1_given_mu4_sign%d",isign+1),";#eta_{sublead};#eta_{lead}",neta_bins_trig_effcy,-2.4,2.4, neta_bins_trig_effcy,-2.4,2.4);
            h_eta_avg_Deta_mu4_mu4noL1_given_mu4[isign] = new TH2D(Form("h_eta_avg_Deta_mu4_mu4noL1_given_mu4_sign%d",isign+1),";#bar{#eta};#Delta#eta", neta_bins_trig_effcy,-2.4,2.4,nDeta_bins_trig_effcy,-4.8,4.8);
            h_eta_avg_Dphi_mu4_mu4noL1_given_mu4[isign] = new TH2D(Form("h_eta_avg_Dphi_mu4_mu4noL1_given_mu4_sign%d",isign+1),";#Delta#phi;#bar{#eta}", nDphi_bins_trig_effcy,-pms.PI,pms.PI,neta_bins_trig_effcy,-2.4,2.4);
            h_minv_pair_pt_log_mu4_mu4noL1_given_mu4[isign] = new TH2D(Form("h_minv_pair_pt_log_mu4_mu4noL1_given_mu4_sign%d",isign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),nminv_bins_log,minv_bins_log[isign]);
            // h_minv_pair_pt_log_mu4_mu4noL1_given_mu4[isign] = new TH2D(Form("h_minv_pair_pt_log_mu4_mu4noL1_given_mu4_sign%d",isign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",pms.npairPT_bins,pms.pairPTBins[isign][2],nminv_bins_log,minv_bins_log[isign]);

            h_Deta_2mu4_given_mu4[isign] = new TH1D(Form("h_Deta_2mu4_given_mu4_sign%d",isign+1),";#Delta#eta;",nDeta_bins_trig_effcy,-4.8,4.8);
            h_Deta_zoomin_2mu4_given_mu4[isign] = new TH1D(Form("h_Deta_zoomin_2mu4_given_mu4_sign%d",isign+1),";#Delta#eta;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
            h_Dphi_2mu4_given_mu4[isign] = new TH1D(Form("h_Dphi_2mu4_given_mu4_sign%d",isign+1),";#Delta#phi;",nDphi_bins_trig_effcy,-pms.PI,pms.PI);
            h_Dphi_zoomin_2mu4_given_mu4[isign] = new TH1D(Form("h_Dphi_zoomin_2mu4_given_mu4_sign%d",isign+1),";#Delta#phi;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
            h_DR_2mu4_given_mu4[isign] = new TH1D(Form("h_DR_2mu4_given_mu4_sign%d",isign+1),";#DeltaR;",nDR_bins_trig_effcy,0,pms.deltaR_thrsh[2]);
            h_DR_zoomin_2mu4_given_mu4[isign] = new TH1D(Form("h_DR_zoomin_2mu4_given_mu4_sign%d",isign+1),";#DeltaR;",nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
            h_minv_zoomin_2mu4_given_mu4[isign] = new TH1D(Form("h_minv_zoomin_2mu4_given_mu4_sign%d",isign+1),";m_{#mu#mu} [GeV];",nminv_bins_trig_effcy,0,3.);
            h_pair_pt_log_2mu4_given_mu4[isign] = new TH1D(Form("h_pair_pt_log_2mu4_given_mu4_sign%d",isign+1),";p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
            h_pt2nd_2mu4_given_mu4[isign] = new TH1D(Form("h_pt2nd_2mu4_given_mu4_sign%d",isign+1),";p_{T,2nd} [GeV];",int(pms.pT_bins_40.size() - 1), pms.pT_bins_40.data());

            h_pt2nd_vs_q_eta_2nd_2mu4_given_mu4[isign] = new TH2D(Form("h_pt2nd_vs_q_eta_2nd_2mu4_given_mu4_sign%d",isign+1),";q*#eta;p_{T} [GeV]",neta_bins_trig_effcy,-2.4,2.4, int(pms.pT_bins_40.size() - 1), pms.pT_bins_40.data());
            h_pair_eta_vs_pair_pT_2mu4_given_mu4[isign] = new TH2D(Form("h_pair_eta_vs_pair_pT_2mu4_given_mu4_sign%d",isign+1),";p_{T}^{pair} [GeV];#eta^{pair}",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),neta_bins_trig_effcy,-2.4,2.4);
            h_Deta_Dphi_2mu4_given_mu4[isign] = new TH2D(Form("h_Deta_Dphi_2mu4_given_mu4_sign%d",isign+1),";#Delta#phi;#Delta#eta", nDphi_bins_trig_effcy,-pms.PI,pms.PI,nDeta_bins_trig_effcy,-4.8,4.8);
            h_eta1_eta2_2mu4_given_mu4[isign] = new TH2D(Form("h_eta1_eta2_2mu4_given_mu4_sign%d",isign+1),";#eta_{sublead};#eta_{lead}",neta_bins_trig_effcy,-2.4,2.4, neta_bins_trig_effcy,-2.4,2.4);
            h_eta_avg_Deta_2mu4_given_mu4[isign] = new TH2D(Form("h_eta_avg_Deta_2mu4_given_mu4_sign%d",isign+1),";#bar{#eta};#Delta#eta", neta_bins_trig_effcy,-2.4,2.4,nDeta_bins_trig_effcy,-4.8,4.8);
            h_eta_avg_Dphi_2mu4_given_mu4[isign] = new TH2D(Form("h_eta_avg_Dphi_2mu4_given_mu4_sign%d",isign+1),";#Delta#phi;#bar{#eta}", nDphi_bins_trig_effcy,-pms.PI,pms.PI,neta_bins_trig_effcy,-2.4,2.4);
            h_minv_pair_pt_log_2mu4_given_mu4[isign] = new TH2D(Form("h_minv_pair_pt_log_2mu4_given_mu4_sign%d",isign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),nminv_bins_log,minv_bins_log[isign]);
            // h_minv_pair_pt_log_2mu4_given_mu4[isign] = new TH2D(Form("h_minv_pair_pt_log_2mu4_given_mu4_sign%d",isign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",pms.npairPT_bins,pms.pairPTBins[isign][2],nminv_bins_log,minv_bins_log[isign]);

            h_Deta_mu4_mu4noL1_given_mu4_w_sig_sel[isign] = new TH1D(Form("h_Deta_mu4_mu4noL1_given_mu4_w_sig_sel_sign%d",isign+1),";#Delta#eta;",nDeta_bins_trig_effcy,-4.8,4.8);
            h_Deta_zoomin_mu4_mu4noL1_given_mu4_w_sig_sel[isign] = new TH1D(Form("h_Deta_zoomin_mu4_mu4noL1_given_mu4_w_sig_sel_sign%d",isign+1),";#Delta#eta;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
            h_Dphi_mu4_mu4noL1_given_mu4_w_sig_sel[isign] = new TH1D(Form("h_Dphi_mu4_mu4noL1_given_mu4_w_sig_sel_sign%d",isign+1),";#Delta#phi;",nDphi_bins_trig_effcy,-pms.PI,pms.PI);
            h_Dphi_zoomin_mu4_mu4noL1_given_mu4_w_sig_sel[isign] = new TH1D(Form("h_Dphi_zoomin_mu4_mu4noL1_given_mu4_w_sig_sel_sign%d",isign+1),";#Delta#phi;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
            h_DR_mu4_mu4noL1_given_mu4_w_sig_sel[isign] = new TH1D(Form("h_DR_mu4_mu4noL1_given_mu4_w_sig_sel_sign%d",isign+1),";#DeltaR;",nDR_bins_trig_effcy,0,pms.deltaR_thrsh[2]);
            h_DR_zoomin_mu4_mu4noL1_given_mu4_w_sig_sel[isign] = new TH1D(Form("h_DR_zoomin_mu4_mu4noL1_given_mu4_w_sig_sel_sign%d",isign+1),";#DeltaR;",nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
            h_minv_zoomin_mu4_mu4noL1_given_mu4_w_sig_sel[isign] = new TH1D(Form("h_minv_zoomin_mu4_mu4noL1_given_mu4_w_sig_sel_sign%d",isign+1),";m_{#mu#mu} [GeV];",nminv_bins_trig_effcy,0,3.);
            h_pair_pt_log_mu4_mu4noL1_given_mu4_w_sig_sel[isign] = new TH1D(Form("h_pair_pt_log_mu4_mu4noL1_given_mu4_w_sig_sel_sign%d",isign+1),";p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
            h_pt2nd_mu4_mu4noL1_given_mu4_w_sig_sel[isign] = new TH1D(Form("h_pt2nd_mu4_mu4noL1_given_mu4_w_sig_sel_sign%d",isign+1),";p_{T,2nd} [GeV];",int(pms.pT_bins_40.size() - 1), pms.pT_bins_40.data());

            h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_given_mu4_w_sig_sel[isign] = new TH2D(Form("h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_given_mu4_w_sig_sel_sign%d",isign+1),";q*#eta;p_{T} [GeV]",neta_bins_trig_effcy,-2.4,2.4, int(pms.pT_bins_40.size() - 1), pms.pT_bins_40.data());
            h_pair_eta_vs_pair_pT_mu4_mu4noL1_given_mu4_w_sig_sel[isign] = new TH2D(Form("h_pair_eta_vs_pair_pT_mu4_mu4noL1_given_mu4_w_sig_sel_sign%d",isign+1),";p_{T}^{pair} [GeV];#eta^{pair}",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),neta_bins_trig_effcy,-2.4,2.4);
            h_Deta_Dphi_mu4_mu4noL1_given_mu4_w_sig_sel[isign] = new TH2D(Form("h_Deta_Dphi_mu4_mu4noL1_given_mu4_w_sig_sel_sign%d",isign+1),";#Delta#phi;#Delta#eta", nDphi_bins_trig_effcy,-pms.PI,pms.PI,nDeta_bins_trig_effcy,-4.8,4.8);
            h_eta1_eta2_mu4_mu4noL1_given_mu4_w_sig_sel[isign] = new TH2D(Form("h_eta1_eta2_mu4_mu4noL1_given_mu4_w_sig_sel_sign%d",isign+1),";#eta_{sublead};#eta_{lead}",neta_bins_trig_effcy,-2.4,2.4, neta_bins_trig_effcy,-2.4,2.4);
            h_eta_avg_Deta_mu4_mu4noL1_given_mu4_w_sig_sel[isign] = new TH2D(Form("h_eta_avg_Deta_mu4_mu4noL1_given_mu4_w_sig_sel_sign%d",isign+1),";#bar{#eta};#Delta#eta", neta_bins_trig_effcy,-2.4,2.4,nDeta_bins_trig_effcy,-4.8,4.8);
            h_eta_avg_Dphi_mu4_mu4noL1_given_mu4_w_sig_sel[isign] = new TH2D(Form("h_eta_avg_Dphi_mu4_mu4noL1_given_mu4_w_sig_sel_sign%d",isign+1),";#Delta#phi;#bar{#eta}", nDphi_bins_trig_effcy,-pms.PI,pms.PI,neta_bins_trig_effcy,-2.4,2.4);
            h_minv_pair_pt_log_mu4_mu4noL1_given_mu4_w_sig_sel[isign] = new TH2D(Form("h_minv_pair_pt_log_mu4_mu4noL1_given_mu4_w_sig_sel_sign%d",isign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),nminv_bins_log,minv_bins_log[isign]);
            // h_minv_pair_pt_log_mu4_mu4noL1_given_mu4_w_sig_sel[isign] = new TH2D(Form("h_minv_pair_pt_log_mu4_mu4noL1_given_mu4_w_sig_sel_sign%d",isign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",pms.npairPT_bins,pms.pairPTBins[isign][2],nminv_bins_log,minv_bins_log[isign]);

            h_Deta_2mu4_given_mu4_w_sig_sel[isign] = new TH1D(Form("h_Deta_2mu4_given_mu4_w_sig_sel_sign%d",isign+1),";#Delta#eta;",nDeta_bins_trig_effcy,-4.8,4.8);
            h_Deta_zoomin_2mu4_given_mu4_w_sig_sel[isign] = new TH1D(Form("h_Deta_zoomin_2mu4_given_mu4_w_sig_sel_sign%d",isign+1),";#Delta#eta;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
            h_Dphi_2mu4_given_mu4_w_sig_sel[isign] = new TH1D(Form("h_Dphi_2mu4_given_mu4_w_sig_sel_sign%d",isign+1),";#Delta#phi;",nDphi_bins_trig_effcy,-pms.PI,pms.PI);
            h_Dphi_zoomin_2mu4_given_mu4_w_sig_sel[isign] = new TH1D(Form("h_Dphi_zoomin_2mu4_given_mu4_w_sig_sel_sign%d",isign+1),";#Delta#phi;",nDR_deta_dphi_zoomin_bins_trig_effcy,-0.8,0.8);
            h_DR_2mu4_given_mu4_w_sig_sel[isign] = new TH1D(Form("h_DR_2mu4_given_mu4_w_sig_sel_sign%d",isign+1),";#DeltaR;",nDR_bins_trig_effcy,0,pms.deltaR_thrsh[2]);
            h_DR_zoomin_2mu4_given_mu4_w_sig_sel[isign] = new TH1D(Form("h_DR_zoomin_2mu4_given_mu4_w_sig_sel_sign%d",isign+1),";#DeltaR;",nDR_deta_dphi_zoomin_bins_trig_effcy,0,0.8);
            h_minv_zoomin_2mu4_given_mu4_w_sig_sel[isign] = new TH1D(Form("h_minv_zoomin_2mu4_given_mu4_w_sig_sel_sign%d",isign+1),";m_{#mu#mu} [GeV];",nminv_bins_trig_effcy,0,3.);
            h_pair_pt_log_2mu4_given_mu4_w_sig_sel[isign] = new TH1D(Form("h_pair_pt_log_2mu4_given_mu4_w_sig_sel_sign%d",isign+1),";p_{T}^{pair} [GeV];",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data());
            h_pt2nd_2mu4_given_mu4_w_sig_sel[isign] = new TH1D(Form("h_pt2nd_2mu4_given_mu4_w_sig_sel_sign%d",isign+1),";p_{T,2nd} [GeV];",int(pms.pT_bins_40.size() - 1), pms.pT_bins_40.data());

            h_pt2nd_vs_q_eta_2nd_2mu4_given_mu4_w_sig_sel[isign] = new TH2D(Form("h_pt2nd_vs_q_eta_2nd_2mu4_given_mu4_w_sig_sel_sign%d",isign+1),";q*#eta;p_{T} [GeV]",neta_bins_trig_effcy,-2.4,2.4, int(pms.pT_bins_40.size() - 1), pms.pT_bins_40.data());
            h_pair_eta_vs_pair_pT_2mu4_given_mu4_w_sig_sel[isign] = new TH2D(Form("h_pair_eta_vs_pair_pT_2mu4_given_mu4_w_sig_sel_sign%d",isign+1),";p_{T}^{pair} [GeV];#eta^{pair}",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),neta_bins_trig_effcy,-2.4,2.4);
            h_Deta_Dphi_2mu4_given_mu4_w_sig_sel[isign] = new TH2D(Form("h_Deta_Dphi_2mu4_given_mu4_w_sig_sel_sign%d",isign+1),";#Delta#phi;#Delta#eta", nDphi_bins_trig_effcy,-pms.PI,pms.PI,nDeta_bins_trig_effcy,-4.8,4.8);
            h_eta1_eta2_2mu4_given_mu4_w_sig_sel[isign] = new TH2D(Form("h_eta1_eta2_2mu4_given_mu4_w_sig_sel_sign%d",isign+1),";#eta_{sublead};#eta_{lead}",neta_bins_trig_effcy,-2.4,2.4, neta_bins_trig_effcy,-2.4,2.4);
            h_eta_avg_Deta_2mu4_given_mu4_w_sig_sel[isign] = new TH2D(Form("h_eta_avg_Deta_2mu4_given_mu4_w_sig_sel_sign%d",isign+1),";#bar{#eta};#Delta#eta", neta_bins_trig_effcy,-2.4,2.4,nDeta_bins_trig_effcy,-4.8,4.8);
            h_eta_avg_Dphi_2mu4_given_mu4_w_sig_sel[isign] = new TH2D(Form("h_eta_avg_Dphi_2mu4_given_mu4_w_sig_sel_sign%d",isign+1),";#Delta#phi;#bar{#eta}", nDphi_bins_trig_effcy,-pms.PI,pms.PI,neta_bins_trig_effcy,-2.4,2.4);
            h_minv_pair_pt_log_2mu4_given_mu4_w_sig_sel[isign] = new TH2D(Form("h_minv_pair_pt_log_2mu4_given_mu4_w_sig_sel_sign%d",isign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",int(pms.pT_bins_80.size() - 1), pms.pT_bins_80.data(),nminv_bins_log,minv_bins_log[isign]);
            // h_minv_pair_pt_log_2mu4_given_mu4_w_sig_sel[isign] = new TH2D(Form("h_minv_pair_pt_log_2mu4_given_mu4_w_sig_sel_sign%d",isign+1),";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]",pms.npairPT_bins,pms.pairPTBins[isign][2],nminv_bins_log,minv_bins_log[isign]);

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
}

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

    bool pass_single_b_signal_selection = (nsign == 1 && minv[nsign] > 1.08 && minv[nsign] < 2.9 && pair_pt[nsign] > 8);
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
    if (mu1PassSingle[nsign]){
        h_pt2nd_mu4[nsign]->Fill(m2pt[nsign],weight[nsign]);
        h_pt2nd_vs_q_eta_2nd_mu4[nsign]->Fill(m2charge[nsign] * m2eta[nsign], m2pt[nsign],weight[nsign]);
    }
    if (mu2PassSingle[nsign]){
        h_pt2nd_mu4[nsign]->Fill(m1pt[nsign],weight[nsign]);
        h_pt2nd_vs_q_eta_2nd_mu4[nsign]->Fill(m1charge[nsign] * m1eta[nsign], m1pt[nsign],weight[nsign]);
    }

    h_pair_eta_vs_pair_pT_mu4[nsign]->Fill(pair_pt[nsign],pair_eta[nsign],weight[nsign]);
    h_Deta_Dphi_mu4[nsign]->Fill(dphi[nsign],deta[nsign],weight[nsign]);
    h_eta1_eta2_mu4[nsign]->Fill(m2eta[nsign],m1eta[nsign],weight[nsign]);
    h_eta_avg_Deta_mu4[nsign]->Fill(deta[nsign],etaavg[nsign],weight[nsign]);
    h_eta_avg_Dphi_mu4[nsign]->Fill(dphi[nsign],etaavg[nsign],weight[nsign]);
    h_minv_zoomin_mu4[nsign]->Fill(minv[nsign],weight[nsign]);
    h_pair_pt_log_mu4[nsign]->Fill(pair_pt[nsign],weight[nsign]);
    h_minv_pair_pt_log_mu4[nsign]->Fill(pair_pt[nsign],minv[nsign],weight[nsign]);

    if (pass_single_b_signal_selection){
        h_Deta_mu4_w_sig_sel[nsign]->Fill(deta[nsign],weight[nsign]);
        h_Deta_zoomin_mu4_w_sig_sel[nsign]->Fill(deta[nsign],weight[nsign]);
        h_Dphi_mu4_w_sig_sel[nsign]->Fill(dphi[nsign],weight[nsign]);
        h_Dphi_zoomin_mu4_w_sig_sel[nsign]->Fill(dphi[nsign],weight[nsign]);
        h_DR_mu4_w_sig_sel[nsign]->Fill(dr[nsign],weight[nsign]);
        h_DR_zoomin_mu4_w_sig_sel[nsign]->Fill(dr[nsign],weight[nsign]);
        if (mu1PassSingle[nsign]){
            h_pt2nd_mu4_w_sig_sel[nsign]->Fill(m2pt[nsign],weight[nsign]);
            h_pt2nd_vs_q_eta_2nd_mu4_w_sig_sel[nsign]->Fill(m2charge[nsign] * m2eta[nsign], m2pt[nsign],weight[nsign]);
        }
        if (mu2PassSingle[nsign]){
            h_pt2nd_mu4_w_sig_sel[nsign]->Fill(m1pt[nsign],weight[nsign]);
            h_pt2nd_vs_q_eta_2nd_mu4_w_sig_sel[nsign]->Fill(m1charge[nsign] * m1eta[nsign], m1pt[nsign],weight[nsign]);
        }

        h_pair_eta_vs_pair_pT_mu4_w_sig_sel[nsign]->Fill(pair_pt[nsign],pair_eta[nsign],weight[nsign]);
        h_Deta_Dphi_mu4_w_sig_sel[nsign]->Fill(dphi[nsign],deta[nsign],weight[nsign]);
        h_eta1_eta2_mu4_w_sig_sel[nsign]->Fill(m2eta[nsign],m1eta[nsign],weight[nsign]);
        h_eta_avg_Deta_mu4_w_sig_sel[nsign]->Fill(deta[nsign],etaavg[nsign],weight[nsign]);
        h_eta_avg_Dphi_mu4_w_sig_sel[nsign]->Fill(dphi[nsign],etaavg[nsign],weight[nsign]);
        h_minv_zoomin_mu4_w_sig_sel[nsign]->Fill(minv[nsign],weight[nsign]);
        h_pair_pt_log_mu4_w_sig_sel[nsign]->Fill(pair_pt[nsign],weight[nsign]);
        h_minv_pair_pt_log_mu4_w_sig_sel[nsign]->Fill(pair_pt[nsign],minv[nsign],weight[nsign]);
    }

    if (passmu4mu4noL1[nsign]){
        h_Deta_mu4_mu4noL1_given_mu4[nsign]->Fill(deta[nsign],weight[nsign]);
        h_Deta_zoomin_mu4_mu4noL1_given_mu4[nsign]->Fill(deta[nsign],weight[nsign]);
        h_Dphi_mu4_mu4noL1_given_mu4[nsign]->Fill(dphi[nsign],weight[nsign]);
        h_Dphi_zoomin_mu4_mu4noL1_given_mu4[nsign]->Fill(dphi[nsign],weight[nsign]);
        h_DR_mu4_mu4noL1_given_mu4[nsign]->Fill(dr[nsign],weight[nsign]);
        h_DR_zoomin_mu4_mu4noL1_given_mu4[nsign]->Fill(dr[nsign],weight[nsign]);
        if (mu1PassSingle[nsign]){
            h_pt2nd_mu4_mu4noL1_given_mu4[nsign]->Fill(m2pt[nsign],weight[nsign]);
            h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_given_mu4[nsign]->Fill(m2charge[nsign] * m2eta[nsign], m2pt[nsign],weight[nsign]);
        }
        if (mu2PassSingle[nsign]){
            h_pt2nd_mu4_mu4noL1_given_mu4[nsign]->Fill(m1pt[nsign],weight[nsign]);
            h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_given_mu4[nsign]->Fill(m1charge[nsign] * m1eta[nsign], m1pt[nsign],weight[nsign]);
        }

        h_pair_eta_vs_pair_pT_mu4_mu4noL1_given_mu4[nsign]->Fill(pair_pt[nsign],pair_eta[nsign],weight[nsign]);
        h_Deta_Dphi_mu4_mu4noL1_given_mu4[nsign]->Fill(dphi[nsign],deta[nsign],weight[nsign]);
        h_eta1_eta2_mu4_mu4noL1_given_mu4[nsign]->Fill(m2eta[nsign],m1eta[nsign],weight[nsign]);
        h_eta_avg_Deta_mu4_mu4noL1_given_mu4[nsign]->Fill(deta[nsign],etaavg[nsign],weight[nsign]);
        h_eta_avg_Dphi_mu4_mu4noL1_given_mu4[nsign]->Fill(dphi[nsign],etaavg[nsign],weight[nsign]);
        h_minv_zoomin_mu4_mu4noL1_given_mu4[nsign]->Fill(minv[nsign],weight[nsign]);
        h_pair_pt_log_mu4_mu4noL1_given_mu4[nsign]->Fill(pair_pt[nsign],weight[nsign]);
        h_minv_pair_pt_log_mu4_mu4noL1_given_mu4[nsign]->Fill(pair_pt[nsign],minv[nsign],weight[nsign]);

        if (pass_single_b_signal_selection){
            h_Deta_mu4_mu4noL1_given_mu4_w_sig_sel[nsign]->Fill(deta[nsign],weight[nsign]);
            h_Deta_zoomin_mu4_mu4noL1_given_mu4_w_sig_sel[nsign]->Fill(deta[nsign],weight[nsign]);
            h_Dphi_mu4_mu4noL1_given_mu4_w_sig_sel[nsign]->Fill(dphi[nsign],weight[nsign]);
            h_Dphi_zoomin_mu4_mu4noL1_given_mu4_w_sig_sel[nsign]->Fill(dphi[nsign],weight[nsign]);
            h_DR_mu4_mu4noL1_given_mu4_w_sig_sel[nsign]->Fill(dr[nsign],weight[nsign]);
            h_DR_zoomin_mu4_mu4noL1_given_mu4_w_sig_sel[nsign]->Fill(dr[nsign],weight[nsign]);
            if (mu1PassSingle[nsign]){
                h_pt2nd_mu4_mu4noL1_given_mu4_w_sig_sel[nsign]->Fill(m2pt[nsign],weight[nsign]);
                h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_given_mu4_w_sig_sel[nsign]->Fill(m2charge[nsign] * m2eta[nsign], m2pt[nsign],weight[nsign]);
            }
            if (mu2PassSingle[nsign]){
                h_pt2nd_mu4_mu4noL1_given_mu4_w_sig_sel[nsign]->Fill(m1pt[nsign],weight[nsign]);
                h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_given_mu4_w_sig_sel[nsign]->Fill(m1charge[nsign] * m1eta[nsign], m1pt[nsign],weight[nsign]);
            }

            h_pair_eta_vs_pair_pT_mu4_mu4noL1_given_mu4_w_sig_sel[nsign]->Fill(pair_pt[nsign],pair_eta[nsign],weight[nsign]);
            h_Deta_Dphi_mu4_mu4noL1_given_mu4_w_sig_sel[nsign]->Fill(dphi[nsign],deta[nsign],weight[nsign]);
            h_eta1_eta2_mu4_mu4noL1_given_mu4_w_sig_sel[nsign]->Fill(m2eta[nsign],m1eta[nsign],weight[nsign]);
            h_eta_avg_Deta_mu4_mu4noL1_given_mu4_w_sig_sel[nsign]->Fill(deta[nsign],etaavg[nsign],weight[nsign]);
            h_eta_avg_Dphi_mu4_mu4noL1_given_mu4_w_sig_sel[nsign]->Fill(dphi[nsign],etaavg[nsign],weight[nsign]);
            h_minv_zoomin_mu4_mu4noL1_given_mu4_w_sig_sel[nsign]->Fill(minv[nsign],weight[nsign]);
            h_pair_pt_log_mu4_mu4noL1_given_mu4_w_sig_sel[nsign]->Fill(pair_pt[nsign],weight[nsign]);
            h_minv_pair_pt_log_mu4_mu4noL1_given_mu4_w_sig_sel[nsign]->Fill(pair_pt[nsign],minv[nsign],weight[nsign]);

        }
    }
    if (pass2mu4[nsign]){
        h_Deta_2mu4_given_mu4[nsign]->Fill(deta[nsign],weight[nsign]);
        h_Deta_zoomin_2mu4_given_mu4[nsign]->Fill(deta[nsign],weight[nsign]);
        h_Dphi_2mu4_given_mu4[nsign]->Fill(dphi[nsign],weight[nsign]);
        h_Dphi_zoomin_2mu4_given_mu4[nsign]->Fill(dphi[nsign],weight[nsign]);
        h_DR_2mu4_given_mu4[nsign]->Fill(dr[nsign],weight[nsign]);
        h_DR_zoomin_2mu4_given_mu4[nsign]->Fill(dr[nsign],weight[nsign]);
        if (mu1PassSingle[nsign]){
            h_pt2nd_2mu4_given_mu4[nsign]->Fill(m2pt[nsign],weight[nsign]);
            h_pt2nd_vs_q_eta_2nd_2mu4_given_mu4[nsign]->Fill(m2charge[nsign] * m2eta[nsign], m2pt[nsign],weight[nsign]);
        }
        if (mu2PassSingle[nsign]){
            h_pt2nd_2mu4_given_mu4[nsign]->Fill(m1pt[nsign],weight[nsign]);
            h_pt2nd_vs_q_eta_2nd_2mu4_given_mu4[nsign]->Fill(m1charge[nsign] * m1eta[nsign], m1pt[nsign],weight[nsign]);
        }

        h_pair_eta_vs_pair_pT_2mu4_given_mu4[nsign]->Fill(pair_pt[nsign],pair_eta[nsign],weight[nsign]);
        h_Deta_Dphi_2mu4_given_mu4[nsign]->Fill(dphi[nsign],deta[nsign],weight[nsign]);
        h_eta1_eta2_2mu4_given_mu4[nsign]->Fill(m2eta[nsign],m1eta[nsign],weight[nsign]);
        h_eta_avg_Deta_2mu4_given_mu4[nsign]->Fill(deta[nsign],etaavg[nsign],weight[nsign]);
        h_eta_avg_Dphi_2mu4_given_mu4[nsign]->Fill(dphi[nsign],etaavg[nsign],weight[nsign]);
        h_minv_zoomin_2mu4_given_mu4[nsign]->Fill(minv[nsign],weight[nsign]);
        h_pair_pt_log_2mu4_given_mu4[nsign]->Fill(pair_pt[nsign],weight[nsign]);
        h_minv_pair_pt_log_2mu4_given_mu4[nsign]->Fill(pair_pt[nsign],minv[nsign],weight[nsign]);

        if (pass_single_b_signal_selection){
            h_Deta_2mu4_given_mu4_w_sig_sel[nsign]->Fill(deta[nsign],weight[nsign]);
            h_Deta_zoomin_2mu4_given_mu4_w_sig_sel[nsign]->Fill(deta[nsign],weight[nsign]);
            h_Dphi_2mu4_given_mu4_w_sig_sel[nsign]->Fill(dphi[nsign],weight[nsign]);
            h_Dphi_zoomin_2mu4_given_mu4_w_sig_sel[nsign]->Fill(dphi[nsign],weight[nsign]);
            h_DR_2mu4_given_mu4_w_sig_sel[nsign]->Fill(dr[nsign],weight[nsign]);
            h_DR_zoomin_2mu4_given_mu4_w_sig_sel[nsign]->Fill(dr[nsign],weight[nsign]);
            if (mu1PassSingle[nsign]){
                h_pt2nd_2mu4_given_mu4_w_sig_sel[nsign]->Fill(m2pt[nsign],weight[nsign]);
                h_pt2nd_vs_q_eta_2nd_2mu4_given_mu4_w_sig_sel[nsign]->Fill(m2charge[nsign] * m2eta[nsign], m2pt[nsign],weight[nsign]);
            }
            if (mu2PassSingle[nsign]){
                h_pt2nd_2mu4_given_mu4_w_sig_sel[nsign]->Fill(m1pt[nsign],weight[nsign]);
                h_pt2nd_vs_q_eta_2nd_2mu4_given_mu4_w_sig_sel[nsign]->Fill(m1charge[nsign] * m1eta[nsign], m1pt[nsign],weight[nsign]);
            }

            h_pair_eta_vs_pair_pT_2mu4_given_mu4_w_sig_sel[nsign]->Fill(pair_pt[nsign],pair_eta[nsign],weight[nsign]);
            h_Deta_Dphi_2mu4_given_mu4_w_sig_sel[nsign]->Fill(dphi[nsign],deta[nsign],weight[nsign]);
            h_eta1_eta2_2mu4_given_mu4_w_sig_sel[nsign]->Fill(m2eta[nsign],m1eta[nsign],weight[nsign]);
            h_eta_avg_Deta_2mu4_given_mu4_w_sig_sel[nsign]->Fill(deta[nsign],etaavg[nsign],weight[nsign]);
            h_eta_avg_Dphi_2mu4_given_mu4_w_sig_sel[nsign]->Fill(dphi[nsign],etaavg[nsign],weight[nsign]);
            h_minv_zoomin_2mu4_given_mu4_w_sig_sel[nsign]->Fill(minv[nsign],weight[nsign]);
            h_pair_pt_log_2mu4_given_mu4_w_sig_sel[nsign]->Fill(pair_pt[nsign],weight[nsign]);
            h_minv_pair_pt_log_2mu4_given_mu4_w_sig_sel[nsign]->Fill(pair_pt[nsign],minv[nsign],weight[nsign]);
        }
    }

}


void MuonPairPlottingPP::CalculateSingleTrigEffcyRatio(TH1* h1, TH1* h2, TH1* h3) {
    if (!h1 || !h2 || !h3) {
        std::cerr << "Null pointer passed to CalculateSingleTrigEffcyRatio." << std::endl;
        return;
    }

    // Clone h1 and add h2
    TH1* hclone = (TH1*)h1->Clone(Form("%s_clone", h1->GetName()));
    hclone->Add(h2);

    // Divide h3 by hclone in place
    h3->Divide(hclone);

    delete hclone; // Clean up
}


void MuonPairPlottingPP::CalculateTrigEffcyRatio(){
    for (int isign = 0; isign < ParamsSet::nSigns; isign++){
        h_Deta_mu4_mu4noL1_given_mu4[isign]->Divide(h_Deta_mu4[isign]);
        h_Deta_zoomin_mu4_mu4noL1_given_mu4[isign]->Divide(h_Deta_zoomin_mu4[isign]);
        h_Dphi_mu4_mu4noL1_given_mu4[isign]->Divide(h_Dphi_mu4[isign]);
        h_Dphi_zoomin_mu4_mu4noL1_given_mu4[isign]->Divide(h_Dphi_zoomin_mu4[isign]);
        h_DR_mu4_mu4noL1_given_mu4[isign]->Divide(h_DR_mu4[isign]);
        h_DR_zoomin_mu4_mu4noL1_given_mu4[isign]->Divide(h_DR_zoomin_mu4[isign]);
        h_pt2nd_mu4_mu4noL1_given_mu4[isign]->Divide(h_pt2nd_mu4[isign]);
        h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_given_mu4[isign]->Divide(h_pt2nd_vs_q_eta_2nd_mu4[isign]);
        h_pair_eta_vs_pair_pT_mu4_mu4noL1_given_mu4[isign]->Divide(h_pair_eta_vs_pair_pT_mu4[isign]);
        h_Deta_Dphi_mu4_mu4noL1_given_mu4[isign]->Divide(h_Deta_Dphi_mu4[isign]);
        h_eta1_eta2_mu4_mu4noL1_given_mu4[isign]->Divide(h_eta1_eta2_mu4[isign]);
        h_eta_avg_Deta_mu4_mu4noL1_given_mu4[isign]->Divide(h_eta_avg_Deta_mu4[isign]);
        h_eta_avg_Dphi_mu4_mu4noL1_given_mu4[isign]->Divide(h_eta_avg_Dphi_mu4[isign]);
        h_minv_zoomin_mu4_mu4noL1_given_mu4[isign]->Divide(h_minv_zoomin_mu4[isign]);
        h_pair_pt_log_mu4_mu4noL1_given_mu4[isign]->Divide(h_pair_pt_log_mu4[isign]);
        h_minv_pair_pt_log_mu4_mu4noL1_given_mu4[isign]->Divide(h_minv_pair_pt_log_mu4[isign]);

        h_Deta_2mu4_given_mu4[isign]->Divide(h_Deta_mu4[isign]);
        h_Deta_zoomin_2mu4_given_mu4[isign]->Divide(h_Deta_zoomin_mu4[isign]);
        h_Dphi_2mu4_given_mu4[isign]->Divide(h_Dphi_mu4[isign]);
        h_Dphi_zoomin_2mu4_given_mu4[isign]->Divide(h_Dphi_zoomin_mu4[isign]);
        h_DR_2mu4_given_mu4[isign]->Divide(h_DR_mu4[isign]);
        h_DR_zoomin_2mu4_given_mu4[isign]->Divide(h_DR_zoomin_mu4[isign]);
        h_pt2nd_2mu4_given_mu4[isign]->Divide(h_pt2nd_mu4[isign]);
        h_pt2nd_vs_q_eta_2nd_2mu4_given_mu4[isign]->Divide(h_pt2nd_vs_q_eta_2nd_mu4[isign]);
        h_pair_eta_vs_pair_pT_2mu4_given_mu4[isign]->Divide(h_pair_eta_vs_pair_pT_mu4[isign]);
        h_Deta_Dphi_2mu4_given_mu4[isign]->Divide(h_Deta_Dphi_mu4[isign]);
        h_eta1_eta2_2mu4_given_mu4[isign]->Divide(h_eta1_eta2_mu4[isign]);
        h_eta_avg_Deta_2mu4_given_mu4[isign]->Divide(h_eta_avg_Deta_mu4[isign]);
        h_eta_avg_Dphi_2mu4_given_mu4[isign]->Divide(h_eta_avg_Dphi_mu4[isign]);
        h_minv_zoomin_2mu4_given_mu4[isign]->Divide(h_minv_zoomin_mu4[isign]);
        h_pair_pt_log_2mu4_given_mu4[isign]->Divide(h_pair_pt_log_mu4[isign]);
        h_minv_pair_pt_log_2mu4_given_mu4[isign]->Divide(h_minv_pair_pt_log_mu4[isign]);

        h_Deta_mu4_mu4noL1_given_mu4_w_sig_sel[isign]->Divide(h_Deta_mu4_w_sig_sel[isign]);
        h_Deta_zoomin_mu4_mu4noL1_given_mu4_w_sig_sel[isign]->Divide(h_Deta_zoomin_mu4_w_sig_sel[isign]);
        h_Dphi_mu4_mu4noL1_given_mu4_w_sig_sel[isign]->Divide(h_Dphi_mu4_w_sig_sel[isign]);
        h_Dphi_zoomin_mu4_mu4noL1_given_mu4_w_sig_sel[isign]->Divide(h_Dphi_zoomin_mu4_w_sig_sel[isign]);
        h_DR_mu4_mu4noL1_given_mu4_w_sig_sel[isign]->Divide(h_DR_mu4_w_sig_sel[isign]);
        h_DR_zoomin_mu4_mu4noL1_given_mu4_w_sig_sel[isign]->Divide(h_DR_zoomin_mu4_w_sig_sel[isign]);
        h_pt2nd_mu4_mu4noL1_given_mu4_w_sig_sel[isign]->Divide(h_pt2nd_mu4_w_sig_sel[isign]);
        h_pt2nd_vs_q_eta_2nd_mu4_mu4noL1_given_mu4_w_sig_sel[isign]->Divide(h_pt2nd_vs_q_eta_2nd_mu4_w_sig_sel[isign]);
        h_pair_eta_vs_pair_pT_mu4_mu4noL1_given_mu4_w_sig_sel[isign]->Divide(h_pair_eta_vs_pair_pT_mu4_w_sig_sel[isign]);
        h_Deta_Dphi_mu4_mu4noL1_given_mu4_w_sig_sel[isign]->Divide(h_Deta_Dphi_mu4_w_sig_sel[isign]);
        h_eta1_eta2_mu4_mu4noL1_given_mu4_w_sig_sel[isign]->Divide(h_eta1_eta2_mu4_w_sig_sel[isign]);
        h_eta_avg_Deta_mu4_mu4noL1_given_mu4_w_sig_sel[isign]->Divide(h_eta_avg_Deta_mu4_w_sig_sel[isign]);
        h_eta_avg_Dphi_mu4_mu4noL1_given_mu4_w_sig_sel[isign]->Divide(h_eta_avg_Dphi_mu4_w_sig_sel[isign]);
        h_minv_zoomin_mu4_mu4noL1_given_mu4_w_sig_sel[isign]->Divide(h_minv_zoomin_mu4_w_sig_sel[isign]);
        h_pair_pt_log_mu4_mu4noL1_given_mu4_w_sig_sel[isign]->Divide(h_pair_pt_log_mu4_w_sig_sel[isign]);
        h_minv_pair_pt_log_mu4_mu4noL1_given_mu4_w_sig_sel[isign]->Divide(h_minv_pair_pt_log_mu4_w_sig_sel[isign]);

        h_Deta_2mu4_given_mu4_w_sig_sel[isign]->Divide(h_Deta_mu4_w_sig_sel[isign]);
        h_Deta_zoomin_2mu4_given_mu4_w_sig_sel[isign]->Divide(h_Deta_zoomin_mu4_w_sig_sel[isign]);
        h_Dphi_2mu4_given_mu4_w_sig_sel[isign]->Divide(h_Dphi_mu4_w_sig_sel[isign]);
        h_Dphi_zoomin_2mu4_given_mu4_w_sig_sel[isign]->Divide(h_Dphi_zoomin_mu4_w_sig_sel[isign]);
        h_DR_2mu4_given_mu4_w_sig_sel[isign]->Divide(h_DR_mu4_w_sig_sel[isign]);
        h_DR_zoomin_2mu4_given_mu4_w_sig_sel[isign]->Divide(h_DR_zoomin_mu4_w_sig_sel[isign]);
        h_pt2nd_2mu4_given_mu4_w_sig_sel[isign]->Divide(h_pt2nd_mu4_w_sig_sel[isign]);
        h_pt2nd_vs_q_eta_2nd_2mu4_given_mu4_w_sig_sel[isign]->Divide(h_pt2nd_vs_q_eta_2nd_mu4_w_sig_sel[isign]);
        h_pair_eta_vs_pair_pT_2mu4_given_mu4_w_sig_sel[isign]->Divide(h_pair_eta_vs_pair_pT_mu4_w_sig_sel[isign]);
        h_Deta_Dphi_2mu4_given_mu4_w_sig_sel[isign]->Divide(h_Deta_Dphi_mu4_w_sig_sel[isign]);
        h_eta1_eta2_2mu4_given_mu4_w_sig_sel[isign]->Divide(h_eta1_eta2_mu4_w_sig_sel[isign]);
        h_eta_avg_Deta_2mu4_given_mu4_w_sig_sel[isign]->Divide(h_eta_avg_Deta_mu4_w_sig_sel[isign]);
        h_eta_avg_Dphi_2mu4_given_mu4_w_sig_sel[isign]->Divide(h_eta_avg_Dphi_mu4_w_sig_sel[isign]);
        h_minv_zoomin_2mu4_given_mu4_w_sig_sel[isign]->Divide(h_minv_zoomin_mu4_w_sig_sel[isign]);
        h_pair_pt_log_2mu4_given_mu4_w_sig_sel[isign]->Divide(h_pair_pt_log_mu4_w_sig_sel[isign]);
        h_minv_pair_pt_log_2mu4_given_mu4_w_sig_sel[isign]->Divide(h_minv_pair_pt_log_mu4_w_sig_sel[isign]);
    }
}

void MuonPairPlottingPP::WriteOutput(){

    if (mode == 1){
        outFile->cd();
        outFile->Write();
        outFile->Close();

    }else{ // mode = 3

    }
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

    output_non_trig_effcy_hists = !(trigger_mode == 0 || trigger_mode == 1);

  	InitInput();
    InitOutput();
  	InitHists();
  	ProcessData();
    CalculateTrigEffcyRatio();
  	WriteOutput();

    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;;
    std::cout << "CPU time used for running mode " << mode << " is " << cpu_time_used << " seconds" << std::endl;

}