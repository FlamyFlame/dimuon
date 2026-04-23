#include"TrigEffPlotterPbPb.cxx"

void trig_effcy_plot_pbpb(){
    std::vector<std::string> var1Ds = {"Deta", "Deta_zoomin", "Dphi", "Dphi_zoomin", "DR", "DR_zoomin", "DR_0_2", "minv_zoomin", "pair_pt_log"};
    std::vector<std::string> var2Ds = {
                                        //"DR_zoomin_vs_pt2nd", "DR_0_2_vs_pt2nd", "pair_eta_vs_pair_pT", "Deta_Dphi", "eta1_eta2", "eta_avg_Deta", "eta_avg_Dphi", "minv_pair_pt_log",
                                        // "Deta_vs_pT_1st", "Deta_zoomin_vs_pT_1st", "Dphi_vs_pT_1st", "Dphi_zoomin_vs_pT_1st", "DR_vs_pT_1st", "DR_zoomin_vs_pT_1st", "minv_zoomin_vs_pT_1st", "pair_pt_log_vs_pT_1st", "pt2nd_vs_pT_1st", 
                                        // "Deta_vs_pair_pT", "Deta_zoomin_vs_pair_pT", "Dphi_vs_pair_pT", "Dphi_zoomin_vs_pair_pT", "DR_vs_pair_pT", "DR_zoomin_vs_pair_pT", "minv_zoomin_vs_pair_pT", "pair_pt_log_vs_pair_pT", "pt2nd_vs_pair_pT",
                                      };
    std::vector<std::string> var2DsProf = {
       "DR_zoomin_vs_pt2nd", "DR_zoomin_vs_pair_pt_log",
       // "Deta_vs_pT_1st", "Deta_zoomin_vs_pT_1st", "Dphi_vs_pT_1st", "Dphi_zoomin_vs_pT_1st", "DR_vs_pT_1st", "DR_zoomin_vs_pT_1st", "minv_zoomin_vs_pT_1st", "pair_pt_log_vs_pT_1st", "pt2nd_vs_pT_1st", 
       // "Deta_vs_pair_pT", "Deta_zoomin_vs_pair_pT", "Dphi_vs_pair_pT", "Dphi_zoomin_vs_pair_pT", "DR_vs_pair_pT", "minv_zoomin_vs_pair_pT", "pair_pt_log_vs_pair_pT", "pt2nd_vs_pair_pT",
    };

    std::map<std::string,std::pair<bool,bool>> logaxes = {
       {"pt2nd",{true,false}},
       {"pair_pt_log",{true,false}},
       // {"pt2nd_vs_q_eta_2nd",{false,true}},
       {"pt2nd_vs_q_eta2nd",{false,true}},
       {"pair_eta_vs_pair_pT",{true,false}},
       {"minv_pair_pt_log",{true,true}},
       {"Deta_vs_pT_1st",{true,false}},
       {"Deta_zoomin_vs_pT_1st",{true,false}},
       {"Dphi_vs_pT_1st",{true,false}},
       {"Dphi_zoomin_vs_pT_1st",{true,false}},
       {"DR_vs_pT_1st",{true,false}},
       {"DR_zoomin_vs_pT_1st",{true,false}},
       {"minv_zoomin_vs_pT_1st",{true,false}},
       {"pair_pt_log_vs_pT_1st",{true,true}},
       {"pt2nd_vs_pT_1st",{true,true}},
       {"Deta_vs_pair_pT",{true,false}},
       {"Deta_zoomin_vs_pair_pT",{true,false}},
       {"Dphi_vs_pair_pT",{true,false}},
       {"Dphi_zoomin_vs_pair_pT",{true,false}},
       {"DR_vs_pair_pT",{true,false}},
       {"DR_zoomin_vs_pair_pT",{true,false}},
       {"DR_zoomin_vs_pt2nd",{true,false}},
       {"minv_zoomin_vs_pair_pT",{true,false}},
       {"pair_pt_log_vs_pair_pT",{true,true}},
       {"pt2nd_vs_pair_p",{true,true}},
    };

    std::map<std::string,std::pair<bool,bool>> var2DProj = { // pair observables ONLY
    };

    std::map<std::string,std::pair<bool,bool>> var2DSingleMuonEffcyProj = {
        // {"pt2nd_vs_q_eta_2nd", {true,  true}},   // PX (→ q_eta_2nd) and PY (→ pt2nd)
        {"pt2nd_vs_q_eta2nd", {true,  true}},   // PX (→ q_eta2nd) and PY (→ pt2nd)
        {"pt2nd_vs_phi2nd",    {true,  false}}  // PX only (→ phi2nd)
    };

    // TrigEffPlotterPbPb::Rect rSigned  = {0.60,0.15,0.88,0.35};
    // TrigEffPlotterPbPb::Rect rSignal  = {0.20,0.70,0.50,0.88};
    // TrigEffPlotterPbPb::Rect rOpSig   = {0.60,0.15,0.88,0.31};

    // std::map<std::string,TrigEffPlotterPbPb::Rect> legSigned  = { {"Dphi", rSigned } };
    // std::map<std::string,TrigEffPlotterPbPb::Rect> legSignal  = { {"Dphi", rSignal } };

    std::map<std::string,TrigEffPlotterPbPb::Rect> legSigned  = {};
    std::map<std::string,TrigEffPlotterPbPb::Rect> legSignal  = {};
    std::map<std::string,TrigEffPlotterPbPb::Rect> legOpSig   = {};

    std::vector<std::string> filters = {};
    // filters = {"_sepr"};
    // std::map<std::string,TrigEffPlotterPbPb::Rect> legOpSig   = {{"DR_zoomin", {0.5,0.15,0.87,0.33}  }};
    // filters = {"_good_accept"};
    // filters = {"_inv_w_by_single_mu_effcy"};

    // only draw single-b for weighted trigger efficiency, showing dR corrections & its influence on other observables
    // never for single-muon efficiencies (single-b gives a biased sample)
    bool draw_single_b = (std::find(filters.begin(), filters.end(), "_inv_w_by_single_mu_effcy") != filters.end());
    TrigEffPlotterPbPb plotter(23, var1Ds, "include_upc", false, true,
                           filters, logaxes,
                           var2Ds, var2DsProf,
                           legSigned, legSignal, legOpSig,
                           var2DProj, var2DSingleMuonEffcyProj);

    plotter.only_draw_single_muon_effcy = true;
    plotter.Run();
}
