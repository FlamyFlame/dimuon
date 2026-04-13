#include <TSystem.h>
#include <algorithm>
#include "RDFBasedHistFillingData.cxx"

void RDFBasedHistFillingPP::SetIOPathsHook(){
    infile_var1D_json = "var1D_pp.json";

    if (run_year == 24){
        std::vector<std::string> input_candidates = {
            "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/muon_pairs_pp_2024" + base_trig_suffix + input_mindR_suffix + "_res_cut_v2.root",
            "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/muon_pairs_pp_2024" + base_trig_suffix + input_mindR_suffix + ".root",
            "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/muon_pairs_pp_2024" + base_trig_suffix + "_no_res_cut.root",
            "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/muon_pairs_pp_2024" + base_trig_suffix + ".root"
        };

        std::string in_path;
        for (const auto& cand : input_candidates) {
            if (!gSystem->AccessPathName(cand.c_str())) {
                in_path = cand;
                break;
            }
        }
        if (in_path.empty()) {
            throw std::runtime_error("RDFBasedHistFillingPP: input file not found for run_year=24 and base_trig_suffix=" + base_trig_suffix);
        }

        input_files.push_back(in_path);
        output_file = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/histograms_real_pairs_pp_2024" + out_file_suffix + ".root";
    } else if (run_year == 17){
        std::vector<std::string> input_candidates = {
            "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_run2/muon_pairs_pp_2017" + base_trig_suffix + input_mindR_suffix + "_res_cut_v2.root",
            "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_run2/muon_pairs_pp_2017" + base_trig_suffix + input_mindR_suffix + ".root",
            "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_run2/muon_pairs_pp_2017" + base_trig_suffix + "_no_res_cut.root",
            "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_run2/muon_pairs_pp_2017" + base_trig_suffix + ".root"
        };

        std::string in_path;
        for (const auto& cand : input_candidates) {
            if (!gSystem->AccessPathName(cand.c_str())) {
                in_path = cand;
                break;
            }
        }
        if (in_path.empty()) {
            throw std::runtime_error("RDFBasedHistFillingPP: input file not found for run_year=17 and base_trig_suffix=" + base_trig_suffix);
        }

        input_files.push_back(in_path);
        output_file = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_run2/histograms_real_pairs_pp_2017" + out_file_suffix + ".root";
    }
}

void RDFBasedHistFillingPP::InitializePPExtra(){
	if (run_year == 17){
        trigs = {"_mu4", "_2mu4"};
        trigs_pair = {{"_2mu4","_mu4"}};
    }

    // Set pp_crossx_lumi_factor from PPBaseClass map.
    // Only trigger_mode 2 (mu4_mu4noL1) and 3 (2mu4) are crossx modes;
    // trigger_mode 1 never calls FillHistogramsCrossx(), so leave factor at -1.
    {
        std::string trig_for_crossx;
        if      (trigger_mode == 2) trig_for_crossx = "mu4_mu4noL1";
        else if (trigger_mode == 3) trig_for_crossx = "2mu4";
        if (!trig_for_crossx.empty()) {
            const auto& m = PPBaseClass::CrossxFactorMap();
            auto it = m.find({run_year, trig_for_crossx});
            if (it != m.end()) pp_crossx_lumi_factor = it->second;
            // else: leave -1; FillHistogramsCrossx() will throw if called
        }
    }

    levels_trg_effcy_to_be_summed_w_musign_summing = {0,1,2}; // ss/op + mu1/2 + mu+/-

    levels_trg_effcy_filters_1D_pre_sum = {{"_ss", "_op"},
                                        {"_mu1passmu4", "_mu2passmu4"},
                                        {"_sign1", "_sign2"},
                                        trigs,
                                        trg_effcy_biases};

    levels_trg_effcy_filters_2D_3D_pre_sum = {{"_ss", "_op"},
                                        {"_mu1passmu4", "_mu2passmu4"},
                                        {"_sign1", "_sign2"},
                                        trigs,
                                        {"", "_sepr"}};

    categories_essential = pair_signs;
}

void RDFBasedHistFillingPP::FlattenTrigEffcyFiltersExtra(){

    // build to-be-summed levels, with mu-sign summing
    for (int level_ind = 0;
         level_ind < static_cast<int>(levels_trg_effcy_filters_2D_3D_pre_sum.size());
         ++level_ind)
    {
        if (std::find(levels_trg_effcy_to_be_summed_w_musign_summing.begin(),
                      levels_trg_effcy_to_be_summed_w_musign_summing.end(),
                      level_ind) != levels_trg_effcy_to_be_summed_w_musign_summing.end())
        {
            levels_trg_effcy_filters_to_be_summed_w_musign_summing.push_back(
                levels_trg_effcy_filters_2D_3D_pre_sum.at(level_ind)
            );
        }
    }

    // flatten to-be-summed levels, with mu-sign summing
    HistFillUtils::flatten_levels(levels_trg_effcy_filters_to_be_summed_w_musign_summing, trg_effcy_filters_to_be_summed_w_musign_summing);    

    // build post-sum levels, with mu-sign summing
    HistFillUtils::write_post_sum_levels(levels_trg_effcy_filters_1D_pre_sum,
                          levels_trg_effcy_to_be_summed_w_musign_summing,
                          levels_trg_effcy_filters_1D_post_sum_w_musign_summing);

    HistFillUtils::write_post_sum_levels(levels_trg_effcy_filters_2D_3D_pre_sum,
                          levels_trg_effcy_to_be_summed_w_musign_summing,
                          levels_trg_effcy_filters_2D_3D_post_sum_w_musign_summing);

    // flatten post-sum levels, with mu-sign summing
    HistFillUtils::flatten_levels(levels_trg_effcy_filters_1D_post_sum_w_musign_summing, trg_effcy_filters_1D_post_sum_w_musign_summing);
    HistFillUtils::flatten_levels(levels_trg_effcy_filters_2D_3D_post_sum_w_musign_summing, trg_effcy_filters_2D_3D_post_sum_w_musign_summing);
}

void RDFBasedHistFillingPP::FillHistogramsSingleMuonEffcy(){
    if (trigger_mode == 1){ // trigger efficiencies: mu4_mu4noL1 | mu4 & 2mu4 | mu4
        FillHistogramsDimuTrigGivenMu4();
    } else if (trigger_mode == 0){
        FillHistogramsMu4GivenMB();
    }
}

void RDFBasedHistFillingPP::FillHistogramsDimuTrigGivenMu4(){
    try{
        for (std::string pair_sign : pair_signs){
            std::string df_name = "df" + pair_sign;
            ROOT::RDF::RNode& node = map_at_checked(df_map, df_name, Form("FillHistogramsSingleMuonEffcy: df_map.at(%s)", df_name.c_str()));

            for (auto mu4sel : {"_mu1passmu4", "_mu2passmu4"}){ // mu4 selection

                std::string df_name = "df" + pair_sign + mu4sel; // e.g, df_ss_mu1passmu4

                std::string ind1st = (std::string(mu4sel) == "_mu1passmu4")? "1" : "2";
                std::string ind2nd = (std::string(mu4sel) == "_mu1passmu4")? "2" : "1";

                df_map.emplace(df_name, node.Filter("m" + ind1st + ".passmu4"));
                df_map.at(df_name) = df_map.at(df_name).Define("pt1st", "m" + ind1st + ".pt");
                df_map.at(df_name) = df_map.at(df_name).Define("pt2nd", "m" + ind2nd + ".pt");
                df_map.at(df_name) = df_map.at(df_name).Define("charge2nd","m" + ind2nd + ".charge");
                df_map.at(df_name) = df_map.at(df_name).Define("eta2nd",    "m" + ind2nd + ".eta");
                df_map.at(df_name) = df_map.at(df_name).Define("phi2nd",    "m" + ind2nd + ".phi");
                df_map.at(df_name) = df_map.at(df_name).Define("q_eta2nd",  "charge2nd * eta2nd");
                df_map.at(df_name) = df_map.at(df_name).Define("mu2nd_passmu4noL1", "m" + ind2nd + ".passmu4noL1");
                df_map.at(df_name) = df_map.at(df_name).Define("mu2nd_good_acceptance",  "pt2nd >= 6 && ((eta2nd > 1.1 && eta2nd < 2.3) || (eta2nd > -2.3 && eta2nd < -1.2))");

                df_map.emplace(df_name + "_sign1", map_at_checked(df_map, df_name, Form("FillHistogramsSingleMuonEffcy: df_map.at(%s)", df_name.c_str())).Filter("charge2nd > 0")); // e.g, df_ss_mu1passmu4_sign1
                df_map.emplace(df_name + "_sign2", map_at_checked(df_map, df_name, Form("FillHistogramsSingleMuonEffcy: df_map.at(%s)", df_name.c_str())).Filter("charge2nd < 0"));

                for (auto mu_sign : {"_sign1", "_sign2"}){
                    std::string df_name = "df" + pair_sign + mu4sel + mu_sign;

                    for (auto trg : {"_mu4", "_mu4_mu4noL1", "_2mu4", "_2mu4_AND_mu4_mu4noL1"}){
                        std::string trg_filter = map_at_checked(trig_to_filter_str_map, trg, Form("trig_to_filter_str_map.at(%s)", trg));
                        if (trg_filter.empty()) df_map.emplace(df_name + trg, map_at_checked(df_map, df_name, Form("FillHistogramsDimuTrigGivenMu4: df_map.at(%s)", df_name.c_str())));
                        else                    df_map.emplace(df_name + trg, map_at_checked(df_map, df_name, Form("FillHistogramsDimuTrigGivenMu4: df_map.at(%s)", df_name.c_str())).Filter(trg_filter));

                        std::string df_name_new = "df" + pair_sign + mu4sel + mu_sign + trg;
                        df_map.emplace(df_name_new + "_sepr",        map_at_checked(df_map, df_name_new, Form("FillHistogramsDimuTrigGivenMu4: df_map.at(%s)", df_name_new.c_str())).Filter("passSeparated"));
                        df_map.emplace(df_name_new + "_good_accept", map_at_checked(df_map, df_name_new, Form("FillHistogramsDimuTrigGivenMu4: df_map.at(%s)", df_name_new.c_str())).Filter("mu2nd_good_acceptance"));

                        for (auto bias : trg_effcy_biases){ // additional selection / bias in data sample
                            std::string df_name_bias = "df" + pair_sign + mu4sel + mu_sign + trg + bias;
                            std::string filter = df_name_bias.substr(2);
                            FillHistogramsSingleDataFrame(filter, map_at_checked(df_map, df_name_bias, Form("FillHistogramsDimuTrigGivenMu4: df_map.at(%s)", df_name_bias.c_str())), false);
                        }
                    }
                }
            } // end loop mu4 selection
        } // end loop pair sign
    } catch(const std::out_of_range& e){
        std::cerr << "FillHistogramsSingleMuonEffcy:: out_of_range exception caught: " << e.what() << std::endl;
    } catch (const std::runtime_error& e) {
        std::cerr << "FillHistogramsSingleMuonEffcy:: RDF runtime error: " << e.what() << std::endl;
    }
}

void RDFBasedHistFillingPP::FillHistogramsMu4GivenMB(){}

void RDFBasedHistFillingPP::FillTrigEffcyHistsInvWeightedbySingleMuonEffcies(){}

void RDFBasedHistFillingPP::SumSingleMuonTrigEffHistsPP(){

    // 1D, with mu-sign summing
    HistFillUtils::SumTrigEffHistsGeneric<TH1D, std::string>(
        single_muon_trig_effcy_var1Ds,
        trg_effcy_filters_1D_post_sum_w_musign_summing,
        trg_effcy_filters_to_be_summed_w_musign_summing,
        hist1D_map,
        [](const std::string& var) {
            return "h_" + var;
        }
    );

    // 2D, with mu-sign summing
    HistFillUtils::SumTrigEffHistsGeneric<TH2D, std::array<std::string,2>>(
        single_muon_trig_effcy_var2Ds,
        trg_effcy_filters_2D_3D_post_sum_w_musign_summing,
        trg_effcy_filters_to_be_summed_w_musign_summing,
        hist2D_map,
        [](const std::array<std::string,2>& vars) {
            const std::string& varx = vars[0];
            const std::string& vary = vars[1];
            return "h_" + vary + "_vs_" + varx;
        }
    );

    // 3D, with mu-sign summing
    HistFillUtils::SumTrigEffHistsGeneric<TH3D, std::array<std::string,3>>(
        single_muon_trig_effcy_var3Ds,
        trg_effcy_filters_2D_3D_post_sum_w_musign_summing,
        trg_effcy_filters_to_be_summed_w_musign_summing,
        hist3D_map,
        [](const std::array<std::string,3>& vars) {
            const std::string& varx = vars[0];
            const std::string& vary = vars[1];
            const std::string& varz = vars[2];
            return "h_" + varz + "_vs_" + vary + "_vs_" + varx;
        }
    );
}

void RDFBasedHistFillingPP::CalculateSingleMuonTrigEffcyRatios(){
    CalculateSingleMuonTrigEffcyRatiosHelper(musigns);
}

void RDFBasedHistFillingPP::MakeAndWriteSingleMuonTrigEffPtGraphs(){
    if(useCoarseQEtaBin)    MakeAndWriteSingleMuonTrigEffPtGraphsHelper({});
    else                    MakeAndWriteSingleMuonTrigEffPtGraphsHelper(musigns);
}

void RDFBasedHistFillingPP::OpenEffcyPtFitFile(){}

void RDFBasedHistFillingPP::MakeAndWriteDRTrigEffGraphs(){}

// ============================================================================
// FillHistogramsCrossx: Single B → J/psi crossx measurements (trigger_mode==2)
// ============================================================================

void RDFBasedHistFillingPP::FillHistogramsCrossx(){
    if (pp_crossx_lumi_factor < 0.) {
        throw std::runtime_error(
            "[PP] FillHistogramsCrossx: pp_crossx_lumi_factor not set (value " +
            std::to_string(pp_crossx_lumi_factor) + "). "
            "No entry in PPBaseClass::CrossxFactorMap for (run_year=" +
            std::to_string(run_year) + ", trigger_mode=" + std::to_string(trigger_mode) + "). "
            "Valid crossx trigger_modes: 2 (mu4_mu4noL1), 3 (2mu4)."
        );
    }
    std::cout << "[PP] FillHistogramsCrossx: opposite-sign only, signal cuts, "
              << "crossx_weight = weight * " << pp_crossx_lumi_factor
              << " (1/L_int)" << std::endl;

    const std::string signal_cuts = "minv > 1.08 && minv < 2.9 && pair_pt > 8 && m1.charge * m1.eta < 2.2 && m2.charge * m2.eta < 2.2 && dr > 0.05";

    ROOT::RDF::RNode df_op_base = map_at_checked(df_map, "df_op", "FillHistogramsCrossx: df_op");
    ROOT::RDF::RNode df_single_b_crossx = df_op_base.Filter(signal_cuts);
    if (df_map.find("df_single_b_crossx") == df_map.end()) {
        df_map.emplace("df_single_b_crossx", df_single_b_crossx);
    }

    const double lumi_factor = pp_crossx_lumi_factor;
    ROOT::RDF::RNode df_single_b_crossx_weighted =
        df_single_b_crossx.Define("crossx_weight",
            [lumi_factor](double weight){ return weight * lumi_factor; },
            {"weight"});

    if (df_map.find("df_single_b_crossx_weighted") == df_map.end()) {
        df_map.emplace("df_single_b_crossx_weighted", df_single_b_crossx_weighted);
    }

    // Store RResultPtrs (lazy-evaluated) instead of immediate clones.
    // They will be evaluated during HistPostProcess() and converted to raw pointers.
    const int    npt   = (int)(pms.pT_bins_120.size() - 1);
    const double* ptbins = pms.pT_bins_120.data();

    // ROOT 6.34 TH3DModel has no mixed (variable/uniform) ctor; generate edge arrays for uniform axes.
    auto make_unif_edges = [](int n, double lo, double hi) {
        std::vector<double> e(n + 1);
        const double step = (hi - lo) / n;
        for (int i = 0; i <= n; ++i) e[i] = lo + i * step;
        return e;
    };
    const auto eta_edges  = make_unif_edges(44, -2.4,  2.4);
    const auto minv_edges = make_unif_edges(50,  1.0,  3.0);
    const auto dr_edges   = make_unif_edges(50,  0.05, 1.0);

    hist2d_rresultptr_map["h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts"] = df_single_b_crossx_weighted.Histo2D(
        ROOT::RDF::TH2DModel("h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts", ";p_{T}^{pair} [GeV];#eta^{pair}", npt, ptbins, 44, -2.4, 2.4),
        "pair_pt", "pair_eta", "crossx_weight");
    hist2d_rresultptr_map["h2d_crossx_pair_pt_minv_w_signal_cuts"] = df_single_b_crossx_weighted.Histo2D(
        ROOT::RDF::TH2DModel("h2d_crossx_pair_pt_minv_w_signal_cuts", ";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]", npt, ptbins, 50, 1.0, 3.0),
        "pair_pt", "minv", "crossx_weight");
    hist2d_rresultptr_map["h2d_crossx_pair_pt_dr_w_signal_cuts"] = df_single_b_crossx_weighted.Histo2D(
        ROOT::RDF::TH2DModel("h2d_crossx_pair_pt_dr_w_signal_cuts", ";p_{T}^{pair} [GeV];#DeltaR", npt, ptbins, 50, 0.05, 1.0),
        "pair_pt", "dr", "crossx_weight");

    hist3d_rresultptr_map["h3d_crossx_minv_vs_pair_eta_vs_pair_pt_w_signal_cuts"] = df_single_b_crossx_weighted.Histo3D(
        ROOT::RDF::TH3DModel("h3d_crossx_minv_vs_pair_eta_vs_pair_pt_w_signal_cuts", ";p_{T}^{pair} [GeV];#eta^{pair};m_{#mu#mu} [GeV]", npt, ptbins, 44, eta_edges.data(), 50, minv_edges.data()),
        "pair_pt", "pair_eta", "minv", "crossx_weight");
    hist3d_rresultptr_map["h3d_crossx_dr_vs_pair_eta_vs_pair_pt_w_signal_cuts"] = df_single_b_crossx_weighted.Histo3D(
        ROOT::RDF::TH3DModel("h3d_crossx_dr_vs_pair_eta_vs_pair_pt_w_signal_cuts", ";p_{T}^{pair} [GeV];#eta^{pair};#DeltaR", npt, ptbins, 44, eta_edges.data(), 50, dr_edges.data()),
        "pair_pt", "pair_eta", "dr", "crossx_weight");

    // --- pT_bins_150 variants ---
    {
        const int    npt150    = (int)(pms.pT_bins_150.size() - 1);
        const double* ptbins150 = pms.pT_bins_150.data();
        const auto eta_edges150 = make_unif_edges(44, -2.4, 2.4);
        const auto dr_edges150  = make_unif_edges(50, 0.05, 1.0);

        hist2d_rresultptr_map["h2d_crossx_pt_150_pair_eta_binned_w_signal_cuts"] = df_single_b_crossx_weighted.Histo2D(
            ROOT::RDF::TH2DModel("h2d_crossx_pt_150_pair_eta_binned_w_signal_cuts", ";p_{T}^{pair} [GeV];#eta^{pair}", npt150, ptbins150, 44, -2.4, 2.4),
            "pair_pt", "pair_eta", "crossx_weight");
        hist3d_rresultptr_map["h3d_crossx_dr_vs_pair_eta_vs_pt_150_w_signal_cuts"] = df_single_b_crossx_weighted.Histo3D(
            ROOT::RDF::TH3DModel("h3d_crossx_dr_vs_pair_eta_vs_pt_150_w_signal_cuts", ";p_{T}^{pair} [GeV];#eta^{pair};#DeltaR", npt150, ptbins150, 44, eta_edges150.data(), 50, dr_edges150.data()),
            "pair_pt", "pair_eta", "dr", "crossx_weight");
    }

    std::cout << "[PP] FillHistogramsCrossx completed" << std::endl;
}
