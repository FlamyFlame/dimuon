#include "RDFBasedHistFillingData.cxx"

void RDFBasedHistFillingPP::SetIOPathsHook(){
    infile_var1D_json = "var1D_pp.json";

    if (run_year == 24){
        input_files.push_back("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/muon_pairs_pp_2024" + trig_suffix + "_res_cut_v2.root");
        output_file = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/histograms_real_pairs_pp_2024" + out_file_suffix + ".root";
    } else if (run_year == 17){
        input_files.push_back("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_run2/muon_pairs_pp_2017" + trig_suffix + "_res_cut_v2.root");
        output_file = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_run2/histograms_real_pairs_pp_2017" + out_file_suffix + ".root";
    }
}

void RDFBasedHistFillingPP::InitializePPExtra(){
	if (run_year == 17){
        trigs = {"_mu4", "_2mu4"};
        trigs_pair = {{"_2mu4","_mu4"}};
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

void RDFBasedHistFillingPP::TrigEffcyFiltersPrePostSumFlatteningExtra(){

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
    TrigEffcyUtils::flatten_levels(levels_trg_effcy_filters_to_be_summed_w_musign_summing, trg_effcy_filters_to_be_summed_w_musign_summing);    

    // build post-sum levels, with mu-sign summing
    TrigEffcyUtils::write_post_sum_levels(levels_trg_effcy_filters_1D_pre_sum,
                          levels_trg_effcy_to_be_summed_w_musign_summing,
                          levels_trg_effcy_filters_1D_post_sum_w_musign_summing);

    TrigEffcyUtils::write_post_sum_levels(levels_trg_effcy_filters_2D_3D_pre_sum,
                          levels_trg_effcy_to_be_summed_w_musign_summing,
                          levels_trg_effcy_filters_2D_3D_post_sum_w_musign_summing);

    // flatten post-sum levels, with mu-sign summing
    TrigEffcyUtils::flatten_levels(levels_trg_effcy_filters_1D_post_sum_w_musign_summing, trg_effcy_filters_1D_post_sum_w_musign_summing);
    TrigEffcyUtils::flatten_levels(levels_trg_effcy_filters_2D_3D_post_sum_w_musign_summing, trg_effcy_filters_2D_3D_post_sum_w_musign_summing);
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
                df_map.at(df_name) = df_map.at(df_name).Define("second_muon_good_acceptance",  "pt2nd >= 6 && ((eta2nd > 1.1 && eta2nd < 2.3) || (eta2nd > -2.3 && eta2nd < -1.2))");

                df_map.emplace(df_name + "_sign1", map_at_checked(df_map, df_name, Form("FillHistogramsSingleMuonEffcy: df_map.at(%s)", df_name.c_str())).Filter("charge2nd > 0")); // e.g, df_ss_mu1passmu4_sign1
                df_map.emplace(df_name + "_sign2", map_at_checked(df_map, df_name, Form("FillHistogramsSingleMuonEffcy: df_map.at(%s)", df_name.c_str())).Filter("charge2nd < 0"));

                for (auto mu_sign : {"_sign1", "_sign2"}){
                    std::string df_name = "df" + pair_sign + mu4sel + mu_sign;
                    df_map.emplace(df_name + "_mu4", map_at_checked(df_map, df_name, Form("FillHistogramsSingleMuonEffcy: df_map.at(%s)", df_name.c_str())));
                    df_map.emplace(df_name + "_mu4_mu4noL1", map_at_checked(df_map, df_name, Form("FillHistogramsSingleMuonEffcy: df_map.at(%s)", df_name.c_str())).Filter("passmu4mu4noL1"));
                    df_map.emplace(df_name + "_2mu4", map_at_checked(df_map, df_name, Form("FillHistogramsSingleMuonEffcy: df_map.at(%s)", df_name.c_str())).Filter("pass2mu4"));
                    df_map.emplace(df_name + "_2mu4_AND_mu4_mu4noL1", map_at_checked(df_map, df_name, Form("FillHistogramsSingleMuonEffcy: df_map.at(%s)", df_name.c_str())).Filter("pass2mu4 && passmu4mu4noL1"));

                    for (auto trg : {"_mu4", "_mu4_mu4noL1", "_2mu4", "_2mu4_AND_mu4_mu4noL1"}){
                        std::string df_name = "df" + pair_sign + mu4sel + mu_sign + trg;
                        df_map.emplace(df_name + "_sepr", map_at_checked(df_map, df_name, Form("FillHistogramsSingleMuonEffcy: df_map.at(%s)", df_name.c_str())).Filter("passSeparated"));
                        df_map.emplace(df_name + "_good_accept", map_at_checked(df_map, df_name, Form("FillHistogramsSingleMuonEffcy: df_map.at(%s)", df_name.c_str())).Filter("second_muon_good_acceptance"));
                        
                        for (auto bias : trg_effcy_biases){ // additional selection / bias in data sample
                            std::string df_name = "df" + pair_sign + mu4sel + mu_sign + trg + bias;
                            std::string filter = df_name.substr(2);

                            FillHistogramsSingleDataFrame(filter, map_at_checked(df_map, df_name, Form("FillHistogramsSingleMuonEffcy: df_map.at(%s)", df_name.c_str())), false); // do not write the sub-dataframe histograms in output file
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

void RDFBasedHistFillingPP::SumSingleMuonTrigEffHists(){
    RDFBasedHistFillingData::SumSingleMuonTrigEffHists();

    // 1D, with mu-sign summing
    TrigEffcyUtils::SumTrigEffHistsGeneric<TH1D, std::string>(
        single_muon_trig_effcy_var1Ds,
        trg_effcy_filters_1D_post_sum_w_musign_summing,
        trg_effcy_filters_to_be_summed_w_musign_summing,
        hist1D_map,
        [](const std::string& var) {
            return "h_" + var;
        }
    );

    // 2D, with mu-sign summing
    TrigEffcyUtils::SumTrigEffHistsGeneric<TH2D, std::array<std::string,2>>(
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
    TrigEffcyUtils::SumTrigEffHistsGeneric<TH3D, std::array<std::string,3>>(
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
