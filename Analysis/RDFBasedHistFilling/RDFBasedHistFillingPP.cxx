#include "RDFBasedHistFillingData.cxx"

void RDFBasedHistFillingPP::Initialize(){
    TriggerModeSettings();
    
    infile_var1D_json = "var1D_pp.json";

    run_year %= 2000;
    if (run_year == 24){
    	input_files.push_back("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/muon_pairs_pp_2024" + trig_suffix + "_res_cut_v2.root");
    	output_file = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/histograms_real_pairs_pp_2024" + trig_suffix + "_res_cut_v2.root";        
    } else if (run_year == 17){
        input_files.push_back("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_run2/muon_pairs_pp_2017" + trig_suffix + "_res_cut_v2.root");
        output_file = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_run2/histograms_real_pairs_pp_2017" + trig_suffix + "_res_cut_v2.root";
        trigs = {"_mu4", "_2mu4"};
        trigs_pair = {{"_2mu4","_mu4"}};

    }
	
    levels_trg_effcy_filters_1D_pre_sum = {{"_ss", "_op"},
                                        {"_mu1passmu4", "_mu2passmu4"},
                                        {"_sign1", "_sign2"},
                                        trigs,
                                        {"", "_good_accept", "_sepr"}};

    levels_trg_effcy_filters_2D_3D_pre_sum = {{"_ss", "_op"},
                                        {"_mu1passmu4", "_mu2passmu4"},
                                        {"_sign1", "_sign2"},
                                        trigs,
                                        {"", "_sepr"}};

    categories_essential = pair_signs;

	RDFBasedHistFillingBaseClass::Initialize();
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
                        
                        for (auto bias : {"", "_good_accept", "_sepr"}){ // additional selection / bias in data sample
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

void RDFBasedHistFillingPP::CalculateSingleMuonTrigEffcyRatios(){
    CalculateSingleMuonTrigEffcyRatiosHelper(musigns);
}

void RDFBasedHistFillingPP::MakeAndWriteSingleMuonPtTrigEffGraphs(){
    MakeAndWriteSingleMuonPtTrigEffGraphsHelper(musigns);
}

void RDFBasedHistFillingPP::OpenEffcyPtFitFile(){}

void RDFBasedHistFillingPP::MakeAndWriteDRTrigEffGraphs(){}
