#include "RDFBasedHistFillingPowheg.cxx"
#include "../Utilities/GeneralUtils.h"

void RDFBasedHistFillingPowhegFullsim::InitializePowhegFullsimExtra(){
    
    levels_reco_effcy_filters = {
        {"_ss", "_op", "_single_b"}, // all truth-same/op-sign pairs & truth single-b-signal pairs
        {"", "_pair_pass_medium", "_pair_pass_tight"}
    };

    levels_detector_response_filters = {
        {"_single_b"}, // truth single-b-signal pairs only
        {"_pair_pass_medium", "_pair_pass_tight"}
    };
}

void RDFBasedHistFillingPowhegFullsim::FlattenFiltersPowhegFullsim(){
    TrigEffcyUtils::flatten_levels(levels_reco_effcy_filters, reco_effcy_filters);
    TrigEffcyUtils::flatten_levels(levels_detector_response_filters, detector_response_filters);
}

void RDFBasedHistFillingPowhegFullsim::BuildFlattenedFilterToVarListMapPowhegFullsim(){
    reco_effcy_var1Ds = {
        "truth_pair_pt", "truth_pair_eta", "truth_dr"
    };

    reco_effcy_var2Ds = {
        {"truth_pair_pt", "truth_pair_eta"}, 
        {"truth_pair_pt", "truth_dr"}, 
        {"truth_pair_eta", "truth_dr"}
    };

    reco_effcy_var3Ds = {
        {"truth_pair_pt", "truth_pair_eta", "truth_dr"}
    };

    detec_resp_var1Ds = {
        "truth_pair_pt", "pair_pt",
        "truth_minv", "minv",
        "truth_dr", "dr"
    };
    
    detec_resp_var2Ds = {
        {"truth_pair_pt", "pair_pt"},
        {"truth_minv", "minv"},
        {"truth_dr", "dr"}
    };

    for (auto filter : reco_effcy_filters)  InsertOrAppend(df_filter_and_weight_to_var1D_list_map, {filter, ""}, reco_effcy_var1Ds);
    for (auto filter : reco_effcy_filters)  InsertOrAppend(df_filter_and_weight_to_var2D_list_map, {filter, ""}, reco_effcy_var2Ds);
    for (auto filter : reco_effcy_filters)  InsertOrAppend(df_filter_and_weight_to_var3D_list_map, {filter, ""}, reco_effcy_var3Ds);

    for (auto filter : detector_response_filters)   InsertOrAppend(df_filter_and_weight_to_var1D_list_map, {filter, ""}, detec_resp_var1Ds);    
    for (auto filter : detector_response_filters)   InsertOrAppend(df_filter_and_weight_to_var2D_list_map, {filter, ""}, detec_resp_var2Ds);
}

void RDFBasedHistFillingPowhegFullsim::CreateBaseRDFsPowhegFullsimExtra(){
    for (std::string pair_catgr : {"_ss", "_op", "_signal_b"}){
        std::string df_name = "df" + pair_catgr;
        ROOT::RDF::RNode& node = map_at_checked(df_map, df_name + "_weighted", Form("FillHistogramsSingleMuonEffcy: df_map.at(%s)", (df_name + "_weighted").c_str()));

        df_map.emplace(df_name + "_pair_pass_medium_weighted" , node.Filter("pair_pass_medium"));
        df_map.emplace(df_name + "_pair_pass_tight_weighted"  , node.Filter("pair_pass_tight"));
    }
}

void RDFBasedHistFillingPowhegFullsim::FillHistogramsFullSim(){
    FillHistogramsFullSimDetecResp();
    FillHistogramsFullSimRecoEffcies();
}

void RDFBasedHistFillingPowhegFullsim::FillHistogramsFullSimDetecResp(){
    try{
        for (std::string pair_catgr : {"_signal_b"}){
            for (auto quality_catgr : {"_pair_pass_medium", "_pair_pass_tight"}){ // mu4 selection

                std::string filter = pair_catgr + quality_catgr;
                std::string df_name = "df" + filter + "_weighted"; // e.g, df_ss_mu1passmu4
                ROOT::RDF::RNode& node = map_at_checked(df_map, df_name, Form("FillHistogramsSingleMuonEffcy: df_map.at(%s)", df_name.c_str()));

                FillHistogramsSingleDataFrame(filter, "", node); // do not write the sub-dataframe histograms in output file
            }
        }
    } catch(const std::out_of_range& e){
        std::cerr << "FillHistogramsFullSimRecoEffcies:: out_of_range exception caught: " << e.what() << std::endl;
    } catch (const std::runtime_error& e) {
        std::cerr << "FillHistogramsFullSimRecoEffcies:: RDF runtime error: " << e.what() << std::endl;
    }
}

void RDFBasedHistFillingPowhegFullsim::FillHistogramsFullSimRecoEffcies(){
    try{
        for (std::string pair_catgr : {"_ss", "_op", "_signal_b"}){
            for (auto quality_catgr : {"", "_pair_pass_medium", "_pair_pass_tight"}){ // mu4 selection

                std::string filter = pair_catgr + quality_catgr;
                std::string df_name = "df" + filter + "_weighted"; // e.g, df_ss_mu1passmu4
                ROOT::RDF::RNode& node = map_at_checked(df_map, df_name, Form("FillHistogramsSingleMuonEffcy: df_map.at(%s)", df_name.c_str()));

                FillHistogramsSingleDataFrame(filter, "", node); // do not write the sub-dataframe histograms in output file
            }
        }
    } catch(const std::out_of_range& e){
        std::cerr << "FillHistogramsFullSimRecoEffcies:: out_of_range exception caught: " << e.what() << std::endl;
    } catch (const std::runtime_error& e) {
        std::cerr << "FillHistogramsFullSimRecoEffcies:: RDF runtime error: " << e.what() << std::endl;
    }
}
