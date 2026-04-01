#include "RDFBasedHistFillingPowheg.cxx"
#include "../Utilities/GeneralUtils.h"
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

void RDFBasedHistFillingPowhegFullsimSingleMuon::SetIOPathsHook(){
    if (!is_fullsim || is_fullsim_overlay){
        throw std::runtime_error("RDFBasedHistFillingPowhegFullsimSingleMuon: requires non-overlay fullsim mode.");
    }
    if (run_year != 17){
        throw std::runtime_error("RDFBasedHistFillingPowhegFullsimSingleMuon: only run_year=17 is supported. Got: " + std::to_string(run_year));
    }

    // The single-muon file uses a single tree with all muon charges.
    // The base class loads both tree_ss and tree_op; we point both to "muon_tree"
    // and only use df_ss downstream (df_op is a duplicate and ignored).
    tree_ss = "muon_tree";
    tree_op = "muon_tree";

    const std::string powheg_dir = "/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample/";
    input_files.clear();

    for (const std::string& mc_mode : {"bb", "cc"}){
        const std::string input_file =
            powheg_dir
            + "user.yuhang.TrigRates.dimuon.PowhegPythia.fullsim.pp17." + mc_mode + ".Feb2026.v1._MYSTREAM/"
            + "single_muon_trees_powheg_" + mc_mode + "_fullsim_pp17.root";

        std::ifstream fin(input_file);
        if (!fin.good()){
            throw std::runtime_error("RDFBasedHistFillingPowhegFullsimSingleMuon: missing input file: " + input_file);
        }
        input_files.push_back(input_file);
    }

    infile_var1D_json = "var1D_powheg_fullsim_single_muon.json";
    output_file = powheg_dir + "histograms_powheg_fullsim_single_muon_pp17.root";
}

void RDFBasedHistFillingPowhegFullsimSingleMuon::InitializePowhegExtra(){
    levels_reco_effcy_filters = {
        {"_sign1", "_sign2"},
        {"", "_pass_medium", "_pass_tight"}
    };

    levels_detector_response_filters = {
        {"_sign1", "_sign2"},
        {"_pass_medium", "_pass_tight"}
    };
}

void RDFBasedHistFillingPowhegFullsimSingleMuon::FlattenFiltersExtra(){
    HistFillUtils::flatten_levels(levels_reco_effcy_filters,      reco_effcy_filters);
    HistFillUtils::flatten_levels(levels_detector_response_filters, detector_response_filters);
}

void RDFBasedHistFillingPowhegFullsimSingleMuon::BuildFlattenedFilterToVarListMapExtra(){
    reco_effcy_var1Ds = {"truth_pt", "truth_eta", "truth_phi"};

    reco_effcy_var2Ds = {
        {"truth_pt",  "truth_eta"},
        {"truth_pt",  "truth_phi"},
        {"truth_eta", "truth_phi"}
    };

    reco_effcy_var3Ds = {
        {"truth_pt", "truth_eta", "truth_phi"}
    };

    detec_resp_var1Ds = {
        "truth_pt", "pt",
        "truth_eta", "eta",
        "truth_phi", "phi"
    };

    detec_resp_var2Ds = {
        {"truth_pt",  "pt"},
        {"truth_eta", "eta"},
        {"truth_phi", "phi"}
    };

    for (const auto& filter : reco_effcy_filters)
        InsertOrAppend(df_filter_and_weight_to_var1D_list_map, std::make_pair(filter, ""), reco_effcy_var1Ds);
    for (const auto& filter : reco_effcy_filters)
        InsertOrAppend(df_filter_and_weight_to_var2D_list_map, std::make_pair(filter, ""), reco_effcy_var2Ds);
    for (const auto& filter : reco_effcy_filters)
        InsertOrAppend(df_filter_and_weight_to_var3D_list_map, std::make_pair(filter, ""), reco_effcy_var3Ds);

    for (const auto& filter : detector_response_filters)
        InsertOrAppend(df_filter_and_weight_to_var1D_list_map, std::make_pair(filter, ""), detec_resp_var1Ds);
    for (const auto& filter : detector_response_filters)
        InsertOrAppend(df_filter_and_weight_to_var2D_list_map, std::make_pair(filter, ""), detec_resp_var2Ds);
}

void RDFBasedHistFillingPowhegFullsimSingleMuon::CreateBaseRDFsExtra(){
    // Fully replaces the powheg-common weight setup.
    // The single-muon tree has no "weight" column; use MuonObj.ev_weight directly.
    weight_specifier_to_column_map[""] = "weight_norm";

    ROOT::RDF::RNode& df_base = map_at_checked(df_map, "df_ss",
        "RDFBasedHistFillingPowhegFullsimSingleMuon::CreateBaseRDFsExtra: df_map.at(df_ss)");

    auto df_weighted = df_base.Define("weight_norm", "MuonObj.ev_weight");

    auto df_sign1_weighted = df_weighted.Filter("MuonObj.truth_charge == 1");
    auto df_sign2_weighted = df_weighted.Filter("MuonObj.truth_charge == -1");

    df_map.emplace("df_sign1_weighted", df_sign1_weighted);
    df_map.emplace("df_sign2_weighted", df_sign2_weighted);

    df_map.emplace("df_sign1_pass_medium_weighted", df_sign1_weighted.Filter("MuonObj.pass_medium"));
    df_map.emplace("df_sign1_pass_tight_weighted",  df_sign1_weighted.Filter("MuonObj.pass_tight"));
    df_map.emplace("df_sign2_pass_medium_weighted", df_sign2_weighted.Filter("MuonObj.pass_medium"));
    df_map.emplace("df_sign2_pass_tight_weighted",  df_sign2_weighted.Filter("MuonObj.pass_tight"));
}

void RDFBasedHistFillingPowhegFullsimSingleMuon::FillHistogramsFullSim(){
    FillHistogramsFullSimDetecResp();
    FillHistogramsFullSimRecoEffcies();
}

void RDFBasedHistFillingPowhegFullsimSingleMuon::FillHistogramsFullSimDetecResp(){
    try{
        for (const auto& filter : detector_response_filters){
            const std::string df_name = "df" + filter + "_weighted";
            ROOT::RDF::RNode& node = map_at_checked(df_map, df_name,
                Form("FillHistogramsFullSimDetecResp: df_map.at(%s)", df_name.c_str()));
            FillHistogramsSingleDataFrame(filter, "", node);
        }
    } catch (const std::out_of_range& e){
        std::cerr << "FillHistogramsFullSimDetecResp:: out_of_range: " << e.what() << std::endl;
    } catch (const std::runtime_error& e){
        std::cerr << "FillHistogramsFullSimDetecResp:: runtime_error: " << e.what() << std::endl;
    }
}

void RDFBasedHistFillingPowhegFullsimSingleMuon::FillHistogramsFullSimRecoEffcies(){
    try{
        for (const auto& filter : reco_effcy_filters){
            const std::string df_name = "df" + filter + "_weighted";
            ROOT::RDF::RNode& node = map_at_checked(df_map, df_name,
                Form("FillHistogramsFullSimRecoEffcies: df_map.at(%s)", df_name.c_str()));
            FillHistogramsSingleDataFrame(filter, "", node);
        }
    } catch (const std::out_of_range& e){
        std::cerr << "FillHistogramsFullSimRecoEffcies:: out_of_range: " << e.what() << std::endl;
    } catch (const std::runtime_error& e){
        std::cerr << "FillHistogramsFullSimRecoEffcies:: runtime_error: " << e.what() << std::endl;
    }
}
