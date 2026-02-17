#include "RDFBasedHistFillingPowheg.h"

void RDFBasedHistFillingPowheg::InitializePowhegImpl(){
    InitializePowhegCommon();
    InitializePowhegExtra();
}

void RDFBasedHistFillingPowheg::InitializePowhegCommon(){
    weight_specifier_to_column_map[""] = "weight_norm";
}

void RDFBasedHistFillingPowheg::SetIOPathsHook(){
    std::string powheg_dir = "/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample/";
    std::string run_year_str = std::to_string(run_year);

    // e.g, "bb" --> {"fullsim" --> FULLSIM_BB_DIR, "truth" --> TRUTH_BB_DIR}
    std::map<std::string, std::map<std::string, std::string>> mc_mode_to_data_subdir_map;
    std::map<std::string, std::map<std::string, std::string>> mc_mode_to_mpair_infile_name_map;
    
    for (std::string mc_mode : {"bb", "cc"}){
        mc_mode_to_data_subdir_map[mc_mode]["fullsim_overlay_run2"] = "user.yuhang.TrigRates.dimuon.PowhegPythia.fullsimOverlay.PbPb" + run_year_str + "." + mc_mode + ".Feb2026.v1._MYSTREAM/";
        mc_mode_to_data_subdir_map[mc_mode]["fullsim_overlay_run3"] = "user.yuhang.TrigRates.dimuon.PowhegPythia.fullsimOverlay.PbPb" + run_year_str + "." + mc_mode + ".Feb2026.v1._MYSTREAM/";
        mc_mode_to_data_subdir_map[mc_mode]["fullsim_run2"]         = "user.yuhang.TrigRates.dimuon.PowhegPythia.fullsim.pp17" + "." + mc_mode + ".Feb2026.v1._MYSTREAM/";
        mc_mode_to_data_subdir_map[mc_mode]["fullsim_run3"]         = "user.yuhang.TrigRates.dimuon.PowhegPythia.fullsim.pp24" + "." + mc_mode + ".Feb2026.v1._MYSTREAM/";
        mc_mode_to_data_subdir_map[mc_mode]["truth"]                = mc_mode + "_evgen_truth_full_sample/";

        mc_mode_to_mpair_infile_name_map[mc_mode]["fullsim_overlay_run2"]   = "muon_pairs_powheg_" + mc_mode + "_fullsim_overlay_PbPb" + run_year_str + "_w_truth.root";
        mc_mode_to_mpair_infile_name_map[mc_mode]["fullsim_overlay_run3"]   = "muon_pairs_powheg_" + mc_mode + "_fullsim_overlay_PbPb" + run_year_str + "_w_truth.root";
        mc_mode_to_mpair_infile_name_map[mc_mode]["fullsim_run2"]           = "muon_pairs_powheg_" + mc_mode + "_fullsim_pp17_w_truth.root";
        mc_mode_to_mpair_infile_name_map[mc_mode]["fullsim_run3"]           = "muon_pairs_powheg_" + mc_mode + "_fullsim_pp24_w_truth.root";
        mc_mode_to_mpair_infile_name_map[mc_mode]["truth"]                  = "muon_pairs_powheg_" + mc_mode + "_truth.root";
    }

    std::string mctype = "";
    if (is_fullsim_overlay){
        if (run_year == 18){
            mctype = "fullsim_overlay_run2";
        }
        else if (run_year >= 23 && run_year <= 26){
            mctype = "fullsim_overlay_run3";
        }
        else{
            throw std::runtime_error("For Powheg fullsim overlay, Run year must be 18/23/24/25/26! Current input invalid: " + run_year_str);
        }
    } else if (is_fullsim){
        if (run_year == 17){
            mctype = "fullsim_run2";
        }
        else if (run_year == 24){
            mctype = "fullsim_run3";
        }
        else{
            throw std::runtime_error("For Powheg fullsim, Run year must be 17/24! Current input invalid: " + run_year_str);
        }
    } else{ // evgen/truth only
        mctype = "truth";
    }

    for (std::string mc_mode : {"bb", "cc"}){
        input_files.push_back(powheg_dir + mc_mode_to_data_subdir_map[mc_mode][mctype] + mc_mode_to_mpair_infile_name_map[mc_mode][mctype]);
    }

    infile_var1D_json = is_fullsim  ? "var1D_powheg_fullsim.json"
                                    : "var1D_powheg_truth.json";

    std::string dt_suffix = is_fullsim_overlay  ? "_fullsim_overlay_pbpb" + run_year_str
                                                : ( is_fullsim  ? "_fullsim_pp" + run_year_str
                                                                : "_truth");

    output_file = powheg_dir + "histograms_powheg" + dt_suffix + ".root";
}

void RDFBasedHistFillingPowheg::CreateBaseRDFsPowhegImpl(){
    CreateBaseRDFsPowhegCommon();
    CreateBaseRDFsPowhegExtra();
}

void RDFBasedHistFillingPowheg::CreateBaseRDFsPowhegCommon(){

    // Sum over #entries before filter & calculate weight for normalizing histogram intergrals to crossx * filter efficiencies
    // Where filter efficiencies include all cuts applied in analysis 

    const double nentries_before_cuts_sum = SumMetaNentriesBeforeFilter(input_files);
    if (nentries_before_cuts_sum <= 0.0) throw std::runtime_error("meta_tree-summed nentries_before_cuts is non-positive");

    std::cout << "nentries_before_cuts_sum = " << nentries_before_cuts_sum << "\n";

    ROOT::RDF::RNode& df_ss = map_at_checked(df_map, "df_ss", "RDFBasedHistFillingPowheg::FillHistograms: df_map.at(df_ss)");
    auto df_ss_weighted = df_ss.Define("weight_norm", [nentries_before_cuts_sum](double w){ return w / nentries_before_cuts_sum; }, {"weight"});
    df_map.emplace("df_ss_weighted", df_ss_weighted);

    ROOT::RDF::RNode& df_op = map_at_checked(df_map, "df_op", "RDFBasedHistFillingPowheg::FillHistograms: df_map.at(df_op)");
    auto df_op_weighted = df_op.Define("weight_norm", [nentries_before_cuts_sum](double w){ return w / nentries_before_cuts_sum; }, {"weight"});
    df_map.emplace("df_op_weighted", df_op_weighted);

    auto df_single_b_weighted = df_op_weighted.Filter("from_same_b");
    df_map.emplace("df_single_b_weighted", df_single_b_weighted);

}

void RDFBasedHistFillingPowheg::FillHistograms(){
    if (is_fullsim){
        FillHistogramsFullSim();
    } else{
        FillHistogramsTruth();
    }
}

void RDFBasedHistFillingPowheg::BuildFilterToVarListMapPowhegImpl(){
    // simple maps powheg commons
    BuildSimpleFilterToVarListMapPowhegCommon();

    // simple maps extras
    BuildSimpleFilterToVarListMapPowhegExtra();

    // leveled filters
    FlattenFilters(); // flatten reco effcy filters
    BuildFlattenedFilterToVarListMap(); // map flattened reco effcy filters to 1/2/3D variable list
}

void RDFBasedHistFillingPowheg::FlattenFilters(){
    FlattenFiltersPowhegCommon();
    FlattenFiltersExtra();
}

void RDFBasedHistFillingPowheg::BuildFlattenedFilterToVarListMap(){
    BuildFlattenedFilterToVarListMapPowhegCommon();
    BuildFlattenedFilterToVarListMapExtra();
}











