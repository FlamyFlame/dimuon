#include "RDFBasedHistFillingPythia.h"

void RDFBasedHistFillingPythia::InitializePythiaImpl(){
    InitializePythiaCommon();
    InitializePythiaExtra();
}

void RDFBasedHistFillingPythia::InitializePythiaCommon(){
    weight_specifier_to_column_map[""] = "weight";
}

void RDFBasedHistFillingPythia::SetIOPathsHook(){
    if (is_fullsim || is_fullsim_overlay){
        throw std::runtime_error("RDFBasedHistFillingPythia: fullsim/fullsim-overlay path configuration is a placeholder and not implemented yet.");
    }

    const std::string cut_suffix = with_data_resonance_cuts
        ? "_with_data_resonance_cuts"
        : "_no_data_resonance_cuts";

    input_files.clear();

    if (isPrivate){
        input_files.push_back(
            "/usatlas/u/yuhanguo/usatlasdata/pythia_private_sample/muon_pairs_pythia_combined"
            + cut_suffix + ".root"
        );

        output_file =
            "/usatlas/u/yuhanguo/usatlasdata/pythia_private_sample/histograms_pythia_combined"
            + cut_suffix + ".root";
    } else {
        std::string ecom_tag;
        std::string ecom_subdir;
        if (std::abs(E_COM - 5.36) < 0.01) {
            ecom_tag = "5p36TeV";
            ecom_subdir = "pythia_5p36TeV";
        } else if (std::abs(E_COM - 5.02) < 0.01) {
            ecom_tag = "5p02TeV";
            ecom_subdir = "pythia_5TeV";
        } else {
            throw std::runtime_error("RDFBasedHistFillingPythia: E_COM must be 5.02 or 5.36 TeV. Got " + std::to_string(E_COM));
        }

        const std::string base_dir =
            "/usatlas/u/yuhanguo/usatlasdata/pythia_truth_full_sample/" + ecom_subdir + "/";

        input_files.push_back(base_dir + "muon_pairs_pythia_" + ecom_tag + cut_suffix + ".root");
        output_file = base_dir + "histograms_pythia_" + ecom_tag + cut_suffix + ".root";
    }

    infile_var1D_json = "var1D_pythia_truth.json";
}

void RDFBasedHistFillingPythia::CreateBaseRDFsPythiaImpl(){
    CreateBaseRDFsPythiaCommon();
    CreateBaseRDFsPythiaExtra();
}

void RDFBasedHistFillingPythia::CreateBaseRDFsPythiaCommon(){
    ROOT::RDF::RNode& df_ss = map_at_checked(df_map, "df_ss", "RDFBasedHistFillingPythia::CreateBaseRDFsPythiaCommon: df_ss");
    ROOT::RDF::RNode& df_op = map_at_checked(df_map, "df_op", "RDFBasedHistFillingPythia::CreateBaseRDFsPythiaCommon: df_op");

    df_map.emplace("df_ss_weighted", df_ss);
    df_map.emplace("df_op_weighted", df_op);
}

void RDFBasedHistFillingPythia::FillHistograms(){
    if (is_fullsim || is_fullsim_overlay){
        FillHistogramsFullSim();
    } else {
        FillHistogramsTruth();
    }
}
