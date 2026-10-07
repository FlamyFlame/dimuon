#include "RDFBasedHistFillingPythia.h"

#include <memory>
#include "TFile.h"
#include "TTree.h"

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

    if (!input_file_override.empty() || !output_file_override.empty() || muon_fiducial_target){
        if (input_file_override.empty() || output_file_override.empty())
            throw std::runtime_error("RDFBasedHistFillingPythia: input_file_override and output_file_override "
                "must be set TOGETHER (and are required with muon_fiducial_target).");
        if (isPrivate)
            throw std::runtime_error("RDFBasedHistFillingPythia: file overrides are for central-format "
                "(external-mode) NTuple outputs, not the standalone private sample.");
        const bool target_input = input_file_override.find("_muon_fid_target") != std::string::npos;
        if (muon_fiducial_target != target_input)
            throw std::runtime_error("RDFBasedHistFillingPythia: muon_fiducial_target=" + std::to_string(muon_fiducial_target)
                + " but the input name " + (target_input ? "carries" : "lacks") + " '_muon_fid_target': " + input_file_override);
        if (output_file_override == input_file_override)
            throw std::runtime_error("RDFBasedHistFillingPythia: output_file_override == input_file_override");
        // TARGET input must be PRODUCTION-WEIGHTED (sigma*eps_filt/N per slice, PythiaAlgCoreT
        // truth_external_sigma_eff_nb): a hadd of raw-counts slices (weight = 1) would mix the pT-hat
        // slices by event count -- slice-dependent, cancels nowhere in A_sig. Every meta entry (one per
        // hadded slice) must carry external_sigma_eff_nb > 0.
        if (muon_fiducial_target) {
            std::unique_ptr<TFile> fmeta(TFile::Open(input_file_override.c_str(), "READ"));
            if (!fmeta || fmeta->IsZombie())
                throw std::runtime_error("RDFBasedHistFillingPythia: cannot open " + input_file_override);
            TTree* meta = dynamic_cast<TTree*>(fmeta->Get("meta_tree_out"));
            if (!meta || !meta->GetBranch("external_sigma_eff_nb") || meta->GetEntries() == 0)
                throw std::runtime_error("RDFBasedHistFillingPythia: TARGET input " + input_file_override
                    + " has no meta_tree_out/external_sigma_eff_nb -- not a production-weighted NTuple.");
            double sig_eff = -1.;
            meta->SetBranchAddress("external_sigma_eff_nb", &sig_eff);
            for (Long64_t i = 0; i < meta->GetEntries(); ++i) {
                meta->GetEntry(i);
                if (!(sig_eff > 0.))
                    throw std::runtime_error("RDFBasedHistFillingPythia: TARGET input meta entry "
                        + std::to_string(i) + " has external_sigma_eff_nb <= 0 (raw-counts slice?)");
            }
            meta->ResetBranchAddresses();
        }
        input_files.push_back(input_file_override);
        output_file = output_file_override;
        infile_var1D_json = "var1D_pythia_truth.json";
        return;
    }

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
