#include "SingleBAnalysisBase.cxx"

// -------------------------- PP analysis class --------------------------

class SingleBAnalysisPP : public SingleBAnalysisBase{
private:
    std::string infile_relative_path;
    std::string outfile_relative_path_truncated;

public:

    // 1/(integrated luminosity) factor to normalize to crossx (unit: pb)
    double crossx_factor;
    
    SingleBAnalysisPP(std::string infile_relative_path_input, std::string outfile_relative_path_input){
        infile_relative_path = infile_relative_path_input;
        outfile_relative_path_truncated = outfile_relative_path_input;
        Initialize();
    }

    void RunAnalysis();
};

// -------------------------- PP-class member functions --------------------------

void SingleBAnalysisPP::RunAnalysis(){

    for (auto bin : pT_bins_80){
        std::cout << bin << ", ";
    }
    std::cout << std::endl;
    // ----------------------------------------------------------------------------------------

    // Create RDataFrame and apply initial cuts
    muon_pair_input_file_name = dataset_base_dir + infile_relative_path;
    RDataFrame df("muon_pair_tree_sign2", muon_pair_input_file_name.c_str());
    auto rdf_with_signal_cuts = df.Filter(signal_cuts.c_str());


    // Create 2D histogram model for metrics calculation
    auto th2dmodel_signal_cuts = ROOT::RDF::TH2DModel(
        "h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts", "Signal (data-like cuts);p_{T} [GeV];#eta",
        int(pT_bins_80.size() - 1), pT_bins_80.data(),
        int(eta_bins.size() - 1), eta_bins.data()
    );

    // ----------------------------------------------------------------------------------------

    // Book 1D histograms
    auto h_pT_80 = rdf_with_signal_cuts.Histo1D(
        {"h_pT_80", "Pair pT (log bins up to 80 GeV)", int(pT_bins_80.size()-1), pT_bins_80.data()}, 
        "pair_pt", "weight");

    auto h_pT_200 = rdf_with_signal_cuts.Histo1D(
        {"h_pT_200", "Pair pT (log bins up to 200 GeV)", int(pT_bins_200.size()-1), pT_bins_200.data()}, 
        "pair_pt", "weight");
    
    // Book 2D histograms for metric calculations
    auto h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts = rdf_with_signal_cuts.Histo2D(th2dmodel_signal_cuts, "pair_pt", "pair_eta", "weight");

    // ----------------------------------------------------------------------------------------

    // Write results
    std::string outfile_name = dataset_base_dir + outfile_relative_path_truncated + ".root";
    TFile out(outfile_name.c_str(), "RECREATE");
    h_pT_80->Write();
    h_pT_200->Write();
    h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts->SetTitle("");
    h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts->SetStats(0);
    h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts->Scale(crossx_factor);
    h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts->Write();
    out.Close();
}


// -------------------------- running pp analysis -------------------------- 

void simplified_single_b_analysis_pp(){
    //Run PP analyses

    double crossx_factor_pp_run2 = 1/256.8; // 256.793
    double crossx_factor_pp_24_2mu4 = 1/410.815;
    double crossx_factor_pp_24_mu4mu4noL1 = 1/113.999;

    SingleBAnalysisPP pprun2("pp_run2/muon_pairs_pp_run2_old_res_cut.root", "pp_run2/pp_run2_single_b_ana_hists");
    pprun2.crossx_factor = crossx_factor_pp_run2;
    pprun2.RunAnalysis();

    SingleBAnalysisPP pp2024_mu4mu4noL1("pp_2024/muon_pairs_pp_2024_mu4_mu4noL1.root", "pp_2024/pp_2024_single_b_ana_hists_mu4_mu4noL1");
    pp2024_mu4mu4noL1.crossx_factor = crossx_factor_pp_24_mu4mu4noL1;
    pp2024_mu4mu4noL1.RunAnalysis();

    SingleBAnalysisPP pp2024_2mu4("pp_2024/muon_pairs_pp_2024_2mu4.root", "pp_2024/pp_2024_single_b_ana_hists_2mu4");
    pp2024_2mu4.crossx_factor = crossx_factor_pp_24_2mu4;
    pp2024_2mu4.RunAnalysis();
}
