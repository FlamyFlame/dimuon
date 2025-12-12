#include "MuonObjectsParamsAndHelpers/PbPbBaseClass.h"
#include "SingleBAnalysisBase.cxx"

// -------------------------- PbPb analysis class --------------------------

class SingleBAnalysisPbPb : public SingleBAnalysisBase, public PbPbBaseClass{
private:
    std::string infile_relative_path;
    std::string outfile_relative_path_truncated;

public:
    SingleBAnalysisPbPb(std::string infile_relative_path_input, std::string outfile_relative_path_input){
        infile_relative_path = infile_relative_path_input;
        outfile_relative_path_truncated = outfile_relative_path_input;
        Initialize(); // run the PbPb version of Initialize()
    }
    ~SingleBAnalysisPbPb(){}

    void Initialize();
    void RunAnalysis();
};


// -------------------------- PbPb-class member functions --------------------------

void SingleBAnalysisPbPb::Initialize(){
    SingleBAnalysisBase::Initialize();
    PbPbBaseClass::InitializePbPb();
}


void SingleBAnalysisPbPb::RunAnalysis(){

    // Create RDataFrame for the diff Centrality groups & signs, correctly weighted
    muon_pair_input_file_name = dataset_base_dir + infile_relative_path;

    RDataFrame df_ss("muon_pair_tree_sign1", muon_pair_input_file_name.c_str());    
    RDataFrame df_op("muon_pair_tree_sign2", muon_pair_input_file_name.c_str());    

    std::vector<ROOT::RDF::RNode> df_pbpb_run2_ctr_binned_list_ss_w_weight_w_signal_cuts;
    std::vector<ROOT::RDF::RNode> df_pbpb_run2_ctr_binned_list_op_w_weight_w_signal_cuts;

    
	auto df_ss_weight = df_ss.Define("weight_for_RAA",
	    [this](int avg_centrality, double weight) {
	        return this->CalculateWeightForRAA(avg_centrality, weight);
	    },
	    {"avg_centrality", "weight"}
	);

	auto df_op_weight = df_op.Define("weight_for_RAA",
	    [this](int avg_centrality, double weight) {
	        return this->CalculateWeightForRAA(avg_centrality, weight);
	    },
	    {"avg_centrality", "weight"}
	);

    auto df_ss_weight_w_signal_cuts = df_ss_weight.Filter(signal_cuts.c_str());
    auto df_op_weight_w_signal_cuts = df_op_weight.Filter(signal_cuts.c_str());

    // fill centrality-binned RDF lists
    for (int ictr = 0; ictr < nCtrBins; ictr++){

        int ctr_bin_low_edge = ctr_bin_edges.at(ictr);
        int ctr_bin_high_edge = ctr_bin_edges.at(ictr + 1);

        auto df_ctr_binned_ss_weight_w_signal_cuts = df_ss_weight_w_signal_cuts.Filter(
            [ctr_bin_low_edge, ctr_bin_high_edge](int centrality) {return centrality >= ctr_bin_low_edge && centrality < ctr_bin_high_edge;},
            {"avg_centrality"});

        auto df_ctr_binned_op_weight_w_signal_cuts = df_op_weight_w_signal_cuts.Filter(
            [ctr_bin_low_edge, ctr_bin_high_edge](int centrality) {return centrality >= ctr_bin_low_edge && centrality < ctr_bin_high_edge;},
            {"avg_centrality"});

        df_pbpb_run2_ctr_binned_list_ss_w_weight_w_signal_cuts.push_back(df_ctr_binned_ss_weight_w_signal_cuts);
        df_pbpb_run2_ctr_binned_list_op_w_weight_w_signal_cuts.push_back(df_ctr_binned_op_weight_w_signal_cuts);
    }

    // ----------------------------------------------------------------------------------------

    // Create 3D histogram model for crossx dependence on centrality, pair pt & pair eta
    auto th3dmodel_ss_crossx_signal_cuts = ROOT::RDF::TH3DModel(
        "h3d_ss_crossx_w_signal_cuts_vs_centr_vs_pair_eta_vs_pair_pt", ";p_{T} [GeV];#eta;Centrality",
        int(pT_bins_80.size() - 1), pT_bins_80.data(),
        int(eta_bins.size() - 1), eta_bins.data(),
        int(ctr_bin_edges_double.size() - 1), ctr_bin_edges_double.data()
    );
    auto th3dmodel_op_crossx_signal_cuts = ROOT::RDF::TH3DModel(
        "h3d_op_crossx_w_signal_cuts_vs_centr_vs_pair_eta_vs_pair_pt", ";p_{T} [GeV];#eta;Centrality",
        int(pT_bins_80.size() - 1), pT_bins_80.data(),
        int(eta_bins.size() - 1), eta_bins.data(),
        int(ctr_bin_edges_double.size() - 1), ctr_bin_edges_double.data()
    );

    // Book 3D histograms: crossx dependence on centrality, pair pt & pair eta
    auto h3d_ss_crossx_w_signal_cuts_vs_centr_vs_pair_eta_vs_pair_pt = df_ss_weight_w_signal_cuts.Histo3D(th3dmodel_ss_crossx_signal_cuts, "pair_pt", "pair_eta", "avg_centrality", "weight_for_RAA");
    auto h3d_op_crossx_w_signal_cuts_vs_centr_vs_pair_eta_vs_pair_pt = df_op_weight_w_signal_cuts.Histo3D(th3dmodel_op_crossx_signal_cuts, "pair_pt", "pair_eta", "avg_centrality", "weight_for_RAA");

    // ----------------------------------------------------------------------------------------
    // booking pT-pair 1D histograms

    std::vector<ROOT::RDF::RResultPtr<TH1D>> hop_pT_80_ctr_binned;
    std::vector<ROOT::RDF::RResultPtr<TH1D>> hss_pT_80_ctr_binned;
    std::vector<ROOT::RDF::RResultPtr<TH1D>> hop_pT_80_counts_ctr_binned;
    std::vector<ROOT::RDF::RResultPtr<TH1D>> hss_pT_80_counts_ctr_binned;
    std::vector<ROOT::RDF::RResultPtr<TH1D>> hop_pT_150_ctr_binned;
    std::vector<ROOT::RDF::RResultPtr<TH1D>> hss_pT_150_ctr_binned;
    std::vector<ROOT::RDF::RResultPtr<TH1D>> hop_pT_150_counts_ctr_binned;
    std::vector<ROOT::RDF::RResultPtr<TH1D>> hss_pT_150_counts_ctr_binned;
    std::vector<ROOT::RDF::RResultPtr<TH1D>> hop_pT_200_ctr_binned;
    std::vector<ROOT::RDF::RResultPtr<TH1D>> hss_pT_200_ctr_binned;
    std::vector<ROOT::RDF::RResultPtr<TH1D>> hop_pT_200_counts_ctr_binned;
    std::vector<ROOT::RDF::RResultPtr<TH1D>> hss_pT_200_counts_ctr_binned;
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h2d_ss_crossx_w_signal_cuts_ctr_binned;
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h2d_op_crossx_w_signal_cuts_ctr_binned;

    for (int ictr = 0; ictr < nCtrBins; ictr++){

        int ctr_bin_low_edge = ctr_bin_edges.at(ictr);
        int ctr_bin_high_edge = ctr_bin_edges.at(ictr + 1);

        try{
            auto df_ctr_binned_op_weight_w_signal_cuts = df_pbpb_run2_ctr_binned_list_op_w_weight_w_signal_cuts.at(ictr);
            auto df_ctr_binned_ss_weight_w_signal_cuts = df_pbpb_run2_ctr_binned_list_ss_w_weight_w_signal_cuts.at(ictr);
    
            // Book 1D histograms
            auto h_pT_80_op = df_ctr_binned_op_weight_w_signal_cuts.Histo1D(
                {Form("hop_pT_80_weighted_ctr%d_%d_w_signal_cuts", ctr_bin_low_edge, ctr_bin_high_edge), Form("h_pT_80_weighted_ctr%d_%d_w_signal_cuts", ctr_bin_low_edge, ctr_bin_high_edge), int(pT_bins_80.size()-1), pT_bins_80.data()}, 
                "pair_pt", "weight_for_RAA");
            hop_pT_80_ctr_binned.push_back(h_pT_80_op);

            auto h_pT_80_ss = df_ctr_binned_ss_weight_w_signal_cuts.Histo1D(
                {Form("hss_pT_80_weighted_ctr%d_%d_w_signal_cuts", ctr_bin_low_edge, ctr_bin_high_edge), Form("h_pT_80_weighted_ctr%d_%d_w_signal_cuts", ctr_bin_low_edge, ctr_bin_high_edge), int(pT_bins_80.size()-1), pT_bins_80.data()}, 
                "pair_pt", "weight_for_RAA");
            hss_pT_80_ctr_binned.push_back(h_pT_80_ss);

            auto h_pT_80_counts_op = df_ctr_binned_op_weight_w_signal_cuts.Histo1D(
                {Form("hop_pT_80_counts_ctr%d_%d_w_signal_cuts", ctr_bin_low_edge, ctr_bin_high_edge), Form("h_pT_80_counts_ctr%d_%d_w_signal_cuts", ctr_bin_low_edge, ctr_bin_high_edge), int(pT_bins_80.size()-1), pT_bins_80.data()}, 
                "pair_pt");
            hop_pT_80_counts_ctr_binned.push_back(h_pT_80_counts_op);

            auto h_pT_80_counts_ss = df_ctr_binned_ss_weight_w_signal_cuts.Histo1D(
                {Form("hss_pT_80_counts_ctr%d_%d_w_signal_cuts", ctr_bin_low_edge, ctr_bin_high_edge), Form("h_pT_80_counts_ctr%d_%d_w_signal_cuts", ctr_bin_low_edge, ctr_bin_high_edge), int(pT_bins_80.size()-1), pT_bins_80.data()}, 
                "pair_pt");
            hss_pT_80_counts_ctr_binned.push_back(h_pT_80_counts_ss);

            auto h_pT_150_op = df_ctr_binned_op_weight_w_signal_cuts.Histo1D(
                {Form("hop_pT_150_weighted_ctr%d_%d_w_signal_cuts", ctr_bin_low_edge, ctr_bin_high_edge), Form("h_pT_150_weighted_ctr%d_%d_w_signal_cuts", ctr_bin_low_edge, ctr_bin_high_edge), int(pT_bins_150.size()-1), pT_bins_150.data()}, 
                "pair_pt", "weight_for_RAA");
            hop_pT_150_ctr_binned.push_back(h_pT_150_op);

            auto h_pT_150_ss = df_ctr_binned_ss_weight_w_signal_cuts.Histo1D(
                {Form("hss_pT_150_weighted_ctr%d_%d_w_signal_cuts", ctr_bin_low_edge, ctr_bin_high_edge), Form("h_pT_150_weighted_ctr%d_%d_w_signal_cuts", ctr_bin_low_edge, ctr_bin_high_edge), int(pT_bins_150.size()-1), pT_bins_150.data()}, 
                "pair_pt", "weight_for_RAA");
            hss_pT_150_ctr_binned.push_back(h_pT_150_ss);

            auto h_pT_150_counts_op = df_ctr_binned_op_weight_w_signal_cuts.Histo1D(
                {Form("hop_pT_150_counts_ctr%d_%d_w_signal_cuts", ctr_bin_low_edge, ctr_bin_high_edge), Form("h_pT_150_counts_ctr%d_%d_w_signal_cuts", ctr_bin_low_edge, ctr_bin_high_edge), int(pT_bins_150.size()-1), pT_bins_150.data()}, 
                "pair_pt");
            hop_pT_150_counts_ctr_binned.push_back(h_pT_150_counts_op);

            auto h_pT_150_counts_ss = df_ctr_binned_ss_weight_w_signal_cuts.Histo1D(
                {Form("hss_pT_150_counts_ctr%d_%d_w_signal_cuts", ctr_bin_low_edge, ctr_bin_high_edge), Form("h_pT_150_counts_ctr%d_%d_w_signal_cuts", ctr_bin_low_edge, ctr_bin_high_edge), int(pT_bins_150.size()-1), pT_bins_150.data()}, 
                "pair_pt");
            hss_pT_150_counts_ctr_binned.push_back(h_pT_150_counts_ss);

            auto h_pT_200_op = df_ctr_binned_op_weight_w_signal_cuts.Histo1D(
                {Form("hop_pT_200_weighted_ctr%d_%d_w_signal_cuts", ctr_bin_low_edge, ctr_bin_high_edge), Form("h_pT_200_weighted_ctr%d_%d_w_signal_cuts", ctr_bin_low_edge, ctr_bin_high_edge), int(pT_bins_200.size()-1), pT_bins_200.data()}, 
                "pair_pt", "weight_for_RAA");
            hop_pT_200_ctr_binned.push_back(h_pT_200_op);

            auto h_pT_200_ss = df_ctr_binned_ss_weight_w_signal_cuts.Histo1D(
                {Form("hss_pT_200_weighted_ctr%d_%d_w_signal_cuts", ctr_bin_low_edge, ctr_bin_high_edge), Form("h_pT_200_weighted_ctr%d_%d_w_signal_cuts", ctr_bin_low_edge, ctr_bin_high_edge), int(pT_bins_200.size()-1), pT_bins_200.data()}, 
                "pair_pt", "weight_for_RAA");
            hss_pT_200_ctr_binned.push_back(h_pT_200_ss);

            auto h_pT_200_counts_op = df_ctr_binned_op_weight_w_signal_cuts.Histo1D(
                {Form("hop_pT_200_counts_ctr%d_%d_w_signal_cuts", ctr_bin_low_edge, ctr_bin_high_edge), Form("h_pT_200_counts_ctr%d_%d_w_signal_cuts", ctr_bin_low_edge, ctr_bin_high_edge), int(pT_bins_200.size()-1), pT_bins_200.data()}, 
                "pair_pt");
            hop_pT_200_counts_ctr_binned.push_back(h_pT_200_counts_op);

            auto h_pT_200_counts_ss = df_ctr_binned_ss_weight_w_signal_cuts.Histo1D(
                {Form("hss_pT_200_counts_ctr%d_%d_w_signal_cuts", ctr_bin_low_edge, ctr_bin_high_edge), Form("h_pT_200_counts_ctr%d_%d_w_signal_cuts", ctr_bin_low_edge, ctr_bin_high_edge), int(pT_bins_200.size()-1), pT_bins_200.data()}, 
                "pair_pt");
            hss_pT_200_counts_ctr_binned.push_back(h_pT_200_counts_ss);

		    auto th2dmodel_ss_crossx_signal_cuts = ROOT::RDF::TH2DModel(
		        Form("h2d_ss_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt_ctr%d_%d", ctr_bin_low_edge, ctr_bin_high_edge), ";p_{T} [GeV];#eta",
		        int(pT_bins_80.size() - 1), pT_bins_80.data(),
		        int(eta_bins.size() - 1), eta_bins.data()
		    );
		    auto th2dmodel_op_crossx_signal_cuts = ROOT::RDF::TH2DModel(
		        Form("h2d_op_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt_ctr%d_%d", ctr_bin_low_edge, ctr_bin_high_edge), ";p_{T} [GeV];#eta",
		        int(pT_bins_80.size() - 1), pT_bins_80.data(),
		        int(eta_bins.size() - 1), eta_bins.data()
		    );

	        auto h2d_ss_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt = df_ss_weight_w_signal_cuts.Histo2D(th2dmodel_ss_crossx_signal_cuts, "pair_pt", "pair_eta", "weight_for_RAA");
		    auto h2d_op_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt = df_op_weight_w_signal_cuts.Histo2D(th2dmodel_op_crossx_signal_cuts, "pair_pt", "pair_eta", "weight_for_RAA");
			h2d_ss_crossx_w_signal_cuts_ctr_binned.push_back(h2d_ss_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt);
			h2d_op_crossx_w_signal_cuts_ctr_binned.push_back(h2d_op_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt);

        }catch (const std::out_of_range& e){
            std::cerr << "WARNING: out_of_range caught when booking 1D histograms" << std::endl;
            std::cerr << "Skipping the current centrality bin: " << ctr_bin_low_edge << "-" << ctr_bin_high_edge << std::endl;
        }
    }

    // ----------------------------------------------------------------------------------------

    // Write results
    std::string outfile_name = dataset_base_dir + outfile_relative_path_truncated + ctr_suffix + ".root";
    cout << outfile_name << endl;
    TFile out(outfile_name.c_str(), "RECREATE");
    h3d_ss_crossx_w_signal_cuts_vs_centr_vs_pair_eta_vs_pair_pt->Write();
	h3d_op_crossx_w_signal_cuts_vs_centr_vs_pair_eta_vs_pair_pt->Write();
    
    for (int ictr = 0; ictr < nCtrBins; ictr++){
        try{
            hop_pT_80_ctr_binned.at(ictr)->Write();
            hss_pT_80_ctr_binned.at(ictr)->Write();
            hop_pT_80_counts_ctr_binned.at(ictr)->Write();
            hss_pT_80_counts_ctr_binned.at(ictr)->Write();
            hop_pT_150_ctr_binned.at(ictr)->Write();
            hss_pT_150_ctr_binned.at(ictr)->Write();
            hop_pT_150_counts_ctr_binned.at(ictr)->Write();
            hss_pT_150_counts_ctr_binned.at(ictr)->Write();
            hop_pT_200_ctr_binned.at(ictr)->Write();
            hss_pT_200_ctr_binned.at(ictr)->Write();
            hop_pT_200_counts_ctr_binned.at(ictr)->Write();
            hss_pT_200_counts_ctr_binned.at(ictr)->Write();
            h2d_ss_crossx_w_signal_cuts_ctr_binned.at(ictr)->Write();
            h2d_op_crossx_w_signal_cuts_ctr_binned.at(ictr)->Write();
        }catch (const std::out_of_range& e){
            std::cerr << "WARNING: out_of_range caught when writing output histograms" << std::endl;
            std::cerr << "Skipping the current centrality bin: " << ictr << std::endl;
        }
    }
    out.Close();
}


// -------------------------- running PbPb analysis -------------------------- 
void simplified_single_b_analysis_PbPb(){
    //Run PbPb analyses

    // SingleBAnalysisPbPb pbpb_run2 ("pbpb_run2/muon_pairs_pbpb_run2.root", "pbpb_run2/pbpb_run2_single_b_ana_hists");
    // pbpb_run2.crossx_factors_ctr_binned = PbPbBaseClass::crossx_factors_pbpb_run2_ctr_binned;
    // pbpb_run2.RunAnalysis();

    // SingleBAnalysisPbPb pbpb_2023 ("pbpb_2023/muon_pairs_pbpb_2023_mu4_mu4noL1_no_res_cut.root", "pbpb_2023/pbpb_2023_single_b_ana_hists_mu4_mu4noL1_no_res_cut");
    // pbpb_2023.crossx_factors_ctr_binned = PbPbBaseClass::crossx_factors_pbpb_2023_ctr_binned;
    // pbpb_2023.RunAnalysis();

    SingleBAnalysisPbPb pbpb_2024 ("pbpb_2024/muon_pairs_pbpb_2024_single_mu4.root", "pbpb_2024/pbpb_2024_single_b_ana_hists_single_mu4");
    pbpb_2024.crossx_factors_ctr_binned = PbPbBaseClass::crossx_factors_pbpb_2024_ctr_binned;
    // pbpb_2024.BinningPrinting();
    pbpb_2024.RunAnalysis();
}

