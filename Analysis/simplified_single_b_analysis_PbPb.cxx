#include "SingleBAnalysisBase.cxx"

// -------------------------- PbPb analysis class --------------------------

class SingleBAnalysisPbPb : public SingleBAnalysisBase{
private:
    std::string infile_relative_path;
    std::string outfile_relative_path_truncated;
    std::vector<int> ctr_bin_edges;
    int nCtrBins;

	double CalculateWeightForRAA(int centrality, double weight);

public:
    // std::vector<std::vector<std::string>> ctr_bins;
    std::string ctr_binning_verion = "default"; // centrality binning version used in current analysis (if there are several)
    std::map<std::string, std::vector<int>> ctr_binning_map; // map from binning versions to the binnings
    std::map<std::string, std::string> ctr_binning_file_suffix_map; // map from binning versions to the suffices
    std::vector<double> crossx_factors_ctr_binned;

    SingleBAnalysisPbPb(std::string infile_relative_path_input, std::string outfile_relative_path_input){
        infile_relative_path = infile_relative_path_input;
        outfile_relative_path_truncated = outfile_relative_path_input;
        Initialize(); // run the PbPb version of Initialize()
    }
    ~SingleBAnalysisPbPb(){}

    void RunAnalysis();
    void Initialize();
};


// -------------------------- PbPb-class member functions --------------------------

void SingleBAnalysisPbPb::Initialize(){
    SingleBAnalysisBase::Initialize();

    ctr_binning_map["default"] = {0,5,10,20,30,50,80};
    ctr_binning_file_suffix_map["default"] = "";
}


double SingleBAnalysisPbPb::CalculateWeightForRAA(int centrality, double weight){
    for (int ictr = 0; ictr < nCtrBins; ictr++){

    	int ctr_bin_low_edge, ctr_bin_high_edge;
    	double crossx_factor;

    	try{
    		ctr_bin_low_edge = ctr_bin_edges.at(ictr);
			ctr_bin_high_edge = ctr_bin_edges.at(ictr + 1);
            crossx_factor = crossx_factors_ctr_binned.at(ictr);
        }catch (const std::out_of_range& e){
            std::cerr << "WARNING: out_of_range caught when accessing the " << ictr << "-th element of ctr_bin_edges and crossx_factors_ctr_binned elements" << std::endl;
            std::cerr << "ctr_bin_edges has size " << ctr_bin_edges.size() << std::endl;
            std::cerr << "crossx_factors_ctr_binned has size " << crossx_factors_ctr_binned.size() << std::endl;
            std::cerr << "Assigning ZERO weight to current event: will result in WRONG total crossx" << std::endl;
            return 0;
        }

        if (centrality >= ctr_bin_low_edge && centrality < ctr_bin_high_edge){
        	return weight * crossx_factor;
        }
    }
    return 0.; // assign weight = 0 to event if out of centrality-range of interest
}


void SingleBAnalysisPbPb::RunAnalysis(){

    // Create RDataFrame for the diff Centrality groups & signs, correctly weighted
    muon_pair_input_file_name = dataset_base_dir + infile_relative_path;

    RDataFrame df_ss("muon_pair_tree_sign1", muon_pair_input_file_name.c_str());    
    RDataFrame df_op("muon_pair_tree_sign2", muon_pair_input_file_name.c_str());    

    std::vector<ROOT::RDF::RNode> df_pbpb_run2_ctr_binned_list_ss_w_weight_w_signal_cuts;
    std::vector<ROOT::RDF::RNode> df_pbpb_run2_ctr_binned_list_op_w_weight_w_signal_cuts;

	ctr_bin_edges = ctr_binning_map[ctr_binning_verion];
	std::vector<double> ctr_bin_edges_double(ctr_bin_edges.begin(), ctr_bin_edges.end());

    if (crossx_factors_ctr_binned.size() != ctr_bin_edges.size() - 1){
        throw std::runtime_error("crossx_factors_ctr_binned MUST equal ctr_bin_edges size - 1");
    }

    nCtrBins = ctr_bin_edges.size() - 1;

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
                {Form("hop_pT_80_ctr%d_%d_w_signal_cuts", ctr_bin_low_edge, ctr_bin_high_edge), Form("h_pT_80_ctr%d_%d_w_signal_cuts", ctr_bin_low_edge, ctr_bin_high_edge), int(pT_bins_80.size()-1), pT_bins_80.data()}, 
                "pair_pt", "weight_for_RAA");
            hop_pT_80_ctr_binned.push_back(h_pT_80_op);

            auto h_pT_80_ss = df_ctr_binned_ss_weight_w_signal_cuts.Histo1D(
                {Form("hss_pT_80_ctr%d_%d_w_signal_cuts", ctr_bin_low_edge, ctr_bin_high_edge), Form("h_pT_80_ctr%d_%d_w_signal_cuts", ctr_bin_low_edge, ctr_bin_high_edge), int(pT_bins_80.size()-1), pT_bins_80.data()}, 
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

            auto h_pT_200_op = df_ctr_binned_op_weight_w_signal_cuts.Histo1D(
                {Form("hop_pT_200_ctr%d_%d_w_signal_cuts", ctr_bin_low_edge, ctr_bin_high_edge), Form("h_pT_200_ctr%d_%d_w_signal_cuts", ctr_bin_low_edge, ctr_bin_high_edge), int(pT_bins_200.size()-1), pT_bins_200.data()}, 
                "pair_pt", "weight_for_RAA");
            hop_pT_200_ctr_binned.push_back(h_pT_200_op);

            auto h_pT_200_ss = df_ctr_binned_ss_weight_w_signal_cuts.Histo1D(
                {Form("hss_pT_200_ctr%d_%d_w_signal_cuts", ctr_bin_low_edge, ctr_bin_high_edge), Form("h_pT_200_ctr%d_%d_w_signal_cuts", ctr_bin_low_edge, ctr_bin_high_edge), int(pT_bins_200.size()-1), pT_bins_200.data()}, 
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
    std::string ctr_binning_version_suffix = ctr_binning_file_suffix_map[ctr_binning_verion];
    std::string outfile_name = dataset_base_dir + outfile_relative_path_truncated + ctr_binning_version_suffix + ".root";
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

    // 67.6 = sigma_{inel} in unit of mb = <TAA> / <N_{coll}>
    std::vector<double> crossx_factors_pbpb_run2_ctr_binned = {
        1./(0.05 * 7.67 * 1000. * 26.2339 * (0.436882 + 1.3803) / 1000.), // 1/[crossx_{incl AA} [unit: mb] * TAA [unit: mb^{-1}] * (L_{int15} + L_{int18}) (unit: pb^{-1})] 0-5%
        1./(0.05 * 7.67 * 1000. * 20.4707 * (0.436882 + 1.3803) / 1000.), // 1/[crossx_{incl AA} [unit: mb] * TAA [unit: mb^{-1}] * (L_{int15} + L_{int18}) (unit: pb^{-1})] 5-10%
        1./(0.1 * 7.67 * 1000. * 14.3345 * (0.436882 + 1.3803) / 1000.), // 1/[crossx_{incl AA} [unit: mb] * TAA [unit: mb^{-1}] * (L_{int15} + L_{int18}) (unit: pb^{-1})] 10-20%
        1./(0.1 * 7.67 * 1000. * 8.63844 * (0.436882 + 1.3803) / 1000.), // 1/[crossx_{incl AA} [unit: mb] * TAA [unit: mb^{-1}] * (L_{int15} + L_{int18}) (unit: pb^{-1})] 20-30%
        1./(0.2 * 7.67 * 1000. * 3.79023 * (0.436882 + 1.3803) / 1000.), // 1/[crossx_{incl AA} [unit: mb] * TAA [unit: mb^{-1}] * (L_{int15} + L_{int18}) (unit: pb^{-1})] 30-50%
        1./(0.3 * 7.67 * 1000. * 0.689695 * (0.436882 + 1.3803) / 1000.) // 1/[crossx_{incl AA} [unit: mb] * TAA [unit: mb^{-1}] * (L_{int15} + L_{int18}) (unit: pb^{-1})] 50-80%
    };

    std::vector<double> crossx_factors_pbpb_2023_ctr_binned = {
        1./(0.05 * 7.8 * 1000. * 26.1428 * 1.3896 / 1000.), // 1/[crossx_{incl AA} [unit: mb] * TAA [unit: mb^{-1}] * L_{int} (unit: pb^{-1})] 0-5%
        1./(0.05 * 7.8 * 1000. * 20.3241 * 1.3896 / 1000.), // 1/[crossx_{incl AA} [unit: mb] * TAA [unit: mb^{-1}] * L_{int} (unit: pb^{-1})] 5-10%
        1./(0.1 * 7.8 * 1000. * 14.0502 * 1.3896 / 1000.), // 1/[crossx_{incl AA} [unit: mb] * TAA [unit: mb^{-1}] * L_{int} (unit: pb^{-1})] 10-20%
        1./(0.1 * 7.8 * 1000. * 8.5074 * 1.3896 / 1000.), // 1/[crossx_{incl AA} [unit: mb] * TAA [unit: mb^{-1}] * L_{int} (unit: pb^{-1})] 20-30%
        1./(0.2 * 7.8 * 1000. * 3.7733 * 1.3896 / 1000.), // 1/[crossx_{incl AA} [unit: mb] * TAA [unit: mb^{-1}] * L_{int} (unit: pb^{-1})] 30-50%
        1./(0.3 * 7.8 * 1000. * 0.6716 * 1.3896 / 1000.) // 1/[crossx_{incl AA} [unit: mb] * TAA [unit: mb^{-1}] * L_{int} (unit: pb^{-1})] 50-80%
    };

    std::vector<double> crossx_factors_pbpb_2024_ctr_binned = {
        1./(0.05 * 7.8 * 1000. * 26.1428 * 1.5411 / 1000.), // 1/[crossx_{incl AA} [unit: mb] * TAA [unit: mb^{-1}] * L_{int} (unit: pb^{-1})] 0-5%
        1./(0.05 * 7.8 * 1000. * 20.3241 * 1.5411 / 1000.), // 1/[crossx_{incl AA} [unit: mb] * TAA [unit: mb^{-1}] * L_{int} (unit: pb^{-1})] 5-10%
        1./(0.1 * 7.8 * 1000. * 14.0502 * 1.5411 / 1000.), // 1/[crossx_{incl AA} [unit: mb] * TAA [unit: mb^{-1}] * L_{int} (unit: pb^{-1})] 10-20%
        1./(0.1 * 7.8 * 1000. * 8.5074 * 1.5411 / 1000.), // 1/[crossx_{incl AA} [unit: mb] * TAA [unit: mb^{-1}] * L_{int} (unit: pb^{-1})] 20-30%
        1./(0.2 * 7.8 * 1000. * 3.7733 * 1.5411 / 1000.), // 1/[crossx_{incl AA} [unit: mb] * TAA [unit: mb^{-1}] * L_{int} (unit: pb^{-1})] 30-50%
        1./(0.3 * 7.8 * 1000. * 0.6716 * 1.5411 / 1000.) // 1/[crossx_{incl AA} [unit: mb] * TAA [unit: mb^{-1}] * L_{int} (unit: pb^{-1})] 50-80%
    };


    // SingleBAnalysisPbPb pbpb_run2 ("pbpb_run2/muon_pairs_pbpb_run2.root", "pbpb_run2/pbpb_run2_single_b_ana_hists");
    // pbpb_run2.crossx_factors_ctr_binned = crossx_factors_pbpb_run2_ctr_binned;
    // pbpb_run2.BinningPrinting();
    // pbpb_run2.RunAnalysis();

    // SingleBAnalysisPbPb pbpb_2023 ("pbpb_2023/muon_pairs_pbpb_2023.root", "pbpb_2023/pbpb_2023_single_b_ana_hists");
    // pbpb_2023.crossx_factors_ctr_binned = crossx_factors_pbpb_2023_ctr_binned;
    // pbpb_2023.RunAnalysis();

    SingleBAnalysisPbPb pbpb_2024 ("pbpb_2024/muon_pairs_pbpb_2024.root", "pbpb_2024/pbpb_2024_single_b_ana_hists");
    pbpb_2024.crossx_factors_ctr_binned = crossx_factors_pbpb_2024_ctr_binned;
    pbpb_2024.RunAnalysis();
}

