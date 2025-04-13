#include <cmath>
#include <iostream>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TH2D.h>
#include <TFile.h>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <TChain.h>


using namespace ROOT;
using namespace ROOT::VecOps;

// -------------------------- analysis base class --------------------------
class SingleBAnalysisBase{
protected:
    std::string dataset_base_dir = "/Users/yuhanguo/Documents/physics/heavy-ion/dimuon/datasets/";
    std::string muon_pair_input_file_name;

    constexpr static double pi = 3.14159265358979323846;

    // Define binning
    std::vector<double> pT_bins_80;
    std::vector<double> pT_bins_200;
    constexpr static int n_eta_bins = 10;
    std::vector<double> eta_bins;
    constexpr static double eta_min = -2.5, eta_max = 2.5;

    virtual void Initialize();
    void fillLogBinningArray(std::vector<double>& bins, int nBins, double low, double high);

public:
    SingleBAnalysisBase(){
        Initialize();
    }
    ~SingleBAnalysisBase(){}

    virtual void RunAnalysis() = 0;
};

// -------------------------- pythia analysis class --------------------------

class SingleBAnalysisPythia : public SingleBAnalysisBase{
public:
    SingleBAnalysisPythia(){}
    ~SingleBAnalysisPythia(){}
    void RunAnalysis();
};

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
    }

    void RunAnalysis();
};

// -------------------------- PbPb analysis class --------------------------

class SingleBAnalysisPbPb : public SingleBAnalysisBase{
private:
    std::string infile_relative_path;
    std::string outfile_relative_path_truncated;

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



// -------------------------- base-class member functions --------------------------

void SingleBAnalysisBase::fillLogBinningArray(std::vector<double>& bins, int nBins, double low, double high) {
    bins.clear();

    double logLow = std::log10(low);
    double logHigh = std::log10(high);
    double logStep = (logHigh - logLow) / nBins;

    for (int i = 0; i <= nBins; ++i) {
        bins.push_back(std::pow(10, logLow + i * logStep));
    }
}

void SingleBAnalysisBase::Initialize(){
    // // Create logarithmic pair-pT bin edges
    fillLogBinningArray(pT_bins_80, 12, 8.0, 80.0);  // 10 log bins from 8 to 200 GeV
    fillLogBinningArray(pT_bins_200, 18, 8.0, 200.0);  // 10 log bins from 8 to 200 GeV

    // Create eta bin edges
    for (int i = 0; i <= n_eta_bins; ++i)
        eta_bins.push_back(eta_min + (eta_max - eta_min) * i / n_eta_bins);
}

// -------------------------- pythia-class member functions --------------------------

void SingleBAnalysisPythia::RunAnalysis() {

    // Create RDataFrame and apply initial cuts
    muon_pair_input_file_name = dataset_base_dir + "pythia/muon_pairs_pythia_combined_no_data_resonance_cuts.root";
    RDataFrame df("muon_pair_tree_sign2", muon_pair_input_file_name.c_str());
    auto rdf_with_signal_cuts = df.Filter("minv > 1.08 && minv < 2.9 && pair_pt > 8");

    // Create truth-separated datasets
    auto rdf_truth_signal_no_cuts = df.Filter("muon_pair_flavor_category == 2"); // truth signal with no data-like cuts
    auto rdf_truth_signal = rdf_with_signal_cuts.Filter("muon_pair_flavor_category == 2");
    auto rdf_truth_bkg = rdf_with_signal_cuts.Filter("muon_pair_flavor_category != 2");

    // ----------------------------------------------------------------------------------------

    // Create 2D histogram model for metrics calculation
    auto th2dmodel_signal_cuts = ROOT::RDF::TH2DModel(
        "h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts", "Signal (data-like cuts);p_{T} [GeV];#eta",
        int(pT_bins_80.size() - 1), pT_bins_80.data(),
        int(eta_bins.size() - 1), eta_bins.data()
    );

    auto model_truth_sig_all_no_cuts = ROOT::RDF::TH2DModel(
        "h2d_truth_sig_all_no_cuts", "Truth Signal NO data-like cuts;p_{T} [GeV];#eta",
        int(pT_bins_80.size() - 1), pT_bins_80.data(),
        int(eta_bins.size() - 1), eta_bins.data()
    );

    auto model_truth_sig_w_cuts = ROOT::RDF::TH2DModel(
        "h2d_truth_sig_w_cuts", "Truth Signal with data-like cuts;p_{T} [GeV];#eta",
        int(pT_bins_80.size() - 1), pT_bins_80.data(),
        int(eta_bins.size() - 1), eta_bins.data()
    );

    auto model_truth_bkg_w_cuts = ROOT::RDF::TH2DModel(
        "h2d_truth_bkg_w_cuts", "Truth Background with data-like cuts;p_{T} [GeV];#eta",
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
    auto h2d_truth_sig_all_no_cuts = rdf_truth_signal_no_cuts.Histo2D(model_truth_sig_all_no_cuts, "pair_pt", "pair_eta", "weight");
    auto h2d_truth_sig_w_cuts = rdf_truth_signal.Histo2D(model_truth_sig_w_cuts, "pair_pt", "pair_eta", "weight");
    auto h2d_truth_bkg_w_cuts = rdf_truth_bkg.Histo2D(model_truth_bkg_w_cuts, "pair_pt", "pair_eta", "weight");

    // ----------------------------------------------------------------------------------------

    // Create output 2D histograms for truth signal-selection metrics
    TH2D h_acceptance("acceptance", "Signal Acceptance;p_{T} [GeV];#eta", 
                     int(pT_bins_80.size()-1), pT_bins_80.data(),
                     n_eta_bins, eta_min, eta_max);
    TH2D h_significance("significance", "Significance S/#sqrt{B};p_{T} [GeV];#eta",
                       int(pT_bins_80.size()-1), pT_bins_80.data(),
                       n_eta_bins, eta_min, eta_max);
    TH2D h_purity("purity", "Signal Purity S/(S+B);p_{T} [GeV];#eta",
                 int(pT_bins_80.size()-1), pT_bins_80.data(),
                 n_eta_bins, eta_min, eta_max);


    // Calculate metrics for each bin
    for (int xbin = 1; xbin <= h_acceptance.GetNbinsX(); ++xbin) {
        for (int ybin = 1; ybin <= h_acceptance.GetNbinsY(); ++ybin) {
            const double S_all_no_cuts = h2d_truth_sig_all_no_cuts->GetBinContent(xbin, ybin); // crossx (truth signal) in (pair pT, pair eta) bin with NO data-like cuts
            const double S = h2d_truth_sig_w_cuts->GetBinContent(xbin, ybin); // crossx (truth signal) in (pair pT, pair eta) bin with data-like cuts
            const double B = h2d_truth_bkg_w_cuts->GetBinContent(xbin, ybin); // crossx (truth background) in (pair pT, pair eta) bin with data-like cuts

            // Calculate metrics
            const double acceptance = (S_all_no_cuts > 0) ? S / S_all_no_cuts : 0;
            const double significance = (B > 0) ? S / std::sqrt(B) : 0;
            const double purity = (S + B > 0) ? S / (S + B) : 0;

            // Fill result histograms
            h_acceptance.SetBinContent(xbin, ybin, acceptance);
            h_significance.SetBinContent(xbin, ybin, significance);
            h_purity.SetBinContent(xbin, ybin, purity);
        }
    }

    // ----------------------------------------------------------------------------------------

    // Write results
    TFile out("/Users/yuhanguo/Documents/physics/heavy-ion/dimuon/datasets/pythia/pythia_single_b_histograms.root", "RECREATE");
    h_pT_80->Write();
    h_pT_200->Write();
    h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts->Write();
    h2d_truth_sig_all_no_cuts->Write();
    h2d_truth_sig_w_cuts->Write();
    h2d_truth_bkg_w_cuts->Write();
    h_acceptance.Write();
    h_significance.Write();
    h_purity.Write();
    out.Close();
}

// -------------------------- PP-class member functions --------------------------

void SingleBAnalysisPP::RunAnalysis(){

    // ----------------------------------------------------------------------------------------

    // Create RDataFrame and apply initial cuts
    muon_pair_input_file_name = dataset_base_dir + infile_relative_path;
    RDataFrame df("muon_pair_tree_sign2", muon_pair_input_file_name.c_str());
    auto rdf_with_signal_cuts = df.Filter("minv > 1.08 && minv < 2.9 && pair_pt > 8");


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
    h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts->Write();
    out.Close();
}

// -------------------------- PbPb-class member functions --------------------------

void SingleBAnalysisPbPb::Initialize(){
    SingleBAnalysisBase::Initialize();

    ctr_binning_map["default"] = {0,5,10,20,30,50,80};
    ctr_binning_file_suffix_map["default"] = "";
}

void SingleBAnalysisPbPb::RunAnalysis(){

    // Create RDataFrame for the diff Centrality groups & signs, correctly weighted
    muon_pair_input_file_name = dataset_base_dir + infile_relative_path;

    RDataFrame df_ss("muon_pair_tree_sign1", muon_pair_input_file_name.c_str());    
    RDataFrame df_op("muon_pair_tree_sign2", muon_pair_input_file_name.c_str());    

    std::vector<ROOT::RDF::RNode> df_pbpb_run2_ctr_binned_list_ss_w_weight;
    std::vector<ROOT::RDF::RNode> df_pbpb_run2_ctr_binned_list_op_w_weight;
    std::vector<ROOT::RDF::RNode> df_pbpb_run2_ctr_binned_list_ss_w_weight_w_signal_cuts;
    std::vector<ROOT::RDF::RNode> df_pbpb_run2_ctr_binned_list_op_w_weight_w_signal_cuts;

    std::vector<int> ctr_bin_edges = ctr_binning_map[ctr_binning_verion];

    for (int ictr = 0; ictr < ctr_bin_edges.size() - 1; ictr++){

        int ctr_bin_low_edge = ctr_bin_edges.at(ictr);
        int ctr_bin_high_edge = ctr_bin_edges.at(ictr + 1);

        double crossx_factor;

        try{
            crossx_factor = crossx_factors_ctr_binned.at(ictr);
        }catch (const std::out_of_range& e){
            std::cerr << "WARNING: out_of_range caught when accessing crossx_factors_ctr_binned elements" << std::endl;
            std::cerr << "crossx_factors_ctr_binned has size " << crossx_factors_ctr_binned.size() << std::endl;
            std::cerr << "Skipping the current centrality bin: " << ctr_bin_low_edge << "-" << ctr_bin_high_edge << std::endl;
        }

        auto df_ctr_binned_ss = df_ss.Filter(
            [ctr_bin_low_edge, ctr_bin_high_edge](int centrality) {return centrality >= ctr_bin_low_edge && centrality < ctr_bin_high_edge;},
            {"avg_centrality"});
        auto df_ctr_binned_ss_weight = df_ctr_binned_ss.Define("weight_crossx_to_TAA", 
            [crossx_factor](double weight) { return weight * crossx_factor; }, {"weight"});

        auto df_ctr_binned_ss_weight_w_signal_cuts = df_ctr_binned_ss_weight.Filter("minv > 1.08 && minv < 2.9 && pair_pt > 8");

        auto df_ctr_binned_op = df_op.Filter(
            [ctr_bin_low_edge, ctr_bin_high_edge](int centrality) {return centrality >= ctr_bin_low_edge && centrality < ctr_bin_high_edge;},
            {"avg_centrality"});
        auto df_ctr_binned_op_weight = df_ctr_binned_op.Define("weight_crossx_to_TAA", 
            [crossx_factor](double weight) { return weight * crossx_factor; }, {"weight"});

        auto df_ctr_binned_op_weight_w_signal_cuts = df_ctr_binned_op_weight.Filter("minv > 1.08 && minv < 2.9 && pair_pt > 8");

        df_pbpb_run2_ctr_binned_list_ss_w_weight.push_back(df_ctr_binned_ss_weight);
        df_pbpb_run2_ctr_binned_list_ss_w_weight_w_signal_cuts.push_back(df_ctr_binned_ss_weight_w_signal_cuts);
        df_pbpb_run2_ctr_binned_list_op_w_weight.push_back(df_ctr_binned_op_weight);
        df_pbpb_run2_ctr_binned_list_op_w_weight_w_signal_cuts.push_back(df_ctr_binned_op_weight_w_signal_cuts);
    }

    // ----------------------------------------------------------------------------------------

    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_pT_80_ctr_binned;
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_pT_80_counts_ctr_binned;
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_pT_200_ctr_binned;
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_pT_200_counts_ctr_binned;
    // std::vector<ROOT::RDF::RResultPtr<TH2D>> h2d_crossx_pair_pt_pair_eta_binned_op_w_signal_cuts_ctr_binned;
    // std::vector<ROOT::RDF::RResultPtr<TH2D>> h2d_crossx_pair_pt_pair_eta_binned_ss_w_signal_cuts_ctr_binned;

    for (int ictr = 0; ictr < ctr_bin_edges.size() - 1; ictr++){

        int ctr_bin_low_edge = ctr_bin_edges.at(ictr);
        int ctr_bin_high_edge = ctr_bin_edges.at(ictr + 1);

        try{
            auto df_ctr_binned_op_weight_w_signal_cuts = df_pbpb_run2_ctr_binned_list_op_w_weight_w_signal_cuts.at(ictr);
            auto df_ctr_binned_ss_weight_w_signal_cuts = df_pbpb_run2_ctr_binned_list_ss_w_weight_w_signal_cuts.at(ictr);
    
            // Book 1D histograms
            auto h_pT_80 = df_ctr_binned_op_weight_w_signal_cuts.Histo1D(
                {Form("h_pT_80_ctr%d_%d_w_signal_cuts", ctr_bin_low_edge, ctr_bin_high_edge), Form("h_pT_80_ctr%d_%d_w_signal_cuts", ctr_bin_low_edge, ctr_bin_high_edge), int(pT_bins_80.size()-1), pT_bins_80.data()}, 
                "pair_pt", "weight_crossx_to_TAA");
            h_pT_80_ctr_binned.push_back(h_pT_80);

            auto h_pT_80_counts = df_ctr_binned_op_weight_w_signal_cuts.Histo1D(
                {Form("h_pT_80_counts_ctr%d_%d_w_signal_cuts", ctr_bin_low_edge, ctr_bin_high_edge), Form("h_pT_80_counts_ctr%d_%d_w_signal_cuts", ctr_bin_low_edge, ctr_bin_high_edge), int(pT_bins_80.size()-1), pT_bins_80.data()}, 
                "pair_pt");
            h_pT_80_counts_ctr_binned.push_back(h_pT_80_counts);

            auto h_pT_200 = df_ctr_binned_op_weight_w_signal_cuts.Histo1D(
                {Form("h_pT_200_ctr%d_%d_w_signal_cuts", ctr_bin_low_edge, ctr_bin_high_edge), Form("h_pT_200_ctr%d_%d_w_signal_cuts", ctr_bin_low_edge, ctr_bin_high_edge), int(pT_bins_200.size()-1), pT_bins_200.data()}, 
                "pair_pt", "weight_crossx_to_TAA");
            h_pT_200_ctr_binned.push_back(h_pT_200);

            auto h_pT_200_counts = df_ctr_binned_op_weight_w_signal_cuts.Histo1D(
                {Form("h_pT_200_counts_ctr%d_%d_w_signal_cuts", ctr_bin_low_edge, ctr_bin_high_edge), Form("h_pT_200_counts_ctr%d_%d_w_signal_cuts", ctr_bin_low_edge, ctr_bin_high_edge), int(pT_bins_200.size()-1), pT_bins_200.data()}, 
                "pair_pt");
            h_pT_200_counts_ctr_binned.push_back(h_pT_200_counts);

            // Create centrality-binned 2D histogram model for metrics calculation
            auto th2dmodel_signal_cuts = ROOT::RDF::TH2DModel(
                Form("h2d_crossx_pair_pt_pair_eta_binned_ctr%d_%d_w_signal_cuts", ctr_bin_low_edge, ctr_bin_high_edge), "Signal (data-like cuts);p_{T} [GeV];#eta",
                int(pT_bins_80.size() - 1), pT_bins_80.data(),
                int(eta_bins.size() - 1), eta_bins.data()
            );

            // Book 2D histograms for metric calculations
            // auto h2d_crossx_pair_pt_pair_eta_binned_op_w_signal_cuts = df_ctr_binned_op_weight_w_signal_cuts.Histo2D(th2dmodel_signal_cuts, "pair_pt", "pair_eta", "weight_crossx_to_TAA");
            // auto h2d_crossx_pair_pt_pair_eta_binned_ss_w_signal_cuts = df_ctr_binned_ss_weight_w_signal_cuts.Histo2D(th2dmodel_signal_cuts, "pair_pt", "pair_eta", "weight_crossx_to_TAA");
            // h2d_crossx_pair_pt_pair_eta_binned_op_w_signal_cuts_ctr_binned.push_back(h2d_crossx_pair_pt_pair_eta_binned_op_w_signal_cuts);
            // h2d_crossx_pair_pt_pair_eta_binned_ss_w_signal_cuts_ctr_binned.push_back(h2d_crossx_pair_pt_pair_eta_binned_ss_w_signal_cuts);
        }catch (const std::out_of_range& e){
            std::cerr << "WARNING: out_of_range caught when booking 1D & 2D histograms" << std::endl;
            std::cerr << "Skipping the current centrality bin: " << ctr_bin_low_edge << "-" << ctr_bin_high_edge << std::endl;
        }
    }

    // ----------------------------------------------------------------------------------------

    // Write results
    std::string ctr_binning_version_suffix = ctr_binning_file_suffix_map[ctr_binning_verion];
    std::string outfile_name = dataset_base_dir + outfile_relative_path_truncated + ctr_binning_version_suffix + ".root";
    cout << outfile_name << endl;
    TFile out(outfile_name.c_str(), "RECREATE");
    for (int ictr = 0; ictr < ctr_bin_edges.size() - 1; ictr++){
        try{
            h_pT_80_ctr_binned.at(ictr)->Write();
            h_pT_80_counts_ctr_binned.at(ictr)->Write();
            h_pT_200_ctr_binned.at(ictr)->Write();
            h_pT_200_counts_ctr_binned.at(ictr)->Write();
            // h2d_crossx_pair_pt_pair_eta_binned_op_w_signal_cuts_ctr_binned.at(ictr)->Write();
            // h2d_crossx_pair_pt_pair_eta_binned_ss_w_signal_cuts_ctr_binned.at(ictr)->Write();
        }catch (const std::out_of_range& e){
            std::cerr << "WARNING: out_of_range caught when writing output histograms" << std::endl;
            std::cerr << "Skipping the current centrality bin: " << ictr << std::endl;
        }
    }
    out.Close();
}


// -------------------------- functions running analyses -------------------------- 

void simplified_single_b_pythia_analysis(){
    // Run Pythia analysis
    SingleBAnalysisPythia py;
    py.RunAnalysis();
}

void simplified_single_b_pp_analysis(){
    //Run PP analyses

    double crossx_factor_pp_run2 = 1/256.8;
    double crossx_factor_pp_24_2mu4 = 1/410.815;
    double crossx_factor_pp_24_mu4mu4noL1 = 1/113.999;

    SingleBAnalysisPP pprun2("pp_run2/muon_pairs_pp2017.root", "pp_run2/pp_run2_single_b_histograms");
    pprun2.crossx_factor = crossx_factor_pp_run2;
    pprun2.RunAnalysis();

    SingleBAnalysisPP pp2024_2mu4("pp_2024/muon_pairs_pp_2024_mu4_mu4noL1.root", "pp_2024/pp_2024_single_b_histograms");
    pp2024_2mu4.crossx_factor = crossx_factor_pp_24_2mu4;
    pp2024_2mu4.RunAnalysis();

    SingleBAnalysisPP pp2024_mu4mu4noL1("pp_2024/muon_pairs_pp_2024_2mu4.root", "pp_2024/pp_2024_single_b_histograms");
    pp2024_mu4mu4noL1.crossx_factor = crossx_factor_pp_24_mu4mu4noL1;
    pp2024_mu4mu4noL1.RunAnalysis();
}

void simplified_single_b_PbPb_analysis(){
    //Run PbPb analyses

    std::vector<double> crossx_factors_pbpb_run2_ctr_binned = {
        0.05128, 0.06536, 0.04630, 0.07602, 0.08503, 0.30441
    };

    SingleBAnalysisPbPb pbpb_run2 ("pbpb_run2/muon_pairs_small_ctr_intvls.root", "pbpb_run2/pbpb_run2_single_b_histograms");
    pbpb_run2.crossx_factors_ctr_binned = crossx_factors_pbpb_run2_ctr_binned;
    pbpb_run2.RunAnalysis();

    // SingleBAnalysisPbPb pbpb_2023 ("pbpb_2023/muon_pairs_small_ctr_intvls.root", "pbpb_2023/pbpb_2023_single_b_histograms");
    // // pbpb_2023.crossx_factors_ctr_binned = crossx_factors_pbpb_run2_ctr_binned;
    // pbpb_2023.RunAnalysis();

    // SingleBAnalysisPbPb pbpb_2024 ("pbpb_2024/muon_pairs_small_ctr_intvls.root", "pbpb_2024/pbpb_2024_single_b_histograms");
    // // pbpb_2024.crossx_factors_ctr_binned = crossx_factors_pbpb_run2_ctr_binned;
    // pbpb_2024.RunAnalysis();


}
void simplified_single_b_analysis(){
    // simplified_single_b_pythia_analysis();
    // simplified_single_b_pp_analysis();
    simplified_single_b_PbPb_analysis();
}