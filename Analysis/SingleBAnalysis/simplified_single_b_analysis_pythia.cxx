#include "SingleBAnalysisBase.cxx"

// -------------------------- pythia analysis class --------------------------

class SingleBAnalysisPythia : public SingleBAnalysisBase{
public:
    SingleBAnalysisPythia(){
        Initialize();
    }
    ~SingleBAnalysisPythia(){}
    void RunAnalysis();
};


// -------------------------- pythia-class member functions --------------------------

void SingleBAnalysisPythia::RunAnalysis() {

    // Create RDataFrame and apply initial cuts
    muon_pair_input_file_name = dataset_base_dir + "pythia/muon_pairs_pythia_combined_no_data_resonance_cuts.root";
    RDataFrame df("muon_pair_tree_sign2", muon_pair_input_file_name.c_str());
    auto rdf_with_signal_cuts = df.Filter("minv > 1.08 && minv < 2.9 && pair_pt > 8 && !data_resonance_or_reso_contam_tagged_old");

    // Create truth-separated datasets
    auto rdf_truth_signal_no_cuts = df.Filter("muon_pair_flavor_category == 2"); // truth signal with no data-like cuts
    auto rdf_truth_signal = rdf_with_signal_cuts.Filter("muon_pair_flavor_category == 2");
    auto rdf_truth_bkg = rdf_with_signal_cuts.Filter("muon_pair_flavor_category != 2");

    // ----------------------------------------------------------------------------------------

    // Create 2D histogram model for metrics calculation
    auto th2dmodel_signal_cuts = ROOT::RDF::TH2DModel(
        "h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts", ";p_{T} [GeV];#eta",
        // "h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts", ";p_{T} [GeV];#eta;d^2#sigma/dp_{T}^{pair}d#eta^{pair}[pb]",
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
    TFile out("/Users/yuhanguo/Documents/physics/heavy-ion/dimuon/datasets/pythia/pythia_single_b_ana_hists.root", "RECREATE");
    h_pT_80->Write();
    h_pT_200->Write();
    h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts->SetTitle("");
    h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts->SetStats(0);
    h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts->Scale(pow(10,6),"width");
    h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts->Write();
    h2d_truth_sig_all_no_cuts->Write();
    h2d_truth_sig_w_cuts->Write();
    h2d_truth_bkg_w_cuts->Write();
    h_acceptance.Write();
    h_significance.Write();
    h_purity.Write();
    out.Close();
}


// -------------------------- running pythia analysis -------------------------- 

void simplified_single_b_analysis_pythia(){
    // Run Pythia analysis
    SingleBAnalysisPythia py;
    py.RunAnalysis();
}
