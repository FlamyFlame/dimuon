#pragma once

class PbPbBaseClass
{
protected:
	// -------- protected class variables --------

    std::vector<int> ctr_bin_edges;
    std::vector<double> ctr_bin_edges_double;
    std::vector<std::string> ctr_bins;
    int nCtrBins;

    std::map<std::string, std::vector<int>> ctr_binning_map; // map from binning versions to the binnings

    std::map<std::pair<int, std::string>, std::vector<double>> run_yr_and_ctrbin_version_to_crossx_factors_map; // map from (run year, ctr binning version) to crossx factors
    std::map<std::string, std::vector<int>> ctr_rebin_scale_factor_map; // map from binning versions to the binnings
    std::map<std::string, std::string> ctr_binning_file_suffix_map; // map from binning versions to the suffices

    std::vector<int> ctr_rebin_scale_factors;

	std::string ctr_suffix;

	// -------- protected class methods --------

    void SanityCheckPbPb();
    void InitializePbPb();
    double CalculateWeightForRAA(int centrality, double weight);

public:
	// -------- public variables for configuration --------
    std::string ctr_binning_version = "default"; // centrality binning version used in current analysis (if there are several)
	std::vector<double> crossx_factors_ctr_binned;

	// -------- crossx factors for Pb+Pb years --------
    // 67.6 = sigma_{inel} in unit of mb = <TAA> / <N_{coll}>
    static std::vector<double> crossx_factors_pbpb_run2_ctr_binned;
    static std::vector<double> crossx_factors_pbpb_2023_ctr_binned;
    static std::vector<double> crossx_factors_pbpb_2024_ctr_binned;


    // -------- public class methods --------
     PbPbBaseClass(){}
    ~ PbPbBaseClass(){}
	
};

void PbPbBaseClass::InitializePbPb(){
    run_yr_and_ctrbin_version_to_crossx_factors_map[{15, "default"}] = crossx_factors_pbpb_run2_ctr_binned;
    run_yr_and_ctrbin_version_to_crossx_factors_map[{18, "default"}] = crossx_factors_pbpb_run2_ctr_binned;
    run_yr_and_ctrbin_version_to_crossx_factors_map[{23, "default"}] = crossx_factors_pbpb_2023_ctr_binned;
    run_yr_and_ctrbin_version_to_crossx_factors_map[{24, "default"}] = crossx_factors_pbpb_2024_ctr_binned;

    run_yr_and_ctrbin_version_to_crossx_factors_map[{15, "include_upc"}] = crossx_factors_pbpb_run2_ctr_binned; // WRONG!! 50-100% NEED TO BE CORRECTED
    run_yr_and_ctrbin_version_to_crossx_factors_map[{18, "include_upc"}] = crossx_factors_pbpb_run2_ctr_binned; // WRONG!! 50-100% NEED TO BE CORRECTED
    run_yr_and_ctrbin_version_to_crossx_factors_map[{23, "include_upc"}] = crossx_factors_pbpb_2023_ctr_binned; // WRONG!! 50-100% NEED TO BE CORRECTED
    run_yr_and_ctrbin_version_to_crossx_factors_map[{24, "include_upc"}] = crossx_factors_pbpb_2024_ctr_binned; // WRONG!! 50-100% NEED TO BE CORRECTED

    ctr_binning_map["default"] = {0,5,10,20,30,50,80};
    ctr_binning_map["include_upc"] = {0,5,10,20,30,50,100};

    ctr_binning_file_suffix_map["default"] = "";
    ctr_binning_file_suffix_map["include_upc"] = "_include_upc";
    
    ctr_rebin_scale_factor_map["default"] = {1,1,1,1,2,4};
    ctr_rebin_scale_factor_map["include_upc"] = {1,1,1,1,2,4};

    ctr_bin_edges = ctr_binning_map[ctr_binning_version];

	ctr_suffix = ctr_binning_file_suffix_map[ctr_binning_version];
    ctr_rebin_scale_factors = ctr_rebin_scale_factor_map[ctr_binning_version];

    ctr_bin_edges_double = std::vector<double>(ctr_bin_edges.begin(), ctr_bin_edges.end());
    nCtrBins = ctr_bin_edges.size() - 1;

    for (int ictr = 0; ictr < nCtrBins; ictr++){
        int ctr_bin_low_edge = ctr_bin_edges.at(ictr);
        int ctr_bin_high_edge = ctr_bin_edges.at(ictr + 1);
        ctr_bins.push_back("_ctr" + std::to_string(ctr_bin_low_edge) + "_" + std::to_string(ctr_bin_high_edge));
    }

    // SanityCheckPbPb(); // at this point, crossx_factors_ctr_binned is not well-defined
}

void PbPbBaseClass::SanityCheckPbPb(){
	if (crossx_factors_ctr_binned.size() != ctr_bin_edges.size() - 1){
        throw std::runtime_error("crossx_factors_ctr_binned MUST equal ctr_bin_edges size - 1");
    }
    if (ctr_rebin_scale_factors.size() != ctr_bin_edges.size() - 1){
        throw std::runtime_error("ctr_rebin_scale_factors MUST equal ctr_bin_edges size - 1");
    }
}

double PbPbBaseClass::CalculateWeightForRAA(int centrality, double weight){
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

std::vector<double> PbPbBaseClass::crossx_factors_pbpb_run2_ctr_binned = {
    1./(0.05 * 7.67 * 1000. * 26.2339 * (0.436882 + 1.3803) / 1000.), // 1/[crossx_{incl AA} [unit: mb] * TAA [unit: mb^{-1}] * (L_{int15} + L_{int18}) (unit: pb^{-1})] 0-5%
    1./(0.05 * 7.67 * 1000. * 20.4707 * (0.436882 + 1.3803) / 1000.), // 1/[crossx_{incl AA} [unit: mb] * TAA [unit: mb^{-1}] * (L_{int15} + L_{int18}) (unit: pb^{-1})] 5-10%
    1./(0.1 * 7.67 * 1000. * 14.3345 * (0.436882 + 1.3803) / 1000.), // 1/[crossx_{incl AA} [unit: mb] * TAA [unit: mb^{-1}] * (L_{int15} + L_{int18}) (unit: pb^{-1})] 10-20%
    1./(0.1 * 7.67 * 1000. * 8.63844 * (0.436882 + 1.3803) / 1000.), // 1/[crossx_{incl AA} [unit: mb] * TAA [unit: mb^{-1}] * (L_{int15} + L_{int18}) (unit: pb^{-1})] 20-30%
    1./(0.2 * 7.67 * 1000. * 3.79023 * (0.436882 + 1.3803) / 1000.), // 1/[crossx_{incl AA} [unit: mb] * TAA [unit: mb^{-1}] * (L_{int15} + L_{int18}) (unit: pb^{-1})] 30-50%
    1./(0.3 * 7.67 * 1000. * 0.689695 * (0.436882 + 1.3803) / 1000.) // 1/[crossx_{incl AA} [unit: mb] * TAA [unit: mb^{-1}] * (L_{int15} + L_{int18}) (unit: pb^{-1})] 50-80%
};

std::vector<double> PbPbBaseClass::crossx_factors_pbpb_2023_ctr_binned = {
    1./(0.05 * 7.8 * 1000. * 26.1428 * 1.3896 / 1000.), // 1/[crossx_{incl AA} [unit: mb] * TAA [unit: mb^{-1}] * L_{int} (unit: pb^{-1})] 0-5%
    1./(0.05 * 7.8 * 1000. * 20.3241 * 1.3896 / 1000.), // 1/[crossx_{incl AA} [unit: mb] * TAA [unit: mb^{-1}] * L_{int} (unit: pb^{-1})] 5-10%
    1./(0.1 * 7.8 * 1000. * 14.0502 * 1.3896 / 1000.), // 1/[crossx_{incl AA} [unit: mb] * TAA [unit: mb^{-1}] * L_{int} (unit: pb^{-1})] 10-20%
    1./(0.1 * 7.8 * 1000. * 8.5074 * 1.3896 / 1000.), // 1/[crossx_{incl AA} [unit: mb] * TAA [unit: mb^{-1}] * L_{int} (unit: pb^{-1})] 20-30%
    1./(0.2 * 7.8 * 1000. * 3.7733 * 1.3896 / 1000.), // 1/[crossx_{incl AA} [unit: mb] * TAA [unit: mb^{-1}] * L_{int} (unit: pb^{-1})] 30-50%
    1./(0.3 * 7.8 * 1000. * 0.6716 * 1.3896 / 1000.) // 1/[crossx_{incl AA} [unit: mb] * TAA [unit: mb^{-1}] * L_{int} (unit: pb^{-1})] 50-80%
};

std::vector<double> PbPbBaseClass::crossx_factors_pbpb_2024_ctr_binned = {
    1./(0.05 * 7.8 * 1000. * 26.1428 * 1.5411 / 1000.), // 1/[crossx_{incl AA} [unit: mb] * TAA [unit: mb^{-1}] * L_{int} (unit: pb^{-1})] 0-5%
    1./(0.05 * 7.8 * 1000. * 20.3241 * 1.5411 / 1000.), // 1/[crossx_{incl AA} [unit: mb] * TAA [unit: mb^{-1}] * L_{int} (unit: pb^{-1})] 5-10%
    1./(0.1 * 7.8 * 1000. * 14.0502 * 1.5411 / 1000.), // 1/[crossx_{incl AA} [unit: mb] * TAA [unit: mb^{-1}] * L_{int} (unit: pb^{-1})] 10-20%
    1./(0.1 * 7.8 * 1000. * 8.5074 * 1.5411 / 1000.), // 1/[crossx_{incl AA} [unit: mb] * TAA [unit: mb^{-1}] * L_{int} (unit: pb^{-1})] 20-30%
    1./(0.2 * 7.8 * 1000. * 3.7733 * 1.5411 / 1000.), // 1/[crossx_{incl AA} [unit: mb] * TAA [unit: mb^{-1}] * L_{int} (unit: pb^{-1})] 30-50%
    1./(0.3 * 7.8 * 1000. * 0.6716 * 1.5411 / 1000.) // 1/[crossx_{incl AA} [unit: mb] * TAA [unit: mb^{-1}] * L_{int} (unit: pb^{-1})] 50-80%
};
