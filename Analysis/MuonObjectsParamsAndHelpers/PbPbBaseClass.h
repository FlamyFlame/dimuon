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
    std::map<std::string, std::string> ctr_binning_file_suffix_map; // map from binning versions to the suffices

	std::string ctr_binning_version_suffix = ctr_binning_file_suffix_map[ctr_binning_verion];

	// -------- protected class methods --------

    SanityCheckPbPb();
    InitializePbPb();

public:
	// -------- public variables for configuration --------
    std::string ctr_binning_verion = "default"; // centrality binning version used in current analysis (if there are several)
	std::vector<double> crossx_factors_ctr_binned;

	// -------- crossx factors for Pb+Pb years --------
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

	// -------- public class methods --------
	 PbPbBaseClass();
	~ PbPbBaseClass();
	
};

PbPbBaseClass::InitializePbPb(){

    ctr_binning_map["default"] = {0,5,10,20,30,50,80};
    ctr_binning_file_suffix_map["default"] = "";

    ctr_bin_edges = ctr_binning_map[ctr_binning_verion];
	ctr_binning_version_suffix = ctr_binning_file_suffix_map[ctr_binning_verion];

    ctr_bin_edges_double = std::vector<double>(ctr_bin_edges.begin(), ctr_bin_edges.end());
    nCtrBins = ctr_bin_edges.size() - 1;

    for (int ictr = 0; ictr < nCtrBins; ictr++){
        int ctr_bin_low_edge = ctr_bin_edges.at(ictr);
        int ctr_bin_high_edge = ctr_bin_edges.at(ictr + 1);
        ctr_bins.push_back("_ctr" + std::to_string(ctr_bin_low_edge) + "_" + std::to_string(ctr_bin_high_edge));
    }

    SanityCheckPbPb();
}

PbPbBaseClass::SanityCheckPbPb(){
	if (crossx_factors_ctr_binned.size() != ctr_bin_edges.size() - 1){
        throw std::runtime_error("crossx_factors_ctr_binned MUST equal ctr_bin_edges size - 1");
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