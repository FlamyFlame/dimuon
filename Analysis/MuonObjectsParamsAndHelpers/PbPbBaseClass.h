#pragma once
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <type_traits>

template <class Derived>
class PbPbBaseClass
{
protected:
	// -------- protected class variables --------
    Derived& self() { return static_cast<Derived&>(*this); }
    const Derived& self() const { return static_cast<const Derived&>(*this); }

    // --- analysis-specific ctr binnings ---
    // Given run year & ctr binning version, these can/will be auto-filled in InitializePbPb()
    std::vector<int> ctr_bin_edges; // e.g, {0,10,...,80}
    std::vector<double> ctr_bin_edges_double;
    std::vector<std::string> ctr_bins; // e.g, {"_ctr0_10", ..., "_ctr50_80"}
    int nCtrBins;

    // --- list of defined ctr binning versions ---
    std::vector<std::string> ctr_binning_versions_allowed = {"default","include_upc"};

    // --- helper maps ---
    // maps ctr binning version to actual binning
    std::map<std::string, std::vector<int>> ctr_binning_map; // map from binning versions to the binnings

    // these map (run year &) ctr binning version to quantities/settings needed for analyses
    std::map<std::pair<int, std::string>, std::vector<double>> run_yr_and_ctrbin_version_to_crossx_factors_map; // map from (run year, ctr binning version) to crossx factors
    std::map<std::string, std::vector<int>> ctr_rebin_scale_factor_map; // map from binning versions to rebin scales
    std::map<std::string, std::string> ctr_binning_file_suffix_map; // map from binning versions to the suffices

    // --- quantities/settings needed for analysis ---
    std::vector<int> ctr_rebin_scale_factors;
	std::vector<double> crossx_factors_ctr_binned;
	std::string ctr_suffix;

	// -------- protected class methods --------

    void SanityCheckPbPb();
    void BuildPbPbMaps();
    void InitializePbPb();

    void SetCrossxFactorsPbPbCtrBinned(int run_yr, const std::string& ctr_binning_v);
    double CalculateWeightForRAA(int centrality, double weight);

public:
	// -------- public variables for configuration --------
    std::string ctr_binning_version = "default"; // centrality binning version used in current analysis (if there are several)

	// -------- crossx factors for Pb+Pb years --------
    // Returned by private helper functions (not static data members) so that cling does not need to
    // JIT-link template-static guard variables (_ZGVN...) which it cannot materialize.
    static std::vector<double> make_crossx_factors_pbpb_run2() {
        return {
            1./(0.05 * 7.67 * 1000. * 26.2339 * (0.436882 + 1.3803) / 1000.),
            1./(0.05 * 7.67 * 1000. * 20.4707 * (0.436882 + 1.3803) / 1000.),
            1./(0.1  * 7.67 * 1000. * 14.3345 * (0.436882 + 1.3803) / 1000.),
            1./(0.1  * 7.67 * 1000. * 8.63844 * (0.436882 + 1.3803) / 1000.),
            1./(0.2  * 7.67 * 1000. * 3.79023 * (0.436882 + 1.3803) / 1000.),
            1./(0.3  * 7.67 * 1000. * 0.689695* (0.436882 + 1.3803) / 1000.)
        };
    }
    static std::vector<double> make_crossx_factors_pbpb_2023() {
        return {
            1./(0.05 * 7.8 * 1000. * 26.1428 * 1.3896 / 1000.),
            1./(0.05 * 7.8 * 1000. * 20.3241 * 1.3896 / 1000.),
            1./(0.1  * 7.8 * 1000. * 14.0502 * 1.3896 / 1000.),
            1./(0.1  * 7.8 * 1000. *  8.5074 * 1.3896 / 1000.),
            1./(0.2  * 7.8 * 1000. *  3.7733 * 1.3896 / 1000.),
            1./(0.3  * 7.8 * 1000. *  0.6716 * 1.3896 / 1000.)
        };
    }
    static std::vector<double> make_crossx_factors_pbpb_2024() {
        return {
            1./(0.05 * 7.8 * 1000. * 26.1428 * 1.5411 / 1000.),
            1./(0.05 * 7.8 * 1000. * 20.3241 * 1.5411 / 1000.),
            1./(0.1  * 7.8 * 1000. * 14.0502 * 1.5411 / 1000.),
            1./(0.1  * 7.8 * 1000. *  8.5074 * 1.5411 / 1000.),
            1./(0.2  * 7.8 * 1000. *  3.7733 * 1.5411 / 1000.),
            1./(0.3  * 7.8 * 1000. *  0.6716 * 1.5411 / 1000.)
        };
    }
    // TODO: replace with actual 2025 T_AA and luminosity values once available
    static std::vector<double> make_crossx_factors_pbpb_2025() {
        return make_crossx_factors_pbpb_2023();
    }

    // -------- public class methods --------
     PbPbBaseClass(){
        BuildPbPbMaps();
     }
    ~ PbPbBaseClass(){}
	
};

template <class Derived>
void PbPbBaseClass<Derived>::BuildPbPbMaps(){
    // Building maps to map run year & ctr binning version to essential quantities/settings
    ctr_binning_map["default"] = {0,5,10,20,30,50,80};
    ctr_binning_map["include_upc"] = {0,5,10,20,30,50,100};

    ctr_binning_file_suffix_map["default"] = "";
    ctr_binning_file_suffix_map["include_upc"] = "_include_upc";
    
    ctr_rebin_scale_factor_map["default"] = {1,1,1,1,2,4};
    ctr_rebin_scale_factor_map["include_upc"] = {1,1,1,1,2,4};

    // 2015 + 2018 are combined --> use 2018 for Run2-combined crossx factors
    run_yr_and_ctrbin_version_to_crossx_factors_map[{18, "default"}] = make_crossx_factors_pbpb_run2();
    run_yr_and_ctrbin_version_to_crossx_factors_map[{23, "default"}] = make_crossx_factors_pbpb_2023();
    run_yr_and_ctrbin_version_to_crossx_factors_map[{24, "default"}] = make_crossx_factors_pbpb_2024();

    run_yr_and_ctrbin_version_to_crossx_factors_map[{25, "default"}]     = make_crossx_factors_pbpb_2025(); // TODO: update with 2025 values
    run_yr_and_ctrbin_version_to_crossx_factors_map[{18, "include_upc"}] = make_crossx_factors_pbpb_run2(); // WRONG!! 50-100% NEED TO BE CORRECTED
    run_yr_and_ctrbin_version_to_crossx_factors_map[{23, "include_upc"}] = make_crossx_factors_pbpb_2023(); // WRONG!! 50-100% NEED TO BE CORRECTED
    run_yr_and_ctrbin_version_to_crossx_factors_map[{24, "include_upc"}] = make_crossx_factors_pbpb_2024(); // WRONG!! 50-100% NEED TO BE CORRECTED
    run_yr_and_ctrbin_version_to_crossx_factors_map[{25, "include_upc"}] = make_crossx_factors_pbpb_2025(); // TODO: update with 2025 values

}

template <class Derived>
void PbPbBaseClass<Derived>::InitializePbPb(){
    // Initialize PbPb settings, assuming ctr_binning_version has been set by user
    
    // Run-year independent functionalities for Pb+Pb analyses 
    if (ctr_binning_version.empty() || 
        std::find(ctr_binning_versions_allowed.begin(), ctr_binning_versions_allowed.end(), ctr_binning_version) == ctr_binning_versions_allowed.end()){
        throw std::runtime_error("PbPbBaseClass::InitializePbPb: ctr_binning_version(" + ctr_binning_version + ") is NOT defined!");
    }

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
    
    // Optional: if Derived has an int-like member named run_year, set crossx factors & perform sanity checks
    if constexpr (requires (const Derived& d) { d.RunYear(); } &&
                  std::is_integral_v<std::remove_cvref_t<decltype(self().RunYear())>>) {

        const int ry = static_cast<int>(self().RunYear()) %2000;
        SetCrossxFactorsPbPbCtrBinned(ry, ctr_binning_version);
        SanityCheckPbPb();
    } else {
        std::cout   << "PbPbBaseClass::InitializePbPb: [WARNING] RunYear() not found in Derived; skipping crossx factors setting\n"
                    << "Crossx factors must be set manually by user!" << std::endl;
    }
}

template <class Derived>
void PbPbBaseClass<Derived>::SanityCheckPbPb(){
	if (crossx_factors_ctr_binned.size() != ctr_bin_edges.size() - 1){
        throw std::runtime_error("crossx_factors_ctr_binned MUST equal ctr_bin_edges size - 1");
    }
    if (ctr_rebin_scale_factors.size() != ctr_bin_edges.size() - 1){
        throw std::runtime_error("ctr_rebin_scale_factors MUST equal ctr_bin_edges size - 1");
    }
}

template <class Derived>
double PbPbBaseClass<Derived>::CalculateWeightForRAA(int centrality, double weight){
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
    return 0.; // assign weight = 0 to event if out of centrality-range of intere
}

template <class Derived>
void PbPbBaseClass<Derived>::SetCrossxFactorsPbPbCtrBinned(int run_yr, const std::string& ctr_binning_v){
    using Key = std::pair<int, std::string>;
    Key key{run_yr, ctr_binning_v};

    std::vector<double> crossx_factors_null = {};
    crossx_factors_null.assign(ctr_bins.size(), -1.);

    auto it = run_yr_and_ctrbin_version_to_crossx_factors_map.find(key);

    if (it != run_yr_and_ctrbin_version_to_crossx_factors_map.end()) {
        crossx_factors_ctr_binned = it->second;
    }else{        
        std::cerr << "[WARNING] PbPbBaseClass::SetCrossxFactorsPbPbCtrBinned: "
                  << "No cross-section factors found for "
                  << "run_yr=" << run_yr
                  << ", ctr_binning_v=\"" << ctr_binning_v << "\"" << std::endl
                  << "Return a list of -1." << std::endl;

        crossx_factors_ctr_binned = crossx_factors_null;
    }
}


// (out-of-line static definitions removed: crossx_factors now use inline static in class body above)
