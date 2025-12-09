#include "RDFBasedHistFilling.h"
#include <map>
#include <string>
#include <vector>
#include <array>
#include <iostream>
#include <algorithm>

void RDFBasedHistFilling::Run(){

	std::cout << "start the run" << std::endl;
	auto verbosity = ROOT::Experimental::RLogScopedVerbosity(ROOT::Detail::RDF::RDFLogChannel(), ROOT::Experimental::ELogLevel::kInfo);

	Initialize();

	TH1::SetDefaultSumw2(kTRUE); // turn on Sumw2 for all histograms

	ROOT::EnableImplicitMT();

    ProcessData();

	Finalize();
}

void RDFBasedHistFilling::ProcessData(){

    // ------- Create RDataFrame and keep ownership -------
    rdf_store.emplace_back(std::make_unique<ROOT::RDataFrame>(tree_ss, input_files));
    ROOT::RDF::RNode node_ss = *(rdf_store.back());
    df_map.emplace("df_ss", node_ss);

    rdf_store.emplace_back(std::make_unique<ROOT::RDataFrame>(tree_op, input_files));
    ROOT::RDF::RNode node_op = *(rdf_store.back());
    df_map.emplace("df_op", node_op);
    
    // ------- Fill histograms -------
    if (hist_filling_cycle == generic){
        FillHistograms();
    } else if (hist_filling_cycle == inv_weight_by_single_mu_effcy && trigger_mode == 1){
        FillTrigEffcyHistsInvWeightedbySingleMuonEffcies();
    }
    
    HistogramPostProcess();

    // info summary for the current for loop (input file)
    std::cout << "#event loop runs per data frame is: " << df_map.at("df_ss").GetNRuns() << std::endl;
    std::cout << "#slots per data frame is: " << df_map.at("df_ss").GetNSlots() << std::endl;
    std::cout << "# 1D histograms: " << hist1D_map.size() << std::endl;
    std::cout << "# 2D histograms: " << hist2D_map.size() << std::endl;
    std::cout << "# 3D histograms: " << hist2D_map.size() << std::endl;
}

void RDFBasedHistFilling::Initialize(){
	std::cout << "Calling Initialize." << std::endl;
	input_files.push_back("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/muon_pairs_pp_2024_combined_single_mu4.root");
    for (auto inFile : input_files){
	    if (!inFile || inFile->IsZombie()) {
	        std::cerr << "Error opening input file!" << std::endl;
	        throw std::exception();
	    }    	
    }
	output_file = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/histograms_real_pairs_pp_2024_single_mu4_new_RDF.root";

    hist1d_rresultptr_map.clear();
    hist2d_rresultptr_map.clear();
    hist3d_rresultptr_map.clear();
    hist1D_map.clear();
    hist2D_map.clear();
    hist3D_map.clear();

	BuildFilterToVarListMap();
	BuildHistBinningMap();
	ReadVar1DJson();
}

// ---------- ----------
void RDFBasedHistFilling::FillHistograms(){

    if (output_non_trig_effcy_hists){ // non-trigger-efficiency histograms
        for (std::string sign : {"_ss", "_op"}){
            FillHistogramsSingleDataFrame(sign, df_map.at(sign));
            FillHistogramsSingleDataFrame(sign, "_jacobian_corrected", df_map.at(sign));
            // FillHistogramsSingleDataFrame(sign + "_wgapcut", df_map.at(sign));
            // FillHistogramsSingleDataFrame(sign + "_wgapcut", "_jacobian_corrected", df_map.at(sign));
        }
    }

	if (trigger_mode == 1){ // trigger efficiencies: mu4_mu4noL1 | mu4 & 2mu4 | mu4
		try{
		    for (std::string pair_sign : {"_ss", "_op"}){
                std::string df_name = "df" + pair_sign;
				ROOT::RDF::RNode& node = df_map.at(df_name);

                for (auto mu4sel : {"_mu1passmu4", "_mu2passmu4"}){ // mu4 selection

                    df_name += mu4sel; // e.g, df_ss_mu1passmu4

					std::string ind1st = (mu4sel == "_mu1passmu4")? "1" : "2";
                    std::string ind2nd = (mu4sel == "_mu1passmu4")? "2" : "1";

					df_map.emplace(df_name, node.Filter("mu" + ind1st + "PassSingle"));
					df_map.at(df_name) = df_map.at(df_name).Define("pt1st",	"m" + ind1st + ".pt");
					df_map.at(df_name) = df_map.at(df_name).Define("pt2nd",	"m" + ind2nd + ".pt");
					df_map.at(df_name) = df_map.at(df_name).Define("charge2nd","m" + ind2nd + ".charge");
					df_map.at(df_name) = df_map.at(df_name).Define("eta2nd",	"m" + ind2nd + ".eta");
					df_map.at(df_name) = df_map.at(df_name).Define("phi2nd",	"m" + ind2nd + ".phi");
					df_map.at(df_name) = df_map.at(df_name).Define("q_eta2nd",	"charge2nd * eta2nd");
                    df_map.at(df_name) = df_map.at(df_name).Define("2nd_muon_good_acceptance",  "pt2nd >= 6 && ((eta2nd > 1.1 && eta2nd < 2.3) || (eta2nd > -2.3 && eta2nd < -1.2))");

					df_map.emplace(df_name + "_sign1", df_map.at(df_name).Filter("charge2nd > 0")); // e.g, df_ss_mu1passmu4_sign1
					df_map.emplace(df_name + "_sign2", df_map.at(df_name).Filter("charge2nd < 0"));

					for (auto mu_sign : {"_sign1", "_sign2"}){
                        df_name += mu_sign;
                        df_map.emplace(df_name + "_mu4", df_map.at(df_name));
                        df_map.emplace(df_name + "_mu4_mu4noL1", df_map.at(df_name).Filter("passmu4mu4noL1"));
                        df_map.emplace(df_name + "_2mu4", df_map.at(df_name).Filter("pass2mu4"));

                        for (auto trg : {"_mu4", "_mu4_mu4noL1", "_2mu4"}){
                            df_name += trg;
                            df_map.emplace(df_name + "_sepr", df_map.at(df_name).Filter("passSeparated"));
                            df_map.emplace(df_name + "_good_accept", df_map.at(df_name).Filter("2nd_muon_good_acceptance"));
                            
                            for (auto bias : {"", "_good_accept", "_sepr"}){ // additional selection / bias in data sample
                                df_name += bias;
                                std::string filter = df_name.substr(2);

                                FillHistogramsSingleDataFrame(filter, df_map.at(df_name), false); // do not write the sub-dataframe histograms in output file
                            }
                        }
					}
				} // end loop mu4 selection
	        } // end loop pair sign
		} catch(const std::out_of_range& e){
			std::cerr << "FillHistograms:: out_of_range exception caught: " << e.what() << ". Likely " << "df_" << pair_sign << " not defined in df_map!" << std::endl;
        } catch (const std::runtime_error& e) {
            std::cerr << "FillHistograms:: RDF runtime error: " << e.what() << std::endl;
        }
	}
}

// ---------- ----------
void RDFBasedHistFilling::HistogramPostProcess(){
    
    for (auto kv : hist1d_rresultptr_map){
        hist1D_map.emplace(kv.first, kv.second.GetPtr());
    }
    for (auto kv : hist2d_rresultptr_map){
        hist2D_map.emplace(kv.first, kv.second.GetPtr());
    }
    for (auto kv : hist3d_rresultptr_map){
        hist3D_map.emplace(kv.first, kv.second.GetPtr());
    }

    if (trigger_mode == 0 || trigger_mode == 1){
        if (hist_filling_cycle == generic){
            SumSingleMuonTrigEffHists();
            MakeAndWriteSingleMuonPtTrigEffGraphs();
            CalculateSingleMuonTrigEffcyRatios();           
        } else{
            MakeAndWriteDRTrigEffGraphs();
        }
    }
}


//--------- ---------


template <class TH, class Var, class MakeBaseName>
void SumTrigEffHistsGeneric(
    const std::vector<Var>& vars,
    const std::vector<std::string>& filters_post_sum,
    const std::vector<std::string>& filters_to_be_summed,
    std::map<std::string, TH*>& hist_map,
    MakeBaseName makeBaseName)        // only depends on Var
{
    for (const auto& var : vars) {
        const std::string base = makeBaseName(var);

        for (const auto& post_sum_filter : filters_post_sum) {

            std::string hname_post_sum = base + post_sum_filter;

            TH*  h_post_sum     = nullptr;
            bool first_instance = true;

            for (const auto& to_sum_filter : filters_to_be_summed) {
                std::string hname_pre_sum = base + to_sum_filter + post_sum_filter;

                auto it = hist_map.find(hname_pre_sum);
                if (it == hist_map.end()) {
                    std::cerr << "SumTrigEffHistsGeneric: histogram "
                              << hname_pre_sum << " not found!\n";
                    continue;
                }

                TH* h_pre_sum = it->second;

                if (first_instance) {
                    h_post_sum = static_cast<TH*>(
                        h_pre_sum->Clone(hname_post_sum.c_str())
                    );
                    first_instance = false;
                } else {
                    h_post_sum->Add(h_pre_sum);
                }
            }

            if (h_post_sum != nullptr) {
                hist_map.emplace(hname_post_sum, h_post_sum);
            }
        }
    }
}

void RDFBasedHistFilling::SumSingleMuonTrigEffHists(){

    // 1D
    SumTrigEffHistsGeneric<TH1D, std::string>(
        single_muon_trig_effcy_var1Ds,
        trg_effcy_filters_1D_post_sum,
        trg_effcy_filters_to_be_summed,
        hist1D_map,
        [](const std::string& var) {
            return "h_" + var;
        }
    );

    // 2D
    SumTrigEffHistsGeneric<TH2D, std::array<std::string,2>>(
        single_muon_trig_effcy_var2Ds,
        levels_trg_effcy_filters_2D_3D_post_sum,
        trg_effcy_filters_to_be_summed,
        hist2D_map,
        [](const std::array<std::string,2>& vars) {
            const std::string& varx = vars[0];
            const std::string& vary = vars[1];
            return "h_" + vary + "_vs_" + varx;
        }
    );

    // 3D
    SumTrigEffHistsGeneric<TH3D, std::array<std::string,3>>(
        single_muon_trig_effcy_var3Ds,
        levels_trg_effcy_filters_2D_3D_post_sum,
        trg_effcy_filters_to_be_summed,
        hist3D_map,
        [](const std::array<std::string,3>& vars) {
            const std::string& varx = vars[0];
            const std::string& vary = vars[1];
            const std::string& varz = vars[2];
            return "h_" + varz + "_vs_" + vary + "_vs_" + varx;
        }
    );
}

//--------- ---------
void RDFBasedHistFilling::MakeAndWriteSingleMuonPtTrigEffGraphs(){}

//--------- ---------
void RDFBasedHistFilling::CalculateSingleMuonTrigEffcyRatios(){}

//--------- ---------
static float RDFBasedHistFilling::EvaluateSingleMuonEffcyPtFitted(bool charge_sign, std::string trg, float pt_2nd, float q_eta_2nd, float phi_2nd){return 0.;}

//--------- ---------
static float RDFBasedHistFilling::EvaluateSingleMuonEffcy(bool charge_sign, std::string trg, float pt_2nd, float q_eta_2nd, float phi_2nd){return 0.;}

//--------- ---------
void RDFBasedHistFilling::FillTrigEffcyHistsInvWeightedbySingleMuonEffcies(){}

//--------- ---------
void RDFBasedHistFilling::MakeAndWriteDRTrigEffGraphs(){}

// ---------- ----------
void RDFBasedHistFilling::BuildFilterToVarListMap(){

	df_filter_to_var1D_and_weight_list_map[{"_single_b", "_crossx"}] = {"pair_pt", "pair_eta"}; // crossx as a weighting

	for (std::string sign : {"_ss", "_op"}){

		df_filter_to_var1D_list_map[sign + ""]         = {"pair_dP_overP", "Dphi", "DR", "DR_zoomin", "pair_y", "pt_asym", "pair_pt_ptlead_ratio"};
		df_filter_to_var1D_list_map[sign + "_wgapcut"] = {"pair_dP_overP", "Dphi", "DR", "DR_zoomin", "pair_y", "pt_asym", "pair_pt_ptlead_ratio"};

		df_filter_to_var1D_and_weight_list_map[{sign + ""			, "_jacobian_corrected"}] = {"DR", "DR_zoomin"};
		df_filter_to_var1D_and_weight_list_map[{sign + "_wgapcut"	, "_jacobian_corrected"}] = {"DR", "DR_zoomin"};

	}
	
	BuildTrgEffcyFilterToVarListMap();
}

// ---------- ----------
void RDFBasedHistFilling::BuildTrgEffcyFilterToVarListMap(){
    levels_trg_effcy_filters_1D_pre_sum = {{"_ss", "_op"}, // add a level of ctr bins for Pb+Pb
                                        {"_mu1passmu4", "_mu2passmu4"},
                                        {"_sign1", "_sign2"},
                                        {"_mu4", "_mu4_mu4noL1", "_2mu4"}, 
                                        {"", "_good_accept", "_sepr"}};

    levels_trg_effcy_filters_2D_3D_pre_sum = {{"_ss", "_op"}, // add a level of ctr bins for Pb+Pb
                                        {"_mu1passmu4", "_mu2passmu4"},
                                        {"_sign1", "_sign2"},
                                        {"_mu4", "_mu4_mu4noL1", "_2mu4"}, 
                                        {"", "_sepr"}};
    TrigEffcyFiltersPrePostSumFlattening();

    for (auto filter : trg_effcy_filters_1D_pre_sum)    df_filter_to_var1D_list_map[filter] = single_muon_trig_effcy_var1Ds;
    for (auto filter : trg_effcy_filters_2D_3D_pre_sum) df_filter_to_var2D_list_map[filter] = single_muon_trig_effcy_var2Ds;
    for (auto filter : trg_effcy_filters_2D_3D_pre_sum) df_filter_to_var3D_list_map[filter] = single_muon_trig_effcy_var3Ds;

    // ----- ADD INVERSE WEIGHTING ONES!! WITH _mu4_mu4noL1_denom, _2mu4_denom (paired with additional filtering) -----
}

//--------- Helper function   ---------
// Given indices of levels to be summed over, separate pre-sum levels (1D & 2D) into to-be-summed and post-sum levels
// And flatten all three levels
void MuonPairPlottingPP::TrigEffcyFiltersPrePostSumFlattening(){

    // ---- Helper that does the flattening ----
    auto flatten_levels = [](const std::vector<std::vector<std::string>>& levels,
                             std::vector<std::string>& flattened)
    {
        std::vector<int> level_sizes;
        int total_size = 1;

        for (const auto& vlevel : levels) {
            level_sizes.push_back(static_cast<int>(vlevel.size()));
            if (vlevel.size() != 0) total_size *= static_cast<int>(vlevel.size());
        }

        flattened.assign(total_size, "");  // allocate & zero-init

        // convert each flattened index into per-level indices
        for (int flattened_ind = 0; flattened_ind < total_size; ++flattened_ind) {
            int remaining_ind = flattened_ind;

            for (int level_ind = 0; level_ind < static_cast<int>(levels.size()); ++level_ind) {
                if (level_sizes.at(level_ind) == 0) continue; // if current level is empty, skip it

                int dim_remaining_levels = 1;
                for (int level_ind_remain = level_ind + 1;
                     level_ind_remain < static_cast<int>(levels.size());
                     ++level_ind_remain)
                {
                    if (level_sizes.at(level_ind_remain) != 0) dim_remaining_levels *= level_sizes.at(level_ind_remain);
                }

                int cur_level_index = remaining_ind / dim_remaining_levels;
                remaining_ind = remaining_ind % dim_remaining_levels;

                flattened.at(flattened_ind) += levels.at(level_ind).at(cur_level_index);
            }
        }
    };

    // ---- Use it for both cases ----
    flatten_levels(levels_trg_effcy_filters_1D_pre_sum, trg_effcy_filters_1D_pre_sum);
    flatten_levels(levels_trg_effcy_filters_2D_3D_pre_sum, trg_effcy_filters_2D_3D_pre_sum);

    auto write_post_sum_levels = [this](const std::vector<std::vector<std::string>>& levels_pre_sum, std::vector<std::vector<std::string>>& levels_post_sum){
        for (int level_ind = 0; level_ind < static_cast<int>(levels_pre_sum.size()); ++level_ind) {
            if (std::find(levels_trg_effcy_to_be_summed.begin(), levels_trg_effcy_to_be_summed.end(), level_ind) == levels_trg_effcy_to_be_summed.end()){ // level is NOT contracted
                levels_post_sum.push_back(levels_pre_sum.at(level_ind));
            }
        }
    };

    write_post_sum_levels(levels_trg_effcy_filters_1D_pre_sum, levels_trg_effcy_filters_1D_post_sum);
    write_post_sum_levels(levels_trg_effcy_filters_2D_3D_pre_sum, levels_trg_effcy_filters_2D_3D_post_sum);


    for (int level_ind = 0; level_ind < static_cast<int>(levels_trg_effcy_filters_2D_3D_pre_sum.size()); ++level_ind) {
        if (std::find(levels_trg_effcy_to_be_contracted.begin(), levels_trg_effcy_to_be_contracted.end(), level_ind) != levels_trg_effcy_to_be_contracted.end()){ // level is NOT contracted
            levels_trg_effcy_filters_to_be_summed.push_back(levels_trg_effcy_filters_2D_3D_pre_sum.at(level_ind));
        }
    }

    flatten_levels(levels_trg_effcy_filters_1D_post_sum, trg_effcy_filters_1D_post_sum);
    flatten_levels(levels_trg_effcy_filters_2D_3D_post_sum, trg_effcy_filters_2D_3D_post_sum);
    flatten_levels(levels_trg_effcy_filters_to_be_summed, trg_effcy_filters_to_be_summed);
}

// ---------- ----------
void RDFBasedHistFilling::BuildHistBinningMap(){

	// ------- pT binnings -------
    hist_binning_map["pT_bins_80"] = pms.pT_bins_80;

	// ------- pT binning for single-muon trigger efficieny -------

    std::vector<double> pT_bins_single_muon (pms.pT_bins_8); // make a copy of a suitable set of single-muon pT bins (adjustable) --> use the copy for histogram settings

    pT_bins_single_muon.insert(pT_bins_single_muon.end(), pms.pT_bins_60.begin(), pms.pT_bins_60.end());

    hist_binning_map["pT_bins_single_muon"] = pT_bins_single_muon;

    // ------- eta binning for single-muon trigger efficiency -------

    std::vector<double> eta_bins_trig_effcy = ParamsSet::makeEtaTrigEffcyBinning();
    int neta_bins_trig_effcy = eta_bins_trig_effcy.size() - 1;
    hist_binning_map["eta_bins_trig_effcy"] = eta_bins_trig_effcy;
    
    // ------- eta binning for single-muon trigger efficiency -------

    // Build uniform phi edges so we can use the (xbins, ybins, zbins) TH3D ctor
    int nphi_bins_trig_effcy = 128; // phi 2nd muon
    
    std::vector<double> phi2nd_bins(nphi_bins_trig_effcy + 1);
    for (int i = 0; i <= nphi_bins_trig_effcy; ++i) {
        phi2nd_bins[i] = -pms.PI + (2.0 * pms.PI) * (static_cast<double>(i) / nphi_bins_trig_effcy);
    }

    hist_binning_map["phi2nd_bins"] = phi2nd_bins;

	// ------- minv log binning -------

	static const int nminv_bins_log = 40;
    float minv_logpow[ParamsSet::nSigns];
    minv_logpow[0] = 0.062;
    minv_logpow[1] = 0.045;
    float minv_max = 60;
    std::vector<double>minv_bins_log[ParamsSet::nSigns];

    std::string dphi_regions[2] = {"near", "away"};

    for (unsigned int isign = 0; isign < ParamsSet::nSigns; isign++){
        for(int iminv = 0; iminv <= nminv_bins_log; iminv++){
            minv_bins_log[isign].push_back(minv_max * pow(10.0, (static_cast<float>(iminv - nminv_bins_log)) * minv_logpow[isign]));
        }
    }

    hist_binning_map["minv_log_bins_ss"] = minv_bins_log[0];
    hist_binning_map["minv_log_bins_op"] = minv_bins_log[1];
}

// ---------- WRITE OUTPUT & FINALIZE ----------
template <typename H>
void write_hist_map_vector(
    std::map<std::string, std::vector<H>>& m,
    const std::vector<std::string>& hists_to_not_write)
{
    for (auto& kv : m) {
        const auto& name = kv.first;
        if (std::find(hists_to_not_write.begin(), hists_to_not_write.end(), name)
            != hists_to_not_write.end()) continue;

        for (auto& h : kv.second) {
            h.Write();        // THnD in a vector
        }
    }
}

void RDFBasedHistFilling::Finalize(){

    // ----- write output -----

    std::string file_mode = (hist_filling_cycle == generic)? "recreate" : "update";
    TFile * m_outfile=new TFile(output_file.c_str(), file_mode.c_str());
    
    // std::vector<THnD> maps
    write_hist_map_vector(hist1D_map, hists_to_not_write);
    write_hist_map_vector(hist2D_map, hists_to_not_write);
    write_hist_map_vector(hist3D_map, hists_to_not_write);

    // ----- delete pointers -----
    delete m_outfile;

    for (auto kv : var1D_dict){
        delete kv.second;
    }
}

// ---------- THROW MISSING FIELD IN JSON READING HELPER ----------
void JsonTestClass::ThrowMissingField(const std::string& field,
                                const std::string& histName) {
    std::ostringstream oss;
    oss << "ReadVar1DJson: mandatory field '" << field
        << "' is missing or empty";
    if (!histName.empty()) {
        oss << " (hist_name='" << histName << "')";
    }
    throw std::runtime_error(oss.str());
}

// ---------- VAR1D JSON FILE READING ----------
void JsonTestClass::ReadVar1DJson() {
    // Open file
    std::ifstream in(infile_var1D_json);
    if (!in) {
        throw std::runtime_error(
            "ReadVar1DJson: cannot open JSON file '" + infile_var1D_json + "'");
    }

    // Parse JSON
    json j;
    try {
        in >> j;
    } catch (const std::exception& e) {
        throw std::runtime_error(
            "ReadVar1DJson: failed to parse JSON file '" + infile_var1D_json +
            "': " + e.what());
    }

    // Check that "variables1D" exists and is an array
    if (!j.contains("variables1D") || !j["variables1D"].is_array()) {
        throw std::runtime_error(
            "ReadVar1DJson: mandatory top-level field 'variables1D' "
            "is missing or is not an array");
    }

    // Loop over entries
    for (const auto& jv : j["variables1D"]) {
        // hist_name (mandatory, non-empty string)
        std::string histName;
        if (jv.contains("hist_name") && jv["hist_name"].is_string()) {
            histName = jv["hist_name"].get<std::string>();
        }
        if (histName.empty()) {
            ThrowMissingField("hist_name", /*histName*/ "");
        }

        // Create var1D object
        auto* v = new var1D{};
        v->name = histName;

        // var_name: optional, default to name
        if (jv.contains("var_name") && jv["var_name"].is_string()) {
            std::string tmp = jv["var_name"].get<std::string>();
            if (!tmp.empty()) {
                v->var = tmp;
            }
        }
        if (v->var.empty()) {
            v->var = v->name;
        }

        // title: optional, default to name
        if (jv.contains("title") && jv["title"].is_string()) {
            std::string tmp = jv["title"].get<std::string>();
            if (!tmp.empty()) {
                v->title = tmp;
            }
        }
        if (v->title.empty()) {
            v->title = v->name;
        }

        // binning: mandatory
        if (!jv.contains("binning") || jv["binning"].is_null()) {
            delete v;
            ThrowMissingField("binning", histName);
        }

        const auto& jb = jv["binning"];

        if (jb.is_string()) {
            // binning by name: look up in hist_binning_map
            std::string binName = jb.get<std::string>();
            if (binName.empty()) {
                delete v;
                ThrowMissingField("binning", histName);
            }

            auto it = hist_binning_map.find(binName);
            auto it_ss = hist_binning_map.find(binName + "_ss");
            auto it_op = hist_binning_map.find(binName + "_op");
            if (it == hist_binning_map.end() && (it_ss == hist_binning_map.end() || it_op == hist_binning_map.end())) {
                delete v;
                std::ostringstream oss;
                oss << "ReadVar1DJson: binning name '" << binName
                    << "' (referenced by hist_name='" << histName
                    << "') not found in hist_binning_map";
                throw std::runtime_error(oss.str());
            }

            if (it_ss != hist_binning_map.end() && it_op != hist_binning_map.end()){
                v->bins_ss = it_ss->second;
                v->nbins_ss = static_cast<int>(v->bins_ss.size()) - 1;
                v->bins_op = it_op->second;
                v->nbins_op = static_cast<int>(v->bins_op.size()) - 1;
            } else{
                v->bins = it->second;
                v->nbins = static_cast<int>(v->bins.size()) - 1;                
            }

        } else if (jb.is_object()) {
            // binning as {nbins, min, max}
            if (!jb.contains("nbins") || !jb.contains("min") || !jb.contains("max")) {
                delete v;
                ThrowMissingField("binning.nbins/min/max", histName);
            }

            if (!jb["nbins"].is_number_integer() ||
                !jb["min"].is_number() ||
                !jb["max"].is_number()) {
                delete v;
                std::ostringstream oss;
                oss << "ReadVar1DJson: binning object for hist_name='"
                    << histName
                    << "' must have integer 'nbins' and numeric 'min', 'max'";
                throw std::runtime_error(oss.str());
            }

            v->nbins = jb["nbins"].get<int>();
            v->vmin  = jb["min"].get<double>();
            v->vmax  = jb["max"].get<double>();

        } else {
            delete v;
            std::ostringstream oss;
            oss << "ReadVar1DJson: 'binning' for hist_name='" << histName
                << "' must be either a string or an object {nbins,min,max}";
            throw std::runtime_error(oss.str());
        }

        // Final validation using your isValid()
        if (!v->isValid()) {
            std::ostringstream oss;
            oss << "ReadVar1DJson: invalid binning specification for hist_name='"
                << histName << "'";
            delete v;
            throw std::runtime_error(oss.str());
        }

        // All good, store pointer in dictionary
        var1D_dict[v->name] = v;
    }
}

// ---------- VAR1D LIST PRINTINT FOR DEBUG ----------
void JsonTestClass::PrintVar1DList() const {
    std::cout << "===== var1D_dict contents (" 
              << var1D_dict.size() << " entries) =====\n";

    for (auto pair : var1D_dict){
        const var1D* v = pair.second;
        if (!v) continue;

        std::cout << "  name  = " << v->name  << "\n"
                  << "  var   = " << v->var   << "\n"
                  << "  title = " << v->title << "\n";

        if (!v->bins_ss.empty()) {
            std::cout << "Same-sign variable-binning (" << v->nbins_ss
                      << " bins): [";
            for (size_t j = 0; j < v->bins_ss.size(); ++j) {
                std::cout << v->bins_ss[j];
                if (j + 1 < v->bins_ss.size()) std::cout << ", ";
            }
            std::cout << "]\n";

            std::cout << "Opposite-sign variable-binning (" << v->nbins_op
                      << " bins): [";
            for (size_t j = 0; j < v->bins_op.size(); ++j) {
                std::cout << v->bins_op[j];
                if (j + 1 < v->bins_op.size()) std::cout << ", ";
            }
            std::cout << "]\n";
        } else if (!v->bins.empty()) {
            std::cout << "  variable-binning (" << v->nbins
                      << " bins): [";
            for (size_t j = 0; j < v->bins.size(); ++j) {
                std::cout << v->bins[j];
                if (j + 1 < v->bins.size()) std::cout << ", ";
            }
            std::cout << "]\n";
        } else {
            std::cout << "  nbins = " << v->nbins << "\n"
                      << "  vmin  = " << v->vmin  << "\n"
                      << "  vmax  = " << v->vmax  << "\n";
        }

        std::cout << "  isValid = " << (v->isValid() ? "true" : "false") << "\n";
        std::cout << "--------------------------------------------\n";
    }
}

// ---------- HELPER FUNCTIONS ----------

void TestClass::FillHistogramsSingleDataFrame(const std::string& filter,
                                            ROOT::RDF::RNode df,
                                            bool hists_not_write = false,
                                            std::array<bool, 3> hists_1_2_3D_not_write = {0,0,0}) {
    static const std::vector<std::string> empty1D;
    static const std::vector<std::array<std::string, 2>> empty2D;
    static const std::vector<std::array<std::string, 3>> empty3D;

    auto it1D = df_filter_to_var1D_list_map.find(filter);
    auto it2D = df_filter_to_var2D_list_map.find(filter);
    auto it3D = df_filter_to_var3D_list_map.find(filter);

    const auto& vars1D = (it1D != df_filter_to_var1D_list_map.end()) ? it1D->second : empty1D;
    const auto& vars2D = (it2D != df_filter_to_var2D_list_map.end()) ? it2D->second : empty2D;
    const auto& vars3D = (it3D != df_filter_to_var3D_list_map.end()) ? it3D->second : empty3D;

    FillHistogramsSingleDataFrame(filter, df, "", vars1D, vars2D, vars3D, hists_not_write, hists_1_2_3D_not_write);
}

void TestClass::FillHistogramsSingleDataFrame(const std::string& filter,
                                            const std::string& weight,
                                            ROOT::RDF::RNode df,
                                            bool weight_before_filter = false,
                                            bool hists_not_write = false,
                                            std::array<bool, 3> hists_1_2_3D_not_write = {0,0,0}) {
    static const std::vector<std::string> empty1D;
    static const std::vector<std::array<std::string, 2>> empty2D;
    static const std::vector<std::array<std::string, 3>> empty3D;

    std::pair<std::string, std::string> key(filter, weight);

    auto it1D = df_filter_and_weight_to_var1D_list_map.find(key);
    auto it2D = df_filter_and_weight_to_var2D_list_map.find(key);
    auto it3D = df_filter_and_weight_to_var3D_list_map.find(key);

    const auto& vars1D = (it1D != df_filter_and_weight_to_var1D_list_map.end()) ? it1D->second : empty1D;
    const auto& vars2D = (it2D != df_filter_and_weight_to_var2D_list_map.end()) ? it2D->second : empty2D;
    const auto& vars3D = (it3D != df_filter_and_weight_to_var3D_list_map.end()) ? it3D->second : empty3D;

    // resolve weight column
    std::string wCol;
    auto itW = weight_specifier_to_column_map.find(weight);
    if (itW != weight_specifier_to_column_map.end()) {
        wCol = itW->second;
    } else {
        wCol = weight.substr(1); // trim off the "_" at the beginning
    }

    std::string suffix = (weight_before_filter)? weight + filter : filter + weight;
    FillHistogramsSingleDataFrame(suffix, df, wCol, vars1D, vars2D, vars3D, hists_not_write, hists_1_2_3D_not_write);
}

void TestClass::FillHistogramsSingleDataFrame(const std::string& suffix, // filter or filter & weight concatenated with custom order
                                            ROOT::RDF::RNode df,
                                            const std::string& weight_col, // weight column, "" if unweighted
                                            const std::vector<std::string>& vars1D,
                                            const std::vector<std::array<std::string, 2>>& vars2D, 
                                            const std::vector<std::array<std::string, 3>>& vars3D,
                                            bool hists_not_write = false,
                                            std::array<bool, 3> hists_1_2_3D_not_write = {0,0,0}) {


    if (vars1D.empty() && vars2D.empty() && vars3D.empty()) {
        std::cerr << "Warning: FillHistogramsSingleDataFrame:: no variables configured "
                  << "for suffix '" << suffix << "'\n";
        return;
    }

    // 1D
    for (const auto& vName : vars1D) {
        var1D* v = nullptr;
        try {
            v = Var1DSearch(vName);
        } catch (const std::exception& e) {
            std::cerr << e.what() << "\n";
            throw; // as you requested: error if not found
        }

        AxisInfo b = GetAxisInfo(*v, suffix);
        std::string hname = "h_" + v->name + suffix;
        std::string htitle = ";" + v->title + ";";

        if (hists_not_write || hists_1_2_3D_not_write.at(0)){
            hists_to_not_write.push_back(hname);
        }
        
        try {
            if (b.bin_edges) {
                ROOT::RDF::TH1DModel model(hname.c_str(), htitle.c_str(),
                                           b.nbins, b.bin_edges);
                if (weight_col.empty()) hist1d_rresultptr_map[hname] = df.Histo1D(model, v->var);
                else                    hist1d_rresultptr_map[hname] = df.Histo1D(model, v->var, weight_col);
            } else {
                ROOT::RDF::TH1DModel model(hname.c_str(), htitle.c_str(),
                                           b.nbins, b.min, b.max);
                if (weight_col.empty()) hist1d_rresultptr_map[hname] = df.Histo1D(model, v->var);
                else                    hist1d_rresultptr_map[hname] = df.Histo1D(model, v->var, weight_col);
            }
        } catch (const std::exception& e) {
            std::cerr << "FillHistogramsSingleDataFrame (1D) exception for suffix '"
                      << suffix << "', column '" << v->var
                      << "': " << e.what() << "\n";
            continue;
        }
    }

    // 2D
    for (const auto& arr : vars2D) {
        const std::string& xName = arr[0];
        const std::string& yName = arr[1];

        var1D *vx = nullptr, *vy = nullptr;
        try {
            vx = Var1DSearch(xName);
            vy = Var1DSearch(yName);
        } catch (const std::exception& e) {
            std::cerr << e.what() << "\n";
            throw;
        }

        AxisInfo bx = GetAxisInfo(*vx, suffix);
        AxisInfo by = GetAxisInfo(*vy, suffix);

        std::string hname  = "h_" + vy->name + "_vs_" + vx->name + suffix;
        ROOT::RDF::TH2DModel model = MakeTH2DModel(hname, vx->title, vy->title, bx, by);

        if (hists_not_write || hists_1_2_3D_not_write.at(1)){
            hists_to_not_write.push_back(hname);
        }

        try {
            if (weight_col.empty()) hist2d_rresultptr_map[hname] = df.Histo2D(model, vx->var, vy->var);
            else                    hist2d_rresultptr_map[hname] = df.Histo2D(model, vx->var, vy->var, weight_col);
        } catch (const std::exception& e) {
            std::cerr << "FillHistogramsSingleDataFrame (2D) exception for suffix '"
                      << suffix << "', columns ('" << vx->var << "', '"
                      << vy->var << "'): " << e.what() << "\n";
            continue;
        }
    }

    // 3D
    for (const auto& arr : vars3D) {
        const std::string& xName = arr[0];
        const std::string& yName = arr[1];
        const std::string& zName = arr[2];

        var1D *vx = nullptr, *vy = nullptr, *vz = nullptr;
        try {
            vx = Var1DSearch(xName);
            vy = Var1DSearch(yName);
            vz = Var1DSearch(zName);
        } catch (const std::exception& e) {
            std::cerr << e.what() << "\n";
            throw;
        }

        AxisInfo bx = GetAxisInfo(*vx, suffix);
        AxisInfo by = GetAxisInfo(*vy, suffix);
        AxisInfo bz = GetAxisInfo(*vz, suffix);

        std::string hname  = "h_" + vz->name + "_vs_" + vy->name + "_vs_" + vx->name + suffix;
        ROOT::RDF::TH3DModel model = MakeTH3DModel(hname, vx->title, vy->title, vz->title,
                                                   bx, by, bz);
        
        if (hists_not_write || hists_1_2_3D_not_write.at(2)){
            hists_to_not_write.push_back(hname);
        }

        try {
            if (weight_col.empty()) hist3d_rresultptr_map[hname] = df.Histo3D(model, vx->var, vy->var, vz->var);
            else                    hist3d_rresultptr_map[hname] = df.Histo3D(model, vx->var, vy->var, vz->var, weight_col);

        } catch (const std::exception& e) {
            std::cerr << "FillHistogramsSingleDataFrame (3D) exception for suffix '"
                      << suffix << "', columns ('" << vx->var << "', '"
                      << vy->var << "', '" << vz->var << "'): " << e.what() << "\n";
            continue;
        }
    }
}

var1D* TestClass::Var1DSearch(const std::string& var1DName) const {
    auto it = var1D_dict.find(var1DName);
    if (it == var1D_dict.end() || !(it->second)) {
        std::ostringstream oss;
        oss << "Var1DSearch: variable '" << var1DName << "' not found in var1D_dict";
        throw std::runtime_error(oss.str());
    }
    return it->second;
}

AxisInfo TestClass::GetAxisInfo(const var1D& v, const std::string& filter) const {
    bool hasSS = (filter.find("ss") != std::string::npos);
    bool hasOP = (filter.find("op") != std::string::npos);

    if (hasSS && hasOP) {
        std::ostringstream oss;
        oss << "GetAxisInfo: filter '" << filter
            << "' contains both 'ss' and 'op'; ambiguous for variable '" << v.name << "'";
        throw std::runtime_error(oss.str());
    }

    AxisInfo b;

    // default: use global nbins/bins
    int nbins = v.nbins;
    const std::vector<double>* binsVec = nullptr;
    if (!v.bins.empty()) binsVec = &v.bins;
    double vmin = v.vmin;
    double vmax = v.vmax;

    // ss/op override
    if (hasSS && v.nbins_ss > 0 && !v.bins_ss.empty()) {
        nbins = v.nbins_ss;
        binsVec = &v.bins_ss;
    } else if (hasOP && v.nbins_op > 0 && !v.bins_op.empty()) {
        nbins = v.nbins_op;
        binsVec = &v.bins_op;
    }

    b.nbins = nbins;

    if (binsVec && !binsVec->empty()) {
        b.bin_edges = binsVec->data();
        b.min  = 0.;
        b.max  = 0.;
    } else {
        b.bin_edges = nullptr;
        b.min  = vmin;
        b.max  = vmax;
    }

    return b;
}

ROOT::RDF::TH2DModel TestClass::MakeTH2DModel(const std::string& hname,
                                            const std::string& xtitle,
                                            const std::string& ytitle,
                                            const AxisInfo& bx,
                                            const AxisInfo& by) const {
    std::string htitle = ";" + xtitle + ";" + ytitle;

    if (bx.bin_edges && by.bin_edges) {
        return ROOT::RDF::TH2DModel(hname.c_str(), htitle.c_str(),
                                    bx.nbins, bx.bin_edges,
                                    by.nbins, by.bin_edges);
    } else if (bx.bin_edges && !by.bin_edges) {
        return ROOT::RDF::TH2DModel(hname.c_str(), htitle.c_str(),
                                    bx.nbins, bx.bin_edges,
                                    by.nbins, by.min, by.max);
    } else if (!bx.bin_edges && by.bin_edges) {
        return ROOT::RDF::TH2DModel(hname.c_str(), htitle.c_str(),
                                    bx.nbins, bx.min, bx.max,
                                    by.nbins, by.bin_edges);
    } else {
        return ROOT::RDF::TH2DModel(hname.c_str(), htitle.c_str(),
                                    bx.nbins, bx.min, bx.max,
                                    by.nbins, by.min, by.max);
    }
}

ROOT::RDF::TH3DModel TestClass::MakeTH3DModel(const std::string& hname,
                                              const std::string& xtitle,
                                              const std::string& ytitle,
                                              const std::string& ztitle,
                                              const AxisInfo& bx,
                                              const AxisInfo& by,
                                              const AxisInfo& bz) const
{
    std::string htitle = ";" + xtitle + ";" + ytitle + ";" + ztitle;

    const bool ex = (bx.bin_edges != nullptr);
    const bool ey = (by.bin_edges != nullptr);
    const bool ez = (bz.bin_edges != nullptr);

    // 1) NONE has bin_edges → use uniform on all
    if (!ex && !ey && !ez) {
        return ROOT::RDF::TH3DModel(hname.c_str(), htitle.c_str(),
                                    bx.nbins, bx.min, bx.max,
                                    by.nbins, by.min, by.max,
                                    bz.nbins, bz.min, bz.max);
    }

    // 2) At least one has bin_edges → all-variable constructor must be used
    // Define storage for any axis that lacks explicit bin_edges
    std::vector<double> xStore, yStore, zStore;

    auto makeEdges = [](const AxisInfo& a, std::vector<double>& store) -> const double* {
        if (a.bin_edges) return a.bin_edges;
        store.resize(a.nbins + 1);
        double step = (a.max - a.min) / a.nbins;
        for (int i = 0; i <= a.nbins; ++i)
            store[i] = a.min + i * step;
        return store.data();
    };

    const double* xEdges = makeEdges(bx, xStore);
    const double* yEdges = makeEdges(by, yStore);
    const double* zEdges = makeEdges(bz, zStore);

    return ROOT::RDF::TH3DModel(hname.c_str(), htitle.c_str(),
                                bx.nbins, xEdges,
                                by.nbins, yEdges,
                                bz.nbins, zEdges);
}








