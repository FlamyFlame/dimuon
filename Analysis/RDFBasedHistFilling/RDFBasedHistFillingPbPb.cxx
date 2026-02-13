#include "RDFBasedHistFillingData.cxx"

void RDFBasedHistFillingPbPb::SetIOPathsHook(){
	infile_var1D_json = "var1D_pbpb.json";

	std::string run_year_str = std::to_string(run_year);

	if (run_year == 23 || run_year == 24 || run_year == 25){
		input_files.push_back("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_20" + run_year_str + "/muon_pairs_pbpb_20" + run_year_str + trig_suffix + "_res_cut_v2.root");
		output_file = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_20" + run_year_str + "/histograms_real_pairs_pbpb_20" + run_year_str + out_file_suffix + ".root";
	} else if (run_year == 15 || run_year == 18){
		input_files.push_back("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_run2/muon_pairs_pbpb_20" + run_year_str + trig_suffix + "_res_cut_v2.root");
		output_file = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_run2/histograms_real_pairs_pbpb_20" + run_year_str + out_file_suffix + ".root";
	} else{
	    throw std::runtime_error("Run year must be 15/18/23/24/25! Current input invalid: " + run_year_str);
	}
}

void RDFBasedHistFillingPbPb::InitializePbPbExtra(){

	PbPbBaseClass::InitializePbPb();

    levels_trg_effcy_filters_1D_pre_sum = {{"_ss", "_op"}, // add a level of ctr bins for Pb+Pb
                                        {"_mu1passmu4", "_mu2passmu4"},
                                        ctr_bins,
                                        musigns,
                                        trigs,
                                        trg_effcy_biases};

    levels_trg_effcy_filters_2D_3D_pre_sum = {{"_ss", "_op"}, // add a level of ctr bins for Pb+Pb
                                        {"_mu1passmu4", "_mu2passmu4"},
                                        ctr_bins,
                                        musigns,
                                        trigs,
                                        trg_effcy_biases};

	levels_trg_effcy_filters_1D_pre_sum_w_musign_ctr_summing 
									=  {{"_ss", "_op"}, // add a level of ctr bins for Pb+Pb
                                        {"_mu1passmu4", "_mu2passmu4"},
                                        trigs,
                                        trg_effcy_biases};

	levels_trg_effcy_filters_2D_3D_pre_sum_w_musign_ctr_summing
									=  {{"_ss", "_op"}, // add a level of ctr bins for Pb+Pb
                                        {"_mu1passmu4", "_mu2passmu4"},
                                        trigs,
                                        trg_effcy_biases};

	levels_trg_effcy_filters_1D_pre_sum_w_ctr_summing
									=  {{"_ss", "_op"}, // add a level of ctr bins for Pb+Pb
                                        {"_mu1passmu4", "_mu2passmu4"},
                                        musigns,
                                        trigs,
                                        trg_effcy_biases};

	levels_trg_effcy_filters_2D_3D_pre_sum_w_ctr_summing
									=  {{"_ss", "_op"}, // add a level of ctr bins for Pb+Pb
                                        {"_mu1passmu4", "_mu2passmu4"},
                                        musigns,
                                        trigs,
                                        trg_effcy_biases};


    for (int ipt = 0; ipt < pT_bins_edges_for_trg_effcy_ctr_dep.size() - 1; ipt++){
        int pt_bin_low_edge = pT_bins_edges_for_trg_effcy_ctr_dep.at(ipt);
        int pt_bin_high_edge = pT_bins_edges_for_trg_effcy_ctr_dep.at(ipt + 1);
        pT_bins_for_trg_effcy_ctr_dep.push_back("_pt" + std::to_string(pt_bin_low_edge) + "_" + std::to_string(pt_bin_high_edge));
    }

    levels_trg_effcy_filters_ctr_dep_pre_sum = {{"_ss", "_op"}, // add a level of ctr bins for Pb+Pb
                                        {"_mu1passmu4", "_mu2passmu4"},
                                        pT_bins_for_trg_effcy_ctr_dep,
										q_eta_ranges_str,
                                        trigs,
                                        trg_effcy_biases};

	for (std::string ctr : ctr_bins){
    	for (std::string sign : {"_sign1", "_sign2"}){
    		ctr_musign_categories.push_back(ctr + sign);
	    }
	}

    for (std::string pair_sign : pair_signs){
		for (std::string ctr : ctr_bins){
    		categories_essential.push_back(pair_sign + ctr);
	    }
	}
}

// ---------- ----------
void RDFBasedHistFillingPbPb::BuildHistBinningMapPbPbExtra(){

	hist_binning_map["ctr_bins"] = ctr_bin_edges_double;

	for (int ictr = 0; ictr < ctr_bins.size(); ictr++){
		auto ctr = ctr_bins.at(ictr);
		int ctr_rebin_factor = ctr_rebin_scale_factors.at(ictr);

	    // ------- ctr-dep eta binning for single-muon trigger efficiency -------

	    std::vector<double> eta_bins_trig_effcy = ParamsSet::makeEtaTrigEffcyBinning(ctr_rebin_factor);
	    hist_binning_map["eta_bins_trig_effcy" + ctr] = eta_bins_trig_effcy;
	    
	    // ------- ctr-dep phi binning for single-muon trigger efficiency -------

	    // Build uniform phi edges so we can use the (xbins, ybins, zbins) TH3D ctor
	    int nphi_bins_trig_effcy = 128 / ctr_rebin_factor;
	    
	    std::vector<double> phi2nd_bins(nphi_bins_trig_effcy + 1);
	    for (int i = 0; i <= nphi_bins_trig_effcy; ++i) {
	        phi2nd_bins[i] = -pms.PI + (2.0 * pms.PI) * (static_cast<double>(i) / nphi_bins_trig_effcy);
	    }

	    hist_binning_map["phi2nd_bins" + ctr] = phi2nd_bins;
	}
}

// ---------- ----------
void RDFBasedHistFillingPbPb::BuildFlattenedTrgEffcyFilterToVarListMapExtra(){
    for (auto filter : trg_effcy_filters_1D_pre_sum_w_ctr_summing)    df_filter_to_var1D_list_map[filter] = single_muon_trig_effcy_var1Ds;
    for (auto filter : trg_effcy_filters_2D_3D_pre_sum_w_ctr_summing) df_filter_to_var2D_list_map[filter] = single_muon_trig_effcy_var2Ds;
    for (auto filter : trg_effcy_filters_2D_3D_pre_sum_w_ctr_summing) df_filter_to_var3D_list_map[filter] = single_muon_trig_effcy_var3Ds;

    for (auto filter : trg_effcy_filters_1D_pre_sum_w_musign_ctr_summing)    df_filter_to_var1D_list_map[filter] = single_muon_trig_effcy_var1Ds;
    for (auto filter : trg_effcy_filters_2D_3D_pre_sum_w_musign_ctr_summing) df_filter_to_var2D_list_map[filter] = single_muon_trig_effcy_var2Ds;
    for (auto filter : trg_effcy_filters_2D_3D_pre_sum_w_musign_ctr_summing) df_filter_to_var3D_list_map[filter] = single_muon_trig_effcy_var3Ds;

    for (auto filter : trg_effcy_filters_ctr_dep_1D_pre_sum) df_filter_to_var1D_list_map[filter] = std::vector<std::string>({"ctr"});
}

//--------- BUILD & FLATTERN PRE-SUM, POST-SUM, TO-BE-SUMMED LEVELS ---------
void RDFBasedHistFillingPbPb::TrigEffcyFiltersPrePostSumFlatteningExtra()
{
    // build post-sum levels, ctr dep
    TrigEffcyUtils::write_post_sum_levels(levels_trg_effcy_filters_ctr_dep_pre_sum,
                          levels_trg_effcy_to_be_summed,
                          levels_trg_effcy_filters_ctr_dep_post_sum);

    // flatten pre-sum & post-sum levels, ctr dep
    TrigEffcyUtils::flatten_levels(levels_trg_effcy_filters_ctr_dep_pre_sum,	trg_effcy_filters_ctr_dep_1D_pre_sum);
    TrigEffcyUtils::flatten_levels(levels_trg_effcy_filters_ctr_dep_post_sum, 	trg_effcy_filters_ctr_dep_1D_post_sum);

    // flatten pre-sum levels, with ctr (+ musign) summing
    TrigEffcyUtils::flatten_levels(levels_trg_effcy_filters_1D_pre_sum_w_musign_ctr_summing, trg_effcy_filters_1D_pre_sum_w_musign_ctr_summing);
    TrigEffcyUtils::flatten_levels(levels_trg_effcy_filters_2D_3D_pre_sum_w_musign_ctr_summing, trg_effcy_filters_2D_3D_pre_sum_w_musign_ctr_summing);

    TrigEffcyUtils::flatten_levels(levels_trg_effcy_filters_1D_pre_sum_w_ctr_summing, trg_effcy_filters_1D_pre_sum_w_ctr_summing);
    TrigEffcyUtils::flatten_levels(levels_trg_effcy_filters_2D_3D_pre_sum_w_ctr_summing, trg_effcy_filters_2D_3D_pre_sum_w_ctr_summing);

    // build post-sum levels, with ctr (+ musign) summing
    TrigEffcyUtils::write_post_sum_levels(levels_trg_effcy_filters_1D_pre_sum_w_musign_ctr_summing,
                          levels_trg_effcy_to_be_summed,
                          levels_trg_effcy_filters_1D_post_sum_w_musign_ctr_summing);

    TrigEffcyUtils::write_post_sum_levels(levels_trg_effcy_filters_2D_3D_pre_sum_w_musign_ctr_summing,
                          levels_trg_effcy_to_be_summed,
                          levels_trg_effcy_filters_2D_3D_post_sum_w_musign_ctr_summing);

    TrigEffcyUtils::write_post_sum_levels(levels_trg_effcy_filters_1D_pre_sum_w_ctr_summing,
                          levels_trg_effcy_to_be_summed,
                          levels_trg_effcy_filters_1D_post_sum_w_ctr_summing);

    TrigEffcyUtils::write_post_sum_levels(levels_trg_effcy_filters_2D_3D_pre_sum_w_ctr_summing,
                          levels_trg_effcy_to_be_summed,
                          levels_trg_effcy_filters_2D_3D_post_sum_w_ctr_summing);

    // flatten post-sum levels, with ctr (+ musign) summing
    TrigEffcyUtils::flatten_levels(levels_trg_effcy_filters_1D_post_sum_w_musign_ctr_summing, trg_effcy_filters_1D_post_sum_w_musign_ctr_summing);
    TrigEffcyUtils::flatten_levels(levels_trg_effcy_filters_2D_3D_post_sum_w_musign_ctr_summing, trg_effcy_filters_2D_3D_post_sum_w_musign_ctr_summing);

    TrigEffcyUtils::flatten_levels(levels_trg_effcy_filters_1D_post_sum_w_ctr_summing, trg_effcy_filters_1D_post_sum_w_ctr_summing);
    TrigEffcyUtils::flatten_levels(levels_trg_effcy_filters_2D_3D_post_sum_w_ctr_summing, trg_effcy_filters_2D_3D_post_sum_w_ctr_summing);
}

// ---------- create essential RDFs ----------
void RDFBasedHistFillingPbPb::CreateRDFs(){
	RDFBasedHistFillingBaseClass::CreateRDFs();
	
	// create ctr-divided RDFs
	for (std::string pair_sign : {"_ss", "_op"}){
		std::string df_name = "df" + pair_sign;
		ROOT::RDF::RNode& node = map_at_checked(df_map, df_name, Form("CreateRDFs: df_map.at(%s)", df_name.c_str()));
	    
	    for (int ictr = 0; ictr < nCtrBins; ictr++){

	        int ctr_bin_low_edge = ctr_bin_edges.at(ictr);
	        int ctr_bin_high_edge = ctr_bin_edges.at(ictr + 1);

	        auto ctr_in_current_bin = [this, ctr_bin_low_edge, ctr_bin_high_edge](int centrality) 
	        	{
	        		if (ctr_binning_version == "include_upc" && ctr_bin_high_edge == 100){ // most peripheral bin
	        			return centrality >= ctr_bin_low_edge || centrality == -1;
	        		} else{
	        			return centrality >= ctr_bin_low_edge && centrality < ctr_bin_high_edge;
	        		}
	        		
	        	};
			std::string ctr = "_ctr" + std::to_string(ctr_bin_low_edge) + "_" + std::to_string(ctr_bin_high_edge);

			df_map.emplace (df_name + ctr, 
							node.Filter(ctr_in_current_bin, {"avg_centrality"})
						    );
	    }
	}
}

void RDFBasedHistFillingPbPb::FillHistogramsSingleMuonEffcy(){
	if (trigger_mode == 1){ // trigger efficiencies: mu4_mu4noL1 | mu4 & 2mu4 | mu4
		FillHistogramsDimuTrigGivenMu4();
		FillHistogramsDimuTrigGivenMu4CtrDep();
	} else if (trigger_mode == 0){
		FillHistogramsMu4GivenMB();
		FillHistogramsMu4GivenMBCtrDep();
	}
}

void RDFBasedHistFillingPbPb::FillHistogramsDimuTrigGivenMu4(){
	try{
	   	for (std::string pair_sign : {"_ss", "_op"}){
	       	std::string df_name = "df" + pair_sign; // e.g, df_ss_mu1passmu4
			ROOT::RDF::RNode& node = map_at_checked(df_map, df_name, Form("FillHistogramsDimuTrigGivenMu4: df_map.at(%s)", df_name.c_str()));

	       	for (auto mu4sel : {"_mu1passmu4", "_mu2passmu4"}){ // mu4 selection
	           	std::string df_name = "df" + pair_sign + mu4sel; // e.g, df_ss_mu1passmu4

				std::string ind1st = (std::string(mu4sel) == "_mu1passmu4")? "1" : "2";
	           	std::string ind2nd = (std::string(mu4sel) == "_mu1passmu4")? "2" : "1";
				
				df_map.emplace(df_name, node.Filter("m" + ind1st + ".passmu4"));
				df_map.at(df_name) = df_map.at(df_name).Define("pt1st",	"m" + ind1st + ".pt");
				df_map.at(df_name) = df_map.at(df_name).Define("pt2nd",	"m" + ind2nd + ".pt");
				df_map.at(df_name) = df_map.at(df_name).Define("charge2nd","m" + ind2nd + ".charge");
				df_map.at(df_name) = df_map.at(df_name).Define("eta2nd",	"m" + ind2nd + ".eta");
				df_map.at(df_name) = df_map.at(df_name).Define("phi2nd",	"m" + ind2nd + ".phi");
				df_map.at(df_name) = df_map.at(df_name).Define("q_eta2nd",	"charge2nd * eta2nd");
	           	df_map.at(df_name) = df_map.at(df_name).Define("second_muon_good_acceptance",  "pt2nd >= 6 && ((eta2nd > 1.1 && eta2nd < 2.3) || (eta2nd > -2.3 && eta2nd < -1.2))");
				
				// ----------- no musign / ctr ----------- 

	            for (auto trg : {"_mu4", "_mu4_mu4noL1", "_2mu4", "_2mu4_AND_mu4_mu4noL1"}){
	            	std::string trg_filter = map_at_checked(trig_to_filter_str_map, trg, Form("trig_to_filter_str_map.at(%s)", trg));
	            	if (trg_filter.empty()) df_map.emplace(df_name + trg, map_at_checked(df_map, df_name, Form("FillHistogramsDimuTrigGivenMu4: df_map.at(%s)", df_name.c_str())));
	            	else 					df_map.emplace(df_name + trg, map_at_checked(df_map, df_name, Form("FillHistogramsDimuTrigGivenMu4: df_map.at(%s)", df_name.c_str())).Filter(trg_filter));
	               	
	               	std::string df_name_new = "df" + pair_sign + mu4sel + trg;
	               	df_map.emplace(df_name_new + "_sepr", map_at_checked(df_map, df_name_new, Form("FillHistogramsDimuTrigGivenMu4: df_map.at(%s)", df_name_new.c_str())).Filter("passSeparated"));
	                
	               	for (auto bias : trg_effcy_biases){ // additional selection / bias in data sample
	                   	std::string df_name_new = "df" + pair_sign + mu4sel + trg + bias;
	                   	std::string filter = df_name_new.substr(2);

	                   	FillHistogramsSingleDataFrame(filter, map_at_checked(df_map, df_name_new, Form("FillHistogramsDimuTrigGivenMu4: df_map.at(%s)", df_name.c_str())), false); // do not write the sub-dataframe histograms in output file
	                }
	            }
				// ----------- with musign, no ctr ----------- 
				df_map.emplace(df_name + "_sign1", map_at_checked(df_map, df_name, Form("FillHistogramsDimuTrigGivenMu4: df_map.at(%s)", df_name.c_str())).Filter("charge2nd > 0")); // e.g, df_ss_mu1passmu4_sign1
				df_map.emplace(df_name + "_sign2", map_at_checked(df_map, df_name, Form("FillHistogramsDimuTrigGivenMu4: df_map.at(%s)", df_name.c_str())).Filter("charge2nd < 0"));

				for (auto mu_sign : {"_sign1", "_sign2"}){
	               	std::string df_name = "df" + pair_sign + mu4sel + mu_sign;

	               	for (auto trg : {"_mu4", "_mu4_mu4noL1", "_2mu4", "_2mu4_AND_mu4_mu4noL1"}){
		            	std::string trg_filter = map_at_checked(trig_to_filter_str_map, trg, Form("trig_to_filter_str_map.at(%s)", trg));
	            		if (trg_filter.empty()) df_map.emplace(df_name + trg, map_at_checked(df_map, df_name, Form("FillHistogramsDimuTrigGivenMu4: df_map.at(%s)", df_name.c_str())));
		            	else 					df_map.emplace(df_name + trg, map_at_checked(df_map, df_name, Form("FillHistogramsDimuTrigGivenMu4: df_map.at(%s)", df_name.c_str())).Filter(trg_filter));

	                   	std::string df_name_new = "df" + pair_sign + mu4sel + mu_sign + trg;
	                   	df_map.emplace(df_name_new + "_sepr", map_at_checked(df_map, df_name_new, Form("FillHistogramsDimuTrigGivenMu4: df_map.at(%s)", df_name_new.c_str())).Filter("passSeparated"));
	                    
	                   	for (auto bias : trg_effcy_biases){ // additional selection / bias in data sample
	                       	std::string df_name_new = "df" + pair_sign + mu4sel + mu_sign + trg + bias;
	                       	std::string filter = df_name_new.substr(2);

	                       	FillHistogramsSingleDataFrame(filter, map_at_checked(df_map, df_name_new, Form("FillHistogramsDimuTrigGivenMu4: df_map.at(%s)", df_name_new.c_str())), false); // do not write the sub-dataframe histograms in output file
	                    }
	                }
				}

				// ----------- with musign & ctr ----------- 
			   	for (int ictr = 0; ictr < nCtrBins; ictr++){

			       	int ctr_bin_low_edge = ctr_bin_edges.at(ictr);
			       	int ctr_bin_high_edge = ctr_bin_edges.at(ictr + 1);

					std::string ctr = "_ctr" + std::to_string(ctr_bin_low_edge) + "_" + std::to_string(ctr_bin_high_edge);

					df_map.emplace (df_name + ctr, 
									map_at_checked(df_map, df_name, Form("FillHistogramsDimuTrigGivenMu4: df_map.at(%s)", df_name.c_str())).Filter(
			            					[ctr_bin_low_edge, ctr_bin_high_edge](int centrality) {return centrality >= ctr_bin_low_edge && centrality < ctr_bin_high_edge;},
			            					{"avg_centrality"})
								    );

					df_map.emplace(df_name + ctr + "_sign1", map_at_checked(df_map, df_name + ctr, Form("FillHistogramsDimuTrigGivenMu4: df_map.at(%s)", (df_name + ctr).c_str())).Filter("charge2nd > 0")); // e.g, df_ss_mu1passmu4_sign1
					df_map.emplace(df_name + ctr + "_sign2", map_at_checked(df_map, df_name + ctr, Form("FillHistogramsDimuTrigGivenMu4: df_map.at(%s)", (df_name + ctr).c_str())).Filter("charge2nd < 0"));

					for (auto mu_sign : {"_sign1", "_sign2"}){
	                   	std::string df_name = "df" + pair_sign + mu4sel + ctr + mu_sign;

	                   	for (auto trg : {"_mu4", "_mu4_mu4noL1", "_2mu4", "_2mu4_AND_mu4_mu4noL1"}){
	                    	std::string trg_filter = map_at_checked(trig_to_filter_str_map, trg, Form("trig_to_filter_str_map.at(%s)", trg));
	            			if (trg_filter.empty()) df_map.emplace(df_name + trg, map_at_checked(df_map, df_name, Form("FillHistogramsDimuTrigGivenMu4: df_map.at(%s)", df_name.c_str())));
		            		else 		            df_map.emplace(df_name + trg, map_at_checked(df_map, df_name, Form("FillHistogramsDimuTrigGivenMu4: df_map.at(%s)", df_name.c_str())).Filter(trg_filter));

	                       	std::string df_name_new = "df" + pair_sign + mu4sel + ctr + mu_sign + trg;
	                       	df_map.emplace(df_name_new + "_sepr", map_at_checked(df_map, df_name_new, Form("FillHistogramsDimuTrigGivenMu4: df_map.at(%s)", df_name_new.c_str())).Filter("passSeparated"));
	                        
	                       	for (auto bias : trg_effcy_biases){ // additional selection / bias in data sample
	                           	std::string df_name_new = "df" + pair_sign + mu4sel + ctr + mu_sign + trg + bias;
	                           	std::string filter = df_name_new.substr(2);

	                           	FillHistogramsSingleDataFrame(filter, map_at_checked(df_map, df_name_new, Form("FillHistogramsDimuTrigGivenMu4: df_map.at(%s)", df_name_new.c_str())), false); // do not write the sub-dataframe histograms in output file
	                        }
	                    }
					}
				} // end loop ctr
			} // end loop mu4 selection
	   	} // end loop pair sign
	} catch(const std::out_of_range& e){
		std::cerr << "FillHistogramsDimuTrigGivenMu4:: out_of_range exception caught: " << e.what() << std::endl;
    } catch (const std::runtime_error& e) {
       	std::cerr << "FillHistogramsDimuTrigGivenMu4:: RDF runtime error: " << e.what() << std::endl;
    }
}

void RDFBasedHistFillingPbPb::FillHistogramsDimuTrigGivenMu4CtrDep(){
	try{
	   	for (std::string pair_sign : {"_ss", "_op"}){

	       	for (auto mu4sel : {"_mu1passmu4", "_mu2passmu4"}){ // mu4 selection
	           	std::string df_name = "df" + pair_sign + mu4sel; // e.g, df_ss_mu1passmu4
				ROOT::RDF::RNode& node = map_at_checked(df_map, df_name, Form("FillHistogramsDimuTrigGivenMu4CtrDep: df_map.at(%s)", df_name.c_str()));
				
				for (int ipt = 0; ipt < pT_bins_edges_for_trg_effcy_ctr_dep.size()-1; ipt++){

				   	int pt_bin_low_edge = pT_bins_edges_for_trg_effcy_ctr_dep.at(ipt);
				   	int pt_bin_high_edge = pT_bins_edges_for_trg_effcy_ctr_dep.at(ipt + 1);

					std::string pt_bin = "_pt" + std::to_string(pt_bin_low_edge) + "_" + std::to_string(pt_bin_high_edge);
					
					df_map.emplace (df_name + pt_bin, 
									map_at_checked(df_map, df_name, Form("FillHistogramsDimuTrigGivenMu4CtrDep: df_map.at(%s)", df_name.c_str()))
									.Filter(
				        					[pt_bin_low_edge, pt_bin_high_edge](float pt2nd) 
				        					{return pt2nd >= pt_bin_low_edge && pt2nd < pt_bin_high_edge;},
				        					{"pt2nd"})
								    );

					for (auto q_eta_pair : q_eta_proj_ranges){
		               	
						std::string q_eta_bin = "_q_eta_" + pairToSuffix(q_eta_pair);

						df_map.emplace(	df_name + pt_bin + q_eta_bin, 
										map_at_checked(df_map, df_name + pt_bin, 
											Form("FillHistogramsDimuTrigGivenMu4CtrDep: df_map.at(%s)", 
												(df_name + pt_bin).c_str()))
										.Filter([q_eta_pair](float q_eta2nd) 
											{return q_eta2nd >= q_eta_pair.first && q_eta2nd < q_eta_pair.second;}
											, {"q_eta2nd"})
										);
						
		               	std::string df_name_new = "df" + pair_sign + mu4sel + pt_bin + q_eta_bin;

		               	for (auto trg : {"_mu4", "_mu4_mu4noL1", "_2mu4", "_2mu4_AND_mu4_mu4noL1"}){
			            	std::string trg_filter = map_at_checked(trig_to_filter_str_map, trg, Form("trig_to_filter_str_map.at(%s)", trg));
	            			if (trg_filter.empty()) df_map.emplace(df_name_new + trg, map_at_checked(df_map, df_name_new, Form("FillHistogramsDimuTrigGivenMu4CtrDep: df_map.at(%s)", df_name_new.c_str())));
			            	else 			        df_map.emplace(df_name_new + trg, map_at_checked(df_map, df_name_new, Form("FillHistogramsDimuTrigGivenMu4CtrDep: df_map.at(%s)", df_name_new.c_str())).Filter(trg_filter));

		                   	std::string df_name_final = "df" + pair_sign + mu4sel + pt_bin + q_eta_bin + trg;
		                   	df_map.emplace(df_name_final + "_sepr", map_at_checked(df_map, df_name_final, Form("FillHistogramsDimuTrigGivenMu4CtrDep: df_map.at(%s)", df_name_final.c_str())).Filter("passSeparated"));
		                    
		                   	for (auto bias : trg_effcy_biases){ // additional selection / bias in data sample
		                       	std::string df_name_final = "df" + pair_sign + mu4sel + pt_bin + q_eta_bin + trg + bias;
		                       	std::string filter = df_name_final.substr(2);

		                       	FillHistogramsSingleDataFrame(filter, map_at_checked(df_map, df_name_final, Form("FillHistogramsDimuTrigGivenMu4CtrDep: df_map.at(%s)", df_name_final.c_str())), false); // do not write the sub-dataframe histograms in output file
		                    }
		                }
					} // end loop 2nd-muon q * eta
				} // end loop 2nd-muon pt
			} // end loop mu4 selection
	   	} // end loop pair sign
	} catch(const std::out_of_range& e){
		std::cerr << "FillHistogramsDimuTrigGivenMu4CtrDep:: out_of_range exception caught: " << e.what() << std::endl;
    } catch (const std::runtime_error& e) {
       	std::cerr << "FillHistogramsDimuTrigGivenMu4CtrDep:: RDF runtime error: " << e.what() << std::endl;
    }
}

void RDFBasedHistFillingPbPb::FillHistogramsMu4GivenMB(){}

void RDFBasedHistFillingPbPb::FillHistogramsMu4GivenMBCtrDep(){}

void RDFBasedHistFillingPbPb::FillTrigEffcyHistsInvWeightedbySingleMuonEffcies(){}

//--------- DATA HIST POST PROCESSING ---------
void RDFBasedHistFillingPbPb::HistPostProcessPbPb(){
    if ((trigger_mode == 0 || trigger_mode == 1) && hist_filling_cycle == generic){
        MakeAndWriteSingleMuonTrigEffCtrGraphs();
    }
}

void RDFBasedHistFillingPbPb::SumSingleMuonTrigEffHistsPbPb(){
    // sum centrality-dependent trigger-efficiency histograms
    TrigEffcyUtils::SumTrigEffHistsGeneric<TH1D, std::string>(
        std::vector<std::string>({"ctr"}),
        trg_effcy_filters_ctr_dep_1D_post_sum,
        trg_effcy_filters_to_be_summed,
        hist1D_map,
        [](const std::string& var) {
            return "h_" + var;
        }
    );

    // 1D, with musign + ctr binning
    TrigEffcyUtils::SumTrigEffHistsGeneric<TH1D, std::string>(
        single_muon_trig_effcy_var1Ds,
        trg_effcy_filters_1D_post_sum_w_musign_ctr_summing,
        trg_effcy_filters_to_be_summed,
        hist1D_map,
        [](const std::string& var) {
            return "h_" + var;
        }
    );

    // 1D, with ctr binning
    TrigEffcyUtils::SumTrigEffHistsGeneric<TH1D, std::string>(
        single_muon_trig_effcy_var1Ds,
        trg_effcy_filters_1D_post_sum_w_ctr_summing,
        trg_effcy_filters_to_be_summed,
        hist1D_map,
        [](const std::string& var) {
            return "h_" + var;
        }
    );

    // 2D, with musign + ctr binning
    TrigEffcyUtils::SumTrigEffHistsGeneric<TH2D, std::array<std::string,2>>(
        single_muon_trig_effcy_var2Ds,
        trg_effcy_filters_2D_3D_post_sum_w_musign_ctr_summing,
        trg_effcy_filters_to_be_summed,
        hist2D_map,
        [](const std::array<std::string,2>& vars) {
            const std::string& varx = vars[0];
            const std::string& vary = vars[1];
            return "h_" + vary + "_vs_" + varx;
        }
    );

    // 2D, with ctr binning
    TrigEffcyUtils::SumTrigEffHistsGeneric<TH2D, std::array<std::string,2>>(
        single_muon_trig_effcy_var2Ds,
        trg_effcy_filters_2D_3D_post_sum_w_ctr_summing,
        trg_effcy_filters_to_be_summed,
        hist2D_map,
        [](const std::array<std::string,2>& vars) {
            const std::string& varx = vars[0];
            const std::string& vary = vars[1];
            return "h_" + vary + "_vs_" + varx;
        }
    );

    // 3D, with musign + ctr binning
    TrigEffcyUtils::SumTrigEffHistsGeneric<TH3D, std::array<std::string,3>>(
        single_muon_trig_effcy_var3Ds,
        trg_effcy_filters_2D_3D_post_sum_w_musign_ctr_summing,
        trg_effcy_filters_to_be_summed,
        hist3D_map,
        [](const std::array<std::string,3>& vars) {
            const std::string& varx = vars[0];
            const std::string& vary = vars[1];
            const std::string& varz = vars[2];
            return "h_" + varz + "_vs_" + vary + "_vs_" + varx;
        }
    );

    // 3D, with ctr binning
    TrigEffcyUtils::SumTrigEffHistsGeneric<TH3D, std::array<std::string,3>>(
        single_muon_trig_effcy_var3Ds,
        trg_effcy_filters_2D_3D_post_sum_w_ctr_summing,
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

void RDFBasedHistFillingPbPb::CalculateSingleMuonTrigEffcyRatios(){
    CalculateSingleMuonTrigEffcyRatiosHelper(ctr_musign_categories);
}

void RDFBasedHistFillingPbPb::MakeAndWriteSingleMuonTrigEffPtGraphs(){
	if (useCoarseQEtaBin)	MakeAndWriteSingleMuonTrigEffPtGraphsHelper({});
	else{
		MakeAndWriteSingleMuonTrigEffPtGraphsHelper(ctr_musign_categories);
		MakeAndWriteSingleMuonTrigEffPtGraphsHelper(musigns);
	} 					
}

void RDFBasedHistFillingPbPb::MakeAndWriteSingleMuonTrigEffCtrGraphs(){

	for (auto pT_bin : pT_bins_for_trg_effcy_ctr_dep){
		for (auto q_eta_bin : q_eta_ranges_str){
	        TH1D* h_ctr_mu4_mu4noL1_sepr  = map_at_checked(hist1D_map, 	Form("h_ctr%s%s_mu4_mu4noL1_sepr", pT_bin.c_str(), q_eta_bin.c_str()),    
	        															Form("MakeAndWriteSingleMuonTrigEffCtrGraphs: hist1D_map.at(h_ctr%s%s_mu4_mu4noL1_sepr)", pT_bin.c_str(), q_eta_bin.c_str()));
	        TH1D* h_ctr_mu4_sepr          = map_at_checked(hist1D_map, 	Form("h_ctr%s%s_mu4_sepr", pT_bin.c_str(), q_eta_bin.c_str()),            
	        															Form("MakeAndWriteSingleMuonTrigEffCtrGraphs: hist1D_map.at(h_ctr%s%s_mu4_sepr)", pT_bin.c_str(), q_eta_bin.c_str()));
    	    TH1D* h_ctr_2mu4_sepr         = map_at_checked(hist1D_map, 	Form("h_ctr%s%s_2mu4_sepr", pT_bin.c_str(), q_eta_bin.c_str()),           
    	    															Form("MakeAndWriteSingleMuonTrigEffCtrGraphs: hist1D_map.at(h_ctr%s%s_2mu4_sepr)", pT_bin.c_str(), q_eta_bin.c_str()));
        	TH1D* h_ctr_2mu4_AND_mu4_mu4noL1_sepr   = map_at_checked(hist1D_map, 	Form("h_ctr%s%s_2mu4_AND_mu4_mu4noL1_sepr", pT_bin.c_str(), q_eta_bin.c_str()),           
        																			Form("MakeAndWriteSingleMuonTrigEffCtrGraphs: hist1D_map.at(h_ctr%s%s_2mu4_AND_mu4_mu4noL1_sepr)", pT_bin.c_str(), q_eta_bin.c_str()));

		    TrigEffcyUtils::divide_and_write(h_ctr_mu4_mu4noL1_sepr,           h_ctr_mu4_sepr, &graph_map);
		    TrigEffcyUtils::divide_and_write(h_ctr_2mu4_sepr,                  h_ctr_mu4_sepr, &graph_map);
		    TrigEffcyUtils::divide_and_write(h_ctr_2mu4_AND_mu4_mu4noL1_sepr,  h_ctr_mu4_sepr, &graph_map);
		}
	}
}

void RDFBasedHistFillingPbPb::OpenEffcyPtFitFile(){}

void RDFBasedHistFillingPbPb::MakeAndWriteDRTrigEffGraphs(){}


AxisInfo RDFBasedHistFillingPbPb::GetAxisInfo(const var1D& v,
                                              const std::string& filter) const {
    // Start from base-class behavior (global / ss / op)
    AxisInfo b = RDFBasedHistFillingBaseClass::GetAxisInfo(v, filter);

    // Now overlay PbPb centrality binning if applicable
    for (int ictr = 0; ictr < ctr_bins.size(); ictr++) {
    	auto ctr = ctr_bins.at(ictr);
        if (filter.find(ctr) == std::string::npos) continue; // current ctr not substring of filter

        auto it_bins = v.bins_ctr.find(ctr);
        if (it_bins == v.bins_ctr.end()) break;          // no special binning for this ctr

        const auto& edges = it_bins->second;
        if (edges.empty()) break;                        // treat empty as "no override"

        auto it_nbins = v.nbins_ctr.find(ctr);
        int nb_ctr = 0;
        if (it_nbins != v.nbins_ctr.end())
            nb_ctr = it_nbins->second;
        else
            nb_ctr = static_cast<int>(edges.size()) - 1;    // fallback if you didn't fill nbins_ctr

        if (nb_ctr <= 0 || static_cast<int>(edges.size()) != nb_ctr + 1) {
            std::ostringstream oss;
            oss << "GetAxisInfo (PbPb): inconsistent centrality binning for variable '"
                << v.name << "', ctr='" << ctr << "'";
            throw std::runtime_error(oss.str());
        }

        // override axis info with centrality-dependent binning
        b.nbins    = nb_ctr;
        b.bin_edges = edges.data();
        b.min      = 0.;
        b.max      = 0.;

        // we assume only one ctr suffix per filter, so we can break
        break;
    }

    // second loop - rebinning for the (nbins, min, max) case
    for (int ictr = 0; ictr < ctr_bins.size(); ictr++) {
    	auto ctr = ctr_bins.at(ictr);
        if (filter.find(ctr) == std::string::npos) continue; // current ctr not substring of filter

	    if (!(b.min == 0 && b.max == 0)){ // no special ctr-dependent binning found, (nbins, min, max) format used --> apply rebin
	        if (std::find(var1D_list_no_ctr_rebin.begin(), var1D_list_no_ctr_rebin.end(), v.name)
	        		 != var1D_list_no_ctr_rebin.end()) 	break; // do not apply rebin

		    if (var1D_excep_ctr_rebin_factor_map.find(v.name)!= var1D_excep_ctr_rebin_factor_map.end()) { // special rebin required
		        if (var1D_excep_ctr_rebin_factor_map.at(v.name).size() == ctr_bins.size()
		        	&&  var1D_excep_ctr_rebin_factor_map.at(v.name).at(ictr) > 0) { // sanity check
		    		b.nbins /= var1D_excep_ctr_rebin_factor_map.at(v.name).at(ictr);
		    	}
			} else{ // apply automatic rebinning
				b.nbins /= ctr_rebin_scale_factors.at(ictr);
			}
	    }
	}

    return b;
}

void RDFBasedHistFillingPbPb::ReadVar1DJson() {
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
            // ---------------- PbPb-specific string binning logic ----------------
            std::string binName = jb.get<std::string>();
            if (binName.empty()) {
                delete v;
                ThrowMissingField("binning", histName);
            }

            auto it     = hist_binning_map.find(binName);
            auto it_ss  = hist_binning_map.find(binName + "_ss");
            auto it_op  = hist_binning_map.find(binName + "_op");

			bool has_inclusive = (it != hist_binning_map.end());
            bool has_ssop      = (it_ss != hist_binning_map.end() &&
                                  it_op != hist_binning_map.end());

            // ss/op override (same logic as base)
            if (has_ssop) {
                v->bins_ss  = it_ss->second;
                v->nbins_ss = static_cast<int>(v->bins_ss.size()) - 1;
                v->bins_op  = it_op->second;
                v->nbins_op = static_cast<int>(v->bins_op.size()) - 1;
            } else if (has_inclusive) {
                v->bins  = it->second;
                v->nbins = static_cast<int>(v->bins.size()) - 1;
            }

            // Centrality-dependent binning: look for binName + ctr
            bool has_ctr_binning = false;
            for (const auto& ctr : ctr_bins) { // e.g. "_ctr0_5"
                auto it_ctr = hist_binning_map.find(binName + ctr);
                if (it_ctr == hist_binning_map.end()) continue;

                const auto& edges = it_ctr->second;
                if (edges.empty()) {
                    std::ostringstream oss;
                    oss << "ReadVar1DJson (PbPb): centrality binning for '"
                        << binName << "', ctr='" << ctr
                        << "' is empty (hist_name='" << histName << "')";
                    delete v;
                    throw std::runtime_error(oss.str());
                }

                v->bins_ctr[ctr] = edges;
                v->nbins_ctr[ctr] = static_cast<int>(edges.size()) - 1;
                has_ctr_binning = true;
            }

            if (!has_inclusive && !has_ssop && !has_ctr_binning) {
                delete v;
                std::ostringstream oss;
                oss << "ReadVar1DJson (PbPb): binning name '" << binName
                    << "' (referenced by hist_name='" << histName
                    << "') not found in hist_binning_map (including centralities)";
                throw std::runtime_error(oss.str());
            }

        } else if (jb.is_object()) {
            // binning as {nbins, min, max} --- same as base
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

        // Final validation using updated isValid()
        if (!v->isValid()) {
            std::ostringstream oss;
            oss << "ReadVar1DJson (PbPb): invalid binning specification for hist_name='"
                << histName << "'";
            delete v;
            throw std::runtime_error(oss.str());
        }

        // All good, store pointer in dictionary
        var1D_dict[v->name] = v;
    }
}
