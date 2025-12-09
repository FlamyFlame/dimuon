#include "RDFBasedHistFilling.h"

void RDFBasedHistFillingPbPb::Initialize(){
	// input_files.push_back("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_2024/muon_pairs_pbpb_2024_single_mu4.root");
	// output_file = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_2024/histograms_real_pairs_pbpb_2024_single_mu4.root";
	input_files.push_back("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_2023/muon_pairs_pbpb_2023_single_mu4.root");
	output_file = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_2023/histograms_real_pairs_pbpb_2023_single_mu4.root";
	
    ctr_binning_map["default"] = {0,5,10,20,30,50,80};
    ctr_binning_file_suffix_map["default"] = "";

	ctr_bin_edges = ctr_binning_map[ctr_binning_verion];
	ctr_bin_edges_double = std::vector<double>(ctr_bin_edges.begin(), ctr_bin_edges.end());
	nCtrBins = ctr_bin_edges.size() - 1;

    for (int ictr = 0; ictr < nCtrBins; ictr++){
        int ctr_bin_low_edge = ctr_bin_edges.at(ictr);
        int ctr_bin_high_edge = ctr_bin_edges.at(ictr + 1);
        ctr_bins.push_back("_ctr" + std::to_string(ctr_bin_low_edge) + "_" + std::to_string(ctr_bin_high_edge));
	}

    levels_trg_effcy_filters_1D_pre_sum = {{"_ss", "_op"}, // add a level of ctr bins for Pb+Pb
                                        {"_mu1passmu4", "_mu2passmu4"},
                                        ctr_bins,
                                        {"_sign1", "_sign2"},
                                        {"_mu4", "_mu4_mu4noL1", "_2mu4"}, 
                                        {"", "_sepr"}};

    levels_trg_effcy_filters_2D_3D_pre_sum = {{"_ss", "_op"}, // add a level of ctr bins for Pb+Pb
                                        {"_mu1passmu4", "_mu2passmu4"},
                                        ctr_bins,
                                        {"_sign1", "_sign2"},
                                        {"_mu4", "_mu4_mu4noL1", "_2mu4"}, 
                                        {"", "_sepr"}};
	RDFBasedHistFilling::Initialize();
}

void RDFBasedHistFillingPbPb::CreateRDFs(){
	RDFBasedHistFilling::CreateRDFs();
	
	// create ctr-divided RDFs
	for (std::string pair_sign : {"_ss", "_op"}){
		std::string df_name = "df" + pair_sign;
		ROOT::RDF::RNode& node = map_at_checked(df_map, df_name, Form("FillHistograms: df_map.at(%s)", df_name.c_str()));
	    
	    for (int ictr = 0; ictr < nCtrBins; ictr++){

	        int ctr_bin_low_edge = ctr_bin_edges.at(ictr);
	        int ctr_bin_high_edge = ctr_bin_edges.at(ictr + 1);

			std::string ctr = "_ctr" + std::to_string(ctr_bin_low_edge) + "_" + std::to_string(ctr_bin_high_edge);

			df_map.emplace (df_name + ctr, 
							node.Filter(
	            					[ctr_bin_low_edge, ctr_bin_high_edge](int centrality) {return centrality >= ctr_bin_low_edge && centrality < ctr_bin_high_edge;},
	            					{"avg_centrality"})
						    );
	    }
	}
}

void RDFBasedHistFillingPbPb::FillHistograms(){

    std::cout << "Calling FillHistograms" << std::endl;

    if (output_non_trig_effcy_hists){ // non-trigger-efficiency histograms
        for (std::string sign : {"_ss", "_op"}){
	        for (std::string ctr : ctr_bins){
	            FillHistogramsSingleDataFrame(sign, df_map.at(sign));
	            FillHistogramsSingleDataFrame(sign, "_jacobian_corrected", df_map.at(sign));
	            // FillHistogramsSingleDataFrame(sign + "_wgapcut", df_map.at(sign));
	            // FillHistogramsSingleDataFrame(sign + "_wgapcut", "_jacobian_corrected", df_map.at(sign));
	        }
	    }
    }

	if (trigger_mode == 1){ // trigger efficiencies: mu4_mu4noL1 | mu4 & 2mu4 | mu4
		try{
		    for (std::string pair_sign : {"_ss", "_op"}){
	            std::string df_name = "df" + pair_sign; // e.g, df_ss_mu1passmu4
				ROOT::RDF::RNode& node = map_at_checked(df_map, df_name, Form("FillHistograms: df_map.at(%s)", df_name.c_str()));

	            for (auto mu4sel : {"_mu1passmu4", "_mu2passmu4"}){ // mu4 selection
	                std::string df_name = "df" + pair_sign + mu4sel; // e.g, df_ss_mu1passmu4

					std::string ind1st = (std::string(mu4sel) == "_mu1passmu4")? "1" : "2";
	                std::string ind2nd = (std::string(mu4sel) == "_mu1passmu4")? "2" : "1";
					
					df_map.emplace(df_name, node.Filter("mu" + ind1st + "PassSingle"));
					df_map.at(df_name) = df_map.at(df_name).Define("pt1st",	"m" + ind1st + ".pt");
					df_map.at(df_name) = df_map.at(df_name).Define("pt2nd",	"m" + ind2nd + ".pt");
					df_map.at(df_name) = df_map.at(df_name).Define("charge2nd","m" + ind2nd + ".charge");
					df_map.at(df_name) = df_map.at(df_name).Define("eta2nd",	"m" + ind2nd + ".eta");
					df_map.at(df_name) = df_map.at(df_name).Define("phi2nd",	"m" + ind2nd + ".phi");
					df_map.at(df_name) = df_map.at(df_name).Define("q_eta2nd",	"charge2nd * eta2nd");
	                df_map.at(df_name) = df_map.at(df_name).Define("second_muon_good_acceptance",  "pt2nd >= 6 && ((eta2nd > 1.1 && eta2nd < 2.3) || (eta2nd > -2.3 && eta2nd < -1.2))");
										
				    for (int ictr = 0; ictr < nCtrBins; ictr++){

				        int ctr_bin_low_edge = ctr_bin_edges.at(ictr);
				        int ctr_bin_high_edge = ctr_bin_edges.at(ictr + 1);

						std::string ctr = "_ctr" + std::to_string(ctr_bin_low_edge) + "_" + std::to_string(ctr_bin_high_edge);

						df_map.emplace (df_name + ctr, 
										map_at_checked(df_map, df_name, Form("FillHistograms: df_map.at(%s)", df_name.c_str())).Filter(
				            					[ctr_bin_low_edge, ctr_bin_high_edge](int centrality) {return centrality >= ctr_bin_low_edge && centrality < ctr_bin_high_edge;},
				            					{"avg_centrality"})
									    );

						df_map.emplace(df_name + ctr + "_sign1", map_at_checked(df_map, df_name + ctr, Form("FillHistograms: df_map.at(%s)", (df_name + ctr).c_str())).Filter("charge2nd > 0")); // e.g, df_ss_mu1passmu4_sign1
						df_map.emplace(df_name + ctr + "_sign2", map_at_checked(df_map, df_name + ctr, Form("FillHistograms: df_map.at(%s)", (df_name + ctr).c_str())).Filter("charge2nd < 0"));

						for (auto mu_sign : {"_sign1", "_sign2"}){
	                        std::string df_name = "df" + pair_sign + mu4sel + ctr + mu_sign;
	                        df_map.emplace(df_name + "_mu4", map_at_checked(df_map, df_name, Form("FillHistograms: df_map.at(%s)", df_name.c_str())));
	                        df_map.emplace(df_name + "_mu4_mu4noL1", map_at_checked(df_map, df_name, Form("FillHistograms: df_map.at(%s)", df_name.c_str())).Filter("passmu4mu4noL1"));
	                        df_map.emplace(df_name + "_2mu4", map_at_checked(df_map, df_name, Form("FillHistograms: df_map.at(%s)", df_name.c_str())).Filter("pass2mu4"));

	                        for (auto trg : {"_mu4", "_mu4_mu4noL1", "_2mu4"}){
	                            std::string df_name = "df" + pair_sign + mu4sel + ctr + mu_sign + trg;
	                            df_map.emplace(df_name + "_sepr", map_at_checked(df_map, df_name, Form("FillHistograms: df_map.at(%s)", df_name.c_str())).Filter("passSeparated"));
	                            
	                            for (auto bias : {"", "_sepr"}){ // additional selection / bias in data sample
	                                std::string df_name = "df" + pair_sign + mu4sel + ctr + mu_sign + trg + bias;
	                                std::string filter = df_name.substr(2);

	                                FillHistogramsSingleDataFrame(filter, map_at_checked(df_map, df_name, Form("FillHistograms: df_map.at(%s)", df_name.c_str())), false); // do not write the sub-dataframe histograms in output file
	                            }
	                        }
						}
					} // end loop ctr
				} // end loop mu4 selection
	        } // end loop pair sign
		} catch(const std::out_of_range& e){
			std::cerr << "FillHistograms:: out_of_range exception caught: " << e.what() << std::endl;
        } catch (const std::runtime_error& e) {
            std::cerr << "FillHistograms:: RDF runtime error: " << e.what() << std::endl;
        }
	}
}


void RDFBasedHistFillingPbPb::HistPostProcess(){
	// RDFBasedHistFilling::HistPostProcess();

    // if (trigger_mode == 0 || trigger_mode == 1){
    //     if (hist_filling_cycle == generic){
    //         SumSingleMuonTrigEffHists();
    //         MakeAndWriteSingleMuonPtTrigEffGraphs();
    //         CalculateSingleMuonTrigEffcyRatios();           
    //     } else{
    //         MakeAndWriteDRTrigEffGraphs();
    //     }
    // }
}

void RDFBasedHistFillingPbPb::MakeAndWriteSingleMuonPtTrigEffGraphs(){
	std::cout << "Calling MakeAndWriteSingleMuonPtTrigEffGraphs" << std::endl;

    // helper for projection, making & writing of TEfficiency graphs
    auto proj_make_and_write = [&](TH2* hNum2D, TH2* hDen2D, bool projy = true, int firstbin = 1, int lastbin = -1, std::string proj_range_str = ""){
        if (!hNum2D || !hDen2D) return (TGraphAsymmErrors*)nullptr;

        // suffix that captures projection axis & range
        std::string proj_suffix = projy? "_py" : "_px";
        if (proj_range_str != "") proj_suffix += "_" + proj_range_str;

        std::unique_ptr<TH1> hNum1D(
            projy ? hNum2D->ProjectionY(Form("%s%s", hNum2D->GetName(), proj_suffix.c_str()), firstbin, lastbin, "e")
                  : hNum2D->ProjectionX(Form("%s%s", hNum2D->GetName(), proj_suffix.c_str()), firstbin, lastbin, "e")
        );

        std::unique_ptr<TH1> hDen1D(
            projy ? hDen2D->ProjectionY(Form("%s%s", hDen2D->GetName(), proj_suffix.c_str()), firstbin, lastbin, "e")
                  : hDen2D->ProjectionX(Form("%s%s", hDen2D->GetName(), proj_suffix.c_str()), firstbin, lastbin, "e")
        );

        auto g = new TGraphAsymmErrors();
        g->BayesDivide(hNum1D.get(), hDen1D.get());

        // Name: "h_<...>" → "g_<...>" from the *numerator* histo name
        std::string n = hNum1D->GetName() ? hNum1D->GetName() : "graph";
        
        n += "_divided";
        if (n.rfind("h_", 0) == 0) n.replace(0, 2, "g_"); else n = "g_" + n;
        g->SetName(n.c_str());

        g->Write(g->GetName(), TObject::kOverwrite);
        return g;
    };

    // q-eta bins
    std::vector<double> eta_bins_trig_effcy = ParamsSet::makeEtaTrigEffcyBinning();
    
    for (std::string sign : {"_sign1", "_sign2"}){
	    for (std::string ctr : ctr_bins){
	        TH2D* h_pt2nd_vs_q_eta2nd_mu4_mu4noL1_sepr  = map_at_checked(hist2D_map, "h_pt2nd_vs_q_eta2nd" + ctr + sign + "_mu4_mu4noL1_sepr",    Form("MakeAndWriteSingleMuonPtTrigEffGraphs: hist2D_map.at(%s)", ("h_pt2nd_vs_q_eta2nd" + ctr + sign + "_mu4_mu4noL1_sepr").c_str()));
	        TH2D* h_pt2nd_vs_q_eta2nd_mu4_sepr          = map_at_checked(hist2D_map, "h_pt2nd_vs_q_eta2nd" + ctr + sign + "_mu4_sepr",            Form("MakeAndWriteSingleMuonPtTrigEffGraphs: hist2D_map.at(%s)", ("h_pt2nd_vs_q_eta2nd" + ctr + sign + "_mu4_sepr").c_str()));
	        TH2D* h_pt2nd_vs_q_eta2nd_2mu4_sepr         = map_at_checked(hist2D_map, "h_pt2nd_vs_q_eta2nd" + ctr + sign + "_2mu4_sepr",           Form("MakeAndWriteSingleMuonPtTrigEffGraphs: hist2D_map.at(%s)", ("h_pt2nd_vs_q_eta2nd" + ctr + sign + "_2mu4_sepr").c_str()));

	        proj_make_and_write(h_pt2nd_vs_q_eta2nd_mu4_mu4noL1_sepr,     h_pt2nd_vs_q_eta2nd_mu4_sepr);
	        proj_make_and_write(h_pt2nd_vs_q_eta2nd_2mu4_sepr,            h_pt2nd_vs_q_eta2nd_mu4_sepr);

	        for (auto range : q_eta_proj_ranges_for_single_muon_effcy_pT_fitting){
	            int bin_first = bin_number(range.first, eta_bins_trig_effcy) + 1; // + 1 pushes into next bin (bin lower end agree with range edge if range edge founded)
	            int bin_last = bin_number(range.second, eta_bins_trig_effcy); // bin higher end agree with range edge if range edge founded
	            std::string proj_suffix = pairToSuffix(range);

	            proj_make_and_write(h_pt2nd_vs_q_eta2nd_mu4_mu4noL1_sepr,     h_pt2nd_vs_q_eta2nd_mu4_sepr, true, bin_first, bin_last, proj_suffix);
	            proj_make_and_write(h_pt2nd_vs_q_eta2nd_2mu4_sepr,            h_pt2nd_vs_q_eta2nd_mu4_sepr, true, bin_first, bin_last, proj_suffix);
	        }
	    }
	}
}







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

void RDFBasedHistFillingPbPb::SumSingleMuonTrigEffHists(){
	std::cout << "Calling SumSingleMuonTrigEffHists" << std::endl;
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
        trg_effcy_filters_2D_3D_post_sum,
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
        trg_effcy_filters_2D_3D_post_sum,
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
void RDFBasedHistFillingPbPb::CalculateSingleMuonTrigEffcyRatios(){}

//--------- ---------
void RDFBasedHistFillingPbPb::MakeAndWriteDRTrigEffGraphs(){}
