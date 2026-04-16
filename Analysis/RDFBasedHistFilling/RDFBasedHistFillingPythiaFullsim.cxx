#include "RDFBasedHistFillingPythia.cxx"
#include "../Utilities/GeneralUtils.h"
#include <cmath>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>

void RDFBasedHistFillingPythiaFullsim::SetIOPathsHook(){
    const std::string data_dir = "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_test_sample/";
    const std::string cut_suffix = with_data_resonance_cuts
        ? "_with_data_resonance_cuts"
        : "_no_data_resonance_cuts";

    input_files.clear();
    input_files.push_back(data_dir + "muon_pairs_pythia_fullsim_pp24" + cut_suffix + ".root");
    output_file  = data_dir + "histograms_pythia_fullsim_pp24" + cut_suffix + ".root";
    infile_var1D_json = "var1D_pythia_fullsim.json";
}

void RDFBasedHistFillingPythiaFullsim::InitializePythiaFullsimExtra(){
    levels_reco_effcy_filters = {
        {"_ss", "_op", "_single_b"},
        {"", "_pass_medium", "_pass_tight",
         "_pass_signal_truth",
         "_pass_medium_and_signal_truth_and_reco",
         "_pass_tight_and_signal_truth_and_reco"}
    };

    levels_detector_response_filters = {
        {"_single_b"},
        {"_pass_medium", "_pass_tight"}
    };

    auto make_ranges_from_edges = [](const std::vector<float>& edges){
        std::vector<std::pair<float,float>> ranges;
        if (edges.size() < 2) return ranges;
        ranges.reserve(edges.size() - 1);
        for (std::size_t i = 0; i + 1 < edges.size(); ++i)
            ranges.emplace_back(edges[i], edges[i+1]);
        return ranges;
    };

    dr_ranges_for_reco_effcy        = make_ranges_from_edges(dr_bins_edges_for_reco_effcy);
    pair_pT_ranges_for_reco_effcy_dR = make_ranges_from_edges(pair_pT_bins_edges_for_reco_effcy_dR);

    mu_pair_reco_eff_proj_cfgs = {
        {"truth_pair_pt",  "truth_dr_zoomin", false, &dr_ranges_for_reco_effcy},
        {"truth_pair_eta", "truth_dr_zoomin", false, &dr_ranges_for_reco_effcy},
        {"truth_pair_pt",  "truth_dr_zoomin", true,  &pair_pT_ranges_for_reco_effcy_dR}
    };

    reco_eff_num_denom_suffix_pairs = {
        {"_pass_medium", ""},
        {"_pass_tight",  ""}
    };
}

void RDFBasedHistFillingPythiaFullsim::BuildFilterToVarListMapExtra(){
    HistFillUtils::flatten_levels(levels_reco_effcy_filters,        reco_effcy_filters);
    HistFillUtils::flatten_levels(levels_detector_response_filters, detector_response_filters);

    reco_effcy_var1Ds = {
        "truth_pair_pt", "truth_pair_eta",
        "truth_dr_zoomin", "truth_dr_2_0",
        "truth_minv_zoomin", "truth_deta_zoomin", "truth_dphi_zoomin"
    };
    reco_effcy_var2Ds = {
        {"truth_pair_pt",   "truth_pair_eta"},
        {"truth_pair_pt",   "truth_dr_zoomin"},
        {"truth_pair_eta",  "truth_dr_zoomin"},
        {"truth_deta_zoomin", "truth_dphi_zoomin"},
        {"truth_pair_pt",   "truth_minv_zoomin"},
        {"truth_minv_zoomin", "truth_dr_zoomin"}
    };
    reco_effcy_var3Ds = {
        {"truth_pair_pt", "truth_pair_eta", "truth_dr_zoomin"}
    };

    detec_resp_var1Ds = {
        "truth_pair_pt", "pair_pt",
        "truth_minv_zoomin", "minv_zoomin",
        "truth_dr_zoomin", "dr_zoomin"
    };
    detec_resp_var2Ds = {
        {"truth_pair_pt",   "pair_pt"},
        {"truth_minv_zoomin", "minv_zoomin"},
        {"truth_dr_zoomin", "dr_zoomin"}
    };

    for (auto& filter : reco_effcy_filters)
        InsertOrAppend(df_filter_and_weight_to_var1D_list_map, std::make_pair(filter, ""), reco_effcy_var1Ds);
    for (auto& filter : reco_effcy_filters)
        InsertOrAppend(df_filter_and_weight_to_var2D_list_map, std::make_pair(filter, ""), reco_effcy_var2Ds);
    for (auto& filter : reco_effcy_filters)
        InsertOrAppend(df_filter_and_weight_to_var3D_list_map, std::make_pair(filter, ""), reco_effcy_var3Ds);

    for (auto& filter : detector_response_filters)
        InsertOrAppend(df_filter_and_weight_to_var1D_list_map, std::make_pair(filter, ""), detec_resp_var1Ds);
    for (auto& filter : detector_response_filters)
        InsertOrAppend(df_filter_and_weight_to_var2D_list_map, std::make_pair(filter, ""), detec_resp_var2Ds);
}

void RDFBasedHistFillingPythiaFullsim::BuildHistBinningMapPythiaFullsimExtra(){
    hist_binning_map["eta_bins_reco_effcy"] = ParamsSet::makeEtaTrigEffcyBinning(1);
}

void RDFBasedHistFillingPythiaFullsim::CreateBaseRDFsPythiaFullsimExtra(){
    {
        var1D* v = map_at_checked(var1D_dict, std::string("truth_dr_2_0"),
            "CreateBaseRDFsPythiaFullsimExtra: var1D_dict.at(truth_dr_2_0)");
        if (v->var != "truth_dr")
            throw std::runtime_error(
                "Invalid var1D mapping: truth_dr_2_0 must map to column 'truth_dr', got '"
                + v->var + "'. Check var1D_pythia_fullsim.json.");
    }

    ROOT::RDF::RNode& df_op_weighted = map_at_checked(df_map, "df_op_weighted",
        "CreateBaseRDFsPythiaFullsimExtra: df_op_weighted");
    df_map.emplace("df_single_b_weighted", df_op_weighted.Filter("from_same_b"));

    for (const std::string& pair_catgr : {"_ss", "_op", "_single_b"}){
        const std::string df_name = "df" + pair_catgr;
        ROOT::RDF::RNode& node = map_at_checked(df_map, df_name + "_weighted",
            Form("CreateBaseRDFsPythiaFullsimExtra: df_map.at(%s)", (df_name + "_weighted").c_str()));

        auto node_sig = node
            .Define("pass_signal_truth",
                "truth_minv > 1.08 && truth_minv < 2.9 && truth_pair_pt > 8 && truth_pair_eta < 2.2 && truth_dr > 0.05")
            .Define("pass_signal_reco",
                "(m1.reco_match && m2.reco_match) ? (minv > 1.08 && minv < 2.9 && pair_pt > 8 && pair_eta < 2.2 && dr > 0.05) : false");

        df_map.emplace(df_name + "_pass_medium_weighted",   node_sig.Filter("pair_pass_medium"));
        df_map.emplace(df_name + "_pass_tight_weighted",    node_sig.Filter("pair_pass_tight"));
        df_map.emplace(df_name + "_pass_signal_truth_weighted", node_sig.Filter("pass_signal_truth"));
        df_map.emplace(df_name + "_pass_medium_and_signal_truth_and_reco_weighted",
            node_sig.Filter("pair_pass_medium && pass_signal_truth && pass_signal_reco"));
        df_map.emplace(df_name + "_pass_tight_and_signal_truth_and_reco_weighted",
            node_sig.Filter("pair_pass_tight && pass_signal_truth && pass_signal_reco"));
    }
}

void RDFBasedHistFillingPythiaFullsim::FillHistogramsFullSim(){
    FillHistogramsFullSimDetecResp();
    FillHistogramsFullSimRecoEffcies();
}

void RDFBasedHistFillingPythiaFullsim::FillHistogramsFullSimDetecResp(){
    try {
        for (const char* pair_catgr : {"_single_b"}){
            for (const char* quality_catgr : {"_pass_medium", "_pass_tight"}){
                const std::string filter  = std::string(pair_catgr) + quality_catgr;
                const std::string df_name = "df" + filter + "_weighted";
                ROOT::RDF::RNode& node = map_at_checked(df_map, df_name,
                    Form("FillHistogramsFullSimDetecResp: df_map.at(%s)", df_name.c_str()));
                FillHistogramsSingleDataFrame(filter, "", node);
            }
        }
    } catch (const std::out_of_range& e){
        std::cerr << "FillHistogramsFullSimDetecResp: out_of_range: " << e.what() << "\n";
    } catch (const std::runtime_error& e){
        std::cerr << "FillHistogramsFullSimDetecResp: runtime_error: " << e.what() << "\n";
    }
}

void RDFBasedHistFillingPythiaFullsim::FillHistogramsFullSimRecoEffcies(){
    try {
        const std::vector<std::string> quality_cats = {
            "", "_pass_medium", "_pass_tight",
            "_pass_signal_truth",
            "_pass_medium_and_signal_truth_and_reco",
            "_pass_tight_and_signal_truth_and_reco"
        };

        for (const std::string& pair_catgr : {"_ss", "_op", "_single_b"}){
            for (const std::string& quality_catgr : quality_cats){
                const std::string filter  = pair_catgr + quality_catgr;
                const std::string df_name = "df" + filter + "_weighted";
                ROOT::RDF::RNode& node = map_at_checked(df_map, df_name,
                    Form("FillHistogramsFullSimRecoEffcies: df_map.at(%s)", df_name.c_str()));
                FillHistogramsSingleDataFrame(filter, "", node);
            }
        }
    } catch (const std::out_of_range& e){
        std::cerr << "FillHistogramsFullSimRecoEffcies: out_of_range: " << e.what() << "\n";
    } catch (const std::runtime_error& e){
        std::cerr << "FillHistogramsFullSimRecoEffcies: runtime_error: " << e.what() << "\n";
    }
}

void RDFBasedHistFillingPythiaFullsim::MakeAndWriteMuPairRecoEffProjGraphsHelper(
    const std::vector<std::string>& categories,
    bool use_TH_divide,
    bool require_signal_cuts)
{
    const std::vector<std::string> default_cats = {"_ss", "_op", "_single_b"};
    const std::vector<std::string>& cats = categories.empty() ? default_cats : categories;

    std::set<std::string> skipped_empty;
    std::set<std::string> xcheck_mismatch;

    auto axis_edges = [](const TAxis* axis){
        std::vector<double> edges;
        if (!axis) return edges;
        const int nbins = axis->GetNbins();
        if (nbins <= 0) return edges;
        if (axis->GetXbins() && axis->GetXbins()->GetSize() > 0){
            const double* arr = axis->GetXbins()->GetArray();
            edges.assign(arr, arr + nbins + 1);
            return edges;
        }
        edges.resize(nbins + 1);
        const double xmin = axis->GetXmin(), xmax = axis->GetXmax();
        for (int i = 0; i <= nbins; ++i)
            edges[i] = xmin + (xmax - xmin) * (static_cast<double>(i) / nbins);
        return edges;
    };

    auto range_to_suffix = [](const std::pair<float,float>& range){
        const bool upper_is_max = !std::isfinite(range.second)
            || range.second >= std::numeric_limits<float>::max() * 0.5f;
        if (upper_is_max){
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(2) << range.first;
            std::string low = oss.str();
            for (auto& c : low) if (c == '.') c = '_';
            while (low.find('-') != std::string::npos) low.replace(low.find('-'), 1, "minus");
            return low + "_TO_MAX";
        }
        return pairToSuffix(range);
    };

    try {
        for (const auto& cfg : mu_pair_reco_eff_proj_cfgs){
            const std::string& varx    = std::get<0>(cfg);
            const std::string& vary    = std::get<1>(cfg);
            const bool project_y       = std::get<2>(cfg);
            const auto* proj_ranges    = std::get<3>(cfg);
            if (!proj_ranges || proj_ranges->empty()) continue;

            for (const std::string& catgr : cats){
                for (const auto& num_denom : reco_eff_num_denom_suffix_pairs){
                    const std::string num_suffix   = require_signal_cuts
                        ? (num_denom.first + "_and_signal_truth_and_reco") : num_denom.first;
                    const std::string denom_suffix = require_signal_cuts
                        ? "_pass_signal_truth" : num_denom.second;

                    const std::string hname_num = "h_" + vary + "_vs_" + varx + catgr + num_suffix;
                    const std::string hname_den = "h_" + vary + "_vs_" + varx + catgr + denom_suffix;

                    TH2D* h_num = map_at_checked(hist2D_map, hname_num,
                        Form("MakeAndWriteMuPairRecoEffProjGraphsHelper: hist2D_map.at(%s)", hname_num.c_str()));
                    TH2D* h_den = map_at_checked(hist2D_map, hname_den,
                        Form("MakeAndWriteMuPairRecoEffProjGraphsHelper: hist2D_map.at(%s)", hname_den.c_str()));

                    const TAxis* axis_for_ranges = project_y ? h_num->GetXaxis() : h_num->GetYaxis();
                    const std::vector<double> edges = axis_edges(axis_for_ranges);
                    if (edges.size() < 2) continue;

                    const int nbins_axis    = axis_for_ranges->GetNbins();
                    const double axis_max   = edges.back();

                    for (const auto& range : *proj_ranges){
                        int bin_first = bin_number(range.first, edges) + 1;
                        const bool upper_is_max = !std::isfinite(range.second)
                            || range.second >= std::numeric_limits<float>::max() * 0.5f
                            || range.second >= axis_max;
                        int bin_last = upper_is_max ? nbins_axis : bin_number(range.second, edges);

                        if (bin_first < 1) bin_first = 1;
                        if (bin_first > nbins_axis) bin_first = nbins_axis;
                        if (bin_last  < 1) bin_last  = 1;
                        if (bin_last  > nbins_axis) bin_last  = nbins_axis;
                        if (bin_last < bin_first) continue;

                        const std::string proj_suffix     = range_to_suffix(range);
                        const std::string proj_axis_suffix = project_y ? "_py" : "_px";
                        const std::string proj_full_suffix = proj_suffix.empty()
                            ? proj_axis_suffix : (proj_axis_suffix + "_" + proj_suffix);

                        std::unique_ptr<TH1D> h_num_proj(
                            project_y
                            ? h_num->ProjectionY(Form("%s%s", h_num->GetName(), proj_full_suffix.c_str()), bin_first, bin_last, "e")
                            : h_num->ProjectionX(Form("%s%s", h_num->GetName(), proj_full_suffix.c_str()), bin_first, bin_last, "e")
                        );
                        std::unique_ptr<TH1D> h_den_proj(
                            project_y
                            ? h_den->ProjectionY(Form("%s%s", h_den->GetName(), proj_full_suffix.c_str()), bin_first, bin_last, "e")
                            : h_den->ProjectionX(Form("%s%s", h_den->GetName(), proj_full_suffix.c_str()), bin_first, bin_last, "e")
                        );
                        if (!h_num_proj || !h_den_proj) continue;

                        const double den_integral = h_den_proj->Integral(1, h_den_proj->GetNbinsX());
                        const bool denom_empty = (h_den_proj->GetEntries() <= 0.0 || den_integral <= 0.0);

                        if (use_TH_divide){
                            TH1D* h_divided = dynamic_cast<TH1D*>(h_num_proj->Clone(
                                (std::string(h_num_proj->GetName()) + "_divided").c_str()));
                            if (!h_divided) continue;
                            h_divided->SetDirectory(nullptr);
                            h_divided->SetStats(0);
                            if (denom_empty){
                                skipped_empty.insert(h_den_proj->GetName());
                                h_divided->Reset("ICES");
                                std::cout << "[RecoEffProjZeroFill] "
                                    << "varx=" << varx << ", vary=" << vary
                                    << ", catgr=" << catgr
                                    << ", num=" << num_suffix << ", den=" << denom_suffix
                                    << ", axis=" << (project_y ? "Y" : "X")
                                    << ", range=" << proj_suffix
                                    << " -> denominator empty, zero-filled\n";
                            } else {
                                h_divided->Divide(h_den_proj.get());
                            }
                            mu_pair_reco_eff_proj_hist_map[h_divided->GetName()] = h_divided;
                        } else {
                            if (denom_empty){
                                skipped_empty.insert(h_den_proj->GetName());
                                continue;
                            }
                            TGraphAsymmErrors* g = HistFillUtils::divide_and_write(
                                h_num_proj.get(), h_den_proj.get(), &mu_pair_reco_eff_proj_graph_map);
                            if (g){
                                g->GetXaxis()->SetTitle(h_num_proj->GetXaxis()->GetTitle());
                            }
                        }
                    }
                }
            }
        }

        std::cout << "[RecoEffProjSummary] skipped-empty-denom=" << skipped_empty.size()
                  << ", xcheck-mismatch=" << xcheck_mismatch.size() << "\n";
    } catch (const std::out_of_range& e){
        std::cerr << "MakeAndWriteMuPairRecoEffProjGraphsHelper: out_of_range: " << e.what() << "\n";
    } catch (const std::runtime_error& e){
        std::cerr << "MakeAndWriteMuPairRecoEffProjGraphsHelper: runtime_error: " << e.what() << "\n";
    }
}

void RDFBasedHistFillingPythiaFullsim::MakeAndWriteMuPairRecoEffProjGraphs(){
    MakeAndWriteMuPairRecoEffProjGraphsHelper({"_ss", "_op", "_single_b"}, true);
    MakeAndWriteMuPairRecoEffProjGraphsHelper({"_ss", "_op", "_single_b"}, true, true);
}

void RDFBasedHistFillingPythiaFullsim::WriteOutputExtra(){
    HistFillUtils::write_hist_map_vector(mu_pair_reco_eff_proj_graph_map, mu_pair_reco_eff_proj_graphs_to_not_write);
    HistFillUtils::write_hist_map_vector(mu_pair_reco_eff_proj_hist_map,  mu_pair_reco_eff_proj_hists_to_not_write);
}

void RDFBasedHistFillingPythiaFullsim::CleanupExtra(){
    for (auto& kv : mu_pair_reco_eff_proj_graph_map) delete kv.second;
    for (auto& kv : mu_pair_reco_eff_proj_hist_map)  delete kv.second;
}
