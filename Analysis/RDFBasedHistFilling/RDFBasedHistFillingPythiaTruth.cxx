#include "RDFBasedHistFillingPythia.cxx"
#include "../MuonObjectsParamsAndHelpers/muon_pair_enums_MC_utils.h"
#include <algorithm>
#include <set>
#include <stdexcept>
#include <sstream>
#include <type_traits>

namespace {
bool replace_first(std::string& s, const std::string& from, const std::string& to){
    const auto pos = s.find(from);
    if (pos == std::string::npos) return false;
    s.replace(pos, from.size(), to);
    return true;
}
}

void RDFBasedHistFillingPythiaTruth::InitializePythiaExtra(){
    vars1D_general = {
       "DR", "DR_zoomin", "Deta_zoomin", "Dphi_zoomin", "Dphi", "minv", "minv_zoomin", "pt_lead", "pair_pt",
        "pair_y", "pt_asym", "psrapidity_ordered_pt_asym", "pair_pt_ptlead_ratio"
    };

    vars1D_general_over_dr = {
        "DR_jacobian_corrected", "DR_zoomin_jacobian_corrected"
    };

    vars1D_general_over_pair_pt = {
        "pair_pt_jacobian_corrected"
    };

    vars2D_general = {
        {"Dphi", "eta_avg"},
        {"Dphi", "Deta"},
        {"eta2", "eta1"},
        {"Deta", "eta_avg"},
        {"pt2", "pt1"},
        {"pair_pt", "pt_lead"},
        {"pair_pt", "pt_lead_zoomin"},
        {"pair_pt_log", "pt_lead"},
        {"pair_pt", "minv"},
        {"pair_pt", "minv_zoomin"},
        {"pair_pt_log", "minv_log"}
    };

    vars2D_general_over_dr = {
        {"pair_pt", "minv_jacobian_corrected"},
        {"pair_pt", "minv_zoomin_jacobian_corrected"}
    };

    vars1D_flavor_origin = {
        "DR", "DR_zoomin", "Deta_zoomin", "Dphi_zoomin", "pt_asym", "psrapidity_ordered_pt_asym", "pair_pt_ptlead_ratio"
    };

    vars1D_flavor_origin_over_dr = {
        "DR_jacobian_corrected", "DR_zoomin_jacobian_corrected"
    };

    vars1D_flavor_origin_over_pair_pt = {
        "pair_pt_jacobian_corrected"
    };

    vars2D_flavor_origin = {
        {"Dphi", "Deta"},
        {"pair_pt", "minv"},
        {"pair_pt", "minv_10GeV"},
        {"pair_pt", "minv_zoomin"},
        {"pair_pt_log", "minv_log"},
        {"pair_pt", "pt_lead"},
        {"pair_pt", "pt_lead_zoomin"},
        {"pair_pt_log", "pt_lead"}
    };

    vars2D_flavor_origin_over_dr = {
        {"pair_pt", "minv_jacobian_corrected"},
        {"pair_pt", "minv_zoomin_jacobian_corrected"}
    };

    flavor_grp_map_build(flavor_suffix_map);
    origin_grp_map_build(origin_suffix_map);
}

void RDFBasedHistFillingPythiaTruth::BuildHistBinningMapExtra(){
    std::vector<double> pair_pt_bins;
    pair_pt_bins.reserve(61);
    for (int i = 0; i <= 60; ++i) pair_pt_bins.push_back(static_cast<double>(i) * 0.5);

    hist_binning_map["pair_pt_bins_linear"] = pair_pt_bins;
}

void RDFBasedHistFillingPythiaTruth::CreateBaseRDFsPythiaExtra(){
    auto define_truth_extras = [](ROOT::RDF::RNode node){
        return node
            .Define("truth_pair_dPoverP", "truth_pair_pt != 0 ? truth_dpt / truth_pair_pt : 0.0")
            .Define("truth_pair_pt_ptlead_ratio", "truth_pt_lead > 0 ? truth_pair_pt / truth_pt_lead : -1.0")
            .Define("weight_over_dr", "truth_dr > 0 ? weight / truth_dr : 0.0")
            .Define("weight_over_pair_pt", "truth_pair_pt > 0 ? weight / truth_pair_pt : 0.0")
            .Define("Qsplit_pTHat_ratio", "pTHat != 0 ? Qsplit / pTHat : -10.0")
            .Define("Qsplit_mHat_ratio", "mHard_relevant != 0 ? Qsplit / mHard_relevant : -10.0")
            .Define("truth_psrapidity_ordered_pt_asym",
                    "(std::abs(m1.truth_eta) > std::abs(m2.truth_eta)) ? "
                    "((m1.truth_pt - m2.truth_pt) / (m1.truth_pt + m2.truth_pt)) : "
                    "((m2.truth_pt - m1.truth_pt) / (m2.truth_pt + m1.truth_pt))");
    };

    ROOT::RDF::RNode& df_ss_weighted = map_at_checked(df_map, "df_ss_weighted", "CreateBaseRDFsPythiaExtra: df_ss_weighted");
    ROOT::RDF::RNode& df_op_weighted = map_at_checked(df_map, "df_op_weighted", "CreateBaseRDFsPythiaExtra: df_op_weighted");

    df_ss_weighted = define_truth_extras(df_ss_weighted);
    df_op_weighted = define_truth_extras(df_op_weighted);

    ValidatePythiaTruthSchemaAndCoverage();
}

void RDFBasedHistFillingPythiaTruth::ValidatePythiaTruthSchemaAndCoverage(){
    Long64_t total_ss_entries = 0;
    Long64_t total_op_entries = 0;

    for (const auto& path : input_files){
        TFile f(path.c_str(), "READ");
        if (f.IsZombie()) {
            throw std::runtime_error("ValidatePythiaTruthSchemaAndCoverage: failed to open input file " + path);
        }

        TTree* tss = dynamic_cast<TTree*>(f.Get(tree_ss.c_str()));
        TTree* top = dynamic_cast<TTree*>(f.Get(tree_op.c_str()));

        total_ss_entries += (tss ? tss->GetEntries() : 0);
        total_op_entries += (top ? top->GetEntries() : 0);
    }

    if (total_ss_entries == 0 && total_op_entries == 0){
        std::ostringstream oss;
        oss << "ValidatePythiaTruthSchemaAndCoverage: input trees are empty for all files. "
            << "RDF cannot produce plotting histograms from empty trees. "
            << "sign1 entries=" << total_ss_entries << ", sign2 entries=" << total_op_entries << ".";
        throw std::runtime_error(oss.str());
    }

    auto& df_ss_weighted = map_at_checked(df_map, "df_ss_weighted", "ValidatePythiaTruthSchemaAndCoverage");
    const auto cols_vec = df_ss_weighted.GetColumnNames();
    std::set<std::string> cols(cols_vec.begin(), cols_vec.end());

    auto require_var_defined = [this](const std::string& vname){
        if (var1D_dict.find(vname) == var1D_dict.end()) {
            throw std::runtime_error("ValidatePythiaTruthSchemaAndCoverage: missing var1D json definition for '" + vname + "'.");
        }
    };

    auto require_column_exists = [this, &cols](const std::string& vname){
        const auto it = var1D_dict.find(vname);
        if (it == var1D_dict.end()) {
            throw std::runtime_error("ValidatePythiaTruthSchemaAndCoverage: missing var1D json definition for '" + vname + "'.");
        }
        const std::string& col = it->second->var;
        if (cols.find(col) == cols.end()) {
            throw std::runtime_error(
                "ValidatePythiaTruthSchemaAndCoverage: JSON var_name '" + col
                + "' (hist_name='" + vname + "') is not an RDF column. "
                + "Check MuonPairPythiaTruth output branch names and Define(...) columns.");
        }
    };

    const std::vector<std::string> required_plot_1d = {
        "DR", "DR_zoomin", "Deta_zoomin", "Dphi_zoomin",
        "DR_jacobian_corrected", "DR_zoomin_jacobian_corrected",
        "pt_asym", "psrapidity_ordered_pt_asym", "pair_pt_ptlead_ratio"
    };

    const std::vector<std::array<std::string, 2>> required_plot_2d = {
        {"pair_pt", "pt_lead"},
        {"pair_pt", "minv"},
        {"pair_pt", "minv_zoomin"},
        {"pair_pt", "minv_jacobian_corrected"},
        {"pair_pt", "minv_zoomin_jacobian_corrected"},
        {"Dphi", "Deta"}
    };

    std::set<std::string> configured_1d;
    configured_1d.insert(vars1D_flavor_origin.begin(), vars1D_flavor_origin.end());
    configured_1d.insert(vars1D_flavor_origin_over_dr.begin(), vars1D_flavor_origin_over_dr.end());
    configured_1d.insert(vars1D_flavor_origin_over_pair_pt.begin(), vars1D_flavor_origin_over_pair_pt.end());

    std::set<std::string> configured_2d_keys;
    auto pack_pair = [](const std::array<std::string, 2>& p){ return p[0] + "|" + p[1]; };
    for (const auto& p : vars2D_flavor_origin) configured_2d_keys.insert(pack_pair(p));
    for (const auto& p : vars2D_flavor_origin_over_dr) configured_2d_keys.insert(pack_pair(p));

    for (const auto& v : required_plot_1d){
        require_var_defined(v);
        require_column_exists(v);
        if (configured_1d.find(v) == configured_1d.end()) {
            throw std::runtime_error(
                "ValidatePythiaTruthSchemaAndCoverage: RDF flavor/origin fills do not include required plotting 1D variable '" + v + "'.");
        }
    }

    for (const auto& p : required_plot_2d){
        require_var_defined(p[0]);
        require_var_defined(p[1]);
        require_column_exists(p[0]);
        require_column_exists(p[1]);
        if (configured_2d_keys.find(pack_pair(p)) == configured_2d_keys.end()) {
            throw std::runtime_error(
                "ValidatePythiaTruthSchemaAndCoverage: RDF flavor/origin fills do not include required plotting 2D variable pair ('"
                + p[0] + "', '" + p[1] + "').");
        }
    }
}

void RDFBasedHistFillingPythiaTruth::FillHistogramsTruth(){
    FillHistogramsGeneral();
    FillHistogramsFlavorBinned();
    FillHistogramsOriginBinned();
    FillHistogramsResonanceStudy();
    FillHistogramsCrossxAndSpecialEta();
    FillHistogramsSignalAcceptance();
}

void RDFBasedHistFillingPythiaTruth::FillHistogramsGeneral(){
    try {
        const std::vector<std::pair<std::string, std::string>> sign_df_map = {
            {"_sign1", "df_ss_weighted"},
            {"_sign2", "df_op_weighted"}
        };

        const std::string near_filter = "std::abs(truth_dphi) < 1.5707963267948966";
        const std::string away_filter = "std::abs(truth_dphi) >= 1.5707963267948966";

        for (const auto& [sign_suffix, df_name] : sign_df_map){
            ROOT::RDF::RNode& node = map_at_checked(df_map, df_name, "FillHistogramsGeneral");

            auto fill_one_region = [this, &sign_suffix](auto region_df, const std::string& region){
                const std::string suffix = "_" + region + sign_suffix;

                FillHistogramsSingleDataFrame(suffix, region_df, "weight", vars1D_general, vars2D_general, {});
                FillHistogramsSingleDataFrame(suffix, region_df, "weight_over_dr", vars1D_general_over_dr, vars2D_general_over_dr, {});
                FillHistogramsSingleDataFrame(suffix, region_df, "weight_over_pair_pt", vars1D_general_over_pair_pt, {}, {});
            };

            auto node_near = node.Filter(near_filter);
            auto node_away = node.Filter(away_filter);

            fill_one_region(node_near, "near");
            fill_one_region(node_away, "away");
        }
    }
    catch (const std::exception& e) {
        std::cerr << "FillHistogramsGeneral: " << e.what() << std::endl;
        throw;
    }
}

void RDFBasedHistFillingPythiaTruth::FillHistogramsFlavorBinned(){
    try {
        const std::vector<std::pair<std::string, std::string>> sign_df_map = {
            {"_sign1", "df_ss_weighted"},
            {"_sign2", "df_op_weighted"}
        };

        for (const auto& [sign_suffix, df_name] : sign_df_map){
            ROOT::RDF::RNode& node = map_at_checked(df_map, df_name, "FillHistogramsFlavorBinned");

            for (const auto& [iflavor, flavor_suffix] : flavor_suffix_map){
                const std::string filter_expr = "muon_pair_flavor_category == " + std::to_string(iflavor);
                const std::string suffix = sign_suffix + flavor_suffix;
                auto node_flavor = node.Filter(filter_expr);

                FillHistogramsSingleDataFrame(suffix, node_flavor, "weight", vars1D_flavor_origin, vars2D_flavor_origin, {});
                FillHistogramsSingleDataFrame(suffix, node_flavor, "weight_over_dr", vars1D_flavor_origin_over_dr, vars2D_flavor_origin_over_dr, {});
                FillHistogramsSingleDataFrame(suffix, node_flavor, "weight_over_pair_pt", vars1D_flavor_origin_over_pair_pt, {}, {});
            }
        }
    }
    catch (const std::exception& e) {
        std::cerr << "FillHistogramsFlavorBinned: " << e.what() << std::endl;
        throw;
    }
}

void RDFBasedHistFillingPythiaTruth::FillHistogramsOriginBinned(){
    try {
        const std::vector<std::pair<std::string, std::string>> sign_df_map = {
            {"_sign1", "df_ss_weighted"},
            {"_sign2", "df_op_weighted"}
        };

        for (const auto& [sign_suffix, df_name] : sign_df_map){
            ROOT::RDF::RNode& node = map_at_checked(df_map, df_name, "FillHistogramsOriginBinned");

            for (const auto& [iorigin, origin_suffix] : origin_suffix_map){
                const std::string filter_expr = "muon_pair_origin_category == " + std::to_string(iorigin);
                const std::string suffix = sign_suffix + origin_suffix;
                auto node_origin = node.Filter(filter_expr);

                FillHistogramsSingleDataFrame(suffix, node_origin, "weight", vars1D_flavor_origin, vars2D_flavor_origin, {});
                FillHistogramsSingleDataFrame(suffix, node_origin, "weight_over_dr", vars1D_flavor_origin_over_dr, vars2D_flavor_origin_over_dr, {});
                FillHistogramsSingleDataFrame(suffix, node_origin, "weight_over_pair_pt", vars1D_flavor_origin_over_pair_pt, {}, {});

                FillHistogramsSingleDataFrame(suffix, node_origin, "weight",
                                              {"Qsplit", "Qsplit_pTHat_ratio", "Qsplit_mHat_ratio"}, {}, {});
            }
        }
    }
    catch (const std::exception& e) {
        std::cerr << "FillHistogramsOriginBinned: " << e.what() << std::endl;
        throw;
    }
}

void RDFBasedHistFillingPythiaTruth::FillHistogramsResonanceStudy(){
    try {
        const std::vector<std::pair<std::string, std::string>> sign_df_map = {
            {"_sign1", "df_ss_weighted"},
            {"_sign2", "df_op_weighted"}
        };

        const std::vector<std::string> signal_vars = {"minv_sub_GeV_signal", "minv_single_b_region_signal"};
        const std::vector<std::string> bkg_vars = {"minv_sub_GeV_resonance_and_res_contam_bkg", "minv_single_b_region_resonance_and_res_contam_bkg"};

        for (const auto& [sign_suffix, df_name] : sign_df_map){
            ROOT::RDF::RNode& node = map_at_checked(df_map, df_name, "FillHistogramsResonanceStudy");

            auto signal_node = node.Filter("from_same_b");
            auto bkg_node = node.Filter("from_same_resonance || resonance_contaminated");

            FillHistogramsSingleDataFrame("_no_res_cut" + sign_suffix, signal_node, "weight", signal_vars, {}, {});
            FillHistogramsSingleDataFrame("_no_res_cut" + sign_suffix, signal_node, "weight_over_dr",
                                          {"minv_sub_GeV_jacobian_corrected_signal", "minv_single_b_region_jacobian_corrected_signal"}, {}, {});

            FillHistogramsSingleDataFrame("_old_res_cut" + sign_suffix,
                                          signal_node.Filter("!Reco_resonance_or_reso_contam_tagged_old"),
                                          "weight", signal_vars, {}, {});
            FillHistogramsSingleDataFrame("_old_res_cut" + sign_suffix,
                                          signal_node.Filter("!Reco_resonance_or_reso_contam_tagged_old"),
                                          "weight_over_dr",
                                          {"minv_sub_GeV_jacobian_corrected_signal", "minv_single_b_region_jacobian_corrected_signal"}, {}, {});

            FillHistogramsSingleDataFrame("_new_res_cut" + sign_suffix,
                                          signal_node.Filter("!Reco_resonance_or_reso_contam_tagged_new"),
                                          "weight", signal_vars, {}, {});
            FillHistogramsSingleDataFrame("_new_res_cut" + sign_suffix,
                                          signal_node.Filter("!Reco_resonance_or_reso_contam_tagged_new"),
                                          "weight_over_dr",
                                          {"minv_sub_GeV_jacobian_corrected_signal", "minv_single_b_region_jacobian_corrected_signal"}, {}, {});

            FillHistogramsSingleDataFrame("_no_res_cut" + sign_suffix, bkg_node, "weight", bkg_vars, {}, {});
            FillHistogramsSingleDataFrame("_no_res_cut" + sign_suffix, bkg_node, "weight_over_dr",
                                          {"minv_sub_GeV_jacobian_corrected_resonance_and_res_contam_bkg",
                                           "minv_single_b_region_jacobian_corrected_resonance_and_res_contam_bkg"}, {}, {});

            FillHistogramsSingleDataFrame("_old_res_cut" + sign_suffix,
                                          bkg_node.Filter("!Reco_resonance_or_reso_contam_tagged_old"),
                                          "weight", bkg_vars, {}, {});
            FillHistogramsSingleDataFrame("_old_res_cut" + sign_suffix,
                                          bkg_node.Filter("!Reco_resonance_or_reso_contam_tagged_old"),
                                          "weight_over_dr",
                                          {"minv_sub_GeV_jacobian_corrected_resonance_and_res_contam_bkg",
                                           "minv_single_b_region_jacobian_corrected_resonance_and_res_contam_bkg"}, {}, {});

            FillHistogramsSingleDataFrame("_new_res_cut" + sign_suffix,
                                          bkg_node.Filter("!Reco_resonance_or_reso_contam_tagged_new"),
                                          "weight", bkg_vars, {}, {});
            FillHistogramsSingleDataFrame("_new_res_cut" + sign_suffix,
                                          bkg_node.Filter("!Reco_resonance_or_reso_contam_tagged_new"),
                                          "weight_over_dr",
                                          {"minv_sub_GeV_jacobian_corrected_resonance_and_res_contam_bkg",
                                           "minv_single_b_region_jacobian_corrected_resonance_and_res_contam_bkg"}, {}, {});
        }
    }
    catch (const std::exception& e) {
        std::cerr << "FillHistogramsResonanceStudy: " << e.what() << std::endl;
        throw;
    }
}

void RDFBasedHistFillingPythiaTruth::FillHistogramsCrossxAndSpecialEta(){
    try {
        const std::vector<std::pair<std::string, std::string>> sign_df_map = {
            {"_sign1", "df_ss_weighted"},
            {"_sign2", "df_op_weighted"}
        };

        for (const auto& [sign_suffix, df_name] : sign_df_map){
            ROOT::RDF::RNode& node = map_at_checked(df_map, df_name, "FillHistogramsCrossxAndSpecialEta");

            auto node_single_b = node.Filter("from_same_b");
            FillHistogramsSingleDataFrame("_truth_from_single_b" + sign_suffix, node_single_b, "weight",
                                          {}, {{"pair_eta", "pair_pt_crossx"}}, {});

            FillHistogramsSingleDataFrame("_DPHI1" + sign_suffix,
                                          node.Filter("std::abs(truth_dphi) < 1.0"),
                                          "weight", {}, {{"eta2", "eta1"}}, {});

            FillHistogramsSingleDataFrame("_DPHI2" + sign_suffix,
                                          node.Filter("std::abs(truth_dphi) > (3.141592653589793 - 1.0)"),
                                          "weight", {}, {{"eta2", "eta1"}}, {});
        }
    }
    catch (const std::exception& e) {
        std::cerr << "FillHistogramsCrossxAndSpecialEta: " << e.what() << std::endl;
        throw;
    }
}

void RDFBasedHistFillingPythiaTruth::FillHistogramsSignalAcceptance(){
    try {
        ROOT::RDF::RNode& df_op = map_at_checked(df_map, "df_op_weighted", "FillHistogramsSignalAcceptance: df_op_weighted");

        const std::string signal_cuts =
            "from_same_b && truth_minv > 1.08 && truth_minv < 2.9 "
            "&& truth_pair_pt > 8 && m1.truth_charge * m1.truth_eta < 2.2 && m2.truth_charge * m2.truth_eta < 2.2 && truth_dr > 0.05";

        auto df_denom = df_op.Filter("from_same_b");
        auto df_num   = df_op.Filter(signal_cuts);

        const int     npt    = static_cast<int>(pms.pT_bins_120.size() - 1);
        const double* ptbins = pms.pT_bins_120.data();

        hist2d_rresultptr_map["h2d_sig_accept_num_pt_eta"] = df_num.Histo2D(
            ROOT::RDF::TH2DModel("h2d_sig_accept_num_pt_eta",
                ";p_{T}^{pair} [GeV];#eta^{pair}", npt, ptbins, 44, -2.4, 2.4),
            "truth_pair_pt", "truth_pair_eta", "weight");

        hist2d_rresultptr_map["h2d_sig_accept_denom_pt_eta"] = df_denom.Histo2D(
            ROOT::RDF::TH2DModel("h2d_sig_accept_denom_pt_eta",
                ";p_{T}^{pair} [GeV];#eta^{pair}", npt, ptbins, 44, -2.4, 2.4),
            "truth_pair_pt", "truth_pair_eta", "weight");

        // pT_bins_150 variants
        const int     npt150    = static_cast<int>(pms.pT_bins_150.size() - 1);
        const double* ptbins150 = pms.pT_bins_150.data();
        hist2d_rresultptr_map["h2d_sig_accept_num_pt_150_eta"] = df_num.Histo2D(
            ROOT::RDF::TH2DModel("h2d_sig_accept_num_pt_150_eta",
                ";p_{T}^{pair} [GeV];#eta^{pair}", npt150, ptbins150, 44, -2.4, 2.4),
            "truth_pair_pt", "truth_pair_eta", "weight");
        hist2d_rresultptr_map["h2d_sig_accept_denom_pt_150_eta"] = df_denom.Histo2D(
            ROOT::RDF::TH2DModel("h2d_sig_accept_denom_pt_150_eta",
                ";p_{T}^{pair} [GeV];#eta^{pair}", npt150, ptbins150, 44, -2.4, 2.4),
            "truth_pair_pt", "truth_pair_eta", "weight");
    }
    catch (const std::exception& e) {
        std::cerr << "FillHistogramsSignalAcceptance (Pythia): " << e.what() << std::endl;
        throw;
    }
}

void RDFBasedHistFillingPythiaTruth::HistPostProcessExtra(){
    BuildAndStoreNearAwaySummedHistograms();

    auto sum_sign_hists = [this](const std::string& base_name, const std::string& out_name){
        const std::string h1_name = base_name + "_sign1";
        const std::string h2_name = base_name + "_sign2";

        auto it1 = hist1D_map.find(h1_name);
        auto it2 = hist1D_map.find(h2_name);
        if (it1 != hist1D_map.end() && it2 != hist1D_map.end()){
            auto* hsum = static_cast<TH1D*>(it1->second->Clone(out_name.c_str()));
            hsum->Add(it2->second);
            hist1D_map[out_name] = hsum;
        }

        auto it21 = hist2D_map.find(h1_name);
        auto it22 = hist2D_map.find(h2_name);
        if (it21 != hist2D_map.end() && it22 != hist2D_map.end()){
            auto* hsum2 = static_cast<TH2D*>(it21->second->Clone(out_name.c_str()));
            hsum2->Add(it22->second);
            hist2D_map[out_name] = hsum2;
        }
    };

    sum_sign_hists("h_pair_pt_crossx_vs_pair_eta_truth_from_single_b", "h_crossx_truth_from_single_b_vs_pair_pt_pair_eta");

    const std::vector<std::string> combined_resonance_hists = {
        "h_minv_sub_GeV_signal_no_res_cut",
        "h_minv_sub_GeV_signal_old_res_cut",
        "h_minv_sub_GeV_signal_new_res_cut",
        "h_minv_sub_GeV_jacobian_corrected_signal_no_res_cut",
        "h_minv_sub_GeV_jacobian_corrected_signal_old_res_cut",
        "h_minv_sub_GeV_jacobian_corrected_signal_new_res_cut",
        "h_minv_sub_GeV_resonance_and_res_contam_bkg_no_res_cut",
        "h_minv_sub_GeV_resonance_and_res_contam_bkg_old_res_cut",
        "h_minv_sub_GeV_resonance_and_res_contam_bkg_new_res_cut",
        "h_minv_sub_GeV_jacobian_corrected_resonance_and_res_contam_bkg_no_res_cut",
        "h_minv_sub_GeV_jacobian_corrected_resonance_and_res_contam_bkg_old_res_cut",
        "h_minv_sub_GeV_jacobian_corrected_resonance_and_res_contam_bkg_new_res_cut",
        "h_minv_single_b_region_signal_no_res_cut",
        "h_minv_single_b_region_signal_old_res_cut",
        "h_minv_single_b_region_signal_new_res_cut",
        "h_minv_single_b_region_jacobian_corrected_signal_no_res_cut",
        "h_minv_single_b_region_jacobian_corrected_signal_old_res_cut",
        "h_minv_single_b_region_jacobian_corrected_signal_new_res_cut",
        "h_minv_single_b_region_resonance_and_res_contam_bkg_no_res_cut",
        "h_minv_single_b_region_resonance_and_res_contam_bkg_old_res_cut",
        "h_minv_single_b_region_resonance_and_res_contam_bkg_new_res_cut",
        "h_minv_single_b_region_jacobian_corrected_resonance_and_res_contam_bkg_no_res_cut",
        "h_minv_single_b_region_jacobian_corrected_resonance_and_res_contam_bkg_old_res_cut",
        "h_minv_single_b_region_jacobian_corrected_resonance_and_res_contam_bkg_new_res_cut"
    };

    for (const auto& hname : combined_resonance_hists){
        sum_sign_hists(hname, hname);
    }

    // Signal acceptance ratio: num / denom (both already in hist2D_map from base HistPostProcessBaseCommon)
    {
        auto it_num = hist2D_map.find("h2d_sig_accept_num_pt_eta");
        auto it_den = hist2D_map.find("h2d_sig_accept_denom_pt_eta");
        if (it_num != hist2D_map.end() && it_den != hist2D_map.end()) {
            TH2D* hratio = static_cast<TH2D*>(it_num->second->Clone("h2d_sig_accept_pt_eta"));
            hratio->Divide(it_den->second);
            hist2D_map["h2d_sig_accept_pt_eta"] = hratio;
        } else {
            std::cerr << "[WARN] HistPostProcessExtra (Pythia): acceptance num/denom not found in hist2D_map." << std::endl;
        }
        auto it_num150 = hist2D_map.find("h2d_sig_accept_num_pt_150_eta");
        auto it_den150 = hist2D_map.find("h2d_sig_accept_denom_pt_150_eta");
        if (it_num150 != hist2D_map.end() && it_den150 != hist2D_map.end()) {
            TH2D* hratio150 = static_cast<TH2D*>(it_num150->second->Clone("h2d_sig_accept_pt_150_eta"));
            hratio150->Divide(it_den150->second);
            hist2D_map["h2d_sig_accept_pt_150_eta"] = hratio150;
        }
    }

    MarkNearAwayHistogramsForNominalExclusion();
}

bool RDFBasedHistFillingPythiaTruth::IsNearAwayDividedHistName(const std::string& hname) const{
    return hname.find("_near_sign1") != std::string::npos
        || hname.find("_near_sign2") != std::string::npos
        || hname.find("_away_sign1") != std::string::npos
        || hname.find("_away_sign2") != std::string::npos;
}

std::string RDFBasedHistFillingPythiaTruth::NearAwayToSummedName(const std::string& hname) const{
    std::string out = hname;

    if (replace_first(out, "_near_sign1", "_sign1")) return out;
    if (replace_first(out, "_near_sign2", "_sign2")) return out;
    if (replace_first(out, "_away_sign1", "_sign1")) return out;
    if (replace_first(out, "_away_sign2", "_sign2")) return out;

    return "";
}

void RDFBasedHistFillingPythiaTruth::BuildAndStoreNearAwaySummedHistograms(){
    auto build_for_map = [this](auto& hmap){
        using HistPtr = typename std::decay_t<decltype(hmap)>::mapped_type;
        using HistT = std::remove_pointer_t<HistPtr>;

        std::vector<std::pair<std::string, HistT*>> to_add;

        for (const auto& [hname, hnear] : hmap){
            if (hname.find("_near_sign1") == std::string::npos
             && hname.find("_near_sign2") == std::string::npos) {
                continue;
            }

            std::string away_name = hname;
            if (!replace_first(away_name, "_near_sign1", "_away_sign1")
             && !replace_first(away_name, "_near_sign2", "_away_sign2")) {
                continue;
            }

            auto it_away = hmap.find(away_name);
            if (it_away == hmap.end()) continue;

            const std::string summed_name = NearAwayToSummedName(hname);
            if (summed_name.empty() || hmap.find(summed_name) != hmap.end()) continue;

            HistT* hsum = static_cast<HistT*>(hnear->Clone(summed_name.c_str()));
            hsum->Add(it_away->second);
            to_add.emplace_back(summed_name, hsum);
        }

        for (auto& kv : to_add){
            hmap.emplace(kv.first, kv.second);
        }
    };

    build_for_map(hist1D_map);
    build_for_map(hist2D_map);
    build_for_map(hist3D_map);
}

void RDFBasedHistFillingPythiaTruth::MarkNearAwayHistogramsForNominalExclusion(){
    auto mark_map = [this](const auto& hmap){
        for (const auto& [hname, _] : hmap){
            if (!IsNearAwayDividedHistName(hname)) continue;
            if (std::find(hists_to_not_write.begin(), hists_to_not_write.end(), hname) == hists_to_not_write.end()) {
                hists_to_not_write.push_back(hname);
            }
        }
    };

    mark_map(hist1D_map);
    mark_map(hist2D_map);
    mark_map(hist3D_map);
}

void RDFBasedHistFillingPythiaTruth::WriteOutputExtra(){
    std::string divided_output_file = output_file;
    const std::string suffix = "_near_away_divided";

    const auto pos = divided_output_file.rfind(".root");
    if (pos == std::string::npos) {
        divided_output_file += suffix;
    } else {
        divided_output_file.insert(pos, suffix);
    }

    TFile divided_file(divided_output_file.c_str(), "recreate");
    if (divided_file.IsZombie()) {
        throw std::runtime_error("RDFBasedHistFillingPythiaTruth::WriteOutputExtra: failed to create output file " + divided_output_file);
    }

    divided_file.cd();

    for (const auto& [hname, h] : hist1D_map){
        if (IsNearAwayDividedHistName(hname)) h->Write();
    }
    for (const auto& [hname, h] : hist2D_map){
        if (IsNearAwayDividedHistName(hname)) h->Write();
    }
    for (const auto& [hname, h] : hist3D_map){
        if (IsNearAwayDividedHistName(hname)) h->Write();
    }
}

void RunRDFBasedHistFillingPythiaTruth(bool isPrivate = false,
                                       double ecom = 5.36,
                                       bool with_data_resonance_cuts = false)
{
    RDFBasedHistFillingPythiaTruth runner(isPrivate, ecom, with_data_resonance_cuts);
    runner.Run();
}
