#include "../MuonObjectsParamsAndHelpers/muon_pair_enums_MC.h"

#include "RDFBasedHistFillingPowheg.cxx"

// ---------------------------------------------------------------------------------------------
// pp24 MC-vs-DATA comparison, GENERIC family (2026-08-25)
// ---------------------------------------------------------------------------------------------
// Binding rule of the generic family: EVERY sample in it carries the SAME cuts as the data --
// all opposite-sign (resp. same-sign) pairs plus the fiducial detector-gap cut on BOTH muons,
// and NOTHING else: no mass window, no pair-pT threshold, no flavour/origin requirement. On the
// MC side the gap cut is taken on TRUTH q*eta.
//
// This block is STRICTLY ADDITIVE. It does NOT touch `FillHistogramsSignalAcceptance`, whose
// one-sided per-muon `q*eta < 2.2` is the POWHEG truth SIGNAL selection feeding the NLO
// template fit (docs/tracking/low_mass_dimuon_template_fit.md): changing that selection in
// place would silently bend the template-fit input. The gap cut below lives BESIDE it, on its
// own dataframes and its own histogram keys.
//
// Everything is declared file-locally rather than as class members because the class layout
// lives in RDFBasedHistFillingPowheg.h, which this change deliberately leaves untouched. Same
// pattern as RDFBasedHistFillingPythiaFullsim.cxx.
namespace {

// The six comparison observables. Their axes are requested BY NAME from `hist_binning_map`
// (see var1D_powheg_truth.json), so they are literally the SAME edge vectors the pp24 data and
// the Pythia fullsim are booked on -- a bin-for-bin data/MC ratio is only meaningful if the two
// sides share the edge vector, and the previous POWHEG file was on a 100 x [0, 5.75] dR axis
// against the data's 40 x [0, 5.75] (ratio 2.5, not even rebinnable).
// NOTE the naming: hist_name has NO `truth_` prefix here, because EVERY variable in
// var1D_powheg_truth.json is a truth variable and none of them carries one.
const std::vector<std::string>& McVsDataGenericVar1Ds(){
    static const std::vector<std::string> v = {
        "dr_ppbin",           // dr_bins_1d           (40 x [0, 5.75])
        "dr_zoomin_ppbin",    // dr_zoomin_bins_1d    (20 x [0, 0.8])
        "dphi_ppbin",         // dphi_bins_1d         (64 x [-pi, pi])
        "dphi_zoomin_ppbin",  // dphi_zoomin_bins_1d  (20 x [-0.8, 0.8])
        "deta_zoomin_ppbin",  // deta_zoomin_bins_1d  (20 x [-0.8, 0.8])
        "minv_zoomin_ppbin"   // minv_zoomin_bins_1d  (40 x [0, 3])
    };
    return v;
}

// The gap-cut FLAVOUR partition, identical in definition to the one
// FillHistogramsTruthFlavorBinned already applies (`:296-306`), so the gap-cut histograms carry
// the very same classification axis as the rest of the file. It IS a partition -- mutually
// exclusive and exhaustive -- so the four members sum to the inclusive gap-cut yield, which is
// what lets the plot side rebuild an inclusive POWHEG curve without the 2.88x over-count that
// summing the three OVERLAPPING axes (mechanism / origin / flavour) produced.
const std::vector<std::pair<std::string, std::string>>& McVsDataGenericFlavourSplit(){
    static const std::vector<std::pair<std::string, std::string>> v = {
        {"_flavor_binned_single_b",    "from_same_b"},
        {"_flavor_binned_both_from_b", "!from_same_b && both_from_b"},
        {"_flavor_binned_both_from_c", "!from_same_b && !both_from_b && both_from_c"},
        {"_flavor_binned_others",      "!from_same_b && !both_from_b && !both_from_c"}
    };
    return v;
}

// The two pair categories of the generic family, and the dataframe key each is built on.
const std::vector<std::string>& McVsDataGenericSigns(){
    static const std::vector<std::string> v = {"_ss", "_op"};
    return v;
}

const char* kGapCutTag = "_gapcut_truth";

// Every filter suffix the generic family books: the INCLUSIVE gap-cut category plus its four
// flavour members, for both signs.
const std::vector<std::string>& McVsDataGenericFilters(){
    static const std::vector<std::string> v = []{
        std::vector<std::string> out;
        for (const std::string& sign : McVsDataGenericSigns()){
            out.push_back(sign + kGapCutTag);
            for (const auto& [flav_suffix, flav_expr] : McVsDataGenericFlavourSplit())
                out.push_back(sign + kGapCutTag + flav_suffix);
        }
        return out;
    }();
    return v;
}

} // namespace

void RDFBasedHistFillingPowhegTruth::InitializePowhegExtra(){
	weight_specifier_to_column_map["_jacobian_corrected"] = "weight_norm_over_truth_dr";

	levels_truth_general_filters = {
		{"_near", "_away"},
		{"_ss", "_op"}
	};

	levels_truth_origin_filters = {
		{"_ss", "_op"},
		{
			"_origin_binned_gg",
			"_origin_binned_gq",
			"_origin_binned_single_gluon",
			"_origin_binned_qqbar",
			"_origin_binned_incoming",
			"_origin_binned_single_b",
			"_origin_binned_all"
		}
	};

	levels_truth_flavor_filters = {
		{"_ss", "_op"},
		{
			"_flavor_binned_single_b",
			"_flavor_binned_both_from_b",
			"_flavor_binned_both_from_c",
			"_flavor_binned_others"
		}
	};

	truth_general_var1Ds = {
		"dphi", "dr", "dr_zoomin",
		"pair_y", "pt_asym", "pair_pt_ptlead_ratio",
		"mQQ", "mQQ_Q_ratio", "mQQ_mHard_ratio"
	};
	truth_general_var2Ds = {
		{"dphi", "eta_avg"},
		{"dphi", "deta"},
		{"eta2", "eta1"},
		{"deta", "eta_avg"},
		{"pt2", "pt1"},
		{"pair_pt", "pt_lead"},
		{"pair_pt", "minv"},
		{"pair_pt", "minv_zoomin"},
		{"pair_pt_log", "minv_log"}
	};
	truth_general_var1Ds_jacobian = {"dr", "dr_zoomin"};
	truth_general_var2Ds_jacobian = {
		{"pair_pt", "minv"},
		{"pair_pt", "minv_zoomin"}
	};

	truth_origin_var1Ds = {
		"dr", "dr_zoomin", "pt_asym", "pair_pt_ptlead_ratio",
		"mQQ", "mQQ_Q_ratio", "mQQ_mHard_ratio"
	};
	truth_origin_var2Ds = {
		{"pair_pt", "pt_lead"},
		{"dphi", "deta"},
		{"pair_pt", "minv"},
		{"pair_pt", "minv_zoomin"}
	};
	truth_origin_var1Ds_jacobian = {"dr", "dr_zoomin"};
	truth_origin_var2Ds_jacobian = {
		{"pair_pt", "minv"},
		{"pair_pt", "minv_zoomin"}
	};

	truth_flavor_var1Ds = {
		"dr", "dr_zoomin", "pt_asym", "pair_pt_ptlead_ratio"
	};
	truth_flavor_var2Ds = {
		{"pair_pt", "pt_lead"},
		{"dphi", "deta"},
		{"pair_pt", "minv"},
		{"pair_pt", "minv_zoomin"}
	};
	truth_flavor_var1Ds_jacobian = {"dr", "dr_zoomin"};
	truth_flavor_var2Ds_jacobian = {
		{"pair_pt", "minv"},
		{"pair_pt", "minv_zoomin"}
	};
}

void RDFBasedHistFillingPowhegTruth::BuildSimpleFilterToVarListMapPowhegExtra(){
	for (const auto& filter : truth_general_filters){
		InsertOrAppend(df_filter_and_weight_to_var1D_list_map, std::make_pair(filter, ""), truth_general_var1Ds);
		InsertOrAppend(df_filter_and_weight_to_var2D_list_map, std::make_pair(filter, ""), truth_general_var2Ds);

		InsertOrAppend(df_filter_and_weight_to_var1D_list_map, std::make_pair(filter, "_jacobian_corrected"), truth_general_var1Ds_jacobian);
		InsertOrAppend(df_filter_and_weight_to_var2D_list_map, std::make_pair(filter, "_jacobian_corrected"), truth_general_var2Ds_jacobian);
	}

	for (const auto& filter : truth_origin_filters){
		InsertOrAppend(df_filter_and_weight_to_var1D_list_map, std::make_pair(filter, ""), truth_origin_var1Ds);
		InsertOrAppend(df_filter_and_weight_to_var2D_list_map, std::make_pair(filter, ""), truth_origin_var2Ds);

		InsertOrAppend(df_filter_and_weight_to_var1D_list_map, std::make_pair(filter, "_jacobian_corrected"), truth_origin_var1Ds_jacobian);
		InsertOrAppend(df_filter_and_weight_to_var2D_list_map, std::make_pair(filter, "_jacobian_corrected"), truth_origin_var2Ds_jacobian);
	}

	for (const auto& filter : truth_flavor_filters){
		InsertOrAppend(df_filter_and_weight_to_var1D_list_map, std::make_pair(filter, ""), truth_flavor_var1Ds);
		InsertOrAppend(df_filter_and_weight_to_var2D_list_map, std::make_pair(filter, ""), truth_flavor_var2Ds);

		InsertOrAppend(df_filter_and_weight_to_var1D_list_map, std::make_pair(filter, "_jacobian_corrected"), truth_flavor_var1Ds_jacobian);
		InsertOrAppend(df_filter_and_weight_to_var2D_list_map, std::make_pair(filter, "_jacobian_corrected"), truth_flavor_var2Ds_jacobian);
	}

	// --- pp24 MC-vs-data comparison, GENERIC family (see the anonymous namespace at the top) ---
	// APPENDED on filters of its own; no existing (filter, weight) entry is altered. Nominal
	// weight only ("" -> weight_norm): the data side has no 1/dR-weighted generic histogram, so
	// a jacobian variant here would have nothing to be compared against.
	for (const std::string& filter : McVsDataGenericFilters())
		InsertOrAppend(df_filter_and_weight_to_var1D_list_map,
		               std::make_pair(filter, std::string("")), McVsDataGenericVar1Ds());
}

void RDFBasedHistFillingPowhegTruth::FlattenFiltersExtra(){
	HistFillUtils::flatten_levels(levels_truth_general_filters, truth_general_filters);
	HistFillUtils::flatten_levels(levels_truth_origin_filters, truth_origin_filters);
	HistFillUtils::flatten_levels(levels_truth_flavor_filters, truth_flavor_filters);
}

void RDFBasedHistFillingPowhegTruth::BuildFlattenedFilterToVarListMapExtra(){
	BuildSimpleFilterToVarListMapPowhegExtra();
}

void RDFBasedHistFillingPowhegTruth::CreateBaseRDFsPowhegExtra(){
	auto define_truth_extras = [](ROOT::RDF::RNode node){
		return node
			.Define("weight_norm_over_truth_dr", "truth_dr > 0 ? weight_norm / truth_dr : 0.0")
			.Define("pair_pt_ptlead_ratio", "truth_pt_lead > 0 ? truth_pair_pt / truth_pt_lead : -1.0")
			.Define("mQQ_Q_ratio", "Q != 0 ? mQQ / Q : -10.0")
			.Define("mQQ_mHard_ratio", "mHard_relevant != 0 ? mQQ / mHard_relevant : -10.0");
	};

	ROOT::RDF::RNode& df_ss_weighted = map_at_checked(df_map, "df_ss_weighted", "CreateBaseRDFsPowhegTruthExtra: df_ss_weighted");
	ROOT::RDF::RNode& df_op_weighted = map_at_checked(df_map, "df_op_weighted", "CreateBaseRDFsPowhegTruthExtra: df_op_weighted");

	df_ss_weighted = define_truth_extras(df_ss_weighted);
	df_op_weighted = define_truth_extras(df_op_weighted);

	auto df_single_b_weighted = df_op_weighted.Filter("from_same_b && truth_dr < 1.0");
	df_map.emplace("df_single_b_weighted", df_single_b_weighted);

	// --- GENERIC family (2026-08-25): the MC mirror of the data generic selection ------------
	// ALL truth OS (resp. SS) pairs plus the fiducial gap cut on BOTH muons in TRUTH q*eta, and
	// NOTHING else. The gap windows are READ from ParamsSet::single_mu_fiducial_gap_cuts through
	// FiducialGapCutExpr and are never retyped -- the same call the pp24 data crossx and the
	// Pythia fullsim make, so the three sides cut on identical windows by construction.
	// Deliberately NOT `q*eta < 2.2`: that one-sided cut belongs to the signal selection at
	// FillHistogramsSignalAcceptance and is left exactly as it is.
	const std::string gap_truth = ParamsSet::FiducialGapCutExpr("m1.truth_charge * m1.truth_eta")
	                            + " && " + ParamsSet::FiducialGapCutExpr("m2.truth_charge * m2.truth_eta");

	for (const std::string& sign : McVsDataGenericSigns()){
		const std::string src_df = "df" + sign + "_weighted";
		ROOT::RDF::RNode& node = map_at_checked(df_map, src_df,
			("CreateBaseRDFsPowhegTruthExtra (generic family): " + src_df).c_str());
		df_map.emplace("df" + sign + kGapCutTag + "_weighted", node.Filter(gap_truth));
	}
}

void RDFBasedHistFillingPowhegTruth::FillHistogramsTruth(){
	FillHistogramsTruthGeneral();
	FillHistogramsTruthOriginBinned();
	FillHistogramsTruthFlavorBinned();
	FillHistogramsSignalAcceptance();

	// --- GENERIC family of the pp24 MC-vs-data comparison (2026-08-25) ----------------------
	// Booked on its own dataframes and its own histogram keys, so nothing above changes.
	// Inclusive gap-cut category FIRST, then the four members of the flavour partition on top
	// of the SAME gap-cut node -- so the four sum to the inclusive by construction and either
	// can be used on the plot side.
	// Written inline rather than as its own member function because the class layout lives in
	// RDFBasedHistFillingPowheg.h, which this change deliberately leaves untouched.
	try {
		for (const std::string& sign : McVsDataGenericSigns()){
			const std::string gap_filter = sign + kGapCutTag;
			const std::string gap_df = "df" + gap_filter + "_weighted";
			ROOT::RDF::RNode& node = map_at_checked(df_map, gap_df,
				("FillHistogramsTruth (generic family): " + gap_df).c_str());

			FillHistogramsSingleDataFrame(gap_filter, "", node);

			for (const auto& [flav_suffix, flav_expr] : McVsDataGenericFlavourSplit())
				FillHistogramsSingleDataFrame(gap_filter + flav_suffix, "", node.Filter(flav_expr));
		}
	}
	catch (const std::exception& e) {
		std::cerr << "FillHistogramsTruth (generic family): " << e.what() << std::endl;
		throw;
	}
}

void RDFBasedHistFillingPowhegTruth::FillHistogramsSignalAcceptance(){
	try {
		ROOT::RDF::RNode& df_op = map_at_checked(df_map, "df_op_weighted", "FillHistogramsSignalAcceptance (Powheg): df_op_weighted");

		const std::string signal_cuts =
			"from_same_b && truth_minv > 1.08 && truth_minv < 2.9 "
			"&& truth_pair_pt > 8 && m1.truth_charge * m1.truth_eta < 2.2 && m2.truth_charge * m2.truth_eta < 2.2";

		auto df_denom = df_op.Filter("from_same_b");
		auto df_num   = df_op.Filter(signal_cuts);

		const int     npt    = static_cast<int>(pms.pT_bins_120.size() - 1);
		const double* ptbins = pms.pT_bins_120.data();

		hist2d_rresultptr_map["h2d_sig_accept_num_pt_eta"] = df_num.Histo2D(
			ROOT::RDF::TH2DModel("h2d_sig_accept_num_pt_eta",
				";p_{T}^{pair} [GeV];#eta^{pair}", npt, ptbins, 44, -2.4, 2.4),
			"truth_pair_pt", "truth_pair_eta", "weight_norm");

		hist2d_rresultptr_map["h2d_sig_accept_denom_pt_eta"] = df_denom.Histo2D(
			ROOT::RDF::TH2DModel("h2d_sig_accept_denom_pt_eta",
				";p_{T}^{pair} [GeV];#eta^{pair}", npt, ptbins, 44, -2.4, 2.4),
			"truth_pair_pt", "truth_pair_eta", "weight_norm");

		// pT_bins_150 variants
		const int     npt150    = static_cast<int>(pms.pT_bins_150.size() - 1);
		const double* ptbins150 = pms.pT_bins_150.data();
		hist2d_rresultptr_map["h2d_sig_accept_num_pt_150_eta"] = df_num.Histo2D(
			ROOT::RDF::TH2DModel("h2d_sig_accept_num_pt_150_eta",
				";p_{T}^{pair} [GeV];#eta^{pair}", npt150, ptbins150, 44, -2.4, 2.4),
			"truth_pair_pt", "truth_pair_eta", "weight_norm");
		hist2d_rresultptr_map["h2d_sig_accept_denom_pt_150_eta"] = df_denom.Histo2D(
			ROOT::RDF::TH2DModel("h2d_sig_accept_denom_pt_150_eta",
				";p_{T}^{pair} [GeV];#eta^{pair}", npt150, ptbins150, 44, -2.4, 2.4),
			"truth_pair_pt", "truth_pair_eta", "weight_norm");
	}
	catch (const std::exception& e) {
		std::cerr << "FillHistogramsSignalAcceptance (Powheg): " << e.what() << std::endl;
		throw;
	}
}

void RDFBasedHistFillingPowhegTruth::HistPostProcessExtra(){
	// Signal acceptance ratio: num / denom (both in hist2D_map after base HistPostProcessBaseCommon)
	auto it_num = hist2D_map.find("h2d_sig_accept_num_pt_eta");
	auto it_den = hist2D_map.find("h2d_sig_accept_denom_pt_eta");
	if (it_num != hist2D_map.end() && it_den != hist2D_map.end()) {
		TH2D* hratio = static_cast<TH2D*>(it_num->second->Clone("h2d_sig_accept_pt_eta"));
		hratio->Divide(it_den->second);
		hist2D_map["h2d_sig_accept_pt_eta"] = hratio;
	} else {
		std::cerr << "[WARN] HistPostProcessExtra (Powheg): acceptance num/denom not found in hist2D_map." << std::endl;
	}

	auto it_num150 = hist2D_map.find("h2d_sig_accept_num_pt_150_eta");
	auto it_den150 = hist2D_map.find("h2d_sig_accept_denom_pt_150_eta");
	if (it_num150 != hist2D_map.end() && it_den150 != hist2D_map.end()) {
		TH2D* hratio150 = static_cast<TH2D*>(it_num150->second->Clone("h2d_sig_accept_pt_150_eta"));
		hratio150->Divide(it_den150->second);
		hist2D_map["h2d_sig_accept_pt_150_eta"] = hratio150;
	}
}

void RDFBasedHistFillingPowhegTruth::FillHistogramsTruthGeneral(){
	try {
		const std::string near_filter = "std::abs(truth_dphi) < 1.5707963267948966";
		const std::string away_filter = "std::abs(truth_dphi) >= 1.5707963267948966";

		ROOT::RDF::RNode& df_ss_weighted = map_at_checked(df_map, "df_ss_weighted", "FillHistogramsTruthGeneral: df_ss_weighted");
		ROOT::RDF::RNode& df_op_weighted = map_at_checked(df_map, "df_op_weighted", "FillHistogramsTruthGeneral: df_op_weighted");

		FillHistogramsSingleDataFrame("_near_ss", "", df_ss_weighted.Filter(near_filter));
		FillHistogramsSingleDataFrame("_away_ss", "", df_ss_weighted.Filter(away_filter));
		FillHistogramsSingleDataFrame("_near_op", "", df_op_weighted.Filter(near_filter));
		FillHistogramsSingleDataFrame("_away_op", "", df_op_weighted.Filter(away_filter));

		FillHistogramsSingleDataFrame("_near_ss", "_jacobian_corrected", df_ss_weighted.Filter(near_filter));
		FillHistogramsSingleDataFrame("_away_ss", "_jacobian_corrected", df_ss_weighted.Filter(away_filter));
		FillHistogramsSingleDataFrame("_near_op", "_jacobian_corrected", df_op_weighted.Filter(near_filter));
		FillHistogramsSingleDataFrame("_away_op", "_jacobian_corrected", df_op_weighted.Filter(away_filter));
	}
	catch (const std::exception& e) {
		std::cerr << "FillHistogramsTruthGeneral: " << e.what() << std::endl;
		throw;
	}
}

void RDFBasedHistFillingPowhegTruth::FillHistogramsTruthOriginBinned(){
	try {
		ROOT::RDF::RNode& df_ss_weighted = map_at_checked(df_map, "df_ss_weighted", "FillHistogramsTruthOriginBinned: df_ss_weighted");
		ROOT::RDF::RNode& df_op_weighted = map_at_checked(df_map, "df_op_weighted", "FillHistogramsTruthOriginBinned: df_op_weighted");

		const std::vector<std::pair<int, std::string>> powheg_ancestor_groups = {
			{powheg_ancestor_categories::gg, "gg"},
			{powheg_ancestor_categories::gq, "gq"},
			{powheg_ancestor_categories::single_gluon, "single_gluon"},
			{powheg_ancestor_categories::qqbar, "qqbar"},
			{powheg_ancestor_categories::incoming, "incoming"}
		};

		auto fill_one_sign = [this, &powheg_ancestor_groups](ROOT::RDF::RNode& df_weighted, const std::string& sign_suffix){
			FillHistogramsSingleDataFrame(sign_suffix + "_origin_binned_single_b", "", df_weighted.Filter("from_same_b"));
			FillHistogramsSingleDataFrame(sign_suffix + "_origin_binned_single_b", "_jacobian_corrected", df_weighted.Filter("from_same_b"));

			std::string matched_ancestor_expr = "(from_same_ancestors && (";
			for (std::size_t idx = 0; idx < powheg_ancestor_groups.size(); ++idx){
				if (idx > 0) matched_ancestor_expr += " || ";
				matched_ancestor_expr += "m1_ancestor_category == " + std::to_string(powheg_ancestor_groups[idx].first);
			}
			matched_ancestor_expr += "))";

			const std::string origin_fallback_expr = "!from_same_b && !(" + matched_ancestor_expr + ")";
			FillHistogramsSingleDataFrame(sign_suffix + "_origin_binned_all", "", df_weighted.Filter(origin_fallback_expr));
			FillHistogramsSingleDataFrame(sign_suffix + "_origin_binned_all", "_jacobian_corrected", df_weighted.Filter(origin_fallback_expr));

			for (const auto& [ancestor_catgr, ancestor_suffix] : powheg_ancestor_groups){
				const std::string filter_expr = "!from_same_b && from_same_ancestors && m1_ancestor_category == " + std::to_string(ancestor_catgr);
				const std::string filter_suffix = sign_suffix + "_origin_binned_" + ancestor_suffix;

				FillHistogramsSingleDataFrame(filter_suffix, "", df_weighted.Filter(filter_expr));
				FillHistogramsSingleDataFrame(filter_suffix, "_jacobian_corrected", df_weighted.Filter(filter_expr));
			}
		};

		fill_one_sign(df_ss_weighted, "_ss");
		fill_one_sign(df_op_weighted, "_op");
	}
	catch (const std::exception& e) {
		std::cerr << "FillHistogramsTruthOriginBinned: " << e.what() << std::endl;
		throw;
	}
}

void RunRDFBasedHistFillingPowhegTruth(){
	RDFBasedHistFillingPowhegTruth runner;
	runner.Run();
}

void RDFBasedHistFillingPowhegTruth::FillHistogramsTruthFlavorBinned(){
	try {
		ROOT::RDF::RNode& df_ss_weighted = map_at_checked(df_map, "df_ss_weighted", "FillHistogramsTruthFlavorBinned: df_ss_weighted");
		ROOT::RDF::RNode& df_op_weighted = map_at_checked(df_map, "df_op_weighted", "FillHistogramsTruthFlavorBinned: df_op_weighted");

		auto fill_one_sign = [this](ROOT::RDF::RNode& df_weighted, const std::string& sign_suffix){
			FillHistogramsSingleDataFrame(sign_suffix + "_flavor_binned_single_b", "", df_weighted.Filter("from_same_b"));
			FillHistogramsSingleDataFrame(sign_suffix + "_flavor_binned_single_b", "_jacobian_corrected", df_weighted.Filter("from_same_b"));

			FillHistogramsSingleDataFrame(sign_suffix + "_flavor_binned_both_from_b", "", df_weighted.Filter("!from_same_b && both_from_b"));
			FillHistogramsSingleDataFrame(sign_suffix + "_flavor_binned_both_from_b", "_jacobian_corrected", df_weighted.Filter("!from_same_b && both_from_b"));

			FillHistogramsSingleDataFrame(sign_suffix + "_flavor_binned_both_from_c", "", df_weighted.Filter("!from_same_b && !both_from_b && both_from_c"));
			FillHistogramsSingleDataFrame(sign_suffix + "_flavor_binned_both_from_c", "_jacobian_corrected", df_weighted.Filter("!from_same_b && !both_from_b && both_from_c"));

			FillHistogramsSingleDataFrame(sign_suffix + "_flavor_binned_others", "", df_weighted.Filter("!from_same_b && !both_from_b && !both_from_c"));
			FillHistogramsSingleDataFrame(sign_suffix + "_flavor_binned_others", "_jacobian_corrected", df_weighted.Filter("!from_same_b && !both_from_b && !both_from_c"));
		};

		fill_one_sign(df_ss_weighted, "_ss");
		fill_one_sign(df_op_weighted, "_op");
	}
	catch (const std::exception& e) {
		std::cerr << "FillHistogramsTruthFlavorBinned: " << e.what() << std::endl;
		throw;
	}
}
