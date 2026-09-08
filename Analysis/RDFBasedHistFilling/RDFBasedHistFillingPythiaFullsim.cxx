// CommonEffcyConfig.h MUST come FIRST: RDFBasedHistFillingPythia.h (pulled in by
// RDFBasedHistFillingPythia.cxx) default-initialises a data member from `CommonEffcyConfig{}`
// (RDFBasedHistFillingPythia.h:134) but does not include the header itself, so with the includes
// the other way round this translation unit does not compile ("use of undeclared identifier
// 'CommonEffcyConfig'"). Pre-existing since the 2026-08-18 header edit; found 2026-08-25.
#include "CommonEffcyConfig.h"
#include "RDFBasedHistFillingPythia.cxx"
#include "../Utilities/GeneralUtils.h"
#include "../MuonObjectsParamsAndHelpers/FullSimSampleType.h"
#include <cmath>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>

// ---------------------------------------------------------------------------------------------
// pp24 MC-vs-DATA comparison families (2026-08-25)
// ---------------------------------------------------------------------------------------------
// The comparison plots are built in TWO families, and within one family EVERY sample carries the
// SAME cuts -- the MC in TRUTH quantities:
//   SIGNAL  : the data signal region. On the MC side that is exactly the pre-existing
//             `_single_b_pass_signal_truth` filter (from_same_b + `pass_signal_truth`, which is
//             already bit-for-bit the data signal region: m_uu in (1.08, 2.9), pair pT > 8 GeV,
//             fiducial gap cut on BOTH muons). No new cut is introduced here.
//   GENERIC : all OS (resp. SS) pairs + the fiducial gap cut on both muons, and NOTHING else --
//             no mass window, no pair-pT threshold, no `from_same_b`. The MC mirror of the data
//             generic histograms, with the gap cut taken on TRUTH q*eta.
// The axes are requested BY NAME from `hist_binning_map` so they are the SAME vectors the data
// crossx histograms are booked on (.claude/CLAUDE.md §Binnings) -- a bin-for-bin ratio is only
// meaningful if the two sides literally share the edge vector.
namespace {

const std::vector<std::string>& McVsDataVar1Ds(){
    static const std::vector<std::string> v = {
        "truth_dr_zoomin_ppbin",    // dr_zoomin_bins_1d
        "truth_dphi_zoomin_ppbin",  // dphi_zoomin_bins_1d
        "truth_deta_zoomin_ppbin",  // deta_zoomin_bins_1d
        "truth_minv_zoomin_ppbin",  // minv_zoomin_bins_1d
        "truth_pair_eta_crossx",    // pair_eta_crossx
        "truth_pair_pt_log_150"     // pT_bins_150
    };
    return v;
}

// FULL-RANGE dR and dphi. GENERIC family ONLY: inside the signal region dR is kinematically
// bounded (m < 2.9 GeV with pair pT > 8 GeV forces dR ~< 2m/pT = 0.725), so a full-range view
// there carries nothing the zoom-in does not -- but the GENERIC family spans the away-side peak
// at dR ~ pi, which is exactly what the generic DR and Dphi panels exist to show, and without
// these two the Pythia curve was simply absent from those panels.
// Same axes as the data's `h_DR_<sign>` / `h_Dphi_<sign>`: dr_bins_1d, dphi_bins_1d.
const std::vector<std::string>& McVsDataGenericOnlyVar1Ds(){
    static const std::vector<std::string> v = {
        "truth_dr_ppbin",   // dr_bins_1d   (40 x [0, 5.75])
        "truth_dphi_ppbin"  // dphi_bins_1d (64 x [-pi, pi])
    };
    return v;
}

// The MC partner of the data's `h2d_crossx_pt_150_pair_eta_binned_w_signal_cuts`: pair pT on
// `pT_bins_150` x pair eta on `pair_eta_crossx`. SIGNAL family only.
const std::vector<std::array<std::string,2>>& McVsDataVar2Ds(){
    static const std::vector<std::array<std::string,2>> v = {
        {"truth_pair_pt_log_150", "truth_pair_eta_crossx"}
    };
    return v;
}

// SIGNAL-family filters: the pre-existing truth signal region, on all three pair categories.
//   _single_b_pass_signal_truth : the physics signal, OS pairs from the same b -- the OS partner
//                                 of the data.
//   _op_pass_signal_truth       : all truth OS pairs in the signal region ("MC including its own
//                                 background") reference.
//   _ss_pass_signal_truth       : all truth SS pairs in the signal region -- the partner of the
//                                 data SS panel. `from_same_b` has no same-sign counterpart by
//                                 construction (df_single_b_weighted = df_op_weighted.Filter(
//                                 "from_same_b")), so the SS partner is the plain SS category.
const std::vector<std::string>& McVsDataSignalFilters(){
    static const std::vector<std::string> v = {
        "_single_b_pass_signal_truth", "_op_pass_signal_truth", "_ss_pass_signal_truth"
    };
    return v;
}

// The 2D pair-pT x pair-eta view is needed for the single-b signal only.
const std::string& McVsDataSignal2DFilter(){
    static const std::string f = "_single_b_pass_signal_truth";
    return f;
}

// GENERIC-family filters (new; created in CreateBaseRDFsPythiaFullsimExtra).
const std::vector<std::string>& McVsDataGenericFilters(){
    static const std::vector<std::string> v = {"_op_gapcut_truth", "_ss_gapcut_truth"};
    return v;
}

} // namespace

void RDFBasedHistFillingPythiaFullsim::SetIOPathsHook(){
    // is_test_sample selects which fullsim production's NTuple-processing output to read, and
    // MUST match the isTestSample used to produce it (PythiaAlgCoreT.h / FullSimSampleType.h):
    //   true  -> pp24 TEST sample  (4 isospin beams; its cross-section carries the Pb 4:6:6:9
    //            average and is therefore NOT a physical pp cross-section -- label it honestly)
    //   false -> pp FULL sample    (pp beam only, isospin weight 1; suffix "_full")
    // TRUE today because only the TEST sample has NTP output; flip when the FULL sample lands.
    const std::string data_dir = FullSimSampleInputDir(FullSimSampleType::pp, is_test_sample);
    const std::string sample_suffix = is_test_sample ? "" : "_full";
    const std::string cut_suffix = with_data_resonance_cuts
        ? "_with_data_resonance_cuts"
        : "_no_data_resonance_cuts";

    input_files.clear();
    input_files.push_back(data_dir + "muon_pairs_pythia_fullsim_pp24" + cut_suffix + sample_suffix + ".root");
    output_file  = data_dir + "histograms_pythia_fullsim_pp24" + cut_suffix + sample_suffix + ".root";
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
        // DIAGNOSTIC view (fine axes; pre-existing).
        {"truth_pair_pt", "truth_pair_eta", "truth_dr_zoomin"},
        // The APPLIED eps_reco cells: canonical coarse pair pT x coarse pair eta x the
        // reco-efficiency dR edges. Consumed by
        // plotting_codes/reco_effcy/build_pp24_fullsim_pair_reco_eff.C.
        {"truth_pair_pt_coarse", "truth_pair_eta_coarse", "truth_dr_effcy"}
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

    // --- pp24 MC-vs-data comparison families (see the anonymous namespace at the top) ---------
    // APPENDED to whatever the reco-efficiency / detector-response booking already put on these
    // filters; no existing histogram set is altered.
    for (const std::string& filter : McVsDataSignalFilters())
        InsertOrAppend(df_filter_and_weight_to_var1D_list_map,
                       std::make_pair(filter, std::string("")), McVsDataVar1Ds());
    InsertOrAppend(df_filter_and_weight_to_var2D_list_map,
                   std::make_pair(McVsDataSignal2DFilter(), std::string("")), McVsDataVar2Ds());

    for (const std::string& filter : McVsDataGenericFilters()){
        InsertOrAppend(df_filter_and_weight_to_var1D_list_map,
                       std::make_pair(filter, std::string("")), McVsDataVar1Ds());
        InsertOrAppend(df_filter_and_weight_to_var1D_list_map,
                       std::make_pair(filter, std::string("")), McVsDataGenericOnlyVar1Ds());
    }
}

void RDFBasedHistFillingPythiaFullsim::BuildHistBinningMapPythiaFullsimExtra(){
    hist_binning_map["eta_bins_reco_effcy"] = ParamsSet::makeEtaTrigEffcyBinning(1);

    // --- axes of the 3D PAIR reco efficiency the cross-section applies (2026-08-17) -------------
    // eps_reco(pair pT, pair eta, dR) is a per-pair WEIGHT, so its cells must be the CANONICAL
    // coarse binnings and nothing else (.claude/CLAUDE.md §Binnings): the same pair-pT vector the
    // MC trigger dR correction is celled on, and the same pair-eta ranges every pair-eta view of
    // this analysis uses. Read from ParamsSet / CommonEffcyConfig, never retyped.
    // The fine `pT_bins_80` / 48-bin eta axes of the pre-existing reco-eff VIEWS are untouched --
    // those are diagnostic projections, not the applied correction.
    hist_binning_map["pair_pt_coarse_bins"] = pms.pair_pt_coarse_bins;

    {
        static const CommonEffcyConfig eff_cfg{};
        std::vector<double> eta_edges;
        eta_edges.push_back(eff_cfg.pair_eta_proj_ranges_coarse_incl_gap.front().first);
        for (const auto& r : eff_cfg.pair_eta_proj_ranges_coarse_incl_gap)
            eta_edges.push_back(r.second);
        hist_binning_map["pair_eta_coarse_bins"] = eta_edges;
    }

    // dR axis = the projection edges the reco-efficiency views have used since the pipeline was
    // written (`dr_bins_edges_for_reco_effcy`), promoted to a histogram axis. It is NOT extended
    // beyond 1.0 and does not need to be: inside the signal region dR is KINEMATICALLY bounded,
    // dR ~< 2 m_uu / pT^pair <= 2 * 2.9 / 8 = 0.725, and the measured maximum over the pp24
    // fullsim single-b truth signal region is 0.701. Nothing falls above the top edge.
    {
        std::vector<double> dr_edges(dr_bins_edges_for_reco_effcy.begin(),
                                     dr_bins_edges_for_reco_effcy.end());
        hist_binning_map["dr_bins_reco_effcy"] = dr_edges;
    }
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

    // SIGNAL REGION, truth and reco legs. Since 2026-08-17 the per-muon one-sided
    // `q*eta < 2.2` is REPLACED by the detector-gap FIDUCIAL cut required of BOTH muons,
    // in LOCKSTEP with the data crossx (RDFBasedHistFillingPP.cxx `signal_cuts`).
    // Windows are READ from ParamsSet::single_mu_fiducial_gap_cuts and never retyped.
    // The cut is applied on TRUTH q*eta in the denominator leg and on RECO q*eta in the
    // numerator leg, so eps_reco is a FIDUCIAL pair efficiency: it does NOT contain the
    // truth-level gap acceptance eps_acc (muon_gap_cuts_acceptance.md F12), which stays a
    // separate, not-yet-built factor.
    // Hoisted out of the category loop (2026-08-25) so the GENERIC family below reuses the
    // very SAME `gap_truth` string as the signal region -- one expression, one definition.
    // The PAIR-LEVEL window |eta^pair| < ParamsSet::pair_eta_fiducial_max = 2.2 (user,
    // 2026-09-07) is part of the same gap definition and travels with the single-muon windows:
    // truth pair eta on the denominator leg, reco pair eta on the numerator leg, so eps_reco
    // stays the efficiency of exactly the region the data crossx selects.
    const std::string gap_truth = ParamsSet::FiducialGapCutExpr("m1.truth_charge * m1.truth_eta")
                                + " && " + ParamsSet::FiducialGapCutExpr("m2.truth_charge * m2.truth_eta")
                                + " && " + ParamsSet::PairFiducialEtaCutExpr("truth_pair_eta");
    const std::string gap_reco  = ParamsSet::FiducialGapCutExpr("m1.charge * m1.eta")
                                + " && " + ParamsSet::FiducialGapCutExpr("m2.charge * m2.eta")
                                + " && " + ParamsSet::PairFiducialEtaCutExpr("pair_eta");

    for (const std::string& pair_catgr : {"_ss", "_op", "_single_b"}){
        const std::string df_name = "df" + pair_catgr;
        ROOT::RDF::RNode& node = map_at_checked(df_map, df_name + "_weighted",
            Form("CreateBaseRDFsPythiaFullsimExtra: df_map.at(%s)", (df_name + "_weighted").c_str()));

        auto node_sig = node
            .Define("pass_signal_truth",
                "truth_minv > 1.08 && truth_minv < 2.9 && truth_pair_pt > 8 && " + gap_truth)
            .Define("pass_signal_reco",
                "(m1.reco_match && m2.reco_match) ? (minv > 1.08 && minv < 2.9 && pair_pt > 8 && "
                + gap_reco + ") : false");

        df_map.emplace(df_name + "_pass_medium_weighted",   node_sig.Filter("pair_pass_medium"));
        df_map.emplace(df_name + "_pass_tight_weighted",    node_sig.Filter("pair_pass_tight"));
        df_map.emplace(df_name + "_pass_signal_truth_weighted", node_sig.Filter("pass_signal_truth"));
        df_map.emplace(df_name + "_pass_medium_and_signal_truth_and_reco_weighted",
            node_sig.Filter("pair_pass_medium && pass_signal_truth && pass_signal_reco"));
        df_map.emplace(df_name + "_pass_tight_and_signal_truth_and_reco_weighted",
            node_sig.Filter("pair_pass_tight && pass_signal_truth && pass_signal_reco"));
    }

    // --- GENERIC family (2026-08-25) ---------------------------------------------------------
    // The MC mirror of the data generic histograms: ALL truth OS (resp. SS) pairs plus the
    // fiducial gap cut on BOTH muons in TRUTH q*eta, and NOTHING else -- deliberately no mass
    // window, no pair-pT threshold and no `from_same_b`, because the data generic selection has
    // none of those either. Same `gap_truth` expression as the signal region above.
    for (const std::string& pair_catgr : {"_ss", "_op"}){
        const std::string df_name = "df" + pair_catgr + "_weighted";
        ROOT::RDF::RNode& node = map_at_checked(df_map, df_name,
            Form("CreateBaseRDFsPythiaFullsimExtra: df_map.at(%s)", df_name.c_str()));
        df_map.emplace("df" + pair_catgr + "_gapcut_truth_weighted", node.Filter(gap_truth));
    }
}

void RDFBasedHistFillingPythiaFullsim::FillHistogramsFullSim(){
    FillHistogramsFullSimDetecResp();
    FillHistogramsFullSimRecoEffcies();

    // GENERIC family of the pp24 MC-vs-data comparison. The SIGNAL family needs no separate call:
    // its filters `_{single_b,op,ss}_pass_signal_truth` are already booked by
    // FillHistogramsFullSimRecoEffcies, and BuildFilterToVarListMapExtra APPENDED the comparison
    // variables to those filters' var lists.
    // DELIBERATELY NOT wrapped in the try/catch the sibling methods below use. If one of these
    // dataframes or variables is missing, the MC curve simply VANISHES from the generic panels
    // and the plotter prints a `[SKIP]` that nobody reads -- a silently incomplete output, the
    // failure mode `.claude/kb/.../reference_root_swallows_rdf_exceptions` exists to prevent.
    // Let it throw: `map_at_checked` names the missing key.
    for (const std::string& filter : McVsDataGenericFilters()){
        const std::string df_name = "df" + filter + "_weighted";
        ROOT::RDF::RNode& node = map_at_checked(df_map, df_name,
            Form("FillHistogramsFullSim (generic family): df_map.at(%s)", df_name.c_str()));
        FillHistogramsSingleDataFrame(filter, "", node);
    }
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
