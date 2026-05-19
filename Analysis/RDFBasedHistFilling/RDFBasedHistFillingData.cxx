#include "RDFBasedHistFillingData.h"
#include <TF1.h>
#include <TH2D.h>
#include <type_traits>

static TFile* s_effcy_pT_fit_file = nullptr;
static std::map<std::string, TF1*> s_effcy_pT_fit_map;
static TFile* s_effcy_2D_hist_file = nullptr;
static std::map<std::string, TH2D*> s_effcy_2D_hist_map;

RDFBasedHistFillingData::RDFBasedHistFillingData(int run_year_input, bool isForSoumya_input)
: run_year (run_year_input), isForSoumya (isForSoumya_input){
    std::cout << " Histogram filling for data:" << std::endl << std::endl;

    std::cout << "The following public variable(s) **MUST** be checked:" << std::endl;
    std::cout << "--> trigger_mode: [INT] default 1 (single_mu4)" << std::endl;
    std::cout << "--> hist_filling_cycle: [INT] default generic" << std::endl;
    std::cout << "--> mu4_nominal_pbpb_NO_trig_calc: [BOOL] default false (set true for PbPb nominal pipeline)" << std::endl << std::endl;
    std::cout << "--> isForSoumya: [BOOL] default false" << std::endl;
    std::cout << "--> useCoarseQEtaBin: [BOOL] default true" << std::endl << std::endl;

    std::cout << "The following public variable(s) **SHOULD** be checked:" << std::endl;
    std::cout << "--> isScram: [BOOL] default false" << std::endl;
    std::cout << "--> isTight: [BOOL] default false" << std::endl;
    std::cout << "--> output_generic_hists: [BOOL] default false" << std::endl;
    std::cout << "--> save_non_sepr_trg_hists: [BOOL] default false" << std::endl;
    std::cout << "--> save_good_accept_trg_hists: [BOOL] default false" << std::endl;
    std::cout << std::endl;
}

void RDFBasedHistFillingData::InitializeDataImpl(){
    InitializeDataCommon();
    InitializeDataExtra();
}

void RDFBasedHistFillingData::InitializeDataCommon(){
    run_year %= 2000;

    isForSoumya_suffix = isForSoumya? "_for_soumya" : "";

    useCoarseQEtaBin |= isForSoumya;

    qEtaBin_suffix = isForSoumya    ? ""
                                    : (useCoarseQEtaBin  ? "_coarse_q_eta_bin"
                                                        : "_fine_q_eta_bin");

    SetQEtaProjRanges(run_year, q_eta_proj_ranges, q_eta_ranges_str, useCoarseQEtaBin);

    if (isForSoumya){
        single_muon_trig_effcy_var1Ds = {};
        single_muon_trig_effcy_var2Ds = {{"q_eta2nd","pt2nd"}};
        single_muon_trig_effcy_var3Ds = {};
    }

    // Compute input_mindR_suffix: > 0 means use new files with mindR suffix; <= 0 means old files (no suffix)
    if (mindR_trig > 0) {
        std::string mindR_str = (mindR_trig == 0.01) ? "0_01" : "0_02";
        input_mindR_suffix = "_mindR_" + mindR_str;
    } else {
        input_mindR_suffix = "";
    }
    std::cout << "RDFBasedHistFillingData: mindR_trig=" << mindR_trig
              << ", input_mindR_suffix='" << input_mindR_suffix << "'"
              << ", useMu4NoL1Leg=" << useMu4NoL1Leg << std::endl;

    trig_to_filter_str_map["_mu4"] = "";
    if (useMu4NoL1Leg) {
        trig_to_filter_str_map["_mu4_mu4noL1"]           = "passmu4mu4noL1 && mu2nd_passmu4noL1";
        trig_to_filter_str_map["_2mu4_AND_mu4_mu4noL1"]  = "pass2mu4 && passmu4mu4noL1 && mu2nd_passmu4noL1";
    } else {
        trig_to_filter_str_map["_mu4_mu4noL1"]           = "passmu4mu4noL1";
        trig_to_filter_str_map["_2mu4_AND_mu4_mu4noL1"]  = "pass2mu4 && passmu4mu4noL1";
    }
    trig_to_filter_str_map["_2mu4"] = "pass2mu4";

    trg_effcy_biases = {"_sepr"};
    if (save_non_sepr_trg_hists) trg_effcy_biases.push_back("");
    if (save_good_accept_trg_hists) trg_effcy_biases.push_back("_good_accept");

    TriggerModeSettings();

    out_file_suffix = trig_suffix + isForSoumya_suffix + qEtaBin_suffix;
    if (save_non_sepr_trg_hists) out_file_suffix += "_w_nonsepr";
    if (save_good_accept_trg_hists) out_file_suffix += "_w_good_accept";
}

void RDFBasedHistFillingData::FillHistograms(){
    std::cout << "Calling FillHistograms (cycle=" << hist_filling_cycle
              << ", trigger_effcy_calc=" << trigger_effcy_calc << ")" << std::endl;

    if (hist_filling_cycle == generic) {
        if (output_generic_hists) {
            FillHistogramsGeneric();
        }

        if (!trigger_effcy_calc) {
            // Pipeline 1 (nominal): crossx for signal pairs
            FillHistogramsCrossx();
        }

        if (trigger_effcy_calc) {
            // Pipeline 2 (trig effcy loop 1): single-muon no-correlation efficiency
            FillHistogramsSingleMuonEffcy();
        }
    } else if (hist_filling_cycle == inv_weight_by_single_mu_effcy && trigger_effcy_calc) {
        // Pipeline 3 (trig effcy loop 2): inverse-weighted dR corrections
        FillTrigEffcyHistsInvWeightedbySingleMuonEffcies();
    }
}

bool PassSingleMuonGapCut(float meta, float mpt, int mcharge){
    // parameters: eta and pT of the same muon
    if (fabs(meta) < ParamsSet::eta_gap_cut1) return false;
    if (mpt < 6){
        for (array<float,2> charge_eta_gap_cut : ParamsSet::charge_eta_gap_cuts){
            if (mcharge * meta > charge_eta_gap_cut[0] && mcharge * meta < charge_eta_gap_cut[1]) return false;
        }
    }
    // if (mpt < 6 && fabs(meta) > ParamsSet::eta_gap_cut2[0] && fabs(meta) < ParamsSet::eta_gap_cut2[1]) return false;
    return true;
}

bool MuPairPassGapCut(float m1eta, float m1pt, int m1charge, float m2eta, float m2pt, int m2charge){
    return PassSingleMuonGapCut(m1eta, m1pt, m1charge) && PassSingleMuonGapCut(m2eta, m2pt, m2charge);
}

void RDFBasedHistFillingData::FillHistogramsGeneric(){
    for (std::string category : categories_essential){
        std::string df_name = "df" + category;
        ROOT::RDF::RNode& df = map_at_checked(df_map, df_name, Form("FillHistogramsGeneric: df_map.at(%s)", df_name.c_str()));

        FillHistogramsSingleDataFrame(category, df);
        FillHistogramsSingleDataFrame(category, "_jacobian_corrected", df);

        if (output_gapcut_hists){
            ROOT::RDF::RNode df_wgapcut = df.Filter(MuPairPassGapCut, {"m1.eta", "m1.pt", "m1.charge", "m2.eta", "m2.pt", "m2.charge"});
            df_map.emplace(df_name + "_wgapcut", df_wgapcut);
            FillHistogramsSingleDataFrame(category + "_wgapcut", df_wgapcut);
            FillHistogramsSingleDataFrame(category + "_wgapcut", "_jacobian_corrected", df_wgapcut);
        }
    }
}

//--------- TRIGGER SETTINGS FROM TRIGGER MODE VARIABLE ---------
void RDFBasedHistFillingData::TriggerModeSettings(){

    switch(trigger_mode){
    case 0:
        base_trig_suffix = "_min_bias";
        trigs = {"_MB", "_mu4"};
        trigs_pair = {{"_mu4","_MB"}};
        break;
    case 1:
        trigs = {"_mu4", "_mu4_mu4noL1", "_2mu4", "_2mu4_AND_mu4_mu4noL1"};
        trigs_pair = {{"_mu4_mu4noL1","_mu4"},{"_2mu4","_mu4"},{"_2mu4_AND_mu4_mu4noL1","_mu4"}};
        break;
    case 2:
        base_trig_suffix = "_mu4_mu4noL1";
        break;
    case 3:
        base_trig_suffix = "_2mu4";
        break;
    default:
        std::cerr << "Trigger mode INVALID: must be 0 / 1 / 2 / 3!" << std::endl;
        throw std::exception();
    }

    // --- Derive trigger_effcy_calc, use_mu6_for_trg_eff ---
    const bool isRun3 = (run_year >= 22);
    trigger_effcy_calc = (trigger_mode == 0 || trigger_mode == 1)
        && !(IsPbPb() && isRun3 && mu4_nominal_pbpb_NO_trig_calc);
    use_mu6_for_trg_eff = trigger_effcy_calc && (trigger_mode == 1);
    doTrigEffcy = trigger_effcy_calc;

    // --- Set mu4 suffix AFTER mu4_nominal_pbpb_NO_trig_calc is finalized ---
    if (trigger_mode == 1) {
        base_trig_suffix = mu4_nominal_pbpb_NO_trig_calc ? "_mu4_nominal" : "_single_mu4";
    }

    // base_trig_suffix: trigger type used for input ntuple paths
    // trig_suffix: output-only, adds markers when trig effcy is off
    trig_suffix = base_trig_suffix;
    if (!doTrigEffcy && (trigger_mode == 0 || trigger_mode == 1)) trig_suffix += "_no_trg_plots";

    if (!filter_out_photo_resn_for_trig_effcy) trig_suffix += "_no_photo_resn_cuts";

    std::cout << "TriggerModeSettings: trigger_mode=" << trigger_mode
              << ", mu4_nominal_pbpb_NO_trig_calc=" << mu4_nominal_pbpb_NO_trig_calc
              << ", trigger_effcy_calc=" << trigger_effcy_calc
              << ", use_mu6_for_trg_eff=" << use_mu6_for_trg_eff
              << ", base_trig_suffix='" << base_trig_suffix << "'"
              << ", trig_suffix='" << trig_suffix << "'" << std::endl;

}

// ---------- ----------
void RDFBasedHistFillingData::InitOutput(){
    std::string file_mode = (hist_filling_cycle == generic)? "recreate" : "update";
    m_outfile = new TFile(output_file.c_str(), file_mode.c_str());
}

// ---------- ----------
void RDFBasedHistFillingData::BuildHistBinningMapDataImpl(){
    BuildHistBinningMapDataCommon();
    BuildHistBinningMapDataExtraHook();
}

void RDFBasedHistFillingData::BuildHistBinningMapDataCommon(){

    // ------- pT binning for single-muon trigger efficieny -------

    std::vector<double> pT_bins_single_muon (pms.pT_bins_8); // make a copy of a suitable set of single-muon pT bins (adjustable) --> use the copy for histogram settings

    pT_bins_single_muon.insert(pT_bins_single_muon.end(), pms.pT_bins_60.begin(), pms.pT_bins_60.end());

    hist_binning_map["pT_bins_single_muon"] = pT_bins_single_muon;

    // ------- eta binning for single-muon trigger efficiency -------

    std::vector<double> eta_bins_trig_effcy = ParamsSet::makeEtaTrigEffcyBinning(1);
    hist_binning_map["eta_bins_trig_effcy"] = eta_bins_trig_effcy;
    
    // ------- phi binning for single-muon trigger efficiency -------

    // Build uniform phi edges so we can use the (xbins, ybins, zbins) TH3D ctor
    int nphi_bins_trig_effcy = 128; // phi 2nd muon
    
    std::vector<double> phi2nd_bins(nphi_bins_trig_effcy + 1);
    for (int i = 0; i <= nphi_bins_trig_effcy; ++i) {
        phi2nd_bins[i] = -pms.PI + (2.0 * pms.PI) * (static_cast<double>(i) / nphi_bins_trig_effcy);
    }

    hist_binning_map["phi2nd_bins"] = phi2nd_bins;
}

// ---------- ----------
void RDFBasedHistFillingData::BuildFilterToVarListMapDataImpl(){
    // data-common filter to variable list maps
    BuildFilterToVarListMapDataCommon();
    // extras
    BuildFilterToVarListMapDataExtraHook();

    // trigger-efficiency specific filter to variable list map
    BuildTrgEffcyFilterToVarListMap();
}

void RDFBasedHistFillingData::BuildFilterToVarListMapDataCommon(){
    for (std::string sign : pair_signs){

        df_filter_to_var1D_list_map[sign + ""]         = {"pair_dPoverP", "Dphi", "DR", "DR_zoomin", "pair_y", "pt_asym", "pair_pt_ptlead_ratio"};
        df_filter_to_var1D_list_map[sign + "_wgapcut"] = {"pair_dPoverP", "Dphi", "DR", "DR_zoomin", "pair_y", "pt_asym", "pair_pt_ptlead_ratio"};

        df_filter_and_weight_to_var1D_list_map[{sign + ""           , "_jacobian_corrected"}] = {"DR", "DR_zoomin"};
        df_filter_and_weight_to_var1D_list_map[{sign + "_wgapcut"   , "_jacobian_corrected"}] = {"DR", "DR_zoomin"};

    }
}

void RDFBasedHistFillingData::BuildTrgEffcyFilterToVarListMap(){
    FlattenTrigEffcyFilters();
    BuildFlattenedTrgEffcyFilterToVarListMap();
}

void RDFBasedHistFillingData::BuildFlattenedTrgEffcyFilterToVarListMap(){
    BuildFlattenedTrgEffcyFilterToVarListMapDataCommon();
    BuildFlattenedTrgEffcyFilterToVarListMapExtra();
}

void RDFBasedHistFillingData::BuildFlattenedTrgEffcyFilterToVarListMapDataCommon(){
    for (auto filter : trg_effcy_filters_1D_pre_sum)    df_filter_to_var1D_list_map[filter] = single_muon_trig_effcy_var1Ds;
    for (auto filter : trg_effcy_filters_2D_3D_pre_sum) df_filter_to_var2D_list_map[filter] = single_muon_trig_effcy_var2Ds;
    for (auto filter : trg_effcy_filters_2D_3D_pre_sum) df_filter_to_var3D_list_map[filter] = single_muon_trig_effcy_var3Ds;
    // ----- ADD INVERSE WEIGHTING ONES!! WITH _mu4_mu4noL1_denom, _2mu4_denom (paired with additional filtering) -----

}

// Default no-op: all concrete derived classes override this.
// Needed explicitly by cling's JIT linker even though a compiled linker elides it
// for abstract classes whose virtual functions are all overridden.
void RDFBasedHistFillingData::FlattenTrigEffcyFiltersExtra() {}

void RDFBasedHistFillingData::FlattenTrigEffcyFilters()
{
    FlattenTrigEffcyFiltersDataCommon();
    FlattenTrigEffcyFiltersExtra();
}

void RDFBasedHistFillingData::FlattenTrigEffcyFiltersDataCommon()
{
    // flatten pre-sum levels
    HistFillUtils::flatten_levels(levels_trg_effcy_filters_1D_pre_sum, trg_effcy_filters_1D_pre_sum);
    HistFillUtils::flatten_levels(levels_trg_effcy_filters_2D_3D_pre_sum, trg_effcy_filters_2D_3D_pre_sum);

    // build post-sum levels
    HistFillUtils::write_post_sum_levels(levels_trg_effcy_filters_1D_pre_sum,
                          levels_trg_effcy_to_be_summed,
                          levels_trg_effcy_filters_1D_post_sum);

    HistFillUtils::write_post_sum_levels(levels_trg_effcy_filters_2D_3D_pre_sum,
                          levels_trg_effcy_to_be_summed,
                          levels_trg_effcy_filters_2D_3D_post_sum);

    // flatten post-sum levels
    HistFillUtils::flatten_levels(levels_trg_effcy_filters_1D_post_sum, trg_effcy_filters_1D_post_sum);
    HistFillUtils::flatten_levels(levels_trg_effcy_filters_2D_3D_post_sum, trg_effcy_filters_2D_3D_post_sum);

    // build to-be-summed levels
    for (int level_ind = 0;
         level_ind < static_cast<int>(levels_trg_effcy_filters_2D_3D_pre_sum.size());
         ++level_ind)
    {
        if (std::find(levels_trg_effcy_to_be_summed.begin(),
                      levels_trg_effcy_to_be_summed.end(),
                      level_ind) != levels_trg_effcy_to_be_summed.end())
        {
            levels_trg_effcy_filters_to_be_summed.push_back(
                levels_trg_effcy_filters_2D_3D_pre_sum.at(level_ind)
            );
        }
    }

    // flatten to-be-summed levels
    HistFillUtils::flatten_levels(levels_trg_effcy_filters_to_be_summed, trg_effcy_filters_to_be_summed);
}

//--------- DATA HIST POST PROCESSING ---------
void RDFBasedHistFillingData::HistPostProcessDataImpl(){
    HistPostProcessDataCommon();
    HistPostProcessDataExtra();
}

void RDFBasedHistFillingData::HistPostProcessDataCommon(){
    if (!trigger_effcy_calc) return; // Pipeline 1 (nominal): no trigger efficiency post-processing

    if (hist_filling_cycle == generic) {
        // Pipeline 2: project & divide trigger efficiency histograms
        SumSingleMuonTrigEffHists();
        MakeAndWriteSingleMuonTrigEffPtGraphs();
        if (!isForSoumya) CalculateSingleMuonTrigEffcyRatios();
    } else if (hist_filling_cycle == inv_weight_by_single_mu_effcy) {
        // Pipeline 3: dR correction graphs from inverse-weighted histograms
        MakeAndWriteDRTrigEffGraphs();
    }
}

//--------- SUMMING 1,2,3D TRIGGER EFFICIENCY HISTOGRAMS ---------
void RDFBasedHistFillingData::SumSingleMuonTrigEffHists(){
    SumSingleMuonTrigEffHistsDataCommon();
    SumSingleMuonTrigEffHistsExtra();
}

void RDFBasedHistFillingData::SumSingleMuonTrigEffHistsDataCommon(){

    // 1D
    HistFillUtils::SumTrigEffHistsGeneric<TH1D, std::string>(
        single_muon_trig_effcy_var1Ds,
        trg_effcy_filters_1D_post_sum,
        trg_effcy_filters_to_be_summed,
        hist1D_map,
        [](const std::string& var) {
            return "h_" + var;
        }
    );

    // 2D
    HistFillUtils::SumTrigEffHistsGeneric<TH2D, std::array<std::string,2>>(
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
    HistFillUtils::SumTrigEffHistsGeneric<TH3D, std::array<std::string,3>>(
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

//--------- MAKE & WRITE SINGLE-MUON TRIG EFFICIENCY GRAPHS AS FUNCTION OF PT ---------
void RDFBasedHistFillingData::MakeAndWriteSingleMuonTrigEffPtGraphsHelper(const std::vector<std::string>& categories){
	std::cout << "Calling MakeAndWriteSingleMuonTrigEffPtGraphs" << std::endl;

    if (categories.empty()){

        TH2D* h_pt2nd_vs_q_eta2nd_mu4_mu4noL1_sepr  = map_at_checked(hist2D_map, "h_pt2nd_vs_q_eta2nd_mu4_mu4noL1_sepr",    "MakeAndWriteSingleMuonTrigEffPtGraphs: hist2D_map.at(h_pt2nd_vs_q_eta2nd_mu4_mu4noL1_sepr)");
        TH2D* h_pt2nd_vs_q_eta2nd_mu4_sepr          = map_at_checked(hist2D_map, "h_pt2nd_vs_q_eta2nd_mu4_sepr",            "MakeAndWriteSingleMuonTrigEffPtGraphs: hist2D_map.at(h_pt2nd_vs_q_eta2nd_mu4_sepr)");
        TH2D* h_pt2nd_vs_q_eta2nd_2mu4_sepr         = map_at_checked(hist2D_map, "h_pt2nd_vs_q_eta2nd_2mu4_sepr",           "MakeAndWriteSingleMuonTrigEffPtGraphs: hist2D_map.at(h_pt2nd_vs_q_eta2nd_2mu4_sepr)");
        TH2D* h_pt2nd_vs_q_eta2nd_2mu4_AND_mu4_mu4noL1_sepr         = map_at_checked(hist2D_map, "h_pt2nd_vs_q_eta2nd_2mu4_AND_mu4_mu4noL1_sepr",           "MakeAndWriteSingleMuonTrigEffPtGraphs: hist2D_map.at(h_pt2nd_vs_q_eta2nd_2mu4_AND_mu4_mu4noL1_sepr)");

        if (!h_pt2nd_vs_q_eta2nd_mu4_mu4noL1_sepr || !h_pt2nd_vs_q_eta2nd_mu4_sepr || !h_pt2nd_vs_q_eta2nd_2mu4_sepr || !h_pt2nd_vs_q_eta2nd_2mu4_AND_mu4_mu4noL1_sepr){
            std::cout << "MakeAndWriteSingleMuonTrigEffPtGraphsHelper: Histogram pointers are nullptr! Return without making graphs!" << std::endl;
            return;
        }

        // ----- extract q-eta bins (might be category dependent) -----
        std::vector<double> eta_bins_trig_effcy = {};

        const TAxis* xaxis = h_pt2nd_vs_q_eta2nd_mu4_sepr->GetXaxis();
        int nbins = xaxis->GetNbins();
        const double* arr = xaxis->GetXbins()->GetArray();
        eta_bins_trig_effcy.assign(arr, arr + nbins + 1);

        for (auto range : q_eta_proj_ranges){
            int bin_first = bin_number(range.first, eta_bins_trig_effcy) + 1; // + 1 pushes into next bin (bin lower end agree with range edge if range edge founded)
            int bin_last = bin_number(range.second, eta_bins_trig_effcy); // bin higher end agree with range edge if range edge founded
            
            std::string proj_suffix = pairToSuffix(range);

            if (isForSoumya){
                HistFillUtils::proj_and_write(h_pt2nd_vs_q_eta2nd_mu4_sepr,                    true, bin_first, bin_last, proj_suffix, &hist1D_map);
                HistFillUtils::proj_and_write(h_pt2nd_vs_q_eta2nd_mu4_mu4noL1_sepr,            true, bin_first, bin_last, proj_suffix, &hist1D_map);
                HistFillUtils::proj_and_write(h_pt2nd_vs_q_eta2nd_2mu4_sepr,                   true, bin_first, bin_last, proj_suffix, &hist1D_map);
                HistFillUtils::proj_and_write(h_pt2nd_vs_q_eta2nd_2mu4_AND_mu4_mu4noL1_sepr,   true, bin_first, bin_last, proj_suffix, &hist1D_map);
            } else{

                HistFillUtils::proj_divide_and_write(h_pt2nd_vs_q_eta2nd_mu4_mu4noL1_sepr,           h_pt2nd_vs_q_eta2nd_mu4_sepr, true, bin_first, bin_last, proj_suffix, &graph_map);
                HistFillUtils::proj_divide_and_write(h_pt2nd_vs_q_eta2nd_2mu4_sepr,                  h_pt2nd_vs_q_eta2nd_mu4_sepr, true, bin_first, bin_last, proj_suffix, &graph_map);
                HistFillUtils::proj_divide_and_write(h_pt2nd_vs_q_eta2nd_2mu4_AND_mu4_mu4noL1_sepr,  h_pt2nd_vs_q_eta2nd_mu4_sepr, true, bin_first, bin_last, proj_suffix, &graph_map);                
            }
        }

        return;
    }

    for (std::string category : categories){
        TH2D* h_pt2nd_vs_q_eta2nd_mu4_mu4noL1_sepr  = map_at_checked(hist2D_map, "h_pt2nd_vs_q_eta2nd" + category + "_mu4_mu4noL1_sepr",    Form("MakeAndWriteSingleMuonTrigEffPtGraphs: hist2D_map.at(%s)", ("h_pt2nd_vs_q_eta2nd" + category + "_mu4_mu4noL1_sepr").c_str()));
        TH2D* h_pt2nd_vs_q_eta2nd_mu4_sepr          = map_at_checked(hist2D_map, "h_pt2nd_vs_q_eta2nd" + category + "_mu4_sepr",            Form("MakeAndWriteSingleMuonTrigEffPtGraphs: hist2D_map.at(%s)", ("h_pt2nd_vs_q_eta2nd" + category + "_mu4_sepr").c_str()));
        TH2D* h_pt2nd_vs_q_eta2nd_2mu4_sepr         = map_at_checked(hist2D_map, "h_pt2nd_vs_q_eta2nd" + category + "_2mu4_sepr",           Form("MakeAndWriteSingleMuonTrigEffPtGraphs: hist2D_map.at(%s)", ("h_pt2nd_vs_q_eta2nd" + category + "_2mu4_sepr").c_str()));
        TH2D* h_pt2nd_vs_q_eta2nd_2mu4_AND_mu4_mu4noL1_sepr         = map_at_checked(hist2D_map, "h_pt2nd_vs_q_eta2nd" + category + "_2mu4_AND_mu4_mu4noL1_sepr",           Form("MakeAndWriteSingleMuonTrigEffPtGraphs: hist2D_map.at(%s)", ("h_pt2nd_vs_q_eta2nd" + category + "_2mu4_AND_mu4_mu4noL1_sepr").c_str()));

        if (!h_pt2nd_vs_q_eta2nd_mu4_mu4noL1_sepr || !h_pt2nd_vs_q_eta2nd_mu4_sepr || !h_pt2nd_vs_q_eta2nd_2mu4_sepr || !h_pt2nd_vs_q_eta2nd_2mu4_AND_mu4_mu4noL1_sepr) continue;

        // HistFillUtils::proj_divide_and_write(h_pt2nd_vs_q_eta2nd_mu4_mu4noL1_sepr,           h_pt2nd_vs_q_eta2nd_mu4_sepr, &graph_map);
        // HistFillUtils::proj_divide_and_write(h_pt2nd_vs_q_eta2nd_2mu4_sepr,                  h_pt2nd_vs_q_eta2nd_mu4_sepr, &graph_map);
        // HistFillUtils::proj_divide_and_write(h_pt2nd_vs_q_eta2nd_2mu4_AND_mu4_mu4noL1_sepr,  h_pt2nd_vs_q_eta2nd_mu4_sepr, &graph_map);

        // ----- extract q-eta bins (might be category dependent) -----
        std::vector<double> eta_bins_trig_effcy = {};

        const TAxis* xaxis = h_pt2nd_vs_q_eta2nd_mu4_sepr->GetXaxis();
        int nbins = xaxis->GetNbins();
        const double* arr = xaxis->GetXbins()->GetArray();
        eta_bins_trig_effcy.assign(arr, arr + nbins + 1);

        for (auto range : q_eta_proj_ranges){
            int bin_first = bin_number(range.first, eta_bins_trig_effcy) + 1; // + 1 pushes into next bin (bin lower end agree with range edge if range edge founded)
            int bin_last = bin_number(range.second, eta_bins_trig_effcy); // bin higher end agree with range edge if range edge founded
            std::string proj_suffix = pairToSuffix(range);

            if (isForSoumya){
                HistFillUtils::proj_and_write(h_pt2nd_vs_q_eta2nd_mu4_sepr,                    true, bin_first, bin_last, proj_suffix, &hist1D_map);
                HistFillUtils::proj_and_write(h_pt2nd_vs_q_eta2nd_mu4_mu4noL1_sepr,            true, bin_first, bin_last, proj_suffix, &hist1D_map);
                HistFillUtils::proj_and_write(h_pt2nd_vs_q_eta2nd_2mu4_sepr,                   true, bin_first, bin_last, proj_suffix, &hist1D_map);
                HistFillUtils::proj_and_write(h_pt2nd_vs_q_eta2nd_2mu4_AND_mu4_mu4noL1_sepr,   true, bin_first, bin_last, proj_suffix, &hist1D_map);
            } else{
                HistFillUtils::proj_divide_and_write(h_pt2nd_vs_q_eta2nd_mu4_mu4noL1_sepr,           h_pt2nd_vs_q_eta2nd_mu4_sepr, true, bin_first, bin_last, proj_suffix, &graph_map);
                HistFillUtils::proj_divide_and_write(h_pt2nd_vs_q_eta2nd_2mu4_sepr,                  h_pt2nd_vs_q_eta2nd_mu4_sepr, true, bin_first, bin_last, proj_suffix, &graph_map);
                HistFillUtils::proj_divide_and_write(h_pt2nd_vs_q_eta2nd_2mu4_AND_mu4_mu4noL1_sepr,  h_pt2nd_vs_q_eta2nd_mu4_sepr, true, bin_first, bin_last, proj_suffix, &graph_map);
            }
        }
    }
}

void RDFBasedHistFillingData::CalculateSingleMuonTrigEffcyRatiosHelper(
    const std::vector<std::string>& categories)
{
    auto makeTrigEffRatios = [this, &categories](const std::string& base,
                                                 auto& hist_map)
    {
        using MapT = std::decay_t<decltype(hist_map)>;           // remove & and cv
        using PtrT = typename MapT::mapped_type;                 // TH2D* or TH3D*
        using TH   = std::remove_pointer_t<PtrT>;                // TH2D or TH3D

        for (const auto& trg_pair : trigs_pair) {
            for (const auto& category : categories) {
                std::string hname_num = base + category + trg_pair.first  + "_sepr";
                std::string hname_dom = base + category + trg_pair.second + "_sepr";

                TH* h_num = map_at_checked(hist_map, hname_num,
                    Form("CalculateSingleMuonTrigEffcyRatios: hist_map.at(%s)",
                         hname_num.c_str()));

                TH* h_dom = map_at_checked(hist_map, hname_dom,
                    Form("CalculateSingleMuonTrigEffcyRatios: hist_map.at(%s)",
                         hname_dom.c_str()));

                TH* h_divided = static_cast<TH*>(
                    h_num->Clone(Form("%s_divided", hname_num.c_str()))
                );
                h_divided->Divide(h_dom);

                hist_map[h_divided->GetName()] = h_divided;
            }
        }
    };

    std::string h2d_base = "h_pt2nd_vs_q_eta2nd";
    std::string h3d_base = "h_pt2nd_vs_q_eta2nd_vs_phi2nd";

    makeTrigEffRatios(h2d_base, hist2D_map);  // MapT::mapped_type = TH2D*
    makeTrigEffRatios(h3d_base, hist3D_map);  // MapT::mapped_type = TH3D*
}


std::string RDFBasedHistFillingData::FindBinReturnStr(float number, const std::vector<std::pair<float, float>>& ranges) {
    for (const auto& range : ranges) {
        if (number >= range.first && number < range.second) {
            return pairToSuffix(range);
        }
    }
    return "";
}

std::string RDFBasedHistFillingData::FindCtrSuffix(int centrality) {
    if (centrality < 0) return "";
    static const std::vector<std::pair<int, int>> ctr_edges = {
        {0, 5}, {5, 10}, {10, 20}, {20, 30}, {30, 50}, {50, 80}
    };
    for (const auto& [lo, hi] : ctr_edges) {
        if (centrality >= lo && centrality < hi)
            return "_ctr" + std::to_string(lo) + "_" + std::to_string(hi);
    }
    return "";
}

float RDFBasedHistFillingData::EvaluateSingleMuonEffcyPtFitted(const std::string& ctr_suffix, bool charge_positive, float pt, float q_eta) {
    static const CommonEffcyConfig cfg{};
    std::string musign = charge_positive ? "_sign1" : "_sign2";
    std::string q_eta_suffix = FindBinReturnStr(q_eta, cfg.q_eta_proj_ranges_fine_excl_gap);

    if (!q_eta_suffix.empty()) {
        std::string key = "f_pt2nd_vs_q_eta2nd" + ctr_suffix + musign + "_2mu4_sepr_py_" + q_eta_suffix + "_divided";
        auto it = s_effcy_pT_fit_map.find(key);
        if (it != s_effcy_pT_fit_map.end()) {
            double val = it->second->Eval(pt);
            if (val < 0.01) val = 0.01;
            if (val > 1.0) val = 1.0;
            return static_cast<float>(val);
        }
    }

    // Fallback: unfitted 2D histogram for gap regions or missing TF1
    std::string h2d_key = "h_pt2nd_vs_q_eta2nd" + ctr_suffix + musign + "_2mu4_sepr_divided";
    auto it2d = s_effcy_2D_hist_map.find(h2d_key);
    if (it2d != s_effcy_2D_hist_map.end()) {
        int bin = it2d->second->FindBin(q_eta, pt);
        double val = it2d->second->GetBinContent(bin);
        if (val > 0.01 && val <= 1.0) return static_cast<float>(val);
    }
    return -1.0f;
}

float RDFBasedHistFillingData::EvaluateSingleMuonEffcy(const std::string& ctr_suffix, bool charge_positive, float pt, float q_eta) {
    return EvaluateSingleMuonEffcyPtFitted(ctr_suffix, charge_positive, pt, q_eta);
}

void RDFBasedHistFillingData::MakeAndWriteDRTrigEffGraphsHelper(const std::vector<std::string>& categories) {
    std::vector<std::string> invw_var1Ds = {
        "DR", "DR_zoomin", "DR_0_2", "Deta", "Deta_zoomin",
        "Dphi", "Dphi_zoomin", "pair_pt_log", "minv_zoomin"
    };

    std::string pair_suffix = IsPbPb() ? "_pair_mu4" : "_2mu4";

    std::vector<std::string> cats = categories.empty() ? std::vector<std::string>{""} : categories;

    for (const std::string& cat : cats) {
        for (const std::string& pair_sign : {"_ss", "_op"}) {
            // 1D ratio graphs
            for (const auto& var : invw_var1Ds) {
                std::string hname_num = "h_" + var + pair_sign + pair_suffix + "_invw_num" + cat;
                std::string hname_den = "h_" + var + pair_sign + pair_suffix + "_denom" + cat;
                TH1D* h_num = map_at_checked(hist1D_map, hname_num, Form("MakeAndWriteDRTrigEffGraphsHelper: %s", hname_num.c_str()));
                TH1D* h_den = map_at_checked(hist1D_map, hname_den, Form("MakeAndWriteDRTrigEffGraphsHelper: %s", hname_den.c_str()));
                if (h_num && h_den)
                    HistFillUtils::ratio_divide_and_write(h_num, h_den, &graph_map);
            }

            // Pair-pT-binned dR corrections from 2D histograms
            struct PtSlice { double lo, hi; std::string label; };
            std::vector<PtSlice> pt_slices = {
                {8.0, 12.0, "_pt8_12"}, {12.0, 20.0, "_pt12_20"},
                {20.0, 40.0, "_pt20_40"}, {40.0, 120.0, "_pt40_120"}
            };
            for (const std::string& dr_var : {"DR_zoomin", "DR"}) {
                std::string h2d_num_name = "h_pair_pt_log_vs_" + dr_var + pair_sign + pair_suffix + "_invw_num" + cat;
                std::string h2d_den_name = "h_pair_pt_log_vs_" + dr_var + pair_sign + pair_suffix + "_denom" + cat;
                auto it_num = hist2D_map.find(h2d_num_name);
                auto it_den = hist2D_map.find(h2d_den_name);
                if (it_num == hist2D_map.end() || it_den == hist2D_map.end()) continue;
                TH2D* h2_num = it_num->second;
                TH2D* h2_den = it_den->second;
                for (const auto& sl : pt_slices) {
                    int ybin_lo = h2_num->GetYaxis()->FindBin(sl.lo + 0.001);
                    int ybin_hi = h2_num->GetYaxis()->FindBin(sl.hi - 0.001);
                    std::string proj_num_name = "h_" + dr_var + pair_sign + pair_suffix + "_invw_num" + cat + sl.label;
                    std::string proj_den_name = "h_" + dr_var + pair_sign + pair_suffix + "_denom" + cat + sl.label;
                    TH1D* h_proj_num = h2_num->ProjectionX(proj_num_name.c_str(), ybin_lo, ybin_hi);
                    TH1D* h_proj_den = h2_den->ProjectionX(proj_den_name.c_str(), ybin_lo, ybin_hi);
                    hist1D_map[proj_num_name] = h_proj_num;
                    hist1D_map[proj_den_name] = h_proj_den;
                    HistFillUtils::ratio_divide_and_write(h_proj_num, h_proj_den, &graph_map);
                }
            }
        }
    }
}

void RDFBasedHistFillingData::WriteOutputExtra(){
    HistFillUtils::write_hist_map_vector(graph_map, graphs_to_not_write);
}

void RDFBasedHistFillingData::CleanupExtra(){
    for (auto g : graph_map) delete g.second;
}