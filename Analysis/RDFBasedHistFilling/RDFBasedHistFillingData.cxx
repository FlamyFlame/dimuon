#include "RDFBasedHistFillingData.h"
#include <type_traits>


RDFBasedHistFillingData::RDFBasedHistFillingData(int run_year_input, bool isForSoumya_input)
: run_year (run_year_input), isForSoumya (isForSoumya_input){
    std::cout << " Histogram filling for data:" << std::endl << std::endl;

    std::cout << "The following public variable(s) **MUST** be checked:" << std::endl;
    std::cout << "--> trigger_mode: [INT] default 1 (single_mu4)" << std::endl;
    std::cout << "--> hist_filling_cycle: [INT] default generic" << std::endl << std::endl;
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

    q_eta_proj_ranges = (run_year > 20)
        ? (useCoarseQEtaBin ? q_eta_proj_ranges_coarse_incl_gap      : q_eta_proj_ranges_fine_excl_gap)
        : (useCoarseQEtaBin ? q_eta_proj_ranges_coarse_incl_gap_run2 : q_eta_proj_ranges_fine_excl_gap_run2);

    // fill q_eta_ranges_str needed for child class Initialize() overriden part
    q_eta_ranges_str.clear();
    for (auto pair : q_eta_proj_ranges){
        q_eta_ranges_str.push_back("_q_eta_" + pairToSuffix(pair));
    }

    if (isForSoumya){
        single_muon_trig_effcy_var1Ds = {};
        single_muon_trig_effcy_var2Ds = {{"q_eta2nd","pt2nd"}};
        single_muon_trig_effcy_var3Ds = {};
    }

    trig_to_filter_str_map["_mu4"] = "";
    trig_to_filter_str_map["_mu4_mu4noL1"] = "passmu4mu4noL1";
    trig_to_filter_str_map["_2mu4"] = "pass2mu4";
    trig_to_filter_str_map["_2mu4_AND_mu4_mu4noL1"] = "pass2mu4 && passmu4mu4noL1";

    trg_effcy_biases = {"_sepr"};
    if (save_non_sepr_trg_hists) trg_effcy_biases.push_back("");
    if (save_good_accept_trg_hists) trg_effcy_biases.push_back("_good_accept");

    TriggerModeSettings();

    out_file_suffix = trig_suffix + isForSoumya_suffix + qEtaBin_suffix;
    if (save_non_sepr_trg_hists) out_file_suffix += "_w_nonsepr";
    if (save_good_accept_trg_hists) out_file_suffix += "_w_good_accept";
}

void RDFBasedHistFillingData::FillHistograms(){    
    std::cout << "Calling FillHistograms" << std::endl;
    
    // ------- Fill histograms -------
    if (hist_filling_cycle == generic){
        if (output_generic_hists){
            FillHistogramsGeneric();
        }
        if (doTrigEffcy){
            FillHistogramsSingleMuonEffcy();
        }
    } else if (hist_filling_cycle == inv_weight_by_single_mu_effcy && trigger_mode == 1){
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
        trig_suffix = "_min_bias";
        if (!doTrigEffcy) trig_suffix += "_no_trg_plots";
        trigs = {"_MB", "_mu4"};
        trigs_pair = {{"_mu4","_MB"}};
        break;
    case 1:
        trig_suffix = "_single_mu4";
        if (!doTrigEffcy) trig_suffix += "_no_trg_plots";
        trigs = {"_mu4", "_mu4_mu4noL1", "_2mu4", "_2mu4_AND_mu4_mu4noL1"}; // can be overwritten if, e.g, mu4_mu4noL1 doesn't exist
        trigs_pair = {{"_mu4_mu4noL1","_mu4"},{"_2mu4","_mu4"},{"_2mu4_AND_mu4_mu4noL1","_mu4"}};
        break;
    case 2:
        trig_suffix = "_mu4_mu4noL1";
        break;
    case 3:
        trig_suffix = "_2mu4";
        break;        
    default:
        std::cerr << "Trigger mode INVALID: must be 0 / 1 / 2 / 3!" << std::endl;
        throw std::exception();
    }

    doTrigEffcy = (doTrigEffcy && (trigger_mode == 0 || trigger_mode == 1)); // turn off trig effcy plotting if trigger_mode isn't 0 or 1

    if (!filter_out_photo_resn_for_trig_effcy) trig_suffix += "_no_photo_resn_cuts"; // if not filter out photoprod/resn pairs for trigger efficiency study

}

// ---------- ----------
void RDFBasedHistFillingData::InitOutput(){
    std::string file_mode = (hist_filling_cycle == generic)? "recreate" : "update";
    m_outfile = new TFile(output_file.c_str(), file_mode.c_str());
}

// ---------- ----------
void RDFBasedHistFillingData::BuildHistBinningMap(){
    
    RDFBasedHistFillingBaseClass::BuildHistBinningMap();

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

    // trigger-efficiency specific filter to variable list map
    BuildTrgEffcyFilterToVarListMap();
}

void RDFBasedHistFillingData::BuildFilterToVarListMapDataCommon(){
    for (std::string sign : pair_signs){

        df_filter_to_var1D_list_map[sign + ""]         = {"pair_dP_overP", "Dphi", "DR", "DR_zoomin", "pair_y", "pt_asym", "pair_pt_ptlead_ratio"};
        df_filter_to_var1D_list_map[sign + "_wgapcut"] = {"pair_dP_overP", "Dphi", "DR", "DR_zoomin", "pair_y", "pt_asym", "pair_pt_ptlead_ratio"};

        df_filter_and_weight_to_var1D_list_map[{sign + ""           , "_jacobian_corrected"}] = {"DR", "DR_zoomin"};
        df_filter_and_weight_to_var1D_list_map[{sign + "_wgapcut"   , "_jacobian_corrected"}] = {"DR", "DR_zoomin"};

    }
}

void RDFBasedHistFillingData::BuildTrgEffcyFilterToVarListMap(){
    TrigEffcyFiltersPrePostSumFlattening();
    BuildFlattenedTrgEffcyFilterToVarListMap();
}

void RDFBasedHistFillingData::BuildFlattenedTrgEffcyFilterToVarListMap(){
    BuildFlattenedTrgEffcyFilterToVarListMapDataCommon();
    BuildFlattenedTrgEffcyFilterToVarListMapExtra();
}

void BuildFlattenedTrgEffcyFilterToVarListMapDataCommon(){
    for (auto filter : trg_effcy_filters_1D_pre_sum)    df_filter_to_var1D_list_map[filter] = single_muon_trig_effcy_var1Ds;
    for (auto filter : trg_effcy_filters_2D_3D_pre_sum) df_filter_to_var2D_list_map[filter] = single_muon_trig_effcy_var2Ds;
    for (auto filter : trg_effcy_filters_2D_3D_pre_sum) df_filter_to_var3D_list_map[filter] = single_muon_trig_effcy_var3Ds;
    // ----- ADD INVERSE WEIGHTING ONES!! WITH _mu4_mu4noL1_denom, _2mu4_denom (paired with additional filtering) -----

}

void RDFBasedHistFillingData::TrigEffcyFiltersPrePostSumFlattening()
{
    TrigEffcyFiltersPrePostSumFlatteningDataCommon();
    TrigEffcyFiltersPrePostSumFlatteningExtra();
}

void RDFBasedHistFillingData::TrigEffcyFiltersPrePostSumFlatteningDataCommon()
{
    // flatten pre-sum levels
    TrigEffcyUtils::flatten_levels(levels_trg_effcy_filters_1D_pre_sum, trg_effcy_filters_1D_pre_sum);
    TrigEffcyUtils::flatten_levels(levels_trg_effcy_filters_2D_3D_pre_sum, trg_effcy_filters_2D_3D_pre_sum);

    // build post-sum levels
    TrigEffcyUtils::write_post_sum_levels(levels_trg_effcy_filters_1D_pre_sum,
                          levels_trg_effcy_to_be_summed,
                          levels_trg_effcy_filters_1D_post_sum);

    TrigEffcyUtils::write_post_sum_levels(levels_trg_effcy_filters_2D_3D_pre_sum,
                          levels_trg_effcy_to_be_summed,
                          levels_trg_effcy_filters_2D_3D_post_sum);

    // flatten post-sum levels
    TrigEffcyUtils::flatten_levels(levels_trg_effcy_filters_1D_post_sum, trg_effcy_filters_1D_post_sum);
    TrigEffcyUtils::flatten_levels(levels_trg_effcy_filters_2D_3D_post_sum, trg_effcy_filters_2D_3D_post_sum);

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
    TrigEffcyUtils::flatten_levels(levels_trg_effcy_filters_to_be_summed, trg_effcy_filters_to_be_summed);
}

//--------- DATA HIST POST PROCESSING ---------
void RDFBasedHistFillingData::HistPostProcessDataImpl(){
    HistPostProcessDataCommon();
    HistPostProcessDataExtra();
}

void RDFBasedHistFillingData::HistPostProcessDataCommon(){
    if (trigger_mode == 0 || trigger_mode == 1){
        if (hist_filling_cycle == generic){
            SumSingleMuonTrigEffHists();
            MakeAndWriteSingleMuonTrigEffPtGraphs();
            if (!isForSoumya) CalculateSingleMuonTrigEffcyRatios();
        } else if (hist_filling_cycle == inv_weight_by_single_mu_effcy){
            MakeAndWriteDRTrigEffGraphs();
        }
    }
}

//--------- SUMMING 1,2,3D TRIGGER EFFICIENCY HISTOGRAMS ---------
void RDFBasedHistFillingData::SumSingleMuonTrigEffHists(){

    // 1D
    TrigEffcyUtils::SumTrigEffHistsGeneric<TH1D, std::string>(
        single_muon_trig_effcy_var1Ds,
        trg_effcy_filters_1D_post_sum,
        trg_effcy_filters_to_be_summed,
        hist1D_map,
        [](const std::string& var) {
            return "h_" + var;
        }
    );

    // 2D
    TrigEffcyUtils::SumTrigEffHistsGeneric<TH2D, std::array<std::string,2>>(
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
    TrigEffcyUtils::SumTrigEffHistsGeneric<TH3D, std::array<std::string,3>>(
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
                TrigEffcyUtils::proj_and_write(h_pt2nd_vs_q_eta2nd_mu4_sepr,                    true, bin_first, bin_last, proj_suffix, &hist1D_map);
                TrigEffcyUtils::proj_and_write(h_pt2nd_vs_q_eta2nd_mu4_mu4noL1_sepr,            true, bin_first, bin_last, proj_suffix, &hist1D_map);
                TrigEffcyUtils::proj_and_write(h_pt2nd_vs_q_eta2nd_2mu4_sepr,                   true, bin_first, bin_last, proj_suffix, &hist1D_map);
                TrigEffcyUtils::proj_and_write(h_pt2nd_vs_q_eta2nd_2mu4_AND_mu4_mu4noL1_sepr,   true, bin_first, bin_last, proj_suffix, &hist1D_map);
            } else{

                TrigEffcyUtils::proj_divide_and_write(h_pt2nd_vs_q_eta2nd_mu4_mu4noL1_sepr,           h_pt2nd_vs_q_eta2nd_mu4_sepr, true, bin_first, bin_last, proj_suffix, &graph_map);
                TrigEffcyUtils::proj_divide_and_write(h_pt2nd_vs_q_eta2nd_2mu4_sepr,                  h_pt2nd_vs_q_eta2nd_mu4_sepr, true, bin_first, bin_last, proj_suffix, &graph_map);
                TrigEffcyUtils::proj_divide_and_write(h_pt2nd_vs_q_eta2nd_2mu4_AND_mu4_mu4noL1_sepr,  h_pt2nd_vs_q_eta2nd_mu4_sepr, true, bin_first, bin_last, proj_suffix, &graph_map);                
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

        // TrigEffcyUtils::proj_divide_and_write(h_pt2nd_vs_q_eta2nd_mu4_mu4noL1_sepr,           h_pt2nd_vs_q_eta2nd_mu4_sepr, &graph_map);
        // TrigEffcyUtils::proj_divide_and_write(h_pt2nd_vs_q_eta2nd_2mu4_sepr,                  h_pt2nd_vs_q_eta2nd_mu4_sepr, &graph_map);
        // TrigEffcyUtils::proj_divide_and_write(h_pt2nd_vs_q_eta2nd_2mu4_AND_mu4_mu4noL1_sepr,  h_pt2nd_vs_q_eta2nd_mu4_sepr, &graph_map);

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
                TrigEffcyUtils::proj_and_write(h_pt2nd_vs_q_eta2nd_mu4_sepr,                    true, bin_first, bin_last, proj_suffix, &hist1D_map);
                TrigEffcyUtils::proj_and_write(h_pt2nd_vs_q_eta2nd_mu4_mu4noL1_sepr,            true, bin_first, bin_last, proj_suffix, &hist1D_map);
                TrigEffcyUtils::proj_and_write(h_pt2nd_vs_q_eta2nd_2mu4_sepr,                   true, bin_first, bin_last, proj_suffix, &hist1D_map);
                TrigEffcyUtils::proj_and_write(h_pt2nd_vs_q_eta2nd_2mu4_AND_mu4_mu4noL1_sepr,   true, bin_first, bin_last, proj_suffix, &hist1D_map);
            } else{
                TrigEffcyUtils::proj_divide_and_write(h_pt2nd_vs_q_eta2nd_mu4_mu4noL1_sepr,           h_pt2nd_vs_q_eta2nd_mu4_sepr, true, bin_first, bin_last, proj_suffix, &graph_map);
                TrigEffcyUtils::proj_divide_and_write(h_pt2nd_vs_q_eta2nd_2mu4_sepr,                  h_pt2nd_vs_q_eta2nd_mu4_sepr, true, bin_first, bin_last, proj_suffix, &graph_map);
                TrigEffcyUtils::proj_divide_and_write(h_pt2nd_vs_q_eta2nd_2mu4_AND_mu4_mu4noL1_sepr,  h_pt2nd_vs_q_eta2nd_mu4_sepr, true, bin_first, bin_last, proj_suffix, &graph_map);
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


//--------- PLACEHOLDER ---------
void RDFBasedHistFillingData::MakeAndWriteDRTrigEffGraphsHelper(const std::vector<std::string>& categories){}

std::string RDFBasedHistFillingData::FindBinReturnStr(float number, const std::vector<std::pair<float, float>>& ranges){return "";}

float RDFBasedHistFillingData::EvaluateSingleMuonEffcyPtFitted(bool charge_sign, std::string trg, float pt_2nd, float q_eta_2nd, float phi_2nd){return 0;}

float RDFBasedHistFillingData::EvaluateSingleMuonEffcy(bool charge_sign, std::string trg, float pt_2nd, float q_eta_2nd, float phi_2nd){return 0;}

void RDFBasedHistFillingData::WriteOutputExtra(){
    TrigEffcyUtils::write_hist_map_vector(graph_map, graphs_to_not_write);
}

void RDFBasedHistFillingData::CleanupExtra(){
    for (auto g : graph_map) delete g.second;
}