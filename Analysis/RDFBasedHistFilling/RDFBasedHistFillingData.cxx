#include "RDFBasedHistFillingData.h"
#include <type_traits>

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
void RDFBasedHistFillingData::BuildFilterToVarListMap(){

    RDFBasedHistFillingBaseClass::BuildFilterToVarListMap();

    for (std::string sign : pair_signs){

        df_filter_to_var1D_list_map[sign + ""]         = {"pair_dP_overP", "Dphi", "DR", "DR_zoomin", "pair_y", "pt_asym", "pair_pt_ptlead_ratio"};
        df_filter_to_var1D_list_map[sign + "_wgapcut"] = {"pair_dP_overP", "Dphi", "DR", "DR_zoomin", "pair_y", "pt_asym", "pair_pt_ptlead_ratio"};

        df_filter_and_weight_to_var1D_list_map[{sign + ""           , "_jacobian_corrected"}] = {"DR", "DR_zoomin"};
        df_filter_and_weight_to_var1D_list_map[{sign + "_wgapcut"   , "_jacobian_corrected"}] = {"DR", "DR_zoomin"};

    }
    
    BuildTrgEffcyFilterToVarListMap();
}


// ---------- ----------
void RDFBasedHistFillingData::BuildTrgEffcyFilterToVarListMap(){

    TrigEffcyFiltersPrePostSumFlattening();

    for (auto filter : trg_effcy_filters_1D_pre_sum)    df_filter_to_var1D_list_map[filter] = single_muon_trig_effcy_var1Ds;
    for (auto filter : trg_effcy_filters_2D_3D_pre_sum) df_filter_to_var2D_list_map[filter] = single_muon_trig_effcy_var2Ds;
    for (auto filter : trg_effcy_filters_2D_3D_pre_sum) df_filter_to_var3D_list_map[filter] = single_muon_trig_effcy_var3Ds;

    // ----- ADD INVERSE WEIGHTING ONES!! WITH _mu4_mu4noL1_denom, _2mu4_denom (paired with additional filtering) -----
}

//--------- BUILD & FLATTERN PRE-SUM, POST-SUM, TO-BE-SUMMED LEVELS ---------
void RDFBasedHistFillingData::TrigEffcyFiltersPrePostSumFlattening()
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
void RDFBasedHistFillingData::HistPostProcess(){
	RDFBasedHistFillingBaseClass::HistPostProcess();

    if (trigger_mode == 0 || trigger_mode == 1){
        if (hist_filling_cycle == generic){
            SumSingleMuonTrigEffHists();
            MakeAndWriteSingleMuonPtTrigEffGraphs();
            CalculateSingleMuonTrigEffcyRatios();           
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
void RDFBasedHistFillingData::MakeAndWriteSingleMuonPtTrigEffGraphsHelper(const std::vector<std::string>& categories){
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

    // ----- draw for each category -----
    for (std::string category : categories){
        TH2D* h_pt2nd_vs_q_eta2nd_mu4_mu4noL1_sepr  = map_at_checked(hist2D_map, "h_pt2nd_vs_q_eta2nd" + category + "_mu4_mu4noL1_sepr",    Form("MakeAndWriteSingleMuonPtTrigEffGraphs: hist2D_map.at(%s)", ("h_pt2nd_vs_q_eta2nd" + category + "_mu4_mu4noL1_sepr").c_str()));
        TH2D* h_pt2nd_vs_q_eta2nd_mu4_sepr          = map_at_checked(hist2D_map, "h_pt2nd_vs_q_eta2nd" + category + "_mu4_sepr",            Form("MakeAndWriteSingleMuonPtTrigEffGraphs: hist2D_map.at(%s)", ("h_pt2nd_vs_q_eta2nd" + category + "_mu4_sepr").c_str()));
        TH2D* h_pt2nd_vs_q_eta2nd_2mu4_sepr         = map_at_checked(hist2D_map, "h_pt2nd_vs_q_eta2nd" + category + "_2mu4_sepr",           Form("MakeAndWriteSingleMuonPtTrigEffGraphs: hist2D_map.at(%s)", ("h_pt2nd_vs_q_eta2nd" + category + "_2mu4_sepr").c_str()));

        if (!h_pt2nd_vs_q_eta2nd_mu4_mu4noL1_sepr || !h_pt2nd_vs_q_eta2nd_mu4_sepr || !h_pt2nd_vs_q_eta2nd_2mu4_sepr) continue;

        proj_make_and_write(h_pt2nd_vs_q_eta2nd_mu4_mu4noL1_sepr,     h_pt2nd_vs_q_eta2nd_mu4_sepr);
        proj_make_and_write(h_pt2nd_vs_q_eta2nd_2mu4_sepr,            h_pt2nd_vs_q_eta2nd_mu4_sepr);

        // ----- extract q-eta bins (might be category dependent) -----
        std::vector<double> eta_bins_trig_effcy = {};

        const TAxis* xaxis = h_pt2nd_vs_q_eta2nd_mu4_sepr->GetXaxis();
        int nbins = xaxis->GetNbins();
        const double* arr = xaxis->GetXbins()->GetArray();
        eta_bins_trig_effcy.assign(arr, arr + nbins + 1);

        for (auto range : q_eta_proj_ranges_for_single_muon_effcy_pT_fitting){
            int bin_first = bin_number(range.first, eta_bins_trig_effcy) + 1; // + 1 pushes into next bin (bin lower end agree with range edge if range edge founded)
            int bin_last = bin_number(range.second, eta_bins_trig_effcy); // bin higher end agree with range edge if range edge founded
            std::string proj_suffix = pairToSuffix(range);

            proj_make_and_write(h_pt2nd_vs_q_eta2nd_mu4_mu4noL1_sepr,     h_pt2nd_vs_q_eta2nd_mu4_sepr, true, bin_first, bin_last, proj_suffix);
            proj_make_and_write(h_pt2nd_vs_q_eta2nd_2mu4_sepr,            h_pt2nd_vs_q_eta2nd_mu4_sepr, true, bin_first, bin_last, proj_suffix);
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
