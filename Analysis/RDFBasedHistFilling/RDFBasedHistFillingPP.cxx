#include <TSystem.h>
#include <TKey.h>
#include <algorithm>
#include "RDFBasedHistFillingData.cxx"

void RDFBasedHistFillingPP::SetIOPathsHook(){
    infile_var1D_json = "var1D_pp.json";

    // Mirror the PbPb SetIOPathsHook ordering (resonance-cut convention, docs/data_analysis.md):
    // NOMINAL/crossx (2mu4, trigger_mode==3) prefers the V1 file (_mindR_0_02, empty res suffix);
    // TRIGGER-EFFICIENCY (single mu4) prefers _res_cut_v2 (V2). V1==V2 after the signal selection;
    // they differ only for generic / low-mass (0–4 GeV) histograms.
    const std::string pp_base =
        "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/muon_pairs_pp_2024" + base_trig_suffix;
    std::string in_path;
    if (low_mass_template_calc) { // low-mass template-fit pass
        // mixed_event_template -> scrambled muon_pairs (ScrambGen mixed-event combinatoric T_mix);
        // else -> _no_res_cut (data D_OS/D_SS with resonances present).
        const std::string tmpl_path = mixed_event_template
            ? pp_base + "_scrambled.root"
            : pp_base + "_no_res_cut.root";
        if (gSystem->AccessPathName(tmpl_path.c_str())) {
            throw std::runtime_error(
                "RDFBasedHistFillingPP: low-mass template-fit input not found: " + tmpl_path +
                " (base_trig_suffix=" + base_trig_suffix + "). " +
                (mixed_event_template ? "mixed_event_template requires the *_scrambled.root (ScrambGen output)."
                                      : "the template-fit pass requires the _no_res_cut ntuples (resonances present)."));
        }
        in_path = tmpl_path;
    } else if (!trigger_effcy_calc) { // pp nominal / crossx (2mu4=mode3, mu4_mu4noL1=mode2) -> V1 ONLY
        // NOMINAL requires V1 _mindR_0_02 (no fallback). The old _res_cut_v2 / _no_res_cut / bare
        // fallback was an OLD-skim crutch (pre-mindR branches) and is removed: a silent fallback to
        // the wrong resonance-cut variant would corrupt crossx/R_AA. (docs/data_analysis.md
        // resonance-cut convention; low_mass_dimuon_template_fit.md Design Decisions 2026-06-23.)
        const std::string v1_path = pp_base + input_mindR_suffix + ".root";
        if (gSystem->AccessPathName(v1_path.c_str())) {
            throw std::runtime_error(
                "RDFBasedHistFillingPP: nominal/crossx V1 input not found: " + v1_path +
                " (base_trig_suffix=" + base_trig_suffix + "). Nominal requires the V1 _mindR_0_02"
                " file; the obsolete _res_cut_v2 / _no_res_cut fallback has been removed.");
        }
        in_path = v1_path;
    } else { // pp trigger-efficiency (single mu4) -> V2 first
        std::vector<std::string> input_candidates = {
            pp_base + input_mindR_suffix + "_res_cut_v2.root",
            pp_base + input_mindR_suffix + ".root",
            pp_base + "_no_res_cut.root",
            pp_base + ".root"
        };
        for (const auto& cand : input_candidates) {
            if (!gSystem->AccessPathName(cand.c_str())) {
                in_path = cand;
                break;
            }
        }
        if (in_path.empty()) {
            throw std::runtime_error("RDFBasedHistFillingPP: input file not found for run_year=24 and base_trig_suffix=" + base_trig_suffix);
        }
    }

    input_files.push_back(in_path);
    output_file = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/histograms_real_pairs_pp_2024" + out_file_suffix + ".root";
}

void RDFBasedHistFillingPP::InitializePPExtra(){
    // Set pp_crossx_lumi_factor from PPBaseClass map.
    // Only trigger_mode 2 (mu4_mu4noL1) and 3 (2mu4) are crossx modes;
    // trigger_mode 1 never calls FillHistogramsCrossx(), so leave factor at -1.
    {
        std::string trig_for_crossx;
        if      (trigger_mode == 2) trig_for_crossx = "mu4_mu4noL1";
        else if (trigger_mode == 3) trig_for_crossx = "2mu4";
        if (!trig_for_crossx.empty()) {
            const auto& m = PPBaseClass::CrossxFactorMap();
            auto it = m.find({run_year, trig_for_crossx});
            if (it != m.end()) pp_crossx_lumi_factor = it->second;
            // else: leave -1; FillHistogramsCrossx() will throw if called
        }
    }

    levels_trg_effcy_to_be_summed_w_musign_summing = {0,1,2}; // ss/op + mu1/2 + mu+/-

    levels_trg_effcy_filters_1D_pre_sum = {{"_ss", "_op"},
                                        {"_mu1passmu4", "_mu2passmu4"},
                                        {"_sign1", "_sign2"},
                                        trigs,
                                        trg_effcy_biases};

    levels_trg_effcy_filters_2D_3D_pre_sum = {{"_ss", "_op"},
                                        {"_mu1passmu4", "_mu2passmu4"},
                                        {"_sign1", "_sign2"},
                                        trigs,
                                        {"", "_sepr"}};

    categories_essential = pair_signs;
}

void RDFBasedHistFillingPP::FlattenTrigEffcyFiltersExtra(){

    // build to-be-summed levels, with mu-sign summing
    for (int level_ind = 0;
         level_ind < static_cast<int>(levels_trg_effcy_filters_2D_3D_pre_sum.size());
         ++level_ind)
    {
        if (std::find(levels_trg_effcy_to_be_summed_w_musign_summing.begin(),
                      levels_trg_effcy_to_be_summed_w_musign_summing.end(),
                      level_ind) != levels_trg_effcy_to_be_summed_w_musign_summing.end())
        {
            levels_trg_effcy_filters_to_be_summed_w_musign_summing.push_back(
                levels_trg_effcy_filters_2D_3D_pre_sum.at(level_ind)
            );
        }
    }

    // flatten to-be-summed levels, with mu-sign summing
    HistFillUtils::flatten_levels(levels_trg_effcy_filters_to_be_summed_w_musign_summing, trg_effcy_filters_to_be_summed_w_musign_summing);    

    // build post-sum levels, with mu-sign summing
    HistFillUtils::write_post_sum_levels(levels_trg_effcy_filters_1D_pre_sum,
                          levels_trg_effcy_to_be_summed_w_musign_summing,
                          levels_trg_effcy_filters_1D_post_sum_w_musign_summing);

    HistFillUtils::write_post_sum_levels(levels_trg_effcy_filters_2D_3D_pre_sum,
                          levels_trg_effcy_to_be_summed_w_musign_summing,
                          levels_trg_effcy_filters_2D_3D_post_sum_w_musign_summing);

    // flatten post-sum levels, with mu-sign summing
    HistFillUtils::flatten_levels(levels_trg_effcy_filters_1D_post_sum_w_musign_summing, trg_effcy_filters_1D_post_sum_w_musign_summing);
    HistFillUtils::flatten_levels(levels_trg_effcy_filters_2D_3D_post_sum_w_musign_summing, trg_effcy_filters_2D_3D_post_sum_w_musign_summing);
}

void RDFBasedHistFillingPP::FillHistogramsSingleMuonEffcy(){
    if (trigger_mode == 1){ // trigger efficiencies: mu4_mu4noL1 | mu4 & 2mu4 | mu4
        FillHistogramsDimuTrigGivenMu4();
    } else if (trigger_mode == 0){
        FillHistogramsMu4GivenMB();
    }
}

void RDFBasedHistFillingPP::FillHistogramsDimuTrigGivenMu4(){
    try{
        for (std::string pair_sign : pair_signs){
            std::string df_name = "df" + pair_sign;
            ROOT::RDF::RNode& node = map_at_checked(df_map, df_name, Form("FillHistogramsSingleMuonEffcy: df_map.at(%s)", df_name.c_str()));

            for (auto mu4sel : {"_mu1passmu4", "_mu2passmu4"}){ // mu4 selection

                std::string df_name = "df" + pair_sign + mu4sel; // e.g, df_ss_mu1passmu4

                std::string ind1st = (std::string(mu4sel) == "_mu1passmu4")? "1" : "2";
                std::string ind2nd = (std::string(mu4sel) == "_mu1passmu4")? "2" : "1";

                df_map.emplace(df_name, node.Filter("m" + ind1st + ".passmu4"));
                df_map.at(df_name) = df_map.at(df_name).Define("pt1st", "m" + ind1st + ".pt");
                df_map.at(df_name) = df_map.at(df_name).Define("pt2nd", "m" + ind2nd + ".pt");
                df_map.at(df_name) = df_map.at(df_name).Define("charge2nd","m" + ind2nd + ".charge");
                df_map.at(df_name) = df_map.at(df_name).Define("eta2nd",    "m" + ind2nd + ".eta");
                df_map.at(df_name) = df_map.at(df_name).Define("phi2nd",    "m" + ind2nd + ".phi");
                df_map.at(df_name) = df_map.at(df_name).Define("q_eta2nd",  "charge2nd * eta2nd");
                df_map.at(df_name) = df_map.at(df_name).Define("mu2nd_passmu4noL1", "m" + ind2nd + ".passmu4noL1");
                df_map.at(df_name) = df_map.at(df_name).Define("mu2nd_good_acceptance",  "pt2nd >= 6 && ((eta2nd > 1.1 && eta2nd < 2.3) || (eta2nd > -2.3 && eta2nd < -1.2))");

                df_map.emplace(df_name + "_sign1", map_at_checked(df_map, df_name, Form("FillHistogramsSingleMuonEffcy: df_map.at(%s)", df_name.c_str())).Filter("charge2nd > 0")); // e.g, df_ss_mu1passmu4_sign1
                df_map.emplace(df_name + "_sign2", map_at_checked(df_map, df_name, Form("FillHistogramsSingleMuonEffcy: df_map.at(%s)", df_name.c_str())).Filter("charge2nd < 0"));

                for (auto mu_sign : {"_sign1", "_sign2"}){
                    std::string df_name = "df" + pair_sign + mu4sel + mu_sign;

                    for (auto trg : {"_mu4", "_mu4_mu4noL1", "_2mu4", "_2mu4_AND_mu4_mu4noL1"}){
                        std::string trg_filter = map_at_checked(trig_to_filter_str_map, trg, Form("trig_to_filter_str_map.at(%s)", trg));
                        if (trg_filter.empty()) df_map.emplace(df_name + trg, map_at_checked(df_map, df_name, Form("FillHistogramsDimuTrigGivenMu4: df_map.at(%s)", df_name.c_str())));
                        else                    df_map.emplace(df_name + trg, map_at_checked(df_map, df_name, Form("FillHistogramsDimuTrigGivenMu4: df_map.at(%s)", df_name.c_str())).Filter(trg_filter));

                        std::string df_name_new = "df" + pair_sign + mu4sel + mu_sign + trg;
                        df_map.emplace(df_name_new + "_sepr",        map_at_checked(df_map, df_name_new, Form("FillHistogramsDimuTrigGivenMu4: df_map.at(%s)", df_name_new.c_str())).Filter("passSeparated"));
                        df_map.emplace(df_name_new + "_good_accept", map_at_checked(df_map, df_name_new, Form("FillHistogramsDimuTrigGivenMu4: df_map.at(%s)", df_name_new.c_str())).Filter("mu2nd_good_acceptance"));

                        for (auto bias : trg_effcy_biases){ // additional selection / bias in data sample
                            std::string df_name_bias = "df" + pair_sign + mu4sel + mu_sign + trg + bias;
                            std::string filter = df_name_bias.substr(2);
                            FillHistogramsSingleDataFrame(filter, map_at_checked(df_map, df_name_bias, Form("FillHistogramsDimuTrigGivenMu4: df_map.at(%s)", df_name_bias.c_str())), false);
                        }
                    }
                }
            } // end loop mu4 selection
        } // end loop pair sign
    } catch(const std::out_of_range& e){
        std::cerr << "FillHistogramsSingleMuonEffcy:: out_of_range exception caught: " << e.what() << std::endl;
    } catch (const std::runtime_error& e) {
        std::cerr << "FillHistogramsSingleMuonEffcy:: RDF runtime error: " << e.what() << std::endl;
    }
}

void RDFBasedHistFillingPP::FillHistogramsMu4GivenMB(){}

void RDFBasedHistFillingPP::FillTrigEffcyHistsInvWeightedbySingleMuonEffcies(){
    OpenEffcyPtFitFile();

    std::vector<std::string> invw_var1Ds = {
        "DR", "DR_zoomin", "DR_0_2", "Deta", "Deta_zoomin",
        "Dphi", "Dphi_zoomin", "pair_pt_log", "minv_zoomin"
    };
    std::vector<std::array<std::string, 2>> invw_var2Ds = {
        {{"DR_zoomin", "pair_pt_log"}},
        {{"DR", "pair_pt_log"}}
    };
    static const std::vector<std::array<std::string, 3>> empty3D;

    try {
        for (std::string pair_sign : {"_ss", "_op"}) {
            std::string df_name = "df" + pair_sign;
            ROOT::RDF::RNode& base_node = map_at_checked(df_map, df_name,
                Form("FillTrigEffcyHistsInvWeighted: df_map.at(%s)", df_name.c_str()));

            auto node = base_node
                .Define("q_eta1", "(float)(m1.charge * m1.eta)")
                .Define("q_eta2", "(float)(m2.charge * m2.eta)")
                .Define("effcy1", [](int q, float pt, float qe) {
                    return RDFBasedHistFillingData::EvaluateSingleMuonEffcy("", q > 0, pt, qe);
                }, {"m1.charge", "m1.pt", "q_eta1"})
                .Define("effcy2", [](int q, float pt, float qe) {
                    return RDFBasedHistFillingData::EvaluateSingleMuonEffcy("", q > 0, pt, qe);
                }, {"m2.charge", "m2.pt", "q_eta2"})
                .Define("valid1", "effcy1 > 0")
                .Define("valid2", "effcy2 > 0")
                .Define("valid_both", "valid1 && valid2")
                .Define("invw_cross", "valid_both ? (double)(1.0f / (effcy1 * effcy2)) : 0.0");

            // §3c: PP 2mu4 — one term only, no tag/probe, no cross-term
            auto df_denom = node.Filter("valid_both");
            FillHistogramsSingleDataFrame(pair_sign + "_2mu4_denom", df_denom, "",
                invw_var1Ds, invw_var2Ds, empty3D, true, {true, true, false});

            auto df_num = df_denom.Filter("pass2mu4");
            FillHistogramsSingleDataFrame(pair_sign + "_2mu4_invw_num", df_num, "invw_cross",
                invw_var1Ds, invw_var2Ds, empty3D, true, {true, true, false});
        }
    } catch (const std::out_of_range& e) {
        std::cerr << "FillTrigEffcyHistsInvWeighted (PP):: out_of_range: " << e.what() << std::endl;
    } catch (const std::runtime_error& e) {
        std::cerr << "FillTrigEffcyHistsInvWeighted (PP):: runtime error: " << e.what() << std::endl;
    }
}

void RDFBasedHistFillingPP::SumSingleMuonTrigEffHistsPP(){

    // 1D, with mu-sign summing
    HistFillUtils::SumTrigEffHistsGeneric<TH1D, std::string>(
        single_muon_trig_effcy_var1Ds,
        trg_effcy_filters_1D_post_sum_w_musign_summing,
        trg_effcy_filters_to_be_summed_w_musign_summing,
        hist1D_map,
        [](const std::string& var) {
            return "h_" + var;
        }
    );

    // 2D, with mu-sign summing
    HistFillUtils::SumTrigEffHistsGeneric<TH2D, std::array<std::string,2>>(
        single_muon_trig_effcy_var2Ds,
        trg_effcy_filters_2D_3D_post_sum_w_musign_summing,
        trg_effcy_filters_to_be_summed_w_musign_summing,
        hist2D_map,
        [](const std::array<std::string,2>& vars) {
            const std::string& varx = vars[0];
            const std::string& vary = vars[1];
            return "h_" + vary + "_vs_" + varx;
        }
    );

    // 3D, with mu-sign summing
    HistFillUtils::SumTrigEffHistsGeneric<TH3D, std::array<std::string,3>>(
        single_muon_trig_effcy_var3Ds,
        trg_effcy_filters_2D_3D_post_sum_w_musign_summing,
        trg_effcy_filters_to_be_summed_w_musign_summing,
        hist3D_map,
        [](const std::array<std::string,3>& vars) {
            const std::string& varx = vars[0];
            const std::string& vary = vars[1];
            const std::string& varz = vars[2];
            return "h_" + varz + "_vs_" + vary + "_vs_" + varx;
        }
    );
}

void RDFBasedHistFillingPP::CalculateSingleMuonTrigEffcyRatios(){
    CalculateSingleMuonTrigEffcyRatiosHelper(musigns);
}

void RDFBasedHistFillingPP::MakeAndWriteSingleMuonTrigEffPtGraphs(){
    if(useCoarseQEtaBin)    MakeAndWriteSingleMuonTrigEffPtGraphsHelper({});
    else                    MakeAndWriteSingleMuonTrigEffPtGraphsHelper(musigns);
}

void RDFBasedHistFillingPP::OpenEffcyPtFitFile() {
    if (!s_effcy_pT_fit_map.empty()) {
        std::cout << "OpenEffcyPtFitFile: TF1 map already loaded (" << s_effcy_pT_fit_map.size() << " entries)" << std::endl;
        return;
    }
    std::string fit_path = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_20" + std::to_string(run_year) + "/trg_effcy_pT_fitting_to_erf_plus_log/single_mu_effcy_pT_fit.root";
    s_effcy_pT_fit_file = TFile::Open(fit_path.c_str(), "READ");
    if (!s_effcy_pT_fit_file || s_effcy_pT_fit_file->IsZombie()) {
        std::cerr << "OpenEffcyPtFitFile: FAILED to open " << fit_path << std::endl;
        return;
    }
    TIter next(s_effcy_pT_fit_file->GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)next())) {
        if (std::string(key->GetClassName()) == "TF1") {
            TF1* func = (TF1*)key->ReadObj();
            s_effcy_pT_fit_map[func->GetName()] = func;
        }
    }
    std::cout << "OpenEffcyPtFitFile: loaded " << s_effcy_pT_fit_map.size() << " TF1s from " << fit_path << std::endl;

    std::string base_dir = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_20" + std::to_string(run_year);
    std::string hist_path = base_dir + "/histograms_real_pairs_pp_20" + std::to_string(run_year) + "_single_mu4_fine_q_eta_bin.root";
    s_effcy_2D_hist_file = TFile::Open(hist_path.c_str(), "READ");
    if (!s_effcy_2D_hist_file || s_effcy_2D_hist_file->IsZombie()) {
        std::cerr << "OpenEffcyPtFitFile: WARNING — 2D hist file not found: " << hist_path << " (gap fallback disabled)" << std::endl;
        return;
    }
    TIter next_h2d(s_effcy_2D_hist_file->GetListOfKeys());
    TKey* key2;
    int n2d = 0;
    while ((key2 = (TKey*)next_h2d())) {
        std::string name = key2->GetName();
        if (std::string(key2->GetClassName()) == "TH2D" && name.find("_divided") != std::string::npos
            && name.find("pt2nd_vs_q_eta2nd") != std::string::npos) {
            TH2D* h = (TH2D*)key2->ReadObj();
            s_effcy_2D_hist_map[h->GetName()] = h;
            n2d++;
        }
    }
    std::cout << "OpenEffcyPtFitFile: loaded " << n2d << " 2D efficiency histograms from " << hist_path << std::endl;
}

void RDFBasedHistFillingPP::MakeAndWriteDRTrigEffGraphs() {
    MakeAndWriteDRTrigEffGraphsHelper({});
}

// ============================================================================
// FillHistogramsGeneric: override to add 2mu4 trigger efficiency weighting
// ============================================================================

void RDFBasedHistFillingPP::FillHistogramsGeneric(){
    if (!trigger_effcy_calc) {
        OpenEffcyPtFitFile();
        OpenRecoEffPlaceholderFile();  // reco-eff PLACEHOLDER (so generic analysis histos are reco-corrected too)

        for (const std::string& category : categories_essential) {
            std::string df_name = "df" + category;
            ROOT::RDF::RNode& df = map_at_checked(df_map, df_name,
                Form("PP::FillHistogramsGeneric: df_map.at(%s)", df_name.c_str()));

            // Generic analysis histograms (incl. the gapcut histos read by the
            // MC-data comparison) are weighted by the SAME per-pair efficiency
            // correction as the crossx: w_reco_trig = w_reco * w_trig. This keeps
            // the MC-data comparison consistent with the reco-eff placeholder
            // (and any future efficiency/det-response/unfolding change to pp crossx).
            // w_reco is the eps1*eps2 reco PLACEHOLDER (pp -> centrality = -1).
            ROOT::RDF::RNode df_with_trig = df
                .Define("q_eta1", "(float)(m1.charge * m1.eta)")
                .Define("q_eta2", "(float)(m2.charge * m2.eta)")
                .Define("effcy1", [](int q, float pt, float qe) {
                    return RDFBasedHistFillingData::EvaluateSingleMuonEffcy("", q > 0, pt, qe);
                }, {"m1.charge", "m1.pt", "q_eta1"})
                .Define("effcy2", [](int q, float pt, float qe) {
                    return RDFBasedHistFillingData::EvaluateSingleMuonEffcy("", q > 0, pt, qe);
                }, {"m2.charge", "m2.pt", "q_eta2"})
                .Define("effcy_pair", "effcy1 > 0 && effcy2 > 0 ? (double)(effcy1 * effcy2) : -1.0")
                .Define("w_trig", "effcy_pair > 0 ? 1.0 / effcy_pair : 0.0")
                .Define("effcy_reco1", [](float pt, float qe) {
                    return RDFBasedHistFillingData::EvaluateSingleMuonRecoEffPlaceholder(-1, pt, qe);
                }, {"m1.pt", "q_eta1"})
                .Define("effcy_reco2", [](float pt, float qe) {
                    return RDFBasedHistFillingData::EvaluateSingleMuonRecoEffPlaceholder(-1, pt, qe);
                }, {"m2.pt", "q_eta2"})
                .Define("effcy_reco_pair", "effcy_reco1 > 0 && effcy_reco2 > 0 ? (double)(effcy_reco1 * effcy_reco2) : -1.0")
                .Define("w_reco", "effcy_reco_pair > 0 ? 1.0 / (effcy_reco_pair < 0.05 ? 0.05 : effcy_reco_pair) : 1.0")
                .Define("w_reco_trig", "w_reco * w_trig");

            df_map.erase(df_name);
            df_map.emplace(df_name, df_with_trig);
        }
        generic_weight_col = "w_reco_trig";  // reco+trig corrected (was "w_trig")
        std::cout << "[PP] FillHistogramsGeneric: w_reco_trig (reco+trig) columns added to "
                  << categories_essential.size() << " dataframes" << std::endl;
    }

    RDFBasedHistFillingData::FillHistogramsGeneric();
}

// ============================================================================
// FillHistogramsCrossx: Single B → J/psi crossx measurements (trigger_mode==2)
// ============================================================================

void RDFBasedHistFillingPP::FillHistogramsCrossx(){
    if (pp_crossx_lumi_factor < 0.) {
        throw std::runtime_error(
            "[PP] FillHistogramsCrossx: pp_crossx_lumi_factor not set (value " +
            std::to_string(pp_crossx_lumi_factor) + "). "
            "No entry in PPBaseClass::CrossxFactorMap for (run_year=" +
            std::to_string(run_year) + ", trigger_mode=" + std::to_string(trigger_mode) + "). "
            "Valid crossx trigger_modes: 2 (mu4_mu4noL1), 3 (2mu4)."
        );
    }
    OpenEffcyPtFitFile();
    OpenRecoEffPlaceholderFile();  // Run 2 reco-eff PLACEHOLDER (eps1*eps2 proxy)

    std::cout << "[PP] FillHistogramsCrossx: opposite-sign only, signal cuts, "
              << "crossx_weight = weight * " << pp_crossx_lumi_factor
              << " (1/L_int), with 2mu4 trig eff correction" << std::endl;

    const std::string signal_cuts = "minv > 1.08 && minv < 2.9 && pair_pt > 8 && m1.charge * m1.eta < 2.2 && m2.charge * m2.eta < 2.2";

    ROOT::RDF::RNode df_op_base = map_at_checked(df_map, "df_op", "FillHistogramsCrossx: df_op");
    ROOT::RDF::RNode df_single_b_crossx = df_op_base.Filter(signal_cuts);
    if (df_map.find("df_single_b_crossx") == df_map.end()) {
        df_map.emplace("df_single_b_crossx", df_single_b_crossx);
    }

    // Per-pair 2mu4 trigger efficiency: ε_pair = ε₁ · ε₂ (AND logic)
    // If FillHistogramsGeneric already added trigger columns to df_op, they
    // propagate through the Filter to df_single_b_crossx — skip re-Define.
    // Per-pair trigger (w_trig) AND reco-eff PLACEHOLDER (w_reco = 1/(eps1*eps2),
    // pp -> centrality=-1, pair eff floored at 0.05, lookup-fail -> 1) columns.
    // If FillHistogramsGeneric already added them to df_op (generic ran; it now
    // defines the full trig+reco chain), they propagate through the signal-cut
    // Filter -> reuse them (skip re-Define to avoid an RDF column collision).
    ROOT::RDF::RNode df_with_trig = generic_weight_col.empty()
        ? df_single_b_crossx
            .Define("q_eta1", "(float)(m1.charge * m1.eta)")
            .Define("q_eta2", "(float)(m2.charge * m2.eta)")
            .Define("effcy1", [](int q, float pt, float qe) {
                return RDFBasedHistFillingData::EvaluateSingleMuonEffcy("", q > 0, pt, qe);
            }, {"m1.charge", "m1.pt", "q_eta1"})
            .Define("effcy2", [](int q, float pt, float qe) {
                return RDFBasedHistFillingData::EvaluateSingleMuonEffcy("", q > 0, pt, qe);
            }, {"m2.charge", "m2.pt", "q_eta2"})
            .Define("effcy_pair", "effcy1 > 0 && effcy2 > 0 ? (double)(effcy1 * effcy2) : -1.0")
            .Define("w_trig", "effcy_pair > 0 ? 1.0 / effcy_pair : 0.0")
            .Define("effcy_reco1", [](float pt, float qe) {
                return RDFBasedHistFillingData::EvaluateSingleMuonRecoEffPlaceholder(-1, pt, qe);
            }, {"m1.pt", "q_eta1"})
            .Define("effcy_reco2", [](float pt, float qe) {
                return RDFBasedHistFillingData::EvaluateSingleMuonRecoEffPlaceholder(-1, pt, qe);
            }, {"m2.pt", "q_eta2"})
            .Define("effcy_reco_pair", "effcy_reco1 > 0 && effcy_reco2 > 0 ? (double)(effcy_reco1 * effcy_reco2) : -1.0")
            .Define("w_reco", "effcy_reco_pair > 0 ? 1.0 / (effcy_reco_pair < 0.05 ? 0.05 : effcy_reco_pair) : 1.0")
        : df_single_b_crossx;

    const double lumi_factor = pp_crossx_lumi_factor;
    ROOT::RDF::RNode df_single_b_crossx_weighted =
        df_with_trig
        // PLACEHOLDER: unfolding identity until det-response unfolding lands (roadmap Q4).
        // NOTE: real unfolding is NOT a per-pair weight — it is a SPECTRUM-LEVEL
        // operation (response matrix / iterative Bayes on the histogram). Implementing
        // it requires a STRUCTURAL change (unfold the trigger-corrected reco spectrum,
        // THEN apply reco-eff binned in TRUTH kinematics), not just replacing this 1.0.
        .Define("w_unfold", "1.0")
        // Base crossx weight (lumi-scaled, no efficiency) and the sequential
        // correction-stage weights (see CorrectionStages.h).
        .Define("crossx_weight",
            [lumi_factor](double weight){ return weight * lumi_factor; },
            {"weight"})
        // NOMINAL corrected weight now INCLUDES the reco-eff PLACEHOLDER (w_reco):
        // every crossx histogram filled with crossx_weight_trig_corr is reco+trig
        // corrected. Equals the _corr_unfolded_reco_trig stage (w_unfold==1). To
        // revert to trig-only, drop "w_reco *". (reco_eff_placeholder_run2.md)
        .Define("crossx_weight_trig_corr", "crossx_weight * w_reco * w_trig")
        .Define("cw_raw",                "crossx_weight")
        .Define("cw_unfolded",           "cw_raw * w_unfold")
        .Define("cw_unfolded_reco",      "cw_unfolded * w_reco")
        .Define("cw_unfolded_reco_trig", "cw_unfolded_reco * w_trig");

    if (df_map.find("df_single_b_crossx_weighted") == df_map.end()) {
        df_map.emplace("df_single_b_crossx_weighted", df_single_b_crossx_weighted);
    }

    // Store RResultPtrs (lazy-evaluated) instead of immediate clones.
    // They will be evaluated during HistPostProcess() and converted to raw pointers.
    const int    npt   = (int)(pms.pT_bins_120.size() - 1);
    const double* ptbins = pms.pT_bins_120.data();

    // --- LOW-MASS TEMPLATE-FIT pass (low_mass_template_calc) ---
    // Reads _no_res_cut (resonances PRESENT; selected in SetIOPathsHook), distinct output.
    // Produces ONLY the 0-4 GeV minv template spectra D_OS/D_SS (1D + 2D vs pair pT/eta)
    // for the low-mass dimuon template fit (docs/tracking/low_mass_dimuon_template_fit.md
    // 3a). Returns early so the signal-region crossx (minv in [1.08,2.9], WRONG from
    // _no_res_cut because of resonance leakage) is NOT produced. Same selection (signal_cuts
    // MINUS minv window, NO dR) and dsigma weight as the nominal h1d_crossx_minv_0_4 block.
    if (low_mass_template_calc) {
        const std::string signal_cuts_no_minv =
            "pair_pt > 8 && m1.charge * m1.eta < 2.2 && m2.charge * m2.eta < 2.2";
        const double lumi_factor_tmpl = pp_crossx_lumi_factor;
        auto attach_crossx_weight = [&](ROOT::RDF::RNode node) -> ROOT::RDF::RNode {
            ROOT::RDF::RNode n = generic_weight_col.empty()
                ? node
                    .Define("q_eta1", "(float)(m1.charge * m1.eta)")
                    .Define("q_eta2", "(float)(m2.charge * m2.eta)")
                    .Define("effcy1", [](int q, float pt, float qe){ return RDFBasedHistFillingData::EvaluateSingleMuonEffcy("", q > 0, pt, qe); }, {"m1.charge","m1.pt","q_eta1"})
                    .Define("effcy2", [](int q, float pt, float qe){ return RDFBasedHistFillingData::EvaluateSingleMuonEffcy("", q > 0, pt, qe); }, {"m2.charge","m2.pt","q_eta2"})
                    .Define("effcy_pair", "effcy1 > 0 && effcy2 > 0 ? (double)(effcy1 * effcy2) : -1.0")
                    .Define("w_trig", "effcy_pair > 0 ? 1.0 / effcy_pair : 0.0")
                    .Define("effcy_reco1", [](float pt, float qe){ return RDFBasedHistFillingData::EvaluateSingleMuonRecoEffPlaceholder(-1, pt, qe); }, {"m1.pt","q_eta1"})
                    .Define("effcy_reco2", [](float pt, float qe){ return RDFBasedHistFillingData::EvaluateSingleMuonRecoEffPlaceholder(-1, pt, qe); }, {"m2.pt","q_eta2"})
                    .Define("effcy_reco_pair", "effcy_reco1 > 0 && effcy_reco2 > 0 ? (double)(effcy_reco1 * effcy_reco2) : -1.0")
                    .Define("w_reco", "effcy_reco_pair > 0 ? 1.0 / (effcy_reco_pair < 0.05 ? 0.05 : effcy_reco_pair) : 1.0")
                : node;
            return n
                .Define("crossx_weight", [lumi_factor_tmpl](double weight){ return weight * lumi_factor_tmpl; }, {"weight"})
                .Define("crossx_weight_trig_corr", "crossx_weight * w_reco * w_trig");
        };
        ROOT::RDF::RNode df_op_t = attach_crossx_weight(map_at_checked(df_map, "df_op", "FillHistogramsCrossx PP: df_op (template)").Filter(signal_cuts_no_minv));
        ROOT::RDF::RNode df_ss_t = attach_crossx_weight(map_at_checked(df_map, "df_ss", "FillHistogramsCrossx PP: df_ss (template)").Filter(signal_cuts_no_minv));
        const int npt150 = (int)(pms.pT_bins_150.size() - 1);
        const double* ptb150 = pms.pT_bins_150.data();
        hist1d_rresultptr_map["h1d_crossx_minv_0_4_op_dsigma"] = df_op_t.Histo1D(ROOT::RDF::TH1DModel("h1d_crossx_minv_0_4_op_dsigma", ";m_{#mu#mu} [GeV];d#sigma/dm_{#mu#mu} [nb GeV^{-1}]", 50, 0.0, 4.0), "minv", "crossx_weight_trig_corr");
        hist1d_rresultptr_map["h1d_crossx_minv_0_4_ss_dsigma"] = df_ss_t.Histo1D(ROOT::RDF::TH1DModel("h1d_crossx_minv_0_4_ss_dsigma", ";m_{#mu#mu} [GeV];d#sigma/dm_{#mu#mu} [nb GeV^{-1}]", 50, 0.0, 4.0), "minv", "crossx_weight_trig_corr");
        hist2d_rresultptr_map["h2d_crossx_minv_0_4_vs_pair_pt_log_150_op_dsigma"] = df_op_t.Histo2D(ROOT::RDF::TH2DModel("h2d_crossx_minv_0_4_vs_pair_pt_log_150_op_dsigma", ";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]", npt150, ptb150, 50, 0.0, 4.0), "pair_pt", "minv", "crossx_weight_trig_corr");
        hist2d_rresultptr_map["h2d_crossx_minv_0_4_vs_pair_pt_log_150_ss_dsigma"] = df_ss_t.Histo2D(ROOT::RDF::TH2DModel("h2d_crossx_minv_0_4_vs_pair_pt_log_150_ss_dsigma", ";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]", npt150, ptb150, 50, 0.0, 4.0), "pair_pt", "minv", "crossx_weight_trig_corr");
        hist2d_rresultptr_map["h2d_crossx_minv_0_4_vs_pair_eta_op_dsigma"] = df_op_t.Histo2D(ROOT::RDF::TH2DModel("h2d_crossx_minv_0_4_vs_pair_eta_op_dsigma", ";#eta^{pair};m_{#mu#mu} [GeV]", 24, -2.4, 2.4, 50, 0.0, 4.0), "pair_eta", "minv", "crossx_weight_trig_corr");
        hist2d_rresultptr_map["h2d_crossx_minv_0_4_vs_pair_eta_ss_dsigma"] = df_ss_t.Histo2D(ROOT::RDF::TH2DModel("h2d_crossx_minv_0_4_vs_pair_eta_ss_dsigma", ";#eta^{pair};m_{#mu#mu} [GeV]", 24, -2.4, 2.4, 50, 0.0, 4.0), "pair_eta", "minv", "crossx_weight_trig_corr");
        std::cout << "[PP] FillHistogramsCrossx (low-mass template mode, "
                  << (mixed_event_template ? "_scrambled/mixed-event" : "_no_res_cut") << ") completed" << std::endl;
        return;
    }

    // ROOT 6.34 TH3DModel has no mixed (variable/uniform) ctor; generate edge arrays for uniform axes.
    auto make_unif_edges = [](int n, double lo, double hi) {
        std::vector<double> e(n + 1);
        const double step = (hi - lo) / n;
        for (int i = 0; i <= n; ++i) e[i] = lo + i * step;
        return e;
    };
    const auto eta_edges  = make_unif_edges(44, -2.4,  2.4);
    const auto minv_edges = make_unif_edges(50,  1.0,  3.0);
    const auto dr_edges   = make_unif_edges(50,  0.0, 1.0);

    hist2d_rresultptr_map["h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts"] = df_single_b_crossx_weighted.Histo2D(
        ROOT::RDF::TH2DModel("h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts", ";p_{T}^{pair} [GeV];#eta^{pair}", npt, ptbins, 44, -2.4, 2.4),
        "pair_pt", "pair_eta", "crossx_weight_trig_corr");
    hist2d_rresultptr_map["h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts_no_trig_corr"] = df_single_b_crossx_weighted.Histo2D(
        ROOT::RDF::TH2DModel("h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts_no_trig_corr", ";p_{T}^{pair} [GeV];#eta^{pair}", npt, ptbins, 44, -2.4, 2.4),
        "pair_pt", "pair_eta", "crossx_weight");

    // Correction-stage histograms (raw -> unfolded -> +reco -> +reco+trig) for the
    // primary pair_pt x pair_eta differential, so each correction's impact is
    // visible at plotting time. See CorrectionStages.h.
    for (const auto& st : CrossxCorrectionStages()) {
        const std::string nm = std::string("h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts") + st.suffix;
        hist2d_rresultptr_map[nm] = df_single_b_crossx_weighted.Histo2D(
            ROOT::RDF::TH2DModel(nm.c_str(), ";p_{T}^{pair} [GeV];#eta^{pair}", npt, ptbins, 44, -2.4, 2.4),
            "pair_pt", "pair_eta", st.weight_col);
    }

    // --- Same-sign (SS) signal-region 2D yield, for the OS-SS combinatorial
    // subtraction in R_AA (analysis_overview.md §4a). Mirrors the OS path with
    // identical cuts and corrected weight (crossx_weight · w_reco · w_trig).
    // Trigger columns come from FillHistogramsGeneric (added to df_ss too, since
    // categories_essential = pair_signs) when generic ran; else define inline.
    {
        ROOT::RDF::RNode df_ss_crossx = map_at_checked(df_map, "df_ss", "FillHistogramsCrossx PP: df_ss").Filter(signal_cuts);
        // Full trig+reco efficiency chain — defined inline only if generic did NOT
        // run (else reuse the columns FillHistogramsGeneric added to df_ss).
        ROOT::RDF::RNode df_ss_with_trig = generic_weight_col.empty()
            ? df_ss_crossx
                .Define("q_eta1", "(float)(m1.charge * m1.eta)")
                .Define("q_eta2", "(float)(m2.charge * m2.eta)")
                .Define("effcy1", [](int q, float pt, float qe) {
                    return RDFBasedHistFillingData::EvaluateSingleMuonEffcy("", q > 0, pt, qe);
                }, {"m1.charge", "m1.pt", "q_eta1"})
                .Define("effcy2", [](int q, float pt, float qe) {
                    return RDFBasedHistFillingData::EvaluateSingleMuonEffcy("", q > 0, pt, qe);
                }, {"m2.charge", "m2.pt", "q_eta2"})
                .Define("effcy_pair", "effcy1 > 0 && effcy2 > 0 ? (double)(effcy1 * effcy2) : -1.0")
                .Define("w_trig", "effcy_pair > 0 ? 1.0 / effcy_pair : 0.0")
                .Define("effcy_reco1", [](float pt, float qe) {
                    return RDFBasedHistFillingData::EvaluateSingleMuonRecoEffPlaceholder(-1, pt, qe);
                }, {"m1.pt", "q_eta1"})
                .Define("effcy_reco2", [](float pt, float qe) {
                    return RDFBasedHistFillingData::EvaluateSingleMuonRecoEffPlaceholder(-1, pt, qe);
                }, {"m2.pt", "q_eta2"})
                .Define("effcy_reco_pair", "effcy_reco1 > 0 && effcy_reco2 > 0 ? (double)(effcy_reco1 * effcy_reco2) : -1.0")
                .Define("w_reco", "effcy_reco_pair > 0 ? 1.0 / (effcy_reco_pair < 0.05 ? 0.05 : effcy_reco_pair) : 1.0")
            : df_ss_crossx;
        ROOT::RDF::RNode df_ss_weighted = df_ss_with_trig
            .Define("crossx_weight",
                [lumi_factor](double weight){ return weight * lumi_factor; }, {"weight"})
            .Define("crossx_weight_trig_corr", "crossx_weight * w_reco * w_trig");
        hist2d_rresultptr_map["h2d_ss_crossx_pair_pt_pair_eta_binned_w_signal_cuts"] = df_ss_weighted.Histo2D(
            ROOT::RDF::TH2DModel("h2d_ss_crossx_pair_pt_pair_eta_binned_w_signal_cuts", ";p_{T}^{pair} [GeV];#eta^{pair}", npt, ptbins, 44, -2.4, 2.4),
            "pair_pt", "pair_eta", "crossx_weight_trig_corr");
    }

    hist2d_rresultptr_map["h2d_crossx_pair_pt_minv_w_signal_cuts"] = df_single_b_crossx_weighted.Histo2D(
        ROOT::RDF::TH2DModel("h2d_crossx_pair_pt_minv_w_signal_cuts", ";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]", npt, ptbins, 50, 1.0, 3.0),
        "pair_pt", "minv", "crossx_weight_trig_corr");
    hist2d_rresultptr_map["h2d_crossx_pair_pt_dr_w_signal_cuts"] = df_single_b_crossx_weighted.Histo2D(
        ROOT::RDF::TH2DModel("h2d_crossx_pair_pt_dr_w_signal_cuts", ";p_{T}^{pair} [GeV];#DeltaR", npt, ptbins, 50, 0.0, 1.0),
        "pair_pt", "dr", "crossx_weight_trig_corr");

    hist3d_rresultptr_map["h3d_crossx_minv_vs_pair_eta_vs_pair_pt_w_signal_cuts"] = df_single_b_crossx_weighted.Histo3D(
        ROOT::RDF::TH3DModel("h3d_crossx_minv_vs_pair_eta_vs_pair_pt_w_signal_cuts", ";p_{T}^{pair} [GeV];#eta^{pair};m_{#mu#mu} [GeV]", npt, ptbins, 44, eta_edges.data(), 50, minv_edges.data()),
        "pair_pt", "pair_eta", "minv", "crossx_weight_trig_corr");
    hist3d_rresultptr_map["h3d_crossx_dr_vs_pair_eta_vs_pair_pt_w_signal_cuts"] = df_single_b_crossx_weighted.Histo3D(
        ROOT::RDF::TH3DModel("h3d_crossx_dr_vs_pair_eta_vs_pair_pt_w_signal_cuts", ";p_{T}^{pair} [GeV];#eta^{pair};#DeltaR", npt, ptbins, 44, eta_edges.data(), 50, dr_edges.data()),
        "pair_pt", "pair_eta", "dr", "crossx_weight_trig_corr");

    // --- pT_bins_150 variants ---
    {
        const int    npt150    = (int)(pms.pT_bins_150.size() - 1);
        const double* ptbins150 = pms.pT_bins_150.data();
        const auto eta_edges150 = make_unif_edges(44, -2.4, 2.4);
        const auto dr_edges150  = make_unif_edges(50, 0.0, 1.0);

        hist2d_rresultptr_map["h2d_crossx_pt_150_pair_eta_binned_w_signal_cuts"] = df_single_b_crossx_weighted.Histo2D(
            ROOT::RDF::TH2DModel("h2d_crossx_pt_150_pair_eta_binned_w_signal_cuts", ";p_{T}^{pair} [GeV];#eta^{pair}", npt150, ptbins150, 44, -2.4, 2.4),
            "pair_pt", "pair_eta", "crossx_weight_trig_corr");
        hist3d_rresultptr_map["h3d_crossx_dr_vs_pair_eta_vs_pt_150_w_signal_cuts"] = df_single_b_crossx_weighted.Histo3D(
            ROOT::RDF::TH3DModel("h3d_crossx_dr_vs_pair_eta_vs_pt_150_w_signal_cuts", ";p_{T}^{pair} [GeV];#eta^{pair};#DeltaR", npt150, ptbins150, 44, eta_edges150.data(), 50, dr_edges150.data()),
            "pair_pt", "pair_eta", "dr", "crossx_weight_trig_corr");
    }

    // --- Low-mass dimuon template-fit inputs (docs/tracking/low_mass_dimuon_template_fit.md,
    // Physics Procedure §2,§3a,§4): OS and SS minv over 0–4 GeV with the dimuon
    // mass window REMOVED (all other single-b cuts kept). Crossx-normalized (1/L_pp)
    // and reco+trig corrected (crossx_weight_trig_corr), IDENTICAL selection/binning/
    // weight for OS and SS so D_OS - D_SS is a clean combinatoric subtraction.
    {
        // signal_cuts MINUS the two minv-window conditions; everything else identical.
        const std::string signal_cuts_no_minv =
            "pair_pt > 8 && m1.charge * m1.eta < 2.2 && m2.charge * m2.eta < 2.2";

        // Attach the trig+reco+crossx weight chain to a node filtered WITHOUT the
        // minv window. If FillHistogramsGeneric ran, it added w_trig/w_reco to
        // df_op/df_ss before any filter, so they propagate through this filter —
        // reuse them; else define inline exactly as the OS/SS crossx blocks above do.
        auto attach_crossx_weight = [&](ROOT::RDF::RNode node) -> ROOT::RDF::RNode {
            ROOT::RDF::RNode n = generic_weight_col.empty()
                ? node
                    .Define("q_eta1", "(float)(m1.charge * m1.eta)")
                    .Define("q_eta2", "(float)(m2.charge * m2.eta)")
                    .Define("effcy1", [](int q, float pt, float qe) {
                        return RDFBasedHistFillingData::EvaluateSingleMuonEffcy("", q > 0, pt, qe);
                    }, {"m1.charge", "m1.pt", "q_eta1"})
                    .Define("effcy2", [](int q, float pt, float qe) {
                        return RDFBasedHistFillingData::EvaluateSingleMuonEffcy("", q > 0, pt, qe);
                    }, {"m2.charge", "m2.pt", "q_eta2"})
                    .Define("effcy_pair", "effcy1 > 0 && effcy2 > 0 ? (double)(effcy1 * effcy2) : -1.0")
                    .Define("w_trig", "effcy_pair > 0 ? 1.0 / effcy_pair : 0.0")
                    .Define("effcy_reco1", [](float pt, float qe) {
                        return RDFBasedHistFillingData::EvaluateSingleMuonRecoEffPlaceholder(-1, pt, qe);
                    }, {"m1.pt", "q_eta1"})
                    .Define("effcy_reco2", [](float pt, float qe) {
                        return RDFBasedHistFillingData::EvaluateSingleMuonRecoEffPlaceholder(-1, pt, qe);
                    }, {"m2.pt", "q_eta2"})
                    .Define("effcy_reco_pair", "effcy_reco1 > 0 && effcy_reco2 > 0 ? (double)(effcy_reco1 * effcy_reco2) : -1.0")
                    .Define("w_reco", "effcy_reco_pair > 0 ? 1.0 / (effcy_reco_pair < 0.05 ? 0.05 : effcy_reco_pair) : 1.0")
                : node;
            return n
                .Define("crossx_weight",
                    [lumi_factor](double weight){ return weight * lumi_factor; }, {"weight"})
                .Define("crossx_weight_trig_corr", "crossx_weight * w_reco * w_trig");
        };

        ROOT::RDF::RNode df_op_no_minv = attach_crossx_weight(
            map_at_checked(df_map, "df_op", "FillHistogramsCrossx PP: df_op (no-minv)").Filter(signal_cuts_no_minv));
        ROOT::RDF::RNode df_ss_no_minv = attach_crossx_weight(
            map_at_checked(df_map, "df_ss", "FillHistogramsCrossx PP: df_ss (no-minv)").Filter(signal_cuts_no_minv));

        hist1d_rresultptr_map["h1d_crossx_minv_0_4_op_dsigma"] = df_op_no_minv.Histo1D(
            ROOT::RDF::TH1DModel("h1d_crossx_minv_0_4_op_dsigma", ";m_{#mu#mu} [GeV];d#sigma/dm_{#mu#mu} [nb GeV^{-1}]", 50, 0.0, 4.0),
            "minv", "crossx_weight_trig_corr");
        hist1d_rresultptr_map["h1d_crossx_minv_0_4_ss_dsigma"] = df_ss_no_minv.Histo1D(
            ROOT::RDF::TH1DModel("h1d_crossx_minv_0_4_ss_dsigma", ";m_{#mu#mu} [GeV];d#sigma/dm_{#mu#mu} [nb GeV^{-1}]", 50, 0.0, 4.0),
            "minv", "crossx_weight_trig_corr");
    }

    std::cout << "[PP] FillHistogramsCrossx completed" << std::endl;
}
