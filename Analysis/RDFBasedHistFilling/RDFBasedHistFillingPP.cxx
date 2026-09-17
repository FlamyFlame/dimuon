#include <TSystem.h>
#include <TKey.h>
#include <algorithm>
#include "RDFBasedHistFillingData.cxx"
#include "../Utilities/PairTrigEffCrossxEvaluator.h"
#include "../Utilities/PairRecoEffEvaluator.h"

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

            // WP config var: nominal = TIGHT (isTight default true). Measure the single-muon trigger
            // turn-on on TIGHT tag+probe pairs so the efficiency matches the tight signal selection
            // (crossx selects pair_pass_tight); Medium (systematic, isTight=false) uses all medium pairs.
            if (isTight) { df_map.at(df_name) = df_map.at(df_name).Filter("pair_pass_tight", "tight_tag_and_probe"); }

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
                // PROBE-side fiducial gap cut (round 8, user decision): the data single-muon
                // efficiency is measured on the SAME fiducial region the analysis applies it to.
                // PROBE ONLY -- the tag is deliberately left uncut: eps^nc is a per-muon
                // efficiency, so the tag's q*eta does not enter its definition, and cutting the
                // tag too would only cost statistics. Windows come from ParamsSet, never retyped.
                if (apply_fiducial_gap_cut)
                    df_map.at(df_name) = df_map.at(df_name)
                        .Filter(ParamsSet::FiducialGapCutExpr("q_eta2nd"), "probe fiducial gap cut");
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
    // The charge-integrated form ({}) exists for the obsolete isForSoumya output. The trigger
    // efficiency itself is ALWAYS per charge -- the toroid bends mu+ and mu- oppositely, which is
    // why the axis is q*eta at all, and the lookup key carries _sign1/_sign2. Before round 8 this
    // branched on useCoarseQEtaBin, so making the coarse binning nominal would have silently
    // produced NO per-charge graphs and left every efficiency lookup missing.
    if (isForSoumya) MakeAndWriteSingleMuonTrigEffPtGraphsHelper({});
    else             MakeAndWriteSingleMuonTrigEffPtGraphsHelper(musigns);
}

// =================================================================================================
// The PAIR corrections the pp24 cross-section applies (2026-09-17;
// docs/tracking/pp24_trig_eff_hybrid_application.md Physics Procedure §2):
//
//   eps_trig^pair  -- the per-pair 2mu4 trigger efficiency, a HYBRID over the canonical coarse
//     pair-pT axis (ParamsSet::pair_pt_coarse_bins): in bins 1..N-2 ([9, 74.24) GeV) the DATA
//     single-muon turn-ons times the MC dR correlation correction eps_dR (opposite sign, 3-group
//     |eta^pair| fold, no a-priori plateau; expo primary, polynomial primary in three user-named
//     forward cells, interpolation fallback, NO raw-bin tier); in the LAST TWO bins the
//     SINGLE-VALUE MC pair efficiency in the signal mass window times the product of the two
//     single-muon data/MC scale factors. Utilities/PairTrigEffCrossxEvaluator.h -- and its four
//     stated, user-decided limitations (D1 refused-cell fallback TEMPORARY, D2 same-sign pairs
//     on the opposite-sign numbers TEMPORARY, D3 signal-window number for every pair, D4 no raw
//     tier). Configuration named ONCE in dr_correction_sample_cfg.h.
//
//   eps_reco^pair(pair pT, pair eta, dR)  -- the pp24-fullsim PAIR reconstruction efficiency,
//     REPLACING the Run 2 single-muon eps_1*eps_2 placeholder (which had no dR dependence at all).
//     Utilities/PairRecoEffEvaluator.h.
//
// Both live on the heap for the lifetime of the process: RDF Defines are LAZY, so a stack-scoped
// evaluator would be destroyed before the event loop runs.
// =================================================================================================
static PairTrigEffCrossxEvaluator* s_pair_trig_eff = nullptr;
static PairRecoEffEvaluator*       s_pair_reco_eff = nullptr;

static int s_pair_eff_loaded_wp = -1;   // -1 = nothing loaded; else the isTight the maps carry

void RDFBasedHistFillingPP::OpenPairEfficiencyInputs()
{
    // Both corrections come from the SAME pp24 fullsim FULL production, so the directory is taken
    // from the one place that names it (dr_correction_sample_cfg.h) instead of being retyped.
    const DrCorrSample mc = GetDrCorrSample("pp_full", isTight);

    // The maps are WP-SPECIFIC (Tight and Medium are different files). Caching them on "already
    // loaded" alone would let a second Run() in the same process, at the other working point,
    // silently reuse the first WP's efficiencies -- every histogram would still fill.
    if (s_pair_eff_loaded_wp >= 0 && s_pair_eff_loaded_wp != (isTight ? 1 : 0))
        throw std::runtime_error("OpenPairEfficiencyInputs: the pair-efficiency maps are already "
                                 "loaded for the OTHER working point. Run one WP per process.");

    if (!s_pair_trig_eff) {
        s_pair_trig_eff = new PairTrigEffCrossxEvaluator();
        s_pair_trig_eff->Load(mc, isTight);
    }
    if (!s_pair_reco_eff) {
        s_pair_reco_eff = new PairRecoEffEvaluator();
        s_pair_reco_eff->Load(DrCorrPairRecoEffFile(mc), isTight);
    }
    s_pair_eff_loaded_wp = isTight ? 1 : 0;
}

// The per-pair census of both MC corrections, AFTER the event loop. Neither header hides anything
// -- how many pairs took a fallback, were clamped past the measured domain, or fell outside the
// eps_dR cell grid (which must be zero) -- but none of it reaches the log unless this is called.
void RDFBasedHistFillingPP::PrintPairEfficiencyStats()
{
    if (s_pair_trig_eff) s_pair_trig_eff->PrintStats();
    if (s_pair_reco_eff) s_pair_reco_eff->PrintStats();
}

// The ONE definition of the per-pair efficiency weight columns. It used to be written out five
// times, verbatim (crossx OS, crossx SS, the two no-minv template passes, and the generic
// histograms); any correction added to one and not the others silently produced two different
// cross-sections in the same file.
ROOT::RDF::RNode RDFBasedHistFillingPP::AddPairEfficiencyWeightColumns(ROOT::RDF::RNode df)
{
    return df
        .Define("q_eta1", "(float)(m1.charge * m1.eta)")
        .Define("q_eta2", "(float)(m2.charge * m2.eta)")
        // DATA tag-and-probe single-muon mu4 turn-on, per leg.
        .Define("effcy1", [](int q, float pt, float qe) {
            return RDFBasedHistFillingData::EvaluateSingleMuonEffcy("", q > 0, pt, qe);
        }, {"m1.charge", "m1.pt", "q_eta1"})
        .Define("effcy2", [](int q, float pt, float qe) {
            return RDFBasedHistFillingData::EvaluateSingleMuonEffcy("", q > 0, pt, qe);
        }, {"m2.charge", "m2.pt", "q_eta2"})
        // eps_trig^pair, the HYBRID of the header block: below the split edge
        // eps^nc_1 * eps^nc_2 * eps_dR(dR; cell) (2mu4 is an AND of the two legs, dressed by the MC
        // dR correlation, exactly 1 for dR >= 1); in the last two coarse pair-pT bins the
        // single-value MC pair efficiency x SF_1 x SF_2, for which the evaluator needs each leg's
        // (pT, eta, charge) to form the MC twin of effcy1/effcy2. ONE column, ONE definition for
        // every pull (crossx OS, SS, template-fit pass, generic).
        .Define("effcy_pair", [](float dr, float pair_pt, float pair_eta, float e1, float e2,
                                 float pt1, float eta1, int q1, float pt2, float eta2, int q2) {
            return (e1 > 0 && e2 > 0)
                 ? s_pair_trig_eff->Eval(dr, pair_pt, pair_eta, e1, e2, pt1, eta1, q1, pt2, eta2, q2)
                 : -1.0;
        }, {"dr", "pair_pt", "pair_eta", "effcy1", "effcy2",
            "m1.pt", "m1.eta", "m1.charge", "m2.pt", "m2.eta", "m2.charge"})
        .Define("w_trig", "effcy_pair > 0 ? 1.0 / effcy_pair : 0.0")
        // pp24-fullsim PAIR reco efficiency. Already floored inside the evaluator; a -1 means no
        // level of the map carries a measurement, which the caller must treat as "no correction".
        .Define("effcy_reco_pair", [](float pair_pt, float pair_eta, float dr) {
            return s_pair_reco_eff->Eval(pair_pt, pair_eta, dr);
        }, {"pair_pt", "pair_eta", "dr"})
        .Define("w_reco", "effcy_reco_pair > 0 ? 1.0 / effcy_reco_pair : 1.0");
}

void RDFBasedHistFillingPP::OpenEffcyPtFitFile() {
    if (!s_effcy_pT_fit_map.empty()) {
        std::cout << "OpenEffcyPtFitFile: TF1 map already loaded (" << s_effcy_pT_fit_map.size() << " entries)" << std::endl;
        return;
    }
    std::string wpsuf = isTight ? "" : "_medium_wp";   // WP-matched trig-eff fit (nominal tight unsuffixed)
    std::string fit_path = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_20" + std::to_string(run_year) + "/trg_effcy_pT_fitting_to_erf_plus_log/single_mu_effcy_pT_fit" + wpsuf + ".root";
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
    std::string hist_path = base_dir + "/histograms_real_pairs_pp_20" + std::to_string(run_year) + "_single_mu4_coarse_q_eta_bin_qeta_fid" + wpsuf + ".root";
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
// ApplyMuonWorkingPointFilter: the NOMINAL Tight muon working point on df_op/df_ss
// ============================================================================
//
// --- Muon working-point (WP) selection for the DATA spectra ---
// NOMINAL WP = TIGHT (isTight=true). Tight is a subset of Medium and the pair-level Tight flag
// (both muons quality&16) is already serialized in the data pair tree, so we select it here with
// a Filter -- NO ntuple reprocessing (verified: OP frac ~0.91). Setting isTight=false recovers
// the Medium spectrum (WP systematic; distinct _medium_wp output).
//
// WHY IT LIVES HERE AND NOT INSIDE FillHistogramsCrossx (fixed 2026-08-25): the driver
// RDFBasedHistFillingData::FillHistograms() calls FillHistogramsGeneric() BEFORE
// FillHistogramsCrossx(). While the Tight Filter was applied inside FillHistogramsCrossx, every
// GENERIC histogram was therefore a MEDIUM-WP yield -- yet it was divided by the TIGHT
// efficiencies (Tight pair eps_reco, Tight single-muon turn-on fits) that
// AddPairEfficiencyWeightColumns builds. A yield and a correction from two different working
// points is a physics error: the WP cancels in neither eps_reco nor eps_trig, so the
// MC-vs-data comparison read a spectrum that was ~10% too large and mis-corrected on top.
// Applying the WP here, above BOTH families, puts the yield and its correction on one WP.
//
// Applied to BOTH df_op and df_ss so every downstream pull (generic, gapcut, crossx, SS crossx,
// template-fit) inherits it; only the quality bit (8->16) changes.
// (docs/muon_wp_registry.md 3; tight_wp_default_change.md S2.3)
void RDFBasedHistFillingPP::ApplyMuonWorkingPointFilter(){
    if (!isTight) {
        // MEDIUM escape hatch, unchanged: no WP Filter at all (the pair tree is already Medium).
        std::cout << "[PP] ApplyMuonWorkingPointFilter: isTight=false -> Medium WP, no filter"
                  << std::endl;
        return;
    }
    if (wp_filter_applied) return;   // exactly once: a second Filter would be a redundant node

    df_map.at("df_op") = map_at_checked(df_map, "df_op", "ApplyMuonWorkingPointFilter: df_op").Filter("pair_pass_tight", "tight WP (pair)");
    df_map.at("df_ss") = map_at_checked(df_map, "df_ss", "ApplyMuonWorkingPointFilter: df_ss").Filter("pair_pass_tight", "tight WP (pair)");
    wp_filter_applied = true;
    std::cout << "[PP] ApplyMuonWorkingPointFilter: TIGHT WP filter applied to df_op and df_ss"
              << std::endl;
}

// ============================================================================
// FillHistogramsGeneric: override to add 2mu4 trigger efficiency weighting
// ============================================================================

void RDFBasedHistFillingPP::FillHistogramsGeneric(){
    if (!trigger_effcy_calc) {
        OpenEffcyPtFitFile();
        OpenPairEfficiencyInputs();           // MC eps_dR (2mu4) + pp24-fullsim PAIR eps_reco

        // WORKING POINT FIRST -- before the gap cut, before the efficiency weight columns and
        // before ANY histogram is booked. The generic histograms are corrected by the Tight
        // efficiencies below, so the yield they correct must be the Tight yield.
        ApplyMuonWorkingPointFilter();

        // THE FIDUCIAL GAP CUT IS MANDATORY HERE, not optional (2026-08-17). These generic
        // dataframes are efficiency-corrected below, and the single-muon turn-on is fitted on a
        // CONTIGUOUS q*eta binning that stops at 2.20 -- the top edge of
        // CommonEffcyConfig::q_eta_proj_ranges_coarse_incl_gap, chosen to match the gap cut's
        // forward window. A muon at q*eta in [2.20, 2.40] therefore has NO fitted turn-on, and
        // EvaluateSingleMuonEffcyPtFitted deliberately THROWS rather than returning a sentinel
        // that would silently delete the pair (pp_trig_eff_highpt_jump.md). Applying the same cut
        // the crossx signal region applies is what makes the efficiency defined for every pair
        // that survives -- and it puts the generic/MC-data-comparison histograms in the SAME
        // fiducial region as the cross-section. Idempotent w.r.t. the later signal-cut filters,
        // which contain the same expression.
        // The PAIR-LEVEL gap cut |eta^pair| < 2.2 travels with the single-muon windows
        // everywhere they are applied to BOTH legs (user, 2026-09-07; ParamsSet.h).
        const std::string fiducial_gap_cut =
            ParamsSet::FiducialGapCutExpr("m1.charge * m1.eta") + " && "
          + ParamsSet::FiducialGapCutExpr("m2.charge * m2.eta") + " && "
          + ParamsSet::PairFiducialEtaCutExpr("pair_eta");

        for (const std::string& category : categories_essential) {
            std::string df_name = "df" + category;
            ROOT::RDF::RNode df = map_at_checked(df_map, df_name,
                Form("PP::FillHistogramsGeneric: df_map.at(%s)", df_name.c_str()))
                .Filter(fiducial_gap_cut, "fiducial gap cut (both muons + pair eta)");

            // Generic analysis histograms (incl. the gapcut histos read by the
            // MC-data comparison) are weighted by the SAME per-pair efficiency
            // correction as the crossx: w_reco_trig = w_reco * w_trig, both built by
            // AddPairEfficiencyWeightColumns. That is what keeps the MC-data comparison
            // and the cross-section on one and the same correction.
            // ..._over_dr is the SAME weight divided by dR: the dR Jacobian that turns dN/ddR
            // into the (1/dR) dN/ddR density behind the "_jacobian_corrected" histograms.
            // Form copied from the truth precedent RDFBasedHistFillingPythiaTruth.cxx:96
            // (`weight_over_dr`), including the dr>0 guard. Before 2026-08-25 the
            // "_jacobian_corrected" data histograms were booked with w_reco_trig itself, so they
            // were bit-identical clones of their un-corrected twins -- a name with no weighting
            // behind it.
            ROOT::RDF::RNode df_with_trig = AddPairEfficiencyWeightColumns(df)
                .Define("w_reco_trig", "w_reco * w_trig")
                .Define("w_reco_trig_over_dr", "dr > 0 ? (w_reco * w_trig) / dr : 0.0");

            df_map.erase(df_name);
            df_map.emplace(df_name, df_with_trig);
        }
        generic_weight_col = "w_reco_trig";  // reco+trig corrected (was "w_trig")
        generic_jacobian_weight_col = "w_reco_trig_over_dr";  // = generic_weight_col / dR
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
    OpenPairEfficiencyInputs();           // MC eps_dR (2mu4) + pp24-fullsim PAIR eps_reco

    std::cout << "[PP] FillHistogramsCrossx: opposite-sign only, signal cuts, "
              << "crossx_weight = weight * " << pp_crossx_lumi_factor
              << " (1/L_int), with 2mu4 trig eff correction" << std::endl;

    // SIGNAL REGION (docs/analysis_overview.md §2). Since 2026-08-17 the per-muon one-sided
    // `q*eta < 2.2` is REPLACED by the detector-gap FIDUCIAL cut, required of BOTH muons
    // (user instruction 2026-08-17; docs/tracking/pp24_crossx_rerun_2026_08.md). The windows are
    // READ from ParamsSet::single_mu_fiducial_gap_cuts -- the single source of truth; the values
    // are NEVER retyped here (muon_gap_cuts_acceptance.md F10/F11a/F17).
    // 2026-09-07 the PAIR-LEVEL gap cut |eta^pair| < ParamsSet::pair_eta_fiducial_max = 2.2 was
    // ADDED alongside them (user): a pair can reach |eta^pair| > 2.2 with both muons passing the
    // one-sided q*eta windows, but only through a charge-dependent corner of phase space whose
    // pair efficiency is shaped by the single-muon cut itself. Removed rather than modelled.
    // The forward window {2.20, 2.40} together with the ntuple-level |eta| < 2.4 makes the
    // effective forward edge 2.20, which is exactly the top edge of the CONTIGUOUS coarse q*eta
    // turn-on binning (CommonEffcyConfig::q_eta_proj_ranges_coarse_incl_gap) -- so every
    // surviving muon has a fitted trigger efficiency and the "no fitted turn-on" throw in
    // EvaluateSingleMuonEffcyPtFitted stays unreachable. Blast radius:
    // docs/signal_selection_change_impact.md.
    const std::string signal_cuts =
        ParamsSet::SignalMinvCutExpr("minv") + " && " + ParamsSet::SignalPairPtCutExpr("pair_pt") + " && "
        + ParamsSet::FiducialGapCutExpr("m1.charge * m1.eta") + " && "
        + ParamsSet::FiducialGapCutExpr("m2.charge * m2.eta") + " && "
        + ParamsSet::PairFiducialEtaCutExpr("pair_eta");

    // --- Muon working-point (WP) selection for the DATA crossx spectrum ---
    // Moved into ApplyMuonWorkingPointFilter() (2026-08-25) so that FillHistogramsGeneric --
    // which the driver runs FIRST -- gets the same Tight yield it is corrected with; see the
    // physics comment on that function. The call is idempotent, so:
    //   * output_generic_hists == true  -> generic already applied it; this is a no-op;
    //   * output_generic_hists == false -> generic never ran; this applies it here.
    // Either way the crossx spectrum is Tight, exactly as before, and the isTight==false
    // (Medium) escape hatch is untouched. FillHistogramsCrossx is only ever reached with
    // trigger_effcy_calc == false (see RDFBasedHistFillingData::FillHistograms), which is
    // exactly the branch of FillHistogramsGeneric that applies the WP -- so the two callers
    // can never disagree about whether the WP is on.
    ApplyMuonWorkingPointFilter();

    ROOT::RDF::RNode df_op_base = map_at_checked(df_map, "df_op", "FillHistogramsCrossx: df_op");
    ROOT::RDF::RNode df_single_b_crossx = df_op_base.Filter(signal_cuts);
    if (df_map.find("df_single_b_crossx") == df_map.end()) {
        df_map.emplace("df_single_b_crossx", df_single_b_crossx);
    }

    // Per-pair 2mu4 trigger efficiency: the HYBRID of PairTrigEffCrossxEvaluator (eps_1 * eps_2 *
    // eps_dR below 74.24 GeV, the single-value MC pair efficiency x SF_1 x SF_2 above).
    // If FillHistogramsGeneric already added trigger columns to df_op, they
    // propagate through the Filter to df_single_b_crossx — skip re-Define.
    // Per-pair trigger (w_trig) AND pp24-fullsim PAIR reco-eff (w_reco) columns.
    // If FillHistogramsGeneric already added them to df_op (generic ran; it now
    // defines the full trig+reco chain), they propagate through the signal-cut
    // Filter -> reuse them (skip re-Define to avoid an RDF column collision).
    ROOT::RDF::RNode df_with_trig = generic_weight_col.empty()
        ? AddPairEfficiencyWeightColumns(df_single_b_crossx)
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
        // NOMINAL corrected weight = trigger (incl. the MC eps_dR correlation correction) x the
        // pp24-fullsim PAIR reco efficiency. Every crossx histogram filled with
        // crossx_weight_trig_corr is reco+trig corrected; it equals the _corr_unfolded_reco_trig
        // stage (w_unfold == 1). To revert to trig-only, drop "w_reco *".
        // (docs/tracking/pp24_crossx_rerun_2026_08.md; the Run-2 single-muon placeholder this
        // replaced is docs/tracking/reco_eff_placeholder_run2.md)
        .Define("crossx_weight_trig_corr", "crossx_weight * w_reco * w_trig")
        .Define("cw_raw",                "crossx_weight")
        .Define("cw_unfolded",           "cw_raw * w_unfold")
        .Define("cw_unfolded_reco",      "cw_unfolded * w_reco")
        .Define("cw_unfolded_reco_trig", "cw_unfolded_reco * w_trig")
        // TRIGGER-FIRST intermediate stage (user, 2026-09-17): the comparison figure shows
        // uncorrected -> + trigger -> + trigger + reco, so the trigger correction's own impact is
        // visible. Same final weight (the two corrections commute); CorrectionStages.h.
        .Define("cw_unfolded_trig",      "cw_unfolded * w_trig");

    if (df_map.find("df_single_b_crossx_weighted") == df_map.end()) {
        df_map.emplace("df_single_b_crossx_weighted", df_single_b_crossx_weighted);
    }

    // Store RResultPtrs (lazy-evaluated) instead of immediate clones.
    // They will be evaluated during HistPostProcess() and converted to raw pointers.
    // DEFAULT fine pair-pT crossx axis = pT_bins_150 (16 log bins 9 -> 150 GeV), user
    // decision D5 (docs/tracking/mu_pt45_gap125_pairpt9_adoption.md). Chosen so every
    // pair_pt_coarse_bins edge IS a fine edge (coarse k = fine 2k), which the 9 -> 120
    // alternative below deliberately does NOT satisfy.
    // The UNSUFFIXED family is the complete nominal set (minv, dR, no_trig_corr, the
    // correction stages, same-sign, the 3Ds and the R_AA globals), so putting the default
    // here -- rather than adding _pt_150 twins for each -- is what keeps R_AA and the
    // cross-section on ONE pair-pT binning.
    const int    npt   = (int)(pms.pT_bins_150.size() - 1);
    const double* ptbins = pms.pT_bins_150.data();

    // --- LOW-MASS TEMPLATE-FIT pass (low_mass_template_calc) ---
    // Reads _no_res_cut (resonances PRESENT; selected in SetIOPathsHook), distinct output.
    // Produces ONLY the 0-4 GeV minv template spectra D_OS/D_SS (1D + 2D vs pair pT/eta)
    // for the low-mass dimuon template fit (docs/tracking/low_mass_dimuon_template_fit.md
    // 3a). Returns early so the signal-region crossx (minv in [1.08,2.9], WRONG from
    // _no_res_cut because of resonance leakage) is NOT produced. Selection = signal_cuts
    // MINUS minv window, NO dR. Weight = TRIGGER-ONLY reco-level dsigma (1/L * w_trig, NO
    // reco-eff): the minv fit runs at reco, BEFORE reco-eff/unfolding (3e ordering reversal).
    if (low_mass_template_calc) {
        // signal_cuts MINUS the minv window; the gap cuts are kept IDENTICAL to signal_cuts
        // -- both muons' q*eta windows AND the pair-level |eta^pair| < 2.2, all read from
        // ParamsSet.
        const std::string signal_cuts_no_minv =
            ParamsSet::SignalPairPtCutExpr("pair_pt") + " && "
            + ParamsSet::FiducialGapCutExpr("m1.charge * m1.eta") + " && "
            + ParamsSet::FiducialGapCutExpr("m2.charge * m2.eta") + " && "
            + ParamsSet::PairFiducialEtaCutExpr("pair_eta");
        const double lumi_factor_tmpl = pp_crossx_lumi_factor;
        // TEMPLATE-FIT DATA INPUT = TRIGGER-ONLY, RECONSTRUCTED LEVEL (correction-ordering
        // reversal 2026-07-01, low_mass_dimuon_template_fit.md 3e). The minv template fit is
        // performed at RECO level, AFTER trigger-efficiency correction (the trigger fires on reco
        // objects) but BEFORE reconstruction-efficiency correction and unfolding. Reco-eff +
        // unfolding are applied to the EXTRACTED signal yield AFTER the fit (they are properties of
        // the detector, origin-blind), NOT to the fit input. Rationale: the mixed-event
        // combinatoric template uses REAL-DATA muons (reco quantities) and the fake/hadronic
        // background has NO truth match, so the fit lives at reco. Hence weight is trigger-only
        // (NO w_reco). The nominal signal-region crossx (above, reco+trig) is a SEPARATE object.
        auto attach_crossx_weight = [&](ROOT::RDF::RNode node) -> ROOT::RDF::RNode {
            // The template-fit input is TRIGGER-ONLY by design (it uses w_trig, not w_reco), but
            // the columns come from the SAME helper so its trigger correction can never differ
            // from the crossx's.
            ROOT::RDF::RNode n = generic_weight_col.empty()
                ? AddPairEfficiencyWeightColumns(node)
                : node;
            return n
                .Define("crossx_weight", [lumi_factor_tmpl](double weight){ return weight * lumi_factor_tmpl; }, {"weight"})
                .Define("crossx_weight_trig_only", "crossx_weight * w_trig");
        };
        ROOT::RDF::RNode df_op_t = attach_crossx_weight(map_at_checked(df_map, "df_op", "FillHistogramsCrossx PP: df_op (template)").Filter(signal_cuts_no_minv));
        ROOT::RDF::RNode df_ss_t = attach_crossx_weight(map_at_checked(df_map, "df_ss", "FillHistogramsCrossx PP: df_ss (template)").Filter(signal_cuts_no_minv));
        const int npt150 = (int)(pms.pT_bins_150.size() - 1);
        const double* ptb150 = pms.pT_bins_150.data();
        hist1d_rresultptr_map["h1d_crossx_minv_0_4_op_dsigma"] = df_op_t.Histo1D(ROOT::RDF::TH1DModel("h1d_crossx_minv_0_4_op_dsigma", ";m_{#mu#mu} [GeV];d#sigma/dm_{#mu#mu} [pb GeV^{-1}]", 50, 0.0, 4.0), "minv", "crossx_weight_trig_only");
        hist1d_rresultptr_map["h1d_crossx_minv_0_4_ss_dsigma"] = df_ss_t.Histo1D(ROOT::RDF::TH1DModel("h1d_crossx_minv_0_4_ss_dsigma", ";m_{#mu#mu} [GeV];d#sigma/dm_{#mu#mu} [pb GeV^{-1}]", 50, 0.0, 4.0), "minv", "crossx_weight_trig_only");
        hist2d_rresultptr_map["h2d_crossx_minv_0_4_vs_pair_pt_log_150_op_dsigma"] = df_op_t.Histo2D(ROOT::RDF::TH2DModel("h2d_crossx_minv_0_4_vs_pair_pt_log_150_op_dsigma", ";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]", npt150, ptb150, 50, 0.0, 4.0), "pair_pt", "minv", "crossx_weight_trig_only");
        hist2d_rresultptr_map["h2d_crossx_minv_0_4_vs_pair_pt_log_150_ss_dsigma"] = df_ss_t.Histo2D(ROOT::RDF::TH2DModel("h2d_crossx_minv_0_4_vs_pair_pt_log_150_ss_dsigma", ";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]", npt150, ptb150, 50, 0.0, 4.0), "pair_pt", "minv", "crossx_weight_trig_only");
        hist2d_rresultptr_map["h2d_crossx_minv_0_4_vs_pair_eta_op_dsigma"] = df_op_t.Histo2D(ROOT::RDF::TH2DModel("h2d_crossx_minv_0_4_vs_pair_eta_op_dsigma", ";#eta^{pair};m_{#mu#mu} [GeV]", 24, -2.4, 2.4, 50, 0.0, 4.0), "pair_eta", "minv", "crossx_weight_trig_only");
        hist2d_rresultptr_map["h2d_crossx_minv_0_4_vs_pair_eta_ss_dsigma"] = df_ss_t.Histo2D(ROOT::RDF::TH2DModel("h2d_crossx_minv_0_4_vs_pair_eta_ss_dsigma", ";#eta^{pair};m_{#mu#mu} [GeV]", 24, -2.4, 2.4, 50, 0.0, 4.0), "pair_eta", "minv", "crossx_weight_trig_only");

        // EXTENDED-MASS 0-20 GeV (40 bins, matching the extended_mass_control_region MC by-origin
        // histos) for the THStack data-vs-MC + control-region study: 1D pair-pT-integrated (pair_pt > ParamsSet::signal_pair_pt_min)
        // and 2D vs the NOMINAL COARSE pair-pT bins (pms.pair_pt_coarse_bins). Same signal_cuts_no_minv
        // selection + trigger-only reco-level weight. In mixed_event_template mode these become the
        // extended mixed-event combinatoric T_mix (from the scrambled input).
        const int    nptc  = ParamsSet::N_COARSE_PAIR_PT_BINS;
        const double* ptc  = pms.pair_pt_coarse_bins.data();
        hist1d_rresultptr_map["h1d_crossx_minv_0_20_op_dsigma"] = df_op_t.Histo1D(ROOT::RDF::TH1DModel("h1d_crossx_minv_0_20_op_dsigma", ";m_{#mu#mu} [GeV];d#sigma/dm_{#mu#mu} [pb GeV^{-1}]", 40, 0.0, 20.0), "minv", "crossx_weight_trig_only");
        hist1d_rresultptr_map["h1d_crossx_minv_0_20_ss_dsigma"] = df_ss_t.Histo1D(ROOT::RDF::TH1DModel("h1d_crossx_minv_0_20_ss_dsigma", ";m_{#mu#mu} [GeV];d#sigma/dm_{#mu#mu} [pb GeV^{-1}]", 40, 0.0, 20.0), "minv", "crossx_weight_trig_only");
        hist2d_rresultptr_map["h2d_crossx_minv_0_20_vs_pair_pt_coarse_op_dsigma"] = df_op_t.Histo2D(ROOT::RDF::TH2DModel("h2d_crossx_minv_0_20_vs_pair_pt_coarse_op_dsigma", ";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]", nptc, ptc, 40, 0.0, 20.0), "pair_pt", "minv", "crossx_weight_trig_only");
        hist2d_rresultptr_map["h2d_crossx_minv_0_20_vs_pair_pt_coarse_ss_dsigma"] = df_ss_t.Histo2D(ROOT::RDF::TH2DModel("h2d_crossx_minv_0_20_vs_pair_pt_coarse_ss_dsigma", ";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]", nptc, ptc, 40, 0.0, 20.0), "pair_pt", "minv", "crossx_weight_trig_only");

        // NO-PAIR-SELECTION trigger-corrected 0-4 GeV data dsigma (muon-level selection only,
        // inherent in the ntuple: pt>4.5, |eta|<2.4, quality, dp/p). NO pair_pt / q*eta cut. Same
        // TRIGGER-ONLY reco-level weight (crossx_weight_trig_only = 1/L * w_trig, NO reco-eff) as
        // above, for the bkg_mc_provenance data-vs-fullsim comparison: the Pythia fullsim MC uses
        // RECONSTRUCTED quantities and therefore already carries the same reconstruction efficiency,
        // so correcting the data for reco-eff would invalidate the data-vs-MC comparison (2026-07-01).
        // CAVEAT: 1/eff_trig is only well defined on the trigger plateau (pair pT >~ 8); for soft
        // muons near the pT>4.5 threshold the trig-eff placeholder clamps, so the corrected soft
        // (low-mass) spectrum carries turn-on/clamp artifacts -- approximate there.
        // "_nosel" = no SIGNAL-REGION cuts (no minv window, no pair-pT threshold). It is NOT "no
        // cuts at all": these histograms are trigger-corrected, and the single-muon turn-on is
        // only fitted inside the fiducial q*eta region, so a muon in a gap window has no
        // efficiency and EvaluateSingleMuonEffcyPtFitted throws by design. The fiducial gap cut is
        // therefore a PREREQUISITE of the weight, not part of the selection being switched off.
        const std::string fiducial_gap_cut_nosel =
            ParamsSet::FiducialGapCutExpr("m1.charge * m1.eta") + " && "
          + ParamsSet::FiducialGapCutExpr("m2.charge * m2.eta") + " && "
          + ParamsSet::PairFiducialEtaCutExpr("pair_eta");
        ROOT::RDF::RNode df_op_nosel = attach_crossx_weight(map_at_checked(df_map, "df_op", "FillHistogramsCrossx PP: df_op (template nosel)").Filter(fiducial_gap_cut_nosel));
        ROOT::RDF::RNode df_ss_nosel = attach_crossx_weight(map_at_checked(df_map, "df_ss", "FillHistogramsCrossx PP: df_ss (template nosel)").Filter(fiducial_gap_cut_nosel));
        hist1d_rresultptr_map["h1d_crossx_minv_0_4_op_dsigma_nosel"] = df_op_nosel.Histo1D(ROOT::RDF::TH1DModel("h1d_crossx_minv_0_4_op_dsigma_nosel", ";m_{#mu#mu} [GeV];d#sigma/dm_{#mu#mu} [pb GeV^{-1}]", 50, 0.0, 4.0), "minv", "crossx_weight_trig_only");
        hist1d_rresultptr_map["h1d_crossx_minv_0_4_ss_dsigma_nosel"] = df_ss_nosel.Histo1D(ROOT::RDF::TH1DModel("h1d_crossx_minv_0_4_ss_dsigma_nosel", ";m_{#mu#mu} [GeV];d#sigma/dm_{#mu#mu} [pb GeV^{-1}]", 50, 0.0, 4.0), "minv", "crossx_weight_trig_only");
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
    // FINE pair-eta axis: SINGLE SOURCE = ParamsSet (N_PAIR_ETA_CROSSX_BINS / PAIR_ETA_CROSSX_MIN
    // / PAIR_ETA_CROSSX_MAX and the generated pair_eta_crossx_bins edge vector, which is also
    // what hist_binning_map["pair_eta_crossx"] serves to the 1D views). The literal
    // "44, -2.4, 2.4" is never retyped here again. The fixed-bin TH2 axes below keep the
    // FIXED-bin constructor (nbins, min, max) so their on-disk representation is byte-for-byte
    // what it was; only the TH3 axes, which always needed an explicit edge array, take the
    // vector (its edges agree with the fixed-bin axis to <= 8.9e-16, i.e. to double rounding).
    const int     n_eta_crossx = ParamsSet::N_PAIR_ETA_CROSSX_BINS;
    const double* eta_edges    = pms.pair_eta_crossx_bins.data();
    const auto minv_edges = make_unif_edges(50,  1.0,  3.0);
    const auto dr_edges   = make_unif_edges(50,  0.0, 1.0);

    // ------------------------------------------------------------------------------------------
    // SIGNAL-REGION 1D DIFFERENTIAL CROSS SECTIONS (both signs)
    // ------------------------------------------------------------------------------------------
    // These are 1D views of EXACTLY the object the 2D/3D crossx histograms above describe: the
    // SAME node (OS: df_single_b_crossx_weighted, i.e. df_op after the Tight WP filter and the
    // single-b signal cuts; SS: the mirror node built in the SS block below) and the SAME weight
    // crossx_weight_trig_corr = weight * (1/L_int) * w_reco * w_trig. Nothing is re-selected and
    // nothing is re-weighted here, so a 1D panel and the 2D/3D panel beside it can never
    // describe different cells.
    //
    // Every axis is requested BY NAME from hist_binning_map (registered once in
    // RDFBasedHistFillingBaseClass::BuildHistBinningMapBaseCommon) -- the pair-eta axis is
    // "pair_eta_crossx", the very edge vector the 2D/3D views use, and the pair-pT axis is
    // "pT_bins_150", the axis of the UNSUFFIXED default family
    // (h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts). No binning is
    // retyped in this file.
    //
    // NOT scaled by bin width at fill time: like h1d_crossx_minv_0_4_*_dsigma, these hold
    // sum(w) per bin and the plotter applies Scale(norm, "width"). The y titles therefore name
    // the density the plotter produces, with the unit of the observable: dimensionless
    // observables (dR, dphi, deta, pair eta) -> [pb], dimensionful ones (minv, pair pT) ->
    // [pb GeV^{-1}].
    //   {hist tag, RDF column, hist_binning_map name, x title, y title}
    const std::vector<std::array<std::string, 5>> signal_region_1d_vars = {
        {"DR_zoomin",   "dr",       "dr_zoomin_bins_1d",   "#DeltaR",            "d#sigma/d#DeltaR [pb]"},
        {"Dphi_zoomin", "dphi",     "dphi_zoomin_bins_1d", "#Delta#phi",         "d#sigma/d#Delta#phi [pb]"},
        {"Deta_zoomin", "deta",     "deta_zoomin_bins_1d", "#Delta#eta",         "d#sigma/d#Delta#eta [pb]"},
        {"minv_zoomin", "minv",     "minv_zoomin_bins_1d", "m_{#mu#mu} [GeV]",   "d#sigma/dm_{#mu#mu} [pb GeV^{-1}]"},
        {"pair_eta",    "pair_eta", "pair_eta_crossx",     "#eta^{pair}",        "d#sigma/d#eta^{pair} [pb]"},
        {"pair_pt_150", "pair_pt",  "pT_bins_150",         "p_{T}^{pair} [GeV]", "d#sigma/dp_{T}^{pair} [pb GeV^{-1}]"}
    };

    auto book_signal_region_1d = [&](ROOT::RDF::RNode node,
                                     const std::string& weight_col,
                                     const std::string& sign_suffix) {
        for (const auto& v : signal_region_1d_vars) {
            const std::string hname =
                "h1d_crossx_" + v[0] + "_w_signal_cuts" + sign_suffix + "_dsigma";
            const std::vector<double>& edges = map_at_checked(
                hist_binning_map, v[2],
                "FillHistogramsCrossx PP: hist_binning_map (signal-region 1D)");
            const std::string htitle = ";" + v[3] + ";" + v[4];
            hist1d_rresultptr_map[hname] = node.Histo1D(
                ROOT::RDF::TH1DModel(hname.c_str(), htitle.c_str(),
                                     static_cast<int>(edges.size()) - 1, edges.data()),
                v[1], weight_col);
        }
    };

    book_signal_region_1d(df_single_b_crossx_weighted, "crossx_weight_trig_corr", "_op");

    hist2d_rresultptr_map["h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts"] = df_single_b_crossx_weighted.Histo2D(
        ROOT::RDF::TH2DModel("h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts", ";p_{T}^{pair} [GeV];#eta^{pair}", npt, ptbins, ParamsSet::N_PAIR_ETA_CROSSX_BINS, ParamsSet::PAIR_ETA_CROSSX_MIN, ParamsSet::PAIR_ETA_CROSSX_MAX),
        "pair_pt", "pair_eta", "crossx_weight_trig_corr");
    hist2d_rresultptr_map["h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts_no_trig_corr"] = df_single_b_crossx_weighted.Histo2D(
        ROOT::RDF::TH2DModel("h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts_no_trig_corr", ";p_{T}^{pair} [GeV];#eta^{pair}", npt, ptbins, ParamsSet::N_PAIR_ETA_CROSSX_BINS, ParamsSet::PAIR_ETA_CROSSX_MIN, ParamsSet::PAIR_ETA_CROSSX_MAX),
        "pair_pt", "pair_eta", "crossx_weight");

    // --- RAW PAIR COUNTS in the signal region -------------------------------------------------
    // The UNWEIGHTED twin of h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts: SAME node
    // (df_single_b_crossx_weighted = Tight WP + signal_cuts on the OS tree) and SAME axes, but
    // NO weight column, so the bin content is an integer number of pairs and the error is
    // sqrt(N). This is the statistical reach of the measurement, which the corrected spectrum
    // hides -- a cell can sit at a healthy dsigma/dpT and rest on three pairs.
    //
    // Deliberately NOT derived from the "_corr_raw" stage histogram: that one is
    // crossx_weight = weight * (1/L_int), so it is the count divided by the luminosity AND
    // carries a weighted (not Poisson) error. Counting needs an unweighted fill.
    //
    // ONE way in which N and dsigma do NOT describe the same pair set, and it is deliberate:
    // a pair whose trigger efficiency could not be evaluated gets the sentinel w_trig = 0
    // (:392-393 above, and the ACTIVE docs/tracking/pp_trig_eff_highpt_jump.md), so it is
    // COUNTED here but contributes exactly 0 to every crossx histogram. N is therefore the raw
    // statistical reach of the selection; dsigma is built from the subset that could be
    // corrected. Measured 2026-09-03: no cell is emptied by this (0 cells with counts > 0 and
    // crossx == 0), so the effect is partial, never gross.
    // Consumed by plotting_codes/single_b_analysis/plot_pp_counts_pair_pt_in_eta.cxx.
    hist2d_rresultptr_map["h2d_counts_pair_pt_pair_eta_binned_w_signal_cuts"] = df_single_b_crossx_weighted.Histo2D(
        ROOT::RDF::TH2DModel("h2d_counts_pair_pt_pair_eta_binned_w_signal_cuts", ";p_{T}^{pair} [GeV];#eta^{pair}", npt, ptbins, ParamsSet::N_PAIR_ETA_CROSSX_BINS, ParamsSet::PAIR_ETA_CROSSX_MIN, ParamsSet::PAIR_ETA_CROSSX_MAX),
        "pair_pt", "pair_eta");


    // Correction-stage histograms (raw -> unfolded -> +reco -> +reco+trig) for the
    // primary pair_pt x pair_eta differential, so each correction's impact is
    // visible at plotting time. See CorrectionStages.h.
    for (const auto& st : CrossxCorrectionStages()) {
        const std::string nm = std::string("h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts") + st.suffix;
        hist2d_rresultptr_map[nm] = df_single_b_crossx_weighted.Histo2D(
            ROOT::RDF::TH2DModel(nm.c_str(), ";p_{T}^{pair} [GeV];#eta^{pair}", npt, ptbins, ParamsSet::N_PAIR_ETA_CROSSX_BINS, ParamsSet::PAIR_ETA_CROSSX_MIN, ParamsSet::PAIR_ETA_CROSSX_MAX),
            "pair_pt", "pair_eta", st.weight_col);
    }

    // --- Same-sign (SS) signal-region 2D yield, for the OS-SS combinatorial
    // subtraction in R_AA (analysis_overview.md §4a). Mirrors the OS path with
    // identical cuts and corrected weight (crossx_weight · w_reco · w_trig).
    // Trigger columns come from FillHistogramsGeneric (added to df_ss too, since
    // categories_essential = pair_signs) when generic ran; else define inline.
    {
        ROOT::RDF::RNode df_ss_crossx = map_at_checked(df_map, "df_ss", "FillHistogramsCrossx PP: df_ss").Filter(signal_cuts);
        // Full trig+reco efficiency chain -- added only if generic did NOT run (else reuse the
        // columns FillHistogramsGeneric already added to df_ss).
        ROOT::RDF::RNode df_ss_with_trig = generic_weight_col.empty()
            ? AddPairEfficiencyWeightColumns(df_ss_crossx)
            : df_ss_crossx;
        ROOT::RDF::RNode df_ss_weighted = df_ss_with_trig
            .Define("crossx_weight",
                [lumi_factor](double weight){ return weight * lumi_factor; }, {"weight"})
            .Define("crossx_weight_trig_corr", "crossx_weight * w_reco * w_trig");
        hist2d_rresultptr_map["h2d_ss_crossx_pair_pt_pair_eta_binned_w_signal_cuts"] = df_ss_weighted.Histo2D(
            ROOT::RDF::TH2DModel("h2d_ss_crossx_pair_pt_pair_eta_binned_w_signal_cuts", ";p_{T}^{pair} [GeV];#eta^{pair}", npt, ptbins, ParamsSet::N_PAIR_ETA_CROSSX_BINS, ParamsSet::PAIR_ETA_CROSSX_MIN, ParamsSet::PAIR_ETA_CROSSX_MAX),
            "pair_pt", "pair_eta", "crossx_weight_trig_corr");

        // SS half of the signal-region 1D set: SAME node/weight as the SS 2D yield above, and the
        // SAME axes as the OS half, so D_OS - D_SS is a clean bin-by-bin combinatoric
        // subtraction (analysis_overview.md 4a).
        book_signal_region_1d(df_ss_weighted, "crossx_weight_trig_corr", "_ss");
    }

    hist2d_rresultptr_map["h2d_crossx_pair_pt_minv_w_signal_cuts"] = df_single_b_crossx_weighted.Histo2D(
        ROOT::RDF::TH2DModel("h2d_crossx_pair_pt_minv_w_signal_cuts", ";p_{T}^{pair} [GeV];m_{#mu#mu} [GeV]", npt, ptbins, 50, 1.0, 3.0),
        "pair_pt", "minv", "crossx_weight_trig_corr");
    hist2d_rresultptr_map["h2d_crossx_pair_pt_dr_w_signal_cuts"] = df_single_b_crossx_weighted.Histo2D(
        ROOT::RDF::TH2DModel("h2d_crossx_pair_pt_dr_w_signal_cuts", ";p_{T}^{pair} [GeV];#DeltaR", npt, ptbins, 50, 0.0, 1.0),
        "pair_pt", "dr", "crossx_weight_trig_corr");

    hist3d_rresultptr_map["h3d_crossx_minv_vs_pair_eta_vs_pair_pt_w_signal_cuts"] = df_single_b_crossx_weighted.Histo3D(
        ROOT::RDF::TH3DModel("h3d_crossx_minv_vs_pair_eta_vs_pair_pt_w_signal_cuts", ";p_{T}^{pair} [GeV];#eta^{pair};m_{#mu#mu} [GeV]", npt, ptbins, n_eta_crossx, eta_edges, 50, minv_edges.data()),
        "pair_pt", "pair_eta", "minv", "crossx_weight_trig_corr");
    hist3d_rresultptr_map["h3d_crossx_dr_vs_pair_eta_vs_pair_pt_w_signal_cuts"] = df_single_b_crossx_weighted.Histo3D(
        ROOT::RDF::TH3DModel("h3d_crossx_dr_vs_pair_eta_vs_pair_pt_w_signal_cuts", ";p_{T}^{pair} [GeV];#eta^{pair};#DeltaR", npt, ptbins, n_eta_crossx, eta_edges, 50, dr_edges.data()),
        "pair_pt", "pair_eta", "dr", "crossx_weight_trig_corr");

    // --- OPT-IN `_pt_120` ALTERNATIVE VIEW (pT_bins_120 = 16 log bins 9 -> 120 GeV) ---
        // Display variant only: it does NOT nest inside pair_pt_coarse_bins, so it must
        // never bin a correction. Partial coverage is intentional -- the default
        // unsuffixed family above carries the complete set.
    {
        const int    npt120    = (int)(pms.pT_bins_120.size() - 1);
        const double* ptbins120 = pms.pT_bins_120.data();
        const auto dr_edges120  = make_unif_edges(50, 0.0, 1.0);

        hist2d_rresultptr_map["h2d_crossx_pt_120_pair_eta_binned_w_signal_cuts"] = df_single_b_crossx_weighted.Histo2D(
            ROOT::RDF::TH2DModel("h2d_crossx_pt_120_pair_eta_binned_w_signal_cuts", ";p_{T}^{pair} [GeV];#eta^{pair}", npt120, ptbins120, ParamsSet::N_PAIR_ETA_CROSSX_BINS, ParamsSet::PAIR_ETA_CROSSX_MIN, ParamsSet::PAIR_ETA_CROSSX_MAX),
            "pair_pt", "pair_eta", "crossx_weight_trig_corr");
        // The UNWEIGHTED twin of h2d_crossx_pt_120_..., same pattern as
        // h2d_counts_pair_pt_pair_eta_binned_w_signal_cuts above but on the pT_bins_120 axis:
        // SAME node, no weight column, so the alternative view carries its own direct Poisson
        // raw counts rather than being inferred from the default axis's bins.
        // (Before 2026-09-08 this comment argued the opposite way round, because the 150 axis
        // was then the OPT-IN one and reached beyond the default's 120 GeV ceiling. The roles
        // are now swapped: the default reaches 150 and this alternative stops at 120.)
        hist2d_rresultptr_map["h2d_counts_pt_120_pair_eta_binned_w_signal_cuts"] = df_single_b_crossx_weighted.Histo2D(
            ROOT::RDF::TH2DModel("h2d_counts_pt_120_pair_eta_binned_w_signal_cuts", ";p_{T}^{pair} [GeV];#eta^{pair}", npt120, ptbins120, ParamsSet::N_PAIR_ETA_CROSSX_BINS, ParamsSet::PAIR_ETA_CROSSX_MIN, ParamsSet::PAIR_ETA_CROSSX_MAX),
            "pair_pt", "pair_eta");
        hist3d_rresultptr_map["h3d_crossx_dr_vs_pair_eta_vs_pt_120_w_signal_cuts"] = df_single_b_crossx_weighted.Histo3D(
            ROOT::RDF::TH3DModel("h3d_crossx_dr_vs_pair_eta_vs_pt_120_w_signal_cuts", ";p_{T}^{pair} [GeV];#eta^{pair};#DeltaR", npt120, ptbins120, ParamsSet::N_PAIR_ETA_CROSSX_BINS, pms.pair_eta_crossx_bins.data(), 50, dr_edges120.data()),
            "pair_pt", "pair_eta", "dr", "crossx_weight_trig_corr");
    }

    // --- Low-mass dimuon template-fit inputs (docs/tracking/low_mass_dimuon_template_fit.md,
    // Physics Procedure §2,§3a,§4): OS and SS minv over 0–4 GeV with the dimuon
    // mass window REMOVED (all other single-b cuts kept). Crossx-normalized (1/L_pp)
    // and reco+trig corrected (crossx_weight_trig_corr), IDENTICAL selection/binning/
    // weight for OS and SS so D_OS - D_SS is a clean combinatoric subtraction.
    {
        // signal_cuts MINUS the minv window; the gap cuts are kept IDENTICAL to signal_cuts
        // -- both muons' q*eta windows AND the pair-level |eta^pair| < 2.2, all read from
        // ParamsSet.
        const std::string signal_cuts_no_minv =
            ParamsSet::SignalPairPtCutExpr("pair_pt") + " && "
            + ParamsSet::FiducialGapCutExpr("m1.charge * m1.eta") + " && "
            + ParamsSet::FiducialGapCutExpr("m2.charge * m2.eta") + " && "
            + ParamsSet::PairFiducialEtaCutExpr("pair_eta");

        // Attach the trig+reco+crossx weight chain to a node filtered WITHOUT the
        // minv window. If FillHistogramsGeneric ran, it added w_trig/w_reco to
        // df_op/df_ss before any filter, so they propagate through this filter —
        // reuse them; else define inline exactly as the OS/SS crossx blocks above do.
        auto attach_crossx_weight = [&](ROOT::RDF::RNode node) -> ROOT::RDF::RNode {
            ROOT::RDF::RNode n = generic_weight_col.empty()
                ? AddPairEfficiencyWeightColumns(node)
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
            ROOT::RDF::TH1DModel("h1d_crossx_minv_0_4_op_dsigma", ";m_{#mu#mu} [GeV];d#sigma/dm_{#mu#mu} [pb GeV^{-1}]", 50, 0.0, 4.0),
            "minv", "crossx_weight_trig_corr");
        hist1d_rresultptr_map["h1d_crossx_minv_0_4_ss_dsigma"] = df_ss_no_minv.Histo1D(
            ROOT::RDF::TH1DModel("h1d_crossx_minv_0_4_ss_dsigma", ";m_{#mu#mu} [GeV];d#sigma/dm_{#mu#mu} [pb GeV^{-1}]", 50, 0.0, 4.0),
            "minv", "crossx_weight_trig_corr");
    }

    std::cout << "[PP] FillHistogramsCrossx completed" << std::endl;
}
