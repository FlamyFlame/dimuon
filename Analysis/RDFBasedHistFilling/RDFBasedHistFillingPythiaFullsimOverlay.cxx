#include "RDFBasedHistFillingPythiaFullsim.cxx"
#include "../MuonObjectsParamsAndHelpers/PbPbBaseClass.h"
#include "../MuonObjectsParamsAndHelpers/FullSimSampleType.h"
#include "../Utilities/GeneralUtils.h"
#include <cmath>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>

class RDFBasedHistFillingPythiaFullsimOverlay
  : public virtual RDFBasedHistFillingPythiaFullsim
  , public PbPbBaseClass<RDFBasedHistFillingPythiaFullsimOverlay>
{
    friend class PbPbBaseClass<RDFBasedHistFillingPythiaFullsimOverlay>;
    int RunYear() const { return 24; }

protected:
    void SetIOPathsHook() override {
        const std::string data_dir = FullSimSampleInputDir(fullsim_sample_type);
        const std::string label    = FullSimSampleLabel(fullsim_sample_type);
        const std::string cut_suffix = with_data_resonance_cuts
            ? "_with_data_resonance_cuts"
            : "_no_data_resonance_cuts";

        input_files.clear();
        input_files.push_back(data_dir + "muon_pairs_pythia_fullsim_" + label + cut_suffix + ".root");
        output_file = data_dir + "histograms_pythia_fullsim_" + label + cut_suffix + ".root";
        infile_var1D_json = "var1D_pythia_fullsim.json";
    }

    void InitializePythiaExtra() override {
        RDFBasedHistFillingPythiaFullsim::InitializePythiaFullsimExtra();
        InitializePbPb();
    }

    void CreateBaseRDFsPythiaExtra() override {
        RDFBasedHistFillingPythiaFullsim::CreateBaseRDFsPythiaFullsimExtra();

        ROOT::RDF::RNode& df_op_weighted = map_at_checked(df_map, "df_op_weighted",
            "CreateBaseRDFsPythiaExtra (overlay): df_op_weighted");
        ROOT::RDF::RNode& df_ss_weighted = map_at_checked(df_map, "df_ss_weighted",
            "CreateBaseRDFsPythiaExtra (overlay): df_ss_weighted");

        for (int ictr = 0; ictr < nCtrBins; ictr++) {
            int lo = ctr_bin_edges[ictr];
            int hi = ctr_bin_edges[ictr + 1];
            const std::string& ctr_tag = ctr_bins[ictr];

            auto ctr_filter = [lo, hi](int ctr) { return ctr >= lo && ctr < hi; };

            df_map.emplace("df_op" + ctr_tag + "_weighted",
                df_op_weighted.Filter(ctr_filter, {"avg_centrality"}));
            df_map.emplace("df_ss" + ctr_tag + "_weighted",
                df_ss_weighted.Filter(ctr_filter, {"avg_centrality"}));

            ROOT::RDF::RNode& df_sb = map_at_checked(df_map, "df_single_b_weighted",
                "CreateBaseRDFsPythiaExtra (overlay): df_single_b_weighted");
            df_map.emplace("df_single_b" + ctr_tag + "_weighted",
                df_sb.Filter(ctr_filter, {"avg_centrality"}));

            for (const std::string& pair_catgr : {"_ss", "_op", "_single_b"}) {
                const std::string base = "df" + pair_catgr + ctr_tag;
                ROOT::RDF::RNode& node = map_at_checked(df_map, base + "_weighted",
                    Form("CreateBaseRDFsPythiaExtra (overlay): %s_weighted", base.c_str()));

                auto node_sig = node
                    .Define("pass_signal_truth" + ctr_tag,
                        "truth_minv > 1.08 && truth_minv < 2.9 && truth_pair_pt > 8 && truth_pair_eta < 2.2 && truth_dr > 0.05")
                    .Define("pass_signal_reco" + ctr_tag,
                        "(m1.reco_match && m2.reco_match) ? (minv > 1.08 && minv < 2.9 && pair_pt > 8 && pair_eta < 2.2 && dr > 0.05) : false");

                df_map.emplace(base + "_pass_medium_weighted",
                    node_sig.Filter("pair_pass_medium"));
                df_map.emplace(base + "_pass_tight_weighted",
                    node_sig.Filter("pair_pass_tight"));
                df_map.emplace(base + "_pass_signal_truth_weighted",
                    node_sig.Filter("pass_signal_truth" + ctr_tag));
                df_map.emplace(base + "_pass_medium_and_signal_truth_and_reco_weighted",
                    node_sig.Filter("pair_pass_medium && pass_signal_truth" + ctr_tag + " && pass_signal_reco" + ctr_tag));
                df_map.emplace(base + "_pass_tight_and_signal_truth_and_reco_weighted",
                    node_sig.Filter("pair_pass_tight && pass_signal_truth" + ctr_tag + " && pass_signal_reco" + ctr_tag));
            }
        }
    }

    void FillHistogramsFullSim() override {
        RDFBasedHistFillingPythiaFullsim::FillHistogramsFullSim();

        for (int ictr = 0; ictr < nCtrBins; ictr++) {
            const std::string& ctr_tag = ctr_bins[ictr];

            for (const std::string& pair_catgr : {"_ss", "_op", "_single_b"}) {
                const std::vector<std::string> quality_cats = {
                    "", "_pass_medium", "_pass_tight",
                    "_pass_signal_truth",
                    "_pass_medium_and_signal_truth_and_reco",
                    "_pass_tight_and_signal_truth_and_reco"
                };

                for (const std::string& qc : quality_cats) {
                    const std::string filter = pair_catgr + ctr_tag + qc;
                    const std::string df_name = "df" + filter + "_weighted";
                    auto it = df_map.find(df_name);
                    if (it == df_map.end()) continue;
                    FillHistogramsSingleDataFrame(filter, "", it->second);
                }
            }
        }
    }

    void HistPostProcessExtra() override {
        RDFBasedHistFillingPythiaFullsim::MakeAndWriteMuPairRecoEffProjGraphs();

        for (int ictr = 0; ictr < nCtrBins; ictr++) {
            const std::string& ctr_tag = ctr_bins[ictr];
            MakeAndWriteMuPairRecoEffProjGraphsHelper(
                {"_ss" + ctr_tag, "_op" + ctr_tag, "_single_b" + ctr_tag}, true);
            MakeAndWriteMuPairRecoEffProjGraphsHelper(
                {"_ss" + ctr_tag, "_op" + ctr_tag, "_single_b" + ctr_tag}, true, true);
        }
    }

    void BuildFilterToVarListMapExtra() override {
        RDFBasedHistFillingPythiaFullsim::BuildFilterToVarListMapExtra();

        for (int ictr = 0; ictr < nCtrBins; ictr++) {
            const std::string& ctr_tag = ctr_bins[ictr];
            std::vector<std::string> ctr_reco_effcy_filters;
            std::vector<std::string> ctr_detector_response_filters;

            for (const std::string& pair_catgr : {"_ss", "_op", "_single_b"}) {
                const std::vector<std::string> quality_cats_eff = {
                    "", "_pass_medium", "_pass_tight",
                    "_pass_signal_truth",
                    "_pass_medium_and_signal_truth_and_reco",
                    "_pass_tight_and_signal_truth_and_reco"
                };
                for (const std::string& qc : quality_cats_eff)
                    ctr_reco_effcy_filters.push_back(pair_catgr + ctr_tag + qc);
            }

            for (const std::string& pair_catgr : {"_single_b"}) {
                for (const std::string& qc : {"_pass_medium", "_pass_tight"})
                    ctr_detector_response_filters.push_back(pair_catgr + ctr_tag + qc);
            }

            for (auto& f : ctr_reco_effcy_filters)
                InsertOrAppend(df_filter_and_weight_to_var1D_list_map, std::make_pair(f, ""), reco_effcy_var1Ds);
            for (auto& f : ctr_reco_effcy_filters)
                InsertOrAppend(df_filter_and_weight_to_var2D_list_map, std::make_pair(f, ""), reco_effcy_var2Ds);
            for (auto& f : ctr_reco_effcy_filters)
                InsertOrAppend(df_filter_and_weight_to_var3D_list_map, std::make_pair(f, ""), reco_effcy_var3Ds);

            for (auto& f : ctr_detector_response_filters)
                InsertOrAppend(df_filter_and_weight_to_var1D_list_map, std::make_pair(f, ""), detec_resp_var1Ds);
            for (auto& f : ctr_detector_response_filters)
                InsertOrAppend(df_filter_and_weight_to_var2D_list_map, std::make_pair(f, ""), detec_resp_var2Ds);
        }
    }

public:
    FullSimSampleType fullsim_sample_type = FullSimSampleType::hijing;

    explicit RDFBasedHistFillingPythiaFullsimOverlay(
        FullSimSampleType sample_type = FullSimSampleType::hijing)
        : fullsim_sample_type(sample_type)
    {
        is_fullsim = true;
        is_fullsim_overlay = true;
    }
    ~RDFBasedHistFillingPythiaFullsimOverlay(){}
};
