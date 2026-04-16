#pragma once

#include "RDFBasedHistFillingBaseClass.cxx"
#include "../MuonObjectsParamsAndHelpers/muon_pair_enums_MC.h"
#include "../Utilities/HistFillUtils.h"
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <vector>

class RDFBasedHistFillingPythia : public virtual RDFBasedHistFillingBaseClass{
protected:
    bool isPrivate{false};
    bool with_data_resonance_cuts{false};
    double E_COM{5.36};

    bool is_fullsim{false};
    bool is_fullsim_overlay{false};

    void            InitializePythiaImpl();
    void            InitializePythiaCommon();
    virtual void    InitializePythiaExtra(){}

    virtual void    InitAnalysisSettingsHook() override{ return InitializePythiaImpl(); }
    virtual void    SetIOPathsHook() override;

    void            CreateBaseRDFsPythiaImpl();
    void            CreateBaseRDFsPythiaCommon();
    virtual void    CreateBaseRDFsPythiaExtra(){}
    virtual void    CreateBaseRDFsExtra() override{ return CreateBaseRDFsPythiaImpl(); }

    virtual void    FillHistograms() override;
    virtual void    FillHistogramsTruth(){}
    virtual void    FillHistogramsFullSim(){} // Placeholder only for now

public:
    explicit RDFBasedHistFillingPythia(bool isPrivate_input = false,
                                       double ecom_input = 5.36,
                                       bool with_data_resonance_cuts_input = false)
    : isPrivate(isPrivate_input),
      with_data_resonance_cuts(with_data_resonance_cuts_input),
      E_COM(ecom_input)
    {}

    ~RDFBasedHistFillingPythia(){}

    void SetIsPrivate(bool v){ isPrivate = v; }
    bool GetIsPrivate() const { return isPrivate; }

    void SetECom(double v){ E_COM = v; }
    double GetECom() const { return E_COM; }

    void SetWithDataResonanceCuts(bool v){ with_data_resonance_cuts = v; }
    bool GetWithDataResonanceCuts() const { return with_data_resonance_cuts; }
};

class RDFBasedHistFillingPythiaTruth : public virtual RDFBasedHistFillingPythia{
protected:
    std::vector<std::string>                vars1D_general;
    std::vector<std::string>                vars1D_general_over_dr;
    std::vector<std::string>                vars1D_general_over_pair_pt;
    std::vector<std::array<std::string, 2>> vars2D_general;
    std::vector<std::array<std::string, 2>> vars2D_general_over_dr;

    std::vector<std::string>                vars1D_flavor_origin;
    std::vector<std::string>                vars1D_flavor_origin_over_dr;
    std::vector<std::string>                vars1D_flavor_origin_over_pair_pt;
    std::vector<std::array<std::string, 2>> vars2D_flavor_origin;
    std::vector<std::array<std::string, 2>> vars2D_flavor_origin_over_dr;

    std::map<int, std::string> flavor_suffix_map;
    std::map<int, std::string> origin_suffix_map;

    virtual void    InitializePythiaExtra() override;
    virtual void    BuildHistBinningMapExtra() override;
    virtual void    CreateBaseRDFsPythiaExtra() override;
    virtual void    FillHistogramsTruth() override;
    virtual void    HistPostProcessExtra() override;
    virtual void    WriteOutputExtra() override;

    void FillHistogramsGeneral();
    void FillHistogramsFlavorBinned();
    void FillHistogramsOriginBinned();
    void FillHistogramsResonanceStudy();
    void FillHistogramsCrossxAndSpecialEta();
    void FillHistogramsSignalAcceptance();
    void ValidatePythiaTruthSchemaAndCoverage();
    void BuildAndStoreNearAwaySummedHistograms();
    void MarkNearAwayHistogramsForNominalExclusion();
    bool IsNearAwayDividedHistName(const std::string& hname) const;
    std::string NearAwayToSummedName(const std::string& hname) const;

public:
    explicit RDFBasedHistFillingPythiaTruth(bool isPrivate_input = false,
                                            double ecom_input = 5.36,
                                            bool with_data_resonance_cuts_input = false)
    : RDFBasedHistFillingPythia(isPrivate_input, ecom_input, with_data_resonance_cuts_input)
    {}

    ~RDFBasedHistFillingPythiaTruth(){}
};

class RDFBasedHistFillingPythiaFullsim : public virtual RDFBasedHistFillingPythia {
protected:
    // --- filter levels ---
    std::vector<std::vector<std::string>> levels_reco_effcy_filters;
    std::vector<std::vector<std::string>> levels_detector_response_filters;

    std::vector<std::string> reco_effcy_filters;
    std::vector<std::string> detector_response_filters;

    std::vector<std::string>               reco_effcy_var1Ds;
    std::vector<std::array<std::string,2>> reco_effcy_var2Ds;
    std::vector<std::array<std::string,3>> reco_effcy_var3Ds;

    std::vector<std::string>               detec_resp_var1Ds;
    std::vector<std::array<std::string,2>> detec_resp_var2Ds;

    // --- projection-graph configs ---
    std::vector<float> dr_bins_edges_for_reco_effcy          = {0.0f, 0.2f, 0.4f, 0.6f, 1.0f};
    std::vector<float> pair_pT_bins_edges_for_reco_effcy_dR  = {8.0f, 12.0f, 20.0f, std::numeric_limits<float>::max()};

    std::vector<std::pair<float,float>> dr_ranges_for_reco_effcy;
    std::vector<std::pair<float,float>> pair_pT_ranges_for_reco_effcy_dR;

    std::vector<std::tuple<std::string, std::string, bool, const std::vector<std::pair<float,float>>*>>
        mu_pair_reco_eff_proj_cfgs;

    std::vector<std::pair<std::string,std::string>> reco_eff_num_denom_suffix_pairs;

    std::map<std::string, TGraphAsymmErrors*> mu_pair_reco_eff_proj_graph_map;
    std::vector<std::string>                  mu_pair_reco_eff_proj_graphs_to_not_write {};
    std::map<std::string, TH1D*>              mu_pair_reco_eff_proj_hist_map;
    std::vector<std::string>                  mu_pair_reco_eff_proj_hists_to_not_write {};

// ----- hooks -----

    void             InitializePythiaFullsimExtra();
    virtual void     InitializePythiaExtra() override { return InitializePythiaFullsimExtra(); }

    virtual void     SetIOPathsHook() override;

    void             BuildHistBinningMapPythiaFullsimExtra();
    virtual void     BuildHistBinningMapExtra() override { return BuildHistBinningMapPythiaFullsimExtra(); }

    virtual void     BuildFilterToVarListMapExtra() override;

    void             CreateBaseRDFsPythiaFullsimExtra();
    virtual void     CreateBaseRDFsPythiaExtra() override { return CreateBaseRDFsPythiaFullsimExtra(); }

    virtual void     FillHistogramsFullSim() override;
    virtual void     FillHistogramsFullSimDetecResp();
    virtual void     FillHistogramsFullSimRecoEffcies();

    virtual void     HistPostProcessExtra() override { return MakeAndWriteMuPairRecoEffProjGraphs(); }
    virtual void     WriteOutputExtra() override;
    virtual void     CleanupExtra() override;

    void             MakeAndWriteMuPairRecoEffProjGraphsHelper(
                         const std::vector<std::string>& categories,
                         bool use_TH_divide = true,
                         bool require_signal_cuts = false);
    void             MakeAndWriteMuPairRecoEffProjGraphs();

public:
    explicit RDFBasedHistFillingPythiaFullsim() { is_fullsim = true; }
    ~RDFBasedHistFillingPythiaFullsim(){}
};
