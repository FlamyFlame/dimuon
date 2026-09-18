#pragma once

#include "RDFBasedHistFillingBaseClass.cxx"
// Needed HERE, not only in the .cxx: the in-class initialiser of
// `dr_bins_edges_for_reco_effcy` below reads CommonEffcyConfig{}. Missing since the
// 2026-08-18 edit that introduced that initialiser, which left this header (and therefore
// RDFBasedHistFillingPythiaFullsim / ...Truth / ...FullsimOverlay) uncompilable with
// "use of undeclared identifier 'CommonEffcyConfig'". Found 2026-08-25.
#include "CommonEffcyConfig.h"
#include "../MuonObjectsParamsAndHelpers/muon_pair_enums_MC.h"
#include "../Utilities/HistFillUtils.h"
#include "../Utilities/PairRecoEffDefinition.h"
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
    void FillHistogramsTemplateMinvSignalRegion();
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
public:
    // Which fullsim production's NTuple-processing output to read. MUST match the isTestSample
    // used to produce it (PythiaAlgCoreT.h). TRUE today because only the TEST sample has NTP
    // output; flip to false once the FULL sample lands (it writes the "_full" suffix).
    // The TEST sample's 4-beam combination carries the Pb 4:6:6:9 isospin average, so any
    // cross-section from it is NOT a physical pp cross-section -- label it honestly.
    bool is_test_sample = true;

    // --- the detector-gap cuts, TRUTH and RECO leg, as ONE expression each -------------------
    // Single-muon fiducial windows on q*eta of BOTH muons + the pair-level |eta^pair| window, all
    // read from ParamsSet (never retyped). Declared here so the HIJING-overlay per-centrality
    // nodes and the pp24 class build the SAME string -- one definition, two call sites.
    // Since 2026-09-17 (docs/tracking/pair_reco_eff_gap_acceptance.md) the gap cuts are an
    // ACCEPTANCE folded into the pair reco efficiency: the RECO expression is part of the reco
    // numerator selection (= the data signal_cuts); the TRUTH expression is NOT part of the
    // truth denominator any more and survives only as the opt-in `_gapcut` truth family used by
    // the MC-vs-data shape comparison.
    static std::string TruthGapCutExpr() {
        return ParamsSet::FiducialGapCutExpr("m1.truth_charge * m1.truth_eta")
             + " && " + ParamsSet::FiducialGapCutExpr("m2.truth_charge * m2.truth_eta")
             + " && " + ParamsSet::PairFiducialEtaCutExpr("truth_pair_eta");
    }
    static std::string RecoGapCutExpr() {
        return ParamsSet::FiducialGapCutExpr("m1.charge * m1.eta")
             + " && " + ParamsSet::FiducialGapCutExpr("m2.charge * m2.eta")
             + " && " + ParamsSet::PairFiducialEtaCutExpr("pair_eta");
    }
    // The definition marker (Utilities/PairRecoEffDefinition.h) is written into the histogram
    // output by WriteOutputExtra; the builder refuses an input file without it.

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
    // READ from CommonEffcyConfig, the single source shared with the crossx plotter. Never
    // retype these edges (.claude/CLAUDE.md §Binnings).
    std::vector<float> dr_bins_edges_for_reco_effcy = CommonEffcyConfig{}.dr_bins_edges_for_reco_effcy;
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
