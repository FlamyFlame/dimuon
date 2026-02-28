#pragma once

#include "RDFBasedHistFillingBaseClass.cxx"
#include "../Utilities/MC_helpers.h"
#include "../Utilities/HistFillUtils.h"
#include <limits>
#include <tuple>

class RDFBasedHistFillingPowheg : public virtual RDFBasedHistFillingBaseClass{
protected:
// --------------------- protected class variables ---------------------------

    int     run_year;
    bool    isRun3; // auto-determined from run_year
    bool    is_fullsim{false}; // if is fullsim, load reco quantities
    bool    is_fullsim_overlay{false}; // if is fullsim overlay, load centrality info
    bool    perform_truth{false};

// --------------------- protected class methods ---------------------------

    // ----- initialize powheg-specific settings -----
    
    virtual void    SetIOPathsHook() override;
    void            InitializePowhegImpl(); // initializations common for all Powheg types
    void            InitializePowhegCommon(); // initializations common for all Powheg types
    virtual void    InitializePowhegExtra(){} // initializations for specific child-class Powheg types

    virtual void    InitAnalysisSettingsHook() override{
        return InitializePowhegImpl();
    }

    void            BuildFilterToVarListMapPowhegImpl();
    virtual void    BuildFilterToVarListMapExtra() override{ 
        return BuildFilterToVarListMapPowhegImpl();
    }

    void            BuildSimpleFilterToVarListMapPowhegCommon(){}
    virtual void    BuildSimpleFilterToVarListMapPowhegExtra(){}

    void            FlattenFilters();
    void            FlattenFiltersPowhegCommon(){}
    virtual void    FlattenFiltersExtra(){}

    void            BuildFlattenedFilterToVarListMap();
    void            BuildFlattenedFilterToVarListMapPowhegCommon(){}
    virtual void    BuildFlattenedFilterToVarListMapExtra(){}

    // ----- fill histograms -----

    void            CreateBaseRDFsPowhegImpl();
    void            CreateBaseRDFsPowhegCommon();
    virtual void    CreateBaseRDFsPowhegExtra(){}
    virtual void    CreateBaseRDFsExtra() override{ return CreateBaseRDFsPowhegImpl(); }

    virtual void    FillHistograms() override;
    virtual void    FillHistogramsTruth(){}
    virtual void    FillHistogramsFullSim(){}

    explicit RDFBasedHistFillingPowheg(bool is_fullsim_input, bool is_fullsim_overlay_input, int run_year_input = 0)
    : is_fullsim(is_fullsim_input), is_fullsim_overlay(is_fullsim_overlay_input), run_year(run_year_input % 2000){
        if (is_fullsim && run_year < 15){
            throw std::runtime_error("For Powheg fullsim / fullsim overlay, Run year must be between 15 - 26! Current input invalid: " + std::to_string(run_year));
        }
        isRun3 = !(run_year <= 18);
    }

public:
    ~RDFBasedHistFillingPowheg(){}

};

class RDFBasedHistFillingPowhegTruth : public virtual RDFBasedHistFillingPowheg{
protected:
    std::vector<std::vector<std::string>> levels_truth_general_filters;
    std::vector<std::vector<std::string>> levels_truth_origin_filters;
    std::vector<std::vector<std::string>> levels_truth_flavor_filters;

    std::vector<std::string> truth_general_filters;
    std::vector<std::string> truth_origin_filters;
    std::vector<std::string> truth_flavor_filters;

    std::vector<std::string> truth_general_var1Ds;
    std::vector<std::array<std::string,2>> truth_general_var2Ds;
    std::vector<std::string> truth_general_var1Ds_jacobian;
    std::vector<std::array<std::string,2>> truth_general_var2Ds_jacobian;

    std::vector<std::string> truth_origin_var1Ds;
    std::vector<std::array<std::string,2>> truth_origin_var2Ds;
    std::vector<std::string> truth_origin_var1Ds_jacobian;
    std::vector<std::array<std::string,2>> truth_origin_var2Ds_jacobian;

    std::vector<std::string> truth_flavor_var1Ds;
    std::vector<std::array<std::string,2>> truth_flavor_var2Ds;
    std::vector<std::string> truth_flavor_var1Ds_jacobian;
    std::vector<std::array<std::string,2>> truth_flavor_var2Ds_jacobian;

    virtual void        InitializePowhegExtra() override;
    virtual void        BuildSimpleFilterToVarListMapPowhegExtra() override;
    virtual void        FlattenFiltersExtra() override;
    virtual void        BuildFlattenedFilterToVarListMapExtra() override;
    virtual void        CreateBaseRDFsPowhegExtra() override;

    virtual void        FillHistogramsTruth() override;
    void                FillHistogramsTruthGeneral();
    void                FillHistogramsTruthOriginBinned();
    void                FillHistogramsTruthFlavorBinned();

public:
    explicit RDFBasedHistFillingPowhegTruth()
    : RDFBasedHistFillingPowheg(false, false){}
    ~RDFBasedHistFillingPowhegTruth(){}
};

class RDFBasedHistFillingPowhegFullsim : public virtual RDFBasedHistFillingPowheg{
protected:
    bool useMixed{false};

    // --- levels of reco efficiency & detector response filters ---
    std::vector<std::vector<std::string>> levels_reco_effcy_filters;
    std::vector<std::vector<std::string>> levels_detector_response_filters;

    std::vector<std::string> reco_effcy_filters;
    std::vector<std::string> detector_response_filters;

    std::vector<std::string>                reco_effcy_var1Ds;
    std::vector<std::array<std::string,2>>  reco_effcy_var2Ds;
    std::vector<std::array<std::string,3>>  reco_effcy_var3Ds;

    std::vector<std::string>                detec_resp_var1Ds;
    std::vector<std::array<std::string,2>>  detec_resp_var2Ds;

    // --- mu-pair reco efficiency projection-graph configurations ---
    std::vector<float> dr_bins_edges_for_reco_effcy = {0.0f, 0.2f, 0.4f, 0.6f, 1.0f};
    std::vector<float> pair_pT_bins_edges_for_reco_effcy_dR = {8.0f, 12.0f, 20.0f, std::numeric_limits<float>::max()};

    std::vector<std::pair<float, float>> dr_ranges_for_reco_effcy;
    std::vector<std::pair<float, float>> pair_pT_ranges_for_reco_effcy_dR;

    // tuple: (varx, vary, project_y_axis, projection_ranges)
    std::vector<std::tuple<std::string, std::string, bool, const std::vector<std::pair<float, float>>*>>
        mu_pair_reco_eff_proj_cfgs;

    // {numerator, denominator}
    std::vector<std::pair<std::string, std::string>> reco_eff_num_denom_suffix_pairs;

    std::map<std::string, TGraphAsymmErrors*> mu_pair_reco_eff_proj_graph_map;
    std::vector<std::string> mu_pair_reco_eff_proj_graphs_to_not_write {};

// --------------------- protected class methods ---------------------------

    void                InitializePowhegFullsimExtra();
    virtual void        InitializePowhegExtra() override{ return InitializePowhegFullsimExtra(); }

    virtual void        SetIOPathsHook() override;
    
    void                BuildHistBinningMapPowhegFullsimExtra();
    virtual void        BuildHistBinningMapExtra() override{ return BuildHistBinningMapPowhegFullsimExtra();}

    void                FlattenFiltersPowhegFullsim();
    void                BuildFlattenedFilterToVarListMapPowhegFullsim();

    virtual void        FlattenFiltersExtra() override{ return FlattenFiltersPowhegFullsim(); }
    virtual void        BuildFlattenedFilterToVarListMapExtra() override{ return BuildFlattenedFilterToVarListMapPowhegFullsim(); }

    void                CreateBaseRDFsPowhegFullsimExtra();
    virtual void        CreateBaseRDFsPowhegExtra() override{ return CreateBaseRDFsPowhegFullsimExtra(); }

    virtual void        FillHistogramsFullSim() override;
    virtual void        FillHistogramsFullSimDetecResp();
    virtual void        FillHistogramsFullSimRecoEffcies();

    virtual void        HistPostProcessExtra() override{ return MakeAndWriteMuPairRecoEffProjGraphs(); }
    virtual void        WriteOutputExtra() override;
    virtual void        CleanupExtra() override;

    void                MakeAndWriteMuPairRecoEffProjGraphsHelper(const std::vector<std::string>& categories);
    void                MakeAndWriteMuPairRecoEffProjGraphs();

public:
    explicit RDFBasedHistFillingPowhegFullsim(int run_year_input = 17, bool useMixed_input = false)
    : RDFBasedHistFillingPowheg(true, false, run_year_input), useMixed(useMixed_input){}
    ~RDFBasedHistFillingPowhegFullsim(){}

    void SetUseMixed(bool useMixed_input){ useMixed = useMixed_input; }
    bool GetUseMixed() const { return useMixed; }

};

class RDFBasedHistFillingPowhegFullsimOverlay : public RDFBasedHistFillingPowhegFullsim{
public:
    explicit RDFBasedHistFillingPowhegFullsimOverlay(int run_year_input)
    : RDFBasedHistFillingPowheg(true, true, run_year_input){}
    ~RDFBasedHistFillingPowhegFullsimOverlay(){}

};