#pragma once

#include "RDFBasedHistFillingBaseClass.cxx"
#include "../Utilities/MC_helpers.h"
#include "../Utilities/HistFillUtils.h"

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
    // virtual void        FillHistogramsTruth() override();

public:
    explicit RDFBasedHistFillingPowhegTruth()
    : RDFBasedHistFillingPowheg(false, false){}
    ~RDFBasedHistFillingPowhegTruth(){}
};

class RDFBasedHistFillingPowhegFullsim : public virtual RDFBasedHistFillingPowheg{
protected:

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

// --------------------- protected class methods ---------------------------

    void                InitializePowhegFullsimExtra();
    virtual void        InitializePowhegExtra() override{ return InitializePowhegFullsimExtra(); }
    
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

public:
    explicit RDFBasedHistFillingPowhegFullsim(int run_year_input = 17)
    : RDFBasedHistFillingPowheg(true, false, run_year_input){}
    ~RDFBasedHistFillingPowhegFullsim(){}

};

class RDFBasedHistFillingPowhegFullsimOverlay : public RDFBasedHistFillingPowhegFullsim{
public:
    explicit RDFBasedHistFillingPowhegFullsimOverlay(int run_year_input)
    : RDFBasedHistFillingPowheg(true, true, run_year_input){}
    ~RDFBasedHistFillingPowhegFullsimOverlay(){}

};