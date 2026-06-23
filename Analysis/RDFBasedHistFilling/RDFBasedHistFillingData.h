#pragma once
#include "../MuonObjectsParamsAndHelpers/PbPbBaseClass.h"
#include "../MuonObjectsParamsAndHelpers/PPBaseClass.h"
#include "../MuonObjectsParamsAndHelpers/DatasetTriggerMap.h"
#include "../Utilities/HistFillUtils.h"
#include "CommonEffcyConfig.h"
#include "CorrectionStages.h"
#include "RDFBasedHistFillingBaseClass.cxx"

class RDFBasedHistFillingData : public virtual RDFBasedHistFillingBaseClass{
protected:
// --------------------- class variables ---------------------------

    int run_year;
    bool isForSoumya;

    std::string isForSoumya_suffix;
    std::string qEtaBin_suffix;
    std::string trig_suffix;       // includes _no_trg_plots suffix when doTrigEffcy=false; used for output path
    std::string base_trig_suffix;  // trigger type only, no _no_trg_plots; used for input ntuple paths
    std::string input_mindR_suffix; // e.g. "_mindR_0_02" appended between trig_suffix and _res_cut_v2 in input paths; "" for old files

    std::string out_file_suffix;

    bool trigger_effcy_calc = false;  // derived in TriggerModeSettings(): true iff trigger_mode ∈ {0,1} AND NOT (isPbPb && isRun3 && mu4_nominal_pbpb_NO_trig_calc)
    bool use_mu6_for_trg_eff = false; // derived: true when trigger_effcy_calc && trigger_mode == 1

    std::vector<std::string> categories_essential;
    std::vector<std::string> trigs;
    std::map<std::string, std::string> trig_to_filter_str_map;
    std::vector<std::pair<std::string, std::string>> trigs_pair;

    std::vector<std::string> trg_effcy_biases;

    std::map<std::string, TGraphAsymmErrors*> graph_map;
    std::vector<std::string> graphs_to_not_write {};

    std::vector<std::string> single_muon_trig_effcy_var1Ds = {"Deta", "Deta_zoomin", "Dphi", "Dphi_zoomin", "DR", "DR_zoomin", "DR_0_2", "minv_zoomin", "pair_pt_log"};
    std::vector<std::array<std::string, 2>> single_muon_trig_effcy_var2Ds = {{"phi2nd","pt2nd"}, {"q_eta2nd","pt2nd"}, {"q_eta2nd","phi2nd"}, {"pt2nd", "DR_zoomin"}, {"pair_pt_log", "DR_zoomin"}};
    std::vector<std::array<std::string, 3>> single_muon_trig_effcy_var3Ds = {{"phi2nd","q_eta2nd","pt2nd"}};

    // --- levels of trigger efficiency filter to be summed, pre- & post sum ---
    std::vector<std::vector<std::string>> levels_trg_effcy_filters_1D_pre_sum;
    std::vector<std::vector<std::string>> levels_trg_effcy_filters_2D_3D_pre_sum;

    std::vector<std::string> trg_effcy_filters_1D_pre_sum;
    std::vector<std::string> trg_effcy_filters_2D_3D_pre_sum;

    std::vector<int> levels_trg_effcy_to_be_summed = {0,1}; // HARD-CODED for now: arbitrary leveling is very complicated; for Pb+Pb, meaning mu1/mu2 pass mu4 MUST precede centrality binning
    std::vector<std::vector<std::string>> levels_trg_effcy_filters_to_be_summed;
    std::vector<std::string> trg_effcy_filters_to_be_summed;

    std::vector<std::vector<std::string>> levels_trg_effcy_filters_1D_post_sum;
    std::vector<std::vector<std::string>> levels_trg_effcy_filters_2D_3D_post_sum;

    std::vector<std::string> trg_effcy_filters_1D_post_sum;
    std::vector<std::string> trg_effcy_filters_2D_3D_post_sum;

// --------------------- kinematics configuration ---------------------------

    using QEtaBinning = std::vector<std::pair<float, float>>;
    QEtaBinning q_eta_proj_ranges;

    std::vector<std::string> q_eta_ranges_str;

    std::vector<std::pair<float, float>> pair_pT_ranges_for_weighted_effcy_dR_fitting = {
        {8, 12},
        {12, 24},
        {24, 80}
    };


// --------------------- protected class methods ---------------------------

    RDFBasedHistFillingData(int run_year_input, bool isForSoumya_input);

    virtual void    InitOutput() override;
    virtual void    InitAnalysisSettingsHook() override{
        return InitializeDataImpl();
    }

    void    InitializeDataImpl(); // initializations common for all data types
    void    InitializeDataCommon(); // initializations common for all data types
    virtual void    InitializeDataExtra(){} // initializations for specific child-class data types
    virtual void    TriggerModeSettings();

    void            BuildHistBinningMapDataImpl();
    void            BuildHistBinningMapDataCommon();
    virtual void    BuildHistBinningMapDataExtraHook(){}
    virtual void    BuildHistBinningMapExtra() override{
        return BuildHistBinningMapDataImpl();
    }

    void            BuildFilterToVarListMapDataImpl();
    virtual void    BuildFilterToVarListMapExtra() override{ 
        return BuildFilterToVarListMapDataImpl();
    }
    void            BuildFilterToVarListMapDataCommon();
    void            BuildTrgEffcyFilterToVarListMap();
    virtual void    BuildFilterToVarListMapDataExtraHook(){}

    void            FlattenTrigEffcyFilters();
    void            FlattenTrigEffcyFiltersDataCommon();
    virtual void    FlattenTrigEffcyFiltersExtra();

    void            BuildFlattenedTrgEffcyFilterToVarListMap();
    void            BuildFlattenedTrgEffcyFilterToVarListMapDataCommon();
    virtual void    BuildFlattenedTrgEffcyFilterToVarListMapExtra(){}

    virtual bool    IsPbPb() const { return false; } // avoids dynamic_cast<PbPb*> cross-dependency when compiling PP and PbPb separately
    std::string     generic_weight_col;
    virtual void    FillHistograms() override;
    virtual void    FillHistogramsGeneric();
    virtual void    FillHistogramsCrossx() = 0; // trigger_mode == 2 crossx filling (opposite-sign only, with signal cuts)
    virtual void    FillHistogramsSingleMuonEffcy() = 0;
    virtual void    FillHistogramsDimuTrigGivenMu4() = 0;
    virtual void    FillHistogramsMu4GivenMB() = 0;
    virtual void    FillTrigEffcyHistsInvWeightedbySingleMuonEffcies() = 0;
    // bool            PassSingleMuonGapCut(float meta, float mpt, int mcharge);

    virtual void    OpenEffcyPtFitFile() = 0;

    void            HistPostProcessDataImpl();
    virtual void    HistPostProcessExtra() override{
        HistPostProcessDataImpl();
    }
    void            HistPostProcessDataCommon();
    virtual void    HistPostProcessDataExtra(){}

    void            SumSingleMuonTrigEffHists();
    void            SumSingleMuonTrigEffHistsDataCommon();
    virtual void    SumSingleMuonTrigEffHistsExtra(){}

    virtual void    CalculateSingleMuonTrigEffcyRatios() = 0;
    void            CalculateSingleMuonTrigEffcyRatiosHelper(const std::vector<std::string>& categories);
    
    virtual void    MakeAndWriteSingleMuonTrigEffPtGraphs() = 0;
    void            MakeAndWriteSingleMuonTrigEffPtGraphsHelper(const std::vector<std::string>& categories);
    virtual void    MakeAndWriteDRTrigEffGraphs() = 0;
    void            MakeAndWriteDRTrigEffGraphsHelper(const std::vector<std::string>& categories);

    static std::string      FindBinReturnStr(float number, const std::vector<std::pair<float, float>>& ranges);
    static std::string      FindCtrSuffix(int centrality);
    static float            EvaluateSingleMuonEffcyPtFitted(const std::string& ctr_suffix, bool charge_positive, float pt, float q_eta);
    static float            EvaluateSingleMuonEffcy(const std::string& ctr_suffix, bool charge_positive, float pt, float q_eta);

    // --- Run 2 single-muon reconstruction-efficiency PLACEHOLDER (eps1*eps2 proxy) ---
    // Loads EfficiencyCorrs/EffFiles/run2_reco_eff_placeholder.root. PbPb: the
    // colleague's Run 2 Medium-mu reco-eff logistic TF1 fits (eps_reco vs pT) per
    // q*eta slice and centrality, evaluated at the exact muon pT. pp: barrel/endcap
    // TGraphs (HF R_AA Fig.31). This is a temporary stand-in until the proper 3D
    // pair efficiency exists (task_05).
    // See docs/tracking/reco_eff_placeholder_run2.md.
    static void             OpenRecoEffPlaceholderFile();
    // centrality < 0 => pp (barrel/endcap by |q_eta|); otherwise PbPb (F.2 ctr interval).
    static float            EvaluateSingleMuonRecoEffPlaceholder(int centrality, float pt, float q_eta);

    virtual void     WriteOutputExtra() override;
    virtual void     CleanupExtra() override;

public:
enum HistFillingCycle{
    generic = 1,
    inv_weight_by_single_mu_effcy = 2,
    inv_weight_by_dR_effcy_corr = 3
};
    // Public method to run histogram filling
    void RunFillHistograms() { FillHistograms(); }

    bool useCoarseQEtaBin = false;
    int hist_filling_cycle = generic;
    int trigger_mode = 1;
    bool useMu4NoL1Leg = true;  // if true, require the probe muon to pass the mu4noL1 (unseeded) leg for _mu4_mu4noL1 filters
    double mindR_trig = 0.02;   // > 0: search input files with _mindR_X_XX suffix; <= 0: use old input files (no suffix)
    bool isScram = false;
    bool isTight = false;
    bool doTrigEffcy = true; // default true; derived from trigger_effcy_calc in TriggerModeSettings()
    bool mu4_nominal_pbpb_NO_trig_calc = false; // PbPb nominal pipeline: use mu4 trigger for event selection only (no trig effcy derivation)
    // PUBLIC (set from run_template_fit_* macros). Low-mass dimuon TEMPLATE-FIT pass: read the
    // _no_res_cut ntuples (resonances PRESENT) and produce ONLY the 0-4 GeV minv template spectra
    // D_OS/D_SS (1D + 2D vs pair pT/eta), to a DISTINCT output (out_file_suffix += "_template_fit").
    // Use the SAME trigger_mode/mu4_nominal flags as nominal crossx. MUST stay public so the cling
    // macro assignment compiles (a protected member silently fails -> runs in nominal mode and
    // overwrites the nominal output; cf. the 2026-06-22 output_single_muon_tree incident). Not to be
    // combined with trigger_effcy_calc. See docs/tracking/low_mass_dimuon_template_fit.md.
    bool low_mass_template_calc = false;

    bool save_non_sepr_trg_hists = false; // save trigger efficiency histograms without separation requirement
    bool save_good_accept_trg_hists = false; // save trigger efficiency histograms with good acceptance requirement

    bool output_generic_hists;
    bool output_gapcut_hists;

    bool filter_out_photo_resn_for_trig_effcy = true;
    bool use_3D_2nd_muon = false; // if true, use 3D kinematics (phi, q*eta, pT) for single (2nd) muon trigger efficiencies
    bool use_pT_fitting_single_muon_effcy = true;


    ~RDFBasedHistFillingData(){}
    
};


// ----------------------------------------------------------------------------------------------------------------
class RDFBasedHistFillingPP : public RDFBasedHistFillingData{

protected:
// --------------------- protected class variables ---------------------------
    double pp_crossx_lumi_factor = -1.; // set in InitializePPExtra(); -1 = not set / not applicable

    std::vector<int> levels_trg_effcy_to_be_summed_w_musign_summing;
    std::vector<std::vector<std::string>> levels_trg_effcy_filters_to_be_summed_w_musign_summing;
    std::vector<std::string> trg_effcy_filters_to_be_summed_w_musign_summing;
    
    std::vector<std::vector<std::string>> levels_trg_effcy_filters_1D_post_sum_w_musign_summing;
    std::vector<std::vector<std::string>> levels_trg_effcy_filters_2D_3D_post_sum_w_musign_summing;

    std::vector<std::string> trg_effcy_filters_1D_post_sum_w_musign_summing;
    std::vector<std::string> trg_effcy_filters_2D_3D_post_sum_w_musign_summing;

// --------------------- protected class methods ---------------------------

    virtual void        SetIOPathsHook() override;
    void                InitializePPExtra();
    virtual void        InitializeDataExtra() override{ return InitializePPExtra(); }
    virtual void        FlattenTrigEffcyFiltersExtra() override;
    
    virtual void        FillHistogramsSingleMuonEffcy() override;
    virtual void        FillHistogramsDimuTrigGivenMu4() override;
    virtual void        FillHistogramsMu4GivenMB() override;
    virtual void        FillTrigEffcyHistsInvWeightedbySingleMuonEffcies() override;
    virtual void        FillHistogramsCrossx() override;
    virtual void        FillHistogramsGeneric() override;

    virtual void        OpenEffcyPtFitFile() override;

    void                SumSingleMuonTrigEffHistsPP();
    virtual void        SumSingleMuonTrigEffHistsExtra() override{ return SumSingleMuonTrigEffHistsPP();}

    virtual void        CalculateSingleMuonTrigEffcyRatios() override;
    virtual void        MakeAndWriteSingleMuonTrigEffPtGraphs() override;
    void                MakeAndWriteDRTrigEffGraphs() override;

public:

    int RunYear() const { return run_year; }

    explicit RDFBasedHistFillingPP(int run_year_input, bool isForSoumya_input = false)
    : RDFBasedHistFillingData (run_year_input, isForSoumya_input){
        std::cout << "constructor for pp called" << std::endl; 
    }
    ~RDFBasedHistFillingPP(){
        std::cout << "destructor for pp called" << std::endl;
    }

};

// ----------------------------------------------------------------------------------------------------------------
class RDFBasedHistFillingPbPb : public RDFBasedHistFillingData, public PbPbBaseClass<RDFBasedHistFillingPbPb>{
protected:

// --------------------- protected class variables ---------------------------

    std::vector<std::string> ctr_musign_categories;
    
    std::vector<std::string>                var1D_list_no_ctr_rebin; // list of 1D variables for which we apply NO ctr rebinning
    std::map<std::string, std::vector<int>> var1D_excep_ctr_rebin_factor_map; // map of 1D variables to exception ctr rebinning

    std::vector<std::vector<std::string>> levels_trg_effcy_filters_ctr_dep_pre_sum;
    std::vector<std::string> trg_effcy_filters_ctr_dep_1D_pre_sum;
    std::vector<std::vector<std::string>> levels_trg_effcy_filters_ctr_dep_post_sum;
    std::vector<std::string> trg_effcy_filters_ctr_dep_1D_post_sum;
    
    std::vector<float> pT_bins_edges_for_trg_effcy_ctr_dep = {4, 6, 8, 12, 30}; // pT binning for centrality dependence of single-muon trigger efficiencies
    std::vector<std::string> pT_bins_for_trg_effcy_ctr_dep;
    
    //-------------- levels of trigger efficiency filter pre- & post sum, with mu sign & ctr summing --------------
    std::vector<std::vector<std::string>> levels_trg_effcy_filters_1D_pre_sum_w_musign_ctr_summing;
    std::vector<std::vector<std::string>> levels_trg_effcy_filters_2D_3D_pre_sum_w_musign_ctr_summing;

    std::vector<std::string> trg_effcy_filters_1D_pre_sum_w_musign_ctr_summing;
    std::vector<std::string> trg_effcy_filters_2D_3D_pre_sum_w_musign_ctr_summing;

    std::vector<std::vector<std::string>> levels_trg_effcy_filters_1D_post_sum_w_musign_ctr_summing;
    std::vector<std::vector<std::string>> levels_trg_effcy_filters_2D_3D_post_sum_w_musign_ctr_summing;

    std::vector<std::string> trg_effcy_filters_1D_post_sum_w_musign_ctr_summing;
    std::vector<std::string> trg_effcy_filters_2D_3D_post_sum_w_musign_ctr_summing;

    //-------------- levels of trigger efficiency filter pre- & post sum, with ctr summing --------------
    std::vector<std::vector<std::string>> levels_trg_effcy_filters_1D_pre_sum_w_ctr_summing;
    std::vector<std::vector<std::string>> levels_trg_effcy_filters_2D_3D_pre_sum_w_ctr_summing;

    std::vector<std::string> trg_effcy_filters_1D_pre_sum_w_ctr_summing;
    std::vector<std::string> trg_effcy_filters_2D_3D_pre_sum_w_ctr_summing;

    std::vector<std::vector<std::string>> levels_trg_effcy_filters_1D_post_sum_w_ctr_summing;
    std::vector<std::vector<std::string>> levels_trg_effcy_filters_2D_3D_post_sum_w_ctr_summing;

    std::vector<std::string> trg_effcy_filters_1D_post_sum_w_ctr_summing;
    std::vector<std::string> trg_effcy_filters_2D_3D_post_sum_w_ctr_summing;

// --------------------- protected class methods ---------------------------

    bool                IsPbPb() const override { return true; }
    virtual void        SetIOPathsHook() override;
    void                InitializePbPbExtra();
    virtual void        InitializeDataExtra() override{ return InitializePbPbExtra(); }
    void                BuildHistBinningMapPbPbExtra();
    virtual void        BuildHistBinningMapDataExtraHook() override{ return BuildHistBinningMapPbPbExtra();}
    virtual void        BuildFlattenedTrgEffcyFilterToVarListMapExtra() override;
    virtual void        FlattenTrigEffcyFiltersExtra() override;

    void                CreateBaseRDFsPbPbExtra();
    virtual void        CreateBaseRDFsExtra() override{ return CreateBaseRDFsPbPbExtra(); }
    virtual void        FillHistogramsSingleMuonEffcy() override;
    virtual void        FillHistogramsDimuTrigGivenMu4() override;
    virtual void        FillHistogramsDimuTrigGivenMu4CtrDep();
    virtual void        FillHistogramsMu4GivenMB() override;
    virtual void        FillHistogramsMu4GivenMBCtrDep();
    virtual void        FillTrigEffcyHistsInvWeightedbySingleMuonEffcies() override;
    virtual void        FillHistogramsCrossx() override;

    virtual void        OpenEffcyPtFitFile() override;

    void                HistPostProcessPbPb();
    virtual void        HistPostProcessDataExtra() override{ return HistPostProcessPbPb();}
    void                SumSingleMuonTrigEffHistsPbPb();
    virtual void        SumSingleMuonTrigEffHistsExtra() override{ return SumSingleMuonTrigEffHistsPbPb();}

    virtual void        CalculateSingleMuonTrigEffcyRatios() override;
    virtual void        MakeAndWriteSingleMuonTrigEffPtGraphs() override;
    void                MakeAndWriteDRTrigEffGraphs() override;
    void                MakeAndWriteSingleMuonTrigEffCtrGraphs();

    AxisInfo            GetAxisInfo(const var1D& v, const std::string& filter) const override;
    void                ReadVar1DJson() override;

public:
    int RunYear() const { return run_year; }
    explicit RDFBasedHistFillingPbPb(int run_year_input, std::string ctr_binning_version_input = "default")
    : RDFBasedHistFillingData (run_year_input, false){
        isForSoumya = false;
        ctr_binning_version = ctr_binning_version_input;
        std::cout << "constructor for PbPb called, run year: " << run_year << ", isForSoumya: " << isForSoumya << ", ctr_binning_version: " << ctr_binning_version << std::endl; 
    }

    explicit RDFBasedHistFillingPbPb(int run_year_input, bool isForSoumya_input)
    : RDFBasedHistFillingData (run_year_input, isForSoumya_input){
        std::cout << "constructor for PbPb called, run year: " << run_year << ", isForSoumya: " << isForSoumya << ", ctr_binning_version: " << ctr_binning_version << std::endl; 
    }

    ~RDFBasedHistFillingPbPb(){
        std::cout << "destructor for PbPb called" << std::endl;
    }
};
