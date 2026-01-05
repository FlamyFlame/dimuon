#pragma once
#include "../MuonObjectsParamsAndHelpers/PbPbBaseClass.h"
#include "../MuonObjectsParamsAndHelpers/TrigEffcyUtils.h"
#include "RDFBasedHistFillingBaseClass.cxx"

class RDFBasedHistFillingData : public virtual RDFBasedHistFillingBaseClass{
protected:
// --------------------- class variables ---------------------------

    int run_year;
    bool isForSoumya;

    std::string isForSoumya_suffix;
    std::string qEtaBin_suffix;
    std::string trig_suffix;

    std::string out_file_suffix;

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

    // maps of q_eta bins to pT projection ranges serving single-muon effciency fitting
    QEtaBinning q_eta_proj_ranges_fine_excl_gap = {
        {-2.4f, -2.0f}, 
        {-2.0f, -1.6f}, 
        {-1.6f, -1.3f}, 
        {-0.9f, -0.5f}, 
        {-0.5f, -0.1f}, 
        {0.1f, 0.5f}, 
        {0.5f, 1.0f}, 
        {1.3f, 1.6f}, 
        {1.6f, 2.0f}, 
        {2.0f, 2.2f}
    };

    QEtaBinning q_eta_proj_ranges_fine_excl_gap_run2 = {
        {-2.4f, -2.0f}, 
        {-2.0f, -1.6f}, 
        {-1.6f, -1.3f}, 
        {-0.9f, -0.5f}, 
        {-0.5f, -0.1f}, 
        {0.1f, 0.5f}, 
        {0.5f, 1.0f}, 
        {1.3f, 1.6f}, 
        {1.6f, 2.0f}, 
        {2.0f, 2.4f}
    };

    QEtaBinning q_eta_proj_ranges_coarse_incl_gap = { // coarse bins including gaps
        {-2.4f, -2.0f}, 
        {-2.0f, -1.5}, 
        {-1.5, -1.0f}, 
        {-1.0f, -0.5f}, 
        {-0.5f, 0.5f}, 
        {0.5f, 1.0f}, 
        {1.0f, 1.5f}, 
        {1.5f, 2.0f}, 
        {2.0f, 2.2f}
    };

    QEtaBinning q_eta_proj_ranges_coarse_incl_gap_run2 = { // including [2.2, 2.4]
        {-2.4f, -2.0f}, 
        {-2.0f, -1.5}, 
        {-1.5, -1.0f}, 
        {-1.0f, -0.5f}, 
        {-0.5f, 0.5f}, 
        {0.5f, 1.0f}, 
        {1.0f, 1.5f}, 
        {1.5f, 2.0f}, 
        {2.0f, 2.4f}
    };

    std::vector<std::pair<float, float>> pair_pT_ranges_for_weighted_effcy_dR_fitting = {
        {8, 12},
        {12, 24},
        {24, 80}
    };


// --------------------- protected class methods ---------------------------

    RDFBasedHistFillingData(int run_year_input, bool isForSoumya_input);

    virtual void    InitOutput() override;
    virtual void    InitializeData(); // initializations common for all data types
    virtual void    TriggerModeSettings();

    virtual void    BuildHistBinningMap() override;

    virtual void    BuildFilterToVarListMap() override;
    virtual void    TrigEffcyFiltersPrePostSumFlattening();
    virtual void    BuildTrgEffcyFilterToVarListMap();

    virtual void    FillHistograms() override;
    virtual void    FillHistogramsGeneric();
    virtual void    FillHistogramsSingleMuonEffcy() = 0;
    virtual void    FillHistogramsDimuTrigGivenMu4() = 0;
    virtual void    FillHistogramsMu4GivenMB() = 0;
    virtual void    FillTrigEffcyHistsInvWeightedbySingleMuonEffcies() = 0;
    // bool            PassSingleMuonGapCut(float meta, float mpt, int mcharge);

    virtual void    OpenEffcyPtFitFile() = 0;

    virtual void    HistPostProcess() override;
    virtual void    SumSingleMuonTrigEffHists();

    virtual void    CalculateSingleMuonTrigEffcyRatios() = 0;
    void            CalculateSingleMuonTrigEffcyRatiosHelper(const std::vector<std::string>& categories);
    
    virtual void    MakeAndWriteSingleMuonTrigEffPtGraphs() = 0;
    void            MakeAndWriteSingleMuonTrigEffPtGraphsHelper(const std::vector<std::string>& categories);
    virtual void    MakeAndWriteDRTrigEffGraphs() = 0;
    void            MakeAndWriteDRTrigEffGraphsHelper(const std::vector<std::string>& categories);

    static std::string      FindBinReturnStr(float number, const std::vector<std::pair<float, float>>& ranges);
    static float            EvaluateSingleMuonEffcyPtFitted(bool charge_sign, std::string trg, float pt_2nd, float q_eta_2nd, float phi_2nd);
    static float            EvaluateSingleMuonEffcy(bool charge_sign, std::string trg, float pt_2nd, float q_eta_2nd, float phi_2nd);

    virtual void     Finalize() override;

public:
enum HistFillingCycle{
    generic = 1,
    inv_weight_by_single_mu_effcy = 2,
    inv_weight_by_dR_effcy_corr = 3
};

    bool useCoarseQEtaBin = true;
    int hist_filling_cycle = generic;
    int trigger_mode = 1;
    bool isScram = false;
    bool isTight = false;
    bool doTrigEffcy = true; // default true; can turn off manually; automatically turn off if trigger_mode != 0 or 1

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
    std::vector<int> levels_trg_effcy_to_be_summed_w_musign_summing;
    std::vector<std::vector<std::string>> levels_trg_effcy_filters_to_be_summed_w_musign_summing;
    std::vector<std::string> trg_effcy_filters_to_be_summed_w_musign_summing;
    
    std::vector<std::vector<std::string>> levels_trg_effcy_filters_1D_post_sum_w_musign_summing;
    std::vector<std::vector<std::string>> levels_trg_effcy_filters_2D_3D_post_sum_w_musign_summing;

    std::vector<std::string> trg_effcy_filters_1D_post_sum_w_musign_summing;
    std::vector<std::string> trg_effcy_filters_2D_3D_post_sum_w_musign_summing;

// --------------------- protected class methods ---------------------------

    virtual void        Initialize() override;
    virtual void        TrigEffcyFiltersPrePostSumFlattening() override;
    
    virtual void        FillHistogramsSingleMuonEffcy() override;
    virtual void        FillHistogramsDimuTrigGivenMu4() override;
    virtual void        FillHistogramsMu4GivenMB() override;
    virtual void        FillTrigEffcyHistsInvWeightedbySingleMuonEffcies() override;

    virtual void        OpenEffcyPtFitFile() override;

    virtual void        SumSingleMuonTrigEffHists() override;

    virtual void        CalculateSingleMuonTrigEffcyRatios() override;
    virtual void        MakeAndWriteSingleMuonTrigEffPtGraphs() override;
    void                MakeAndWriteDRTrigEffGraphs() override;

public:

    explicit RDFBasedHistFillingPP(int run_year_input, bool isForSoumya_input = false)
    : RDFBasedHistFillingData (run_year_input, isForSoumya_input){
        std::cout << "constructor for pp called" << std::endl; 
    }
    ~RDFBasedHistFillingPP(){
        std::cout << "destructor for pp called" << std::endl;
    }

};

// ----------------------------------------------------------------------------------------------------------------
class RDFBasedHistFillingPbPb : public RDFBasedHistFillingData, public PbPbBaseClass{
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

    virtual void        Initialize() override;
    virtual void        BuildHistBinningMap() override;
    virtual void        BuildTrgEffcyFilterToVarListMap() override;
    virtual void        TrigEffcyFiltersPrePostSumFlattening() override;

    virtual void        CreateRDFs() override;
    virtual void        FillHistogramsSingleMuonEffcy() override;
    virtual void        FillHistogramsDimuTrigGivenMu4() override;
    virtual void        FillHistogramsDimuTrigGivenMu4CtrDep();
    virtual void        FillHistogramsMu4GivenMB() override;
    virtual void        FillHistogramsMu4GivenMBCtrDep();
    virtual void        FillTrigEffcyHistsInvWeightedbySingleMuonEffcies() override;

    virtual void        OpenEffcyPtFitFile() override;

    virtual void        HistPostProcess() override;
    virtual void        SumSingleMuonTrigEffHists() override;

    virtual void        CalculateSingleMuonTrigEffcyRatios() override;
    virtual void        MakeAndWriteSingleMuonTrigEffPtGraphs() override;
    void                MakeAndWriteDRTrigEffGraphs() override;
    void                MakeAndWriteSingleMuonTrigEffCtrGraphs();

    AxisInfo            GetAxisInfo(const var1D& v, const std::string& filter) const override;
    void                ReadVar1DJson() override;

public:
    explicit RDFBasedHistFillingPbPb(int run_year_input, std::string ctr_binning_version_input = "include_upc")
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
