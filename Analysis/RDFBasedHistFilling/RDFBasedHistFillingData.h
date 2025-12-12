#pragma once
#include "../MuonObjectsParamsAndHelpers/PbPbBaseClass.h"
#include "../MuonObjectsParamsAndHelpers/TrigEffcyUtils.h"
#include "RDFBasedHistFillingBaseClass.cxx"

class RDFBasedHistFillingData : public virtual RDFBasedHistFillingBaseClass{
protected:
// --------------------- class variables ---------------------------

    std::string trig_suffix;

    std::vector<std::string> categories_essential;
    std::vector<std::string> trigs;
    std::vector<std::pair<std::string, std::string>> trigs_pair;

    std::vector<std::string> single_muon_trig_effcy_var1Ds = {"Deta", "Deta_zoomin", "Dphi", "Dphi_zoomin", "DR", "DR_zoomin", "DR_0_2", "minv_zoomin", "pair_pt_log"};
    std::vector<std::array<std::string, 2>> single_muon_trig_effcy_var2Ds = {{"phi2nd","pt2nd"}, {"q_eta2nd","pt2nd"}, {"q_eta2nd","phi2nd"}, {"pt2nd", "DR_zoomin"}, {"pair_pt_log", "DR_zoomin"}};
    std::vector<std::array<std::string, 3>> single_muon_trig_effcy_var3Ds = {{"phi2nd","q_eta2nd","pt2nd"}};

    std::vector<std::vector<std::string>> levels_trg_effcy_filters_1D_pre_sum;
    std::vector<std::vector<std::string>> levels_trg_effcy_filters_2D_3D_pre_sum;
    std::vector<std::vector<std::string>> levels_trg_effcy_filters_to_be_summed;

    std::vector<std::string> trg_effcy_filters_1D_pre_sum;
    std::vector<std::string> trg_effcy_filters_2D_3D_pre_sum;

    std::vector<int> levels_trg_effcy_to_be_summed = {0,1}; // HARD-CODED for now: arbitrary leveling is very complicated; for Pb+Pb, meaning mu1/mu2 pass mu4 MUST precede centrality binning
    std::vector<std::string> trg_effcy_filters_to_be_summed;

    std::vector<std::vector<std::string>> levels_trg_effcy_filters_1D_post_sum;
    std::vector<std::vector<std::string>> levels_trg_effcy_filters_2D_3D_post_sum;

    std::vector<std::string> trg_effcy_filters_1D_post_sum;
    std::vector<std::string> trg_effcy_filters_2D_3D_post_sum;

// --------------------- kinematics configuration ---------------------------

    // maps of q_eta bins to pT projection ranges serving single-muon effciency fitting
    std::vector<std::pair<float, float>> q_eta_proj_ranges_for_single_muon_effcy_pT_fitting = {
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

    std::vector<std::pair<float, float>> pair_pT_ranges_for_weighted_effcy_dR_fitting = {
        {8, 12},
        {12, 24},
        {24, 80}
    };

// --------------------- protected class methods ---------------------------

    virtual void    InitOutput() override;
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
    virtual void    MakeAndWriteSingleMuonPtTrigEffGraphs() = 0;
    void            MakeAndWriteSingleMuonPtTrigEffGraphsHelper(const std::vector<std::string>& categories);
    virtual void    MakeAndWriteDRTrigEffGraphs() = 0;
    void            MakeAndWriteDRTrigEffGraphsHelper(const std::vector<std::string>& categories);

    static std::string      FindBinReturnStr(float number, const std::vector<std::pair<float, float>>& ranges);
    static float            EvaluateSingleMuonEffcyPtFitted(bool charge_sign, std::string trg, float pt_2nd, float q_eta_2nd, float phi_2nd);
    static float            EvaluateSingleMuonEffcy(bool charge_sign, std::string trg, float pt_2nd, float q_eta_2nd, float phi_2nd);

public:
enum HistFillingCycle{
    generic = 1,
    inv_weight_by_single_mu_effcy = 2,
    inv_weight_by_dR_effcy_corr = 3
};

    int hist_filling_cycle = generic;
    int trigger_mode = 1;
    int run_year;
    bool isScram = false;
    bool isTight = false;
    bool doTrigEffcy = true; // default true; can turn off manually; automatically turn off if trigger_mode != 0 or 1

    bool output_generic_hists;
    bool output_gapcut_hists;

    bool filter_out_photo_resn_for_trig_effcy = true;
    bool use_3D_2nd_muon = false; // if true, use 3D kinematics (phi, q*eta, pT) for single (2nd) muon trigger efficiencies
    bool use_pT_fitting_single_muon_effcy = true;


    RDFBasedHistFillingData(){
        std::cout << " Histogram filling for data:" << std::endl;
        std::cout << "The following public variable(s) **MUST** be set:" << std::endl;
        std::cout << "--> run_year: [INT]" << std::endl << std::endl;
        std::cout << "The following public variable(s) **MUST** be checked:" << std::endl;
        std::cout << "--> trigger_mode: [INT]" << std::endl;
        std::cout << "--> hist_filling_cycle: [INT]" << std::endl << std::endl;
        std::cout << "The following public variable(s) **SHOULD** be checked:" << std::endl;
        std::cout << "--> isScram: [BOOL]" << std::endl;
        std::cout << "--> isTight: [BOOL]" << std::endl;
        std::cout << "--> output_generic_hists: [BOOL]" << std::endl;

    }
    ~RDFBasedHistFillingData(){}
    
};


// ----------------------------------------------------------------------------------------------------------------
class RDFBasedHistFillingPP : public RDFBasedHistFillingData{

protected:
    virtual void        Initialize() override;
    
    virtual void        FillHistogramsSingleMuonEffcy() override;
    virtual void        FillHistogramsDimuTrigGivenMu4() override;
    virtual void        FillHistogramsMu4GivenMB() override;
    virtual void        FillTrigEffcyHistsInvWeightedbySingleMuonEffcies() override;

    virtual void        OpenEffcyPtFitFile() override;

    virtual void        CalculateSingleMuonTrigEffcyRatios() override;
    virtual void        MakeAndWriteSingleMuonPtTrigEffGraphs() override;
    void                MakeAndWriteDRTrigEffGraphs() override;

public:

    RDFBasedHistFillingPP(){
        std::cout << "constructor for pp called" << std::endl; 
    }
    ~RDFBasedHistFillingPP(){
        std::cout << "destructor for pp called" << std::endl;
    }

};

// ----------------------------------------------------------------------------------------------------------------
class RDFBasedHistFillingPbPb : public RDFBasedHistFillingData, public PbPbBaseClass{
protected:
    std::vector<std::string> ctr_musign_categories;
    
    std::vector<std::string>                var1D_list_no_ctr_rebin; // list of 1D variables for which we apply NO ctr rebinning
    std::map<std::string, std::vector<int>> var1D_excep_ctr_rebin_factor_map; // map of 1D variables to exception ctr rebinning

    std::vector<std::vector<std::string>> levels_trg_effcy_filters_ctr_dep_1D_pre_sum;
    std::vector<std::string> trg_effcy_filters_ctr_dep_1D_pre_sum;
    std::vector<std::vector<std::string>> levels_trg_effcy_filters_ctr_dep_1D_post_sum;
    std::vector<std::string> trg_effcy_filters_ctr_dep_1D_post_sum;
    
    std::vector<float> pT_bins_edges_for_trg_effcy_ctr_dep = {4, 6, 8, 12, 30}; // pT binning for centrality dependence of single-muon trigger efficiencies
    std::vector<std::string> pT_bins_for_trg_effcy_ctr_dep;

    std::pair<float, float> mid_rapidity_range = {-0.5, 0.5}; // mid-rapidity range for pT-binned centrality dependence of single-muon trigger efficiencies
    
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

    virtual void        SumSingleMuonTrigEffHists() override;

    virtual void        CalculateSingleMuonTrigEffcyRatios() override;
    virtual void        MakeAndWriteSingleMuonPtTrigEffGraphs() override;
    void                MakeAndWriteDRTrigEffGraphs() override;

    AxisInfo            GetAxisInfo(const var1D& v, const std::string& filter) const override;
    void                ReadVar1DJson() override;

public:
    RDFBasedHistFillingPbPb(){
        std::cout << "constructor for PbPb called" << std::endl; 
    }
    ~RDFBasedHistFillingPbPb(){
        std::cout << "destructor for PbPb called" << std::endl;
    }
};
