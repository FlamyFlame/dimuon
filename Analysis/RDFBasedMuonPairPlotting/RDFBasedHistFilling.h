#ifndef __RDFBasedHistFilling_h__
#define __RDFBasedHistFilling_h__

#include <map>
#include <vector>
#include <string>
#include <utility> 
#include <numeric>
#include "../MuonObjectsParamsAndHelpers/ParamsSet.h"
#include "time.h"
#include <ROOT/RLogger.hxx>
#include <ROOT/RDF/InterfaceUtils.hxx>

#include <nlohmann/json.hpp>
using nlohmann::json;



/* ----------------------- simple POD ----------------------- */
struct var1D {
    std::string var{""}; // e.g, dR / minv / m1.pt
    std::string name{""}; // e.g, dR_zoomin / minv_log / pt_lead
    std::string title{""}; // e.g, #Delta R / m_{#mu#mu}
    int         nbins{0};
    int         nbins_ss{0};
    int         nbins_op{0};
    double      vmin{0.};
    double      vmax{0.};
    std::vector<double> bins{};          // optional variable binning
    std::vector<double> bins_ss{};          // optional variable binning
    std::vector<double> bins_op{};          // optional variable binning

    bool isValid() const {
        if (nbins_ss > 0){ // ss/op-separate configurations in Json file
            if (nbins_op <= 0) return false;
            if (bins_ss.empty() || static_cast<int>(bins_ss.size()) != nbins_ss + 1) return false;
            if (bins_op.empty() || static_cast<int>(bins_op.size()) != nbins_op + 1) return false;
            return true;
        }

        // no ss/op distinction at Json file level
        if (nbins <= 0) return false;
        if (!bins.empty()) return (static_cast<int>(bins.size()) == nbins + 1);
        return vmin < vmax;
    }
};


class RDFBasedHistFilling{
protected:
	ParamsSet pms;

	std::string tree_ss = "muon_pair_tree_sign1";
	std::string tree_op = "muon_pair_tree_sign2";

    std::vector<std::string> input_files {};
    std::string output_file;

    std::vector<std::unique_ptr<ROOT::RDataFrame>> rdf_store; // pointers for the original dataframes

    std::map<std::string, ROOT::RDF::RNode> df_map; // Map of dataframe label to correspondant RNode

	std::map<std::string, TH1D*> hist1D_map;
	std::map<std::string, TH2D*> hist2D_map;
    std::map<std::string, TH3D*> hist3D_map;

	std::map<std::string, ROOT::RDF::RResultPtr< ::TH1D>> hist1d_rresultptr_map;
	std::map<std::string, ROOT::RDF::RResultPtr< ::TH2D>> hist2d_rresultptr_map;
    std::map<std::string, ROOT::RDF::RResultPtr< ::TH3D>> hist3d_rresultptr_map;

	std::vector<std::string> hists_to_not_write {}; // list of 1, 2, 3D histogram names to not write (e.g, histograms for sub-DF's where only the added histograms should be written)
   
    std::map<std::string, std::vector<double>> hist_binning_map;
    std::string infile_var1D_json = "var1D.json";
    std::map<std::string, var1D*> var1D_dict; // dictionary for 1D variables (to be used in 1D & 2D histograms)

    std::vector<std::string> single_muon_trig_effcy_var1Ds = {"Deta", "Deta_zoomin", "Dphi", "Dphi_zoomin", "DR", "DR_zoomin", "DR_0_2", "minv_zoomin", "pair_pt_log"};
    std::vector<std::array<std::string, 2>> single_muon_trig_effcy_var2Ds = {{"phi2nd","pt2nd"}, {"q_eta2nd","pt2nd"}, {"pt2nd", "DR_zoomin"}, {"pair_pt", "DR_zoomin"}};
    std::vector<std::array<std::string, 3>> single_muon_trig_effcy_var3Ds = {{"phi2nd","q_eta2nd","pt2nd"}};

    // ----- map of df filters or <filter, weight> pairs to variables to be plotted -----
    std::map<std::string, std::vector<std::string>> df_filter_to_var1D_list_map; // map of a dataframe filter to a list of 1D variables to plot for the filtered sample
    std::map<std::pair<std::string, std::string>, std::vector<std::string>> df_filter_and_weight_to_var1D_list_map; // map of a <dataframe filter, weight> pair to a list of 1D variables to plot for the filtered sample

    std::map<std::string, std::vector<std::array<std::string, 2>>> df_filter_to_var2D_list_map; // map of a dataframe filter to a list of 2D variables to plot for the filtered sample
    std::map<std::pair<std::string, std::string>, std::vector<std::array<std::string, 2>>> df_filter_and_weight_to_var2D_list_map; // map of a <dataframe filter, weight> pair to a list of 2D variables to plot for the filtered sample

    std::map<std::string, std::vector<std::array<std::string, 3>>> df_filter_to_var3D_list_map; // map of a dataframe filter to a list of 2D variables to plot for the filtered sample
    std::map<std::pair<std::string, std::string>, std::vector<std::array<std::string, 3>>> df_filter_and_weight_to_var3D_list_map; // map of a <dataframe filter, weight> pair to a list of 2D variables to plot for the filtered sample
    
    // ----- map of weight specifier to column -----
    std::map<std::string, std::string> weight_specifier_to_column_map; // map of a string specifying a weight (e.g, "_jacobian_corrected" or "_inv_w_by_single_mu_effcy") to the column used for weight

    std::vector<std::vector<std::string>> levels_trg_effcy_filters_1D_pre_sum;
    std::vector<std::vector<std::string>> levels_trg_effcy_filters_2D_3D_pre_sum;
    std::vector<std::vector<std::string>> levels_trg_effcy_filters_to_be_summed;

    std::vector<std::string> trg_effcy_filters_1D_pre_sum;
    std::vector<std::string> trg_effcy_filters_2D_3D_pre_sum;
    std::vector<std::string> trg_effcy_filters_to_be_summed;

    std::vector<std::vector<std::string>> levels_trg_effcy_filters_1D_post_sum;
    std::vector<std::vector<std::string>> levels_trg_effcy_filters_2D_3D_post_sum;

    std::vector<std::string> trg_effcy_filters_1D_post_sum;
    std::vector<std::string> trg_effcy_filters_2D_3D_post_sum;

    std::vector<int> levels_trg_effcy_to_be_summed = {0,1}; // HARD-CODED for now: arbitrary leveling is very complicated; for Pb+Pb, meaning mu1/mu2 pass mu4 MUST precede centrality binning



// --------------------- class methods ---------------------------

    // ----- main functions in workflow -----
    virtual void    Initialize();
    virtual void    ProcessData();
    virtual void    Finalize();

    // ----- histogram configuration helpers -----
    virtual void    BuildHistBinningMap();
    
    // ----- json reading -----
    static void ThrowMissingField(const std::string& field,
                                  const std::string& histName); // helper function for ReadVar1DJson()
    virtual void    ReadVar1DJson();
    void            PrintVar1DList() const; //print list of 1D variables for testing
    
    // ----- filter & filter-to-variable mapping helpers -----
    virtual void    TrigEffcyFiltersPrePostSumFlattening();
    virtual void    BuildTrgEffcyFilterToVarListMap();
    virtual void    BuildFilterToVarListMap();

    // ----- data processing functions -----
    void            OpenEffcyPtFitFile();
    bool            PassSingleMuonGapCut(float meta, float mpt, int mcharge);
    std::string     FindBinReturnStr(float number, const std::vector<std::pair<float, float>>& ranges);
    
    void            FillHistogramsSingleDataFrame(  const std::string& filter,
                                                    ROOT::RDF::RNode df,
                                                    bool hists_not_write = false,
                                                    std::array<bool, 3> hists_1_2_3D_not_write = {0,0,0});

    void            FillHistogramsSingleDataFrame(  const std::string& filter,
                                                    const std::string& weight,
                                                    ROOT::RDF::RNode df,
                                                    bool weight_before_filter = false,
                                                    bool hists_not_write = false,
                                                    std::array<bool, 3> hists_1_2_3D_not_write = {0,0,0});

    void             FillHistogramsSingleDataFrame( const std::string& suffix, // filter or filter & weight concatenated with custom order
                                                    ROOT::RDF::RNode df,
                                                    const std::string& weight_col, // weight column, "" if unweighted
                                                    const std::vector<std::string>& vars1D,
                                                    const std::vector<std::array<std::string, 2>>& vars2D, 
                                                    const std::vector<std::array<std::string, 3>>& vars3D,
                                                    bool hists_not_write = false,
                                                    std::array<bool, 3> hists_1_2_3D_not_write = {0,0,0});
    virtual void    FillHistograms();
    virtual void    HistogramPostProcess();

    void            SumSingleMuonTrigEffHists();
    void            MakeAndWriteSingleMuonPtTrigEffGraphs();
    void            CalculateSingleMuonTrigEffcyRatios();
    
    static float    EvaluateSingleMuonEffcyPtFitted(bool charge_sign, std::string trg, float pt_2nd, float q_eta_2nd, float phi_2nd);
    static float    EvaluateSingleMuonEffcy(bool charge_sign, std::string trg, float pt_2nd, float q_eta_2nd, float phi_2nd);

    void            FillTrigEffcyHistsInvWeightedbySingleMuonEffcies();

    void            MakeAndWriteDRTrigEffGraphs();


    var1D*          Var1DSearch(const std::string& var1DName) const;
    AxisInfo        GetAxisInfo(const var1D& v, const std::string& filter) const;

    ROOT::RDF::TH2DModel MakeTH2DModel(const std::string& hname,
                                       const std::string& xtitle,
                                       const std::string& ytitle,
                                       const AxisInfo& bx,
                                       const AxisInfo& by) const;

    ROOT::RDF::TH3DModel MakeTH3DModel(const std::string& hname,
                                       const std::string& xtitle,
                                       const std::string& ytitle,
                                       const std::string& ztitle,
                                       const AxisInfo& bx,
                                       const AxisInfo& by,
                                       const AxisInfo& bz) const;

public:
enum HistFillingCycle{
    generic = 1,
    inv_weight_by_single_mu_effcy = 2,
    inv_weight_by_dR_effcy_corr = 3
};

    bool output_non_trig_effcy_hists;
    bool isScram;
    bool isTight;
    bool isRun3;
    // bool isMu4mu4noL1;
    int trigger_mode = 1;
    int hist_filling_cycle = generic;
    // bool require_exclusive_mu4_for_mu4_mu4noL1 = false; // if true, require that only one muon fire mu4, to ensure unambiguous knowledge of which muon fires mu4noL1
    bool doTrigEffcy = true;
    bool filter_out_photo_resn_for_trig_effcy = true;

    bool use_3D_2nd_muon = false; // if true, use 3D kinematics (phi, q*eta, pT) for single (2nd) muon trigger efficiencies
  	
    bool use_pT_fitting_single_muon_effcy = true;

    RDFBasedHistFilling(){}
  	~RDFBasedHistFilling(){}
  	void Run();

};

// ----------------------------------------------------------------------------------------------------------------
class RDFBasedHistFillingPP : public RDFBasedHistFilling{
public:

    RDFBasedHistFillingPP(){
        std::cout << "constructor for pp called" << std::endl; 
    }
    ~RDFBasedHistFillingPP(){
        std::cout << "destructor for pp called" << std::endl;
    }
    
    virtual void        Initialize() override{ // placeholder for now
        RDFBasedHistFilling::Initialize();
    }
};

#endif