#pragma once

#include <map>
#include <vector>
#include <string>
#include <utility> 
#include <numeric>
#include "../MuonObjectsParamsAndHelpers/ParamsSet.h"
#include "../Utilities/bin_number.cxx"
#include "../Utilities/proj_range_to_suffix.cxx"
#include "../Utilities/GeneralUtils.h"
#include "../Utilities/HistFillUtils.h"
#include "time.h"
#include <ROOT/RLogger.hxx>
#include <ROOT/RDF/InterfaceUtils.hxx>

#include <nlohmann/json.hpp>
using nlohmann::json;


/* ----------------------- simple POD ----------------------- */
struct var1D {
    std::string var{""};   // e.g. dR / minv / m1.pt
    std::string name{""};  // e.g. dR_zoomin / minv_log / pt_lead
    std::string title{""}; // e.g. #Delta R / m_{#mu#mu}

    int         nbins{0};
    int         nbins_ss{0};
    int         nbins_op{0};
    double      vmin{0.};
    double      vmax{0.};

    std::vector<double> bins{};     // default binning
    std::vector<double> bins_ss{};  // ss-specific binning
    std::vector<double> bins_op{};  // op-specific binning

    // PbPb-only centrality-dependent binning
    std::map<std::string, int>                nbins_ctr; // key: "_ctr0_5", etc.
    std::map<std::string, std::vector<double>> bins_ctr; // per-ctr bin edges

    bool isValid() const {
        // ss/op-separate configurations in Json file
        if (nbins_ss > 0) {
            if (nbins_op <= 0) return false;
            if (bins_ss.empty() || static_cast<int>(bins_ss.size()) != nbins_ss + 1) return false;
            if (bins_op.empty() || static_cast<int>(bins_op.size()) != nbins_op + 1) return false;
            return true;
        }

        // PbPb centrality-specific configurations: allow var to be valid
        // even if only centrality binning is provided.
        if (!bins_ctr.empty()) {
            for (const auto& kv : bins_ctr) {
                const std::string& ctr = kv.first;
                const auto& edges = kv.second;

                auto it_n = nbins_ctr.find(ctr);
                if (it_n == nbins_ctr.end() || it_n->second <= 0) return false;
                if (static_cast<int>(edges.size()) != it_n->second + 1) return false;
            }
            return true;
        }

        // no ss/op/centrality distinction at Json file level
        if (nbins <= 0) return false;
        if (!bins.empty()) return (static_cast<int>(bins.size()) == nbins + 1);
        return vmin < vmax;
    }
};



class RDFBasedHistFillingBaseClass{
protected:
	ParamsSet pms;

	std::string tree_ss = "muon_pair_tree_sign1";
	std::string tree_op = "muon_pair_tree_sign2";

    std::vector<std::string> musigns = {"_sign1", "_sign2"};
    std::vector<std::string> pair_signs = {"_ss", "_op"};

    std::vector<std::string> input_files {};
    std::string output_file;
    TFile * m_outfile;

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
    std::string infile_var1D_json;
    std::map<std::string, var1D*> var1D_dict; // dictionary for 1D variables (to be used in 1D & 2D histograms)

    // ----- map of df filters or <filter, weight> pairs to variables to be plotted -----
    std::map<std::string, std::vector<std::string>> df_filter_to_var1D_list_map; // map of a dataframe filter to a list of 1D variables to plot for the filtered sample
    std::map<std::pair<std::string, std::string>, std::vector<std::string>> df_filter_and_weight_to_var1D_list_map; // map of a <dataframe filter, weight> pair to a list of 1D variables to plot for the filtered sample

    std::map<std::string, std::vector<std::array<std::string, 2>>> df_filter_to_var2D_list_map; // map of a dataframe filter to a list of 2D variables to plot for the filtered sample
    std::map<std::pair<std::string, std::string>, std::vector<std::array<std::string, 2>>> df_filter_and_weight_to_var2D_list_map; // map of a <dataframe filter, weight> pair to a list of 2D variables to plot for the filtered sample

    std::map<std::string, std::vector<std::array<std::string, 3>>> df_filter_to_var3D_list_map; // map of a dataframe filter to a list of 2D variables to plot for the filtered sample
    std::map<std::pair<std::string, std::string>, std::vector<std::array<std::string, 3>>> df_filter_and_weight_to_var3D_list_map; // map of a <dataframe filter, weight> pair to a list of 2D variables to plot for the filtered sample
    
    // ----- map of weight specifier to column -----
    std::map<std::string, std::string> weight_specifier_to_column_map; // map of a string specifying a weight (e.g, "_jacobian_corrected" or "_inv_w_by_single_mu_effcy") to the column used for weight


// --------------------- class methods ---------------------------


    // ----- main functions in workflow -----
    void            Initialize();
    void            ProcessData();
    void            Finalize();

    // ----- data processing -----
    void            CreateBaseRDFs();
    void            CreateBaseRDFsBaseCommon();
    virtual void    CreateBaseRDFsExtra(){}

    virtual void    FillHistograms() = 0;

    // ----- data post processing -----
    void            HistPostProcess();
    void            HistPostProcessBaseCommon();
    virtual void    HistPostProcessExtra(){}

    // ----- initialize base & analysis settings -----
    void            InitializeBaseCommon();
    virtual void    InitAnalysisSettingsHook(){}

    // ----- initialize I/O -----
    virtual void    SetIOPathsHook() = 0;
    void    InitOutput();

    // ----- histogram binning name to vector map -----
    void            BuildHistBinningMap();
    void            BuildHistBinningMapBaseCommon();
    virtual void    BuildHistBinningMapExtra(){}
    
    // ----- json reading -----
    static void     ThrowMissingField(  const std::string& field,
                                        const std::string& histName); // helper function for ReadVar1DJson()
    virtual void    ReadVar1DJson();
    void            PrintVar1DList() const; //print list of 1D variables for testing
    
    // ----- filter & filter-to-variable mapping helpers -----
    void            BuildFilterToVarListMap();
    void            BuildFilterToVarListMapBaseCommon();
    virtual void    BuildFilterToVarListMapExtra(){}

    // ----- write output -----
    void            WriteOutput();
    void            WriteOutputBaseCommon();
    virtual void    WriteOutputExtra(){}
    
    // ----- cleanup -----
    void            Cleanup();
    void            CleanupBaseCommon();
    virtual void    CleanupExtra(){}

    // ----- hist filling helper functions -----
    
    void            FillHistogramsSingleDataFrame(  const std::string& filter,
                                                    ROOT::RDF::RNode df,
                                                    bool hists_write = true,
                                                    std::array<bool, 3> hists_1_2_3D_write = {1,1,1});

    void            FillHistogramsSingleDataFrame(  const std::string& filter,
                                                    const std::string& weight,
                                                    ROOT::RDF::RNode df,
                                                    bool weight_before_filter = false,
                                                    bool hists_write = true,
                                                    std::array<bool, 3> hists_1_2_3D_write = {1,1,1});

    void             FillHistogramsSingleDataFrame( const std::string& suffix, // filter or filter & weight concatenated with custom order
                                                    ROOT::RDF::RNode df,
                                                    const std::string& weight_col, // weight column, "" if unweighted
                                                    const std::vector<std::string>& vars1D,
                                                    const std::vector<std::array<std::string, 2>>& vars2D, 
                                                    const std::vector<std::array<std::string, 3>>& vars3D,
                                                    bool hists_write = true,
                                                    std::array<bool, 3> hists_1_2_3D_write = {1,1,1});
    

    var1D*                  Var1DSearch(const std::string& var1DName) const;
    virtual AxisInfo        GetAxisInfo(const var1D& v, const std::string& filter) const;

    ROOT::RDF::TH2DModel    MakeTH2DModel(const std::string& hname,
                                       const std::string& xtitle,
                                       const std::string& ytitle,
                                       const AxisInfo& bx,
                                       const AxisInfo& by) const;

    ROOT::RDF::TH3DModel    MakeTH3DModel(const std::string& hname,
                                       const std::string& xtitle,
                                       const std::string& ytitle,
                                       const std::string& ztitle,
                                       const AxisInfo& bx,
                                       const AxisInfo& by,
                                       const AxisInfo& bz) const;

public:

    bool debug_mode = false;
    RDFBasedHistFillingBaseClass(){}
  	~RDFBasedHistFillingBaseClass(){}
  	void Run();

};
