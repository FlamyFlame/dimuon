#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include "../MuonObjectsParamsAndHelpers/ParamsSet.h"

struct var1D {
    std::string var{""};   // column name
    std::string name{""};  // histogram base name
    std::string title{""}; // axis title
    int         nbins{0};
    int         nbins_ss{0};
    int         nbins_op{0};
    double      vmin{0.};
    double      vmax{0.};
    std::vector<double> bins{};
    std::vector<double> bins_ss{};
    std::vector<double> bins_op{};

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

class TestClass {
public:
    void Run();

private:
    ParamsSet pms;
    
    std::string tree_ss = "muon_pair_tree_sign1";
    std::string tree_op = "muon_pair_tree_sign2";

    std::vector<std::string> input_files {};
    std::string output_file;

    std::map<std::string, std::vector<double>> hist_binning_map;
    std::string infile_var1D_json = "var1D.json";
    std::map<std::string, var1D*> var1D_dict;

    std::vector<std::string> hists_to_not_write {}; // list of 1, 2, 3D histogram names to not write (e.g, histograms for sub-DF's where only the added histograms should be written)

    std::map<std::string, std::vector<std::string>> df_filter_to_var1D_list_map;
    std::map<std::pair<std::string, std::string>, std::vector<std::string>> df_filter_and_weight_to_var1D_list_map;

    std::map<std::string, std::vector<std::array<std::string, 2>>> df_filter_to_var2D_list_map;
    std::map<std::pair<std::string, std::string>, std::vector<std::array<std::string, 2>>> df_filter_and_weight_to_var2D_list_map;

    std::map<std::string, std::vector<std::array<std::string, 3>>> df_filter_to_var3D_list_map;
    std::map<std::pair<std::string, std::string>, std::vector<std::array<std::string, 3>>> df_filter_and_weight_to_var3D_list_map;

    std::map<std::string, std::string> weight_specifier_to_column_map;

    std::map<std::string, ROOT::RDF::RResultPtr< ::TH1D>> hist1d_rresultptr_map;
    std::map<std::string, ROOT::RDF::RResultPtr< ::TH2D>> hist2d_rresultptr_map;
    std::map<std::string, ROOT::RDF::RResultPtr< ::TH3D>> hist3d_rresultptr_map;


    virtual void    Initialize();
    virtual void    ProcessData();
    virtual void    Finalize();

    void FillHistogramsSingleDataFrame(const std::string& filter,
                                            ROOT::RDF::RNode df,
                                            bool hists_not_write = false,
                                            std::array<bool, 3> hists_1_2_3D_not_write = {0,0,0});

    void FillHistogramsSingleDataFrame(const std::string& filter,
                                            const std::string& weight,
                                            ROOT::RDF::RNode df,
                                            bool weight_before_filter = false,
                                            bool hists_not_write = false,
                                            std::array<bool, 3> hists_1_2_3D_not_write = {0,0,0});
    
    void FillHistogramsSingleDataFrame(const std::string& suffix, // filter or filter & weight concatenated with custom order
                                            ROOT::RDF::RNode df,
                                            const std::string& weight_col, // weight column, "" if unweighted
                                            const std::vector<std::string>& vars1D,
                                            const std::vector<std::array<std::string, 2>>& vars2D, 
                                            const std::vector<std::array<std::string, 3>>& vars3D,
                                            bool hists_not_write = false,
                                            std::array<bool, 3> hists_1_2_3D_not_write = {0,0,0});

    void BuildHistBinningMap();
    static void ThrowMissingField(const std::string& field,
                                  const std::string& histName);

    void ReadVar1DJson();
    void PrintVar1DList() const;

    var1D* Var1DSearch(const std::string& var1DName) const;
    AxisInfo GetAxisInfo(const var1D& v, const std::string& filter) const;

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
};
