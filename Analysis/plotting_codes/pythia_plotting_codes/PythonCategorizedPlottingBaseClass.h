#pragma once

#include <TROOT.h>
// #include <TChain.h>
#include <TFile.h>
#include <string.h>
#include  <stdlib.h>
#include "vector"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include<algorithm>
#include "../helper_functions.c"
#include "../../MuonObjectsParamsAndHelpers/ParamsSet.h"
#include "../../MuonObjectsParamsAndHelpers/muon_pair_enums_MC.h"
// #include "string"
// #include <vector>
// #include "time.h"
// #include "struct_hist.h"


struct line{ // attributes to each histogram appear as a line in a subplot
    std::string line_category;
    std::string line_label;
    Color_t line_color;
    Color_t fill_color;
};

enum sign{
    same_sign = 0,
    op_sign = 1,
    nSigns
};

class PythonCategorizedPlottingBaseClass{
protected:

    double min_integral_thrsh = 1e-6;

    std::map<int, line*> line_map;
    std::map<int, std::vector<int>> thstack_order_map;
    std::map<std::pair<int, int>, TH1D*> hist_map;

    TH2D* h2d;
    TFile* f;

    std::string with_data_resonance_cuts_dir;
    std::string with_data_resonance_cuts_suffix;

    std::string pythia_path = "/usatlas/u/yuhanguo/usatlasdata/pythia_private_sample/";
    std::string fname;
    std::map<int, std::string> signs;
    std::map<int, std::string> subpl_titles;

    std::vector<std::string> single_b_analysis_observables = {"DR_zoomin", "Deta_zoomin", "Dphi_zoomin", "DR_zoomin_jacobian_corrected", "minv_zoomin", "minv_zoomin_jacobian_corrected"};

    virtual void initialize();
    virtual void fill_line_map() = 0;
    virtual void fill_thstack_order_map() = 0;
public:
    bool with_data_resonance_cuts = false;

    std::string kin;
    bool projx_2d;
    bool projy_2d;
    bool staggered;
    bool norm_unity;
    std::string kin1d;
    std::string kin_title;
    std::vector<std::array<float,2>> cuts;
    bool logx=false;

    std::string subdir_name = "";
    std::string optional_suffix = ""; // allow child classes to have an optional suffix if needed

    PythonCategorizedPlottingBaseClass();
	PythonCategorizedPlottingBaseClass(std::string kin_in, bool projx_2d_in, bool projy_2d_in, bool staggered_in, bool norm_unity_in, std::string kin1d_in, std::string kin_title_in): projx_2d (projx_2d_in), projy_2d (projy_2d_in), staggered (staggered_in), norm_unity (norm_unity_in){
	    kin = kin_in;
	    kin1d = kin1d_in;
	    kin_title = kin_title_in;
	}

	PythonCategorizedPlottingBaseClass(std::string kin_in, bool projx_2d_in, bool projy_2d_in, bool staggered_in, bool norm_unity_in, std::string kin1d_in, std::string kin_title_in, std::vector<std::array<float,2>> cuts_in): projx_2d (projx_2d_in), projy_2d (projy_2d_in), staggered (staggered_in), norm_unity (norm_unity_in){
	    kin = kin_in;
	    kin1d = kin1d_in;
	    kin_title = kin_title_in;
	    cuts = cuts_in;
	}
    ~PythonCategorizedPlottingBaseClass();
    void Run();
};

