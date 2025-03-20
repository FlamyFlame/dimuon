#ifndef __PlotMCDataComprBaseClass_h__
#define __PlotMCDataComprBaseClass_h__


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
#include "../DimuonPlottingBaseClass.cxx" // bad practice but keep it for now
#include "../../MuonObjectsParamsAndHelpers/ParamsSet.h"


class PlotMCDataComprBaseClass : public DimuonPlottingBaseClass{
protected:
    // private class variables
    static const int s_nDtTypes = 6;

    enum DataType{
        powheg_bb,
        powheg_cc,
        pythia,
        pp_2017,
        pbpb_2024_2mu4,
        pbpb_2024_mu4mu4
    };

    double pi = acos(-1.0);
    // NEED TO CHANGE TO 2024 LUMINOSITY
    std::array<float, s_nDtTypes> norm_factor = {1., 1., 1., 1/256.8, 1./410.815, 1/113.999}; // normalizing to differential crossx; unit is pb
    std::array<Color_t, s_nDtTypes> colors = {kRed, kRed, kBlue, kBlack, kGreen+2, kMagenta};

    std::string pythia_with_data_resonance_cuts_dir;
    std::string pythia_with_data_resonance_cuts_suffix;

    std::array<TFile*, s_nDtTypes> f;
    std::array<bool, s_nDtTypes> is_data = {false, false, false, true, true, true}; //for data files, the histogram names need to have _gapcut1/2 suffix
    std::array<std::string, s_nDtTypes> dt_paths;
    std::array<std::string, s_nDtTypes> fnames;
    std::array<std::string, s_nDtTypes> dtTitles = {"POWHEG", "", "Pythia", "pp data 2017", "pp data 2024 2mu4", "pp data 2024 mu4mu4noL1"};

    // private class methods
    virtual void initialize();
    virtual void set_legend_position();
    virtual void hist_helper(TH1* h, float norm, bool norm_unity, std::string title, std::string ytitle="");

public:

    // public class variables

    bool pythia_with_data_resonance_cuts = true;

    std::array<float,4> legend_position_same_sign = {s_NULL, s_NULL, s_NULL, s_NULL}; // xmin, xmax, ymin, ymax
    std::array<float,4> legend_position_opp_sign = {s_NULL, s_NULL, s_NULL, s_NULL}; // xmin, xmax, ymin, ymax
};

#endif