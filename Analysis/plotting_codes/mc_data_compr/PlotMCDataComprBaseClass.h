#ifndef __MuonPairPlotting_h__
#define __MuonPairPlotting_h__


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

class PlotMCDataComprBaseClass{
protected:
    // private class variables
    static const int s_nDtTypes = 6;
    static const int s_nSigns = 2;
    static const int s_nDphis = 2;

    static constexpr float s_NULL = -1000.;

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
    float norm_factor[s_nDtTypes] = {1., 1., 1., 1/256.8, 1./410.815, 1/113.999}; // normalizing to differential crossx; unit is pb
    Color_t colors[s_nDtTypes] = {kRed, kRed, kBlue, kBlack, kGreen+2, kMagenta};

    std::string signs[s_nSigns] = {"_sign1", "_sign2"};
    std::string signTitles[s_nSigns] = {"same sign", "opposite sign"};

    std::string dphis[s_nDphis] = {"_near", "_away"};
    std::string dphiTitles[s_nDphis] = {"#Delta #phi < #pi/2","#Delta #phi #geq #pi/2"};

    TFile* f[s_nDtTypes];
    bool is_data[s_nDtTypes] = {false, false, false, true, true, true}; //for data files, the histogram names need to have _gapcut1/2 suffix
    std::string dt_paths[s_nDtTypes] = {"/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample/bb_full_sample/","/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample/cc_full_sample/","/usatlas/u/yuhanguo/usatlasdata/pythia/","/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_run2/","/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/","/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/"};
    std::string fnames[s_nDtTypes] = {"histograms_mc_truth_bb_combined.root","histograms_mc_truth_cc_combined.root","histograms_pythia_combined.root","histograms_real_pairs_pp_2017.root","histograms_real_pairs_pp_2024_2mu4.root","histograms_real_pairs_pp_2024_mu4_mu4noL1.root"};
    std::string dtTitles[s_nDtTypes] = {"POWHEG", "", "Pythia", "pp data 2017", "pp data 2024 2mu4", "pp data 2024 mu4mu4noL1"};

    // private class methods
    virtual void initialize();
    virtual void set_legend_position();
    virtual void hist_helper(TH1* h, float norm, bool norm_unity, std::string title, std::string ytitle="");

public:

    // public class variables
    std::string kin;
    std::string kin_title;

    std::array<float,4> legend_position_same_sign = {s_NULL, s_NULL, s_NULL, s_NULL}; // xmin, xmax, ymin, ymax
    std::array<float,4> legend_position_opp_sign = {s_NULL, s_NULL, s_NULL, s_NULL}; // xmin, xmax, ymin, ymax

    virtual void Run() = 0;
};

#endif