#ifndef __PlotMCDataComprBaseClass_h__
#define __PlotMCDataComprBaseClass_h__

#include <TROOT.h>
#include <TFile.h>
#include <string.h>
#include <stdlib.h>
#include "vector"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include <algorithm>
#include "../helper_functions.c"
#include "../DimuonPlottingBaseClass.cxx"
#include "../../MuonObjectsParamsAndHelpers/ParamsSet.h"

class PlotMCDataComprBaseClass : public DimuonPlottingBaseClass{
protected:
    static const int s_nDtTypes = 4;

    enum DataType{
        powheg_bb,
        powheg_cc,
        pythia,
        pp_2024_2mu4
    };

    double pi = acos(-1.0);
    std::array<float, s_nDtTypes> norm_factor = {1., 1., 1., 1./410.815};
    std::array<Color_t, s_nDtTypes> colors = {kRed, kRed, kBlue, kGreen+2};

    std::string pythia_with_data_resonance_cuts_dir;
    std::string pythia_with_data_resonance_cuts_suffix;

    std::array<TFile*, s_nDtTypes> f;
    std::array<bool, s_nDtTypes> is_data = {false, false, false, true};
    std::array<std::string, s_nDtTypes> dt_paths;
    std::array<std::string, s_nDtTypes> fnames;
    std::array<std::string, s_nDtTypes> dtTitles = {"POWHEG", "", "Pythia", "pp data 2024"};

    std::vector<std::string> powheg_mechanisms = {
        "gg", "qg", "single_g", "qq", "from_same_b", "others",
        "single_b", "bb", "cc", "other_flavors"
    };

    // Build histogram name for each data type.
    //   kin: e.g. "DR", "Dphi"
    //   isign: 0=SS, 1=OS
    //   jacobian: if true, append jacobian_corrected
    //   gapcut: if true, append gapcut suffix (data only)
    // Naming conventions differ:
    //   Data:   h_{kin}_{ss|op}[_jacobian_corrected][_wgapcut]
    //   Pythia: h_{kin}[_jacobian_corrected]_{sign1|sign2}
    //   POWHEG: h_{kin}_{sign1|sign2}_{mechanism}  (no jacobian)
    std::string BuildHistName(int idt, const std::string& kin_name, int isign,
                              bool jacobian, bool gapcut) const {
        if (is_data[idt]) {
            std::string sign_s = (isign == 0) ? "_ss" : "_op";
            std::string gap_s = gapcut ? "_wgapcut" : "";
            std::string jac_s = jacobian ? "_jacobian_corrected" : "";
            return "h_" + kin_name + sign_s + gap_s + jac_s;
        }
        if (idt == DataType::pythia) {
            std::string jac_s = jacobian ? "_jacobian_corrected" : "";
            std::string sign_s = (isign == 0) ? "_sign1" : "_sign2";
            return "h_" + kin_name + jac_s + sign_s;
        }
        // POWHEG: no jacobian, mechanism suffix handled by GetHist1D
        std::string sign_s = (isign == 0) ? "_sign1" : "_sign2";
        return "h_" + kin_name + sign_s;
    }

    TH1D* GetHist1D(int idt, const std::string& base_name) const;

    virtual void initialize();
    virtual void set_legend_position();

public:
    bool pythia_with_data_resonance_cuts = false;

    std::array<float,4> legend_position_same_sign = {s_NULL, s_NULL, s_NULL, s_NULL};
    std::array<float,4> legend_position_opp_sign = {s_NULL, s_NULL, s_NULL, s_NULL};
};

#endif
