#include "PlotMCDataComprBaseClass.h"

void PlotMCDataComprBaseClass::initialize(){
    std::cout << "Base class initialization function" << std::endl;

    pythia_with_data_resonance_cuts_dir = pythia_with_data_resonance_cuts? "with_data_resonance_cuts/" : "no_data_resonance_cuts/";
    pythia_with_data_resonance_cuts_suffix = pythia_with_data_resonance_cuts? "_with_data_resonance_cuts" : "_no_data_resonance_cuts";

    dt_paths = {
        "/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample/bb_evgen_truth_full_sample/",
        "/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample/cc_evgen_truth_full_sample/",
        "/usatlas/u/yuhanguo/usatlasdata/pythia_private_sample/",
        "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/"
    };
    fnames = {
        "histograms_mc_truth_bb_combined.root",
        "histograms_mc_truth_cc_combined.root",
        "histograms_pythia_combined" + pythia_with_data_resonance_cuts_suffix + ".root",
        "histograms_real_pairs_pp_2024_2mu4_nominal.root"
    };

    for (int idt = 0; idt < s_nDtTypes; idt++){
        std::cout << "opening file: " << dt_paths[idt] << fnames[idt] << std::endl;
        f[idt] = TFile::Open((dt_paths[idt] + fnames[idt]).c_str());
        if (!f[idt] || f[idt]->IsZombie()){
            throw std::runtime_error("Cannot open file: " + dt_paths[idt] + fnames[idt]);
        }
    }

    set_legend_position();
}

void PlotMCDataComprBaseClass::set_legend_position(){
    if (legend_position_same_sign[0] < -999.){
        legend_position_same_sign = {0.5,0.6,0.95,0.89};
    }
    if (legend_position_opp_sign[0] < -999.){
        legend_position_opp_sign = {0.5,0.6,0.93,0.89};
    }
}

TH1D* PlotMCDataComprBaseClass::GetHist1D(int idt, const std::string& base_name) const {
    if (idt == DataType::powheg_bb || idt == DataType::powheg_cc) {
        // POWHEG: no inclusive histogram; sum over production mechanisms
        TH1D* h_sum = nullptr;
        for (const auto& mech : powheg_mechanisms) {
            std::string hname = base_name + "_" + mech;
            TH1D* hm = (TH1D*)f[idt]->Get(hname.c_str());
            if (!hm) {
                std::cerr << "[WARN] Missing POWHEG mechanism hist: " << hname
                          << " in " << fnames[idt] << std::endl;
                continue;
            }
            if (!h_sum) {
                h_sum = (TH1D*)hm->Clone(Form("hsum_%s_%d", base_name.c_str(), idt));
            } else {
                h_sum->Add(hm);
            }
        }
        if (!h_sum) {
            throw std::runtime_error("GetHist1D: no POWHEG mechanism hists found for " + base_name);
        }
        return h_sum;
    }

    // Pythia and data: direct lookup
    TH1D* h = (TH1D*)f[idt]->Get(base_name.c_str());
    if (!h) {
        throw std::runtime_error("GetHist1D: histogram '" + base_name + "' not found in " + fnames[idt]);
    }
    return (TH1D*)h->Clone(Form("h_%s_%d", base_name.c_str(), idt));
}
