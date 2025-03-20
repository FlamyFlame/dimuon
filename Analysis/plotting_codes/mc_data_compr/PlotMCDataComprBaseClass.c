#include "PlotMCDataComprBaseClass.h"

void PlotMCDataComprBaseClass::initialize(){
    std::cout << "Base class initialization function" << std::endl;

    pythia_with_data_resonance_cuts_dir = pythia_with_data_resonance_cuts? "with_data_resonance_cuts/" : "no_data_resonance_cuts/";
    pythia_with_data_resonance_cuts_suffix = pythia_with_data_resonance_cuts? "_with_data_resonance_cuts" : "_no_data_resonance_cuts";

    dt_paths = {"/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample/bb_full_sample/","/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample/cc_full_sample/","/usatlas/u/yuhanguo/usatlasdata/pythia/","/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_run2/","/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/","/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/"};
    fnames = {"histograms_mc_truth_bb_combined.root","histograms_mc_truth_cc_combined.root","histograms_pythia_combined" + pythia_with_data_resonance_cuts_suffix + ".root","histograms_real_pairs_pp_2017.root","histograms_real_pairs_pp_2024_2mu4.root","histograms_real_pairs_pp_2024_mu4_mu4noL1.root"};

    for (int idt = 0; idt < s_nDtTypes; idt++){
        std::cout << "opening file with the name " << dt_paths[idt] << fnames[idt] << endl;
        f[idt] = TFile::Open((dt_paths[idt] + fnames[idt]).c_str());
        if (!f[idt]){
            std::cout << "File with the name " << dt_paths[idt] << fnames[idt] << "does not exist or cannot be opened";
            throw std::exception();
        }
    }

    set_legend_position();
}

void PlotMCDataComprBaseClass::set_legend_position(){
    std::cout << "Base class function to set default legend positions" << std::endl;
    if (legend_position_same_sign[0] < -999.){ // no user input --> set default values
        legend_position_same_sign = {0.5,0.6,0.95,0.89};
    }
    if (legend_position_opp_sign[0] < -999.){ // no user input --> set default values
        legend_position_opp_sign = {0.5,0.6,0.93,0.89};
    }
}

void PlotMCDataComprBaseClass::hist_helper(TH1* h, float norm, bool norm_unity, std::string title, std::string ytitle=""){

  h->SetStats(0);

    if (norm_unity){
        h->Scale(1.,"width");
        try{
            if (h->Integral("width") == 0) throw std::runtime_error("Histogram interval is ZERO! Cannot normalize to unity!");
            h->Scale(1./h->Integral("width"));
        }catch (const std::runtime_error& e){
            std::cout << "runtime_error caught: " << e.what() << std::endl;
            std::cout << "Proceed without normalizing to unity!" << std::endl; 
        }
        h->GetYaxis()->SetTitle("pdf");
    }else{
        h->Scale(norm,"width");
        if (ytitle.length() != 0){
            h->GetYaxis()->SetTitle(ytitle.c_str());
        }
    }

    // h->SetTitle(title.c_str());
    // h->SetTitleSize(35);
    h->GetYaxis()->SetLabelFont(43);
    h->GetYaxis()->SetLabelSize(28);
    h->GetYaxis()->SetLabelOffset(0.01);
    h->GetYaxis()->SetTitleFont(43);
    h->GetYaxis()->SetTitleSize(28);
    h->GetYaxis()->SetTitleOffset(2.1);
    h->GetXaxis()->SetLabelFont(43);
    h->GetXaxis()->SetLabelSize(28);
    h->GetXaxis()->SetLabelOffset(0.01);
    h->GetXaxis()->SetTitleFont(43);  
    h->GetXaxis()->SetTitleSize(28);
    h->GetXaxis()->SetTitleOffset(1);
    h->SetMarkerStyle(20);
    h->SetMarkerSize(0.9);
}

