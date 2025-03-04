#include "PlotMCDataComprBaseClass.h"

void PlotMCDataComprBaseClass::initialize(){
    std::cout << "Base class initialization function" << std::endl;
    for (int idt = 0; idt < s_nDtTypes; idt++){
        f[idt] = TFile::Open((dt_paths[idt] + fnames[idt]).c_str());
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

