#include "PythonCategorizedPlottingBaseClass.h"

void PythonCategorizedPlottingBaseClass::initialize(){

    with_data_resonance_cuts_dir = with_data_resonance_cuts? "with_data_resonance_cuts/" : "no_data_resonance_cuts/";
    with_data_resonance_cuts_suffix = with_data_resonance_cuts? "_with_data_resonance_cuts" : "_no_data_resonance_cuts";
    // with_data_resonance_cuts_suffix = with_data_resonance_cuts? "_with_data_resonance_cuts" : "";
    
    fname = "histograms_pythia_combined" + with_data_resonance_cuts_suffix + ".root";

    signs[sign::same_sign]  = "_sign1";
    signs[sign::op_sign]    = "_sign2";

    subpl_titles[sign::same_sign]   = "Pythia8, same sign";
    subpl_titles[sign::op_sign]     = "Pythia8, opposite sign";

    fill_line_map();
    fill_thstack_order_map();  

    legend_position_default = {0.65,0.3,0.89,0.7};

    legend_ss_position_map["default"] = {0.65,0.3,0.89,0.7};
    legend_op_position_map["default"] = {0.65,0.3,0.89,0.7};
    
    legend_op_position_map["minv"] =        {0.65,0.2,0.89,0.7};

    legend_ss_position_map["minv_10GeV"] =  {0.25,0.3,0.51,0.7};
    legend_op_position_map["minv_10GeV"] =  {0.5,0.25,0.89,0.7};

    legend_ss_position_map["minv_zoomin"] = {0.2,0.3,0.46,0.7};
    legend_op_position_map["minv_zoomin"] = {0.25,0.3,0.51,0.7};

}

PythonCategorizedPlottingBaseClass::~PythonCategorizedPlottingBaseClass(){
    for (auto& [grp_id, grp_line_ptr] : line_map){
        delete grp_line_ptr;
    }
}

void PythonCategorizedPlottingBaseClass::Run(){

    initialize();

    THStack *hs[sign::nSigns];
    if (projx_2d && projy_2d){
        std::cout << "For single kinematic plotting, can only project either X or Y variable." << std::endl;
        throw std::exception();
    }

    if (staggered && norm_unity){
        std::cout << "Cannot both stagger and normalize to unity." << std::endl;
        throw std::exception();
    }

    TCanvas* c = new TCanvas("c1","c1",1400,500);
    c->Divide(2,1);

    f = TFile::Open((pythia_path + fname).c_str());

    if (!f){
        std::cout << "File with the name " << pythia_path << fname << "does not exist or cannot be opened";
        throw std::exception();
    }
    for (unsigned int ksign = 0; ksign < sign::nSigns; ksign++){

        c->cd(ksign + 1);
        gPad->SetLeftMargin(0.16);
        // gPad->SetRightMargin(0.16);
        gPad->SetBottomMargin(0.135);
        gPad->SetLogx(logx);

        // ---------- legend box ----------
        std::map<std::string, std::array<double, 4>>* legend_position_map = (ksign == 0)? &legend_ss_position_map : &legend_op_position_map;
        std::array<double, 4> legend_position;

        try {
            if (legend_position_map->find(kin1d) != legend_position_map->end()) { // if the kinematic variable has specialized legend position for the same/opposite-sign pairs
                legend_position = legend_position_map->at(kin1d);
            } else { // use default same/opposite-sign-pair legend positions
                legend_position = legend_position_map->at("default");
            }            
        } catch(const std::out_of_range& e){
            std::cerr << "Out of range error caught for kinematic " << kin << " for legend positions! Use default legend position." << std::endl;
            legend_position = legend_position_default;
        }

        TLegend* l = new TLegend(legend_position[0], legend_position[1], legend_position[2], legend_position[3]);
        l->SetBorderSize(0);
        l->SetFillStyle(0);
        l->SetTextFont(43);
        l->SetMargin(0.2);
        l->SetTextColor(1);

        // ---------- text box ----------
        // draw text box with the same x position + sits on top of y position of the legend box
        TPaveText* tbox = new TPaveText(legend_position[0], legend_position[3], legend_position[2], legend_position[3]+0.2, "NDC");

        tbox->SetFillColor(0);       // transparent background
        tbox->SetFillStyle(0);       // fully transparent
        tbox->SetBorderSize(0);
        tbox->SetTextAlign(12);
        
        std::string panel_descr_text = subpl_titles[ksign];
        if(norm_unity)  panel_descr_text += ", unity";

        tbox->AddText(panel_descr_text.c_str());
        tbox->AddText("nCTEQ15npFullNuc_208_82");
        tbox->AddText("p_{T}^{#mu} > 4GeV, |#eta^{#mu}| < 2.4");

        // ---------- THStack ----------
        if (staggered){
            hs[ksign] = new THStack(("hs" + signs[ksign]).c_str(), ("hs" + signs[ksign]).c_str());
        }

        // ---------- loop over the categories ----------
        for (int kgrp : thstack_order_map[ksign]){

            if (!line_map[kgrp]){
                std::cerr << "Line map for the current flavor/origin category enum " << kgrp << " gives nullptr!" << std::endl;
                std::cerr << "Skip the current flavor/origin category" << std::endl;
                continue;
            }

            if (projx_2d){
                h2d = (TH2D*) f->Get(("h_" + kin + signs[ksign] + line_map[kgrp]->line_category).c_str());
                if (norm_unity){
                    hist_map[{ksign, kgrp}] = (TH1D*) h2d->ProjectionX(("h_unity_" + kin1d + signs[ksign] + line_map[kgrp]->line_category).c_str());
                }else{
                    hist_map[{ksign, kgrp}] = (TH1D*) h2d->ProjectionX(("h_" + kin1d + signs[ksign] + line_map[kgrp]->line_category).c_str());
                }
            }
            else if (projy_2d){
                h2d = (TH2D*) f->Get(("h_" + kin + signs[ksign] + line_map[kgrp]->line_category).c_str());
                if (norm_unity){
                    hist_map[{ksign, kgrp}] = (TH1D*) h2d->ProjectionY(("h_unity_" + kin1d + signs[ksign] + line_map[kgrp]->line_category).c_str());
                }else{
                    hist_map[{ksign, kgrp}] = (TH1D*) h2d->ProjectionY(("h_" + kin1d + signs[ksign] + line_map[kgrp]->line_category).c_str());
                }
            }else{ // 1D
                hist_map[{ksign, kgrp}] = (TH1D*) f->Get(("h_" + kin + signs[ksign] + line_map[kgrp]->line_category).c_str());
            }

            if (!hist_map[{ksign, kgrp}]){
                std::cout << "file h_" << kin << line_map[kgrp]->line_category << signs[ksign] << "does not exist." << std::endl;
                throw std::exception();
            }

            // if the argument cuts is not empty, do not draw the bins overlapping with the cuts
            if (!cuts.empty()) ApplyCutsTo1DHistogram(hist_map[{ksign, kgrp}], cuts);

            if (hist_map[{ksign, kgrp}]->Integral() > min_integral_thrsh){ // skipping categories with too low or zero crossx
                float norm = norm_unity? 1. : pow(10,6);
                hist_helper(hist_map[{ksign, kgrp}], norm, norm_unity, "");

                if (!staggered){
                    hist_map[{ksign, kgrp}]->SetLineColor(line_map[kgrp]->line_color);
                    hist_map[{ksign, kgrp}]->SetMarkerColor(line_map[kgrp]->line_color);          
                }else{
                    hs[ksign]->Add(hist_map[{ksign, kgrp}]);
                    hist_map[{ksign, kgrp}]->SetFillColor(line_map[kgrp]->fill_color);          
                }

                if (staggered){
                    l->AddEntry(hist_map[{ksign, kgrp}], (line_map[kgrp]->line_label).c_str(),"f");
                }else{
                    l->AddEntry(hist_map[{ksign, kgrp}], (line_map[kgrp]->line_label).c_str(),"lp");
                }
            }
        } // end loop over flavor/origin categories

        float ymax = hist_map[{ksign, thstack_order_map[ksign].at(0)}]->GetMaximum();
        if (!staggered){

            for (int kgrp : thstack_order_map[ksign]){
                ymax = (ymax > hist_map[{ksign, kgrp}]->GetMaximum())? ymax : hist_map[{ksign, kgrp}]->GetMaximum();
            }

            // find the same index where the category has nonnegligible integral
            int first_ind = 0;
            while (first_ind < thstack_order_map[ksign].size() && hist_map[{ksign, thstack_order_map[ksign].at(first_ind)}]->Integral("width") < min_integral_thrsh){
                first_ind++;
            }

            if (first_ind >= thstack_order_map[ksign].size()){
                std::cerr << "Trying to find first index of category with integral over minimum threshold, but NONE found!" << std::endl;
                std::cerr << "Return without plotting or saving histograms!" << std::endl;
            }

            hist_map[{ksign, thstack_order_map[ksign].at(first_ind)}]->GetYaxis()->SetRangeUser(0., ymax * 1.1);
    
            hist_map[{ksign, thstack_order_map[ksign].at(first_ind)}]->Draw("E");

            for (int grp_ind = first_ind+1; grp_ind < thstack_order_map[ksign].size(); grp_ind++){
                int kgrp = thstack_order_map[ksign].at(grp_ind);
                if (hist_map[{ksign, kgrp}]->Integral("width") > min_integral_thrsh){ // skipping categories with too low or zero crossx
                    hist_map[{ksign, kgrp}]->Draw("E,same");
                }
            }
        }else{
            hs[ksign]->Draw("hist");
            thstack_helper(hs[ksign], kin_title, "");
        }

        tbox->Draw();
        l->Draw();
    }

    // determine which directory to output png files to
    // if not needed for single-b analysis, output into others/ sub-directory
    auto it_obs = std::find(single_b_analysis_observables.begin(), single_b_analysis_observables.end(), kin1d);
    std::string sub_dir = (it_obs != single_b_analysis_observables.end())? "" : "others/";

    std::string outdir_name = pythia_path + "plots/" + subdir_name + with_data_resonance_cuts_dir + sub_dir;

    std::string drawing_option_suffix = "";
    if (staggered)  drawing_option_suffix += "_staggered";
    else if (norm_unity)    drawing_option_suffix += "_unity";

    std::string outfile_name = outdir_name + "pythia_" + kin1d + drawing_option_suffix + optional_suffix + with_data_resonance_cuts_suffix;

    c->SaveAs(Form("%s.png", outfile_name.c_str()));
    if (saveAsC) c->SaveAs(Form("%s.c", outfile_name.c_str()));

    c->Close();
    delete c;
}
