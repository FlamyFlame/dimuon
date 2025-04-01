#include "PythonCategorizedPlottingBaseClass.h"

void PythonCategorizedPlottingBaseClass::initialize(){

    with_data_resonance_cuts_dir = with_data_resonance_cuts? "with_data_resonance_cuts/" : "no_data_resonance_cuts/";
    with_data_resonance_cuts_suffix = with_data_resonance_cuts? "_with_data_resonance_cuts" : "_no_data_resonance_cuts";
    // with_data_resonance_cuts_suffix = with_data_resonance_cuts? "_with_data_resonance_cuts" : "";
    
    fname = "histograms_pythia_combined" + with_data_resonance_cuts_suffix + ".root";

    signs[sign::same_sign]  = "_sign1";
    signs[sign::op_sign]    = "_sign2";

    subpl_titles[sign::same_sign]   = "pythia, same sign";
    subpl_titles[sign::op_sign]     = "pythia, opposite sign";

    fill_line_map();
    fill_thstack_order_map();  
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

    TCanvas* c = new TCanvas("c1","c1",1250,500);
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

        TLegend* l = new TLegend(0.63,0.58,0.89,0.9);
        l->SetBorderSize(0);
        l->SetFillStyle(0);
        l->SetTextFont(43);
        l->SetMargin(0.2);
        l->SetTextColor(1);

        if (staggered){
            hs[ksign] = new THStack(("hs" + signs[ksign]).c_str(), ("hs" + signs[ksign]).c_str());
        }

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
                if (staggered){
                    hist_helper(hist_map[{ksign, kgrp}], pow(10,6), norm_unity, subpl_titles[ksign] + ", accumulative");
                }else if(norm_unity){
                    hist_helper(hist_map[{ksign, kgrp}], 1., norm_unity, subpl_titles[ksign] + ", unity");
                }else{
                    hist_helper(hist_map[{ksign, kgrp}], pow(10,6), norm_unity, subpl_titles[ksign]);
                }

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
        }

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
            thstack_helper(hs[ksign], kin_title, subpl_titles[ksign]);
            // hs[ksign]->GetXaxis()->SetTitle(kin1d.c_str());
            // hs[ksign]->GetYaxis()->SetTitle(("d#sigma/d" + kin1d).c_str());
            // hs[ksign]->SetTitle(subpl_titles[ksign].c_str());          
        }

        l->Draw();
    }

    // determine which directory to output png files to
    // if not needed for single-b analysis, output into others/ sub-directory
    auto it_obs = std::find(single_b_analysis_observables.begin(), single_b_analysis_observables.end(), kin1d);
    std::string sub_dir = (it_obs != single_b_analysis_observables.end())? "" : "others/";

    if (staggered){
        c->SaveAs(Form("%splots/%s%s%spythia_%s_staggered%s%s.png", pythia_path.c_str(), subdir_name.c_str(), with_data_resonance_cuts_dir.c_str(), sub_dir.c_str(), kin1d.c_str(), optional_suffix.c_str(), with_data_resonance_cuts_suffix.c_str()));
    }else if (norm_unity){
        c->SaveAs(Form("%splots/%s%s%spythia_%s_unity%s%s.png", pythia_path.c_str(), subdir_name.c_str(), with_data_resonance_cuts_dir.c_str(), sub_dir.c_str(), kin1d.c_str(), optional_suffix.c_str(), with_data_resonance_cuts_suffix.c_str()));
    }else{
        c->SaveAs(Form("%splots/%s%s%spythia_%s%s%s.png", pythia_path.c_str(), subdir_name.c_str(), with_data_resonance_cuts_dir.c_str(), sub_dir.c_str(), kin1d.c_str(), optional_suffix.c_str(), with_data_resonance_cuts_suffix.c_str()));
    }

    c->Close();
    delete c;

}
