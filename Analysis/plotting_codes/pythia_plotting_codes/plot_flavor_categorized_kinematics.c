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
// #include "string"
// #include <vector>
// #include "time.h"
// #include "struct_hist.h"


const int nMCmodes = 2;
const int nSigns = 2;
// const int nDphi = 2;
const int nflavorCategories = 5;
// const int nflavorCategories = 3;

// TH1D* h[nSigns][nDphi][nflavorCategories];
TH2D* h2d;
TFile* f;

bool with_data_resonance_cuts = false;
std::string with_data_resonance_cuts_dir = with_data_resonance_cuts? "with_data_resonance_cuts/" : "no_data_resonance_cuts/";
std::string with_data_resonance_cuts_suffix = with_data_resonance_cuts? "_with_data_resonance_cuts" : "_no_data_resonance_cuts";

std::string pythia_path = "/usatlas/u/yuhanguo/usatlasdata/pythia/";
std::string fname;
std::vector<std::string> signs = {"_sign1", "_sign2"};
std::vector<std::string> flavor_categories = {"_cc", "_other_flavors", "_bb", "_single_b", "_resonance"};
std::vector<std::string> flavor_catg_labels = {"cc", "other flavors", "bb", "single b", "resonances"};
// std::vector<std::string> flavor_catg_labels = {"cc", "bb", "single b"};
std::vector<Color_t> line_colors = {kRed, kGray, kBlue, kOrange, kCyan+1};
std::vector<Color_t> fill_colors = {kRed, kGray, kBlue, kYellow, kCyan+1};

// std::string subpl_titles[nSigns] = {{"bb, same sign", "bb, opposite sign"}, {"cc, same sign", "cc, opposite sign"}};
std::string subpl_titles[nSigns] = {"pythia, same sign", "pythia, opposite sign"};

void initialize(){
  fname = "histograms_pythia_combined" + with_data_resonance_cuts_suffix + ".root";
}

void hist_helper(TH1* h, std::string title, bool norm_unity){
  h->SetStats(0);
  if (norm_unity){
    h->Scale(1./h->Integral());
    h->GetYaxis()->SetTitle("pdf");
  }else{
    h->Scale(1000.,"width");
  }
  h->SetTitle(title.c_str());
  h->GetYaxis()->SetLabelFont(43);
  h->GetYaxis()->SetLabelSize(23);
  h->GetYaxis()->SetLabelOffset(0.01);
  h->GetYaxis()->SetTitleFont(43);
  h->GetYaxis()->SetTitleSize(25);
  h->GetYaxis()->SetTitleOffset(2);
  h->GetXaxis()->SetLabelFont(43);
  h->GetXaxis()->SetLabelSize(23);
  h->GetXaxis()->SetLabelOffset(0.01);
  h->SetMarkerStyle(20);
  h->GetYaxis()->SetRangeUser(0.,h->GetMaximum() * 1.1);
  h->GetXaxis()->SetTitleFont(43);  
  h->GetXaxis()->SetTitleSize(25);
  h->GetXaxis()->SetTitleOffset(1);
}

void thstack_helper(THStack* h, std::string kin_name, std::string title){
  h->SetTitle(title.c_str());
  h->GetXaxis()->SetTitle(kin_name.c_str());
  h->GetYaxis()->SetTitle(("d#sigma/d" + kin_name + " [nb]").c_str());

  h->GetYaxis()->SetLabelFont(43);
  h->GetYaxis()->SetLabelSize(23);
  h->GetYaxis()->SetLabelOffset(0.01);
  h->GetYaxis()->SetTitleFont(43);
  h->GetYaxis()->SetTitleSize(25);
  h->GetYaxis()->SetTitleOffset(2);

  h->GetXaxis()->SetLabelFont(43);
  h->GetXaxis()->SetLabelSize(23);
  h->GetXaxis()->SetLabelOffset(0.01);
  h->GetXaxis()->SetTitleFont(43);  
  h->GetXaxis()->SetTitleSize(25);
  h->GetXaxis()->SetTitleOffset(1);
}

void plot_flavor_categorized_kinematic_single(std::string kin, bool projx_2d, bool projy_2d, bool staggered, bool norm_unity, std::string kin1d, std::string kin_title, std::vector<std::array<float,2>> cuts = {}, bool logx=false){

  initialize();

  TH1D* h[nSigns][nflavorCategories];
  THStack *hs[nSigns];
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
  for (unsigned int ksign = 0; ksign < nSigns; ksign++){

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

    for (int mgrp = 0; mgrp < nflavorCategories; mgrp++){
      if (projx_2d){
        h2d = (TH2D*) f->Get(("h_" + kin + signs[ksign] + flavor_categories[mgrp]).c_str());
        if (norm_unity){
          h[ksign][mgrp] = (TH1D*) h2d->ProjectionX(("h_unity_" + kin1d + signs[ksign] + flavor_categories[mgrp]).c_str());
        }else{
          h[ksign][mgrp] = (TH1D*) h2d->ProjectionX(("h_" + kin1d + signs[ksign] + flavor_categories[mgrp]).c_str());
        }
      }
      else if (projy_2d){
        h2d = (TH2D*) f->Get(("h_" + kin + signs[ksign] + flavor_categories[mgrp]).c_str());
        if (norm_unity){
          h[ksign][mgrp] = (TH1D*) h2d->ProjectionY(("h_unity_" + kin1d + signs[ksign] + flavor_categories[mgrp]).c_str());
        }else{
          h[ksign][mgrp] = (TH1D*) h2d->ProjectionY(("h_" + kin1d + signs[ksign] + flavor_categories[mgrp]).c_str());
        }
      }else{ // 1D
        h[ksign][mgrp] = (TH1D*) f->Get(("h_" + kin + signs[ksign] + flavor_categories[mgrp]).c_str());
        // std::cout << kin1d << ", " << subpl_titles[ksign] << ", " << flavor_catg_labels[mgrp] << ", dphi " << ldphi << ", integral: " << h[ksign][mgrp]->Integral() << std::endl;
      }

      // std::cout << kin1d << ", " << subpl_titles[ksign] << ", " << flavor_catg_labels[mgrp] << ", integral: " << h[ksign][mgrp]->Integral() << std::endl;

      if (!h[ksign][mgrp]){
        std::cout << "file h_" << kin << flavor_categories[mgrp] << signs[ksign] << "does not exist." << std::endl;
        throw std::exception();
      }

      // if the argument cuts is not empty, do not draw the bins overlapping with the cuts
      if (!cuts.empty()) ApplyCutsTo1DHistogram(h[ksign][mgrp], cuts);

      if (staggered){
        hist_helper(h[ksign][mgrp], subpl_titles[ksign] + ", accumulative", norm_unity);
      }else if(norm_unity){
        hist_helper(h[ksign][mgrp], subpl_titles[ksign] + ", unity", norm_unity);
      }else{
        hist_helper(h[ksign][mgrp], subpl_titles[ksign], norm_unity);
      }
      if (!staggered){
        h[ksign][mgrp]->SetLineColor(line_colors[mgrp]);
        h[ksign][mgrp]->SetMarkerColor(line_colors[mgrp]);          
      }else{
        hs[ksign]->Add(h[ksign][mgrp]);
        h[ksign][mgrp]->SetFillColor(fill_colors[mgrp]);
      }
      // l->AddEntry(h[ksign][mgrp], "hello","lp");
      if (staggered){
        l->AddEntry(h[ksign][mgrp], flavor_catg_labels[mgrp].c_str(),"f");
      }else{
        l->AddEntry(h[ksign][mgrp], flavor_catg_labels[mgrp].c_str(),"lp");
      }
    }

    float ymax = h[ksign][0]->GetMaximum();
    if (!staggered){
      for (int i = 1; i < nflavorCategories; i++){
        ymax = (ymax > h[ksign][i]->GetMaximum())? ymax : h[ksign][i]->GetMaximum();
      }
      h[ksign][0]->GetYaxis()->SetRangeUser(0., ymax * 1.1);

      h[ksign][0]->Draw("E");
      for (int i = 1; i < nflavorCategories; i++){
        h[ksign][i]->Draw("E,same");
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

  if (staggered){
    c->SaveAs(Form("%splots/flavor_categoried/%spythia_%s_flavor_staggered%s.png", pythia_path.c_str(), with_data_resonance_cuts_dir.c_str(), kin1d.c_str(), with_data_resonance_cuts_suffix.c_str()));
  }else if (norm_unity){
    c->SaveAs(Form("%splots/flavor_categoried/%spythia_%s_flavor_unity%s.png", pythia_path.c_str(), with_data_resonance_cuts_dir.c_str(), kin1d.c_str(), with_data_resonance_cuts_suffix.c_str()));
  }else{
    c->SaveAs(Form("%splots/flavor_categoried/%spythia_%s_flavor%s.png", pythia_path.c_str(), with_data_resonance_cuts_dir.c_str(), kin1d.c_str(), with_data_resonance_cuts_suffix.c_str()));
  }
  c->Close();
  delete c;

}


void plot_flavor_categorized_kinematics(){
  ParamsSet pms;

  plot_flavor_categorized_kinematic_single("DR", false, false, false, false, "DR", "#Delta R");
  plot_flavor_categorized_kinematic_single("DR", false, false, true, false, "DR", "#Delta R"); // accumulative
  // // plot_flavor_categorized_kinematic_single("DR", false, false, false, true, "DR", "#Delta R"); // norm to unity

  plot_flavor_categorized_kinematic_single("DR_zoomin", false, false, false, false, "DR_zoomin", "#Delta R");
  plot_flavor_categorized_kinematic_single("DR_zoomin", false, false, true, false, "DR_zoomin", "#Delta R"); // accumulative
  // // plot_flavor_categorized_kinematic_single("DR_zoomin", false, false, false, true, "DR_zoomin", "#Delta R"); // norm to unity

  plot_flavor_categorized_kinematic_single("DR_jacobian_corrected", false, false, false, false, "DR_jacobian_corrected", "#Delta R");
  plot_flavor_categorized_kinematic_single("DR_jacobian_corrected", false, false, true, false, "DR_jacobian_corrected", "#Delta R"); // accumulative
  // // plot_flavor_categorized_kinematic_single("DR_jacobian_corrected", false, false, false, true, "DR_jacobian_corrected", "#Delta R"); // norm to unity

  plot_flavor_categorized_kinematic_single("DR_zoomin_jacobian_corrected", false, false, false, false, "DR_zoomin_jacobian_corrected", "#Delta R");
  plot_flavor_categorized_kinematic_single("DR_zoomin_jacobian_corrected", false, false, true, false, "DR_zoomin_jacobian_corrected", "#Delta R"); // accumulative
  // // plot_flavor_categorized_kinematic_single("DR_zoomin_jacobian_corrected", false, false, false, true, "DR_zoomin_jacobian_corrected", "#Delta R"); // norm to unity

  plot_flavor_categorized_kinematic_single("pt_asym", false, false, false, false, "pt_asym", "A");
  plot_flavor_categorized_kinematic_single("pt_asym", false, false, true, false, "pt_asym", "A"); // accumulative
  // plot_flavor_categorized_kinematic_single("pt_asym", false, false, false, true, "pt_asym", "A"); // norm to unity
  
  plot_flavor_categorized_kinematic_single("psrapidity_ordered_pt_asym", false, false, false, false, "psrapidity_ordered_pt_asym", "#Tilde{A}");
  plot_flavor_categorized_kinematic_single("psrapidity_ordered_pt_asym", false, false, true, false, "psrapidity_ordered_pt_asym", "#Tilde{A}"); // accumulative
  // plot_flavor_categorized_kinematic_single("psrapidity_ordered_pt_asym", false, false, false, true, "psrapidity_ordered_pt_asym", "#Tilde{A}"); // norm to unity
  
  plot_flavor_categorized_kinematic_single("pair_pt_ptlead_ratio", false, false, false, false, "pair_pt_ptlead_ratio", "p_{T}^{pair} / p_{T}^{lead}");
  plot_flavor_categorized_kinematic_single("pair_pt_ptlead_ratio", false, false, true, false, "pair_pt_ptlead_ratio", "p_{T}^{pair} / p_{T}^{lead}");
  // plot_flavor_categorized_kinematic_single("pair_pt_ptlead_ratio", false, false, false, true, "pair_pt_ptlead_ratio", "p_{T}^{pair} / p_{T}^{lead}"); // norm to unity
  
  // // plot_flavor_categorized_kinematic_single("Qsplit", false, false, false, false, "Qsplit", "Q_{split}");
  // // plot_flavor_categorized_kinematic_single("Qsplit", false, false, true, false, "Qsplit", "Q_{split}");
  // // // plot_flavor_categorized_kinematic_single("Qsplit", false, false, false, true, "Qsplit", "Q_{split}"); // norm to unity
  
  // // plot_flavor_categorized_kinematic_single("Qsplit_pTHat_ratio", false, false, false, false, "Qsplit_pTHat_ratio", "Q_{split} / #hat{p}_T");
  // // plot_flavor_categorized_kinematic_single("Qsplit_pTHat_ratio", false, false, true, false, "Qsplit_pTHat_ratio", "Q_{split} / #hat{p}_T");
  // // // plot_flavor_categorized_kinematic_single("Qsplit_pTHat_ratio", false, false, false, true, "Qsplit_pTHat_ratio", "Q_{split} / #hat{p}_T"); // norm to unity
  
  // // plot_flavor_categorized_kinematic_single("Qsplit_mHat_ratio", false, false, false, false, "Qsplit_mHat_ratio", "Q_{split} / #sqrt{#hat{s}}");
  // // plot_flavor_categorized_kinematic_single("Qsplit_mHat_ratio", false, false, true, false, "Qsplit_mHat_ratio", "Q_{split} / #sqrt{#hat{s}}");
  // // // plot_flavor_categorized_kinematic_single("Qsplit_mHat_ratio", false, false, false, true, "Qsplit_mHat_ratio", "Q_{split} / #sqrt{#hat{s}}"); // norm to unity

  plot_flavor_categorized_kinematic_single("ptlead_pair_pt", true, false, false, false, "pair_pt", "p_{T}^{pair}");
  plot_flavor_categorized_kinematic_single("ptlead_pair_pt", true, false, true, false, "pair_pt", "p_{T}^{pair}");
  // plot_flavor_categorized_kinematic_single("ptlead_pair_pt", true, false, false, true, "pair_pt", "p_{T}^{pair}"); // norm to unity

  plot_flavor_categorized_kinematic_single("ptlead_pair_pt", false, true, false, false, "ptlead", "p_{T}^{lead}");
  plot_flavor_categorized_kinematic_single("ptlead_pair_pt", false, true, true, false, "ptlead", "p_{T}^{lead}");
  // plot_flavor_categorized_kinematic_single("ptlead_pair_pt", false, true, false, true, "ptlead", "p_{T}^{lead}"); // norm to unity
  
  std::vector<std::array<float,2>> minv_cuts;
  if (with_data_resonance_cuts) minv_cuts = pms.minv_cuts;

  plot_flavor_categorized_kinematic_single("minv_pair_pt", false, true, false, false, "minv", "m_{#mu#mu}", minv_cuts);
  plot_flavor_categorized_kinematic_single("minv_pair_pt", false, true, true, false, "minv", "m_{#mu#mu}", minv_cuts);
  // plot_flavor_categorized_kinematic_single("minv_pair_pt", false, true, false, true, "minv", "m_{#mu#mu}", minv_cuts); // norm to unity

  plot_flavor_categorized_kinematic_single("minv_pair_pt_zoomin", false, true, false, false, "minv_zoomin", "m_{#mu#mu}", minv_cuts);
  plot_flavor_categorized_kinematic_single("minv_pair_pt_zoomin", false, true, true, false, "minv_zoomin", "m_{#mu#mu}", minv_cuts);
  // plot_flavor_categorized_kinematic_single("minv_pair_pt_zoomin", false, true, false, true, "minv_zoomin", "m_{#mu#mu}", minv_cuts); // norm to unity

  plot_flavor_categorized_kinematic_single("minv_pair_pt_jacobian_corrected", false, true, false, false, "minv_jacobian_corrected", "m_{#mu#mu}", minv_cuts);
  plot_flavor_categorized_kinematic_single("minv_pair_pt_jacobian_corrected", false, true, true, false, "minv_jacobian_corrected", "m_{#mu#mu}", minv_cuts);
  // plot_flavor_categorized_kinematic_single("minv_pair_pt_jacobian_corrected", false, true, false, true, "minv_jacobian_corrected", "m_{#mu#mu}", minv_cuts); // norm to unity

  plot_flavor_categorized_kinematic_single("minv_pair_pt_zoomin_jacobian_corrected", false, true, false, false, "minv_zoomin_jacobian_corrected", "m_{#mu#mu}", minv_cuts);
  plot_flavor_categorized_kinematic_single("minv_pair_pt_zoomin_jacobian_corrected", false, true, true, false, "minv_zoomin_jacobian_corrected", "m_{#mu#mu}", minv_cuts);
  // plot_flavor_categorized_kinematic_single("minv_pair_pt_zoomin_jacobian_corrected", false, true, false, true, "minv_zoomin_jacobian_corrected", "m_{#mu#mu}", minv_cuts); // norm to unity
  
  plot_flavor_categorized_kinematic_single("Deta_Dphi", true, false, false, false, "Dphi", "#Delta #phi");
  plot_flavor_categorized_kinematic_single("Deta_Dphi", true, false, true, false, "Dphi", "#Delta #phi");
  // plot_flavor_categorized_kinematic_single("Deta_Dphi", true, false, false, true, "Dphi", "#Delta #phi"); // norm to unity
  
  plot_flavor_categorized_kinematic_single("Deta_Dphi", false, true, false, false, "Deta", "#Delta #eta");
  plot_flavor_categorized_kinematic_single("Deta_Dphi", false, true, true, false, "Deta", "#Delta #eta");
  // plot_flavor_categorized_kinematic_single("Deta_Dphi", false, true, false, true, "Deta", "#Delta #eta"); // norm to unity
}



