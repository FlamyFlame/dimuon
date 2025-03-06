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
// #include "string"
// #include <vector>
// #include "time.h"
// #include "struct_hist.h"


const int nMCmodes = 2;
const int nSigns = 2;
// const int nDphi = 2;
// const int nflavorCategories = 4;
const int nflavorCategories = 3;

// TH1D* h[nMCmodes][nSigns][nDphi][nflavorCategories];
TH2D* h2d;
TFile* f[nMCmodes];
// TFile* f;

  // KEY: TH1D h_pair_pt_ptlead_ratio_gs_ISR_one_hard_scatt_sign2;1  
  // KEY: TH1D h_pair_pt_ptlead_ratio_gs_FSR_sign2;1 
  // KEY: TH1D h_pair_pt_ptlead_ratio_FC_sign2;1 
  // KEY: TH1D h_pair_pt_ptlead_ratio_from_same_b_sign2;1


std::string mcdir = "/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample/";
std::vector<std::string> sub_dirs = {"bb_full_sample/", "cc_full_sample/"};

std::string fnames[nMCmodes] = {"histograms_mc_truth_bb_combined.root", "histograms_mc_truth_cc_combined.root"};
std::vector<std::string> mcmodes = {"_bb","_cc"};
std::vector<std::string> signs = {"_sign1", "_sign2"};
// std::vector<std::string> dphis = {"_near", "_away"};
// std::vector<std::string> flavor_categories = {"_cc", "_other_flavors", "_bb", "_single_b"};
std::vector<std::string> flavor_categories = {"_cc", "_bb", "_single_b"};
std::vector<std::string> flavor_catg_labels = {"cc", "bb", "single b"};
std::vector<Color_t> line_colors = {kRed, kBlue, kBlack};
std::vector<Color_t> fill_colors = {kRed, kBlue, kYellow};

std::string subpl_titles[nMCmodes][nSigns] = {{"POWHEG, bb, same sign", "POWHEG, bb, opposite sign"}, {"POWHEG, cc, same sign", "POWHEG, cc, opposite sign"}};
std::string subpl_titles_sum_bb_cc[nSigns] = {"POWHEG bb+cc, same sign", "POWHEG bb+cc, opposite sign"};

void hist_helper(TH1* h, std::string title, bool norm_unity){
  h->SetStats(0);
  if (norm_unity){
    h->Scale(1./h->Integral());
    h->GetYaxis()->SetTitle("pdf");
  }else{
    h->Scale(1./1000.,"width");
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
  // h->SetMarkerSize(1.2);
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

void plot_flavor_categorized_kinematic_single(std::string kin, bool projx_2d, bool projy_2d, bool staggered, bool norm_unity, std::string kin1d, std::string kin_title, bool logx=false){

  TH1D* h[nMCmodes][nSigns][nflavorCategories];
  THStack *hs[nMCmodes][nSigns];
  if (projx_2d && projy_2d){
    std::cout << "For single kinematic plotting, can only project either X or Y variable." << std::endl;
    throw std::exception();
  }

  if (staggered && norm_unity){
    std::cout << "Cannot both staggered and normalize to unity." << std::endl;
    throw std::exception();
  }

  TCanvas* c = new TCanvas("c1","c1",2900,2000);
  c->Divide(2,2);

  for (int imc = 0; imc < nMCmodes; imc++){
    f[imc] = TFile::Open((mcdir + sub_dirs[imc] + fnames[imc]).c_str());

    for (unsigned int ksign = 0; ksign < nSigns; ksign++){

      c->cd(imc * nMCmodes + ksign + 1);
      gPad->SetLeftMargin(0.16);
      // gPad->SetRightMargin(0.16);
      gPad->SetBottomMargin(0.135);
      gPad->SetLogx(logx);

      TLegend* l = new TLegend(0.72,0.68,0.9,0.9);
      // l->SetMarkerSize(2);
      l->SetBorderSize(0);
      l->SetFillStyle(0);
      l->SetTextFont(43);
      l->SetMargin(0.2);
      l->SetTextColor(1);
      // l->AddEntry("",Form("%d",imc * nMCmodes + ksign + 1),"");
      // l->Draw();

      if (staggered){
        hs[imc][ksign] = new THStack(("hs" + mcmodes[imc] + signs[ksign]).c_str(), ("hs" + mcmodes[imc] + signs[ksign]).c_str());
      }

      for (int mgrp = 0; mgrp < nflavorCategories; mgrp++){
        if (projx_2d){
          h2d = (TH2D*) f[imc]->Get(("h_" + kin + signs[ksign] + flavor_categories[mgrp]).c_str());
          if (norm_unity){
            h[imc][ksign][mgrp] = (TH1D*) h2d->ProjectionX(("h_unity_" + kin1d + mcmodes[imc] + signs[ksign] + flavor_categories[mgrp]).c_str());
          }else{
            h[imc][ksign][mgrp] = (TH1D*) h2d->ProjectionX(("h_" + kin1d + mcmodes[imc] + signs[ksign] + flavor_categories[mgrp]).c_str());
          }
        }
        else if (projy_2d){
          h2d = (TH2D*) f[imc]->Get(("h_" + kin + signs[ksign] + flavor_categories[mgrp]).c_str());
          if (norm_unity){
            h[imc][ksign][mgrp] = (TH1D*) h2d->ProjectionY(("h_unity_" + kin1d + mcmodes[imc] + signs[ksign] + flavor_categories[mgrp]).c_str());
          }else{
            h[imc][ksign][mgrp] = (TH1D*) h2d->ProjectionY(("h_" + kin1d + mcmodes[imc] + signs[ksign] + flavor_categories[mgrp]).c_str());
          }
        }else{ // 1D
          h[imc][ksign][mgrp] = (TH1D*) f[imc]->Get(("h_" + kin + signs[ksign] + flavor_categories[mgrp]).c_str());
          // std::cout << kin1d << ", " << subpl_titles[imc][ksign] << ", " << flavor_catg_labels[mgrp] << ", dphi " << ldphi << ", integral: " << h[imc][ksign][mgrp]->Integral() << std::endl;
        }

        if (!h[imc][ksign][mgrp]){
          std::cout << "file h_" << kin << flavor_categories[mgrp] << signs[ksign] << "does not exist." << std::endl;
          throw std::exception();
        }

        if (staggered){
          hist_helper(h[imc][ksign][mgrp], subpl_titles[imc][ksign] + ", accumulative", norm_unity);
        }else if(norm_unity){
          hist_helper(h[imc][ksign][mgrp], subpl_titles[imc][ksign] + ", unity", norm_unity);
        }else{
          hist_helper(h[imc][ksign][mgrp], subpl_titles[imc][ksign], norm_unity);
        }
        if (!staggered){
          h[imc][ksign][mgrp]->SetLineColor(line_colors[mgrp]);
          h[imc][ksign][mgrp]->SetMarkerColor(line_colors[mgrp]);          
        }else{
          hs[imc][ksign]->Add(h[imc][ksign][mgrp]);
          h[imc][ksign][mgrp]->SetFillColor(fill_colors[mgrp]);
        }
        // l->AddEntry(h[imc][ksign][mgrp], "hello","lp");
        if (staggered){
          l->AddEntry(h[imc][ksign][mgrp], flavor_catg_labels[mgrp].c_str(),"f");
        }else{
          l->AddEntry(h[imc][ksign][mgrp], flavor_catg_labels[mgrp].c_str(),"lp");
        }
      }

      float ymax = h[imc][ksign][0]->GetMaximum();
      if (!staggered){
        for (int i = 1; i < nflavorCategories; i++){
          ymax = (ymax > h[imc][ksign][i]->GetMaximum())? ymax : h[imc][ksign][i]->GetMaximum();
        }
        h[imc][ksign][0]->GetYaxis()->SetRangeUser(0., ymax * 1.1);

        h[imc][ksign][0]->Draw("E");
        for (int i = 1; i < nflavorCategories; i++){
          h[imc][ksign][i]->Draw("E,same");
        }
      }else{
        hs[imc][ksign]->Draw("hist");
        thstack_helper(hs[imc][ksign], kin_title, subpl_titles[imc][ksign]);
      }
      l->Draw();
    }
  }

  if (staggered){
    c->SaveAs(Form("/usatlas/u/yuhanguo/usatlasdata/plots/powheg/flavor_categoried/sum_bb_cc/powheg_%s_flavor_staggered.png", kin1d.c_str()));
  }else if (norm_unity){
    c->SaveAs(Form("/usatlas/u/yuhanguo/usatlasdata/plots/powheg/flavor_categoried/sum_bb_cc/powheg_%s_flavor_unity.png", kin1d.c_str()));
  }else{
    c->SaveAs(Form("/usatlas/u/yuhanguo/usatlasdata/plots/powheg/flavor_categoried/sum_bb_cc/powheg_%s_flavor.png", kin1d.c_str()));
  }
  c->Close();
  delete c;
}


void plot_flavor_categorized_kinematic_single_sum_bb_cc(std::string kin, bool projx_2d, bool projy_2d, bool staggered, bool norm_unity, std::string kin1d, std::string kin_title, bool logx=false){
  
  TH1D* h[nMCmodes][nSigns][nflavorCategories];
  THStack *hs_sum_bb_cc[nSigns];
  if (projx_2d && projy_2d){
    std::cout << "For single kinematic plotting, can only project either X or Y variable." << std::endl;
    throw std::exception();
  }

  if (staggered && norm_unity){
    std::cout << "Cannot both staggered and normalize to unity." << std::endl;
    throw std::exception();
  }

  TCanvas* c = new TCanvas("c1","c1",1250,500);
  c->Divide(2,1);

  for (int imc = 0; imc < nMCmodes; imc++){
    f[imc] = TFile::Open((mcdir + sub_dirs[imc] + fnames[imc]).c_str());
  }

  for (unsigned int ksign = 0; ksign < nSigns; ksign++){
    c->cd(ksign + 1);

    gPad->SetLeftMargin(0.16);
    // gPad->SetRightMargin(0.16);
    gPad->SetBottomMargin(0.135);
    gPad->SetLogx(logx);

    TLegend* l = new TLegend(0.72,0.68,0.9,0.9);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetTextFont(43);
    l->SetMargin(0.2);
    l->SetTextColor(1);
    // l->AddEntry("",Form("%d",imc * nMCmodes + ksign + 1),"");
    // l->Draw();

    if (staggered){
      hs_sum_bb_cc[ksign] = new THStack(("hs_sum_bb_cc" + signs[ksign]).c_str(), ("hs_sum_bb_cc" + signs[ksign]).c_str());
    }

    for (int mgrp = 0; mgrp < nflavorCategories; mgrp++){
      for (int imc = 0; imc < nMCmodes; imc++){
        if (projx_2d){
          h2d = (TH2D*) f[imc]->Get(("h_" + kin + signs[ksign] + flavor_categories[mgrp]).c_str());
          if (norm_unity){
            h[imc][ksign][mgrp] = (TH1D*) h2d->ProjectionX(("h_unity_" + kin1d + mcmodes[imc] + signs[ksign] + flavor_categories[mgrp]).c_str());
          }else{
            h[imc][ksign][mgrp] = (TH1D*) h2d->ProjectionX(("h_" + kin1d + mcmodes[imc] + signs[ksign] + flavor_categories[mgrp]).c_str());
          }
        }
        else if (projy_2d){
          h2d = (TH2D*) f[imc]->Get(("h_" + kin + signs[ksign] + flavor_categories[mgrp]).c_str());
          if (norm_unity){
            h[imc][ksign][mgrp] = (TH1D*) h2d->ProjectionY(("h_unity_" + kin1d + mcmodes[imc] + signs[ksign] + flavor_categories[mgrp]).c_str());
          }else{
            h[imc][ksign][mgrp] = (TH1D*) h2d->ProjectionY(("h_" + kin1d + mcmodes[imc] + signs[ksign] + flavor_categories[mgrp]).c_str());
          }
        }else{ // 1D
          h[imc][ksign][mgrp] = (TH1D*) f[imc]->Get(("h_" + kin + signs[ksign] + flavor_categories[mgrp]).c_str());
        }

        if (!h[imc][ksign][mgrp]){
          std::cout << "file h_" << kin << flavor_categories[mgrp] << signs[ksign] << "does not exist." << std::endl;
          throw std::exception();
        }
      }

      h[0][ksign][mgrp]->Add(h[1][ksign][mgrp]);
      if(norm_unity){
        hist_helper(h[0][ksign][mgrp], subpl_titles_sum_bb_cc[ksign] + ", unity", norm_unity);
      }else{
        hist_helper(h[0][ksign][mgrp], subpl_titles_sum_bb_cc[ksign], norm_unity);
      }
      if (!staggered){
        h[0][ksign][mgrp]->SetLineColor(line_colors[mgrp]);
        h[0][ksign][mgrp]->SetMarkerColor(line_colors[mgrp]);          
      }else{
        hs_sum_bb_cc[ksign]->Add(h[0][ksign][mgrp]);
        h[0][ksign][mgrp]->SetFillColor(fill_colors[mgrp]);
      }
      if (staggered){
        l->AddEntry(h[0][ksign][mgrp], flavor_catg_labels[mgrp].c_str(),"f");
      }else{
        l->AddEntry(h[0][ksign][mgrp], flavor_catg_labels[mgrp].c_str(),"lp");
      }
    }

    float ymax = h[0][ksign][0]->GetMaximum();
    if (!staggered){
      for (int i = 1; i < nflavorCategories; i++){
        ymax = (ymax > h[0][ksign][i]->GetMaximum())? ymax : h[0][ksign][i]->GetMaximum();
      }
      h[0][ksign][0]->GetYaxis()->SetRangeUser(0., ymax * 1.1);

      h[0][ksign][0]->Draw("E");
      for (int i = 1; i < nflavorCategories; i++){
        h[0][ksign][i]->Draw("E,same");
      }
    }else{
      hs_sum_bb_cc[ksign]->Draw("hist");
      thstack_helper(hs_sum_bb_cc[ksign], kin_title, subpl_titles_sum_bb_cc[ksign]);
    }
    l->Draw();
  }

  if (staggered){
    c->SaveAs(Form("/usatlas/u/yuhanguo/usatlasdata/plots/powheg/flavor_categoried/sum_bb_cc/powheg_%s_flavor_staggered_sum_bb_cc.png", kin1d.c_str()));
  }else if (norm_unity){
    c->SaveAs(Form("/usatlas/u/yuhanguo/usatlasdata/plots/powheg/flavor_categoried/sum_bb_cc/powheg_%s_flavor_unity_sum_bb_cc.png", kin1d.c_str()));
  }else{
    c->SaveAs(Form("/usatlas/u/yuhanguo/usatlasdata/plots/powheg/flavor_categoried/sum_bb_cc/powheg_%s_flavor_sum_bb_cc.png", kin1d.c_str()));
  }
  c->Close();
  delete c;
}


void plot_flavor_categorized_kinematics_powheg(){
  // plot_flavor_categorized_kinematic_single("DR", false, false, false, false, "DR", "#Delta R");
  // plot_flavor_categorized_kinematic_single("DR", false, false, true, false, "DR", "#Delta R"); // accumulative
  // plot_flavor_categorized_kinematic_single("DR", false, false, false, true, "DR", "#Delta R"); // norm to unity
  plot_flavor_categorized_kinematic_single_sum_bb_cc("DR", false, false, false, false, "DR", "#Delta R");
  plot_flavor_categorized_kinematic_single_sum_bb_cc("DR", false, false, true, false, "DR", "#Delta R"); // accumulative
  plot_flavor_categorized_kinematic_single_sum_bb_cc("DR", false, false, false, true, "DR", "#Delta R");

  // plot_flavor_categorized_kinematic_single("pt_asym", false, false, false, false, "pt_asym","A");
  // plot_flavor_categorized_kinematic_single("pt_asym", false, false, true, false, "pt_asym","A"); // accumulative
  // plot_flavor_categorized_kinematic_single("pt_asym", false, false, false, true, "pt_asym","A"); // norm to unity
  plot_flavor_categorized_kinematic_single_sum_bb_cc("pt_asym", false, false, false, false, "pt_asym","A");
  plot_flavor_categorized_kinematic_single_sum_bb_cc("pt_asym", false, false, true, false, "pt_asym","A"); // accumulative
  plot_flavor_categorized_kinematic_single_sum_bb_cc("pt_asym", false, false, false, true, "pt_asym","A");
    
  // plot_flavor_categorized_kinematic_single("pair_pt_ptlead_ratio", false, false, false, false, "pair_pt_ptlead_ratio", "p_{T}^{pair}/p_{T}^{lead}");
  // plot_flavor_categorized_kinematic_single("pair_pt_ptlead_ratio", false, false, true, false, "pair_pt_ptlead_ratio", "p_{T}^{pair}/p_{T}^{lead}");
  // plot_flavor_categorized_kinematic_single("pair_pt_ptlead_ratio", false, false, false, true, "pair_pt_ptlead_ratio", "p_{T}^{pair}/p_{T}^{lead}"); // norm to unity
  plot_flavor_categorized_kinematic_single_sum_bb_cc("pair_pt_ptlead_ratio", false, false, false, false, "pair_pt_ptlead_ratio", "p_{T}^{pair}/p_{T}^{lead}");
  plot_flavor_categorized_kinematic_single_sum_bb_cc("pair_pt_ptlead_ratio", false, false, true, false, "pair_pt_ptlead_ratio", "p_{T}^{pair}/p_{T}^{lead}");
  plot_flavor_categorized_kinematic_single_sum_bb_cc("pair_pt_ptlead_ratio", false, false, false, true, "pair_pt_ptlead_ratio", "p_{T}^{pair}/p_{T}^{lead}");
  
  // // plot_flavor_categorized_kinematic_single("mQQ", false, false, false, false, "mQQ","m_{QQ}");
  // // plot_flavor_categorized_kinematic_single("mQQ", false, false, true, false, "mQQ","m_{QQ}");
  // // plot_flavor_categorized_kinematic_single("mQQ", false, false, false, true, "mQQ","m_{QQ}"); // norm to unity
  // // plot_flavor_categorized_kinematic_single_sum_bb_cc("mQQ", false, false, false, false, "mQQ","m_{QQ}");
  // // plot_flavor_categorized_kinematic_single_sum_bb_cc("mQQ", false, false, true, false, "mQQ","m_{QQ}");
  // // plot_flavor_categorized_kinematic_single_sum_bb_cc("mQQ", false, false, false, true, "mQQ","m_{QQ}"); // norm to unity

  // // plot_flavor_categorized_kinematic_single("mQQ_Q_ratio", false, false, false, false, "mQQ_Q_ratio", "m_{QQ} / Q");
  // // plot_flavor_categorized_kinematic_single("mQQ_Q_ratio", false, false, true, false, "mQQ_Q_ratio", "m_{QQ} / Q");
  // // plot_flavor_categorized_kinematic_single("mQQ_Q_ratio", false, false, false, true, "mQQ_Q_ratio", "m_{QQ} / Q"); // norm to unity
  // // plot_flavor_categorized_kinematic_single_sum_bb_cc("mQQ_Q_ratio", false, false, false, false, "mQQ_Q_ratio", "m_{QQ} / Q");
  // // plot_flavor_categorized_kinematic_single_sum_bb_cc("mQQ_Q_ratio", false, false, true, false, "mQQ_Q_ratio", "m_{QQ} / Q");
  // // plot_flavor_categorized_kinematic_single_sum_bb_cc("mQQ_Q_ratio", false, false, false, true, "mQQ_Q_ratio", "m_{QQ} / Q"); // norm to unity
  
  // // plot_flavor_categorized_kinematic_single("mQQ_mHard_ratio", false, false, false, false, "mQQ_mHard_ratio","m_{QQ} / #sqrt{#hat{s}}");
  // // plot_flavor_categorized_kinematic_single("mQQ_mHard_ratio", false, false, true, false, "mQQ_mHard_ratio","m_{QQ} / #sqrt{#hat{s}}");
  // // plot_flavor_categorized_kinematic_single("mQQ_mHard_ratio", false, false, false, true, "mQQ_mHard_ratio","m_{QQ} / #sqrt{#hat{s}}"); // norm to unity
  // // plot_flavor_categorized_kinematic_single_sum_bb_cc("mQQ_mHard_ratio", false, false, false, false, "mQQ_mHard_ratio","m_{QQ} / #sqrt{#hat{s}}");
  // // plot_flavor_categorized_kinematic_single_sum_bb_cc("mQQ_mHard_ratio", false, false, true, false, "mQQ_mHard_ratio","m_{QQ} / #sqrt{#hat{s}}");
  // // plot_flavor_categorized_kinematic_single_sum_bb_cc("mQQ_mHard_ratio", false, false, false, true, "mQQ_mHard_ratio","m_{QQ} / #sqrt{#hat{s}}"); // norm to unity

  // plot_flavor_categorized_kinematic_single("ptlead_pair_pt", true, false, false, false, "pair_pt", "p_{T}^{pair}");
  // plot_flavor_categorized_kinematic_single("ptlead_pair_pt", true, false, true, false, "pair_pt", "p_{T}^{pair}");
  // plot_flavor_categorized_kinematic_single("ptlead_pair_pt", true, false, false, true, "pair_pt", "p_{T}^{pair}"); // norm to unity
  plot_flavor_categorized_kinematic_single_sum_bb_cc("ptlead_pair_pt", true, false, false, false, "pair_pt", "p_{T}^{pair}");
  plot_flavor_categorized_kinematic_single_sum_bb_cc("ptlead_pair_pt", true, false, true, false, "pair_pt", "p_{T}^{pair}");
  plot_flavor_categorized_kinematic_single_sum_bb_cc("ptlead_pair_pt", true, false, false, true, "pair_pt", "p_{T}^{pair}");

  // plot_flavor_categorized_kinematic_single("ptlead_pair_pt", false, true, false, false, "ptlead", "p_{T}^{lead}");
  // plot_flavor_categorized_kinematic_single("ptlead_pair_pt", false, true, true, false, "ptlead", "p_{T}^{lead}");
  // plot_flavor_categorized_kinematic_single("ptlead_pair_pt", false, true, false, true, "ptlead", "p_{T}^{lead}"); // norm to unity
  plot_flavor_categorized_kinematic_single_sum_bb_cc("ptlead_pair_pt", false, true, false, false, "ptlead", "p_{T}^{lead}");
  plot_flavor_categorized_kinematic_single_sum_bb_cc("ptlead_pair_pt", false, true, true, false, "ptlead", "p_{T}^{lead}");
  plot_flavor_categorized_kinematic_single_sum_bb_cc("ptlead_pair_pt", false, true, false, true, "ptlead", "p_{T}^{lead}");
  
  // plot_flavor_categorized_kinematic_single ("minv_pair_pt", false, true, false, false, "minv", "m_{#mu#mu}");
  // plot_flavor_categorized_kinematic_single("minv_pair_pt", false, true, true, false, "minv", "m_{#mu#mu}");
  // plot_flavor_categorized_kinematic_single("minv_pair_pt", false, true, false, true, "minv", "m_{#mu#mu}"); // norm to unity
  plot_flavor_categorized_kinematic_single_sum_bb_cc ("minv_pair_pt", false, true, false, false, "minv", "m_{#mu#mu}");
  plot_flavor_categorized_kinematic_single_sum_bb_cc("minv_pair_pt", false, true, true, false, "minv", "m_{#mu#mu}");
  plot_flavor_categorized_kinematic_single_sum_bb_cc ("minv_pair_pt", false, true, false, true, "minv", "m_{#mu#mu}");
  
  // plot_flavor_categorized_kinematic_single("Deta_Dphi", true, false, false, false, "Dphi", "#Delta #phi");
  // plot_flavor_categorized_kinematic_single("Deta_Dphi", true, false, true, false, "Dphi", "#Delta #phi");
  // plot_flavor_categorized_kinematic_single("Deta_Dphi", true, false, false, true, "Dphi", "#Delta #phi"); // norm to unity
  plot_flavor_categorized_kinematic_single_sum_bb_cc("Deta_Dphi", true, false, false, false, "Dphi", "#Delta #phi");
  plot_flavor_categorized_kinematic_single_sum_bb_cc("Deta_Dphi", true, false, true, false, "Dphi", "#Delta #phi");
  plot_flavor_categorized_kinematic_single_sum_bb_cc("Deta_Dphi", true, false, false, true, "Dphi", "#Delta #phi");
  
  // plot_flavor_categorized_kinematic_single("Deta_Dphi", false, true, false, false, "Deta", "#Delta #eta");
  // plot_flavor_categorized_kinematic_single("Deta_Dphi", false, true, true, false, "Deta", "#Delta #eta");
  // plot_flavor_categorized_kinematic_single("Deta_Dphi", false, true, false, true, "Deta", "#Delta #eta"); // norm to unity
  plot_flavor_categorized_kinematic_single_sum_bb_cc("Deta_Dphi", false, true, false, false, "Deta", "#Delta #eta");
  plot_flavor_categorized_kinematic_single_sum_bb_cc("Deta_Dphi", false, true, true, false, "Deta", "#Delta #eta");
  plot_flavor_categorized_kinematic_single_sum_bb_cc("Deta_Dphi", false, true, false, true, "Deta", "#Delta #eta");
}

