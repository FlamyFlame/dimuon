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
const int nDphi = 2;
const int nAncestorGroups = 4;

// TH1D* h[nMCmodes][nSigns][nDphi][nAncestorGroups];
TH2D* h2d;
TFile* f[nMCmodes];

std::string mc_path = "/usatlas/u/yuhanguo/usatlasdata/athena/runMCV2/";
std::string fnames[nMCmodes] = {"muon_pairs_mc_truth_bb.root", "muon_pairs_mc_truth_cc.root"};
std::vector<std::string> mcmodes = {"_bb","_cc"};
std::vector<std::string> signs = {"_sign1", "_sign2"};
std::vector<std::string> dphis = {"_near", "_away"};
std::vector<std::string> ancestor_grps = {"_gg", "_qg","_single_g","_qq"};
std::vector<std::string> ancestor_grp_labels = {"gg", "qg","single g","qq"};
std::vector<Color_t> line_colors = {kBlack, kRed, kGreen + 2, kBlue};
std::vector<Color_t> fill_colors = {kYellow, kRed, kGreen, kBlue};

std::string subpl_titles[nMCmodes][nSigns] = {{"bb, same sign", "bb, opposite sign"}, {"cc, same sign", "cc, opposite sign"}};

void hist_helper(TH1* h, std::string title, bool norm_unity){
  h->SetStats(0);
  if (norm_unity){
    h->Scale(1./h->Integral());
    h->GetYaxis()->SetTitle("pdf");
  }else{
    h->Scale(1.,"width");
  }
  h->SetTitle(title.c_str());
  h->GetYaxis()->SetLabelFont(43);
  h->GetYaxis()->SetLabelSize(38);
  h->GetYaxis()->SetLabelOffset(0.01);
  h->GetYaxis()->SetTitleFont(43);
  h->GetYaxis()->SetTitleSize(32);
  h->GetYaxis()->SetTitleOffset(2);
  h->GetXaxis()->SetLabelFont(43);
  h->GetXaxis()->SetLabelSize(38);
  h->GetXaxis()->SetLabelOffset(0.01);
  h->SetMarkerStyle(20);
  h->GetYaxis()->SetRangeUser(0.,h->GetMaximum() * 1.1);
  h->GetXaxis()->SetTitleFont(43);  
  h->GetXaxis()->SetTitleSize(32);
  h->GetXaxis()->SetTitleOffset(1);
}

void thstack_helper(THStack* h, std::string kin_name, std::string title){
  h->SetTitle(title.c_str());
  h->GetXaxis()->SetTitle(kin_name.c_str());
  h->GetYaxis()->SetTitle(("d#sigma/d" + kin_name).c_str());

  h->GetYaxis()->SetLabelFont(43);
  h->GetYaxis()->SetLabelSize(38);
  h->GetYaxis()->SetLabelOffset(0.01);
  h->GetYaxis()->SetTitleFont(43);
  h->GetYaxis()->SetTitleSize(45);
  h->GetYaxis()->SetTitleOffset(2);

  h->GetXaxis()->SetLabelFont(43);
  h->GetXaxis()->SetLabelSize(38);
  h->GetXaxis()->SetLabelOffset(0.01);
  h->GetXaxis()->SetTitleFont(43);  
  h->GetXaxis()->SetTitleSize(45);
  h->GetXaxis()->SetTitleOffset(1);
}

void plot_hard_scatt_single_kinematic(std::string kin, bool projx_2d, bool projy_2d, bool add, bool norm_unity, std::string kin1d, bool logx=false){

  TH1D* h[nMCmodes][nSigns][nDphi][nAncestorGroups];
  THStack *hs[nMCmodes][nSigns];
  if (projx_2d && projy_2d){
    std::cout << "For single kinematic plotting, can only project either X or Y variable." << std::endl;
    throw std::exception();
  }

  if (add && norm_unity){
    std::cout << "Cannot both add and normalize to unity." << std::endl;
    throw std::exception();
  }

  TCanvas* c = new TCanvas("c1","c1",2900,2000);
  c->Divide(2,2);

  for (int imc = 0; imc < nMCmodes; imc++){
    f[imc] = TFile::Open((mc_path + fnames[imc]).c_str());

    for (unsigned int ksign = 0; ksign < nSigns; ksign++){

      c->cd(imc * nMCmodes + ksign + 1);
      gPad->SetLeftMargin(0.16);
      // gPad->SetRightMargin(0.16);
      gPad->SetBottomMargin(0.135);
      gPad->SetLogx(logx);

      TLegend* l = new TLegend(0.75,0.68,0.89,0.9);
      l->SetBorderSize(0);
      l->SetFillStyle(0);
      l->SetTextFont(43);
      l->SetMargin(0.02);
      l->SetTextColor(1);
      // l->AddEntry("",Form("%d",imc * nMCmodes + ksign + 1),"");
      // l->Draw();

      if (add){
        hs[imc][ksign] = new THStack(("hs" + mcmodes[imc] + signs[ksign]).c_str(), ("hs" + mcmodes[imc] + signs[ksign]).c_str());
      }

      for (int mgrp = 0; mgrp < nAncestorGroups; mgrp++){
        for (int ldphi = 0; ldphi < nDphi; ldphi++){
          if (projx_2d){
            h2d = (TH2D*) f[imc]->Get(("h_QQ_" + kin + signs[ksign] + dphis[ldphi] + ancestor_grps[mgrp]).c_str());
            if (norm_unity){
              h[imc][ksign][ldphi][mgrp] = (TH1D*) h2d->ProjectionX(("h_QQ_unity_" + kin1d + signs[ksign] + dphis[ldphi] + ancestor_grps[mgrp]).c_str());
            }else{
              h[imc][ksign][ldphi][mgrp] = (TH1D*) h2d->ProjectionX(("h_QQ_" + kin1d + signs[ksign] + dphis[ldphi] + ancestor_grps[mgrp]).c_str());
            }
          }
          else if (projy_2d){
            h2d = (TH2D*) f[imc]->Get(("h_QQ_" + kin + signs[ksign] + dphis[ldphi] + ancestor_grps[mgrp]).c_str());
            if (norm_unity){
              h[imc][ksign][ldphi][mgrp] = (TH1D*) h2d->ProjectionY(("h_QQ_unity_" + kin1d + signs[ksign] + dphis[ldphi] + ancestor_grps[mgrp]).c_str());
            }else{
              h[imc][ksign][ldphi][mgrp] = (TH1D*) h2d->ProjectionY(("h_QQ_" + kin1d + signs[ksign] + dphis[ldphi] + ancestor_grps[mgrp]).c_str());
            }
          }else{ // 1D
            h[imc][ksign][ldphi][mgrp] = (TH1D*) f[imc]->Get(("h_QQ_" + kin + signs[ksign] + dphis[ldphi] + ancestor_grps[mgrp]).c_str());
            // std::cout << kin1d << ", " << subpl_titles[imc][ksign] << ", " << ancestor_grp_labels[mgrp] << ", dphi " << ldphi << ", integral: " << h[imc][ksign][ldphi][mgrp]->Integral() << std::endl;
          }
        }

        h[imc][ksign][0][mgrp]->Add(h[imc][ksign][1][mgrp]);
        // std::cout << kin1d << ", " << subpl_titles[imc][ksign] << ", " << ancestor_grp_labels[mgrp] << ", integral: " << h[imc][ksign][0][mgrp]->Integral() << std::endl;

        if (add){
          hist_helper(h[imc][ksign][0][mgrp], subpl_titles[imc][ksign] + ", accumulative", norm_unity);
        }else if(norm_unity){
          hist_helper(h[imc][ksign][0][mgrp], subpl_titles[imc][ksign] + ", unity", norm_unity);
        }else{
          hist_helper(h[imc][ksign][0][mgrp], subpl_titles[imc][ksign], norm_unity);
        }
        if (!add){
          h[imc][ksign][0][mgrp]->SetLineColor(line_colors[mgrp]);
          h[imc][ksign][0][mgrp]->SetMarkerColor(line_colors[mgrp]);          
        }else{
          hs[imc][ksign]->Add(h[imc][ksign][0][mgrp]);
          h[imc][ksign][0][mgrp]->SetFillColor(fill_colors[mgrp]);
        }
        // l->AddEntry(h[imc][ksign][0][mgrp], "hello","lp");
        if (add){
          l->AddEntry(h[imc][ksign][0][mgrp], ancestor_grp_labels[mgrp].c_str(),"f");
        }else{
          l->AddEntry(h[imc][ksign][0][mgrp], ancestor_grp_labels[mgrp].c_str(),"lp");
        }
      }

      float ymax;
      // if (add){
      //   hs[imc][ksign]->Add(h[imc][ksign][0][3]);
      //   hs[imc][ksign]->Add(h[imc][ksign][0][2]);
      //   hs[imc][ksign]->Add(h[imc][ksign][0][1]);
      //   hs[imc][ksign]->Add(h[imc][ksign][0][0]);
      //   // h[imc][ksign][0][2]->Add(h[imc][ksign][0][3]);
      //   // h[imc][ksign][0][1]->Add(h[imc][ksign][0][2]);
      //   // h[imc][ksign][0][0]->Add(h[imc][ksign][0][1]);
      // }else{
      if (!add){
        ymax = (h[imc][ksign][0][0]->GetMaximum() > h[imc][ksign][0][1]->GetMaximum())? h[imc][ksign][0][0]->GetMaximum() : h[imc][ksign][0][1]->GetMaximum();
        ymax = (ymax > h[imc][ksign][0][2]->GetMaximum())? ymax : h[imc][ksign][0][2]->GetMaximum();
        ymax = (ymax > h[imc][ksign][0][3]->GetMaximum())? ymax : h[imc][ksign][0][3]->GetMaximum();
        h[imc][ksign][0][0]->GetYaxis()->SetRangeUser(0., ymax * 1.1);
        h[imc][ksign][0][0]->Draw("E");
        h[imc][ksign][0][1]->Draw("E,same");
        h[imc][ksign][0][2]->Draw("E,same");
        h[imc][ksign][0][3]->Draw("E,same");
      }else{
        hs[imc][ksign]->Draw("hist");
        thstack_helper(hs[imc][ksign], kin1d, subpl_titles[imc][ksign]);
        // hs[imc][ksign]->GetXaxis()->SetTitle(kin1d.c_str());
        // hs[imc][ksign]->GetYaxis()->SetTitle(("d#sigma/d" + kin1d).c_str());
        // hs[imc][ksign]->SetTitle(subpl_titles[imc][ksign].c_str());          
      }
      l->Draw();
    }
  }

  if (add){
    c->SaveAs(Form("plots/powheg/hard_scatt/accumulative/hard_scatt_%s_accumulative.png", kin1d.c_str()));
  }else if (norm_unity){
    c->SaveAs(Form("plots/powheg/hard_scatt/hard_scatt_%s_unity.png", kin1d.c_str()));
  }else{
    c->SaveAs(Form("plots/powheg/hard_scatt/hard_scatt_%s.png", kin1d.c_str()));
  }
  c->Close();
  delete c;

}


void plot_hard_scatt(){
  // plot_hard_scatt_single_kinematic("DR", false, false, false, false, "DR");
  plot_hard_scatt_single_kinematic("DR", false, false, true, false, "DR"); // accumulative
  // plot_hard_scatt_single_kinematic("DR", false, false, false, true, "DR"); // norm to unity
  // plot_hard_scatt_single_kinematic("Dphi", false, false, false, false, "Dphi");
  plot_hard_scatt_single_kinematic("Dphi", false, false, true, false, "Dphi"); // accumulative
  // plot_hard_scatt_single_kinematic("Dphi", false, false, false, true, "Dphi"); // norm to unity
  // plot_hard_scatt_single_kinematic("minv", false, false, false, false, "minv", true);
  plot_hard_scatt_single_kinematic("minv", false, false, true, false, "minv", true); // accumulative
  // plot_hard_scatt_single_kinematic("minv", false, false, false, true, "minv", true); // norm to unity
  
  // plot_hard_scatt_single_kinematic("pair_pt_ptlead_ratio", false, false, false, false, "pair_pt_ptlead_ratio");
  plot_hard_scatt_single_kinematic("pair_pt_ptlead_ratio", false, false, true, false, "pair_pt_ptlead_ratio");
  // plot_hard_scatt_single_kinematic("pair_pt_ptlead_ratio", false, false, false, true, "pair_pt_ptlead_ratio"); // norm to unity
  // plot_hard_scatt_single_kinematic("asym", false, false, false, false, "asym");
  plot_hard_scatt_single_kinematic("asym", false, false, true, false, "asym");
  // plot_hard_scatt_single_kinematic("asym", false, false, false, true, "asym"); // norm to unity
  // plot_hard_scatt_single_kinematic("minv_s_cm_ratio", false, false, false, false, "minv_s_cm_ratio");
  plot_hard_scatt_single_kinematic("minv_s_cm_ratio", false, false, true, false, "minv_s_cm_ratio");
  // plot_hard_scatt_single_kinematic("minv_s_cm_ratio", false, false, false, true, "minv_s_cm_ratio"); // norm to unity

  // plot_hard_scatt_single_kinematic("pt_avg", false, false, false, false, "pt_avg", true);
  plot_hard_scatt_single_kinematic("pt_avg", false, false, true, false, "pt_avg", true);
  // plot_hard_scatt_single_kinematic("pt_avg", false, false, false, true, "pt_avg", true); // norm to unity

  // plot_hard_scatt_single_kinematic("ptlead_pair_pt", true, false, false, false, "pair_pt",true);
  plot_hard_scatt_single_kinematic("ptlead_pair_pt", true, false, true, false, "pair_pt",true);
  // plot_hard_scatt_single_kinematic("ptlead_pair_pt", true, false, false, true, "pair_pt",true); // norm to unity
  // plot_hard_scatt_single_kinematic("ptlead_pair_pt", false, true, false, false, "ptlead",true);
  plot_hard_scatt_single_kinematic("ptlead_pair_pt", false, true, true, false, "ptlead",true);
  // plot_hard_scatt_single_kinematic("ptlead_pair_pt", false, true, false, true, "ptlead",true); // norm to unity
  // plot_hard_scatt_single_kinematic("Deta_Dphi", false, true, false, false, "Deta");
  plot_hard_scatt_single_kinematic("Deta_Dphi", false, true, true, false, "Deta");
  // plot_hard_scatt_single_kinematic("Deta_Dphi", false, true, false, true, "Deta"); // norm to unity
}



