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

std::string subpl_titles[nMCmodes][nSigns][nDphi] = {
  {{"bb, same sign, #Dphi < #pi/2","bb, same sign, #Dphi #geq #pi/2"}, {"bb, opposite sign, #Dphi < #pi/2","bb, opposite sign, #Dphi #geq #pi/2"}}, 
  {{"cc, same sign, #Dphi < #pi/2","cc, same sign, #Dphi #geq #pi/2"}, {"cc, opposite sign, #Dphi < #pi/2","cc, opposite sign, #Dphi #geq #pi/2"}}
};

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


void plot_hard_scatt_single_kinematic(std::string kin, bool projx_2d, bool projy_2d, bool add, bool norm_unity, std::string kin1d, bool logx=false){

  TH1D* h[nMCmodes][nSigns][nDphi][nAncestorGroups];

  if (projx_2d && projy_2d){
    std::cout << "For single kinematic plotting, can only project either X or Y variable." << std::endl;
    throw std::exception();
  }

  if (add && norm_unity){
    std::cout << "Cannot both add and normalize to unity." << std::endl;
    throw std::exception();
  }

  for (int imc = 0; imc < nMCmodes; imc++){
    TCanvas* c = new TCanvas(("c"+mcmodes[imc]).c_str(),("c"+mcmodes[imc]).c_str(),2900,2000);
    c->Divide(2,2);
    f[imc] = TFile::Open((mc_path + fnames[imc]).c_str());

    for (unsigned int ksign = 0; ksign < nSigns; ksign++){
      for (int ldphi = 0; ldphi < nDphi; ldphi++){

        c->cd(ksign * nSigns + ldphi + 1);
        gPad->SetLeftMargin(0.16);
        // gPad->SetRightMargin(0.16);
        gPad->SetBottomMargin(0.135);
        gPad->SetLogx(logx);

        TLegend* l = new TLegend(0.55,0.74,0.87,0.87);
        l->SetBorderSize(0);
        l->SetFillStyle(0);
        l->SetTextFont(43);
        l->SetMargin(0.02);
        l->SetTextColor(1);

        for (int mgrp = 0; mgrp < nAncestorGroups; mgrp++){
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
          }

          if (add){
            hist_helper(h[imc][ksign][ldphi][mgrp], subpl_titles[imc][ksign][ldphi] + ", accumulative", norm_unity);
          }else if(norm_unity){
            hist_helper(h[imc][ksign][ldphi][mgrp], subpl_titles[imc][ksign][ldphi] + ", unity", norm_unity);
          }else{
            hist_helper(h[imc][ksign][ldphi][mgrp], subpl_titles[imc][ksign][ldphi], norm_unity);
          }
          h[imc][ksign][ldphi][mgrp]->SetLineColor(line_colors[mgrp]);
          h[imc][ksign][ldphi][mgrp]->SetMarkerColor(line_colors[mgrp]);
          l->AddEntry(h[imc][ksign][ldphi][mgrp], ancestor_grp_labels[mgrp].c_str(),"lp");
          // l->AddEntry(h[imc][ksign][ldphi][mgrp], "hello","lp");
        }

        float ymax;
        if (add){
          h[imc][ksign][ldphi][2]->Add(h[imc][ksign][ldphi][3]);
          h[imc][ksign][ldphi][1]->Add(h[imc][ksign][ldphi][2]);
          h[imc][ksign][ldphi][0]->Add(h[imc][ksign][ldphi][1]);
        }else{
          ymax = (h[imc][ksign][ldphi][0]->GetMaximum() > h[imc][ksign][ldphi][1]->GetMaximum())? h[imc][ksign][ldphi][0]->GetMaximum() : h[imc][ksign][ldphi][1]->GetMaximum();
          ymax = (ymax > h[imc][ksign][ldphi][2]->GetMaximum())? ymax : h[imc][ksign][ldphi][2]->GetMaximum();
          ymax = (ymax > h[imc][ksign][ldphi][3]->GetMaximum())? ymax : h[imc][ksign][ldphi][3]->GetMaximum();
          h[imc][ksign][ldphi][0]->GetYaxis()->SetRangeUser(0., ymax * 1.1);
        }

        h[imc][ksign][ldphi][0]->Draw("E");
        h[imc][ksign][ldphi][1]->Draw("E,same");
        h[imc][ksign][ldphi][2]->Draw("E,same");
        h[imc][ksign][ldphi][3]->Draw("E,same");
        l->Draw();
      }
    }

    if (add){
      c->SaveAs(Form("plots/mc_truth/hard_scatt/hard_scatt_%s%s_accumulative.png", kin1d.c_str(), mcmodes[imc].c_str()));
    }else if (norm_unity){
      c->SaveAs(Form("plots/mc_truth/hard_scatt/hard_scatt_%s%s_unity.png", kin1d.c_str(), mcmodes[imc].c_str()));
    }else{
      c->SaveAs(Form("plots/mc_truth/hard_scatt/hard_scatt_%s%s.png", kin1d.c_str(), mcmodes[imc].c_str()));
    }
    c->Close();
    delete c;


  }
}


void plot_hard_scatt(){
  plot_hard_scatt_single_kinematic("DR", false, false, false, false, "DR");
  // plot_hard_scatt_single_kinematic("DR", false, false, true, false, "DR"); // accumulative
  plot_hard_scatt_single_kinematic("DR", false, false, false, true, "DR"); // norm to unity
  plot_hard_scatt_single_kinematic("Dphi", false, false, false, false, "Dphi");
  // plot_hard_scatt_single_kinematic("Dphi", false, false, true, false, "Dphi"); // accumulative
  plot_hard_scatt_single_kinematic("Dphi", false, false, false, true, "Dphi"); // norm to unity
  plot_hard_scatt_single_kinematic("minv", false, false, false, false, "minv");
  // plot_hard_scatt_single_kinematic("minv", false, false, true, false, "minv"); // accumulative
  plot_hard_scatt_single_kinematic("minv", false, false, false, true, "minv"); // norm to unity

  plot_hard_scatt_single_kinematic("ptlead_pair_pt", true, false, false, false, "pair_pt",true);
  plot_hard_scatt_single_kinematic("ptlead_pair_pt", true, false, false, true, "pair_pt",true); // norm to unity
  plot_hard_scatt_single_kinematic("ptlead_pair_pt", false, true, false, false, "ptlead",true);
  plot_hard_scatt_single_kinematic("ptlead_pair_pt", false, true, false, true, "ptlead",true); // norm to unity
  plot_hard_scatt_single_kinematic("Deta_Dphi", false, true, false, false, "Deta");
  plot_hard_scatt_single_kinematic("Deta_Dphi", false, true, false, true, "Deta"); // norm to unity
}



