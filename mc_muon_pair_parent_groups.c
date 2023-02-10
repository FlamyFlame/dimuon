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


const int nSigns = 2;
const int nMCmodes = 2;
const int nDphi = 2;

std::string mcmodes[nMCmodes] = {"bb","cc"};
std::string mc_path = "/usatlas/u/yuhanguo/usatlasdata/athena/runMCV2/";
std::string fnames[nMCmodes] = {"muon_pairs_mc_truth_bb.root", "muon_pairs_mc_truth_cc.root"};
std::vector<std::vector<std::string>> hist_names;

std::string labels[nMCmodes][nSigns][nDphi] = {
  {{"bb, same sign, #Dphi < #pi/2","bb, same sign, #Dphi #geq #pi/2"}, {"bb, opposite sign, #Dphi < #pi/2","bb, opposite sign, #Dphi #geq #pi/2"}}, 
  {{"cc, same sign, #Dphi < #pi/2","cc, same sign, #Dphi #geq #pi/2"}, {"cc, opposite sign, #Dphi < #pi/2","cc, opposite sign, #Dphi #geq #pi/2"}}
};

void hist_helper(TH1* h, std::string title){
  h->SetStats(0);
  h->SetTitle(title.c_str());
  h->GetYaxis()->SetLabelFont(43);
  h->GetYaxis()->SetLabelSize(38);
  h->GetYaxis()->SetLabelOffset(0.01);
  // h->GetYaxis()->SetTitleFont(43);
  // h->GetYaxis()->SetTitleSize(32);
  // h->GetYaxis()->SetTitleOffset(2);
  h->GetXaxis()->SetLabelFont(43);
  h->GetXaxis()->SetLabelSize(38);
  h->GetXaxis()->SetLabelOffset(0.01);
  // h->GetXaxis()->SetTitleFont(43);  
  // h->GetXaxis()->SetTitleSize(32);
  // h->GetXaxis()->SetTitleOffset(1);
}

void plot_one_mode(int mode){

  TH2D* h;
  TFile* f[nMCmodes];

  if (mode == 1) hist_names = {{"h_parent_groups_sign1_near","h_parent_groups_sign1_away"},{"h_parent_groups_sign2_near","h_parent_groups_sign2_away"}};
  else           hist_names = {{"h_parent_groups_weighted_sign1_near","h_parent_groups_weighted_sign1_away"},{"h_parent_groups_weighted_sign2_near","h_parent_groups_weighted_sign2_away"}};

  for (int imc = 0; imc < nMCmodes; imc++){

    f[imc] = TFile::Open((mc_path + fnames[imc]).c_str());
    TCanvas* c = new TCanvas(Form("c-%s",mcmodes[imc].c_str()),Form("c-%s",mcmodes[imc].c_str()),2900,2000);
    c->Divide(2,2);

    for (unsigned int ksign = 0; ksign < nSigns; ksign++){
      for (int lphi = 0; lphi < nDphi; lphi++){

        c->cd(ksign * nSigns + lphi + 1);
        h = (TH2D*) f[imc]->Get(hist_names[ksign][lphi].c_str());
        hist_helper(h,labels[imc][ksign][lphi] + " (" + std::to_string(static_cast<int>(h->GetEntries())) + ")");
        gPad->SetLeftMargin(0.16);
        gPad->SetRightMargin(0.16);
        gPad->SetBottomMargin(0.135);

        // TLegend* l = new TLegend(0.65,0.74,0.87,0.87);
        // l->SetBorderSize(0);
        // l->SetFillStyle(0);
        // l->SetTextFont(38);
        // l->SetTextSize(l->GetTextSize()*3);
        // l->SetMargin(0.02);
        // l->SetTextColor(1);
        // l->AddEntry(h, labels[imc][ksign][lphi].c_str(),"lp");
        h->Draw("colz");
        // l->Draw(); 
      }
    }

    if (mode == 1) c->SaveAs(Form("plots/mc_truth/prt_grouping/muon_pair_parent_groups_%s.png",mcmodes[imc].c_str()));
    else           c->SaveAs(Form("plots/mc_truth/prt_grouping/muon_pair_parent_groups_%s_weighted.png",mcmodes[imc].c_str()));
    c->Close();
    delete c;
  }
}

void mc_muon_pair_parent_groups(){
  plot_one_mode(1);
  plot_one_mode(2);
}




