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
#include "string"
// #include "time.h"
// #include "struct_hist.h"


const int nSigns = 2;
const int nMCmodes = 2;

std::string mc_path = "/usatlas/u/yuhanguo/usatlasdata/athena/runMCV2/";
std::string fnames[nMCmodes] = {"muon_pairs_mc_truth_bb.root", "muon_pairs_mc_truth_cc.root"};
std::string hist_names[nSigns] = {"h_MuonPairParentGroups_sign1","h_MuonPairParentGroups_sign2"};
std::string labels[nMCmodes][nSigns] = {{"bb, same sign", "bb, opposite sign"}, {"cc, same sign", "cc, opposite sign"}};

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

void mc_muon_pair_parent_groups(){

  TH2D* h;
  TFile* f[nMCmodes];

  TCanvas* c = new TCanvas("c","c",2900,2000);
  c->Divide(2,2);

  for (int imc = 0; imc < nMCmodes; imc++){
    f[imc] = TFile::Open((mc_path+fnames[imc]).c_str());

    for (unsigned int ksign = 0; ksign < nSigns; ksign++){

      c->cd(imc * nMCmodes + ksign + 1);
      h = (TH2D*) f[imc]->Get(hist_names[ksign].c_str());
      hist_helper(h,labels[imc][ksign]);
      gPad->SetLeftMargin(0.16);
      gPad->SetRightMargin(0.16);
      gPad->SetBottomMargin(0.135);

      TLegend* l = new TLegend(0.65,0.74,0.87,0.87);
      l->SetBorderSize(0);
      l->SetFillStyle(0);
      l->SetTextFont(38);
      l->SetTextSize(l->GetTextSize()*3);
      l->SetMargin(0.02);
      l->SetTextColor(1);
      l->AddEntry(h, labels[imc][ksign].c_str(),"lp");
      h->Draw("colz");
      // l->Draw();
        
    }
  }
  c->SaveAs("plots/mc_truth/prt_grouping/muon_pair_parent_groups.png");
  c->Close();
  delete c;

}




