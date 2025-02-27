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
// const int nMCmodes = 2;
const int nDphi = 2;

// std::string mcmodes[nMCmodes] = {"bb","cc"};
std::string py_path = "/usatlas/u/yuhanguo/usatlasdata/pythia/";
// std::string fname = "muon_pairs_pythia.root";
std::string fname = "muon_pairs_pythia_0429.root";
std::vector<std::vector<std::string>> hist_names;

// std::string labels[nMCmodes][nSigns][nDphi] = {
//   {{"bb, same sign, #Delta #phi < #pi/2","bb, same sign, #Delta #phi #geq #pi/2"}, {"bb, opposite sign, #Delta #phi < #pi/2","bb, opposite sign, #Delta #phi #geq #pi/2"}}, 
//   {{"cc, same sign, #Delta #phi < #pi/2","cc, same sign, #Delta #phi #geq #pi/2"}, {"cc, opposite sign, #Delta #phi < #pi/2","cc, opposite sign, #Delta #phi #geq #pi/2"}}
// };
std::string labels[nSigns][nDphi] = {{"same sign, #Delta #phi < #pi/2","same sign, #Delta #phi #geq #pi/2"}, {"opposite sign, #Delta #phi < #pi/2","opposite sign, #Delta #phi #geq #pi/2"}};

void hist_helper(TH1* h, char * title){
  h->SetStats(0);
  h->SetTitle(title);
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

void pythia_muon_pair_parent_groups(){

  TH2D* h;
  TFile* f = TFile::Open((py_path + fname).c_str());

  hist_names = {{"h_parent_groups_sign1_near","h_parent_groups_sign1_away"},{"h_parent_groups_sign2_near","h_parent_groups_sign2_away"}};

  // for (int imc = 0; imc < nMCmodes; imc++){

  TCanvas* c = new TCanvas("c","c",2900,2000);
  c->Divide(2,2);

  for (unsigned int ksign = 0; ksign < nSigns; ksign++){
    for (int lphi = 0; lphi < nDphi; lphi++){

      c->cd(ksign * nSigns + lphi + 1);
      h = (TH2D*) f->Get(hist_names[ksign][lphi].c_str());
      h->Scale(1000);
      hist_helper(h, Form("%s (%.2f nb)",(labels[ksign][lphi]).c_str(), h->Integral()));

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

  c->SaveAs("/usatlas/u/yuhanguo/usatlasdata/pythia/plots/ancestor_tracing/pythia_muon_pair_parent_groups.png");
  c->Close();
  delete c;
  // }
}




