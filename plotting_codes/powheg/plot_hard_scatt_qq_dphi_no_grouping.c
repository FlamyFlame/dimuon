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
std::vector<std::string> signs = {"_sign1", "_sign2"};
std::vector<std::string> dphis = {"_near", "_away"};

std::string labels[nMCmodes][nSigns][nDphi] = {
  {{"bb, same sign, all","bb, same sign, #Dphi #geq #pi/2 only"}, {"bb, opposite sign, all","bb, opposite sign, #Dphi #geq #pi/2 only"}}, 
  {{"cc, same sign, all","cc, same sign, #Dphi #geq #pi/2 only"}, {"cc, opposite sign, all","cc, opposite sign, #Dphi #geq #pi/2 only"}}
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
  h->SetMarkerStyle(20);
  h->GetYaxis()->SetRangeUser(0.,h->GetMaximum() * 1.1);
  // h->GetXaxis()->SetTitleFont(43);  
  // h->GetXaxis()->SetTitleSize(32);
  // h->GetXaxis()->SetTitleOffset(1);
}


void plot_hard_scatt_qq_dphi_separate(){ // mode = 1: unweighted; mode = 2: weighted

  TH1D* h;
  TFile* f[nMCmodes];

  for (int imc = 0; imc < nMCmodes; imc++){

    f[imc] = TFile::Open((mc_path + fnames[imc]).c_str());
    TCanvas* c = new TCanvas(Form("c-%s",mcmodes[imc].c_str()),Form("c-%s",mcmodes[imc].c_str()),2900,2000);
    c->Divide(2,2);

    for (unsigned int ksign = 0; ksign < nSigns; ksign++){
      for (int lphi = 0; lphi < nDphi; lphi++){

        c->cd(ksign * nSigns + lphi + 1);
        h = (TH1D*) f[imc]->Get(("h_QQ_dphi" + signs[ksign] + dphis[lphi]).c_str());
        hist_helper(h,labels[imc][ksign][lphi]);
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
        h->Draw();
        // l->Draw(); 
      }
    }

    c->SaveAs(Form("plots/mc_truth/hard_scatt/hard_scatt_qq_dphi_%s.png",mcmodes[imc].c_str()));
    c->Close();
    delete c;
  }
}

void plot_hard_scatt_qq_dphi_no_grouping(){ // mode = 1: unweighted; mode = 2: weighted

  TH1D* h[nDphi];
  TFile* f[nMCmodes];

  TCanvas* c = new TCanvas("c","c",2900,2000);
  c->Divide(2,2);

  for (int imc = 0; imc < nMCmodes; imc++){
    f[imc] = TFile::Open((mc_path + fnames[imc]).c_str());

    for (unsigned int ksign = 0; ksign < nSigns; ksign++){

      c->cd(imc * nMCmodes + ksign + 1);
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

      for (int lphi = 0; lphi < nDphi; lphi++){
        h[lphi] = (TH1D*) f[imc]->Get(("h_QQ_dphi" + signs[ksign] + dphis[lphi]).c_str());
        hist_helper(h[lphi],labels[imc][ksign][lphi]);
        l->AddEntry(h[lphi], labels[imc][ksign][lphi].c_str(),"lp");
      }

      h[0]->Add(h[1]);
      h[0]->SetLineColor(kRed);
      h[0]->SetMarkerColor(kRed);
      h[1]->SetLineColor(kBlack);
      h[1]->SetMarkerColor(kBlack);
      h[0]->Draw("E");
      h[1]->Draw("E,same");
    }
  }

  c->SaveAs("plots/mc_truth/hard_scatt/hard_scatt_qq_dphi.png");
  c->Close();
  delete c;


}



