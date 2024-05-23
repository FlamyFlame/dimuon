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
  {{"bb, same sign, #Dphi < #pi/2","bb, same sign, #Dphi #geq #pi/2"}, {"bb, opposite sign, #Dphi < #pi/2","bb, opposite sign, #Dphi #geq #pi/2"}}, 
  {{"cc, same sign, #Dphi < #pi/2","cc, same sign, #Dphi #geq #pi/2"}, {"cc, opposite sign, #Dphi < #pi/2","cc, opposite sign, #Dphi #geq #pi/2"}}
};

void hist_helper(TH1* h, std::string title){
  h->SetStats(0);
  h->SetTitle(title.c_str());
  h->GetYaxis()->SetLabelFont(43);
  h->GetYaxis()->SetLabelSize(38);
  h->GetYaxis()->SetLabelOffset(0.01);
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetLogz();
  // h->GetYaxis()->SetTitleFont(43);
  // h->GetYaxis()->SetTitleSize(32);
  // h->GetYaxis()->SetTitleOffset(2);
  h->GetXaxis()->SetLabelFont(43);
  h->GetXaxis()->SetLabelSize(38);
  h->GetXaxis()->SetLabelOffset(0.01);
  h->SetMarkerStyle(20);
  // h->GetXaxis()->SetTitleFont(43);  
  // h->GetXaxis()->SetTitleSize(32);
  // h->GetXaxis()->SetTitleOffset(1);
}


void plot_one_mode(int mode){ // mode = 1: 2D; mode = 2: pt_lead 1D; mode = 3: pair pT 1D

  TH2D* h;
  TH1D* hx[nMCmodes][nSigns][nDphi];
  TH1D* hy[nMCmodes][nSigns][nDphi];

  TFile* f[nMCmodes];

  if (mode == 1){
    for (int imc = 0; imc < nMCmodes; imc++){

      f[imc] = TFile::Open((mc_path + fnames[imc]).c_str());
      TCanvas* c = new TCanvas(Form("c-%s",mcmodes[imc].c_str()),Form("c-%s",mcmodes[imc].c_str()),2900,2000);
      c->Divide(2,2);

      for (unsigned int ksign = 0; ksign < nSigns; ksign++){
        for (int lphi = 0; lphi < nDphi; lphi++){

          c->cd(ksign * nSigns + lphi + 1);
          h = (TH2D*) f[imc]->Get(("h_ptlead_pair_pt" + signs[ksign] + dphis[lphi]).c_str());
          // hx = (TH1D*) h->ProjectionX("hx");
          // hy = (TH1D*) h->ProjectionY("hy");

          hist_helper(h,labels[imc][ksign][lphi]);
          gPad->SetLeftMargin(0.16);
          gPad->SetRightMargin(0.16);
          gPad->SetBottomMargin(0.135);

          h->Draw();
        }
      }

      c->SaveAs(Form("plots/mc_truth/ptlead_pair_pt_%s.png",mcmodes[imc].c_str()));
      c->Close();
      delete c;
    }
  }else{
    TCanvas* c = new TCanvas("c","c",2900,2000);

    c->Divide(2,2);

    for (unsigned int ksign = 0; ksign < nSigns; ksign++){
      for (int lphi = 0; lphi < nDphi; lphi++){
        c->cd(ksign * nSigns + lphi + 1);

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

        for (int imc = 0; imc < nMCmodes; imc++){

          f[imc] = TFile::Open((mc_path + fnames[imc]).c_str());
          h = (TH2D*) f[imc]->Get(("h_ptlead_pair_pt" + signs[ksign] + dphis[lphi]).c_str());
          hx[imc][ksign][lphi] = (TH1D*) h->ProjectionX(Form("hx%u%u",ksign+1,lphi+1));
          hy[imc][ksign][lphi] = (TH1D*) h->ProjectionY(Form("hy%u%u",ksign+1,lphi+1));

          if (mode == 2){
            hist_helper(hx[imc][ksign][lphi],labels[imc][ksign][lphi]);
            l->AddEntry(hx[imc][ksign][lphi], labels[imc][ksign][lphi].c_str(),"lp");
          }else{
            hist_helper(hy[imc][ksign][lphi],labels[imc][ksign][lphi]);
            l->AddEntry(hy[imc][ksign][lphi], labels[imc][ksign][lphi].c_str(),"lp");
          }
        }

        hx[0][ksign][lphi]->Draw();
        hx[1][ksign][lphi]->Draw("same");
        l->Draw(); 
      }
    }

    if (mode == 2) c->SaveAs("plots/mc_truth/pair_pt.png");
    else c->SaveAs("plots/mc_truth/ptlead.png");
    c->Close();
    delete c;
  }
}


void plot_pt_lead(){
  plot_one_mode(1);
  plot_one_mode(2);
  plot_one_mode(3);
}






