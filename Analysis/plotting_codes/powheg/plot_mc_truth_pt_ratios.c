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
const int nRatios = 3;

std::string mcmodes[nMCmodes] = {"bb","cc"};
std::string mc_path = "/usatlas/u/yuhanguo/usatlasdata/athena/runMCV2/";
std::string fnames[nMCmodes] = {"muon_pairs_mc_truth_bb.root", "muon_pairs_mc_truth_cc.root"};
std::vector<std::string> signs = {"_sign1", "_sign2"};
std::vector<std::string> dphis = {"_near", "_away"};
std::vector<std::string> hratios = {"h_pt_muon_pt_closest_hadr_ratio", "h_pt_closest_hadr_pt_furthest_hadr_ratio", "h_pt_hadr_hq_ratio"};
std::vector<std::string> hratio_titles = {"p_{T}^{#mu} / p_{T}^{closest hadron}", "p_{T}^{closest hadron} / p_{T}^{furthest hadron}", "p_{T}^{closest hadron} / p_{T}^{Q}"};
std::vector<Color_t> line_colors = {kBlack, kRed, kGreen + 2, kBlue};

std::string labels[nSigns][nDphi] ={{"same sign, near","same sign, away"}, {"opp sign, near","opp sign, away"}};

void hist_helper(TH1* h, std::string xtitle, std::string title, bool norm_unity){
  h->SetStats(0);
  h->Rebin(2);
  if (norm_unity){
    h->Scale(1./h->Integral());
    h->GetYaxis()->SetTitle("pdf");
  }else{
    h->Scale(1.,"width");
    h->GetYaxis()->SetTitle("d#sigma/dratio");
  }
  h->SetTitle(title.c_str());
  h->GetXaxis()->SetTitle(xtitle.c_str());
  h->GetYaxis()->SetLabelFont(43);
  h->GetYaxis()->SetLabelSize(32);
  h->GetYaxis()->SetLabelOffset(0.01);
  h->GetYaxis()->SetTitleFont(43);
  h->GetYaxis()->SetTitleSize(36);
  h->GetYaxis()->SetTitleOffset(1.95);
  h->GetXaxis()->SetLabelFont(43);
  h->GetXaxis()->SetLabelSize(32);
  h->GetXaxis()->SetLabelOffset(0.01);
  h->SetMarkerStyle(20);
  h->GetYaxis()->SetRangeUser(0.,h->GetMaximum() * 1.1);
  h->GetXaxis()->SetTitleFont(43);  
  h->GetXaxis()->SetTitleSize(36);
  h->GetXaxis()->SetTitleOffset(1.35);
}


void plot_mc_truth_pt_ratios_one_norm_mode(bool norm_unity){ // normalized to unity or absolute crossx

  TH1D* h[nMCmodes][nRatios][nSigns][nDphi];
  TFile* f[nMCmodes];

  TCanvas* c = new TCanvas("c","c",3000,1800);
  c->Divide(3,2);

  for (int imc = 0; imc < nMCmodes; imc++){
  // for (int imc = 0; imc < 1; imc++){
    f[imc] = TFile::Open((mc_path + fnames[imc]).c_str());

    for (int jratio = 0; jratio < nRatios; jratio++){

      c->cd(imc * nRatios + jratio + 1);
      gPad->SetLeftMargin(0.2);
      // gPad->SetRightMargin(0.16);
      gPad->SetBottomMargin(0.2);

      // TLegend* l = new TLegend(0.65,0.74,0.87,0.87);
      TLegend* l = new TLegend(0.27,0.74,0.62,0.88);
      l->SetBorderSize(0);
      l->SetFillStyle(0);
      l->SetTextFont(43);
      l->SetMargin(0.2);
      l->SetTextColor(1);

      for (unsigned int ksign = 0; ksign < nSigns; ksign++){
        for (int lphi = 0; lphi < nDphi; lphi++){
          h[imc][jratio][ksign][lphi] = (TH1D*) f[imc]->Get((hratios[jratio] + signs[ksign] + dphis[lphi]).c_str());
          hist_helper(h[imc][jratio][ksign][lphi],hratio_titles[jratio], mcmodes[imc] + ", " + hratio_titles[jratio], norm_unity);
          h[imc][jratio][ksign][lphi]->SetLineColor(line_colors[ksign * nSigns + lphi]);
          h[imc][jratio][ksign][lphi]->SetMarkerColor(line_colors[ksign * nSigns + lphi]);
          l->AddEntry(h[imc][jratio][ksign][lphi], labels[ksign][lphi].c_str(),"lp");
        }
      }
      l->AddEntry("",mcmodes[imc].c_str(),"");

      h[imc][jratio][1][1]->Draw("E");
      h[imc][jratio][1][0]->Draw("E,same");
      h[imc][jratio][0][1]->Draw("E,same");
      h[imc][jratio][0][0]->Draw("E,same");
      l->Draw();
    }
  }

  if (norm_unity) c->SaveAs("plots/mc_truth/hard_scatt/pt_ratios_unity.png");
  else c->SaveAs("plots/mc_truth/hard_scatt/pt_ratios.png");
  c->Close();
  delete c;
}


void plot_mc_truth_pt_ratios(){
  plot_mc_truth_pt_ratios_one_norm_mode(false);
  plot_mc_truth_pt_ratios_one_norm_mode(true);
}



