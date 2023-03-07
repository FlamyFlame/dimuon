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
// #include "time.h"
// #include "struct_hist.h"


const int nDtTypes = 3;

double pi = acos(-1.0);
float norm_factor[nDtTypes] = {1/256.8, 1/256.8, 1/256.8}; // normalizing to differential crossx; unit is pb
Color_t colors[nDtTypes] = {kBlack, kBlue, kRed};

TFile* f[nDtTypes];
std::string dt_path = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/";
std::string fnames[nDtTypes] = {"histograms_real_pairs_pp_no_resn_cut.root","histograms_real_pairs_pp_1s_only.root","histograms_real_pairs_pp.root"};
std::string dtTitles[nDtTypes] = {"pp - no resonance cut", "pp - Upsilon 1s cut only", "pp - Upsilon 1s, 2s, 3s cuts"};

TCanvas* c;
TH1D* h[nDtTypes];

void initialize(){
  for (int jdt = 0; jdt < nDtTypes; jdt++){
    f[jdt] = TFile::Open((dt_path + fnames[jdt]).c_str());
  }
}

void hist_helper(TH1* h, float norm, std::string ytitle=""){

  h->SetStats(0);
  h->Scale(norm,"width");
  if (ytitle.length() != 0){
    h->GetYaxis()->SetTitle(ytitle.c_str());
  }
  h->GetYaxis()->SetLabelFont(43);
  h->GetYaxis()->SetLabelSize(32);
  h->GetYaxis()->SetLabelOffset(0.01);
  h->GetYaxis()->SetTitleFont(43);
  h->GetYaxis()->SetTitleSize(32);
  h->GetYaxis()->SetTitleOffset(2.1);
  h->GetXaxis()->SetLabelFont(43);
  h->GetXaxis()->SetLabelSize(32);
  h->GetXaxis()->SetLabelOffset(0.02);
  h->GetXaxis()->SetTitleFont(43);  
  h->GetXaxis()->SetTitleSize(32);
  h->GetXaxis()->SetTitleOffset(1);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.9);
}

void dphi_resn_cuts_one_mode(bool ratio){

  initialize();
  
  TCanvas* c = new TCanvas("c","c",1200,1000);

  c->cd(1);
  gPad->SetLeftMargin(0.16);
  gPad->SetBottomMargin(0.135);

  TLegend* l = new TLegend(0.24,0.6,0.5,0.89);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->SetTextFont(42);
  l->SetTextSize(l->GetTextSize()*3);
  l->SetMargin(0.02);
  l->SetTextColor(1);

  for (unsigned int jdt = 0; jdt < nDtTypes; jdt++){
    TH2D* h2d = (TH2D*) f[jdt]->Get("h_Deta_Dphi_dr3_sign2_gapcut1");
    h[jdt] = (TH1D*) h2d->ProjectionX(Form("h%d",jdt+1));

    h[jdt]->SetMarkerColor(colors[jdt]);
    h[jdt]->SetLineColor(colors[jdt]);
    if (ratio)
      hist_helper(h[jdt], norm_factor[jdt], "ratio");
    else
      hist_helper(h[jdt], norm_factor[jdt], "#frac{d#sigma}{d #Delta #phi} [pb]");
    l->AddEntry(h[jdt],dtTitles[jdt].c_str(),"lp");
  }

  if (ratio){
    h[1]->Divide(h[0]);
    h[2]->Divide(h[0]);
    h[0]->Divide(h[0]);
  }

  float ylim = 1.1 * h[0]->GetMaximum();
  h[0]->GetYaxis()->SetRangeUser(0,ylim);

  h[0]->Draw("E");
  h[1]->Draw("E,same");
  h[2]->Draw("E,same");
  l->Draw();

  if (ratio)  c->SaveAs("plots/pp_data/resn_cuts/dphi_resn_cuts_ratio.png");
  else        c->SaveAs("plots/pp_data/resn_cuts/dphi_resn_cuts.png");
  c->Close();
  delete c;

}

void dphi_resn_cuts(){
  dphi_resn_cuts_one_mode(false);
  dphi_resn_cuts_one_mode(true);
}

