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


const int nKinRanges = 5;
const int nSigns = 2;
const int nDphis = 2;

double pi = acos(-1.0);
Color_t colors[nKinRanges] = {kRed, kBlack, kBlue, kGreen+2, kViolet};

std::string signs[nSigns] = {"_sign1", "_sign2"};
std::string signTitles[nSigns] = {"same sign", "opposite sign"};

std::string dphis[nDphis] = {"_near", "_away"};
std::string dphiTitles[nDphis] = {"#Delta #phi < #pi/2","#Delta #phi #geq #pi/2"};

std::string kins[nKinRanges] = {"_k0","_k1", "_k2","_k3","_k4"};
std::string kinTitles[nKinRanges] = {"pTHat: 5-10 GeV", "pTHat: 10-25 GeV", "pTHat: 25-60 GeV", "pTHat: 60-120 GeV", "pTHat: 120-3200 GeV"};

TH1D* h[nKinRanges][nSigns][nDphis];

void hist_helper(TH1* h, bool norm_unity, bool accumulate, std::string title, std::string ytitle=""){

  h->SetStats(0);

  if (norm_unity){
    h->Scale(1.,"width");
    h->Scale(1./h->Integral("width"));
    h->GetYaxis()->SetTitle("pdf");
  }else{
    h->Scale(1.,"width");
    if (ytitle.length() != 0){
      h->GetYaxis()->SetTitle(ytitle.c_str());
    }
  }

  // h->SetTitle(title.c_str());
  // h->SetTitleSize(35);
  h->GetYaxis()->SetLabelFont(43);
  h->GetYaxis()->SetLabelSize(32);
  h->GetYaxis()->SetLabelOffset(0.01);
  h->GetYaxis()->SetTitleFont(43);
  h->GetYaxis()->SetTitleSize(32);
  h->GetYaxis()->SetTitleOffset(2.1);
  h->GetXaxis()->SetLabelFont(43);
  h->GetXaxis()->SetLabelSize(32);
  h->GetXaxis()->SetLabelOffset(0.01);
  h->GetXaxis()->SetTitleFont(43);  
  h->GetXaxis()->SetTitleSize(32);
  h->GetXaxis()->SetTitleOffset(1);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.9);
}

void plot_pythia_kn_single_kinematic(std::string kin, bool projx_2d, bool projy_2d, bool norm_unity, bool accum, std::string kin1d, std::string kin_title, bool logx=false){

  if (projx_2d && projy_2d){
    std::cout << "For single kinematic plotting, can only project either X or Y variable." << std::endl;
    throw std::exception();
  }

  if (norm_unity && accum){
    std::cout << "Cannot both normalize to unity and add the histograms." << std::endl;
    throw std::exception();
  }

  TFile* f = TFile::Open("/usatlas/u/yuhanguo/usatlasdata/pythia/histograms_pythia.root");
  
  TCanvas* c = new TCanvas("c","c",1800,1200);
  c->Divide(nSigns,2);

  for (unsigned int ksign = 0; ksign < nSigns; ksign++){
    for (unsigned int jdphi = 0; jdphi < 2; jdphi++){

      c->cd(ksign * 2 + jdphi + 1);
      gPad->SetLeftMargin(0.2);
      gPad->SetBottomMargin(0.135);
      // gPad->SetTopMargin(0.18);
      gPad->SetLogx(logx);

      // TLegend* l = new TLegend(0.24,0.6,0.5,0.89);
      TLegend* l;
      l = new TLegend(0.64,0.55,0.89,0.89);

      l->SetBorderSize(0);
      l->SetFillStyle(0);
      l->SetTextFont(42);
      // l->SetTextSize(20);
      l->SetMargin(0.02);
      l->SetTextColor(1);

      for (int ikin = 0; ikin < nKinRanges; ikin++){

        std::string hist_name = "h_" + kin + kins[ikin] + dphis[jdphi] + signs[ksign];
        
        if (projx_2d){
          TH2D* h2d = (TH2D*) f->Get(hist_name.c_str());
          if (norm_unity){
            h[ikin][ksign][jdphi] = (TH1D*) h2d->ProjectionX(Form("hx_unity_%d%d%d",ikin+1,ksign+1,jdphi+1));
          }else{
            h[ikin][ksign][jdphi] = (TH1D*) h2d->ProjectionX(Form("hx_%d%d%d",ikin+1,ksign+1,jdphi+1));
          }
        }else if (projy_2d){
          TH2D* h2d = (TH2D*) f->Get(hist_name.c_str());
          if (norm_unity){
            h[ikin][ksign][jdphi] = (TH1D*) h2d->ProjectionY(Form("hy_unity_%d%d%d",ikin+1,ksign+1,jdphi+1));
          }else{
            h[ikin][ksign][jdphi] = (TH1D*) h2d->ProjectionY(Form("hy_%d%d%d",ikin+1,ksign+1,jdphi+1));
          }
        }else{ // 1D
          h[ikin][ksign][jdphi] = (TH1D*) f->Get(hist_name.c_str());
        }

        h[ikin][ksign][jdphi]->SetMarkerColor(colors[ikin]);
        h[ikin][ksign][jdphi]->SetLineColor(colors[ikin]);

        // h[0][ksign][jdphi]->Add(h[1][ksign][jdphi]);

        h[ikin][ksign][jdphi]->Scale(pow(10,6));
        hist_helper(h[ikin][ksign][jdphi], norm_unity, accum, (signTitles[ksign] + ", " + dphiTitles[jdphi]).c_str(), Form("#frac{d#sigma}{d %s} [pb]", kin_title.c_str()));
        l->AddEntry(h[ikin][ksign][jdphi],kinTitles[ikin].c_str(),"lp");
      }

      float ylim = 0;
      if (accum){
        for (int ikin = 1; ikin < nKinRanges; ikin++){
          h[ikin][ksign][jdphi]->Add(h[ikin-1][ksign][jdphi]);
        }
        ylim = h[nKinRanges-1][ksign][jdphi]->GetMaximum();
      }else{
        for (int ikin = 0; ikin < nKinRanges; ikin++){
          ylim = (ylim > h[ikin][ksign][jdphi]->GetMaximum())? ylim : h[ikin][ksign][jdphi]->GetMaximum();
        }
      }

      ylim *= 1.1;
      h[0][ksign][jdphi]->GetYaxis()->SetRangeUser(0,ylim);
      h[0][ksign][jdphi]->Draw("E");
      if (ksign == 1 && jdphi == 1) std::cout << h[0][ksign][jdphi]->Integral("width") << std::endl;

      for (int ikin = 1; ikin < nKinRanges; ikin++){
        h[ikin][ksign][jdphi]->Draw("E,same");
        if (ksign == 1 && jdphi == 1) std::cout << h[ikin][ksign][jdphi]->Integral("width") << std::endl;
      }

      // if (!norm_unity) std::cout << h[0][ksign][jdphi]->Integral("width") << " " << h[2][ksign][jdphi]->Integral("width") << " " << h[3][ksign][jdphi]->Integral("width") << std::endl;

      l->AddEntry("",(signTitles[ksign]).c_str(),"");
      l->AddEntry("",(dphiTitles[jdphi]).c_str(),"");
      l->Draw();
    }
  }
  if (norm_unity){
    c->SaveAs(Form("plots/pythia/%s_pythia_kn_unity.png", kin1d.c_str()));
  }else if(accum){
    c->SaveAs(Form("plots/pythia/%s_pythia_kn_accum.png", kin1d.c_str()));
  }else{
    c->SaveAs(Form("plots/pythia/%s_pythia_kn.png", kin1d.c_str()));
  }

  c->Close();
  delete c;
  
}

void plot_pythia_kn(){
  // kin, projx_2d, projy_2d, norm_unity, kin1d, logx=false
  // plot_pythia_kn_single_kinematic("DR", false, false, false, false, "DR", "#Delta R");
  // plot_pythia_kn_single_kinematic("DR", false, false, true, false, "DR", "#Delta R"); // norm to unity
  // plot_pythia_kn_single_kinematic("Dphi", false, false, false, false, "Dphi", "#Delta #phi");
  // plot_pythia_kn_single_kinematic("Dphi", false, false, true, false, "Dphi", "#Delta #phi"); // norm to unity
  // plot_pythia_kn_single_kinematic("pt_asym", false, false, false, false, "pt_asym", "A");
  // plot_pythia_kn_single_kinematic("pt_asym", false, false, true, false, "pt_asym", "A"); // norm to unity
  // plot_pythia_kn_single_kinematic("pair_pt_ptlead_ratio", false, false, false, false, "pair_pt_ptlead_ratio", "#frac{p_{T}^{pair}}{p_{T}^{lead}}");
  // plot_pythia_kn_single_kinematic("pair_pt_ptlead_ratio", false, false, true, false, "pair_pt_ptlead_ratio", "#frac{p_{T}^{pair}}{p_{T}^{lead}}"); // norm to unity
  
  // plot_pythia_kn_single_kinematic("ptlead_pair_pt", true, false, false, false, "pair_pt", "p_T^{pair}", true);
  // plot_pythia_kn_single_kinematic("ptlead_pair_pt", true, false, true, false, "pair_pt", "p_T^{pair}", true); // norm to unity
  // plot_pythia_kn_single_kinematic("ptlead_pair_pt", true, false, false, true, "pair_pt", "p_T^{pair}", true); // norm to unity
  plot_pythia_kn_single_kinematic("ptlead_pair_pt", false, true, false, false, "ptlead", "p_T^{lead}", true);
  plot_pythia_kn_single_kinematic("ptlead_pair_pt", false, true, true, false, "ptlead", "p_T^{lead}", true); // norm to unity
  plot_pythia_kn_single_kinematic("ptlead_pair_pt", false, true, false, true, "ptlead", "p_T^{lead}", true); // accumulative
  // plot_pythia_kn_single_kinematic("minv_pair_pt", false, true, false, false, "minv","m_{#mu#mu}",true);
  // plot_pythia_kn_single_kinematic("minv_pair_pt", false, true, true, false, "minv","m_{#mu#mu}",true); // norm to unity
  // plot_pythia_kn_single_kinematic("Deta_Dphi", false, true, false, false, "Deta", "#Delta #eta");
  // plot_pythia_kn_single_kinematic("Deta_Dphi", false, true, true, false, "Deta", "#Delta #eta"); // norm to unity

}


