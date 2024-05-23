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


const int nDtTypes = 5;
const int nSigns = 2;
const int nDphis = 2;

double pi = acos(-1.0);
float norm_factor[nDtTypes] = {1., 1., 1., 1/256.8, 1/256.8}; // normalizing to differential crossx; unit is pb
Color_t colors[nDtTypes] = {kRed, kRed, kBlue, kBlack, kGreen+2};

std::string signs[nSigns] = {"_sign1", "_sign2"};
std::string signTitles[nSigns] = {"same sign", "opposite sign"};

std::string dphis[nDphis] = {"_near", "_away"};
std::string dphiTitles[nDphis] = {"#Delta #phi < #pi/2","#Delta #phi #geq #pi/2"};

TFile* f[nDtTypes];
std::string dt_paths[nDtTypes] = {"/usatlas/u/yuhanguo/usatlasdata/athena/runMCV2/","/usatlas/u/yuhanguo/usatlasdata/athena/runMCV2/","/usatlas/u/yuhanguo/usatlasdata/pythia/","/usatlas/u/yuhanguo/usatlasdata/dimuon_data/","/usatlas/u/yuhanguo/usatlasdata/dimuon_data/"};
std::string fnames[nDtTypes] = {"histograms_mc_truth_bb.root","histograms_mc_truth_cc.root","histograms_pythia.root","histograms_real_pairs_pp.root","histograms_real_pairs_pp.root"};
// std::string dtTitles[nDtTypes] = {"MC truth bb", "MC truth cc", "pp data"};
std::string dtTitles[nDtTypes] = {"POWHEG", "", "Pythia", "pp data - no gap cut", "pp data - with gap cut"};

TCanvas* c;
// TH1D* h[nDtTypes][nSigns][nGapCuts];
TH1D* h[nDtTypes][nSigns][nDphis];


void initialize(){
  for (int idt = 0; idt < nDtTypes; idt++){
    f[idt] = TFile::Open((dt_paths[idt] + fnames[idt]).c_str());
  }
}

void hist_helper(TH1* h, float norm, bool norm_unity, std::string title, std::string ytitle=""){

  h->SetStats(0);

  if (norm_unity){
    h->Scale(1.,"width");
    h->Scale(1./h->Integral("width"));
    h->GetYaxis()->SetTitle("pdf");
  }else{
    h->Scale(norm,"width");
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

void plot_mc_data_compr_w_gapcut_single_kinematic(std::string kin, bool projx_2d, bool projy_2d, bool norm_unity, std::string kin1d, std::string kin_title, bool ratio, bool logx=false){

  initialize();

  if (projx_2d && projy_2d){
    std::cout << "For single kinematic plotting, can only project either X or Y variable." << std::endl;
    throw std::exception();
  }

  if (norm_unity && ratio){
    std::cout << "Cannot both normalize to unity and take ratio." << std::endl;
    throw std::exception();
  }
  
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
      if (ksign == 0) l = new TLegend(0.73,0.64,0.89,0.89);
      else            l = new TLegend(0.7,0.64,0.91,0.89);

      l->SetBorderSize(0);
      l->SetFillStyle(0);
      l->SetTextFont(42);
      // l->SetTextSize(20);
      l->SetMargin(0.2);
      l->SetTextColor(1);

      for (unsigned int idt = 0; idt < nDtTypes; idt++){
        std::string hist_gapcut_postfix = "_gapcut1";
        if (idt == 2) hist_gapcut_postfix = "";
        else if (idt == 4) hist_gapcut_postfix = "_gapcut2";
        std::string hist_name = "h_" + kin + dphis[jdphi] + signs[ksign] + hist_gapcut_postfix;
        
        if (projx_2d){
          TH2D* h2d = (TH2D*) f[idt]->Get(hist_name.c_str());
          if (norm_unity){
            h[idt][ksign][jdphi] = (TH1D*) h2d->ProjectionX(Form("hx_unity_%d%d%d",idt+1,ksign+1,jdphi+1));
          }else{
            h[idt][ksign][jdphi] = (TH1D*) h2d->ProjectionX(Form("hx_%d%d%d",idt+1,ksign+1,jdphi+1));
          }
        }else if (projy_2d){
          TH2D* h2d = (TH2D*) f[idt]->Get(hist_name.c_str());
          if (norm_unity){
            h[idt][ksign][jdphi] = (TH1D*) h2d->ProjectionY(Form("hy_unity_%d%d%d",idt+1,ksign+1,jdphi+1));
          }else{
            h[idt][ksign][jdphi] = (TH1D*) h2d->ProjectionY(Form("hy_%d%d%d",idt+1,ksign+1,jdphi+1));
          }
        }else{ // 1D
          h[idt][ksign][jdphi] = (TH1D*) f[idt]->Get(hist_name.c_str());
        }

        h[idt][ksign][jdphi]->SetMarkerColor(colors[idt]);
        h[idt][ksign][jdphi]->SetLineColor(colors[idt]);
      }

      h[0][ksign][jdphi]->Add(h[1][ksign][jdphi]);
      // if (kin1d == "DR" || kin1d == "Dphi" || kin1d == "Deta") h[2][ksign][jdphi]->Rebin(2);
      // h[0][ksign][jdphi]->Rebin(2); // must first rebin then call hist_helper: the Scale(1,"width") depends on the # of bins

      if (ratio)
        hist_helper(h[0][ksign][jdphi], norm_factor[0], norm_unity, (signTitles[ksign] + ", " + dphiTitles[jdphi]).c_str(), "ratio to pp (no gap cut)");
      else
        hist_helper(h[0][ksign][jdphi], norm_factor[0], norm_unity, (signTitles[ksign] + ", " + dphiTitles[jdphi]).c_str(), Form("#frac{d#sigma}{d %s} [pb]", kin_title.c_str()));
      l->AddEntry(h[0][ksign][jdphi],dtTitles[0].c_str(),"lp");

      h[2][ksign][jdphi]->Scale(pow(10,6));
      for (unsigned int idt = 2; idt < nDtTypes; idt++){
        if (ratio)
          hist_helper(h[idt][ksign][jdphi], norm_factor[idt], norm_unity, (signTitles[ksign] + ", " + dphiTitles[jdphi]).c_str(), "ratio to pp (no gap cut)");
        else
          hist_helper(h[idt][ksign][jdphi], norm_factor[idt], norm_unity, (signTitles[ksign] + ", " + dphiTitles[jdphi]).c_str(), Form("#frac{d#sigma}{d %s} [pb]", kin_title.c_str()));
        // hist_helper(h[idt][ksign][jdphi], norm_factor[idt], norm_unity, (signTitles[ksign] + ", " + dphiTitles[jdphi]).c_str(), Form("#frac{d#sigma}{d %s} [pb]", kin_title.c_str()));
        l->AddEntry(h[idt][ksign][jdphi],dtTitles[idt].c_str(),"lp");
      }

      if (ratio){
        h[0][ksign][jdphi]->Divide(h[3][ksign][jdphi]);
        h[2][ksign][jdphi]->Divide(h[3][ksign][jdphi]);
        h[4][ksign][jdphi]->Divide(h[3][ksign][jdphi]);
        h[3][ksign][jdphi]->Divide(h[3][ksign][jdphi]);
      }

      float ylim = (h[0][ksign][jdphi]->GetMaximum() > h[2][ksign][jdphi]->GetMaximum())? h[0][ksign][jdphi]->GetMaximum() : h[2][ksign][jdphi]->GetMaximum();
      ylim = (ylim > h[3][ksign][jdphi]->GetMaximum())? ylim : h[3][ksign][jdphi]->GetMaximum();
      ylim = (ylim > h[4][ksign][jdphi]->GetMaximum())? ylim : h[4][ksign][jdphi]->GetMaximum();
      ylim *= 1.1;
      h[0][ksign][jdphi]->GetYaxis()->SetRangeUser(0,ylim);
      h[0][ksign][jdphi]->Draw("E");
      for (unsigned int idt = 2; idt < nDtTypes; idt++){
        h[idt][ksign][jdphi]->Draw("E,same");
      }
      // h[3][ksign][jdphi]->Draw("E,same");
      // h[4][ksign][jdphi]->Draw("E,same");

      // if (!norm_unity) std::cout << h[0][ksign][jdphi]->Integral("width") << " " << h[2][ksign][jdphi]->Integral("width") << " " << h[3][ksign][jdphi]->Integral("width") << std::endl;

      l->AddEntry("",(signTitles[ksign]).c_str(),"");
      l->AddEntry("",(dphiTitles[jdphi]).c_str(),"");
      l->Draw();
    }
  }
  if (norm_unity){
    c->SaveAs(Form("plots/mc_data_compr/%s_mc_data_compr_w_gapcut_unity.png", kin1d.c_str()));
  }else if (ratio){
    c->SaveAs(Form("plots/mc_data_compr/%s_mc_data_compr_w_gapcut_ratio.png", kin1d.c_str()));
  }else{
    c->SaveAs(Form("plots/mc_data_compr/%s_mc_data_compr_w_gapcut.png", kin1d.c_str()));
  }

  c->Close();
  delete c;

}

void plot_mc_data_compr_w_gapcut(){
  // // kin, projx_2d, projy_2d, norm_unity, kin1d, logx=false
  // plot_mc_data_compr_w_gapcut_single_kinematic("DR", false, false, false, "DR", "#Delta R", false);
  // plot_mc_data_compr_w_gapcut_single_kinematic("DR", false, false, true, "DR", "#Delta R", false); // norm to unity
  // plot_mc_data_compr_w_gapcut_single_kinematic("Dphi", false, false, false, "Dphi", "#Delta #phi", false);
  // plot_mc_data_compr_w_gapcut_single_kinematic("Dphi", false, false, true, "Dphi", "#Delta #phi", false); // norm to unity
  // plot_mc_data_compr_w_gapcut_single_kinematic("pt_asym", false, false, false, "pt_asym", "A", false);
  // plot_mc_data_compr_w_gapcut_single_kinematic("pt_asym", false, false, true, "pt_asym", "A", false); // norm to unity
  // plot_mc_data_compr_w_gapcut_single_kinematic("pair_pt_ptlead_ratio", false, false, false, "pair_pt_ptlead_ratio", "#frac{p_{T}^{pair}}{p_{T}^{lead}}", false);
  // plot_mc_data_compr_w_gapcut_single_kinematic("pair_pt_ptlead_ratio", false, false, true, "pair_pt_ptlead_ratio", "#frac{p_{T}^{pair}}{p_{T}^{lead}}", false); // norm to unity

  // plot_mc_data_compr_w_gapcut_single_kinematic("ptlead_pair_pt", true, false, false, "pair_pt", "p_T^{pair}", false, true);
  // plot_mc_data_compr_w_gapcut_single_kinematic("ptlead_pair_pt", true, false, true, "pair_pt", "p_T^{pair}", false, true); // norm to unity
  // plot_mc_data_compr_w_gapcut_single_kinematic("ptlead_pair_pt", false, true, false, "ptlead", "p_T^{lead}", false, true);
  // plot_mc_data_compr_w_gapcut_single_kinematic("ptlead_pair_pt", false, true, true, "ptlead", "p_T^{lead}", false, true); // norm to unity
  // plot_mc_data_compr_w_gapcut_single_kinematic("minv_pair_pt", false, true, false, "minv","m_{#mu#mu}", false, true);
  // plot_mc_data_compr_w_gapcut_single_kinematic("minv_pair_pt", false, true, true, "minv","m_{#mu#mu}", false, true); // norm to unity
  // plot_mc_data_compr_w_gapcut_single_kinematic("minv_pair_pt", false, true, false, "minv","m_{#mu#mu}", true, true); // ratio
  // plot_mc_data_compr_w_gapcut_single_kinematic("Deta_Dphi", false, true, false, "Deta", "#Delta #eta", false);
  // plot_mc_data_compr_w_gapcut_single_kinematic("Deta_Dphi", false, true, true, "Deta", "#Delta #eta", false); // norm to unity
  // plot_mc_data_compr_w_gapcut_single_kinematic("Deta_Dphi", false, true, false, "Deta", "#Delta #eta", false);
  // plot_mc_data_compr_w_gapcut_single_kinematic("Deta_Dphi", false, true, true, "Deta", "#Delta #eta", false); // norm to unity
  // plot_mc_data_compr_w_gapcut_single_kinematic("eta1_eta2", false, true, false, "eta1", "#eta_{lead}", false);
  plot_mc_data_compr_w_gapcut_single_kinematic("eta1_eta2", false, true, false, "eta1", "#eta_{lead}", false, true); // ratio
  // plot_mc_data_compr_w_gapcut_single_kinematic("eta1_eta2", false, true, true, "eta1", "#eta_{lead}", false); // norm to unity

}


