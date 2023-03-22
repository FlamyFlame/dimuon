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
const int nSigns = 2;

double pi = acos(-1.0);
float norm_factor[nDtTypes] = {1., 1., 1/256.8}; // normalizing to differential crossx; unit is pb
Color_t colors[nDtTypes] = {kRed, kRed, kBlack};

// std::string dRs[ndRcuts] = {"_dr1","_dr2","_dr3"};
// std::string dRTitles[ndRcuts] = {"#Delta R < 0.8","#Delta R < 1.2", "no #Delta R cut"};

std::string signs[nSigns] = {"_sign1", "_sign2"};
std::string signTitles[nSigns] = {"same sign", "opposite sign"};

TFile* f[nDtTypes];
std::string dt_paths[nDtTypes] = {"/usatlas/u/yuhanguo/usatlasdata/athena/runMCV2/","/usatlas/u/yuhanguo/usatlasdata/athena/runMCV2/","/usatlas/u/yuhanguo/usatlasdata/dimuon_data/"};
std::string fnames[nDtTypes] = {"histograms_mc_truth_bb_old.root","histograms_mc_truth_cc_old.root","histograms_real_pairs_pp_old.root"};
// std::string dtTitles[nDtTypes] = {"MC truth bb", "MC truth cc", "pp data"};
std::string dtTitles[nDtTypes] = {"MC truth", "", "pp data"};

TCanvas* c;
// TH1D* h[nDtTypes][nSigns][nGapCuts];
TH1D* h[nDtTypes][nSigns];


void initialize(){
  for (int jdt = 0; jdt < nDtTypes; jdt++){
    f[jdt] = TFile::Open((dt_paths[jdt] + fnames[jdt]).c_str());
  }
}

void hist_helper(TH1* h, float norm, bool norm_unity, std::string ytitle=""){

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

void plot_powheg_data_compr_noDPHI_div_single_kinematic(std::string kin, bool projx_2d, bool projy_2d, bool norm_unity, std::string kin1d, std::string kin_title, bool logx=false){

  initialize();

  if (projx_2d && projy_2d){
    std::cout << "For single kinematic plotting, can only project either X or Y variable." << std::endl;
    throw std::exception();
  }
  
  TCanvas* c = new TCanvas("c","c",2900,1000);
  c->Divide(nSigns,1);

  for (unsigned int ksign = 0; ksign < nSigns; ksign++){
    c->cd(ksign + 1);
    gPad->SetLeftMargin(0.2);
    gPad->SetBottomMargin(0.135);
    gPad->SetLogx(logx);


    // TLegend* l = new TLegend(0.24,0.6,0.5,0.89);
    TLegend* l = new TLegend(0.76,0.72,0.87,0.87);

    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetTextFont(42);
    // l->SetTextSize(32);
    l->SetMargin(0.02);
    l->SetTextColor(1);

    for (unsigned int jdt = 0; jdt < nDtTypes; jdt++){
      if (projx_2d){
        TH2D* h2d = (TH2D*) f[jdt]->Get(("h_" + kin + "_dr3" + signs[ksign] + "_gapcut1").c_str());
        if (norm_unity){
          h[jdt][ksign] = (TH1D*) h2d->ProjectionX(Form("hx_unity_%d%d",jdt+1,ksign+1));
        }else{
          h[jdt][ksign] = (TH1D*) h2d->ProjectionX(Form("hx_%d%d",jdt+1,ksign+1));
        }
      }else if (projy_2d){
        TH2D* h2d = (TH2D*) f[jdt]->Get(("h_" + kin + "_dr3" + signs[ksign] + "_gapcut1").c_str());
        if (norm_unity){
          h[jdt][ksign] = (TH1D*) h2d->ProjectionY(Form("hy_unity_%d%d",jdt+1,ksign+1));
        }else{
          h[jdt][ksign] = (TH1D*) h2d->ProjectionY(Form("hy_%d%d",jdt+1,ksign+1));
        }
      }else{ // 1D
        h[jdt][ksign] = (TH1D*) f[jdt]->Get(("h_" + kin + "_dr3" + signs[ksign] + "_gapcut1").c_str());
      }

      h[jdt][ksign]->SetMarkerColor(colors[jdt]);
      h[jdt][ksign]->SetLineColor(colors[jdt]);
    }

    h[0][ksign]->Add(h[1][ksign]);
    // h[0][ksign]->Rebin(2); // must first rebin then call hist_helper: the Scale(1,"width") depends on the # of bins

    hist_helper(h[0][ksign], norm_factor[0], norm_unity, Form("#frac{d#sigma}{d %s} [pb]", kin_title.c_str()));
    // if (norm_unity) cout << h[0][ksign]->Integral() << endl;
    l->AddEntry(h[0][ksign],dtTitles[0].c_str(),"lp");

    hist_helper(h[2][ksign], norm_factor[2], norm_unity, Form("#frac{d#sigma}{d %s} [pb]", kin_title.c_str()));
    l->AddEntry(h[2][ksign],dtTitles[2].c_str(),"lp");

    float ylim = 1.1 * ((h[0][ksign]->GetMaximum() > h[2][ksign]->GetMaximum())? h[0][ksign]->GetMaximum() : h[2][ksign]->GetMaximum());
    h[0][ksign]->GetYaxis()->SetRangeUser(0,ylim);
    h[0][ksign]->Draw("E");
    h[2][ksign]->Draw("E,same");

    l->AddEntry("",signTitles[ksign].c_str(),"");
    l->Draw();
  }

  if (norm_unity){
    c->SaveAs(Form("plots/mc_data_compr/noDPHI_div/%s_powheg_data_compr_unity_noDPHI_div.png", kin1d.c_str()));
  }else{
    c->SaveAs(Form("plots/mc_data_compr/noDPHI_div/%s_powheg_data_compr_noDPHI_div.png", kin1d.c_str()));
  }


  c->Close();
  delete c;

}


void plot_powheg_data_compr_noDPHI_div(){
  // kin, projx_2d, projy_2d, norm_unity, kin1d, logx=false
  // plot_powheg_data_compr_noDPHI_div_single_kinematic("DR", false, false, false, "DR", "#Delta R");
  // plot_powheg_data_compr_noDPHI_div_single_kinematic("DR", false, false, true, "DR", "#Delta R"); // norm to unity
  // plot_powheg_data_compr_noDPHI_div_single_kinematic("Dphi", false, false, false, "Dphi", "#Delta #phi");
  // plot_powheg_data_compr_noDPHI_div_single_kinematic("Dphi", false, false, true, "Dphi", "#Delta #phi"); // norm to unity

  // plot_powheg_data_compr_noDPHI_div_single_kinematic("ptlead_pair_pt", true, false, false, "pair_pt", "p_T^{pair}", true);
  // plot_powheg_data_compr_noDPHI_div_single_kinematic("ptlead_pair_pt", true, false, true, "pair_pt", "p_T^{pair}", true); // norm to unity
  // plot_powheg_data_compr_noDPHI_div_single_kinematic("ptlead_pair_pt", false, true, false, "ptlead", "p_T^{lead}", true);
  // plot_powheg_data_compr_noDPHI_div_single_kinematic("ptlead_pair_pt", false, true, true, "ptlead", "p_T^{lead}", true); // norm to unity
  plot_powheg_data_compr_noDPHI_div_single_kinematic("minv_pair_pt", false, true, false, "minv","m_{#mu#mu}",true);
  plot_powheg_data_compr_noDPHI_div_single_kinematic("minv_pair_pt", false, true, true, "minv","m_{#mu#mu}",true); // norm to unity
  // plot_powheg_data_compr_noDPHI_div_single_kinematic("Deta_Dphi", false, true, false, "Deta", "#Delta #eta");
  // plot_powheg_data_compr_noDPHI_div_single_kinematic("Deta_Dphi", false, true, true, "Deta", "#Delta #eta"); // norm to unity

}


