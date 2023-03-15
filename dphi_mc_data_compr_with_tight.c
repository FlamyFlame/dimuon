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


// const int ndRs = 3;
const int nGapCuts = 2;
const int nDtTypes = 4;
const int nSigns = 2;

double pi = acos(-1.0);
int nbins[nDtTypes] = {64,64,128,128};
// int nbins[nDtTypes] = {64,64,64,64};
// float plateau_list[nGapCuts][nDtTypes][nSigns] = {{{7.295, 11.441}, {0.033, 0.067}, {3587.4, 7614.5}},{{0.,0.},{0.,0.},{0.,0.}}};
float plateau_list[nDtTypes][nSigns] = {{7.295, 11.441}, {0.033, 0.067}, {3587.4, 7614.5},{3154.12, 5867.67}};
float norm_factor[nDtTypes] = {1., 1., 1/256.8, 1/256.8}; // normalizing to differential crossx; unit is pb
Color_t colors[nDtTypes] = {kRed, kRed, kBlack, kViolet};

// std::string dRs[ndRcuts] = {"_dr1","_dr2","_dr3"};
// std::string dRTitles[ndRcuts] = {"#Delta R < 0.8","#Delta R < 1.2", "no #Delta R cut"};

std::string signs[nSigns] = {"_sign1", "_sign2"};
std::string signTitles[nSigns] = {"same sign", "opposite sign"};

std::string gapcuts[nGapCuts] = {"_gapcut1", "_gapcut2"};
std::string gapcutTitles[nGapCuts] = {"no gap cut", "with gap cut"};

TFile* f[nDtTypes];
std::string dt_paths[nDtTypes] = {"/usatlas/u/yuhanguo/usatlasdata/athena/runMCV2/","/usatlas/u/yuhanguo/usatlasdata/athena/runMCV2/","/usatlas/u/yuhanguo/usatlasdata/dimuon_data/","/usatlas/u/yuhanguo/usatlasdata/dimuon_data/"};
std::string fnames[nDtTypes] = {"histograms_mc_truth_bb_old.root","histograms_mc_truth_cc_old.root","histograms_real_pairs_pp_old.root","histograms_real_pairs_pp_tight_old.root"};
// std::string dtTitles[nDtTypes] = {"MC truth bb", "MC truth cc", "pp data"};
std::string dtTitles[nDtTypes] = {"MC truth", "", "pp medium", "pp tight"};

TCanvas* c;
// TH1D* h[nDtTypes][nSigns][nGapCuts];
TH1D* h[nDtTypes][nSigns];

// std::string hist_names[nDtTypes][nSigns][nGapCuts];
std::string hist_names[nSigns];


void initialize(){
  for (int jdt = 0; jdt < nDtTypes; jdt++){
    f[jdt] = TFile::Open((dt_paths[jdt] + fnames[jdt]).c_str());
  }

  for (int ksign = 0; ksign < nSigns; ksign++){
    hist_names[ksign] = "h_Deta_Dphi_dr3" + signs[ksign] + "_gapcut1";
    // for (int lgapcut = 0; lgapcut < nGapCuts; lgapcut++){
      // hist_names[ksign][lgapcut] = "h_Deta_Dphi_dr3" + signs[ksign] + gapcuts[lgapcut];
    // }
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
  h->GetYaxis()->SetLabelOffset(0.036);
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
  // h->GetYaxis()->SetTitleOffset(h->GetYaxis()->GetTitleOffset()*0.8);
    // h->GetXaxis()->SetTitleOffset(h->GetXaxis()->GetTitleOffset()*0.8);
  // h->GetXaxis()->SetTitleSize(h->GetXaxis()->GetTitleSize()*1.8);
    // h->GetXaxis()->SetLabelSize(h->GetXaxis()->GetLabelSize()*1.35);
    // h->GetYaxis()->SetTitleSize(h->GetYaxis()->GetTitleSize()*1.8);
    // h->GetYaxis()->SetLabelSize(h->GetYaxis()->GetLabelSize()*1.35);
    // h->GetZaxis()->SetLabelSize(h->GetZaxis()->GetLabelSize()*1.35);
}

void dphi_mc_data_compr_with_tight_one_mode(bool substract_plateau, bool ratio, float scaling_factor = 1.){

  initialize();

  if (substract_plateau && ratio){
    std::cout << "Cannot both substract plateau and take ratio." << std::endl;
    throw std::exception();
  }

  
  TCanvas* c = new TCanvas("c","c",2900,1000);
  c->Divide(nSigns,1);

  for (unsigned int ksign = 0; ksign < nSigns; ksign++){
    c->cd(ksign + 1);
    gPad->SetLeftMargin(0.16);
    gPad->SetBottomMargin(0.135);

    TLegend* l = new TLegend(0.24,0.65,0.6,0.89);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetTextFont(42);
    l->SetTextSize(l->GetTextSize()*3);
    l->SetMargin(0.02);
    l->SetTextColor(1);

    for (unsigned int jdt = 0; jdt < nDtTypes; jdt++){
      TH2D* h2d = (TH2D*) f[jdt]->Get(hist_names[ksign].c_str());
      h[jdt][ksign] = (TH1D*) h2d->ProjectionX(Form("h%d%d",jdt+1,ksign+1));
      // if (!substract_plateau) std::cout << h[jdt][ksign]->Integral() << std::endl;

      if (substract_plateau){
        for (int i = 1; i <= nbins[jdt]; i++){ //important: bin number goes from 1 to 128, not 0 to 127
          float ni = plateau_list[jdt][ksign];
          // float ni = float(h[jdt][ksign]->GetBinContent(i));
          // if (ni > plateau_list[jdt][ksign]){
          //   ni = plateau_list[jdt][ksign];
          // }
          h[jdt][ksign]->AddBinContent(i,-ni);
        }
      }

      h[jdt][ksign]->SetMarkerColor(colors[jdt]);
      h[jdt][ksign]->SetLineColor(colors[jdt]);

    }

    h[0][ksign]->Add(h[1][ksign]);
    if (ratio){
      h[2][ksign]->Rebin(2);
      h[3][ksign]->Rebin(2);
    }else{
      h[0][ksign]->Rebin(2); // must first rebin then call hist_helper: the Scale(1,"width") depends on the # of bins
    }

    std::string y_title = (ratio)? "ratio" : "#frac{d#sigma}{d #Delta #phi} [pb]";
    hist_helper(h[0][ksign], norm_factor[0], y_title);
    l->AddEntry(h[0][ksign],dtTitles[0].c_str(),"lp");

    hist_helper(h[2][ksign], norm_factor[2], y_title);
    l->AddEntry(h[2][ksign],dtTitles[2].c_str(),"lp");

    hist_helper(h[3][ksign], norm_factor[3], y_title);
    l->AddEntry(h[3][ksign],dtTitles[3].c_str(),"lp");

    // std::cout << h[0][ksign]->Integral("width") << std::endl;

    // std::cout << h[0][ksign]->Integral("width") << ", " << h[2][ksign]->Integral("width") * 256.8 << std::endl;
        

    if (ratio){
      h[0][ksign]->Divide(h[2][ksign]);
      h[3][ksign]->Divide(h[2][ksign]);
      h[2][ksign]->Divide(h[2][ksign]);
    }else{
      float ylim = 1.1 * ((h[0][ksign]->GetMaximum() > h[2][ksign]->GetMaximum())? h[0][ksign]->GetMaximum() : h[2][ksign]->GetMaximum());
      h[2][ksign]->GetYaxis()->SetRangeUser(0,ylim);
    }

    h[2][ksign]->Draw("E");
    h[3][ksign]->Draw("E,same");
    h[0][ksign]->Draw("E,same");
    // h[ymax_ind][ksign]->Draw("E");
    // for (int jdt = 0; jdt < nDtTypes; jdt++){
    //   if (jdt != ymax_ind) h[jdt][ksign]->Draw("E,same");
    // }

    l->AddEntry("",signTitles[ksign].c_str(),"");
    if (substract_plateau)    l->AddEntry("","subtracting plateau","");
    l->Draw();
  }

  if (substract_plateau) c->SaveAs("plots/mc_data_compr/dphi_mc_pp_tight_compr_substract_plateau.png");
  else if (scaling_factor != 1.){
    if (ratio)  c->SaveAs(Form("plots/mc_data_compr/dphi_SCALED_%.1f_mc_pp_tight_compr_ratio.png",scaling_factor));
    else        c->SaveAs(Form("plots/mc_data_compr/dphi_SCALED_%.1f_mc_pp_tight_compr.png",scaling_factor));
  }
  else{
    if (ratio)  c->SaveAs("plots/mc_data_compr/dphi_mc_pp_tight_compr_ratio.png");
    else        c->SaveAs("plots/mc_data_compr/dphi_mc_pp_tight_compr.png");
  }                   
  c->Close();
  delete c;

}


void dphi_mc_data_compr_with_tight(){
  // dphi_mc_data_compr_with_tight_one_mode(true, false);
  dphi_mc_data_compr_with_tight_one_mode(false, true);
  dphi_mc_data_compr_with_tight_one_mode(false, false);
  // dphi_mc_data_compr_with_tight_one_mode(false, true, 5.);
  // dphi_mc_data_compr_with_tight_one_mode(false, false, 5.);
}


