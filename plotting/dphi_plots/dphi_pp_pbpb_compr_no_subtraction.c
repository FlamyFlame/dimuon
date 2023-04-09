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
const int nCtrBins = 3;
const int nSigns = 2;

double pi = acos(-1.0);
int nbins[nCtrBins] = {128,128,128};
int plateau_list[nGapCuts][nCtrBins][nSigns] = {{{4774,4776},{29,47},{3659,7680}},{{4458,4470},{27,44},{3446,7264}}};
float norm_factor[nCtrBins] = {1./18.8, 1./0.42, 1/256.8};
  // normalizing to differential yield / TAA / N_coll for PbPb (in each centrality bin)
  // normalizing to differential crossx for pp
  // unit is pb

// std::string dRs[ndRcuts] = {"_dr1","_dr2","_dr3"};
// std::string dRTitles[ndRcuts] = {"#Delta R < 0.8","#Delta R < 1.2", "no #Delta R cut"};

std::string ctrs[nCtrBins] = {"_ctr1","_ctr4",""};
std::string ctr_dirs[nCtrBins] = {"ctr-binned/","ctr-binned/",""};
std::string ctrTitles[nCtrBins] = {"Centrality 0-20", "Centrality 60-80", "pp"};

std::string signs[nSigns] = {"_sign1", "_sign2"};
std::string signTitles[nSigns] = {"same sign", "opposite sign"};

std::string gapcuts[nSigns] = {"_gapcut1", "_gapcut2"};
std::string gapcutTitles[nSigns] = {"no gap cut", "with gap cut"};

TFile* f[nCtrBins];
std::string dt_path = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/";
std::string fnames[nCtrBins] = {"histograms_real_pairs.root","histograms_real_pairs_no_rebinning.root","histograms_real_pairs_pp.root"};


TCanvas* c;
TH1D* h[nCtrBins][nSigns][nGapCuts];

std::string hist_names[nCtrBins][nSigns][nGapCuts];


void initialize(){
  for (int jctr = 0; jctr < nCtrBins; jctr++){
    f[jctr] = TFile::Open((dt_path+fnames[jctr]).c_str());
    for (int ksign = 0; ksign < nSigns; ksign++){
      for (int lgapcut = 0; lgapcut < nGapCuts; lgapcut++){
        hist_names[jctr][ksign][lgapcut] = ctr_dirs[jctr] + "h_eta_avg_Dphi_dr3" + ctrs[jctr] + signs[ksign] + gapcuts[lgapcut];
      }
    }
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
  h->GetYaxis()->SetTitleOffset(1.8);
  h->GetXaxis()->SetLabelFont(43);
  h->GetXaxis()->SetLabelSize(32);
  h->GetXaxis()->SetLabelOffset(0.02);
  h->GetXaxis()->SetTitleFont(43);  
  h->GetXaxis()->SetTitleSize(32);
  h->GetXaxis()->SetTitleOffset(1);
  // h->GetYaxis()->SetTitleOffset(h->GetYaxis()->GetTitleOffset()*0.8);
    // h->GetXaxis()->SetTitleOffset(h->GetXaxis()->GetTitleOffset()*0.8);
  // h->GetXaxis()->SetTitleSize(h->GetXaxis()->GetTitleSize()*1.8);
    // h->GetXaxis()->SetLabelSize(h->GetXaxis()->GetLabelSize()*1.35);
    // h->GetYaxis()->SetTitleSize(h->GetYaxis()->GetTitleSize()*1.8);
    // h->GetYaxis()->SetLabelSize(h->GetYaxis()->GetLabelSize()*1.35);
    // h->GetZaxis()->SetLabelSize(h->GetZaxis()->GetLabelSize()*1.35);
}

void dphi_pp_pbpb_compr_no_subtraction(){

  initialize();

  for (int lgap = 0; lgap < nGapCuts; lgap++){
    TCanvas* c = new TCanvas(Form("c%d",lgap),Form("c%d",lgap),2900,1000);
    c->Divide(nSigns,1);

    for (unsigned int ksign = 0; ksign < nSigns; ksign++){
      c->cd(ksign + 1);
      gPad->SetLeftMargin(0.16);
      gPad->SetBottomMargin(0.135);

      TLegend* l = new TLegend(0.2,0.7,0.5,0.87);
      l->SetBorderSize(0);
      l->SetFillStyle(0);
      l->SetTextFont(42);
      l->SetTextSize(l->GetTextSize()*3);
      l->SetMargin(0.02);
      l->SetTextColor(1);

      for (unsigned int jctr = 0; jctr < nCtrBins; jctr++){
        TH2D* h2d = (TH2D*) f[jctr]->Get(hist_names[jctr][ksign][lgap].c_str());
        h[jctr][ksign][lgap] = (TH1D*) h2d->ProjectionX();

        hist_helper(h[jctr][ksign][lgap], norm_factor[jctr], "#frac{1}{T_{AA}} #frac{1}{N_{coll}} #frac{dN}{d #Delta #phi} [pb]");
        l->AddEntry(h[jctr][ksign][lgap],ctrTitles[jctr].c_str(),"lp");
        h[jctr][ksign][lgap]->SetMarkerStyle(20);
        h[jctr][ksign][lgap]->SetMarkerSize(0.6);
      }
                
      h[0][ksign][lgap]->SetMarkerColor(kRed);
      h[0][ksign][lgap]->SetLineColor(kRed);
      h[1][ksign][lgap]->SetMarkerColor(kBlue);
      h[1][ksign][lgap]->SetLineColor(kBlue);
      // first plot whichever with largest normalized peak at Delta phi ~ pi
      // since this sets the maximum y value

      float ylim_arr[nCtrBins];
      for (int ictr = 0; ictr < nCtrBins; ictr++){
        ylim_arr[ictr] = 1.1 * std::max(h[ictr][ksign][lgap]->GetBinContent(128), h[ictr][ksign][lgap]->GetBinContent(1));
        // std::cout << "ctr-group "<< ictr << ", max y value: " << ylim_arr[ictr] << std::endl;
        std::cout << "ctr-group "<< ictr << ", total integral (normalized): " << h[ictr][ksign][lgap]->Integral("width") << std::endl;
      }
      int ymax_ind = std::max_element(ylim_arr, ylim_arr + nCtrBins) - ylim_arr;
      float ylim = *std::max_element(ylim_arr, ylim_arr + nCtrBins);

      h[ymax_ind][ksign][lgap]->GetYaxis()->SetRangeUser(0,ylim);
      h[ymax_ind][ksign][lgap]->Draw("E");
      for (int ictr = 0; ictr < nCtrBins; ictr++){
        if (ictr != ymax_ind) {
          h[ictr][ksign][lgap]->Draw("E,same");
          // std::cout << h[ictr][ksign][lgap]->GetEntries() << std::endl;
        }
      }

      l->AddEntry("",(gapcutTitles[lgap]+", " + signTitles[ksign]).c_str(),"");
      l->Draw();
    }

    c->SaveAs(Form("plots/pbpb_pp_compr/dphi_pp_pbpb_compr_gapcut%d_no_subtraction.png",lgap+1));
    c->Close();
    delete c;
  }

}


