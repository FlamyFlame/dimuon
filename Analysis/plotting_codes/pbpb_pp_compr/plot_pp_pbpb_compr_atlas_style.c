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
#include "/usatlas/u/yuhanguo/plotstyle/AtlasStyle.C"
#include "/usatlas/u/yuhanguo/plotstyle/AtlasUtils.C"
// #include "string"
// #include "time.h"
// #include "struct_hist.h"


// const int ndRs = 3;
const int nGapCuts = 2;
const int nCtrBins = 7;
const int nLines = 4;
const int nSigns = 2;

double pi = acos(-1.0);

float norm_factor[nCtrBins] = {0.05128, 0.06536, 0.04630, 0.07602, 0.08503, 0.30441, 1/256.8};
  // normalizing to differential yield / TAA / N_coll for PbPb (in each centrality bin)
  // normalizing to differential crossx for pp
  // unit is pb

// std::string dRs[ndRcuts] = {"_dr1","_dr2","_dr3"};
// std::string dRTitles[ndRcuts] = {"#Delta R < 0.8","#Delta R < 1.2", "no #Delta R cut"};

std::string ctrs[nCtrBins] = {"_ctr1","_ctr2","_ctr3","_ctr4","_ctr5","_ctr6",""};
std::string ctr_dirs[nCtrBins] = {"ctr-binned/","ctr-binned/","ctr-binned/","ctr-binned/","ctr-binned/","ctr-binned/",""};
std::string ctrTitles[nCtrBins] = {"Centrality 0-5", "Centrality 5-10", "Centrality 10-20", "Centrality 20-30", "Centrality 30-50", "Centrality 50-80", "pp"};

std::string signs[nSigns] = {"_sign1", "_sign2"};
std::string signTitles[nSigns] = {"same sign", "opposite sign"};

std::string gapcuts[nSigns] = {"_gapcut1", "_gapcut2"};
std::string gapcutTitles[nSigns] = {"no gap cut", "with gap cut"};

TFile* f[nCtrBins];
std::string dt_path = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/";
std::string fnames[nCtrBins] = {"histograms_real_pairs_5_10_20_30_50_80.root",
"histograms_real_pairs_5_10_20_30_50_80.root",
"histograms_real_pairs_5_10_20_30_50_80.root",
"histograms_real_pairs_5_10_20_30_50_80.root",
"histograms_real_pairs_5_10_20_30_50_80.root",
"histograms_real_pairs_5_10_20_30_50_80.root",
"histograms_real_pairs_pp.root"};


TCanvas* c;
TH1D* h[nLines][nSigns][nGapCuts];

std::string hist_names[nCtrBins][nSigns][nGapCuts];


void initialize(std::string kin){
  for (int jctr = 0; jctr < nCtrBins; jctr++){
    f[jctr] = TFile::Open((dt_path+fnames[jctr]).c_str());
    // for (int ksign = 0; ksign < nSigns; ksign++){
    //   for (int lgapcut = 0; lgapcut < nGapCuts; lgapcut++){
    //     hist_names[jctr][ksign][lgapcut] = ctr_dirs[jctr] + "h_" + kin + "_dr3" + ctrs[jctr] + signs[ksign] + gapcuts[lgapcut];
    //   }
    // }
  }
}

void hist_helper(TH1* h, float norm, std::string ytitle=""){

  h->SetStats(0);
  h->Scale(norm,"width");
  if (ytitle.length() != 0){
    h->GetYaxis()->SetTitle(ytitle.c_str());
  }
  // h->SetMarkerStyle(20);
  // h->SetMarkerSize(4);
  // h->SetLineWidth(2);
  // h->GetYaxis()->SetLabelFont(43);
  // h->GetYaxis()->SetLabelSize(32);
  // h->GetYaxis()->SetLabelOffset(0.01);
  // h->GetYaxis()->SetTitleFont(43);
  // h->GetYaxis()->SetTitleSize(32);
  // h->GetYaxis()->SetTitleOffset(2);
  // h->GetXaxis()->SetLabelFont(43);
  // h->GetXaxis()->SetLabelSize(32);
  // h->GetXaxis()->SetLabelOffset(0.02);
  // h->GetXaxis()->SetTitleFont(43);  
  // h->GetXaxis()->SetTitleSize(32);
  // h->GetXaxis()->SetTitleOffset(1);
}

void plot_pp_pbpb_compr_single_kin(std::string kin, bool projx_2d, bool projy_2d, std::string kin1d, std::string kin_title, bool logx=false){

  SetAtlasStyle();
  initialize(kin);
  TH2D* h2d;

  if (projx_2d && projy_2d){
    std::cout << "For single kinematic plotting, can only project either X or Y variable." << std::endl;
    throw std::exception();
  }


  for (int lgap = 0; lgap < nGapCuts; lgap++){
    TCanvas* c = new TCanvas(Form("c%d",lgap),Form("c%d",lgap),2900,2000);
    c->Divide(nSigns,2);

    for (unsigned int ksign = 0; ksign < nSigns; ksign++){
      ////////////////////////////////////////// pp //////////////////////////////////////////
      std::string hist_name_pp = ctr_dirs[nCtrBins - 1] + "h_" + kin + "_dr3" + ctrs[nCtrBins-1] + signs[ksign] + gapcuts[lgap];
      // h2d = (TH2D*) f[nCtrBins-1]->Get(hist_names[nCtrBins - 1][ksign][lgap].c_str());
      h2d = (TH2D*) f[nCtrBins-1]->Get(hist_name_pp.c_str());
      // cout << ctrs[nCtrBins-1] << hist_name_pp << endl;
      // h[nLines-1][ksign][lgap] = (TH1D*) h2d->ProjectionY();
      if (projx_2d){
        h2d = (TH2D*) f[nCtrBins-1]->Get(hist_name_pp.c_str());
        h[nLines-1][ksign][lgap] = (TH1D*) h2d->ProjectionX(("h_" + kin1d + "_pp" + signs[ksign] + gapcuts[lgap]).c_str());
      }
      else if (projy_2d){
        h2d = (TH2D*) f[nCtrBins-1]->Get(hist_name_pp.c_str());
        h[nLines-1][ksign][lgap] = (TH1D*) h2d->ProjectionY(("h_" + kin1d + "_pp" + signs[ksign] + gapcuts[lgap]).c_str());
      }else{ // 1D
        h[nLines-1][ksign][lgap] = (TH1D*) f[nCtrBins-1]->Get(hist_name_pp.c_str());
        // std::cout << kin1d << ", " << subpl_titles[imc][ksign] << ", " << origin_catg_labels[mgrp] << ", dphi " << ldphi << ", integral: " << h[imc][ksign][mgrp]->Integral() << std::endl;
      }

      h[nLines-1][ksign][lgap]->Scale(1./1000.);
      // hist_helper(h[nLines-1][ksign][lgap], norm_factor[nCtrBins-1], "#frac{1}{T_{AA}} #frac{1}{N_{coll}} #frac{dN}{d " + kin_title + "} [nb]");


      for (unsigned int subpl = 0; subpl < 2; subpl++){
        c->cd(ksign * 2 + subpl + 1);
        // gPad->SetLeftMargin(0.16);
        // gPad->SetBottomMargin(0.135);

        // TLegend* l = new TLegend(0.18,0.64,0.53,0.87);
        TLegend* l = new TLegend(0.6,0.67,0.9,0.9);
        l->SetBorderSize(0);
        l->SetFillStyle(0);
        l->SetTextFont(42);
        l->SetTextSize(l->GetTextSize()*3.5);
        l->SetMargin(0.2);
        l->SetTextColor(1);

        int ctrfirst, ctrlast;
        if (subpl == 0){ // first 3 lines (0-5, 5-10, 10-20)
          ctrfirst = 0;
          ctrlast = 2;
        }else{ // last 3 lines (20-30, 30-50, 50-80)
          ctrfirst = 3;
          ctrlast = 5;
        }
        for (unsigned int jctr = ctrfirst; jctr <= ctrlast; jctr++){
          // h2d = (TH2D*) f[jctr]->Get(hist_names[jctr][ksign][lgap].c_str());
          // h[jctr-ctrfirst][ksign][lgap] = (TH1D*) h2d->ProjectionY(); // Dphi distribution with no |Deta| restriction
          
          std::string hist_name = ctr_dirs[jctr] + "h_" + kin + "_dr3" + ctrs[jctr] + signs[ksign] + gapcuts[lgap];

          if (projx_2d){
            h2d = (TH2D*) f[jctr]->Get(hist_name.c_str());
            h[jctr-ctrfirst][ksign][lgap] = (TH1D*) h2d->ProjectionX(("h_" + kin1d + ctrs[jctr] + signs[ksign] + gapcuts[lgap]).c_str());
          }
          else if (projy_2d){
            h2d = (TH2D*) f[jctr]->Get(hist_name.c_str());
            h[jctr-ctrfirst][ksign][lgap] = (TH1D*) h2d->ProjectionY(("h_" + kin1d + ctrs[jctr] + signs[ksign] + gapcuts[lgap]).c_str());
          }else{ // 1D
            h[jctr-ctrfirst][ksign][lgap] = (TH1D*) f[jctr]->Get(hist_name.c_str());
            // std::cout << kin1d << ", " << subpl_titles[imc][ksign] << ", " << origin_catg_labels[mgrp] << ", dphi " << ldphi << ", integral: " << h[imc][ksign][mgrp]->Integral() << std::endl;
          }

          // h[jctr-ctrfirst][ksign][lgap] = (TH1D*) h2d->ProjectionY("",81,120); // Dphi distribution for |Deta| < 0.96
        
          // if (jctr == 2)
          //   std::cout << "Bin counts in Dphi ~ pi and Dphi ~ 0: " << h[2][0][lgap]->GetBinContent(128) << ", " << h[2][0][lgap]->GetBinContent(64) << std::endl;

          h[jctr-ctrfirst][ksign][lgap]->Scale(1./1000.);
          hist_helper(h[jctr-ctrfirst][ksign][lgap], norm_factor[jctr], "#frac{1}{T_{AA}} #frac{1}{N_{coll}} #frac{dN}{d " + kin_title + "} [nb]");
          l->AddEntry(h[jctr-ctrfirst][ksign][lgap], ctrTitles[jctr].c_str(),"lp");
          h[jctr-ctrfirst][ksign][lgap]->SetMarkerStyle(20);
          h[jctr-ctrfirst][ksign][lgap]->SetMarkerSize(0.6);
        }

        l->AddEntry(h[3][ksign][lgap], ctrTitles[nCtrBins-1].c_str(),"lp");
        h[3][ksign][lgap]->SetMarkerStyle(20);
        h[3][ksign][lgap]->SetMarkerSize(0.6);
                  
        h[0][ksign][lgap]->SetMarkerColor(kRed);
        h[0][ksign][lgap]->SetLineColor(kRed);
        h[1][ksign][lgap]->SetMarkerColor(kBlue);
        h[1][ksign][lgap]->SetLineColor(kBlue);
        h[2][ksign][lgap]->SetMarkerColor(kGreen+2);
        h[2][ksign][lgap]->SetLineColor(kGreen+2);
        // first plot whichever with largest normalized peak at Delta phi ~ pi
        // since this sets the maximum y value

        float ymax = h[0][ksign][lgap]->GetMaximum();
        for (int iline = 1; iline < nLines; iline++){
          ymax = (ymax > h[iline][ksign][lgap]->GetMaximum())? ymax : h[iline][ksign][lgap]->GetMaximum();
          std::cout << "gapcut" << lgap+1 << ", sign" << ksign+1 << "subplot" << subpl+1 << ", line" << iline+1 << ", total integral: " << h[iline][ksign][lgap]->Integral("width") << std::endl;
        }
        h[0][ksign][lgap]->GetYaxis()->SetRangeUser(0., ymax * 1.1);

        h[0][ksign][lgap]->Draw("E");
        for (int iline = 1; iline < nLines; iline++){
          h[iline][ksign][lgap]->Draw("E,same");
        }

        l->AddEntry("",(gapcutTitles[lgap] + ", " + signTitles[ksign]).c_str(),"");
        l->Draw();
        // vl.push_back(l);
      }
    }

    c->SaveAs(Form("/usatlas/u/yuhanguo/usatlasdata/plots/pbpb_pp_compr/%s_pp_pbpb_compr_gapcut%d_atlasstyle.png", kin1d.c_str(), lgap+1));
    c->Close();
    delete c;
  }
}

void plot_pp_pbpb_compr(){
  plot_pp_pbpb_compr_single_kin("minv_pair_pt", false, true, "minv", "m_{#mu#mu}");
  // plot_pp_pbpb_compr_single_kin("Deta_Dphi", true, false, "DDDphi", "#Delta #phi");
  // plot_origin_categorized_kinematic_single("Deta_Dphi", false, true, false, false, "Deta", "#Delta #eta");


}


