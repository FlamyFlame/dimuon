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
const int nCtrBins = 7;
const int nLines = 4;
const int nSigns = 2;

double pi = acos(-1.0);
int nbins[nCtrBins] = {128,128,128,128,128,128,128};
int plateau_list[nGapCuts][nCtrBins * nSigns] = 
{{1649, 1614, 1313, 1312, 1811, 1825, 975, 1004, 668, 733, 104, 139, 3587, 7615}, 
 {1268, 1247, 1006, 1010, 1394, 1413, 755, 777, 517, 570, 81, 108, 3415, 7200}};

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
  h->Rebin(2);
  h->Scale(norm / 1000.,"width"); // change unit from pb to nb
  if (ytitle.length() != 0){
    h->GetYaxis()->SetTitle(ytitle.c_str());
  }
  h->SetMarkerStyle(20);
  h->SetMarkerSize(4);
  h->SetLineWidth(4);
  h->GetYaxis()->SetLabelFont(43);
  h->GetYaxis()->SetLabelSize(38);
  h->GetYaxis()->SetLabelOffset(0.01);
  h->GetYaxis()->SetTitleFont(43);
  h->GetYaxis()->SetTitleSize(40);
  h->GetYaxis()->SetTitleOffset(1.5);
  h->GetXaxis()->SetLabelFont(43);
  h->GetXaxis()->SetLabelSize(38);
  h->GetXaxis()->SetLabelOffset(0.01);
  h->GetXaxis()->SetTitleFont(43);  
  h->GetXaxis()->SetTitleSize(40);
  h->GetXaxis()->SetTitleOffset(1);
  // h->GetYaxis()->SetTitleOffset(h->GetYaxis()->GetTitleOffset()*0.8);
    // h->GetXaxis()->SetTitleOffset(h->GetXaxis()->GetTitleOffset()*0.8);
  // h->GetXaxis()->SetTitleSize(h->GetXaxis()->GetTitleSize()*1.8);
    // h->GetXaxis()->SetLabelSize(h->GetXaxis()->GetLabelSize()*1.35);
    // h->GetYaxis()->SetTitleSize(h->GetYaxis()->GetTitleSize()*1.8);
    // h->GetYaxis()->SetLabelSize(h->GetYaxis()->GetLabelSize()*1.35);
    // h->GetZaxis()->SetLabelSize(h->GetZaxis()->GetLabelSize()*1.35);
}

void dphi_pp_pbpb_compr(){

  initialize();
  TH2D* h2d;

  for (int lgap = 0; lgap < nGapCuts; lgap++){
    TCanvas* c = new TCanvas(Form("c%d",lgap),Form("c%d",lgap),2500,2000);
    c->Divide(nSigns,2);

    for (unsigned int ksign = 0; ksign < nSigns; ksign++){
      ////////////////////////////////////////// pp //////////////////////////////////////////
      h2d = (TH2D*) f[nCtrBins-1]->Get(hist_names[nCtrBins - 1][ksign][lgap].c_str());
      // h[nLines-1][ksign][lgap] = (TH1D*) h2d->ProjectionX("",81,120);
      h[nLines-1][ksign][lgap] = (TH1D*) h2d->ProjectionX();

      for (int i = 1; i <= nbins[nCtrBins-1]; i++){ //important: bin number goes from 1 to 128, not 0 to 127
        int ni = h[nLines-1][ksign][lgap]->GetBinContent(i);
        if (ni > plateau_list[lgap][(nCtrBins-1) * nSigns + ksign]){ // no buffer region
          ni = plateau_list[lgap][(nCtrBins-1) * nSigns + ksign];
        }
        h[nLines-1][ksign][lgap]->AddBinContent(i,-ni);
      }

      hist_helper(h[nLines-1][ksign][lgap], norm_factor[nCtrBins-1], "#frac{1}{T_{AA}} #frac{1}{N_{coll}} #frac{dN}{d #Delta #phi} [nb]");


      for (unsigned int subpl = 0; subpl < 2; subpl++){
        c->cd(ksign * 2 + subpl + 1);
        gPad->SetLeftMargin(0.17);
        gPad->SetBottomMargin(0.135);

        // TLegend* l = new TLegend(0.18,0.64,0.53,0.87);
        TLegend* l = new TLegend(0.58,0.64,0.88,0.87);
        l->SetBorderSize(0);
        l->SetFillStyle(0);
        l->SetTextFont(42);
        l->SetMargin(0.15);
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
          h2d = (TH2D*) f[jctr]->Get(hist_names[jctr][ksign][lgap].c_str());
          h[jctr-ctrfirst][ksign][lgap] = (TH1D*) h2d->ProjectionX(); // Dphi distribution with no |Deta| restriction
          // h[jctr-ctrfirst][ksign][lgap] = (TH1D*) h2d->ProjectionX("",81,120); // Dphi distribution for |Deta| < 0.96
        
          for (int i = 1; i <= nbins[jctr]; i++){ //important: bin number goes from 1 to 128, not 0 to 127
            int ni = h[jctr-ctrfirst][ksign][lgap]->GetBinContent(i);
            // if (ni > plateau_list[lgap][jctr * nSigns + ksign] * 1.01){ // 1% buffer region
            if (ni > plateau_list[lgap][jctr * nSigns + ksign]){ // no buffer region
              ni = plateau_list[lgap][jctr * nSigns + ksign];
            }
            h[jctr-ctrfirst][ksign][lgap]->AddBinContent(i,-ni);
          }

          // if (jctr == 2)
          //   std::cout << "Bin counts in Dphi ~ pi and Dphi ~ 0: " << h[2][0][lgap]->GetBinContent(128) << ", " << h[2][0][lgap]->GetBinContent(64) << std::endl;

          hist_helper(h[jctr-ctrfirst][ksign][lgap], norm_factor[jctr], "#frac{1}{T_{AA}} #frac{1}{N_{coll}} #frac{dN}{d #Delta #phi} [nb]");
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

        // float ylim_arr[nCtrBins];
        // for (int iline = 0; iline < nLines; iline++){
        //   ylim_arr[iline] = 1.1 * std::max(h[iline][ksign][lgap]->GetBinContent(128), h[iline][ksign][lgap]->GetBinContent(1));
        //   // std::cout << "ctr-group "<< iline << ", max y value: " << ylim_arr[iline] << std::endl;
        //   std::cout << "gapcut" << lgap+1 << ", sign" << ksign+1 << "subplot" << subpl+1 << ", line" << iline+1 << ", total integral: " << h[iline][ksign][lgap]->Integral("width") << std::endl;
        // }
        // // float ylim = *std::max_element(ylim_arr, ylim_arr + nLines);
        // int ymax_ind = std::max_element(ylim_arr, ylim_arr + nLines) - ylim_arr;
        // // std::cout << "y-max = " << ylim << std::endl;

        // // h[2][ksign]->GetYaxis()->SetRange(0,ylim);
        // h[ymax_ind][ksign][lgap]->Draw("E");
        // for (int iline = 0; iline < nLines; iline++){
        //   if (iline != ymax_ind) h[iline][ksign][lgap]->Draw("E,same");
        // }

        // l->AddEntry("",(gapcutTitles[lgap] + ", " + signTitles[ksign]).c_str(),"");
        l->AddEntry("",(signTitles[ksign]).c_str(),"");
        l->Draw();
        // vl.push_back(l);
      }
    }

    c->SaveAs(Form("/usatlas/u/yuhanguo/usatlasdata/plots/pbpb_pp_compr/dphi_pp_pbpb_compr_gapcut%d.png",lgap+1));
    c->Close();
    delete c;
  }

}


