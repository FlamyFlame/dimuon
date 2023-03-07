#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <string.h>
#include  <stdlib.h>
#include <algorithm>
#include "vector"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "string"
#include "time.h"
#include "struct_hist.h"


// can consider putting the global variables in a class or struct

const int nSigns = 2; //3*2 canvas
const int nWeights = 2; //3*2 canvas
const int nMCs = 2; //bb, cc
const int nGapCuts = 2; //with, without gap cut
// const int ndRCuts = 3; //with, without gap cut

std::string weight_prefixs[nWeights] = {"h_unweighted_", "h_"};
// std::string dRs[ndRcuts] = {"_dr1","_dr2","_dr3"};
// std::string dphis[ndphicuts] = {"_dphi1","_dphi2","_dphi3"};
std::string signs[nSigns] = {"_sign1", "_sign2"};
std::string gapcuts[nGapCuts] = {"_gapcut1", "_gapcut2"};

// std::string dRTitles[ndRcuts] = {"#Delta R < 0.8","#Delta R < 1.2", "no #Delta R cut"};
// std::string dphiTitles[ndphicuts] = {"#Delta #phi < 1","#Delta #phi > #pi - 1", "no #Delta #phi cut"};
std::string weightTitles[nSigns] = {"unweighted", "weighted"};
std::string signTitles[nSigns] = {"same sign", "opposite sign"};
std::string mcTitles[nMCs] = {"MC truth bb", "MC truth cc"};
std::string gapCutTitles[nGapCuts] = {"no gap cut", "with gap cut"};

std::string png_title_1D[nGapCuts] = {"_mc_truth_no_gap_cut.png", "_mc_truth_with_gap_cut.png"};


void hist_helper(TH1D& h, std::string title, int scale_to_unity, bool enhance_by_100 = false, std::string ytitle=""){

	h.SetTitle(title.c_str());
	h.SetStats(0);
	if (enhance_by_100) h.Scale(100,"width");
	else 				h.Scale(1,"width");
	if(scale_to_unity){
		h.Scale(1./h.Integral("width"));
	}
	if (ytitle.length() != 0){
		h.GetYaxis()->SetTitle(ytitle.c_str());
	}
	h.GetXaxis()->SetNdivisions(505);
	h.GetYaxis()->SetNdivisions(505);
	h.GetYaxis()->SetLabelFont(43);
    h.GetYaxis()->SetLabelSize(36);
    h.GetYaxis()->SetLabelOffset(0.02);
    h.GetYaxis()->SetTitleFont(43);
    h.GetYaxis()->SetTitleSize(36);
    h.GetYaxis()->SetTitleOffset(1.5);
	h.GetXaxis()->SetLabelFont(43);
    h.GetXaxis()->SetLabelSize(36);
    h.GetXaxis()->SetLabelOffset(0.02);
    h.GetXaxis()->SetTitleFont(43);
    h.GetXaxis()->SetTitleSize(36);
    h.GetXaxis()->SetTitleOffset(1);

    h.SetMarkerStyle(20);
    h.SetMarkerSize(0.6);
}


void plot_one_scaling_mode(bool scale_to_unity){

	TCanvas *c[nGapCuts];
	TH2D* h2[nMCs];
	// TH1D* hx[nMCs];
	// TH1D hx[nGapCuts][nWeights][nSigns][nMCs];
	TH1D (*hx)[nWeights][nSigns][nMCs] = new TH1D[nGapCuts][nWeights][nSigns][nMCs];

	TFile* f[nMCs];
	f[0] = TFile::Open("/usatlas/u/yuhanguo/usatlasdata/athena/runMCV2/histograms_mc_truth_bb.root"); 
	f[1] = TFile::Open("/usatlas/u/yuhanguo/usatlasdata/athena/runMCV2/histograms_mc_truth_cc.root"); 


	// for (unsigned int mgapcut = 0; mgapcut < nGapCuts; mgapcut++){
	for (unsigned int mgapcut = 0; mgapcut < 1; mgapcut++){
		char name[100];
		sprintf(name, "c_dphi_gapcut%d", mgapcut);
		c[mgapcut] = new TCanvas(name,name,1500,1200);
		c[mgapcut]->Divide(nSigns,2);

		// std::vector<TLegend*> vl;
		// std::vector<TH1*> vh;
		for (unsigned int jweight = 0; jweight < nWeights; jweight++){
			for (unsigned int ksign = 0; ksign < nSigns; ksign++){
				std::string histName;
				histName = weight_prefixs[jweight] + "Deta_Dphi_dr3" + signs[ksign] + gapcuts[mgapcut];
				
				c[mgapcut]->cd(jweight * 2 + ksign + 1);
				gPad->SetLeftMargin(0.16);
				gPad->SetBottomMargin(0.155);
				// gPad->SetRightMargin(0.135);

				TLegend* l = new TLegend(0.18,0.75,0.5,0.88);
		    	l->SetBorderSize(0);
		    	l->SetFillStyle(0);
		    	l->SetTextFont(43);
		    	l->SetMargin(0.02);
		    	l->SetTextColor(1);


				// TH2D* h2_bb;
				// TH2D* h2_cc;
				// TH1D* hdphi_bb;
				// TH1D* hdphi_cc;
				for (int lmc = 0; lmc < nMCs; lmc++){

					h2[lmc] = (TH2D*) f[lmc]->Get(histName.c_str());
					if (h2[lmc]->GetEntries() == 0){
					    std::cout<<"Error:: Read in histogram with name " << histName << " is empty,  quitting"<<std::endl;
		                throw std::exception();
		            }
					// hx[lmc] = h2[lmc]->ProjectionX(Form("hx%d",lmc+1));
					hx[mgapcut][jweight][ksign][lmc] = *(h2[lmc]->ProjectionX(Form("hx%d%d%d%d",mgapcut,jweight,ksign,lmc)));
					hx[mgapcut][jweight][ksign][lmc].Rebin(2);
					if (jweight == 1 && lmc == 1 && !scale_to_unity) hist_helper(hx[mgapcut][jweight][ksign][lmc], weightTitles[jweight] + ", " + signTitles[ksign], scale_to_unity, true);
					else 											 hist_helper(hx[mgapcut][jweight][ksign][lmc], weightTitles[jweight] + ", " + signTitles[ksign], scale_to_unity);
				}


				// h2_bb = (TH2D*) f[0]->Get(histName.c_str());
				// h2_cc = (TH2D*) f[1]->Get(histName.c_str());

				// hdphi_bb = h2_bb->ProjectionX("hdphi_bb");
				// hdphi_cc = h2_cc->ProjectionX("hdphi_cc");
				// hist_helper(hdphi_bb, weightTitles[jweight] + ", " + signTitles[ksign], scale_to_unity);
				// hist_helper(hdphi_cc, weightTitles[jweight] + ", " + signTitles[ksign], scale_to_unity);

				// cout << h2[0]->GetEntries() << ", "	<< h2[0]->Integral() << ", " << hx[mgapcut][jweight][ksign][0].GetEntries() << ", " << hx[mgapcut][jweight][ksign][0].Integral("width") << std::endl;
				// cout << h2[1]->GetEntries() << ", "	<< h2[1]->Integral() << ", " << hx[mgapcut][jweight][ksign][1].GetEntries() << ", " << hx[mgapcut][jweight][ksign][1].Integral("width") << std::endl;

				int ymax = 1.1 * (std::max(hx[mgapcut][jweight][ksign][0].GetMaximum(), hx[mgapcut][jweight][ksign][1].GetMaximum()));
				int ymin = 0.9 * (std::min(hx[mgapcut][jweight][ksign][0].GetMinimum(), hx[mgapcut][jweight][ksign][1].GetMinimum()));
				hx[mgapcut][jweight][ksign][0].GetYaxis()->SetRangeUser(0,ymax);
				hx[mgapcut][jweight][ksign][1].GetYaxis()->SetRangeUser(0,ymax);
				hx[mgapcut][jweight][ksign][0].SetLineColor(kRed);
				hx[mgapcut][jweight][ksign][0].SetMarkerColor(kRed);
				hx[mgapcut][jweight][ksign][1].SetLineColor(kBlue);
				hx[mgapcut][jweight][ksign][1].SetMarkerColor(kBlue);

			    l->AddEntry(&hx[mgapcut][jweight][ksign][0],mcTitles[0].c_str(),"lp");
			    if (jweight == 1 && !scale_to_unity) l->AddEntry(&hx[mgapcut][jweight][ksign][1],(mcTitles[1] + " * 100").c_str(),"lp");
			    else 								 l->AddEntry(&hx[mgapcut][jweight][ksign][1],mcTitles[1].c_str(),"lp");
				
			    // l->AddEntry("",gapCutTitles[mgapcut].c_str(),"");
			 //    hdphi_bb->Draw("E");
				// hdphi_cc->Draw("E,same");
			    hx[mgapcut][jweight][ksign][0].Draw("E");
				hx[mgapcut][jweight][ksign][1].Draw("E,same");
				l->Draw();
							
				// vl.push_back(l);
				// vh.push_back(h2[0]);
				// vh.push_back(h2[1]);
				// vh.push_back(hx[0]);
				// vh.push_back(hx[1]);
			}
		}

		if (scale_to_unity) c[mgapcut]->SaveAs(("plots/mc_truth/no_prt_grouping/h_DPHI_unity" + png_title_1D[mgapcut]).c_str());
		else 				c[mgapcut]->SaveAs(("plots/mc_truth/no_prt_grouping/h_DPHI_absolute" + png_title_1D[mgapcut]).c_str());
  		c[mgapcut]->Close();
  		delete c[mgapcut];
  		// for (auto i = vl.begin(); i != vl.end(); i++){
  		// 	if (*i != nullptr) delete *i;
  		// }
  		// for (auto i = vh.begin(); i != vh.end(); i++){
  		// 	if (*i != nullptr) delete *i;
  		// }
  		// vl.clear();
  		// vh.clear();
	}

	delete f[0];
  	delete f[1];

}


void plot_dphi_mc_truth(){

	plot_one_scaling_mode(true);
	plot_one_scaling_mode(false);

}


