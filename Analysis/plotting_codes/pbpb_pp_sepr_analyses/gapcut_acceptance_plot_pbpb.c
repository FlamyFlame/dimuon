#include <string.h>
#include  <stdlib.h>
#include "vector"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include<algorithm>
#include "struct_hist.h"
#define NPLOTS 4

const int ndRs = 3;
const int nSigns = 2;
const int nCtrBins = 5;
const int nGapCuts = 2;


void hist_helper(TH1* h, bool scale, std::string ytitle=""){

	h->SetStats(0);
	if(scale){
		h->Scale(1,"width");
		h->Scale(1./h->Integral("width"));
	}
	if (ytitle.length() != 0){
		h->GetYaxis()->SetTitle(ytitle.c_str());
	}
	h->GetYaxis()->SetLabelFont(43);
    h->GetYaxis()->SetLabelSize(32);
    h->GetYaxis()->SetLabelOffset(0.02);
    h->GetYaxis()->SetTitleFont(43);
    h->GetYaxis()->SetTitleSize(32);
    h->GetYaxis()->SetTitleOffset(1.5);
	h->GetXaxis()->SetLabelFont(43);
    h->GetXaxis()->SetLabelSize(32);
    h->GetXaxis()->SetLabelOffset(0.02);
    h->GetXaxis()->SetTitleFont(43);
    h->GetXaxis()->SetTitleSize(32);
    h->GetXaxis()->SetTitleOffset(1);
}

void gapcut_acceptance_plot_pbpb(){
	std::string dpath = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/";
	std::string scr_file = "histograms_scrambled_pairs.root";
	TFile* f = TFile::Open((dpath + scr_file).c_str());

	std::string dRs[ndRs] = {"_dr1","_dr2","_dr3"};
	std::string dRTitles[ndRs] = {"#Delta R < 0.8","#Delta R < 1.2", "no #Delta R cut"};

	std::string signs[nSigns] = {"_sign1", "_sign2"};
	std::string signTitles[nSigns] = {"same sign", "opposite sign"};

	std::string ctrs[nCtrBins] = {"_ctr1","_ctr2","_ctr3","_ctr4","_ctr5"};
	std::string ctrTitles[nCtrBins] = {"Centrality 0-20", "Centrality 20-40", "Centrality 40-60", "Centrality 60-80", "Centrality 80-100"};

	std::string gapcuts[nGapCuts] = {"_gapcut1", "_gapcut2"};

	std::vector<std::string> png_title_list[nCtrBins];
	for (int ictr = 0; ictr < nCtrBins; ictr++){
	png_title_list[ictr].push_back(("acceptance_dR" + ctrs[ictr] + ".png").c_str());
	png_title_list[ictr].push_back(("acceptance_Deta" + ctrs[ictr] + ".png").c_str());
	png_title_list[ictr].push_back(("acceptance_eta_avg" + ctrs[ictr] + ".png").c_str());
	png_title_list[ictr].push_back(("acceptance_Dphi" + ctrs[ictr] + ".png").c_str());
	}

	std::vector<TH1D*> h_list[ndRs][nCtrBins][nSigns];

	for (unsigned int jdr = 0; jdr < ndRs; jdr++){
		for (unsigned int ictr = 0; ictr < nCtrBins; ictr++){
			for (unsigned int ksign = 0; ksign < nSigns; ksign++){
				TH1D* h_dR = (TH1D*) f->Get(("ctr-binned/h_DR" + dRs[jdr] + ctrs[ictr] + signs[ksign] + "_gapcut2").c_str()); // after gap cut
				TH1D* h_dR_before = (TH1D*) f->Get(("ctr-binned/h_DR" + dRs[jdr] + ctrs[ictr] + signs[ksign] + "_gapcut1").c_str()); // after gap cut
				h_dR->Divide(h_dR_before);
				h_list[jdr][ictr][ksign].push_back(h_dR);

				TH2D* h_eta_avg_Deta = (TH2D*) f->Get(("ctr-binned/h_eta_avg_Deta" + dRs[jdr] + ctrs[ictr] + signs[ksign] + "_gapcut2").c_str());
				TH2D* h_eta_avg_Deta_before = (TH2D*) f->Get(("ctr-binned/h_eta_avg_Deta" + dRs[jdr] + ctrs[ictr] + signs[ksign] + "_gapcut1").c_str());
				TH1D* h_Deta = h_eta_avg_Deta->ProjectionX();
				TH1D* h_eta_avg = h_eta_avg_Deta->ProjectionY();
				TH1D* h_Deta_before = h_eta_avg_Deta_before->ProjectionX();
				TH1D* h_eta_avg_before = h_eta_avg_Deta_before->ProjectionY();
				h_Deta->Divide(h_Deta_before);
				h_eta_avg->Divide(h_eta_avg_before);
				h_list[jdr][ictr][ksign].push_back(h_Deta);
				h_list[jdr][ictr][ksign].push_back(h_eta_avg);

				TH2D* h_eta_avg_Dphi = (TH2D*) f->Get(("ctr-binned/h_eta_avg_Dphi" + dRs[jdr] + ctrs[ictr] + signs[ksign] + "_gapcut2").c_str());
				TH2D* h_eta_avg_Dphi_before = (TH2D*) f->Get(("ctr-binned/h_eta_avg_Dphi" + dRs[jdr] + ctrs[ictr] + signs[ksign] + "_gapcut1").c_str());
				TH1D* h_Dphi = h_eta_avg_Dphi->ProjectionX();
				TH1D* h_Dphi_before = h_eta_avg_Dphi_before->ProjectionX();
				h_Dphi->Divide(h_Dphi_before);
				h_list[jdr][ictr][ksign].push_back(h_Dphi);

				assert(h_list[jdr][ictr][ksign].size == NPLOTS);
			}
		}
	}
	

	for (int iplt = 0; iplt < NPLOTS; iplt++){
		for (unsigned int ictr = 0; ictr < nCtrBins; ictr++){

			char name[100];
			sprintf(name, "ratio_%s", png_title_list[ictr][iplt].c_str());
			TCanvas *c = new TCanvas(name,name,4200,2000);
			c->Divide(ndRs,nSigns);

			std::vector<TLegend*> vl;

			for (unsigned int ksign = 0; ksign < nSigns; ksign++){
				for (unsigned int jdr = 0; jdr < ndRs; jdr++){
							
					c->cd(ksign * ndRs + jdr + 1);
					gPad->SetLeftMargin(0.16);
					gPad->SetBottomMargin(0.135);

					TLegend* l = new TLegend(0.22,0.7,0.5,0.87);
		    		l->SetBorderSize(0);
		    		l->SetFillStyle(0);
		    		l->SetTextFont(42);
		    		l->SetMargin(0.2);
		    		l->SetTextColor(1);
		    		vl.push_back(l);

					// std::string subplot_title = cutTitles[jdr] + ", " + signTitles[ksign];
					hist_helper(h_list[jdr][ictr][ksign][iplt], false, "Ratio of events passing gapcut");
					l->AddEntry(h_list[jdr][ictr][ksign][iplt],(dRTitles[jdr] + ", " + signTitles[ksign]).c_str(),"lp");
					l->AddEntry("",ctrTitles[ictr].c_str(),"");
					h_list[jdr][ictr][ksign][iplt]->SetMarkerStyle(20);
					h_list[jdr][ictr][ksign][iplt]->SetMarkerSize(0.5);
					h_list[jdr][ictr][ksign][iplt]->Draw();
					l->Draw();

					// for (auto i = h_list[jdr][ictr][ksign].begin(); i != h_list[jdr][ictr][ksign].end(); i++){
	  		// 			if (*i != nullptr) delete *i;
	  		// 		}
	  		// 		h_list[jdr][ictr][ksign].clear();
				}
			}
			
			c->SaveAs(("plots/real_scramb_comparison/gapcut_acceptance/" + png_title_list[ictr][iplt]).c_str());
		  	c->Close();
		  	delete c;
					
			for (auto i = vl.begin(); i != vl.end(); i++){
	  			if (*i != nullptr) delete *i;
	  		}
	  		vl.clear();	
	  	}
  	}
 }

