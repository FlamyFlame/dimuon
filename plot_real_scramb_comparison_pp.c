#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <string.h>
#include  <stdlib.h>
#include "vector"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "string"
#include "time.h"
#include "struct_hist.h"


// can consider putting the global variables in a class or struct



const int ndRcuts = 3;
const int ndphicuts = 3;
const int nSigns = 2; //3*2 canvas
const int nGroups = 2; //real, scrambled pairs
const int nGapCuts = 2; //with, without gap cut

std::vector<std::string> hist1DNames = {"h_pair_dP_overP", "h_pair_y","h_DR"};
std::vector<bool> logx1Ds = {false,false,false};
std::vector<bool> logy1Ds = {false,false,false};

std::vector<std::string> hist2DNames = {"h_eta_avg_Dphi", "h_eta1_eta2", "h_eta_avg_Deta", "h_pt1_pt2","h_ptlead_pair_pt", "h_minv_pair_pt"};
std::vector<bool> projxs = {true, true, true, true, true, false};
std::vector<bool> projys = {false, true, true, true, true, true};
std::vector<bool> logxs = {false, false, false, true, true, true};
std::vector<bool> logys = {false, false, false, true, true, true};
std::vector<bool> logzs = {false, false, false, true, true, true};
std::vector<std::string> hxNames = {"h_Dphi", "h_eta2", "h_Deta", "h_pt2", "h_pair_pt", ""};
std::vector<std::string> hyNames = {"", "h_eta1", "h_eta_avg", "h_pt1", "h_ptlead", "h_minv"};

const int n2Ds = hist2DNames.size();
const int n1Ds = hist1DNames.size();

std::string dRs[ndRcuts] = {"_dr1","_dr2","_dr3"};
std::string dphis[ndphicuts] = {"_dphi1","_dphi2","_dphi3"};
std::string signs[nSigns] = {"_sign1", "_sign2"};
std::string gapcuts[nGapCuts] = {"_gapcut1", "_gapcut2"};

std::string dRTitles[ndRcuts] = {"#Delta R < 0.8","#Delta R < 1.2", "no #Delta R cut"};
std::string dphiTitles[ndphicuts] = {"#Delta #phi < 1","#Delta #phi > #pi - 1", "no #Delta #phi cut"};
std::string signTitles[nSigns] = {"same sign", "opposite sign"};
std::string grpTitles[nGroups] = {"real muon pairs", "scrambled pairs"};
std::string gapCutTitles[nGroups] = {"no gap cut", "with gap cut"};

std::string png_title_1D[nGapCuts];
std::string png_title_2D[nGroups][nGapCuts];
std::string subplot_titles[ndRcuts][nSigns];
std::string subplot_dphi_titles[ndRcuts][nSigns];
// std::string subplot_titles_1D[ndRcuts][nSigns];
// std::string subplot_titles_2D[ndRcuts][nSigns][nGroups];

TFile* f[nGroups];

void initialize(){
	f[0] = TFile::Open("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/histograms_real_pairs_pp.root"); //real pairs
	f[1] = TFile::Open("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/histograms_scrambled_pairs_pp.root"); //scrambled pairs

	for (unsigned int ksign = 0; ksign < nSigns; ksign++){
		for (unsigned int jdr = 0; jdr < ndRcuts; jdr++){  	
			subplot_titles[jdr][ksign] = dRTitles[jdr] + ", " + signTitles[ksign];
		}
		for (unsigned int jdphi = 0; jdphi < ndRcuts; jdphi++){  	
			subplot_dphi_titles[jdphi][ksign] = dRTitles[jdphi] + ", " + signTitles[ksign];
		}

	}

	png_title_1D[0] = "_real_scr_no_gap_cut_pp.png";
	png_title_1D[1] = "_real_scr_with_gap_cut_pp.png";
	png_title_2D[0][0] = "_real_no_gap_cut_pp.png";
	png_title_2D[0][1] = "_real_with_gap_cut_pp.png";
	png_title_2D[1][0] = "_scr_no_gap_cut_pp.png";
	png_title_2D[1][1] = "_scr_with_gap_cut_pp.png";
}

void hist_helper(TH1* h, std::string title, bool scale, std::string ytitle=""){

	h->SetTitle(title.c_str());
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

void plot1D(Hist1D& h, const int nCuts, std::string cuts[], std::string cutTitles[], std::string name_specf=""){
	h.name_specifier = name_specf;
	
	TCanvas *c;

	for (int mgapcut = 0; mgapcut < nGapCuts; mgapcut++){

		char name[100];
		sprintf(name, "c_%s_gapcut%d", h.name.c_str(), mgapcut);
		c = new TCanvas(name,name,4200,2000);
		c->Divide(nCuts,nSigns);

		std::vector<TLegend*> vl;
		std::vector<TH1D*> vh;

		for (unsigned int ksign = 0; ksign < nSigns; ksign++){
			for (unsigned int jcut = 0; jcut < nCuts; jcut++){
					
				c->cd(ksign * nCuts + jcut + 1);
				gPad->SetLogx(h.logx);
				gPad->SetLogy(h.logy);
				gPad->SetLeftMargin(0.16);
				gPad->SetBottomMargin(0.135);
				TLegend* l = new TLegend(0.2,0.7,0.5,0.87);
	   			l->SetBorderSize(0);
	   			l->SetFillStyle(0);
	   			l->SetTextFont(42);
	   			l->SetMargin(0.02);
	   			l->SetTextColor(1);

				std::string histName = h.name + cuts[jcut]+ signs[ksign] + gapcuts[mgapcut];
	    		TH1D* h1[nGroups];
				for (int igrp = 0; igrp < nGroups; igrp++){
					h1[igrp] = (TH1D*) f[igrp]->Get(histName.c_str());
					std::string subplot_title = cutTitles[jcut] + ", " + signTitles[ksign];
					hist_helper(h1[igrp],subplot_title, true);
	    			l->AddEntry(h1[igrp],grpTitles[igrp].c_str(),"lp");
					h1[igrp]->SetMarkerStyle(20);
					h1[igrp]->SetMarkerSize(0.5);
				}
	    		l->AddEntry("pp","");
		    	l->AddEntry("",gapCutTitles[mgapcut].c_str(),"");
						
				h1[0]->SetMarkerColor(kRed);
				h1[0]->SetLineColor(kRed);
				h1[0]->Draw("E");
				h1[1]->Draw("E,same");
				l->Draw();
				
				vl.push_back(l);
				vh.push_back(h1[0]);
				vh.push_back(h1[1]);
			}
		}

		c->SaveAs(("plots/real_scramb_comparison_pp/" + h.name + h.name_specifier + png_title_1D[mgapcut]).c_str());
	  	c->Close();
	  	delete c;

	  	for (auto i = vl.begin(); i != vl.end(); i++){
  			if (*i != nullptr) delete *i;
  		}
  		for (auto i = vh.begin(); i != vh.end(); i++){
  			if (*i != nullptr) delete *i;
  		}
  		vl.clear();
  		vh.clear();
	}
}

// void plot1D(TH1D* h[nGapCuts][ndRcuts][nSigns][nGroups], std::string hist1DName, bool logx, bool logy, std::string name_specifier = ""){
void plot1D(TH1D (*h)[nGapCuts][nSigns][nGroups], Hist1D& h_struct, const int nCuts, std::string cuts[], std::string cutTitles[]){
	TCanvas *c;

		for (int mgapcut = 0; mgapcut < nGapCuts; mgapcut++){

			char name[100];
			sprintf(name, "c_%s_gapcut%d", h_struct.name.c_str(), mgapcut);
			c = new TCanvas(name,name,4200,2000);
			c->Divide(nCuts,nSigns);

			std::vector<TLegend*> vl;

			for (unsigned int ksign = 0; ksign < nSigns; ksign++){
				for (unsigned int jcut = 0; jcut < nCuts; jcut++){
						
					c->cd(ksign * nCuts + jcut + 1);
					gPad->SetLogx(h_struct.logx);
					gPad->SetLogy(h_struct.logy);
					gPad->SetLeftMargin(0.16);
					gPad->SetBottomMargin(0.135);

					TLegend* l = new TLegend(0.22,0.7,0.5,0.87);
	    			l->SetBorderSize(0);
	    			l->SetFillStyle(0);
	    			l->SetTextFont(42);
	    			l->SetMargin(0.02);
	    			l->SetTextColor(1);

	    			// the gap cut situation already handled
	    			TH1D* h1[nGroups] = {&h[jcut][mgapcut][ksign][0], &h[jcut][mgapcut][ksign][1]};
					for (int igrp = 0; igrp < nGroups; igrp++){
						std::string subplot_title = cutTitles[jcut] + ", " + signTitles[ksign];
						hist_helper(h1[igrp],subplot_title, true);
	    				l->AddEntry(h1[igrp],grpTitles[igrp].c_str(),"lp");
						h1[igrp]->SetMarkerStyle(20);
						h1[igrp]->SetMarkerSize(0.5);
					}
	    			l->AddEntry("","pp","");
		    		l->AddEntry("",gapCutTitles[mgapcut].c_str(),"");
							
					h1[0]->SetMarkerColor(kRed);
					h1[0]->SetLineColor(kRed);
					h1[0]->Draw("E");
					h1[1]->SetMarkerColor(kBlue);
					h1[1]->SetLineColor(kBlue);
					h1[1]->Draw("E,same");
					l->Draw();
					vl.push_back(l);
				}
			}

			c->SaveAs(("plots/real_scramb_comparison_pp/" + h_struct.name + h_struct.name_specifier + png_title_1D[mgapcut]).c_str());
	  		c->Close();
	  		delete c;

	  		for (auto i = vl.begin(); i != vl.end(); i++){
  				if (*i != nullptr) delete *i;
  			}
  			vl.clear();

		}	
}

void Process2D(Hist2D h, const int nCuts, std::string cuts[], std::string cutTitles[], std::string name_specf=""){
	
	h.name_specifier = name_specf;
	// TH1D (*hx)[nGapCuts][nCuts][nSigns][nGroups] = new TH1D[nGapCuts][nCuts][nSigns][nGroups];
	// TH1D (*hy)[nGapCuts][nCuts][nSigns][nGroups] = new TH1D[nGapCuts][nCuts][nSigns][nGroups];
	TH1D (*hx)[nGapCuts][nSigns][nGroups] = new TH1D[nCuts][nGapCuts][nSigns][nGroups];
	TH1D (*hy)[nGapCuts][nSigns][nGroups] = new TH1D[nCuts][nGapCuts][nSigns][nGroups];

	TCanvas *c[nGapCuts][nGroups];

	for (unsigned int mgapcut = 0; mgapcut < nGapCuts; mgapcut++){
		for (int lgrp = 0; lgrp < nGroups; lgrp++){
			char name[100];
			sprintf(name, "c_%s_gapcut%d_group%d", h.name.c_str(), mgapcut,lgrp);
			c[mgapcut][lgrp] = new TCanvas(name,name,4200,2000);
			c[mgapcut][lgrp]->Divide(nCuts,nSigns);

			std::vector<TLegend*> vl;
			std::vector<TH2D*> vh;

			for (unsigned int ksign = 0; ksign < nSigns; ksign++){
				for (unsigned int jcut = 0; jcut < nCuts; jcut++){
					c[mgapcut][lgrp]->cd(ksign * nCuts + jcut + 1);
					gPad->SetLogx(h.logx);
					gPad->SetLogy(h.logy);
					gPad->SetLeftMargin(0.16);
					gPad->SetBottomMargin(0.135);
					gPad->SetLogz(h.logz);
					gPad->SetRightMargin(0.135);

					TLegend* l = new TLegend(0.22,0.82,0.46,0.87);
	    			l->SetBorderSize(0);
	    			l->SetFillStyle(0);
	    			l->SetTextFont(43);
	    			l->SetMargin(0.02);
	    			l->SetTextColor(1);
	    			l->AddEntry("","pp","");
		    		l->AddEntry("",gapCutTitles[mgapcut].c_str(),"");
		    		l->AddEntry("",grpTitles[lgrp].c_str(),"");
					l->Draw();
						
					std::string histName = h.name + cuts[jcut] + signs[ksign] + gapcuts[mgapcut];
					TH2D* h2 = (TH2D*) f[lgrp]->Get(histName.c_str());
					if (h2->GetEntries() == 0){
					    std::cout<<"Error:: Read in histogram with name " << histName << " is empty,  quitting"<<std::endl;
                   		throw std::exception();
                   	}

					hx[jcut][mgapcut][ksign][lgrp] = *(h2->ProjectionX());
					hy[jcut][mgapcut][ksign][lgrp] = *(h2->ProjectionY());

					// hist_helper(h2,subplot_titles[jcut][ksign], true);
					std::string subplot_title = cutTitles[jcut] + ", " + signTitles[ksign];
					hist_helper(h2,subplot_title, true);
					h2->Draw("COLZ");

					vl.push_back(l);
					vh.push_back(h2);
				}
			}
			c[mgapcut][lgrp]->SaveAs(("plots/real_scramb_comparison_pp/" + h.name + h.name_specifier + png_title_2D[lgrp][mgapcut]).c_str());
  			c[mgapcut][lgrp]->Close();
  			delete c[mgapcut][lgrp];

  			for (auto i = vl.begin(); i != vl.end(); i++){
  				if (*i != nullptr) delete *i;
  			}
  			for (auto i = vh.begin(); i != vh.end(); i++){
  				if (*i != nullptr) delete *i;
  			}
  			vl.clear();
  			vh.clear();
		}
	}

	if (h.projx){
		Hist1D hx_struct = {h.hxName, h.logx, h.logz, h.name_specifier};
		plot1D(hx, hx_struct, nCuts, cuts, cutTitles);
	}
	if (h.projy){
		Hist1D hy_struct = {h.hyName, h.logy, h.logz, h.name_specifier};
		plot1D(hy, hy_struct, nCuts, cuts, cutTitles);
	}

	delete[] hx;
	delete[] hy;
}

void plot_real_scramb_comparison_pp(){

	clock_t start, end;
    double cpu_time_used;
    start = clock();

	initialize();


	// for (unsigned int i2d = 0; i2d < n2Ds; i2d++){
	// 	Hist2D h2d = {hist2DNames[i2d], projxs[i2d], projys[i2d], hxNames[i2d], hyNames[i2d], logxs[i2d], logys[i2d], logzs[i2d]};
	// 	Process2D(h2d, ndRcuts, dRs, dRTitles);
	// }

	for (unsigned int i1d = 0; i1d < n1Ds; i1d++){
		Hist1D h1d = {hist1DNames[i1d], logx1Ds[i1d], logy1Ds[i1d]};
		plot1D(h1d, ndRcuts, dRs, dRTitles);
	}

	Hist2D h2d_minvnolog = {"h_minv_pair_pt",false,true,"h_pair_pt","h_minv",true,false,false};
	Process2D(h2d_minvnolog,ndRcuts, dRs, dRTitles,"_MINVNOLOG_");

	Hist2D h2d_phi = {"h_eta1_eta2", true, true, "h_eta1", "h_eta2", false, false, false};
	Process2D(h2d_phi, ndphicuts, dphis, dphiTitles, "_DPHI_");

	delete f[0];
  	delete f[1];

  	end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;;
    std::cout << "CPU time used is " << cpu_time_used << " seconds" << std::endl;

}

