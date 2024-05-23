#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <string.h>
#include "vector"
#include "TH1.h"
#include "TH2.h"

void hist_helper(TH1* h, std::string title, bool scale, std::string ytitle=""){
	h->SetTitle(title.c_str());
	h->SetStats(0);
	if(scale){
		h->Scale(1./h->Integral("width"));
	}
	if (ytitle.length() != 0){
		h->GetYaxis()->SetTitle(ytitle.c_str());
	}
	h->GetYaxis()->SetTitleOffset(h->GetYaxis()->GetTitleOffset()*0.8);
    h->GetXaxis()->SetTitleOffset(h->GetXaxis()->GetTitleOffset()*0.8);
	h->GetXaxis()->SetTitleSize(h->GetXaxis()->GetTitleSize()*1.8);
   	h->GetXaxis()->SetLabelSize(h->GetXaxis()->GetLabelSize()*1.35);
    h->GetYaxis()->SetTitleSize(h->GetYaxis()->GetTitleSize()*1.8);
    h->GetYaxis()->SetLabelSize(h->GetYaxis()->GetLabelSize()*1.35);
    h->GetZaxis()->SetLabelSize(h->GetZaxis()->GetLabelSize()*1.35);
}

void plot_raw_killres_comparison(){
	const int ndRcuts = 3;
	const int nSigns = 2; //3*2 canvas
	const int nGroups = 2; //raw, excluding-resonances
	const int nKins = 16; //total number of canvases we need
	const int n1Ds = 10;
	const int n2Ds = nKins - n1Ds;

	// std::string kin1Ds[n1Ds] = {"h_pair_dP_overP", "h_Dphi", "h_Deta", "h_DR", "h_Minv", "h_pt_lead", "h_eta_avg", "h_pair_pt", "h_pair_eta", "h_pair_y"};
	// std::string kin2Ds[n1Ds] = {"h_eta_phi", "h_eta1_eta2", "h_eta_avg_Deta", "h_eta_avg_pair_eta", "h_pt1_pt2","h_ptlead_pair_pt"};
	std::string kinNames[nKins] = {"h_pair_dP_overP", "h_Dphi", "h_Deta", "h_DR", "h_Minv", "h_pt_lead", "h_eta_avg", "h_pair_pt", "h_pair_eta", "h_pair_y", "h_eta_phi", "h_eta1_eta2", "h_eta_avg_Deta", "h_eta_avg_pair_eta", "h_pt1_pt2","h_ptlead_pair_pt"};
  	std::string histName[nKins][ndRcuts][nSigns];
  	bool logx[nKins] = {false,false,false,false,false,true,false,true,false,false,false,false,false,false,true,true};
  	bool logy[nKins] = {false,false,false,false,false,true,false,true,false,false,false,false,false,false,true,true};
  	bool logz[n2Ds] = {false,false,false,false,true,true};

	std::string dRs[ndRcuts] = {"_dr1","_dr2","_dr3"};
	std::string signs[nSigns] = {"_sign1", "_sign2"};

	for (unsigned int ikin = 0; ikin < nKins; ikin++){
		for (unsigned int jdr = 0; jdr < ndRcuts; jdr++){
			for (unsigned int ksign = 0; ksign < nSigns; ksign++){
				histName[ikin][jdr][ksign] = kinNames[ikin] + dRs[jdr] + signs[ksign];
			}
		}
	}

	std::string dRTitles[ndRcuts] = {"#Delta R < 0.8","#Delta R < 1.2", "no #Delta R cut"};
	std::string signTitles[nSigns] = {"same sign", "opposite sign"};
	std::string grpTitles[nGroups] = {"raw", "excluding resonances"};

	std::string png_title_2D[nGroups] = {"_raw.png", "_killres.png"};

	std::string subplot_titles_1D[ndRcuts][nSigns][nGroups];
	std::string subplot_titles_2D[ndRcuts][nSigns];
	for (unsigned int jdr = 0; jdr < ndRcuts; jdr++){  	
		for (unsigned int ksign = 0; ksign < nSigns; ksign++){
			
			subplot_titles_2D[jdr][ksign] = dRTitles[jdr] + ", " + signTitles[ksign];
			
			for (unsigned int lres = 0; lres < nGroups; lres++){
				subplot_titles_1D[jdr][ksign][lres] = dRTitles[jdr] + ", " + signTitles[ksign] + ", " + grpTitles[lres];
			}
		}
	}


	// TFile* f_raw = TFile::Open("analysis_raw.root");
	// TFile* f_new = TFile::Open("analysis_kill_resonances.root");
	TFile* f[nGroups];
	f[0] = TFile::Open("analysis_raw.root");
	f[1] = TFile::Open("analysis_kill_resonances.root");

	for (unsigned int ikin = 0; ikin < nKins; ikin++){
		bool is1D = (ikin < n1Ds);

		int numCanvas = (is1D)? 1:nGroups;
		TCanvas *c[numCanvas];
		for (int icanv = 0; icanv < numCanvas; icanv++){
			c[icanv] = new TCanvas();
			c[icanv]->SetWindowSize(1600,1450);
			c[icanv]->Divide(ndRcuts,nSigns);

			for (unsigned int ksign = 0; ksign < nSigns; ksign++){
				for (unsigned int jdr = 0; jdr < ndRcuts; jdr++){
					
					c[icanv]->cd(ksign * ndRcuts + jdr + 1);
					gPad->SetLogx(logx[ikin]);
					gPad->SetLogy(logy[ikin]);
					gPad->SetLeftMargin(0.16);
					gPad->SetBottomMargin(0.135);
	
					if (is1D){
						TH1* h1[nGroups];

						TLegend* l = new TLegend(0.22,0.77,0.46,0.87);
    					l->SetBorderSize(0);
    					l->SetFillStyle(0);
    					l->SetTextFont(42);
    					l->SetMargin(0.2);
    					l->SetTextColor(1);

						for (int igrp = 0; igrp < nGroups; igrp++){
							h1[igrp] = (TH1*) f[igrp]->Get((histName[ikin][jdr][ksign]).c_str());
							hist_helper(h1[igrp],subplot_titles_1D[jdr][ksign][igrp], false,"Events");
    						l->AddEntry(h1[igrp],grpTitles[igrp].c_str(),"lp");
							h1[igrp]->SetMarkerStyle(20);
							h1[igrp]->SetMarkerSize(0.7);
						}

						h1[0]->Draw("E");
						h1[1]->SetLineColor(kRed);
						h1[1]->SetMarkerColor(kRed);
						h1[1]->Draw("E,same");
						l->Draw();
					}else{   //2D
						gPad->SetLogz(logz[ikin-n1Ds]);
						gPad->SetRightMargin(0.135);
	
						// for (int icanv = 0; icanv < numCanvas; icanv++){
						TH2* h2 = (TH2*) f[icanv]->Get((histName[ikin][jdr][ksign]).c_str());
						hist_helper(h2,subplot_titles_2D[jdr][ksign], false);
						h2->Draw("COLZ");
						// }
					}
				}
			}
			// c[icanv]->SaveAs(("plots/nobinning/" + kinNames[ikin] + "_raw_killres.pdf").c_str());
			if(is1D){
				c[icanv]->SaveAs(("plots/nobinning/raw_killres_comparison/" + kinNames[ikin] + "_raw_killres.png").c_str());
			}else{
				c[icanv]->SaveAs(("plots/nobinning/raw_killres_comparison/" + kinNames[ikin] + png_title_2D[icanv]).c_str());
			}
	  		c[icanv]->Close();
	  		delete c[icanv];
	  	}
	}
}