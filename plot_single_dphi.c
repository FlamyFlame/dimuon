#include <iostream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <string.h>
#include  <stdlib.h>
#include "TH1.h"
#include "TH2.h"
#include <assert.h>

void plot_single_dphi(unsigned int ctr, unsigned int sign, unsigned int gapcut=1){
	assert (ctr > 0 && ctr <= 6);
	assert (sign == 1 || sign == 2);
	assert (gapcut == 1 || gapcut == 2);

	TFile *f;
	TH2D* h2d;
	if (ctr != 0){
		f = TFile::Open("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/histograms_real_pairs_5_10_20_30_50_80.root");
		h2d = (TH2D*) f->Get(Form("ctr-binned/h_eta_avg_Dphi_dr3_ctr%u_sign%u_gapcut%u",ctr,sign,gapcut));
	}else{ //pp
		f = TFile::Open("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/histograms_real_pairs_pp.root");
		h2d = (TH2D*) f->Get(Form("h_eta_avg_Dphi_dr3_sign%u_gapcut%u",sign,gapcut));
	}

	TH1D* h = nullptr;
	if (sign == 1){
		// h = (TH1D*) h2d->ProjectionX();
		h = (TH1D*) h2d->ProjectionX("",81,120);
	}else{
		// h = (TH1D*) h2d->ProjectionX();
		h = (TH1D*) h2d->ProjectionX("",81,120);
	}
	// h_c1_s1_g1->SetLineColor(kRed);
	// h_c1_s1_g2->SetLineColor(kGreen);
	h->Draw();
}