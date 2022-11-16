#include <iostream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <string.h>
#include  <stdlib.h>
#include "TH1.h"
#include "TH2.h"
#include <assert.h> 

using namespace std;

void find_plateau_avg_single_hist(string fname, string hist_name, int a1, int a2, int a3, int a4){
	// fname: name of the root file containing histogram of interest
	// hist_name: name of histogram of interest: should start with h_eta_avg_Dphi, and without dR cut
	// a1, a2: first and last bin in the left (dphi < 0) plateau region
	// a3, a4: first and last bin in the right (dphi > 0) plateau region

	int nbins_plateau = 0; // total number of bins in the plateau region
	int total_counts_plateau = 0; // total number of bin counts in the plateau region

	std::string path = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/";
	TFile* f = TFile::Open((path + fname).c_str());
    TH2D* h2d = (TH2D*) f->Get(hist_name.c_str());
    TH1D* h = (TH1D*) h2d->ProjectionX("",81,120);

	for (int i = a1; i <= a2; i++){
		nbins_plateau += 1;
		total_counts_plateau += h->GetBinContent(i);
	}

	for (int i = a3; i <= a4; i++){
		nbins_plateau += 1;
		total_counts_plateau += h->GetBinContent(i);
	}

	assert(nbins_plateau > 0);
	double avg_bin_count = float(total_counts_plateau) / nbins_plateau;
	int floor = static_cast<int>(avg_bin_count);
	// std::cout << "FILE NAME: " << fname << endl;
	// std::cout << "HISTOGRAM NAME: " << hist_name << endl;
	// std::cout << "The average bin count in the plateau regions is: " << avg_bin_count << endl;
	if (avg_bin_count - floor >= 0.5)
		std::cout << floor + 1 << ", ";
	else
		std::cout << floor << ", ";
}


void find_plateau_avg(){
	find_plateau_avg_single_hist("histograms_real_pairs_5_10_20_30_50_80.root","ctr-binned/h_eta_avg_Dphi_dr3_ctr1_sign1_gapcut1",18,56,80,108);
	find_plateau_avg_single_hist("histograms_real_pairs_5_10_20_30_50_80.root","ctr-binned/h_eta_avg_Dphi_dr3_ctr1_sign2_gapcut1",19,51,76,106);
	find_plateau_avg_single_hist("histograms_real_pairs_5_10_20_30_50_80.root","ctr-binned/h_eta_avg_Dphi_dr3_ctr2_sign1_gapcut1",30,45,85,108);
	find_plateau_avg_single_hist("histograms_real_pairs_5_10_20_30_50_80.root","ctr-binned/h_eta_avg_Dphi_dr3_ctr2_sign2_gapcut1",24,51,81,109);
	find_plateau_avg_single_hist("histograms_real_pairs_5_10_20_30_50_80.root","ctr-binned/h_eta_avg_Dphi_dr3_ctr3_sign1_gapcut1",31,50,79,105);
	find_plateau_avg_single_hist("histograms_real_pairs_5_10_20_30_50_80.root","ctr-binned/h_eta_avg_Dphi_dr3_ctr3_sign2_gapcut1",26,51,77,106);
	find_plateau_avg_single_hist("histograms_real_pairs_5_10_20_30_50_80.root","ctr-binned/h_eta_avg_Dphi_dr3_ctr4_sign1_gapcut1",30,51,78,107);
	find_plateau_avg_single_hist("histograms_real_pairs_5_10_20_30_50_80.root","ctr-binned/h_eta_avg_Dphi_dr3_ctr4_sign2_gapcut1",23,52,78,107);
	find_plateau_avg_single_hist("histograms_real_pairs_5_10_20_30_50_80.root","ctr-binned/h_eta_avg_Dphi_dr3_ctr5_sign1_gapcut1",27,51,81,100);
	find_plateau_avg_single_hist("histograms_real_pairs_5_10_20_30_50_80.root","ctr-binned/h_eta_avg_Dphi_dr3_ctr5_sign2_gapcut1",27,52,80,104);
	find_plateau_avg_single_hist("histograms_real_pairs_5_10_20_30_50_80.root","ctr-binned/h_eta_avg_Dphi_dr3_ctr6_sign1_gapcut1",27,53,81,102);
	find_plateau_avg_single_hist("histograms_real_pairs_5_10_20_30_50_80.root","ctr-binned/h_eta_avg_Dphi_dr3_ctr6_sign2_gapcut1",27,51,80,102);
	find_plateau_avg_single_hist("histograms_real_pairs_pp.root","h_eta_avg_Dphi_dr3_sign1_gapcut1",42,42,87,90);
	find_plateau_avg_single_hist("histograms_real_pairs_pp.root","h_eta_avg_Dphi_dr3_sign2_gapcut1",37,45,88,92);

	std::cout << std::endl;
	find_plateau_avg_single_hist("histograms_real_pairs_5_10_20_30_50_80.root","ctr-binned/h_eta_avg_Dphi_dr3_ctr1_sign1_gapcut2",18,56,80,108);
	find_plateau_avg_single_hist("histograms_real_pairs_5_10_20_30_50_80.root","ctr-binned/h_eta_avg_Dphi_dr3_ctr1_sign2_gapcut2",19,51,76,106);
	find_plateau_avg_single_hist("histograms_real_pairs_5_10_20_30_50_80.root","ctr-binned/h_eta_avg_Dphi_dr3_ctr2_sign1_gapcut2",30,45,85,108);
	find_plateau_avg_single_hist("histograms_real_pairs_5_10_20_30_50_80.root","ctr-binned/h_eta_avg_Dphi_dr3_ctr2_sign2_gapcut2",24,51,81,109);
	find_plateau_avg_single_hist("histograms_real_pairs_5_10_20_30_50_80.root","ctr-binned/h_eta_avg_Dphi_dr3_ctr3_sign1_gapcut2",31,50,79,105);
	find_plateau_avg_single_hist("histograms_real_pairs_5_10_20_30_50_80.root","ctr-binned/h_eta_avg_Dphi_dr3_ctr3_sign2_gapcut2",26,51,77,106);
	find_plateau_avg_single_hist("histograms_real_pairs_5_10_20_30_50_80.root","ctr-binned/h_eta_avg_Dphi_dr3_ctr4_sign1_gapcut2",30,51,78,107);
	find_plateau_avg_single_hist("histograms_real_pairs_5_10_20_30_50_80.root","ctr-binned/h_eta_avg_Dphi_dr3_ctr4_sign2_gapcut2",23,52,78,107);
	find_plateau_avg_single_hist("histograms_real_pairs_5_10_20_30_50_80.root","ctr-binned/h_eta_avg_Dphi_dr3_ctr5_sign1_gapcut2",27,51,81,100);
	find_plateau_avg_single_hist("histograms_real_pairs_5_10_20_30_50_80.root","ctr-binned/h_eta_avg_Dphi_dr3_ctr5_sign2_gapcut2",27,52,80,104);
	find_plateau_avg_single_hist("histograms_real_pairs_5_10_20_30_50_80.root","ctr-binned/h_eta_avg_Dphi_dr3_ctr6_sign1_gapcut2",27,53,81,102);
	find_plateau_avg_single_hist("histograms_real_pairs_5_10_20_30_50_80.root","ctr-binned/h_eta_avg_Dphi_dr3_ctr6_sign2_gapcut2",27,51,80,102);
	find_plateau_avg_single_hist("histograms_real_pairs_pp.root","h_eta_avg_Dphi_dr3_sign1_gapcut2",40,45,87,90);
	find_plateau_avg_single_hist("histograms_real_pairs_pp.root","h_eta_avg_Dphi_dr3_sign2_gapcut2",37,45,88,92);

}