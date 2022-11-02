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
    TH1D* h = (TH1D*) h2d->ProjectionX();

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
	std::cout << "FILE NAME: " << fname << endl;
	std::cout << "HISTOGRAM NAME: " << hist_name << endl;
	std::cout << "The average bin count in the plateau regions is: " << avg_bin_count << endl;
}


void find_plateau_avg(){
	find_plateau_avg_single_hist("histograms_real_pairs.root","ctr-binned/h_eta_avg_Dphi_dr3_ctr1_sign1_gapcut1",26,46,83,108);
	find_plateau_avg_single_hist("histograms_real_pairs.root","ctr-binned/h_eta_avg_Dphi_dr3_ctr1_sign2_gapcut1",26,52,76,106);
	find_plateau_avg_single_hist("histograms_real_pairs_no_rebinning.root","ctr-binned/h_eta_avg_Dphi_dr3_ctr4_sign1_gapcut1",33,54,77,100);
	find_plateau_avg_single_hist("histograms_real_pairs_no_rebinning.root","ctr-binned/h_eta_avg_Dphi_dr3_ctr4_sign2_gapcut1",24,52,78,103);
	find_plateau_avg_single_hist("histograms_real_pairs_pp.root","h_eta_avg_Dphi_dr3_sign1_gapcut1",35,45,83,90);
	find_plateau_avg_single_hist("histograms_real_pairs_pp.root","h_eta_avg_Dphi_dr3_sign2_gapcut1",37,45,84,94);

	find_plateau_avg_single_hist("histograms_real_pairs.root","ctr-binned/h_eta_avg_Dphi_dr3_ctr1_sign1_gapcut2",26,46,83,108);
	find_plateau_avg_single_hist("histograms_real_pairs.root","ctr-binned/h_eta_avg_Dphi_dr3_ctr1_sign2_gapcut2",26,52,76,106);
	find_plateau_avg_single_hist("histograms_real_pairs_no_rebinning.root","ctr-binned/h_eta_avg_Dphi_dr3_ctr4_sign1_gapcut2",33,54,77,100);
	find_plateau_avg_single_hist("histograms_real_pairs_no_rebinning.root","ctr-binned/h_eta_avg_Dphi_dr3_ctr4_sign2_gapcut2",24,52,78,103);
	find_plateau_avg_single_hist("histograms_real_pairs_pp.root","h_eta_avg_Dphi_dr3_sign1_gapcut2",35,45,83,90);
	find_plateau_avg_single_hist("histograms_real_pairs_pp.root","h_eta_avg_Dphi_dr3_sign2_gapcut2",37,45,84,94);

}