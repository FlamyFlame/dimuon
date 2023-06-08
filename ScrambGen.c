#include "ScrambGen.h"
#include <time.h>


void ScrambSampleGen::GenerateRandPair(int nctr, int mode, bool opsign_only){
	// generate a random muon pair - only satisfies the centrality and pair-sign requirement

	int n1st = -1;
	int n2nd = -1;
	int num_muon = muon_pt[nctr]->size();
	bool same_sign;
	int ctr1;
	int ctr2;
	float avg_ctr; // useful for the left/right spill modes
	float ctr_diff; // useful for the left/right spill modes

	int ctr_bin1, ctr_bin2;

	switch(mode){
		case 0:
			ctr_bin1 = nctr;
			ctr_bin2 = nctr;
			break;
		case 1:
			ctr_bin1 = nctr - 1;
			ctr_bin2 = nctr;
			break;
		case 2:
			ctr_bin1 = nctr;
			ctr_bin2 = nctr + 1;
			break;
	}

	if (mode == 0){ // both muons in the same centrality bin
		
		n1st = rand()%num_muon;
		do{
			n2nd = rand()%num_muon;
			same_sign = (muon_charge[ctr_bin2]->at(n2nd) == muon_charge[ctr_bin1]->at(n1st));
		}while((ev_num[ctr_bin2]->at(n2nd) == ev_num[ctr_bin1]->at(n1st)) || (opsign_only && same_sign));
	
	}else if (mode == 1){ // left spill: the second muon is in the curent centrality bin; the first muon is in the left adjacent bin
		if (nctr < 1){
			std::cout << "For the left spill mode, the centrality bin CANNOT be the first (0)!" << std::endl;
			throw std::exception();
		}

		// cout << "ctr bin 1st: " << ctr_bin1 << ", ctr bin 2nd: " << ctr_bin2 << endl;
		int num_muon_left = muon_pt[nctr - 1]->size();
		// cout << "num muon left" << num_muon_left << ", num muon" << num_muon << endl;

		do {
			n1st = rand()%num_muon_left;
			// cout << "trial n1st: " << n1st << endl;
			ctr1 = ev_centrality[ctr_bin1]->at(n1st);
		}while(ctr1 <= ctr_step * (nctr - 0.5));

		// cout << "final choice: n1st " << n1st << ", ctr1 " << ctr1 << endl;
		
		do{
			n2nd = rand()%num_muon;
			// cout << "trial n2nd: " << n2nd << endl;
			same_sign = (muon_charge[ctr_bin2]->at(n2nd) == muon_charge[ctr_bin1]->at(n1st));
			ctr2 = ev_centrality[ctr_bin2]->at(n2nd);
			avg_ctr = (ctr1 + ctr2) / 2.;
			ctr_diff = ctr2 - ctr1;
		}while((avg_ctr < ctr_step * nctr || ctr_diff > ctr_step) || (opsign_only && same_sign));

		// cout << "final choice: n2nd " << n2nd << ", ctr2 " << ctr2 << endl;
	
	}else{ // mode == 2, right spill: the first muon is in the curent centrality bin; the second muon is in the right adjacent bin
		if (nctr > nctr_intvls - 2){
			std::cout << "For the right spill mode, the centrality bin CANNOT be the last (nctr_intvls-1)!" << std::endl;
			throw std::exception();
		}

		int num_muon_right = muon_pt[nctr + 1]->size();

		do{
			n2nd = rand()%num_muon_right;
			ctr2 = ev_centrality[ctr_bin2]->at(n2nd);
		}while(ctr2 >= ctr_step * (nctr + 1.5));

		do{
			n1st = rand()%num_muon;
			same_sign = (muon_charge[ctr_bin2]->at(n2nd) == muon_charge[ctr_bin1]->at(n1st));
			ctr1 = ev_centrality[ctr_bin1]->at(n1st);
			avg_ctr = (ctr2 + ctr1) / 2.;
			ctr_diff = ctr2 - ctr1;
		}while((avg_ctr > ctr_step * (nctr+1) || ctr_diff > ctr_step) || (opsign_only && same_sign));
	}

	mpair->m1.pt = muon_pt[ctr_bin1]->at(n1st);
	mpair->m1.eta = muon_eta[ctr_bin1]->at(n1st);
	mpair->m1.phi = muon_phi[ctr_bin1]->at(n1st);
	mpair->m1.dP_overP = muon_dP_overP[ctr_bin1]->at(n1st);
	mpair->m1.d0 = muon_d0[ctr_bin1]->at(n1st);
	mpair->m1.z0 = muon_z0[ctr_bin1]->at(n1st);
	mpair->m1.charge = muon_charge[ctr_bin1]->at(n1st);
	mpair->m1.quality = muon_quality[ctr_bin1]->at(n1st);
	mpair->m1.ev_num = ev_num[ctr_bin1]->at(n1st);
	mpair->m1.ev_centrality = ev_centrality[ctr_bin1]->at(n1st);
	mpair->m1.ev_FCal_Et = ev_FCal_Et[ctr_bin1]->at(n1st);

	mpair->m2.pt = muon_pt[ctr_bin2]->at(n2nd);
	mpair->m2.eta = muon_eta[ctr_bin2]->at(n2nd);
	mpair->m2.phi = muon_phi[ctr_bin2]->at(n2nd);
	mpair->m2.dP_overP = muon_dP_overP[ctr_bin2]->at(n2nd);
	mpair->m2.d0 = muon_d0[ctr_bin2]->at(n2nd);
	mpair->m2.z0 = muon_z0[ctr_bin2]->at(n2nd);
	mpair->m2.charge = muon_charge[ctr_bin2]->at(n2nd);
	mpair->m2.quality = muon_quality[ctr_bin2]->at(n2nd);
	mpair->m2.ev_num = ev_num[ctr_bin2]->at(n2nd);
	mpair->m2.ev_centrality = ev_centrality[ctr_bin2]->at(n2nd);
	mpair->m2.ev_FCal_Et = ev_FCal_Et[ctr_bin2]->at(n2nd);
}

bool ScrambSampleGen::ResonanceCut(){ //assumes opposite sign; check if minv falls into any resonance range

	// IMPORTANT: here we are assuming we have already filled up the muon pair (mpair)
	if (mpair->same_sign) return false;
	
	if (mpair->minv > pms.minv_upper){ // upper cut at 60 GeV
		return true;
	}
	
	for (array<float,2> ires : pms.minv_cuts){
		if (mpair->minv > ires[0] && mpair->minv < ires[1]){
			return true;
		}
	}

	return false;

}

bool ScrambSampleGen::PhotoProductionCut(){ 
	return (!(mpair->same_sign) && mpair->asym < 0.05 && mpair->acop < 0.01);
}

void ScrambSampleGen::ImplementOneScramPair(int nctr, int mode, bool opsign_only = false){

	// call the function GenerateRandPair until a pair passing both resonance and photo-production cuts is generated
	// then update the pair and fill in the corresponding TTrees

	// std::cout << "Start implementing one scrambled pair for centrality bin " << nctr << ", mode " << mode << std::endl;
	bool resonance_cut = false;
	bool photoprod_cut = false;

	do{
		GenerateRandPair(nctr, mode, opsign_only);
		mpair->Update(); // IMPORTANT: first update then cut if fail the resonance and photoproduction cuts

		resonance_cut = ResonanceCut();
		photoprod_cut = PhotoProductionCut();

	}while (resonance_cut || photoprod_cut);

	if (mpair->same_sign) 	n_ss_scr_pairs++;
	else 					n_op_scr_pairs++;
	
	unsigned int ksign = (mpair->same_sign)? 0:1;

	// for (unsigned int idr = 0; idr < ParamsSet::ndRselcs; idr++){
	for (unsigned int idr = ParamsSet::ndRselcs-1; idr < ParamsSet::ndRselcs; idr++){
		if (mpair->dr < pms.deltaR_thrsh[idr]){
			// FillHistograms(idr, jctr, ksign);
			outTree[idr][nctr][ksign]->Fill();
		}
	}
}


void ScrambSampleGen::Run(){
	clock_t start, end;
	double cpu_time_used;

	// start = clock();
	InitInput();
   	InitOutput();
	InitTempVariables();
	// end = clock();
	// cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	// std::cout << "CPU time used (in seconds) for initializaing and reading the input trees is " << cpu_time_used << std::endl;

	std::cout << "Finished initialization and reading the input trees. Starting to generate scrambled pairs." << std::endl;
	start = clock();

	mpair = new MuonPair();

	nmodes = (keep_spill)? 3 : 1;

	for (int jctr = 0; jctr < nctr_intvls; jctr++){

		int ctr_start;

		switch (jctr){
		case 0:
			ctr_start = 0;
			break;
		case nctr_intvls - 1:
			ctr_start = nctr_intvls + 1;
			break;
		default:
			ctr_start = jctr + 1;
		}
		cout << "ctr intvl: " << jctr << ", jjctr: ";
		for (int jjctr = ctr_start; jjctr <= jctr + 1; jjctr++){
			cout << jjctr << " ";
	        Long64_t nentries = inTree[jjctr]->GetEntries(); //#muon pairs
			// Long64_t nentries = 10000;
		  	for (Long64_t kentry=0; kentry<nentries; kentry++) {//loop over the events
				int num_bytes = inTree[jjctr]->GetEntry(kentry);//read in an event
		    	if(num_bytes==0){
		      		std::cout<<"Error:: Read in event has size of zero bytes,  quitting"<<std::endl;
		      		throw std::exception();
		    	}
		    	muon_pt[jjctr]->push_back(pt[jjctr]);
		    	muon_eta[jjctr]->push_back(eta[jjctr]);
		    	muon_phi[jjctr]->push_back(phi[jjctr]);
		    	muon_dP_overP[jjctr]->push_back(dP_overP[jjctr]);
		    	muon_d0[jjctr]->push_back(d0[jjctr]);
		    	muon_z0[jjctr]->push_back(z0[jjctr]);
		    	muon_charge[jjctr]->push_back(charge[jjctr]);
		    	muon_quality[jjctr]->push_back(quality[jjctr]);
		    	ev_num[jjctr]->push_back(event_num[jjctr]);
		    	ev_centrality[jjctr]->push_back(centrality[jjctr]);
			   	ev_FCal_Et[jjctr]->push_back(FCal_Et[jjctr]);
		   	}
		}
		// cout << endl;

		cout << "centrality " << jctr << ", #single muons " << muon_pt[jctr]->size() << endl;
		for (int mode = 0; mode < nmodes; mode++){ // perform all 3 modes if keep spill; only mode 0 if not

			if (jctr == 0 && mode == 1) continue; // the first bin cannot have left spill
			if (jctr == nctr_intvls - 1 && mode == 2) continue; // the last bin cannot have right spill
			n_ss_scr_pairs = 0; // restart the count for the current centrality bin & mode
			n_op_scr_pairs = 0;

			cout << "Current mode: " << mode << endl;
			do{
				ImplementOneScramPair(jctr, mode);
				if (n_op_scr_pairs + n_ss_scr_pairs % 1000000 == 0){
					std::cout << "Generating the " << n_op_scr_pairs + n_ss_scr_pairs << "-th scrambled pair, for centrality bin " << jctr+1 << std::endl;
				}
			}while (n_ss_scr_pairs < nScramb_ss[mode][jctr]);

			do{
				ImplementOneScramPair(jctr, mode, true);
				if (n_op_scr_pairs + n_ss_scr_pairs % 1000000 == 0){
					std::cout << "Generating the " << n_op_scr_pairs + n_ss_scr_pairs << "-th scrambled pair, for centrality bin " << jctr+1 << std::endl;
				}
			}while (n_op_scr_pairs < nScramb_op[mode][jctr]);

		}
		
		cout << "#scrambled pairs [same sign, current centrality bin] should be: " << nScramb_ss[0][jctr] + nScramb_ss[1][jctr] + nScramb_ss[2][jctr] << std::endl;
		cout << "the actual #scrambled pairs [same sign, current centrality bin] is: " << outTree[2][jctr][0]->GetEntries() << std::endl;
		cout << "#scrambled pairs [opposite sign, current centrality bin] should be: " << nScramb_op[0][jctr] + nScramb_op[1][jctr] + nScramb_op[2][jctr] << std::endl;
		cout << "the actual #scrambled pairs [opposite sign, current centrality bin] is: " << outTree[2][jctr][1]->GetEntries() << std::endl;

		if (jctr >= 1){ 
		// delete muon info in not the current, but the previous centrality bin
		// so that the next centrality bin can perform left-spill correctly
			delete muon_pt[jctr-1];
	   		delete muon_eta[jctr-1];
	   		delete muon_phi[jctr-1];
	   		delete muon_dP_overP[jctr-1];
	   		delete muon_d0[jctr-1];
	   		delete muon_z0[jctr-1];
	   		delete muon_charge[jctr-1];
	   		delete muon_quality[jctr-1];
	   		delete ev_num[jctr-1];
	   		delete ev_centrality[jctr-1];
	   		delete ev_FCal_Et[jctr-1];
		}
	}

	cout << "Have finished the for loop over the centrality bins.";
	cout << "Ready to write to the output file" << endl;
	outFile->Write();
	delete mpair;

	// for (unsigned int idr = 0; idr < ParamsSet::ndRselcs; idr++){
 //      for (unsigned int jctr = 0; jctr < nctr_intvls; jctr++){
 //         for (unsigned int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
 //            outTree[idr][jctr][ksign]->Write();
 //         }
 //      }
 //   }

	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	std::cout << "CPU time used (in seconds) is " << cpu_time_used << std::endl;
}

