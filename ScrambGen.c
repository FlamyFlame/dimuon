#include "ScrambGen.h"
#include <time.h>

//#define ParamsSet::ndRselcs 3

void ScrambSampleGen::ReadData(){
	
	int cur_ctr_intvl = 0;
	//for (int i = 0; i < ParamsSet::ndRselcs; i++){
	for (int jctr = 0; jctr < nCtrBins; jctr++){
		muon_pt[jctr] =  new std::vector<float>();
		muon_eta[jctr] =  new std::vector<float>();
		muon_phi[jctr] =  new std::vector<float>();
		muon_dP_overP[jctr] =  new std::vector<float>();
		muon_d0[jctr] =  new std::vector<float>();
		muon_z0[jctr] =  new std::vector<float>();
		muon_charge[jctr] =  new std::vector<int>();
		muon_quality[jctr] =  new std::vector<int>();
		ev_num[jctr] =  new std::vector<int>();
		ev_centrality[jctr] =  new std::vector<int>();
		ev_FCal_Et[jctr] =  new std::vector<float>();

		for (unsigned int jj = 0; jj < ctrBins[jctr].size(); jj++){
        	Long64_t nentries = inTree[cur_ctr_intvl]->GetEntries(); //#muon pairs
			// Long64_t nentries = 10000;
	  		for (Long64_t kentry=0; kentry<nentries; kentry++) {//loop over the events
				int num_bytes = inTree[cur_ctr_intvl]->GetEntry(kentry);//read in an event
	    		if(num_bytes==0){
	    	  		std::cout<<"Error:: Read in event has size of zero bytes,  quitting"<<std::endl;
	    	  		throw std::exception();
	    		}
	    		muon_pt[jctr]->push_back(pt[cur_ctr_intvl]);
	    		muon_eta[jctr]->push_back(eta[cur_ctr_intvl]);
	    		muon_phi[jctr]->push_back(phi[cur_ctr_intvl]);
	    		muon_dP_overP[jctr]->push_back(dP_overP[cur_ctr_intvl]);
	    		muon_d0[jctr]->push_back(d0[cur_ctr_intvl]);
	    		muon_z0[jctr]->push_back(z0[cur_ctr_intvl]);
	    		muon_charge[jctr]->push_back(charge[cur_ctr_intvl]);
	    		muon_quality[jctr]->push_back(quality[cur_ctr_intvl]);
	    		ev_num[jctr]->push_back(event_num[cur_ctr_intvl]);
	    		ev_centrality[jctr]->push_back(centrality[cur_ctr_intvl]);
			   ev_FCal_Et[jctr]->push_back(FCal_Et[cur_ctr_intvl]);
	   		}
	   		cur_ctr_intvl++;
	   	}
    }
}


void ScrambSampleGen::GenerateRandPair(int num_muon, int nctr, bool opsign_only){
	// generate a random muon pair
	int nfirst = rand()%num_muon;
	int nsecond;
	do{
		nsecond = rand()%num_muon;
	}while((ev_num[nctr]->at(nsecond) == ev_num[nctr]->at(nfirst)) || (fabs(ev_centrality[nctr]->at(nsecond) - ev_centrality[nctr]->at(nfirst)) >= 10) || (opsign_only && muon_charge[nctr]->at(nsecond) == muon_charge[nctr]->at(nfirst)));

	mpair = new MuonPair();
	mpair->m1.pt = muon_pt[nctr]->at(nfirst);
	mpair->m1.eta = muon_eta[nctr]->at(nfirst);
	mpair->m1.phi = muon_phi[nctr]->at(nfirst);
	mpair->m1.dP_overP = muon_dP_overP[nctr]->at(nfirst);
	mpair->m1.d0 = muon_d0[nctr]->at(nfirst);
	mpair->m1.z0 = muon_z0[nctr]->at(nfirst);
	mpair->m1.charge = muon_charge[nctr]->at(nfirst);
	mpair->m1.quality = muon_quality[nctr]->at(nfirst);
	mpair->m1.ev_num = ev_num[nctr]->at(nfirst);
	mpair->m1.ev_centrality = ev_centrality[nctr]->at(nfirst);
	mpair->m1.ev_FCal_Et = ev_FCal_Et[nctr]->at(nfirst);

	mpair->m2.pt = muon_pt[nctr]->at(nsecond);
	mpair->m2.eta = muon_eta[nctr]->at(nsecond);
	mpair->m2.phi = muon_phi[nctr]->at(nsecond);
	mpair->m2.dP_overP = muon_dP_overP[nctr]->at(nsecond);
	mpair->m2.d0 = muon_d0[nctr]->at(nsecond);
	mpair->m2.z0 = muon_z0[nctr]->at(nsecond);
	mpair->m2.charge = muon_charge[nctr]->at(nsecond);
	mpair->m2.quality = muon_quality[nctr]->at(nsecond);
	mpair->m2.ev_num = ev_num[nctr]->at(nsecond);
	mpair->m2.ev_centrality = ev_centrality[nctr]->at(nsecond);
	mpair->m2.ev_FCal_Et = ev_FCal_Et[nctr]->at(nsecond);
}

bool ScrambSampleGen::IsResonance(){ //assumes opposite sign; check if minv falls into any resonance range
    if (mpair->minv > pms.minv_upper) continue; // upper cut at 80 GeV
	
	bool isresonance = false;
	for (array<float,2> ires : pms.minv_cuts){
	   if (mpair->minv > ires[0] && mpair->minv < ires[1]) isresonance = true;
	}

	return isresonance;
}

bool ScrambSampleGen::IsPhotoProduction(){ 

void ScrambSampleGen::ImplementOneScramPair(int num_muon, int nctr, bool opsign_only = false){

	bool isresonance = false;

	do{
		GenerateRandPair(num_muon,nctr,opsign_only);
		mpair->Update();

		if (opsign_only){
			assert (!(mpair->same_sign));
			isresonance = IsResonance();
		}
		else if (!(mpair->same_sign)) isresonance = IsResonance();

	}while (isresonance);

	if (mpair->same_sign) 	n_ss_scr_pairs++;
	else 					n_op_scr_pairs++;
	
	unsigned int ksign = (mpair->same_sign)? 0:1;

	for (unsigned int idr = 0; idr < ParamsSet::ndRselcs; idr++){
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
	ReadData();
	// end = clock();
	// cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	// std::cout << "CPU time used (in seconds) for initializaing and reading the input trees is " << cpu_time_used << std::endl;

	std::cout << "Finished initialization and reading the input trees. Starting to generate scrambled pairs." << std::endl;
	start = clock();

	for (int jctr = 0; jctr < nCtrBins; jctr++){
		int nmuon = muon_pt[jctr]->size(); //#muons in the centrality bin
		n_ss_scr_pairs = 0; // restart the count for the current centrality bin
		n_op_scr_pairs = 0;

		do{
			ImplementOneScramPair(nmuon, jctr);
			if (n_op_scr_pairs + n_ss_scr_pairs % 1000000 == 0){
				std::cout << "Generating the " << n_op_scr_pairs + n_ss_scr_pairs << "-th scrambled pair, for centrality bin " << jctr+1 << std::endl;
			}
		}while (n_ss_scr_pairs < nScramb[jctr][0]);

		do{
			ImplementOneScramPair(nmuon, jctr, true);
			if (n_op_scr_pairs + n_ss_scr_pairs % 1000000 == 0){
				std::cout << "Generating the " << n_op_scr_pairs + n_ss_scr_pairs << "-th scrambled pair, for centrality bin " << jctr+1 << std::endl;
			}
		}while (n_op_scr_pairs < nScramb[jctr][1]);

		delete muon_pt[jctr];
   		delete muon_eta[jctr];
   		delete muon_phi[jctr];
   		delete muon_dP_overP[jctr];
   		delete muon_d0[jctr];
   		delete muon_z0[jctr];
   		delete muon_charge[jctr];
   		delete muon_quality[jctr];
   		delete ev_num[jctr];
   		delete ev_centrality[jctr];
   		delete ev_FCal_Et[jctr];
	}

	outFile->Write();

	// for (unsigned int idr = 0; idr < ParamsSet::ndRselcs; idr++){
 //      for (unsigned int jctr = 0; jctr < nCtrBins; jctr++){
 //         for (unsigned int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
 //            outTree[idr][jctr][ksign]->Write();
 //         }
 //      }
 //   }

	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	std::cout << "CPU time used (in seconds) is " << cpu_time_used << std::endl;
}

