#include "ScrambGenPP.h"
#include <time.h>

void ScrambSampleGenPP::ReadData(){
	
	muon_pt =  new std::vector<float>();
	muon_eta =  new std::vector<float>();
	muon_phi =  new std::vector<float>();
	muon_dP_overP =  new std::vector<float>();
	muon_d0 =  new std::vector<float>();
	muon_z0 =  new std::vector<float>();
	muon_charge =  new std::vector<int>();
	muon_quality =  new std::vector<int>();
	ev_num =  new std::vector<int>();

	Long64_t nentries = inTree->GetEntries();//number of events
	// Long64_t nentries = 10000;
  	for (Long64_t kentry=0; kentry<nentries; kentry++) {//loop over the events
		int num_bytes = inTree->GetEntry(kentry);//read in an event
    	if(num_bytes==0){
      		std::cout<<"Error:: Read in event has size of zero bytes,  quitting"<<std::endl;
      		throw std::exception();
    	}
    	muon_pt->push_back(pt);
    	muon_eta->push_back(eta);
    	muon_phi->push_back(phi);
    	muon_dP_overP->push_back(dP_overP);
    	muon_d0->push_back(d0);
    	muon_z0->push_back(z0);
    	muon_charge->push_back(charge);
    	muon_quality->push_back(quality);
    	ev_num->push_back(event_num);
   }
}


void ScrambSampleGenPP::GenerateRandPair(int num_muon, bool opsign_only){
	// generate a random muon pair
	int nfirst = rand()%num_muon;
	int nsecond = -1;
	do{
		nsecond = rand()%num_muon;
	}while((ev_num->at(nsecond) == ev_num->at(nfirst)) || (opsign_only && muon_charge->at(nsecond) == muon_charge->at(nfirst)));

	mpair = new MuonPair();
	mpair->m1.pt = muon_pt->at(nfirst);
	mpair->m1.eta = muon_eta->at(nfirst);
	mpair->m1.phi = muon_phi->at(nfirst);
	mpair->m1.dP_overP = muon_dP_overP->at(nfirst);
	mpair->m1.d0 = muon_d0->at(nfirst);
	mpair->m1.z0 = muon_z0->at(nfirst);
	mpair->m1.charge = muon_charge->at(nfirst);
	mpair->m1.quality = muon_quality->at(nfirst);
	mpair->m1.ev_num = ev_num->at(nfirst);

	mpair->m2.pt = muon_pt->at(nsecond);
	mpair->m2.eta = muon_eta->at(nsecond);
	mpair->m2.phi = muon_phi->at(nsecond);
	mpair->m2.dP_overP = muon_dP_overP->at(nsecond);
	mpair->m2.d0 = muon_d0->at(nsecond);
	mpair->m2.z0 = muon_z0->at(nsecond);
	mpair->m2.charge = muon_charge->at(nsecond);
	mpair->m2.quality = muon_quality->at(nsecond);
	mpair->m2.ev_num = ev_num->at(nsecond);
}

bool ScrambSampleGenPP::CheckResonance(){ //assumes opposite sign; check if minv falls into any resonance range
	bool isresonance = false;
	for (array<float,2> ires : pms.minv_cuts){
	   if (mpair->minv > ires[0] && mpair->minv < ires[1]) isresonance = true;
	}

	return isresonance;
}

void ScrambSampleGenPP::ImplementOneScramPair(int num_muon, bool opsign_only = false){

	bool isresonance = false;

	do{
		GenerateRandPair(num_muon, opsign_only);
		mpair->Update();

		if (opsign_only){
			assert (!(mpair->same_sign));
			isresonance = CheckResonance();
		}
		else if (!(mpair->same_sign)) isresonance = CheckResonance();

	}while (isresonance); // generate a new pair until we get a non-resonance pair

	if (mpair->same_sign) 	n_ss_scr_pairs++;
	else 					n_op_scr_pairs++;
	
	unsigned int ksign = (mpair->same_sign)? 0:1;

	for (unsigned int idr = 0; idr < ParamsSet::ndRselcs; idr++){
		if (mpair->dr < pms.deltaR_thrsh[idr]){
			outTree[idr][ksign]->Fill();
		}
	}
}


void ScrambSampleGenPP::Run(){
	clock_t start, end;
	double cpu_time_used;

	InitInput();
   	InitOutput();
	ReadData();

	std::cout << "Finished initialization and reading the input trees. Starting to generate scrambled pairs." << std::endl;
	start = clock();

	int nmuon = muon_pt->size();
	n_ss_scr_pairs = 0;
	n_op_scr_pairs = 0;

	do{
		ImplementOneScramPair(nmuon);
		if (n_op_scr_pairs + n_ss_scr_pairs % 1000000 == 0){
			std::cout << "Generating the " << n_op_scr_pairs + n_ss_scr_pairs << "-th scrambled pair" << std::endl;
		}
	}while (n_ss_scr_pairs < nScramb[0]);

	do{
		ImplementOneScramPair(nmuon, true); //only generate opposite-sign pairs now
		if (n_op_scr_pairs + n_ss_scr_pairs % 1000000 == 0){
			std::cout << "Generating the " << n_op_scr_pairs + n_ss_scr_pairs << "-th scrambled pair" << std::endl;
		}
	}while (n_op_scr_pairs < nScramb[1]);

	delete muon_pt;
   	delete muon_eta;
   	delete muon_phi;
   	delete muon_dP_overP;
   	delete muon_d0;
   	delete muon_z0;
   	delete muon_charge;
   	delete muon_quality;
   	delete ev_num;

	outFile->Write();

	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	std::cout << "CPU time used (in seconds) is " << cpu_time_used << std::endl;
}

