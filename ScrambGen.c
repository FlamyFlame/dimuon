#include "ScrambGen.h"
#include <time.h>

//#define ParamsSet::ndRselcs 3

void ScrambSampleGen::ReadData(){
	
	//for (int i = 0; i < ParamsSet::ndRselcs; i++){
	for (int jctr = 0; jctr < ParamsSet::nCtrBins; jctr++){
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

		Long64_t nentries = inTree[jctr]->GetEntries();//number of events
		// Long64_t nentries = 10000;
  		for (Long64_t kentry=0; kentry<nentries; kentry++) {//loop over the events
			int num_bytes = inTree[jctr]->GetEntry(kentry);//read in an event
    		if(num_bytes==0){
    	  		std::cout<<"Error:: Read in event has size of zero bytes,  quitting"<<std::endl;
    	  		throw std::exception();
    		}
    		muon_pt[jctr]->push_back(pt[jctr]);
    		muon_eta[jctr]->push_back(eta[jctr]);
    		muon_phi[jctr]->push_back(phi[jctr]);
    		muon_dP_overP[jctr]->push_back(dP_overP[jctr]);
    		muon_d0[jctr]->push_back(d0[jctr]);
    		muon_z0[jctr]->push_back(z0[jctr]);
    		muon_charge[jctr]->push_back(charge[jctr]);
    		muon_quality[jctr]->push_back(quality[jctr]);
    		ev_num[jctr]->push_back(event_num[jctr]);
    		ev_centrality[jctr]->push_back(centrality[jctr]);
		   ev_FCal_Et[jctr]->push_back(FCal_Et[jctr]);
   	}
   }
}


// void ScrambSampleGen::FillHistograms(unsigned int ndr, unsigned int nctr, unsigned int nsign){
//   assert(ndr < ParamsSet::ndRselcs && nctr < ParamsSet::nCtrBins && nsign < ParamsSet::nSigns);
//   h_pair_dP_overP[ndr][nctr][nsign]->Fill(mpair->pair_dPoverP);
//   h_Minv[ndr][nctr][nsign]   ->Fill(mpair->minv);
//   h_Dphi[ndr][nctr][nsign]   ->Fill(mpair->dphi);
//   h_Deta[ndr][nctr][nsign]   ->Fill(mpair->deta);
//   h_DR[ndr][nctr][nsign]     ->Fill(mpair->dr);
//   h_pt_lead[ndr][nctr][nsign] ->Fill(mpair->m1.pt);
//   h_eta_avg[ndr][nctr][nsign] ->Fill(mpair->etaavg);
//   h_pair_pt[ndr][nctr][nsign]->Fill(mpair->pair_pt);
//   h_pair_eta[ndr][nctr][nsign]->Fill(mpair->pair_eta);
//   h_pair_y[ndr][nctr][nsign]->Fill(mpair->pair_y);
//   h_eta_phi[ndr][nctr][nsign]->Fill(mpair->phiavg,mpair->etaavg);
//   h_eta1_eta2[ndr][nctr][nsign]->Fill(mpair->m2.eta,mpair->m1.eta);
//   h_pt1_pt2[ndr][nctr][nsign]->Fill(mpair->m2.pt,mpair->m1.pt);
//   h_eta_avg_Deta[ndr][nctr][nsign]->Fill(mpair->deta,mpair->etaavg);        
//   h_eta_avg_pair_eta[ndr][nctr][nsign]->Fill(mpair->pair_eta,mpair->etaavg);
//   h_ptlead_pair_pt[ndr][nctr][nsign]->Fill(mpair->pair_pt,mpair->m1.pt);
// }


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

bool ScrambSampleGen::CheckResonance(){ //assumes opposite sign; check if minv falls into any resonance range
    if (mpair->minv > pms.minv_upper) continue; // upper cut at 80 GeV
	
	bool isresonance = false;
	for (array<float,2> ires : pms.minv_cuts){
	   if (mpair->minv > ires[0] && mpair->minv < ires[1]) isresonance = true;
	}

	return isresonance;
}

void ScrambSampleGen::ImplementOneScramPair(int num_muon, int nctr, bool opsign_only = false){

	bool isresonance = false;

	do{
		GenerateRandPair(num_muon,nctr,opsign_only);
		mpair->Update();

		if (opsign_only){
			assert (!(mpair->same_sign));
			isresonance = CheckResonance();
		}
		else if (!(mpair->same_sign)) isresonance = CheckResonance();

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

	for (int jctr = 0; jctr < ParamsSet::nCtrBins; jctr++){
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
 //      for (unsigned int jctr = 0; jctr < ParamsSet::nCtrBins; jctr++){
 //         for (unsigned int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
 //            outTree[idr][jctr][ksign]->Write();
 //         }
 //      }
 //   }

	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	std::cout << "CPU time used (in seconds) is " << cpu_time_used << std::endl;
}

