#include "DimuonDataAnalysisBaseClass.h"

void DimuonDataAnalysisBaseClass::PrintInstructions(){
}

void DimuonDataAnalysisBaseClass::InitInput(){
	TChainFill();
	if (!fChain){
		std::cerr << "FATAL:: TChain for analysis is nullptr! Throwing exception." << std::endl;
    	throw std::exception();
	}

    fChain->SetBranchAddress("RunNumber"                   , &RunNumber);
    fChain->SetBranchAddress("lbn"                         , &lbn);
    fChain->SetBranchAddress("bcid"                        , &bcid);
    
    fChain->SetBranchAddress("muon_deltaP_overP"           , &muon_deltaP_overP);
  
    fChain->SetBranchAddress("muon_pair_muon1_index"       , &muon_pair_muon1_index);
    fChain->SetBranchAddress("muon_pair_muon1_pt"          , &muon_pair_muon1_pt);
    fChain->SetBranchAddress("muon_pair_muon1_eta"         , &muon_pair_muon1_eta);
    fChain->SetBranchAddress("muon_pair_muon1_phi"         , &muon_pair_muon1_phi);
    fChain->SetBranchAddress("muon_pair_muon1_quality"     , &muon_pair_muon1_quality);
    fChain->SetBranchAddress("muon_pair_muon1_d0"          , &muon_pair_muon1_d0);
    fChain->SetBranchAddress("muon_pair_muon1_z0"          , &muon_pair_muon1_z0);
    fChain->SetBranchAddress("muon_pair_muon1_trk_pt"      , &muon_pair_muon1_trk_pt);
    fChain->SetBranchAddress("muon_pair_muon1_trk_eta"     , &muon_pair_muon1_trk_eta);
    fChain->SetBranchAddress("muon_pair_muon1_trk_phi"     , &muon_pair_muon1_trk_phi);
  
    fChain->SetBranchAddress("muon_pair_muon2_index"       , &muon_pair_muon2_index);
    fChain->SetBranchAddress("muon_pair_muon2_pt"          , &muon_pair_muon2_pt);
    fChain->SetBranchAddress("muon_pair_muon2_eta"         , &muon_pair_muon2_eta);
    fChain->SetBranchAddress("muon_pair_muon2_phi"         , &muon_pair_muon2_phi);
    fChain->SetBranchAddress("muon_pair_muon2_quality"     , &muon_pair_muon2_quality);
    fChain->SetBranchAddress("muon_pair_muon2_d0"          , &muon_pair_muon2_d0);
    fChain->SetBranchAddress("muon_pair_muon2_z0"          , &muon_pair_muon2_z0);
    fChain->SetBranchAddress("muon_pair_muon2_trk_pt"      , &muon_pair_muon2_trk_pt);
    fChain->SetBranchAddress("muon_pair_muon2_trk_eta"     , &muon_pair_muon2_trk_eta);
    fChain->SetBranchAddress("muon_pair_muon2_trk_phi"     , &muon_pair_muon2_trk_phi);
    
    std::string mu4_trigger_name = (isRun3)? "HLT_mu4_L1MU3V" : "HLT_mu4";
    std::string mu4_trigger_branch = "b_" + mu4_trigger_name;
    std::string mu4_trigger_match_branch = "muon_b_" + mu4_trigger_name;

    std::string mu4_mu4noL1_trigger_name = (isRun3)? "HLT_mu4_mu4noL1_L1MU3V" : "HLT_mu4_mu4noL1";
    std::string mu4_mu4noL1_trigger_branch = "b_" + mu4_mu4noL1_trigger_name;
    std::string mu4_mu4noL1_trigger_match_branch = "dimuon_b_" + mu4_mu4noL1_trigger_name;

    std::string twomu4_trigger_name = (isRun3)? "HLT_2mu4_L12MU3V" : "HLT_2mu4";
    std::string twomu4_trigger_branch = "b_" + twomu4_trigger_name;
    std::string twomu4_trigger_match_branch = "dimuon_b_" + twomu4_trigger_name;

    fChain->SetBranchAddress(mu4_trigger_branch.c_str()                       , &b_HLT_mu4);
    fChain->SetBranchAddress(twomu4_trigger_branch.c_str()                    , &b_HLT_2mu4);
    fChain->SetBranchAddress(mu4_trigger_match_branch.c_str()                 , &muon_b_HLT_mu4);
    fChain->SetBranchAddress(twomu4_trigger_match_branch.c_str()              , &dimuon_b_HLT_2mu4);        
    
    if (isRun3){ // run3: also have mu4_mu4_noL1 trigger
        fChain->SetBranchAddress(mu4_mu4noL1_trigger_branch.c_str()           , &b_HLT_mu4_mu4noL1);
        fChain->SetBranchAddress(mu4_mu4noL1_trigger_match_branch.c_str()     , &dimuon_b_HLT_mu4_mu4noL1);
    }

    //SetBranch Status
    fChain->SetBranchStatus("*"                               ,0);//switch off all branches, then enable just the ones that we need
    fChain->SetBranchStatus("RunNumber"                       ,1);
    fChain->SetBranchStatus("muon_deltaP_overP"               ,1);
  
    fChain->SetBranchStatus("muon_pair_muon1_index"           ,1);
    fChain->SetBranchStatus("muon_pair_muon1_pt"              ,1);
    fChain->SetBranchStatus("muon_pair_muon1_eta"             ,1);
    fChain->SetBranchStatus("muon_pair_muon1_phi"             ,1);
    fChain->SetBranchStatus("muon_pair_muon1_trk_pt"              ,1);
    fChain->SetBranchStatus("muon_pair_muon1_trk_eta"             ,1);
    fChain->SetBranchStatus("muon_pair_muon1_trk_phi"             ,1);
    fChain->SetBranchStatus("muon_pair_muon1_quality"         ,1);
    fChain->SetBranchStatus("muon_pair_muon1_d0"              ,1);
    fChain->SetBranchStatus("muon_pair_muon1_z0"              ,1);
  
    fChain->SetBranchStatus("muon_pair_muon2_index"           ,1);
    fChain->SetBranchStatus("muon_pair_muon2_pt"              ,1);
    fChain->SetBranchStatus("muon_pair_muon2_eta"             ,1);
    fChain->SetBranchStatus("muon_pair_muon2_phi"             ,1);
    fChain->SetBranchStatus("muon_pair_muon2_trk_pt"              ,1);
    fChain->SetBranchStatus("muon_pair_muon2_trk_eta"             ,1);
    fChain->SetBranchStatus("muon_pair_muon2_trk_phi"             ,1);
    fChain->SetBranchStatus("muon_pair_muon2_quality"         ,1);
    fChain->SetBranchStatus("muon_pair_muon2_d0"              ,1);
    fChain->SetBranchStatus("muon_pair_muon2_z0"              ,1);
  
    fChain->SetBranchStatus(mu4_trigger_branch.c_str()                          ,1);
    fChain->SetBranchStatus(twomu4_trigger_branch.c_str()                       ,1);
    fChain->SetBranchStatus(mu4_trigger_match_branch.c_str()                    ,1);
    fChain->SetBranchStatus(twomu4_trigger_match_branch.c_str()                 ,1);

    if (isRun3){
        fChain->SetBranchStatus(mu4_mu4noL1_trigger_branch.c_str()              ,1);
        fChain->SetBranchStatus(mu4_mu4noL1_trigger_match_branch.c_str()        ,1);
    }
}

void DimuonDataAnalysisBaseClass::TrigModeToSuffixMap(){
    trig_suffix = "";

    switch(trigger_mode){
    case 1:
        trig_suffix = "_single_mu4";
        break;
    case 2:
        trig_suffix = "_mu4_mu4noL1";
        break;
    case 3:
        trig_suffix = "_2mu4";
        break;        
    default:
    	std::cerr << "Trigger mode INVALID: must be 1 / 2 / 3!" << std::endl;
    	throw std::exception();
    }
}

bool DimuonDataAnalysisBaseClass::EventTrigFilter(){
    switch(trigger_mode){
    case 1:
		return b_HLT_mu4;
    case 2:
		return b_HLT_mu4_mu4noL1;
    case 3:
		return b_HLT_2mu4;
    }

    return false;
}

bool DimuonDataAnalysisBaseClass::TrigMatch(int pair_ind, int m1_ind, int m2_ind){
	//Trigger match for muon pair or single muon

	try{
		switch (trigger_mode){
		case 1:
			if (!b_HLT_mu4) return false;
			
			if (!muon_b_HLT_mu4){ // nullptr
			  	throw std::runtime_error("The pointer muon_b_HLT_mu4 (to vector of bool) is a null pointer!");
			}
			return muon_b_HLT_mu4->at(m1_ind) || muon_b_HLT_mu4->at(m2_ind);
		case 2:
			if (!b_HLT_mu4_mu4noL1) return false;
			
			if (!dimuon_b_HLT_mu4_mu4noL1){ // nullptr
			  	throw std::runtime_error("The pointer dimuon_b_HLT_mu4_mu4noL1 (to vector of bool) is a null pointer!");
			}
			return dimuon_b_HLT_mu4_mu4noL1->at(pair_ind);
		case 3:
			if (!b_HLT_2mu4) return false;
			
			if (!dimuon_b_HLT_2mu4){ // nullptr
			  	throw std::runtime_error("The pointer dimuon_b_HLT_2mu4 (to vector of bool) is a null pointer!");
			}
			return dimuon_b_HLT_2mu4->at(pair_ind);
		default:
	    	std::cerr << "Trigger mode INVALID: must be 1 / 2 / 3!" << std::endl;
	    	throw std::exception();
		}
	}catch (const std::out_of_range& e){
	 	std::cerr << "out-of-range error occured: " << e.what() << endl;
	 	std::cout << "By default, event will not be saved." << endl;
	}catch (const std::runtime_error& e){
	 	std::cerr << "Run time error occured: " << e.what() << endl;
	 	std::cout << "By default, event will not be saved." << endl;
	}catch(...){
	 	std::cerr << "Some other unknown error occured at the trigger matching stage" << endl;
	 	std::cout << "By default, event will not be saved." << endl;
	}

	return false;
}

void DimuonDataAnalysisBaseClass::FillMuonPair(int pair_ind, std::shared_ptr<MuonPairData> const& mpair){
  	if (debug_mode) std::cout << "DEBUG: Calling DimuonDataAnalysisBaseClass::FillMuonPair" << std::endl;
  	
  	mpair->run_number = RunNumber;
  	mpair->lb         = lbn;
  	mpair->bcid       = bcid;

  	mpair->m1.ind     = muon_pair_muon1_index->at(pair_ind);
  	mpair->m2.ind     = muon_pair_muon2_index->at(pair_ind);
	
  	mpair->m1.pt    = fabs(muon_pair_muon1_pt->at(pair_ind))/1000.0;//pt of the first muon in the pair
  	mpair->m2.pt    = fabs(muon_pair_muon2_pt->at(pair_ind))/1000.0;//pt of the second muon in the pair
  	mpair->m1.eta   = muon_pair_muon1_eta->at(pair_ind);
  	mpair->m2.eta   = muon_pair_muon2_eta->at(pair_ind);
  	mpair->m1.phi   = muon_pair_muon1_phi->at(pair_ind);
  	mpair->m2.phi   = muon_pair_muon2_phi->at(pair_ind);

  	mpair->m1.d0    = muon_pair_muon1_d0 ->at(pair_ind);
  	mpair->m2.d0    = muon_pair_muon2_d0 ->at(pair_ind);
  	mpair->m1.z0    = muon_pair_muon1_z0 ->at(pair_ind);
  	mpair->m2.z0    = muon_pair_muon2_z0 ->at(pair_ind);
  	mpair->m1.charge  =(muon_pair_muon1_pt ->at(pair_ind) > 0)? 1:-1;//sign of pt stores charge
  	mpair->m2.charge  =(muon_pair_muon2_pt ->at(pair_ind) > 0)? 1:-1;//sign of pt stores charge
  	mpair->m1.quality = muon_pair_muon1_quality->at(pair_ind);
  	mpair->m2.quality = muon_pair_muon2_quality->at(pair_ind);
  	mpair->m1.dP_overP = muon_deltaP_overP->at(mpair->m1.ind);
  	mpair->m2.dP_overP = muon_deltaP_overP->at(mpair->m2.ind);
	
  	mpair->m1_trk_pt      = fabs(muon_pair_muon1_trk_pt->at(pair_ind))/1000.0;//pt of the first muon in the pair
  	mpair->m2_trk_pt      = fabs(muon_pair_muon2_trk_pt->at(pair_ind))/1000.0;//pt of the second muon in the pair
  	mpair->m1_trk_eta     = muon_pair_muon1_trk_eta->at(pair_ind);
  	mpair->m2_trk_eta     = muon_pair_muon2_trk_eta->at(pair_ind);
  	mpair->m1_trk_phi     = muon_pair_muon1_trk_phi->at(pair_ind);
  	mpair->m2_trk_phi     = muon_pair_muon2_trk_phi->at(pair_ind);

  	if (turn_on_track_charge){
  	    mpair->m1_trk_charge  =(muon_pair_muon1_trk_pt ->at(pair_ind) > 0)? 1:-1;//sign of pt stores charge
  	    mpair->m2_trk_charge  =(muon_pair_muon2_trk_pt ->at(pair_ind) > 0)? 1:-1;//sign of pt stores charge    
  	}else{ // do not turn on track charge: set to nonsense
  	  	mpair->m1_trk_charge  = 0;
  	  	mpair->m2_trk_charge  = 0;
  	}
}



bool DimuonDataAnalysisBaseClass::PassCuts(const std::shared_ptr<MuonPair>& mpair){ // by default only require medium + tag pairs passing tight WP
    return PassCuts(std::static_pointer_cast<MuonPairData>(mpair), false);
}

bool DimuonDataAnalysisBaseClass::PassCuts(const std::shared_ptr<MuonPairData>& mpair, bool requireTight){

	//require some quality cuts on the muons
	if((mpair->m1.quality&mpair->m2.quality&1  )==0) return false;//combined muon

	if (requireTight){
		if ((mpair->m1.quality&mpair->m2.quality&16)==0) return false;//tight muon
	}else{
		if ((mpair->m1.quality&mpair->m2.quality&8  )==0) return false;//Medium muon
	}

	if((mpair->m1.quality&mpair->m2.quality&32 )==0) return false;//IDCuts
	if((mpair->m1.quality&mpair->m2.quality&256)==0) return false;//MuonCuts
	h_cutAcceptance[mpair->m1.charge != mpair->m2.charge]->Fill(static_cast<int>(CutsCommon::pass_muon_quality) + 0.5, mpair->weight); // if same sign: fill the h_cutAcceptance[0] histogram; if opposite sign, fill the h_cutAcceptance[1] histogram

	if (fabs(mpair->m1.eta) > 2.4 || fabs(mpair->m2.eta) > 2.4) return false;
	h_cutAcceptance[mpair->m1.charge != mpair->m2.charge]->Fill(static_cast<int>(CutsCommon::pass_muon_eta) + 0.5, mpair->weight);

	if (mpair->m1.pt < 4 || mpair->m2.pt < 4) return false;
	h_cutAcceptance[mpair->m1.charge != mpair->m2.charge]->Fill(static_cast<int>(CutsCommon::pass_muon_pt) + 0.5, mpair->weight);
	
	if( fabs(mpair->m1.dP_overP) > pms.deltaP_overP_thrsh || fabs(mpair->m2.dP_overP) > pms.deltaP_overP_thrsh ) return false;
	h_cutAcceptance[mpair->m1.charge != mpair->m2.charge]->Fill(static_cast<int>(CutsCommon::pass_muon_dP_overP) + 0.5, mpair->weight);
	
	//cut on d0 & z0 sin(theta)
	double m1z0sinTheta = fabs(mpair->m1.z0 * sin(2.0*atan(exp(-mpair->m1.eta))));
	double m2z0sinTheta = fabs(mpair->m2.z0 * sin(2.0*atan(exp(-mpair->m2.eta))));
	bool pass_d0_z0_cuts = (fabs(mpair->m1.d0) < pms.d0cut && fabs(mpair->m2.d0) < pms.d0cut && m1z0sinTheta < pms.z0cut && m2z0sinTheta < pms.z0cut);
	if (!pass_d0_z0_cuts) return false;
	h_cutAcceptance[mpair->m1.charge != mpair->m2.charge]->Fill(static_cast<int>(CutsCommon::pass_muon_d0_z0) + 0.5, mpair->weight);

	// for both muons, require muon + track charge to agree
	if (turn_on_track_charge){
		if (mpair->m1_trk_charge != mpair->m1.charge || mpair->m2_trk_charge != mpair->m2.charge) return false;
	}
	h_cutAcceptance[mpair->m1.charge != mpair->m2.charge]->Fill(static_cast<int>(CutsCommon::pass_muon_trk_charge) + 0.5, mpair->weight); // always fill if NOT turn on track charge
	
	return true;
}


void DimuonDataAnalysisBaseClass::FillSingleMuonTree(){
  muonOutTree->Fill();
}


void DimuonDataAnalysisBaseClass::FillMuonPairTree(){
  // NECESSARY step: ALWAYS update the raw pointer with the current content of the shared pointer BEFORE FILLING OUTPUT TREES
  try{
    if (!mpair) throw std::runtime_error("FillMuonPairTree: Muon Pair doesn't exist!");
    mpair_raw_ptr = mpair.get();
    cout << "mpair pt_lead: " << mpair->pt_lead << std::endl;
    cout << "mpair_raw_ptr pt_lead: " << mpair_raw_ptr->pt_lead << std::endl;
  }catch(const std::runtime_error& e){
    std::cout << "Runtime error caught in function FillMuonPairTree: " << e.what() << std::endl;
    std::cout << "Return without filling the muon pair in the output trees!" << std::endl;
    return;
  }

  int nsign = (mpair->same_sign)? 0:1;
  muonPairOutTree[nsign]->Fill();
}


void DimuonDataAnalysisBaseClass::ProcessData(){

	Long64_t nentries = fChain->GetEntries();//number of events
	// for (Long64_t jentry=0; jentry<nentries;jentry++) {//loop over the events
	for (Long64_t jentry=0; jentry<10000;jentry++) {//loop over the events

		if(jentry%100==0) cout<<"Processing "<<jentry<<" event out of "<<nentries<<" events"<<std::endl;
		// if(jentry%100000==0) cout<<"Processing "<<jentry<<" event out of "<<nentries<<" events"<<std::endl;

		int num_bytes = fChain->GetEntry(jentry);//read in an event
		if(num_bytes==0){
		  	std::cout<<"Error:: Read in event has size of zero bytes,  quitting"<<std::endl;
		  	throw std::exception();
		}
		
		muon_pair_list_cur_event_pre_resonance_cut.clear();
		resonance_tagged_muon_index_list.clear(); // MUST CLEAR for each event!!
        resonance_tagged_muon_index_list_v2.clear(); // MUST CLEAR for each event!!

		//trigger requirement for event
		EventTrigFilter();

		std::vector<int> muon_index_list = {};
		std::vector<int>::iterator it;

		int NPairs=muon_pair_muon1_pt->size();//number of muon pairs in the event

		for(int pair_ind=0;pair_ind<NPairs;pair_ind++){//first loop over all muon-pairs in the event

			mpair = MakeMuonPair();

			FillMuonPair(pair_ind, mpair);
			mpair->m1.ev_num = jentry;
			mpair->m2.ev_num = jentry;

			//------------------------------------------------------------

			h_cutAcceptance[mpair->m1.charge != mpair->m2.charge]->Fill(static_cast<int>(CutsCommon::nocut) + 0.5, mpair->weight);

			//Trigger match for muon pair
			if(!TrigMatch(pair_ind, mpair->m1.ind, mpair->m2.ind)) continue;
			h_cutAcceptance[mpair->m1.charge != mpair->m2.charge]->Fill(static_cast<int>(CutsCommon::pass_trigger_match) + 0.5, mpair->weight);

			if (!PassCuts(mpair))continue;
			
			mpair->passTight = (mpair->m1.quality&mpair->m2.quality&16); //tag tight muon pairs
		
			//------------------------------------------------------------

			//Two things at this step: 
			//1) sort pt, eta, phi by pt
			//2) update the muon-pair values
			mpair->Update();
			mpair->passmu4mu4noL1 = dimuon_b_HLT_mu4_mu4noL1->at(pair_ind);
			mpair->pass2mu4 = dimuon_b_HLT_2mu4->at(pair_ind);
			mpair->mu1PassSingle = muon_b_HLT_mu4->at(mpair->m1.ind);
			mpair->mu2PassSingle = muon_b_HLT_mu4->at(mpair->m2.ind);

			// resonance tag
			ResonanceTagging(mpair);
            ResonanceTaggingV2(mpair);

			// photo-production cut
			if (isPbPb){
			  	if (IsPhotoProduction(mpair)) continue;
			  	h_cutAcceptance[mpair->m1.charge != mpair->m2.charge]->Fill(static_cast<int>(CutsPbPb::pass_photoprod) + 0.5, mpair->weight);      	
			}

			muon_pair_list_cur_event_pre_resonance_cut.push_back(std::move(mpair));

		} // finish first loop over all muon pairs

		for(int pair_ind = 0; pair_ind < muon_pair_list_cur_event_pre_resonance_cut.size(); pair_ind++){//second loop over all muon-pairs in the event
			// discard the pair if either muon is resonance-tagged

			mpair = std::move(muon_pair_list_cur_event_pre_resonance_cut.at(pair_ind));

			if (!mpair){
			  	std::cerr << "mpair at second muon-pair loop NOT found! Return without resonance-tag checking or pair analysis." << std::endl;
			  	continue;
			}

			std::vector<int>::iterator itres_m1;
			std::vector<int>::iterator itres_m2;

            // apply resonance cuts if resonance_cut_mode != 0
            if (resonance_cut_mode == 1){ // apply old cuts
                itres_m1 = std::find(resonance_tagged_muon_index_list.begin(),resonance_tagged_muon_index_list.end(),mpair->m1.ind);
                if(itres_m1 != resonance_tagged_muon_index_list.end())  continue;

                itres_m2 = std::find(resonance_tagged_muon_index_list.begin(),resonance_tagged_muon_index_list.end(),mpair->m2.ind);
                if(itres_m2 != resonance_tagged_muon_index_list.end())  continue;
            } else if (resonance_cut_mode == 2){ // apply new cuts
                itres_m1 = std::find(resonance_tagged_muon_index_list_v2.begin(),resonance_tagged_muon_index_list_v2.end(),mpair->m1.ind);
                if(itres_m1 != resonance_tagged_muon_index_list_v2.end())  continue;

                itres_m2 = std::find(resonance_tagged_muon_index_list_v2.begin(),resonance_tagged_muon_index_list_v2.end(),mpair->m2.ind);
                if(itres_m2 != resonance_tagged_muon_index_list_v2.end())  continue;
            }

			int pass_resonance_ind = isPbPb? static_cast<int>(CutsPbPb::pass_resonance) : static_cast<int>(CutsPP::pass_resonance);
			h_cutAcceptance[mpair->m1.charge != mpair->m2.charge]->Fill( + 0.5, mpair->weight);

			//------------------------------------------------------------
			// perform additional pair-level analysis if needed
			PerformAdditionalPairAnalysis();

			if(output_single_muon_tree){
			  it = std::find(muon_index_list.begin(),muon_index_list.end(),mpair->m1.ind);
			  if(it == muon_index_list.end()){ //muon1 index NOT found
			    	muon_index_list.push_back(mpair->m1.ind);
			    	tempmuon = &(mpair->m1);
			    	FillSingleMuonTree();
			  }
			  it = std::find(muon_index_list.begin(),muon_index_list.end(),mpair->m2.ind);
			  if(it == muon_index_list.end()){ //muon1 index NOT found
			    	muon_index_list.push_back(mpair->m2.ind);
			    	tempmuon = &(mpair->m2);
			    	FillSingleMuonTree();
			  }
			}
			else{ // fill muon pair trees
			 	FillMuonPairTree();
			}
		}

		mpair.reset(); // make sure the reference to the last muon pair gets reset
	}//loop over events
}