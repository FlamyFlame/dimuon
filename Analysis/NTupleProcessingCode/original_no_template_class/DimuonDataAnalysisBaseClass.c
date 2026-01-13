#include "DimuonDataAnalysisBaseClass.h"

void DimuonDataAnalysisBaseClass::PrintInstructions(){
	std::string datatype = isPbPb? "PbPb" : "PP";
	std::cout << datatype << " Data Ntuple processing script:" << std::endl;
    std::cout << "The following variable(s) are required by constructor:" << std::endl;
    std::cout << "--> file_batch: [INT] Decides which run3-file batch to process, only has effect when isRun3 is true" << std::endl;
    std::cout << "--> 					PbPb: file batch: 1-6 for 2023 data, 1-9 for 2024 data, 1-7 for 2015/2018 data" << std::endl;
    std::cout << "--> 					PP:   file batch: 1-4 for 2024 data, 1-3 for 2017 data" << std::endl;
    std::cout << "--> run_year: [INT] year of data taking" << std::endl;
    std::cout << std::endl;

    std::cout << "The following public variable(s) **MUST** be checked:" << std::endl;
    std::cout << "--> trigger_mode: INT, value = 0,1,2,3" << std::endl;
    std::cout << "                       value = 0 (set true if isMinBias): require MinBias trigger: tag if either muon files single-muon trigger" << std::endl;
    std::cout << "                       value = 1 (DEFAULT): require single-muon trigger: tag which muon files trigger & if pair passes mu4_mu4_noL1 & 2mu4" << std::endl;
    std::cout << "                       value = 2: require mu4_mu4_noL1" << std::endl;
    std::cout << "                       value = 3: require 2mu4" << std::endl;
    std::cout << std::endl;

    std::cout << "The following public variable(s) should be checked:" << std::endl;
    std::cout << "--> isMinBias: [BOOL] if true, use MinBias data + perform single-muon analysis + trigger mode set to 0" << std::endl;
    std::cout << "--> output_single_muon_tree: [BOOL] if true, output single-muon tree" << std::endl;
    std::cout << "                                    if false (DEFAULT), output muon-pair tree" << std::endl;
    if (isPbPb) std::cout << "--> turn_on_ctr_binned_tree_writing (only for PbPb): [BOOL] if true, output centrality-binned muon-pair trees (DEFAULT FALSE)" << std::endl;
    std::cout << "--> requireTight: boolean, default false - if true: require tight WP; false; require medium WP" << std::endl;
    std::cout << "--> resonance_cut_mode: integer that determines which set of resonant cuts to apply" << std::endl;
    std::cout << "        * resonance_cut_mode = 0: NO resonance cut" << std::endl;
    std::cout << "        * resonance_cut_mode = 1: old resonance cut (default)" << std::endl;
    std::cout << "        * resonance_cut_mode = 2: new resonance cut" << std::endl;
    std::cout << "        If resonance_cut_mode value is outside {0,1,2}: assume default option" << std::endl;
    std::cout << "--> pbpb24_mu4_NO_trig_calc: boolean, default false - if true: use Pb+Pb 24 single-mu4 data for nominal analysis, not trigger efficiency evaluation" << std::endl;
    std::cout << "--> filter_out_photo_resn_for_trig_effcy: boolean, default false - if true: do NOT filter out photoproduction / resonance pairs for trigger efficiency study" << std::endl;

    std::cout << std::endl;

    if (isPbPb){
	    std::cout << "if run_year == 23, output files will be written to /usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_2023" << std::endl;
	    std::cout << "if run_year == 24, output files will be written to /usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_2024" << std::endl;
	    std::cout << "else,      output files will be written to /usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_run2" << std::endl;
    }else{
	    std::cout << "if isRun3, output files will be written to /usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024" << std::endl;
	    std::cout << "else,      output files will be written to /usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_run2" << std::endl;
    }
}

void DimuonDataAnalysisBaseClass::InitParams(){
    run_year = run_year % 2000; // e.g, 2023 --> 23
    if (run_year < 15) throw std::runtime_error("Run year invalid: must be >= 2015!");

    isRun3 = (run_year > 20);

    if (isMinBias) trigger_mode = 0; // if use MB data, do not require any muon trigger
    
    use_mu6_for_trg_eff = ((isPbPb && run_year == 23) || (!isPbPb && run_year == 24)); // Pb+Pb 23 or pp24
    use_mu8_for_trg_eff = (!isPbPb && run_year == 24); // pp24
    if (use_mu6_for_trg_eff) std::cout << "Using HLT_mu6_L1MU3V as a support trigger!" << std::endl;
    if (use_mu8_for_trg_eff) std::cout << "Using HLT_mu8_L1MU5VF as a support trigger!" << std::endl;
    
    if (!(isPbPb && run_year == 24 && pbpb24_mu4_NO_trig_calc)) trigger_effcy_calc = (trigger_mode == 0 || trigger_mode == 1);
    if (trigger_effcy_calc) resonance_cut_mode = 2; // for sub-GeV mass region, apply narrow-window cuts for each individual resonance peak

    if (pbpb24_mu4_NO_trig_calc && trigger_mode == 1) std::cout << "For Pb+Pb 2024 data, using single_mu4 for NOMINAL analysis (not trigger efficiency)." << std::endl;
}

void DimuonDataAnalysisBaseClass::InitInput(){
    std::string base_dir = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/";
    std::string dt_str = isPbPb? "pbpb" : "pp";
    std::string data_subdir = isRun3? dt_str + "_20" + std::to_string(run_year) + "/" : dt_str + "_run2/";

    data_dir = base_dir + data_subdir;

	TChainFill();
	if (!fChain){
		std::cerr << "FATAL:: TChain for analysis is nullptr! Throwing exception." << std::endl;
    	throw std::exception();
	}

    if (isMinBias){
        InitInputBranchesSingleMuonAnalysis();
    }else{
        InitInputBranchesDimuonAnalysis();
    }

}

void DimuonDataAnalysisBaseClass::InitInputBranchesSingleMuonAnalysis(){
    fChain->SetBranchAddress("muon_pt"              , &muon_pt);
    fChain->SetBranchAddress("muon_eta"             , &muon_eta);
    fChain->SetBranchAddress("muon_phi"             , &muon_phi);
    fChain->SetBranchAddress("muon_quality"         , &muon_quality);
    fChain->SetBranchAddress("muon_deltaP_overP"          , &muon_deltaP_overP);
    fChain->SetBranchAddress("muon_d0"              , &muon_d0);
    fChain->SetBranchAddress("muon_z0"              , &muon_z0);

    std::string mu4_trigger_name = (isRun3)? "HLT_mu4_L1MU3V" : "HLT_mu4";
    std::string mu4_trigger_branch = "b_" + mu4_trigger_name;
    std::string mu4_trigger_match_branch = "muon_b_" + mu4_trigger_name;

    fChain->SetBranchAddress(mu4_trigger_branch.c_str()                       , &b_HLT_mu4);
    fChain->SetBranchAddress(mu4_trigger_match_branch.c_str()                 , &muon_b_HLT_mu4);

    fChain->SetBranchStatus("*"                     ,0);//switch off all branches, then enable just the ones that we need
    fChain->SetBranchStatus("muon_pt"               ,1);
    fChain->SetBranchStatus("muon_eta"              ,1);
    fChain->SetBranchStatus("muon_phi"              ,1);
    fChain->SetBranchStatus("muon_quality"          ,1);
    fChain->SetBranchStatus("muon_deltaP_overP"     ,1);
    fChain->SetBranchStatus("muon_d0"               ,1);
    fChain->SetBranchStatus("muon_z0"               ,1);

    fChain->SetBranchStatus(mu4_trigger_branch.c_str()                          ,1);
    fChain->SetBranchStatus(mu4_trigger_match_branch.c_str()                    ,1);
}

void DimuonDataAnalysisBaseClass::InitInputBranchesDimuonAnalysis(){

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
    
    if (isRun3 || isPbPb){
        fChain->SetBranchAddress(mu4_mu4noL1_trigger_branch.c_str()           , &b_HLT_mu4_mu4noL1);
        fChain->SetBranchAddress(mu4_mu4noL1_trigger_match_branch.c_str()     , &dimuon_b_HLT_mu4_mu4noL1);
    }

    if (use_mu6_for_trg_eff){ // only use mu6 for Pb+Pb 23 for now
        fChain->SetBranchAddress("b_HLT_mu6_L1MU3V"                      , &b_HLT_mu6_L1MU3V);
        fChain->SetBranchAddress("muon_b_HLT_mu6_L1MU3V"                 , &muon_b_HLT_mu6_L1MU3V);
    }

    if (use_mu8_for_trg_eff){ // only use mu6 for Pb+Pb 23 for now
        fChain->SetBranchAddress("b_HLT_mu8_L1MU5VF"                      , &b_HLT_mu8_L1MU5VF);
        fChain->SetBranchAddress("muon_b_HLT_mu8_L1MU5VF"                 , &muon_b_HLT_mu8_L1MU5VF);
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

    if (isRun3 || isPbPb){
        fChain->SetBranchStatus(mu4_mu4noL1_trigger_branch.c_str()              ,1);
        fChain->SetBranchStatus(mu4_mu4noL1_trigger_match_branch.c_str()        ,1);
    }
    
    if (use_mu6_for_trg_eff){ // only use mu6 for Pb+Pb 23 for now
        fChain->SetBranchStatus("b_HLT_mu6_L1MU3V",         1);
        fChain->SetBranchStatus("muon_b_HLT_mu6_L1MU3V",    1);
    }

    if (use_mu8_for_trg_eff){ // only use mu6 for Pb+Pb 23 for now
        fChain->SetBranchStatus("b_HLT_mu8_L1MU5VF",         1);
        fChain->SetBranchStatus("muon_b_HLT_mu8_L1MU5VF",    1);
    }
}

void DimuonDataAnalysisBaseClass::TrigModeToSuffixMap(){
    trig_suffix = "";

    switch(trigger_mode){
    case 0:
        trig_suffix = "_min_bias";
        break;
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
    	std::cerr << "Trigger mode INVALID: must be 0 / 1 / 2 / 3!" << std::endl;
    	throw std::exception();
    }
    
    if (trigger_effcy_calc && !filter_out_photo_resn_for_trig_effcy) trig_suffix += "_no_photo_resn_cuts"; // if not filter out photoprod/resn pairs for trigger efficiency study
}


void DimuonDataAnalysisBaseClass::InitOutput() {
    // Common suffix setup
    std::string tight_suffix = (requireTight) ? "_tight" : "";
    TrigModeToSuffixMap(); // Sets trig_suffix based on trigger_mode

    // Determine resonance cut suffix
    std::string resonance_cut_suffix;
    // if (!isPbPb) { // commented out PP requirement 9/25/2025
    switch (resonance_cut_mode) {
        case 0: resonance_cut_suffix = "_no_res_cut"; break;
        case 1: resonance_cut_suffix = ""; break;
        case 2: resonance_cut_suffix = "_res_cut_v2"; break;
        default: 
			resonance_cut_mode = 1;
        	resonance_cut_suffix = "";
        	break;
    }
    // }

    // Construct output paths
    std::string file_name_base = output_single_muon_tree ? "single_muon_trees" : "muon_pairs";
    std::string run_suffix;

    if (!isPbPb) { // PP
        run_suffix = isRun3 ? "_pp_2024" : "_pp_2017";
    } else { // PbPb
        run_suffix = "_pbpb_20" + std::to_string(run_year);
    }
    
	std::string file_batch_suffix = "_part" + std::to_string(file_batch);
    std::string outfile_ending = run_suffix + file_batch_suffix + trig_suffix + tight_suffix + resonance_cut_suffix + ".root";
    output_file_path = data_dir + file_name_base + outfile_ending;
    output_hist_file_path = data_dir + "hists_cut_acceptance" + outfile_ending;

    // Histograms file
    m_outHistFile = new TFile(output_hist_file_path.c_str(), "recreate");
    DimuonAnalysisBaseClass::InitOutput(); // Base histogram setup

    // Create muon/muon-pair output file
    m_outfile = new TFile(output_file_path.c_str(), "recreate");
    
    // Create main trees
    if (output_single_muon_tree) {
        muonOutTree = new TTree("muon_tree", "All single muons");
        muonOutTree->Branch("MuonObj", &tempmuon);
    } else {
        for (unsigned int ksign = 0; ksign < ParamsSet::nSigns; ksign++) {
            muonPairOutTree[ksign] = new TTree(Form("muon_pair_tree_sign%u", ksign+1), Form("All muon pairs, sign%u", ksign+1));
            muonPairOutTree[ksign]->Branch("MuonPairObj", &mpair_raw_ptr);
        }
    }
}

bool DimuonDataAnalysisBaseClass::TrigMatch(int pair_ind, int m1_ind, int m2_ind){
	//Trigger match for muon pair or single muon

	try{
		switch (trigger_mode){
        case 0:
            return true;
		case 1:{
            // Default case: first consider mu4
            bool pass_mu4 = b_HLT_mu4;

            if (!muon_b_HLT_mu4){ // nullptr
                throw std::runtime_error("TrigMatch:: muon_b_HLT_mu4 is a null pointer!");
            }

            pass_mu4 &= muon_b_HLT_mu4->at(m1_ind) || muon_b_HLT_mu4->at(m2_ind);

            if (pass_mu4 || !use_mu6_for_trg_eff) return pass_mu4; // if pass mu4 or not using higher pT triggers, return mu4 result

            // Case II: Not passing mu4 && using mu6 (mu4 prescaled) --> look at mu6
            bool pass_mu6 = b_HLT_mu6_L1MU3V;

            if (!muon_b_HLT_mu6_L1MU3V){ // nullptr
                throw std::runtime_error("TrigMatch:: muon_b_HLT_mu6_L1MU3V is a null pointer!");
            }
            if (!mpair) throw std::runtime_error("TrigMatch: Muon Pair doesn't exist!");

            pass_mu6 &= ((muon_b_HLT_mu6_L1MU3V->at(m1_ind) && mpair->m1.pt >= 6) || (muon_b_HLT_mu6_L1MU3V->at(m2_ind) && mpair->m2.pt >= 6));

            if (pass_mu6 || !use_mu8_for_trg_eff) return pass_mu6; // if pass mu4 or not using higher 

            // Case III: Not passing either mu4 or mu6 && using mu8 (mu4 & mu6 prescaled) --> look at mu8
            
            if (!b_HLT_mu8_L1MU5VF) return false;

            if (!muon_b_HLT_mu8_L1MU5VF){ // nullptr
                throw std::runtime_error("TrigMatch:: muon_b_HLT_mu8_L1MU5VF is a null pointer!");
            }

            return ((muon_b_HLT_mu8_L1MU5VF->at(m1_ind) && mpair->m1.pt >= 8) || (muon_b_HLT_mu8_L1MU5VF->at(m2_ind) && mpair->m2.pt >= 8));
        }
		case 2:
            if (!isPbPb && !isRun3){ // pp Run2: no mu4_mu4noL1 trigger
                std::cerr << "PP Run2 data: no mu4_mu4_noL1 trigger!" << std::endl;
                throw std::exception();
            }

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

void DimuonDataAnalysisBaseClass::FillMuonPairTree(){
  // NECESSARY step: ALWAYS update the raw pointer with the current content of the shared pointer BEFORE FILLING OUTPUT TREES
  try{
    if (!mpair) throw std::runtime_error("FillMuonPairTree: Muon Pair doesn't exist!");
    mpair_raw_ptr = mpair.get();
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
	for (Long64_t jentry=0; jentry<nentries;jentry++) {//loop over the events
	// for (Long64_t jentry=0; jentry<1000;jentry++) {//loop over the events

		if(jentry%100000==0) cout<<"Processing "<<jentry<<" event out of "<<nentries<<" events"<<std::endl;

		int num_bytes = fChain->GetEntry(jentry);//read in an event
		if(num_bytes==0){
		  	std::cout<<"Error:: Read in event has size of zero bytes,  quitting"<<std::endl;
		  	throw std::exception();
		}
		
		muon_pair_list_cur_event_pre_resonance_cut.clear();
		resonance_tagged_muon_index_list.clear(); // MUST CLEAR for each event!!
        resonance_tagged_muon_index_list_v2.clear(); // MUST CLEAR for each event!!

		std::vector<int> muon_index_list = {};
		std::vector<int>::iterator it;

        if (isMinBias){ // MinBias data --> perform offline-single-muon analysis

        }else{ // HardProbe data --> perform offline-dimuon analysis
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
    			mpair->passmu4mu4noL1 = (!isPbPb && !isRun3)? false : dimuon_b_HLT_mu4_mu4noL1->at(pair_ind); // pp run2: no mu4_mu4noL1 trigger
    			mpair->pass2mu4 = dimuon_b_HLT_2mu4->at(pair_ind);
    	
        		mpair->m1.passmu4 = muon_b_HLT_mu4->at(mpair->m1.ind);
                if (use_mu6_for_trg_eff) mpair->m1.passmu4 |= (muon_b_HLT_mu6_L1MU3V->at(mpair->m1.ind) && mpair->m1.pt > 6);
                if (use_mu8_for_trg_eff) mpair->m1.passmu4 |= (muon_b_HLT_mu8_L1MU5VF->at(mpair->m1.ind) && mpair->m1.pt > 8);
         
                mpair->m2.passmu4 = muon_b_HLT_mu4->at(mpair->m2.ind);
                if (use_mu6_for_trg_eff) mpair->m2.passmu4 |= (muon_b_HLT_mu6_L1MU3V->at(mpair->m2.ind) && mpair->m2.pt > 6);
                if (use_mu8_for_trg_eff) mpair->m2.passmu4 |= (muon_b_HLT_mu8_L1MU5VF->at(mpair->m2.ind) && mpair->m2.pt > 8);
                

                mpair->passSeparated = (mpair->dr > 0.8);
                mpair->passSeparatedDeta = (abs(mpair->deta) > 0.8);

    			// resonance tag
    			ResonanceTagging(mpair);
                ResonanceTaggingV2(mpair);

    			// photo-production cut
    			if (isPbPb){
    			  	if (!(trigger_effcy_calc && !filter_out_photo_resn_for_trig_effcy) && IsPhotoProduction(mpair)) continue;
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
                if (!(trigger_effcy_calc && !filter_out_photo_resn_for_trig_effcy)){ // perform resn cuts
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
        }
	}//loop over events
}


void DimuonDataAnalysisBaseClass::Run(){

    clock_t start, end;
    double cpu_time_used;
    start = clock();

    InitParams();
    InitInput();
    InitOutput();
    ProcessData();
    HistAdjust();

    m_outfile->Write();
    m_outHistFile->Write();
    m_outfile->Close();
    m_outHistFile->Close();
    std::cout << "Output muon-pair trees have been written to: " << output_file_path << std::endl;
    std::cout << "Output histograms have been written to: " << output_hist_file_path << std::endl;

    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    std::cout << "CPU time used (in seconds) is " << cpu_time_used << std::endl;

}
