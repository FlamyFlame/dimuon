#include "DimuonDataAlgCoreT.h"

template <class PairT, class MuonT, class Derived, class... Extras>
void DimuonDataAlgCoreT<PairT, MuonT, Derived, Extras...>::PrintInstructions(){
    PrintInstructionsHook();
}

template <class PairT, class MuonT, class Derived, class... Extras>
void DimuonDataAlgCoreT<PairT, MuonT, Derived, Extras...>::PrintInstructions_DataCore(){
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
    std::cout << "--> filter_out_photo_resn_for_trig_effcy: boolean, default true - if true: filter out photoproduction / resonance pairs even in trigger efficiency study mode" << std::endl;

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

template <class PairT, class MuonT, class Derived, class... Extras>
void DimuonDataAlgCoreT<PairT, MuonT, Derived, Extras...>::InitParams_DataCore(){
    this->isMC = false; // is MC: regardless of truth only or fullsim
    
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

    if (mindR_trig != 0.01 && mindR_trig != 0.02)
        throw std::runtime_error("mindR_trig must be 0.01 or 0.02! Got: " + std::to_string(mindR_trig));
    std::cout << "mindR_trig = " << mindR_trig << " (will use matching branches if available)" << std::endl;
}

template <class PairT, class MuonT, class Derived, class... Extras>
void DimuonDataAlgCoreT<PairT, MuonT, Derived, Extras...>::InitInput_DataCore(){
    std::string base_dir = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/";
    std::string dt_str = isPbPb? "pbpb" : "pp";
    std::string data_subdir = isRun3? dt_str + "_20" + std::to_string(run_year) + "/" : dt_str + "_run2/";

    data_dir = base_dir + data_subdir;

	TChainFill();
	if (!fChainRef()){
		std::cerr << "FATAL:: TChain for analysis is nullptr! Throwing exception." << std::endl;
    	throw std::exception();
	}

    if (isMinBias){
        self().InitInputBranchesSingleMuonAnalysisHook();
    }else{
        self().InitInputBranchesDimuonAnalysisHook();
    }
}

template <class PairT, class MuonT, class Derived, class... Extras>
void DimuonDataAlgCoreT<PairT, MuonT, Derived, Extras...>::InitInputBranchesSingleMuonAnalysis_DataCore(){
    fChainRef()->SetBranchAddress("muon_pt"              , &muon_pt);
    fChainRef()->SetBranchAddress("muon_eta"             , &muon_eta);
    fChainRef()->SetBranchAddress("muon_phi"             , &muon_phi);
    fChainRef()->SetBranchAddress("muon_quality"         , &muon_quality);
    fChainRef()->SetBranchAddress("muon_deltaP_overP"          , &muon_deltaP_overP);
    fChainRef()->SetBranchAddress("muon_d0"              , &muon_d0);
    fChainRef()->SetBranchAddress("muon_z0"              , &muon_z0);

    std::string mu4_trigger_name = (isRun3)? "HLT_mu4_L1MU3V" : "HLT_mu4";
    std::string mu4_trigger_branch = "b_" + mu4_trigger_name;
    std::string mu4_trigger_match_branch = "muon_b_" + mu4_trigger_name;

    fChainRef()->SetBranchAddress(mu4_trigger_branch.c_str()                       , &b_HLT_mu4);
    fChainRef()->SetBranchAddress(mu4_trigger_match_branch.c_str()                 , &muon_b_HLT_mu4);

    fChainRef()->SetBranchStatus("*"                     ,0);//switch off all branches, then enable just the ones that we need
    fChainRef()->SetBranchStatus("muon_pt"               ,1);
    fChainRef()->SetBranchStatus("muon_eta"              ,1);
    fChainRef()->SetBranchStatus("muon_phi"              ,1);
    fChainRef()->SetBranchStatus("muon_quality"          ,1);
    fChainRef()->SetBranchStatus("muon_deltaP_overP"     ,1);
    fChainRef()->SetBranchStatus("muon_d0"               ,1);
    fChainRef()->SetBranchStatus("muon_z0"               ,1);

    fChainRef()->SetBranchStatus(mu4_trigger_branch.c_str()                          ,1);
    fChainRef()->SetBranchStatus(mu4_trigger_match_branch.c_str()                    ,1);
}

template <class PairT, class MuonT, class Derived, class... Extras>
void DimuonDataAlgCoreT<PairT, MuonT, Derived, Extras...>::InitInputBranchesDimuonAnalysis_DataCore(){

    fChainRef()->SetBranchAddress("RunNumber"                   , &RunNumber);
    fChainRef()->SetBranchAddress("lbn"                         , &lbn);
    fChainRef()->SetBranchAddress("bcid"                        , &bcid);
    
    fChainRef()->SetBranchAddress("muon_deltaP_overP"           , &muon_deltaP_overP);
  
    fChainRef()->SetBranchAddress("muon_pair_muon1_index"       , &muon_pair_muon1_index);
    fChainRef()->SetBranchAddress("muon_pair_muon1_pt"          , &muon_pair_muon1_pt);
    fChainRef()->SetBranchAddress("muon_pair_muon1_eta"         , &muon_pair_muon1_eta);
    fChainRef()->SetBranchAddress("muon_pair_muon1_phi"         , &muon_pair_muon1_phi);
    fChainRef()->SetBranchAddress("muon_pair_muon1_quality"     , &muon_pair_muon1_quality);
    fChainRef()->SetBranchAddress("muon_pair_muon1_d0"          , &muon_pair_muon1_d0);
    fChainRef()->SetBranchAddress("muon_pair_muon1_z0"          , &muon_pair_muon1_z0);
    fChainRef()->SetBranchAddress("muon_pair_muon1_trk_pt"      , &muon_pair_muon1_trk_pt);
    fChainRef()->SetBranchAddress("muon_pair_muon1_trk_eta"     , &muon_pair_muon1_trk_eta);
    fChainRef()->SetBranchAddress("muon_pair_muon1_trk_phi"     , &muon_pair_muon1_trk_phi);
  
    fChainRef()->SetBranchAddress("muon_pair_muon2_index"       , &muon_pair_muon2_index);
    fChainRef()->SetBranchAddress("muon_pair_muon2_pt"          , &muon_pair_muon2_pt);
    fChainRef()->SetBranchAddress("muon_pair_muon2_eta"         , &muon_pair_muon2_eta);
    fChainRef()->SetBranchAddress("muon_pair_muon2_phi"         , &muon_pair_muon2_phi);
    fChainRef()->SetBranchAddress("muon_pair_muon2_quality"     , &muon_pair_muon2_quality);
    fChainRef()->SetBranchAddress("muon_pair_muon2_d0"          , &muon_pair_muon2_d0);
    fChainRef()->SetBranchAddress("muon_pair_muon2_z0"          , &muon_pair_muon2_z0);
    fChainRef()->SetBranchAddress("muon_pair_muon2_trk_pt"      , &muon_pair_muon2_trk_pt);
    fChainRef()->SetBranchAddress("muon_pair_muon2_trk_eta"     , &muon_pair_muon2_trk_eta);
    fChainRef()->SetBranchAddress("muon_pair_muon2_trk_phi"     , &muon_pair_muon2_trk_phi);
    
    std::string mu4_trigger_name = (isRun3)? "HLT_mu4_L1MU3V" : "HLT_mu4";
    std::string mu4_trigger_branch = "b_" + mu4_trigger_name;
    std::string mu4_trigger_match_branch = "muon_b_" + mu4_trigger_name;

    std::string mu4_mu4noL1_trigger_name = (isRun3)? "HLT_mu4_mu4noL1_L1MU3V" : "HLT_mu4_mu4noL1";
    std::string mu4_mu4noL1_trigger_branch = "b_" + mu4_mu4noL1_trigger_name;
    std::string mu4_mu4noL1_trigger_match_branch = "dimuon_b_" + mu4_mu4noL1_trigger_name;

    std::string twomu4_trigger_name = (isRun3)? "HLT_2mu4_L12MU3V" : "HLT_2mu4";
    std::string twomu4_trigger_branch = "b_" + twomu4_trigger_name;
    std::string twomu4_trigger_match_branch = "dimuon_b_" + twomu4_trigger_name;

    // mindR suffix strings: "_dR_0_02" for leg branches, "_0_02" for order-insensitive mindR branch
    std::string mindR_str = (mindR_trig == 0.01) ? "0_01" : "0_02";
    std::string leg_dR_suffix  = "_dR_" + mindR_str;       // e.g. "_dR_0_02"
    std::string mindR_suffix   = "_" + mindR_str;           // e.g. "_0_02"

    fChainRef()->SetBranchAddress(mu4_trigger_branch.c_str()                       , &b_HLT_mu4);
    fChainRef()->SetBranchAddress(twomu4_trigger_branch.c_str()                    , &b_HLT_2mu4);
    fChainRef()->SetBranchAddress(mu4_trigger_match_branch.c_str()                 , &muon_b_HLT_mu4);

    // --- 2mu4: prefer mindR branch, fall back to old branch ---
    if (isRun3 || isPbPb) {
        std::string twomu4_mindR_branch = twomu4_trigger_match_branch + mindR_suffix;
        if (fChainRef()->GetBranch(twomu4_mindR_branch.c_str())) {
            fChainRef()->SetBranchAddress(twomu4_mindR_branch.c_str(), &dimuon_b_2mu4_mindR);
            use_mindR_branch_2mu4 = true;
            use_mindR_suffix_in_output = true;
            std::cout << "INFO: 2mu4 mindR branch found: " << twomu4_mindR_branch << std::endl;
        } else {
            fChainRef()->SetBranchAddress(twomu4_trigger_match_branch.c_str(), &dimuon_b_HLT_2mu4);
            std::cerr << "SERIOUS WARNING: 2mu4 mindR branch '" << twomu4_mindR_branch
                      << "' not found; falling back to old branch (likely outdated skimmed data)" << std::endl;
        }
    } else {
        fChainRef()->SetBranchAddress(twomu4_trigger_match_branch.c_str(), &dimuon_b_HLT_2mu4);
    }

    // --- mu4_mu4noL1: try leg branches first, then mindR-only, then old branch ---
    if (isRun3 || isPbPb){
        fChainRef()->SetBranchAddress(mu4_mu4noL1_trigger_branch.c_str(), &b_HLT_mu4_mu4noL1);

        std::string mu1Leg1_branch = mu4_mu4noL1_trigger_match_branch + "_mu1passLeg1" + leg_dR_suffix;
        std::string mu1Leg2_branch = mu4_mu4noL1_trigger_match_branch + "_mu1passLeg2" + leg_dR_suffix;
        std::string mu2Leg1_branch = mu4_mu4noL1_trigger_match_branch + "_mu2passLeg1" + leg_dR_suffix;
        std::string mu2Leg2_branch = mu4_mu4noL1_trigger_match_branch + "_mu2passLeg2" + leg_dR_suffix;

        if (fChainRef()->GetBranch(mu1Leg2_branch.c_str()) && use_per_leg_matching) {
            // Per-leg branches found and use_per_leg_matching enabled - use per-leg matching
            fChainRef()->SetBranchAddress(mu1Leg1_branch.c_str(), &dimuon_b_mu4_mu4noL1_mu1passLeg1);
            fChainRef()->SetBranchAddress(mu1Leg2_branch.c_str(), &dimuon_b_mu4_mu4noL1_mu1passLeg2);
            fChainRef()->SetBranchAddress(mu2Leg1_branch.c_str(), &dimuon_b_mu4_mu4noL1_mu2passLeg1);
            fChainRef()->SetBranchAddress(mu2Leg2_branch.c_str(), &dimuon_b_mu4_mu4noL1_mu2passLeg2);
            use_leg_branches_mu4_mu4noL1 = true;
            use_mindR_suffix_in_output = true;
            std::cout << "INFO: mu4_mu4noL1 per-leg branches enabled (mindR=" << mindR_trig << ")" << std::endl;
        } else {
            if (fChainRef()->GetBranch(mu1Leg2_branch.c_str()))
                std::cout << "INFO: mu4_mu4noL1 per-leg branches exist but use_per_leg_matching=false; using pair-level branch" << std::endl;
            // No leg branches - try mindR-only (order-insensitive) branch
            std::string mindR_only_branch = mu4_mu4noL1_trigger_match_branch + mindR_suffix;
            if (fChainRef()->GetBranch(mindR_only_branch.c_str())) {
                fChainRef()->SetBranchAddress(mindR_only_branch.c_str(), &dimuon_b_mu4_mu4noL1_mindR_only);
                use_mindR_only_branch_mu4_mu4noL1 = true;
                use_mindR_suffix_in_output = true;
                std::cout << "WARNING: mu4_mu4noL1 leg branches not found; using mindR-only branch '"
                          << mindR_only_branch << "' (no per-muon leg info available)" << std::endl;
            } else {
                // Full fallback: old branch without mindR
                fChainRef()->SetBranchAddress(mu4_mu4noL1_trigger_match_branch.c_str(), &dimuon_b_HLT_mu4_mu4noL1);
                std::cerr << "SERIOUS WARNING: mu4_mu4noL1 mindR branches not found; falling back to old branch '"
                          << mu4_mu4noL1_trigger_match_branch << "' (likely outdated skimmed data). "
                          << "passmu4noL1 per-muon flags will be false for all pairs." << std::endl;
            }
        }
    }

    if (use_mu6_for_trg_eff){ // only use mu6 for Pb+Pb 23 for now
        fChainRef()->SetBranchAddress("b_HLT_mu6_L1MU3V"                      , &b_HLT_mu6_L1MU3V);
        fChainRef()->SetBranchAddress("muon_b_HLT_mu6_L1MU3V"                 , &muon_b_HLT_mu6_L1MU3V);
    }

    if (use_mu8_for_trg_eff){ // only use mu6 for Pb+Pb 23 for now
        fChainRef()->SetBranchAddress("b_HLT_mu8_L1MU5VF"                      , &b_HLT_mu8_L1MU5VF);
        fChainRef()->SetBranchAddress("muon_b_HLT_mu8_L1MU5VF"                 , &muon_b_HLT_mu8_L1MU5VF);
    }

    //SetBranch Status
    fChainRef()->SetBranchStatus("*"                               ,0);//switch off all branches, then enable just the ones that we need
    fChainRef()->SetBranchStatus("RunNumber"                       ,1);
    fChainRef()->SetBranchStatus("muon_deltaP_overP"               ,1);
  
    fChainRef()->SetBranchStatus("muon_pair_muon1_index"           ,1);
    fChainRef()->SetBranchStatus("muon_pair_muon1_pt"              ,1);
    fChainRef()->SetBranchStatus("muon_pair_muon1_eta"             ,1);
    fChainRef()->SetBranchStatus("muon_pair_muon1_phi"             ,1);
    fChainRef()->SetBranchStatus("muon_pair_muon1_trk_pt"              ,1);
    fChainRef()->SetBranchStatus("muon_pair_muon1_trk_eta"             ,1);
    fChainRef()->SetBranchStatus("muon_pair_muon1_trk_phi"             ,1);
    fChainRef()->SetBranchStatus("muon_pair_muon1_quality"         ,1);
    fChainRef()->SetBranchStatus("muon_pair_muon1_d0"              ,1);
    fChainRef()->SetBranchStatus("muon_pair_muon1_z0"              ,1);
  
    fChainRef()->SetBranchStatus("muon_pair_muon2_index"           ,1);
    fChainRef()->SetBranchStatus("muon_pair_muon2_pt"              ,1);
    fChainRef()->SetBranchStatus("muon_pair_muon2_eta"             ,1);
    fChainRef()->SetBranchStatus("muon_pair_muon2_phi"             ,1);
    fChainRef()->SetBranchStatus("muon_pair_muon2_trk_pt"              ,1);
    fChainRef()->SetBranchStatus("muon_pair_muon2_trk_eta"             ,1);
    fChainRef()->SetBranchStatus("muon_pair_muon2_trk_phi"             ,1);
    fChainRef()->SetBranchStatus("muon_pair_muon2_quality"         ,1);
    fChainRef()->SetBranchStatus("muon_pair_muon2_d0"              ,1);
    fChainRef()->SetBranchStatus("muon_pair_muon2_z0"              ,1);
  
    fChainRef()->SetBranchStatus(mu4_trigger_branch.c_str()                          ,1);
    fChainRef()->SetBranchStatus(twomu4_trigger_branch.c_str()                       ,1);
    fChainRef()->SetBranchStatus(mu4_trigger_match_branch.c_str()                    ,1);

    if (isRun3 || isPbPb){
        fChainRef()->SetBranchStatus(mu4_mu4noL1_trigger_branch.c_str()              ,1);

        // 2mu4: enable whichever branch is in use
        if (use_mindR_branch_2mu4) {
            fChainRef()->SetBranchStatus((twomu4_trigger_match_branch + mindR_suffix).c_str(), 1);
        } else {
            fChainRef()->SetBranchStatus(twomu4_trigger_match_branch.c_str(), 1);
        }

        // mu4_mu4noL1: enable whichever branch(es) are in use
        if (use_leg_branches_mu4_mu4noL1) {
            fChainRef()->SetBranchStatus((mu4_mu4noL1_trigger_match_branch + "_mu1passLeg1" + leg_dR_suffix).c_str(), 1);
            fChainRef()->SetBranchStatus((mu4_mu4noL1_trigger_match_branch + "_mu1passLeg2" + leg_dR_suffix).c_str(), 1);
            fChainRef()->SetBranchStatus((mu4_mu4noL1_trigger_match_branch + "_mu2passLeg1" + leg_dR_suffix).c_str(), 1);
            fChainRef()->SetBranchStatus((mu4_mu4noL1_trigger_match_branch + "_mu2passLeg2" + leg_dR_suffix).c_str(), 1);
        } else if (use_mindR_only_branch_mu4_mu4noL1) {
            fChainRef()->SetBranchStatus((mu4_mu4noL1_trigger_match_branch + mindR_suffix).c_str(), 1);
        } else {
            fChainRef()->SetBranchStatus(mu4_mu4noL1_trigger_match_branch.c_str()        ,1);
        }
    } else {
        fChainRef()->SetBranchStatus(twomu4_trigger_match_branch.c_str()                 ,1);
    }
    
    if (use_mu6_for_trg_eff){ // only use mu6 for Pb+Pb 23 for now
        fChainRef()->SetBranchStatus("b_HLT_mu6_L1MU3V",         1);
        fChainRef()->SetBranchStatus("muon_b_HLT_mu6_L1MU3V",    1);
    }

    if (use_mu8_for_trg_eff){ // only use mu6 for Pb+Pb 23 for now
        fChainRef()->SetBranchStatus("b_HLT_mu8_L1MU5VF",         1);
        fChainRef()->SetBranchStatus("muon_b_HLT_mu8_L1MU5VF",    1);
    }
}

template <class PairT, class MuonT, class Derived, class... Extras>
void DimuonDataAlgCoreT<PairT, MuonT, Derived, Extras...>::TrigModeToSuffixMap(){
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


template <class PairT, class MuonT, class Derived, class... Extras>
void DimuonDataAlgCoreT<PairT, MuonT, Derived, Extras...>::InitOutputSettings_DataCore() {
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
    std::string mindR_suffix_output = use_mindR_suffix_in_output
        ? ("_mindR_" + (mindR_trig == 0.01 ? std::string("0_01") : std::string("0_02")))
        : "";
    std::string test_suffix = this->is_test_run ? "_test" : "";
    std::string outfile_ending = run_suffix + file_batch_suffix + trig_suffix + mindR_suffix_output + tight_suffix + resonance_cut_suffix + test_suffix + ".root";
    output_file_path = data_dir + file_name_base + outfile_ending;
    output_hist_file_path = data_dir + "hists_cut_acceptance" + outfile_ending;

}

template <class PairT, class MuonT, class Derived, class... Extras>
bool DimuonDataAlgCoreT<PairT, MuonT, Derived, Extras...>::TrigMatch(int pair_ind, int m1_ind, int m2_ind){
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
            if (!mpairRef()) throw std::runtime_error("TrigMatch: Muon Pair doesn't exist!");

            pass_mu6 &= ((muon_b_HLT_mu6_L1MU3V->at(m1_ind) && mpairRef()->m1.pt >= 6) || (muon_b_HLT_mu6_L1MU3V->at(m2_ind) && mpairRef()->m2.pt >= 6));

            if (pass_mu6 || !use_mu8_for_trg_eff) return pass_mu6; // if pass mu4 or not using higher 

            // Case III: Not passing either mu4 or mu6 && using mu8 (mu4 & mu6 prescaled) --> look at mu8
            
            if (!b_HLT_mu8_L1MU5VF) return false;

            if (!muon_b_HLT_mu8_L1MU5VF){ // nullptr
                throw std::runtime_error("TrigMatch:: muon_b_HLT_mu8_L1MU5VF is a null pointer!");
            }

            return ((muon_b_HLT_mu8_L1MU5VF->at(m1_ind) && mpairRef()->m1.pt >= 8) || (muon_b_HLT_mu8_L1MU5VF->at(m2_ind) && mpairRef()->m2.pt >= 8));
        }
		case 2:
            if (!isPbPb && !isRun3){ // pp Run2: no mu4_mu4noL1 trigger
                std::cerr << "PP Run2 data: no mu4_mu4_noL1 trigger!" << std::endl;
                throw std::exception();
            }
			if (!b_HLT_mu4_mu4noL1) return false;

            if (use_leg_branches_mu4_mu4noL1) {
                if (!dimuon_b_mu4_mu4noL1_mu1passLeg1 || !dimuon_b_mu4_mu4noL1_mu2passLeg2 ||
                    !dimuon_b_mu4_mu4noL1_mu2passLeg1 || !dimuon_b_mu4_mu4noL1_mu1passLeg2)
                    throw std::runtime_error("TrigMatch case 2: leg branch pointer(s) are null!");
                bool r1 = dimuon_b_mu4_mu4noL1_mu1passLeg1->at(pair_ind) && dimuon_b_mu4_mu4noL1_mu2passLeg2->at(pair_ind);
                bool r2 = dimuon_b_mu4_mu4noL1_mu2passLeg1->at(pair_ind) && dimuon_b_mu4_mu4noL1_mu1passLeg2->at(pair_ind);
                return r1 || r2;
            } else if (use_mindR_only_branch_mu4_mu4noL1) {
                if (!dimuon_b_mu4_mu4noL1_mindR_only)
                    throw std::runtime_error("TrigMatch case 2: dimuon_b_mu4_mu4noL1_mindR_only is null!");
                return dimuon_b_mu4_mu4noL1_mindR_only->at(pair_ind);
            } else {
                if (!dimuon_b_HLT_mu4_mu4noL1)
                    throw std::runtime_error("TrigMatch case 2: dimuon_b_HLT_mu4_mu4noL1 is null!");
                return dimuon_b_HLT_mu4_mu4noL1->at(pair_ind);
            }
		case 3:
			if (!b_HLT_2mu4) return false;

            if (use_mindR_branch_2mu4) {
                if (!dimuon_b_2mu4_mindR)
                    throw std::runtime_error("TrigMatch case 3: dimuon_b_2mu4_mindR is null!");
                return dimuon_b_2mu4_mindR->at(pair_ind);
            } else {
                if (!dimuon_b_HLT_2mu4)
                    throw std::runtime_error("TrigMatch case 3: dimuon_b_HLT_2mu4 is null!");
                return dimuon_b_HLT_2mu4->at(pair_ind);
            }
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

template <class PairT, class MuonT, class Derived, class... Extras>
void DimuonDataAlgCoreT<PairT, MuonT, Derived, Extras...>::FillMuonPair_DataCore(int pair_ind){
  	if (debug_mode) std::cout << "DEBUG: Calling DimuonDataAlgCoreT<PairT, MuonT, Derived, Extras>::FillMuonPair" << std::endl;
  	
  	mpairRef()->run_number = RunNumber;
  	mpairRef()->lb         = lbn;
  	mpairRef()->bcid       = bcid;

  	mpairRef()->m1.ind     = muon_pair_muon1_index->at(pair_ind);
  	mpairRef()->m2.ind     = muon_pair_muon2_index->at(pair_ind);
	
  	mpairRef()->m1.pt    = fabs(muon_pair_muon1_pt->at(pair_ind))/1000.0;//pt of the first muon in the pair
  	mpairRef()->m2.pt    = fabs(muon_pair_muon2_pt->at(pair_ind))/1000.0;//pt of the second muon in the pair
  	mpairRef()->m1.eta   = muon_pair_muon1_eta->at(pair_ind);
  	mpairRef()->m2.eta   = muon_pair_muon2_eta->at(pair_ind);
  	mpairRef()->m1.phi   = muon_pair_muon1_phi->at(pair_ind);
  	mpairRef()->m2.phi   = muon_pair_muon2_phi->at(pair_ind);

  	mpairRef()->m1.d0    = muon_pair_muon1_d0 ->at(pair_ind);
  	mpairRef()->m2.d0    = muon_pair_muon2_d0 ->at(pair_ind);
  	mpairRef()->m1.z0    = muon_pair_muon1_z0 ->at(pair_ind);
  	mpairRef()->m2.z0    = muon_pair_muon2_z0 ->at(pair_ind);
  	mpairRef()->m1.charge  =(muon_pair_muon1_pt ->at(pair_ind) > 0)? 1:-1;//sign of pt stores charge
  	mpairRef()->m2.charge  =(muon_pair_muon2_pt ->at(pair_ind) > 0)? 1:-1;//sign of pt stores charge
  	mpairRef()->m1.quality = muon_pair_muon1_quality->at(pair_ind);
  	mpairRef()->m2.quality = muon_pair_muon2_quality->at(pair_ind);
  	mpairRef()->m1.dP_overP = muon_deltaP_overP->at(mpairRef()->m1.ind);
  	mpairRef()->m2.dP_overP = muon_deltaP_overP->at(mpairRef()->m2.ind);
	
  	mpairRef()->m1.trk_pt      = fabs(muon_pair_muon1_trk_pt->at(pair_ind))/1000.0;//pt of the first muon in the pair
  	mpairRef()->m2.trk_pt      = fabs(muon_pair_muon2_trk_pt->at(pair_ind))/1000.0;//pt of the second muon in the pair
  	mpairRef()->m1.trk_eta     = muon_pair_muon1_trk_eta->at(pair_ind);
  	mpairRef()->m2.trk_eta     = muon_pair_muon2_trk_eta->at(pair_ind);
  	mpairRef()->m1.trk_phi     = muon_pair_muon1_trk_phi->at(pair_ind);
  	mpairRef()->m2.trk_phi     = muon_pair_muon2_trk_phi->at(pair_ind);

  	if (turn_on_track_charge){
  	    mpairRef()->m1.trk_charge  =(muon_pair_muon1_trk_pt ->at(pair_ind) > 0)? 1:-1;//sign of pt stores charge
  	    mpairRef()->m2.trk_charge  =(muon_pair_muon2_trk_pt ->at(pair_ind) > 0)? 1:-1;//sign of pt stores charge    
  	}else{ // do not turn on track charge: set to nonsense
  	  	mpairRef()->m1.trk_charge  = 0;
  	  	mpairRef()->m2.trk_charge  = 0;
  	}
}


template <class PairT, class MuonT, class Derived, class... Extras>
bool DimuonDataAlgCoreT<PairT, MuonT, Derived, Extras...>::PassCuts_DataCore(bool requireTight){

	//require some quality cuts on the muons
	if((mpairRef()->m1.quality&mpairRef()->m2.quality&1  )==0) return false;//combined muon

	if (requireTight){
		if ((mpairRef()->m1.quality&mpairRef()->m2.quality&16)==0) return false;//tight muon
	}else{
		if ((mpairRef()->m1.quality&mpairRef()->m2.quality&8  )==0) return false;//Medium muon
	}

	if((mpairRef()->m1.quality&mpairRef()->m2.quality&32 )==0) return false;//IDCuts
	if((mpairRef()->m1.quality&mpairRef()->m2.quality&256)==0) return false;//MuonCuts
	h_cutAcceptanceRef()[mpairRef()->m1.charge != mpairRef()->m2.charge]->Fill(static_cast<int>(CutsCommon::pass_muon_quality) + 0.5, mpairRef()->weight); // if same sign: fill the h_cutAcceptanceRef()[0] histogram; if opposite sign, fill the h_cutAcceptanceRef()[1] histogram

	if (fabs(mpairRef()->m1.eta) > 2.4 || fabs(mpairRef()->m2.eta) > 2.4) return false;
	h_cutAcceptanceRef()[mpairRef()->m1.charge != mpairRef()->m2.charge]->Fill(static_cast<int>(CutsCommon::pass_muon_eta) + 0.5, mpairRef()->weight);

	if (mpairRef()->m1.pt < 4 || mpairRef()->m2.pt < 4) return false;
	h_cutAcceptanceRef()[mpairRef()->m1.charge != mpairRef()->m2.charge]->Fill(static_cast<int>(CutsCommon::pass_muon_pt) + 0.5, mpairRef()->weight);
	
	if( fabs(mpairRef()->m1.dP_overP) > pms.deltaP_overP_thrsh || fabs(mpairRef()->m2.dP_overP) > pms.deltaP_overP_thrsh ) return false;
	h_cutAcceptanceRef()[mpairRef()->m1.charge != mpairRef()->m2.charge]->Fill(static_cast<int>(CutsCommon::pass_muon_dP_overP) + 0.5, mpairRef()->weight);
	
	//cut on d0 & z0 sin(theta)
	double m1z0sinTheta = fabs(mpairRef()->m1.z0 * sin(2.0*atan(exp(-mpairRef()->m1.eta))));
	double m2z0sinTheta = fabs(mpairRef()->m2.z0 * sin(2.0*atan(exp(-mpairRef()->m2.eta))));
	bool pass_d0_z0_cuts = (fabs(mpairRef()->m1.d0) < pms.d0cut && fabs(mpairRef()->m2.d0) < pms.d0cut && m1z0sinTheta < pms.z0cut && m2z0sinTheta < pms.z0cut);
	if (!pass_d0_z0_cuts) return false;
	h_cutAcceptanceRef()[mpairRef()->m1.charge != mpairRef()->m2.charge]->Fill(static_cast<int>(CutsCommon::pass_muon_d0_z0) + 0.5, mpairRef()->weight);

	// for both muons, require muon + track charge to agree
	if (turn_on_track_charge){
        if (mpairRef()->m1.trk_charge != mpairRef()->m1.charge || mpairRef()->m2.trk_charge != mpairRef()->m2.charge) return false;
	}
	h_cutAcceptanceRef()[mpairRef()->m1.charge != mpairRef()->m2.charge]->Fill(static_cast<int>(CutsCommon::pass_muon_trk_charge) + 0.5, mpairRef()->weight); // always fill if NOT turn on track charge
	
	return true;
}


template <class PairT, class MuonT, class Derived, class... Extras>
bool DimuonDataAlgCoreT<PairT, MuonT, Derived, Extras...>::IsPhotoProduction(){
    return (!(mpairRef()->same_sign) && mpairRef()->asym < 0.05 && mpairRef()->acop < 0.01);
}


template <class PairT, class MuonT, class Derived, class... Extras>
void DimuonDataAlgCoreT<PairT, MuonT, Derived, Extras...>::ProcessDataHook(){

	nentries = fChainRef()->GetEntries();//number of events
    const Long64_t nentries_to_process = (this->nevents_max <= 0)
        ? nentries
        : std::min<Long64_t>(nentries, this->nevents_max);

    for (Long64_t jentry=0; jentry<nentries_to_process;jentry++) {//loop over the events
	// for (Long64_t jentry=0; jentry<1000;jentry++) {//loop over the events

        if(jentry%100000==0) cout<<"Processing "<<jentry<<" event out of "<<nentries_to_process<<" events"<<std::endl;

		int num_bytes = fChainRef()->GetEntry(jentry);//read in an event
		if(num_bytes==0){
		  	std::cout<<"Error:: Read in event has size of zero bytes,  quitting"<<std::endl;
		  	throw std::exception();
		}
		
		muon_pair_list_cur_event_pre_resonance_cut.clear();
		resonance_tagged_muon_index_list.clear(); // MUST CLEAR for each event!!
        resonance_tagged_muon_index_list_v2.clear(); // MUST CLEAR for each event!!

        if (!this->PassEventSelHook()) continue;

		std::vector<int> muon_index_list = {};
		std::vector<int>::iterator it;

        if (isMinBias){ // MinBias data --> perform offline-single-muon analysis

        }else{ // HardProbe data --> perform offline-dimuon analysis
    		int NPairs=muon_pair_muon1_pt->size();//number of muon pairs in the event

    		for(int pair_ind=0;pair_ind<NPairs;pair_ind++){//first loop over all muon-pairs in the event
                mpairRef() = std::make_shared<pair_t>();

    			self().FillMuonPairHook(pair_ind);
    			mpairRef()->m1.ev_num = jentry;
    			mpairRef()->m2.ev_num = jentry;

    			//------------------------------------------------------------

    			h_cutAcceptanceRef()[mpairRef()->m1.charge != mpairRef()->m2.charge]->Fill(static_cast<int>(CutsCommon::nocut) + 0.5, mpairRef()->weight);

    			//Trigger match for muon pair
    			if(!TrigMatch(pair_ind, mpairRef()->m1.ind, mpairRef()->m2.ind)) continue;
    			h_cutAcceptanceRef()[mpairRef()->m1.charge != mpairRef()->m2.charge]->Fill(static_cast<int>(CutsCommon::pass_trigger_match) + 0.5, mpairRef()->weight);

                if (!this->PassCuts())continue;
    			
    			mpairRef()->pair_pass_tight = (mpairRef()->m1.quality&mpairRef()->m2.quality&16); //tag tight muon pairs
    		
    			//------------------------------------------------------------

    			// Set per-muon passmu4noL1 BEFORE Update() so that the flag follows
    			// the correct muon if pt-ordering swaps m1 and m2.
    			// mu1passLeg2 = mu1 passes the unseeded (noL1) leg (Leg2) of the trigger.
    			if (isRun3 || isPbPb) {
    			    mpairRef()->m1.passmu4noL1 = use_leg_branches_mu4_mu4noL1
    			        ? dimuon_b_mu4_mu4noL1_mu1passLeg2->at(pair_ind) : false;
    			    mpairRef()->m2.passmu4noL1 = use_leg_branches_mu4_mu4noL1
    			        ? dimuon_b_mu4_mu4noL1_mu2passLeg2->at(pair_ind) : false;
    			} else {
    			    mpairRef()->m1.passmu4noL1 = false;
    			    mpairRef()->m2.passmu4noL1 = false;
    			}

    			//Two things at this step: 
    			//1) sort pt, eta, phi by pt (also swaps passmu4noL1 correctly)
    			//2) update the muon-pair values
    			mpairRef()->Update();

    			// pair-level trigger flags (after Update so m1/m2 are pt-sorted)
    			if (!isPbPb && !isRun3) {
    			    mpairRef()->passmu4mu4noL1 = false;
    			    mpairRef()->passmu4noL1    = false;
    			} else if (use_leg_branches_mu4_mu4noL1) {
    			    bool r1 = dimuon_b_mu4_mu4noL1_mu1passLeg1->at(pair_ind) && dimuon_b_mu4_mu4noL1_mu2passLeg2->at(pair_ind);
    			    bool r2 = dimuon_b_mu4_mu4noL1_mu2passLeg1->at(pair_ind) && dimuon_b_mu4_mu4noL1_mu1passLeg2->at(pair_ind);
    			    mpairRef()->passmu4mu4noL1 = r1 || r2;
    			    mpairRef()->passmu4noL1    = mpairRef()->m1.passmu4noL1 || mpairRef()->m2.passmu4noL1;
    			} else if (use_mindR_only_branch_mu4_mu4noL1) {
    			    mpairRef()->passmu4mu4noL1 = dimuon_b_mu4_mu4noL1_mindR_only->at(pair_ind);
    			    mpairRef()->passmu4noL1    = false; // no per-leg info available
    			} else {
    			    mpairRef()->passmu4mu4noL1 = dimuon_b_HLT_mu4_mu4noL1->at(pair_ind);
    			    mpairRef()->passmu4noL1    = false; // fallback: no per-leg info
    			}

    			if (use_mindR_branch_2mu4) {
    			    mpairRef()->pass2mu4 = dimuon_b_2mu4_mindR->at(pair_ind);
    			} else {
    			    mpairRef()->pass2mu4 = dimuon_b_HLT_2mu4->at(pair_ind);
    			}
    	
        		mpairRef()->m1.passmu4 = muon_b_HLT_mu4->at(mpairRef()->m1.ind);
                if (use_mu6_for_trg_eff) mpairRef()->m1.passmu4 |= (muon_b_HLT_mu6_L1MU3V->at(mpairRef()->m1.ind) && mpairRef()->m1.pt > 6);
                if (use_mu8_for_trg_eff) mpairRef()->m1.passmu4 |= (muon_b_HLT_mu8_L1MU5VF->at(mpairRef()->m1.ind) && mpairRef()->m1.pt > 8);
         
                mpairRef()->m2.passmu4 = muon_b_HLT_mu4->at(mpairRef()->m2.ind);
                if (use_mu6_for_trg_eff) mpairRef()->m2.passmu4 |= (muon_b_HLT_mu6_L1MU3V->at(mpairRef()->m2.ind) && mpairRef()->m2.pt > 6);
                if (use_mu8_for_trg_eff) mpairRef()->m2.passmu4 |= (muon_b_HLT_mu8_L1MU5VF->at(mpairRef()->m2.ind) && mpairRef()->m2.pt > 8);
                

                mpairRef()->passSeparated = (mpairRef()->dr > 0.8);
                mpairRef()->passSeparatedDeta = (abs(mpairRef()->deta) > 0.8);

    			// resonance tag
    			ResonanceTagging();
                ResonanceTaggingV2();

    			// photo-production cut
    			if (isPbPb){
              	if (!(trigger_effcy_calc && !filter_out_photo_resn_for_trig_effcy) && IsPhotoProduction()) continue;
    			  	h_cutAcceptanceRef()[mpairRef()->m1.charge != mpairRef()->m2.charge]->Fill(static_cast<int>(CutsPbPb::pass_photoprod) + 0.5, mpairRef()->weight);      	
    			}

    			muon_pair_list_cur_event_pre_resonance_cut.push_back(std::move(mpairRef()));

    		} // finish first loop over all muon pairs

    		for(int pair_ind = 0; pair_ind < muon_pair_list_cur_event_pre_resonance_cut.size(); pair_ind++){//second loop over all muon-pairs in the event
    			// discard the pair if either muon is resonance-tagged

    			mpairRef() = std::move(muon_pair_list_cur_event_pre_resonance_cut.at(pair_ind));

    			if (!mpairRef()){
    			  	std::cerr << "mpairRef() at second muon-pair loop NOT found! Return without resonance-tag checking or pair analysis." << std::endl;
    			  	continue;
    			}

    			std::vector<int>::iterator itres_m1;
    			std::vector<int>::iterator itres_m2;

                // apply resonance cuts if resonance_cut_mode != 0
                if (!(trigger_effcy_calc && !filter_out_photo_resn_for_trig_effcy)){ // perform resn cuts
                    if (resonance_cut_mode == 1){ // apply old cuts
                        itres_m1 = std::find(resonance_tagged_muon_index_list.begin(),resonance_tagged_muon_index_list.end(),mpairRef()->m1.ind);
                        if(itres_m1 != resonance_tagged_muon_index_list.end())  continue;

                        itres_m2 = std::find(resonance_tagged_muon_index_list.begin(),resonance_tagged_muon_index_list.end(),mpairRef()->m2.ind);
                        if(itres_m2 != resonance_tagged_muon_index_list.end())  continue;
                    } else if (resonance_cut_mode == 2){ // apply new cuts
                        itres_m1 = std::find(resonance_tagged_muon_index_list_v2.begin(),resonance_tagged_muon_index_list_v2.end(),mpairRef()->m1.ind);
                        if(itres_m1 != resonance_tagged_muon_index_list_v2.end())  continue;

                        itres_m2 = std::find(resonance_tagged_muon_index_list_v2.begin(),resonance_tagged_muon_index_list_v2.end(),mpairRef()->m2.ind);
                        if(itres_m2 != resonance_tagged_muon_index_list_v2.end())  continue;
                    }                    
                }

    			int pass_resonance_ind = isPbPb? static_cast<int>(CutsPbPb::pass_resonance) : static_cast<int>(CutsPP::pass_resonance);
    			h_cutAcceptanceRef()[mpairRef()->m1.charge != mpairRef()->m2.charge]->Fill(pass_resonance_ind + 0.5, mpairRef()->weight);

    			//------------------------------------------------------------
    			// perform additional pair-level analysis if needed
                PerformAdditionalPairAnalysisHook();

    			if(output_single_muon_tree){
    			  it = std::find(muon_index_list.begin(),muon_index_list.end(),mpairRef()->m1.ind);
    			  if(it == muon_index_list.end()){ //muon1 index NOT found
    			    	muon_index_list.push_back(mpairRef()->m1.ind);
    			    	muon_raw_ptr = &(mpairRef()->m1);
    			    	FillSingleMuonTree();
    			  }
    			  it = std::find(muon_index_list.begin(),muon_index_list.end(),mpairRef()->m2.ind);
    			  if(it == muon_index_list.end()){ //muon1 index NOT found
    			    	muon_index_list.push_back(mpairRef()->m2.ind);
    			    	muon_raw_ptr = &(mpairRef()->m2);
    			    	FillSingleMuonTree();
    			  }
    			}
    			else{ // fill muon pair trees
    			 	FillMuonPairTree();
    			}
    		}

            mpairRef().reset(); // make sure the reference to the last muon pair gets reset
        }
	}//loop over events
}


