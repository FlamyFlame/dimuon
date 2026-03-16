#ifndef __PythiaTreeParser_C__
#define __PythiaTreeParser_C__

#include "PythiaNTupleFirstPass.h"
#include "Riostream.h"
#include "TLorentzVector.h"
#include <math.h>

// #include<bits/stdc++.h>

PythiaNTupleFirstPass::PythiaNTupleFirstPass(){
  cutLabels = cutLabels_MC;
  numCuts = cutLabels_MC.size();

  std::cout << "Pythia Ntuple processing script:" << std::endl;
  std::cout << std::endl;

  std::cout << "The following public variable(s) **MUST** be set:" << std::endl;
  std::cout << "--> batch_num: Integer with value 1 or 2" << std::endl;
  std::cout << "  * batch_num = 1: runs pythia samples run on or before 03/18/2023 (1st batch)" << std::endl;
  std::cout << "  * batch_num = 2: runs pythia samples run on or after 03/22/2023 (2nd batch)" << std::endl;
  std::cout << std::endl;

  std::cout << "The following public variable(s) should be checked:" << std::endl;
  std::cout << "--> turn_data_resonance_cuts_on:  boolean to turn on/off minv-based resonance cuts (same as in data)" << std::endl;
  std::cout << "                                  default false" << std::endl;
  std::cout << std::endl;

  // std::cout << "The following public variable(s) can be (optionally) set:" << std::endl;
  // std::cout << "print_prt_history:            print parent history for muon pairs from open HF's" << std::endl;
  // std::cout << "                              default false" << std::endl;
  // std::cout << "print_HF_pair_origin_others_history:         print parent history for muon pairs with origin tag \"others\"" << std::endl;
  // std::cout << "                              default false" << std::endl;
  // std::cout << "print_unspecified_parent:     print history in the middle of origin tracing for various places the parents are \"unspecified\"" << std::endl;
  // std::cout << "                              default false" << std::endl;
  // std::cout << "print_FE:                     print history for pairs from flavor excitation" << std::endl;
  // std::cout << "                              default false" << std::endl;
  // std::cout << "print_bad_warnings:           print bad warnings" << std::endl;
  // std::cout << "                              default false" << std::endl;
  // std::cout << "print_low_minv_resonances:          print history for pairs from resonances with very low invariant mass (< low_minv_threshold)" << std::endl;
  // std::cout << "                              default false" << std::endl;
  // std::cout << std::endl;

  std::cout << "Output files will be written to /usatlas/u/yuhanguo/usatlasdata/pythia" << std::endl;
  std::cout << "    Including a root file for muon pairs, a root file for histograms, and output txt files" << std::endl;
}

void PythiaNTupleFirstPass::InputSanityCheck(){
  if (job_dirs.size() != kn_in_job.size() || job_dirs.size() != nfiles_factor.size()){
    std::cout << "The vectors job_dirs, kn_in_job and nfiles_factor must have the same length!" << std::endl;
    std::cout << "The length should be the number of runs (directories) to include." << std::endl;
    throw std::exception();
  }
  for (auto v : kn_in_job){
    if (v.size() != nKinRanges){
      std::cout << "Each vector in kn_in_job must have length 5 (mKinRange)!" << std::endl;
      throw std::exception();
    }
  }
    for (auto v : nfiles_factor){
    if (v.size() != nKinRanges){
      std::cout << "Each vector in nfiles_factor must have length 5 (mKinRange)!" << std::endl;
      throw std::exception();
    }
  }
}

void PythiaNTupleFirstPass::SetInputOutputFilesFromBatch(){
  switch(batch_num){
  case 1:
    new_run = false;
    batch_suffix = "_allto0318";
    outfile_name = "muon_pairs_pythia_allto0318";
    outhistfile_name = "hists_pythia_ntuple_processing_allto0318";
    job_dirs = {"0317_all_k/", "0318_all_k/", "0318_k0/"};
    kn_in_job = {{true, true, true, true, true},{true, true, true, true, true},{true,false,false,false,false}};
    nfiles_factor = {{20, 10, 4, 1, 1}, {20, 10, 4, 1, 1},{40,0,0,0,0}};
    break;

  case 2:
    new_run = true;
    batch_suffix = "_after0322";
    outfile_name = "muon_pairs_pythia_after0322";
    outhistfile_name = "hists_pythia_ntuple_processing_after0322";
    job_dirs = {"0322_k0_k1/", "0323_k0_k1/", "0325_all_k/", "0401_all_k/", "0429_all_k/"};
    kn_in_job = {{true,true,false,false,false},{true,true,false,false,false},{true, true, true, true, true},{true, true, true, true, true},{true, true, true, true, true}};
    nfiles_factor = {{80,20,0,0,0},{80,20,0,0,0},{80,40,4,1,1},{80,40,4,1,1},{80,40,4,1,1}};
    break;

  default:
    std::cerr << "ERROR: The integer batch number must be specified and set between 1 and 2!" << std::endl;
    throw std::exception();
  }
}

void PythiaNTupleFirstPass::ResonanceNameMap(){
  // reference: https://pdg.lbl.gov/2020/reviews/rpp2020-rev-monte-carlo-numbering.pdf
  // besides X --> mu+ mu-, we added X \in {eta, eta'}, which decays mainly into mu+ mu- photon

  // resonance_ids = {113, 223, 333, 443, 553}; // we do %10000 in tracing --> do not list higher resonances
  resonance_ids = {113, 221, 223, 331, 333, 443, 553}; // we do %10000 in tracing --> do not list higher resonances
  resonance_id_to_name_and_crossx_map[113] = {"rho0", 0.};
  resonance_id_to_name_and_crossx_map[221] = {"eta", 0.};
  resonance_id_to_name_and_crossx_map[223] = {"omega", 0.};
  resonance_id_to_name_and_crossx_map[331] = {"eta'", 0.};
  resonance_id_to_name_and_crossx_map[333] = {"phi", 0.};
  resonance_id_to_name_and_crossx_map[443] = {"J/psi", 0.};
  resonance_id_to_name_and_crossx_map[553] = {"Upsilon", 0.};
}


//initialize the TChain
void PythiaNTupleFirstPass::InitInput(){
    clock_t start, end;
    double cpu_time_used;

    evChain = new TChain("PyTree","PyTree");
    metaChain = new TChain("meta_tree","meta_tree");
    evChain->SetMakeClass(1);
    metaChain->SetMakeClass(1);

  for (int ikin = 0; ikin < nKinRanges; ikin++){
    std::cout << "Currently loading root files for the " << ikin << "-th kinematic range into the TChains." << std::endl;

    // for (int ikin = 3; ikin <= 3; ikin++){
    // evChain.push_back(new TChain("PyTree","PyTree"));
      // metaChain.push_back(new TChain("meta_tree","meta_tree"));

      for (int jjob = 0; jjob < job_dirs.size(); jjob++){
            std::cout << "For kin-range" << ikin << ", currently loading root files from the directory " << job_dirs[jjob] << "." << std::endl;
            start = clock();

            if (kn_in_job[jjob][ikin]){ // we have run the i-th kinematic range in the j-th job

          for (int kbeam = 0; kbeam < nBeamTypes; kbeam++){
            std::string job_path = py_dir + job_dirs[jjob] + kin_dirs[ikin] + beam_dirs[kbeam];
            int nfiles = nfiles_base[kbeam] * nfiles_factor[jjob][ikin];
            njobs[ikin] += nfiles;
            nevents[ikin] += nfiles * nevents_per_file[ikin];
            // std::cout << "Dir: " << job_dirs[jjob] << kin_dirs[ikin] << beam_dirs[kbeam] << std::endl;
            // std::cout << "Number of root files: " << nfiles << std::endl;

            for (int lfile = 1; lfile <= nfiles; lfile++){ // note: pytree_N.root starts with 1 not 0
                        // const char* filename = (job_path + "pytree_" + std::to_string(lfile) + ".root").c_str();
                        std::ifstream infile((job_path + "pytree_" + std::to_string(lfile) + ".root").c_str());
                        if (!infile.good()){
                            std::cout << "job path " << job_path << ".\n";
                            // std::cout << "file name " << filename << ".\n";
                            std::cout << "file name should be " << (job_path + "pytree_" + std::to_string(lfile) + ".root").c_str() << ".\n";
                            std::cout << "Warning: File " << (job_path + "pytree_" + std::to_string(lfile) + ".root").c_str() << " not found. Skip.\n";
                            njobs[ikin] -= 1;
                            nevents[ikin] -= nevents_per_file[ikin];
                            continue; // skip this file
                        }
                        evChain->Add((job_path + "pytree_" + std::to_string(lfile) + ".root?#PyTree").c_str());
                        metaChain->Add((job_path + "pytree_" + std::to_string(lfile) + ".root?#meta_tree").c_str());
                        // evChain->Add((job_path + "pytree_" + std::to_string(lfile) + ".root?#PyTree").c_str(), 0);
                        // metaChain->Add((job_path + "pytree_" + std::to_string(lfile) + ".root?#meta_tree").c_str(), 0);
                    }
          }
        }

            end = clock();
            cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
            std::cout << "Took " << cpu_time_used << " seconds to load the events in kin-range" << ikin << ", job directory " << job_dirs[jjob] << std::endl;
      }

        std::cout << "Finished loading all roots files from all job directories for the " << ikin << "-th kinematic range." << std::endl;

        nentries_k0 = nevents[0];
        nentries_k1 = nevents[1];
        nentries_k2 = nevents[2];
        nentries_k3 = nevents[3];
        nentries_k4 = nevents[4];

        std::cout << "#events in kin-range k0: " << nevents[0] << std::endl;
        std::cout << "#events in kin-range k1: " << nevents[1] << std::endl;
        std::cout << "#events in kin-range k2: " << nevents[2] << std::endl;
        std::cout << "#events in kin-range k3: " << nevents[3] << std::endl;
        std::cout << "#events in kin-range k4: " << nevents[4] << std::endl;
        // std::cout << "Before performming GetEntries() for the two TChains." << std::endl;
        // start = clock();

        nevents_accum[ikin] = (ikin == 0)? nevents[ikin] : nevents_accum[ikin - 1] + nevents[ikin];
        njobs_accum[ikin] = (ikin == 0)? njobs[ikin] : njobs_accum[ikin - 1] + njobs[ikin];

        // end = clock();
        // cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
        // std::cout << "Took " << cpu_time_used << " seconds to perform GetEntries() for the two TChains." << std::endl;

        // int nevents_prev = (ikin == 0)? 0 : nevents_accum[ikin - 1];
        // int njobs_prev = (ikin == 0)? 0 : njobs_accum[ikin - 1];

        std::cout << "Number of events in kinematic range " << ikin << ": " << nevents[ikin] << std::endl;
        std::cout << "Number of jobs (files) in kinematic range " << ikin << ": " << njobs[ikin] << std::endl;
      // if (njobs[ikin] != njobs_accum[ikin] - njobs_prev){
    //  std::cout << "The " << ikin << "-th kinematic-range contributions to metaChain doesn't have the correct number of entries" << std::endl;
    //  std::cout << "The correct number of entries should be " << njobs[ikin] << std::endl;
    //  std::cout << "The actual number of entries is " << njobs_accum[ikin] - njobs_prev << std::endl;
    //  throw std::exception();
      // }
      // if (nevents[ikin] != nevents_accum[ikin] - nevents_prev){
    //  std::cout << "The " << ikin << "-th kinematic-range contributions to evChain doesn't have the correct number of entries" << std::endl;
    //  std::cout << "The correct number of entries should be " << nevents[ikin] << std::endl;
    //  std::cout << "The actual number of entries is " << nevents_accum[ikin] - nevents_prev << std::endl;
    //  throw std::exception();
      // }

      metaChain->SetBranchAddress("efficiency"                , &efficiency);
        if (new_run){
            metaChain->SetBranchAddress("eventWeight"                , &ev_weight);
        }else{
            metaChain->SetBranchAddress("totalSigma"                , &ev_weight);
        }

        metaChain->SetBranchStatus("*"                          ,0);
        metaChain->SetBranchStatus("efficiency"                 ,1);
        if (new_run){
            metaChain->SetBranchStatus("eventWeight"                ,1);
        }else{
            metaChain->SetBranchStatus("totalSigma"                ,1);
        }

        // metaChain->SetBranchStatus("eventWeight"                ,1);

        evChain->SetBranchAddress("QHard"                     , &QHard);
        evChain->SetBranchAddress("pTHat"                   , &pTHat);
        evChain->SetBranchAddress("mHat"                    , &mHat);
        evChain->SetBranchAddress("truth_id"                  , &truth_id);
        evChain->SetBranchAddress("truth_barcode"               , &truth_barcode);
        evChain->SetBranchAddress("truth_status"               , &truth_status);
        // evChain->SetBranchAddress("truth_qual"                , &truth_qual);
        evChain->SetBranchAddress("truth_m"                   , &truth_m);
        evChain->SetBranchAddress("truth_pt"                    , &truth_pt);
        evChain->SetBranchAddress("truth_eta"                   , &truth_eta);
        evChain->SetBranchAddress("truth_phi"                   , &truth_phi);
        evChain->SetBranchAddress("truth_mother1"               , &truth_mother1);
        evChain->SetBranchAddress("truth_mother2"               , &truth_mother2);
        // if (new_run){
        //     evChain->SetBranchAddress("truth_daughter1"               , &truth_daughter1);
        //     evChain->SetBranchAddress("truth_daughter2"               , &truth_daughter2);
        // }

        evChain->SetBranchAddress("truth_mupair_pt1"           , &muon_pair_muon1_pt);
        evChain->SetBranchAddress("truth_mupair_eta1"          , &muon_pair_muon1_eta);
        evChain->SetBranchAddress("truth_mupair_phi1"          , &muon_pair_muon1_phi);
        evChain->SetBranchAddress("truth_mupair_ch1"           , &muon_pair_muon1_ch);
        evChain->SetBranchAddress("truth_mupair_bar1"          , &muon_pair_muon1_bar);

        evChain->SetBranchAddress("truth_mupair_pt2"           , &muon_pair_muon2_pt);
        evChain->SetBranchAddress("truth_mupair_eta2"          , &muon_pair_muon2_eta);
        evChain->SetBranchAddress("truth_mupair_phi2"          , &muon_pair_muon2_phi);
        evChain->SetBranchAddress("truth_mupair_ch2"           , &muon_pair_muon2_ch);
        evChain->SetBranchAddress("truth_mupair_bar2"          , &muon_pair_muon2_bar);

        //SetBranch Status
        evChain->SetBranchStatus("*"                               ,0);//switch off all branches, then enable just the ones that we need
      
        evChain->SetBranchStatus("QHard"              ,1);
        evChain->SetBranchStatus("pTHat"                ,1);
        evChain->SetBranchStatus("mHat"                 ,1);
        evChain->SetBranchStatus("truth_id"             ,1);
        evChain->SetBranchStatus("truth_barcode"         ,1);
        evChain->SetBranchStatus("truth_status"         ,1);
        // evChain->SetBranchStatus("truth_qual"             ,1);
        evChain->SetBranchStatus("truth_m"             ,1);
        evChain->SetBranchStatus("truth_pt"             ,1);
        evChain->SetBranchStatus("truth_eta"             ,1);
        evChain->SetBranchStatus("truth_phi"             ,1);
        evChain->SetBranchStatus("truth_mother1"             ,1);
        evChain->SetBranchStatus("truth_mother2"             ,1);
        // if (new_run){
        //     evChain->SetBranchStatus("truth_daughter1"             ,1);
        //     evChain->SetBranchStatus("truth_daughter2"             ,1);
        // }

        evChain->SetBranchStatus("truth_mupair_pt1"           ,1);
        evChain->SetBranchStatus("truth_mupair_eta1"              ,1);
        evChain->SetBranchStatus("truth_mupair_phi1"             ,1);
        evChain->SetBranchStatus("truth_mupair_ch1"             ,1);
        evChain->SetBranchStatus("truth_mupair_bar1"             ,1);

        evChain->SetBranchStatus("truth_mupair_pt2"             ,1);
        evChain->SetBranchStatus("truth_mupair_eta2"         ,1);
        evChain->SetBranchStatus("truth_mupair_phi2"              ,1);
        evChain->SetBranchStatus("truth_mupair_ch2"              ,1);
        evChain->SetBranchStatus("truth_mupair_bar2"              ,1);

    }
}

void PythiaNTupleFirstPass::InitTempVariables(){
    m1_history = new std::vector<std::vector<int>>();
    m2_history = new std::vector<std::vector<int>>();
    m1_history_particle = new std::vector<std::vector<Particle>>();
    m2_history_particle = new std::vector<std::vector<Particle>>();
    // m1_single_gluon_history = new std::vector<std::vector<int>>();
    // m2_single_gluon_history = new std::vector<std::vector<int>>();
    // m1_single_gluon_history_particle = new std::vector<std::vector<Particle>>();
    // m2_single_gluon_history_particle = new std::vector<std::vector<Particle>>();
    // m1_last_hf_hadron_prt_pt_eta_phi_m = new std::vector<float>();
    // m2_last_hf_hadron_prt_pt_eta_phi_m = new std::vector<float>();
    m1_first_hf_hadron_prt_pt_eta_phi_m = new std::vector<float>();
    m2_first_hf_hadron_prt_pt_eta_phi_m = new std::vector<float>();
    // m1_hq_ancestor_pt_eta_phi_m = new std::vector<float>();
    // m2_hq_ancestor_pt_eta_phi_m = new std::vector<float>();
}

void PythiaNTupleFirstPass::InitOutput(){

    std::string sign_labels[ParamsSet::nSigns] = {"same_sign", "op_sign"};
    std::string signs[ParamsSet::nSigns] = {"_ss", "_op"};
    std::string dphis[ParamsSet::nSigns] = {"_near", "_away"};
    std::string ancestor_grps[nAncestorGroups] = {"_gg", "_qg","_single_g","_qq"};

    // ---------------------------------------------------------------------------------------------------------------------------

    std::string py_dir_output_diagnostic_files = py_dir + "ancestor_tracing_output_diagnostic_files/";

    if (print_bad_warnings){
      m_very_bad_warning_file = new std::ofstream(py_dir_output_diagnostic_files + "very_bad_warnings_pythia" + batch_suffix + ".txt");
      *m_very_bad_warning_file << "Test writing into warning files [begining]" << std::endl;
      m_hard_scattering_warning_file = new std::ofstream(py_dir_output_diagnostic_files + "hard_scattering_warnings_pythia" + batch_suffix + ".txt");      
    }
    m_crossx_summary_file = new std::ofstream(py_dir_output_diagnostic_files + "crossx_summary_pythia" + batch_suffix + ".txt");

    if (print_low_minv_resonances){
      m_very_low_minv_resonance_file = new std::ofstream(py_dir_output_diagnostic_files + "very_low_minv_resonances_pythia" + batch_suffix + ".txt");
    }

    if (print_unspecified_parent){
      m_unspecified_parent_file = new std::ofstream(py_dir_output_diagnostic_files + "unspecified_parents_pythia" + batch_suffix + ".txt");
    }

    if (print_FE){
      m_FE_file = new std::ofstream(py_dir_output_diagnostic_files + "FE_pythia" + batch_suffix + ".txt");
    }
    
    if (print_HF_pair_origin_others_history){ // debug mode for pairs in the "others" category
      m_HF_pair_origin_others_category_file = new std::ofstream(py_dir_output_diagnostic_files + "HF_pair_origin_others_category" + batch_suffix + ".txt");
    }

    if (print_other_flavor_history){ // debug mode for pairs in the "others" category
      m_other_flavor_category_file = new std::ofstream(py_dir_output_diagnostic_files + "other_flavor_category" + batch_suffix + ".txt");
    }

    if (print_prt_history){
        for (int isign = 0; isign < ParamsSet::nSigns; isign++){
            for (int jdphi = 0; jdphi < 2; jdphi++){
                m_b_parent_file[isign][jdphi] = new std::ofstream(Form("%sb_parents_%s%s%s.txt", py_dir_output_diagnostic_files.c_str(), sign_labels[isign].c_str(), dphis[jdphi].c_str(), batch_suffix.c_str()));
            }
        }
    }
    if (print_prt_history){
        for (int isign = 0; isign < ParamsSet::nSigns; isign++){
            for (int jdphi = 0; jdphi < 2; jdphi++){
                m_c_parent_file[isign][jdphi] = new std::ofstream(Form("%sc_parents_%s%s%s.txt", py_dir_output_diagnostic_files.c_str(), sign_labels[isign].c_str(), dphis[jdphi].c_str(), batch_suffix.c_str()));
            }
        }
    }

    // ---------------------------------------------------------------------------------------------------------------------------
    // --------------------------------------------- output file for tree writing ------------------------------------------------

    std::string apply_data_cuts_suffix = turn_data_resonance_cuts_on? "_with_data_resonance_cuts" : "_no_data_resonance_cuts";
    output_file_path = py_dir + outfile_name + apply_data_cuts_suffix + ".root";
    output_hist_file_path = py_dir + outhistfile_name + apply_data_cuts_suffix + ".root";
    
    m_outfile=new TFile(output_file_path.c_str(), "recreate");
    
    meta_tree_out = new TTree("meta_tree_out","meta_tree_out");
    meta_tree_out->Branch("nentries_k0_with_3.7GeV_cuts", &nentries_k0, "nentries_k0_with_3.7GeV_cuts/i");
    meta_tree_out->Branch("nentries_k1_with_3.7GeV_cuts", &nentries_k1, "nentries_k1_with_3.7GeV_cuts/i");
    meta_tree_out->Branch("nentries_k2_with_3.7GeV_cuts", &nentries_k2, "nentries_k2_with_3.7GeV_cuts/i");
    meta_tree_out->Branch("nentries_k3_with_3.7GeV_cuts", &nentries_k3, "nentries_k3_with_3.7GeV_cuts/i");
    meta_tree_out->Branch("nentries_k4_with_3.7GeV_cuts", &nentries_k4, "nentries_k4_with_3.7GeV_cuts/i");

    for (unsigned int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
        muonPairOutTree[ksign] = new TTree(Form("muon_pair_tree_sign%u",ksign+1),Form("all muon pairs, sign%u",ksign+1));
        muonPairOutTree[ksign]->Branch("MuonPairObj",&mpair_raw_ptr);
        for (int ikin = 0; ikin < nKinRanges; ikin++){
            muonPairOutTreeKinRange[ikin][ksign] = new TTree(Form("muon_pair_tree_kin%d_sign%u",ikin,ksign+1),Form("all muon pairs, kin range%u, sign%u",ikin,ksign+1));
            muonPairOutTreeKinRange[ikin][ksign]->Branch("MuonPairObj",&mpair_raw_ptr);
        }
    }

    // ---------------------------------------------------------------------------------------------------------------------------
    // --------------------------------------------- output file for hist writing ------------------------------------------------

    // for (unsigned int isign = 0; isign < ParamsSet::nSigns; isign++){
    //     for (int jdphi = 0; jdphi < 2; jdphi++){
    //         for (int kgrp = 0; kgrp < nAncestorGroups; kgrp++){
    //             QQPairOutTree[isign][jdphi][kgrp] = new TTree(Form("QQ_pair_tree%s%s%s",signs[isign].c_str(),dphis[jdphi].c_str(),ancestor_grps[kgrp].c_str()),Form("QQ pairs, muon %s%s%s",signs[isign].c_str(),dphis[jdphi].c_str(),ancestor_grps[kgrp].c_str()));
    //             QQPairOutTree[isign][jdphi][kgrp]->Branch("QQPairObj",&qqpair);
    //         }
    //     }
    // }

    m_outHistFile=new TFile(output_hist_file_path.c_str(), "recreate");

    h_numMuonPairs = new TH1D("h_numMuonPairs","h_numMuonPairs",6,0,6);

    DimuonAnalysisBaseClass::InitOutput();

    for (int isign = 0; isign < ParamsSet::nSigns; isign++){
        for (int jdphi = 0; jdphi < 2; jdphi++){
            // h_unweighted_parent_groups[isign][jdphi] = new TH2D(Form("h_unweighted_parent_groups_sign%d%s",isign+1, dphis[jdphi].c_str()),Form("h_unweighted_parent_groups_sign%d%s",isign+1, dphis[jdphi].c_str()),nParentGroups,0,nParentGroups,nParentGroups,0,nParentGroups);
            h_parent_groups[isign][jdphi] = new TH2D(Form("h_parent_groups_sign%d%s",isign+1, dphis[jdphi].c_str()),Form("h_parent_groups_sign%d%s",isign+1, dphis[jdphi].c_str()),nParentGroups,0,nParentGroups,nParentGroups,0,nParentGroups);
            
            h_hard_scatt_categr_vs_origin_categr[isign][jdphi] = new TH2D(Form("h_hard_scatt_categr_vs_origin_categr_sign%d%s",isign+1, dphis[jdphi].c_str()),Form("h_hard_scatt_categr_vs_origin_categr_sign%d%s",isign+1, dphis[jdphi].c_str()),nHFMuonPairOriginCategories,0,nHFMuonPairOriginCategories,nHardScattCategories,0,nHardScattCategories);
            h_both_from_b_hard_scatt_categr_vs_origin_categr[isign][jdphi] = new TH2D(Form("h_both_from_b_hard_scatt_categr_vs_origin_categr_sign%d%s",isign+1, dphis[jdphi].c_str()),Form("h_both_from_b_hard_scatt_categr_vs_origin_categr_sign%d%s",isign+1, dphis[jdphi].c_str()),nHFMuonPairOriginCategories,0,nHFMuonPairOriginCategories,nHardScattCategories,0,nHardScattCategories);
            h_both_from_c_hard_scatt_categr_vs_origin_categr[isign][jdphi] = new TH2D(Form("h_both_from_c_hard_scatt_categr_vs_origin_categr_sign%d%s",isign+1, dphis[jdphi].c_str()),Form("h_both_from_c_hard_scatt_categr_vs_origin_categr_sign%d%s",isign+1, dphis[jdphi].c_str()),nHFMuonPairOriginCategories,0,nHFMuonPairOriginCategories,nHardScattCategories,0,nHardScattCategories);
            h_one_from_b_one_from_c_hard_scatt_categr_vs_origin_categr[isign][jdphi] = new TH2D(Form("h_one_from_b_one_from_c_hard_scatt_categr_vs_origin_categr_sign%d%s",isign+1, dphis[jdphi].c_str()),Form("h_one_from_b_one_from_c_hard_scatt_categr_vs_origin_categr_sign%d%s",isign+1, dphis[jdphi].c_str()),nHFMuonPairOriginCategories,0,nHFMuonPairOriginCategories,nHardScattCategories,0,nHardScattCategories);
            
            h_muon_pair_origin_categr[isign][jdphi] = new TH1D(Form("h_muon_pair_origin_categr_sign%d%s",isign+1, dphis[jdphi].c_str()),Form("h_muon_pair_origin_categr_sign%d%s",isign+1, dphis[jdphi].c_str()),nHFMuonPairOriginCategories,0,nHFMuonPairOriginCategories);
            h_both_from_b_muon_pair_origin_categr[isign][jdphi] = new TH1D(Form("h_both_from_b_muon_pair_origin_categr_sign%d%s",isign+1, dphis[jdphi].c_str()),Form("h_both_from_b_muon_pair_origin_categr_sign%d%s",isign+1, dphis[jdphi].c_str()),nHFMuonPairOriginCategories,0,nHFMuonPairOriginCategories);
            h_both_from_c_muon_pair_origin_categr[isign][jdphi] = new TH1D(Form("h_both_from_c_muon_pair_origin_categr_sign%d%s",isign+1, dphis[jdphi].c_str()),Form("h_both_from_c_muon_pair_origin_categr_sign%d%s",isign+1, dphis[jdphi].c_str()),nHFMuonPairOriginCategories,0,nHFMuonPairOriginCategories);
            h_one_from_b_one_from_c_muon_pair_origin_categr[isign][jdphi] = new TH1D(Form("h_one_from_b_one_from_c_muon_pair_origin_categr_sign%d%s",isign+1, dphis[jdphi].c_str()),Form("h_one_from_b_one_from_c_muon_pair_origin_categr_sign%d%s",isign+1, dphis[jdphi].c_str()),nHFMuonPairOriginCategories,0,nHFMuonPairOriginCategories);
            
            // bins to be changed later --> get a crude feeling for now
            h_Qsplit_gs_ISR_one_hard_scatt[isign][jdphi] = new TH1D(Form("h_Qsplit_gs_ISR_one_hard_scatt_sign%d%s",isign+1, dphis[jdphi].c_str()),Form("h_Qsplit_gs_ISR_one_hard_scatt_sign%d%s",isign+1, dphis[jdphi].c_str()),80,0,80);
            h_Qsplit_gs_FSR[isign][jdphi] = new TH1D(Form("h_Qsplit_gs_FSR_sign%d%s",isign+1, dphis[jdphi].c_str()),Form("h_Qsplit_gs_FSR_sign%d%s",isign+1, dphis[jdphi].c_str()),80,0,80);
            h_Qsplit_to_mHard_gs_ISR_one_hard_scatt[isign][jdphi] = new TH1D(Form("h_Qsplit_to_mHard_gs_ISR_one_hard_scatt_sign%d%s",isign+1, dphis[jdphi].c_str()),Form("h_Qsplit_to_mHard_gs_ISR_one_hard_scatt_sign%d%s",isign+1, dphis[jdphi].c_str()),50,0,1.);
            h_Qsplit_to_mHard_gs_FSR[isign][jdphi] = new TH1D(Form("h_Qsplit_to_mHard_gs_FSR_sign%d%s",isign+1, dphis[jdphi].c_str()),Form("h_Qsplit_to_mHard_gs_FSR_sign%d%s",isign+1, dphis[jdphi].c_str()),50,0,1.);

            // h_ptlead_pair_pt[isign][jdphi] = new TH2D(Form("h_ptlead_pair_pt_sign%d%s",isign+1, dphis[jdphi].c_str()),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",pms.npairPT_bins,pms.pairPTBins[isign][2],pms.npt_bins,pms.pTBins);
            // h_num_hard_scatt_out[isign][jdphi] = new TH1D(Form("h_num_hard_scatt_out_sign%d%s",isign+1, dphis[jdphi].c_str()),Form("h_num_hard_scatt_out_sign%d%s",isign+1, dphis[jdphi].c_str()),3,2,5);
            // h_pt_muon_pt_closest_hadr_ratio[isign][jdphi] = new TH1D(Form("h_pt_muon_pt_closest_hadr_ratio_sign%d%s",isign+1, dphis[jdphi].c_str()),Form("h_pt_muon_pt_closest_hadr_ratio_sign%d%s",isign+1, dphis[jdphi].c_str()),40,0,1.);
            // h_pt_closest_hadr_pt_furthest_hadr_ratio[isign][jdphi] = new TH1D(Form("h_pt_closest_hadr_pt_furthest_hadr_ratio_sign%d%s",isign+1, dphis[jdphi].c_str()),Form("h_pt_closest_hadr_pt_furthest_hadr_ratio_sign%d%s",isign+1, dphis[jdphi].c_str()),40,0,1.);
            // h_pt_hadr_hq_ratio[isign][jdphi] = new TH1D(Form("h_pt_hadr_hq_ratio_sign%d%s",isign+1, dphis[jdphi].c_str()),Form("h_pt_hadr_hq_ratio_sign%d%s",isign+1, dphis[jdphi].c_str()),40,0,1.);
            // h_dphi_muon_closest_hadr[isign][jdphi] = new TH1D(Form("h_dphi_muon_closest_hadr_sign%d%s",isign+1, dphis[jdphi].c_str()),Form("h_dphi_muon_closest_hadr_sign%d%s",isign+1, dphis[jdphi].c_str()),32,-pms.PI,pms.PI);

            h_parent_groups[isign][jdphi]->Sumw2();
            h_hard_scatt_categr_vs_origin_categr[isign][jdphi]->Sumw2();
            h_both_from_b_hard_scatt_categr_vs_origin_categr[isign][jdphi]->Sumw2();
            h_both_from_c_hard_scatt_categr_vs_origin_categr[isign][jdphi]->Sumw2();
            h_one_from_b_one_from_c_hard_scatt_categr_vs_origin_categr[isign][jdphi]->Sumw2();
            h_muon_pair_origin_categr[isign][jdphi]->Sumw2();
            h_both_from_b_muon_pair_origin_categr[isign][jdphi]->Sumw2();
            h_both_from_c_muon_pair_origin_categr[isign][jdphi]->Sumw2();
            h_one_from_b_one_from_c_muon_pair_origin_categr[isign][jdphi]->Sumw2();

            h_Qsplit_gs_ISR_one_hard_scatt[isign][jdphi]->Sumw2();
            h_Qsplit_gs_FSR[isign][jdphi]->Sumw2();
            h_Qsplit_to_mHard_gs_ISR_one_hard_scatt[isign][jdphi]->Sumw2();
            h_Qsplit_to_mHard_gs_FSR[isign][jdphi]->Sumw2();

            // h_ptlead_pair_pt[isign][jdphi]->Sumw2();
            // h_num_hard_scatt_out[isign][jdphi]->Sumw2();
            // h_pt_muon_pt_closest_hadr_ratio[isign][jdphi]->Sumw2();
            // h_pt_closest_hadr_pt_furthest_hadr_ratio[isign][jdphi]->Sumw2();
            // h_pt_hadr_hq_ratio[isign][jdphi]->Sumw2();
            // h_dphi_muon_closest_hadr[isign][jdphi]->Sumw2();

            // for (int kgrp = 0; kgrp < nAncestorGroups; kgrp++){
            //     h_QQ_DR[isign][jdphi][kgrp] = new TH1D(Form("h_QQ_DR_sign%d%s%s",isign+1, dphis[jdphi].c_str(), ancestor_grps[kgrp].c_str()),";#Delta R;d#sigma/d#Delta R", 50,0,5.75);
            //     h_QQ_Dphi[isign][jdphi][kgrp] = new TH1D(Form("h_QQ_Dphi_sign%d%s%s",isign+1, dphis[jdphi].c_str(), ancestor_grps[kgrp].c_str()),";#Delta#phi;d#sigma/d#Delta#phi", 32,-pms.PI,pms.PI);
            //     h_QQ_minv[isign][jdphi][kgrp] = new TH1D(Form("h_QQ_minv_sign%d%s%s",isign+1, dphis[jdphi].c_str(), ancestor_grps[kgrp].c_str()),";m_{QQ} [GeV]; d#sigma/dm_{QQ}",pms.n_hq_minv_bins,pms.hq_minvBins);
            //     h_QQ_pair_pt_ptlead_ratio[isign][jdphi][kgrp] = new TH1D(Form("h_QQ_pair_pt_ptlead_ratio_sign%d%s%s",isign+1, dphis[jdphi].c_str(), ancestor_grps[kgrp].c_str()),";#frac{p_{T}^{pair}}{p_{T}^{lead}};d#sigma/d#frac{p_{T}^{pair}}{p_{T}^{lead}}", 25,0,2.);
            //     h_QQ_pt_avg[isign][jdphi][kgrp] = new TH1D(Form("h_QQ_pt_avg_sign%d%s%s",isign+1, dphis[jdphi].c_str(), ancestor_grps[kgrp].c_str()),";p_{T}^{avg};d#sigma/dp_{T}^{avg}", pms.n_hq_pt_bins,pms.hq_pTBins);
            //     h_QQ_asym[isign][jdphi][kgrp] = new TH1D(Form("h_QQ_asym_sign%d%s%s",isign+1, dphis[jdphi].c_str(), ancestor_grps[kgrp].c_str()),";A = (pT1 - pT2)/(pT1 + pT2);d#sigma/dA", 25,0,1);
            //     // if (kgrp != 2) // do not initialize for single gluon case
            //     h_QQ_minv_s_cm_ratio[isign][jdphi][kgrp] = new TH1D(Form("h_QQ_minv_s_cm_ratio_sign%d%s%s",isign+1, dphis[jdphi].c_str(), ancestor_grps[kgrp].c_str()),";#frac{m_{QQ}}{#hat{s}};d#sigma/d#frac{m_{QQ}}{#hat{s}}", 25,0,1);
                
            //     // h_QQ_ptlead_pair_pt[isign][jdphi] = new TH2D(Form("h_QQ_ptlead_pair_pt_sign%d%s",isign+1, dphis[jdphi].c_str()),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",pms.npairPT_bins,pms.pairPTBins[isign][2],pms.n_hq_pt_bins,pms.hq_pTBins);
            //     h_QQ_ptlead_pair_pt[isign][jdphi][kgrp] = new TH2D(Form("h_QQ_ptlead_pair_pt_sign%d%s%s",isign+1, dphis[jdphi].c_str(), ancestor_grps[kgrp].c_str()),";p_{T}^{pair} [GeV];p_{T}^{lead} [GeV]",pms.n_hq_pt_bins,pms.hq_pTBins,pms.n_hq_pt_bins,pms.hq_pTBins);
            //     h_QQ_pt1_pt2[isign][jdphi][kgrp] = new TH2D(Form("h_QQ_pt1_pt2_sign%d%s%s",isign+1, dphis[jdphi].c_str(), ancestor_grps[kgrp].c_str()),";p_{T}^{sublead} [GeV];p_{T}^{lead} [GeV]",pms.n_hq_pt_bins,pms.hq_pTBins,pms.n_hq_pt_bins,pms.hq_pTBins);
            //     h_QQ_Deta_Dphi[isign][jdphi][kgrp] = new TH2D(Form("h_QQ_Deta_Dphi_sign%d%s%s",isign+1, dphis[jdphi].c_str(), ancestor_grps[kgrp].c_str()),";#Delta#phi;#Delta#eta", 32,-pms.PI,pms.PI,40,-4.8,4.8);
            //     h_QQ_eta1_eta2[isign][jdphi][kgrp] = new TH2D(Form("h_QQ_eta1_eta2_sign%d%s%s",isign+1, dphis[jdphi].c_str(), ancestor_grps[kgrp].c_str()),";#eta_{sublead};#eta_{lead}",40,-2.4,2.4, 40,-2.4,2.4);
            //     h_QQ_minv_pair_pt[isign][jdphi][kgrp] = new TH2D(Form("h_QQ_minv_pair_pt_sign%d%s%s",isign+1, dphis[jdphi].c_str(), ancestor_grps[kgrp].c_str()),";p_{T}^{pair} [GeV];m_{QQ} [GeV]",pms.n_hq_pt_bins,pms.hq_pTBins,pms.n_hq_minv_bins,pms.hq_minvBins);
            //     h_QQ_minv_Dphi[isign][jdphi][kgrp] = new TH2D(Form("h_QQ_minv_Dphi_sign%d%s%s",isign+1, dphis[jdphi].c_str(), ancestor_grps[kgrp].c_str()),";#Delta#phi;m_{QQ} [GeV]",pms.npt_bins,pms.pTBins,pms.n_hq_minv_bins,pms.hq_minvBins);
            // }
        }
    }
}


void PythiaNTupleFirstPass::Finalize(){
    delete m1_history;
    delete m2_history;
    delete m1_history_particle;
    delete m2_history_particle;

    delete m1_first_hf_hadron_prt_pt_eta_phi_m;
    delete m2_first_hf_hadron_prt_pt_eta_phi_m;

    m_crossx_summary_file->close();
    delete m_crossx_summary_file;

    if (print_bad_warnings){
      *m_very_bad_warning_file << "Test writing into warning files [end]" << std::endl;
      m_very_bad_warning_file->close();
      delete m_very_bad_warning_file;

      m_hard_scattering_warning_file->close();
      delete m_hard_scattering_warning_file;      
    }

    if (print_low_minv_resonances){
      m_very_low_minv_resonance_file->close();
      delete m_very_low_minv_resonance_file;      
    }

    if (print_unspecified_parent && m_unspecified_parent_file){
        m_unspecified_parent_file->close();
        delete m_unspecified_parent_file;      
    }

    if (print_FE && m_FE_file){
        m_FE_file->close();
        delete m_FE_file;
    }

    if (print_other_flavor_history && m_other_flavor_category_file){
        m_other_flavor_category_file->close();
        delete m_other_flavor_category_file;
    }

    if (print_prt_history){
        for (int isign = 0; isign < ParamsSet::nSigns; isign++){
            if (m_b_parent_file[isign][0]){
                m_b_parent_file[isign][0]->close();
                delete m_b_parent_file[isign][0];
            }
            if (m_b_parent_file[isign][1]){
                m_b_parent_file[isign][1]->close();
                delete m_b_parent_file[isign][1];
            }
        }
        for (int isign = 0; isign < ParamsSet::nSigns; isign++){
            if (m_c_parent_file[isign][0]){
                m_c_parent_file[isign][0]->close();
                delete m_c_parent_file[isign][0];
            }
            if (m_c_parent_file[isign][1]){
                m_c_parent_file[isign][1]->close();
                delete m_c_parent_file[isign][1];
            }
        }
    }


    for (int jdphi = 0; jdphi < 2; jdphi++){
        for (int isign = 0; isign < ParamsSet::nSigns; isign++){
            // delete h_ptlead_pair_pt[isign][jdphi];
            delete h_parent_groups[isign][jdphi];
            delete h_hard_scatt_categr_vs_origin_categr[isign][jdphi];
            delete h_both_from_b_hard_scatt_categr_vs_origin_categr[isign][jdphi];
            delete h_both_from_c_hard_scatt_categr_vs_origin_categr[isign][jdphi];
            delete h_one_from_b_one_from_c_hard_scatt_categr_vs_origin_categr[isign][jdphi];
            delete h_muon_pair_origin_categr[isign][jdphi];
            delete h_both_from_b_muon_pair_origin_categr[isign][jdphi];
            delete h_both_from_c_muon_pair_origin_categr[isign][jdphi];
            delete h_one_from_b_one_from_c_muon_pair_origin_categr[isign][jdphi];

            delete h_Qsplit_gs_ISR_one_hard_scatt[isign][jdphi];
            delete h_Qsplit_gs_FSR[isign][jdphi];
            delete h_Qsplit_to_mHard_gs_ISR_one_hard_scatt[isign][jdphi];
            delete h_Qsplit_to_mHard_gs_FSR[isign][jdphi];

            // for (int kgrp = 0; kgrp < nAncestorGroups; kgrp++){
            //     delete h_QQ_DR[isign][jdphi][kgrp];
            //     delete h_QQ_Dphi[isign][jdphi][kgrp];
            //     delete h_QQ_minv[isign][jdphi][kgrp];

            //     delete h_QQ_ptlead_pair_pt[isign][jdphi][kgrp];
            //     delete h_QQ_pt1_pt2[isign][jdphi][kgrp];
            //     delete h_QQ_Deta_Dphi[isign][jdphi][kgrp];
            //     delete h_QQ_eta1_eta2[isign][jdphi][kgrp];
            //     delete h_QQ_minv_pair_pt[isign][jdphi][kgrp];
            //     delete h_QQ_minv_Dphi[isign][jdphi][kgrp];
            // }
        }
    }
}

void PythiaNTupleFirstPass::FillMuonPair(int pair_ind, std::shared_ptr<MuonPairPythia> const& mpair){
  // // use ind to record barcode instead
  mpair->m1.ind     = muon_pair_muon1_bar->at(pair_ind);
  mpair->m2.ind     = muon_pair_muon2_bar->at(pair_ind);
  mpair->m1.charge  = muon_pair_muon1_ch->at(pair_ind);//sign of pt stores charge
  mpair->m2.charge  = muon_pair_muon2_ch->at(pair_ind);//sign of pt stores charge

  mpair->m1.pt    = static_cast<float>(muon_pair_muon1_pt->at(pair_ind));//pt of the first muon in the pair
  mpair->m2.pt    = static_cast<float>(muon_pair_muon2_pt->at(pair_ind));//pt of the second muon in the pair
  mpair->m1.eta   = static_cast<float>(muon_pair_muon1_eta->at(pair_ind));
  mpair->m2.eta   = static_cast<float>(muon_pair_muon2_eta->at(pair_ind));
  mpair->m1.phi   = static_cast<float>(muon_pair_muon1_phi->at(pair_ind));
  mpair->m2.phi   = static_cast<float>(muon_pair_muon2_phi->at(pair_ind));

  mpair->effcy      = efficiency;
  mpair->QHard      = QHard;
  mpair->pTHat      = pTHat;
  mpair->mHat       = mHat;
}

bool PythiaNTupleFirstPass::PassCuts(const std::shared_ptr<MuonPair>& mpair){
  //Apply ALL CUTS but for resonances
  // NO quality cuts or dP over P cut (s + light hadron decay turned OFF) on the MC truth muons

  // passed muon eta cut
  if (fabs(mpair->m1.eta) > 2.4 || fabs(mpair->m2.eta) > 2.4) return false;
  h_cutAcceptance[mpair->m1.charge != mpair->m2.charge]->Fill(pass_muon_eta + 0.5, mpair->weight); // if same sign: fill the h_cutAcceptance[0] histogram; if opposite sign, fill the h_cutAcceptance[1] histogram

  // passed muon pt cut
  if (mpair->m1.pt < 4 || mpair->m2.pt < 4) return false;
  h_cutAcceptance[mpair->m1.charge != mpair->m2.charge]->Fill(pass_muon_pt + 0.5, mpair->weight);
   
  return true;
}


int PythiaNTupleFirstPass::ParentGrouping(std::vector<int>& parent_ids, bool c_tag, bool prev_is_lepton){

  if (parent_ids.size() == 2 && abs(parent_ids[0]) <= 5 && parent_ids[1] == (-1) * parent_ids[0] && prev_is_lepton){
    return prt_drell_yan;
  } 

  int parent_id = abs(parent_ids[0]) % 10000;

  auto it_truth_resonance = std::find(resonance_ids.begin(), resonance_ids.end(), parent_id);

  if (it_truth_resonance != resonance_ids.end()) { // current muon comes from resonance
    return resonance_decay;
  }

  if ((parent_id >= 500 && parent_id < 600) || (parent_id >= 5000 && parent_id < 6000)){
    if (! c_tag) return single_muon_parent_group::direct_b; // direct b
    else         return single_muon_parent_group::b_to_c; // b -> c
  }

  if (c_tag) return single_muon_parent_group::direct_c; // c not from b
  
  if ((parent_id >= 100 && parent_id < 400) || (parent_id >= 1000 && parent_id < 4000)) // light and s-hadrons
    return s_light; // light and s-flavored hadrons

  if (parent_id == 22) // photons
    return single_photon;

  std::cout << "Not in the given set of parent groups. Parent ids: " << std::endl;
  for (auto id : parent_ids) std::cout << id << " ";
  std::cout << std::endl;
  return -1; // others
}

void PythiaNTupleFirstPass::GetPtEtaPhiMFromBarcode(int barcode, std::vector<float>* pt_eta_phi_m){
  pt_eta_phi_m->clear();
  pt_eta_phi_m->push_back(truth_pt->at(barcode));
  pt_eta_phi_m->push_back(truth_eta->at(barcode));
  pt_eta_phi_m->push_back(truth_phi->at(barcode));
}

std::pair<int,int> PythiaNTupleFirstPass::UpdateCurParents(bool isMuon1, std::vector<int>& cur_prt_bars, std::vector<int>& cur_prt_ids, bool before_gs, bool& prev_out_hard_scatt, int hf_quark_index = -1000000){
  
  // function that traces one step back and updates the current parent barcodes + ids
  // it returns a pair of integer that is {previous first parent barcode, previous first parent id} if the previous first parent was a HF hadron
  // and returns {0,0} otherwise
  // besides tracing, it also handles many specific cases, such as resonance tagging and c tagging, as well as hadron to quark transition

  bool cur_is_quark_gluon = (abs(cur_prt_ids[0]) <= 5 || cur_prt_ids[0] == 21);
  if (cur_prt_bars.size() > 1 && !(cur_is_quark_gluon)){
  // if (cur_prt_bars.size() > 1 && hf_quark_index < 0){
    // LATER REMEMBER TO ADD [EXCLUDE 21 21 / 21 2]
    if (print_unspecified_parent && m_unspecified_parent_file){
      *m_unspecified_parent_file << "Event#: " << mpair->m1.ev_num << std::endl;
      *m_unspecified_parent_file << "More than one parent at hadronic level:" << std::endl;
      for (int parent_id : cur_prt_ids) *m_unspecified_parent_file << parent_id << " ";
      *m_unspecified_parent_file << std::endl << std::endl;      
    }
  }

  int prev_first_prt_id = cur_prt_ids[0];
  int prev_first_prt_bar = cur_prt_bars[0];
  cur_prt_bars.clear();
  cur_prt_ids.clear();

  int index = (hf_quark_index >= 0)? hf_quark_index : 0;
  int mother1 = truth_mother1->at(cur_prt_bars[index]);
  int mother2 = truth_mother2->at(cur_prt_bars[index]);
  int status = abs(truth_status->at(cur_prt_bars[index]));

  // std::cout << "Event: " << mpair->m1.ev_num << ", barcode " << cur_prt_bars[index] << ", id " << truth_id->at(cur_prt_bars[index]) << ", status " << status << ", mother1 " << mother1 << ", mother2 " << mother2 << std::endl;

  if (mother1 == 0){
    std::cout << "Error:: Mother1 cannot be 0, quitting" << std::endl;
    PrintHistory(&std::cout, true, isMuon1);
    cout << "bar: " << prev_first_prt_bar << ", mother1 " << mother1 << ", mother2 " << mother2 << endl;
  }

  if (mother2 == mother1 || mother2 == 0){ // only one mother
    cur_prt_bars.push_back(mother1);
    cur_prt_ids.push_back(truth_id->at(mother1) % 10000);
    if ((mother1 == 1 || mother1 == 2) && before_gs){
      if (isMuon1) m1_ancestor_is_incoming = true;
      else         m2_ancestor_is_incoming = true;
    }
  }else if(((status >= 81 && status <= 86) || (status >= 101 && status <= 106)) && mother1 < mother2){ // all particles between mother1 and mother2 are parents
  // }else if(((status >= 81 && status <= 86) || (status >= 101 && status <= 106))){ // all particles between mother1 and mother2 are parents
    // std::cout << "Event: " << mpair->m1.ev_num << ", status " << status << ", mother1 " << mother1 << ", mother2 " << mother2 << std::endl;
    for (int cur_mother = mother1; cur_mother <= mother2; cur_mother++){
      cur_prt_bars.push_back(cur_mother);
      cur_prt_ids.push_back(truth_id->at(cur_mother) % 10000);
    }
  }else{ // two truly different parents 
    cur_prt_bars.push_back(mother1);
    cur_prt_ids.push_back(truth_id->at(mother1) % 10000);

    cur_prt_bars.push_back(mother2);
    cur_prt_ids.push_back(truth_id->at(mother2) % 10000);
  }

  if (cur_prt_bars.size() == 0 || cur_prt_ids.size() == 0 && (print_bad_warnings && m_very_bad_warning_file && m_hard_scattering_warning_file)){
    *m_very_bad_warning_file << "Error:: current parent barcodes & ids have zero length, quitting" << std::endl;
    *m_very_bad_warning_file << "Event: " << mpair->m1.ev_num << ", previous first parent barcode: " << prev_first_prt_bar << std::endl;
    *m_very_bad_warning_file << "Mother1: " << mother1 << ", Mother2: " << mother2 << std::endl;
  }

  // fill in vector of the profiles of the current parents
  std::vector<Particle> cur_prt_profiles {};
  for (int bar : cur_prt_bars){
    float pt = static_cast<float>(truth_pt->at(bar));
    float eta = static_cast<float>(truth_eta->at(bar));
    float phi = static_cast<float>(truth_phi->at(bar));
    int id = truth_id->at(bar);
    int status = truth_status->at(bar);
    Particle p {pt, eta, phi, bar, id, status};
    cur_prt_profiles.push_back(p);
  }

  // save the vector of the updated parent ids into the muon history chain
  if (isMuon1){
    m1_history->push_back(cur_prt_ids);
    m1_history_particle->push_back(cur_prt_profiles);
  }
  else{
    m2_history->push_back(cur_prt_ids);
    m2_history_particle->push_back(cur_prt_profiles);
  }

  if (prev_out_hard_scatt){
    if (cur_prt_bars.size() != 2 || (abs(truth_status->at(cur_prt_bars[0])) != 21 && abs(truth_status->at(cur_prt_bars[0])) != 31) && (print_bad_warnings && m_very_bad_warning_file && m_hard_scattering_warning_file)){
      *m_very_bad_warning_file << "The current parents should be the incoming particles of a hard scattering." << std::endl;
      PrintHistory(m_very_bad_warning_file, true, isMuon1);
    }

    int catgr = HardAnalysisCategr(cur_prt_bars[0], cur_prt_bars[1]);

    if ((abs(truth_status->at(cur_prt_bars[0]+2)) != 23 && abs(truth_status->at(cur_prt_bars[0]+2)) != 33) || (abs(truth_status->at(cur_prt_bars[1]+2)) != 23 && abs(truth_status->at(cur_prt_bars[1]+2)) != 33) && (print_bad_warnings && m_very_bad_warning_file && m_hard_scattering_warning_file)){
      *m_very_bad_warning_file << "Outgoing particles of a HS don't have barcode as 2 + incoming particles' barcodes" << std::endl;
      PrintHistory(m_very_bad_warning_file, true, isMuon1);
    }

    // std::cout << "Event# " << mpair->m1.ev_num << ", current barcode: " << cur_prt_bars[0] << std::endl;

    if (isMuon1){
      vm1out1.SetPtEtaPhiM(truth_pt->at(cur_prt_bars[0]+2), truth_eta->at(cur_prt_bars[0]+2), truth_phi->at(cur_prt_bars[0]+2), truth_m->at(cur_prt_bars[0]+2));
      vm1out2.SetPtEtaPhiM(truth_pt->at(cur_prt_bars[1]+2), truth_eta->at(cur_prt_bars[1]+2), truth_phi->at(cur_prt_bars[1]+2), truth_m->at(cur_prt_bars[1]+2));
    }else{
      vm2out1.SetPtEtaPhiM(truth_pt->at(cur_prt_bars[0]+2), truth_eta->at(cur_prt_bars[0]+2), truth_phi->at(cur_prt_bars[0]+2), truth_m->at(cur_prt_bars[0]+2));
      vm2out2.SetPtEtaPhiM(truth_pt->at(cur_prt_bars[1]+2), truth_eta->at(cur_prt_bars[1]+2), truth_phi->at(cur_prt_bars[1]+2), truth_m->at(cur_prt_bars[1]+2));
    }

    if (isMuon1){
      if (before_gs)  m1_from_hard_scatt_before_gs = true;
      else            m1_from_hard_scatt_after_gs = true;
      m1_hard_scatt_in_bar1 = cur_prt_bars[0];
      mpair->m1_hard_scatt_category = catgr;
      TLorentzVector vhard = vm1out1 + vm1out2;

      m1_hard_scatt_Q = vhard.M();
    }  
    else{
      if (before_gs)  m2_from_hard_scatt_before_gs = true;
      else            m2_from_hard_scatt_after_gs = true;
      m2_hard_scatt_in_bar1 = cur_prt_bars[0];
      mpair->m2_hard_scatt_category = catgr;
      TLorentzVector vhard = vm2out1 + vm2out2;
      m2_hard_scatt_Q = vhard.M();
    }

    prev_out_hard_scatt = false;
  }

  if (abs(truth_status->at(cur_prt_bars[0])) == 23 || abs(truth_status->at(cur_prt_bars[0])) == 33){
    if ((isMuon1 && m1_hard_scatt_in_bar1 > 0) || (!isMuon1 && m2_hard_scatt_in_bar1 > 0) && (print_bad_warnings && m_very_bad_warning_file && m_hard_scattering_warning_file)){
      *m_very_bad_warning_file << "There should be at most one hard scattering in a muon's history chain." << std::endl;
      PrintHistory(m_very_bad_warning_file, true, isMuon1);
    }

    if (cur_prt_bars.size() != 1 && (print_bad_warnings && m_very_bad_warning_file && m_hard_scattering_warning_file)){
      *m_very_bad_warning_file << "There should only be one [outgoing particle from a hard scattering] as part of a parent history chain." << std::endl;
      PrintHistory(m_very_bad_warning_file, true, isMuon1);
    }
    prev_out_hard_scatt = true;
  }

  // record if there is oscillation
  if (cur_prt_ids[0] == (-1) * prev_first_prt_id && abs(prev_first_prt_id) != 4 && abs(prev_first_prt_id) != 5){
    if (isMuon1) m1_osc = true;
    else         m2_osc = true;
  }

  auto it_truth_resonance = std::find(resonance_ids.begin(), resonance_ids.end(), abs(cur_prt_ids[0]) % 10000);

  if (it_truth_resonance != resonance_ids.end()) { // current muon comes from resonance
    if (isMuon1){
      m1_resonance_barcode = cur_prt_bars[0];
      mpair->m1_parent_group = resonance_decay; 
    }
    else{
      m2_resonance_barcode = cur_prt_bars[0];
      mpair->m2_parent_group = resonance_decay; 
    }
    return std::pair<int,int>(0,0);
  }

  //c-tag
  if ((abs(cur_prt_ids[0]) >= 400 && abs(cur_prt_ids[0]) < 500) || (abs(cur_prt_ids[0]) >= 4000 && abs(cur_prt_ids[0]) < 5000)){ // c-hadrons
    if (isMuon1) m1_c_tag = true;
    else         m2_c_tag = true;
  }


  // from muon to hadron stage
  if (abs(prev_first_prt_id) == 13 || abs(prev_first_prt_id) == 15){
    if ((abs(cur_prt_ids[0]) >= 400 && abs(cur_prt_ids[0]) < 600) || (abs(cur_prt_ids[0]) >= 4000 && abs(cur_prt_ids[0]) < 6000)){ // c-hadrons
      if (isMuon1)  GetPtEtaPhiMFromBarcode(cur_prt_bars[0], &mpair->m1_last_hf_hadron_prt_pt_eta_phi_m);
      else          GetPtEtaPhiMFromBarcode(cur_prt_bars[0], &mpair->m2_last_hf_hadron_prt_pt_eta_phi_m);
    }
  }

  // last b-flavored hadron
  bool prev_is_b_hadron = ((abs(prev_first_prt_id) >= 500 && abs(prev_first_prt_id) < 600) || (abs(prev_first_prt_id) >= 5000 && abs(prev_first_prt_id) < 6000));
  bool current_is_b_hadron = ((abs(cur_prt_ids[0]) >= 500 && abs(cur_prt_ids[0]) < 600) || (abs(cur_prt_ids[0]) >= 5000 && abs(cur_prt_ids[0]) < 6000));
  if ((!prev_is_b_hadron) && current_is_b_hadron){ // the last b hadron
    if (isMuon1)  GetPtEtaPhiMFromBarcode(cur_prt_bars[0], &mpair->m1_last_b_hadron_prt_pt_eta_phi_m);
    else          GetPtEtaPhiMFromBarcode(cur_prt_bars[0], &mpair->m2_last_b_hadron_prt_pt_eta_phi_m);
  }

  // from hadron to quark/gluon stage
  bool prev_is_hf_hadron = ((abs(prev_first_prt_id) >= 400 && abs(prev_first_prt_id) < 600) || (abs(prev_first_prt_id) >= 4000 && abs(prev_first_prt_id) < 6000));
  if (prev_is_hf_hadron && ((abs(cur_prt_ids[0]) < 4000 && abs(cur_prt_ids[0]) % 100 <= 5) || cur_prt_ids[0] == 21)){ // including light diquarks
    if (isMuon1)  GetPtEtaPhiMFromBarcode(prev_first_prt_bar, m1_first_hf_hadron_prt_pt_eta_phi_m);
    else          GetPtEtaPhiMFromBarcode(prev_first_prt_bar, m2_first_hf_hadron_prt_pt_eta_phi_m);
    return std::pair<int,int>(prev_first_prt_bar, prev_first_prt_id);
  }
  return std::pair<int,int>(0,0);
}


int PythiaNTupleFirstPass::FindHeavyQuarks(std::vector<int>& cur_prt_ids, std::vector<int>& cur_prt_bars, int quark_type, bool isMuon1, int prev_hq_bar, int hadron_child_id = 0){
  if ((quark_type != 4 && quark_type != 5) && (print_bad_warnings && m_very_bad_warning_file && m_hard_scattering_warning_file)){
    *m_very_bad_warning_file << "Error:: the parameter quark_type must take value of 4 (c) or 5 (b), quitting" << std::endl;
    *m_very_bad_warning_file << "Event #: " << mpair->m1.ev_num << std::endl;
  }
  int status = truth_status->at(cur_prt_bars[0]);

  if (cur_prt_ids.size() == 0 && (print_bad_warnings && m_very_bad_warning_file && m_hard_scattering_warning_file)){
    *m_very_bad_warning_file << "Error:: parent list is empty, quitting" << std::endl;
    *m_very_bad_warning_file << "Event #: " << mpair->m1.ev_num << std::endl;
  }

  if (cur_prt_ids.size() == 1 && (cur_prt_ids[0] == 2212 || cur_prt_ids[0] == 2112) ){
    return -2; // case II
  }

  int index_candidate = -1; // default value if -1 (only to be reassigned if found either Q or Qbar)

  std::vector<int>::iterator it_q = std::find(cur_prt_ids.begin(),cur_prt_ids.end(),quark_type);
  std::vector<int>::iterator it_qbar = std::find(cur_prt_ids.begin(),cur_prt_ids.end(), (-1) * quark_type);
  if (it_q != cur_prt_ids.end()){ // found q
    index_candidate = it_q - cur_prt_ids.begin();
    if (it_qbar != cur_prt_ids.end()){ // also found qbar
      // int sign_correction = (quark_type == 4)? +1 : -1; // b quark to hadron changes sign
      int sign_correction = (abs(hadron_child_id) >= 500 && abs(hadron_child_id) < 600)? -1 : +1; // b quark to B meson changes sign; c to D meson or b/c(bar) to hadron does not change sign
      if (hadron_child_id * sign_correction > 0)
        index_candidate = it_q - cur_prt_ids.begin();
      else if (hadron_child_id * sign_correction < 0)
        index_candidate = it_qbar - cur_prt_ids.begin();
      else if (abs(status) == 21 || abs(status) == 31){ // hard scattering is Q Qbar (or Qbar Q)
        if (prev_hq_bar - 2 == cur_prt_bars[0]){
          index_candidate = 0;
        }else{
          index_candidate = 1;
        }
        // sanity check: require the quark at prev_hq_bar to have the same id as the corresponding quark in cur_prt_ids
        if (truth_id->at(prev_hq_bar) != cur_prt_ids[index_candidate]){
          std::cout << "hard scattering being Q Qbar or Qbar Q" << std::endl;
          std::cout << "the quark at prev_hq_bar must have the same id as the corresponding quark in cur_prt_ids" << std::endl;
          std::cout << "the previous HQ barcode is " << prev_hq_bar << ", the index candidate is " << index_candidate << std::endl;
          PrintHistory(&std::cout, true, isMuon1);
        }
      }else{ // not the immediate parent of a hadron (hadron_child_id == 0)
        if (print_bad_warnings && m_very_bad_warning_file){
          *m_very_bad_warning_file << "The current parent vector has both Q and Qbar, and is not input to hard scattering, yet is not the immediate parent of a hadron." << std::endl;
          PrintHistory(m_very_bad_warning_file, true, isMuon1);          
        }
      }
      // std::cout << "To explicitly check for the case of q and qbar in the same parent vector: " << std::endl;
      // PrintHistory(&std::cout, true, isMuon1);
      // std::cout << "The determined quark index is " << index_candidate << std::endl << std::endl;
    }

    if (std::find(it_q + 1,cur_prt_ids.end(),quark_type) != cur_prt_ids.end()){ // find two q's 
      // require it's from a hard scattering of c c --> c c
      if (cur_prt_ids.size() != 2 || (abs(status) != 21 && abs(status) != 31) || it_qbar != cur_prt_ids.end()){
        // std::cout << "In the current parents there are more than one same-sign, same-flavor HQ's." << std::endl;
        // std::cout << "Yet the parents are not incoming particles of a hard scattering." << std::endl;
        // PrintHistory(&std::cout, true, isMuon1);
        skip_event_origin_analysis = true;
        return index_candidate;
      }
      // e.g, current may have barcode {3 4}, and previous can either be 5 or 6
      // if (prev_hq_bar - 2 != cur_prt_bars[0] && prev_hq_bar - 2 != cur_prt_bars[1]){
        // std::cout << "Previous barcode (outgoing from a hard scattering minus two) is not an element in the current parents' barcode link." << std::endl;
        // std::cout << prev_hq_bar << " " << it_qbar_second - cur_prt_ids.begin() << " " << cur_prt_ids.end() - cur_prt_ids.begin() << std::endl;
        // PrintHistory(&std::cout, true, isMuon1);
      // }
      if (prev_hq_bar - 2 == cur_prt_bars[0]){
        index_candidate = 0;
      }else{
        index_candidate = 1;
      }
    }
  }
  else if (it_qbar != cur_prt_ids.end()){
    index_candidate = it_qbar - cur_prt_ids.begin();
    std::vector<int>::iterator it_qbar_second = std::find(it_qbar + 1,cur_prt_ids.end(), (-1) * quark_type);
    if (it_qbar_second != cur_prt_ids.end()){ // find two q's 
      // cout << prev_hq_bar << " " << it_qbar_second - cur_prt_ids.begin() << " " << cur_prt_ids.end() - cur_prt_ids.begin() << endl;
      if (cur_prt_ids.size() != 2 || (abs(status) != 21 && abs(status) != 31)){
        std::cout << "In the current parents there are more than one same-sign, same-flavor HQ's." << std::endl;
        std::cout << "Yet the parents are not incoming particles of a hard scattering." << std::endl;
        PrintHistory(&std::cout, true, isMuon1);
        return index_candidate;
      }
      // if (prev_hq_bar - 2 != cur_prt_bars[0] && prev_hq_bar - 2 != cur_prt_bars[1]){
        // std::cout << prev_hq_bar << " " << it_qbar_second - cur_prt_ids.begin() << " " << cur_prt_ids.end() - cur_prt_ids.begin() << std::endl;
        // std::cout << "Previous barcode (outgoing from a hard scattering minus two is not an elment in the current parents' barcode link." << std::endl;
        // PrintHistory(&std::cout, true, isMuon1);
      // }
      if (prev_hq_bar - 2 == cur_prt_bars[0]){
        index_candidate = 0;
      }else{
        index_candidate = 1;
      }
    }
  }

  // if (mpair->m1.ev_num ==58 || mpair->m1.ev_num == 117 || mpair->m1.ev_num == 203 || mpair->m1.ev_num == 293 || mpair->m1.ev_num == 478 || mpair->m1.ev_num == 484 ||  mpair->m1.ev_num ==601){
  //   std::cout << "FindHeavyQuarks - current step parent ids: " << std::endl;
  //   for (int i = 0; i < cur_prt_ids.size(); i++)
  //     std::cout << cur_prt_ids[i] << " ";
  //   std::cout << std::endl << "The heavy quark index is: " << index_candidate << std::endl;
  // }
  return index_candidate;
}

void PythiaNTupleFirstPass::SingleMuonAncestorTracing(bool isMuon1){

  // ----------------------------------------------------------------
  // Single-muon ancestor tracing through leptonic, hadronic and HQ level
  // For a muon from open HQ, traces all the way back to partonic parents of the earliest HQ ancestor
  // ----------------------------------------------------------------
  std::vector<int> parent_bars;
  std::vector<int> parent_ids;
  int first_hadron_id = 0;
  int first_hadron_barcode = 0;
  int prev_first_prt_id = -1;
  std::pair<int, int> first_hadron_bar_id;

  bool prev_step_status_m23_m33 = false; // for sanity check

  // ---------------- add muon itself into history ----------------
  float pt = (isMuon1)? static_cast<float>(mpair->m1.pt) : static_cast<float>(mpair->m2.pt);
  float eta = (isMuon1)? static_cast<float>(mpair->m1.eta) : static_cast<float>(mpair->m2.eta);
  float phi = (isMuon1)? static_cast<float>(mpair->m1.phi) : static_cast<float>(mpair->m2.phi);
  int barcode = (isMuon1)? mpair->m1.ind : mpair->m2.ind;
  int charge = (isMuon1)? mpair->m1.charge : mpair->m2.charge;
  int id = -13 * charge;
  int status = truth_status->at(barcode);
  Particle p {pt, eta, phi, barcode, id, status};

  parent_bars.push_back(barcode);
  parent_ids.push_back(id);
  std::vector<Particle> cur_muon_profile {};
  cur_muon_profile.push_back(p);

  if (isMuon1){
    m1_history->push_back(parent_ids);
    m1_history_particle->push_back(cur_muon_profile);
  }
  else{
    m2_history->push_back(parent_ids);
    m2_history_particle->push_back(cur_muon_profile);
  }

  // ---------------- trace back to first non-lepton, non-resonance & non-c-hadron parent ----------------

  // leptons, resonances & c-flavored hadrons
  while(abs(parent_ids[0]) == 13 || abs(parent_ids[0]) == 15 || (abs(parent_ids[0]) >= 400 && abs(parent_ids[0]) < 500) || (abs(parent_ids[0]) >= 4000 && abs(parent_ids[0]) < 5000)){
    if (abs(parent_ids[0]) == 15){
      if (isMuon1)  m1_from_tau = true;
      else          m2_from_tau = true;
    }
    prev_first_prt_id = parent_ids[0];
    first_hadron_bar_id = UpdateCurParents(isMuon1, parent_bars, parent_ids, true, prev_step_status_m23_m33);

    // if muon comes from resonance (J/Psi in this case): return --> avoid J/Psi from B causing muon to be tagged as b-to-c
    if ((isMuon1 && mpair->m1_parent_group == resonance_decay) || (!isMuon1 && mpair->m2_parent_group == resonance_decay)) return;

    first_hadron_barcode = first_hadron_bar_id.first;
    first_hadron_id = first_hadron_bar_id.second;
  }

  // we have encounted the first non-c-flavored hadronic parent here: if the muon is from a resonance decay, the resonance tag has been turned on, and we return without performing the rest
  if ((isMuon1 && mpair->m1_parent_group == resonance_decay) || (!isMuon1 && mpair->m2_parent_group == resonance_decay)) return;

  // ---------------- determine single-muon parent group for non-resonance case ----------------

  bool prev_is_lepton = (abs(prev_first_prt_id) == 13 || abs(prev_first_prt_id) == 15);

  int cur_parent_group; // necessary to determine at this stage for the heavy-quark finding step
  if (isMuon1){
    m1_earliest_parent_id = parent_ids[0];
    m1_youngest_non_chadron_parent_barcode = parent_bars[0];
    cur_parent_group = ParentGrouping(parent_ids, m1_c_tag, prev_is_lepton);
    mpair->m1_parent_group = cur_parent_group;
  } 
  else{
    m2_earliest_parent_id = parent_ids[0];
    m2_youngest_non_chadron_parent_barcode = parent_bars[0];
    cur_parent_group = ParentGrouping(parent_ids, m2_c_tag, prev_is_lepton);
    mpair->m2_parent_group = cur_parent_group;
  }

  // return at this point if single muon is not from (open) HF
  if (cur_parent_group != single_muon_parent_group::direct_b && cur_parent_group != single_muon_parent_group::b_to_c && cur_parent_group != single_muon_parent_group::direct_c){ // current muon not from HF
    return;
  }

  // trace back through all b-flavored hadrons
  if((abs(parent_ids[0]) >= 500 && abs(parent_ids[0]) < 600) || (abs(parent_ids[0]) >= 5000 && abs(parent_ids[0]) < 6000)){
    while((abs(parent_ids[0]) >= 500 && abs(parent_ids[0]) < 600) || (abs(parent_ids[0]) >= 5000 && abs(parent_ids[0]) < 6000)){
      first_hadron_bar_id = UpdateCurParents(isMuon1, parent_bars, parent_ids, true, prev_step_status_m23_m33);
      first_hadron_barcode = first_hadron_bar_id.first;
      first_hadron_id = first_hadron_bar_id.second;
    }
    if (isMuon1) m1_eldest_bhadron_barcode = first_hadron_barcode; // != -10 only if the current (single) muon is from a b hadron (either direct b or b to c to muon)
    else         m2_eldest_bhadron_barcode = first_hadron_barcode; // != -10 only if the current (single) muon is from a b hadron (either direct b or b to c to muon)
  }

  // ---------------- Heavy-quark level tracing ----------------

  int quark = (cur_parent_group == single_muon_parent_group::direct_c)? 4:5;
  
  int quark_index = FindHeavyQuarks(parent_ids, parent_bars, quark, isMuon1, parent_bars[0], first_hadron_id);

  int prev_hq_bar = -1;
  while(quark_index >= 0){
    prev_hq_bar = parent_bars[quark_index];
    UpdateCurParents(isMuon1, parent_bars, parent_ids, true, prev_step_status_m23_m33, quark_index);
    quark_index = FindHeavyQuarks(parent_ids, parent_bars, quark, isMuon1, prev_hq_bar);
  }

  if (prev_hq_bar < 0){
    std::cout << "Previous HQ barcode is -1: heavy quark never found." << std::endl;
    std::cout << "In this case the ANCESTOR TRACING IS NOT TRUST-WORTHY!!" << std::endl;
    std::cout << "In fact, the parton ancestor information should be left unfilled." << std::endl;
    PrintHistory(&std::cout, true, isMuon1);
    skip_event_origin_analysis = true;
    return;
  }

  // record the kinematics (pt, eta, phi, m) of the first (earliest) HQ parent
  // needed for gluon-splitting kinematic studies
  if (isMuon1){
    m1_parton_ancestor_ids = parent_ids;
    m1_parton_ancestor_bars = parent_bars;
    GetPtEtaPhiMFromBarcode(prev_hq_bar, &mpair->m1_first_hq_ancestor_pt_eta_phi_m);
  }
  else{
    m2_parton_ancestor_ids = parent_ids;
    m2_parton_ancestor_bars = parent_bars;
    GetPtEtaPhiMFromBarcode(prev_hq_bar, &mpair->m2_first_hq_ancestor_pt_eta_phi_m);
  }
}


int PythiaNTupleFirstPass::HardAnalysisCategr(int in_bar1, int in_bar2){
  // determine the category of a given hard scattering
  // the default is the hardest scattering
  // if given input barcodes (two integers), then categorize the 
  // categories:
  // 1. gg --> gg
  // 2. gq --> gq
  // 3. gg --> QQbar, qqbar --> QQbar (FC)
  // 4. gQ --> gQ, qQ --> qQ (FE)
  // 5. QQ' --> QQ'

  int status1 = abs(truth_status->at(in_bar1));
  int status2 = abs(truth_status->at(in_bar2));
  int id1 = abs(truth_id->at(in_bar1));
  int id2 = abs(truth_id->at(in_bar2));
  if ((in_bar2 != in_bar1 + 1 || (status1 != 21 && status1 != 31) || (status2 != 21 && status2 != 31)) && (print_bad_warnings && m_very_bad_warning_file && m_hard_scattering_warning_file)){
    *m_hard_scattering_warning_file << "Input paramters must be valid barcodes of incoming particles of a hard scattering!" << std::endl;
  }

  bool in1_heavy_quark = (id1 == 4 || id1 == 5);
  bool in2_heavy_quark = (id2 == 4 || id2 == 5);
  bool in1_g_lightq = (id1 <= 3 || id1 == 21);
  bool in2_g_lightq = (id2 <= 3 || id2 == 21);
  int out_bar1 = in_bar1 + 2;
  int out_bar2 = in_bar2 + 2;
  int out_id1 = abs(truth_id->at(out_bar1));
  int out_id2 = abs(truth_id->at(out_bar2));

  bool flavor_excitation = ((in1_heavy_quark && in2_g_lightq) || (in2_heavy_quark && in1_g_lightq));
  if (flavor_excitation) return flavor_excit;

  // HQ scattering - require flavor conservation (since we can have, e.g, c cbar --> b bbar, which would be flavor creation)
  if (in1_heavy_quark && in2_heavy_quark && truth_id->at(in_bar1) == truth_id->at(out_bar1) && truth_id->at(in_bar2) == truth_id->at(out_bar2)) return hq_scatt;
  
  bool final_QQbar = ((out_id1 == 4 || out_id1 == 5) && (truth_id->at(out_bar2) == - truth_id->at(out_bar1)));
  if (final_QQbar){ // we've ruled out Q -Q --> Q -Q, so final state being -+Q +-Q indicates flavor creation
    if (!((id1 <= 5 || id1 == 21) && id2 == id1) && (print_bad_warnings && m_very_bad_warning_file && m_hard_scattering_warning_file)){
      *m_hard_scattering_warning_file << "Hard scattering has final state " << truth_id->at(out_bar1) << " " << truth_id->at(out_bar2);
      *m_hard_scattering_warning_file << ", but initial state is " << truth_id->at(in_bar1) << " " << truth_id->at(in_bar2) << std::endl;
    }
    return flavor_creat;
  }

  if (id1 == 21 && id2 == 21){ // at this stage ruled out FC: this has to be  gg --> gg (since gg --> q qbar makes no sense)
    if (truth_id->at(out_bar1) == 21) return gg_gg;
    if(out_id1 > 3 && (print_bad_warnings && m_very_bad_warning_file && m_hard_scattering_warning_file)){
      *m_hard_scattering_warning_file << "Hard scattering should be gg --> qqbar, but it is";
      *m_hard_scattering_warning_file << truth_id->at(in_bar1) << " " << truth_id->at(in_bar2) << " " << truth_id->at(out_bar1) << " " << truth_id->at(out_bar2) << std::endl;
      *m_hard_scattering_warning_file << "Muon pair minv: " << mpair->minv << std::endl;
      *m_hard_scattering_warning_file << "History: " << std::endl;
      PrintHistory(m_hard_scattering_warning_file, true, mpair->from_same_ancestors);
    }
    return gg_qqbar;
  }
  if (out_id1 == 21 && out_id2 == 21){
    if (id1 <= 3) return qqbar_gg;
    if (id1 != 4 && id2 != 5 && (print_bad_warnings && m_very_bad_warning_file && m_hard_scattering_warning_file)){
      *m_hard_scattering_warning_file << "Hard scattering should be QQbar --> gg, but it is";
      *m_hard_scattering_warning_file << truth_id->at(in_bar1) << " " << truth_id->at(in_bar2) << " " << truth_id->at(out_bar1) << " " << truth_id->at(out_bar2) << std::endl;
    }
    // return QQbar_gg;
  }

  if ((id1 == 21 && id2 <=3) || (id2 == 21 && id1 <=3)) return gq_gq;

  if ((out_id1 == 13 || out_id1 == 15) && (out_id2 == 13 || out_id2 == 15)){
    if (id1 != id2 && (print_bad_warnings && m_very_bad_warning_file && m_hard_scattering_warning_file)){
      *m_hard_scattering_warning_file << "Hard scattering should be Drell-Yan, but the incoming particles are unexpected." << std::endl;
    }
    return hard_scatt_drell_yan;
  }
  // at this stage, only light-quark scattering is left
  // qq' --> qq' should NOT appear in the parental chain of a single muon
  // but it's worth doing a double/sanity check
  if ((id1 > 3 || id2 > 3 || out_id1 > 3 || out_id2 > 3) && (print_bad_warnings && m_very_bad_warning_file && m_hard_scattering_warning_file)){
    *m_hard_scattering_warning_file << "Hard scattering should be qq' --> qq', but it is";
    *m_hard_scattering_warning_file << truth_id->at(in_bar1) << " " << truth_id->at(in_bar2) << " " << truth_id->at(out_bar1) << " " << truth_id->at(out_bar2) << std::endl;
  }
  return qqprime_qqprime;

}


void PythiaNTupleFirstPass::PrintHistory(std::ostream* f, bool print_single, bool muon1_sameancestor, bool print_category = false){
  // print_single: if True then print single muon history; else then print both muons' history
  // muon1_sameancestor: meaning depends on print_single
  // if (print_single): True = muon1, False = muon2
  // if (!print_single): True = same partonic ancestors, False = different partonic ancestors
  // if True then print single muon history; else then print both muons' history
  try{
    if (!f) throw std::runtime_error("PrintHistory: the std::ostream pointer is nullptr!");
  }catch(const std::runtime_error& e){
    std::cout << "Runtime error caught in function PrintHistory: " << e.what() << std::endl;
    std::cout << "Return without printing history!" << std::endl;
    return; // return so that code doesn't crash
  }

  *f << "Event #: " << mpair->m1.ev_num << std::endl;

  if (print_single){
    if (print_category){
      int hard_scatt_category = (muon1_sameancestor)? mpair->m1_hard_scatt_category : mpair->m2_hard_scatt_category;
      std::string hard_category_name = (hard_scatt_category > 0)? HardScattCategoryLabels[hard_scatt_category] : "no hard scattering (yet)";
      std::string pair_origin_name;
      if (mpair->muon_pair_origin_category > 0) pair_origin_name = HFMuonPairOriginCategoryLabels[mpair->muon_pair_origin_category];
      else pair_origin_name = "origin not determined (yet)";
      *f << "hard scattering category: " << hard_category_name << std::endl;
      *f << "muon pair origin category: " << pair_origin_name << std::endl;
    }

    std::vector<std::vector<Particle>>* m_history_particle = (muon1_sameancestor)? m1_history_particle : m2_history_particle;
    // *f << "id (barcode, status, pt, eta, phi)" << std::endl;
    *f << "id (barcode, status)" << std::endl;
    for (size_t i = 0; i < (*m_history_particle).size(); i++){
      for (size_t j = 0; j < (*m_history_particle)[i].size(); j++) {
        *f << (*m_history_particle)[i][j].id << " (" 
        << (*m_history_particle)[i][j].barcode  << ", "
        << (*m_history_particle)[i][j].status   << ") ";
        // << (*m_history_particle)[i][j].pt       << ", "
        // << (*m_history_particle)[i][j].eta      << ", "
        // << (*m_history_particle)[i][j].phi      << ") ";
      }
      *f << "<--- ";
    }
    *f << std::endl << std::endl;
  }
  else{ // print out both muons' history
    if (print_category){
      // int hard_scatt_category = (muon1_sameancestor)? mpair->m1_hard_scatt_category : mpair->m2_hard_scatt_category;
      // std::string hard_category_name = (hard_scatt_category > 0)? HardScattCategoryLabels[hard_scatt_category] : "no hard scattering (yet)";
      std::string same_ancestors = (mpair->from_same_ancestors > 0)? "same ancestors" : "either different ancestors or undetermined";
      
      std::string pair_origin_name;
      if (mpair->muon_pair_origin_category > 0) pair_origin_name = HFMuonPairOriginCategoryLabels[mpair->muon_pair_origin_category];
      else pair_origin_name = "origin not determined (yet)";
      *f << same_ancestors << std::endl;
      *f << "muon pair origin category: " << pair_origin_name << std::endl;
    }

    // if (print_if_same_ancestor){
    //   if (muon1_sameancestor) // same ancestor
    //     *f << "Same ancestors." << std::endl;
    //   else
    //     *f << "Different ancestors." << std::endl;
    // }

    // *f << "id (barcode, status, pt, eta, phi)" << std::endl;
    *f << "id (barcode, status)" << std::endl;
    for (size_t i = 0; i < (*m1_history_particle).size(); i++){
      for (size_t j = 0; j < (*m1_history_particle)[i].size(); j++) {
        *f << (*m1_history_particle)[i][j].id << " (" 
        << (*m1_history_particle)[i][j].barcode  << ", "
        << (*m1_history_particle)[i][j].status   << ") ";
        // << (*m1_history_particle)[i][j].pt       << ", "
        // << (*m1_history_particle)[i][j].eta      << ", "
        // << (*m1_history_particle)[i][j].phi      << ") ";
      }
      *f << "<--- ";
    }
    *f << std::endl << std::endl;

    // *f << "id (barcode, status, pt, eta, phi)" << std::endl;
    // *f << "id (barcode, status)" << std::endl;
    for (size_t i = 0; i < (*m2_history_particle).size(); i++){
      for (size_t j = 0; j < (*m2_history_particle)[i].size(); j++) {
        *f << (*m2_history_particle)[i][j].id << " (" 
        << (*m2_history_particle)[i][j].barcode  << ", "
        << (*m2_history_particle)[i][j].status   << ") ";
        // << (*m2_history_particle)[i][j].pt       << ", "
        // << (*m2_history_particle)[i][j].eta      << ", "
        // << (*m2_history_particle)[i][j].phi      << ") ";
      }
      *f << "<--- ";
    }
    *f << std::endl << std::endl;
  }
}


int PythiaNTupleFirstPass::GluonHistoryTracking(int gluon_bar, bool isMuon1){
  // history tracking of a single gluon in the muon history chain
  // determine if a gluon is [from FSR] or [from ISR and of whose Q, Qbar children neither participates in a hard scattering]
  // return: mode (=1 if from FSR, =2 if from ISR)

  // cout << gluon_bar << endl;
  bool is_from_fsr_GUESS = (truth_mother1->at(gluon_bar) < gluon_bar);
  std::vector<int> parent_bars;
  std::vector<int> parent_ids;

  bool prev_step_status_m23_m33 = (abs(truth_status->at(gluon_bar)) == 23 || abs(truth_status->at(gluon_bar)) == 33); // for sanity check

  float pt = truth_pt->at(gluon_bar);
  float eta = truth_eta->at(gluon_bar);
  float phi = truth_phi->at(gluon_bar);
  int id = 21;
  int status = truth_status->at(gluon_bar);
  Particle p {pt, eta, phi, gluon_bar, id, status};

  parent_bars.push_back(gluon_bar);
  parent_ids.push_back(id);
  std::vector<Particle> cur_muon_profile {};
  cur_muon_profile.push_back(p);

  // if (isMuon1){
    // m1_single_gluon_history->push_back(parent_ids);
    // m1_single_gluon_history_particle->push_back(cur_muon_profile);
  // }
  // else{
    // m2_single_gluon_history->push_back(parent_ids);
    // m2_single_gluon_history_particle->push_back(cur_muon_profile);
  // }

  bool terminate = false;
  while(!terminate){
    UpdateCurParents(isMuon1, parent_bars, parent_ids, false, prev_step_status_m23_m33);
    // cout << parent_bars[0] << endl;
    terminate = (abs(truth_status->at(parent_bars[0])) == 21 || truth_mother1->at(parent_bars[0]) <= 2);
  }

  // cout << "terminate" << endl;
  // at this stage, if the gluon is from FSR
  // the hard scattering information (hard_scatt_after_gs & category) has been updated 

  if (abs(truth_status->at(parent_bars[0])) == 21){ // gluon is from FSR
    return 1;
  }

  // gluon from ISR
  return 2;
}


void PythiaNTupleFirstPass::MuonPairTagsReinit(){
  // re-initialize muon pair tags at per-pair level
  // to be called at the beginning of the method MuonPairAncestorTracing()
  m1_c_tag = false;
  m2_c_tag = false;
  m1_osc = false;
  m2_osc = false;
  m1_from_tau = false;
  m2_from_tau = false;

  skip_event_origin_analysis = false;
  // gs_4vec_correctly_set = false;


  from_same_gluon_photon_splitting_or_both_HQ_incoming = false;
  m1_ancestor_is_incoming = false;
  m2_ancestor_is_incoming = false;
  m1_from_hard_scatt_before_gs = false;
  m1_from_hard_scatt_after_gs = false;
  m2_from_hard_scatt_before_gs = false;
  m2_from_hard_scatt_after_gs = false;
  m1_hard_scatt_in_bar1 = -10;
  m2_hard_scatt_in_bar1 = -10;
  m1_eldest_bhadron_barcode = -10;
  m2_eldest_bhadron_barcode = -10;
  m1_hard_scatt_Q = -10.;
  m2_hard_scatt_Q = -10.;

  m1_resonance_barcode = -10;
  m2_resonance_barcode = -10;

  mpair->m1_parent_group = -10;
  mpair->m2_parent_group = -10;

  mpair->Qsplit = -10.;
  mpair->mHard_relevant = -10.;
  // mpair->wrong_4vec_mode_012 = -10;

  mpair->m1_hard_scatt_category = -10;
  mpair->m2_hard_scatt_category = -10;
  mpair->muon_pair_origin_category = -10; // no default "other" origin - need to separate into [non-HF] & [HF uncategorized origin]
  mpair->muon_pair_flavor_category = pair_flavor_index::other_flavor; // initialize to other_flavor if no category found

  mpair->from_same_b = false;
  mpair->from_same_ancestors = false;

  mpair->m1_from_pdf = false;
  mpair->m2_from_pdf = false;

  m1_youngest_non_chadron_parent_barcode = -10;
  m2_youngest_non_chadron_parent_barcode = -10;

  for (std::vector<int> v : *m1_history) v.clear();
  for (std::vector<int> v : *m2_history) v.clear();
  m1_history->clear();
  m2_history->clear();

  for (std::vector<Particle> v : *m1_history_particle) v.clear();
  for (std::vector<Particle> v : *m2_history_particle) v.clear();
  m1_history_particle->clear();
  m2_history_particle->clear();

  // m1_last_hf_hadron_prt_pt_eta_phi_m->clear();
  // m2_last_hf_hadron_prt_pt_eta_phi_m->clear();
  m1_first_hf_hadron_prt_pt_eta_phi_m->clear();
  m2_first_hf_hadron_prt_pt_eta_phi_m->clear();
  // m1_hq_ancestor_pt_eta_phi_m->clear();
  // m2_hq_ancestor_pt_eta_phi_m->clear();

  m1_parton_ancestor_ids.clear();
  m2_parton_ancestor_ids.clear();
  m1_parton_ancestor_bars.clear();
  m2_parton_ancestor_bars.clear();

  m1_multi_hf_quark_ids.clear();
  m2_multi_hf_quark_ids.clear();

  // m1_multi_hadronic_parents_ids.clear();
  // m2_multi_hadronic_parents_ids.clear();
}


void PythiaNTupleFirstPass::FillCategoryHistograms(TH1D* horig, TH2D* hhard_vs_orig){

  horig->Fill(mpair->muon_pair_origin_category, mpair->weight);
  if (m1_hard_scatt_in_bar1 > 0){
    if (m2_hard_scatt_in_bar1 < 0 || m1_hard_scatt_in_bar1 == m2_hard_scatt_in_bar1){ // only one hard or the same hard scattering
      hhard_vs_orig->Fill(mpair->muon_pair_origin_category, mpair->m1_hard_scatt_category, mpair->weight);
    }

  }else if (m2_hard_scatt_in_bar1 > 0){ // need else here: otherwise the [same hard scattering case would get filled twice]
    
    if (m1_hard_scatt_in_bar1 > 0 && (print_bad_warnings && m_very_bad_warning_file && m_hard_scattering_warning_file)){
      *m_very_bad_warning_file << "two cases should be exclusive!!!" << std::endl;
      PrintHistory(m_very_bad_warning_file, false, mpair->from_same_ancestors);
    }

    hhard_vs_orig->Fill(mpair->muon_pair_origin_category, mpair->m2_hard_scatt_category, mpair->weight);
  
  }
}


void PythiaNTupleFirstPass::HFMuonPairAnalysis(){
  // this function finds the origin of a HF muon pair (where both muons come from open HF's)
  // when this function is called, we have already finished tracing of all leptonic, hadronic, and HQ-level history of both open-HF muons
  // for each muon, the vectors m1(2)_parton_ancestor_bars and m1(2)_parton_ancestor_ids record the barcodes + ID's of 
  // the (latest) partonic (non-HQ) parent(s) of the first (earliest) HQ ancestor
  // this method considers (1) the number and identities of these latest non-HQ partonic parents, 
  // (2) whether a hard scattering (HS) is encountered either by the HQ's or by the gluon/photon that splits into HQ's
  // and (3) the order of the partonic-to-HQ transition and the HS
  // Once figured out the pair origin, it then plots this origin category against the relevant hard scattering category

  // ---------------- set pair flavor tag ----------------

  if (m1_youngest_non_chadron_parent_barcode == m2_youngest_non_chadron_parent_barcode){
    int prt_id = abs(truth_id->at(m1_youngest_non_chadron_parent_barcode));
    // necessary to make sure the "hadronic parent" is actually a b-flavored hadron
    // since it's possible that the first not-c-hadron parent vectors is quark level
    // and contains both 4 and -4, where the 4's barcode is recorded as m1_youngest_non_chadron_parent_barcode
    if ((prt_id >= 500 && prt_id < 600) || (prt_id >= 5000 && prt_id < 6000)){

      if (m1_eldest_bhadron_barcode == 10 || m1_eldest_bhadron_barcode != m2_eldest_bhadron_barcode){
        std::cout << "FROM SAME B EVENT WRONG!!! Printing history" << std::endl;
        PrintHistory(&std::cout, false, true);
      }
    }
  }
  
  if (m1_eldest_bhadron_barcode != -10 && m1_eldest_bhadron_barcode == m2_eldest_bhadron_barcode){
    mpair->from_same_b = true;
    mpair->muon_pair_flavor_category = pair_flavor_index::from_single_b;
    mpair->muon_pair_origin_category = muon_pair_both_from_open_HF_origin_catgr::not_both_from_open_and_different_HF_hadrons;
    return;
  }


  // HQ-flavor categories & whether from the same (first non-HQ partonic
  bool both_from_b = ((mpair->m1_parent_group == single_muon_parent_group::direct_b || mpair->m1_parent_group == single_muon_parent_group::b_to_c) && (mpair->m2_parent_group == single_muon_parent_group::direct_b || mpair->m2_parent_group == single_muon_parent_group::b_to_c));
  bool both_from_c = (mpair->m1_parent_group == single_muon_parent_group::direct_c && mpair->m2_parent_group == single_muon_parent_group::direct_c);
  bool one_from_b_one_from_c = (((mpair->m1_parent_group == single_muon_parent_group::direct_b || mpair->m1_parent_group == single_muon_parent_group::b_to_c) && mpair->m2_parent_group == single_muon_parent_group::direct_c) || ((mpair->m2_parent_group == single_muon_parent_group::direct_b || mpair->m2_parent_group == single_muon_parent_group::b_to_c) && mpair->m1_parent_group == single_muon_parent_group::direct_c));
  
  bool same_partonic_ancestors = (m1_parton_ancestor_bars == m2_parton_ancestor_bars);
  mpair->from_same_ancestors = same_partonic_ancestors;

  if (both_from_b) mpair->muon_pair_flavor_category = pair_flavor_index::bb;
  else if (both_from_c) mpair->muon_pair_flavor_category = pair_flavor_index::cc;
  else if (one_from_b_one_from_c) mpair->muon_pair_flavor_category = pair_flavor_index::one_b_one_c;


  bool near_side = (abs(mpair->dphi) < pms.PI / 2.);
  bool small_dphi = (abs(mpair->dphi) < 0.4);
  bool not_near = !(abs(mpair->dphi) < pms.PI / 2.);

  // HF pair origin initialized to "null" - will be overwritten if a specific category is found
  // This way we ensure all pairs both from HF have a valid muon_pair_origin_category given by the enum muon_pair_both_from_open_HF_origin_catgr
  // Avoid confusion + errors in subsequent analysis
  mpair->muon_pair_origin_category = muon_pair_both_from_open_HF_origin_catgr::others;

  // if skip the pair, skip the muon-origin analysis AFTER correctly setting the flavor tag
  // and setting muon_pair_origin_category to (both-from-open-HF) others
  // this is needed since, even if we don't understand the pair, we need to record it to ensure meaning data comparison
  if (skip_event_origin_analysis){
    skipped_event_crossx += mpair->weight;
    return; // return without filling in the ancestor-category histograms
  }

  // ---------------- muon-pair origin categorizing ----------------

  if (same_partonic_ancestors){
    // --------------------------------
    // Case I: the HQ's come from the same partonic ancestors
    // --------------------------------

    if (m1_parton_ancestor_bars.size() > 1){ // from flavor creation (for an LO generator)
      // --------------------------------
      // Case I.1: more than one partonic ancestor --> FC
      // double-check that HS is FC
      // --------------------------------
      if (mpair->m1_hard_scatt_category != flavor_creat && (print_bad_warnings && m_very_bad_warning_file && m_hard_scattering_warning_file)){ // this ensures muon1 has undergone through a hard scattering, and it is a flavor creation
        *m_very_bad_warning_file << "HF muon pair - More than one partonic ancestors, but hard scattering is not flavor creation." << std::endl;
        *m_very_bad_warning_file << "The current parents' barcodes: " << m1_parton_ancestor_bars[0] << " " << m1_parton_ancestor_bars[1] << endl;
        *m_very_bad_warning_file << "The current parents' ids: " << m1_parton_ancestor_ids[0] << " " << m1_parton_ancestor_ids[1] << endl;
        PrintHistory(m_very_bad_warning_file, false, same_partonic_ancestors);
      }
      mpair->muon_pair_origin_category = muon_pair_both_from_open_HF_origin_catgr::fc;
      mpair->mHard_relevant = m1_hard_scatt_Q;
      // mpair->vout1 = vm1out1;
      // mpair->vout2 = vm1out2;

      // print out history for FC events with |Dphi| < 0.4
      if (abs(mpair->dphi) < 0.4 && mpair->m1.ev_num < 10000){ // only print out 1/10
        std::cout << "Flavor creation, not from same b, but |Dphi_{muon pair}| < 0.4" << std::endl;
        std::cout << "How could this be possible?" << std::endl;
        PrintHistory(&std::cout, false, same_partonic_ancestors);
      }

    }else{
      // --------------------------------
      // Case I.2: single parton ancestor --> from gluon/photon splitting / both HQ's incoming 
      // --------------------------------
      if (m1_parton_ancestor_ids[0] != 21 && m1_parton_ancestor_ids[0] != 22 && !m1_ancestor_is_incoming && (print_bad_warnings && m_very_bad_warning_file && m_hard_scattering_warning_file)){
        *m_very_bad_warning_file << "A single parton-level parent that is neither a gluon nor a photon nor incoming." << std::endl;
        PrintHistory(m_very_bad_warning_file, false, same_partonic_ancestors);
      }

      from_same_gluon_photon_splitting_or_both_HQ_incoming = true;

      // sanity checks for gluon/photon splitting
      if (!m1_ancestor_is_incoming){ // excluding the cases where both HQ's are from the PDF
        if (!both_from_b && !both_from_c && (print_bad_warnings && m_very_bad_warning_file && m_hard_scattering_warning_file)){
          *m_very_bad_warning_file << "HF muons from the same gluon splitting, yet neither both from b nor both from c." << std::endl;
          PrintHistory(m_very_bad_warning_file, false, true);
        }

        if ((*m1_history_particle)[m1_history_particle->size() - 2].size() > 2 || (*m2_history_particle)[m2_history_particle->size() - 2].size() > 2 && (print_bad_warnings && m_very_bad_warning_file && m_hard_scattering_warning_file)){
          if (turn_data_resonance_cuts_on){ // if no resonance cuts, pairs from resonances can show up here
            *m_very_bad_warning_file << "The second to last stage in the history of an HF muon (up to the gluon splitting stage) should have at most two particles." << std::endl;
          }
          PrintHistory(m_very_bad_warning_file, false, true);
        }

        int quark_type = (both_from_b)? 5 : 4;
        // float mQ = (both_from_b)? 4.8 : 1.5; // pole mass

        // Particle pg = (*m1_history_particle)[m1_history_particle->size() - 1][0]; // the last element
        int barQ1, barQ2;
        // Particle pQ1;
        // Particle pQ2;

        if ((*m1_history_particle)[m1_history_particle->size() - 2].size() == 1){
          barQ1 = (*m1_history_particle)[m1_history_particle->size() - 2][0].barcode;
          if (abs((*m1_history_particle)[m1_history_particle->size() - 2][0].id) != quark_type){
            if (turn_data_resonance_cuts_on && (print_bad_warnings && m_very_bad_warning_file && m_hard_scattering_warning_file)){  // if no resonance cuts, pairs from resonances can show up here
              *m_very_bad_warning_file << "The second to last stage in the history of an HF muon (up to the gluon splitting stage) has only ONE element" << std::endl;
              *m_very_bad_warning_file << "but that element is NOT a HQ" << std::endl;
              PrintHistory(m_very_bad_warning_file, false, true);              
            }
          }
        }else{ // the second to last stage is the incoming particles to the hard scattering (contains two elements)
          // this should handle all annoyingly error-prone cases, where the incoming-to-hard-scattering particles can be
          // Q + X, Q + Q, or Q + Qbar (X not a HQ of the same flavor)
          if ((*m1_history_particle)[m1_history_particle->size() - 3].size() != 1 || (abs((*m1_history_particle)[m1_history_particle->size() - 3][0].id) != quark_type) && (print_bad_warnings && m_very_bad_warning_file && m_hard_scattering_warning_file)){
            *m_very_bad_warning_file << "The third to last stage in the history of an HF muon (up to the gluon splitting stage) should only have one element, and it should be the HQ of interest" << std::endl;
            PrintHistory(m_very_bad_warning_file, false, true);
          }
          barQ1 = (*m1_history_particle)[m1_history_particle->size() - 3][0].barcode - 2;
          if ((*m1_history_particle)[m1_history_particle->size() - 2][0].barcode != barQ1 && (*m1_history_particle)[m1_history_particle->size() - 2][1].barcode != barQ1 && (print_bad_warnings && m_very_bad_warning_file && m_hard_scattering_warning_file)){
            *m_very_bad_warning_file << "The previous HQ barcode (" << (*m1_history_particle)[m1_history_particle->size() - 3][0].barcode  << ") - 2 method doesn't work." << std::endl;
            PrintHistory(m_very_bad_warning_file, false, true);
          }
        }

        if ((*m2_history_particle)[m2_history_particle->size() - 2].size() == 1){
          barQ2 = (*m2_history_particle)[m2_history_particle->size() - 2][0].barcode;
          if (abs((*m2_history_particle)[m2_history_particle->size() - 2][0].id) != quark_type && (print_bad_warnings && m_very_bad_warning_file && m_hard_scattering_warning_file)){
            *m_very_bad_warning_file << "The second to last stage in the history of an HF muon (up to the gluon splitting stage) has only ONE element" << std::endl;
            *m_very_bad_warning_file << "but that element is NOT a HQ" << std::endl;
            PrintHistory(m_very_bad_warning_file, false, true);
          }
        }else{ // the second to last stage is the incoming particles to the hard scattering (contains two elements)
          // this should handle all annoyingly error-prone cases, where the incoming-to-hard-scattering particles can be
          // Q + X, Q + Q, or Q + Qbar (X not a HQ of the same flavor)
          if ((*m2_history_particle)[m2_history_particle->size() - 3].size() != 1 || (abs((*m2_history_particle)[m2_history_particle->size() - 3][0].id) != quark_type) && (print_bad_warnings && m_very_bad_warning_file && m_hard_scattering_warning_file)){
            *m_very_bad_warning_file << "The third to last stage in the history of an HF muon (up to the gluon splitting stage) should only have one element, and it should be the HQ of interest" << std::endl;
            PrintHistory(m_very_bad_warning_file, false, true);
          }
          barQ2 = (*m2_history_particle)[m2_history_particle->size() - 3][0].barcode - 2;
          if ((*m2_history_particle)[m2_history_particle->size() - 2][0].barcode != barQ2 && (*m2_history_particle)[m2_history_particle->size() - 2][1].barcode != barQ2 && (print_bad_warnings && m_very_bad_warning_file && m_hard_scattering_warning_file)){
            *m_very_bad_warning_file << "The previous HQ barcode (" << (*m2_history_particle)[m2_history_particle->size() - 3][0].barcode  << ") - 2 method doesn't work." << std::endl;
            PrintHistory(m_very_bad_warning_file, false, true);
          }
        }

        if (truth_id->at(barQ1) != -1 * truth_id->at(barQ2)){ // from the same Q but via some complicated mechanism (perhaps at the hadronization step); kinematics likely to be similar to [from same b]
          if (turn_data_resonance_cuts_on && (print_bad_warnings && m_very_bad_warning_file && m_hard_scattering_warning_file)){  // if no resonance cuts, pairs from resonances can show up here
            *m_very_bad_warning_file << "HF muons from the same gluon splitting. Should have quark type +/-" << quark_type;
            *m_very_bad_warning_file << "But actually the id's of the particles are: " << truth_id->at(barQ1) << ", " << truth_id->at(barQ2) << std::endl;
            PrintHistory(m_very_bad_warning_file, false, true); // print out to make sure it's physical not a bug in our code            
          }
        }

        int barg = m1_parton_ancestor_bars[0];
        vg.SetPtEtaPhiM(truth_pt->at(barg), truth_eta->at(barg), truth_phi->at(barg), truth_m->at(barg));
        vQ1.SetPtEtaPhiM(truth_pt->at(barQ1), truth_eta->at(barQ1), truth_phi->at(barQ1), truth_m->at(barQ1));
        vQ2.SetPtEtaPhiM(truth_pt->at(barQ2), truth_eta->at(barQ2), truth_phi->at(barQ2), truth_m->at(barQ2));
        // record the on-shell 4-vectors
        // we'll use an integer-valued flag to indicate which one of these three is wrong
        // mpair->vg = vg;
        // mpair->vQ1 = vQ1;
        // mpair->vQ2 = vQ2;
      }

      // --------------------------------
      // find out if GS/PhS comes from ISR/FSR + if HS is part of the gluon's history
      // --------------------------------
      if (m1_from_hard_scatt_before_gs || m2_from_hard_scatt_before_gs){ // only need one muon to have gone through hard scattering to know the gluon is from ISR
        if (m1_from_hard_scatt_before_gs && m2_from_hard_scatt_before_gs){
          // mpair->muon_pair_origin_category = muon_pair_both_from_open_HF_origin_catgr::same_gs_isr_both_hard_scatt; // if occurs: should be uncorrelated --> classify origin as "others"
          if (mpair->m1.ev_num < nevents_accum[4]/10. || mpair->m1.ev_num > nevents_accum[4]/10. * 9){
            if (turn_data_resonance_cuts_on){ // if no resonance cuts, pairs from resonances can show up here
              std::cout << "For study purpose - same GS ISR both hard scattering" << std::endl;
              std::cout << "Value of pair dphi is " << (mpair->dphi) << std::endl;
              PrintHistory(&std::cout, false, true);              
            }
          }
        }else{
          mpair->muon_pair_origin_category = muon_pair_both_from_open_HF_origin_catgr::same_gs_isr_one_hard_scatt;
        }
      }else if (m1_ancestor_is_incoming){ // same ancestor --> both HQ's are from PDF (of the same nucleus)
        // std::cout << "For study purpose - how two HQ's from PDF could each produce a 4GeV pT muon without either going through a hard scattering." << std::endl;
        // PrintHistory(&std::cout, false, true);
        mpair->muon_pair_origin_category = muon_pair_both_from_open_HF_origin_catgr::same_gs_isr_zero_hard_scatt;
      }else{ // determine history & origin of the gluon
        int gluon_mode = GluonHistoryTracking(m1_parton_ancestor_bars[0], true);
        // at this point, either gluon is from FSR (so "same hard scattering" for m1, m2), or is from ISR yet neither HQ's experiences any hard scattering before producing the two muons
        // in both cases, muon2 have the same hard scattering category as muon1
        mpair->m2_hard_scatt_category = mpair->m1_hard_scatt_category; // m2 hard scattering info is not filled since we only called GluonHistoryTracking with muon1
        if (gluon_mode == 1){ // gluon from FSR
          if (m1_parton_ancestor_ids[0] == 21){
            mpair->muon_pair_origin_category = muon_pair_both_from_open_HF_origin_catgr::same_gs_fsr;
          } else if (m1_parton_ancestor_ids[0] == 22){
            mpair->muon_pair_origin_category = muon_pair_both_from_open_HF_origin_catgr::same_phs_fsr;
          }
        }else{
          mpair->muon_pair_origin_category = muon_pair_both_from_open_HF_origin_catgr::same_gs_isr_zero_hard_scatt;
          if (mpair->m1.ev_num < nevents_accum[4]/100. || mpair->m1.ev_num > nevents_accum[4]/100. * 99){
            // std::cout << "For study purpose - same GS ISR zero hard scattering and NOT from PDF." << std::endl;
            // std::cout << "Value of dphi is " << (mpair->dphi) << std::endl;
            // PrintHistory(&std::cout, false, true);
          }
        }
      }

      // --------------------------------
      // In case of gluon/photon splitting, store the 4-vectors vg, vQ1, vQ2
      // --------------------------------
      if (!m1_ancestor_is_incoming){ // excluding the cases where both HQ's are from the PDF --> gluon/photon splitting
        // the three on-shell 4-vectors have been recorded
        // now find the correct Q_split depending on whether it's ISR or FSR
        // and record which 4 vector is the most off-shell (to be corrected in future code if we need 4-momentum conservation)
        if (mpair->muon_pair_origin_category == muon_pair_both_from_open_HF_origin_catgr::same_gs_fsr || mpair->muon_pair_origin_category == muon_pair_both_from_open_HF_origin_catgr::same_phs_fsr ){ // GS from FSR: Q^2 is the virtuality of the gluon (internal propagator)
          // mpair->wrong_4vec_mode_012 = 0;
          vg = vQ1 + vQ2;
          mpair->Qsplit = vg.M();

          if (m1_hard_scatt_Q < 0 && (print_bad_warnings && m_very_bad_warning_file && m_hard_scattering_warning_file)){
            *m_very_bad_warning_file << "GS from FSR, but Q of hard scattering is unfilled or negative, with value " << m1_hard_scatt_Q << std::endl;
            PrintHistory(m_very_bad_warning_file, false, true);
          }

          mpair->mHard_relevant = m1_hard_scatt_Q;
          // mpair->vout1 = vm1out1;
          // mpair->vout2 = vm1out2;

          h_Qsplit_gs_FSR[!mpair->same_sign][not_near]->Fill(abs(vg.M()), mpair->weight);
          h_Qsplit_to_mHard_gs_FSR[!mpair->same_sign][not_near]->Fill(abs(mpair->Qsplit) / mpair->mHard_relevant, mpair->weight);

        }else if(mpair->muon_pair_origin_category == same_gs_isr_one_hard_scatt){
          // TLorentzVector vQQ_sum = vQ1 + vQ2;
          // mpair->mQQ = vQQ_sum.M(); // this will simply give zero since the gluon is recorded as if on shell

          // now, the correct off-shell treatment
          if (m1_hard_scatt_in_bar1 > 0){ // Q1 is the most off-shell
          // mpair->wrong_4vec_mode_012 = 1;
            vQ1 = vg - vQ2; // explicitly conserve 4-momentum
            mpair->Qsplit = vQ1.M();
            mpair->mHard_relevant = m1_hard_scatt_Q;
            // mpair->vout1 = vm1out1;
            // mpair->vout2 = vm1out2;
          }else{ // Q2 is the most off-shell
          // mpair->wrong_4vec_mode_012 = 2;
            vQ2 = vg - vQ1;
            mpair->Qsplit = vQ2.M();
            mpair->mHard_relevant = m2_hard_scatt_Q;
          }

          h_Qsplit_gs_ISR_one_hard_scatt[!mpair->same_sign][not_near]->Fill(abs(mpair->Qsplit), mpair->weight);
          h_Qsplit_to_mHard_gs_ISR_one_hard_scatt[!mpair->same_sign][not_near]->Fill(abs(mpair->Qsplit) / mpair->mHard_relevant, mpair->weight);
          // h_mQQ_to_Qsplit_gs_ISR_one_hard_scatt->Fill(mpair->mQQ / mpair->Qsplit, mpair->weight);
        }
      }
    }
  }else{
    // --------------------------------
    // Case II: the HQ's come from the same partonic ancestors
    // --------------------------------
    if (m1_from_hard_scatt_before_gs && m1_hard_scatt_in_bar1 == m2_hard_scatt_in_bar1 && mpair->m1_hard_scatt_category != hq_scatt){      
      std::cout << "Not same partonic ancestors. Same hard scattering, and it's not Q Q' scattering." << std::endl;
      // PrintHistory(&std::cout, false, same_partonic_ancestors);
    }
    if (!m1_from_hard_scatt_before_gs){ // m1 not yet experiened hard scattering      
      if (!m1_parton_ancestor_bars.empty()){
        if ((m1_parton_ancestor_bars.size() != 1 || (m1_parton_ancestor_ids[0] != 21 && m1_parton_ancestor_ids[0] != 22 && !m1_ancestor_is_incoming)) && (print_bad_warnings && m_very_bad_warning_file && m_hard_scattering_warning_file)){
          *m_very_bad_warning_file << "Muon hasn't gone through a hard scattering; yet is neither from gluon/photon splitting or from PDF." << std::endl;
          PrintHistory(m_very_bad_warning_file, true, true); // single, m1
        }else if (m1_ancestor_is_incoming){
          // std::cout << "For study purpose - Muon1 is from PDF and have not gone through hard scattering." << std::endl;
          // PrintHistory(&std::cout, true, true); // single, m1
          // mpair->muon_pair_origin_category = muon_pair_both_from_open_HF_origin_catgr::others;
        }else{ // determine if [GS from FSR] or [GS from ISR & no hard scattering]
          GluonHistoryTracking(m1_parton_ancestor_bars[0], true);
        }
      }        
    }

    if (!m2_from_hard_scatt_before_gs){ // m2 not yet experiened hard scattering      
      if (!m2_parton_ancestor_bars.empty()){
        if ((m2_parton_ancestor_bars.size() != 1 || (m2_parton_ancestor_ids[0] != 21 && m2_parton_ancestor_ids[0] != 22 && !m2_ancestor_is_incoming)) && (print_bad_warnings && m_very_bad_warning_file && m_hard_scattering_warning_file)){
          *m_very_bad_warning_file << "Muon hasn't gone through a hard scattering; yet is neither from gluon/photon splitting or from PDF." << std::endl;
          PrintHistory(m_very_bad_warning_file, true, false); // single, m2
        }else if (m2_ancestor_is_incoming){
          // std::cout << "For study purpose - Muon2 is from PDF and have not gone through hard scattering." << std::endl;
          // PrintHistory(&std::cout, true, false); // single, m1
          // mpair->muon_pair_origin_category = muon_pair_both_from_open_HF_origin_catgr::others;
        }else{
          GluonHistoryTracking(m2_parton_ancestor_bars[0], false);
        }
      }
    }

    // now that we've updated the hard scattering information. Let's see if they are from the same hard scattering
    if (m1_hard_scatt_in_bar1 == m2_hard_scatt_in_bar1 && m1_hard_scatt_in_bar1 > 0){      
      mpair->muon_pair_origin_category = diff_gs_same_hard_scatt;
      mpair->mHard_relevant = m1_hard_scatt_Q;
    }
    
    if (((m1_ancestor_is_incoming && !m1_from_hard_scatt_before_gs) || (m2_ancestor_is_incoming && !m2_from_hard_scatt_before_gs)) && (mpair->muon_pair_origin_category != muon_pair_both_from_open_HF_origin_catgr::others)){      
      std::cout << "Different ancestors - at least one satisfies: (1) is incoming (2) never experienced hard scattering." << std::endl;
      std::cout << "The muons should not be correlated - pair origin category should be [others]." << std::endl;
      PrintHistory(&std::cout, false, false);
    }
  }
  
  // print history for others origin category if minv is very low
  if (mpair->muon_pair_origin_category == muon_pair_both_from_open_HF_origin_catgr::others){
    if (print_HF_pair_origin_others_history && m_HF_pair_origin_others_category_file && mpair->minv < low_minv_threshold){
      *m_HF_pair_origin_others_category_file << "pair minv: " << mpair->minv << endl;
      PrintHistory(m_HF_pair_origin_others_category_file, false, mpair->from_same_ancestors, true); // print both muons history + print hard-scattering category
    }
  }

  // ---------------- Fill in origin + HS category histograms ----------------

  // analysis of muon pairs both from b & NOT from the same b
  if (both_from_b){ // have excluded the [from same b] possibility
    
    // fill in histograms
    FillCategoryHistograms(h_both_from_b_muon_pair_origin_categr[!mpair->same_sign][not_near], h_both_from_b_hard_scatt_categr_vs_origin_categr[!mpair->same_sign][not_near]);
    FillCategoryHistograms(h_muon_pair_origin_categr[!mpair->same_sign][not_near], h_hard_scatt_categr_vs_origin_categr[!mpair->same_sign][not_near]);

    // print history
    if (print_prt_history && mpair->m1.ev_num < 5000){
      if (same_partonic_ancestors)   PrintHistory(m_b_parent_file[!mpair->same_sign][not_near], false, true, true);
      else                  PrintHistory(m_b_parent_file[!mpair->same_sign][not_near], false, false, true);
    }
  }
  if (one_from_b_one_from_c){
    FillCategoryHistograms(h_one_from_b_one_from_c_muon_pair_origin_categr[!mpair->same_sign][not_near], h_one_from_b_one_from_c_hard_scatt_categr_vs_origin_categr[!mpair->same_sign][not_near]);
    FillCategoryHistograms(h_muon_pair_origin_categr[!mpair->same_sign][not_near], h_hard_scatt_categr_vs_origin_categr[!mpair->same_sign][not_near]); 
  }
  // case of c-cbar same-sign muon pairs both from c (directly)
  if (both_from_c){
    FillCategoryHistograms(h_both_from_c_muon_pair_origin_categr[!mpair->same_sign][not_near], h_both_from_c_hard_scatt_categr_vs_origin_categr[!mpair->same_sign][not_near]);
    FillCategoryHistograms(h_muon_pair_origin_categr[!mpair->same_sign][not_near], h_hard_scatt_categr_vs_origin_categr[!mpair->same_sign][not_near]);

    // print history for inclusive open-HF pairs if printing option turned on
    if (print_prt_history && mpair->m1.ev_num < 50000){
      if (same_partonic_ancestors)   PrintHistory(m_c_parent_file[!mpair->same_sign][not_near], false, true, true);
      else                  PrintHistory(m_b_parent_file[!mpair->same_sign][not_near], false, false, true);
    }
  }
}

void PythiaNTupleFirstPass::MuonPairAncestorTracing(){
  MuonPairTagsReinit();
  bool not_near = !(abs(mpair->dphi) < pms.PI / 2.);

  // implement mpair->m1_parent_group, mpair->m2_parent_group within SingleMuonAncestorTracing function call
  SingleMuonAncestorTracing(true);
  SingleMuonAncestorTracing(false);
  
  mpair->from_same_resonance = (mpair->m1_parent_group == resonance_decay && mpair->m2_parent_group == resonance_decay && m1_resonance_barcode == m2_resonance_barcode);
  if (mpair->from_same_resonance) mpair->muon_pair_flavor_category = pair_flavor_index::from_resonance;

  if (mpair->from_same_resonance && mpair->minv < low_minv_threshold && print_low_minv_resonances && m_very_low_minv_resonance_file){
    // *m_very_low_minv_resonance_file << "Pair from same resonance but minv = " << mpair->minv << std::endl;
    PrintHistory(m_very_low_minv_resonance_file, false, true);
  }

  mpair->resonance_contaminated = (!mpair->from_same_resonance && (mpair->m1_parent_group == resonance_decay || mpair->m2_parent_group == resonance_decay));
  if (mpair->resonance_contaminated) mpair->muon_pair_flavor_category = pair_flavor_index::resonance_contaminated;

  if (mpair->m1_parent_group == single_muon_parent_group::single_photon && mpair->m2_parent_group == single_muon_parent_group::single_photon){
    mpair->muon_pair_flavor_category = pair_flavor_index::photon_to_dimuon_splitting;
  }

  if (mpair->m1_parent_group == prt_drell_yan && mpair->m2_parent_group == prt_drell_yan){
    mpair->muon_pair_flavor_category = pair_flavor_index::drell_yan;
  }

  mpair->m1_from_pdf = m1_ancestor_is_incoming;
  mpair->m2_from_pdf = m2_ancestor_is_incoming;

  // cout << "find single muon parents fine" << endl;

  // if (mpair->m1_parent_group < 0) *m_unspecified_parent_file << "*** Ungrouped parent: Event#" << mpair->m1.ev_num << ", parent id: " << m1_earliest_parent_id << " ***" << std::endl << std::endl;
  // if (mpair->m2_parent_group < 0) *m_unspecified_parent_file << "*** Ungrouped parent: Event#" << mpair->m1.ev_num << ", parent id: " << m2_earliest_parent_id << " ***" << std::endl << std::endl;
  
  if (mpair->m1_parent_group < 0){
    std::cout << "parents not in any group: " << std::endl;
    PrintHistory(&std::cout, true, true);
  }
  if (mpair->m2_parent_group < 0){
    std::cout << "parents not in any group: " << std::endl;
    PrintHistory(&std::cout, true, true);
  }

  h_parent_groups[!mpair->same_sign][not_near]->Fill(mpair->m1_parent_group + 0.5, mpair->m2_parent_group + 0.5, mpair->weight);

  mpair->pair_origin_analysis_skipped = skip_event_origin_analysis;

  if ((mpair->m1_parent_group == single_muon_parent_group::direct_b || mpair->m1_parent_group == single_muon_parent_group::b_to_c || mpair->m1_parent_group == single_muon_parent_group::direct_c) && (mpair->m2_parent_group == single_muon_parent_group::direct_b || mpair->m2_parent_group == single_muon_parent_group::b_to_c || mpair->m2_parent_group == single_muon_parent_group::direct_c)){    
    HFMuonPairAnalysis();
  }else{
    // 
    mpair->muon_pair_origin_category = muon_pair_both_from_open_HF_origin_catgr::not_both_from_open_and_different_HF_hadrons;
  }

  // print history for others flavor if minv is very low
  if (mpair->muon_pair_flavor_category == pair_flavor_index::other_flavor){
    if (print_other_flavor_history && m_other_flavor_category_file && mpair->minv < low_minv_threshold){
      *m_other_flavor_category_file << "others flavor with low minv" << endl;
      *m_other_flavor_category_file << "pair minv: " << mpair->minv << endl;
      PrintHistory(m_other_flavor_category_file, false, mpair->from_same_ancestors); // print both muons history
    }
  }
}

// void PythiaNTupleFirstPass::KinematicCorrPlots(int isign, int idphi){
//   h_pt_muon_pt_closest_hadr_ratio[isign][idphi]->Fill(mpair->m1.pt / (*m1_last_hf_hadron_prt_pt_eta_phi_m)[0],mpair->weight);
//   h_pt_muon_pt_closest_hadr_ratio[isign][idphi]->Fill(mpair->m2.pt / (*m2_last_hf_hadron_prt_pt_eta_phi_m)[0],mpair->weight);
//   h_pt_closest_hadr_pt_furthest_hadr_ratio[isign][idphi]->Fill((*m1_last_hf_hadron_prt_pt_eta_phi_m)[0] / (*m1_first_hf_hadron_prt_pt_eta_phi_m)[0], mpair->weight);
//   h_pt_closest_hadr_pt_furthest_hadr_ratio[isign][idphi]->Fill((*m2_last_hf_hadron_prt_pt_eta_phi_m)[0] / (*m2_first_hf_hadron_prt_pt_eta_phi_m)[0], mpair->weight);
//   h_pt_hadr_hq_ratio[isign][idphi]->Fill((*m1_last_hf_hadron_prt_pt_eta_phi_m)[0] / (*m1_hq_ancestor_pt_eta_phi_m)[0], mpair->weight);
//   h_pt_hadr_hq_ratio[isign][idphi]->Fill((*m2_last_hf_hadron_prt_pt_eta_phi_m)[0] / (*m2_hq_ancestor_pt_eta_phi_m)[0], mpair->weight);
//   float muon_hadr_dphi = mpair->m1.phi - (*m1_last_hf_hadron_prt_pt_eta_phi_m)[2];
//   muon_hadr_dphi = atan2(sin(muon_hadr_dphi),cos(muon_hadr_dphi));//fold muon_hadr_dphi to [-pi,pi]
//   h_dphi_muon_closest_hadr[isign][idphi]->Fill(abs(muon_hadr_dphi), mpair->weight);
//   muon_hadr_dphi = mpair->m2.phi - (*m2_last_hf_hadron_prt_pt_eta_phi_m)[2];
//   muon_hadr_dphi = atan2(sin(muon_hadr_dphi),cos(muon_hadr_dphi));//fold muon_hadr_dphi to [-pi,pi]
//   h_dphi_muon_closest_hadr[isign][idphi]->Fill(abs(muon_hadr_dphi), mpair->weight);
// }


void PythiaNTupleFirstPass::FillMuonPairTreePythia(int nkin){

  // NECESSARY step: ALWAYS update the raw pointer with the current content of the shared pointer BEFORE FILLING OUTPUT TREES
  try{
    if (!mpair) throw std::runtime_error("When trying to fill current muon pair into output trees, the shared_ptr mpair is nullptr!");
    mpair_raw_ptr = mpair.get();
  }catch(const std::runtime_error& e){
    std::cout << "Runtime error caught in function FillMuonPairTree: " << e.what() << std::endl;
    std::cout << "Return without filling the muon pair in the output trees!" << std::endl;
    return;
  }

  int nsign = (mpair->same_sign)? 0:1;
  muonPairOutTree[nsign]->Fill();
  muonPairOutTreeKinRange[nkin][nsign]->Fill();
}


void PythiaNTupleFirstPass::HistAdjust(){
  DimuonAnalysisBaseClass::HistAdjust();

  // for (int ibin = 0; ibin < 3; ibin++){
  //   for (int isign = 0; isign < ParamsSet::nSigns; isign++){
  //     for (int jdphi = 0; jdphi < 2; jdphi++){
  //       h_num_hard_scatt_out[isign][jdphi]->GetXaxis()->SetBinLabel(ibin+1,num_hard_scatt_out_labels[ibin].c_str());
  //     }
  //   }
  // }

  for (int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
    for (int lphi = 0; lphi < 2; lphi++){
      for (int iprt = 0; iprt < nParentGroups; iprt++){
        h_parent_groups[ksign][lphi]->GetXaxis()->SetBinLabel(iprt+1,parentGroupLabels[iprt].c_str());
        h_parent_groups[ksign][lphi]->GetYaxis()->SetBinLabel(iprt+1,parentGroupLabels[iprt].c_str());
      }

      for (int iorigin = 0; iorigin < nHFMuonPairOriginCategories; iorigin++){
        h_muon_pair_origin_categr[ksign][lphi]->GetXaxis()->SetBinLabel(iorigin+1,HFMuonPairOriginCategoryLabels[iorigin].c_str());
        h_both_from_b_muon_pair_origin_categr[ksign][lphi]->GetXaxis()->SetBinLabel(iorigin+1,HFMuonPairOriginCategoryLabels[iorigin].c_str());
        h_both_from_c_muon_pair_origin_categr[ksign][lphi]->GetXaxis()->SetBinLabel(iorigin+1,HFMuonPairOriginCategoryLabels[iorigin].c_str());
        h_one_from_b_one_from_c_muon_pair_origin_categr[ksign][lphi]->GetXaxis()->SetBinLabel(iorigin+1,HFMuonPairOriginCategoryLabels[iorigin].c_str());
        
        h_hard_scatt_categr_vs_origin_categr[ksign][lphi]->GetXaxis()->SetBinLabel(iorigin+1,HFMuonPairOriginCategoryLabels[iorigin].c_str());
        h_both_from_b_hard_scatt_categr_vs_origin_categr[ksign][lphi]->GetXaxis()->SetBinLabel(iorigin+1,HFMuonPairOriginCategoryLabels[iorigin].c_str());
        h_both_from_c_hard_scatt_categr_vs_origin_categr[ksign][lphi]->GetXaxis()->SetBinLabel(iorigin+1,HFMuonPairOriginCategoryLabels[iorigin].c_str());
        h_one_from_b_one_from_c_hard_scatt_categr_vs_origin_categr[ksign][lphi]->GetXaxis()->SetBinLabel(iorigin+1,HFMuonPairOriginCategoryLabels[iorigin].c_str());
      }
      for (int ihard = 0; ihard < nHardScattCategories; ihard++){
        h_hard_scatt_categr_vs_origin_categr[ksign][lphi]->GetYaxis()->SetBinLabel(ihard+1,HardScattCategoryLabels[ihard].c_str());
        h_both_from_b_hard_scatt_categr_vs_origin_categr[ksign][lphi]->GetYaxis()->SetBinLabel(ihard+1,HardScattCategoryLabels[ihard].c_str());
        h_both_from_c_hard_scatt_categr_vs_origin_categr[ksign][lphi]->GetYaxis()->SetBinLabel(ihard+1,HardScattCategoryLabels[ihard].c_str());
        h_one_from_b_one_from_c_hard_scatt_categr_vs_origin_categr[ksign][lphi]->GetYaxis()->SetBinLabel(ihard+1,HardScattCategoryLabels[ihard].c_str());
      
      }
    }
  }
}


void PythiaNTupleFirstPass::CrossxClear(){
  pair_counter = 0;
  total_crossx = 0.;
  from_same_b_total_crossx = 0.;
  either_from_tau_total_crossx = 0.;
  both_incoming_total_crossx = 0.;
  both_incoming_FE_QQ_from_same_b_total_crossx = 0.;
  FE_total_crossx = 0.;
  from_same_gluon_spitting_total_crossx = 0.;
  hard_QQ_scatt_total_crossx = 0.;
  FE_from_same_b_total_crossx = 0.;
  FE_from_same_GS_total_crossx = 0.;
  FE_from_diff_ancestors_total_crossx = 0.;
  FE_from_same_ancestors_not_same_b_or_gs_total_crossx = 0.;
}

void PythiaNTupleFirstPass::PerPairCrossxUpdate(){
  total_crossx += mpair->weight;

  bool parton1_heavy_quark = (abs(truth_id->at(3)) == 4 || abs(truth_id->at(3)) == 5);
  bool parton2_heavy_quark = (abs(truth_id->at(4)) == 4 || abs(truth_id->at(4)) == 5);
  bool parton1_g_lightq = (abs(truth_id->at(3)) == 1 || abs(truth_id->at(3)) == 2 || abs(truth_id->at(3)) == 3 || abs(truth_id->at(3)) == 21);
  bool parton2_g_lightq = (abs(truth_id->at(4)) == 1 || abs(truth_id->at(4)) == 2 || abs(truth_id->at(4)) == 3 || abs(truth_id->at(4)) == 21);
  bool flavor_excitation = ((parton1_heavy_quark && parton2_g_lightq) || (parton2_heavy_quark && parton1_g_lightq));
  bool hard_QQ_scatt = (parton1_heavy_quark && parton2_heavy_quark);
  bool both_incoming = (m1_ancestor_is_incoming && m2_ancestor_is_incoming);

  if (m1_from_tau || m2_from_tau){
    either_from_tau_total_crossx += mpair->weight;
  }

  if (hard_QQ_scatt){
    hard_QQ_scatt_total_crossx += mpair->weight;
  }

  if (flavor_excitation){
    FE_total_crossx += mpair->weight;
    if (mpair->from_same_b){ // not specify where the b is from, since effects from the origin of the b is kind of washed out by the sequential decay
      FE_from_same_b_total_crossx += mpair->weight;
    }else if (from_same_gluon_photon_splitting_or_both_HQ_incoming){
      FE_from_same_GS_total_crossx += mpair->weight;
    }else if (!mpair->from_same_ancestors){
      FE_from_diff_ancestors_total_crossx += mpair->weight;
    }else{ // same ancestors, not from same b, and not from the same gluon splitting
      FE_from_same_ancestors_not_same_b_or_gs_total_crossx += mpair->weight;
      if (print_FE && m_FE_file){
        *m_FE_file << "Event #: " << mpair->m1.ev_num << ". Hard scattering is FE, but muon pairs are from same ancestors and not same b or gluon splitting." << std::endl;
        *m_FE_file << truth_id->at(3) << " " << truth_id->at(4) << " " << truth_id->at(5) << " " << truth_id->at(6) << " " << mpair->from_same_ancestors << std::endl;
        PrintHistory(m_FE_file, false, mpair->from_same_ancestors);              
      }
    }
  }

  if (mpair->from_same_b){
    from_same_b_total_crossx += mpair->weight;
  }

  if (from_same_gluon_photon_splitting_or_both_HQ_incoming){
    from_same_gluon_spitting_total_crossx += mpair->weight;
  }

  if (both_incoming) {
    both_incoming_total_crossx += mpair->weight;

    if ((flavor_excitation && mpair->from_same_b) || hard_QQ_scatt){
      both_incoming_FE_QQ_from_same_b_total_crossx += mpair->weight;
    }
    else{
      // std::cout << "Event #: " << mpair->m1.ev_num << ". Both muons are from quarks that are incoming partons, and not FE + from the same b." << std::endl;
      // // std::cout << truth_id->at(3) << " " << truth_id->at(4) << " " << truth_id->at(5) << " " << truth_id->at(6) << " " << mpair->from_same_b << std::endl;

      // PrintHistory(&std::cout, false, mpair->from_same_ancestors);
    }
  }
}

void PythiaNTupleFirstPass::WriteCrossxSummary(){
  *m_crossx_summary_file << "#muon pairs: " << pair_counter << std::endl;
  *m_crossx_summary_file << "Total crossx: " << total_crossx << std::endl;
  *m_crossx_summary_file << "Total crossx of skipped HF events (not filled in HF-muon-pair-ancestor-category trees): " << skipped_event_crossx << std::endl;
  *m_crossx_summary_file << "Total crossx [from resonance]: " << from_resonance_total_crossx << std::endl;
  *m_crossx_summary_file << "Total crossx [from same b]: " << from_same_b_total_crossx << std::endl;
  *m_crossx_summary_file << "Total crossx [BOTH INCOMING]: " << both_incoming_total_crossx << std::endl;
  *m_crossx_summary_file << "Total crossx [both incoming] - from FE & same b: " << both_incoming_FE_QQ_from_same_b_total_crossx << std::endl;
  *m_crossx_summary_file << "Total crossx [same GS] (gluon splitting): " << from_same_gluon_spitting_total_crossx << std::endl;

  *m_crossx_summary_file << "Total crossx [QQ']: " << hard_QQ_scatt_total_crossx << std::endl;
  *m_crossx_summary_file << "Total crossx [FE] (hard scattering is flavor excitation): " << FE_total_crossx << std::endl;
  *m_crossx_summary_file << "Out of such:" << std::endl;
  *m_crossx_summary_file << "Total crossx [FE - from same b] (b any origin): " << FE_from_same_b_total_crossx << std::endl;
  *m_crossx_summary_file << "Total crossx [FE - from same GS]: " << FE_from_same_GS_total_crossx << std::endl;
  *m_crossx_summary_file << "Total crossx [FE - from different ancestores]: " << FE_from_diff_ancestors_total_crossx << std::endl;
  *m_crossx_summary_file << "Total crossx [FE - others]: " << FE_from_same_ancestors_not_same_b_or_gs_total_crossx << std::endl;
  *m_crossx_summary_file << "Total crossx [either is from tau]: " << either_from_tau_total_crossx << std::endl;
  
  *m_crossx_summary_file << "For pairs from resonance decay:" << std::endl;
  for (const auto& resonance : resonance_id_to_name_and_crossx_map){
    *m_crossx_summary_file << "crossx [" << resonance.second.first << " decay]: " << resonance.second.second << std::endl;
  }
}

void PythiaNTupleFirstPass::ProcessData(){

  // qqpair = new TruthQQPair(quark);

  CrossxClear();

  for (int ikin = 0; ikin < nKinRanges; ikin++){
  // for (int ikin = 3; ikin <= 3; ikin++){

    int nevent_start = (ikin == 0)? 0 : nevents_accum[ikin - 1];
    int njob_start = (ikin == 0)? 0 : njobs_accum[ikin - 1];
    std::cout << "Processsing Kinematic Range: " << ikin << ", Entries from " << nevent_start << " to " << nevents_accum[ikin]-1 << std::endl;
  
    for (Long64_t jevent = nevent_start; jevent < nevents_accum[ikin]; jevent++) {//loop over the events
    // for (Long64_t jevent = nevent_start; jevent < nevent_start + 10; jevent++) {//loop over the events // debug mode

      // cout << "event" << jevent << "***********" << endl;

      // cout << "processing event" << jevent << endl;
      int kjob = njob_start + static_cast<int>((jevent - nevent_start) / nevents_per_file[ikin]);
      int num_bytes_ev = evChain->GetEntry(jevent);//read in an event
      int num_bytes_job = metaChain->GetEntry(kjob);//read in an event
      if(num_bytes_ev==0){
          std::cout<<"Error:: Read in event has size of zero bytes, quitting"<<std::endl;
          throw std::exception();
      }
      if(num_bytes_job==0){
          std::cout<<"Error:: Read in job has size of zero bytes, quitting"<<std::endl;
          throw std::exception();
      }
      if (jevent == nevents_accum[ikin] - 1 && kjob != njobs_accum[ikin] - 1){
          std::cout<<"Error:: number of jobs not matching expectation, quitting"<<std::endl;
          throw std::exception();
      }

      muon_pair_list_cur_event_pre_resonance_cut.clear();
      resonance_tagged_muon_index_list_v2.clear(); // MUST CLEAR for each event!!
      resonance_tagged_muon_index_list.clear(); // MUST CLEAR for each event!!

      int NPairs = muon_pair_muon1_pt->size();//number of muon pairs in the event
      int NPairsAfter = 0;

      for(int pair_ind=0;pair_ind<NPairs;pair_ind++){//first loop over all muon-pairs in the event

        mpair = std::make_shared<MuonPairPythia>(MuonPairPythia());

        FillMuonPair(pair_ind, mpair);

        mpair->m1.ev_num = jevent;
        mpair->m2.ev_num = jevent;

        // mpair->weight     = ev_weight / njobs[ikin]; // old version: individual muon-pairs output files can be used by themselves
                                                        // but combining them for histogram filling requires reweighting
        mpair->weight     = ev_weight / njobs_all_files_combined[ikin]; // new version: MUST use combined muon-pairs output file
                                                                        // the individual files will have wrong weights
        mpair->m1.ev_weight = mpair->weight;
        mpair->m2.ev_weight = mpair->weight;
        mpair->crossx     = mpair->weight * nevents[ikin] / efficiency;

        // ------------------------------------------------------------

        h_cutAcceptance[mpair->m1.charge != mpair->m2.charge]->Fill(nocut + 0.5, mpair->weight);


        // Apply cuts
        if (!PassCuts(mpair))continue;
              
        //------------------------------------------------------------

        //Two things at this step: 
        //1) sort pt, eta, phi by pt
        //2) update the muon-pair values
        mpair->Update();
          
        // resonance tag
        ResonanceTagging(mpair);
        ResonanceTaggingV2(mpair);

        // photo-production cut - do NOT apply for MC
        // if (IsPhotoProduction()) continue;
              
        muon_pair_list_cur_event_pre_resonance_cut.push_back(std::move(mpair));
        
      } // finish first loop over all muon pairs

      for(int pair_ind=0;pair_ind<muon_pair_list_cur_event_pre_resonance_cut.size();pair_ind++){//second loop over all muon-pairs in the event
        // discard the pair if either muon is resonance-tagged
        mpair = std::move(muon_pair_list_cur_event_pre_resonance_cut.at(pair_ind));

        if (!mpair){
          std::cerr << "mpair at second muon-pair loop NOT found! Return without resonance-tag checking or pair analysis." << std::endl;
          continue;
        }

        mpair->data_resonance_or_reso_contam_tagged_old = false;
        
        // examine data resonance cuts (old + new)
        std::vector<int>::iterator itres_m1;
        std::vector<int>::iterator itres_m2;

        // check for new data-resonance cut
        // if turn data cuts on: throw away pairs tagged as "from resonance" / "contaminated by resonances" under new data-resonance cuts
        // else: tag the pairs + keep them
        itres_m1 = std::find(resonance_tagged_muon_index_list.begin(),resonance_tagged_muon_index_list.end(),mpair->m1.ind);
        itres_m2 = std::find(resonance_tagged_muon_index_list.begin(),resonance_tagged_muon_index_list.end(),mpair->m2.ind);
        
        if(itres_m1 != resonance_tagged_muon_index_list.end() || itres_m2 != resonance_tagged_muon_index_list.end()){ // pair is resonance-tagged
          if (turn_data_resonance_cuts_on){
            continue; // throw away the pair
          }else{
            mpair->data_resonance_or_reso_contam_tagged_old = true;
          }
        }

        // check for data-resonance cut v2 - no pair elimination, just tag
        itres_m1 = std::find(resonance_tagged_muon_index_list_v2.begin(),resonance_tagged_muon_index_list_v2.end(),mpair->m1.ind);
        itres_m2 = std::find(resonance_tagged_muon_index_list_v2.begin(),resonance_tagged_muon_index_list_v2.end(),mpair->m2.ind);

        mpair->data_resonance_or_reso_contam_tagged_new = (itres_m1 != resonance_tagged_muon_index_list_v2.end() || itres_m2 != resonance_tagged_muon_index_list_v2.end());

        h_cutAcceptance[mpair->m1.charge != mpair->m2.charge]->Fill(pass_resonance + 0.5, mpair->weight);
        
        //------------------------------------------------------------
        // the pair survives all cuts --> proceed with analysis + fill the muon pair in the output trees

        // fill in muon parents (now that we have finished applying all the cuts)
        MuonPairAncestorTracing();
        FillMuonPairTreePythia(ikin); // fill all-muon-pair tree & kinematic-range-binned trees


        PerPairCrossxUpdate();


        NPairsAfter++;        
        h_numMuonPairs->Fill(NPairsAfter - 0.5, mpair->weight);


      } // finish second loop over all muon pair candidates in the current event
    } // loop over events in the current kinematic range
  } // loop over all kinematic ranges
  
  // delete qqpair;
  WriteCrossxSummary();
}


void PythiaNTupleFirstPass::Run(){
  clock_t start, end;
  double cpu_time_used;
  start = clock();
  
  ResonanceNameMap();

  SetInputOutputFilesFromBatch(); // must be called before InitInput and InputOutput

  InputSanityCheck();
  InitInput();
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  std::cout << "Took " << cpu_time_used << " seconds to initiate the inputs." << std::endl;

  InitTempVariables();
  InitOutput();
  
  meta_tree_out->Fill();

  start = clock();
  ProcessData();
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  std::cout << "Took " << cpu_time_used << " seconds to run the event loop." << std::endl;
  
  HistAdjust();

  m_outfile->Write();
  m_outHistFile->Write();
  std::cout << "Output muon-pair trees have been written to: " << output_file_path << std::endl;
  std::cout << "Output histograms have been written to: " << output_hist_file_path << std::endl;

  Finalize();
}

#endif
