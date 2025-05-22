#ifndef __PbPbDataNTupleFirstPass_C__
#define __PbPbDataNTupleFirstPass_C__

#include "PbPbDataNTupleFirstPass.h"
#include "Riostream.h"
#include "TLorentzVector.h"
#include "time.h"


PbPbDataNTupleFirstPass::PbPbDataNTupleFirstPass(){
    // cuts (filters) on muon pairs, read from muon_pair_enums_data.h
    cutLabels.assign(cutLabelsPbPb.begin(), cutLabelsPbPb.end());
    numCuts = static_cast<int>(CutsPbPb::nCuts_PbPb_data);

    isPbPb = true;
    PrintInstructions();

    std::cout << "PbPb Data Ntuple processing script:" << std::endl;
    std::cout << "The following public variable(s) should be checked:" << std::endl;
    std::cout << "--> run_year: [INT] (20)23 or (20)24; MUST be set if isRun3 is true" << std::endl;
    std::cout << "--> file_batch: [INT] file batch: 1-6 for 2023 data, 1-9 for 2024 data" << std::endl;
    std::cout << "--> isRun3: [BOOL] if true (DEFAULT), use run3 data" << std::endl;
    std::cout << "                   if false: use run2 data" << std::endl;
    std::cout << std::endl;

    std::cout << "The following public variable(s) should be checked:" << std::endl;
    std::cout << "--> trigger_mode: INT, value = 1,2,3" << std::endl;
    std::cout << "                       value = 1 (DEFAULT): require single-muon trigger: tag which muon files trigger & if pair passes mu4_mu4_noL1 & 2mu4" << std::endl;
    std::cout << "                       value = 2: require mu4_mu4_noL1" << std::endl;
    std::cout << "                       value = 3: require 2mu4" << std::endl;
    std::cout << "--> apply_resonance_cuts: [BOOL] if true (DEFAULT), apply (old) resonance cuts" << std::endl;
    std::cout << "--> output_single_muon_tree: [BOOL] if true, output single-muon tree" << std::endl;
    std::cout << "                                    if false (DEFAULT), output muon-pair tree" << std::endl;
    std::cout << "--> turn_on_ctr_binned_tree_writing: [BOOL] if true, output centrality-binned muon-pair trees (DEFAULT FALSE)" << std::endl;
    std::cout << std::endl;

    std::cout << "if run_year == 23, output files will be written to /usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_2023" << std::endl;
    std::cout << "if run_year == 24, output files will be written to /usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_2024" << std::endl;
    std::cout << "else,      output files will be written to /usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_run2" << std::endl;
}

//initialize the TChain
void PbPbDataNTupleFirstPass::TChainFill(){

  fChain = new TChain("HeavyIonD3PD","HeavyIonD3PD");
  fChain->SetMakeClass(1);

  if (!isRun3){
    fChain->Add("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_run2/All_Data2018_12March2022.root");
    fChain->Add("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_run2/All_Data2015_12March2022.root");     
  } else{
    // add file to TChain if file path exists
    std::string file_path = in_out_file_dir + "data_pbpb" + std::to_string(run_year) + "_part" + std::to_string(file_batch) + ".root";
    if (!gSystem->AccessPathName(file_path.c_str())) {
        fChain->Add(file_path.c_str());
    } else {
        std::cerr << "File does not exist: " << file_path << std::endl;
        throw std::exception();
    }
  }
}

void PbPbDataNTupleFirstPass::InitInput(){

  DimuonDataAnalysisBaseClass::InitInput();
  
  fChain->SetBranchAddress("FCal_Et"                     , &FCal_Et);
  fChain->SetBranchAddress("centrality"                  , &centrality);

  fChain->SetBranchStatus("FCal_Et"                         ,1);
  fChain->SetBranchStatus("centrality"                      ,1);
}

void PbPbDataNTupleFirstPass::InitOutput(){

  std::string file_name_base;
  std::string outfile_path;
  std::string resonance_cuts_suffix = (apply_resonance_cuts)? "" : "_no_resn_cuts";

  TrigModeToSuffixMap();

  std::string outfile_ending = "_pbpb_20" + std::to_string(run_year) + "_part" + std::to_string(file_batch) + trig_suffix + resonance_cuts_suffix +".root";
  
  if (output_single_muon_tree){
    file_name_base = "single_muon_trees";

    outfile_path = in_out_file_dir + file_name_base + outfile_ending;
    m_outfile=new TFile(outfile_path.c_str(),"recreate");

    muonOutTree = new TTree("muon_tree","all single muons");
    muonOutTree->Branch("MuonObj",&tempmuon);

    // by default write centrality-binned trees: need for scrambling
    for (int jctr = 0; jctr < ParamsSet::nCtrIntvls; jctr++){
      muonOutTreeBinned[jctr] = new TTree(Form("muon_tree_ctr%d",jctr+1),Form("all muons, centrality bin %d",jctr+1));
      muonOutTreeBinned[jctr]->Branch("MuonObj",&tempmuon);
    }
  }
  else{
    file_name_base = "muon_pairs";

    outfile_path = in_out_file_dir + file_name_base + outfile_ending;
    m_outfile=new TFile(outfile_path.c_str(),"recreate");


    for (unsigned int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
      muonPairOutTree[ksign] = new TTree(Form("muon_pair_tree_sign%u",ksign+1),Form("all muon pairs, sign%u",ksign+1));
      muonPairOutTree[ksign]->Branch("MuonPairObj",&mpair_raw_ptr);

      if (turn_on_ctr_binned_tree_writing){
        for (unsigned int jctr = 0; jctr < ParamsSet::nCtrIntvls; jctr++){
          muonPairOutTreeBinned[jctr][ksign] = new TTree(Form("muon_pair_tree_ctr%u_sign%u",jctr+1,ksign+1),Form("all muon pairs, ctr%u, sign%u",jctr+1,ksign+1));
          muonPairOutTreeBinned[jctr][ksign]->Branch("MuonPairObj",&mpair_raw_ptr);
        }        
      }
    }
  }
  
  // ---------------------------------------------------------------------------------------------------------------------------

  std::string outhistfile_name_base = "hists_cut_acceptance";
  std::string outhistfile_path = in_out_file_dir + outhistfile_name_base + outfile_ending;

  m_outHistFile=new TFile(outhistfile_path.c_str(),"recreate");
  DimuonAnalysisBaseClass::InitOutput();
}


void PbPbDataNTupleFirstPass::FillMuonPair(int pair_ind, std::shared_ptr<MuonPairData> const& mpair){
  if (debug_mode) std::cout << "DEBUG: Calling PbPbDataNTupleFirstPass::FillMuonPair" << std::endl;
  DimuonDataAnalysisBaseClass::FillMuonPair(pair_ind, mpair);
  FillMuonPairPbPb(pair_ind, std::static_pointer_cast<MuonPairPbPb>(mpair)); 
}


void PbPbDataNTupleFirstPass::FillMuonPairPbPb(int pair_ind, std::shared_ptr<MuonPairPbPb> const& mpair){
  if (debug_mode) std::cout << "DEBUG: Calling PbPbDataNTupleFirstPass::FillMuonPairPbPb" << std::endl;

  mpair->m1.ev_centrality = centrality;
  mpair->m2.ev_centrality = centrality;
  mpair->FCal_Et = FCal_Et;
  mpair->m1.ev_FCal_Et = FCal_Et;
  mpair->m2.ev_FCal_Et = FCal_Et;

  if (isRun3){
    mpair->year = run_year;
  }else{ // run3: no need to update centrality (for now)
    run_year = 2015;
    mpair->year = run_year;
    mpair->UpdateCentrality();
  }
}


void PbPbDataNTupleFirstPass::FillSingleMuonTree(){
  DimuonDataAnalysisBaseClass::FillSingleMuonTree();
  for (int jctr = 0; jctr < pms.nCtrIntvls; jctr++){
    if (tempmuon->ev_centrality >= pms.CtrStep * jctr && tempmuon->ev_centrality < pms.CtrStep *(jctr+1)){
      muonOutTreeBinned[jctr]->Fill();
    }
  }
}

void PbPbDataNTupleFirstPass::FillMuonPairTree(){

  DimuonDataAnalysisBaseClass::FillMuonPairTree();

  if (turn_on_ctr_binned_tree_writing){
    int nsign = (mpair->same_sign)? 0:1;
    for (unsigned int jctr = 0; jctr < pms.nCtrIntvls; jctr++){
      if (mpair->avg_centrality >= pms.CtrStep * jctr && mpair->avg_centrality < pms.CtrStep *(jctr+1)){
        muonPairOutTreeBinned[jctr][nsign]->Fill();
      }
    }      
  }
}


void PbPbDataNTupleFirstPass::ParamCheck(){
  run_year = run_year % 2000; // 2023 --> 23, 2024 --> 24

  if (isRun3){ // run3: run_year & file_batch must be specified
    if (run_year != 23 && run_year != 24){
      std::cerr<<"Error:: Either isRun3 be false, or run_year must be set to (20)23 or (20)24"<<std::endl;
      throw std::exception();
    }
    if (file_batch <= 0){
      std::cerr<<"Error:: run3 file_batch is invalid! Must be in range 1-6 for 2023 data, and 1-9 for 2024 data"<<std::endl;
      throw std::exception();
    }
  }

  // input, output file directory setting
  if (isRun3){
    in_out_file_dir = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_20" + std::to_string(run_year) + "/";
  }else{
    in_out_file_dir = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_run2/";
  }
}

void PbPbDataNTupleFirstPass::Run(){

  clock_t start, end;
  double cpu_time_used;
  start = clock();

  ParamCheck();
  InitInput();
  InitOutput();
  ProcessData();
  HistAdjust();

  m_outfile->Write();
  m_outHistFile->Write();

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  std::cout << "CPU time used (in seconds) is " << cpu_time_used << std::endl;

}

#endif
