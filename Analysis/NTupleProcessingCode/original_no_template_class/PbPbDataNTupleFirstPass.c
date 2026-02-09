#include "PbPbDataNTupleFirstPass.h"
#include "Riostream.h"
#include "TLorentzVector.h"
#include "time.h"


PbPbDataNTupleFirstPass::PbPbDataNTupleFirstPass(int run_year_input, int file_batch_input)
    : DimuonDataAnalysisBaseClass(run_year_input, file_batch_input){
        isPbPb = true;
        cutLabels.assign(cutLabelsPbPb.begin(), cutLabelsPbPb.end()); // cuts (filters) on muon pairs, read from muon_pair_enums_data.h
        numCuts = static_cast<int>(CutsPbPb::nCuts_PbPb_data);
        PrintInstructions();
}

//initialize the TChain
void PbPbDataNTupleFirstPass::TChainFill(){

  fChain = new TChain("HeavyIonD3PD","HeavyIonD3PD");
  fChain->SetMakeClass(1);

  std::string file_path;
  if (!isRun3){
    std::string data_trigeffcy_dir = "user.soumya.TrigRates.physics_HP.PbPb20" + std::to_string(run_year) + ".15April2022._MYSTREAM/";
    file_path = data_dir + data_trigeffcy_dir + "data_pbpb" + std::to_string(run_year) + "_trigeffcy_part" + std::to_string(file_batch) + ".root";
  } else{
    file_path = data_dir + "data_pbpb" + std::to_string(run_year) + "_part" + std::to_string(file_batch) + ".root";
  }

  // add file to TChain if file path exists
  if (!gSystem->AccessPathName(file_path.c_str())) {
      fChain->Add(file_path.c_str());
  } else {
      std::cerr << "File does not exist: " << file_path << std::endl;
      throw std::exception();
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

  DimuonDataAnalysisBaseClass::InitOutput(); // pair-acceptance hist writing + non-centrality-binned muon(-pair) tree writing

  if (turn_on_ctr_binned_tree_writing){
    if (output_single_muon_tree){
      // in addition, write centrality-binned trees: need for scrambling
      for (int jctr = 0; jctr < ParamsSet::nCtrIntvls; jctr++){
        muonOutTreeBinned[jctr] = new TTree(Form("muon_tree_ctr%d",jctr+1),Form("all muons, centrality bin %d",jctr+1));
        muonOutTreeBinned[jctr]->Branch("MuonObj",&tempmuon);
      }
    }
    else{
      for (unsigned int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
        for (unsigned int jctr = 0; jctr < ParamsSet::nCtrIntvls; jctr++){
          muonPairOutTreeBinned[jctr][ksign] = new TTree(Form("muon_pair_tree_ctr%u_sign%u",jctr+1,ksign+1),Form("all muon pairs, ctr%u, sign%u",jctr+1,ksign+1));
          muonPairOutTreeBinned[jctr][ksign]->Branch("MuonPairObj",&mpair_raw_ptr);
        }        
      }      
    }
  }
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

  if (isRun3){  // run3: no need to update centrality (for now)
    mpair->year = run_year;
  }else{
    mpair->year = run_year;
    mpair->UpdateCentrality();
  }
}


void PbPbDataNTupleFirstPass::FillSingleMuonTree(){
  DimuonAnalysisBaseClass::FillSingleMuonTree();

  if (turn_on_ctr_binned_tree_writing){
    for (int jctr = 0; jctr < pms.nCtrIntvls; jctr++){
      if (tempmuon->ev_centrality >= pms.CtrStep * jctr && tempmuon->ev_centrality < pms.CtrStep *(jctr+1)){
        muonOutTreeBinned[jctr]->Fill();
      }
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


void PbPbDataNTupleFirstPass::InitParams(){

  DimuonDataAnalysisBaseClass::InitParams();
  
  std::map<int, int> run_year_to_file_batch_max_map = {
    {23, 6}, {24, 9}, {15, 7}, {18, 7}
  };

  // check for run year
  if (isRun3){
    if (run_year != 23 && run_year != 24){
      std::cerr<<"Error:: If isRun3 is true, run_year must be set to (20)23 or (20)24"<<std::endl;
      throw std::exception();
    }
  }else{
    if (run_year != 15 && run_year != 18){
      std::cerr<<"Error:: If isRun3 be false, run_year must be set to (20)15 or (20)18"<<std::endl;
      throw std::exception();      
    }
  }

  // check for file batch
  if (file_batch <= 0 || file_batch > run_year_to_file_batch_max_map[run_year]){
    std::cerr<<"Error:: run3 file_batch is invalid! Must be in range 1-6 for 2023 data / 1-9 for 2024 data / 1-7 for 2015/2018 data"<<std::endl;
    throw std::exception();
  }
}
