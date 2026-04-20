#include "PbPbExtras.h"
#include "Riostream.h"
#include "TLorentzVector.h"
#include "time.h"

template <class Derived>
void PbPbExtras<Derived>::PerformTChainFill(){

  self().fChainRef() = new TChain("HeavyIonD3PD","HeavyIonD3PD");
  self().fChainRef()->SetMakeClass(1);

  std::string file_path;
  if (!self().isRun3){
    std::string data_trigeffcy_dir = "user.soumya.TrigRates.physics_HP.PbPb20" + std::to_string(self().run_year) + ".15April2022._MYSTREAM/";
    file_path = self().data_dir + data_trigeffcy_dir + "data_pbpb" + std::to_string(self().run_year) + "_trigeffcy_part" + std::to_string(self().file_batch) + ".root";
  } else{
    file_path = self().data_dir + "data_pbpb" + std::to_string(self().run_year) + "_part" + std::to_string(self().file_batch) + ".root";
  }

  // add file to TChain if file path exists
  if (!gSystem->AccessPathName(file_path.c_str())) {
      self().fChainRef()->Add(file_path.c_str());
  } else {
      std::cerr << "File does not exist: " << file_path << std::endl;
      throw std::exception();
  }
}

template <class Derived>
void PbPbExtras<Derived>::InitInputExtra(){
  
  // Enable branches BEFORE SetBranchAddress: with SetMakeClass(1) + prior SetBranchStatus("*",0),
  // calling SetBranchAddress on a disabled branch returns code 3 and doesn't connect the buffer.
  self().fChainRef()->SetBranchStatus("FCal_Et"   , 1);
  self().fChainRef()->SetBranchStatus("FCal_Et_P" , 1);
  self().fChainRef()->SetBranchStatus("FCal_Et_N" , 1);
  self().fChainRef()->SetBranchStatus("centrality", 1);
  self().fChainRef()->SetBranchAddress("FCal_Et"   , &FCal_Et);
  self().fChainRef()->SetBranchAddress("FCal_Et_P" , &FCal_Et_P);  // side A
  self().fChainRef()->SetBranchAddress("FCal_Et_N" , &FCal_Et_N);  // side C
  self().fChainRef()->SetBranchAddress("centrality", &centrality);

  if (self().isRun3) {
    self().fChainRef()->SetBranchStatus("zdc_ZdcEnergy",              1);
    self().fChainRef()->SetBranchStatus("zdc_ZdcTime",                1);
    self().fChainRef()->SetBranchStatus("zdc_ZdcModulePreSampleAmp",  1);
    self().fChainRef()->SetBranchStatus("trk_numqual",                1);
    self().fChainRef()->SetBranchAddress("zdc_ZdcEnergy",             zdc_ZdcEnergy);
    self().fChainRef()->SetBranchAddress("zdc_ZdcTime",               zdc_ZdcTime);
    self().fChainRef()->SetBranchAddress("zdc_ZdcModulePreSampleAmp", zdc_ZdcModulePreSampleAmp);
    self().fChainRef()->SetBranchAddress("trk_numqual",               &trk_numqual);
  }
}

template <class Derived>
void PbPbExtras<Derived>::InitOutputSettingsExtra(){
  if (turn_on_ctr_binned_tree_writing){
    if (self().output_single_muon_tree){
      // in addition, write centrality-binned trees: need for scrambling
      for (int jctr = 0; jctr < ParamsSet::nCtrIntvls; jctr++){
        muonOutTreeBinned[jctr] = new TTree(Form("muon_tree_ctr%d",jctr+1),Form("all muons, centrality bin %d",jctr+1));
        muonOutTreeBinned[jctr]->Branch("MuonObj",&self().muon_raw_ptr);
      }
    }
    else{
      for (unsigned int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
        for (unsigned int jctr = 0; jctr < ParamsSet::nCtrIntvls; jctr++){
          muonPairOutTreeBinned[jctr][ksign] = new TTree(Form("muon_pair_tree_ctr%u_sign%u",jctr+1,ksign+1),Form("all muon pairs, ctr%u, sign%u",jctr+1,ksign+1));
          muonPairOutTreeBinned[jctr][ksign]->Branch("MuonPairObj",&self().mpair_raw_ptr);
        }        
      }      
    }
  }
}


template <class Derived>
void PbPbExtras<Derived>::FillMuonPairExtra(int pair_ind){
  if (self().debug_mode) std::cout << "DEBUG: Calling PbPbExtras::FillMuonPairPbPb" << std::endl;

  if (!self().mpairRef()) {
    throw std::runtime_error("PbPbExtras::FillMuonPairExtra: mpair is nullptr");
  }

  self().mpairRef()->m1.ev_centrality = centrality;
  self().mpairRef()->m2.ev_centrality = centrality;
  self().mpairRef()->FCal_Et   = FCal_Et   * 1e-6f;  // MeV → TeV
  self().mpairRef()->FCal_Et_A = FCal_Et_P * 1e-6f;  // MeV → TeV
  self().mpairRef()->FCal_Et_C = FCal_Et_N * 1e-6f;  // MeV → TeV
  self().mpairRef()->m1.ev_FCal_Et = FCal_Et * 1e-6f;
  self().mpairRef()->m2.ev_FCal_Et = FCal_Et * 1e-6f;

  if (self().isRun3){
    self().mpairRef()->year = self().run_year;

    self().mpairRef()->ZDC_E_tot = zdc_ZdcEnergy[0] + zdc_ZdcEnergy[1];  // [0]=A, [1]=C
    self().mpairRef()->ZDC_t_A   = zdc_ZdcTime[0];   // [0] = A-side
    self().mpairRef()->ZDC_t_C   = zdc_ZdcTime[1];   // [1] = C-side
    float preamp_A = 0.f, preamp_C = 0.f;
    for (int i = 0; i < 4; ++i) preamp_A += zdc_ZdcModulePreSampleAmp[0][i];  // [0]=A
    for (int i = 0; i < 4; ++i) preamp_C += zdc_ZdcModulePreSampleAmp[1][i];  // [1]=C
    self().mpairRef()->ZDC_preamp_A = preamp_A;
    self().mpairRef()->ZDC_preamp_C = preamp_C;

    if (trk_numqual && trk_numqual->size() >= 8) {
      self().mpairRef()->ntrk_HIloose         = (*trk_numqual)[2];
      self().mpairRef()->ntrk_HItight         = (*trk_numqual)[3];
      self().mpairRef()->ntrk_HIloose_noPtCut = (*trk_numqual)[6];
      self().mpairRef()->ntrk_HItight_noPtCut = (*trk_numqual)[7];
    }
  } else {
    self().mpairRef()->year = self().run_year;
    self().mpairRef()->UpdateCentrality();
  }
}

template <class Derived>
void PbPbExtras<Derived>::FillSingleMuonTreeExtra(){
  if (turn_on_ctr_binned_tree_writing){
    for (int jctr = 0; jctr < self().pms.nCtrIntvls; jctr++){
      if (self().muon_raw_ptr->ev_centrality >= self().pms.CtrStep * jctr && self().muon_raw_ptr->ev_centrality < self().pms.CtrStep *(jctr+1)){
        muonOutTreeBinned[jctr]->Fill();
      }
    }    
  }
}

template <class Derived>
void PbPbExtras<Derived>::FillMuonPairTreeExtra(){
  if (turn_on_ctr_binned_tree_writing){
    int nsign = (self().mpairRef()->same_sign)? 0:1;
    for (unsigned int jctr = 0; jctr < self().pms.nCtrIntvls; jctr++){
      if (self().mpairRef()->avg_centrality >= self().pms.CtrStep * jctr && self().mpairRef()->avg_centrality < self().pms.CtrStep *(jctr+1)){
        muonPairOutTreeBinned[jctr][nsign]->Fill();
      }
    }      
  }
}


template <class Derived>
void PbPbExtras<Derived>::InitParamsExtra(){
  self().isPbPb = true;
  self().cutLabels.assign(cutLabelsPbPb.begin(), cutLabelsPbPb.end()); // cuts (filters) on muon pairs, read from muon_pair_enums_data.h
  self().numCuts = static_cast<int>(CutsPbPb::nCuts_PbPb_data);

  int run_year_short = self().run_year % 2000;
  bool is_run3_local = (run_year_short > 20);

  std::map<int, int> run_year_to_file_batch_max_map = {
    {23, 6}, {24, 2}, {15, 7}, {18, 7}
  };

  // check for run year
  if (is_run3_local){
    if (run_year_short != 23 && run_year_short != 24){
      std::cerr<<"Error:: If isRun3 is true, run_year must be set to (20)23 or (20)24"<<std::endl;
      throw std::exception();
    }
  }else{
    if (run_year_short != 15 && run_year_short != 18){
      std::cerr<<"Error:: If isRun3 be false, run_year must be set to (20)15 or (20)18"<<std::endl;
      throw std::exception();      
    }
  }

  // check for file batch
  if (self().file_batch <= 0 || self().file_batch > run_year_to_file_batch_max_map[run_year_short]){
    std::cerr<<"Error:: run3 file_batch is invalid! Must be in range 1-6 for 2023 data / 1-2 for 2024 data / 1-7 for 2015/2018 data"<<std::endl;
    throw std::exception();
  }
}
