#ifndef __PPDataNTupleFirstPass_h__
#define __PPDataNTupleFirstPass_h__

#include "MuonPairData.h"
#include "MuonNTupleFirstPassBaseClass.c"
#include "muon_pair_enums_data.h"

class PPDataNTupleFirstPass : public MuonNTupleFirstPassBaseClass{
    // Read through the N-tuple, apply appropriate cuts
    // Then fill in histograms and/or output trees
    // mode = 3: kill resonances & output single-muon information into a TTree
    // mode = 4: kill resonances & output muon-pair information into a TTree
    // NOW ONLY HAVE MODE = 3, 4


private:
// --------------------- general settings ---------------------------

    int mode = 4;
    ParamsSet pms;
    
// --------------------- input files & trees & data for setting branches---------------------------
  
    TChain          *fChain;   //!pointer to the analyzed TTree or TChain
    std::shared_ptr<MuonPairData> mpair;
    MuonPairData* mpair_raw_ptr = nullptr;

    // Declaration of leaf types
    UInt_t          RunNumber;
    std::vector<float>   *muon_deltaP_overP           =nullptr;
  
    std::vector<int>     *muon_pair_muon1_index       =nullptr;
    std::vector<float>   *muon_pair_muon1_pt          =nullptr;
    std::vector<float>   *muon_pair_muon1_eta         =nullptr;
    std::vector<float>   *muon_pair_muon1_phi         =nullptr;
    std::vector<float>   *muon_pair_muon1_trk_pt          =nullptr;
    std::vector<float>   *muon_pair_muon1_trk_eta         =nullptr;
    std::vector<float>   *muon_pair_muon1_trk_phi         =nullptr;
    std::vector<int>     *muon_pair_muon1_quality     =nullptr;
    std::vector<float>   *muon_pair_muon1_d0          =nullptr;
    std::vector<float>   *muon_pair_muon1_z0          =nullptr;
  
    std::vector<int>     *muon_pair_muon2_index       =nullptr;
    std::vector<float>   *muon_pair_muon2_pt          =nullptr;
    std::vector<float>   *muon_pair_muon2_eta         =nullptr;
    std::vector<float>   *muon_pair_muon2_phi         =nullptr;
    std::vector<float>   *muon_pair_muon2_trk_pt          =nullptr;
    std::vector<float>   *muon_pair_muon2_trk_eta         =nullptr;
    std::vector<float>   *muon_pair_muon2_trk_phi         =nullptr;
    std::vector<int>     *muon_pair_muon2_quality     =nullptr;
    std::vector<float>   *muon_pair_muon2_d0          =nullptr;
    std::vector<float>   *muon_pair_muon2_z0          =nullptr;
  
    Bool_t               b_HLT_2mu4;
    Bool_t               b_HLT_mu4_mu4noL1 = false; // by default set to false: trigger non-existent for run2 data
    std::vector<bool>    *dimuon_b_HLT_2mu4    =nullptr;
    std::vector<bool>    *dimuon_b_HLT_mu4_mu4noL1    =nullptr;
  
// --------------------- temporary muon and muonpair objects ---------------------------
  
    Muon* tempmuon = nullptr;
    // std::vector<std::shared_ptr<MuonPairData>> muon_pair_list_cur_event;
    std::vector<std::shared_ptr<MuonPairData>> muon_pair_list_cur_event_pre_resonance_cut;
  
// --------------------- class methods ---------------------------
  
    void InitInput() override;
    void InitOutput() override;
    void ProcessData() override;
    bool PassCuts(const std::shared_ptr<MuonPair>& mpair) override;
    
    void FillSingleMuonTree() override;
    void FillMuonPairTree();

    void FillMuonPair(int pair_ind, std::shared_ptr<MuonPairData>& mpair);
    bool TrigMatch(int pair_ind);
  
  
// --------------------- output file, histograms & trees ---------------------------
  
    TFile *m_outfile = nullptr;
    
    TTree* muonOutTree;
    TTree* muonPairOutTree[ParamsSet::nSigns];

public :
    bool isTight;
    bool isRun3;
    int run3_file_batch;
    bool run3_use_mu4_mu4_noL1;

    PPDataNTupleFirstPass();
    ~PPDataNTupleFirstPass(){}
    void Run() override;
};


PPDataNTupleFirstPass::PPDataNTupleFirstPass(){
    numCuts = numCuts_data;
    cutLabels = cutLabels_data;

    isTight = false;
    isRun3 = true;
    run3_file_batch = 0;
    run3_use_mu4_mu4_noL1 = true;
}

#endif


