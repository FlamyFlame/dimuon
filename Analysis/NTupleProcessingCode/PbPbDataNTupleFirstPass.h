#ifndef __PbPbDataNTupleFirstPass_h__
#define __PbPbDataNTupleFirstPass_h__

#include "../MuonObjectsParamsAndHelpers/MuonPairData.h"
#include "../MuonObjectsParamsAndHelpers/muon_pair_enums_data.h"
#include "MuonNTupleFirstPassBaseClass.c"

class PbPbDataNTupleFirstPass : public MuonNTupleFirstPassBaseClass{
    // ONLY ALLOW MODE 3 and 4 NOW!!!
    // Read through the N-tuple, apply appropriate cuts
    // Then fill in histograms and/or output trees
    // mode = 1: raw histograms
    // mode = 2: kill resonances & make no-binning histograms
    // mode = 3: kill resonances & output single-muon information into a TTree
    // mode = 4: kill resonances & output muon-pair information into a TTree
    // mode = 5: kill resonances & make centrality-binned histograms
    // mode = 6: kill resonances & make pT-binned histograms


private:
    // --------------------- general settings ---------------------------

    ParamsSet pms;

    // const ParamsSet pms;
    // int scaleFactorCtrs[ParamsSet::nCtrIntvls] = {1,1,3,5,5};
    // int scaleFactorPts[ParamsSet::nPtBins] = {1,1,1,1,1};

      // --------------------- input files & trees & data for setting branches---------------------------

    TChain          *fChain;   //!pointer to the analyzed TTree or TChain
    std::shared_ptr<MuonPairData> mpair;
    MuonPairData* mpair_raw_ptr = nullptr;

    // Declaration of leaf types
    UInt_t          RunNumber;
    Float_t         ZdcEtA;
    Float_t         ZdcEtC;
    Float_t         FCal_Et;
    Int_t           centrality;
    std::vector<float>   *muon_deltaP_overP           =nullptr;

    std::vector<int>     *muon_pair_muon1_index       =nullptr;
    std::vector<float>   *muon_pair_muon1_pt          =nullptr;
    std::vector<float>   *muon_pair_muon1_eta         =nullptr;
    std::vector<float>   *muon_pair_muon1_phi         =nullptr;
    std::vector<int>     *muon_pair_muon1_quality     =nullptr;
    std::vector<float>   *muon_pair_muon1_d0          =nullptr;
    std::vector<float>   *muon_pair_muon1_z0          =nullptr;

    std::vector<int>     *muon_pair_muon2_index       =nullptr;
    std::vector<float>   *muon_pair_muon2_pt          =nullptr;
    std::vector<float>   *muon_pair_muon2_eta         =nullptr;
    std::vector<float>   *muon_pair_muon2_phi         =nullptr;
    std::vector<int>     *muon_pair_muon2_quality     =nullptr;
    std::vector<float>   *muon_pair_muon2_d0          =nullptr;
    std::vector<float>   *muon_pair_muon2_z0          =nullptr;

    Bool_t               b_HLT_mu4_mu4noL1;
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

    // --------------------- output file, histograms & trees ---------------------------

    TTree* muonOutTree;
    TTree* muonOutTreeBinned[ParamsSet::nCtrIntvls];
    TTree* muonPairOutTree[ParamsSet::nSigns];
    TTree* muonPairOutTreeBinned[ParamsSet::nCtrIntvls][ParamsSet::nSigns];


public :
    int mode = 4;
    bool isRun3;
    PbPbDataNTupleFirstPass(){
        numCuts = numCuts_data;
        cutLabels = cutLabels_data;
        std::cout << "PbPb Data Ntuple processing script:" << std::endl;
        std::cout << "The following public variable(s) should be checked:" << std::endl;
        std::cout << "mode: integer that sets the mode (default 4)" << std::endl;
        std::cout << "isRun3: if true, use run3 data; false: use run2 data" << std::endl;
        std::cout << "if isRun3, output files will be written to /usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_2023" << std::endl;
        std::cout << "else,      output files will be written to /usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_run2" << std::endl;
    }
    ~PbPbDataNTupleFirstPass(){}
    void Run() override;
};

#endif


