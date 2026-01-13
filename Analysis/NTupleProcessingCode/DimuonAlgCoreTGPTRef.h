// DimuonAlgCoreT.h
#pragma once
#include <memory>

template <class PairT, class MuonT>
class DimuonAlgCoreT {
public:
    using pair_t = PairT;
    using muon_t = MuonT;

protected:
    std::shared_ptr<pair_t> mpair;
    pair_t* mpair_raw_ptr{nullptr};

    muon_t* muon_raw_ptr{nullptr};

    TTree* pairTree{nullptr};
    TTree* muonTree{nullptr};

protected:
    // hooks called by your Run()/event-loop implementation
    virtual void Initialize_Impl() {}
    virtual void InitTempVariables_Impl() {}
    virtual void InitOutput_Impl() {}
    virtual void HistAdjust_Impl() {}
    virtual void Finalize_Impl() {}

    virtual bool PassCuts_Impl(pair_t const& p) { (void)p; return true; }
    virtual void FillPair_Impl() {}
    virtual void FillSingleMuon_Impl(muon_t const& m) { (void)m; }

    void SetupPairBranch() {
        mpair = std::make_shared<pair_t>();
        mpair_raw_ptr = mpair.get();
        if (pairTree) pairTree->Branch("pair", &mpair_raw_ptr);
    }

    void SetupMuonBranch() {
        if (muonTree) muonTree->Branch("muon", &muon_raw_ptr);
    }

    void FillMuonFromPair(bool first) {
        muon_raw_ptr = first ? &mpair->m1 : &mpair->m2;
        FillSingleMuon_Impl(*muon_raw_ptr);
        if (muonTree) muonTree->Fill();
    }
};
