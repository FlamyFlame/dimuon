// PowhegNTupleFirstPassT.h
#pragma once

#include "DimuonAlgCoreT.h"
#include "PowhegExtrasInterface.h"

// your existing includes:
#include "../MuonObjectsParamsAndHelpers/TruthQQPair.h"
#include "../MuonObjectsParamsAndHelpers/muon_pair_enums_MC.h"
#include "../MuonObjectsParamsAndHelpers/struct_particle.h"

template <class PairT, class MuonT, class ExtraT>
class PowhegNTupleFirstPassT
    : public DimuonAlgCoreT<PairT, MuonT>
{
public:
    using core_t = DimuonAlgCoreT<PairT, MuonT>;
    using pair_t = typename core_t::pair_t;
    using muon_t = typename core_t::muon_t;

protected:
    ExtraT extra;   // owns truth/fullsim/overlay-specific state

    // ---- your current Powheg members (keep them here) ----
    int file_batch{};
    std::string mc_mode;
    bool is_fullsim{false};

    double crossx_cut{0.0};
    double filter_effcy{0.0};

    // input chain, branches, etc (from your current PowhegNTupleFirstPass)
    TChain* fChain{nullptr};
    Long64_t nentries{0};

    // ... all Powheg-specific branch pointers, vectors, etc ...

public:
    // keep your Run() pattern if you like
    void Run() {
        Initialize();
        ProcessData();
        Finalize();
    }

protected:
    // ---- lifecycle wrapper calls ----
    void Initialize() {
        core_t::Initialize_Impl();
        extra.Initialize(*this);
    }

    void InitTempVariables() {
        core_t::InitTempVariables_Impl();
        extra.InitTempVariables(*this);
    }

    void InitOutput() {
        core_t::InitOutput_Impl();
        extra.InitOutput(*this);
    }

    void HistAdjust() {
        core_t::HistAdjust_Impl();
        extra.HistAdjust(*this);
    }

    void Finalize() {
        extra.Finalize(*this);
        core_t::Finalize_Impl();
    }

    // ---- key: one PassCuts that chains base + extra ----
    bool PassCuts(pair_t const& p) {
        if (!core_t::PassCuts_Impl(p)) return false;
        if (!extra.PassCuts(*this, p)) return false;
        return true;
    }

    // ---- your main loop; call extra hooks ----
    void ProcessData() {
        nentries = fChain->GetEntries();

        for (Long64_t jentry = 0; jentry < nentries; ++jentry) {
            if (jentry % 10000 == 0)
                std::cout << "Processing " << jentry << " / " << nentries << "\n";

            int num_bytes = fChain->GetEntry(jentry);
            if (num_bytes == 0) throw std::runtime_error("GetEntry returned 0 bytes");

            extra.OnBeginEvent(*this);

            int NPairs = muon_pair_muon1_pt->size();
            for (int pair_ind = 0; pair_ind < NPairs; ++pair_ind) {
                // allocate the *concrete* pair type
                this->mpair = std::make_shared<pair_t>();
                this->mpair_raw_ptr = this->mpair.get();

                FillMuonPair(pair_ind, *this->mpair);    // common fill
                this->mpair->Update();                    // base + hook in pair object

                // cuts
                if (!PassCuts(*this->mpair)) continue;

                // truth/fullsim/overlay per-pair analysis
                extra.OnPair(*this, *this->mpair);

                // output pair tree
                // (assuming pairTree exists)
                if (this->pairTree) this->pairTree->Fill();

                // optional single muon output
                if (output_single_muon_tree && this->muonTree) {
                    this->FillMuonFromPair(true);
                    extra.OnMuon(*this, this->mpair->m1);
                    this->FillMuonFromPair(false);
                    extra.OnMuon(*this, this->mpair->m2);
                }
            }

            extra.OnEndEvent(*this);
        }
    }

    // ---- your existing FillMuonPair logic, now templated on PairT ----
    void FillMuonPair(int pair_ind, pair_t& p) {
        // port your PowhegNTupleFirstPass::FillMuonPair content here
        // use p.m1, p.m2 and p.<powheg extras> fields
    }

public:
    bool output_single_muon_tree{false};
};
