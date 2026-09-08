#include "PPExtras.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include "TSystem.h"
#include "../Utilities/AllVertexIPSelection.h"

template <class Derived>
void PPExtras<Derived>::InitParamsExtra(){
    self().isPbPb = false;
    self().cutLabels.assign(cutLabelsPP.begin(), cutLabelsPP.end());
    self().numCuts = static_cast<int>(CutsPP::nCuts_pp_data);

    if (self().run_year != 24) {
        std::cerr << "Error:: pp run_year must be 24" << std::endl;
        throw std::exception();
    }

    const int file_batch_max = 12;
    if (self().file_batch <= 0 || self().file_batch > file_batch_max) {
        std::cerr << "Error:: pp file_batch invalid! Must be 1-"
                  << file_batch_max
                  << " for 20" << self().run_year << " data" << std::endl;
        throw std::exception();
    }
}

//initialize the TChain
template <class Derived>
void PPExtras<Derived>::PerformTChainFill(){

    self().fChainRef() = new TChain("HeavyIonD3PD","HeavyIonD3PD");
    self().fChainRef()->SetMakeClass(1);

    std::string file_path = self().data_dir + "data_pp" + std::to_string(self().run_year)
                          + "_part" + std::to_string(self().file_batch) + ".root";

    if (!gSystem->AccessPathName(file_path.c_str())) {
        self().fChainRef()->Add(file_path.c_str());
    } else {
        std::cerr << "File does not exist: " << file_path << std::endl;
        throw std::exception();
    }
}


// =========================================================================================
// ALL-VERTEX impact-parameter selection (pp DATA).
// The algorithm itself lives in Utilities/AllVertexIPSelection.h and is SHARED with the
// pp-conditions fullsim MC, so the data selection and the eps_reco it is corrected by can
// never drift apart. Read that header for the physics, the z0_v identity and the branch-binding
// order this file depends on.
//
// NOTE: these branches are enabled and bound only on the DIMUON path. On the (currently dead)
// isMinBias / single-muon path InitInputBranchesDimuonAnalysisExtra never runs, so vtx_z stays
// null and AllVertexIP::CheckVertexArrays THROWS -- deliberately, rather than silently
// rejecting every pair against an unread vertex list.
// =========================================================================================

template <class Derived>
void PPExtras<Derived>::InitInputBranchesDimuonAnalysisExtra(){
    // Runs AFTER InitInputBranchesDimuonAnalysis_DataCore, i.e. after its
    // SetBranchStatus("*",0) and after the chain's first tree has been loaded (GetBranch()
    // loads it). ENABLE BEFORE BINDING -- that order is what makes the binding work here, and
    // it is not cosmetic: measured on this very chain, SetBranchAddress on an STL-collection
    // branch of a SetMakeClass(1) TChain leaves the pointer NULL when the branch is still
    // disabled AND the first tree is already loaded, which is exactly the situation at this
    // point. Enabled first it allocates normally. Full truth table in
    // Utilities/AllVertexIPSelection.h (BRANCH BINDING).
    auto* ch = self().fChainRef();

    if (!ch->GetBranch("vtx_z") || !ch->GetBranch("vtx_ntrk"))
        throw std::runtime_error("PPExtras: the pp NTUP has no vtx_z/vtx_ntrk branches, so the "
                                 "all-vertex impact-parameter selection cannot be applied. "
                                 "Re-skim with m_store_Vtx enabled.");

    ch->SetBranchStatus ("vtx_z"   , 1);
    ch->SetBranchStatus ("vtx_ntrk", 1);
    ch->SetBranchAddress("vtx_z"   , &vtx_z);
    ch->SetBranchAddress("vtx_ntrk", &vtx_ntrk);

    std::cout << "INFO: pp all-vertex impact-parameter selection ENABLED "
                 "(pairs may come from secondary vertices)." << std::endl;
}


template <class Derived>
bool PPExtras<Derived>::PassD0Z0Extra(){
    auto& pair = *self().mpairRef();
    const auto& pms = self().pmsRef();

    pair.matched_vtx_ind = AllVertexIP::BestCommonVertex(
        pair.m1.d0, pair.m1.z0, pair.m1.eta,
        pair.m2.d0, pair.m2.z0, pair.m2.eta,
        vtx_z, vtx_ntrk, pms.d0cut, pms.z0cut,
        &pair.pass_primary_vtx);

    return (pair.matched_vtx_ind >= 0);
}
