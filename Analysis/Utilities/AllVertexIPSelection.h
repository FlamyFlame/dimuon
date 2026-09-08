#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

// =========================================================================================
// ALL-VERTEX impact-parameter selection -- THE single implementation.
//
// SINGLE SOURCE OF TRUTH. The pp data NTuple stage (NTupleProcessingCode/PPExtras.c) and the
// pp-conditions fullsim MC (PythiaFullSimExtras.c, PowhegFullSimExtras.c) MUST both go through
// this header. If data and MC ever applied two separately-written copies of this algorithm,
// eps_reco would silently correct a selection it was not measured on -- so do not re-implement
// it anywhere; call these functions.
//
// PHYSICS. pp24 has ~4 inelastic collisions per bunch crossing, and the pp luminosity is
// defined over ALL of them. Counting only muon pairs that point at the primary vertex therefore
// puts an under-counted numerator over an all-collision L_int and biases dsigma/dX low. So a
// pair of two good muons is kept when both muons pass the |d0| and |z0 sin(theta)| cuts with
// respect to the SAME reconstructed vertex, primary or secondary.
// Pb+Pb is the opposite case (<#collisions>/crossing ~ 1e-3 AND pile-up rejected at event
// level) and keeps the primary-vertex-only cut. The HIJING overlay has exactly one
// track-bearing vertex in 100% of events, where this is an exact no-op.
// See docs/tracking/pp24_all_vertex_pairs.md.
//
// MECHANISM. The skim (SkimCode/.../HFtrigValidation/src/TrigRates.cxx) stores, per muon,
//     z0_stored = idTrk->z0() + idTrk->vz() - vtx_z[0]                      (ProcessMuons)
// and stores the whole vertex list in the SAME container order              (ProcessVertex),
// so the z0 with respect to any vertex v follows exactly, with no approximation:
//     z0_v = z0_stored + vtx_z[0] - vtx_z[v].
// d0 = idTrk->d0() is referenced to the BEAMLINE, not to a vertex, so it is vertex-independent
// and is tested once. (Vertices in a pp event differ essentially only in z: beam-spot
// transverse size ~10 um against a 2 mm d0 cut.) Keeping d0 as-is also leaves the d0 cut
// byte-identical to the primary-vertex procedure it replaces.
// If the skim's vertex handling ever changes, THIS identity is what has to be re-verified.
//
// BRANCH BINDING (measured, ROOT 6.34.04, on data_pp24_part11.root). These chains run in
// SetMakeClass(1) mode. SetBranchAddress on an STL-collection branch fails to allocate the
// std::vector -- leaving the pointer null -- if and ONLY IF BOTH of these hold: the chain's
// first tree is already loaded AND the branch is currently DISABLED:
//     bind before first load, enabled  -> rc = 5 (kNoCheck), allocated
//     bind before first load, disabled -> rc = 5 (kNoCheck), allocated
//     bind after  first load, enabled  -> rc = 3 (kMakeClass), allocated
//     bind after  first load, disabled -> rc = 3 (kMakeClass), **NULL**
// Note that GetBranch()/GetListOfBranches()/GetEntries() all load the first tree, so "after the
// load" is the normal state by the time any Extra hook runs. The robust rule is therefore
// simply: **SetBranchStatus(name, 1) BEFORE SetBranchAddress(name, ...)** -- which is exactly
// what Utilities/tchain_helpers.h::enable_and_bind does, and why every call site here uses it
// or replicates that order. Ordering against the first tree load does NOT matter.
// =========================================================================================

namespace AllVertexIP {

// The skim dumps PrimaryVertices unfiltered, so the container always ALSO holds exactly one
// dummy beam-spot entry with ntrk == 0 whose z equals the primary's (measured on pp24: always
// the LAST entry, z identical to vtx_z[0] in 100% of events). Without this guard it would act
// as a duplicate primary vertex. Same guard as MuonFullsimExtra::n_vtx.
inline bool IsEligibleVertex(int ntrk) { return ntrk >= 2; }

inline double Z0SinTheta(double z0, double eta) {
    return std::fabs(z0 * std::sin(2.0 * std::atan(std::exp(-eta))));
}

inline void CheckVertexArrays(const std::vector<float>* vtx_z,
                              const std::vector<int>*   vtx_ntrk,
                              const char* caller) {
    if (!vtx_z || !vtx_ntrk)
        throw std::runtime_error(std::string(caller) + ": vtx_z/vtx_ntrk not bound. Enable the "
                                 "branch (SetBranchStatus(name,1)) BEFORE SetBranchAddress -- "
                                 "see the BRANCH BINDING note at the top of this header.");
    if (vtx_z->size() != vtx_ntrk->size())
        throw std::runtime_error(std::string(caller) + ": vtx_z and vtx_ntrk differ in length");
    // An empty vertex list cannot happen: the skim rejects the event unless the PrimaryVertices
    // container holds more than one entry (TrigRates.cxx, m_MaxZvtx cleaning), and it always
    // also contains the dummy beam-spot entry. If it IS empty, the branch is silently unread
    // (disabled, or bound on only some chains of a TChain) -- and every pair would then be
    // rejected with no error at all. Fail loudly instead.
    if (vtx_z->empty())
        throw std::runtime_error(std::string(caller) + ": the vertex list is EMPTY. The skim "
                                 "cannot produce that, so the branch is not actually being "
                                 "read -- check that vtx_z/vtx_ntrk are enabled on EVERY chain.");
}

// SINGLE-MUON form: does this muon point at ANY eligible vertex?
// Used for the per-muon working-point flags, where no partner exists to share a vertex with.
// It is strictly looser than the pair form below, which the pair flags then apply on top.
inline bool PassAnyVertex(double d0, double z0, double eta,
                          const std::vector<float>* vtx_z,
                          const std::vector<int>*   vtx_ntrk,
                          double d0cut, double z0cut) {
    CheckVertexArrays(vtx_z, vtx_ntrk, "AllVertexIP::PassAnyVertex");
    if (std::fabs(d0) >= d0cut) return false;          // vertex-independent

    const double z_primary = static_cast<double>(vtx_z->at(0));
    for (std::size_t iv = 0; iv < vtx_z->size(); ++iv) {
        if (!IsEligibleVertex(vtx_ntrk->at(iv))) continue;
        const double dz = z_primary - static_cast<double>(vtx_z->at(iv));
        if (Z0SinTheta(z0 + dz, eta) < z0cut) return true;
    }
    return false;
}

// PAIR form: index of the eligible vertex with respect to which BOTH muons pass, or -1.
//
// When several vertices qualify, the one returned is the vertex minimising the WORSE of the two
// muons' |z0 sin(theta)| -- the best-matching vertex. Ties break to the lower index, so vertex 0
// wins whenever it is as good as any other and the assignment can only move off the primary
// when a secondary vertex is strictly better. (The ambiguity is small by construction: on
// data_pp24_part11.root, first 300 000 events, eligible vertices only, just 0.60 % of vertex
// pairs sit within 2 mm in z, against an RMS spread of 58.3 mm and a mean |dz| of 67.3 mm.)
//
// `pass_primary_out`, when non-null, receives whether vertex 0 ITSELF qualifies -- i.e. whether
// the pair would also have survived the old primary-vertex-only cut. The two statistics differ
// and must not be confused:
//     "from a secondary vertex"      <=>  returned index != 0
//     "ADDED by the all-vertex rule" <=>  !pass_primary
inline int BestCommonVertex(double d0_1, double z0_1, double eta_1,
                            double d0_2, double z0_2, double eta_2,
                            const std::vector<float>* vtx_z,
                            const std::vector<int>*   vtx_ntrk,
                            double d0cut, double z0cut,
                            bool* pass_primary_out = nullptr) {
    if (pass_primary_out) *pass_primary_out = false;

    CheckVertexArrays(vtx_z, vtx_ntrk, "AllVertexIP::BestCommonVertex");

    // d0 is vertex-independent, so a pair failing it can be rescued by no vertex.
    if (std::fabs(d0_1) >= d0cut || std::fabs(d0_2) >= d0cut) return -1;

    const double z_primary = static_cast<double>(vtx_z->at(0));
    int    best_iv         = -1;
    double best_worst      = 0.0;

    for (std::size_t iv = 0; iv < vtx_z->size(); ++iv) {
        if (!IsEligibleVertex(vtx_ntrk->at(iv))) continue;

        const double dz     = z_primary - static_cast<double>(vtx_z->at(iv));
        const double z0sin1 = Z0SinTheta(z0_1 + dz, eta_1);
        const double z0sin2 = Z0SinTheta(z0_2 + dz, eta_2);
        if (z0sin1 >= z0cut || z0sin2 >= z0cut) continue;

        if (iv == 0 && pass_primary_out) *pass_primary_out = true;

        const double worst = std::max(z0sin1, z0sin2);
        if (best_iv < 0 || worst < best_worst) {
            best_worst = worst;
            best_iv    = static_cast<int>(iv);
        }
    }
    return best_iv;
}

} // namespace AllVertexIP
