#ifndef __MODULE_TRIGMUONMATCHING_H__
#define __MODULE_TRIGMUONMATCHING_H__

#include "HFtrigValidation/AthenaVersion.h"

#include <string>
#include <map>
#include <vector>
#include <utility>
#include <limits>
#include <cmath>

#include "TMath.h"
#include "TLorentzVector.h"
#include "Rtypes.h"

#include "xAODMuon/MuonContainer.h"

#if defined(HF_IS_R21) || defined(HF_IS_R25)
  #include "TriggerMatchingTool/IMatchingTool.h"
  #include "TrigDecisionTool/TrigDecisionTool.h"
#endif

// -------------------------------------------------------------------------
// Helper structs (in global scope to avoid nested-class forward-decl issues)
// -------------------------------------------------------------------------

struct TrigMuonEFmuon {
  bool  valid = false;
  float pt    = 0.f;
  float eta   = 0.f;
  float phi   = 0.f;
};

struct TrigMuonDimuonChainInfo {
  std::string chain;
  std::pair<std::string, std::string> thresholds; // (leg1_chain, leg2_chain)
  bool isSymmetric = false;
  bool isValid     = false;
  explicit TrigMuonDimuonChainInfo(const std::string& c) : chain(c) {}
  TrigMuonDimuonChainInfo() = default;
};

// -------------------------------------------------------------------------
// TrigMuonMatchingModule
//
// Self-contained dimuon trigger matching utility.  Does NOT inherit from any
// Athena/ASG class; it is a plain C++ object owned by TrigRates.
//
// For HF_IS_R25:
//   Uses IMatchingTool::match(particle, chain, dR) per leg.  No explicit
//   trigger-object disambiguation (not exposed by the R3 API).
//
// For HF_IS_R21:
//   Uses TDT feature-access (containerFeature<MuonContainer>) exactly as in
//   the original TrigMuonMatching library, including the disambiguation step
//   that prevents two offline muons from claiming the same trigger object.
// -------------------------------------------------------------------------

class TrigMuonMatchingModule {
public:

#if defined(HF_IS_R21) || defined(HF_IS_R25)
  // matchTool  : the IMatchingTool already retrieved by TrigRates
  // trigDecTool: only used in R21 for matchedTrackDetail feature access
  TrigMuonMatchingModule(Trig::IMatchingTool*        matchTool,
                         Trig::TrigDecisionTool*     trigDecTool = nullptr);
#else
  TrigMuonMatchingModule() = default;
#endif

  ~TrigMuonMatchingModule() = default;

  // Main interface — mirrors the original TrigMuonMatching::matchDimuon.
  //
  // result1 = {mu1_passed_leg1_in_assignment_A,  mu1_passed_leg2_in_assignment_B}
  // result2 = {mu2_passed_leg1_in_assignment_B,  mu2_passed_leg2_in_assignment_A}
  //
  // For symmetric chains (e.g. HLT_2mu4_L12MU3V):
  //   result1.first == result1.second  (same leg for both assignments)
  //   result2.first == result2.second
  //
  // Returns false if the chain cannot be decoded (logs to stderr).
  Bool_t matchDimuon(const xAOD::Muon*         mu1,
                     const xAOD::Muon*         mu2,
                     const std::string&        chain,
                     std::pair<Bool_t, Bool_t>& result1,
                     std::pair<Bool_t, Bool_t>& result2,
                     double                    mindelR);

private:

  // ---- tool pointers (raw, non-owning) -----------------------------------
#if defined(HF_IS_R21) || defined(HF_IS_R25)
  Trig::IMatchingTool*     m_matchTool    = nullptr;
  Trig::TrigDecisionTool*  m_trigDecTool  = nullptr; // R21 only
#endif

  // ---- decode cache -------------------------------------------------------
  std::map<std::string, TrigMuonDimuonChainInfo> m_DimuonChainMap;

  // ---- private helpers ----------------------------------------------------

  // Parse chain name into leg threshold sub-chains + symmetry flag.
  // Behaviour differs between R21 and R25 (see .cxx).
  bool decodeDimuonChain(TrigMuonDimuonChainInfo& chainInfo);

  // Inner ordered matching: (mu1 vs threshold1, mu2 vs threshold2).
  // For R25: single-particle IMatchingTool::match calls.
  // For R21: matchedTrackDetail with disambiguation.
  std::pair<Bool_t, Bool_t>
  matchDimuon_inner(const TLorentzVector& lv1, const xAOD::Muon* xmu1,
                    const TLorentzVector& lv2, const xAOD::Muon* xmu2,
                    const TrigMuonDimuonChainInfo& chainInfo,
                    double mindelR);

#if defined(HF_IS_R21)
  // Find the closest EF muon in the features of 'chainForEventTrigger' to
  // (eta, phi), excluding the trigger object described by usedEFMuonId.
  // Returns the minimum dR found (always <= mindelR on success).
  double matchedTrackDetail(TrigMuonEFmuon&       efMuonId,
                            const TrigMuonEFmuon& usedEFMuonId,
                            double                eta,
                            double                phi,
                            double                mindelR,
                            const std::string&    chainForEventTrigger) const;
#endif

  // Utility: delta-R between two (eta, phi) points
  double dR(double eta1, double phi1,
            double eta2, double phi2) const;

  // Utility: floating-point equality within float epsilon
  bool isEqual(double x, double y) const;

  // Utility: tokenize string by delimiter characters
  void tokenize(const std::string&        str,
                std::vector<std::string>& tokens,
                const std::string&        delimiters) const;
};

#endif // __MODULE_TRIGMUONMATCHING_H__
