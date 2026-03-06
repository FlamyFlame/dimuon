/*
  Self-contained dimuon trigger matching module for TrigRates.
  Replicates the matchDimuon logic from the Athena-only TrigMuonMatching
  library, using only APIs available in AthAnalysis.

  For HF_IS_R25 (Run3 / AthAnalysis 25.x):
    Per-leg matching is done via IMatchingTool::match(particle, chain, dR).
    decodeDimuonChain maps:
      HLT_2muX_<L1seed>           -> symmetric, both legs use the full chain
      HLT_muX_muYnoL1_<L1seed>   -> asymmetric,
                                     leg1 = "HLT_muX_<L1seed>"  (seeded proxy)
                                     leg2 = full chain            (noL1 proxy)

  For HF_IS_R21 (Run2 / AthAnalysis 21.x):
    matchedTrackDetail uses TDT feature access + explicit disambiguation,
    mirroring the original TrigMuonMatching::matchedTrackDetail.
    decodeDimuonChain maps:
      HLT_2muX                  -> symmetric, both legs = "HLT_muX"
      HLT_muX_muYnoL1           -> asymmetric,
                                   leg1 = "HLT_muX"
                                   leg2 = full chain
*/

#include "HFtrigValidation/Module_TrigMuonMatching.h"

#if defined(HF_IS_R21)
  #include "xAODMuon/MuonContainer.h"
  // TrigDefs for alsoDeactivateTEs flag
  #include "TrigSteeringEvent/TrigDefs.h"
#endif

#include <sstream>
#include <iostream>
#include <algorithm>

// muon mass in MeV
static constexpr double MUON_MASS_MEV = 105.6583755;

// ============================================================
// Constructor
// ============================================================

#if defined(HF_IS_R21) || defined(HF_IS_R25)
TrigMuonMatchingModule::TrigMuonMatchingModule(
    Trig::IMatchingTool*    matchTool,
    Trig::TrigDecisionTool* trigDecTool)
  : m_matchTool(matchTool)
  , m_trigDecTool(trigDecTool)
{}
#endif

// ============================================================
// decodeDimuonChain
// ============================================================

bool TrigMuonMatchingModule::decodeDimuonChain(TrigMuonDimuonChainInfo& chainInfo)
{
  chainInfo.isValid = false;

  // --- check the cache first -------------------------------------------
  auto it = m_DimuonChainMap.find(chainInfo.chain);
  if (it != m_DimuonChainMap.end()) {
    chainInfo = it->second;
    return chainInfo.isValid;
  }

  // --- tokenise on '_' --------------------------------------------------
  std::vector<std::string> tokens;
  tokenize(chainInfo.chain, tokens, "_");

  if (tokens.size() < 2)         return false;
  if (tokens[0] != "HLT")        return false;

  chainInfo.isSymmetric = (tokens[1].substr(0, 3) == "2mu");

  if (chainInfo.isSymmetric) {
    // -----------------------------------------------------------------------
    // Symmetric: HLT_2muX[_L1seed...]
    // -----------------------------------------------------------------------
#if defined(HF_IS_R25)
    // For R25, use the full dimuon chain as the threshold for both legs.
    // IMatchingTool::match(particle, chain, dR) returns true if the particle
    // is associated with ANY leg, which is what we want here.
    chainInfo.thresholds.first  = chainInfo.chain;
    chainInfo.thresholds.second = chainInfo.chain;
#else
    // For R21: threshold = "HLT_muX" (drop the "2" prefix)
    chainInfo.thresholds.first  = "HLT_" + tokens[1].substr(1);
    chainInfo.thresholds.second = "HLT_" + tokens[1].substr(1);
#endif
    chainInfo.isValid = true;

  } else {
    // -----------------------------------------------------------------------
    // Asymmetric: HLT_muX_muYnoL1[_L1seed...]
    // -----------------------------------------------------------------------
#if defined(HF_IS_R25)
    // Run3 naming: HLT_mu4_mu4noL1_L1MU3V -> tokens = [HLT,mu4,mu4noL1,L1MU3V]
    // leg1 (seeded) = "HLT_mu4_L1MU3V"  (standalone chain that exists in menu)
    // leg2 (noL1)   = full chain          (proxy for the full-scan leg)
    if (tokens.size() < 3) return false;
    chainInfo.thresholds.first  = "HLT_" + tokens[1] + "_" + tokens.back();
    chainInfo.thresholds.second = chainInfo.chain;
#else
    // Run2 naming: HLT_muX_muYnoL1 -> exactly 3 tokens
    if (tokens.size() != 3) return false;
    chainInfo.thresholds.first  = "HLT_" + tokens[1]; // e.g. "HLT_mu4"
    chainInfo.thresholds.second = chainInfo.chain;     // full chain as leg2 proxy
#endif
    chainInfo.isValid = true;
  }

  m_DimuonChainMap[chainInfo.chain] = chainInfo;
  return chainInfo.isValid;
}

// ============================================================
// matchedTrackDetail  (R21 only)
// ============================================================

#if defined(HF_IS_R21)
double TrigMuonMatchingModule::matchedTrackDetail(
    TrigMuonEFmuon&       efMuonId,
    const TrigMuonEFmuon& usedEFMuonId,
    double                eta,
    double                phi,
    double                mindelR,
    const std::string&    chainForEventTrigger) const
{
  efMuonId.valid = false;
  double drmin   = mindelR;

  auto cg = m_trigDecTool->getChainGroup(chainForEventTrigger);
  auto fc = cg->features();

#if defined(XAOD_STANDALONE) || defined(XAOD_ANALYSIS)
  auto MuFeatureContainers = fc.containerFeature<xAOD::MuonContainer>();
#else
  const std::vector< Trig::Feature<xAOD::MuonContainer> >
      MuFeatureContainers = fc.get<xAOD::MuonContainer>();
#endif

  for (auto mucont : MuFeatureContainers) {
    for (auto mu : *mucont.cptr()) {

      // skip the trigger object already claimed by the other offline muon
      if (isEqual(usedEFMuonId.pt,  mu->pt())  &&
          isEqual(usedEFMuonId.eta, mu->eta()) &&
          isEqual(usedEFMuonId.phi, mu->phi())) continue;

      double dr = dR(eta, phi, mu->eta(), mu->phi());
      if (drmin > dr) {
        drmin          = dr;
        efMuonId.pt    = mu->pt();
        efMuonId.eta   = mu->eta();
        efMuonId.phi   = mu->phi();
        efMuonId.valid = true;
      }
    }
  }
  return drmin;
}
#endif // HF_IS_R21

// ============================================================
// matchDimuon_inner  — ordered assignment (mu1->threshold1, mu2->threshold2)
// ============================================================

std::pair<Bool_t, Bool_t>
TrigMuonMatchingModule::matchDimuon_inner(
    const TLorentzVector& lv1, const xAOD::Muon* xmu1,
    const TLorentzVector& lv2, const xAOD::Muon* xmu2,
    const TrigMuonDimuonChainInfo& chainInfo,
    double mindelR)
{
#if defined(HF_IS_R25)
  // -----------------------------------------------------------------------
  // R25: delegate to IMatchingTool::match(particle, threshold_chain, dR)
  // -----------------------------------------------------------------------
  (void)lv1; (void)lv2; // TLorentzVectors not needed here

  bool mu1_ok = m_matchTool->match(*xmu1, chainInfo.thresholds.first,  mindelR);
  bool mu2_ok = m_matchTool->match(*xmu2, chainInfo.thresholds.second, mindelR);

  return {static_cast<Bool_t>(mu1_ok), static_cast<Bool_t>(mu2_ok)};

#elif defined(HF_IS_R21)
  // -----------------------------------------------------------------------
  // R21: matchedTrackDetail with disambiguation
  // -----------------------------------------------------------------------
  (void)xmu1; (void)xmu2; // xAOD pointers not needed here

  TrigMuonEFmuon trkId1, trkId2, dummy;

  const double dr1 = matchedTrackDetail(trkId1, dummy,
                                        lv1.Eta(), lv1.Phi(),
                                        mindelR,
                                        chainInfo.thresholds.first);

  const double dr2 = matchedTrackDetail(trkId2, dummy,
                                        lv2.Eta(), lv2.Phi(),
                                        mindelR,
                                        chainInfo.thresholds.second);

  // Disambiguation: if both offline muons matched to the same trigger object,
  // re-match the one that was farther away, excluding the closer match.
  if (trkId1.valid && trkId2.valid &&
      isEqual(trkId1.pt,  trkId2.pt)  &&
      isEqual(trkId1.eta, trkId2.eta) &&
      isEqual(trkId1.phi, trkId2.phi))
  {
    if (dr1 > dr2) {
      matchedTrackDetail(trkId1, trkId2,
                         lv1.Eta(), lv1.Phi(),
                         mindelR, chainInfo.thresholds.first);
    } else {
      matchedTrackDetail(trkId2, trkId1,
                         lv2.Eta(), lv2.Phi(),
                         mindelR, chainInfo.thresholds.second);
    }
  }

  return {static_cast<Bool_t>(trkId1.valid),
          static_cast<Bool_t>(trkId2.valid)};

#else
  // Neither R21 nor R25 — should not be reached
  (void)lv1; (void)lv2; (void)xmu1; (void)xmu2; (void)mindelR;
  return {false, false};
#endif
}

// ============================================================
// matchDimuon  (public interface)
// ============================================================

Bool_t TrigMuonMatchingModule::matchDimuon(
    const xAOD::Muon*          mu1,
    const xAOD::Muon*          mu2,
    const std::string&         chain,
    std::pair<Bool_t, Bool_t>& result1,
    std::pair<Bool_t, Bool_t>& result2,
    double                     mindelR)
{
  // Build 4-vectors
  TLorentzVector lv1, lv2;
  lv1.SetPtEtaPhiM(mu1->pt(), mu1->eta(), mu1->phi(), MUON_MASS_MEV);
  lv2.SetPtEtaPhiM(mu2->pt(), mu2->eta(), mu2->phi(), MUON_MASS_MEV);

  // Decode chain
  TrigMuonDimuonChainInfo chainInfo(chain);
  if (!decodeDimuonChain(chainInfo)) {
    std::cerr << "TrigMuonMatchingModule::matchDimuon : "
                 "Failed to decode chain \"" << chain << "\". "
                 "Accepted formats: HLT_2muX[_L1seed] and "
                 "HLT_muX_muYnoL1[_L1seed]." << std::endl;
    result1 = {false, false};
    result2 = {false, false};
    return false;
  }

  // rc12: ordered assignment  mu1 -> threshold1, mu2 -> threshold2
  // rc21: reversed assignment mu2 -> threshold1, mu1 -> threshold2
  std::pair<Bool_t, Bool_t> rc12, rc21;

  rc12 = matchDimuon_inner(lv1, mu1, lv2, mu2, chainInfo, mindelR);

  if (chainInfo.isSymmetric) {
    // For symmetric chains, the reversed assignment is just the swap
    rc21.first  = rc12.second;
    rc21.second = rc12.first;
  } else {
    rc21 = matchDimuon_inner(lv2, mu2, lv1, mu1, chainInfo, mindelR);
  }

  // Combine into the result pairs following the original TrigMuonMatching
  // convention:
  //   result1 = (mu1 passed leg1 in assignment A,  mu1 passed leg2 in assignment B)
  //   result2 = (mu2 passed leg1 in assignment B,  mu2 passed leg2 in assignment A)
  result1.first  = rc12.first;   // mu1 as leg1 in assignment (mu1->1, mu2->2)
  result1.second = rc21.second;  // mu1 as leg2 in assignment (mu2->1, mu1->2)
  result2.first  = rc21.first;   // mu2 as leg1 in assignment (mu2->1, mu1->2)
  result2.second = rc12.second;  // mu2 as leg2 in assignment (mu1->1, mu2->2)

  return true;
}

// ============================================================
// Utilities
// ============================================================

double TrigMuonMatchingModule::dR(double eta1, double phi1,
                                   double eta2, double phi2) const
{
  double deta = std::fabs(eta1 - eta2);
  double dphi = std::fabs(phi1 - phi2);
  if (dphi > TMath::Pi()) dphi = 2.0 * TMath::Pi() - dphi;
  return std::sqrt(deta * deta + dphi * dphi);
}

bool TrigMuonMatchingModule::isEqual(double x, double y) const
{
  return std::fabs(x - y) < std::numeric_limits<float>::epsilon();
}

void TrigMuonMatchingModule::tokenize(const std::string&        str,
                                       std::vector<std::string>& tokens,
                                       const std::string&        delimiters) const
{
  tokens.clear();
  std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  std::string::size_type pos     = str.find_first_of(delimiters, lastPos);

  while (std::string::npos != pos || std::string::npos != lastPos) {
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    lastPos = str.find_first_not_of(delimiters, pos);
    pos     = str.find_first_of(delimiters, lastPos);
  }
}
