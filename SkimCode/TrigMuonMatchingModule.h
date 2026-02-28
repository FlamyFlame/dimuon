#pragma once

// --- ROOT ---
#include "TLorentzVector.h"
#include "TMath.h"

// --- STL ---
#include <cmath>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// --- xAOD ---
#include "xAODMuon/Muon.h"
#include "xAODMuon/MuonContainer.h"

// --- Athena messaging (so you can get ATH_MSG_* style output) ---
#include "AthenaBaseComps/AthMessaging.h"

// --- Trigger ---
#include "TrigDecisionTool/TrigDecisionTool.h"
#include "TrigDecisionTool/ChainGroup.h"
#include "TrigDecisionTool/FeatureContainer.h"
#include "TrigDecisionTool/TrigDefs.h"

class TrigMuonMatchingModule : public AthMessaging {
public:
  // You pass the *retrieved* TDT pointer from your algorithm.
  TrigMuonMatchingModule(const std::string& name,
                         const Trig::TrigDecisionTool* tdt)
    : AthMessaging(name), m_tdt(tdt) {}

  void setTrigDecisionTool(const Trig::TrigDecisionTool* tdt) { m_tdt = tdt; }

  // ---- This is the one you want (xAOD::Muon version) ----
  bool matchDimuon(const xAOD::Muon* mu1,
                   const xAOD::Muon* mu2,
                   const std::string& chain,
                   std::pair<bool,bool>& result1,
                   std::pair<bool,bool>& result2,
                   const double mindelR = 0.1)
  {
    if (!mu1 || !mu2) return false;

    TLorentzVector mu_1, mu_2;
    mu_1.SetPtEtaPhiM(mu1->pt(), mu1->eta(), mu1->phi(), MUONMASS);
    mu_2.SetPtEtaPhiM(mu2->pt(), mu2->eta(), mu2->phi(), MUONMASS);

    return matchDimuon(mu_1, mu_2, chain, result1, result2, mindelR);
  }

  // ---- TLorentzVector overload (mirrors TrigMuonMatching.cxx) ----
  bool matchDimuon(const TLorentzVector& muon1,
                   const TLorentzVector& muon2,
                   const std::string& chain,
                   std::pair<bool,bool>& result1,
                   std::pair<bool,bool>& result2,
                   const double mindelR = 0.1)
  {
    if (!m_tdt) {
      ATH_MSG_ERROR("TrigDecisionTool pointer is null. Did you retrieve it in your algorithm?");
      return false;
    }

    DimuonChainInfo chainInfo(chain);
    if (!decodeDimuonChain(chainInfo)) {
      ATH_MSG_ERROR("Failed to decode chain " << chain
                    << " (only supports names like HLT_2muXX or HLT_muXX_mu8noL1).");
      return false;
    }

    std::pair<bool,bool> rc12, rc21;
    rc12 = matchDimuon(muon1, muon2, chainInfo, mindelR);

    if (chainInfo.isSymmetric) {
      rc21.first  = rc12.second;
      rc21.second = rc12.first;
    } else {
      rc21 = matchDimuon(muon2, muon1, chainInfo, mindelR);
    }

    // Same assignment as upstream:
    result1.first  = rc12.first;   result1.second = rc21.second;
    result2.first  = rc21.first;   result2.second = rc12.second;
    return true;
  }

private:
  // Matches (muon1 -> threshold1) and (muon2 -> threshold2) by EF-feature matching.
  std::pair<bool,bool> matchDimuon(const TLorentzVector& muon1,
                                  const TLorentzVector& muon2,
                                  const DimuonChainInfo& chainInfo,
                                  const double mindelR)
  {
    EFmuon trkId1, trkId2, dummy;

    const double dr1 = matchedTrackDetail(trkId1, dummy,
                                          muon1.Eta(), muon1.Phi(),
                                          mindelR, chainInfo.thresholds.first);

    const double dr2 = matchedTrackDetail(trkId2, dummy,
                                          muon2.Eta(), muon2.Phi(),
                                          mindelR, chainInfo.thresholds.second);

    // Upstream “deduplicate if both matched to identical EF muon”
    if (trkId1.valid && trkId2.valid &&
        isEqual(trkId1.pt,  trkId2.pt)  &&
        isEqual(trkId1.eta, trkId2.eta) &&
        isEqual(trkId1.phi, trkId2.phi))
    {
      if (dr1 > dr2) {
        matchedTrackDetail(trkId1, trkId2,
                           muon1.Eta(), muon1.Phi(),
                           mindelR, chainInfo.thresholds.first);
      } else {
        matchedTrackDetail(trkId2, trkId1,
                           muon2.Eta(), muon2.Phi(),
                           mindelR, chainInfo.thresholds.second);
      }
    }

    return {trkId1.valid, trkId2.valid};
  }

  // Pulls EF muons from TDT features for `chainForEventTrigger` and finds closest in dR,
  // skipping `usedEFMuonId` to prevent double-use.
  double matchedTrackDetail(EFmuon& efMuonId,
                            const EFmuon& usedEFMuonId,
                            const double eta,
                            const double phi,
                            const double mindelR,
                            const std::string& chainForEventTrigger) const
  {
    efMuonId.valid = false;
    double drmin = mindelR;

    auto cg = m_tdt->getChainGroup(chainForEventTrigger);
    auto fc = cg->features();

#if defined(XAOD_STANDALONE) || defined(XAOD_ANALYSIS)
    auto muFeatures = fc.containerFeature<xAOD::MuonContainer>();
    for (const auto& mucont : muFeatures) {
      const xAOD::MuonContainer* cont = mucont.cptr();
      if (!cont) continue;
      for (const xAOD::Muon* mu : *cont) {
        if (!mu) continue;

        if (isEqual(usedEFMuonId.pt,  mu->pt())  &&
            isEqual(usedEFMuonId.eta, mu->eta()) &&
            isEqual(usedEFMuonId.phi, mu->phi())) continue;

        const double dr = dR(eta, phi, mu->eta(), mu->phi());
        if (drmin > dr) {
          drmin = dr;
          efMuonId.pt = mu->pt();
          efMuonId.eta = mu->eta();
          efMuonId.phi = mu->phi();
          efMuonId.valid = true;
        }
      }
    }
#else
    const std::vector< Trig::Feature<xAOD::MuonContainer> > muFeatures =
      fc.get<xAOD::MuonContainer>();

    for (const auto& mucont : muFeatures) {
      const xAOD::MuonContainer* cont = mucont.cptr();
      if (!cont) continue;
      for (const xAOD::Muon* mu : *cont) {
        if (!mu) continue;

        if (isEqual(usedEFMuonId.pt,  mu->pt())  &&
            isEqual(usedEFMuonId.eta, mu->eta()) &&
            isEqual(usedEFMuonId.phi, mu->phi())) continue;

        const double dr = dR(eta, phi, mu->eta(), mu->phi());
        if (drmin > dr) {
          drmin = dr;
          efMuonId.pt = mu->pt();
          efMuonId.eta = mu->eta();
          efMuonId.phi = mu->phi();
          efMuonId.valid = true;
        }
      }
    }
#endif

    return drmin;
  }

  static double dR(const double eta1, const double phi1,
                   const double eta2, const double phi2)
  {
    const double deta = std::fabs(eta1 - eta2);
    const double dphi_raw = std::fabs(phi1 - phi2);
    const double dphi = (dphi_raw < TMath::Pi()) ? dphi_raw : (2.0*TMath::Pi() - dphi_raw);
    return std::sqrt(deta*deta + dphi*dphi);
  }

  static void tokenize(const std::string& str,
                       std::vector<std::string>& tokens,
                       const std::string& delims)
  {
    tokens.clear();
    std::string::size_type lastPos = str.find_first_not_of(delims, 0);
    std::string::size_type pos     = str.find_first_of(delims, lastPos);

    while (pos != std::string::npos || lastPos != std::string::npos) {
      tokens.push_back(str.substr(lastPos, pos - lastPos));
      lastPos = str.find_first_not_of(delims, pos);
      pos     = str.find_first_of(delims, lastPos);
    }
  }

  bool decodeDimuonChain(DimuonChainInfo& chainInfo)
  {
    chainInfo.isValid = false;

    // cache (same as upstream)
    auto it = m_DimuonChainMap.find(chainInfo.chain);
    if (it != m_DimuonChainMap.end()) {
      chainInfo = it->second;
      return chainInfo.isValid;
    }

    std::vector<std::string> tokens;
    tokenize(chainInfo.chain, tokens, "_");
    if (tokens.size() < 2) return false;
    if (tokens[0] != "HLT") return false;

    chainInfo.isSymmetric = (tokens[1].substr(0,3) == "2mu");

    if (chainInfo.isSymmetric) {
      // e.g. HLT_2mu14 -> thresholds are both HLT_mu14
      const std::string thr = std::string("HLT_" + tokens[1].substr(1));
      chainInfo.thresholds.first  = thr;
      chainInfo.thresholds.second = thr;
      chainInfo.isValid = true;
    } else {
      // e.g. HLT_mu22_mu8noL1
      if (tokens.size() != 3) return false;
      const std::string high = std::string("HLT_" + tokens[1]);
      chainInfo.thresholds.first  = high;
      chainInfo.thresholds.second = chainInfo.chain;
      chainInfo.isValid = true;
      // upstream returns immediately here
      return chainInfo.isValid;
    }

    m_DimuonChainMap[chainInfo.chain] = chainInfo;
    return chainInfo.isValid;
  }

  static bool isEqual(const double x, const double y)
  {
    return (std::fabs(x - y) < std::numeric_limits<float>::epsilon());
  }

private:
  // Same constant as upstream file.
  static constexpr double MUONMASS = 105.65837;

  struct EFmuon {
    bool  valid{false};
    float pt{-1.e30f};
    float eta{-1.e30f};
    float phi{-1.e30f};
  };

  struct DimuonChainInfo {
    std::string chain;
    std::pair<std::string,std::string> thresholds;
    bool isEFFS{false};
    bool isSymmetric{false};
    bool isValid{false};

    DimuonChainInfo(const std::string& chain_ = "") : chain(chain_) {}
  };

  const Trig::TrigDecisionTool* m_tdt{nullptr};
  std::map<std::string, DimuonChainInfo> m_DimuonChainMap;
};
