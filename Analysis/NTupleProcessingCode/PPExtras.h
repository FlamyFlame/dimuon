#pragma once

#include <vector>
#include "Riostream.h"
#include "TChain.h"
#include "../MuonObjectsParamsAndHelpers/muon_pair_enums_data.h"

template <class PairT, class MuonT, class Derived, class... Extras>
class DimuonDataAlgCoreT;

template <class Derived>
class PPExtras {
  template <class, class, class, class...> friend class DimuonDataAlgCoreT;

protected:
  Derived& self() { return static_cast<Derived&>(*this); }
  const Derived& self() const { return static_cast<const Derived&>(*this); }

    void InitParamsExtra();
  void PerformTChainFill();

  // ---- ALL-VERTEX impact-parameter selection (pp only) --------------------------------
  // The skim's PrimaryVertices dump, in container order. vtx_z[0] is the primary vertex --
  // the same one the skim used to reference every muon's stored z0.
  std::vector<float>* vtx_z    {nullptr};
  std::vector<int>*   vtx_ntrk {nullptr};

  void InitInputBranchesDimuonAnalysisExtra();
  bool PassD0Z0Extra();

public:
  ~PPExtras(){}
};
