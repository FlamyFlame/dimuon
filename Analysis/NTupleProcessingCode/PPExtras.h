#pragma once

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

public:
  ~PPExtras(){}
};
