#pragma once
#include "Muon.h"
#include "TLorentzVector.h"
#include "RtypesCore.h"
#include <iostream>
#include <cstring>
#include <cmath>
#include <cassert>
#include <algorithm>   // for std::swap
#include <stdexcept>
#include <type_traits>
#include <utility>   // std::swap

template <class Derived, class MuonT>
struct MuonPairBaseT {
  	using muon_t = MuonT;
  	MuonT m1, m2;

  	// base/common pair quantities
	double weight{1.};

  	Derived& self() { return static_cast<Derived&>(*this); }
  	const Derived& self() const { return static_cast<const Derived&>(*this); }

    void Clear() { self() = Derived{}; }

    // --- base-controlled dispatch ---
    void Sort() {
        // Prefer reco sort if available
        if constexpr (requires(Derived& d) { d.SortTruth(); }) {
            self().SortTruth();
        } else if constexpr (requires(Derived& d) { d.SortReco(); }) {
            self().SortReco();
        } else {
            static_assert([]{ return false; }(),
                          "MuonPairBaseT::Sort(): need SortTruth() or SortReco() in Derived");
        }
    }

    void PairValueCalc() {
        // Call whichever exist (possibly both)
        if constexpr (requires(Derived& d) { d.PairValueCalcTruth(); }) {
            self().PairValueCalcTruth();
        }
        if constexpr (requires(Derived& d) { d.PairValueCalcReco(); }) {
            self().PairValueCalcReco();
        }

        // Optional: enforce at least one exists
        static_assert(
            (requires(Derived& d) { d.PairValueCalcReco(); }) ||
            (requires(Derived& d) { d.PairValueCalcTruth(); }),
            "MuonPairBaseT::PairValueCalc(): need PairValueCalcReco() and/or PairValueCalcTruth() in Derived"
        );
        
        self().PairValueCalcHook();
    }

    void Update(){
        Sort();
        PairValueCalc();
    }

protected:
    void PairValueCalcHook() {}  // no-op by default
};

