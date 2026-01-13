// PowhegTruthExtra.h
#pragma once
#include "PowhegExtrasInterface.h"

template <class AlgT>
struct PowhegTruthExtra : PowhegExtrasInterface<AlgT> {
    using pair_t = typename AlgT::pair_t;

    TruthQQPair* qqpair{nullptr};

    float crossx_2_to_2{0.f};
    float crossx_2_to_3{0.f};

    // replace your new/delete vectors with owning members if possible:
    std::vector<std::vector<int>> single_gluon_history;
    std::vector<std::vector<int>> m1_history;
    std::vector<std::vector<int>> m2_history;

    void Initialize(AlgT& a) {
        int quark = (a.mc_mode == "bb") ? 5 : 4;
        qqpair = new TruthQQPair(quark);
        // any other “begin job” setup that your old ProcessData did
    }

    void InitTempVariables(AlgT&) {
        single_gluon_history.clear();
        m1_history.clear();
        m2_history.clear();
        // if you truly need pointer vectors exactly as before, keep pointers here
    }

    void OnPair(AlgT& a, pair_t& p) {
        // this replaces PowhegTruthNTupleFirstPass::PerformTruthPairAnalysis()

        // You can call helper methods you move into this extra:
        // MuonPairAncestorTracing(a, p);
        // HardScatteringAnalysis(a, p);
        // etc.
        // Update crossx integrals like your old code
        (void)a;
        (void)p;
    }

    void HistAdjust(AlgT& a) {
        // replaces PowhegTruthNTupleFirstPass::HistAdjust()
        (void)a;
    }

    void Finalize(AlgT& a) {
        // replaces PowhegTruthNTupleFirstPass::Finalize()
        delete qqpair;
        qqpair = nullptr;

        std::cout << "2->2 integral = " << crossx_2_to_2 << "\n";
        std::cout << "2->3 integral = " << crossx_2_to_3 << "\n";

        (void)a;
    }

    // Move your big helper functions here (as methods),
    // changing signature from PowhegTruthNTupleFirstPass::Foo(...)
    // to Foo(AlgT& a, pair_t& p, ...)
};
