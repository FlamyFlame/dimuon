// PowhegExtrasInterface.h
#pragma once

template <class AlgT>
struct PowhegExtrasInterface {
    using pair_t = typename AlgT::pair_t;
    using muon_t = typename AlgT::muon_t;

    // lifecycle hooks
    void Initialize(AlgT&) {}
    void InitTempVariables(AlgT&) {}
    void InitOutput(AlgT&) {}
    void HistAdjust(AlgT&) {}
    void Finalize(AlgT&) {}

    // analysis hooks
    bool PassCuts(AlgT&, pair_t const&) { return true; }
    void OnBeginEvent(AlgT&) {}
    void OnPair(AlgT&, pair_t&) {}
    void OnEndEvent(AlgT&) {}

    // single-muon hook (optional)
    void OnMuon(AlgT&, muon_t const&) {}
};
