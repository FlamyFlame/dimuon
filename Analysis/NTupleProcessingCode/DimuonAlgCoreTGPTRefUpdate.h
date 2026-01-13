template <class PairT, class MuonT, class Derived>
struct DimuonAlgCoreT {
    using pair_t = PairT;
    using muon_t = MuonT;

    PairT* mpair_raw_ptr = nullptr;
    MuonT* muon_raw_ptr  = nullptr;

    Derived& self() { return static_cast<Derived&>(*this); }

    // final wrappers define the flow
    void Initialize() {
        self().Dimuon_InitializeHook();
    }

    void ProcessPair(pair_t& p) {
        self().Dimuon_BeginPairHook(p);
        if (self().Dimuon_PassCutsHook(p)) {
            self().Dimuon_OnPairHook(p);
        }
        self().Dimuon_EndPairHook(p);
    }

    void Finalize() {
        self().Dimuon_FinalizeHook();
    }

    // default no-op hooks (provided via separate default base or in Derived)
    void Dimuon_InitializeHook() {}
    void Dimuon_BeginPairHook(pair_t&) {}
    bool Dimuon_PassCutsHook(pair_t const&) { return true; }
    void Dimuon_OnPairHook(pair_t&) {}
    void Dimuon_EndPairHook(pair_t&) {}
    void Dimuon_FinalizeHook() {}
};
