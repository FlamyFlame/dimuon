template <class PairT, class MuonT, class Derived, class... Extras>
struct PowhegAlgCoreT
    : DimuonAlgCoreT<PairT, MuonT, Derived>
    , Extras...
{
    using base_t = DimuonAlgCoreT<PairT, MuonT, Derived>;
    using pair_t = typename base_t::pair_t;

    Derived& self() { return static_cast<Derived&>(*this); }

    // Powheg-common state
    bool do_truth   = false;
    bool isFullSim  = false;
    bool isOverlay  = false;

    // ---- Powheg core work ----
    void Powheg_InitializeCore() {
        // set up branches, load configs, etc.
    }

    bool Powheg_PassCutsCore(pair_t const& p) {
        (void)p;
        return true;
    }

    void Powheg_OnPairCore(pair_t& p) {
        (void)p;
        // fill powheg-common fields, etc.
    }

    void Powheg_FinalizeCore() {
        // write summary, close files, etc.
    }

    // ---- Optional-call helpers (C++20 `requires`; can be done in C++17 too) ----
    template <class E>
    void CallInitExtra() {
        if constexpr (requires(E& e){ e.InitializeExtra(); }) {
            static_cast<E&>(*this).InitializeExtra();
        }
    }

    template <class E>
    bool CallPassCutsExtra(pair_t const& p) {
        if constexpr (requires(E& e, pair_t const& pp){ e.PassCutsExtra(pp); }) {
            return static_cast<E&>(*this).PassCutsExtra(p);
        } else {
            return true;
        }
    }

    template <class E>
    void CallOnPairExtra(pair_t& p) {
        if constexpr (requires(E& e, pair_t& pp){ e.OnPairExtra(pp); }) {
            static_cast<E&>(*this).OnPairExtra(p);
        }
    }

    template <class E>
    void CallFinalizeExtra() {
        if constexpr (requires(E& e){ e.FinalizeExtra(); }) {
            static_cast<E&>(*this).FinalizeExtra();
        }
    }

    // ---- Implement Dimuon hooks by chaining Powheg core + extras ----
    void Dimuon_InitializeHook() {
        Powheg_InitializeCore();
        (CallInitExtra<Extras>(), ...);
    }

    bool Dimuon_PassCutsHook(pair_t const& p) {
        if (!Powheg_PassCutsCore(p)) return false;
        bool ok = true;
        (void)std::initializer_list<int>{ (ok = ok && CallPassCutsExtra<Extras>(p), 0)... };
        return ok;
    }

    void Dimuon_OnPairHook(pair_t& p) {
        Powheg_OnPairCore(p);
        (CallOnPairExtra<Extras>(p), ...);
    }

    void Dimuon_FinalizeHook() {
        (CallFinalizeExtra<Extras>(), ...);
        Powheg_FinalizeCore();
    }
};
