#pragma once

// Placeholder for Pythia PbPb fullsim overlay extras.
// Will add centrality / FCal branches when overlay samples are available.
template <class Derived>
class PythiaFullSimOverlayExtras {
    template <class, class, class, class...> friend class PythiaAlgCoreT;

protected:
    Derived& self() { return static_cast<Derived&>(*this); }
    const Derived& self() const { return static_cast<const Derived&>(*this); }

    void InitParamsExtra(){
        self().setIsFullsimOverlay(true);
    }
    void InitInputExtra() {}
};
