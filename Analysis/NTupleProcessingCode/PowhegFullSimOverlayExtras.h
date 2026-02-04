#pragma once

template <class Derived>
class PowhegFullSimOverlayExtras{
    template <class, class, class, class...> friend class PowhegAlgCoreT;
protected:
    Derived& self() { return static_cast<Derived&>(*this); }
    const Derived& self() const { return static_cast<const Derived&>(*this); }

    void InitParamsExtra(){
        self().setIsFullsimOverlay(true);
    }

    void FillMuonPairExtra(int pair_ind);
};