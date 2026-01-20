template <class Derived>
class PowhegFullSimOverlayExtras{
protected:
    Derived& self() { return static_cast<Derived&>(*this); }
    const Derived& self() const { return static_cast<const Derived&>(*this); }

    void InitParamsExtra(){
        self().is_fullsim_overlay = true;
    }

    void FillMuonPairExtra(int pair_ind);
};