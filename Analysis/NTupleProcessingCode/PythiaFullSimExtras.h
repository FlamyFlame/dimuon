#pragma once

template <class PairT, class MuonT, class Derived>
class PythiaFullSimExtras {
    template <class, class, class, class...> friend class PythiaAlgCoreT;
protected:
    Derived& self() { return static_cast<Derived&>(*this); }
    const Derived& self() const { return static_cast<const Derived&>(*this); }
    void InitInputExtra() {}
};
