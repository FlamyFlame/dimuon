#pragma once

#include "../Utilities/tchain_helpers.h"

template <class Derived>
class PythiaFullSimOverlayExtras {
    template <class, class, class, class...> friend class PythiaAlgCoreT;
    template <class, class, class> friend class PythiaFullSimExtras;

protected:
    Derived& self() { return static_cast<Derived&>(*this); }
    const Derived& self() const { return static_cast<const Derived&>(*this); }

    Float_t overlay_FCal_Et{};
    Float_t overlay_FCal_Et_P{}, overlay_FCal_Et_N{};
    Int_t   overlay_centrality{};
    std::vector<int>* overlay_trk_numqual{nullptr};

    void InitParamsExtra(){
        self().setIsFullsimOverlay(true);
    }

    void InitInputExtra() {
        for (int ikin = 0; ikin < self().nKinRanges; ikin++)
            for (int ibeam = 0; ibeam < self().nBeamTypes; ibeam++) {
                auto* ch = self().evChains_kn_beam[ikin][ibeam];
                if (!ch) continue;
                enable_and_bind(ch, "FCal_Et",    &overlay_FCal_Et);
                enable_and_bind(ch, "FCal_Et_P",  &overlay_FCal_Et_P);
                enable_and_bind(ch, "FCal_Et_N",  &overlay_FCal_Et_N);
                enable_and_bind(ch, "centrality",  &overlay_centrality);
                enable_and_bind(ch, "trk_numqual", &overlay_trk_numqual);
            }
    }

    template <class MuonT>
    void FillMuonOverlay(MuonT& muon) {
        if constexpr (requires { muon.ev_centrality; muon.ev_FCal_Et; }) {
            muon.ev_centrality = overlay_centrality;
            muon.ev_FCal_Et    = overlay_FCal_Et * 1e-6f;
        }
    }

    void FillPairOverlay() {
        auto& pair = *self().mpairRef();
        if constexpr (requires { pair.FCal_Et; pair.avg_centrality; }) {
            pair.FCal_Et   = overlay_FCal_Et   * 1e-6f;
            pair.FCal_Et_A = overlay_FCal_Et_P * 1e-6f;
            pair.FCal_Et_C = overlay_FCal_Et_N * 1e-6f;
            pair.year = 24;
            if (overlay_trk_numqual && overlay_trk_numqual->size() >= 8) {
                pair.ntrk_HIloose         = (*overlay_trk_numqual)[2];
                pair.ntrk_HItight         = (*overlay_trk_numqual)[3];
                pair.ntrk_HIloose_noPtCut = (*overlay_trk_numqual)[6];
                pair.ntrk_HItight_noPtCut = (*overlay_trk_numqual)[7];
            }
        }
    }
};
