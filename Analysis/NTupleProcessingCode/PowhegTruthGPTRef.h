struct PowhegTruthExtra {
    template <class Alg, class Pair>
    bool PassCutsExtra(Alg& a, Pair const& p) {
        if (!a.do_truth) return true;
        (void)p;
        return true;
    }

    template <class Alg, class Pair>
    void OnPairExtra(Alg& a, Pair& p) {
        if (!a.do_truth) return;
        (void)a; (void)p;
    }
};

struct PowhegFullSimExtra { /* similar; gated on a.isFullSim */ };
struct PowhegOverlayExtra { /* similar; gated on a.isOverlay */ };

class PowhegFullSimOverlayWithTruth
  : public PowhegAlgCoreT<
        MuonPairPowhegFullSimOverlay, // PairT
        MuonPowhegFullSimOverlay, // MuonT
        PowhegFullSimOverlayWithTruth,   // Derived
        PowhegFullSimExtra<PowhegFullSimOverlayWithTruth>, // Extra1<Derived>
        PowhegOverlayExtra<PowhegFullSimOverlayWithTruth>,  // Extra2<Derived>
        PowhegTruthExtra<PowhegFullSimOverlayWithTruth>  // Extra3<Derived> // ordering chosen here
    >
    , public PowhegTruthExtra<PowhegTruthNTupleFirstPass>
{};

class PowhegTruthNTupleFirstPass
  : public PowhegAlgCoreT<
        MuonPairPowheg, MuonPowheg,
        PowhegTruthNTupleFirstPass,
        PowhegTruthExtra<PowhegTruthNTupleFirstPass>
    >
  , public PowhegTruthExtra<PowhegTruthNTupleFirstPass>
{
public:
    using pair_t = MuonPairPowheg;   // so extras can refer to Derived::pair_t if desired
    bool do_truth{true};
};
