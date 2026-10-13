#pragma once
// Event-level data-quality exclusions for the pp NTuple stage (PPExtras::PassEventSelExtra).
//
// pp 2024 (USER DECISION 2026-10-12): the new pp-ref GRL
// SkimCode/xmls/physics_ppref2024_all_good.xml (DetStatus-v144-pro42-02) differs from the GRL
// the skim was made with (physics_2024ppRef_25ns.xml, v132-pro42-01; now
// physics_2024ppRef_25ns_legend_old_lumi_calib.xml) ONLY in run 488456: it starts at LB 221
// instead of 84.  The skim still contains LB 84-220, so they are dropped here.
// NOTE: anything excluded here must ALSO be absent from the luminosity in
// PPBaseClass.h::CrossxFactorMap() (395.836 pb^-1 = lumitable from the new GRL), or the
// cross-section numerator and denominator cover different LBs.
#include <vector>

struct PPBadLBRange { int run; unsigned lb_lo; unsigned lb_hi; };   // inclusive

inline std::vector<PPBadLBRange> PPBadLBRanges(int run_year) {
    if (run_year % 2000 == 24) return { {488456, 84u, 220u} };
    return {};
}

inline bool PPIsBadRunLB(const std::vector<PPBadLBRange>& bad, int run, unsigned lb) {
    for (const auto& r : bad)
        if (r.run == run && lb >= r.lb_lo && lb <= r.lb_hi) return true;
    return false;
}
