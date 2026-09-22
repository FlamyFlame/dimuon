#ifndef PAIR_RECO_EFF_DEFINITION_H
#define PAIR_RECO_EFF_DEFINITION_H

// =================================================================================================
// PairRecoEffDefinition -- the ONE statement of what the `*_pass_signal_truth*` legs of the pair
// reconstruction efficiency mean, shared by the writer (RDFBasedHistFillingPythiaFullsim), the
// builder (plotting_codes/reco_effcy/build_pp24_fullsim_pair_reco_eff.C) and the consumer
// (Utilities/PairRecoEffEvaluator.h).
//
// WHY. On 2026-09-17 the truth DENOMINATOR lost the detector-gap cuts (single-muon q*eta windows +
// |eta^pair| < 2.2), which moved to the RECO leg only so that the gap ACCEPTANCE is part of
// eps_reco (docs/tracking/pair_reco_eff_gap_acceptance.md). The histogram NAMES did not change.
// A file filled before that date therefore still opens, every number still comes out, and every
// one is the old fiducial definition -- the textual-vs-inherited hazard of
// docs/tracking/mu_pt45_gap125_pairpt9_adoption.md R4. So the writer stamps this marker into its
// output, the builder REFUSES an input without it and copies it into the product, and the
// evaluator refuses a product without it. Change the VALUE whenever the definition changes.
// =================================================================================================
namespace PairRecoEffDefinition {

inline const char* Key() { return "pair_reco_eff_definition"; }

inline const char* Value() {
    return "denominator = truth signal cuts only (mass window, pair pT; NO gap cut); "
           "numerator = + reco-matched + WP + reco signal cuts + gap cuts on RECO q*eta and RECO eta^pair; "
           "gap acceptance INCLUDED (2026-09-17, docs/tracking/pair_reco_eff_gap_acceptance.md)";
}

}  // namespace PairRecoEffDefinition

#endif  // PAIR_RECO_EFF_DEFINITION_H
