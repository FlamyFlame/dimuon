#pragma once
#include <Rtypes.h>

// =================================================================================================
// Colour convention of the MC-vs-data comparison plot sets.
//
// SINGLE SOURCE, at NAMESPACE scope on purpose: the signal family is drawn by two independent
// producers -- `plot_mc_data_compr_signal.cxx`, which derives from PlotMCDataComprBaseClass, and
// `plot_mc_data_pair_pt_in_eta.cxx`, which is a free-function macro that does not include the
// class at all. Keeping these constants as protected class members made the second macro fail to
// compile, and since `pipeline_pp_crossx.sh` runs under `set -Eeuo pipefail` that killed Stage 8
// before the generic family was ever produced (found in review, 2026-08-25).
//
// SIGNAL (single-b) family, user instruction 2026-08-25:
//   pp data          -> RED
//   Pythia (primary) -> BLACK  (single-b on the OS pad, all-SS on the SS pad: the one curve that
//                               is the data's counterpart on that pad)
//   Pythia, all OS   -> BLUE   (secondary reference; deliberately NOT black, so it can never be
//                               mistaken for the single-b signal)
//   POWHEG           -> GREEN  (added 2026-09-03; the fourth colour, so that inside `signal/`
//                               red/black/blue/green each mean ONE thing across every PNG in the
//                               directory. It is deliberately NOT blue, which already means
//                               "Pythia, all OS" on the neighbouring canvas.)
// The GENERIC family keeps its own palette (`PlotMCDataComprBaseClass::colors`), where the same
// three colours carry DIFFERENT meanings -- the two families are separate directories and
// separate conventions.
// =================================================================================================
namespace McDataComprColors {
    constexpr Color_t kSignalData   = kRed;
    constexpr Color_t kSignalMc     = kBlack;
    constexpr Color_t kSignalMcRef  = kBlue;
    constexpr Color_t kSignalPowheg = kGreen + 2;
}
