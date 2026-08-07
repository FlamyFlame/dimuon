#ifndef MC_TRIG_EFF_PAIR_PT_BINNING_H
#define MC_TRIG_EFF_PAIR_PT_BINNING_H

#include <string>
#include <vector>

#include <TSystem.h>

#include "../MuonObjectsParamsAndHelpers/ParamsSet.h"

// =================================================================================================
// COARSE pair-pT binning FOR THE MC TRIGGER-EFFICIENCY Steps 2/3/4 -- selection + suffix in ONE
// place, so the FILL stage and the PLOT stage can never disagree about which binning is in use.
//
// WHY THIS HEADER EXISTS AT ALL. Two coexisting coarse pair-pT binnings once ran side by side for
// two weeks (R17): the plateau tables and the panel beside them described DIFFERENT cells, every
// plot rendered, every number came out. The only robust defence is that the binning AND the
// output name are derived from a single switch that both stages read. Never re-derive either.
//
// NOMINAL (default, no environment variable): `ParamsSet::pair_pt_coarse_bins`
//   = 8 logarithmic bins, 8 -> 150 GeV. This is the advisor's requested default: pair pT is the
//   key observable and 4 bins smear its dependence out of the efficiency correction.
//
// COMPARISON VARIANT (`MCTRIGEFF_PAIRPT_4BIN=1`): `ParamsSet::pair_pt_coarse_bins_4bin`
//   = 4 logarithmic bins over the SAME 8 -> 150 GeV range. Its purpose is a like-for-like
//   comparison against the 8-bin default with EVERYTHING else held fixed (same samples, same gap
//   cut, same q*eta binning, same plateau window, same fit methods), so the only difference is the
//   number of pair-pT cells. It answers whether the finer binning buys resolution or just noise:
//   on 8 bins the top cells run past where pp Pythia has yield (12 of 72 pp Step-3 cells fail the
//   plateau guard), and the 4-bin run is the control that shows what those cells look like when
//   merged.
//
//   It is OPT-IN and SUFFIXED, never a silent second default (.claude/CLAUDE.md §Binnings item 4):
//   histograms go to `..._pt4bin.root` and plots to a `pt4bin/` subdirectory, so a 4-bin run can
//   neither overwrite nor be mistaken for the nominal.
//
// SCOPE: Steps 2, 3 and 4 only -- the pair-pT-binned views. Step 1 (the eps_MC(pT,q*eta) map) and
// the sanity check carry no pair-pT binning and are unaffected.
// =================================================================================================
namespace MCTrigEffPairPt {

// The switch. Presence of the variable is what counts, matching MCTRIGEFF_NO_GAPCUT.
inline bool UseFourBin() {
    return gSystem->Getenv("MCTRIGEFF_PAIRPT_4BIN") != nullptr;
}

// The edges actually used by the Step-3/Step-4 TH3D pair-pT axis.
inline std::vector<double> Edges(const ParamsSet& pms) {
    return UseFourBin() ? pms.pair_pt_coarse_bins_4bin : pms.pair_pt_coarse_bins;
}

inline int NBins() {
    return UseFourBin() ? ParamsSet::N_COARSE_PAIR_PT_BINS_4BIN
                        : ParamsSet::N_COARSE_PAIR_PT_BINS;
}

// Output-file token. Empty for the nominal so existing nominal filenames are untouched.
inline std::string FileSuffix() {
    return UseFourBin() ? "_pt4bin" : "";
}

// Plot subdirectory, appended AFTER the step directory and BEFORE the working-point directory,
// so the 4-bin set mirrors the nominal tree instead of interleaving with it.
inline std::string PlotSubdir() {
    return UseFourBin() ? "pt4bin/" : "";
}

// One line for the log, so a run always states which binning it used.
inline std::string Describe() {
    return UseFourBin() ? "COMPARISON: 4 log pair-pT bins (8-150 GeV), output tagged _pt4bin"
                        : "NOMINAL: 8 log pair-pT bins (8-150 GeV)";
}

}  // namespace MCTrigEffPairPt

#endif  // MC_TRIG_EFF_PAIR_PT_BINNING_H
