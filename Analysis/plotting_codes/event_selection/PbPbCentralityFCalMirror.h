#ifndef PBPB_CENTRALITY_FCAL_MIRROR_H
#define PBPB_CENTRALITY_FCAL_MIRROR_H

// PbPbCentralityFCalMirror.h
//
// Centrality recalculation from FCal ET, for the event-selection plotting macros
// (plot_pbpb_event_sel_cuts.cxx, plot_pbpb_event_sel_event_level.cxx), used for years
// whose 'centrality' branch is unfilled in the skim (pbpb2025, pbpb2026).
//
// ############################################################################
// #  MIRROR — NOT the canonical table.                                        #
// #                                                                           #
// #  The canonical PbPb2023 Glauber FCal->centrality table is                 #
// #      PairPbPbExtras::FCal_ET_Bins_PbPb2023                                #
// #  in MuonObjectsParamsAndHelpers/MuonPairPbPb.h.                           #
// #                                                                           #
// #  The values below were hand-typed from it long ago and 40 of the 85       #
// #  entries are rounded to 5 significant figures, so they DIFFER from the    #
// #  canonical values by up to 8.1e-5 relative (largest at index 82).         #
// #  In particular index 79 — the 0-80% boundary that drives `is_ctr80` —     #
// #  is 0.06321 here versus 0.063208 canonical.                               #
// #                                                                           #
// #  DO NOT "reconcile" these numbers.  Replacing them with the canonical     #
// #  values would perturb the 2023/24/25 event-selection derivation, which    #
// #  is an OPEN USER DECISION, deliberately parked.  It is tracked in         #
// #  docs/tracking/pbpb2026_analysis_support.md and is the user's call, not   #
// #  a cleanup.                                                               #
// #                                                                           #
// #  This header exists solely so there is ONE copy of the mirror instead of  #
// #  two byte-identical hand-typed copies that could drift apart.  The values #
// #  are EXACTLY those that were in both macros before it was created.        #
// ############################################################################

// PbPb2023 Glauber v3.2 table: FCal ET lower bounds (TeV) per integer percentile bin (0-84).
static const float kFCalBinsPbPb2023[85] = {
  4.51272f, 4.32043f, 4.15372f, 3.99602f, 3.84498f, 3.69944f, 3.55802f, 3.42045f, 3.28744f, 3.15972f, // 0-10
  3.03748f, 2.92012f, 2.80723f, 2.69878f, 2.59464f, 2.49406f, 2.39646f, 2.30180f, 2.21028f, 2.12188f, // 10-20
  2.03659f, 1.95428f, 1.87489f, 1.79842f, 1.72484f, 1.65387f, 1.58516f, 1.51853f, 1.45406f, 1.39178f, // 20-30
  1.33168f, 1.27363f, 1.21752f, 1.16336f, 1.11112f, 1.06069f, 1.01210f, 0.96518f, 0.91991f, 0.87632f, // 30-40
  0.83431f, 0.79384f, 0.75493f, 0.71753f, 0.68162f, 0.64714f, 0.61393f, 0.58195f, 0.55126f, 0.52171f, // 40-50
  0.49348f, 0.46640f, 0.44043f, 0.41569f, 0.39200f, 0.36933f, 0.34768f, 0.32701f, 0.30731f, 0.28853f, // 50-60
  0.27064f, 0.25357f, 0.23726f, 0.22178f, 0.20720f, 0.19325f, 0.17983f, 0.16748f, 0.15557f, 0.14435f, // 60-70
  0.13388f, 0.12389f, 0.11468f, 0.10598f, 0.09765f, 0.09015f, 0.08264f, 0.07589f, 0.06955f, 0.06321f, // 70-80
  0.05760f, 0.05273f, 0.04787f, 0.04300f, 0.03885f                                                     // 80-85
};
// Returns integer centrality percentile [0,84], or -1 if FCal ET below the 85th-percentile boundary.
// Identical logic to PairPbPbExtras::GetCentrality() with std::greater binary search.
static int CentralityFromFCal2023(float fcal_AC) {
    const int n = 85;
    // lower_bound with greater<float>: first position where bins[k] <= fcal_AC
    int lo = 0, hi = n;
    while (lo < hi) {
        int mid = (lo + hi) / 2;
        if (kFCalBinsPbPb2023[mid] > fcal_AC) lo = mid + 1;
        else                                   hi = mid;
    }
    return (lo >= n) ? -1 : lo;
}

#endif // PBPB_CENTRALITY_FCAL_MIRROR_H
