# Task 03 — Overlay reco-efficiency Bug #1/#2 fixes + pipeline rerun

**Roadmap:** `analysis_roadmap_2026_06.md` §Q3a, Q4 / chain [4]
**Can start:** immediately. **Reviewer:** `/review-analysis-code` (cite
the implementation specs below in the task prompt), then `/review-plot`.
**Authoritative spec:**
`docs/tracking/hijing_overlay_reco_effcy_investigation.md` Steps 8–10 —
read fully first (INVARIANT). Implementation specs for both fixes are
already written there ("Implementation: Bug #1 fix", "Implementation:
Bug #2 fix").

## Objective

Make the overlay (and pp fullsim) reconstruction efficiency physically
meaningful by fixing:
- **Bug #1 (index mapping):** position in filtered
  `real_muon_truth_barcode_list` used to index full `muon_*` vectors →
  60% of matched muons in overlay get wrong kinematics/quality flags.
  Fix: parallel `real_muon_orig_index` vector. File:
  `NTupleProcessingCode/PythiaFullSimExtras.c` (~lines 86–97, 129).
- **Bug #2 (AOD truth-match failure):** 18.3% of fiducial Pythia truth
  muons have no barcode match (HIJING barcode collision at AOD level);
  69% of those have a nearby reco muon (dR < 0.05–0.12). Fix: dR < 0.05
  fallback matching with claimed-index set, **gated to overlay samples
  only** (pp fullsim's 1.4% unmatched is physical).

Expected after both: single-muon match 81.7% → ~94%; pair ε ~20% →
~80–82% (≈ pp fullsim), with mild physical centrality dependence.

## Steps

1. Read the investigation doc; write plan to its Latest Stage.
2. Implement Bug #1 exactly per spec. Compile with ACLiC
   (`.L PythiaAnalysisClasses.h+` pattern used by the condor script).
3. Implement Bug #2 per spec (overlay-only guard, claimed-index set).
4. `/review-analysis-code` with both spec sections in the prompt.
5. Rerun overlay pipeline on r17618 60k-evt sample
   (`pipelines/pipeline_pythia_fullsim_overlay.sh hijing`; condor per
   user rules — submit, monitor, sanity-check). Also rerun pp fullsim
   test (regression: pp numbers should barely move, ≤1%).
6. Re-plot reco efficiencies; cross-check vs r17662 signalOnlyTruth
   sample (no barcode collisions → Bug #2 path unused there).
7. Optional validation: run the P3-style cross-term dR correction on the
   fixed overlay sample (unbiased MC — mu4 doc Remaining Work item).
8. Update investigation doc (Steps 6'–7'), roadmap ledger. If results
   are physical, this sample becomes the **dummy ε_reco source** for
   task 04.

## Verification

- Pair ε (medium, 0–5% ctr) rises from ~20% to ~55% after Bug #1 alone,
  then to ~80% after Bug #2 (budget table in doc Step 8).
- Single-B vs all-OS efficiencies physically distinct or
  consistently explained.
- pp fullsim regression: match rate stays ≈98.6%, pair ε ≈81.7%.
