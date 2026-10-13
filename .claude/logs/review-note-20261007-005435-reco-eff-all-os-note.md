# Note Review Log
**Task**: IntNote reco-eff section to the 2026-10-06 procedure (all OS, origin test, mass axis, unfolding plan) + trigger section up to date
**Log file**: review-note-20261007-005435-reco-eff-all-os-note.md
**Started**: 2026-10-07T00:54:35-04:00
**Status**: ESCALATED at iteration 5 (max iterations; last WARNING fixed post-review)
**Iterations completed**: 5
**Max iterations**: 5

## Iteration 1 (reco section only)
**Reviewer verdict**: FAIL
**Issues found**: 9 (2 CRITICAL, 6 WARNING, 1 INFO group)
**Details**:
1. CRITICAL [METHOD/FACTUAL] "clamps never on a signal-region pair" — 3 signal pairs with pT>150 GeV are pT-clamped into the top cell.
2. CRITICAL [PHYSICS-RESULTS C1] kinematic-edge ΔR cells (e.g. 0.2–0.4 at 36.7–52.2 GeV 0.49±0.05 vs 0.67; 0.4–0.6 at 12.8–18.2 0.537 vs 0.61) shown but unexplained.
3. WARNING [C2] "close-by inefficiency resolved by ΔR" not shown — no small-ΔR loss; compare with 6–13 % not "factor two".
4. WARNING origin test quotes only the cell-level χ²; projection tests significant (104.7/71, 39.1/16; inclusive 3.7σ).
5. WARNING composition "bound" unsupported.
6. WARNING planned template-fit ordering written in present tense.
7. WARNING purity caveat (truth-seeded tree) missing.
8. WARNING writing: long sentences, 21 em-dashes, "0.817 of the loss".
9. INFO: crossx file mtime 10-07 00:58; weighted-mean wording; s/p bias estimate; collimation claim untested; stage-figure axis title overlap.
**Numerical verification**: all MATCH except the clamp statement.
**C1 disposition (Step 5b)**: EXPECTED PHYSICS — named kinematic threshold. The flagged cells are exactly the last populated ΔR interval of each pT bin, at the documented bound ΔR_max = m/√(pT,min(pT−pT,min)) (0.64 @9, 0.47 @12.8, 0.24 @36.7 GeV), where one muon sits at the 4.5 GeV threshold and the mass at the window top → reco migration out of the signal region; the bin-migration study (pair_bin_migration_bin_by_bin_feasibility.md R2) independently finds the lowest purity / largest outside feed-in in the same edge cells. No investigation needed; justification written into the note and given to the next reviewer.
**Amendments**: clamp statement corrected (3 of 516 649 pairs); edge-cell mechanism paragraph with examples + cross-ref to migration purity; "no close-by loss visible at ΔR<0.2 granularity"; mass vs ΔR comparison 6–13 % (two lowest pT bins) and ×2 incl. edge; projection-level significance added (104.7/71, 39.1/16, 3.7σ), independence stated as cell-level; composition item → Δf × 0.006 estimate; planned-procedure tense + PHbox rewritten (present per-pair weighting of both signs ≡ per-cell correction of OS−SS); purity caveat; 0.817 factorisation wording; em-dashes 21→10; long sentences split; source-comment crossx date note. Rebuilt: 43 pp, 0 undefined, 0 overfull.

## Iteration 2 (reco + trigger)
**Reviewer verdict**: FAIL
**Issues found**: 6 (2 WARNING, 4 INFO)
**Details**:
1. WARNING mass paragraph: inserted ΔR comparison broke logic ("it" w/o subject; "Most of that variation" misreferent; "only ones with three populated ΔR" false; edge already inside "6–13 %"; 6 % MISMATCH (6.7 %); long sentence).
2. WARNING trigger: single-value uncertainty minimum is 5.9 % at 1≤|η^pair|<2 (74–105 GeV), not 7 % barrel.
3. INFO "same efficiency cell by cell" overstated. 4. INFO "not lower" → "not significantly lower". 5. INFO "0.61" referent. 6. INFO near-threshold muon intrinsically less efficient.
**Numerical verification**: all MATCH except the 6–13 % ΔR spread and the 7 % minimum.
**Amendments**: mass paragraph reordered (3.8–5.1 % → top-bin dip/migration; long sentence split); one ΔR comparison with stated (max−min)/mean convention: 7 % across the first three ΔR intervals at 9–12.8 GeV, edge intervals 3.5–47 % below their neighbour (12.2 %, 27.0 %, 47.4 %, 3.5 %, 10.6 % from the value table); trigger 6 % (1<|η^pair|<2, 74–105 GeV) to 22 %; INFO 3–6 applied. Rebuilt 43 pp clean.

## Iteration 3
**Reviewer verdict**: FAIL
**Issues found**: 6 (3 WARNING, 3 INFO)
**Details**:
1. WARNING forward-η mass slope 17–18 % lacks pT qualifier (16.6/18.3 at the two lowest pT; 13.6/6.5/2.1 next).
2. WARNING Results says ΔR loss only in the last populated interval, but the "7 %" at 9–12.8 GeV is the 0.4–0.6 cell, itself edge-cut above 9.7 GeV → contradiction.
3. WARNING 56-word sentence.
4. INFO dip 3.2–3.8 % / "most of" in five of six. 5. INFO 18–26 → 18.2–25.8; 0.61 → 0.611. 6. INFO mass dependence within a cell not in systematics list.
**Numerical verification**: all MATCH except the unqualified 17–18 %.
**Amendments**: qualifier (9–18.2 GeV; 13.6 % and 6.5 % next two); Results: edge = intervals cut by the bound (last populated + 0.4–0.6 at 9–12.8 GeV, bound < 0.6 above 9.7 GeV); mass comparator now non-edge 0.8 %/1.2 % vs edge 3.5–47 %; sentence split; INFO 4–6 applied (systematics: within-cell mass distribution). Rebuilt 43 pp clean.

## Iteration 4
**Reviewer verdict**: FAIL
**Issues found**: 4 (2 WARNING, 2 INFO)
**Details**:
1. WARNING edge-cut enumeration inconsistent: bound = 0.4 at 16.2 GeV also cuts 0.2–0.4 at 12.8–18.2 GeV above 16.2 GeV, yet that cell was the "away from the edge" 1.2 % comparator.
2. WARNING three changed sentences > 40 words.
3. INFO three names for one concept (bound / edge / boundary). 4. INFO 3.5–47 % covers only edge cells with a neighbour.
**Numerical verification**: all MATCH (slopes, comparators, edge table, bound arithmetic incl. 0.6 at 9.69 GeV, 0.4 at 16.2 GeV).
**Amendments**: bound rule stated generally (cuts the last populated interval of every pT interval and, over part of the pT range, the interval below it); only the clean comparator kept (0.8 % at 9–12.8 GeV where the bound stays > 0.4); sentences split; one term "kinematic bound" throughout. Rebuilt clean.

## Iteration 5
**Reviewer verdict**: FAIL
**Issues found**: 2 (1 WARNING, 1 INFO)
**Details**:
1. WARNING "bound-cut intervals lie 3.5–47 % below their neighbour" contradicts the Results rule: the partially cut 0.2–0.4 at 12.8–18.2 GeV is only 1.2 % below; 3.5–47 % holds for the last populated intervals only.
2. INFO "over part of the pT range, the interval below it" true only in the two lowest pT intervals.
**Numerical verification**: all MATCH (bound 0.644/0.475/0.241; 0.4 at 16.18, 0.2 at 51.2 GeV; 0.8 %; examples; 3.5–47 % for last populated intervals).
**Amendment (applied AFTER the last review, not re-reviewed)**: "The last populated ΔR interval of each pair-pT interval, which the kinematic bound cuts most deeply, lies 3.5–47 % below its neighbour; the intervals cut only over part of their range show smaller deficits of 1–7 %" (1.2 %, 6.5 %); Results rule restricted to the two lowest pT intervals. Rebuilt 43 pp, 0 undefined, 0 overfull; SYNC CLEAN.

**Status**: ESCALATED at iteration 5 (max iterations reached)
**Unresolved issues**:
- None outstanding in substance: the single iteration-5 WARNING was fixed after the last review with the reviewer's own suggested wording, but that fix has not been independently re-reviewed. Zero CRITICAL in iterations 2–5.
**Recommended action**: optional one-pass re-review of §Results ΔR paragraph + §mass comparator sentence; otherwise accept.
