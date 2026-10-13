# Note Review Log
**Task**: Update IntNotes/tex/reconstruction_efficiency.tex: (1) 3D bin-by-bin decision (bin-by-bin efficiency + purity as criteria); (2) reco-eff procedure/results with gap cuts and (pT, eta^pair, m) reco-binned eps^reco_j = N[R_j]/N[T_j] incl. feed-in; vs-reco plots only; within-cell eta^pair shape systematic (option A).
**Log file**: review-note-20261013-002420-reco-eff-bin-by-bin.md
**Started**: 2026-10-13T00:24:20-04:00
**Status**: APPROVED
**Iterations completed**: 4
**Max iterations**: 5

## Iteration 1
**Reviewer verdict**: FAIL
**Issues found**: 13 (4 WARNING, 9 INFO)
**Details**:
1. WARNING MC closure "exact" omits the 6 unmeasured cells (105.5–150, |eta|>2.0; 1e-5 integrated) — fixed (text :56, :407, caption).
2. WARNING 1D pT migration attributed to mass resolution — fixed (muon momentum resolution).
3. WARNING mass "flat within stat" — fixed (≤1.5 %, lowest mass slightly higher 13–26 GeV).
4. WARNING same-sign PHbox named "efficiency of this section" — fixed (currently applied eps(pT,eta,dR)).
5. INFO delta ranges 1–4 → 0.1–4 %, 0.4–1.5 → 0.2–1.5 % — fixed (note + tracking R4).
6. INFO "fewer than 30" → "at most 30" — fixed.
7. INFO pT "next three" double count — fixed. 8. INFO pT rise also gap collimation — fixed.
9. INFO Application tense / garbled stages sentence — fixed. 10. INFO data_corr caption y quantity — fixed.
11. INFO long sentences; 12. INFO overfull <10pt; 13. INFO RUN2-CROSSCHECK UNVERIFIED (no comparable Run 2 pair eff) — not changed.
**Numerical verification**: all MATCH except closure (PARTIAL, issue 1), mass flatness (issue 3), delta prose ranges (issue 5) — all fixed.

## Iteration 2
**Reviewer verdict**: FAIL
**Issues found**: 5 (2 WARNING, 3 INFO)
**Details**:
1. WARNING unrecovered yield "1e-5 of integrated" wrong (actual 3.7e-7; 1e-5 = integrated-panel per-interval residual); sentence order broke "it does not" — fixed (4e-7, rewritten).
2. WARNING mass "≤1.5 % at fixed pT" false above 37 GeV — fixed (below 37 GeV; above within MC stat).
3. INFO closure caption: integrated panel 1e-5 (pairs >150 GeV clamped) — fixed.
4. INFO RUN2-CROSSCHECK UNVERIFIED (acceptable). 5. INFO G7 consistent.
**Numerical verification**: all MATCH except the two WARNING numbers (fixed).

## Iteration 3
**Reviewer verdict**: FAIL
**Issues found**: 4 (1 WARNING, 3 INFO)
**Details**:
1. WARNING closure "exact in every pT interval" false for the two outermost eta intervals at 105.5–150 GeV (−60 %/−50 %) — fixed (text + caption: 40 %/50 % recovered, 4e-7 of total).
2. INFO mass >37 GeV "within" → "compatible with" — fixed. 3. INFO placeholder.md item 3 (already carries the NOT-yet-applied sentence). 4. INFO C3 UNVERIFIED.
**Numerical verification**: fixes 1–3 of it.2 verified (3.67e-7, ≤1.5 %, 1.034e-5); no regressions.

## Iteration 4
**Reviewer verdict**: PASS
**Issues found**: 2 INFO (wording of the 4e-7 missing yield — applied; C3 UNVERIFIED acceptable)
**Numerical verification**: all OK (closure 3.55e-15, 40/50 %, 3.7e-7, 1.034e-5, mass ≤1.5 %).

**Status**: APPROVED at iteration 4
**Summary**: reconstruction_efficiency.tex rewritten for the 3D bin-by-bin reco-binned correction in (pT, eta^pair, m) with migration criteria, reco-quantity figures, within-cell eta^pair shape systematic, and honest not-yet-applied disclosure.
