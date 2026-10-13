# Note Review Log
**Task**: Update IntNotes/tex/reconstruction_efficiency.tex for user decision D7 (2026-10-13): numerator = truth-matched && reco cuts && truth in T; feed-in = template-fit fake; regenerated numbers/figures; D7 not yet propagated to crossx.
**Log file**: review-note-20261013-165550-reco-eff-d7-numerator.md
**Started**: 2026-10-13T16:55:50-04:00
**Status**: APPROVED
**Iterations completed**: 3
**Max iterations**: 5

## Iteration 1
**Reviewer verdict**: FAIL
**Issues found**: 9 (0 CRITICAL, 3 WARNING, 6 INFO)
**Details**:
1. WARNING: SF framing — propagation will be D7 × SF1SF2 (expected change +1.3 %, not −7.1 %); planned-application weight lacks SFs; PHtext implies no Run 3 SFs (pp24 SFs are applied).
2. WARNING: |η^pair_gen| < 2.2 is in T in code (cell grid) but the note says it is not part of the truth signal region; define measured region = truth signal region ∩ |η^pair_gen| < 2.2.
3. WARNING: top-mass dip 3.5–5.5 % reco-binned vs 3.3–3.8 % generated (not "as is"); 26–37 GeV "within 3 %" is 3.2 %.
4. INFO: integrated closure outermost intervals 9e-6. 5. INFO: GRL vintage of data figures not stated. 6. INFO: long sentences. 7. INFO: Pb+Pb PHbox stale (out of scope). 8. INFO: RUN2-CROSSCHECK UNVERIFIED. 9. INFO: stages figure dashed line / axis label overlap (pre-existing).
**Numerical verification**: all MATCH except the mass-dip/3.2 % items (Issue 3) and integrated closure partial (Issue 4).
**Amendments (iteration 1)**: W1 — Application states SF-weighted N[R_j] in the crossx; PHbox leads with the +1.3 % propagation effect (0.592→0.584), the 7.1 % / 4–13 % no-SF comparison as reference; stages caption +1.3 %; Simulation-modelling item: pp24 SFs applied, PHtext = Pb+Pb SFs missing. W2 — feed-in = outside the measured region (truth signal region ∩ grid |η^pair|<2.2, pT<150); |η^pair|<2.2 removed from the reco-only gap list and stated to apply on both legs. W3 — 3.2 %; dip 3.5–5.5 % reco-binned vs 3–4 % generated. INFO 4 (9e-6) and 5 (GRL in captions) also applied.

## Iteration 2
**Reviewer verdict**: FAIL
**Issues found**: 7 (0 CRITICAL, 1 WARNING, 6 INFO)
**Details**:
1. WARNING: mass paragraph "Above 37 GeV ... compatible with statistics" contradicts the 3.5–5.5 % dip claim below 52 GeV (36.7–52 dip 5.5 %, ~3σ). Fix: "Above 52 GeV".
2. INFO: integrated closure 9e-6 only in [−2.2,−2.0] (+2.9e-6 in [2.0,2.2]). 3. INFO: two places "outside the truth signal region" → measured region. 4. INFO: 0.40 → 0.39. 5. INFO: SF handling below 5 GeV / SF systematic not described (propagation stream). 6. INFO: Pb+Pb subsection stale (out of scope). 7. INFO: PHtext missing period.
**Numerical verification**: all MATCH.
**Amendments (iteration 2)**: 1 fixed (52 GeV); INFO 2, 3, 4, 7 fixed; 5, 6 left to the propagation work stream. Recompiled rc 0.

## Iteration 3
**Reviewer verdict**: PASS
**Issues found**: 4 (all INFO)
**Details**:
1. INFO: data_corr caption implied feed-in was the only excess (background too) — FIXED.
2. INFO: +4.5(3) half-way rounding of +0.0455 (kept, matches R4). 3. INFO: 0.39 vs 0.3951 (accepted). 4. INFO: RUN2-CROSSCHECK UNVERIFIED.
**Numerical verification**: all MATCH.

**Status**: APPROVED at iteration 3
**Summary**: Reco-eff section updated to the D7 numerator (truth in the measured region; feed-in = template-fit fake), all numbers/figures regenerated, crossx status disclosed in a PHbox (+1.3 % on propagation).
