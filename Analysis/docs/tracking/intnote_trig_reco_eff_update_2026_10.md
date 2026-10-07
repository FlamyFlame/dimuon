# Internal note: trigger + reconstruction efficiency sections — October 2026 update

**Mode:** Implementation (academic writing). **Opened:** 2026-10-07.
**Predecessor:** `intnote_trig_reco_eff_update_2026_09.md` (DONE 2026-09-18).

## Objective
User request 2026-10-07: "update internal notes (methodology + results) for the updated reco
efficiency procedure. ensure trigger and reco efficiency parts of the internal notes are up to date".
1. `IntNotes/tex/reconstruction_efficiency.tex`: the ε_reco procedure of 2026-10-06 — all truth
   opposite-sign pairs (any origin; was single-b), the origin-independence test (method + result),
   the mass-as-closeness-variable diagnostic, the truth-binned definition + Bayesian-unfolding plan
   (bin-by-bin only if justified; bin-migration study status), every number re-read from the current
   product / crossx.
2. `IntNotes/tex/trigger_efficiency.tex`: verified against current code + outputs (no trigger code or
   product changed since 2026-09-18; check every quoted crossx-derived number and status statement).
3. Figures re-synced (applied reco plots regenerated 2026-10-06), new figures where the procedure
   changed (origin test; mass view); `/compile-note` clean; `/review-note` PASS.

## Physics Procedure (authoritative = the CODE + the owning tracking docs)
- ε_reco: `pair_reco_eff_origin_test_and_views.md` §2 + D3/D4 + R1–R3; definition otherwise
  `pair_reco_eff_gap_acceptance.md` §2 (gap cuts on the reco leg only).
- Correction ordering / unfolding: user decisions 2026-10-06 (D-U1, D-U2 of the origin-test doc):
  ε_reco binned in TRUTH, migration by (multidimensional) Bayesian unfolding; bin-by-bin only if
  shown sufficient (`pair_bin_migration_bin_by_bin_feasibility.md`: migration criterion met,
  prior-dependence/closure untested, user decision pending); no same-sign reco efficiency; the
  template fit is trigger-corrected only, ε_reco + unfolding act on the extracted signal yield
  (`low_mass_dimuon_template_fit.md` §3e).
- Trigger: `pp24_trig_eff_hybrid_application.md` §2 (as in the predecessor doc R1).
- Negative constraints: never retype binnings; Pb+Pb reco efficiency stays a PHbox placeholder (G7);
  do not edit sibling tracking docs — record disagreements in §Inconsistencies.

## Implementation Plan
| # | Step | Review | Status |
|---|---|---|---|
| 1 | Ground truth: reco numbers (R1) + trigger audit (R2) | — | DONE |
| 2 | Rewrite reco section | /review-note | DONE |
| 3 | Trigger-section fixes from the audit | /review-note | DONE |
| 4 | Manifest + `/sync-note-figures`; `/compile-note` | — | DONE (SYNC CLEAN; 43 pp, 0 undefined, 0 overfull) |
| 5 | Commit (IntNotes submodule, then parent), INDEX | — | DONE (IntNotes 54b7680, bf6381a, e29da28; not pushed) |

## Progress Log
- 2026-10-07 Step 1 DONE: reco numbers (R1); trigger audit by subagent (R2, merged from
  `_sub_note_trig_audit_1.md`, deleted after merge). Figures: manifest +2 (origin test, mass view),
  `sync_note_figures.py` copied 6 stale + 2 new → SYNC CLEAN. Bib: `DAgostini:1995` added to
  `ANA-HION-2023-07-INT1.bib`. `placeholder.md` item 3 updated (all OS, 0.5841/0.6351).
- 2026-10-07 Steps 2–3 drafted: reco section rewritten (backup in session scratchpad
  `reconstruction_efficiency.bak_20261007.tex`); trigger section fixes from R2 (backup
  `trigger_efficiency.bak_20261007.tex`). Build: 42 pages, 0 undefined, 0 overfull.

- 2026-10-07 Step 2–3 `/review-note` loop, 5 iterations (log
  `.claude/logs/review-note-20261007-005435-reco-eff-all-os-note.md`), closed ESCALATED at the cap:
  zero CRITICAL in iterations 2–5; the single iteration-5 WARNING (scope of the "3.5–47 %" kinematic-
  bound deficit) was fixed after the last review with the reviewer's suggested wording, not
  re-reviewed. Substantive catches, all fixed: clamp statement (3 of 516 649 signal pairs above 150 GeV
  ARE pT-clamped); kinematic-bound ΔR cells explained (C1 → EXPECTED PHYSICS: bound
  m/√(pT,min(pT−pT,min)) = 0.64/0.47/0.24 at 9/12.8/36.7 GeV, near-threshold muon + window-top mass →
  reco migration out; lowest purity in the migration study); no close-by (small-ΔR) loss visible at
  ΔR<0.2 granularity; mass vs ΔR comparison on stated conventions (non-bound 0.8 %, bound-cut last
  intervals 3.5–47 %, partial 1–7 %); origin test: independence is cell-level, projections significant
  (104.7/71, 39.1/16, inclusive 3.7σ, partly composition); composition term Δf × 0.006; planned tense
  for the template-fit ordering + PHbox for the present OS/SS per-pair weighting; purity caveat
  (truth-seeded tree); forward-η mass slope qualified (9–18.2 GeV); trigger: single-value uncertainty
  6 % (1<|η^pair|<2) to 22 %.

## Results & Observations

### R1 — Reco numbers on the 2026-10-06 product (all OS, Tight unless stated)
Sources: `reco_eff/pair_reco_eff_pp24_full.root` + `plots/pp24_reco_effcy_plots/tight/applied/pair_reco_eff_values.txt`
(2026-10-06 11:27); RDF `histograms_pythia_fullsim_pp24_no_data_resonance_cuts_full.root` (11:22);
crossx `histograms_real_pairs_pp_2024_2mu4_nominal.root` (11:28). Scratch macros `rf2.C`, `ov.C`.
- inclusive 0.5841 T / 0.6351 M. vs ΔR 0.6419/0.5842/0.5119/0.2680; at 9–12.79 GeV 0.5409/0.5454/0.5097/0.2680;
  at 18.18–25.85 0.6581/0.6350. vs pT 0.5324/0.6100/0.6474/0.6703/0.6734/0.6572/0.6665/0.6391.
  vs η (−2.2→2.2) 0.6348/0.7895/0.4488/0.6560/0.4758/0.6557/0.4459/0.7871/0.6352 (±η within 0.003).
- 2D map: 1.0<|η|<1.5 rows 0.377–0.479; |η|<0.5 row 0.420 (lowest pT) → 0.59–0.63 above 25.85 GeV;
  other rows 0.574–0.746 at lowest pT, 0.75–0.85 at 25.85–52.23 GeV; outermost rows above 52 GeV
  0.50–0.83 (errors 0.05–0.35); top-pT cell of −2.2<η<−2.0: 0.600 ± 0.219; of 2.0<η<2.2: 0.500 ± 0.354.
- census: 144 of 288 cells no generated pair, 0 measured zeros (was 144 + 1 with single-b).
- top ΔR interval 0.069 % of the in-range denominator, only in the lowest pT bin; no ΔR overflow.
- truth denominator outside the measured range: 3.29 % of the truth signal region; reco-inside /
  truth-outside (η) migration-in 0.015 % (T) / 0.013 % (M) of the numerator.
- crossx application (pp24): reco-only factor Σw_reco·w_unfold/Σw_unfold = 1.6875 (mean ε 0.593); per
  populated (fine pT × 0.1 η) cell 1.171–2.654; per coarse pT 1.849 (9–12.8) then 1.50–1.62; trigger
  factor 1.880; trigger-weighted reco factor 1.718; total 3.230.
- Origin test (origin-test doc R1), mass diagnostic (R2), migration study (bin-migration doc R3) as cited.

### R2 — Trigger-section audit (merged from `_sub_note_trig_audit_1.md`)
- No trigger code/product change since 2026-09-18 (only b6ef0e9 09-17: poly cells OS-only, already
  described). All 7 trigger figures byte-identical to sources. Every trigger number re-extracted
  from the 2026-10-06 crossx and fit files: unchanged (1.8799; per-bin 2.357/1.735/2.839/3.771/4.07–4.59;
  cells 1.271–7.250; region-B 0.0310 %; ε_ΔR, single-value table, SF, seam numbers).
- Stale (fixed in the tex): ε_reco definition wording l.28–30, 37–42, 193–195 (reco now
  generated-pair based, gap cuts on the reco side); "single-b yield" l.588 → OS yield before
  background subtraction; SS weight statement l.582–584; "identical signal regions" l.716–718
  (R_AA mixed-fiducial); MC request "requested" l.613/631 → "in preparation" (no submission recorded
  in `additional_mc_request_jo_test.md`); 2.5× → loss ≈3× (0.780 vs 0.307); 23 % → 22 %; 1.74 at
  15–18 GeV; comment line pointers :386-411, :650-677.

## Inconsistencies (doc vs code)
1. `pair_trig_eff_pp24_full.root` provenance names a moved MC-fit path
   (`mc_trig_eff/singles_fits/single_mu_effcy_pT_fit_mc_pp24_full.root` now) and still says "NOT wired into the
   cross-section" (wired since 2026-09-17).
2. `pp24_trig_eff_hybrid_application.md` D3 quotes the 1–4 GeV wide/sig difference as 10–21 %; current file
   −2.6 % to +15.5 % (the note is right).
3. INDEX scope of `pbpb2026_analysis_support.md` still says step 14 RUNNING; the 2026-inclusive R_AA exists
   (2026-10-06 11:32).
4. SS signal-region yield (`RDFBasedHistFillingPP.cxx:897-913`) is weighted with w_reco from the OS product.
   Not a contradiction of D-U2 in substance: weighting OS and SS with the same per-cell ε is equivalent to
   correcting the OS−SS difference per cell; stated in the note as a provisional PHbox.

## Remaining Work

## Latest Stage
DONE 2026-10-07. Both sections updated, figures synced, note builds clean, IntNotes committed (not
pushed). Open for the user: optional re-review of the post-cap fix; §Inconsistencies items 1–3
(sibling-doc / provenance staleness, not edited); push of the IntNotes submodule.
