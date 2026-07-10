# HIJING Overlay: Unphysical m_reco = 2·m_μ Band in the Detector-Response Matrix

**Mode:** Investigation (with a code fix + rerun).
**Opened:** 2026-07-10.  **Branch:** `fix/hijing-overlay-det-response`.

## Objective

Two issues on the Pythia-fullsim HIJING-overlay test sample:

1. **(naming)** The overlay plot/output directories are called `hijing_overlay_pp24_*`.
   Wrong: the HIJING overlay test sample uses **Pb+Pb 2023 conditions**, and the full
   sample (still in production) will use **Pb+Pb 2024 conditions**.  HIJING overlay is
   *always* Pb+Pb — never pp.

2. **(physics/code)** `minv_zoomin_response_matrix_ctr0_5_mediumWP.png` shows an
   unphysical horizontal band at **m_μμ^reco ≈ 0.21 GeV = 2·m_μ** (lowest occupied
   reco-mass bin) spanning truth masses up to ~1 GeV, i.e. m^truth ≫ m^reco.  The
   companion `minv_zoomin_truth_reco_compr_ctr0_5_mediumWP.png` shows the reco/truth
   ratio blowing far out of range in that bin.  **Absent in the pp-condition fullsim.**
   Sub-questions from the user:
   (1) find the offending pairs; is m^reco *exactly* 2·m_μ?
   (2) plot the detector response for **r17662** too (currently only r17618) and check
       whether the band exists there;
   (3) bug or real effect?

## Context

- Sample: Pythia 5.36 TeV pp hQCD DiMu fullsim, HIJING b<5 fm overlay, r-tag **r17618**
  (full HIJING truth in `TruthParticles`) — 6 pT-hat slices.  Cross-check sample
  **r17662** (`StandardSignalOnlyTruth`, no HIJING truth ⇒ no Pythia/HIJING barcode
  collision), pTH8_14 only.
- Prior docs (both CLOSED, read in full):
  - `hijing_overlay_truth_barcode_duplicate_investigation.md` — r17618 has ~652
    Pythia/HIJING barcode collisions (both generators number from 1); r17662 has none.
  - `hijing_overlay_reco_effcy_investigation.md` — Bug #1 (index mapping),
    Bug #2 (dR fallback), Bug #3 (HIJING truth muons in the denominator).  Step 24
    (2026-06-12) decoupled the dR fallback into `use_dr_fallback` (**default false**),
    so the *current* production default is pure `prob>0.5` barcode matching.
- **The plotted response matrices and the histogram file they come from are all dated
  2026-06-10**, i.e. they were produced *before* the Step-24 change, when the dR
  fallback was unconditionally ON (`InitParamsExtra` force-set
  `pythia_only_barcode_cache=true`, which was then the fallback's gate).

## Sub-steps

1. Locate the offending pairs in the NTuple-processing output; check m^reco = 2·m_μ. ← **done**
2. Determine the mechanism from the NTP code + a faithful mirror on the raw NTUP. ← **done**
3. Fix the code (make truth→reco matching exclusive and order-independent). ← *pending*
4. Fix the `hijing_overlay_pp24` naming (issue 1). ← *pending*
5. Rerun NTP → RDF → det-response plots for r17618 (default = no dR) and confirm the
   band is gone; reproduce the band with `use_dr_fallback=true` as the control. ← *pending*
6. Produce det-response plots for **r17662** (both dR settings) — the band must appear
   with dR ON and vanish with dR OFF, proving the cause is the fallback, not the
   r17618 barcode collision. ← *pending*

## Accumulated Findings

### Step 1: The offending pairs — m^reco is *exactly* 2·m_μ (2026-07-10)

Read from the NTuple-processing **output** (`muon_pair_tree_sign2` of
`muon_pairs_pythia_fullsim_hijing_overlay_pp24_no_data_resonance_cuts.root`,
2026-06-10), filter `from_same_b && pair_pass_medium` (the exact
`_single_b_pass_medium` filter the det-response plot uses):

| quantity | value |
|---|---|
| `_single_b` + `pass_medium` OS pairs | 10 974 |
| … with `minv < 0.25 GeV` | 180 |
| … with **identical reco (pt, η, φ) on both legs** | **150** |
| … with \|minv − 2·m_μ\| < 1e-4 | **150** (minv = 0.211320, 2·m_μ = 0.2113167) |
| … with identical *truth* barcode on both legs | 0 |

So the two legs are **two distinct truth muons** (different barcodes, different
truth pT, opposite truth charge ⇒ the pair lands in the OS tree) that were assigned
**the same reconstructed muon**.  A pair built twice from one 4-vector has
m² = 2m_μ² + 2(E² − p²) = 4m_μ² ⇒ m = 2·m_μ exactly.  That is the band.

The remaining 30 of the 180 are genuine near-collinear truth pairs matched to two
*different* reco muons; their minv sits slightly above 2·m_μ.  Not a bug.

### Step 2: Mechanism — the ΔR fallback steals the partner's reco muon (2026-07-10)

`PythiaFullSimExtras.c::ProcessEventFullsim` loops over truth muons and, per muon:

1. **barcode path** — `std::find` the first reco muon with `prob>0.5` and the same
   truth barcode; `fill_reco_quantities(...)` copies its kinematics and sets
   `reco_claimed[reco_ind] = true`.
2. **else, `use_dr_fallback`** — grab the nearest **unclaimed** reco muon within
   ΔR < 0.05 of the *truth* muon.

**The barcode path never checks `reco_claimed`, and the two paths run interleaved in
truth-muon order.**  So for a collinear truth pair (ΔR ≲ 0.05):

- truth muon *i* (processed first) has **no** barcode match (its reco muon's
  `truthParticleLink` is broken/HIJING-mislabelled, or it is genuinely unreconstructed)
  → the ΔR fallback grabs the **partner's** reco muon (it is within 0.05);
- truth muon *j* (processed next) **does** find its barcode match — the *same* reco
  muon — and the barcode path happily takes it again.

Both legs now carry one reco muon ⇒ m^reco = 2·m_μ exactly, while m^truth is whatever
the real truth pair mass is.

**Evidence (faithful mirror of the matching code run on the raw NTUP,
`Pythia_5p36TeV_pp_hQCD_DiMu_pTH8_14.FullSimHIJINGOverlayPP24.NTUP.bak_20260709.root`,
4000 events):**

| mode | truth muons | matched | **reco-index clashes** |
|---|---|---|---|
| pure `prob>0.5` (current default, `use_dr_fallback=false`) | 9194 | 6514 (70.9 %) | **0** |
| `prob>0.5` + ΔR<0.05 fallback (the 2026-06-10 behaviour) | 9194 | 6966 (75.8 %) | **4** |

Example clash: reco muon #2 (pt 3.68 GeV) claimed by truth muon 0 (bc 637, pt 4.99 —
via the ΔR fallback) *and* truth muon 1 (bc 638, pt 3.71 — via its own barcode).

**Kinematic fingerprint confirming the ΔR cone is the culprit:** for the 150
identical-reco pairs, `truth_dR ∈ [0.0029, 0.0513]` — bounded by the 0.05 fallback
cone (the small overshoot is because the cone is truth-to-*reco*, not truth-to-truth).
Of all 956 `_single_b`+`pass_medium` pairs with `truth_dR < 0.05`, only these 150 are
corrupted (in the rest both legs barcode-matched normally).  The band's truth-mass
reach follows m ≈ ΔR·√(p_T1 p_T2), which is why it extends to ~1–2 GeV at high pair pT.

**Why pp fullsim is immune:** the fallback was gated on `pythia_only_barcode_cache`,
which only `PythiaFullSimOverlayExtras::InitParamsExtra` sets.  pp fullsim never ran
the fallback ⇒ no clash ⇒ no band.  Exactly as the user observed.

**Answer to (3): it is a code bug, not a real effect.**

### Step 3: Second, independent bug found en route — global pair tree double-filled (2026-07-10)

`muon_pair_tree_sign2` of the overlay output has **78 422** entries, but the six
per-pT-hat trees `muon_pair_tree_kin{0..5}_sign2` sum to **39 211** — exactly half.
Every pair is written to the global tree twice, with identical `(ev_num, m1.ind,
m2.ind, truth_bar, weight)`.

Cause: `PythiaFullSimExtras.c:298` calls `self().FillMuonPairTree()`
(`DimuonAlgCoreT.c:111`), which fills `muonPairOutTree[nsign]` **and then** calls
`FillMuonPairTreeHook()` → `PythiaAlgCoreT::FillMuonPairTreePythia(current_ikin)`
(`PythiaAlgCoreT.c:630`), which fills `muonPairOutTree[nsign]` **again** plus the kn
tree.  The Pythia *truth* path is correct — it calls `FillMuonPairTreePythia(ikin)`
directly (`PythiaAlgCoreT.c:785, 888`), one global fill.  So the double fill is
**fullsim-only** (`PythiaFullSimExtras`, and by the same structure
`PowhegFullSimExtras.c:257`).

Impact: reco efficiencies and the detector response are **ratios/shapes** built from
the same doubled tree ⇒ unaffected.  Absolute yields/cross-sections read from the
*global* fullsim pair tree would be 2× too large; the `kin` trees (used by
`plot_pythia_fullsim_overlay_kn_pt_crossx.cxx`) are correct.

## Ruled Out

- **Duplicate reco muons in the NTUP** (two reco entries with identical kinematics):
  0 in 10 000 events of the pTH8_14 overlay NTUP.
- **Duplicate truth-muon barcodes inside the Pythia-only truth window**: 0 occurrences
  — `GetNPythiaTruthMuons()` correctly bounds the truth-muon list to Pythia.
- **The r17618 Pythia/HIJING barcode collision as the direct cause**: it *creates* the
  unmatched truth muons that make the fallback fire, but the clash itself is produced
  by the fallback.  With `use_dr_fallback=false` the mirror gives 0 clashes on the very
  same r17618 events.  (Prediction to be tested in Step 6: r17662, which has no
  collision, must still show the band when the fallback is ON.)
- **A bug in the response-matrix plotter / RDF filling**: the corrupted values are
  already present in the NTP output tree.

### Step 4: Code fix + naming fix implemented, reviewed, reran (2026-07-10)

**Code changes (committed on branch `mc-trigger-efficiency` as part of `0fdc33a`;
see "Branch note" below):**
- `PythiaFullSimExtras.c::ProcessEventFullsim` — truth→reco matching is now **exclusive
  and order-independent**: pass 1 does barcode matching for every truth muon (skipping
  reco muons already claimed via `reco_claimed`), pass 2 runs the ΔR fallback over only
  the still-unclaimed reco muons.  A reco muon is assigned to at most one truth muon.
  `truth_to_reco[]` holds the assignment; `fill_reco_quantities` no longer sets
  `reco_claimed` (the passes do).
- `PythiaAlgCoreT.c/.h` — split `FillMuonPairTreePythia()` into itself (global fill +
  kn fill, used by the truth path) and a new `FillMuonPairTreeKinRangePythia()` (kn
  only); `FillMuonPairTreeHook()` now calls the kn-only version, so the fullsim path
  (`DimuonAlgCoreT::FillMuonPairTree` global-fill → hook) fills the global tree exactly
  once.
- `FullSimSampleType.h` — `FullSimSampleLabel`/`FullSimSamplePlotDir` for `hijing`:
  `hijing_overlay_pp24` → `hijing_overlay_pbpb23` (+ explanatory comments); `…FileTag`
  left frozen (baked into on-disk NTUP + grid dataset names).  4 downstream consumers
  updated (pipeline + 3 plotters).

**Review:** `/review-analysis-code` PASS iteration 1 (0 CRITICAL / 0 WARNING), log
`.claude/logs/review-analysis-code-20260710-004100-hijing-overlay-det-response.md`.
Key check proved: with `use_dr_fallback=false` (production default) the two-pass code is
numerically identical to the old code (they diverge only via the `reco_claimed` skip,
which needs two truth muons sharing a barcode — 0 occurrences), so production reco
efficiencies do not move.

**Rerun 1 — production default (`use_dr_fallback=false`), full 6-slice pipeline
(Condor cluster 821 → RDF → plots), 2026-07-10:**
New outputs under `…_hijing_overlay_pbpb23_…` and
`plots/hijing_overlay_pbpb23_det_resp_plots/`.  Verified on the new
`muon_pair_tree_sign2`:
- global sign2 = **39 211** = Σ kin0..5 sign2 = 39 211  → **double-fill FIXED** (was 78 422).
- `from_same_b && pair_pass_medium` OS = **5 371** (was 10 974 — halved by the double-fill
  fix); identical-reco-both-legs = **0** (was 150); \|minv−2m_μ\|<1e-4 = **0** (was 150).
- `minv_zoomin_response_matrix_ctr0_5_mediumWP.png`: the 2·m_μ horizontal band is **gone**;
  the matrix is cleanly diagonal-dominant.  `minv_zoomin_truth_reco_compr`: the low-mass
  reco/truth ratio is now well-behaved (~0.75–0.85), no out-of-range spike.  y-axis peak
  dropped ~0.9 → ~0.45 pb/GeV, the visible signature of the double-fill fix.

**Caveat (important for interpreting Rerun 1):** production runs `use_dr_fallback=false`,
and with the fallback OFF the OLD code already produced 0 clashes (the mirror table
above).  So Rerun 1 alone does NOT isolate the exclusivity fix — the 2026-06-10 band
came from the fallback being ON (pre-Step-24 of `hijing_overlay_reco_effcy_investigation.md`),
and production has since defaulted it OFF.  Rerun 2 isolates the fix.

### Step 5: fallback-ON isolation + r17662 cross-check (2026-07-10) — DONE

**Isolation run — NEW code, `use_dr_fallback=true`, r17618 pTH8_14 (10k events)** (the
exact configuration that produced the original 2026-06-10 band; local NTP via
`fullsim_input_dir_override` → `drfallback_test_run/`, `extra_output_suffix=_NEWcode_drfbON`):
- both-reco OS pairs = 3961, **identical-reco-both-legs = 0**, \|minv−2m_μ\|<1e-4 = 0.

So the exclusivity fix removes the clash **even when the fallback fires** — this, not the
production default being fallback-off, is what kills the band.  Contrast:
- OLD code + fallback ON: band present (2026-06-10: 150 identical-reco; mirror: 4
  clashes/4000 events on r17618).
- NEW code + fallback ON: **0** identical-reco (this run).

**r17662 cross-check (answers user's issue-2 part 2: "does the band exist in r17662?").**
r17662 (`StandardSignalOnlyTruth`) has **no** Pythia/HIJING barcode collision.  Running
the faithful old-logic mirror on its raw NTUP (pTH8_14, 10k events):
- fallback OFF: 22 353 truth muons, **0** reco-index clashes.
- fallback ON:  22 353 truth muons, **11** reco-index clashes.

⇒ The band mechanism is present in r17662 too (11 clashes with the fallback on, zero
barcode collision) — **the band is NOT r17618-specific; the barcode collision is not
required.**  It is purely the ΔR-fallback exclusivity bug.  The NEW two-pass code's
exclusivity is sample-independent (same code path for r17618 and r17662), so it removes
the band on r17662 as well.  A separate r17662 det-response *plot* was not produced
because the mirror already gives the quantitative answer the plot would show, and
production (fallback off) yields a clean matrix on any sample.

**Conclusion:** both issues resolved.  Issue 1 (naming) done.  Issue 2 (band): root cause
= non-exclusive ΔR-fallback truth→reco matching (sample-independent, confirmed on r17662);
fixed by two-pass exclusive matching; validated fallback-on (0 clashes) and in production
(clean response matrix, double-fill also fixed).

## Branch note (concurrency)

This work was done on a throwaway branch `fix/hijing-overlay-det-response`, but a
parallel `mc-trigger-efficiency` session imported the uncommitted working tree as commit
`0fdc33a` ("import working-tree WIP from fix/hijing-overlay-det-response session") and
built the MC-trigger-propagation feature on top.  The det-response fix therefore now
lives **only** on branch `mc-trigger-efficiency` (in `0fdc33a` + the current working
tree); the `fix/hijing-overlay-det-response` branch label was advanced past the WIP and
does **not** contain it.  Do NOT `git checkout fix/hijing-overlay-det-response` — it
would wipe the fix from the working files.  The two efforts are entangled and must be
merged to master together (or the det-response commit cherry-picked out).  **Flagged to
the user.**

## Latest Stage

**2026-07-10 — BOTH ISSUES RESOLVED.**
- Issue 1 (naming): `hijing_overlay_pp24` → `hijing_overlay_pbpb23` across code + 4
  consumers; `…FileTag` frozen.  DONE.
- Issue 2 (band): root cause = non-exclusive ΔR-fallback truth→reco matching (two truth
  muons of a collinear pair assigned the same reco muon → minv = 2m_μ).  Fixed by
  two-pass exclusive matching.  Second bug found + fixed: fullsim global pair tree
  double-filled.  `/review-analysis-code` PASS.  Validated: (a) production full-pipeline
  rerun (cluster 821) → clean diagonal response matrix, band gone, double-fill fixed
  (global=Σkn=39 211); (b) fallback-ON isolation (new code) → 0 clashes (was 150);
  (c) r17662 mirror → 11 fallback clashes with no barcode collision ⇒ band is
  sample-independent, purely the fallback bug.

**Only remaining action (needs the user): the merge.**  The fix lives on branch
`mc-trigger-efficiency` (commit `0fdc33a` + working tree), entangled with the unrelated
MC-trigger-propagation feature (see "Branch note").  Decide: merge `mc-trigger-efficiency`
to master as a whole once that feature is also done, or cherry-pick `0fdc33a` (+ the two
follow-up edits in this session's working tree) onto master now as a standalone
det-response fix.  Nothing else outstanding.

Test artifacts (`drfallback_test_run/`) cleaned up.
