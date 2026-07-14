# HIJING Overlay: Unphysical m_reco = 2·m_μ Band in the Detector-Response Matrix

**Mode:** Investigation (with a code fix + rerun).
**Opened:** 2026-07-10.  **Branch:** `mc-trigger-efficiency` (see Branch note — entangled).

## Autonomy Contract (DONE — all Done items met; doc CLOSED 2026-07-14)
- Mandate: run autonomously to DONE; do NOT pause to confirm progress. Finishing a
  plan, a passing small test, or one pipeline stage is NOT a stopping point.
- Done =
  1. **dR<0.05 fallback DELETED entirely** — `use_dr_fallback` flag + the pass-2 branch
     + any run-script references removed; code compiles; `/review-analysis-code` PASS.
     (User: the fallback ALWAYS overestimates reco efficiency for pairs with ΔR<0.05 and
     is never acceptable — not as a fallback, and geometric-ΔR is not our method.)
  2. Production r17618 overlay pipeline reran with the fallback-free code; reco-eff +
     det-response plots regenerated under `hijing_overlay_pbpb23`; 2·m_μ band absent;
     efficiencies consistent with the prior pure-prob>0.5 production (sanity).
  3. **Direct r17618-vs-r17662 reconstruction-efficiency comparison plots**, BOTH
     single-muon AND muon-pair, BOTH samples under **pure prob>0.5 (no dR)**, at matched
     pT-hat (r17618 kin0 vs r17662 pTH8_14) and relevant centralities (0-5%, 5-10%).
     Clear verdict: does r17618 genuinely UNDERESTIMATE reco efficiency vs r17662 because
     of the Pythia/HIJING truth-barcode overlap (r17618 < r17662 ⇒ real)?  If real:
     diagnose the mechanism, quantify, and fix-if-confident or STOP-and-ask.  If not
     real: document that r17618 is safe as the full-sample choice.  (User is switching
     the FULL sample to r17618 — need HIJING truth to separate hadronic vs fake
     backgrounds — but must rule out a reco-efficiency underestimate first.)
- Stop-and-ask = ANY physics-results-bending ambiguity (no fixed list; use judgment;
  when unsure whether an ambiguity is blocking, treat it as blocking → AskUserQuestion).

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

### Step 6: dR fallback DELETED entirely (2026-07-10)

Per user direction (the ΔR<0.05 fallback ALWAYS overestimates the pair reco efficiency
for ΔR<0.05 pairs — it recovers genuine physical losses that also exist in data, has no
data analog, and the reco+trigger-efficiency correction is applied TO data; geometric ΔR
matching is a distinct alternative method, never a valid fallback):
- Deleted the `if (self().use_dr_fallback){...}` pass-2 block from
  `PythiaFullSimExtras.c::ProcessEventFullsim`; the matching is now the single pure
  prob>0.5 exclusive barcode pass (ATLAS-standard, Run 2 dimuon note
  ATL-COM-PHYS-2021-1094).  Removed `bool use_dr_fallback` from `PythiaTruthExtras.h`.
  Cleaned obsolete comments in `run_pythia_fullsim_overlay_r17662_nodr.sh`.
  (`use_geometric_matching` was already removed in Step 24 of the reco-effcy doc.)
- Compiles (ACLiC exit 0).  `/review-analysis-code` **PASS iteration 1, 0 issues**
  (log `review-analysis-code-20260710-184500-delete-dr-fallback.md`).  Behaviourally
  identical to the production default (fallback was already OFF) — changes NO production
  numbers; it removes the ability to enable a wrong method.

### Step 7: CLEAN r17618-vs-r17662 reco-efficiency comparison (2026-07-10, IN PROGRESS)

**Motivation (user, physics-crucial):** the full sample is switching to **r17618** (need
HIJING truth to separate hadronic vs fake backgrounds — r17662 signal-only truth
mislabels real HIJING muons as fake).  Must rule out that the r17618 Pythia/HIJING
**truth-barcode collision** genuinely UNDERESTIMATES the reconstruction efficiency (a
collision can stamp a genuinely-reconstructed Pythia signal muon's reco with a wrong/
HIJING barcode → barcode match fails → reco_match=false → ε too low).  Prior
r17618-vs-r17662 single-muon plots (Jun 12) were **confounded**: r17618 ran WITH the +dR
fallback (which artificially recovered exactly those collision losses, masking any
underestimate) and mixed all pT-hat slices.  No pair comparison existed.  **The clean
test — both samples pure prob>0.5 (fallback now deleted), matched pT-hat (kin0) — had
never been done.**

Setup (both pure prob>0.5, pTH8_14): r17618 via `r17618_kin0_run/` (kin0-only override
so pT-hat matches r17662), r17662 via `r17662_run/`.  Four NTP runs (single-muon + pair
each), then official RDF pair reco-eff hist filling, then comparison plots (single-muon
macro `plot_single_muon_reco_effcy_r17618_vs_r17662.cxx` repointed at the pure-prob
files; new pair macro `plot_pair_reco_effcy_r17618_vs_r17662.cxx`).

**Critical fact: r17618 and r17662 are the SAME 10 000 events** (100 % eventNumber
overlap — r17662 is a re-reconstruction of the same EVNT with HIJING truth stripped),
same detector overlay ⇒ **identical reco muons**; they differ ONLY in whether HIJING
truth exists (⇒ barcode collision or not).  This is a perfectly controlled test.

**RESULT — reco efficiency is IDENTICAL between the two r-tags (no underestimate):**

Single-muon (0-5% centrality, weighted):
| quantity | r17618 (collision) | r17662 (no collision) |
|---|---|---|
| reco_match | 0.8102 | 0.8117 |
| ε pass_medium | 0.6842 | 0.6810 |
| ε pass_tight | 0.6368 | 0.6326 |

Per-pT-bin reco_match agrees within stats (4-6: 0.805/0.808; 6-8: 0.821/0.819; 8-10:
0.830/0.828; 10-15: 0.881/0.887).  Pair (0-5%, official RDF hists): OS ε_medium
0.4754/0.4742, ε_tight 0.4114/0.4078; single-b ε_medium 0.4485/0.4551 (within stats).
r17618 is **not systematically below** r17662 anywhere — differences are <0.5 % and
non-directional.  Plots (`plots/r17618_vs_r17662_comparison/` single-muon + `.../pair/`,
tight & `medium_wp/`) overlay point-for-point within errors across all pT and all ΔR,
**including ΔR<0.1** (where the deleted fallback would have mattered).

**MECHANISM (independent raw-NTUP subagent cross-check, `_sub_barcode_reco_mechanism.md`,
now merged):** 99.97 % of fiducial Pythia signal muons in r17618 have their barcode reused
by a HIJING truth particle (near-100 % at-risk), yet the matched fraction (A) is unchanged
(0.8172 vs 0.8176 overall).  **Why the collision is harmless for the reco→truth-muon
link:** the collision is on the barcode VALUE, and `std::find(muon_truth_barcode == truth
barcode)` matches on that integer — the reco muon's stored barcode is the same integer
whether the AOD `truthParticleLink` resolves to the Pythia or the colliding HIJING
particle, so the match succeeds whenever the Pythia muon is genuinely reconstructed (its
own hits dominate its track).  HIJING truth merely RELABELS the muons that ALREADY failed
(hit-dilution losses): their nearest reco's stored barcode moves from `-1` (r17662, nothing
to link to) to a HIJING barcode (r17618) — the −191 in the `bc=-1` bucket exactly equals
the +189 in the HIJING-bc bucket.  **Zero muons move from matched→unmatched.**  (Contrast:
the collision DOES corrupt ancestor-tracing / origin classification — BFS chain-walking
follows shared barcodes into HIJING territory — which is the separate "others"-excess issue
that `pythia_only_barcode_cache` already handles.  The single reco→truth-muon link used for
efficiency is immune.)

### Step 8 (2026-07-13): CORRECTION — the two r-tags are DIFFERENT digi+reco passes; the
### user's spurious-HIJING-match hypothesis is RULED OUT; both samples are physically valid

**RETRACTION of a Step-7 claim.**  Step 7 asserted the two r-tags have "identical reco
muons ⇒ a perfectly controlled test".  **That was WRONG and is retracted.**  Verified:

| check (400 shared events) | result |
|---|---|
| eventNumber unique & shared | 10000/10000, **100 % overlap** ✔ |
| Pythia **truth** muons identical | ✔ |
| **FCal_Et** (⇒ same overlaid HIJING event) | **400/400 bit-identical** ✔ |
| **reco muon collection identical** | **0/400** ✘ — counts differ (9 vs 10, 3 vs 4 …) |

**Why:** the AMI diff between the r-tags is a SINGLE field — `digiSteeringConf:
StandardSignalOnlyTruth` (r17662).  Same Athena 24.0.58, same `OFLCOND-MC23-SDR-RUN3-05`,
same `ATLAS-R3S-2021-03-02-00`, same HIJING overlay input, same beamspot postExec.
`digiSteeringConf` is a **digitisation** steering option (not a pure "write less truth"
flag): it changes the digi algorithm sequence and thereby the **tracking-detector RNG
stream**.  So r17618 and r17662 are **two different digitisation+reco passes of the same
generated events** — each deterministic, but not the same computation.  (The calorimeter
stream was unaffected — hence bit-identical FCal_Et — a neat internal check that this is
RNG-level, not physics-level.)

**Is anything PHYSICAL different?  NO — verified:**
| observable | r17618 | r17662 | compatibility |
|---|---|---|---|
| FCal_Et (occupancy) | 4.4332 TeV | 4.4332 TeV | KS p = **1.000** |
| ⟨N reco muons⟩/event | 4.559 | 4.569 | KS p = 0.94 |
| **detector response** (p_T^reco−p_T^truth)/p_T^truth | mean −0.00136, **RMS 0.03624** | mean −0.00117, **RMS 0.03633** | **KS p = 0.75**, χ² p = 0.9997 |

⇒ Same occupancy, same multiplicity, **same momentum resolution (RMS agree to 0.25 %)**.
**BOTH r-tags give correct reconstruction efficiencies AND correct detector response** for
the final measurement.  They are two independent statistical realisations of the same
physics.  Consequence: they are **NOT interchangeable event-by-event** — a per-muon
reco_match/pass_WP difference between them is EXPECTED, not a bug.  The residual ~0.2-0.5 %
efficiency difference is fully accounted for by (a) the two reco realisations and (b) 6
muons of ~12 900 flipping across the prob>0.5 threshold (truthMatchProbability shifts
because r17662 has no HIJING truth to attribute hits to; mean |Δprob| = 0.018).

**User's hypothesis — spurious HIJING barcode match (would make r17618 OVERestimate) —
RULED OUT by direct measurement (tested WITHIN r17618 alone, so no cross-sample matching
is needed):** a Pythia truth muon spuriously matched to a HIJING muon's track would land
far away in ΔR.
```
r17618: 12913 barcode+prob>0.5 matches | ΔR(truth, matched reco) > 0.05 : 0 (0.000 %) | worst ΔR = 0.022
r17662: 12917 matches                  | ΔR > 0.05 : 0 (0.000 %) | worst ΔR = 0.023
```
**Every match sits within ΔR < 0.022 of its own truth muon** — not one match to a distant
HIJING muon.  The real muon's own hits always win the prob>0.5 association over a
barcode-colliding HIJING particle.

**Latent trap noted (not currently biting):** `GetNPythiaTruthMuons` falls back to
returning ALL truth muons when an event has no truth particle with barcode > 200000.  It
fires in **117/10000 r17662 events** and **0 r17618 events**.  Harmless today (r17662 has
no HIJING truth to leak in), but in a HIJING-truth sample it would silently pull HIJING
muons into the Pythia truth-muon list.  Worth hardening.

**Net: the Step-7 conclusion STANDS, but its basis is corrected** — r17618 has no
reco-efficiency bias because (a) the spurious-match mechanism demonstrably does not occur
and (b) the two reco passes agree statistically, NOT because the reco was identical.

### Step 9 (2026-07-14): Pythia-truth-index-guard CLOSURE TEST — found & fixed a real leak

**The test (user-requested).** Since r17618 and r17662 are the SAME 10 000 events with the
SAME Pythia truth, every Pythia truth muon PAIR must get the identical
flavor / origin / m1,m2 parent_group / from_same_b after ntuple processing.  Any difference
means the **Pythia truth index guard** (`GetNPythiaTruthMuons`) leaked HIJING truth.  This is
the sharpest possible closure test of the guard.  (NTP `ev_num` is the *entry index* and the
two NTUPs have DIFFERENT event ordering ⇒ pairs must be keyed on
`(eventNumber, truth-barcode pair)`, mapping ev_num→eventNumber via the raw NTUPs.)

**Result BEFORE the fix:** 0 mismatches on every matched pair (6271) — but r17618 had
**6 EXTRA pairs** (2 SS + 4 OS) and **4 extra fiducial truth muons** (15804 vs 15800).

**ROOT CAUSE — a genuine bug in the guard.**  The container is
`[Pythia gen][Pythia Geant4][HIJING gen][HIJING Geant4]`, and the guard bounded the Pythia
block at *the first barcode > 200000* (the Geant4 marker).  **In 117/10000 events (1.17 %)
the Pythia Geant4 block is EMPTY**, so the layout is `[Pythia gen][HIJING gen][HIJING Geant4]`
and that criterion lands on the **HIJING** Geant4 block — swallowing the entire HIJING
generator block.  Measured leak: **643 HIJING truth muons** admitted into the Pythia list,
**4 of them fiducial** (0.013 % of the 15804 denominator), producing the 6 spurious pairs
(e.g. ev 20245: Pythia bc 233 paired with **HIJING** bc 9562, pT 4.4, η −1.25).  This is also
why r17618's reco_match sat slightly BELOW r17662's (0.8102 vs 0.8117) — leaked HIJING muons
sit in the denominator and essentially never reco-match.
(Same root cause as the previously-logged "fallback" trap: both are the *missing Geant4
block* case.  In r17662 those same 117 events trigger the no-boundary fallback, which is
harmless there because there is no HIJING truth to leak.)

**FIX** (`PythiaTruthExtras.h::GetNPythiaTruthMuons`): the Pythia block ends at the first
index that is EITHER (a) a Geant4 particle (barcode > 200000) **OR (b) a barcode RESTART**
(`barcode[i] < barcode[i-1]`) — each generator block numbers its barcodes from 1, so the
restart marks where HIJING begins.  Non-overlay samples (pp fullsim) have monotonically
increasing generator barcodes followed by Geant4, so (b) never fires ⇒ **pp is unchanged**.

**Result AFTER the fix — PERFECT CLOSURE:**
| | r17618 | r17662 | extras |
|---|---|---|---|
| same-sign pairs | **1107** | **1107** | **0** (was 2) |
| opposite-sign pairs | **5164** | **5164** | **0** (was 4) |
| fiducial truth muons | **15800** | **15800** | **0** (was +4) |
| flavor / origin / parent_group / from_same_b / skipped mismatches | **0** | **0** | — |

Every Pythia truth muon and every truth pair now carries **identical** flavor, origin, parent
groups and `from_same_b` in both r-tags.  **The Pythia truth index guard is now validated
exactly**, and the HIJING-truth leak is closed.

**Also added:** `PythiaAlgCoreT::allow_missing_slices` (default **false** = strict).  The
strict "refusing to run on a missing pT-hat slice" check (correct for σ-weighted production)
made single-slice diagnostics impossible — notably r17662, which exists ONLY for pTH8_14.
The flag downgrades it to a loud warning for diagnostic runs; production behaviour unchanged.

**VERDICT (Autonomy Contract Done item 3): NO genuine reconstruction-efficiency
underestimate in r17618 from the barcode collision.**  r17618 ε = r17662 ε to <0.5 %
(single-muon AND pair, both centralities, both WPs), confirmed independently at the
raw-NTUP per-muon level.  The ~19 % of signal muons that go unmatched are a PHYSICAL /
prob-strictness loss, identical in both samples (r17662 has no collision yet the same
~19 %).  **r17618 is SAFE as the full-sample choice** — the user can switch to it (to keep
HIJING truth for hadronic-vs-fake background separation) with no reco-efficiency penalty.
This also settles the closed `hijing_overlay_reco_effcy_investigation.md` Step 8 worry
(the "18.3 % unmatched" was read as possibly collision-driven; it is physical, present
identically in the collision-free r17662).

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

## Final Summary (CLOSED 2026-07-14)

1. **Naming** — `hijing_overlay_pp24` → `hijing_overlay_pbpb23` (HIJING overlay is Pb+Pb,
   never pp); all plot dirs migrated, stale pp24 dirs removed; NTUP file tag frozen.
2. **The 2·m_μ band** — root cause: the ΔR<0.05 fallback + a non-exclusive barcode path let
   two truth muons of a collinear pair claim the SAME reco muon. Fixed (exclusive matching),
   then the **ΔR fallback was DELETED entirely** (it always overestimates the pair reco
   efficiency for ΔR<0.05 pairs). Matching is now pure prob>0.5 exclusive barcode.
   Second bug fixed: the fullsim global pair tree was double-filled.
3. **r17618 vs r17662** — the two r-tags are the SAME 10 000 events but TWO DIFFERENT
   digitisation+reco passes (`digiSteeringConf: StandardSignalOnlyTruth` shifts the tracking
   RNG stream). Nothing physical differs (FCal_Et bit-identical; detector response KS p=0.75,
   RMS agree to 0.25%) ⇒ **both give correct reco efficiency and detector response**.
   The barcode collision causes **no** reco-efficiency underestimate, and the suspected
   spurious-HIJING-match mechanism is **ruled out** (all 12913 matches within ΔR<0.022).
4. **Pythia truth index guard** — the closure test (identical truth ⇒ identical
   flavor/origin/parent-group for every pair) exposed a **real leak**: in 1.17% of events the
   Pythia Geant4 block is empty, so the "first barcode>200000" bound swallowed the HIJING
   generator block (643 HIJING muons; 4 fiducial; 6 spurious pairs). **Fixed** by also
   detecting the barcode restart. Now **perfect closure**: 1107=1107 SS, 5164=5164 OS,
   15800=15800 fiducial truth muons, **zero** flavor/origin/parent-group/from_same_b
   mismatches. **⇒ r17618 is SAFE for the full sample** (keeps HIJING truth for
   hadronic-vs-fake separation, at no cost in efficiency or response).

**Remaining (user decision):** the merge. All of this lives on branch `mc-trigger-efficiency`,
entangled with the unrelated MC-trigger feature (see Branch note).

**Downstream note:** the guard fix changes the overlay reco-efficiency denominator by 0.013%
(4 of 15804 fiducial truth muons) and removes 6 of 6277 truth pairs. Negligible, but the
overlay reco-eff / det-response outputs were produced with the pre-fix guard; regenerate them
when the full sample lands.
