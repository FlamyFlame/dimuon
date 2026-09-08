# pp24: muon pairs from ALL vertices (not only the primary), + gap-cut rerun

Mode: **IMPLEMENTATION**. Opened 2026-09-07.

## Objective

Change the **pp data** NTuple-processing stage so that a muon pair is kept if its two good
muons pass the impact-parameter cuts with respect to **the same vertex**, which may be a
**secondary (pile-up) vertex** — not only the primary vertex. Mirror the same procedure in the
**pp-conditions fullsim MC** that supplies epsilon_reco. Record the resulting secondary-vertex
statistics. Then rerun the whole pp24 chain, which also picks up the 2026-09-07 gap-cut change
(`muon_gap_cuts_acceptance.md` F17/F18).

## Autonomy Contract (ACTIVE — re-read on every compaction)
- Mandate: run autonomously to DONE; do NOT pause to confirm progress. Finishing a
  plan, a passing small test, or one pipeline stage is NOT a stopping point.
- Done =
  1. pp data NTuple stage keeps every pair of two good muons that pass |d0| and
     |z0 sin(theta)| w.r.t. the SAME reconstructed vertex (any vertex, not only index 0);
     the matched vertex index is stored on the pair. Pb+Pb NTuple output byte-identical.
  2. The same all-vertex procedure mirrored in the pp-conditions fullsim MC
     (`PythiaFullSimExtras` pp path + `PowhegFullSimExtras`), so epsilon_reco is measured on
     the same selection as the data. HIJING overlay untouched (provable no-op).
  3. Statistics record: fraction of good muon pairs whose matched vertex is NOT the primary,
     for pp24 in total AND in each of the 8 canonical `ParamsSet::pair_pt_coarse_bins`
     (user decision 2026-09-07: the canonical coarse axis, not a 16-bin axis — none exists).
     Written as a CSV + a table in this doc. The MC counterpart comes for free from the
     end-of-job counters added to the fullsim Extras and is recorded alongside.
  4. pp24 NTuple Condor rerun for BOTH the nominal (2mu4) and the trigger-efficiency
     (trigger_mode=1) modes, then the data mu4 turn-on refit on the NEW trees.
  5. pp24 Pythia-fullsim pipeline rerun (reco-eff + MC trig-eff), dR-correction fits rerun
     (including today's MINUIT-limits change, which goes through /review-analysis-code),
     `pair_reco_eff_pp24_full.root` rebuilt, pp24 crossx pipeline rerun to completion.
  6. Crossx Stage 7 sanity plot produced pp-only this round (user decision; Pb+Pb deferred).
  7. Docs updated (this doc, `signal_selection_change_impact.md`, `analysis_overview.md`
     if the signal region text moves, `muon_gap_cuts_acceptance.md` cross-ref, INDEX.md) and
     the work committed.
- Stop-and-ask = ANY physics-results-bending ambiguity (no fixed list; use judgment;
  when unsure whether an ambiguity is blocking, treat it as blocking → AskUserQuestion).

## Physics Procedure

### 1. Motivation

For **Pb+Pb** the mean number of collisions per bunch crossing is O(10^-3) and the analysis
applies pile-up **rejection** cuts, so every surviving event has one collision and the primary
vertex is the collision. For **pp24** the mean number of collisions per crossing is **~4**, and
— decisively — the **pp luminosity is defined for ALL inelastic collisions in the crossing, not
only the one that produced the primary vertex**. A cross-section

    dsigma/dX = (1/L_int) * SUM_pairs w(pair)

whose numerator counts only primary-vertex pairs while its denominator L_int counts all
collisions is therefore **inconsistent, and biased low**. The fix is on the numerator side:
count muon pairs from every vertex.

This is NOT a diagnostic or a systematic — the all-vertex procedure is the CORRECT one for pp.
The secondary-vertex fraction recorded in §3 is a statistics record, not a cross-check.

### 2. Top-level equation

Unchanged in form (see `pp24_crossx_rerun_2026_08.md` §2). What changes is the SUM's domain:

    pair kept  <=>  EXISTS a reconstructed vertex v (track-bearing, vtx_ntrk >= 2) such that
                    BOTH muons of the pair pass |d0| < d0cut AND |z0_v sin(theta)| < z0cut

where `z0_v` is the muon's z0 measured with respect to **that** vertex v. The two muons must
use the SAME v. Every other cut (quality/WP, pT > 4, |eta| < 2.4, one-sided Delta-p/p,
trigger match, resonance veto) is unchanged and still applied per muon / per pair.

### 3. Step-by-step method

**(a) Per-vertex z0.** The skim (`SkimCode/.../TrigRates.cxx`) stores, per muon,
`z0 = idTrk->z0() + idTrk->vz() - z_vtx(0)`, i.e. ALREADY relative to vertex index 0
(`ProcessMuons()`, `priVtx = *(VertexContainer->cbegin())`), and stores the full vertex list in
`vtx_x/vtx_y/vtx_z/vtx_ntrk` **in the same container order** (`ProcessVertex()`). Therefore

    z0 w.r.t. vertex i  =  z0_stored + vtx_z[0] - vtx_z[i]

exactly, with no approximation. This identity is the whole mechanism; it is verified against
the skim source and must be re-verified if the skim's vertex handling ever changes.

**(b) d0 is vertex-independent.** `d0 = idTrk->d0()` is the transverse impact parameter with
respect to the **beamline**, not to a vertex. Vertices in a pp event differ from one another
essentially only in z (beam-spot transverse size ~10 um vs a 2 mm d0 cut), so the d0 cut is the
same for every vertex and is applied once. The vertex-dependent part of the selection is
entirely the z0 sin(theta) term. This is standard ATLAS practice (d0 w.r.t. beamline, z0 w.r.t.
the vertex) and it also keeps the d0 cut byte-identical to nominal.

**(c) Which vertices are eligible.** Entries of the skim's `PrimaryVertices` dump with
`vtx_ntrk >= 2`. The container is stored unfiltered and always contains exactly ONE dummy
beam-spot vertex with `ntrk == 0`; measured on pp24 it is always the LAST entry and its
`vtx_z` equals `vtx_z[0]` exactly in 100% of events, so without the `ntrk >= 2` guard it would
act as a duplicate of the primary. Same guard as `Muon.h::n_vtx` /
`PythiaFullSimExtras.c`.

**(d) Vertex assignment and the "secondary" definition.** A pair is kept if ANY eligible
vertex satisfies (a)+(b) for both muons. The vertex recorded on the pair
(`matched_vtx_ind`) is the eligible vertex minimising max(|z0_v sin(theta)|) over the two
muons — i.e. the best-matching vertex — with ties broken by the lower index, so index 0 wins
whenever it is as good as any other. A pair is counted as **secondary-vertex** when
`matched_vtx_ind != 0`.

Vertex-assignment ambiguity is small by construction. Measured on `data_pp24_part11.root`,
first 300 000 events, **eligible (`vtx_ntrk >= 2`) vertices only**, averaging over vertex
*pairs*: only **0.60 %** of vertex pairs sit within 2 mm in z, against a `vtx_z` RMS of
**58.3 mm** and a mean |z_i - z_j| of **67.3 mm**. (Two reviewers reproduced these three
numbers to three digits; the values 0.73 % / 61.0 / 68.3 quoted here before 2026-09-08 were
not reproducible under any definition and are withdrawn. The sample and the
"eligible-only, over pairs" definition are part of the statement — keeping that definition and
changing only eligible->all gives 10.3 % / 57.9 mm / 60.3 mm instead. Do not confuse the
pair-averaged mean |z_i - z_j| quoted here with the mean distance to the PRIMARY,
|z_i - z_0| over i >= 1, which is a different quantity: 65.4 mm eligible-only, 48.1 mm over all
vertices.)

A related property this design rests on, **empirical rather than structural**: the returned
index can only be the primary if vertex 0 is itself eligible. Were `vtx_ntrk[0] < 2` ever to
occur, a pair the old primary-vertex cut kept could be dropped, breaking the strict-superset
guarantee. Measured over all 1 696 134 events of `data_pp24_part11.root`: **zero** such events,
minimum `vtx_ntrk[0] = 2` — as expected of the highest-sum-pT^2 vertex. Re-check if the skim's
vertex ordering or cleaning changes.

**(e) Statistics record.** Fraction of kept pairs with `matched_vtx_ind != 0`, computed from
the NTuple-processing OUTPUT pair tree (never from the raw NTUPs), for pp24, reported
- in total, and
- in each of the 8 canonical `ParamsSet::pair_pt_coarse_bins` (8-150 GeV, log),
for the full pair-tree population and, separately, for the single-b signal region.

**(e2) A known, accepted asymmetry in the SINGLE-muon population.** In data the single-muon
tree is filled from inside the pair loop, after `PassCuts()`, so a data single muon inherits the
requirement that *its partner* point at the same vertex; in MC the per-muon `pass_medium` /
`pass_tight` flags use only the ANY-vertex form (§3f). This is not new in kind — the old code
put "partner must pass the primary-vertex cut" on the data side too — and it enters the trigger
efficiency as a population definition, not as a ratio bias, since ε(mu4) is essentially flat in
|z0 sinθ| across a 2 mm window. The PAIR quantities, which is what ε_reco is, match exactly.
Recorded so it is not rediscovered as a bug.

**(f) MC mirror.** The pp-conditions fullsim MC (`PythiaFullSimExtras` pp path,
`PowhegFullSimExtras`) has genuine pile-up — track-bearing vertices 1:5.0%, 2:14.7%, 3:21.8%,
4:22.7%, 5:17.1%, ... , compared with data's 1:6.4%, 2:17.9%, 3:25.1%, 4:23.0%, 5:15.3%. Its
reco muons get the identical all-vertex treatment, so `eps_reco` is measured on the same pair
selection the data applies. The HIJING **overlay** has exactly one track-bearing vertex in
100% of events, where the procedure is an exact no-op; it is left on the primary-vertex path
explicitly rather than relying on that.

### 4. Negative constraints (things the code must NOT do)

1. Do **NOT** change Pb+Pb. Pb+Pb rejects pile-up at event level; its pairs must stay
   primary-vertex pairs and its NTuple output must be byte-identical. The d0/z0 decision is
   made dispatchable and only the pp path overrides it.
2. Do **NOT** recompute d0 per vertex (see §3b) — that would be a silent deviation from the
   nominal d0 definition for a sub-micron-scale effect.
3. Do **NOT** treat the `ntrk == 0` dummy beam-spot entry as a vertex (§3c).
4. Do **NOT** re-invent any binning. The statistics table uses
   `ParamsSet::pair_pt_coarse_bins` (8 bins); read it, never retype it.
5. Do **NOT** apply the all-vertex procedure to the HIJING overlay (§3f).
6. Do **NOT** let the pp trigger-efficiency tag-and-probe keep the old trees: its input pair
   trees come from this same NTuple stage, so eps^nc must be refit on the NEW trees.

## Context

- Gap-cut change of the same day: `muon_gap_cuts_acceptance.md` F17/F18. Code edited,
  nothing rerun. This task's rerun therefore also delivers the gap-cut rerun.
- The pp24 crossx physics procedure is owned by `pp24_crossx_rerun_2026_08.md`.
- Rerun blast radius: `docs/signal_selection_change_impact.md` — note this change breaks that
  doc's §1 assertion that "NTuple processing ... UNCHANGED"; the doc must be amended.

## Scope

IN: pp data NTuple stage; pp-conditions fullsim MC IP cut; the statistics record; the pp24
rerun chain (NTuple -> trig-eff refit -> fullsim reco-eff + MC trig-eff -> dR fits ->
pair_reco_eff -> crossx + MC-data comparison plots).

OUT: Pb+Pb anything; R_AA; unfolding; template fits; the Pb+Pb pair-eta axis migration
(deferred, F18).

## Design Decisions

- **D1 (2026-09-07).** The d0/z0 test is extracted from `PassCuts_DataCore` into a
  dispatchable `PassD0Z0Hook`, with the DataCore implementation (primary vertex) as the
  default and a `PPExtras::PassD0Z0Extra()` override. Reason: `PassCutsHook` is
  `PassCuts_DataCore(...) && (CallPassCuts<Extras>(), ...)` — DataCore runs FIRST and
  short-circuits, so a plain `PassCutsExtra` could never rescue a pair the primary-vertex cut
  had already rejected. Pb+Pb implements no override, so its path is unchanged.
- **D2 (2026-09-07, user).** Statistics on the canonical 8-bin `pair_pt_coarse_bins`.
- **D3 (2026-09-07, user).** Mirror the procedure in the pp-conditions fullsim MC.
- **D4 (2026-09-07, user).** Today's uncommitted MINUIT-limits change to the Step-3 `expo`
  dR fit is INCLUDED in the rerun and goes through /review-analysis-code.
- **D5 (2026-09-07, user).** Crossx Stage 7 draws pp only this round; the Pb+Pb panels are
  deferred with the rest of the Pb+Pb migration.

## Implementation Plan

1. NTuple stage: dispatchable d0/z0 hook + `PPExtras::PassD0Z0Extra()` + `vtx_*` branch
   binding + `matched_vtx_ind` on the pp pair. (per §3a-d) -> /review-analysis-code.
2. MC mirror in `PythiaFullSimExtras` (pp path) + `PowhegFullSimExtras`. (per §3f)
   -> /review-analysis-code.
3. Small local test on one pp24 file + one fullsim file before any Condor submission.
4. pp24 NTuple Condor rerun: nominal (2mu4) + trigger-efficiency mode.
5. Statistics record from the new pair trees. (per §3e)
6. pp trigger-efficiency refit on the new trees.
7. Pythia fullsim pp pipeline (reco-eff + MC trig-eff), dR fits, pair_reco_eff rebuild.
8. pp24 crossx pipeline (Stage 7 pp-only) + MC-data comparison plots. -> /review-plot.
9. Docs + commits.

## Progress Log

- 2026-09-08 Step 1 DONE — **pp data NTuple stage**. The d0/z0 decision is now dispatchable:
  `DimuonDataAlgCoreT.h` gains `PassD0Z0_DataCore()` (the verbatim old primary-vertex
  expression, Pb+Pb's default), `CallPassD0Z0<E>` and `PassD0Z0Hook()`;
  `DimuonDataAlgCoreT.c:617-619` replaces the inline test with `if (!PassD0Z0Hook()) return
  false;` and `:578` holds the DataCore body. The algorithm itself lives in the NEW
  `Utilities/AllVertexIPSelection.h` (`PassAnyVertex`, `BestCommonVertex`, `IsEligibleVertex`,
  `Z0SinTheta`) — ONE implementation, shared with the MC so the two can never drift.
  `PPExtras.{h,c}` binds `vtx_z`/`vtx_ntrk` and implements `PassD0Z0Extra()`.
  `MuonPairReco.h` gains a pp-only `PairPPExtras` mixin with `matched_vtx_ind` and
  `pass_primary_vtx` (deliberately NOT in the shared `PairDataExtras`, so Pb+Pb's tree is
  untouched).

  **Non-obvious ROOT mechanism, measured twice (the first two accounts were both wrong).**
  These chains run in `SetMakeClass(1)` mode. `SetBranchAddress` on an STL-collection branch
  leaves the pointer NULL if and ONLY IF **both** of these hold: the chain's first tree is
  already loaded **and** the branch is currently DISABLED. The 2x2 measurement, on the real
  chain (reproduced independently by the iteration-2 reviewer):

  | | branch enabled | branch DISABLED |
  |---|---|---|
  | bound BEFORE first tree load | rc = 5 (kNoCheck), allocated | rc = 5, allocated |
  | bound AFTER first tree load | rc = 3 (kMakeClass), allocated | rc = 3, **NULL** |

  My first diagnosis blamed the load ordering alone (it is not sufficient); the iteration-1
  reviewer, varying only the other axis, concluded the trap did not exist at all (it does).
  **The rule that follows is simply `SetBranchStatus(name, 1)` BEFORE `SetBranchAddress`** —
  which is exactly what `Utilities/tchain_helpers.h::enable_and_bind` already does, and why all
  three call sites now use it. Ordering against the tree load is then irrelevant, which matters
  because `GetBranch()`, `GetListOfBranches()` and `GetEntries()` all load the first tree, so
  "after the load" is the normal state by the time any Extra hook runs. The truth table lives in
  `Utilities/AllVertexIPSelection.h`; `bind_branch` now also accepts rc == 5, which it used to
  reject — so the working order had been unusable through that helper.

- 2026-09-08 Step 2 DONE — **MC mirror**. `PythiaFullSimExtras.{h,c}`: `UseAllVertexIP()` is
  true iff `getFullSimSampleType() == FullSimSampleType::pp`; `vtx_z`/`vtx_ntrk` bound inside
  `bind_reco` on every chain (before any tree load), the `store_mc_trigger` re-bind of
  `vtx_ntrk` now guarded by `!has_vtx_ntrk` so it cannot swap a good deferred address for a
  MakeClass one; `PassMuonMediumCuts` uses the ANY-vertex form; `PassPairAllVertexIP` adds the
  SAME-vertex requirement to `pair_pass_{medium,tight}`; `FinalizeExtra` reports the
  secondary/added fractions. `PowhegFullSimExtras.{h,c}`: identical treatment, gated on
  `!getIsFullsimOverlay()`. All three analysis-class headers compile clean.

  **Why the per-muon flag is loosened to ANY-vertex rather than dropped:** the IP cut has to
  leave `PassMuonMediumCuts` (a per-muon primary-vertex cut inside `pass_medium` would veto
  exactly the secondary-vertex pairs the change exists to keep), but dropping it outright would
  make the single-muon MC reco-efficiency carry no IP cut at all. The ANY-vertex form is the
  single-muon analogue, and since the pair form implies it, `pair_pass_{medium,tight}` are
  exact.

- 2026-09-08 Steps 4b + 5 DONE — **pp24 data Condor rerun + the statistics record**.
  Both pp24 modes resubmitted and finished clean in 18 min: cluster 10570
  (`run_pp_24_nominal.sub`, trigger_mode=3 -> `muon_pairs_pp_2024_part*_2mu4_mindR_0_02.root`)
  and cluster 10571 (`run_pp_24.sub`, trigger_mode=1 -> the `_single_mu4_mindR_0_02_res_cut_v2`
  variant, distinct paths — the two clusters do NOT collide). 24/24 jobs, no errors in any
  `.err` beyond the pre-existing TStreamerInfo warnings, 12+12 outputs present. hadded to
  `muon_pairs_pp_2024_2mu4_mindR_0_02.root` (652 MB) + the cut-acceptance hists. Statistics
  record produced from it — see Results & Observations. Headline: **1.400 ± 0.014 % of
  signal-region OS pairs are ADDED by the all-vertex rule, but with a factor-~9 pair-pT
  dependence (2.73 % -> 0.3 %)**, so the cross-section change is a shape change, not a
  normalisation.

- 2026-09-08 Step 3 (partial) — **local tests**. pp: `PPAnalysis(24,11)`, `trigger_mode=3`,
  `is_test_run=true`, 300 000 events -> 2 641 pairs (SS 777, OS 1 864), of which
  **58 (2.20 %) have a secondary best-matching vertex** and **42 (1.59 %) fail the primary
  vertex outright**, i.e. are ADDED by the change. That 1.59 % agrees with the independent
  raw-NTUP scouting estimate (+1.57 % signal region, +1.9 % nominal), which is the cross-check
  that the `z0_v` identity is right. Refactoring onto the shared header reproduced the numbers
  byte-for-byte. Pb+Pb: `PbPbAnalysis(23,1)` runs unchanged, no `matched_vtx_ind` branch, no
  all-vertex message — the compile-time dispatch keeps it on `PassD0Z0_DataCore`.

- 2026-09-07 Step 0: doc created; Autonomy Contract pinned. Doc triage read
  `INDEX.md`, `pp24_crossx_rerun_2026_08.md`, `muon_gap_cuts_acceptance.md`.
  Two read-only scouting subagents mapped (a) the rerun/staleness chain and (b) the
  NTuple implementation sites; their scratch docs are
  `_sub_pp_rerun_staleness_1.md` and `_sub_pp_allvtx_sites_1.md` (merged below, then
  deleted). Four user decisions taken (D2-D5 above, plus the 8-bin binning).

## Results & Observations

### THE STATISTICS RECORD (Done item 3) — pp24, full sample, Tight WP

Produced 2026-09-08 by `plotting_codes/single_b_analysis/pp24_secondary_vertex_stats.cxx` from
the hadded ntuple-processing output `muon_pairs_pp_2024_2mu4_mindR_0_02.root` (all 12 batches,
rerun with the all-vertex selection AND the 2026-09-07 gap cuts). CSV:
`plots/single_b_analysis/pp24/pp24_secondary_vertex_fraction.csv`.

`f_sec` = the pair's own best-matching vertex is not the primary. `f_new` = the pair fails the
primary vertex outright, i.e. it is ADDED by this change. `f_sec >= f_new` always.

**Integrated**

| population | N | f_sec | f_new |
|---|---|---|---|
| all pairs (OS+SS) | 4 920 638 | 2.314 ± 0.007 % | **1.657 ± 0.006 %** |
| all OS pairs | 3 482 923 | 2.296 ± 0.008 % | 1.659 ± 0.007 % |
| all SS pairs | 1 437 715 | 2.357 ± 0.013 % | 1.651 ± 0.011 % |
| **single-b signal region (OS)** | **661 237** | 2.138 ± 0.018 % | **1.400 ± 0.014 %** |

**Signal region (OS), canonical `ParamsSet::pair_pt_coarse_bins` [GeV]**

| pair pT | N | f_sec | f_new |
|---|---|---|---|
| 8.00-11.54 | 224 389 | 3.506 ± 0.039 % | **2.727 ± 0.034 %** |
| 11.54-16.65 | 275 901 | 1.669 ± 0.024 % | 0.914 ± 0.018 % |
| 16.65-24.01 | 117 877 | 1.079 ± 0.030 % | 0.409 ± 0.019 % |
| 24.01-34.64 | 34 046 | 0.937 ± 0.052 % | 0.303 ± 0.030 % |
| 34.64-49.97 | 7 526 | 0.811 ± 0.103 % | 0.345 ± 0.068 % |
| 49.97-72.08 | 1 300 | 1.077 ± 0.286 % | 0.308 ± 0.154 % |
| 72.08-103.98 | 179 | 0.559 ± 0.557 % | 0.000 + 1.023 % |
| 103.98-150.00 | 16 | 0.000 + 10.869 % | 0.000 + 10.869 % |
| >= 150 (above axis) | 3 | 0.000 + 45.865 % | 0.000 + 45.865 % |

Errors are nested-binomial; at k = 0 (and k = N) the 68 % Clopper-Pearson width is given
instead, because 0/16 is "no events yet", not a measurement of exactly zero. The above-axis row
exists so the per-bin rows exhaust the integrated total — the macro now ASSERTS that
(`check_complete`, throws) rather than claiming it in a comment, which is precisely the claim
that turned out to be false in the first version (3 pairs above 150 GeV were being dropped).

(The all-pairs-vs-pair-pT table, including the 2 591 627 pairs below the 8 GeV axis where
`f_new` reaches 2.116 ± 0.009 %, is in the CSV.)

**PHYSICS — the effect is strongly pair-pT DEPENDENT, and in the direction it must be.**
`f_new` falls monotonically from 2.73 % in the lowest signal-region bin to ~0.3 % above 25 GeV,
a factor ~9. That is the expected behaviour and is a real consistency check on the whole
procedure: the primary vertex is the highest-sum-pT^2 vertex in the event, so a collision hard
enough to make a high-pair-pT dimuon almost always WINS that ranking and is therefore the
primary — it can only be "secondary" if some other collision in the same crossing was harder
still. A soft collision making a low-pair-pT dimuon loses the ranking often. Hence the
correction is largest exactly where the cross-section is largest.

**CONSEQUENCE for the cross-section: this is a SHAPE change, not only a normalisation.**
dsigma/dpair-pT rises by ~2.7 % in the first coarse bin and by ~0.3 % in the highest populated
ones, so the spectrum steepens slightly less than before. It must NOT be quoted as "a 1.4 %
overall increase" — 1.4 % is only the signal-region integral. Any downstream ratio that assumed
a flat normalisation shift (R_AA included, once Pb+Pb is brought over) has to take the
pair-pT dependence.

The one non-monotonic point, signal-region 49.97-72.08 with `f_sec` 1.077 % against 0.811 % in
the bin below, is 14 pairs against an error of 0.286 % — under 1 sigma, not a feature.

### Scouting numbers (pre-implementation, raw-NTUP estimate)

From `data_pp24_part11.root` (1 696 134 events), mirroring `PassCuts_DataCore` (Medium WP,
2mu4 trigger match, mindR 0.02; no resonance veto):

| population | pass primary | rescued by other vertices | fraction |
|---|---|---|---|
| nominal pp24 (Medium) | 20 778 | 398 | **+1.915 +- 0.096 %** |
| nominal pp24 (Tight) | — | — | +1.796 % |
| single-b signal region | 2 099 | — | **+1.572 +- 0.274 %** |

Rescuing vertex index: [1] 378, [2] 19, [3] 1. 147 of the 1 289 primary-failing pairs fail
|d0| and are therefore unrescuable by any vertex (as expected — d0 is vertex-independent).
These are SCOUTING numbers; the authoritative statistics come from the NTuple output (§3e).

## Rerun plan (Steps 4-8), in dependency order

Established by a read-only scouting subagent (`_sub_pp_rerun_staleness_1.md`, merged and
deleted) and amended for THIS change. The amendment matters: that subagent concluded
`SKIP_NTP=1` / `SKIP_CONDOR=1` were safe because "the cuts live entirely at the RDF stage" —
true for the gap-cut change alone, and **false now**. The all-vertex change IS an NTuple-stage
change, so every NTuple stage must actually run, in data AND in MC.

| # | What | Command | Why it must precede the next |
|---|---|---|---|
| 4a | pp24 NTuple, trigger-eff mode (`trigger_mode=1`, `run_pp_24.sub`) | inside `pipelines/pipeline_pp_trig_eff.sh` | the T&P probes come from this tree |
| 4b | pp24 NTuple, nominal mode (`trigger_mode=3`, `run_pp_24_nominal.sub`) | inside `pipelines/pipeline_pp_crossx.sh` | the crossx pairs |
| 6 | data mu4 turn-on refit, Tight + Medium | `pipelines/pipeline_pp_trig_eff.sh`; `SAMPLES="pp" pipelines/run_data_trigeff_medium_wp.sh` | crossx reads eps^nc; ALSO stale on its own from the F17 coarse q·η edge 2.3 -> 2.2 |
| 7a | Pythia pp24 fullsim, NTuple + RDF + MC trig-eff | `ENABLE_MC_TRIG_EFF=1 pipelines/pipeline_pythia_fullsim_pp.sh full` (**NOT** `SKIP_NTP=1`) | supplies eps_reco, the MC trig-eff hists, **the detector response AND the template-fit MC** — all four come from this one production, so all four go stale together |
| 7b | ΔR correction fits | `SAMPLES="pp_full" STEPS="3" METHODS="expo polyu_fixedRp" SIGNS="os" PLATEAU_MODES="nocorr_ptmerge" WPS="tight" pipelines/run_dr_correction_fits.sh` | reads 7a's hists; crossx reads its output |
| 7c | pair reco-efficiency file | `root -l -b -q 'build_pp24_fullsim_pair_reco_eff.C+(true)'` in `plotting_codes/reco_effcy/` | not in any pipeline; crossx reads `pair_reco_eff_pp24_full.root` |
| ~~7d~~ | ~~POWHEG pp17 fullsim rerun~~ | **NOT NEEDED — see below** | its only live product is truth-level |
| 8 | pp24 crossx + plots + MC-data comparison | `pipelines/pipeline_pp_crossx.sh` | last |
| 5 | statistics record | `root -l -b -q 'pp24_secondary_vertex_stats.cxx+()'` | needs 4b's hadded tree |

**Wider MC blast radius than eps_reco alone.** `UseAllVertexIP()` keys off
`FullSimSampleType::pp`, i.e. the *whole* pp24 fullsim production — which also supplies the
**detector response** and the **template-fit MC**, not only ε_reco and the MC trigger
efficiency. Those outputs are regenerated by the same Stage 5 of `pipeline_pythia_fullsim_pp.sh`
and so are covered by row 7a, but they must be named, because a reader looking only at "reco-eff
+ crossx" would not expect the detector response to move. For POWHEG the gate is
`!getIsFullsimOverlay()`, which switches the change ON for the **pp17** POWHEG fullsim as well
(its NTUPs do carry `vtx_z`/`vtx_ntrk`, verified, so it cannot throw).

**Why POWHEG is code-changed but NOT rerun.** `PowhegFullSimExtras` was mirrored for
consistency (so the two fullsim samples cannot drift), but its reco-level products have no live
consumer: the POWHEG entry in the MC-vs-data comparison is the node
`df_single_b_pass_signal_truth_gapcut_weighted`
(`RDFBasedHistFillingPowhegFullsim.cxx:309`), which is **purely truth-level** —
`from_same_b && truth_minv ∈ (1.08,2.9) && truth_pair_pt > 8` plus the truth gap cuts — with no
`reco_match` and no `pair_pass_*`. The nodes the IP cut does reach
(`_pass_{medium,tight}_weighted`, `_pass_*_and_signal_truth_and_reco_weighted`, and the
single-muon trees feeding the Run-2 mixed-pair study) belong to the POWHEG fullsim reco-efficiency
line, which pp24 stopped using on 2026-08-18 when the genuine Pythia-fullsim pair ε_reco replaced
it (`project_mc_sample_roles`: POWHEG fullsim obsolete; POWHEG truth is kept for the NLO
template). So the POWHEG trees on disk keep the primary-vertex selection and the code no longer
matches them — **acceptable only because nothing live reads those nodes. Rerun
`pipeline_powheg_fullsim_single_muon.sh` + `pipeline_powheg_fullsim_mixed_pairs.sh` BEFORE
reviving any POWHEG fullsim RECO product.**

Three stale inputs verified by reading their axes: `dr_correction_fits_pp24_full_step3_{expo,
polyu_fixedRp}_os_nocorr_ptmerge.root` (2026-08-24) and `pair_reco_eff_pp24_full.root`
(2026-08-18) all carry pair-η ±2.4 against the canonical ±2.2, and
`DrCorrectionCrossxEvaluator.h::CheckCanonicalBinning()` / `PairRecoEffEvaluator.h::
CheckCanonicalBinning()` THROW on them — the chain fails loudly, not silently. The eps^nc file
was already refit to `2_00_TO_2_20` on 2026-09-07, but it is stale AGAIN here because its input
tree changes.

**Clobber risks accepted** (all are intended reruns): fullsim Stage 10 overwrites
`plots/pp_trigger_efficiency/mc_based/` (backs up first); `run_dr_correction_fits.sh` writes
into the same tree with NO backup; `build_pp24_fullsim_pair_reco_eff.C` uses RECREATE with no
backup; crossx Stage 4 `rm -f`s the combined files before hadd. Disk: the data fileset is at
~1.37 TB of a ~1.54 TB quota (halved figures), and the NTuple reruns overwrite in place, so the
net change is ~0.

## Remaining Work

Steps 4-8 of the Implementation Plan, per the table above.

## Latest Stage

**DATA HALF DONE; MC HALF HELD ON A CONCURRENT SESSION (2026-09-08).**

Done: Steps 1, 2, 3, 4a, 4b, 5 — code written, reviewed once and amended for all 10 findings,
committed (2bb281d, 99a9efa, ee46d2c, 568542e); both pp24 Condor modes rerun clean and hadded;
the statistics record produced at full statistics (see Results & Observations).

**HELD — do not start without checking:** Steps 6, 7a, 7b, 7c, 8 (eps^nc refit, Pythia fullsim
NTuple+RDF+MC trig-eff, dR fits, pair_reco_eff rebuild, crossx). **A SECOND Claude session is
working in this same checkout on the MC trigger efficiency** — `FillMCTrigEffPairEff.cxx`,
`FillMCTrigEffClosure.cxx`, `Utilities/PairTrigEffEvaluator.h` and a new ACTIVE tracking doc
`mc_trigeff_single_value_pair_eff.md` ("single-value pair 2mu4 efficiency in the top pair-pT
cells"), all written 2026-09-08 00:31-00:41 with their own ACTIVE Autonomy Contract. The
collision is not merely a race: **Step 7a's fullsim NTuple rerun MOVES the fullsim pair tree
their whole measurement is built on** (by ~1.2 %), and Stage 10 + the dR fits overwrite
`mc_trig_eff_hists_*`, `dr_correction_fits_*` and `plots/pp_trigger_efficiency/mc_based/` —
none of which git can see (`.gitignore` drops `*.root *.png`). USER DECISION 2026-09-08: hold
the MC half, finish the data half, resume when that session reports done. Second user decision
the same day: the cross-section keeps the **dR-fit** trigger weight (with today's MINUIT-limits
change, which still needs `/review-analysis-code`); the single-value pair-efficiency method is
NOT adopted into the crossx — consistent with that doc's own Done list ("NOT wired into the
cross-section").

Also in flight: `/review-analysis-code` iteration 2 on the amended code + the statistics macro.
