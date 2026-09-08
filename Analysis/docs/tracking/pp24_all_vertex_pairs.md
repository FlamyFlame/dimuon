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
`matched_vtx_ind != 0`. Vertex-assignment ambiguity is small by construction: only 0.73% of
vertex pairs are within 2 mm in z (`vtx_z` RMS 61.0 mm, mean |z_i - z_0| 68.3 mm).

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

  **Non-obvious mechanism found the hard way (a ROOT trap, documented at both call sites):**
  the input chain runs in `SetMakeClass(1)` mode. If `SetBranchAddress` for an STL-collection
  branch is issued AFTER the chain's first tree has been loaded, ROOT takes the immediate
  MakeClass path (`rc = kMakeClass = 3`), treats the argument as the address of the DATA rather
  than of a pointer-to-vector, never allocates the `std::vector`, and would write leaf bytes
  over the member itself. Bound before the first load it is deferred (`rc = kNoCheck = 5`) and
  `TChain::LoadTree` allocates correctly. `GetBranch()` itself loads the tree — which is why the
  vertex ADDRESSES are set in `PPExtras::PerformTChainFill` (right after `Add`) while the
  existence check and `SetBranchStatus` stay in `InitInputBranchesDimuonAnalysisExtra`.
  This is also why every DataCore branch is bound before the first `GetBranch()` call.

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

**STEPS 1-2 DONE, STEP 3 (local tests) IN PROGRESS (2026-09-07/08).** Data + MC code written
and compiling; pp and Pb+Pb smoke tests pass; the MC smoke test is running. Next: statistics
macro (Step 5 tooling), then the Condor reruns.
