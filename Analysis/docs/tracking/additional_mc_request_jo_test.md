# Additional high-pT / same-sign MC request: JO + GeneratorFilters patch test

**Mode:** Implementation. **Created:** 2026-09-22. **Session:** "JO test".
**Working dir:** `PythiaCentralProdJO/evgen_high_pT_additional/` (all relative paths below are
w.r.t. it unless they start with `Analysis/`, `SkimCode/` or `~`).
**Parent:** `pthat_slice_mass_stats_sample_request.md` (ACTIVE — sized this request; owns the 4
reference 3×3 tables `cells_{pTH70_125,pTH125_300}_{m_lt_10,m_sig}.csv` in
`~/usatlasdata/dimuon_data/plots/pp_trigger_efficiency/mc_pthat_slice_mass_stats/`).

## Objective

Test the GeneratorFilters patch (`xAODDiLeptonMassFilter::RequireSameCharge`) and the 6 job options
of the additional MC request before handing them to the ATLAS production expert:

| group | pT-hat slices [GeV] | filter (truth, both legs pT>3.7 GeV, \|η\|<2.5) | requested events |
|---|---|---|---|
| sign-independent (SI) | 70–125, 125–300 | ∃ μμ pair: m < 10 GeV, pT^μμ > 45 GeV | 1 M / slice |
| same-sign (SS) | 24–40, 40–70, 70–125, 125–300 | ∃ SAME-SIGN μμ pair: m < 10 GeV, pT^μμ > 20 GeV | 0.5 M / slice |

## Autonomy Contract (ACTIVE — re-read on every compaction)
- Mandate: run autonomously to DONE; do NOT pause to confirm progress. Finishing a
  plan, a passing small test, or one pipeline stage is NOT a stopping point.
- Done = (1) patch compiled against AthGeneration 23.6.68; (2) 10-event local Gen_tf test per
  config (6) passes the event-by-event filter check (≥2 muons, pair pT, max mass, same sign when
  imposed), with filter efficiency + CPU time/event recorded and a projected time for 1 000
  events; (3) patch + JOs + scripts committed; (4) grid tests: 1 000 events/config EVNT → truth
  AOD → NTUP (SkimCode truth-only skim) downloaded — OR paused with a warning to the user if any
  config is projected > 1–2 days; (5) 1 000-event sanity check (same 4 conditions) + plots: per
  sample (SI, SS; pT-hat combined) nMuons/event, pair pT (9–150 GeV, log x), pair mass (0–40 GeV),
  pair sign (2 bins); per sample pair pT with pT-hat slices overlaid (2 / 4 lines); (6) projected
  statistics boost in the 4 reference 3×3 tables from 1 M (SI) / 0.5 M (SS) events per slice.
- Stop-and-ask = ANY physics-results-bending ambiguity (no fixed list; use judgment;
  when unsure whether an ambiguity is blocking, treat it as blocking → AskUserQuestion).

## Physics Procedure (AUTHORITATIVE)

**§1 Motivation.** The pp24 MC trigger efficiency (single-value pair eff and ε_ΔR) is
statistics-starved in the top pair-pT cells, and the same-sign (SS) population is poor down to
pair pT ≈ 25 GeV (`pthat_slice_mass_stats_sample_request.md` R5, R7; `pp24_trig_eff_hybrid`
D1/D2 are TEMPORARY pending this MC). A generator-level low-mass (m<10 GeV) dimuon filter removes
the back-to-back population; the SS filter buys SS pairs directly.

**§2 Generator.** Identical to the pp24 FULL sample (DSIDs 803015–803020, e8599 `_pdf`): Pythia8
A14, `NNPDF23_lo_as_0119_qed`, HardQCD:all, pT-hat slice, √s = 5.36 TeV, pp beams,
`xAODMultiMuonFilter` (≥2 muons pT>3.7 GeV, |η|<2.5) — checked against
`/cvmfs/.../MCJobOptions/803xxx/803015/`. NEW: `xAODDiLeptonMassFilter` on `TruthMuons`,
MinPt 3.7 GeV, MaxEta 2.5, 0 < m < 10 GeV, pT^μμ > 45 (SI) / 20 (SS) GeV, AllowSameCharge=True,
RequireSameCharge=True for SS only (patch: skip a pair if pdgId1·pdgId2 < 0).

**§3 Steps.**
- (a) *Local test* (10 events/config, Gen_tf, `--ecmEnergy=5360`): event-by-event check from the
  EVNT truth: ≥2 stable muons pT>3.7, |η|<2.5; ∃ pair (SS pair for SS configs) with m<10 GeV and
  pT^μμ above the cut. Filter efficiency = accepted / generated from the log (the filter's own
  count and the Gen_tf summary); time/accepted-event from the log → projected time for 1 000.
- (b) *Grid test*: 1 000 accepted events/config, EVNT → truth xAOD (the user asked for `Reco_tf.py`;
  if it cannot take EVNT input, `Derivation_tf.py --formats TRUTH0`, the chain the existing truth
  NTUPs were made with) → NTUP via SkimCode truth-only config (`TruthTrigRates.py`, IsEvgen, no
  trigger, no reco).
- (c) *Sanity check + plots* on the 1 000-event NTUPs: truth-level. "Muon" and "pair" for the
  filter check = the FILTER's definition (stable truth muon, pT>3.7 GeV, |η|<2.5) — this is a
  validation of the generator filter, not analysis physics.
- (d) *Statistics boost* of the 4 reference tables (reco-level Step-3 pairs, per slice, per sign,
  m<10 / 1.08<m<2.9, top-3 pair-pT × 3 |η^pair| cells). Method: TBD in Design Decisions before
  implementation (truth-only new samples vs reco-level reference tables).

**§4 Negative constraints.**
- No reconstruction anywhere in this test (truth only).
- AMI weights: these are NEW datasets, but the test outputs are raw counts/shape checks only; no
  σ-weighted combination across slices is made here (would need the new samples' own AMI records).
- The reference binning is read from the reference CSVs / `ParamsSet` (never retyped).
- Counts are raw per slice; never summed across slices as yields.

## Context
- Patch in `source/athena` (branch `dimuon-same-sign-filter`, based on upstream `main`
  `a028797ea07`): `xAODDiLeptonMassFilter.{h,cxx}`, `share/common/xAODDiLeptonMassFilter_Common.py`.
- Base release: AthGeneration 23.6.68 (Pythia 8.317; the original sample used 23.6.47 / 8.308 —
  flagged to the user 2026-09-21).
- References: `../evgen_4isospins` (JO test scripts, grid `prun --useAthenaPackage`),
  `SkimCode/run_pythia_truth` (truth skim), `.claude/worktrees/.../PythiaTruthDerivation/scripts`
  (EVNT→TRUTH0 pathena).

## Implementation Plan
1. [x] Compile patch (`setup_AthGeneration.sh`) — check build + filter property visible.
2. [ ] Local Gen_tf 10 events × 6 configs (parallel, background) → event-by-event check script
   (§3a), filter eff, time/event, projected time for 1 000 events.
3. [ ] Commit (source/athena branch; JOs + scripts in dimuon_codes repo).
4. [ ] Grid: Gen_tf 1 000 events/config (split over jobs if slow; WARN USER + pause if >1–2 days);
   EVNT→AOD; truth skim (new `SkimCode/run_evgen_high_pT_additional/`); download.
5. [ ] Sanity check + plots (§3c).
6. [ ] Statistics boost (§3d).
7. [ ] Wrap-up: results here, INDEX, commit.

## Design Decisions

**D1 — Local build on a 23.6.68-based branch.** The checkout's branch `dimuon-same-sign-filter`
sits on `main` HEAD `a028797ea07` (2026-09-21), whose `GeneratorFilters/src/Common.cxx` needs an
`IMCTruthClassifier::particleTruthClassifier_full` the 23.6.68 release lacks → compile error. The
patch was moved (uncommitted diff, stash/pop) onto new branch `dimuon-same-sign-filter-23.6.68`
from `d525b4dfe72` (= 23.6.68's `ReleaseData` nightly commit). The 6 files the JOs touch
(xAODDiLeptonMassFilter .h/.cxx/_Common.py, xAODMultiMuonFilter .cxx/_Common.py,
CreatexAODSlimContainers.py) are byte-identical between the two bases, so the tested code = the MR
code. The main-based branch stays for the MR (cherry-pick).

**D2 — Package filter.** `build/package_filters.txt` = `+ Generators/GeneratorFilters`, `- .*`.
The previous build dir (configured for `--sourcedir=../source` with only `- athena/.*`) built
nothing; with `../source/athena` and no filter, `Projects/WorkDir` was configured twice. Old build
dirs moved aside to `build_stale_20260922/`, `build_failed_cfg_20260922/` (not deleted).

**D3 — JO renaming (ATLAS naming rule, not physics).** Gen_tf refuses any JO whose physicsShort
(name between `mc.` and `.py`) exceeds 50 characters; the original names were 55–65. Renamed for the
test (USER may override the scheme): `..._pdf_SameSign_maxMass10_DiMupT20` → `..._pdf_SS_m10_pT20`
(47–49 chars), `..._pdf_maxMass10_DiMupT45` → `..._pdf_m10_pT45` (45–46 chars).

**D4 — EVNT→AOD in release 25.** SkimCode `TrigRates.cxx` (since `85de057`) reads truth
`uid()`; a DAOD_TRUTH0 made with the Feb-2026 24.0 nightly throws `ExcBadAuxVar ::uid` in the
current SkimCode build (AthAnalysis 25.2.89), tested on the mc15_5TeV example file. So the EVNT→AOD
step runs in Athena 25.0.72.

**D5 — Statistics-boost method (plan; implement at step 6).** Projected reco count in a reference
cell c (slice s, sign, mass range) = `N_ref,reco(c) × [N_req · ν_new(c)] / [N_ref · ν_ref(c)]`, i.e.
the reference table scaled by the ratio of TRUTH-level pair rates. ν = truth pairs per event in c
(truth pair pT/|η| cell, truth mass window, truth sign), with the truth fiducial of
`PassCuts_PythiaCore` (both legs truth pT>4.5, |η|<2.4) plus truth analogues of the Step-3 leg cuts
(forward low-pT veto, gap windows on truth q·η, |η^pair|<2.2). Assumes the truth→Step-3 reco
efficiency per cell is the same in filtered and unfiltered samples (cells are kinematic).
Sources (ntuple-processing output only): ν_ref from `~/usatlasdata/pythia_fullsim_full_sample/
muon_pairs_pythia_fullsim_pp24_no_data_resonance_cuts_mc_trig_full.root`, `muon_pair_tree_kin<K>_sign<d>`
(truth-seeded: every truth pair written, reco only flags it; `PythiaFullSimExtras.c:366-604`),
N_ref = `meta_tree_out/nproc_kin<K>_beam0`; ν_new from the truth-only processing
(`PythiaTruthAnalysis`) run on the new NTUPs — which needs a NEW mode (input dir / file list /
labels override, raw-counts mode skipping the AMI read + nPDF DSID guard; hard-coded today at
`PythiaAlgCoreT.c:65-67,284,323-357`). SS slices 24–40 / 40–70 have no reference table: their
projected counts are reported against the same-slice reference computed the same way (kin2, kin3).

## Progress Log

- **2026-09-22, step 1.** Compiled on AthGeneration 23.6.68 (x86_64-el9-gcc15-opt) after D1/D2;
  `libGeneratorFilters.so` built, `RequireSameCharge` present in the generated
  `GeneratorFiltersConf.py`. Log: `compile.log`.
- **2026-09-22, step 2 started.** `local_test/run_local_tests.sh 10` (6 configs in parallel,
  seeds 1001–1006, `--ecmEnergy=5360`); control file symlinked into each `run/run_*/` (Gen_tf needs
  includes inside `--jobConfig`; the production DSID dirs will need a real copy, as 803015 has).
- **2026-09-22, step 2 (partial).** First run failed for all 6 on the physicsShort ≤ 50 rule → D3
  rename, relaunched 18:47. SI pTH125_300 finished (wall 1547 s, 6 jobs sharing an 8-core VM):
  filter eff 2.58e-4 (10 / 38 715), 139 s CPU per accepted event; `check_evnt_filter.py` on the
  EVNT: 10/10 events pass (every event ≥2 filter muons and a pair with m<10, pT^μμ>45).
  Printout of the patched filter confirms `RequireSameCharge: False` for SI.
- **2026-09-22, chain test on that EVNT.** `Reco_tf.py` (25.0.72) rejects `--inputEVNTFile`
  (unrecognized argument) → EVNT→AOD = `Derivation_tf.py --formats TRUTH0` (25.0.72, rc 0), then
  SkimCode `TruthTrigRates.py` (AthAnalysis 25.2.89, rc 0) → NTUP, 10 entries, `truth_muon_*` and
  `truth_mupair_*` branches present. Dir: `local_test/aod_test/`.
- **2026-09-22, filter-efficiency estimate (`estimate/emulate_filter.py`).** The new pair filter
  emulated on the existing pp truth-only NTUPs (nPDF e8599, every event already passed the ≥2μ
  filter; `truth_mupair_*` = all status-1 μμ pairs, MeV) × AMI genFiltEff of the same slice's
  `_pdf` DSID. Validation: SI pTH125_300 est 2.75e-4 vs measured 2.58e-4 (10/38 715) ✓.
  Output `estimate/emulate_filter.out`; numbers → R1.

## Results & Observations

### R1 — Filter efficiency and CPU cost (the feasibility question)

Generation cost per ATTEMPTED Pythia event t₀ ≈ 139 s × 2.58e-4 = 0.036 s CPU (measured, SI
pTH125_300, this VM); per accepted event ≈ t₀ / ε.

| config | f (pair filter \| ≥2μ) | ε_total est. | CPU / accepted evt | 1 000 evts | requested N | CPU for requested N |
|---|---|---|---|---|---|---|
| SI pTH70_125 | 3.5e-2 | 4.4e-5 | ~820 s | ~230 CPU-h (9.5 d) | 1 M | ~26 CPU-yr |
| SI pTH125_300 | 1.1e-1 | 2.6e-4 (meas.) | 139 s (meas.) | ~39 CPU-h | 1 M | ~4.4 CPU-yr |
| SS pTH24_40 | 1.0e-3 | 1.9e-7 | ~1.9e5 s | ~53 000 CPU-h (6 CPU-yr) | 0.5 M | ~3 000 CPU-yr |
| SS pTH40_70 | 6.8e-3 | 3.6e-6 | ~1.0e4 s | ~2 800 CPU-h (117 d) | 0.5 M | ~160 CPU-yr |
| SS pTH70_125 | 2.5e-2 | 3.1e-5 | ~1 200 s | ~330 CPU-h (14 d) | 0.5 M | ~19 CPU-yr |
| SS pTH125_300 | 5.3e-2 | 1.3e-4 | ~290 s | ~80 CPU-h (3.3 d) | 0.5 M | ~4.6 CPU-yr |


## Remaining Work

## Latest Stage
Step 2: local 10-event tests running (`local_test/<config>/gen_tf.out`). Next: event-by-event
check script on the EVNT, filter eff + time/event.
