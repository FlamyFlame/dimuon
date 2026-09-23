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
| same-sign (SS) | 24–40, 40–70, 70–125, 125–300 | ∃ SAME-SIGN μμ pair: m < 10 GeV, pT^μμ > 20 GeV (**24–40: NO pT^μμ cut**, D6) | 0.5 M / slice |

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
3. [x] Commit (source/athena branch; JOs + scripts in dimuon_codes repo).
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

**D6 — USER DECISION (2026-09-22, "final"): SS pTH24_40 has NO pair-pT cut.** With pT^μμ > 20
the slice's ε ≈ 1.9e-7 (~37 CPU-h per accepted event, ~2 100 CPU-yr for 0.5 M); without it
ε ≈ 1.2e-5 (R3). JO renamed `mc.Py8EG_A14_pp_hQCD_DiMu_pTH24_40_pdf_SS_m10.py`, `PairPtMin = -1`
(filter condition `dilepPt > MinDilepPt` then always true). Grid test sizes (user): 1 000 events for
SI 70–125/125–300 and SS 70–125/125–300, 100 for SS 40–70 and SS 24–40.

**D7 — Step-5 plots: filter-level definitions, σ·ε-weighted combination.** Macro
`analysis/plot_test_sample_sanity.C` reads the downloaded NTUPs (`truth_muon_*`, `truth_mupair_*`)
with the FILTER's muon/pair definition (pT>3.7 GeV, |η|<2.5) — it validates the generator filter,
so the analysis 4.5/2.4 fiducial would be the wrong definition here (§3c). "pT-hat combined" and the
per-slice overlays are weighted by σ·ε_filt/N per slice, with σ and ε read from THIS sample's own
grid `log.generate` files (these test datasets have no AMI record; no other production's numbers
are used). Sanity counts are raw per slice (`sanity_table.csv`).

**D8 — Grid chain.** `grid/chain.py` (resumable, `grid/chain_state.json`): evgen done →
`Derivation_tf TRUTH0` (pathena, Athena 25.0.72, outDS `user.yuhang.JOtest.<phys>.TRUTH0.v1`) →
`TruthTrigRates.py` (pathena from `SkimCode/run_evgen_high_pT_additional/`, EvtMax −1, outDS
`user.yuhang.NTUP.JOtest.<phys>.v1`) → `rucio download` to
`~/usatlasdata/evgen_high_pT_additional_test/<config>/` (NTUP + evgen `log.generate`). Exits on any
non-`done` terminal task state.

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
- **2026-09-22, SS pTH125_300 local.** ε = 1.32e-4 (10 / 75 956; est 1.26e-4 ✓), 262 s CPU per
  accepted event; `RequireSameCharge: True` printed; 10/10 events hold a passing SAME-SIGN pair,
  the three 3-muon events show OS pairs skipped (patch works as intended).
- **2026-09-22, step 3 DONE.** athena `71ef28ce20d` (branch `dimuon-same-sign-filter-23.6.68`),
  cherry-picked as `4b66ba5e232` onto `dimuon-same-sign-filter` (main-based, for the MR);
  dimuon_codes `9109bde`.
- **2026-09-22, PAUSE (user's own rule: warn if 1 000 events > 1–2 days).** R1: every config
  exceeds ~1.5 CPU-days per 1 000 events; SS pTH24_40 ~6 CPU-years, SS pTH40_70 ~117 CPU-days.
  Grid submission NOT done; question put to the user. Local 10-event jobs for SI pTH70_125,
  SS pTH70_125, SS pTH40_70, SS pTH24_40 left running (projected ~2.5 h, ~4 h, ~1 d, ~3 weeks).
- **2026-09-22, USER: "Stop."** Questions: (1) why CPU-years when a private ~1e-7 sample was
  generated before — how long is one unfiltered event? (2) why is SS pTH24_40 so hard; ε without the
  pair-pT cut? Answered with R2/R3 (`estimate/breakdown_ss24_40.py`). Grid remains NOT submitted.
- **2026-09-22, local 70–125 done.** SI pTH70_125 ε = 5.5e-5 (10/181k, est 4.4e-5), SS pTH70_125
  ε = 4.1e-5 (10/244k, est 3.1e-5) — both within the ±32 % of a 10-event count; EVNT check 10/10
  pass for both. 4 of 6 configs now locally validated; SS 40–70 / 24–40 still running (0 accepted
  after ~630k generated each, consistent with ε ≈ 4e-6 / 2e-7).
- **2026-09-22, USER decision.** Grid test sizes: 100 events for SS 40–70 and SS 24–40 (5–10 per
  job), 1 000 for the other 4. Asked for a time estimate before submitting (see R5).
- **2026-09-23, step 4a — grid evgen submitted** (`grid/submit_evgen.sh`, `prun --useAthenaPackages`
  from the 23.6.68 build, worker `grid/run_gen.sh`, seeds `%RNDM=<base>`; outputs
  `EVNT.pool.root,log.generate`; log `grid/submit_evgen.log`):

  | config | evts/job × jobs | seed base | jediTaskID | outDS |
  |---|---|---|---|---|
  | SI pTH125_300 | 250 × 4 | 120000 | 52691377 | user.yuhang.JOtest.Py8EG_A14_pp_hQCD_DiMu_pTH125_300_pdf_m10_pT45.v1 |
  | SS pTH125_300 | 200 × 5 | 160000 | 52691387 | user.yuhang.JOtest.Py8EG_A14_pp_hQCD_DiMu_pTH125_300_pdf_SS_m10_pT20.v1 |
  | SI pTH70_125 | 100 × 10 | 110000 | 52691402 | user.yuhang.JOtest.Py8EG_A14_pp_hQCD_DiMu_pTH70_125_pdf_m10_pT45.v1 |
  | SS pTH70_125 | 50 × 20 | 150000 | 52691417 | user.yuhang.JOtest.Py8EG_A14_pp_hQCD_DiMu_pTH70_125_pdf_SS_m10_pT20.v1 |
  | SS pTH40_70 | 5 × 20 | 140000 | 52691431 | user.yuhang.JOtest.Py8EG_A14_pp_hQCD_DiMu_pTH40_70_pdf_SS_m10_pT20.v1 |
  | SS pTH24_40 | 10 × 10 | 130000 | 52691447 | user.yuhang.JOtest.Py8EG_A14_pp_hQCD_DiMu_pTH24_40_pdf_SS_m10.v1 |

  Expected CPU/job (25 ms/generated): ~7 / 11 / 13–16 / 9–11 / ~10 / ~6 h. Local 10-event test of the
  revised SS 24–40 JO launched in parallel (`local_test/run_same_sign_pT24_40/`; the pT>20 attempt
  moved to `local_test/superseded_same_sign_pT24_40_pT20/`).
- **2026-09-23, step 6 prerequisite — truth-only external-sample mode (subagent, worktree).**
  `PythiaAlgCoreT.{h,c}` gain `truth_input_override` (glob/list, ONE slice kn), `truth_output_dir_override`,
  `truth_output_label`, `truth_raw_counts` (all default off; override ⇔ raw_counts guarded both ways;
  throws on fullsim/overlay/private). Output `<dir>/muon_pairs_pythia_truth_<label>_kn<N>_no_data_resonance_cuts.root`;
  weight ≡ 1, AMI read + nPDF DSID guard skipped; `meta_tree_out` `nproc_kin<N>_beam0`/`nbeam_kin<N>_beam0`.
  Branch binding factored into `BindTruthPairBranches_PythiaCore` (line-by-line identical). Test on
  the local 10-event kn5 NTUP: 3 SS + 11 OS pairs = hand count, weight 1, nproc 10. Doc section in
  `Analysis/docs/pythia_truth.md`. Commit `bd81e27`; `/review-analysis-code` iteration 1 **PASS**
  (0 C/W, 4 INFO: failed input guard can leave an empty diagnostic dir in the override dir; no
  E_COM tag in the output name → put it in the label; pre-existing zero `nentries_*` meta branches;
  nproc = N_to_process). Merged `26328ad`; merged state recompiled with ACLiC in the main checkout.
  Run recipe (per slice): `.L PythiaAnalysisClasses.h; PythiaTruthAnalysis py(kn,false,5.36);
  py.truth_input_override="<NTUP glob>"; py.truth_output_dir_override="<dir>";
  py.truth_output_label="<label incl. 5p36>"; py.truth_raw_counts=true; py.Run();`
- **2026-09-23, step 6 tooling.** `Analysis/plotting_codes/trig_effcy/mc_based/mc_request_stats_boost.cxx`
  (D5) + `analysis/run_truth_processing.sh` (runs the new mode per config: labels `5p36_SI_m10_pT45`,
  `5p36_SS_m10_pT20`, `5p36_SS_m10` for 24–40). Truth analogue of Step 3 = `TruthFiducial()` +
  `ForwardLowPtVeto()` on truth-aliased legs + the gap terms rebuilt from the same `ParamsSet`
  primitives with the pair term on `truth_pair_eta` (`pair_eta` is a real RECO branch in the trees and
  cannot be aliased). Test on the 10-event stand-in: the macro's reference Step-3 reco counts for
  pTH125_300 m<10 reproduce `cells_pTH125_300_m_lt_10.csv` EXACTLY (SS 527/264/13, 115/47/0, 17/3/0;
  OS 6361/3007/314, 2831/1254/126, 643/257/19). Pending: real inputs, then `/review-analysis-code`
  (macro) and `/review-plot` (step-5 figures).
- **2026-09-23, local SS pTH24_40 (revised, no pT^μμ cut).** ε = 1.28e-5 (10 / 779 815; est
  1.22e-5 ✓), 2 087 s CPU per accepted event locally (26.8 ms/generated); `MinDilepPt: -1`,
  `RequireSameCharge: True` printed; EVNT check 10/10 pass. 5 of 6 configs locally validated; the SS
  pTH40_70 local job is still running.

## Results & Observations

### R2 — Cost of one generated (unfiltered) event

SI pTH125_300 job: 1 419 s CPU total (PerfMonMT; configure/init excluded) for 38 715 generated
events → **36.7 ms CPU per generated event** (Pythia8 + EvtGen + per-event xAOD truth conversion in
`prefiltSeq` + the two filters; component-level PerfMon was disabled, so the split is not
measured). SS pTH125_300: 262 s × 75 956⁻¹·10 → 34.5 ms. Cost per accepted event = 36 ms / ε.

### R3 — SS pTH24_40: what each requirement costs (emulation, 720 000 NTUP events, ε_≥2μ = 1.82e-4)

| requirement on ≥1 pair (legs pT>3.7, \|η\|<2.5) | f | ε_total | CPU / accepted |
|---|---|---|---|
| any pair (= existing filter) | 1 | 1.82e-4 | 198 s |
| OS | 0.81 | 1.47e-4 | 245 s |
| SS | 0.21 | 3.8e-5 | 950 s |
| SS, m<10 | 0.067 | **1.2e-5** | 2 900 s |
| SS, pT^μμ>20 | 3.0e-3 | 5.5e-7 | 6.6e4 s |
| SS, m<10, pT^μμ>10 | 2.6e-2 | 4.8e-6 | 7.6e3 s |
| SS, m<10, pT^μμ>15 | 5.7e-3 | 1.0e-6 | 3.5e4 s |
| SS, m<10, pT^μμ>20 (JO) | 1.0e-3 | 1.9e-7 | 1.9e5 s |
| OS, m<10, pT^μμ>20 | 5.7e-2 | 1.05e-5 | 3.4e3 s |

The pair-pT cut is the killer: given SS m<10, pT^μμ>20 keeps 1.5 % (740 / 48 337). With
pT-hat 24–40 a 20 GeV low-mass pair must carry most of one jet's momentum, and SS needs two muons
from one jet that are NOT the b→μ⁻ / b→c→μ⁺ cascade (which is OS).

### R4 — The EXISTING e8599 `_pdf` production (≥2μ filter only): AMI + PanDA

AMI (`ami_info_PDF/`) genFiltEff, totalEvents, sumOfWeightsNoFilter (= generated events); PanDA
(bigpanda JSON, finished Gen_tf jobs of the evgen task = `prodsysIdentifier_0`, median of ≤300
jobs) CPU per job:

| slice | DSID | genFiltEff | kept events | generated | evts/job | CPU/job | CPU/kept evt | CPU/generated | total CPU |
|---|---|---|---|---|---|---|---|---|---|
| 8–14 | 803020 | 5.93e-6 | 2.04 M | 3.44e11 | 50 | 49.8 h | 3 585 s | 21.2 ms | ~232 CPU-yr |
| 14–24 | 803016 | 4.36e-5 | 3.6 M | 8.25e10 | 200 | 39.1 h | 703 s | 30.7 ms | ~80 CPU-yr |
| 24–40 | 803017 | 1.82e-4 | 2.4 M | 1.32e10 | 1000 | 38.1 h | 137 s | 25.0 ms | ~10 CPU-yr |
| 40–70 | 803018 | 5.26e-4 | 1.2 M | 2.28e9 | 1000 | 13.4 h | 48 s | 25.5 ms | ~1.8 CPU-yr |
| 70–125 | 803019 | 1.24e-3 | 0.32 M | 2.58e8 | 2000 | 13.4 h | 24 s | 29.9 ms | ~0.2 CPU-yr |
| 125–300 | 803015 | 2.40e-3 | 0.32 M | 1.34e8 | 5000 | 13.0 h | 9.3 s | 22.4 ms | ~0.1 CPU-yr |

**The previous request cost ≈ 320 CPU-years**, dominated by pTH8_14 (ε ≈ 6e-6) — so the ATLAS
production system already ran ~1e-6-level filters at this scale. Per GENERATED event the cost is
flat, 21–31 ms (my VM: 35–37 ms). Re-costing the new request at 25 ms/generated:
SI 70–125 ≈ 18, SI 125–300 ≈ 3, SS 125–300 ≈ 3, SS 70–125 ≈ 13, **SS 40–70 ≈ 110, SS 24–40 ≈
2 100 CPU-yr** (vs 320 for the whole previous request).

### R5 — Wall-clock of the existing low-ε productions, and what 0.5 M of the hard SS slices means

PanDA (bigpanda JSON) + AMI:
- **802769** np pTH8_14 (nPDF, ε = 4.96e-6, 472 400 events = 2 tasks × 240 k): task 45576771 =
  2 400 jobs × 100 events, created 2025-07-11, ended 2025-07-29 (**18 days**); per job median
  **55 h CPU** (p10 36 h, p90 77 h), wall max 96 h; ≈ 9.9 ms per generated event;
  ≈ 2 000 s CPU per kept event; ≈ 15 CPU-yr for the 240 k.
- **803020** pp pTH8_14 `_pdf` (ε = 5.9e-6, 2.04 M): ONE task 50094092 = 40 800 jobs × 50 events,
  2026-05-04 → 2026-06-10 (**37 days**); 49.8 h CPU per job; 3 585 s per kept event; ≈ 232 CPU-yr.

New 0.5 M requests (CPU per kept event at 10–25 ms/generated):
- SS 24–40 (no pT^μμ cut, ε ≈ 1.2e-5): 0.8–2.1 ks per kept event ≈ the 8–14 slices' cost;
  0.5 M ≈ 13–33 CPU-yr (≈ 1/7–1/17 of 803020). nEventsPerJob ≈ 100 → ~55 h jobs (as 802769).
- SS 40–70 (pT^μμ > 20, ε ≈ 3.6e-6 est.): 2.8–6.9 ks per kept event (≈ 1–2× 803020);
  0.5 M ≈ 44–110 CPU-yr (≈ 0.2–0.5 × 803020). nEventsPerJob ≈ 25 → ~48 h jobs, ~20 000 jobs.

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
Step 4: 6 grid evgen tasks RUNNING; `grid/chain.py` (background) chains TRUTH0 → skim → download.
Steps 5/6 code ready and tested on stand-ins (`analysis/plot_test_sample_sanity.C`,
`analysis/run_truth_processing.sh`, `mc_request_stats_boost.cxx`). On download: run sanity macro →
`/review-plot`; run truth processing → boost macro → `/review-analysis-code`; results → R6/R7.
