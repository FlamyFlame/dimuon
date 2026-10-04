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

## Autonomy Contract (DONE 2026-09-29 — all Done items met)
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
- AMI weights: these are NEW datasets without AMI records. Quoted numbers (filter checks, ε, boost
  tables) are raw per slice. The ONLY cross-slice combination is the DISPLAY of the "pT-hat combined"
  sanity figures, weighted by σ·ε_filt/N from each test sample's OWN grid `log.generate` files (USER
  APPROVED 2026-09-29, D7); no other production's weights are ever used.
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
4. [x] Grid: Gen_tf 1 000 events/config (split over jobs if slow; WARN USER + pause if >1–2 days);
   EVNT→AOD; truth skim (new `SkimCode/run_evgen_high_pT_additional/`); download.
5. [x] Sanity check + plots (§3c). `/review-plot` PASS it. 2.
6. [x] Statistics boost (§3d).
7. [x] Wrap-up: results here, INDEX, commit.

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

**D5 REVISED (2026-09-25) — pooled reco/truth efficiency.** Old: n_new(c) = n_ref,reco(c) ×
N_req ν_new(c) / (N_ref ν_ref(c)) per slice and sign. New: n_new(c) = N_req · ν_new(c) · ε(c), with
ε(c) = Σ reco / Σ truth over kn2–5 and both signs (truth → Step-3-reco efficiency incl. migration of
the cell). Reason: the per-slice/per-sign ratio is noisy in sparse reference cells and ZERO wherever
the reference has no reco pair — exactly the cells the request is for. Physics: ε(c) is a detector
(kinematic) quantity of the cell; it should not depend on the pT-hat slice or pair sign. Check:
populated cells agree between the two methods (pTH125_300 OS m<10 [52,74)×[0,1): boost 30.3 vs 30.5);
pooled ε = 0.51–0.80 across the 9 cells. Reviewer checks (`/review-analysis-code` it. 1–3): per-sign ε within 8 % of the pooled one in populated cells (OS <1 %, SS 1–8 % higher, e.g. m<10 [52,74)×[1,2): SS 0.851 vs OS 0.798 → pooled ε is conservative for SS); per-slice kn4 vs kn5 differ up to ~13 % (m_sig [74,105)×[0,1) OS, 3.5σ; m<10 same cell 9 %, 3.6σ; pooled ε conservative for kn4 there). The assumption holds to ≲13 %, far inside the ±7–100 % test-sample errors. ε binomial error is not propagated (matters only in the forward top cell, ε = 3/6).

**D9 — Same-sign request mass cut = m < 10 GeV (USER, 2026-09-29).** An uncommitted edit of another
session in `Analysis/Utilities/MCRequestLooseMassFilter.h` adds `kLooseMassMaxSameSign = 100` ("user,
2026-09-17: 100, not 10"); the user confirmed 10 GeV is correct for this request (as the JOs and tests
use). The file is left untouched (another session's working-tree change); that constant is stale.

**D10 — Step-6 output = 4 tables of projected TOTAL counts (USER, 2026-10-02).** Old: one long-format
`stats_boost.csv` (per slice × config, boost factors). New (user): 4 CSVs
`projected_counts_{m_lt_10,m_sig}_{SS,OS}.csv`, 3 pair-pT rows × 3 |η^pair| columns, each cell = projected
RAW pair count n_ref + Σ n_new **summed over the request's slices kn2–5** (user chose counts over a boost
factor and over σ-weighted N_eff). Raw sum = a statistics count, not a yield. 0-test configs enter with
n_new = 0 (cell = lower value; the header lists them with the cell's 68 % CL upper value). Per-slice
boosts stay in `stats_boost_stdout.txt`.

**D11 — Projection of the SS pTH125_300 pT^μμ>50 (0.2 M) config without its own truth test sample (USER,
2026-10-04).** pT dependence: the pT>20 SS test sample's events whose filter-level SS pair (legs ≥3.7 GeV,
|η|≤2.5, 0<m<10) has pT^μμ>50 (90/1000 events; exactly the events the pT50 filter keeps), ν = pairs/90.
|η^pair| dependence: SS cells as they are in 52–74 / 74–105; in 105–150 the SS bin total is split by the OS
|η^pair| fractions of the slice's (kn5) test samples, pooled. An empty OS shape bin → 0-test cell with UL.
OS cells: the subset's OS pairs, direct.

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
- **2026-09-23, local SS pTH40_70 done.** ε = 3.19e-6 (10 / 3 131 262; est 3.6e-6 ✓), 8 729 s CPU
  per accepted event locally (27.9 ms/generated); EVNT check 10/10 SS pass. **All 6 configs locally
  validated (Done item 2 met).**
- **2026-09-23, grid: first two configs through the chain** (`grid/chain.log`): SS pTH24_40 evgen done →
  TRUTH0 → skim 52697315 → downloaded 08:46; SI pTH125_300 evgen 52691377 done 11:19 → TRUTH0
  52701612 → skim 52702229 → downloaded 11:42. Sanity macro: SI pTH125_300 1 000/1 000 events pass,
  ε = 1 000/3 652 565 = 2.74e-4; SS pTH24_40 100/100 pass, ε = 100/7 276 136 = 1.37e-5. Others
  at 4–9 of 5–20 jobs finished (19:33).
- **2026-09-29, `/review-plot` iteration 1 FAIL** (2 WARNING): pair-pT filter cuts inside bins (fixed:
  edges anchored at 20/45 GeV); §4 vs D7 (fixed after user approval, see §4). INFO recorded: out-of-range
  pairs SI 2.0 % below 9 GeV, 0.35 % ≥150 GeV, 3.5 % ≥40 GeV mass; SS 36 % below 9 GeV (nearly all SS
  24–40, no pT^μμ cut) — the SS pT panel covers 64 % of the pairs of its mass/sign panels; the SS
  combined figures are ~94 % the SS 24–40 slice by σ·ε weight.
- **2026-09-29, step 5 DONE.** `/review-plot` iteration 2 PASS (after edges anchored on the cuts and
  §4 approval). Figures `~/usatlasdata/evgen_high_pT_additional_test/plots/{SI,SS}_combined_sanity.png`,
  `{SI,SS}_pair_pt_by_pthat.png`, `sanity_table.csv`. Display pT axis = local piecewise-log (9→20→45→150
  GeV, 8/9/13 bins), not a canonical binning (none applies to a generator-filter spectrum). SS 24–40
  has 0 pairs above 20 GeV in its 100 events (≈1.6 expected, Poisson fluctuation).

- **2026-10-02, follow-up (user).** (a) D10: `mc_request_stats_boost.cxx` writes the 4 tables; rerun
  → `~/usatlasdata/evgen_high_pT_additional_test/stats_boost/projected_counts_*.csv`; `stats_boost.csv`
  deleted. Cells = Σ over slices of (n_ref_reco + n_new_sum) of the old CSV's ALL rows EXACTLY (independent
  re-sum). (b) Pair-pT plots: `SetLogx` was already called (`plot_test_sample_sanity.C:197` gPad, `:219`
  canvas) and the PNGs were log-x; the overlay canvas now uses `gPad->SetLogx()/SetLogy()` too; figures
  regenerated (identical content). Plot-reviewer criteria gained R5 "a log axis is made ONLY with
  SetLogx/y/z" (`.claude/commands/review-plot.md`, `conventions/atlas-plotting.md` R2a,
  `agents/plot-reviewer.md`).

- **2026-10-02, follow-up 2 (user).** User keeps the test-sample projection (D5 REVISED); the reference-scaling
  alternative offered for 0-test cells is NOT adopted. (1) Production bundle
  `~/usatlasdata/evgen_high_pT_additional_test/for_production/<physicsShort>/` = JO + control file (real copy) +
  ONE grid `log.generate` (user, simplified: the job whose ε is closest to the combined one, within 0–5 %; combined ε per config = R6). Staged grid JOs
  byte-identical to `run/`. `evgenConfig.nEventsPerJob = 500` in all 6 JOs is still the placeholder. (2) R8.

- **2026-10-03, follow-up 3 (7th config SS pTH125_300 m<10 pT^μμ>50).** JO
  `run/run_same_sign_pT125_300_pT50/mc.Py8EG_A14_pp_hQCD_DiMu_pTH125_300_pdf_SS_m10_pT50.py` (diff vs pT20 JO:
  PairPtMin 50000, nEventsPerJob 50). Local job seed 1101: ε = 2.20e-5 (10/453 840), 32 ms CPU/generated, wall
  4.2 h; EVNT check 10/10 pass (SS, m 0.7–9.5, pT^μμ 52–82). Expected ε: emulation on 720 k existing nPDF
  pTH125_300 truth events × genFiltEff 2.396e-3 → **1.34e-5** (±2 % stat; same method gives 1.264e-4 for the
  pT20 config vs grid 1.247e-4) → seed 1101 is a +64 % fluctuation (6.1 expected, 10 seen, p≈9 %). Not
  representative for the production log → 6 more 10-event seeds (1102–1107) in parallel; pick the log closest to
  the combined ε (same rule as the other configs). `estimate/emulate_filter_pT50.{py,out}`.

- **2026-10-04, follow-up 3 DONE.** Seeds 1102–1107: ε = 1.72/1.22/1.26/1.05/1.87/1.58 e-5; combined 7 jobs
  **ε = 1.47e-5 ± 12 %** (70 / 4 764 458), consistent with the 1.34e-5 emulation. EVNT check 70/70 pass, no OS
  pass. Production log = seed 1107 (1.58e-5, +8 % vs combined) → `for_production/Py8EG_A14_pp_hQCD_DiMu_pTH125_300_pdf_SS_m10_pT50/`
  (JO + control + log.generate). Local CPU ≈ 35 ms/generated → ≈ 2.4 ks per kept event; at grid 41 ms ≈ 2.8 ks
  → 50 evts/job ≈ 39 h; 0.5 M ≈ 44 CPU-yr (R9's 58 used ε 1.1e-5).

- **2026-10-04, nEventsPerJob set (user).** From R6 grid CPU per kept event, targeting ~30–50 CPU-h jobs (803020
  precedent 50 evts/~50 h): SI 70–125 200 (31 h), SI 125–300 1000 (40 h), SS 24–40 100 (53 h), SS 40–70 20 (47 h),
  SS 70–125 200 (48 h), SS 125–300 500 (46 h), SS 125–300 pT50 50 (39 h). JOs in `run/` and `for_production/`
  identical. The production logs were made with the placeholder (500) / test maxEvents; only the comment+value changed.
  athena branch `dimuon-same-sign-filter` (`4b66ba5e232`) pushed to origin `yuhang/athena` for the MR to atlas/athena main.

- **2026-10-04, follow-up 4 DONE.** (1) D11 in `mc_request_stats_boost.cxx` (`/review-analysis-code` PASS; its
  2 WARNINGs fixed: empty OS-shape bin → UL, stale OS-note; cells unchanged by the fixes). Subset 90/1000;
  `ev_num` → NTUP match validated pair by pair; subset holds all SS top-bin pairs (asserted). SS50 m<10 SS
  contributions = 0.2 M × {43,13,2}/90 × ε (top bin split 25:9:1 OS); m_sig {9,1,1}/90 (OS split 8:5:0).
  Tables (summed over |η|), before → after: m<10 SS 41.2k→117.4k / 12.0k→35.0k / 0.8k→4.1k; m_sig SS
  10.2k→26.2k / 2.0k→3.8k / 0.4k→2.2k; m<10 OS +≤3 % (52–105); m_sig OS unchanged. (2) Helper plot
  `plots/ss_os_pair_pt_45_150.png` (`analysis/plot_ss_os_pair_pt_slope.C`, `/review-plot` PASS; overlap fixed
  by a ±1.2 % x offset, fiducial in headline). SS/OS σ ratio in the canonical cells: pTH125_300 0.055±0.008 /
  0.037±0.010 / 0.020±0.014 (SS steeper, ~2σ); pTH70_125 0.019±0.006 / 0.013±0.013 / 0 (30 SS pairs: too few).
  SS 125–300 first display bin (45–49 GeV) SS/OS = 0.17 is a ~3σ fluctuation (not correlated pairs).

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

### R6 — Grid test (1 000 / 100 events per config): filter check, efficiency, measured CPU

All 6 evgen tasks done, TRUTH0 + skim done, NTUPs + `log.generate` downloaded to
`~/usatlasdata/evgen_high_pT_additional_test/<config>/` (`grid/chain.log`, `grid/chain_state.json`).
Filter check (`analysis/plot_test_sample_sanity.C`, `plots/sanity_table.csv`): **0 failing events in
4 200** (every event has ≥2 filter muons and a qualifying pair, same-sign where required).

| config | evts | ε_filt (grid) | jobs | CPU/job median (max) | total CPU | ms / generated | CPU / kept evt | projected CPU for requested N |
|---|---|---|---|---|---|---|---|---|
| SI pTH70_125 | 1 000 | 4.31e-5 (1000/23.19 M) | 10 | 15.9 h (22.2) | 155 h | 24 | 558 s | 1 M → ~18 CPU-yr |
| SI pTH125_300 | 1 000 | 2.74e-4 (1000/3.65 M) | 4 | 10.0 h (10.9) | 40 h | 39 | 144 s | 1 M → ~4.6 CPU-yr |
| SS pTH24_40 (no pT^μμ cut) | 100 | 1.37e-5 (100/7.28 M) | 10 | 5.6 h (7.4) | 53 h | 26 | 1 908 s | 0.5 M → ~30 CPU-yr |
| SS pTH40_70 | 100 | 4.0e-6 (100/25.0 M) | 20 | 12.6 h (21.2) | 234 h | 34 | 8 424 s | 0.5 M → ~134 CPU-yr |
| SS pTH70_125 | 1 000 | 3.21e-5 (1000/31.15 M) | 20 | 12.0 h (16.5) | 242 h | 28 | 871 s | 0.5 M → ~14 CPU-yr |
| SS pTH125_300 | 1 000 | 1.25e-4 (1000/8.02 M) | 5 | 18.8 h (20.6) | 92 h | 41 | 331 s | 0.5 M → ~5.2 CPU-yr |

Test total 816 CPU-h (estimate was 650–720). Whole request ≈ **205 CPU-yr** at these grid speeds,
comparable to the ≈ 232 CPU-yr of the existing 803020 pTH8_14 sample alone (R4/R5); SS pTH40_70 is
2/3 of it. ε agrees with the local 10-event runs and the emulation (R1) in every config.

### R7 — Projected statistics boost of the 4 reference tables (+ SS 24–40 / 40–70 slices)

`mc_request_stats_boost.cxx` (D5 REVISED; `/review-analysis-code` PASS it. 4); outputs (D10)
`~/usatlasdata/evgen_high_pT_additional_test/stats_boost/projected_counts_{m_lt_10,m_sig}_{SS,OS}.csv`
(projected total counts, slices 24–300 summed) + per-slice boosts in `stats_boost_stdout.txt`.
Boost = (n_ref + Σ n_new)/n_ref per slice, for 1 M (SI) + 0.5 M (SS) events; per-cell test-sample
errors ±7–100 %; cells with 0 test pairs are quoted as upper values, never as "no gain".

- **pTH125_300, opposite sign:** ×23–42 in every cell (m<10 and signal window) — e.g. [52,74)×|η|<1:
  6 361 → ~194 000 (×30.5). Matches the expectation (N_req/N_ref)(ε_old/ε_new) ≈ 27.
- **pTH125_300, same sign, m<10:** ×39–81 in the cells with ≥ 2 test pairs ([52,74)×|η|<1: 527 → ~22 800,
  ×43; single-test-pair cells 23.8 (≤ 76) and 124.5 (≤ 408)); signal window ×20–64; the forward (2.0–2.2) and top-pT cells are < ×100–400 (0 test pairs).
- **pTH70_125, opposite sign:** m<10 ×33–101, signal window ×67–165 (≈ 90 expected; the low values are
  single-test-pair cells); several top-pT / forward cells are 0-test (upper values).
- **pTH70_125, same sign:** n_ref is tiny (34/27/3, 2/0/0, 0/0/0 for m<10); projected ≈ 4 000 pairs in
  each of the two lowest-|η| cells of [52,74) (×118–149); the other cells are 0-test in the 1 000-event
  test sample (upper values ~400–900 per config).
- **SS pTH24_40 / pTH40_70:** contribute NOTHING measurable in the top 3 pair-pT cells (all 0-test; the
  reference is also ~0) — their pairs sit below ~50 GeV pair pT, i.e. they serve the 20–52 GeV region
  of R7 in `pthat_slice_mass_stats_sample_request.md`, which the 4 reference tables do not cover.
- Method checks: reference Step-3 reco counts reproduce `cells_*.csv` exactly; pooled ε 0.50–0.83;
  truth-rate ratio 28.2 (kn5) / 89.6 (kn4) vs filter-efficiency expectation 27.3 / 89.8.

### R8 — Why SS pTH125_300 gives few top-pT pairs (`analysis/ss_pair_pt_fractions.cxx`, 2026-10-02)

Every accepted event has a qualifying SS pair (1000/1000; 3 with two), but its pair-pT spectrum falls
steeply — not flat. Leading qualifying pair (filter level, legs 3.7/2.5), SS pTH125_300 test, 1 000 events:
20–30 54.4 %, 30–40 24.4 %, 40–52 13.6 %, 52–74 5.8 %, 74–105 1.6 %, 105–150 0.2 %, ≥150 0.
(SS pTH70_125: 71.8 / 22.0 / 4.8 / 1.3 / 0.1 / 0 %.) After the truth Step-3 analogue: 43 / 13 / 2 pairs in
the top 3 bins, then × ε_reco ≈ 0.8. So 0.5 M → ≈ 17 k / 5 k / 0.8 k SS pairs (summed over |η|) from this
config. SI pTH125_300 (pT^μμ>45): 7.7 % of events hold an SS pair, but at high pT (Step-3: 19 / 8 / 0).
Physics: an SS low-mass pair needs two muons from separate decay chains in one jet (g→bb̄ with b→μ⁻ +
b̄→c̄→μ⁻, or B⁰ mixing), each carrying a small fraction of the jet pT; the in-slice pT-hat spectrum also
falls from 125 GeV. Doubling the SS pTH125_300 request doubles only its contribution; a higher pair-pT
threshold instead keeps the same evgen CPU per top-bin pair but needs far fewer events to simulate.

### R9 — Extra SS sample with pT^μμ > 50 GeV (estimate, 2026-10-02, `analysis/ss_pair_pt_threshold_fractions.C`)

Fraction of SS (pT^μμ>20) test events whose leading filter SS pair (legs 3.7/2.5, m<10) exceeds 40/45/50/60
GeV: pTH125_300 0.212/0.144/0.090/0.047; pTH70_125 0.062/0.033/0.019/0.004. ⇒ ε(SS, m<10, pT^μμ>50) ≈
1.1e-5 (125–300, ±11 %) and 6.1e-7 (70–125, ±23 %); at the measured grid 41 / 28 ms per generated event:
3.7 ks / 46 ks CPU per kept event → 0.5 M ≈ 58 / 730 CPU-yr. Yield per kept event in the top 3 Step-3 cells
(truth pairs 58/90 resp. 9/19 × ε_reco 0.8) ≈ 0.5 / 0.4 → 0.5 M at 125–300 ≈ 190 k / 58 k / 9 k SS pairs in
52–74 / 74–105 / 105–150 (vs today's projected 41 k / 12 k / 0.8 k, m<10 summed over η). Combination with
the pT^μμ>20 sample needs an event-level overlap treatment on the filter variable (same as already needed
between the original ≥2μ, SI and SS samples) and the new DSID's own AMI σ·ε_filt.

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
- For the production request (user): JO names (D3 scheme) and the SS pTH24_40 no-pT^μμ-cut variant (D6);
  the GeneratorFilters MR (branch `dimuon-same-sign-filter`, commit `4b66ba5e232` on `main`); set
  `evgenConfig.nEventsPerJob` per config from R6 (e.g. SS 40–70 ≈ 25, SS 24–40 ≈ 100 for ~50 h jobs);
  quote Pythia 8.317 (23.6.68) vs 8.308 (23.6.47) of the existing sample.
- Optional: boost projection in the 20–52 GeV pair-pT cells (where SS 24–40 / 40–70 contribute) —
  outside the 4 reference tables, not requested.

## Latest Stage
DONE (2026-09-29; follow-ups 2026-10-02/03/04 DONE: D10 4 tables; production bundle of 7 configs; R8/R9; D11
tables for the final request incl. SS50; SS/OS slope helper plot). Open for the user: submit the MR; send the bundle.
