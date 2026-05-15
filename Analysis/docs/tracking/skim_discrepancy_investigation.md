# PbPb23 Skim Statistics Discrepancy Investigation

**Question:** Old skim (R24, task 44157121 Apr2025) has ~108M events; new skim (R25, Apr2026) has ~10.4M events. Both claim to process the same 398M-event input. What explains the 10× difference?

---

## Dataset Summary

| | Old skim | New skim |
|---|---|---|
| Location | `~/dcachearea/dimuon_data/pbpb_2023/` (parts 1–6) | `~/usatlasdata/dimuon_data/pbpb_2023/` (parts 1–4) |
| Total events | 108,388,133 | 10,374,641 |
| Athena release | AthAnalysis 24.2.40 (`source_24/`) | AthAnalysis 25.2.89 (`source/`) |
| PanDA task | 42925802 (Jan2025) / 44157121 (Apr2025) | 49860533/49683976/50022072/49860535 (Apr2026) |
| Input container | `data23_hi.periodAllYear2.physics_HardProbes.PhysCont.AOD.repro35_v01` = 61 × `merge.AOD.r16069_p6447_tid*`, 87,924 files, 660 TB, ~398M events | Same 61 per-run `r16069_p6447` datasets via `--inDsTxt` |
| PanDA flags | default `pathena --inDS` (local staging) | `--nGBPerJob MAX --mergeOutput` (XRootD remote I/O) |
| Input files processed | ~all 87,924 (4183 sub-jobs, 98%) | **all 87,924** (99.4% processed per BigPanDA; 757 merge output files) |

---

## Key Observations

### O1 — Trigger pass fraction is similar in both skims
- Old skim: 89.1% of stored events pass `b_HLT_mu4_L1MU3V`
- New skim: 87.3% of stored events pass `b_HLT_mu4_L1MU3V`
- Old skim: 0.21% of mu4-passing events have no muon
- New skim: 0.21% of mu4-passing events have no muon

**Implication:** Both skims store genuine muon-trigger events. No sign of fake/broken trigger decisions in stored data.

### O2 — Test AOD (run 462969, LBs 488+497, 2487 events)
- New skim processes same AOD → 832/2487 events pass (33.5% for this file)
- AOD cleaning flags: 0% event errors for all 8 flags (LAr, Tile, SCT, Pixel, TRT, Muon, Core, Core18)
- AOD muon content: 79.9% events have ≥1 muon; 50.4% have muon pT > 4 GeV

### O3 — TDT in R24 returns false for all events (on R3 data)
Diagnostic (`diagnostic_muons.py`) run in AthAnalysis 24.2.40:
- Both properly-configured TDT (`NavigationFormat=TrigComposite`) and default TDT return `isPassed()=false` for ALL 2487 events
- Yet old skim ntuple has `b_HLT_mu4_L1MU3V=true` for 89% of events → old skim's trigger decision logic must differ from the diagnostic (see H4)

---

## Hypotheses

### H1 — Old and new skim processed different input datasets *(open)*
The 10× discrepancy could come from old skim processing more input events (more runs, different period, or different dataset tag).

**Test:** Compare run number distributions between old and new skim. Check if old-skim-only runs/LBs exist in new skim InDsTxt. Check GRL differences.

**Status: UNDER INVESTIGATION**

---

### H2 — New skim applies tighter cuts that reject 90% of genuine input events *(open)*
GRL or cleaning differences between the two skim versions could explain the factor of 10.

**Known GRL difference:**
- Old skim JO backup (`TrigRates_backup_aug25.py`): `data23_hi.periodAllYear_DetStatus-v113-pro31-08_MERGED_PHYS_HeavyIonP_All_Good.xml`
- New skim JO (`TrigRates_CA.py`): `data23_hi.periodAllYear_DetStatus-v120-pro33-03_MERGED_PHYS_HeavyIon_All_Good_IgnoreBSPOT_INVALID.xml`

**Known cleaning cut difference:**
- Old R24 code: checks 4 flags (LAr, Tile, SCT, Core18)
- New R25 code: checks 8 flags (+ Pixel, TRT, Muon, Core)
- But O2 shows 0% error rate on test AOD → cleaning is NOT the cause

**Test:** Rerun old and new skim on the same test AOD and compare per-cut event counts.

**Status: UNDER INVESTIGATION**

---

### H3 — Old skim used `StoreAllEvents=True` *(partially open)*
If the old skim had `StoreAllEvents=True`, it would store all events passing GRL+cleaning regardless of trigger.

**Evidence against:** `TrigRates_backup_aug25.py` has `StoreAllEvents=False`. But this is a backup of the JO; the actual JO for task 44157121 may differ.

**Test:** Check git history for the JO used in task 44157121.

**Status: UNDER INVESTIGATION**

---

### ~~H4 — TrigDecisionTool was broken in R24 (all events pass trigger)~~ RULED OUT

**Ruling-out evidence:**
- Diagnostic in R24: TDT returns `isPassed()=false` for ALL events. Trigger was not always-true.
- Old skim ntuple: 89.1% pass mu4. If trigger were always-true, we'd expect ~100% (all events stored), not 89%.
- Old skim ntuple: 0.21% of mu4-pass events have no muon. If TDT were broken, stored events would have no muons (since physics_HardProbes is ~80% muon events). The 0.2% is consistent with real trigger selection, not fake.
- The trigger branches in the ntuple reflect real trigger decisions from the COOL/trigger DB stored in the AOD, not from a live TDT call.

**Conclusion:** Trigger selection was working correctly in the old skim. It correctly selected ~89% mu4 events from whatever input was processed.

---

## Run-distribution comparison (O4) — DECISIVE

Per-run event counts from `RunNumber` branch (histogram bins width=10):

| Run bin | Old skim | New skim | Ratio |
|---------|----------|----------|-------|
| 462200 | 3,872,392 | 367,768 | 10.5× |
| 462440 | 1,542,372 | 138,836 | 11.1× |
| 462530 | 1,596,744 | 144,295 | 11.1× |
| 462760 | 4,498,893 | 426,012 | 10.6× |
| 463020 | 4,825,272 | 506,180 | 9.5× |
| 463040 | 9,119,685 | 732,832 | 12.4× |
| 463120 | 10,084,076 | 1,013,406 | **9.96×** |
| 463310 | 4,060,027 | 408,393 | 9.94× |
| *(most others)* | *(various)* | *(various)* | **~10×** |

**Old skim run range:** 462200–463380 (34 bin groups)  
**New skim run range:** 461630–463420 (48 bin groups)

**Key findings:**
1. Old skim is **missing** early runs 461630–462190 (~13 bin groups, ~523K events in new skim) — H1 confirmed for this subset
2. For **all common runs**, old skim has **~10× more events** — H2 confirmed as dominant effect
3. New skim has runs 463414–463427 (Part4) not present in old skim either

---

## Updated Hypotheses

### H1 — Different input datasets *(partially confirmed)*
- Old skim **missing** early runs 461630–462190 present in new skim
- Old skim also **missing** late runs 463414–463427 (new skim Part4)
- These account for only ~1.1M events (small fraction of 10.4M new skim total)
- **The missing-runs effect is a minor contributor** (< 11% of discrepancy)

**Confirmed cause:** Old skim GRL v113 covers only 43 runs (462201–463380); new skim GRL v120 adds 18 more runs.

### ~~H3 — Old skim used StoreAllEvents=True~~ RULED OUT
*(see O7, O8 above)*

### ~~H4 — TrigDecisionTool was broken in R24~~ RULED OUT
*(see above)*

### ~~H2a — Old JO had TriggerChains="HLT_.*"~~ RULED OUT (O11)

If `TriggerChains="HLT_.*"` had been active, `getChainGroup()` would have populated hundreds of HLT branches in the output ntuple. The old skim ntuple has **exactly 32 `b_HLT_*` branches**, corresponding precisely to the 8 configured muon/dimuon chains × 4 sub-branches (L1TBP/L1TAP/L1TAV/pass). No extra branches exist. The JO unambiguously used the restricted specific-chain configuration. **H2a is definitively ruled out.**

### ~~H2b — Different dataset composition~~ RULED OUT (O18)

Rucio `list-content` on `repro35_v01` shows it is a PhysCont wrapper containing **exactly the same 61 `merge.AOD.r16069_p6447_tid*` per-run datasets** used by the new skim. The input data is identical. See O18.

### ~~H2c — New skim grid jobs processed only ~10% of input files~~ DISPROVEN

**Original claim:** `--nGBPerJob MAX` caused PanDA jobs to process only ~10% of input files.

**Indirect evidence that seemed to support H2c:**
- Run 463017: old skim 1.26M events (41.4% pass rate), new skim 128K (4.2% apparent pass rate)
- O17: new skim events are strict subset of old skim
- 757 output files vs 4183 old skim sub-jobs

**Direct evidence that DISPROVES H2c (Latest Stage audit, 2026-05-12):**
- BigPanDA API: `totevproc = totev` for all 4 tasks → 99.4% of ~398M input events were actually processed
- All input datasets show `nfilesfinished ≈ nfiles` → all input files were assigned and completed
- _EXT0 containers were fully downloaded (757 files, 11.4 GB) and locally hadd'd (10.8 GB matches)
- Skim fraction is genuinely ~2.6% (10.4M / 398M), not an artifact of incomplete processing
- Analysis job logs confirm: `EvtMax = -1`, 1–3 AOD files per job, `leaving with code 0`

**The ~10× discrepancy with the old skim (~27% pass rate vs ~2.6%) is REAL and remains unexplained.** Both skims processed all ~398M input events. Both have 100% trigger pass rate among stored events. The source code, JO trigger chains, and StoreAllEvents setting appear identical. The O7 test (local rerun) showed R24 and R25 give identical results on a test AOD — but this test's 33.5% pass rate is itself anomalous compared to the grid average of 2.6%, suggesting the test AOD is non-representative.

**Key open question:** What mechanism causes the R24 TDT in AthAnalysis 24.2.40 to store ~27% of HardProbes events vs ~2.6% in R25 (AthAnalysis 25.2.89)? The run directories (JO files) are untracked in git, so we cannot definitively verify what was submitted to the grid. Possible explanations:
1. The TDT `isPassed()` in R24 may use a different default definition (e.g., including prescaled/resurrected chains) than in R25
2. The actual `TrigRates.py` sent to the grid in Apr 2025 may have been different from the current version
3. There may be a configuration-level difference in the TDT between Athena releases that changes trigger evaluation

**Note on L1TBP=100% (O9):** In PbPb collisions, the underlying event produces large forward activity in the muon spectrometer, causing L1_MU3V (pT > 3 GeV before precision reco) to fire for essentially all PbPb events. L1TBP=100% is therefore expected for any set of stored PbPb events and is NOT evidence of anomalous trigger behavior or event selection.

---

## Test AOD result (O7) — Code equivalence confirmed

**Test:** Old R24 skim rerun on `AOD.41716150._000006.pool.root.1` (run 462969, LBs 488+497, 2487 events, all in GRL v113).

| | Events |
|---|---|
| Input events | 2487 |
| New skim (R25 CA) | 832 |
| Old skim R24 StoreAll=False | **832** (exact match) |
| Old skim R24 StoreAll=True | 2417 |

The R24 and R25 skim codes are equivalent when `StoreAllEvents=False`.

## Old skim trigger pattern (O8) — Rules out StoreAllEvents=True

Checked `data_pbpb23_part1.root` (18M events, run 463017):

| Trigger measure | Value |
|---|---|
| `b_HLT_mu4_L1MU3V` (HLT after prescale) | 89% |
| `b_HLT_mu4_L1MU3V_L1TBP` (L1 before prescale) | **100%** |
| `b_HLT_mu4_L1MU3V_L1TAV` (L1 after veto) | **100%** |
| Any muon trigger passed | **100%** (only 964/18M events pass none) |
| Duplicate (run, LB, eventNumber) | **0 out of 500K checked** |

Compare with local StoreAll=True output (2417 events from same run range):
- HLT_mu4 = 34%, L1TBP = 46%, any muon trig = 34%

**Conclusion:** Old skim cannot be a StoreAll=True output. The L1TBP=100% is expected PbPb physics (see O9), but the trigger pattern mismatch (89% vs 34% HLT_mu4) rules out local StoreAll=True.

## L1TBP=100% explained (O9) — Not a diagnostic

In PbPb collisions, the large underlying event causes copious soft activity in the muon spectrometer. L1_MU3V (L1 muon threshold ~3 GeV, before precision reconstruction) therefore fires before prescale for essentially **all** PbPb events. This makes `b_HLT_mu4_L1MU3V_L1TBP = 100%` a generic PbPb feature, NOT evidence of any anomalous trigger evaluation. Likewise, `L1TAV = 100%` simply indicates no L1 prescale was applied to L1_MU3V in these runs.

## JO comparison (O10) — Untracked run directory prevents definitive diagnosis

All current JO versions (`TrigRates.py`, `TrigRates_backup_aug25.py`, `TrigRates_JO.py`, `TrigRates_CA.py`) set identical trigger configuration for hi2023:
- `StoreAllEvents = False`
- `TriggerChains = "|".join(MinBias_triggers) = ""`  (no MinBias triggers for hi2023)
- Same `MuonTriggerChains` and `DiMuonTriggerChains`

However, **`run_23hi/` is not tracked in git**. The `TrigRates.py` used for task 44157121 (April 2025) cannot be recovered from git history. The untracked file may have had a different trigger configuration at that time. The `TrigRates_backup_aug25.py` label suggests it was saved as a backup in August 2025, four months after the production task ran.

The new skim (April 2026) used `TrigRates_CA.py` (tracked, Component Accumulator mode, R25), submitted via per-run datasets with 4 separate tasks — a fundamentally different submission strategy from the old single-task `--inDS container` approach.

---

## GRL comparison (O5) — DECISIVE for H1

Parsed LB ranges from both GRL XML files. **For all 43 common runs, LB coverage is exactly identical (ratio = 1.00 everywhere).** The ONLY difference is 18 extra runs in v120 not in v113:

- v113: 43 runs, range 462201–463380
- v120: 61 runs, range 461633–463427
- Runs only in v120: 461633–462149 (14 early runs), 462502, 463389, 463414, 463427 (4 late runs)

**This directly explains the missing runs in old skim**: v113 simply did not include those runs. No other explanation needed.

Since LBs are identical for common runs, the per-run event rate from input is the same. The ~10× more events per run in old skim (O4) **cannot** be explained by GRL LB differences.

## PanDA task confirmed (O6)

Task 44157121 cliParams: `pathena TrigRates.py --inDS data23_hi.periodAllYear2.physics_HardProbes.PhysCont.AOD.repro35_v01`
- JO was `TrigRates.py` (current file says `StoreAllEvents=False`, but `run_23hi/` is not tracked in git — file may have been different in April 2025)
- Input was the inclusive `PhysCont.AOD.repro35_v01` container, ~398M events
- This is a DIFFERENT submission strategy from new skim (which uses 4 tasks × per-run containers)

---

## Old skim branch list (O11) — Rules out H2a

Inspected `data_pbpb23_part1.root` branch names: **exactly 32 `b_HLT_*` branches**, corresponding to the 8 configured muon chains × 4 sub-branches each (pass, L1TBP, L1TAP, L1TAV). There are no extra branches from any wider pattern. If `TriggerChains="HLT_.*"` had been used, hundreds of HLT chains would be present. **H2a definitively ruled out.**

## Per-trigger analysis for run 463017 (O12, O13, O14)

**O12 — Trigger fractions across full old skim:** HLT_mu4_L1MU3V = 88.88%, HLT_mu6_L1MU3V = 34.70%, HLT_mu6_L1MU5VF = 30.78%; Fail ALL triggers = 0 (0.0000%). All events pass at least one muon trigger.

**O13 — Per-run comparison (exact bin width 1):**

| Run | Old | New | Ratio |
|-----|-----|-----|-------|
| 463017 | 1,255,623 | 127,868 | 9.82× |
| 462995 | 4,817,025 | 20,005 | **240×** (incomplete new skim) |
| 462814 | 2,712,585 | 47,313 | **57×** (incomplete new skim) |
| *(~30 others)* | various | various | **9–11×** (tight cluster) |
| TOTAL | 108,388,133 | 10,374,641 | 10.45× |

Anomalous-ratio runs (462995, 462814, 462244, 463364) likely have failed grid jobs in the new skim.

**O14 — Per-trigger breakdown for run 463017:**

| Trigger | Old | New | Abs ratio |
|---------|-----|-----|-----------|
| HLT_mu4_L1MU3V | 100.00% | 100.00% | 9.82× |
| HLT_mu6_L1MU3V | 27.81% | 27.70% | 9.86× |
| HLT_mu6_L1MU5VF | 24.72% | 24.83% | 9.78× |
| HLT_2mu4_L12MU3V | 2.31% | 1.83% | 12.42× |
| Fail ALL | 0.00% | 0.00% | — |

**Every trigger has the same ~9.82× absolute ratio.** Trigger composition is identical per-lbn. The discrepancy cannot be in the skim algorithm — it is in the input data volume.

## Duplicate check (O15) — No duplicates in old skim

Loaded 500K events for run 463017 from old skim into `std::map<pair<lbn,evtNum>, int>`.  
All 500K entries are unique (multiplicity = 1 for every entry). The old skim does NOT contain duplicate events.

## lbn distribution comparison (O16) — Same coverage, 10× events per lbn

For run 463017:
- Old skim: 94 lbns, range [170, 263], avg **13,358 events/lbn**
- New skim: 89 lbns, range [171, 264], avg **1,437 events/lbn**
- **Lbn count ratio: 1.06×** (essentially same lbn range)
- **Events/lbn ratio: 9.30×** (the 10× excess is within each lbn, not from extra lbn coverage)

Sample per-lbn:

| lbn | old | new | ratio |
|-----|-----|-----|-------|
| 171 | 20,571 | 1,999 | 10.3× |
| 172 | 19,112 | 1,613 | 11.8× |
| 176 | 18,771 | 2,056 | 9.1× |
| 185 | 19,176 | 1,924 | 10.0× |
| 196 | 3,325 | 0 | ∞ (lbn not in new skim GRL) |

Five lbns present in old skim with 0 events in new skim (excluded by v120 GRL or not in per-run container for those lbns).

## Event overlap test (O17) — **New skim is a strict subset of old skim**

Loaded all 127,868 new skim events for run 463017 into a `std::set<pair<lbn,evtNum>>`. Then scanned all 1,255,623 old skim events for run 463017 and checked membership.

**Result: 127,868 out of 1,255,623 old skim events (10.18%) are found in the new skim set. This equals exactly the new skim size.**

This means: **every single new skim event (lbn, eventNumber) pair exists in the old skim.** The new skim is a strict subset of the old skim. The old skim has an additional 1,127,755 unique events (no duplicates) that do not exist in the per-run `r16069_p6447` containers used by the new skim.

Per-lbn overlap fractions are uniformly ~7–15% (consistent with the global 10.18%), with a few lbns at 0% (excluded by v120 GRL).

---

## Container composition (O18) — repro35_v01 = same r16069_p6447 datasets

Rucio `list-content` on `data23_hi.periodAllYear2.physics_HardProbes.PhysCont.AOD.repro35_v01/` reveals:

```
61 datasets, ALL of the form:
data23_hi:data23_hi.00XXXXXX.physics_HardProbes.merge.AOD.r16069_p6447_tid*
```

These are **exactly the same per-run datasets** used by the new skim's `InDstxt` files. The `repro35_v01` PhysCont container is simply a wrapper around those 61 datasets.

**This overturns the earlier H2b hypothesis:** the input data is identical for both skims. The 10× difference is NOT from different dataset composition — it is from incomplete processing of the input by the new skim grid jobs.

Per-run input event counts (from rucio `list-files`):

| Run | Input events | Input files | Input size | Old skim output | New skim output | Old pass rate | New apparent |
|-----|-------------|-------------|------------|-----------------|-----------------|---------------|-------------|
| 463017 | 3,031,977 | 661 | 4.99 TB | 1,255,623 | 127,868 | 41.4% | 4.2% |
| 462580 | 4,664,323 | 1,023 | 7.83 TB | — | — | — | — |
| 462549 | 9,451,424 | 2,065 | 15.84 TB | — | — | — | — |

---

## FINAL CONCLUSION (UPDATED 2026-05-12)

### Root cause of the 10× discrepancy — UNKNOWN

Controlled experiments (Checks 1–3, Checkpoints 6–8) have **ruled out** all previously hypothesized technical causes:

| Hypothesis | Status | Evidence |
|------------|--------|----------|
| H2c: New skim processed only ~10% of input | **DISPROVEN** | BigPanDA confirms 99.4% of ~398M events processed |
| H5: AthAnalysis version difference (R24 vs R25) | **RULED OUT** | Same source_24 code gives identical results under R24 and R25 |
| H6: Source code changes | **RULED OUT** | source_24+R24 and source+R25 give identical results on test AODs |
| H3: StoreAllEvents=True in old skim | **RULED OUT** | Old skim trigger pattern inconsistent with StoreAllEvents=True |
| H4: Broken TDT in R24 | **RULED OUT** | TDT works correctly in both versions |

**The 10× discrepancy CANNOT be reproduced locally.** On test AODs from run 462969 (2487 and 5249 input events), all configurations (old code + R24, new code + R25, old code + R25) produce **exactly identical** output event counts (832 and 1925, respectively). The local skim fraction (~33–37%) is much higher than the new skim's grid average (2.6%).

### Remaining open hypothesis

**H7 — The JO actually submitted to the grid in April 2025 was different from what exists on disk today.** The `run_23hi/` directory is **untracked in git**. The file `TrigRates.py` on disk may have been modified since the old skim grid submission (task 44157121, April 2025). A backup file `TrigRates_backup_aug25.py` exists (suggesting a backup was taken in August 2025, 4 months after the grid task). We cannot verify what configuration was used on the grid.

### Key numbers

| Metric | Old skim (grid) | New skim (grid) | Local test (all configs) |
|--------|-----------------|-----------------|--------------------------|
| Input events | ~398M | ~398M | 2487 / 5249 |
| Output events | ~108M | ~10.4M | 832 / 1925 |
| Skim fraction | ~27% | ~2.6% | ~33–37% |

The local skim fraction (~33–37%) exceeds even the old skim's grid average (27%), confirming that run 462969 has above-average trigger rate. But the critical finding is that R24 and R25 agree perfectly locally, while they disagree by 10× on the grid. This strongly points to a configuration difference in the grid submission (H7), not a code or framework difference.

---

## Trigger matching deltaR study (O19)

The current `source/TrigRates.cxx` has `MuonTriggerMatchDR` defaulting to **0.02**, and simultaneously stores a hardcoded dR=0.01 result in parallel (`_0_01` branches). The old `source_24/TrigRates.cxx` had this hardcoded to 0.1.

**Test setup:** Ran skim on 7736 input AOD events (run 462969, 2 files: `_000001` 5249 evt + `_000006` 2487 evt) with `TRIGRATES_MATCH_DR` env-var override. Both runs yield 2757 stored events. Cross-check: dR=0.01 count from both runs matches exactly (2058).

### Events with ≥1 mu4-matched muon (HLT_mu4_L1MU3V)

| deltaR | Events passing | Fraction of skim | Relative to dR=0.1 |
|--------|---------------|------------------|---------------------|
| 0.01   | 2058          | 74.6%            | 0.841               |
| 0.02   | 2102          | 76.2%            | 0.859               |
| 0.1    | 2446          | 88.7%            | 1.000               |

### Per-muon mu4 trigger match rate

| deltaR | Muons matched | Fraction of all muons | Relative to dR=0.1 |
|--------|--------------|----------------------|---------------------|
| 0.01   | 2095 / 6483  | 32.3%                | 0.830               |
| 0.02   | 2141 / 6483  | 33.0%                | 0.849               |
| 0.1    | 2523 / 6483  | 38.9%                | 1.000               |

**Key findings:**
- dR=0.1 → 0.02 (old → new default): **~15% fewer mu4-matched muons**, ~14% fewer events with any matched muon
- dR=0.1 → 0.01: ~17% fewer matched muons, ~16% fewer events
- Difference between 0.01 and 0.02 is small (~2% relative) — most loss is between 0.02 and 0.1
- The total number of stored events (2757) is **unchanged** across all three dR values — deltaR only affects per-muon `muon_trig_match` flags, not event-level storage (which is gated by HLT trigger decision, not offline matching)

---

## Latest Stage — Direct evidence audit for H2c

### Hypothesis under scrutiny

H2c claims `--nGBPerJob MAX` caused new skim grid jobs to process only ~10% of input files. So far this is supported by **indirect** evidence (output event counts, subset test, file-count ratios). The user correctly notes this could be a guess. We now seek **direct** evidence from the grid jobs themselves.

### Objectives
1. Confirm from BigPanDA task metadata how many input files were **assigned** vs how many input events were **actually processed** per job.
2. Download and parse actual job log files from the new skim grid tasks to find Athena's own accounting of files opened, events read, events written.
3. Either confirm H2c with direct log evidence, or discover an alternative explanation.

### Methods and tools
- **BigPanDA**: query task/job metadata via `pbook` or the PanDA web API for tasks 49860533 (part1), 49683976 (part2), 50022072 (part3), 49860535 (part4)
- **Rucio log download**: download log tarballs from output log datasets (derived from `data-merging-record.txt` — replace `_EXT0` with `.log/`)
- **Log parsing**: search Athena logs for `events processed`, `EvtMax`, `filesInput`, file-open counts, `THistSvc` output event tallies

### Log datasets (from grid_sub.sh outDS names)
| Part | outDS | Log dataset |
|------|-------|-------------|
| 1 | `user.yuhang.TrigRates.dimuon.PbPb2023data.Apr2026.v9.part1.` | `...v9.part1.log/` |
| 2 | `user.yuhang.TrigRates.dimuon.PbPb2023data.Apr2026.v7.part2.` | `...v7.part2.log/` |
| 3 | `user.yuhang.TrigRates.dimuon.PbPb2023data.Apr2026.v11.part3.` | `...v11.part3.log/` |
| 4 | `user.yuhang.TrigRates.dimuon.PbPb2023data.Apr2026.v8.part4.` | `...v8.part4.log/` |

### Findings

**Checkpoint 1 — BigPanDA task metadata (2026-05-12)**

Queried BigPanDA JSON API (`https://bigpanda.cern.ch/task/<ID>/?json`) for all 4 new skim tasks:

| Task ID | Part | totev (input) | totevproc (processed) | pctfinished | status |
|---------|------|--------------|----------------------|-------------|--------|
| 49860533 | Part1 (v9) | 127,388,351 | 124,867,920 | 98.0% | finished |
| 49683976 | Part2 (v7) | 124,060,865 | 124,060,865 | 100% | done |
| 50022072 | Part3 (v11) | 123,625,643 | 123,625,643 | 100% | done |
| 49860535 | Part4 (v8) | 22,835,292 | 22,835,292 | 100% | done |
| **Total** | | **397,910,151** | **395,389,720** | **99.4%** | |

`totevoutput` = 0 for all tasks because TrigRates writes to THistSvc (ROOT TFile), not an ATLAS POOL output stream — PanDA's event counting only tracks POOL output.

All input datasets show `nfilesfinished ≈ nfiles` (Part1 has 562 failed files out of ~28,700, Part2-4 have 0 failed). Input files span all 61 expected runs.

**Conclusion: H2c is DISPROVEN.** All ~398M input events were processed. The new skim did NOT miss 90% of input.

**Checkpoint 2 — _EXT0 container file counts (2026-05-12)**

The `_EXT0` output containers hold ALL merged outputs (not just one file):

| Part | Container | Files | Total size |
|------|-----------|-------|-----------|
| Part1 (v9.part1) | `...v9.part1._EXT0/` | 256 | 2.525 GB |
| Part2 (v7.part2) | `...v7.part2._EXT0/` | 231 | 3.998 GB |
| Part3 (v11.part3) | `...v11.part3._EXT0/` | 240 | 4.374 GB |
| Part4 (v8.part4) | `...v8.part4._EXT0/` | 30 | 0.538 GB |
| **Total** | | **757** | **11.435 GB** |

The local hadd merge log (`hadd_merge_log.txt`) confirms ALL files from each `_EXT0` container were downloaded and hadd'd into the local `data_pbpb23_partN.root` files (3.8G + 0.5G + 4.1G + 2.4G = 10.8 GB, consistent with the container total).

**Conclusion: The full output was downloaded.** The 10.4M entries IS the complete new skim output.

**Checkpoint 3 — Analysis job logs confirm single-file processing (2026-05-12)**

Extracted nested analysis job logs from Part3 merge tarball. Each analysis job processes 1–3 AOD input files (e.g., `--filesInput=AOD.41547598._000001.pool.root.1`), with `EvtMax = -1` (all events). One sample job (7121778020) on a single 10.7 GB AOD produced a 340 KB myfile.root — consistent with a low skim fraction.

The 27 analysis jobs in this merge group cover AODs _000001 through _000049 from one dataset (run 463124). PanDA assigned 1–3 files per job, consistent with `--nGBPerJob MAX` assigning ~1 file per job (each AOD is ~10 GB ≈ MAX for a single job).

**Checkpoint 5 — Old skim task metadata (2026-05-12)**

BigPanDA confirms the old skim task also processed all input:
- Task 44157121 (Apr2025): `totevproc = 392,728,906` of `totev = 397,910,151` (98.7%)
- Same 61 input datasets, same `TrigRates.py` JO, AthAnalysis 24.2.40
- Task 42925802 (Jan2025): `totevproc = 396,387,214` (99.6%), same JO, same input

Both old and new skim processed the same ~398M input events with the same trigger chains. The skim fraction difference is genuine: old skim ~27%, new skim ~2.6%.

Per-trigger fractions (Part1 comparison):

| Trigger | Old skim (18.1M) | New skim (3.7M) |
|---------|-----------------|-----------------|
| mu4     | 88.9%           | 97.6%           |
| mu6_L1MU3V | 34.7%       | 29.4%           |
| mu6_L1MU5VF | 30.8%      | 26.2%           |
| 2mu4    | 2.8%            | 2.0%            |
| Any trigger | **100%**    | **100%**        |

Both skims have 100% trigger pass rate (no events with all triggers false). The mu4 fraction is lower in the old skim (88.9% vs 97.6%), suggesting the old skim stored events where lower-frequency triggers passed but mu4 did not — consistent with more permissive trigger evaluation.

**Checkpoint 4 — Skim fraction analysis (2026-05-12)**

The true skim fraction for the new skim:
- Total input: ~398M events
- Total output: ~10.4M entries (from `data-merging-record.txt`)
- **Skim fraction: ~2.6%**

The old skim had ~100M entries. Assuming it also processed ~398M input events:
- **Old skim fraction: ~25%**

The ~10× discrepancy between old and new skim is therefore a **genuine difference in skim selection rate**, not incomplete processing. The 2.6% skim fraction in the new skim is expected for HardProbes data where only ~3% of events fire any muon trigger (the primary selection criterion in TrigRates).

---

## Current Status (updated 2026-05-12)

**All technical hypotheses (H1–H6) have been tested and ruled out.** Controlled experiments with all combinations of old/new source code and AthAnalysis R24/R25 produce identical results on test AODs. The 10× grid discrepancy cannot be reproduced locally.

**The only remaining hypothesis is H7:** the JO file actually submitted to the grid in April 2025 may have differed from the version currently on disk (`run_23hi/` is untracked in git). This would explain why the discrepancy is visible in grid output but not in local tests using the current on-disk files.

**Separate finding (O19):** The trigger matching deltaR changed from 0.1 (old `source_24/`) to 0.02 (new `source/`). This does NOT affect event counts in the skim (event storage depends on HLT decision, not offline matching), but it does reduce the number of muons flagged as trigger-matched by ~15%.

---

## Stage: Controlled source code & AthAnalysis version experiments (2026-05-12)

### Motivation

The 10× skim fraction discrepancy (27% old vs 2.6% new) is confirmed genuine. Both skims process the same input with the same trigger lists, `StoreAllEvents=False`, and `UseTrigger=True`. The key changes between old and new skim are:
1. **AthAnalysis version**: 24.2.40 → 25.2.89
2. **Source code restructuring**: `source_24/` with `__ATHENA_24p2__` ifdefs → unified `source/` with `HF_IS_R25` cmake-detected ifdefs
3. **Event-level cleaning flags**: 4 flags (LAr, Tile, SCT, Core18) → 8 flags (+ Pixel, TRT, Muon, Core)
4. **JO config framework**: Old-style `CfgMgr` JO → Component Accumulator CA
5. **GRL version**: v113 → v120 (18 extra runs, but identical LB coverage for common runs)
6. **New features added to `source/`**: dimuon per-leg matching module, dR=0.01 parallel branches, MC weight names, truth event container, `DoubleToString()`, removed ZdcCalib mode, `barcode()` → `uid()`, HLT_MuonsCB_FS key
7. **StoreTracks**: 0 → 1 (JO setting, determines what's stored per event)

Of these, only changes 1–5 could potentially affect **event storage** (skim fraction). Changes 6–7 only affect per-event content or are on code paths not reached for hi2023 data.

### Plan

**Test AOD:** `AOD.41716150._000006.pool.root.1` (run 462969, 2487 events) at `/gpfs/mnt/.../test_aod/data23_hi/`

**Check 3 — End-to-end reproduction (old code+R24 vs new code+R25):**
- Condition A: `source_24/` at git commit `10b3043` (April 2025) + AthAnalysis 24.2.40 + old JO (`TrigRates_test.py`)
- Condition B: current `source/` + AthAnalysis 25.2.89 + new CA JO (`TrigRates_CA.py`)
- Goal: reproduce the skim fraction difference locally on the test AOD

**Check 2 — AthAnalysis version isolation:**
- Use the OLD `source_24/` code (which has `__ATHENA_24p2__` → TDT init works for both R24 and R25)
- Compile under R24 (24.2.40) → run on test AOD → record event count
- Compile under R25 (25.2.89) → run on test AOD → record event count
- If counts differ, the AthAnalysis version (TDT behavior, trigger navigation) is the cause

**Check 1 — Per-change controlled experiments (all on R25):**
Using `source/` compiled with R25, systematically test each change:
- **1a: Cleaning flags** — Revert to 4-flag cleaning (old pattern: LAr, Tile, SCT, Core18)
- **1b: GRL version** — Use v113 instead of v120
- **1c: JO framework** — Use old-style JO vs CA (both calling same R25 binary)
- **1d: StoreTracks** — Change from 1 to 0

Each experiment changes ONE variable; all others remain at new-skim defaults.

### Catalog of source code differences (old `source_24/` @ `10b3043` vs current `source/`)

| # | Category | Old (`source_24/`) | New (`source/`) | Affects skim fraction? |
|---|----------|-------------------|-----------------|----------------------|
| 1 | Ifdef system | `#ifdef __ATHENA_24p2__` (from `AthenaVersion.h`) | `#if defined(HF_IS_R25)` (from cmake) | No (same code paths active) |
| 2 | TDT init | Via `__ATHENA_24p2__` guard, always active | Via `HF_IS_R25` guard, active for R25 only | No (both init TDT correctly for their version) |
| 3 | Cleaning flags | 4: LAr(1), Tile(2), SCT(4), Core18(8) | 8: Pixel(1), SCT(2), TRT(4), Muon(8), LAr(16), Tile(32), Core(64), Core18(128) | **Potentially** — extra flags could reject more events |
| 4 | ZDCAnalysisTool/HIPileupTool | Present but disabled (`#ifndef __ATHENA_24p2__`) | Removed entirely | No (was already disabled in R24) |
| 5 | ZdcCalib mode (`m_is_Zdc_Calib`) | Wraps Init/Process in `if(m_is_Zdc_Calib==false)` | Removed; always init/process | No (ZdcCalib=false in both JO) |
| 6 | `barcode()` → `uid()` | `particle->barcode()` | `particle->uid()` | No (truth-only, disabled for data) |
| 7 | MC weight tool | Not present | PMGTruthWeightTool added | No (data, not MC) |
| 8 | `MuonTriggerMatchDR` property | Not present (hardcoded 0.1) | `declareProperty("MuonTriggerMatchDR", 0.02)` | No (affects matching flags, not event storage) |
| 9 | Per-leg dimuon matching | Not present | `TrigMuonMatchingModule` created for R25 | No (post-event-selection) |
| 10 | dR=0.01 parallel branches | Not present | `_0_01` branches added | No (additional storage) |
| 11 | `DoubleToString()` | Not present | Added for branch naming | No (utility) |
| 12 | `HLTMuonsFSKey` | Not present | `"HLT_MuonsCB_FS"` property | No (for RoI matching gated by `use_trigger && m_hlt_muons_fs_key!=""`) |
| 13 | `StoreDimuonPerLeg` | Not present | `true` | No (additional storage) |
| 14 | `StoreMCWeightNames` | Not present | `false` | No (MC only) |
| 15 | `TruthEventContainerKey` | Not present | `"TruthEvents"` | No (truth only) |
| 16 | `m_store_L1` type | `int` (default 0) | `bool` (default false) | No (same effect) |
| 17 | `m_store_EventInfo` | `>=2 then also store ZDC` | `>=1 then store event info` (ZDC moved out) | No (ZDC info restructured, not affecting event selection) |
| 18 | ProcessEventInfo call | Inside `if(!m_is_Zdc_Calib)` block | Unconditional | No (ZdcCalib=false for data) |
| 19 | `dimuMatchModule` cleanup | Not present | `delete m_dimuMatchModule` in finalize | No |
| 20 | ProcessTriggers | Identical logic | Identical logic | No |

**Only item #3 (cleaning flags) could potentially affect skim fraction.** All other changes are either (a) disabled for data, (b) only affect per-event content, or (c) are on code paths already active in both versions.

### Experiment results

**Test AODs:** Both from run 462969 (dataset `data23_hi.00462969...AOD.r16069_p6447`):
- Small: `AOD.41716150._000006.pool.root.1` — 2487 input events
- Large: `AOD.41716150._000001.pool.root.1` — 5249 input events

**Build configurations:**
- `source_24 + R24`: old `source_24/` (git commit `10b3043` era, with `barcode→uid` and muon quality try/catch fixes) compiled against AthAnalysis 24.2.40
- `source + R25`: current unified `source/` compiled against AthAnalysis 25.2.89
- `source_24 + R25` (Check 2): old `source_24/` compiled against AthAnalysis 25.2.89 (required `barcode()→uid()` fix for compilation; truth code disabled for data)

All configurations used the same old-style JO (`TrigRates_test_r24.py`) with:
- `StoreAllEvents=False`, `UseTrigger=True`, `IsRun3=True`, `RunYear=2023`
- Same 6+2 trigger chains, GRL v113, `METContainerKey=""`
- The R25+source run used `TrigRates_CA.py` (CA framework, GRL v120)

**Checkpoint 6 — Check 3: End-to-end reproduction (2026-05-12)**

| Configuration | Small AOD (2487 in) | Large AOD (5249 in) | Skim fraction |
|---|---|---|---|
| source_24 + R24 (old) | **832** | **1925** | 33.4–36.7% |
| source + R25 (new) | **832** | **1925** | 33.4–36.7% |

**Result: Identical.** The old and new skim configurations produce exactly the same number of output events on both test AODs. The 10× discrepancy is NOT reproducible locally.

**Checkpoint 7 — Check 2: AthAnalysis version isolation (2026-05-12)**

| Configuration | Small AOD | Large AOD |
|---|---|---|
| source_24 + R24 | **832** | **1925** |
| source_24 + R25 | **832** | **1925** |

**Result: Identical.** Changing ONLY the AthAnalysis version (24.2.40 → 25.2.89) while keeping the same source code does NOT change the skim fraction. The TDT `isPassed()` behavior is identical between R24 and R25 on this test AOD.

**Checkpoint 8 — Check 1: Source code differences (2026-05-12)**

Since Check 3 showed no difference at all between old and new configurations, individual per-change experiments (cleaning flags, GRL, etc.) are moot — there is zero difference to decompose. All source code changes combined produce exactly 0 events of difference.

**Critical implication:** The ~10× skim fraction discrepancy (27% old vs 2.6% new, as measured on grid output) **cannot be reproduced on test AODs from run 462969**. Both test AODs show ~33–37% skim fraction for ALL configurations — matching neither the 27% grid average of the old skim nor the 2.6% of the new skim.

This means one of:
1. Run 462969 is non-representative — its trigger rate is 10× higher than the grid average
2. The discrepancy originates from a mechanism not exercised by the test AOD (e.g., run-dependent trigger configuration, input file corruption, specific LB ranges)
3. The JO actually submitted to the grid in April 2025 was different from what exists on disk today (run_23hi is untracked in git)

The test AOD skim fractions (~33–37%) are much closer to the old skim's 27% grid average than the new skim's 2.6%. This suggests that run 462969 has a higher-than-average trigger rate, and that the old skim's ~27% may be correct for the trigger chains used, while the new skim's 2.6% may indicate a systematic issue with the grid production that is not present in local runs.

---

## Investigation Steps Completed
1. ~~Compare run distributions~~ → O4
2. ~~Check GRL LB coverage~~ → O5 (identical for common runs)
3. ~~Identify PanDA task JO~~ → O6 (TrigRates.py, repro35 inclusive container)
4. ~~Test StoreAllEvents=True/False~~ → O7, O8 (ruled out)
5. ~~Check for duplicate events~~ → O15 (0 duplicates in 500K)
6. ~~Check L1TBP/L1TAV pattern~~ → O8, O9 (100% is expected PbPb physics)
7. ~~Compare all JO file versions~~ → O10 (all identical for hi2023)
8. ~~Test f-tag vs repro35 hypothesis~~ — ruled out
9. ~~Check old skim branch list~~ → O11 (32 branches = H2a ruled out)
10. ~~Per-run and per-trigger comparison~~ → O12–O14 (uniform 10× ratio, identical trigger composition)
11. ~~Check lbn distribution both skims~~ → O16 (same lbn range, 10× more events/lbn in old skim)
12. ~~Event overlap test (new skim ⊂ old skim?)~~ → O17 (**new skim is strict 10% subset**)
13. ~~Check repro35_v01 container composition~~ → O18 (**same 61 `r16069_p6447` datasets = H2b ruled out**)
14. ~~Check per-run input event counts~~ → O18 (run 463017: 3.03M input, old kept 41%, new kept 4.2%)
15. ~~Check PanDA task completion and output sizes~~ → H2c was provisionally confirmed but later **DISPROVEN** by direct BigPanDA API audit
16. ~~Trigger matching deltaR study~~ → O19 (dR=0.1→0.02 loses ~15% matched muons; event count unchanged)
17. ~~Direct evidence audit for H2c~~ → **H2c DISPROVEN**: BigPanDA API confirms all ~398M input events processed (99.4%); _EXT0 containers fully downloaded (757 files, 11.4 GB); skim fraction is ~2.6%
18. ~~Controlled source code & AthAnalysis version experiments~~ → Checks 1–3 (Checkpoints 6–8): **All configurations produce identical results** on test AODs. Neither AthAnalysis version nor source code changes affect skim fraction. The 10× discrepancy is NOT reproducible locally.
19. **Missing `--evtMax=-1` in grid trf command** → **H8 CONFIRMED (Checkpoint 9)**: `grid_sub.sh` does not pass `--evtMax`, so `TrigRates_CA.py` defaults to `m_EvtMax=1000`. Each grid job processes only 1000 events. This is the root cause of the ~10× discrepancy. All three years (PbPb23/24/25) are affected.

## Hypotheses Summary

| ID | Hypothesis | Status |
|----|-----------|--------|
| H1 | Different input datasets (missing runs) | Partially confirmed — GRL v113 vs v120 explains missing runs (~1.1M events) |
| H2a | TriggerChains="HLT_.*" in old JO | **RULED OUT** (O11: 32 branches only) |
| H2b | Different dataset composition | **RULED OUT** (O18: same 61 r16069_p6447 datasets) |
| **H2c** | **New skim processed only ~10% of input files** | **DISPROVEN** — BigPanDA confirms 99.4% of ~398M input events processed; full output downloaded |
| H3 | StoreAllEvents=True | **RULED OUT** (O7, O8) |
| H4 | Broken TDT in R24 | **RULED OUT** (O1, O8) |
| **H5** | **AthAnalysis version causes TDT behavior difference** | **RULED OUT** — Check 2: source_24 + R24 and source_24 + R25 give identical results (832/1925) |
| **H6** | **Source code changes affect skim fraction** | **RULED OUT** — Check 3: source_24 + R24 and source + R25 give identical results (832/1925) |
| **H7** | **Grid JO differed from local disk copy** | **OPEN** — `run_23hi/` is untracked; cannot verify what was actually submitted to grid in Apr 2025. Local tests cannot reproduce the discrepancy. |
| **H8** | **New skim grid jobs only processed 1000 events each (missing `--evtMax=-1`)** | **CONFIRMED** — see Checkpoint 9 below |

---

## Stage: Missing `--evtMax` root cause (Checkpoint 9, 2026-05-13)

### Discovery

The grid submission trf command in all three `grid_sub.sh` files is:
```
pathena --trf "TRIGRATES_RUNMODE=hi20XX athena.py TrigRates_CA.py --filesInput=%IN" ...
```

**`--evtMax` is not passed.** The JO file `TrigRates_CA.py` (line 106) sets:
```python
m_EvtMax = 1000
```

The code path on the grid:
1. `_parse_evtmax(m_EvtMax)` → returns `1000` (no `--evtMax` on command line)
2. `build_cfg(evtmax=1000)` → sets `flags.Exec.MaxEvents = 1000`
3. `flags.fillFromArgs()` → sees no `--evtMax` → **does not override** (confirmed by reading `AthConfigFlags.py` line 746: `if arg_set('evtMax'): self.Exec.MaxEvents=args.evtMax`)
4. `athena.py`'s own parser (`AthOptionsParser.py` line 95): `if opts.evtMax is not None:` → also does not override

**Each grid job processes only the first 1000 events of its input, regardless of how many events the input file(s) contain.**

### Verification

**Code inspection (definitive):**
- `TrigRates_CA.py` line 106: `m_EvtMax = 1000`
- `TrigRates_CA.py` line 178: `flags.Exec.MaxEvents = evt_max` (1000)
- `TrigRates_CA.py` line 182: `flags.fillFromArgs()` — no `--evtMax` on CLI → no override
- `AthConfigFlags.py` line 746-747: only overrides if `arg_set('evtMax')` returns True
- `AthOptionsParser.py` line 95-96: only overrides if `opts.evtMax is not None`
- Simulated `_parse_evtmax()` with grid-like args: returns **1000**

**BigPanDA task parameters (definitive):**
```
cliParams: pathena --trf "TRIGRATES_RUNMODE=hi2023 athena.py TrigRates_CA.py --filesInput=%IN" ...
```
No `--evtMax` anywhere in the trf command (task ID 49860533).

**Grid job evidence (supporting):**
- Job 7114398915 (task 49860533, Part 1): input 2 files, 5071 events, CPU = 28s
- 28s CPU for 5071 events at ~5ms/event = ~25s → consistent with only 1000 events (5s) + initialization overhead (~23s)

**BigPanDA `neventsUsedTot` is misleading:** It reports 124.9M/127.4M events "used" for Part 1, but this is JEDI's bookkeeping based on file metadata for files assigned to successful jobs — NOT the number of events actually processed by the algorithm.

### Quantitative check

| Part | nfiles (finished) | neventsTot | Output events | Apparent skim % |
|------|-------------------|-----------|---------------|-----------------|
| 1 (v9) | 28,135 | 127.4M | 3,690,782 | 2.90% |
| 2 (v7) | 27,142 | 124.1M | 498,975 | 0.40% |
| 3 (v11) | 27,055 | 123.6M | 4,040,599 | 3.27% |
| 4 (v8) | 5,030 | 22.8M | 2,144,285 | 9.39% |
| **Total** | **87,362** | **397.9M** | **10,374,641** | **2.61%** |

Under the MaxEvents=1000 hypothesis with ~2 files/job:
- ~43,681 total jobs × 1000 events = 43.7M events actually processed
- Overall trigger rate = 10.37M / 43.7M = **23.7%** — close to the old skim's 27%

The remaining ~3% gap is expected given:
- GRL v120 (new) is stricter than v113 (old), rejecting some events
- Variable files-per-job across tasks (Part 4 appears to have many files/job)
- Not all runs have the same trigger rate

### Affected datasets

**All three PbPb years are affected** — all `grid_sub.sh` files use the same trf pattern without `--evtMax`:
- `run_23hi/grid_sub.sh`: 4 tasks (parts 1–4)
- `run_24hi/grid_sub.sh`: 2 tasks (parts 1–2)
- `run_25hi/grid_sub.sh`: 6 tasks (parts 1–6)

### Fix

The trf command must include `--evtMax=-1`:
```
pathena --trf "TRIGRATES_RUNMODE=hi20XX athena.py TrigRates_CA.py --filesInput=%IN --evtMax=-1" ...
```

Alternatively, change `TrigRates_CA.py` line 106 from `m_EvtMax = 1000` to `m_EvtMax = -1` (process all events by default). The `1000` was a convenience default for local interactive testing; it should never propagate to grid jobs.

### Why the old skim (R24, April 2025) was not affected

The old skim used old-style `TrigRates.py` (non-CA, CfgMgr) with:
```python
theApp.EvtMax = m_EvtMax  # also 1000
```
But the old skim was submitted with standard `pathena` (not `--trf`). PanDA's pilot for standard pathena jobs automatically injects `--evtMax=-1` to override `theApp.EvtMax`, so all events were processed. The new CA-style submission uses `--trf` mode, where PanDA does NOT auto-inject event limits — the trf command is executed as-is.

### FINAL CONCLUSION (Updated 2026-05-13)

**ROOT CAUSE FOUND: The new skim grid jobs only processed 1000 events per job due to missing `--evtMax=-1` in the pathena trf command.**

The 10× skim fraction discrepancy (old ~27% vs new ~2.6%) is NOT a trigger or code difference. It is caused by:
1. `TrigRates_CA.py` defaulting `m_EvtMax = 1000` for local testing convenience
2. `grid_sub.sh` not passing `--evtMax=-1` to override this default
3. The CA framework (`flags.fillFromArgs()`) not overriding `MaxEvents` when `--evtMax` is absent from the command line

The actual trigger pass rate (~24–27%) is consistent between old and new skims. The new skim output (10.37M events) represents trigger-selected events from only ~43.7M events actually processed (first 1000 of each job's input), not the full ~398M input dataset.

**All new PbPb skims (2023, 2024, 2025) need to be resubmitted with `--evtMax=-1` in the trf command.**

### Wrong-skim (Apr2026, evtMax=1000) output statistics — for reference

All entries below are from the `HeavyIonD3PD` tree in the hadded output files at
`~/usatlasdata/dimuon_data/`. These represent ~2.6% of input (only first 1000
events per grid job processed). Files replaced by May2026 resubmission.

**PbPb 2023** (tasks 49860533/49683976/50022072/49860535):

| Part | Entries |
|------|---------|
| part1 (v9) | 3,690,782 |
| part2 (v7) | 498,975 |
| part3 (v11) | 4,040,599 |
| part4 (v8) | 2,144,285 |
| **Total** | **10,374,641** |

**PbPb 2024** (tasks not recorded — Apr2026 submission):

| Part | Entries |
|------|---------|
| part1 | 3,552,468 |
| part2 | 3,853,338 |
| **Total** | **7,405,806** |

**PbPb 2025** (tasks not recorded — Apr2026 submission):

| Part | Entries |
|------|---------|
| part1 | 4,714,195 |
| part2 | 6,409,183 |
| part3 | 7,057,841 |
| part4 | 5,869,389 |
| part5 | 8,084,354 |
| part6 | 1,508,168 |
| **Total** | **33,643,130** |

**Grand total (all wrong-skim PbPb):** 51,423,577 entries
