# PbPb 2026 Data Skim + LOCALGROUPDISK Storage Migration

**Created:** 2026-09-10 · **Mode:** Implementation

## Objective

1. Skim the 2026 Pb+Pb `physics_HardProbes` data (`data26_hi`) on the grid with the
   existing `TrigRates` skim, producing NTUPs equivalent to the 2023/24/25 skims.
2. Monitor the grid tasks and download the output with `scripts/grid_monitor.sh`,
   keeping the data-status bookkeeping current; sanity-check the downloaded NTUPs.
3. Free local GPFS storage by migrating the **raw skim NTUPs** (inputs to the NTuple
   processing code) for PbPb 2023/2024/2025 and pp 2024 to
   `BNL-OSG2_LOCALGROUPDISK`, replacing them with a same-path symlink farm so the
   analysis code keeps working unchanged.
4. If local storage stays tight, migrate the full PbPb2026 skim output to
   LOCALGROUPDISK as well and clean up locally.

## Autonomy Contract (ACTIVE — re-read on every compaction)
- Mandate: run autonomously to DONE; do NOT pause to confirm progress. Finishing a
  plan, a passing small test, or one pipeline stage is NOT a stopping point.
- Done =
  1. `hi2026` run mode implemented in `SkimCode/scripts/TrigRates_CA.py` (canonical) and
     synced to `SkimCode/run_26hi/`, using GRL `xmls/physics_HI2026_50ns_noIBL.xml` and the
     same muon/dimuon trigger lists as the other Run-3 PbPb years; skim source recompiled;
     committed.
  2. A local test job over a real `data26_hi` HardProbes AOD runs to completion and its
     output NTUP is branch-sanity-checked (muon kinematics, muon SF branches, trigger and
     trigger-matching branches, ZDC, FCal/centrality all present and filled).
  3. `data26_hi.*.physics_HardProbes.*.AOD.*` datasets enumerated (latest f-tag per run,
     **no x-tags, no duplicates**), partitioned into `InDstxt_PbPb2026_5p36TeV_part*.txt`
     within PanDA per-job / per-task limits, `run_26hi/grid_sub.sh` written, and all grid
     tasks SUBMITTED.
  4. Grid tasks monitored to completion (resubmitting/repartitioning if they fail on job
     count or per-job data volume); output downloaded and hadd-ed via `grid_monitor.sh`
     into `~/usatlasdata/dimuon_data/pbpb_2026/`; entry counts + branch content
     sanity-checked; `data-merging-record.txt` updated.
  5. PbPb23/24/25 + pp24 raw skim NTUPs on `BNL-OSG2_LOCALGROUPDISK` (rule state OK),
     byte-size and TTree-entry verified against the local originals (no data loss),
     same-path symlink farm built, NTuple-processing smoke test passing against the farm,
     and local originals deleted (one small raw file per data-taking period may be kept).
  6. If storage still tight: PbPb2026 raw skim NTUPs likewise on LOCALGROUPDISK with a
     verified symlink farm and local cleanup done.
- Stop-and-ask = ANY physics-results-bending ambiguity (no fixed list; use judgment;
  when unsure whether an ambiguity is blocking, treat it as blocking → AskUserQuestion).

## Physics Procedure

### 1. Motivation

The 2026 Pb+Pb run adds a fourth Run-3 heavy-ion data-taking period to the low-mass
dimuon measurement. Its skim must be **procedurally identical** to the 2023/24/25
skims so the years can be combined (and compared) without a selection-level offset:
same `TrigRates` algorithm, same muon selection / calibration / scale-factor tools,
same trigger chain lists, same stored branches. The only per-year inputs that legitimately
change are (a) the Good Runs List and (b) the run/dataset list.

### 2. What the skim produces

For every event passing the GRL, the vertex cut and the OR of the configured trigger
chains, one entry in the `HeavyIonD3PD` tree carrying: muon kinematics (calibrated),
Medium+Tight muon selection flags and efficiency scale factors, L1 + HLT trigger
decisions and prescales, single-muon and dimuon trigger matching, track multiplicities,
ZDC energies/times/PreSampleAmp, and HI event-shape quantities (FCal E_T A/C, centrality).

### 3. Step-by-step method

a. **Run mode `hi2026`** — `RunYear = 2026`, `is_HION = True`, `is_Run3 = True`,
   `is_MC = False`, `use_trigger = True`, `StoreAllEvents = False` (the data skim keeps
   only events firing the trigger OR — this is correct for data and must NOT be changed;
   only MC keeps every event, because MC needs unbiased efficiency denominators).
b. **GRL** — `xmls/physics_HI2026_50ns_noIBL.xml` (35 runs, periods J and K,
   defect tag `DetStatus-v144-pro58-01`, ignoring PIXEL/IBL and BSPOT-invalid defects).
c. **Trigger chains** — the same Run-3 Pb+Pb muon and dimuon lists used for 2023/24/25.
   The single-muon analysis trigger in Pb+Pb is `HLT_mu4_L1MU3V`; the dimuon list is
   `HLT_2mu4_L12MU3V` + `HLT_mu4_mu4noL1_L1MU3V`.
d. **Centrality** — `Module_EventShape.cxx` maps Run-3 years 2023/24/25/**2026** onto the
   PbPb2023 FCal-E_T thresholds (per-year calibrations not finalised). The downstream
   FCal rescaling to 2023 conditions is an NTuple-processing-stage correction and is
   NOT part of this skim.
e. **Dataset selection** — `data26_hi.*.physics_HardProbes.*.AOD.*`, **f-tags only**
   (production reprocessing), never x-tags (express/debug-stream style tags); where one
   run has more than one f-tag, keep only the **latest** so no run is double-counted.
   Double-counting a run would double its luminosity in every yield.

### 4. Negative constraints

- Do NOT set `StoreAllEvents = True` for data — that would blow up the output size and
  change what a "skimmed event" means relative to 2023/24/25.
- Do NOT include the same run twice under two f-tags.
- Do NOT invent a new trigger list for 2026; a different chain set would give the
  2026 year a different effective trigger efficiency from the other Run-3 years.
- Do NOT change any muon tool setting (selection WP, calibration mode, SF calibration
  release) relative to the other Run-3 years.

## Context

- Skim code: `SkimCode/source/HFtrigValidation/`; config `SkimCode/scripts/TrigRates_CA.py`
  (canonical) copied into each `run_*/`; build in `SkimCode/build_25/` (AthAnalysis 25.2.89).
- `SkimCode/README.md` documents run modes, submission, `grid_monitor.sh`, and the 2025
  dataset inventory / partition plan used here as a sizing reference.
- Storage (2026-09-10, GPFS `atlasgpfs01`, **numbers halved for real usage**):
  `usatlast3-data` 1.38 TB used / 1.54 TB quota → **~157 GB real free (90 % full)**.
- Local raw skim NTUPs (real sizes): pbpb_2023 125 GB (4 files), pbpb_2024 100 GB (2),
  pbpb_2025 245 GB (6), pp_2024 437 GB (12) → **~907 GB** recoverable.
- LOCALGROUPDISK procedure and its gotchas are recorded in
  `localgroupdisk_migration.md` (CLOSED) — upload to `BNL-OSG2_SCRATCHDISK`, then
  `rucio add-rule ... 1 BNL-OSG2_LOCALGROUPDISK`, then build a symlink farm from
  `rucio list-file-replicas --pfns`. Quota there is 50 TB.

### External constraint — BNL dCache maintenance

**BNL US ATLAS dCache storage element is DOWN for a major software upgrade on
2026-09-14, 09:00–17:00 EDT (13:00–21:00 UTC).** Consequences:
- Uploading to LOCALGROUPDISK via Rucio/FTS should be unaffected in principle, but the
  **pnfs symlink farm** (`/pnfs/usatlas.bnl.gov/LOCALGROUPDISK/...`) will not be readable
  during that window, and freshly created rules may not transfer.
- Therefore: **try to complete all uploads and symlink-farm construction before
  2026-09-14**. If a step lands in the window and misbehaves, pause it and resume after
  the service is back rather than fighting transient errors.
- The current VOMS proxy expires ~2026-09-14 16:30 UTC — i.e. inside the window. A proxy
  renewal by the user may be needed to finish late steps.

### Coordination with the analysis-code agent

Analysis-code support for PbPb 2026 (NTuple processing, event selection, centrality,
FCal scaling) is assigned to a **different agent in a different session**. The paths and
names this work will produce, for that agent to consume:

| Item | Value |
|---|---|
| Grid outDS prefix | `user.yuhang.TrigRates.dimuon.PbPb2026data.Sep2026.v1.part<N>.` |
| Download dir | `~/usatlasdata/dimuon_data/pbpb_2026/` |
| Merged NTUP names | `data_pbpb26_part<N>.root` |
| Tree name | `HeavyIonD3PD` |
| Skim run mode | `hi2026` (env `TRIGRATES_RUNMODE=hi2026`) |
| GRL | `SkimCode/xmls/physics_HI2026_50ns_noIBL.xml` |
| grid_monitor mapping | `PbPb2026data...partN._EXT0` → `pbpb_2026/data_pbpb26_partN.root` |

## Scope

IN: skim config, compile, test job, dataset enumeration + partitioning, grid submission,
monitoring, download, NTUP sanity checks, LOCALGROUPDISK migration of raw skim NTUPs
(pbpb23/24/25, pp24, and pbpb26 if needed), symlink farms, storage cleanup.

OUT: any change to the NTuple-processing / RDF / plotting analysis code to *use* the 2026
data (assigned to another agent/session); physics results from the 2026 data.

## Design Decisions

*(appended as taken)*

### D1 — the six 2026 AOD runs absent from the GRL are EXCLUDED from the submission

**Runs 522011, 522089, 522112, 522145, 522185 (the pre-GRL block at the start of the
year) and 523168 (an *interior* run, between GRL runs 523138 and 523188) have
`merge.AOD` datasets but do not appear in `physics_HI2026_50ns_noIBL.xml`.**

Decision: **exclude them from the numbered partition files**, but write them to
`run_26hi/InDstxt_PbPb2026_5p36TeV_nonGRL.txt` so a reissued GRL can be served by
re-running just those.

Rationale — this is *not* a physics choice. The skim applies the GRL inside the job
(`alg.UseGRL = True`), so every event of a non-GRL run is dropped by `TrigRates` whether
or not the dataset is submitted. Submitting them only burns grid CPU to produce empty
output — precisely what happened to the 2024 `part1` task, which was resubmitted (v6) with
its pre-GRL runs removed. Excluding them therefore changes no yield, no luminosity and no
efficiency; it changes only the CPU bill. 523168 is flagged separately because an interior
omission is more likely to be a certification lag than a bad run.

### D2 — the 2026 skim needs AthAnalysis **25.2.90**, not the 25.2.89 used for 2023/24/25

**Forced, not chosen.** The first local test job on a real 2026 AOD (run 522474, lb0122)
under 25.2.89 aborted on the **first event**:

```
ERROR: problem when building L1 menu structure (algorithms).
       Flavour gRISTRETTO for EnergyThreshold algorithm not recongnised!
IncidentProcAlg1 FATAL Standard std::exception is caught in sysExecute
Py:Athena INFO leaving with code 65: "failure in an algorithm execute"
```
The 2026 L1 menu contains a gFEX `EnergyThreshold` algorithm flavour, `gRISTRETTO`, that
25.2.89's `TrigConfData` L1-menu parser does not know. This is not a configuration
mistake and cannot be worked around in `TrigRates_CA.py` — the failure happens inside
`TrigDecisionTool`'s menu decoding, which the skim needs for `StoreL1Decision`.

**Release chosen: 25.2.90** — established by `strings` on
`libTrigConfData.so` across the release series that **25.2.90 is the FIRST release
containing `gRISTRETTO`** (25.2.89 no; 25.2.90…25.2.110 yes). It also keeps the same
platform as 25.2.89 (`x86_64-el9-gcc14-opt`, LCG_108a_ATLAS_9), whereas 25.2.110 has
moved to gcc15. So 25.2.90 is the smallest possible step away from the release the other
Run-3 years were skimmed with. New setup script `SkimCode/setup_26.sh`, build dir
`build_26/`; `setup_25.sh` / `build_25/` are left untouched so 2023/24/25 remain
reproducible.

**The step is proven physics-neutral — verified, not assumed.** The muon libraries do
differ between 25.2.89 and 25.2.90 (`libMuonEfficiencyCorrectionsLib.so`,
`libMuonMomentumCorrectionsLib.so`, `libMuonSelectorToolsLib.so` all change md5 and size;
25.2.90 and 25.2.91 are identical to each other), so a binary-level "nothing changed"
argument was **not** available. Instead an explicit A/B was run: the **same** `data25_hi`
AOD (`00512049 ... f1671_m2272._lb0166._0007.1`, 2.553 GB, lb 166 inside GRL range
147–168) was skimmed under `TRIGRATES_RUNMODE=hi2025` with **both** releases.

> **Result: 488 events out of 1311 in both; output files byte-identical in size; and all
> 169 branches agree exactly** (entry count, sum and sum-of-squares of every element of
> every branch, compared by `/tmp/.../cmp_ntup.C`) — **0 branches differ.**

Consequence: the 2023/24/25 skims do **not** need to be regenerated, and the 2026 skim is
directly combinable with them. No muon calibration, selection or scale-factor value moves.

## Implementation Plan

| # | Step | Status |
|---|------|--------|
| 1 | Tracking doc + INDEX registration | DONE |
| 2 | Enumerate `data26_hi` HardProbes AOD datasets (latest f-tag/run), build partition plan + `InDstxt_*` files | **DONE** |
| 3 | Add `hi2026` run mode to `scripts/TrigRates_CA.py`, sync to `run_26hi/`, check `Module_EventShape.cxx` year switch | **DONE** |
| 4 | Recompile `SkimCode` (AthAnalysis 25.2.89, `build_25/`) | **DONE** |
| 5 | Local test job on one 2026 AOD + branch sanity check | **DONE** |
| 6 | Write `run_26hi/grid_sub.sh`, submit all parts | **DONE** (5 tasks) |
| 7 | LGD pre-flight (quota, DID conflicts, pnfs mount, plugin install) | **DONE** |
| 8 | Migrate pbpb23/24/25 + pp24 raw NTUPs → LGD, symlink farm, verify, smoke test, delete originals | pending |
| 9 | Monitor grid tasks → download → hadd → validate (`grid_monitor.sh`) | pending |
| 10 | Sanity-check downloaded pbpb26 NTUPs | pending |
| 11 | If tight: migrate pbpb26 NTUPs → LGD + cleanup | pending |
| 12 | Final: commit, close doc | pending |

## Progress Log

### 2026-09-10 — Step 1

Doc created. Environment surveyed:
- `SkimCode/run_26hi/` exists but is **empty**.
- GRL `xmls/physics_HI2026_50ns_noIBL.xml` present: 35 runs (522041 … 523437),
  periods `data26_hi.periodJ` + `periodK`, defect tag `DetStatus-v144-pro58-01`,
  ignoring `PIXEL_PERFORMANCE_INTOLERABLE, TRIG_HLT_IDT_BSPOT_INVALID_STATUS,
  ID_IBL_TRACKCOVERAGE_SEVERE, PIXEL_IBL_DISABLED`.
- `scripts/TrigRates_CA.py` (381 lines) has modes hi2023/hi2024/hi2025/pp2024/
  ppmcfullsim2024; `hi2026` must be added in 4 places (flag decl, `_set_run_mode`,
  the three dispatch branches, the `is_HION`/`is_Run3`/guard conjunctions, and the
  per-mode config block).
- VOMS proxy valid, 95 h 57 m left at 16:30 UTC 2026-09-10 (→ expires ~2026-09-14 16:30 UTC).
- Rucio account `yuhang` ACTIVE.
- GPFS quota (halved): 1.38 TB / 1.54 TB used → ~157 GB real free.

### 2026-09-10 — Steps 3 + 4 (config + compile) DONE

**`hi2026` run mode added to `SkimCode/scripts/TrigRates_CA.py`** (canonical copy) in
seven places: flag declaration, `_set_run_mode()` global + assignment, the
`TRIGRATES_RUNMODE` whitelist tuple, the `--filesInput` substring dispatch
(`data26_hi` → `hi2026`), the run-directory dispatch (`run_26hi` → `hi2026`), the
`is_HION` and `is_Run3` conjunctions, the CA-only guard, and a new `elif do_hi2026:`
configuration block.

Configuration chosen (per §3 of the Physics Procedure):
- `RunYear = 2026`
- `GRL = ["physics_HI2026_50ns_noIBL.xml"]`
- `Muon_triggers` and `DiMuon_triggers` **byte-identical to the `do_hi2025` block**:
  `HLT_mu4_L1MU3V, HLT_mu6_L1MU3V, HLT_mu6_L1MU5VF, HLT_mu8_L1MU5VF, HLT_mu10_L1MU8F,
  HLT_mu10_L1MU5VF` and `HLT_2mu4_L12MU3V, HLT_mu4_mu4noL1_L1MU3V`.
  *(2024 additionally carried five HI-specific `mu4noL1_*` UCC/ZDC chains that 2023 and
  2025 do not; the 2023 ∩ 2025 list is the stable "Run-3 Pb+Pb" list and is what 2026
  inherits. The local test job must confirm these chains resolve in the 2026 HLT menu.)*
- `InputFile` is still the literal `PLACEHOLDER_DATA26_AOD` — to be replaced with a real
  disk-resident 2026 AOD before the local test. It is a local-test default only; the grid
  always overrides it via `--filesInput=%IN`.
- Everything else (muon selection/calibration/SF tools, `StoreAllEvents=False`,
  `StoreZdc=1`, `StoreTracks=1`, `StoreEventInfo=1`) is inherited unchanged from the
  shared code path — no per-year divergence.

**No C++ change was needed.** `Module_EventShape.cxx:121-131` already has `case 2026:`
falling through to `GetCentralityPbPb2023()`, and `m_year` is used *only* for centrality
(`grep -rn m_year source/` → 6 hits, all in the EventShape path). Recompiled anyway with
`source setup_25.sh` (AthAnalysis 25.2.89, `build_25/`): all targets **already up to
date**, build clean, exit 0.

**Config synced** `scripts/TrigRates_CA.py` → `run_23hi/ run_24hi/ run_24pp/ run_25hi/
run_26hi/`. Note `run_24pp/TrigRates_CA.py` had been **stale** (76 diff lines — it
pre-dated the July-2026 `mc_has_trigger_sim` / `TRIGRATES_OUTPUT_FILE` work); the sync
brings it current. `run_23hi/ run_24hi/ run_25hi/` were byte-identical to the committed
`scripts/` copy before the sync.

**`scripts/grid_monitor.sh` extended for 2026** in three places:
`*PbPb2026*` → `pbpb_2026/` + `data_pbpb26_part<N>.root` (outDS→dir map, line 199), the
chunked-hadd `file_prefix` case (line 409), and `get_code_update_info()` (line 351) which
deliberately returns **empty** for `pbpb_2026` — `PbPbExtras.c` has no `{26, ...}` entry
and there is no `run_pbpb_26.sub`, so the NTuple-processing auto-update is skipped rather
than failing. The agent adding 2026 analysis support must set `file_batch_max` from the
part count actually on disk. `bash -n` clean.

### 2026-09-10 — Step 5 tooling: skim-NTUP sanity checker

Wrote `SkimCode/scripts/check_skim_output.C`. It (a) diffs the branch list of a new NTUP
against a validated reference NTUP and (b) reports, per branch, the **fraction of entries
in which the branch is actually filled** (scalar != 0 / array any != 0 / vector size > 0),
grouped as event-info, muon kinematics, muon SF, trigger, ZDC, event-shape/FCal, tracks,
plus a full list of always-empty branches. This targets the specific silent failure mode
of a new year: the branch exists (the skim still writes it) but is never *set*, because a
container key or a tool did not resolve.

**Validated on known-good NTUPs** (`data_pbpb25_part6.root` vs `data_pbpb23_part4.root`):
- both have **168 branches**; the diff is *exactly* the 2023-vs-2025 trigger-list
  difference (2023-only `*_VTE50` chains missing; 2025-only `mu8_L1MU5VF`, `mu10_L1MU8F`,
  `mu10_L1MU5VF` chains extra) — nothing structural.
- always-empty in 2025: `L1TE`, `L1TE24`, `b_HLT_mu4_mu4noL1_L1MU3V` (3/168). The
  `mu4_mu4noL1` chain is known not to be maintained in Run 3 Pb+Pb.
- filled fractions on the good branches: ZDC / FCal / `trk_numqual` = 1.0000,
  `centrality` = 0.9515, single-muon branches (incl. `muon_eff_SF_medium`,
  `muon_eff_SF_tight`, `muon_eff_corr_{medium,tight}_{data,mc}`,
  `muon_b_HLT_*` trigger matching) = 0.9970, pair-level branches = 0.5635.

**Acceptance criterion for the 2026 test job (derived from the above):** because `hi2026`
uses the *same* chain lists as `hi2025`, the 2026 NTUP must have a branch list **identical
to `data_pbpb25_*`** (0 missing, 0 extra, 168 branches), the same always-empty set, and
non-zero fill fractions in every group above.

### 2026-09-10 — Step 2 (dataset survey, in progress) + Step 8 tooling

**Dataset survey (subagent).** `rucio list-dids "data26_hi:data26_hi.*physics_HardProbes*AOD*"`
returns **124 DATASET DIDs, of which 42 are `merge.AOD`**, and 43 CONTAINER DIDs which are
*all* `deriv.DAOD_HION5.*_p7386` (plus one PhysCont and one OpenEnded AOD) — i.e. **no
`merge.AOD` container exists**, so the DATASET query is the complete enumeration.

- **1 x-tag rejected:** `data26_hi.00522056.physics_HardProbes.merge.AOD.x973_m2281`
  (run 522056 also has f-tag `f1714_m2281`, which is kept).
- **41 f-tag datasets kept; NO run carries two f-tags**, so the "keep the latest f-tag"
  rule had nothing to resolve and no run can be double-counted. f-tags present:
  f1714 (13 runs), f1717 (2), f1720 (13), f1723 (12), f1727 (1).
- **GRL cross-check:** all **35** GRL runs have an f-tag `merge.AOD` dataset (none missing).
  Six AOD runs are *not* in the GRL → see Design Decision **D1**.

**Symlink-farm smoke test recipe (for step 8).** The NTuple processing reads its raw input
as `data_dir + "data_pbpb" + run_year + "_part" + file_batch + ".root"`
(`PbPbExtras.c:19`) from `/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_20YY/`
(`PbPbEventSelConfig.h:50`), and writes its *outputs* into the same directory — so a
per-file symlink farm is compatible: symlinked inputs and real output files coexist.
`DimuonAlgCoreT.h:122-123` exposes `nevents_max` and `is_test_run`, so the smoke test is

```
PbPbAnalysis pbpb(25, 1);
pbpb.is_test_run = true;    // appends "_test" -> cannot clobber production outputs
pbpb.nevents_max = 20000;
pbpb.Run();
```

run through `run_pbpb_25.sh`'s LCG view (`views LCG_107a_ATLAS_2 x86_64-el9-gcc13-opt`).

**Migration driver written:** `~/usatlasdata/lgd_migration/rawskim_2026/lgd_migrate_rawskim.sh`
— staged (`upload | dataset | rule | wait | farm | verify | all`) per group
(`pbpb23 pbpb24 pbpb25 pp24 pbpb26`), one Rucio dataset per group
(`user.yuhang:dimuon_rawskim_<key>`). It **never deletes anything** and **never calls
`rucio update-rule`** (the `--lifetime 0` trap that purged replicas in the 2026-06
migration). `farm` refuses to run unless the LGD PFN count equals the local file count,
parks the originals in `<dir>_orig_<key>/` rather than deleting them, and `verify`
compares byte size **and** `HeavyIonD3PD` entry count farm-vs-original per file.

### 2026-09-10 — Step 2 COMPLETE: 2026 dataset inventory and partition plan

**Enumeration.** `rucio list-dids "data26_hi:data26_hi.*physics_HardProbes*AOD*"
--filter type=DATASET` → 124 DIDs, **42 `merge.AOD`**. The CONTAINER query (43 DIDs) is
entirely `deriv.DAOD_HION5.*_p7386` plus one PhysCont and one `OpenEnded.AOD` — **no
`merge.AOD` container exists**, so the DATASET query is the complete enumeration.

- **41 f-tag datasets kept**: f1714 (13 runs), f1717 (2), f1720 (13), f1723 (12), f1727 (1).
- **1 x-tag rejected**: `data26_hi.00522056...merge.AOD.x973_m2281` (run 522056 keeps its
  f1714 dataset).
- **No run carries two f-tags** → the "keep the latest f-tag" rule had nothing to resolve,
  and no run can be double-counted.
- **6 non-GRL runs excluded** (Design Decision D1): 522011, 522089, 522112, 522145, 522185
  (pre-GRL block) and **523168** (interior *and* near-empty: 3 files / 17 905 events — the
  2026 analogue of 2025's excluded run 511013). Cost of excluding them: **977 files,
  8.7 TB, 4 562 969 events = 0.87 % of files, 0.63 % of events.** Parked in
  `run_26hi/InDstxt_PbPb2026_5p36TeV_nonGRL.txt` with a dated rationale header.

**Partition plan (35 GRL runs).** Contiguous-by-run grouping chosen to minimise the largest
part. K=4 → 486 jobs in the biggest part; **K=5 → 398 (chosen** — same per-part scale as
2025 with ~100 jobs of headroom below the ~500-jobs/task guideline**)**; K=6 → 345 (more
tasks than needed).

| Part | Runs | #ds | Files | Jobs @60 f/job | TB | Events |
|---|---|---:|---:|---:|---:|---:|
| 1 | 522041–522327 | 8 | 23 879 | 398 | 273.8 | 142 143 968 |
| 2 | 522336–522518 | 7 | 23 559 | 393 | 298.1 | 157 557 726 |
| 3 | 522541–522780 | 8 | 22 542 | 376 | 288.4 | 149 502 017 |
| 4 | 522830–522998 | 6 | 20 965 | 350 | 266.3 | 138 161 913 |
| 5 | 523023–523437 | 6 | 20 334 | 339 | 257.6 | 134 128 557 |
| **Total** | **522041–523437** | **35** | **111 279** | **1 856** | **1 384.2** | **721 494 181** |

**Average file: 12.41 GB / 6 468 events** — *larger* than 2025's 11.3 GB, so
`--nGBPerJob MAX` is even more necessary than last year (local staging would need
~745 GB/job).

**Coverage assertion:** `sorted(runs in part1..5) == sorted(35 GRL runs)` → **True**,
35/35, none extra; `concat(parts) == combined file`; zero overlap between the part files
and the nonGRL file. Independently re-verified by the orchestrator.

**Disk vs tape: 100 % DISK, ZERO TAPE.** All 41 datasets have a complete disk replica; all
41 are complete at `CERN-PROD_DERIVED`; 13 are additionally complete at
`BNL-OSG2_DATADISK` (522089, 522125, 522327, 522474, 522492, 522518, 522628, 522658,
522837, 523138, 523372, 523418, 523437). `grep -i tape` over the full replica sweep: no
matches. **No staging rule is needed** — a real difference from the older `data23_hi`
reprocessings.

**Scale comparison.** 2026 (35 GRL runs, 111 279 files, 1 384 TB, 721 M events) is ~87 % of
2025 (47 runs, 143 528 files, 1 587 TB, 828 M events) in events, in 5 tasks instead of 6.

**Test file for the local job.** Run 522089 (the original candidate) is one of the excluded
non-GRL runs — with `alg.UseGRL = True` the test would have produced an **empty** tree and
the branch-fill check would have been meaningless. Replaced with:
```
data26_hi.00522474.physics_HardProbes.merge.AOD.f1720_m2281._lb0122._0014.1
1.487 GB / 731 events, complete at BNL-OSG2_DATADISK
root://dcgftp.usatlas.bnl.gov:1094//pnfs/usatlas.bnl.gov/BNLT0D1/rucio/data26_hi/b6/6b/<that file>
```
Orchestrator verified independently that **lb 122 falls inside run 522474's only GRL
LumiBlock range, 117–224**, so the test job must produce a non-empty tree.

### 2026-09-10 — Step 5 COMPLETE: local test job on a real 2026 AOD

First attempt under 25.2.89 **failed on event 1** with the gFEX `gRISTRETTO` L1-menu
error — see Design Decision **D2** for the diagnosis, the release survey, and the A/B
proof that moving to 25.2.90 changes no physics.

Rebuilt against 25.2.90 (`setup_26.sh` → `build_26/`) with **no source change**, then:

```
TRIGRATES_RUNMODE=hi2026 athena.py TrigRates_CA.py --evtMax=731
  -> Py:Athena INFO leaving with code 0: "successful run"
```

Output `myfile.root`: **317 of 731 events kept** (43 %, i.e. the trigger-OR skim is
working — `StoreAllEvents=False`), tree `HeavyIonD3PD`, **169 branches**.

**Branch-list diff vs the validated 2025 production NTUP: 0 missing, 1 extra.** The extra
branch is `muon_match_L1MU3V`, added to the skim source *after* the 2025 data was
produced (it is the branch the L1/HLT trigger-efficiency split needs). Its presence in a
fresh skim is expected and correct; nothing is missing.

**Fill fractions** — every group healthy and matching the 2025 rates:

| group | 2026 test | 2025 reference |
|---|---|---|
| ZDC (`zdc_*`, 8 branches) | 1.0000 | 1.0000 |
| `FCal_Et`, `FCal_Et_P/N`, `trk_numqual` | 1.0000 | 1.0000 |
| `centrality` | filled | 0.9515 |
| muon kinematics `muon_pt/eta/phi` + `*_precorr` | 0.9968 | 0.9970 |
| **muon SF** `muon_eff_SF_{medium,tight}`, `muon_eff_corr_{medium,tight}_{data,mc}` | 0.9968 | 0.9970 |
| single-muon trigger matching `muon_b_HLT_*` | 0.9968 | 0.9970 |
| pair-level branches | 0.6183 | 0.5635 |
| `b_HLT_2mu4_L12MU3V` | 0.0410 | — |

**Always-empty branches: exactly 3 — `L1TE`, `L1TE24`, `b_HLT_mu4_mu4noL1_L1MU3V` — the
identical set found in the 2025 reference.** `mu4_mu4noL1` is known not to be maintained
in Run-3 Pb+Pb. So every configured 2026 chain (`mu4`, `mu6` ×2, `mu8`, `mu10` ×2,
`2mu4`) resolves in the 2026 HLT menu; the inherited chain list is valid for 2026.

### 2026-09-10 — Step 6: grid submission

`run_26hi/make_grid_sub.sh` generates `grid_sub.sh` from whatever `InDstxt_*_part*.txt`
files are on disk (regenerating rather than hand-editing keeps the part count, outDS tags
and per-part run ranges in lock-step with the partition). 5 parts, campaign tag
`Sep2026.v1`, outDS `user.yuhang.TrigRates.dimuon.PbPb2026data.Sep2026.v1.part<N>.`,
`--nGBPerJob MAX --mergeOutput --extFile physics_HI2026_50ns_noIBL.xml`.
**Submitted from a shell set up with `setup_26.sh` so the tasks carry 25.2.90.**

### 2026-09-10 — Step 6 COMPLETE: all 5 grid tasks SUBMITTED

Submitted from `run_26hi/` in a shell set up with `../setup_26.sh` (AthAnalysis 25.2.90),
`lsetup panda`. `pathena` confirmed `CMTCONFIG=x86_64-el9-gcc14-opt` and build directory
`SkimCode/build_26` for every task. Exit 0, no errors.

| part | jediTaskID | runs | #ds | files | jobs @60 | TB | events |
|---|---|---|---:|---:|---:|---:|---:|
| 1 | **52488079** | 522041–522327 | 8 | 23 879 | 398 | 273.8 | 142.1 M |
| 2 | **52488080** | 522336–522518 | 7 | 23 559 | 393 | 298.1 | 157.6 M |
| 3 | **52488081** | 522541–522780 | 8 | 22 542 | 376 | 288.4 | 149.5 M |
| 4 | **52488083** | 522830–522998 | 6 | 20 965 | 350 | 266.3 | 138.2 M |
| 5 | **52488085** | 523023–523437 | 6 | 20 334 | 339 | 257.6 | 134.1 M |

outDS: `user.yuhang.TrigRates.dimuon.PbPb2026data.Sep2026.v1.part<N>.`
Task-ID file for `grid_monitor.sh`:
`~/usatlasdata/dimuon_data/sep2026_pbpb26_skim.txt`.
BigPanDA: <https://bigpanda.cern.ch/?user=yuhang>.

**`grid_monitor.sh` is deliberately NOT started yet.** It downloads as soon as a task
reports done, and there are currently only ~157 GB real free on the GPFS data fileset —
not enough for a 2026 part. It will be started once the LOCALGROUPDISK migration has freed
space. Freshly submitted tasks have nothing to download for many hours, so nothing is lost
by waiting.

### 2026-09-10 — Step 8 in progress: uploads to SCRATCHDISK

`lgd_migrate_rawskim.sh pp24 upload` running; 9 of 12 pp24 files on
`BNL-OSG2_SCRATCHDISK` after ~30 min (~17 GB/min sustained, davs). Order is
pp24 → pbpb25 → pbpb23 → pbpb24 (largest first).

Pre-flight (subagent) established, before any upload:
- **LOCALGROUPDISK: 50 TB limit, 1.720 TB used, 48.28 TB free.** SCRATCHDISK was at 0 B
  (the June staging aged out).
- **23 existing rules, all `OK`** — zero REPLICATING, zero STUCK.
- **The raw skims are NOT already uploaded.** The only `dimuon_data`-ish DID under
  `user.yuhang` is `dimuon_data_test_aod_and_backups` (analysis-output backups + 4 test
  AODs). This is a from-scratch migration.
- **DID conflicts: 0 of 24.** Every one of the 24 target filenames is free, no
  cross-year collision, and none begins with `user.yuhang.` — so the server-500
  name-stripping workaround from the June migration is not needed.
- **The symlink-farm approach still works today**: `/pnfs/usatlas.bnl.gov/LOCALGROUPDISK/
  rucio/user/yuhang/` lists 256 hash dirs, and a spot-checked symlink under
  `~/dcachearea/pythia_truth_full_sample/pythia_5p36TeV/` resolved to a real
  8 881 262 188-byte file whose `HeavyIonD3PD` tree ROOT opened (480 000 entries)
  **with `X509_USER_PROXY` deliberately broken** — no grid proxy needed for reads.
- **No-data-loss baseline recorded** in `_sub_lgd_preflight_2_baseline.txt`: 24 rows of
  `relpath bytes GiB tree entries mtime`. All 24 open cleanly, every tree is
  `HeavyIonD3PD`, **1 080 739 586 entries / 974 421 871 775 bytes = 907.50 GiB** total.

**⚠ The plugin's `migrate` skill must NOT be pointed at the period directories.** It
uploads `$source_dir/*.root` wholesale and, in same-path-swap mode, `mv`s the whole
directory aside. Those directories are not clean: pbpb_2023 holds 90 `.root` of which only
4 are targets, pbpb_2024 79/2, pbpb_2025 81/6, pp_2024 191/12 — the other ~23 GiB are live
pipeline outputs (`muon_pairs_*`, `histograms_real_pairs_*`, `hists_cut_acceptance_*`),
plus 11 non-`.root` files and 15 subdirectories. A blind same-path swap would register
analysis outputs as immutable replicas, turn them into read-only symlinks the pipelines
could no longer rewrite, and orphan the rest. Hence the purpose-built
`lgd_migrate_rawskim.sh`, which touches **only** the `data_*_part*.root` files by name.

**Plugin status:** `bnl-localgroupdisk@usatlas-marketplace` v1.1.0 installed (user scope,
enabled) — but from the local fork `~/workarea/marketplace-fork` (branch
`add-bnl-localgroupdisk-plugin`, HEAD `f1d6be3`), because (a) `lsetup` puts LCG git 2.29.2
first on `PATH` and it has no `remote-https` helper, so any `marketplace add <repo>` fails
with `git: 'remote-https' is not a git command` (fix: `export PATH=/usr/bin:/bin:$PATH`),
and (b) `FlamyFlame/claude-bnl-localgroupdisk` is a *plugin* repo, not a marketplace, while
upstream `usatlas/marketplace` still has no entry — **PR #61 is not merged**. The skill
will only be visible in a *new* session. Its full procedure is transcribed in
`_sub_lgd_preflight_2.md`, so nothing depends on the skill being loadable here.

### 2026-09-10 — Step 8: uploads done for pp24; two driver bugs found and fixed

**pp24 complete through the rule stage.** 12/12 files on `BNL-OSG2_SCRATCHDISK`
(~17 GB/min sustained over davs, 437 GB in ~40 min), dataset
`user.yuhang:dimuon_rawskim_pp24` holds 12/12, LGD rule
**`c2df30b0707d497588778c26712b4b78`** created. pbpb25 → pbpb23 → pbpb24 following.

Two bugs in `lgd_migrate_rawskim.sh`, both found by running it, both fixed:

1. **`set -u` silently killed the script inside `source ~/setup.sh`.** The ATLAS setup
   scripts read unbound variables, so with `set -u` the shell died mid-source — before
   any log line was written, so the failure looked like "nothing happened at all"
   (`upload_all.log` showed four instant `UPLOAD FAILED` lines with no per-group log).
   Fix: `set +u` around the sourcing block, `set -u` after.
2. **`rucio list-dids <scope>:<exact-file> --short` returns NOTHING for a file DID** —
   `--filter type=FILE` is required. Without it the "already uploaded?" check always
   answered "free", so every completed upload was retried and failed on the duplicate,
   which aborted the group before its dataset/rule stages. Fix: add the type filter and
   anchor the match.

A third, subtler failure: the progress logger was `echo ... | tee -a "$LOG"`. When the
harness moved the foreground run into the background, the parent's stdout pipe closed and
every subsequent `tee` died on SIGPIPE — so the uploads kept working but stopped being
logged, and the log went silent after file 3 of 12 while all 12 actually succeeded. Fixed
by appending to the log with a plain redirect and echoing to stdout separately.

**Stages added to the driver:** `keep` restores one small *real* file per data-taking
period over its symlink — `data_pbpb23_part4` (8.8 GB), `data_pbpb24_part1` (48.3 GB),
`data_pbpb25_part6` (9.7 GB), `data_pp24_part11` (1.2 GB), 68 GB total — so local tests
(and any work during the 09-14 dCache outage) do not depend on `/pnfs` being up. Those
files remain on LOCALGROUPDISK too; the redundancy is deliberate. `purge` is the only
stage that deletes, and it re-runs `verify` first.

**Deletion plan (staged hedge).** After farm+verify passes: purge `pp_2024` and
`pbpb_2025` originals (~683 GB → ~829 GB real free, ample for the 2026 download), but
**keep `pbpb_2023` + `pbpb_2024` originals parked until after the 09-14 dCache upgrade**
and a post-outage re-verify. Holding them costs nothing (that space is not needed yet) and
hedges the only real risk in this migration — a bad outcome from the dCache upgrade, since
LOCALGROUPDISK is disk-only with no tape backup.

### 2026-09-10 — Step 9: grid tasks are running

First status sweep ~25 min after submission: 52488079 **running** (41 input files done),
52488080 **running** (35 done), 52488081/83/85 still **scouting**. Zero failed files.

52488080 carries a brokerage note — *"no candidates. brokerage failed for 3 input datasets
when trying 5 datasets"* — while nonetheless progressing. This is JEDI's transient
site-brokerage message and it retries; it is only a repartitioning signal if a task stalls
on it with no file progress. Being watched.

## Results & Observations

*(to be filled)*

## Remaining Work

Everything from step 2 onward.

## Latest Stage

**Step 2 + 3 + 7 launched in parallel.** Dataset enumeration delegated to a subagent
(owns `run_26hi/InDstxt_*`), LGD pre-flight delegated to a second subagent (owns no repo
files), while the orchestrator writes the `hi2026` run mode and compiles.
