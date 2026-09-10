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

## Implementation Plan

| # | Step | Status |
|---|------|--------|
| 1 | Tracking doc + INDEX registration | DONE |
| 2 | Enumerate `data26_hi` HardProbes AOD datasets (latest f-tag/run), build partition plan + `InDstxt_*` files | pending |
| 3 | Add `hi2026` run mode to `scripts/TrigRates_CA.py`, sync to `run_26hi/`, check `Module_EventShape.cxx` year switch | pending |
| 4 | Recompile `SkimCode` (AthAnalysis 25.2.89, `build_25/`) | pending |
| 5 | Local test job on one 2026 AOD + branch sanity check | pending |
| 6 | Write `run_26hi/grid_sub.sh`, submit all parts | pending |
| 7 | LGD pre-flight (quota, DID conflicts, pnfs mount, plugin install) | pending |
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

## Results & Observations

*(to be filled)*

## Remaining Work

Everything from step 2 onward.

## Latest Stage

**Step 2 + 3 + 7 launched in parallel.** Dataset enumeration delegated to a subagent
(owns `run_26hi/InDstxt_*`), LGD pre-flight delegated to a second subagent (owns no repo
files), while the orchestrator writes the `hi2026` run mode and compiles.
