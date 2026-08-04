# Pythia fullsim pp24 FULL sample — grid skim + storage headroom

**Mode:** Implementation. **Created:** 2026-07-13. **Session:** "pp24 fullsim full sample".
**Reviewer rules:** SkimCode/Athena config changes → `/review-analysis-code`.
**Siblings:** `mc_trigger_info_skim.md` (CLOSED — the pp24-fullsim *test*-sample skim, trigger
enablement, `grid_sub.sh`, `grid_monitor.sh --mode fullsim_pp`; **its open question "role of the
`_pdf` production 803015–803020" is answered here: that IS the full sample**),
`grid_monitor_data_infra.md` (grid_monitor internals), `localgroupdisk_migration.md`
(the 2026-06 pnfs→`BNL-OSG2_LOCALGROUPDISK` migration recipe), `mc_trigger_efficiency.md`
(the downstream consumer of these NTUPs).

---

## Objective

**Task I.** Take the newly-produced **full** Pythia fullsim pp24-conditions sample (6 DSIDs,
`_pdf` production 803015–803020, r16578) from AOD to skimmed NTUP on disk:
(1) a local **test run** proving the job runs and the AOD carries everything the analysis needs
(truth + reco kinematics, truth-origin branches, scale factors, trigger simulation);
(2) grid submission of all 6 datasets; (3) `grid_monitor.sh` monitoring → download → hadd →
validate → bookkeeping.
Directory reorganisation: `SkimCode/run_pythia_fullsim` → `run_pythia_fullsim_test_sample`;
new `run_pythia_fullsim_full_sample` holds everything for this task. Skimmed NTUPs +
bookkeeping txt go to `~/usatlasdata/pythia_fullsim_full_sample/`.

**Task II.** GPFS data-area quota is nearly exhausted. Census the data area, identify the
datasets worth relocating to `BNL-OSG2_LOCALGROUPDISK`, present the list, **wait for the user's
decision**, then delegate the migration to a subagent using the migration skill.

**Task II is a prerequisite for Task I's download stage** — there is not enough free space for
the full-sample NTUPs (see Physics/Operational Procedure §4).

## Autonomy Contract (ACTIVE — re-read on every compaction)
- Mandate: run autonomously to DONE; do NOT pause to confirm progress. Finishing a
  plan, a passing small test, or one pipeline stage is NOT a stopping point.
- Done =
  1. `SkimCode/run_pythia_fullsim_test_sample/` (renamed) and `run_pythia_fullsim_full_sample/`
     (new, holding the full-sample `grid_sub.sh` + `TrigRates_CA.py` + test-run artefacts) exist.
  2. A local test skim over ≥1 AOD file of the full sample ran to `rc=0`, and its NTUP was
     **explicitly verified** to contain: reco muon kinematics, truth muon kinematics, truth-origin/
     provenance branches (`muon_truth_*`, `truth_*`, parent/barcode navigation), MC event weights
     + scale factors, and the trigger branches (`b_HLT_*`, per-muon and per-leg dimuon matching).
  3. All 6 full-sample DSIDs submitted to the grid (task IDs recorded in this doc).
  4. `grid_monitor.sh` running against them in a mode that writes to
     `~/usatlasdata/pythia_fullsim_full_sample/` (NTUPs + `merging-record.txt`), monitored until
     all tasks resolve; outputs sanity-checked (non-empty trees, entry counts match the grid).
  5. Task II: the storage census + relocation candidate list delivered to the user; **user decision
     obtained**; if the user approves, the migration is executed by a subagent and the freed space
     confirmed via `mmlsquota`.
- Stop-and-ask = ANY physics-results-bending ambiguity (no fixed list; use judgment;
  when unsure whether an ambiguity is blocking, treat it as blocking → AskUserQuestion).

## Physics / Operational Procedure (AUTHORITATIVE)

### 1. What this sample is

6 datasets, one per pT-hat slice, **pp beam only**:

| DSID | slice | files | events | AOD size |
|---|---|---|---|---|
| 803015 | pTH125_300 | 32 | 320 000 | 275.9 GB |
| 803016 | pTH14_24 | 360 | 3 600 000 | 2.571 TB |
| 803017 | pTH24_40 | 240 | 2 400 000 | 1.793 TB |
| 803018 | pTH40_70 | 120 | 1 200 000 | 941.3 GB |
| 803019 | pTH70_125 | 32 | 320 000 | 263.0 GB |
| 803020 | pTH8_14 | 203 | 2 030 000 | 1.413 TB |
| **total** | | **987** | **9 870 000** | **7.26 TB** |

Full DID pattern:
`mc23_5p36TeV:mc23_5p36TeV.<DSID>.Py8EG_A14_pp_hQCD_DiMu_<slice>_pdf.recon.AOD.e8599_e8586_s4521_s4483_r16578`

Reco tag **r16578** — the SAME tag as the test sample ⇒ per `mc_trigger_info_skim.md` §R1 it
carries `doRDO_TRIG` + `doTRIGtoALL`, HLT menu `PhysicsP1_pp_lowMu_run3_v1`, and `mu4` / `2mu4` /
`mu4_mu4noL1` at prescale 1. The test run must nevertheless **confirm this in the produced file**
(an AMI flag says what was requested; the file says what is there).

**⚠ Difference from the test sample (flag, do not silently absorb):** the test sample was
802758–802781 = 4 isospin beams {pp,pn,np,nn} × 6 slices, combined downstream with the Pb-like
isospin ratio 4:6:6:9 (`pythia_fullsim_test_sample/isospin_weight.md`). The `_pdf` full production
has **only the `pp` beam** (verified: `rucio list-dids mc23_5p36TeV.803*.Py8EG_A14_*hQCD_DiMu*`
returns `pp` names only). The `_pdf` tag and the extra e-tag `e8586` suggest the isospin/nuclear
content is now carried by the PDF choice instead of by 4 explicit beams. **Consequence:** the
downstream isospin weighting must NOT be applied to this sample, and its AMI cross-sections
(`crossSection` × `genFiltEff`, in **nb**) must be re-read for the new DSIDs. This is a
downstream (NTuple-processing) concern, out of scope for the skim itself, but it is recorded here
because it changes how the sample is combined.

### 2. Skim configuration (must be IDENTICAL to the test sample)

Run mode `ppmcfullsim2024` in `TrigRates_CA.py`:
- `mc_has_trigger_sim = True` ⇒ `UseTrigger`, `StoreL1Decision`, `TrigDecisionTool`,
  `R3MatchingTool`, `HLTMuonsKey = HLT_MuonsCB_RoI`, `HLTMuonsFSKey = HLT_MuonsCB_FS`.
- **`StoreAllEvents = True` for MC** — CRITICAL. With `UseTrigger=True` and
  `StoreAllEvents=False` the skim silently becomes a **trigger-OR-filtered** sample, destroying
  the reco-efficiency and trigger-efficiency denominators (`mc_trigger_info_skim.md` §R1b).
  The test run MUST verify NTUP entries == AOD events read.
- No change to any muon selection, calibration, WP, or kinematic setting.

### 3. What the test run must prove (Task I item 2)

Open the test NTUP and confirm the presence AND non-emptiness of:
- (a) **reco muon kinematics** — `muon_pt` (SIGNED: charge is the sign — see
  `mc_trigger_info_skim.md` 8e), `muon_eta`, `muon_phi`, `muon_quality`, `muon_trk_pt`,
  `muon_d0`, `muon_z0`;
- (b) **truth muon kinematics** — `truth_muon_*` (truth record) and `muon_truth_*`
  (reco→truth match);
- (c) **truth-origin / provenance** — the branches the provenance classifier needs
  (`reference_muon_truth_provenance_nav`): truth barcodes, PDG ids, parent links, match
  probability, `IsPrimary`-type flags;
- (d) **MC weights / scale factors** — event weight(s) and muon SFs, whatever the test-sample
  NTUP carries (compare branch lists: the full-sample NTUP must be a superset-or-equal);
- (e) **trigger simulation** — `b_HLT_mu4_L1MU3V`, `b_HLT_2mu4_L12MU3V`,
  `b_HLT_mu4_mu4noL1_L1MU3V`, per-muon `muon_b_HLT_*`, and the 8 per-leg dimuon matching branches
  `dimuon_b_HLT_2mu4_L12MU3V_mu{1,2}pass Leg{1,2}_dR_*` (the ΔR-correlation ingredient).
- (f) **entry count == events read** (the `StoreAllEvents` check).

The decisive comparison is against the **test-sample** NTUP branch list: any branch present there
and absent here is a red flag → STOP and report.

### 4. Storage (Task II) — why it blocks Task I

`mmlsquota atlasgpfs01`, fileset `usatlast3-data`, **numbers halved** (GPFS reports allocated
blocks = 2× real bytes; `du -sb` = real bytes — memory `reference_storage_quota`):

| | GiB |
|---|---|
| used | **1471** |
| soft quota | 1536 |
| hard limit | 1792 |
| **headroom to soft** | **65** |

Expected full-sample NTUP size: the test-sample pp24 NTUPs are ~180 MB / 10 000 events
(≈18 kB/event) ⇒ 9.87 M events ⇒ **≈ 180 GB**, plus a transient download+hadd peak of up to
~65 GB on the largest slice (pTH14_24). **Required headroom ≈ 250 GB; available 65 GiB.**
⇒ space must be freed BEFORE the download stage. Submission and monitoring can start immediately;
the download is what must not begin under-provisioned. (Precedent: the 2026-07-09 storage
incident in `mc_trigger_info_skim.md` — a quota-exhausted `rucio download` failed on every RSE in
<1 s, and `grid_monitor` had already renamed the canonical NTUP to `.bak`, briefly leaving the
`.bak` as the only copy.)

### 5. Negative constraints

- Do **NOT** delete or move any `.bak_*` NTUP until its canonical replacement is verified on disk
  (memory `project_grid_monitor_bak_rename`).
- Do **NOT** write the full-sample NTUPs into `pythia_fullsim_test_sample/` — they go to
  `~/usatlasdata/pythia_fullsim_full_sample/`. The test-sample NTUPs stay where they are (the
  `mc_trigger_efficiency.md` results depend on them).
- Do **NOT** change any selection/WP/calibration while skimming the full sample.
- Do **NOT** apply the 4:6:6:9 isospin weight to this sample (see §1).

## Scope

**In:** run-dir reorganisation; full-sample `grid_sub.sh`; local test job + NTUP content
verification; grid submission of the 6 DSIDs; a `grid_monitor.sh` mode writing to
`pythia_fullsim_full_sample/`; download/hadd/validate/bookkeeping; the storage census +
LOCALGROUPDISK migration (user-approved).

**Out:** downstream NTuple-processing of the new sample (new DSID cross-sections, dropping the
isospin weight, `FullSimSampleType.h` wiring) — a follow-up; the HIJING-overlay full sample
(still in production).

## Design Decisions

### D1 (2026-07-13). The 6 full-sample NTUPs live on LOCALGROUPDISK, not on GPFS.

**Old approach (Steps 1–5 as originally planned):** `grid_monitor --mode fullsim_pp_full`
downloads each task, hadds it into one canonical NTUP per slice on GPFS, validates entry counts.

**New approach (USER DECISION 2026-07-13):** replicate the 6 grid output datasets **directly to
`BNL-OSG2_LOCALGROUPDISK`** (`rucio add-rule` — no resubmission needed, the outputs already exist
on the grid), build a **symlink farm**, and point the NTuple processing at the farm. **Download
only ONE small slice locally** — `pTH70_125` (~11 GB) — as a dev/test sample.

**Reason (two, both decisive):**
1. **Space.** The measured NTUP size is **30.6 kB/event** (task 51419511: 36.723 GB for 1.2 M
   events), not the ~18 kB/event extrapolated from the 10 k-event test-sample files. The 6 slices
   therefore total **≈ 280 GB**, not ~180 GB. The LGD migration frees ~292 GB — downloading would
   consume essentially all of it and put us straight back at the quota wall.
2. **The overlay sample forces this architecture anyway.** r17618 overlay NTUPs are 22 GB per
   10 000 events; a full overlay sample (~10 M events) would skim to **~20 TB**, an order of
   magnitude beyond the entire GPFS quota. There is no version of the future in which that sample
   lives on GPFS. Build the read-from-LGD path now, on the smaller sample, where it is cheap to
   get right.

**Why this is cheap — three things already exist in the code:**
- `fullsim_input_dir_override` (`PythiaAlgCoreT.h:258`, used at `PythiaAlgCoreT.c:27-28`) ⇒
  pointing the reader at the farm needs **no code change**, just a constructor argument.
- Reading ROOT over the dCache/pnfs POSIX mount is **already proven in this analysis**:
  `PythiaAlgCoreT.c:50` has a `pythia_pnfs_dir` branch behind `getUseLocal()`, and the June-2026
  migration (`localgroupdisk_migration.md`) moved 1.4 TB — including `dimuon_data` — to LGD and
  repointed `~/dcachearea`.
- `only_pp_isospin` (`PythiaAlgCoreT.c:382-389`) already exists ⇒ the pp-beam-only full sample is
  supported without new flags. (The AMI cross-sections for DSIDs 803015–803020 still must be
  re-read; see §1.)

**The ONE code change required.** `PythiaAlgCoreT.c:391-404` (`InitInputFullsim`) assumes
**exactly one file per (slice, beam)**: it builds the literal name
`Pythia_5p36TeV_<beam>_hQCD_DiMu_pTH<lo>_<hi>.FullSimPP24.NTUP.root`, tests it with
`std::ifstream::good()`, then `ch->Add(fname)`. Each grid task emits **12 files** (measured on
803018), and with no local download there is no hadd to collapse them into one. Fix: the reader is
**already a `TChain`**, so name the farm symlinks `...NTUP.part01.root`, `...part02.root`, … and
switch to `ch->Add(dir + "...NTUP.part*.root")` (`TChain::Add` accepts globs), replacing the
`ifstream` existence test with a glob-count test. Small, local, and must go through
`/review-analysis-code` (it changes how MC input is assembled — a provenance-critical path).

**Accepted cost:** every NTuple-processing pass streams ~280 GB over dCache rather than GPFS, so
each pass is slower. Judged worth it against a permanently-clear quota and the 20 TB overlay
problem. The local `pTH70_125` slice exists so that iteration/debugging does not pay that cost.

## Implementation Plan

1. [ ] **Storage census + candidate list (Task II)** — `mmlsquota` + `du -sb`; deliver the
   relocation list; **wait for user decision**. (§4)
2. [ ] **Run-dir reorganisation** — `run_pythia_fullsim` → `run_pythia_fullsim_test_sample`;
   create `run_pythia_fullsim_full_sample`.
3. [ ] **Local test run** on one full-sample AOD; verify every item of §3. → `/review-analysis-code`
   if any config changed.
4. [ ] **Grid submission** of the 6 DSIDs (`grid_sub.sh`, new VER_TAG).
5. [ ] **grid_monitor mode** for the full sample (`--mode fullsim_pp_full` → DATA_BASE
   `pythia_fullsim_full_sample`); launch monitoring.
6. [ ] **Migration (Task II)** — subagent, migration skill, only if the user approves.
7. [ ] **Download + validate + bookkeeping**; post-skim sanity check.

## Progress Log

*(append-only)*

- 2026-07-13 — **Doc created.** Doc triage: read `INDEX.md` + CLOSED `mc_trigger_info_skim.md`
  (the direct predecessor). Full-sample DIDs resolved from rucio (§1 table): 6 DSIDs
  803015–803020, `_pdf`, **pp beam only**, 987 files / 9.87 M events / 7.26 TB AOD.
  Storage census done (§4 + Results R1): **1471 GiB used of a 1536 GiB soft quota — 65 GiB
  headroom, vs ~250 GB needed** ⇒ Task II blocks Task I's download stage.

- 2026-07-13 — **User decision on Task II (all 4 candidates approved, MOVE semantics).**
  Approved: (A) overlay trigger-off `.bak_20260709` × 6 = 133 GB; (B) `pythia_private_sample`
  raw `*_all_k` = 113 GB; (C) POWHEG fullsim; (D) `dimuon_data` test AODs + old backups = 20 GB.
  **User instruction on (C) partially falsified by the facts — recorded here so it is not lost:**
  the user asked to check that the two `user.yuhang.*._MYSTREAM/` dirs hold *one file each*, and
  if so (and if >10 GB) delete them and drop the `useLocal` option from the POWHEG fullsim
  NTuple-processing code. **They do NOT.** `bb` = 28 files, `cc` = 35 files; each holds the raw
  grid NTUP (`user.yuhang.4859*.MYSTREAM._000001.root`, 2.19 / 2.20 GB → **4.4 GB total, not
  >10 GB**) **plus all the derived NTuple-processing output** (`muon_pairs_powheg_*`,
  `single_muon_trees_powheg_*`, `hists_powheg_*`, ~9 GB each dir). Both legs of the condition
  fail ⇒ **the `_MYSTREAM` dirs are OUT OF SCOPE: not deleted, and `useLocal` is NOT touched.**
  Only `powheg_full_sample/mixed` (13.8 GB) + `mixed_mass_1_3GeV` (13.0 GB) are migrated.
  Revised total to free: **~293 GB** (still > the ~250 GB needed).

- 2026-07-13 — **Migration tooling: the `bnl-localgroupdisk` plugin (user-supplied).**
  `https://github.com/FlamyFlame/claude-bnl-localgroupdisk` — installed to
  `dimuon_codes/.claude/plugins/bnl-localgroupdisk` (skills: `preflight`, `migrate`,
  `check-rule`, `build-symlinks`). **This supersedes the hand-rolled recipe** and changes the
  design: the correct procedure is upload → `BNL-OSG2_SCRATCHDISK` → `rucio add-rule` →
  `BNL-OSG2_LOCALGROUPDISK` → **symlink farm with a same-path swap** (original dir renamed
  `<dir>_orig`, symlinks placed at the ORIGINAL path pointing into
  `/pnfs/usatlas.bnl.gov/LOCALGROUPDISK/rucio/user/...`). Analysis code therefore needs **no path
  edits**; the quota is freed by deleting `<dir>_orig` *after* the farm is verified. FTS
  replication can take 1–12 h. Delegated to a subagent (scratch doc `_sub_lgd_migration_1.md`).

- 2026-07-13 — **Steps 2–3 DONE: run-dir reorg + local test run PASSES all §3 checks.**
  Reorg: `SkimCode/run_pythia_fullsim` → `run_pythia_fullsim_test_sample`;
  `run_pythia_fullsim_full_sample` created (copy of the build + `TrigRates_CA.py`, new
  `grid_sub.sh`, old test-sample artefacts removed). Run dirs are git-excluded
  (`.git/info/exclude: SkimCode/run*/**`), so this is not a tracked move.
  Test AOD staged **outside** the exhausted data fileset, in `/usatlas/scratch/yuhanguo/`
  (`usatlast3-scratch`, ~2 TB free): `AOD.50094090._000110.pool.root.1` (8.21 GB, DSID 803019
  pTH70_125), via `rucio download --nrandom 1`.
  Job: `TRIGRATES_RUNMODE=ppmcfullsim2024 athena.py --no-excabort TrigRates_CA.py --evtMax=300`
  → **rc=0**, output `test_fullsample_pp24.root` (9.77 MB).

  | §3 item | result |
  |---|---|
  | (f) entries == events read | **300 == 300** ⇒ `StoreAllEvents=True` reached the job; NO trigger filtering |
  | branch list vs test-sample NTUP | **260 vs 260, ZERO difference either way** (the decisive check) |
  | (a) reco kinematics | 687 muons; `muon_pt` correctly SIGNED (336 neg / 351 pos); `muon_eta` ∈ [−2.47, 2.55]; `muon_quality` ∈ [66, 511]; `muon_trk_pt`, `muon_d0` populated |
  | (b) truth kinematics | `truth_muon_pt` (n=668), `muon_truth_pt` (n=687) populated |
  | (c) truth origin / provenance | 337 279 truth particles; `truth_barcode` max 7.2e6; `truth_parents` 211–2491 links/event (mean 1124); `muon_truth_prob` mean 0.971, **670/687 muons matched at prob>0.5**; `muon_truth_IsPrimary` populated |
  | (d) scale factors / weights | `muon_eff_SF_medium` mean 0.678, `muon_eff_SF_tight` mean 0.639, both ∈ [0, 1.036] (SF=0 = muon fails that WP, expected); `EventWeights` ≡ 1 (Pythia unweighted; the per-slice σ·ε_filt weighting is applied downstream from AMI, in **nb**) |
  | (e) trigger simulation | `b_HLT_mu4_L1MU3V` 281/300, `b_HLT_mu4_mu4noL1_L1MU3V` 232/300, `b_HLT_2mu4_L12MU3V` 165/300; **2mu4 ⊆ mu4 EXACTLY** (0 violations); per-leg `dimuon_b_HLT_2mu4_..._mu1passLeg1_dR_0_02` populated in 160 events ⇒ the ΔR-correlation ingredient is present |

  Trigger rates (94 % / 77 % / 55 %) are exactly where the pT-hat trend predicts for the
  `pTH70_125` slice (cf. the overlay slice table in `mc_trigger_info_skim.md`: rates rise
  monotonically with pT-hat; the pp `pTH8_14` slice is 80 / 52 / 30 %).

  **Benign warning noted, not a regression:** `MuonEfficiencyTool_{Medium,Tight}` emits
  *"Failed to find the RandomRunNumber decoration … You'll receive SFs from the most recent
  period."* The config is byte-identical to the one that produced the validated test sample, so
  this is pre-existing behaviour, not something introduced here: the muon SFs are
  period-inclusive rather than period-weighted. Flagged for the SF systematic; not blocking.

  ⇒ **The full sample carries everything the analysis needs.** Nothing in the test-sample NTUP is
  missing from it.

- 2026-07-13 — **Step 4 DONE: all 6 grid tasks SUBMITTED, rc=0, zero submission errors.**
  `run_pythia_fullsim_full_sample/grid_sub.sh`, mode `ppmcfullsim2024`, `VER_TAG=FullJuly2026.v1`
  (distinct from the test sample's `July2026.v1` — pathena task names are immutable).
  outDS pattern `user.yuhang.NTUP.Pythia_5p36TeV_pp_hQCD_DiMu_<slice>.FullSimPP24.FullJuly2026.v1.`

  | DSID | slice | jediTaskID | files → jobs (`nFilesPerJob`) |
  |---|---|---|---|
  | 803020 | pTH8_14 | **51419497** | 203 → 41 (5) |
  | 803016 | pTH14_24 | **51419500** | 360 → 72 (5) |
  | 803017 | pTH24_40 | **51419506** | 240 → 48 (5) |
  | 803018 | pTH40_70 | **51419511** | 120 → 24 (5) |
  | 803019 | pTH70_125 | **51419515** | 32 → 8 (4) |
  | 803015 | pTH125_300 | **51419520** | 32 → 8 (4) |

  `nFilesPerJob` is set EXPLICITLY — never `--nGBPerJob MAX`, which with `--mergeOutput` + remote
  XRootD I/O silently lost events in the Apr-2026 data skim
  (`skim_discrepancy_investigation.md`).

  **Gotcha for the future:** `SkimCode/setup_25.sh` does a *relative* `cd build_25` and then
  `acm compile`. It MUST be sourced with cwd = `SkimCode/`, else it silently does nothing and
  `pathena` is never on PATH. (Cost a failed submission attempt.)

- 2026-07-13 — **Step 5 DEFERRED ON PURPOSE: `grid_monitor` NOT launched yet.**
  `--mode fullsim_pp_full` added to `SkimCode/scripts/grid_monitor.sh` (DATA_BASE =
  `~/usatlasdata/pythia_fullsim_full_sample`, `IS_MC_FLAT=1`, its own
  `merging-record.txt`/state/log). But `grid_monitor` **downloads as soon as a task completes**,
  and with only 65 GiB headroom that download would fail on the quota — precisely the 2026-07-09
  incident (`mc_trigger_info_skim.md`): a quota-blocked `rucio download` fails on every RSE in
  <1 s, `grid_monitor` marks the healthy task `failed`, and it has ALREADY renamed the canonical
  NTUP to `.bak`. Launch the monitor only once the migration has freed the space
  (target: ≥250 GB headroom).

### D2 (2026-07-13/14). Isospin: `isTestSample` is the single switch. **A REAL BUG WAS FOUND.**

**USER-SUPPLIED PHYSICS RULE (authoritative):**
- **pp-CONDITIONS fullsim simulates pp COLLISIONS** ⇒ **ONE** isospin beam (pp), isospin weight
  **1**. There is nothing to isospin-average. The pp24 **TEST** sample having 4 beams was a
  **PRODUCTION MISTAKE**; the FULL `_pdf` sample correctly has only the pp beam.
- **PbPb-CONDITIONS HIJING overlay simulates Pb+Pb**, whose nucleons are a p/n mix ⇒ **FOUR**
  beams {pp,pn,np,nn} combined with the Pb ratio **4:6:6:9** (Z=82, N=126 ⇒ (2:3)⊗(2:3), Σ=25/25).
- The two TEST samples are exceptions in **opposite** directions ⇒ the rule collapses to a XOR:
  `four_beams = FullSimSampleIsOverlay(t) != isTestSample` (`FullSimSampleType.h`).

**Implementation:** ONE public switch `isTestSample` (`PythiaAlgCoreT.h`, **default `false` = the
FULL production**) drives **input dir + AMI cross-section dir + isospin treatment** together, so
they can never disagree. `setIsospinBeams()` remains only as an escape hatch (used by the
`pp_only` cross-check). All 13 run scripts declare `isTestSample`; the full-sample script sets
nothing (the default IS the full sample). Verified by compiled test macro — all 5 cases correct.

**THE BUG (pre-existing, now fixed):** `fullsim_weight_factor = ami_w * nom_ratio / N_beam`
(`PythiaAlgCoreT.c`) used `nominal_beam_ratio["pp"] = 4/25` **even when only the pp beam was
read** — which was EVERY overlay run (`FullSimSampleIsOverlay` forced pp-only) and every
`pp_only` run. A spurious factor **4/25 = 0.16**.
- **Where it CANCELS (and why):** it is **slice-independent**, and it sits in **both numerator and
  denominator**. ⇒ reco efficiency, detector response, MC trigger efficiency (all steps),
  area-normalized templates: **unaffected**.
- ⚠ **Correction to an earlier claim in this doc/session:** reco-eff and det-response are **NOT
  filled unweighted**. An empty weight *specifier* resolves to the `weight` **column**
  (`RDFBasedHistFillingPythia.cxx:9` maps `""` → `"weight"`). They ARE weighted; the factor
  cancels because it is slice-independent, not because it is absent. **This distinction is
  load-bearing:** an AMI cross-section error is *slice-dependent* and therefore cancels
  **nowhere**.
- **Where it does NOT cancel:** any **absolute cross-section**. Overlay NTP outputs and
  `plot_pythia_fullsim_overlay_kn_pt_crossx.cxx` are **6.25× low and STALE** ⇒ **rerun pending**
  (see Remaining Work).
- The pp24 **test**-sample outputs are **UNCHANGED** (4 beams + 4:6:6:9 before and after), so
  `plot_pythia_fullsim_kn_pt_crossx.cxx` does **not** need a rerun for this fix.

**USER DECISION (2026-07-14): "ANY crossx plot must be honest — this is a key physics
observable."** A cross-section from the 4-beam pp24 **TEST** sample carries the **Pb isospin
average** and is therefore **NOT a physical pp cross-section**. It must be labelled as such.
The honest absolute pp σ comes from the FULL sample (pp-only, weight 1). **Plot labelling is
still TODO** (see Remaining Work).

### D3 (2026-07-14). AMI provenance guard — **the landmine the reviewer caught.**

AMI files are named by **beam + slice ONLY**, so they do **not** identify the production. The
default AMI dir is the *truth* sample's `ami_info/`, holding DSIDs **802758–802781**. The full
sample is **803015–803020** with **different, slice-dependent** cross-sections:

| slice | TEST σ·ε_filt (nb) | FULL σ·ε_filt (nb) | ratio |
|---|---|---|---|
| pTH8_14 | 23.78 (802781) | **36.7432** (803020) | **1.55** |
| pTH14_24 | 31.89 | **35.6365** (803016) | 1.12 |
| pTH24_40 | 19.34 | **18.5418** (803017) | 0.96 |
| pTH40_70 | — | **6.9018** (803018) | — |
| pTH70_125 | — | **1.4588** (803019) | — |
| pTH125_300 | — | **0.1978** (803015) | — |

The ratio spans **0.96–1.55** ⇒ it is **slice-dependent** ⇒ it does **NOT cancel in any ratio**.
Running the full sample against the old AMI would have silently corrupted the MC trigger
efficiency **and** every cross-section, with **no error message** (the file exists under the
right name).
**Fixed:** AMI written to `~/usatlasdata/pythia_fullsim_full_sample/ami_info/`; new
`ami_info_dir_override` + `expected_ami_dsids`; `InitInputFullsim` parses `datasetNumber` and
**THROWS** on mismatch. **Guard proven to fire** (pointed at the test sample while declaring the
full-sample DSIDs → threw `AMI PROVENANCE MISMATCH`). All 13 run scripts now declare their DSIDs.

Also made **fatal** (were silent): ambiguous input (a hadded file *and* `.part*` for one slice —
a stale local file would shadow the farm and give a wrong `N_beam` with an unchanged σ); a
**missing pT-hat slice** (biases the σ-weighted combination — unlike a *partial* farm, which is
self-correcting because `N_beam` is measured from the files actually chained); zero entries.
Single-slice **diagnostic** runs (the r-tag dirs hold `pTH8_14` alone) opt out with
`allow_missing_slices = true`.

- 2026-07-14 — **AMI weights fetched + verified; the rule written into permanent space (user request).**
  `pyami` re-queried for **all 24** test DSIDs (802758–802781) **and all 6** full-sample DSIDs
  (803015–803020): **zero drift** vs what is on disk. New authoritative registry
  **`Analysis/docs/ami_weights.md`** (per-DSID σ, genFiltEff, σ·ε_filt; the rule; the code guard;
  the rerun blast radius). **BLOCKING rule added to `.claude/CLAUDE.md`** and to
  `.claude/conventions/ntuple-provenance.md` (reviewer enforcement). Commit `be9b77a`.
  > **THE RULE:** AMI weights are a HARD BLOCK on every MC dataset. Moving to a NEW dataset
  > (test→full, new overlay, new generator, new tag, any new DSID) ⇒ fetch that dataset's OWN AMI
  > weights from `pyami` FIRST. **NEVER reuse the old ones — it is a silent failure that
  > propagates to final results.**
  Why silent: AMI files are keyed by **beam+slice only**, so the filename is byte-identical across
  productions — swap the dataset, leave the old `ami_info/` in place, and everything runs and every
  number is wrong. Why it cancels nowhere: the error is **slice-dependent** (pp24 FULL/TEST σ·ε
  ratio spans **0.879–1.545**), so it reweights the pT-hat mixture and survives every ratio,
  including the MC trigger efficiency (a σ-weighted average over slices).

- 2026-07-14 — **Overlay crossx rerun (minimal) + TWO honesty bugs found and fixed. Commit `92d9f06`.**
  - **NTuple processing did NOT need a rerun.** The on-disk overlay NTP output (mtime 07-14 01:17)
    was verified against the AMI table for **all 6 slices**: mean `weight` = σ·ε/N_beam with
    isospin ratio **1.0000** in every slice ⇒ it is already the post-fix (D2) output. Only the
    **plot** was stale. (This is the "only rerun what is needed" answer.)
  - ⚠ **BUG 1 — the y-axis was mislabelled by 1000×.** Both crossx plots printed
    `d#sigma/dp_{T} [#mub/GeV]` while applying **no unit conversion anywhere**. The weight carries
    the AMI cross-section in **nb**, so the values are **nb/GeV**. Verified numerically: overlay
    single-b total σ = **13.12 nb** (kin0 1.741, kin1 4.867, kin2 4.163, kin3 1.886, kin4 0.413,
    kin5 0.049 nb), and the plotted kin0 integral matches **1.74 nb**, not 1.74 μb. Axis corrected
    to **[nb/GeV]** on both plots (values unchanged) + a "do not fix this back" units comment.
  - ⚠ **BUG 2 — the pp crossx plot was not honest** (USER: *"ANY crossx plot must be honest — this
    is a key physics observable"*). It reads the pp24 **TEST** sample, which was produced with 4
    isospin beams **by mistake** and is therefore combined with the **Pb 4:6:6:9 isospin average**.
    A pp-conditions sample has nothing to isospin-average ⇒ that absolute σ is **NOT a physical pp
    cross-section**. Both panels now carry an explicit red caption saying so. The honest pp σ will
    come from the FULL sample (pp beam only, weight 1).
  - Regenerated: 16 overlay plots (`.../pythia_fullsim_hijing_overlay_test_sample/plots/`) + 12 pp
    plots (`.../pythia_fullsim_test_sample/plots/`), all `rc=0`. Shape check: smoothly falling
    dσ/dp_T, each pT-hat slice peaking in its own range — physically sane.
  - Backup of the pre-rerun overlay NTP kept:
    `muon_pairs_...pbpb23_no_data_resonance_cuts.bak_pre_isospin_fix_20260714.root`.

- 2026-07-14 — **LGD farm LIVE for the first slices — the read-from-LGD architecture (D1) WORKS.**
  `fullsim_pp24_full_to_lgd.sh` running detached. Grid: 3/6 tasks `done` with **zero failed files**
  (pTH8_14 203/203, pTH24_40 240/240, pTH40_70 120/120); the other 3 still scouting/running.
  Rules created (persisted in `pythia_fullsim_full_sample/lgd_rules.txt`) and **replication was
  fast, not the feared 1–12 h** (BNL→BNL is local):

  | slice | rule | parts | entries read **through the symlink farm** |
  |---|---|---|---|
  | pTH24_40 | `d2fa6f50473f4c37…` | 17 | **2 399 995** |
  | pTH40_70 | `450f1755c68a45bd…` | 12 | **1 199 997** |
  | pTH8_14 | `33bd318b47dd4bd2…` | — | replicating |

  ⇒ ROOT reads the unmerged grid outputs through the pnfs symlink farm exactly as designed; nothing
  is hadded onto GPFS.

- 2026-07-14 — **⚠ Small event deficit in the skim: ~2–3 events lost per MILLION (2e-6).**
  pTH24_40: NTUP 2 399 995 vs AOD 2 400 000 (**−5**). pTH40_70: 1 199 997 vs 1 200 000 (**−3**).
  **Not a grid failure** — BigPanDA reports 240/240 and 120/120 files finished, `nfilesfailed=0`,
  `neventsTot` = the full 2.4 M / 1.2 M. So the events are dropped **inside the skim**, despite
  `StoreAllEvents=True`; most likely a handful of events taking an early return in `TrigRates.cxx`
  (e.g. no primary vertex). **Not root-caused.**
  **Why the test sample never showed it:** at 10 000 events the expected loss is 0.02 events — the
  10 000/10 000 match in `mc_trigger_info_skim.md` is fully consistent with a 2e-6 rate.
  **Physics impact: negligible, and it does not bias.** `N_beam` is measured from the files
  actually chained, so `w = σ·ε_filt / N_stored` stays self-consistent, and the dropped events are
  an unbiased subset. 2e-6 is orders of magnitude below any systematic here. **Recorded, not
  chased.** Revisit only if a slice ever shows a loss ≫ 1e-5 (that would indicate a real problem,
  e.g. a truncated/failed merge, not this).

---

## TASK III (added 2026-07-14) — push the FULL sample through the whole pp-fullsim chain

**Request:** after the grid jobs finish, run the pp24-fullsim FULL sample through the entire pp
fullsim chain (smoke-test first; back up anything that would be overwritten) and reproduce every pp
fullsim plot — **reco efficiency, detector response, the new MC-based trigger efficiency** — plus a
**statistics plot: pair-pT differential crossx for each pT-hat slice on one canvas**.

### III.0 Code exploration (DONE 2026-07-14) — the chain, and what is missing

**The chain** (there is **NO pp-fullsim pipeline script**; the overlay one,
`pipelines/pipeline_pythia_fullsim_overlay.sh`, is the template):

| stage | code | note |
|---|---|---|
| NTP nominal | `run_pythia_fullsim_full_sample.sh` | **exists** ✅ → `..._no_data_resonance_cuts_full.root` |
| NTP `_single_muon` | — | ❌ **MISSING for full sample** |
| NTP `_mc_trig` | — | ❌ **MISSING** |
| NTP `_mc_trig_single_muon` | — | ❌ **MISSING** |
| RDF hists | `RDFBasedHistFillingPythiaFullsim` | needs `is_test_sample=false` |
| reco-eff + det-response | `PythiaFullsimRecoEffPlotter` + `plot_reco_effcy_pythia_fullsim_pp24.cxx` | needs `is_test_sample=false` |
| single-muon reco-eff | `plot_single_muon_reco_effcy.cxx` | needs the `_single_muon` NTP |
| MC trig-eff | `FillMCTrigEffHists` → `FitMCSinglesEffcy` → `FillMCTrigEffHists(step3)` → `plot_mc_trig_eff.cxx` | **dir HARD-CODED to test sample**; no knob |
| crossx per pT-hat slice | `plot_pythia_fullsim_kn_pt_crossx.cxx` | **paths hard-coded**; see below |
| kn contributor table | `make_kn_contributor_table.cxx` | hard-coded paths |
| reco distr single-b vs OS | `plot_reco_distr_singleb_vs_op_pp24.C` | hard-coded paths |

**⚠ The user's "statistics plot" ALREADY EXISTS** — it is the LEFT panel of
`plot_pythia_fullsim_kn_pt_crossx.cxx` (`plot_impl`): dσ/dp_T markers, one colour per pT-hat slice,
all on one canvas. What is new with the full sample is that its **error bars become real** instead
of a forecast. (Plus the right-panel stack, and the `err_fraction` / `err_ratio` diagnostics.)

### III.1 BLOCKERS found by the exploration (must be fixed before any full-sample run)

1. **3 NTP run scripts missing** for the full sample (`_single_muon`, `_mc_trig`,
   `_mc_trig_single_muon`). Without them the full sample produces **no MC-trigger-efficiency inputs
   at all** and no single-muon tree.
2. **`RDFBasedHistFillingPythia.h:114` `is_test_sample = true`** (default) — must be set false, else
   the RDF stage silently reads the TEST sample.
3. **`PythiaFullsimRecoEffPlotter.cxx:111-114` filename bug:** the base `GetInputFilePath()` builds
   `histograms_..._no_data_resonance_cuts.root` **without the `_full` suffix**, while the RDF filler
   writes `..._full.root`. Flipping `is_test_sample=false` changes the *directory* but not the
   *filename* ⇒ open fails (or silently picks a stale file).
4. **`FillMCTrigEffHists.cxx:71-74`, `FitMCSinglesEffcy.cxx:101`, `plot_mc_trig_eff.cxx:211`** —
   the sample dir is a **hard-coded string** to `pythia_fullsim_test_sample/`. No knob exists; one
   must be added.
5. **`plot_pythia_fullsim_kn_pt_crossx.cxx`** — input + output paths hard-coded to the test sample
   in 4 functions.
6. ⚠ **`kScaleFactors` / `scale_factors` would DOUBLE-COUNT the statistics.** They are
   `N_full/N_test` forecast factors applied as `SetBinError(err/sqrt(sf))` (`:288`, `:487`, `:675`;
   overlay `:195`, `:306`, `:407`). **Central values are untouched** (no `SetBinContent`), so with
   the real full sample the *error bars* would be divided by √sf a SECOND time — understated by a
   further **7.1× / 9.5× / 7.7× / 5.5× / 2.8×** (kn0–kn4), **slice-dependently** ⇒ does not cancel,
   and it fails **silently** (the σ curve still looks right). The `err_fraction` / `err_ratio`
   "which slice is the statistical bottleneck" maps would rank slices by a fictitious ratio.
   ⇒ **Delete the forecast factors** (keep the functions — they are genuinely useful once fed real
   errors).
7. **The honesty caption must be REMOVED for the full sample.** `plot_pythia_fullsim_kn_pt_crossx`
   currently prints *"TEST sample … Pb isospin avg (4:6:6:9) … NOT a physical pp σ"*. On the FULL
   sample (pp beam only, isospin weight 1) that warning becomes **false** — the full-sample σ IS a
   physical pp cross-section.
8. **WP: the crossx macro hard-codes `pair_pass_medium`** (`:64, 258, 456, 644`) while the analysis
   nominal is **TIGHT** and every other pp-fullsim stage already defaults to Tight. Violates
   `feedback_plots_wp_config_var` / `muon_wp_registry.md`.
9. **`docs/muon_wp_registry.md` is STALE:** §4 claims the reco-eff plotters default to *Medium*;
   the code defaults to **Tight** (`PythiaFullsimRecoEffPlotter.cxx:28,627`). §1 line anchors drifted.
   The crossx macro is **absent from the registry** despite hard-coding a WP on a key observable.

### III.2 USER DECISIONS (2026-07-14)

- **Scope = produce + validate the plots ONLY.** Do **NOT** wire the full-sample results into
  crossx/R_AA: the pp reco-eff stays the Run-2 placeholder (`project_pp_reco_eff_placeholder`) and
  ε_ΔR stays the ≡1 dummy until the user has seen the new plots. Keeps the blast radius small.
- **⚠ CONCURRENCY — a SIBLING SESSION owns the MC-trig-eff chain right now.**
  `mc_trigger_efficiency.md` is being actively edited (mtime 2026-07-14 03:04); it has already
  rerun both samples on all 24 pp slices and has an **open physics decision** (§2 assumption (i)
  VIOLATED — the single-muon efficiency IS ΔR-dependent; pp/2mu4 product SAFE, **PbPb mu4 UNION
  weight AT RISK**). Its merge is HELD.
  ⇒ **Do NOT run the MC-trig-eff chain, and do NOT edit `FillMCTrigEffHists.cxx` /
  `FitMCSinglesEffcy.cxx` / `plot_mc_trig_eff.cxx`, until that session has committed.** Those are
  its files; concurrent writes clobber.
  **After the sibling lands AND all grid jobs finish:** back up the TEST-sample MC-trig-eff plots,
  then remake them from the FULL sample.
- **Stale-warning CORRECTED:** the `INDEX.md` note that *"`8c917b4` changed the pp isospin default
  without its companion runner fix ⇒ committed pp runners silently drop 3/4 of pp stats"* was true
  at `8c917b4` but is **STALE at HEAD**. Verified 2026-07-14: all four pp TEST runners set
  `isTestSample = true` (⇒ 4 beams), added in `216f750`/`1b927d2`, which land *after* `8c917b4`.
  No statistics loss at HEAD. (Also confirmed by a compiled test: pp + `isTestSample=true` →
  `four_beams=1`.)
- **Everything else: proceed.** Fix the silent failures, add config variables, add the missing NTP
  scripts, and write a **pp-fullsim pipeline** (with the MC-trig-eff stage present but
  **switchable off**, default OFF until the sibling lands). Only genuinely
  physics-ambiguous items get escalated.

## Autonomy Contract — TASK III **DONE 2026-07-21** (both tasks, both WPs, validated)

**COMPLETE.** All pp24-fullsim results reproduced on the FULL sample (9.87 M events) at BOTH
working points and physics-validated:
- **Task 1 (non-trig-eff):** reco efficiency, detector response, single-muon reco efficiency,
  reco distributions, differential crossx per pT-hat slice (the "statistics" plot) + kn table —
  Tight (47 plots) + Medium (reco-eff 25, det-resp 6, single-mu 3). Crossx/kn/reco-distr are
  single-WP nominal (Tight) by design (no tight/medium output separation).
- **Task 2 (trig-eff):** MC-based trigger efficiency (step1 singles data/MC, step2 ΔR-binned incl.
  L1/HLT/full_chain, step3 ε_ΔR, step4 singles ε_ΔR) — Tight + Medium. ε_ΔR^2mu4 plateau
  0.988±0.001 (Tight) / 0.9875±0.0012 (Medium), ≈1 with the expected small-ΔR suppression.
  Test-sample plots backed up first (`mc_based_TESTSAMPLE_backup_20260720` + timestamped `.bak`s).
- **Bug found + fixed mid-run:** `store_mc_trigger` threw on the multi-file farm (ROOT exits 0 on
  the swallowed throw ⇒ silently skipped the mc_trig NTP); fixed with a per-file trigger-branch
  uniformity check (`31c2724`), pipeline Stage-4 now validates the mc_trig NTP.
- **Not adopted into crossx/R_AA** (per the scope decision): pp reco-eff stays the Run-2
  placeholder, ε_ΔR stays the ≡1 dummy until the user reviews these plots.
- **Concurrency:** the parallel r17663 session's edits to `plot_mc_trig_eff.cxx` were insulated
  (pp_full block byte-identical; r17663 code path is `noovl`-only). Its WIP left uncommitted for it.

Commits: `c700b0c a4462db e01e374 dcbc434 31c2724 206f127 e3cfbef` (+ doc commits).

## Autonomy Contract (SUPERSEDED — TASK III RE-SCOPED 2026-07-20; see DONE block above)

**Context at re-scope:** the concurrent `mc_trigger_efficiency` session has **merged to master**
(`8edc4fb`) and been stopped; the farm is **COMPLETE (all 6 slices, 9 869 980 events)**. The user
now directs: **(1) reproduce ALL non-trigger-efficiency pp-fullsim results on the full sample;
(2) back up the pp trigger-efficiency results made from the TEST sample and reproduce them from the
FULL sample.** Explicit negative constraint: **do NOT touch the PbPb (overlay) results or the
r17663 no-overlay results.** Proceed autonomously; stop on ambiguity.
- Done =
  1. **Non-trig-eff (task 1):** full-sample NTP (nominal + single_muon) → RDF hists → reco-eff +
     detector-response + single-muon reco-eff + crossx/statistics + kn table + reco-distr, all
     regenerated from the FULL sample; `/review-plot`.
  2. **pp trig-eff (task 2):** add a pp-FULL-sample option to the MC-trig-eff chain
     (`FillMCTrigEffHists`/`FitMCSinglesEffcy`/`plot_mc_trig_eff`) WITHOUT changing the overlay or
     r17663 configs; run the 2 full-sample `_mc_trig` NTP scripts; **back up the TEST-sample pp
     trig-eff plots**; reproduce from the FULL sample; `/review-plot`.
  3. Orphaned WIP from the stopped session committed; `muon_wp_registry.md` +
     `docs/pythia_fullsim_pp.md` updated; `/review-analysis-code` on the code.
- Stop-and-ask = ANY physics-results-bending ambiguity (no fixed list; use judgment;
  when unsure whether an ambiguity is blocking, treat it as blocking → AskUserQuestion).

## Autonomy Contract (SUPERSEDED 2026-07-20 — original Task III below, kept for the record)
- Mandate: run autonomously to DONE; do NOT pause to confirm progress. Finishing a plan, a passing
  smoke test, or one pipeline stage is NOT a stopping point.
- Done =
  1. Blockers 1,2,3,5,6,7,8,9 fixed (NOT 4 — that is the sibling's file; deferred).
  2. 3 new full-sample NTP run scripts (`_single_muon`, `_mc_trig`, `_mc_trig_single_muon`).
  3. Every pp-fullsim stage/plot parameterized so it can read EITHER sample (no hard-coded
     test-sample paths), with a WP config var defaulting to TIGHT.
  4. `pipelines/pipeline_pythia_fullsim_pp.sh` written (MC-trig-eff stage switchable, default OFF).
  5. Smoke test of each stage on the farm before the full run.
  6. FULL-sample run: NTP → RDF hists → reco-eff + det-response + single-muon reco-eff →
     crossx/statistics plots + kn table. **Back up anything that would be overwritten.**
  7. `/review-analysis-code` on the code; `/review-plot` on the regenerated plots.
  8. `muon_wp_registry.md` + `docs/pythia_fullsim_pp.md` updated.
  9. **AFTER the sibling session commits AND all 6 grid tasks finish:** back up the TEST-sample
     MC-trig-eff plots and remake them from the FULL sample (blocker 4 + the MC-trig-eff run).
- Stop-and-ask = ANY physics-results-bending ambiguity (no fixed list; use judgment;
  when unsure whether an ambiguity is blocking, treat it as blocking → AskUserQuestion).

### III.3 Plan

1. **Code prep (no data needed — do NOW, in parallel with the grid):** blockers 1,2,3,5,6,7,8;
   the 3 NTP scripts; the pipeline; the WP config var; registry/doc updates.
2. **Smoke test** each stage on the farm with a small `nevents_max`.
3. **Full run** once the grid + farm are complete (back up first where a clobber is possible).
4. **MC-trig-eff (LAST):** only after the sibling session commits — blocker 4 + rerun + backup.

- 2026-07-20 — **TASK III execution (re-scoped): non-trig-eff chain + pp trig-eff on the FULL sample.**
  Context: sibling `mc_trigger_efficiency` MERGED to master (`8edc4fb`) and stopped; farm COMPLETE
  (6/6 slices, **9 869 980** events). Committed work:
  - `c700b0c` — finished the 3 orphaned plot-macro parameterizations left by the stopped session
    (reco-eff driver, kn table, reco-distr); fixed two it left broken (missing
    `FullSimSampleType.h` include; `.C` is interpreted not ACLiC).
  - `a4462db` — **`pp_full` sample knob** added to the MC-trig-eff chain (`FillMCTrigEffHists`,
    `FitMCSinglesEffcy`, `plot_mc_trig_eff`): same physics as `pp` (2mu4 PRODUCT weight, same pp24
    data reference), reads `_full` intermediate hists (label `pp24_full`), writes final plots to
    the canonical `pp_trigger_efficiency/mc_based/`. Overlay + r17663 configs untouched.
    Also **`set -e`/`set -u` hardening** of the pipeline (ALRB setup and ROOT both exit non-zero
    under `set -Eeuo pipefail` and were killing it silently — the same trap as the migration script).
  - `e01e374` — **BUG the smoke test caught:** the crossx macro's `g_is_test_sample`/`g_use_tight_wp`
    were `static` (internal linkage) ⇒ the ACLiC-compiled globals were INVISIBLE to the ROOT
    interpreter, so the pipeline's `g_is_test_sample=false` silently no-op'd and the FULL-sample
    crossx wrote to the TEST dir with the TEST caption. Fixed: dropped `static` + entry-point args.
    Pipeline: Stage 9 now actually runs `make_kn_contributor_table` (was logged, never invoked);
    Stage 8 adds reco-distr; **Stage 10 fully wired** for pp_full.
  - **Pipeline `pipelines/pipeline_pythia_fullsim_pp.sh`** (NEW) validated end-to-end by a
    2000-event smoke test: Stages 0–9 all ✅ (preflight → NTP×2 → validate → RDF → validate →
    reco-eff/det-resp → single-mu → crossx/kn). All 3 trig-eff files compile with the pp_full edits.
  - **Test-sample pp trig-eff plots backed up** → `pp_trigger_efficiency/mc_based_TESTSAMPLE_backup_20260720`
    (24 pngs) before any full-sample overwrite (user requirement).
  - **Full production run LAUNCHED** (detached, `USE_TIGHT_WP=1 ENABLE_MC_TRIG_EFF=1`): task 1
    (non-trig-eff, Stages 1–9) + task 2 (trig-eff, Stage 10) in one pass. Log
    `pythia_fullsim_full_sample/pipeline_run_tight_20260720_181419.log`.
  - **Still to do:** medium-WP pass (Stages 7–10, NTP/RDF are WP-agnostic → skip); `/review-plot`
    on the regenerated plots; `muon_wp_registry.md` + `docs/pythia_fullsim_pp.md` updates.

- 2026-07-20 — **CONCURRENCY (user-confirmed intentional): a parallel r17663 session is live.**
  PID 2654515 (`claude --resume`) is NOT stopped — it edited `plot_mc_trig_eff.cxx` at 18:35,
  adding an r17663 comparison overlay (HIJING r17618 on the r17663 step-1 panels). **User
  directive: it is intentional r17663 work, and r17663 must not affect pp.** Verified:
  - it touched ONLY `plot_mc_trig_eff.cxx` (my `FillMCTrigEffHists`/`FitMCSinglesEffcy` pp_full
    versions are untouched);
  - the edit only adds a `noovl`-only `cmp_fit_file` path (guarded `if(fcmp)`), so for
    `pp`/`pp_full` (`cmp_fit_file` empty → `fcmp=nullptr`) the new code is skipped;
  - the `pp_full` block is **byte-identical** to my committed `a4462db`; the file compiles.
  **Insulation:** my pp_full trig-eff code is committed clean (`a4462db`). The only residual risk
  is compile-time coupling — if that session leaves the file non-compiling at the instant the
  production run's Stage 10 compiles it, the plot step fails. Recovery is cheap: `git checkout
  a4462db -- plot_mc_trig_eff.cxx` (or its committed blob) and re-run only Stage 10d (the
  mc_trig NTP + hists + fits from 10a–c persist). Do NOT commit the r17663 WIP — the parallel
  session owns it.

- 2026-07-20 — **BUG in the production run (task 2 blocker, FIXED `31c2724`): `store_mc_trigger`
  vs the multi-file farm.** Stage 3 (the two mc_trig NTP passes) THREW
  `store_mc_trigger: expected one file per fullsim chain` — the trigger-branch presence check in
  `PythiaFullSimExtras.c` inspected one file and guarded against multi-file chains, but the LGD
  farm gives each pT-hat slice a 12–25-part TChain. **ROOT exits 0 on the swallowed C++ exception,
  so the pipeline `|| fail` did NOT catch it** and the run continued WITHOUT the mc_trig NTP
  (Stage 3 finished in 77 s — impossible for two full passes; that was the tell).
  Fix: per-file trigger-branch check requiring uniformity (all parts of a full-sample slice come
  from the same trigger-enabled grid task ⇒ uniform); all-with → bind, all-without → drop
  (trigger-off skim), MIXED → throw. Recompiles clean.
  **Consequence:** task 1 (non-trig-eff, Stages 1–9) is UNAFFECTED and completed. Task 2 needs a
  re-run of the 2 mc_trig NTP passes + Stage 10. **RE-RUN PLAN (run the mc_trig NTP DIRECTLY, not
  the whole pipeline, to avoid redoing the nominal+single NTP already done):**
  ```
  cd NTupleProcessingCode
  bash run_pythia_fullsim_mc_trig_full_sample.sh              # ~2-3h streaming
  bash run_pythia_fullsim_single_muon_mc_trig_full_sample.sh  # ~2-3h streaming
  # then Stage 10 for pp_full (Fill step1 -> Fit -> Fill step3 -> plot_mc_trig_eff)
  ```
  **Lesson for the pipeline:** a ROOT stage that throws but exits 0 slips past `|| fail`. The NTP
  scripts should propagate ROOT exceptions as a non-zero exit (e.g. wrap the `.q` / check output
  existence). TODO: harden the NTP run scripts + the pipeline's post-NTP validation to detect a
  missing/empty mc_trig NTP.

- 2026-07-21 — **TASK 1 (non-trig-eff) COMPLETE + VALIDATED on the FULL sample (Tight WP).**
  47 plots in `pythia_fullsim_full_sample/plots/`. Physics sanity (visual, C1–C3):
  - **crossx / statistics plot** (`reco_pair_pt_kn.png`, `truth_pair_pt_kn.png`): honest FULL-sample
    dσ/dp_T in **nb/GeV**, title says "FULL sample", the "NOT a physical pp σ" caption correctly
    OMITTED (pp-only, weight 1). Each pT-hat slice peaks in its own range, smoothly falling ~6
    decades; error bars far tighter than the test sample. Reco panel uses **Tight** (the WP fix).
  - **single-muon reco eff** (`single_muon_reco_effcy_vs_pt.png`): clean turn-on to a 0.90–0.95
    plateau by ~20 GeV; central |q·η|<0.5 bin correctly lower (~0.86, barrel crack); tight errors
    (plateau 0.9245±0.0006 at pT[20,30]). Titled "(FULL)", "tight WP".
  - **detector response** (`pair_pt_response_matrix_tightWP.png`): strong diagonal, modest
    resolution smearing, NO anomalous band (pp; the overlay 2·m_μ band is a closed, overlay-only
    issue). 
  ⇒ full-statistics improvement evident throughout. Task 1 physics is sound.
- 2026-07-21 — **TASK 2 (trig-eff) re-run in progress** (`run_pp_fullsim_trigeff_fullsample.sh`,
  detached, survived a session restart). Step [1] mc_trig pair NTP DONE + validated (1 652 223
  entries, 4.5 GB — the store_mc_trigger fix works on the full farm). Step [2] mc_trig single-muon
  NTP streaming; then Stage 10 (Fill→Fit→Fill-step3→`plot_mc_trig_eff` for pp_full).
  **Remaining after this:** the Medium-WP pass (reco-eff + trig-eff systematic variant; NTP/RDF are
  WP-agnostic so it skips them), and `/review-plot` sign-off.

## Results & Observations

### R1. Disk census of `~/usatlasdata` (real bytes, `du -sb`)

| dir | size | status |
|---|---|---|
| `dimuon_data/` | **1067 GB** | ACTIVE — pp_2024 511 G, pbpb_2025 271 G, pbpb_2023 140 G, pbpb_2024 111 G |
| `pythia_fullsim_hijing_overlay_test_sample/` | **277 GB** | 6 canonical r17618 NTUPs (22.1 G each = 133 G, ACTIVE) + 6 `.bak_20260709` trigger-off NTUPs (133 G, **OBSOLETE**) + r17662 |
| `pythia_private_sample/` | **127 GB** | `0317/0318/0325/0401/0429_all_k` ≈ 113 G |
| `powheg_full_sample/` | **61 GB** | fullsim parts (`mixed` 14 G, `mixed_mass_1_3GeV` 13 G, 2 × `_MYSTREAM` 14 G) — POWHEG fullsim is **obsolete** (memory `project_mc_sample_roles`); truth (bb/cc evgen 17 G) must stay |
| `pythia_truth_full_sample/` | **33 GB** | ACTIVE |
| `pythia_fullsim_test_sample/` | **9.5 GB** | ACTIVE (24 NTUPs + 24 `.bak`) |
| others | <1 GB | |

**The single biggest, safest win: the 6 overlay `.bak_20260709.root` files = 133 GB.** They are
the trigger-OFF NTUPs superseded by the July2026 trigger-on ones, which were validated as
identical outside the trigger branches (`mc_trigger_info_skim.md` step 8a) and carry the same
10 000 entries each.

**⚠ Forward-looking hazard:** the r17618 overlay NTUPs are **22 GB per 10 000 events** because
r17618 keeps the full HIJING truth record (r17662 signal-only-truth is 333 MB, ×65 smaller). If
the HIJING-overlay **full** sample has comparable statistics to this pp full sample (~10 M
events), its skim would be **~20 TB** — an order of magnitude beyond the entire quota. This must
be resolved (truth slimming in the skim, or r17662-style truth) before that sample lands.

## Remaining Work

**RESUME HERE (2026-07-14). Two things run in tmux with NO Claude session needed:**

```bash
tmux new -s lgd   ; bash ~/usatlasdata/lgd_migration/lgd_migration_resume.sh
tmux new -s farm  ; ~/workarea/dimuon_codes/SkimCode/scripts/fullsim_pp24_full_to_lgd.sh
```

1. **Finish the LGD migration** (~160 GB still to free). All 292 GB is ON LOCALGROUPDISK (5 rules
   `State: OK`); only group A was deleted (132.8 GB freed). The subagent stopped deleting because
   the LGD **POSIX door wedged** after the 14.7k-file write burst (ROOT readers in D-state;
   `stat` and **xrootd** both fine ⇒ transient, NOT data loss). All B/B2/C/D sources are intact
   on disk. Durable state (rule IDs, the irreplaceable `local|pnfs` maps): `~/usatlasdata/lgd_migration/`.
   The resume script probes the door first and **deletes nothing** while it is sick.
2. **`fullsim_pp24_full_to_lgd.sh`** — grid task done → `rucio add-rule` to LGD → symlink farm
   (`...NTUP.partNN.root`) → verify → optional dev slice. Idempotent; `--status` to inspect.
3. **Then:** `run_pythia_fullsim_full_sample.sh` (needs the farm).

**Open items (NOT done):**
- [x] ~~Overlay NTP rerun + overlay crossx replot~~ **DONE 2026-07-14** (`92d9f06`). NTP needed no
  rerun — verified already post-fix; only the plot was stale. Two honesty bugs fixed along the way
  (1000× ub/nb axis mislabel; pp plot now says it is NOT a physical pp σ).
- [x] ~~AMI weights + the hard rule~~ **DONE 2026-07-14** (`be9b77a`) — `Analysis/docs/ami_weights.md`.
- **`docs/pythia_fullsim_pp.md`: add a FULL-sample section** (farm, the `isTestSample` switch, the
  AMI dir + DSIDs, the `_full` output suffix).
- **Downstream `is_test_sample` defaults:** the NTuple-processing `isTestSample` defaults to
  **false** (full sample) as requested, but the *consumers* (`RDFBasedHistFillingPythiaFullsim`,
  `PythiaFullsimRecoEffPlotter`, `plot_single_muon_reco_effcy`,
  `RDFBasedHistFillingPythiaFullsimOverlay`) default to **`is_test_sample = true`**, because only
  the TEST sample has NTP output today. **Flip them when the full sample lands.** Flagged to the
  user.
- **Review loop: iteration 2 returned FAIL** (2 CRITICAL + 4 WARNING, all now fixed + verified;
  the reruns above are the outstanding WARNING). A 3rd `/review-analysis-code` pass is owed.
- Truth path (`InitInputCentrProd`) intentionally keeps 4:6:6:9 — it is a genuinely 4-beam
  Pb-intent sample. Left as-is; noted so a future agent does not "fix" it.

## Latest Stage

**Steps 1–5 DONE.** Storage census + user-approved migration (all 292 GB on LGD, 132.8 GB freed so
far, resume script written). Run-dir reorg done. Local test run PASSED every content check. All 6
grid tasks submitted (**51419497/500/506/511/515/520**); `51419497` and `51419511` already
`done` with zero failures. Design pivoted (D1) to **LGD-only, no 280 GB download**. Code:
multi-part farm reader + `isTestSample` isospin switch + AMI provenance guard, all compiled and
verified; 3 commits (`216f750`, `1b927d2`, `6af24dd`).

**Next:** run the two tmux scripts above (migration finish + farm build), then the overlay rerun
and the crossx-plot honesty labelling.
