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

*(to be filled as work proceeds)*

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

Everything in the Implementation Plan.

## Latest Stage

**Step 1 (storage census) — census DONE, candidate list being prepared for the user.**
Next, in parallel: Step 2 (run-dir reorganisation) and Step 3 (local test run) — neither needs
much disk if the test AOD is staged outside the `usatlast3-data` fileset.
