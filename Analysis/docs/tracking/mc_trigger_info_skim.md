# Trigger information in Run-3 fullsim MC skims (pp24 fullsim + HIJING overlay)

**Mode:** Implementation. **Created:** 2026-07-09. **Session:** "fullsim MC trigger info".
**Orchestrated** (subagents for discovery; orchestrator owns git + this doc).
**Reviewer rules:** C++/Athena/python skim-config changes → `/review-analysis-code`.
**Siblings:** `hijing_overlay_r17618_grid_reprocessing.md` (overlay grid workflow),
`grid_monitor_data_infra.md` (grid_monitor.sh), `mu4_trig_effcy_implementation.md`
(data trigger-efficiency method), `analysis_roadmap_2026_06.md` Q4 (dR trigger correction
blocked on "unbiased trigger decision in MC").

---

## Objective

Determine whether the Run-3 Pythia fullsim MC samples (pp24 conditions, and HIJING-overlay
PbPb conditions) contain **trigger simulation** output in their AODs; if so, turn the trigger
branches ON by default in the SkimCode configuration for **all Run-3 fullsim/overlay MC**
(everything except EvGen/truth-only), re-skim on the grid with a `_July2026` tag, and update
bookkeeping.

## Physics Procedure (AUTHORITATIVE)

### 1. Motivation

The analysis corrects data yields by a per-pair trigger efficiency ε_trig^pair. Today that
efficiency is measured **from data** (tag-and-probe / `FillHistogramsDimuTrigGivenMu4`), and
the **ΔR correlation correction** between the two muons' trigger decisions is a **dummy
(ε_dR ≡ 1)** because a data-derived cross-term is biased on a mu4-triggered sample
(roadmap D9 / Q4). The roadmap explicitly names the missing prerequisite:

> "dR trigger-correlation correction (final) — **Unbiased trigger decision in MC** → needs
> fullsim overlay with trigger sim."  (`analysis_roadmap_2026_06.md` Q4)

Trigger simulation in the fullsim/overlay MC therefore unblocks:
- (a) an **unbiased** MC-based trigger efficiency, including the two-muon ΔR correlation;
- (b) applying the **same trigger selection** to MC as to data (pp24 `HLT_2mu4_L12MU3V`;
  PbPb `HLT_mu4_L1MU3V`), so MC-based reco efficiency and detector response are evaluated
  on a trigger-selected sample, matching the data selection;
- (c) a full **2mu4 trigger-efficiency evaluation** in MC (single-leg turn-on × correlation).

### 2. What must be present in the AOD / NTUP

For the skim to reproduce the data trigger content, the AOD must carry the HLT/L1 trigger
decision and navigation, so that `TrigDecisionTool` and `Trig::R3MatchingTool` work:
- **Trigger decisions** for the analysis chains: `HLT_2mu4_L12MU3V` (pp analysis trigger),
  `HLT_mu4_L1MU3V` (PbPb analysis trigger), `HLT_mu4_mu4noL1_L1MU3V` (support chain).
- **Per-muon trigger matching** to those chains (needed for tag-and-probe and for the
  per-leg turn-on), i.e. HLT muon feature navigation.
- Optionally L1 items / L1 TE (data-only in the current config; not required for MC).

The convener states the reco tag **r16578** AMI configuration has `doRDO_TRIG` and
`doTRIGtoALL` in `steering`, i.e. trigger simulation ran before reconstruction and trigger
content was written to all output formats. This must be **verified two ways**:
1. **AMI tag configuration** (`pyami` / `AMIGetDatasetInfo` / `ami show`) for the reco tag.
2. **Explicitly in the produced sample** — open the AOD and confirm the trigger containers
   (`xTrigDecision`, `TrigConfKeys`, HLT navigation summary, HLT muon containers) exist and
   that the chains above are in the menu.

Both checks are required: an AMI steering flag says what was *requested*; the file says what
is *there*.

### 3. Step-by-step method

- **(a) AMI check** for pp24 fullsim reco tag `r16578` (example DSID
  `mc23_5p36TeV.803018.Py8EG_A14_pp_hQCD_DiMu_pTH40_70_pdf.recon.AOD.e8599_e8586_s4521_s4483_r16578`)
  and for the overlay reco tags `r17618` and `r17662` (example DSID
  `mc23_5p36TeV.802781.Py8EG_A14_pp_hQCD_DiMu_pTH8_14.merge.AOD.e8599_s4614_r17618_r15970`).
  For the overlay the trigger content could equally come from the **digitization tag r15970**
  — check that tag too. If trigger info is in **neither** r17618/r17662 nor r15970 → the
  overlay has no trigger sim (report; we NEED it).
- **(b) File-level check.** Open one AOD per configuration (via rucio / xrootd) and list the
  trigger-related containers + the chain names in the menu. Confirm `HLT_2mu4_L12MU3V`,
  `HLT_mu4_mu4noL1_L1MU3V`, `HLT_mu4_L1MU3V` are present.
- **(c) Enable trigger for MC in SkimCode.** In `scripts/TrigRates_CA.py` and
  `run_pythia_fullsim_HIJING_overlay/TrigRates_CA.py`, the trigger is currently disabled
  **solely because the sample is MC**:
  `alg.UseTrigger = (not is_MC)`, `alg.StoreL1Decision = (not is_MC)`,
  `alg.StoreL1TE = (not is_MC)`, `TrigDecisionTool`/`R3MatchingTool` built only `if not is_MC`,
  `alg.HLTMuonsKey = ""` for MC.
  The change is to gate trigger on a **new, explicit** `has_trigger_sim` flag (true for the
  Run-3 fullsim + overlay MC modes, false for EvGen/truth-only and for legacy Run-2 MC),
  not on `is_MC`, and to populate `Muon_triggers` / `DiMuon_triggers` for the MC modes with
  the SAME chain lists as the corresponding data mode (pp24 fullsim ← `pp2024` list;
  HIJING overlay ← `hi2023` list, since the overlay uses PbPb-2023-like conditions).
- **(d) Test job.** Run a local `athena.py TrigRates_CA.py --evtMax=N` on one AOD per
  configuration.
  - If it crashes with an **obvious non-physics bug** (missing container key, wrong tool
    config, GRL applied to MC, etc.) → fix and re-test.
  - If the error is **physics-affecting or ambiguous** → STOP and escalate to the user.
- **(e) Re-skim on the grid** with `_July2026` suffix for all DSIDs of each sample family
  (pp24 fullsim: 4 beam combos × 6 pTHat slices = 24; overlay: 6 pTHat slices), monitor with
  `SkimCode/scripts/grid_monitor.sh`, download + hadd + validate entry counts, update
  bookkeeping (`merging-record.txt`, `SkimCode/README.md` run-mode table).

### 4. Negative constraints

- Do **NOT** enable trigger for EvGen/truth-only skims (`TruthTrigRates.py`, `run_pythia_truth`)
  or for legacy Run-2 samples (`do_pp_MC_fullsim_17`, `TrigRates_JO.py`) — no trigger sim there.
- Do **NOT** apply the GRL to MC (`UseGRL` must stay consistent with current behavior for MC).
- Do **NOT** change any muon selection, calibration, WP, or kinematic setting while enabling
  trigger. This task adds trigger branches only. The current default muon WP is **Tight**
  (`tight_wp_default_change.md`) — the skim packs the full quality bitmask regardless.
- Do **NOT** clobber the existing (no-trigger) NTUPs — new grid outputs carry `_July2026` and
  the local NTUPs must be written to a distinct name/backup dir.
- The skim writes trigger *decisions and matching*; it does **not** apply a trigger cut.
  Any trigger requirement is a downstream (NTupleProcessing/RDF) decision.
- If the overlay AODs turn out to have **no** trigger sim in any of r17618 / r17662 / r15970,
  do **not** fabricate a workaround — report to the user (we need it).

## Context

- SkimCode run modes and trigger status (`SkimCode/README.md`): `ppmcfullsim2024` and
  `ppmcfullsim_hioverlay24` both have **Triggers: off**.
- Trigger gating lives in `scripts/TrigRates_CA.py` (canonical, copied into `run_23hi`,
  `run_24hi`, `run_24pp`, `run_25hi`) and in the specialized
  `run_pythia_fullsim_HIJING_overlay/TrigRates_CA.py`.
- `README.md` "Non-obvious details": *"`HLT_MuonsCB_RoI` / `HLT_MuonsCB_FS` keys are data-only.
  For MC modes the HLT muon container keys are left empty; setting them on a MC AOD causes
  retrieval failures."* — this note was written when MC had no trigger sim; it must be
  re-evaluated (and the README updated) if the fullsim AODs do have HLT muon containers.
- pp24 fullsim NTUPs: 24 files `Pythia_5p36TeV_{pp,pn,np,nn}_hQCD_DiMu_pTH{8_14,14_24,24_40,
  40_70,70_125,125_300}.FullSimPP24.NTUP.root` in `~/usatlasdata/pythia_fullsim_test_sample/`
  (10k events each). **No grid_sub.sh exists in SkimCode for this family** — it must be
  reconstructed (DSID list to be obtained from AMI).
- Overlay NTUPs: 6 files `Pythia_5p36TeV_pp_hQCD_DiMu_pTH*.FullSimHIJINGOverlayPP24.NTUP.root`
  in `~/usatlasdata/pythia_fullsim_hijing_overlay_test_sample/` (10k each);
  `run_pythia_fullsim_HIJING_overlay/grid_sub.sh` (r17618) and
  `grid_sub_r17662_signalonly.sh` (r17662, pTH8_14 only).
- Overlay r-tag comparison: `hijing_overlay_truth_barcode_duplicate_investigation.md`
  (r17618 = full HIJING truth, baseline; r17662 = StandardSignalOnlyTruth).

## Scope

**In:** AMI + file-level trigger verification; SkimCode config change to enable trigger for
Run-3 fullsim/overlay MC; local test jobs; grid re-skim (`_July2026`); monitoring, download,
hadd, validation; bookkeeping (merging records, SkimCode README, this doc, roadmap Q4 note).

**Out:** downstream NTupleProcessing/RDF consumption of the new trigger branches; the MC-based
ΔR trigger-correlation correction itself (a follow-up, unblocked by this task).

## Design Decisions

*(to be filled as work proceeds)*

## Implementation Plan

1. [x] **Discovery** — SkimCode trigger code map (TrigRates.cxx branches + gating);
   pp24-fullsim DSID list; grid_monitor usage for MC modes.
2. [x] **AMI check** — reco tags r16578 (pp fullsim), r17618 / r17662 (+ digi r15970, overlay):
   `steering` field, `doRDO_TRIG` / `doTRIGtoALL`. (§3a) → all three have both.
3. [x] **File check** — open one AOD per config; confirm trigger containers + the 3 chains. (§3b)
   → all present, unprescaled.
4. [x] **Enable trigger for MC** in both `TrigRates_CA.py` copies (new `mc_has_trigger_sim` gate,
   MC chain lists). → `/review-analysis-code` **PASS** iter 1. Commit `434de17`. (§3c)
   Includes the `StoreAllEvents` fix (R1b).
5. [x] **Test job** on one AOD per config → all rc=0, trigger branches present, rates sane. (§3d)
6. [~] **Grid re-skim** `_July2026` — 31 tasks submitted 2026-07-09, monitors running.
   Awaiting completion → download/hadd/validate. (§3e)
7. [ ] **Bookkeeping** — merging records (auto by grid_monitor), SkimCode README run-mode table
   (done, `434de17`), roadmap Q4 note, this doc.
8. [ ] **Post-skim sanity check** — each new NTUP must have **10 000 entries** (same as the
   trigger-off NTUPs) and non-empty `b_HLT_*` / `dimuon_b_HLT_2mu4_*` branches. A count < 10 000
   would mean the `StoreAllEvents` fix did not reach the grid job.

## Progress Log

*(append-only)*

- 2026-07-09 — **Step 1 (Discovery) DONE.** Read-only Explore agent mapped the trigger flow.
  Key facts (file:line):
  - Trigger is gated ONLY on `is_MC` in both CA configs: `alg.UseTrigger = (not is_MC)`
    (`scripts/TrigRates_CA.py:255`), `StoreL1Decision` (`:256`), `StoreL1TE` (`:269`);
    `TrigDecisionTool` built only `if not is_MC` (`:229`), `R3MatchingTool` (`:242`);
    `HLTMuonsKey=""` for MC (`:283-286`; overlay copy `:278-286` also clears `HLTMuonsFSKey`).
  - `TrigRates.cxx`: `AddTriggerBranches()` (458-629) + `ProcessTriggers()` (631-689) run only
    `if(m_use_trigger && !m_is_evgen)` (284-288, 301-303). Branches written:
    `b_<chain>` (decision, always); `b_<chain>_L1TBP/_L1TAP/_L1TAV` only if `StoreL1Decision`;
    prescale/isPrescaled/isRerun branches are dead (local `save=false`, line 488).
    Muon matching: `muon_b_<chain>`, `_0_01`, `_V2`, `_V3` (541-566, needs `StoreSingleMuon`).
    Dimuon matching: `dimuon_b_<chain>_<DR>`, `_0_01`, and 8 per-leg branches
    `dimuon_b_<chain>_mu{1,2}pass Leg{1,2}_dR_*` (568-617, needs `StoreAcoplanarMuon` +
    `StoreDimuonPerLeg`). Also `muon_trig_match`, `muon_match_mu4roi`, `muon_match_mu6roi`
    (918-983, guarded `if(m_use_trigger)`).
  - **Hazard:** `ProcessMuons` (1156-1162) does `CHECK(evtStore()->retrieve(l_muons_trig_roi,
    m_hlt_muons_key))` when `m_use_trigger` — hard failure if `HLTMuonsKey` is set but the
    container is absent (this is exactly why MC had it blank). So enabling trigger REQUIRES
    `HLT_MuonsCB_RoI` / `HLT_MuonsCB_FS` to exist in the MC AOD.
  - Chain lists are TDT chain-group expansions of the `"|"`-joined property strings
    (`TrigRates.cxx:459-461`) → every requested chain must exist in the sample's HLT menu.
  - GRL already skipped for MC at runtime (`CleaningCuts`, `:385` `if(m_use_GRL && !isMC)`),
    so `alg.UseGRL=True` is harmless for MC.
  - `grid_monitor.sh` supports only `--mode data` and `--mode overlay` (38-51). **No mode for
    the non-overlay `FullSimPP24` family.**
  - **No pathena submission script for `FullSimPP24` exists** anywhere in the repo or in git
    history — it must be written from scratch. (Overlay has `grid_sub.sh` + `grid_sub_r17662_signalonly.sh`.)
  - Stale copies: `run_24pp/TrigRates_CA.py`, `run_24{hi,pp}_with_dimuon_trig_match/TrigRates_CA.py`
    differ from canonical; `run_23hi`, `run_24hi`, `run_25hi` are byte-identical to canonical.

- 2026-07-09 — **Step 2 (AMI) DONE — all three reco tags run trigger simulation.**
  `ami cmd AMIGetAMITagInfo -amiTag=<t>`:
  | tag | step | release | steering | trigger menu (preExec) | note |
  |---|---|---|---|---|---|
  | `r16578` | recon | Athena_24.0.95 | `"doRDO_TRIG" "doTRIGtoALL"` | (none in preExec; `preInclude="all:Campaigns.MC23ppReferenceRun2024"`) | "HI pp tag with pileup, clone of r16554" |
  | `r17618` | recon | Athena_24.0.58 | `"doRDO_TRIG" "doTRIGtoALL"` | `flags.Trigger.AODEDMSet='AODFULL'; flags.Trigger.triggerMenuSetup='Dev_HI_run3_v1_TriggerValidation_prescale'` | overlay baseline (ATLHI-576 vtx fix) |
  | `r17662` | recon | Athena_24.0.58 | `"doRDO_TRIG" "doTRIGtoALL"` | same as r17618 | clone of r17618 + `StandardSignalOnlyTruth` |
  | `r15970` | **merge** | Athena_24.0.52 | — (`AODMerge_tf.py`) | — | merge step only; trigger comes from the recon tag, not here |
  Also `r17663` = "clone r17618 w/o Hijing overlay", same steering (this is the tag of the AOD
  that happened to be cached locally at `pythia_fullsim_hijing_overlay_test_sample/mc23_5p36TeV/`).
  ⇒ Convener's statement CONFIRMED for r16578 **and** independently for r17618/r17662.

- 2026-07-09 — **Step 3 (file-level check) DONE — trigger content present, all 3 chains, PS=1.**
  Opened one AOD per configuration and dumped `CollectionTree` branches + the
  `MetaData` `TriggerMenuJson_{HLT,L1,HLTPS}` payloads.
  - **pp fullsim r16578** (`mc23_5p36TeV.802781...recon.AOD.e8599_s4521_s4483_r16578`,
    `AOD.49223786._000003.pool.root.1`, 1000 evt, 3206 branches):
    present → `xTrigDecision`, `TrigConfKeys`, `HLTNav_Summary_AODSlimmed`,
    **`HLT_MuonsCB_RoI`**, **`HLT_MuonsCB_FS`**, `HLT_Muons_RoI`, `LVL1MuonRoIs`.
    HLT menu = **`PhysicsP1_pp_lowMu_run3_v1`** (610 chains); L1 menu `Physics_HI_run3_v1`.
    Contains ALL pp2024 data chains: `HLT_mu4_L1MU3V`, `HLT_mu6_L1MU3V`, `HLT_mu8_L1MU5VF`,
    `HLT_mu10_L1MU8F`, `HLT_mu12_L1MU8F`, `HLT_mu15_L1MU8F`, `HLT_mu15_L1MU14FCH`,
    `HLT_2mu3_L12MU3V`, `HLT_2mu4_L12MU3V`, `HLT_mu4_mu6_L12MU3V`, `HLT_mu4_mu4noL1_L1MU3V`.
    Prescale set `PhysicsP1_pp_lowMu_run3_v1_TriggerValidation_prescale`:
    `mu4`, `2mu4`, `mu4_mu4noL1` all **prescale=1, enabled=true**.
  - **overlay (HI menu)** — checked on the locally cached r17663 AOD, which shares r17618's
    trigger config (clone, identical `preExec`/steering): `xTrigDecision`, `TrigConfKeys`,
    `HLTNav_Summary_AODSlimmed`, `HLT_MuonsCB_RoI`, `HLT_MuonsCB_FS`, `LVL1MuonRoIs` present.
    HLT menu = **`Dev_HI_run3_v1`** (690 chains); L1 menu `MC_HI_run3_v1`.
    Contains ALL hi2023 data chains: `HLT_mu4_L1MU3V`, `HLT_mu6_L1MU3V`, `HLT_mu6_L1MU5VF`,
    `HLT_mu4_L1MU3V_VTE50`, `HLT_mu6_L1MU3V_VTE50`, `HLT_mu8_L1MU5VF_VTE50`,
    `HLT_2mu4_L12MU3V`, `HLT_mu4_mu4noL1_L1MU3V`. Prescale set
    `Dev_HI_run3_v1_TriggerValidation_prescale`: `mu4`, `2mu4`, `mu4_mu4noL1` all
    **prescale=1, enabled=true**. (r17618/r17662 AODs downloading for direct confirmation.)
  ⇒ **Everything needed for skimming (2mu4, mu4_mu4noL1, mu4) and a full 2mu4 trigger-efficiency
  evaluation (per-leg matching via `HLT_MuonsCB_RoI`/`HLT_MuonsCB_FS`) is present.**

- 2026-07-09 — **Step 4 (enable trigger for MC) CODE WRITTEN.** Both CA configs changed
  identically. New explicit gate replaces the `is_MC` proxy:
  ```python
  mc_has_trigger_sim = do_pp_MC_fullsim_24 [or do_pp_MC_fullsim_hioverlay24]
  use_trigger = (not is_MC) or mc_has_trigger_sim
  ```
  Applied to: `TrigDecisionTool` + `R3MatchingTool` construction, `alg.TriggerMatchTool`,
  `alg.UseTrigger`, `alg.StoreL1Decision`, `alg.StoreL1TE`, and the
  `HLTMuonsKey`/`HLTMuonsFSKey` block (now `HLT_MuonsCB_RoI`/`HLT_MuonsCB_FS` whenever
  `use_trigger`, verified present in all three MC AODs).
  MC chain lists added: `ppmcfullsim2024` ← the pp2024 data lists (7 single-μ + 4 dimuon);
  `ppmcfullsim_hioverlay24` ← the hi2023 data lists (6 single-μ + 2 dimuon). Every chain was
  verified present in the corresponding simulated HLT menu (Step 3).
  Legacy Run-2 fullsim (`do_pp_MC_fullsim_17`) and truth-only skims keep trigger OFF.
  Also added `TRIGRATES_OUTPUT_FILE` support to the canonical script (default `myfile.root`,
  so data run dirs are unchanged) — needed for per-sample grid output names.
  Synced canonical → `run_23hi`, `run_24hi`, `run_25hi`; created **new run dir
  `run_pythia_fullsim/`** (none existed for the pp24-fullsim family).
  Both files `ast.parse` clean; no `(not is_MC)` trigger gates remain.

  **Note (repo hygiene, pre-existing):** `.git/info/exclude` has `SkimCode/run*/**`, so
  `run_pythia_fullsim_HIJING_overlay/TrigRates_CA.py` — the overlay's *specialized* config —
  is **not tracked by git**. Only `SkimCode/scripts/TrigRates_CA.py` is.

- 2026-07-09 — **Step 5 (test jobs) DONE — all three configurations succeed.**
  Invocation (the working recipe, from `run_test.sh`): source `build_25/.asetup.save` +
  `build_25/x86_64-el9-gcc14-opt/setup.sh`, then
  `athena.py --no-excabort TrigRates_CA.py --filesInput=<AOD> --evtMax=N`.
  (Sourcing `setup_25.sh` from inside a run dir puts athena.py in *legacy* mode, where
  `flags.fillFromArgs()` chokes on the script name in `sys.argv` — an environment artifact,
  not a code bug. `athena.py:101` only strips the script path in CA mode.)

  | config | mode | AOD | evts read | entries | trig branches | rc |
  |---|---|---|---|---|---|---|
  | pp fullsim | `ppmcfullsim2024` | 802781 r16578 `AOD.49223786._000003` | 200 | 200 | 112 | 0 |
  | overlay r17618 | `ppmcfullsim_hioverlay24` | `AOD.50035726._000013` | 100 | 100 | 76 | 0 |
  | overlay r17662 | `ppmcfullsim_hioverlay24` | `AOD.50427580._000029` | 100 | 100 | 76 | 0 |

  Trigger rates (unbiased, `StoreAllEvents=True`): pp fullsim mu4 166/200, mu4_mu4noL1 114/200,
  2mu4 68/200. Overlay r17618: 75 / 41 / 27. Overlay r17662: 71 / 43 / 20. Ordering
  `mu4 > mu4_mu4noL1 > 2mu4` holds everywhere (a mu4_mu4noL1 event must fire mu4; 2mu4 is a
  subset). Per-leg 2mu4 matching branches (8) present in all — this is what the ΔR
  trigger-correlation correction needs.

- 2026-07-09 — **Step 4/5 REVIEWED: `/review-analysis-code` PASS, iter 1, 0 CRITICAL / 0 WARNING**
  (log `.claude/logs/review-analysis-code-20260709-223736-mc-trigger-info-skim.md`).
  Reviewer independently reproduced every number and confirmed data-mode behavior is
  bit-identical. 5 INFO notes: (1) `HLTMuonsFSKey` is a **dead property** — never read;
  `Module_TrigMuonMatching.cxx:256` hardcodes the `HLT_MuonsCB_FS` string. (2) `StoreL1TE`
  hard-retrieves `LVL1EnergySumRoI` independently of `m_use_trigger`; container verified
  present in all three MC AODs. (3) `grid_monitor` auto-backs-up the old NTUP to
  `.bak_YYYYMMDD.root` before writing, so the canonical filename is reused safely (and the
  canonical name is *required* by `FullSimSampleType.h`). (4) overlay `RunYear=2024` +
  hi2023 chain list is self-consistent (`RunYear` only drives the `Module_EventShape`
  centrality calibration). (5) overlay CA config is not git-tracked.
  Committed `434de17` (+ `SkimCode/README.md` run-mode table / HLT-key note corrected).

- 2026-07-09 — **Step 6 (grid re-skim `_July2026`) SUBMITTED — 31 tasks, 0 submission errors.**
  User decisions (2026-07-09): (1) MC keeps all events, decision recorded only
  (`StoreAllEvents` fix confirmed); (2) re-skim the **802758–802781** family (the 24 DSIDs
  behind the current NTUPs), **not** the `_pdf` 8030xx production.
  Version tag `July2026.v1` (matches the April2026.v1 / June2026.v1 convention; the local
  hadded NTUP filename stays canonical because `FullSimSampleType.h` requires it —
  `grid_monitor` renames the old file to `.bak_YYYYMMDD.root`).

  **pp fullsim** (`run_pythia_fullsim/grid_sub.sh`, new file), mode `ppmcfullsim2024`:

  | beam | pTH8_14 | pTH14_24 | pTH24_40 | pTH40_70 | pTH70_125 | pTH125_300 |
  |---|---|---|---|---|---|---|
  | nn | 51360068 | 51360048 | 51360056 | 51360060 | 51360064 | 51360052 |
  | np | 51360093 | 51360077 | 51360081 | 51360085 | 51360089 | 51360073 |
  | pn | 51360116 | 51360100 | 51360104 | 51360109 | 51360113 | 51360096 |
  | pp | 51360141 | 51360124 | 51360128 | 51360133 | 51360137 | 51360120 |

  **HIJING overlay** (`run_pythia_fullsim_HIJING_overlay/grid_sub.sh`, r17618), mode
  `ppmcfullsim_hioverlay24`: pTH125_300 51360153, pTH14_24 51360157, pTH24_40 51360162,
  pTH40_70 51360166, pTH70_125 51360170, pTH8_14 51360174.
  **Overlay r17662 signal-only** (`grid_sub_r17662_signalonly.sh`): pTH8_14 **51360178**.

  Monitoring: `grid_monitor.sh --mode fullsim_pp -i 10 <24 ids>` and
  `grid_monitor.sh --mode overlay -i 10 <7 ids>` (both running; download → hadd → validate →
  append to each family's `merging-record.txt`).

- 2026-07-09 19:25 — **First task landed & Step-8 sanity check PASSES on it.**
  Task `51360178` (overlay r17662, pTH8_14) validated OK — **10 000 entries**, i.e. *identical*
  to the trigger-off NTUP ⇒ the `StoreAllEvents` fix reached the grid job (a trigger-OR-filtered
  output would have been ~7–8k). `grid_monitor` renamed the old file to
  `...FullSimHIJINGOverlayPP24_r17662.NTUP.bak_20260709.root` (June-11 file preserved).

  | | entries | branches | trigger branches |
  |---|---|---|---|
  | old (trigger-off, Jun 11) | 10 000 | 143 | 0 |
  | new (trigger-on, Jul 9) | 10 000 | 224 | 76 |

  Full-stat trigger rates: `mu4` 7838/10000 (78.4%), `mu4_mu4noL1` 5258 (52.6%),
  `2mu4` 2959 (29.6%) — ordering `mu4 > mu4_mu4noL1 > 2mu4` holds and the values track the
  100-event test (75/41/27%). Per-leg matching is populated:
  `dimuon_b_HLT_2mu4_L12MU3V_mu1passLeg1_dR_0_02` has ≥1 matched muon in 2557 events
  (≤ the 2959 that fire 2mu4, as it must be). **The ΔR trigger-correlation ingredient is real.**

- 2026-07-09 — **STORAGE INCIDENT during Step 6 download: GPFS quota exhausted → 2 tasks
  falsely "failed".** `usatlast3-data` hit its soft limit with the grace period **expired**
  (`used 3223186944 ≥ soft 3221225472` blocks), so GPFS refused new writes. `rucio download`
  then failed for tasks `51360153` and `51360174` — every protocol/RSE failed in <1 s, the
  signature of a local write refusal, not a transfer problem. The grid tasks themselves were
  fine.
  - **Danger discovered:** `grid_monitor` renames the old NTUP → `.bak_YYYYMMDD.root` **before**
    downloading. For those two slices the rename succeeded and the download did not, leaving
    `pTH8_14` and `pTH125_300` with **no canonical NTUP** and the `.bak` as the *only* copy.
    ⇒ Never move/delete a `.bak_*` until its canonical replacement has landed.
  - **Response:** stopped both monitors (so healthy tasks stopped being marked failed), reset
    the two tasks `failed`→`pending`, freed space, restarted the monitors on the outstanding IDs.
  - **Space freed (user ran the deletes; agent `rm` on `usatlasdata` is permission-blocked):**
    (1) `aod_trigger_check/` — the AODs downloaded for the Step-3 file check (28 G blocks);
    (2) **`backup_wrong_vtx_z/` moved to pnfs** — the 6 obsolete r17044 wrong-beamspot-z NTUPs
    (`hijing_overlay_r17618_grid_reprocessing.md`), 132 GiB real. Copied to
    `/pnfs/usatlas.bnl.gov/users/yuhanguo/pythia_fullsim_hijing_overlay_test_sample/backup_wrong_vtx_z/`
    (= `~/dcachearea/...`), each file size-verified **and** re-opened in ROOT (all 6 → 10 000
    entries) before the source was removed; (3) the leftover `user.yuhang.NTUP.*_EXT0` partial
    download dir. Quota 99% → 90.3% real (1.35/1.50 TiB).
  - **Quota accounting note (confirms memory `reference_storage_quota`):** `mmlsquota` and
    `du -sh` report GPFS **allocated blocks = 2× real bytes** (data replication). `du -sb`
    reports real bytes. So the mmlsquota numbers MUST be halved: usable = 1.50 TiB, not 3.0 TiB.
  - **Disk census (real bytes):** `dimuon_data` 995 G (pp_2024 476 G, pbpb_2025 253 G,
    pbpb_2023 131 G, pbpb_2024 104 G) · overlay test sample 311 G · `pythia_private_sample`
    119 G · `powheg_full_sample` 57 G · `pythia_truth_full_sample` 31 G · `pythia_fullsim_test_sample`
    7.4 G. **The overlay NTUPs are 21 GiB per 10 000 events** because r17618 keeps the FULL
    HIJING truth record; the r17662 signal-only-truth equivalent is 333 MB (×65 smaller).
    Relevant if overlay statistics ever grow.

- 2026-07-09 — **✅ OVERLAY FAMILY COMPLETE & VALIDATED (7/7 tasks).** All six r17618 slices +
  r17662: **10 000 entries each** (= the trigger-off originals ⇒ no trigger bias) and **76
  trigger branches** each. Trigger rates rise **monotonically with pT-hat**, and the ordering
  `mu4 > mu4_mu4noL1 > 2mu4` holds in every slice:

  | slice | mu4 | mu4_mu4noL1 | 2mu4 |
  |---|---|---|---|
  | pTH8_14 | 7833 | 5234 | 2943 |
  | pTH14_24 | 8559 | 6238 | 3838 |
  | pTH24_40 | 8895 | 6838 | 4639 |
  | pTH40_70 | 9202 | 7404 | 5144 |
  | pTH70_125 | 9346 | 7739 | 5629 |
  | pTH125_300 | 9371 | 7996 | 5618 |

  Harder scatters → stiffer muons → higher trigger rates. This is a physics-level validation,
  not just a plumbing check.

- 2026-07-09 — **Step 8 SANITY CHECKS on the downloaded pp-fullsim files (user-requested).**

  **(8a) Re-skim identical to backup outside trigger branches — PASS (22/22 file pairs).**
  ⚠ *Method trap:* a first, positional (event-index) comparison reported huge differences in
  round multiples of 2000. Root cause: each NTUP is `hadd`ed from 5 grid jobs × 2000 events and
  the **chunk order differs between skims**. The *sets* of `eventNumber` are identical
  (10 000 = 10 000, zero one-sided). Re-run keyed on `eventNumber`: **every one of
  `muon_pt`, `muon_eta`, `muon_quality`, `muon_trk_pt`, `muon_d0`, `muon_truth_pt`,
  `truth_muon_pt` is bit-identical for all 10 000 events in all 22 pairs.** ⇒ enabling the
  trigger added branches and changed nothing else. (Always key MC-vs-MC comparisons on
  `eventNumber`, never on entry index.)

  **(8b) Cross-section-weighted trigger fractions, per isospin beam.** Complete beams `nn`,
  `np`, `pn` (6/6 pTHat slices each; `pp` was 4/6 at the time). Weight per event
  `w_s = σ_s·ε_filt,s / N_s` (AMI, from the processing run log; `pp pTH40_70` cross-checked
  against AMI directly).
  `f(X|mu4) = Σ_s w_s·n_{X∧mu4,s} / Σ_s w_s·n_{mu4,s}`

  | beam | f(2mu4 \| mu4) | f(mu4_mu4noL1 \| mu4) |
  |---|---|---|
  | nn | 0.4545 | 0.7203 |
  | np | 0.4557 | 0.7233 |
  | pn | 0.4552 | 0.7276 |

  Beams agree to <1% ⇒ the trigger response is isospin-blind, as it must be. Both fractions
  rise monotonically with pTHat. **In MC, `2mu4 ⊆ mu4` and `mu4_mu4noL1 ⊆ mu4` EXACTLY** (zero
  events fire the dimuon chain without mu4) — a strong internal-consistency check.

  **(8c) Data comparison (pp24, 603 M events) — and the prescale structure it reveals.**
  In DATA the subset relation FAILS, and that is the prescale:
  - `N(2mu4)/N(2mu4∧mu4) = 11 292 618 / 1 212 963 = 9.310` — 2mu4 is unprescaled and seeds a
    *different* L1 item (`L1_2MU3V`), so this measures the **total** mu4 prescale.
  - `N(noL1)/N(noL1∧mu4) = 13 985 355 / 5 403 586 = 2.588` — `mu4_mu4noL1` shares the
    `L1_MU3V` seed, so the L1 draw is common and this measures the **HLT-only** mu4 prescale.
  - ⇒ implied `PS(L1_MU3V) = 9.310/2.588 = 3.597`.

  Independent bookkeeping (`IntNotes/data/luminosity/pp_2024/*.csv`, `LAr Corrected` ÷
  `Prescale Corrected`): `PS(mu4)=9.283`, `PS(2mu4)=1.000`, `PS(mu4_mu4noL1)=3.657`.
  **Agreement: 0.3% on the total mu4 prescale, 1.6% on the L1_MU3V prescale.** Two completely
  independent sources (event counts vs lumi bookkeeping) — and it proves the conditional ratios
  `f(X|mu4)` are **prescale-free** (the mu4 prescale cancels between numerator and denominator;
  for noL1 the shared L1 draw cancels too, and its HLT prescale is 1).

  **(8d) MC vs data — the comparison only works at fixed muon quality.** Inclusive fractions are
  meaningless across samples (MC is DiMu-filtered; data mu4 is dominated by single muons):
  MC 0.455 vs data 0.0041. Requiring ≥2 fiducial reco muons (`pT>4 GeV`, `|η|<2.4`) is NOT
  enough either (MC/data = 13.0). **The driver is muon quality.** Requiring both muons **Tight**
  (`quality & 305 == 305`) — the analysis nominal WP — collapses the discrepancy, differentially
  in sub-leading-muon pT:

  | sublead μ pT [GeV] | f(2mu4\|mu4) MC | data | MC/data | f(noL1\|mu4) MC | data | MC/data |
  |---|---|---|---|---|---|---|
  | 4–5   | 0.549 | 0.380 | 1.44 | 0.874 | 0.792 | 1.10 |
  | 5–6   | 0.696 | 0.530 | 1.31 | 0.953 | 0.887 | 1.07 |
  | 6–8   | 0.769 | 0.620 | 1.24 | 0.978 | 0.946 | 1.03 |
  | 8–10  | 0.832 | 0.661 | 1.26 | 0.994 | 0.965 | 1.03 |
  | 10–15 | 0.820 | 0.662 | 1.24 | 0.986 | 0.973 | 1.01 |
  | 15–25 | 0.702 | 0.701 | 1.00 | 0.942 | 0.979 | 0.96 |
  | 25–50 | 0.788 | 0.635 | 1.24 | 0.947 | 0.981 | 0.97 |
  | **integrated** | **0.638** | **0.459** | **1.39** | **0.918** | **0.842** | **1.09** |

  Only **250 517** of the 4 068 784 data "≥2 fiducial muon, mu4-fired" events survive the Tight
  requirement (6%) — i.e. **94% of the loose data dimuon sample has a non-tight second muon**
  that never produces a trigger object. MC (DiMu-filtered, two real muons) has almost none.
  This independently corroborates the Medium→Tight WP default (`tight_wp_default_change.md`).

  **Physics reading of the residual (this is the important part):** after the Tight requirement,
  **2mu4 (needs TWO L1 muon RoIs) has MC/data ≈ 1.24–1.44**, while **mu4_mu4noL1 (needs only ONE
  L1 RoI, the second leg is full-scan) has MC/data ≈ 0.96–1.10**. A per-leg L1 muon efficiency
  ratio `r = ε_L1^MC/ε_L1^data ≈ 1.13` reproduces BOTH: `r² ≈ 1.28` (2mu4) and `r ≈ 1.13` →
  ~1.0–1.1 (noL1). That is the known MC over-efficiency of the L1 muon trigger (simulation lacks
  real RPC/TGC chamber inefficiency), and the discrepancy is largest in the turn-on
  (4–5 GeV: 1.44) and vanishes on the plateau (15–25 GeV: 1.00) — exactly the expected shape.
  **Consequence for this analysis:** benign. The trigger efficiency is measured *from data*;
  MC supplies only the ΔR *correlation* correction, which is a ratio in which a per-leg
  efficiency normalization largely cancels. The residual is a systematic to quantify, not a bug.

- 2026-07-10 — **Step 8e: comparison REDONE with the EXACT NTuple-processing single-muon
  selection (user request). Supersedes the 8d numbers.**

  ⚠ **BUG IN THE 8d NUMBERS (mine, now fixed): `muon_pt` in the NTUP is SIGNED — the charge is
  carried in the sign** (`PythiaFullSimExtras.c:117` `m.pt = fabs(muon_pt->at(i))/1000`,
  `m.charge = (muon_pt->at(i) > 0)`). The 8d cut `muon_pt>4000` therefore kept only
  **positive-charge** muons, i.e. it was silently a same-sign-(++) dimuon selection.
  Verified: MC 10 142 muons with `muon_pt<0` vs 10 044 with `>0`; data 939 513 vs 985 264.
  **All 8d fiducial/turn-on numbers are void.** Use the table below.
  (`muon_trk_pt` is unsigned; `muon_pair_muon{1,2}_pt` is signed.)

  **Mirrored selection** — `DimuonDataAlgCoreT::PassCuts_DataCore` (data, `requireTight=true`)
  and `PythiaFullSimExtras::PassMuonMediumCuts` + `pass_tight` (MC), applied per muon:
  `quality & 305 == 305` (1 combined | 16 Tight | 32 IDCuts | 256 MuonCuts) · `|η| ≤ 2.4` ·
  `|pT|/1000 ≥ 4 GeV` · `|Δp/p| ≤ 0.12` · `|d0| < 2 mm` · `|z0·sinθ| < 2 mm`
  (`ParamsSet.h:370,375,376`; `turn_on_track_charge=false` in all three Extras headers).
  Event requires ≥2 such muons. Verified `muon_pair_muon{1,2}_*` == flat `muon_*` at
  `muon_pair_muonN_index` (0 mismatches / 3152 pairs), so the flat vectors are a faithful proxy.

  **Two code asymmetries found between the data and MC processing (both checked, both benign):**
  1. **Tight bitmask.** data `1|16|32|256` = 305; MC `pass_medium && (q&16)` = `1|8|16|32|256`
     = 313. **Equivalent:** `TrigRates.cxx:1196-1197` sets bit 8 for `quality<=Medium` and bit 16
     for `quality<=Tight`, and Tight ⊂ Medium ⇒ bit16 ⟹ bit8. Measured: **0** muons with bit16
     and not bit8, out of 13 874 (MC) and 1 576 769 (data).
  2. **Δp/p `fabs`.** data `fabs(dP_overP) > 0.12` → reject; **MC `dP_overP > 0.12` → reject
     (no `fabs`)** — MC keeps muons with `Δp/p < −0.12`. Affects 0.66% of otherwise-tight MC
     muons and 1.30% of data muons. Ran BOTH definitions applied to both samples: results agree
     to the 3rd decimal (integrated MC/data 1.23 either way). **Flagged as a probable bug in
     `PythiaFullSimExtras.c:41` (and `PowhegFullSimExtras`), not affecting this comparison.**

  **Result** (MC = isospin-weighted 6:6:9 over the complete `pn`,`np`,`nn` beams — the `pp` beam
  was still re-skimming; cross-section-weighted within each beam; data = pp24, 603 M events,
  8 907 591 events with ≥2 selected muons, 1 325 245 mu4-fired):

  | sublead μ pT [GeV] | f(2mu4\|mu4) MC | data | MC/data | f(noL1\|mu4) MC | data | MC/data |
  |---|---|---|---|---|---|---|
  | 4–5   | 0.5236 | 0.4199 | 1.25 | 0.8695 | 0.8514 | 1.02 |
  | 5–6   | 0.6920 | 0.5703 | 1.21 | 0.9466 | 0.9291 | 1.02 |
  | 6–8   | 0.7602 | 0.6406 | 1.19 | 0.9699 | 0.9554 | 1.02 |
  | 8–10  | 0.7822 | 0.6572 | 1.19 | 0.9753 | 0.9645 | 1.01 |
  | 10–15 | 0.7702 | 0.6513 | 1.18 | 0.9756 | 0.9662 | 1.01 |
  | 15–25 | 0.7078 | 0.6544 | 1.08 | 0.9728 | 0.9708 | 1.00 |
  | 25–50 | 0.5740 | 0.6882 | 0.83 | 0.9314 | 0.9812 | 0.95 |
  | **integrated** | **0.6240** | **0.5075** | **1.23** | **0.9134** | **0.8946** | **1.02** |

  MC per-beam integrated f(2mu4|mu4): pn 0.6235, np 0.6218, nn 0.6257 (spread <0.7% ⇒ still
  isospin-blind under the full selection).

  **Physics reading (cleaner than 8d, and it corrects 8d's interpretation):**
  - `f(mu4_mu4noL1 | mu4)` **cannot test L1**: `HLT_mu4_mu4noL1_L1MU3V` and `HLT_mu4_L1MU3V`
    seed the **same** `L1_MU3V` item, so the single L1 RoI is required in both numerator and
    denominator and cancels exactly. Observed MC/data = **1.02**, i.e. consistent with 1 — the
    HLT full-scan second leg is well modelled. (The 8d claim that this probes "one L1 leg"
    was wrong.)
  - `f(2mu4 | mu4)` isolates the **second muon leg**: `L1_2MU3V` needs a *second* L1 RoI on top
    of the denominator's one, plus a second HLT leg. Observed MC/data = **1.18–1.25** across the
    turn-on and plateau ⇒ **MC over-estimates the second-leg (L1 RoI × HLT) efficiency by
    ~20%.** This is the classic missing RPC/TGC chamber inefficiency in L1 muon simulation.
  - The 25–50 GeV bin flips to 0.83, but that bin has the least MC weight (high-pTHat slices
    carry `w ≈ 0.21 nb` vs `≈ 24–32 nb` at low pTHat) — **treat as MC statistics, not a trend**,
    pending the `pp` beam and a proper uncertainty.
  - **Consequence for the analysis: benign but quantified.** ε_trig is measured *from data*; MC
    supplies only the **ΔR correlation** correction, i.e. a *ratio* of two-leg efficiencies in
    which a per-leg normalisation largely cancels. The ~20% second-leg over-efficiency is
    exactly the quantity that must NOT be taken from MC in absolute terms, and it is the
    systematic to carry on the MC-based ΔR correction.

- 2026-07-10 — **R2 data-side menu names VERIFIED (and corrected).** User asked where the
  PbPb-data menu claim came from; it had been asserted without a file check. Read
  `TriggerMenuJson_{HLT,L1,HLTPS}` from the local PbPb23 test AOD
  (`dimuon_data/test_aod/data23_hi/AOD.41716150._000001.pool.root.1`, run 462240,
  physics_HardProbes; AthAnalysis env + `xAOD::MakeTransientMetaTree`): data HLT menu =
  `PhysicsP1_HI_run3_v1` (R2 had wrongly said `Physics_HI_run3_v1`, which is the data **L1**
  menu name), HLT PS sets `Physics_HI_{2e-05,1.6e-05}e32_880b`. R2 updated. The substantive
  systematic note stands: overlay MC simulates a Dev menu + validation prescales; data runs a
  PhysicsP1 menu + physics prescales.

## Results & Observations

*(organized, mutable)*

### R1. Trigger simulation — VERIFIED PRESENT in both configurations
See Progress Log Steps 2–3. AMI steering (`doRDO_TRIG` + `doTRIGtoALL`) **and** the produced
files agree. The `README.md` claim *"`HLT_MuonsCB_RoI`/`HLT_MuonsCB_FS` keys are data-only"* is
**FALSE for these Run-3 fullsim/overlay AODs** and must be corrected.

### R1b. ⚠ PHYSICS-AFFECTING SIDE EFFECT FOUND: `StoreAllEvents=False` + `UseTrigger=True`
`TrigRates.cxx:301-307`:
```cpp
if(m_use_trigger && m_is_evgen==false){
  bool pass_trigger_cuts=false;
  CHECK(ProcessTriggers(pass_trigger_cuts));   // true if ANY configured chain fired
  if(m_StoreAllEvents==false && pass_trigger_cuts==false) return StatusCode::SUCCESS;  // event DROPPED
}
```
Both CA configs set `alg.StoreAllEvents = False`. That is correct for the **data** skim (keep
only trigger-selected events), but with trigger enabled it turns the MC skim into a
**trigger-OR-filtered** sample. Measured on the pp-fullsim test (200 events read):
**166 stored, 34 dropped (17%)**, and `b_HLT_mu4_L1MU3V` = 100% by construction. The existing
(trigger-off) NTUP `Pythia_5p36TeV_pp_hQCD_DiMu_pTH8_14.FullSimPP24.NTUP.root` has all
**10 000** events.

This would silently bias exactly the two things this task is meant to unblock:
the **reco-efficiency denominator** and the **MC trigger-efficiency denominator** (you cannot
measure a trigger efficiency on a trigger-selected sample).

**Fix applied (pending user confirmation): `alg.StoreAllEvents = (is_MC and mc_has_trigger_sim)`**
— MC keeps every event and merely *records* the decision; data behavior is bit-identical
(`is_MC` False ⇒ `False`, as before). Consistent with §4 negative constraint: "the skim writes
trigger decisions and matching; it does not apply a trigger cut."

### R2. Trigger menus differ between MC and data (note for systematics)
- pp fullsim r16578 → HLT `PhysicsP1_pp_lowMu_run3_v1`; pp24 data → the pp-reference physics menu.
- overlay r17618/r17662 → HLT **`Dev_HI_run3_v1`** (a *development* menu), L1 `MC_HI_run3_v1`,
  prescale set `Dev_HI_run3_v1_TriggerValidation_prescale`; PbPb23 data → HLT
  **`PhysicsP1_HI_run3_v1`**, L1 `Physics_HI_run3_v1`, HLT prescale sets
  `Physics_HI_{2e-05,1.6e-05}e32_880b`. *(Data side verified 2026-07-10 by reading
  `MetaData` `TriggerMenuJson_{HLT,L1,HLTPS}` from the local test AOD
  `dimuon_data/test_aod/data23_hi/AOD.41716150._000001.pool.root.1`, run 462240
  physics_HardProbes — an earlier version of this entry asserted the data HLT menu was
  `Physics_HI_run3_v1` without a file check; that is actually the data L1 menu name.)*
The `mu4`/`2mu4`/`mu4_mu4noL1` chain *definitions* are the same and all run **unprescaled**
(PS=1) in MC, but the Dev-vs-Physics menu difference should be recorded as a possible
trigger-efficiency systematic when the MC-based ΔR correction is derived.

### R3. Dataset census (mc23_5p36TeV, `Py8EG_A14_*_hQCD_DiMu_pTH*`)
Two distinct families exist; **this changes what "all DSIDs" means** (see Open Question O1).

**(a) Existing FullSimPP24 family — DSIDs 802758–802781 (24 = 4 beams × 6 pTHat), r16578.**
`recon.AOD.e8599_s4521_r16578` and `recon.AOD.e8599_s4521_s4483_r16578` are the SAME files
(identical GUIDs; two AMI names). **10 files, 6.5 GB, 10 000 events each** → 24 × 10k = 240k
events, ~156 GB total. These are the AODs behind the 24 NTUPs in
`~/usatlasdata/pythia_fullsim_test_sample/`. Beams `{pp,pn,np,nn}` are combined with the
isospin ratio 4:6:6:9 (`isospin_weight.md`, `PythiaAlgCoreT.c`).

**(b) NEW `_pdf` family — DSIDs 803015–803020 (6 = pp beam only × 6 pTHat), also r16578.**
`recon.AOD.e8599_e8586_s4521_s4483_r16578`. **1 200 000 events per slice** (803018: 120 files,
941 GB). recon.AOD exists for **803015, 803018, 803019, 803020** (pTH125_300, 40_70, 70_125,
8_14); **803016 (pTH14_24) and 803017 (pTH24_40) are still in production** (only
`NA.recon.AOD.<taskid>_sub*` sub-datasets so far).
Different cross-section from the 8027xx pp-beam sample (803018: σ=13114 pb, ε=5.2629e-4 vs
802779: σ=13860 pb, ε=5.5694e-4) ⇒ a **different generation setup**, not a statistics extension
of 802779. Only the `pp` beam exists → the 4:6:6:9 isospin combination **cannot** be formed from
this family alone.

### R4. Storage constraint
`mmlsquota atlasgpfs01` → `usatlast3-data` at **~1.47 TB used of a 1.5 TB soft quota (~96%)**
(hard limit ~1.79 TB). Re-skimming family (a) reproduces ~5 GB of NTUPs (fine). Re-skimming
family (b) would download O(4–5 TB) of AOD on the grid (fine, grid-side) but produce NTUPs
scaled by 120× events → O(100+ GB), which does **not** fit comfortably. Must be planned.

### O1. OPEN QUESTION — which DSIDs does "all DSIDs" mean?
The user's example DSID (`803018 ..._pdf`) belongs to family (b), but the 24 NTUPs we actually
use come from family (a). These are different productions (different σ, different beams,
120× statistics, family (b) incomplete). **Blocking for Step 6 (grid re-skim) only** —
Steps 4–5 (code + test jobs) are DSID-independent and proceed.

## Remaining Work

- Step 6: 31 grid tasks in flight → download/hadd/validate (grid_monitor, both modes).
- Step 8: post-skim sanity check (10 000 entries per NTUP; trigger branches non-empty).
- Step 7 leftovers: roadmap Q4 row ("dR trigger-correlation correction ... needs fullsim overlay
  with trigger sim") → mark the prerequisite satisfied once the NTUPs land.
- **Follow-up (out of scope here, now unblocked):** derive the MC-based ΔR trigger-correlation
  correction from `dimuon_b_HLT_2mu4_L12MU3V_mu{1,2}passLeg{1,2}_dR_*`; consume the new trigger
  branches in `NTupleProcessingCode` / RDF.
- **Housekeeping:** `/usatlas/u/yuhanguo/usatlasdata/aod_trigger_check/` holds ~29 GB of AODs
  downloaded for the file-level verification. They are re-downloadable and can be deleted
  (`usatlast3-data` quota was at ~96%).
- **Open (not blocking):** the new `_pdf` production 803015–803020 (1.2M events/slice, pp beam
  only, 4/6 slices reconstructed) is a *different* generation setup (σ 13114 vs 13860 pb for
  pTH40_70) and cannot form the 4:6:6:9 isospin combination alone. Its role vs the 24-DSID
  family needs a physics decision before it can be used (roadmap Q4 "Full Pythia fullsim pp24
  sample").

## Latest Stage

**2026-07-09 — Steps 1–5 DONE, Step 6 IN FLIGHT.**

Answer to the original question: **YES — trigger simulation is present in every Run-3 fullsim
sample we use.** Confirmed twice over (AMI steering `doRDO_TRIG`+`doTRIGtoALL` on r16578,
r17618 and r17662; and the trigger containers + menu read straight out of one AOD per tag).
Everything needed for skimming (`mu4`, `2mu4`, `mu4_mu4noL1`) **and** a full 2mu4
trigger-efficiency evaluation (per-leg matching via `HLT_MuonsCB_RoI` / `HLT_MuonsCB_FS`)
is there, all unprescaled. `r15970` is only the AODMerge tag — the trigger comes from the
recon tag, so the "trigger might be in r15970" worry is moot.

Trigger is now ON by default for all Run-3 fullsim/overlay MC (code `434de17`,
`/review-analysis-code` PASS). The one real trap — `StoreAllEvents=False` silently converting
the MC skim into a trigger-OR-filtered sample — was caught, escalated, and fixed (R1b).

**PAUSE POINT 2026-07-10.** Overlay 7/7 complete + validated. pp fullsim **23/24** downloaded &
validated (10 000 entries each); **only task `51360141` (pp beam, pTH8_14) is still running on
the grid at 0%**. Step 8 sanity checks (a)–(e) all done and committed.

**To resume:**
1. `ps -eo cmd | grep grid_monitor` — one `--mode fullsim_pp` instance should still be polling
   task `51360141` (launched with `nohup`, survives the shell). If it died:
   `cd SkimCode/scripts && ./grid_monitor.sh --mode fullsim_pp -i 10 51360141`
2. When it lands: confirm 10 000 entries + trigger branches, then re-run the Step-8e comparison
   **including the `pp` beam** so the isospin 4:6:6:9 combination is complete
   (scripts kept: `/tmp/claude-101379/{sel.h,sel_mc.py,sel_data.py}`; note `/tmp` may be wiped —
   the selection is fully specified in the 8e Progress-Log entry above and can be rebuilt).
3. Then: mark roadmap Q4 prerequisite ("unbiased trigger decision in MC") satisfied; update
   `merging-record.txt` bookkeeping (auto by grid_monitor); consider moving the now-superseded
   `.bak_20260709.root` files (84 GiB overlay + ~9 GiB pp) to pnfs — **safe now**, every
   canonical file exists (but re-verify per `project_grid_monitor_bak_rename`).
4. Open follow-up (separate task): fix the missing `fabs` on `Δp/p` in
   `PythiaFullSimExtras.c:41` / `PowhegFullSimExtras` so MC matches
   `DimuonDataAlgCoreT.c:597`. Affects 0.66% of tight MC muons; does not change any result here.

**Next physics action (unblocked):** derive the MC-based ΔR trigger-correlation correction from
`dimuon_b_HLT_2mu4_L12MU3V_mu{1,2}passLeg{1,2}_dR_*`, carrying the ~20% second-leg MC
over-efficiency (Step 8e) as a systematic on the correction.

Monitor state files:
`~/usatlasdata/pythia_fullsim_test_sample/grid_monitor_{status.log,state.txt}` and
`~/usatlasdata/pythia_fullsim_hijing_overlay_test_sample/grid_monitor_{status.log,state.txt}`.
