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
- overlay r17618/r17662 → HLT **`Dev_HI_run3_v1`** (a *development* menu) with the
  `..._TriggerValidation_prescale` prescale set; PbPb data → `Physics_HI_run3_v1`.
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

**Next action when the grid finishes:** run Step 8 (sanity check: every NTUP must have 10 000
entries and non-empty `b_HLT_*`), let `grid_monitor` update both `merging-record.txt` files,
then mark the roadmap Q4 prerequisite ("unbiased trigger decision in MC") satisfied.

Monitor state files:
`~/usatlasdata/pythia_fullsim_test_sample/grid_monitor_{status.log,state.txt}` and
`~/usatlasdata/pythia_fullsim_hijing_overlay_test_sample/grid_monitor_{status.log,state.txt}`.
