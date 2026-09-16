# HIJING overlay test sample, Pb+Pb 2024 conditions (r17864): skim, download, verify

## Objective

Skim the new **Pb+Pb-2024-conditions** HIJING-overlay test sample (one DSID, pTH125_300,
r17864) with the existing overlay `TrigRates` skim, download the NTUP into
`~/usatlasdata/pythia_fullsim_hijing_overlay_test_sample/` (the old pbpb23-conditions
sample now lives in `..._test_sample_pbpb23/`), and define the checks that decide whether
the r-tag is acceptable as the recipe for the full 4-beam × 6-slice production.

## Autonomy Contract (DONE 2026-09-16 — all Done items met)
- Mandate: run autonomously to DONE; do NOT pause to confirm progress. Finishing a
  plan, a passing small test, or one pipeline stage is NOT a stopping point.
- Done = (1) skim code compiled in `build_25`; (2) 100-event local test on a BNL replica
  passes `check_skim_output.C` vs the pbpb23 pTH125_300 NTUP; (3) grid task submitted from
  a NEW `grid_sub_pbpb24_test_sample.sh`; (4) `grid_monitor.sh --mode overlay` has
  downloaded + hadd-ed the NTUP to
  `~/usatlasdata/pythia_fullsim_hijing_overlay_test_sample/Pythia_5p36TeV_pp_hQCD_DiMu_pTH125_300.FullSimHIJINGOverlayPP24.NTUP.root`
  with 10 000 entries; (5) the dataset's OWN AMI info + r-tag info saved under that dir;
  (6) a written verification proposal for the user (the verification itself is NOT in
  scope until the user picks).
- Stop-and-ask = ANY physics-results-bending ambiguity (no fixed list; use judgment;
  when unsure whether an ambiguity is blocking, treat it as blocking → AskUserQuestion).

## Context

- Previous overlay test sample: 6 slices, **Pb+Pb 2023 conditions**, r17618 over r15970,
  evgen e8599 (DSIDs 802776-802781), skimmed June/July 2026
  (`hijing_overlay_r17618_grid_reprocessing.md`; `grid_sub_pbpb23_test_sample.sh`). Now
  moved to `pythia_fullsim_hijing_overlay_test_sample_pbpb23/`.
- New sample (one DSID only):
  `mc23_5p36TeV.802776.Py8EG_A14_pp_hQCD_DiMu_pTH125_300.merge.AOD.e8613_e8586_s4684_r17864_r17855_tid52480369_00`
  — 10 files, 10 000 events, 54.9 GB (pbpb23 one was ~84 GB / 10k). Full replica only at
  UKI-LT2-QMUL; 1 file (`_000026`) at BNL-OSG2_DATADISK.
  Requested by the production contact: vertex z = -71.2 mm, impact parameter 0-5 fm,
  2024 run conditions.
- **AMI (fetched 2026-09-15, this dataset):** σ = 89.541 nb, genFiltEff = 2.342314e-3,
  σ·ε = 0.20973 nb. The evgen tag is NEW (`e8613_e8586`, was `e8599`) while the DSID is the
  same 802776 → `expected_ami_dsids` would NOT distinguish the productions; the AMI file
  for this sample MUST live in the new dir's `ami_info/` and be re-fetched (done, step 2).
  Numerically it differs from registry A (89.529 / 2.342579e-3) by 1e-4 relative.
- Uncommitted `TrigRates.cxx/.h` ZDC-policy edits in the working tree belong to the
  pbpb26 data skim session. They only act when `StoreZdc>0`; the overlay mode sets
  `StoreZdc = 0`, so they are inert for this skim. Not committed here.

## Procedure

1. AMI: `ami show dataset info <ds>` + `ami show tag r17864` / `r17855` → save to
   `~/usatlasdata/pythia_fullsim_hijing_overlay_test_sample/{ami_info/,r17864.txt,r17855.txt}`.
2. Compile: `source setup_25.sh` (runs `acm compile` in `build_25`).
3. Local test: `run_test.sh r17864 orig root://<BNL PFN>` (100 events,
   `TRIGRATES_RUNMODE=ppmcfullsim_hioverlay24`), then
   `check_skim_output.C(test, pbpb23 pTH125_300 NTUP)`.
4. Grid: new `grid_sub_pbpb24_test_sample.sh` — same trf, OUT_TAG `FullSimHIJINGOverlayPP24`
   (frozen file tag), VER_TAG bumped to `Sep2026.v1`, inDS = the r17864 container.
5. Monitor: `nohup scripts/grid_monitor.sh --mode overlay <taskid>` (writes flat into the
   new test-sample dir) + an artifact-polling waiter.
6. Verification proposal → user.

## Implementation Plan

1. [x] AMI + tag info saved (Procedure 1)
2. [x] compile (2)
3. [x] local 100-event test + branch check (3)
4. [x] grid submission (4) — jediTaskID **52568863**
5. [x] monitor/download/hadd (5) — done 2026-09-16 06:20, 10 000 entries, 22.4 GB
6. [x] verification proposal (6) — delivered 2026-09-15; awaiting user decision

## Progress Log

### 2026-09-15 steps 1-4
- Step 1: `~/usatlasdata/pythia_fullsim_hijing_overlay_test_sample/ami_info/ami_info_mc23_5p36TeV_Py8EG_A14_pp_hQCD_DiMu_pTH125_300.txt`
  (σ 89.541 nb, genFiltEff 2.342314e-3, DSID 802776) + `r17864.txt`, `r17855.txt`, `s4684.txt`, `e8613.txt`.
  r17864: beam spot (-0.7,-0.6,-71.2) mm, σ_z = 0.01 mm (pinned), `Campaigns.MC23HeavyIons2024NoPileUp`,
  bunch structure `BunchStructureHeavyIon2022`, HIJING HITS `860250.Hijing_PbPb_UCC_Flow_JJFV6_ip0_5.e8613_s4684_s4688`,
  r17855 = "clone r15970 in 24.0.90", steering doRDO_TRIG+doTRIGtoALL. AOD run number **488600** (was 460000 for r17618).
- Step 2: `setup_25.sh` → `acm compile` in `build_25`, up to date (working-tree ZDC edits compiled in; inert, StoreZdc=0).
- Step 3: `run_test.sh r17864 orig root://dcgftp.usatlas.bnl.gov:1094//pnfs/usatlas.bnl.gov/BNLT0D1/rucio/mc23_5p36TeV/62/c2/AOD.52480369._000026.pool.root.1`
  → `test_r17864/test_r17864_orig.root`, 100 events, exit 0. Same `meTrk Not Found for Combined muon` ERROR lines as r17618 (212 vs 369/248 ev).
  `check_skim_output.C` vs pbpb23 pTH125_300 NTUP: only diff = the 3 `_VTE50` muon chains (2023 menu) are **not in the 2024 HI menu**
  (it has `_VjTE50`); analysis chains mu4 / 2mu4 / mu4_mu4noL1 present and firing (fill 0.93 / 0.47 / 0.73 in 100 ev).
  **Config change:** `run_pythia_fullsim_HIJING_overlay/TrigRates_CA.py` overlay `Muon_triggers` list switched from the hi2023 list
  to the hi2024 data list (adds mu8/mu10 and the 2024 `mu4noL1_hi_ucc*` / `mu4noL1_L1ZDC_HELT*` chains; the L1ZDC ones are always
  empty in the overlay, no ZDC sim). Nothing downstream reads the VTE50 branches (grep of Analysis/: only the Run-2 mb chain in TrigEff).
  Re-test: 0 "Is Not Configured" warnings; 24 missing / 64 extra branches vs pbpb23, all trigger-list bookkeeping.
- Step 4: `grid_sub_pbpb24_test_sample.sh` (VER_TAG `Sep2026.v1`, single comma-separated `--excludeFile`, root/log/test dirs excluded)
  → `user.yuhang.NTUP.Pythia_5p36TeV_pp_hQCD_DiMu_pTH125_300.FullSimHIJINGOverlayPP24.Sep2026.v1.`, **jediTaskID 52568863**, 14:16.
- Step 5 started 14:17: `nohup scripts/grid_monitor.sh --mode overlay -i 10 52568863` (PID in `$D/grid_monitor_r17864.pid`,
  stdout `$D/grid_monitor_r17864.out`, status `$D/grid_monitor_status.log`); waiter polls the status log for
  "All tasks resolved. Worker exiting.".

## Results & Observations

### R0. Full 10k NTUP (2026-09-16)
`~/usatlasdata/pythia_fullsim_hijing_overlay_test_sample/Pythia_5p36TeV_pp_hQCD_DiMu_pTH125_300.FullSimHIJINGOverlayPP24.NTUP.root`
— 10 000 entries, 22.4 GB (pbpb23 analogue 22.1 GB). `check_skim_output.C` vs pbpb23: 24 missing / 64 extra branches, all the
trigger-list bookkeeping of step 3; 17 always-empty branches = the 5 known MC ones + the 12 L1ZDC-seeded `mu4noL1` chain
branches (no ZDC sim). All muon / truth / vertex / FCal branches filled.

R1 numbers on the FULL sample (10k vs 10k): FCal ΣE_T quantiles 16/50/84 = **-0.77 / 0.41 / 2.07 TeV**, **35.2 %** of events
negative, 87.6 % below 2.5 TeV (pbpb23: 3.79 / 4.40 / 5.11, 0 %, 0 %). Truth-matched muons/event (pT>4, |η|<2.4, prob>0.5,
|id|=13) **0.674 vs 0.751 (-10 %, n = 6741 vs 7510, stat ±1.2 %)**; non-matched 1.50 vs 0.755 (×2, Δz(PV) 7.0 ± 21.2 mm vs
0.40 ± 1.07 mm); MS η/φ hits of matched muons 2.79/3.04 vs 2.85/3.10; pT/pT_truth rms 0.090 vs 0.074. The muon system is
degraded as well, not only the calorimeter.

### R1. r17864 calorimeter response is broken (found 2026-09-15 on the 100-event local test; confirmed in the AOD)

Evidence (skim NTUP `test_r17864/test_r17864_orig.root`, 100 ev, vs pbpb23 r17618 pTH125_300 NTUP, 2000 ev):

| quantity | r17864 (2024 cond.) | r17618 (2023 cond.) | expectation for b < 5 fm |
|---|---|---|---|
| FCal ΣE_T quantiles 16/50/84 [TeV] | **-0.66 / 0.49 / 1.70** | 3.79 / 4.39 / 5.11 | ≳ 3 TeV (0-10 %) |
| fraction of events with FCal ΣE_T < 0 | **33 %** | 0 | 0 |
| `centrality` (2023 calib.) | mean 23 %, up to 77 %, 33 % undefined (-1) | mean 2 %, max 14 % | ≤ ~10 % |
| `trk_numqual` (ID tracks) | 3546 ± 1686 | 3576 ± 1413 | same ✓ |
| `vtx_ntrk[0]` | 4787 ± 1240 | 4609 ± 1118 | same ✓ |
| truth particles / event | 7.0e4 | 6.9e4 | same ✓ |
| reco muons / ev, pT>4, |η|<2.4, truth-matched μ (prob>0.5) | 0.65 (n=65) | 0.74 | ~same (1σ) |
| reco muons / ev, pT>4, |η|<2.4, NOT truth-matched | **1.55** | 0.79 | — |
| `muon_trk_vz - vtx_z[0]` of the not-matched muons [mm] | 6.4 ± 20.4 | 0.38 ± 1.06 | — |
| MS η-hits per truth-matched muon | 2.66 ± 1.30 | 2.87 ± 1.11 | — |
| `muon_pt/muon_truth_pt` (matched) | 0.997 ± 0.035 | 0.999 ± 0.135 | — |

Directly in the AOD (`HIEventShape`, script `scratchpad/hies.py`, 8 events of `AOD.52480369._000026`):
FCal 0.42-3.37 TeV, all-calo ΣE_T **-2.4 … +20.3 TeV**, 22-359 of the 562 (layer × η) slices negative per event;
`CaloSums` ev0: ALL 2.4, EMCal 1.56, HCal 0.42, FCal 0.42 TeV. Reference 2023-cond. overlay AOD
(`..._pbpb23/test_AOD_pythia_fullsim_hijing.root`): FCal 3.35-4.52 TeV, all-calo 13-18 TeV, 5-14 negative slices;
`CaloSums` ev0: ALL 13.4, EMCal 8.0, HCal 2.0, FCal 3.35 TeV. Per-layer ev0: every LAr AND Tile layer is 3-10× lower in r17864
(Tile 12: -3.5 vs 130 GeV; EMB2 (2): 185 vs 1388; FCAL0 (21): 312 vs 2701).

Interpretation: the ID sees the full HIJING event (same track multiplicity), so the generator/simulation input is right and the
problem is in digitisation/reconstruction of the calorimeters (LAr + Tile) — energy far too low AND large negative sums.
Negative reconstructed E_T across LAr is the signature of sampling the bipolar pulse off-peak (undershoot), i.e. a
signal/readout **timing (BCID) offset**. The only stated change from r17835 → r17864 is
`HITtoRDO: Digitization.PU.BunchStructureConfig = 'RunDependentSimData.BunchStructureHeavyIon2022'` ("bunch structure added"),
which is the natural suspect. Not proven here — needs the production contact / an r17835-style clone without it.
The doubled non-truth-matched muon rate with |Δz| ~ 20 mm from the PV and the slightly lower MS hit count are consistent with
the same timing problem reaching the muon system, but that is a weaker, 100-event statement.

Consequence for us: FCal-based centrality and any calorimeter-dependent muon quantity (calo-tagged muons, isolation, E-loss)
are unusable in this r-tag; reco-efficiency / det-response derived from it would be wrong in the 2024 conditions.

## Remaining Work

- User decision on the verification plan / on reporting R1 to the production contact (message draft offered).
- The verification itself (items 1-6 of the proposal) is NOT started; item 1 (calorimeter) is already answered negatively by R1.

## Latest Stage

2026-09-16: skim + download DONE (task 52568863). R0/R1 recorded. Waiting for the user's decision on verification scope and on reporting the r17864 calorimeter/muon degradation to the production contact.