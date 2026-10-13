# HIJING overlay test sample, Pb+Pb 2024 conditions (r17864): skim, download, verify

## Objective

Skim the new **Pb+Pb-2024-conditions** HIJING-overlay test sample (one DSID, pTH125_300)
with the existing overlay `TrigRates` skim, download the NTUP into
`~/usatlasdata/pythia_fullsim_hijing_overlay_test_sample/` (the old pbpb23-conditions
sample now lives in `..._test_sample_pbpb23/`), and run the checks that decide whether
the r-tag is acceptable as the recipe for the full 4-beam × 6-slice production.

**R-tag succession of the same (vtx z = −71.2 mm, b 0–5 fm) configuration** — the NTUP name is
configuration-keyed, so all of them share it and are told apart only by directory / VER_TAG:

| r-tag | AOD merge tid | skim VER_TAG / task | verdict |
|---|---|---|---|
| r17864 | 52480369 | Sep2026.v1 / 52568863 | **FAIL** — calorimeter broken (R1, R4, R9) |
| r17864 | 52629267 ("reproduction") | Sep2026.v3 / 52687659 | **FAIL, bit-identical re-run** (R7, R8) |
| **r17913** | **52733789** | **Sep2026.v4 / 52785715** | session 6, in progress (R10–) |

## Autonomy Contract (DONE 2026-09-29 — all six Done items met; session 6, r17913)
- Mandate: run autonomously to DONE; do NOT pause to confirm progress. Finishing a
  plan, a passing small test, or one pipeline stage is NOT a stopping point.
- Done = (1) r17913 provenance saved under
  `~/usatlasdata/pythia_fullsim_hijing_overlay_test_sample/` (AMI r-tag record + the
  r17864→r17913 preExec diff + rucio dataset record), and the r17864 NTUPs archived so the
  configuration-keyed NTUP name is free; (2) 100-event local test on a downloaded AOD file →
  first FCal / L1TE / centrality verdict; (3) grid skim submitted, downloaded and hadd-ed by
  `grid_monitor.sh --mode overlay` to the canonical
  `Pythia_5p36TeV_pp_hQCD_DiMu_pTH125_300.FullSimHIJINGOverlayPbPb24.vtxz-71_2mm_b0_5fm.NTUP.root`
  with 10 000 entries, branch-checked vs the pbpb23 pTH125_300 NTUP; (4) the FULL event-level
  check set of `plots/r17864_rtag_sanity/` regenerated for r17913 into
  `plots/r17913_rtag_sanity/`: `fcal`, `vtx`, `ntrk`, `l1te`, `centr`, `bcid`, `trainpos`;
  (5) the efficiency-chain half of that same check set (single-µ reco ε, mu4 ε, pT/η response,
  `from_same_b`) vs the pbpb23 pTH125_300 slice — i.e. the R2 redo, whose only blocker was the
  broken FCal that r17913 targets — OR, if the FCal is still broken, explicitly skipped per the
  2026-09-17 user ruling and said so; (6) tracking doc + INDEX scope updated.
- Stop-and-ask = ANY physics-results-bending ambiguity (no fixed list; use judgment;
  when unsure whether an ambiguity is blocking, treat it as blocking → AskUserQuestion).

*(Session 1's contract, DONE 2026-09-16, is preserved in the Progress Log.)*

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
- **AMI — CORRECTED 2026-09-17 (user ruling; see §"AMI-weight correction" below).** The
  2026-09-15 fetch on the AOD dataset (σ = 89.541 nb, genFiltEff = 2.342314e-3) was WRONG in
  principle: `e8613` is the HIJING evgen tag, not a new Pythia evgen. The Pythia evgen of this
  sample is 802776 e8599 (nPDF, registry table A: 89.529 nb / 2.342579e-3), whose weight it now
  carries. The AOD record is kept as `ami_info_AOD_record_20260915/`, read by nothing.
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

### Session 6 (2026-09-29, r17913)

7. [x] AMI/provenance for r17913 + r17864→r17913 preExec diff (§R10)
8. [x] archive the two r17864 NTUPs so the configuration-keyed name is free (§R10)
9. [x] `grid_sub_pbpb24_test_sample.sh` → VER_TAG `Sep2026.v4`, inDS `…r17913…tid52733789_00`;
       submitted **jediTaskID 52785715**; `grid_monitor.sh --mode overlay -i 10` launched
10. [x] 100-event local test on one downloaded AOD → first verdict (§R11) — **FCal FIXED**
11. [x] grid download + hadd (10 000 entries, 21.66 GB) + branch check vs pbpb23 (§R12)
12. [x] full event-level check set → `plots/r17913_rtag_sanity/` (§R12) — **R9 signature GONE**; /review-plot
13. [x] efficiency-chain comparison vs the pbpb23 pTH125_300 slice (R2 redo) (§R13) — **deficit GONE**

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

### 2026-09-29 session 6 (r17913) — final stage summary, moved here from Latest Stage on 2026-10-06

**2026-09-29 (session 6) DONE — r17913 PASSES every check. The r-tag is acceptable.**

Sequence: provenance (R10 — the r17864 → r17913 `HITtoRDO` preExec diff is the whole change and is a
direct answer to R9) → r17864 NTUPs archived so the configuration-keyed name is free → skim
`jediTaskID 52785715` (VER_TAG Sep2026.v4) → 100-event local verdict (R11) → download + hadd
(10 000 entries, 21.66 GB) + branch check → full event-level set (R12) → NTP + RDF chain and the
R2 redo (R13).

**Both halves pass.**
- **R12, calorimeter:** FCal median 4.29 TeV with **0.0 % negative** (r17864: 0.41 TeV, 36.2 %),
  L1 ΣE_T 15.70 vs 16.14 TeV, centrality population **identical to the reference to two decimals**
  (99.32 % / 0.68 %, 0 % undefined). **The R9 signature is gone**: the mod-8 train-position swing is
  0.074 TeV (r17864: 3.58 TeV), flat to the same ±0.02 TeV as r17618, and the 2 isolated bunches —
  the only healthy crossings under r17864 — now agree with the in-train ones.
- **R13, reconstruction:** Tight reco ε centrality-blind **0.6518 vs 0.6595 (ratio 0.988**, was
  0.885), mu4 ε 1.017 (was 0.975), p_T-resolution RMS +1.3 % (was +8 %), and the Tight OS pair count
  **2997 vs 2998** (was 2389, −20 %). R2's deficits were conditions artefacts and all three closed.

Two residuals, both understood and neither a defect of the r-tag: a uniform ~2 % FCal/L1TE offset vs
the 2023 reference (cross-year comparison, vertex −71.2 vs −3.3 mm, percent-level by R2's own Δη
estimate), and a −6 % HITight track count. The 5–10 % efficiency ratio of 0.955 follows from the
first as bin migration; the centrality-blind rows are immune to it. **A 2024-conditions FCal
centrality calibration is the prerequisite before pbpb24 efficiencies are used for physics.**

Open for the user: (i) the R5 evgen decisions; (ii) whether to send the production contact a
confirmation that the six-flag preExec recipe works; (iii) green light for the full
4-beam × 6-slice × 5-vertex × 4-b production.

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

### R1b. Notes from the FCal / vertex plot review (2026-09-16)
- FCal: sides A (η>0) and C (η<0) are affected identically in r17864 (per-side curves overlay; min -1.56 TeV on both) — a symmetric, not single-side, readout failure.
- The r17618 REFERENCE FCal distribution has a ×1.5 step between 3.9 and 4.0 TeV (+4.6σ in the 4.0-4.1 TeV bin, present on both sides, not from HIJING event reuse: 10 000/10 000 distinct truth fingerprints). A ~1 % localized feature of the r17618 production (HIJING b-sampling / HITS composition); irrelevant to R1 but worth knowing before r17618 FCal is ever used for a centrality calibration.
- Vertex: r17864 reco PV residual mean -5.4 µm (σ_mean 0.1 µm → highly significant; absent in r17618) and RMS 11.1 vs 9.1 µm (22 % worse). Negligible for muon-pair physics; a genuine 2024-conditions difference (z = -71 mm and/or the same conditions problem reaching the ID).

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

## Verification round (user request 2026-09-16)

User-selected checks, each vs the pbpb23 (r17618) sample restricted to the SAME pT-hat 125-300 slice:
1. [x] FCal ΣE_T distributions, both samples → `plots/r17864_rtag_sanity/fcal_sum_et_r17864_vs_r17618_pTH125_300.png`
2. [x] reco vs truth vertex z, r17864 → `plots/r17864_rtag_sanity/vertex_z_reco_vs_truth_r17864.png`
   (reco PV mean -71.2054 mm on a truth point of -71.2000; residual -5.4 ± 11.1 µm vs -0.2 ± 9.1 µm for r17618)
3. [x] muon reco efficiency + detector response — full chain run (`plots/hijing_overlay_pbpb24_{reco_effcy_plots,det_resp_plots,reco_effcy_plots_require_signal_cuts}/`, 300+72+300 PNGs) + slice-matched comparison `plots/r17864_rtag_sanity/single_muon_reco_effcy_vs_pt_ctr{0_5,5_10}_r17864_vs_r17618.png`, `single_muon_pt_eta_response_ctr0_5_r17864_vs_r17618.png`
4. [x] MC single-muon mu4 efficiency (Step-1 mirror, 0-5 %) — `single_muon_mu4_effcy_vs_pt{,_coarse}_ctr0_5_r17864_vs_r17618.png`
5. [x] `from_same_b` fraction — `from_same_b_fraction_r17864_vs_r17618.csv`
Code: `plotting_codes/overlay_rtag_checks/plot_r17864_event_level.cxx` (items 1-2; raw skim event-level branches, no
processing procedure involved). Items 3-5: the sample-dir-layout refactor (`fullsim_sample_dir_layout.md`) landed
2026-09-16 ~18:00 (commits 2d78e97 / edd66ba / 14d6067) with the `overlay_pbpb_year` knob, so the chain runs as-is.

**User decisions 2026-09-16 (centrality):** FCal distribution without centrality cut; reco AND trigger efficiency keep
the nominal FCal-centrality selection (0-5 % for the MC trig-eff Step 1, per-centrality bins for reco-eff) — NO code
change. Consequence for r17864 (broken FCal): only 4.2 % of events are "0-5 %" and 37.4 % have centrality = -1.

Centrality-bin population (skim `centrality`, 2023 FCal calibration; pTH125_300, 10 000 events each):

| bin | pbpb23 (r17618) | pbpb24 (r17864) |
|---|---|---|
| 0-10 % | 99.32 % | 7.59 % |
| 10-20 % | 0.68 % | 7.89 % |
| 20-30 % | 0 | 8.19 % |
| 30-40 % | 0 | 9.22 % |
| 40-50 % | 0 | 12.57 % |
| 50-80 % | 0 | 16.31 % |
| 80-100 % | 0 | 0.80 % |
| undefined (-1: below the calibration's 100 % edge, incl. all negative FCal) | 0 | 37.43 % |
| (0-5 % / 5-10 % sub-split) | 81.54 / 17.78 % | 4.20 / 3.39 % |

**Chain runs (all with the post-refactor code, same day):**
- r17864: `OVERLAY_YEAR=24 pipeline_pythia_fullsim_overlay.sh hijing` (Condor NTP → RDF → reco-eff/det-resp plots;
  5 min, single slice) + `run_pythia_fullsim_overlay_{mc_trig,single_muon_mc_trig,single_muon}.sh` locally. Outputs flat in
  the sample dir (`muon_pairs_..._pbpb24_no_data_resonance_cuts{,_mc_trig,_mc_trig_single_muon,_single_muon}.root`,
  `histograms_..._pbpb24_...root`, `plots/hijing_overlay_pbpb24_{reco_effcy_plots,reco_effcy_plots_require_signal_cuts,det_resp_plots}/`).
  Sanity: OS pairs 6689, SS 2052 (kin5 only, kin0-4 empty); mc_trig single-muon tree 13 596 reco-matched truth muons
  (pbpb23 slice: 15 205 → -10.6 %, the reco-eff deficit again). AMI read from the sample's own `ami_info/`
  (w_factor 2.09733e-05 = 0.20973 nb / 10 000). The 5 missing slices are SKIPPED via `allow_missing_slices`, which the four
  overlay run scripts + the Condor worker now set for `overlay_year == 24` ONLY (TEMPORARY — delete when the 6-slice pbpb24
  production exists).
- pbpb23 reference restricted to the SAME slice: `..._pbpb23/r17618_pTH125_300_run/` (symlink to the pTH125_300 NTUP,
  `fullsim_input_dir_override`, outputs suffixed `_r17618pTH125_300`), driven by the new
  `NTupleProcessingCode/run_pythia_fullsim_overlay_pbpb23_pTH125_300_slice.sh {pair,single_muon,mc_trig,mc_trig_single_muon}`.
  Needed because the existing pbpb23 products are 6-slice AND of mixed vintage (Jul 22 / Aug 3 / Aug 13 / Sep 9).
- Comparison macro: `plotting_codes/overlay_rtag_checks/plot_r17864_vs_r17618_pTH125_300.cxx` (recoeff | trigeff | detresp | samb),
  selections mirrored from `plot_single_muon_reco_effcy_r17618_vs_r17662.cxx` and `FillMCTrigEffHists.cxx` Step 1.

### R2. Slice-matched comparison r17864 vs r17618 (pTH125_300, same NTP code, 2026-09-16; reviewer-verified numbers)

**USER RULING 2026-09-17 — R2 is NOT studiable until the FCal is fixed.** The reco/trigger efficiencies must be compared in the
SAME centrality bin, but the r17864 centrality is derived from a mis-reconstructed FCal, so its "0-5 %" / "5-10 %" bins are not the
same event populations as r17618's; the bin-to-bin deficits below (in particular the 5-10 % ones) may be centrality mis-assignment,
not a reconstruction-efficiency problem. Keep the table as a record of what was run; draw no r-tag conclusion from the
centrality-binned rows. The centrality-blind rows do not depend on the FCal but still describe a sample whose conditions are
broken. Redo R2 once a fixed-calorimeter production exists.

**Vertex-z clarification (user, 2026-09-17):** the production requests FIVE vertex positions, x = -0.7, y = -0.6 mm,
z = -71.2 / -28.5 / -4.8 / 18.9 / 61.6 mm, sampled from μ = -4.8 mm, σ = 47.4 mm; only the z = -71.2 mm slice (r17864) is
finished. So "−71.2 mm ≠ data" is not a mismatch of the design, just one of the five points. The -5.8 ± 47 mm quoted in R2 is the
per-event RECONSTRUCTED primary vertex (`vtx_z[0]`) of the pbpb24 DATA skim, first 20k events of `data_pbpb24_part1.root` — a
distribution, not the conditions-DB beam spot; per-run values over the full dataset: see the table below (R3) once filled.

| quantity (Tight WP) | r17864 (2024 cond.) | r17618 (2023 cond.) | ratio |
|---|---|---|---|
| single-μ reco ε, FCal 0-5 %, bins 4.5-8 / 8-14 / 14-25 / 25-100 GeV | 0.563 / 0.624 / 0.660 / 0.640 (N = 252/213/153/139) | 0.588 / 0.661 / 0.675 / 0.709 | 0.96 / 0.94 / 0.98 / 0.90 (±0.05-0.06) |
| single-μ reco ε, FCal 5-10 % | 0.54 / 0.55 / 0.56 / 0.63 | 0.65 / 0.70 / 0.73 / 0.75 | 0.83 / 0.79 / 0.76 / 0.83 (> 3σ each; r17864 "5-10 %" is not a comparable population) |
| single-μ reco ε, **centrality-blind** (all events, same trees) | 10681 / 18304 = **0.584** | 12058 / 18284 = **0.660** | **0.885** |
| mu4 ε, Step-1 selection, 0-5 % | 314 / 424 = 0.741 | 6882 / 8902 = 0.773 | 0.958 ± 0.03 |
| mu4 ε, coarse bins 4.5-8 / 8-14 / 14-25 / 25-100 | 0.64 / 0.85 / 0.71 / 0.78 | 0.70 / 0.80 / 0.80 / 0.81 | 0.92 / 1.07 / 0.89 / 0.96 (±0.04-0.06) |
| mu4 ε, centrality-blind | 7388 / 9710 = 0.761 | 8654 / 11084 = 0.781 | 0.975 |
| (pT_reco − pT_truth)/pT_truth, 0-5 %: mean / RMS | −0.25 % / **3.24 %** | −0.13 % / 2.99 % | RMS +8 % |
| η_reco − η_truth RMS | 0.74e-3 | 0.74e-3 | = |
| from_same_b, all OS reco pairs | 0.3135 ± 0.0057 (6689) | 0.3200 ± 0.0057 (6593) | consistent |
| from_same_b, Tight OS pairs (pT>4.5, |η|<2.4 both) | 0.2947 (2389 pairs) | 0.2982 (2998 pairs) | consistent; pair count −20 % ≈ 0.885² |

Truth structure is unaffected by the r-tag (from_same_b, barcode layout, truth-muon multiplicity 1.88 vs 1.84 per event all
agree); the r17864 deficits are all reconstruction-level: −11 % Tight reco efficiency (centrality-blind), −2..−4 % mu4
efficiency, +8 % pT-resolution width, on top of the broken calorimeter (R1). Nominal pipeline plots for r17864 exist but
their per-centrality panels are statistics-starved (only 4.2 % of events in "0-5 %").

**Vertex z (user question):** the real Pb+Pb 2024 PV is at z = −5.8 mm (RMS 47 mm; pbpb23 −5.7 ± 43, pbpb25 −2.0 ± 49,
pp24 −5.1 ± 52 mm, from the data skims), so the pinned −71.2 mm of s4684/r17864 does NOT reproduce the 2024 beam spot
(the pbpb23-conditions sample's −3.3 mm was close to reality). A 7 cm shift changes the FCal η-acceptance by Δη ≈ 0.015
(percent-level ΣE_T change, A/C asymmetric) and the calorimeter time-of-flight by 0.23 ns (≪ 25 ns) — it cannot produce
negative sums, a 10× drop, or the observed A = C symmetry; the vertex is also reconstructed correctly (residual −5 µm).
It is therefore not the cause of R1; it could at most shift η-edge acceptances by ~0.03-0.07 and should be corrected to the
data value (with the data's ~47 mm spread) in the full production regardless.

### R3. Pb+Pb 2024 data primary-vertex position, per run (2026-09-17; for the vertex-sampling μ/σ discussion)

Source: `~/usatlasdata/dimuon_data/pbpb_2024/data_pbpb24_part*.root`, `HeavyIonD3PD`, `vtx_{x,y,z}[0]` with `vtx_ntrk[0] >= 10`
(reconstructed PV of HardProbes-triggered events; skim keeps |z| < 250 mm). NOT the conditions-DB beam spot (`/Indet/Beampos`).

| run | N | ⟨x⟩ | ⟨y⟩ | ⟨z⟩ [mm] | RMS z [mm] |
|---|---|---|---|---|---|
| 489703 | 8 977 701 | −0.732 | −0.648 | −5.6 | 46.8 |
| 489718 | 3 499 471 | −0.732 | −0.643 | −5.5 | 48.1 |
| 489749 | 5 412 686 | −0.734 | −0.644 | −5.4 | 47.8 |
| 489764 | 8 761 704 | −0.731 | −0.646 | −5.6 | 46.8 |
| 489801 | 7 038 785 | −0.731 | −0.648 | −5.6 | 48.3 |
| 489895 | 1 904 486 | −0.724 | −0.646 | −6.2 | 48.5 |
| 489909 | 8 237 955 | −0.729 | −0.649 | −7.6 | 47.7 |
| 489938 | 1 002 470 | −0.739 | −0.648 | −6.9 | 48.2 |
| 489961 | 7 293 753 | −0.735 | −0.649 | −6.3 | 47.7 |
| 490085 | 9 413 078 | −0.737 | −0.655 | −5.7 | 47.8 |
| 490145 | 9 466 199 | −0.738 | −0.648 | −5.7 | 47.8 |
| 490156 | 9 192 522 | −0.734 | −0.651 | −4.6 | 48.1 |
| 490182 | 6 667 419 | −0.735 | −0.652 | −4.3 | 48.7 |
| 490223 | 5 629 562 | −0.735 | −0.648 | −4.6 | 48.8 |
| **all** | 92 497 791 | −0.73 | −0.65 | **−5.61** | **47.8** |

vs the production request (x −0.7, y −0.6, μ_z −4.8, σ_z 47.4 mm): x, y within 0.05 mm; μ_z differs by 0.8 mm (inside the
run-to-run range −4.3 … −7.6), σ_z within 1 %. No change to μ/σ is warranted on this basis (0.8 mm = 0.02 σ, Δη ~ 1e-4).

### R4. ID track multiplicity and L1 total E_T, r17864 vs r17618 (pTH125_300, no centrality cut, 2026-09-17; /review-plot PASS)

`plots/r17864_rtag_sanity/trk_multiplicity_r17864_vs_r17618_pTH125_300.png` (4 panels, `trk_numqual[0,3,4,7]`) and
`l1te_r17864_vs_r17618_pTH125_300.png` (`plot_r17864_event_level.cxx` modes `ntrk` / `l1te`; raw skim counters, unweighted —
`EventWeights[0]` = 1 in every event of both single-slice samples).

| quantity (per event, 10 000 ev each) | r17864 (2024 cond.) | r17618 (2023 cond.) | ratio |
|---|---|---|---|
| N_trk, p_T > 400 MeV (`trk_numqual[0]`) | 5830 ± 1740 (RMS) | 5435 ± 1460 | 1.073 |
| N_trk, HITight (`trk_numqual[3]`) | 2183 ± 359 | 2400 ± 362 | 0.910 |
| L1 ΣE_T median (mean) [TeV] | **3.21 (3.47)** | 16.14 (16.26) | 0.20 |

- **Tracker is essentially unaffected** — the HIJING event is fully there (+7 % untagged tracks, −9 % HITight tracks). The
  2024 sample has a high-multiplicity shoulder at 8.5–11.5 k untagged tracks that is ABSENT in HITight (2024 HITight max 3540
  < 2023 max 3718): the 2024-conditions reconstruction produces extra LOW-QUALITY tracks and slightly fewer tight ones —
  consistent with the −11 % Tight muon reco efficiency of R2 (ID-track quality), a second, milder symptom next to the
  calorimeter.
- **L1Calo sees the same deficit as the offline calorimeter**: L1TE ×5 low (3.2 vs 16.1 TeV; no negative L1 sums since the
  L1 ET sum is unsigned). So the problem is upstream of both the offline cell energies and the L1Calo trigger-tower path —
  digitisation-level (consistent with the R1 timing/BCID hypothesis), not offline reconstruction.
- Side notes: one r17618 event has L1TE = 32.767 TeV (16-bit saturation of the L1 sum); the r17618 L1TE distribution is
  double-humped (15 and 19 TeV) like its FCal ×1.5 step (R1b) — a feature of that production's HITS composition.

### NTUP naming change (2026-09-17, user decision; /review-analysis-code)

Overlay NTUPs renamed `…FullSimHIJINGOverlayPP24.NTUP.root` → `…FullSimHIJINGOverlayPbPb<yy>.vtxz<z>mm_b<lo>_<hi>fm.NTUP.root`
("PP24" was wrong; the tag now encodes the production configuration — Pb+Pb 2024 requests 5 vertex-z points × 4 b intervals
on top of 4 beams × 6 slices). On disk: pbpb24 `FullSimHIJINGOverlayPbPb24.vtxz-71_2mm_b0_5fm` (1 file), pbpb23
`FullSimHIJINGOverlayPbPb23.vtxz-3_3mm_b0_5fm` (6 files) + `_r17662` variant; diagnostic run-dir symlinks re-pointed;
`.bak_20260709` LGD symlinks and grid dataset names keep the legacy tag. Code: `FullSimSampleType.h`
(`FullSimOverlayConfigTag` — ONE config per year, never a glob; `FullSimSampleFileTag(t, year)`), `PythiaAlgCoreT.c`
(year-aware tag + new zero-input guard: a dir/year mismatch now throws instead of exiting 0 with empty outputs under
`allow_missing_slices`), `grid_monitor.sh` (tag delimited by the `<Month><Year>.v<n>` token; `-` allowed), the three overlay
`grid_sub*.sh` (`_cfg_tag <z> <ip>`, one `_submit` per (beam, slice, z, b)), `plot_r17864_event_level.cxx` (paths from the
header). Verified: pbpb24 + pbpb23-slice single-muon NTP reruns reproduce N = 10000, w_factor 2.09733e-05 / 2.09729e-05 and
byte-identical output sizes. Rule + skim side documented in `SkimCode/README.md` "Overlay NTUP naming",
`docs/pythia_fullsim_overlay.md`.

**Observed during the review (not this session's change):** `~/usatlasdata/pythia_truth_full_sample/pythia_5p36TeV/ami_info/`
was renamed to `ami_info_PDF/` (empty) + `ami_info_nPDF/` at 19:17 by another session. `PythiaAlgCoreT.c:41` falls back to that
`ami_info/` for the pbpb23 overlay and pp24 TEST samples → their NTP now throws "missing AMI file" until that session finishes.

### AMI-weight correction (2026-09-17, USER RULING — the physics rule, recorded here for good)

**The AMI weight of every event is the weight of its PYTHIA evgen dataset, not of anything in the
AOD chain.** `e8613` in `e8613_e8586_s4684_r17864_r17855` is the evgen used for **HIJING**
(`860250.Hijing_PbPb_UCC_Flow_JJFV6_ip0_5.e8613_...`), NOT a Pythia truth evgen. The Pythia evgen
e-tag is not shown in the overlay's AMI tag chain because the Pythia events were already generated
and are used as INPUT to the simulation, not re-run. For each overlay event we use the AMI weight
from Pythia truth, never from HIJING — and that is exactly why, in the HIJING-overlay samples, only
truth muons from Pythia enter the reconstruction- and trigger-efficiency evaluation, never truth
muons from HIJING: HIJING is background / environment / underlying event for the Pythia truth
particles; a HIJING truth muon would carry the wrong (Pythia) weight.

What was wrong (session 2, commit 2d78e97, 2026-09-16): `PythiaAlgCoreT` was changed to prefer "a
sample's own `ami_info/`" whenever the directory existed, and this sample's `ami_info/` held the AOD
dataset's AMI numbers (registered as `ami_weights.md` B2, Δ 1e-4 vs table A — numerically harmless
here, wrong in principle, and a trap for the next production). Undone 2026-09-17, in the same
commit series:

- `FullSimSampleType.h`: `enum PythiaEvgen { nPDF, PDF }` + `FullSimSampleEvgen(type, isTestSample)`
  (pp24 FULL → PDF, everything else → nPDF), `PythiaEvgenAmiDir/FileName/Dsids`. AMI files live
  with the evgen in ONE place each: `pythia_truth_full_sample/pythia_5p36TeV/ami_info_nPDF/`
  (e8599, 802758–802781, 4 isospins 4:6:6:9, nuclear PDF nNNPDF30 — Pb+Pb-suitable; truth-only
  analysis, pp24 TEST, every HIJING overlay incl. this one) and `ami_info_PDF/` (`_pdf`,
  803015–803020, pp only, proton PDF NNPDF23LO — pp-suitable; pp24 fullsim FULL). Both dirs +
  `fetch_ami_info_{nPDF,PDF}.sh` were laid out by the user on 2026-09-17.
- `PythiaAlgCoreT.c` InitParams: `ami_evgen` derived from the sample switch; the own-`ami_info/`
  probe deleted; `expected_ami_dsids` defaults to the evgen's DSIDs. InitInputFullsim: file name
  from `PythiaEvgenAmiFileName` (the PDF files carry the `_pdf` suffix). **InitInputCentrProd
  (truth-only reader): was silently skipping** (WARNING + `ami_weight = 0`) since the user's
  rename of `ami_info/` → `ami_info_nPDF/` at 19:17; now reads `ami_info_nPDF/`, missing file
  fatal, DSID checked against 802758–802781.
- The four `run_pythia_fullsim_*_full_sample.sh` no longer set `ami_info_dir_override`;
  `pipeline_pythia_fullsim_pp.sh` preflight and `mc_pthat_slice_mass_statistics.cxx` read
  `ami_info_PDF/`. On disk: `pythia_fullsim_full_sample/ami_info/` →
  `ami_info_AOD_record_20260713/`, `pythia_fullsim_hijing_overlay_test_sample/ami_info/` →
  `ami_info_AOD_record_20260915/` (renames, read by nothing).
- **pp24 fullsim FULL sample: no rerun needed.** Its `ami_info/` copy (fetched 2026-07-13 on the
  `_pdf` AOD) had σ and ε_filt byte-identical to the evgen's `ami_info_PDF/` files; the nominal
  `muon_pairs_pythia_fullsim_pp24_no_data_resonance_cuts_full.root` (2026-09-10) carries
  kin0 weight 1.81001e-5 nb = 36.7432/2029996 and kin5 6.1807e-7 = 0.19779/319999 — the PDF
  evgen's numbers. Every pp24 fullsim result, including `plots/mc_data_compr/` (which reads
  `histograms_pythia_fullsim_pp24_*_full.root`, NOT the truth sample, since 2026-08-25), was
  already built on the pp-suitable evgen with its own weights.
- Docs: `docs/ami_weights.md` rewritten around the evgen rule (B2 retired), `docs/pythia_truth.md`,
  `docs/pythia_fullsim_pp.md`, root `README.md`, `pipelines/fullsim_sample_layout.sh`, SkimCode
  comments (`grid_sub_pbpb24_test_sample.sh`, `fullsim_pp24_full_to_lgd.sh`).

### R5. Where the Pythia TRUTH-only (nPDF evgen) quantities are consumed — survey 2026-09-17 (two independent subagents, cross-checked; awaiting user decisions)

Facts established: (a) the pp24 fullsim NTUP is a strict superset of the truth-only TRUTH0 NTUP — all 53
truth branches (`truth_id/barcode/status/pt/eta/phi/m/parents`, `truth_mupair_*`, `Q`, …) are present, all
events stored — and `PythiaFullSimAnalysis` composes the SAME `PythiaTruthExtras` ancestry tracing
(`from_same_b`, flavor/origin categories, last-b-hadron vector) into truth-seeded pair trees
`muon_pair_tree_sign{1,2}` / `_kin<k>_sign<s>` of `muon_pairs_pythia_fullsim_pp24_..._full.root`; (b) the
PDF (pp) evgen has the SAME event count per slice as the 4-beam nPDF total in slices 1-4 but 3.7x / 6.75x
fewer in pTH70_125 / pTH125_300 (320 k vs 1.195 M / 2.16 M), and the fullsim FULL sample already IS the
whole PDF evgen (2 029 996 / 2 039 700 …) → a truth-only skim of the PDF evgen would add NO events;
(c) the truth-only RDF fills (`RDFBasedHistFillingPythiaTruth.cxx`: flavor/origin-binned kinematics,
`_sigsel` templates, `h2d_sig_accept_*`) have no counterpart in `RDFBasedHistFillingPythiaFullsim.cxx`
(which has `pass_signal_truth` numerators but no `from_same_b`-only denominator, no category split).

| # | result | code / output | truth role | today | pp-suitable needed? |
|---|---|---|---|---|---|
| 1 | **pp24 MC-vs-data comparison** `plots/mc_data_compr/{signal,generic}` | `plotting_codes/mc_data_compr/*`, `McDataComprConfig.h` → `histograms_pythia_fullsim_pp24_*_full.root` | truth-quantity hists (`*_pass_signal_truth`) from the FULLSIM truth block | **already pp24 fullsim + PDF weights since 2026-08-25** (truth sample rejected then for its 4:6:6:9 isospin average, `mc_data_compr_signal_generic_split.md`) | **no remake needed** |
| 2 | signal acceptance ε_acc(pair pT, pair η) `plots/single_b_analysis/pythia/pythia_sig_accept_*` | `plot_signal_acceptance_pythia.cxx` ← `h2d_sig_accept_{num,denom}_pt_eta` in `histograms_pythia_5p36TeV_*.root` | acceptance of the pp signal region on truth single-b pairs; the α reported with the pp crossx (not applied in crossx code) | nPDF truth sample | YES if reported for pp — needs the denominator fill added to the fullsim RDF (or the truth pipeline run on the fullsim NTUP) |
| 3 | acceptance cutflow >60 GeV, ΔR-vs-pair-pT diagnostic, single-b pT / last-b-hadron ratio | `plot_sig_accept_cutflow_above_60GeV.cxx`, `plot_dr_vs_pair_pt_diagnostic.cxx`, `plot_single_b_pt_ratio_to_last_b_hadron.cxx` (RDF on `muon_pair_tree_sign2`) | settled pp-signal diagnostics (closed docs) | nPDF truth sample | path edit only (fullsim OS tree has the same branches); low priority |
| 4 | flavor / origin categorized generator-level plot sets `pythia_5p36TeV/plots/{flavor,origin}_categoried/` | `plot_{flavor,origin}_categorized_kinematics.cxx` ← truth RDF; labelled nNNPDF30 | generator-level composition studies | nPDF truth sample | USER DECISION: pp-reference composition → PDF (needs a fullsim-NTUP input mode of the truth pipeline or category fills in the fullsim RDF); Pb+Pb-signal composition → nPDF is right |
| 5 | low-mass template-fit truth `_sigsel` templates; k = G_SS/G_OS studies (data-area macros `template_fitting/*/code/*.C`) | `FillHistogramsTemplateMinvSignalRegion` | S/G shapes, k ratio | nPDF truth; superseded for the nominal reco-level fit (`tf_nominal_fit_build.md`) | no (k is a shape ratio; nominal templates come from fullsim reco) |
| 6 | 5.36 vs 5.02 TeV truth comparison | `plot_pythia_compr_priv_nonpriv_5_02TeV.cxx` | √s dependence, both central productions nuclear-PDF | nPDF | no (like-for-like) |
| 7 | HIJING overlay reco-eff / det-response / MC trig-eff / closure; overlay truth reference | fullsim NTUP truth block of each overlay | Pb+Pb | nPDF (= the overlay's own evgen) | no — nPDF is correct |
| 8 | crossx_summary / kn-slice tables from the fullsim runs | `PythiaTruthExtras.c` diagnostics written under the truth dir | bookkeeping | evgen of each run | no |

Nothing uses the truth sample as a same-event reference for the pp24 fullsim (all efficiency denominators
come from the fullsim NTUP's own truth block), and no truth-sample figure is in the note yet.

Also: the truth-only `_pp_only` variant (nPDF, pp beam only) has no consumer.

**Answer to "skim the pp-suitable evgen?"** — not needed: the fullsim NTUP already carries the complete
PDF-evgen truth block with identical branch content, and the ntuple processing already produces the truth
pair trees from it. What is missing is only on the RDF/plot side (item 2 denominator, item 4 category
fills) or a `fullsim_input_dir`-style input mode for `PythiaTruthAnalysis`. The statistics argument goes the
other way: the nPDF truth sample has 3.7-6.75x more events in the two highest slices (the same deficit
`pthat_slice_mass_stats_sample_request.md` is addressing).

### R6. Provenance of the REPRODUCED slice (2026-09-22): what actually differs from the broken production

New dataset (user-provided, now the only content of the merge container):
`…merge.AOD.e8613_e8586_s4684_r17864_r17855_tid52629267_00`, 10 files / 10 000 events / 54.891 GB.

| | broken production | reproduction |
|---|---|---|
| merge task / dataset | 52480369 (`tid52480369_00`) | **52629267** (`tid52629267_00`) |
| recon ("pile") task | 52480366, created 2026-09-10, done 09-11 | **52629263**, created 2026-09-18, done 09-21 |
| production request | 71862 | **71862 — the same request** |
| r-tag / s-tag / e-tag | e8613_e8586_s4684_r17864_r17855 | **identical** |
| release / transform | Atlas-24.0.90, `Reco_tf.py`, x86_64-el9-gcc13-opt | **identical** |
| input files / events | 1000 / 10 000, 0 failed | **identical** |
| AMI `ami show ami-tag r17864` | saved 2026-09-15 as `r17864.txt` | re-fetched 2026-09-22 as `r17864_refetched_20260922.txt` — **byte-identical (`diff` empty)** |
| AMI dataset record (AOD) | σ 89.541 nb, genFiltEff 2.342314e-3 | **identical**; the AMI container now lists BOTH tids (`prodsysIdentifier_0` = 52480369, `_1` = 52629267) |

The old `tid52480369_00` dataset has been **unregistered from Rucio** (`rucio stat` finds nothing; the
recon container likewise holds only `…recon.AOD…_tid52629263_00`). Rucio replicas of the new one:
UKI-LT2-QMUL 10/10, CERN-PROD 1/10, RAL-LCG2-ECHO 1/10 — **no BNL replica**.

**Nothing visible in the bookkeeping changed.** The user flagged this before the work started: the
`preExec` list of r17864 in AMI is the same as before the reproduction (and still differs from the
r17618 list used for the pbpb23 test samples). PanDA does not expose job-level `jobParameters`
through its JSON API without authentication, so a task-level `preExec`/`preInclude` override applied
at submission — the one thing that would NOT show in the AMI tag — can be neither confirmed nor
excluded from the outside. **The data therefore decide**, which is why the verification below re-runs
the identical check chain rather than trusting the tag.

Archived, so the two productions never get confused (the NTUP file name is configuration-keyed and
therefore IDENTICAL between them): `…/pythia_fullsim_hijing_overlay_test_sample/broken_calo_tid52480369/`
holds the broken production's `muon_pairs_*`, `histograms_*`, `hists_*`, its three plot dirs and a copy
of the old `plots/r17864_rtag_sanity/` (verified byte-identical before the move), plus a README naming
the dataset and the R1/R4 evidence. The broken NTUP itself stays in the parent dir until
`grid_monitor.sh` renames it to `*.bak_<date>` at download time.

Skim resubmitted with the new input: `grid_sub_pbpb24_test_sample.sh` VER_TAG **Sep2026.v3**,
**jediTaskID 52687659** (2026-09-22 15:2x), outDS
`user.yuhang.NTUP.Pythia_5p36TeV_pp_hQCD_DiMu_pTH125_300.FullSimHIJINGOverlayPbPb24.vtxz-71_2mm_b0_5fm.Sep2026.v3.`

### R7. FIRST VERDICT on the reproduction: the calorimeter is STILL BROKEN (100-event local test, 2026-09-22)

Ahead of the grid skim, one AOD file of the new dataset (`AOD.52629267._000020.pool.root.1`, 5.44 GB,
downloaded from UKI-LT2-QMUL to `…/test_aod_tid52629267/`) was skimmed locally with the standard
100-event test (`run_test.sh r17864_tid52629267 orig <local file>`, `TRIGRATES_RUNMODE=ppmcfullsim_hioverlay24`,
athena exit 0). Raw per-event skim branches only (`FCal_Et`, `L1TE`, `centrality`, `trk_numqual`) —
same provenance justification as `plot_r17864_event_level.cxx`.

| | FCal SE_T q16/50/84 [TeV] | negative | L1TE median | centrality undefined | 0-10 % | N_trk | HITight |
|---|---|---|---|---|---|---|---|
| **NEW pbpb24 tid52629267** (100 ev) | **-0.81 / 0.28 / 1.94** | **44 %** | **2.86 TeV** | **45 %** | 7 % | 5256 | 2068 |
| OLD pbpb24 tid52480369 (100 ev) | -0.63 / 0.53 / 1.78 | 34 % | 3.09 TeV | 36 % | 4 % | 5895 | 2194 |
| REF pbpb23 r17618 (10 000 ev) | 3.79 / **4.40** / 5.11 | **0 %** | **16.14 TeV** | **0 %** | 99.3 % | 5435 | 2400 |

**The reproduction reproduces the failure, not the fix.** Every symptom of R1/R4 is present with the same
magnitude: FCal an order of magnitude low with ~40 % of events reconstructing NEGATIVE energy, L1TE a
factor ~5 low (so the deficit is again upstream of both offline cells and the L1Calo towers, i.e.
digitisation-level), ~45 % of events off the end of the centrality calibration, while the ID track
multiplicity is normal (the HIJING event is fully there). Old vs new differ by ~1.4 sigma on 100 events
(binomial sigma ~5 % on the negative fraction) — statistically indistinguishable.

Combined with R6 (identical AMI tag, identical release, identical task configuration as far as the
outside can see), the conclusion is that the slice was re-run with the SAME configuration, so the
suspect setting — `HITtoRDO: Digitization.PU.BunchStructureConfig='RunDependentSimData.BunchStructureHeavyIon2022'`,
the only stated r17835 -> r17864 change — was never removed. Nothing here points to a NEW failure mode.

Consequences, unchanged from the first production: FCal-based centrality and every calorimeter-dependent
muon quantity are unusable in this r-tag; R2 stays not studiable; the r-tag is NOT acceptable as the
recipe for the full 4-beam x 6-slice x 5-vertex x 4-b production.

**BCID cross-check (negative so far).** Since the suspect setting is a *bunch-structure* option, the
bunch-crossing ID was checked as a possible handle: `bcid` spans 33-3365 (95 distinct in 100 events,
55 % odd, median gap 20) in the new sample vs 11-3413 (97 distinct, 64 % odd) in the old one and
1-3415 (1010 distinct in 10 000 events, median gap 2) in r17618 — the two pbpb24 productions are
indistinguishable from each other, and nothing separates their filled-bunch pattern from the reference
at this statistics. With ~1 event per bcid, FCal-vs-bcid carries no information yet. **Redo on the full
10 000-event NTUP**: (a) the SET of filled bcids, 2024 vs 2023 structure; (b) <FCal> vs bcid, which
would localise an off-peak-sampling explanation.

### R8. FULL-STATISTICS VERDICT: the "reproduction" is a BIT-IDENTICAL re-run (2026-09-23)

Skim `jediTaskID 52687659` (VER_TAG Sep2026.v3) completed 2026-09-23 16:16: 5 grid NTUPs hadd-ed to
10 000 entries, 22.44 GB, validated (pre- and post-merge counts agree). `grid_monitor.sh` renamed the
first production's NTUP to `…vtxz-71_2mm_b0_5fm.NTUP.bak_20260923.root` (22.44 GB, kept), so BOTH
productions are on disk under the same configuration-keyed name.

**The two productions contain the SAME EVENTS.** Matched on `eventNumber` (10 000 of 10 000 common;
the tree ORDER differs because the grid job/hadd order differs):

| quantity, per event | events differing (of 10 000) |
|---|---|
| `FCal_Et` | **0** — bit-identical, max relative difference exactly 0 |
| `L1TE` | **0** |
| `centrality` | **0** |
| `trk_numqual[0]` | 3 |

Summed over the sample: `truth_pt` multiplicity identical to the unit (698 526 429), `vtx_z` sum
identical, N(vtx) identical; reconstructed muons 95 776 vs 95 771 (+5, 5e-5 relative) and
sum(`muon_pt`) 8.314348e9 vs 8.315876e9 (1.8e-4 relative) — job-to-job non-determinism in
track/muon reconstruction, nothing more. **The calorimeter output is literally the same numbers.**

Full-statistics event-level checks, regenerated into `plots/r17864_rtag_sanity/` (the first
production's copies are archived in `broken_calo_tid52480369/plots/r17864_rtag_sanity/`):

| check (10 000 events each) | pbpb24 re-run (tid52629267) | pbpb24 first (tid52480369) | pbpb23 r17618 | verdict |
|---|---|---|---|---|
| FCal SE_T q16/50/84 [TeV] | -0.77 / **0.41** / 2.08 | -0.77 / 0.41 / 2.08 | 3.79 / **4.40** / 5.11 | **FAIL** |
| FCal SE_T < 0 | **36.2 %** | 36.2 % | **0 %** | **FAIL** |
| L1 SE_T median [TeV] | **3.20** | 3.20 | **16.14** | **FAIL** |
| centrality undefined / 0-10 % | **37.4 % / 7.6 %** | 37.4 % / 7.6 % | 0 % / 99.3 % | **FAIL** |
| N_trk (p_T>400 MeV) / HITight | 5830 +- 1740 / 2183 +- 359 | same | 5435 +- 1460 / 2400 +- 362 | tracker OK (+7 % / -9 %) |
| reco PV residual [um] | -5.43, RMS 11.09 | -5.4, RMS 11.1 | -0.2, RMS 9.1 | unchanged |

Files: `fcal_sum_et_…png`, `l1te_…png`, `trk_multiplicity_…png`, `vertex_z_reco_vs_truth_r17864.png`,
`centrality_population_…txt`, `bcid_structure_…png`.

**BCID check.** Both samples fill exactly **1010 of 3564** crossings (identical pattern). Binned 12 bcid
per point, <FCal SE_T> in the 2024 re-run scatters over -0.54 … 5.15 TeV about a mean of 0.64 TeV while
the reference sits at 4.43 TeV — and that scatter is NOT noise: it resolves into a strict structure once
the crossings are ordered by position within the bunch train. **See R9, which identifies the failure
mechanism.**

**Conclusion: the r-tag is NOT fixed and was NOT re-configured.** Bit-identical calorimeter output under
an unchanged AMI tag, unchanged release and unchanged task configuration (R6) means the slice was
re-run with exactly the same recipe. Nothing was learned about a fix; the r17864 recipe remains
unusable for the full 4-beam x 6-slice x 5-vertex x 4-b production.

**The efficiency chain was deliberately NOT re-run.** It is gated on the event-level checks passing
(user ruling 2026-09-17: R2 is not studiable while the FCal is wrong), and it is now also pointless:
with an input NTUP that differs in 3 track counts and 5 muons out of ~96 000, the chain would reproduce
the R2 table digit for digit, far inside its 5 % statistical errors. R2 stands as written.

**Correction to R1's bookkeeping:** R1 records "35.2 % of events negative" for the first production.
Re-measured today on that very file (the `.bak`) with the same macro: **36.2 %**. The quantiles
(-0.77 / 0.41 / 2.07-2.08 TeV) reproduce; the 35.2 % does not and should be read as 36.2 %.

### R9. ROOT CAUSE LOCALISED: the deficit is a periodic function of the position in the bunch train (2026-09-23; /review-plot FAIL -> amended -> re-reviewed)

`plot_r17864_event_level.cxx+("bcid")` and `+("trainpos")`, full 10 000-event samples, raw per-event
skim branches.

**The filled-bunch pattern, measured from the data and IDENTICAL in both samples** (1010 of 3564
crossings; gap spectrum 882 gaps of 2 bcid + 108 of 4): **18 trains of 7 x 8 bunches** — 8 bunches at
2-bcid spacing, then one crossing skipped, so a sub-burst repeats every **18 bcid (450 ns)** and the
4-bcid gaps fall after intra-train positions 7, 15, 23, 31, 39, 47 in every train — **plus 2 ISOLATED
single bunches (bcid 826 and 2611)**. (Train boundary = a bcid DIFFERENCE of more than 5, i.e. more than 4 empty crossings; the
result is unchanged for thresholds 5 / 10 / 20 / 31.)

**The failure is bcid-deterministic, not event-by-event.** Classifying bcids with >=5 events by whether
their events have FCal SE_T > 3 TeV (a b < 5 fm event cannot legitimately be below that):

| | events with FCal > 3 TeV | bcids all-good | all-bad | mixed |
|---|---|---|---|---|
| pbpb24 r17864 | **8.71 %** | 39 | 813 | 121 |
| pbpb23 r17618 | 99.70 % | 947 | **0** | 30 |

**38 of those 39** are exactly the first two bunches of each of the 18 trains plus the two isolated
bunches; the remaining one is bcid 1993 (position 25, 8 events, all in 3.04-4.24 TeV) — a chance pass
of the threshold in the recovery phase, not a second pathology. The classification is robust to the
threshold in shape though not in count: all-good / all-bad = 55/694, 42/757, 39/813, 35/845, 16/904 at
2.0 / 2.5 / 3.0 / 3.5 / 4.0 TeV.

**Only the ISOLATED bunches are healthy — every bunch inside a train is wrong, the head included.**
Splitting the position-0 population (this is the correction that the /review-plot round forced; merging
them had hidden the sample's only clean control):

| position-0 sub-population | N | <FCal SE_T> [TeV] | <L1 SE_T> [TeV] |
|---|---|---|---|
| **2 isolated bunches (no train)** | 23 | **4.564 +- 0.127** | **16.671 +- 0.532** |
| 18 train heads | 184 | **5.104 +- 0.054** (+15 % vs the 2023 reference) | **10.421 +- 0.188** (-36 % vs it) |
| pbpb23 reference, all positions | 10 000 | 4.43 | 15.9-16.7 (all-event mean 16.27) |
| pbpb23 isolated bunches | 14 | 4.349 | 15.902 |

The isolated bunches agree with the 2023 reference in BOTH quantities; the train heads agree in
neither. A crossing with no predecessors is reconstructed correctly; a crossing with predecessors is
not, and already the FIRST bunch of a train is wrong.

**<FCal SE_T> vs position in the train** (`train_position_…png` left panel; right panel L1TE),
in-train bunches only, 2024 re-run / 2023 reference [TeV]:

| pos | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|---|---|---|---|
| 2024 | +5.10 | +4.96 | +1.91 | +0.44 | -0.16 | -0.38 | -0.38 | -0.79 | +0.41 | +2.80 |
| 2023 | +4.44 | +4.51 | +4.47 | +4.46 | +4.53 | +4.40 | +4.44 | +4.53 | +4.41 | +4.42 |

Beyond the head it is **not a monotonic decay but strictly PERIODIC with period 8 filled bunches —
i.e. exactly the 8-bunch sub-burst of the filled pattern above**. Folding positions 8-55 modulo 8
(2024 re-run / 2023 reference, TeV):

| phase | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|---|---|
| **2024 re-run** | +0.371 | **+2.836** | +1.388 | +0.458 | +0.011 | -0.086 | -0.160 | **-0.743** |
| uncertainty | +-0.018 | +-0.021 | +-0.026 | +-0.025 | +-0.024 | +-0.025 | +-0.024 | +-0.022 |
| 2023 r17618 | +4.443 | +4.437 | +4.419 | +4.419 | +4.435 | +4.443 | +4.405 | +4.409 |

The 2023 reference is flat to +-0.02 TeV in every phase; the 2024 sample swings by 3.6 TeV with the
same precision. **The modulation is locked to the sub-burst, not to absolute bcid** — verified two
ways: <FCal> vs absolute `bcid % 16` is featureless (0.56-0.78 TeV), and the mod-8 argmax is phase 1
in all 18 trains although the train starts take 14 distinct values of `start % 16`.

**Interpretation (hypothesis for the production contact, NOT established here).** Everything fits a
mismatch between the bunch pattern assumed by the digitisation and the one the events were actually
placed in: a crossing with no predecessors comes out right, a crossing with predecessors sits on the
accumulated undershoot of their bipolar LAr pulses and is dragged low and then negative, recovering
once per sub-burst. That is the failure mode the ONE stated r17835 -> r17864 change would produce if
its pattern does not match this sample:
`HITtoRDO: Digitization.PU.BunchStructureConfig='RunDependentSimData.BunchStructureHeavyIon2022'`
("bunch structure added"). Caveat, and it is a real one: a pure undershoot-accumulation picture does
not explain the train head being 15 % HIGH in FCal while its L1TE is 36 % LOW — there is a
baseline/pedestal component too. Confirming or refuting this needs the production side (an r17835-style
clone without the setting, or the digitisation configuration actually used); it cannot be settled from
the NTUP.

**This is the actionable content of the report to the production contact** — a specific, reproducible,
periodic signature with a named suspect setting, not "the calorimeter looks wrong".

### R10. r17913 — the FIRST genuinely re-configured production (2026-09-29)

AMI `ami show ami-tag r17913` (saved as `r17913.txt`): **"clone r17864 with additional preExec
parameters"**, created 2026-09-25 by the same production contact, release Atlas-24.0.90,
`Reco_tf.py`. Everything r17864 had is unchanged — `postExec` beam spot (−0.7, −0.6, −71.2) mm with
σ = 0.01 mm, `preInclude Campaigns.MC23HeavyIons2024NoPileUp` + `HIRecConfig.HIModeFlags.HImode`,
`inputCavernHitsFile mc23_5p36TeV.860250.Hijing_PbPb_UCC_Flow_JJFV6_ip0_5` (b 0–5 fm),
`conditionsTag OFLCOND-MC23-SDR-RUN3-09`, `geometryVersion ATLAS-R3S-2021-03-02-00`,
steering `doRDO_TRIG`+`doTRIGtoALL`, merge tag r17855 — so the configuration tag
`vtxz-71_2mm_b0_5fm` and hence the NTUP file name are UNCHANGED.

**The whole difference is the `HITtoRDO` preExec:**

| | r17864 (broken) | r17913 |
|---|---|---|
| `Digitization.PU.BunchStructureConfig` | `'RunDependentSimData.BunchStructureHeavyIon2022'` | same |
| `Digitization.PileUp` | — | `True` |
| `Digitization.DoXingByXingPileUp` | — | `True` |
| `Digitization.PU.InitialBunchCrossing` | — | **`0`** |
| `Digitization.PU.FinalBunchCrossing` | — | **`0`** |
| `LAr.ROD.nSamples` | — | **`4`** |
| `LAr.ROD.NumberOfCollisions` | — | **`0`** |
| `Digitization.HighGainFCal` | — | **`True`** |

This is a targeted answer to R9. `Initial/FinalBunchCrossing = 0` collapses the out-of-time pile-up
window to the in-time crossing alone, which removes exactly the mechanism R9 localised — every
crossing WITH predecessors reconstructed low and then negative, only the 2 isolated bunches were
healthy, and the deficit was periodic with the 8-bunch sub-burst. `LAr.ROD.NumberOfCollisions = 0`
sets the LAr optimal-filtering coefficients / pedestal for zero pile-up rather than for an assumed
collision rate (the baseline/pedestal component R9 flagged as not explained by undershoot alone),
`nSamples = 4` matches the 2024 HI 4-sample readout, and `HighGainFCal = True` the FCal gain used in
HI running. **Prediction if the hypothesis is right: FCal ΣE_T back to ≈ 4.4 TeV with 0 % negative
events, L1 ΣE_T back to ≈ 16 TeV, centrality 0–10 % for ~all events, and the train-position
modulation flat.**

Dataset (user-provided; `rucio stat`, 2026-09-29):
`mc23_5p36TeV:mc23_5p36TeV.802776.Py8EG_A14_pp_hQCD_DiMu_pTH125_300.merge.AOD.e8613_e8586_s4684_r17913_r17855_tid52733789_00`
— 10 files / **10 000 events** / **50.572 GB** (r17864 tid52629267: 54.891 GB — the AODs genuinely
differ in size, unlike the bit-identical re-run). Replicas: **DESY-ZN_DATADISK 10/10**,
RAL-LCG2-ECHO 1/10; **no BNL replica** (as for r17864). Not yet in the AMI *dataset* catalogue
(`ami show dataset info` → not found), which is bookkeeping only: the AMI **weight** is the Pythia
evgen's (802776 e8599, `ami_info_nPDF/`, registry table A), never the AOD chain's — see the
2026-09-17 ruling above and `docs/ami_weights.md`.

**Archiving, so three productions cannot share one file name.** The NTUP name encodes the
configuration, not the r-tag (2026-09-17 user decision), so before the r17913 download:

| moved to | file |
|---|---|
| `broken_calo_tid52480369/` (joins the first production's derived products) | the first production's NTUP, was `…NTUP.bak_20260923.root` in the parent dir |
| `broken_calo_r17864_rerun_tid52629267/` (new, with README) | the re-run's NTUP, was the canonical `…NTUP.root` |

Both under the canonical base name inside their archive dir; nothing globs the sample dir
(`PythiaAlgCoreT`/`FullSimSampleFileTag` build the exact file name), so the archived copies are
invisible to the chain. `plots/r17864_rtag_sanity/` is left in place — it is r-tag-keyed and does not
collide with the new `plots/r17913_rtag_sanity/`. Nothing was deleted.

Skim: `grid_sub_pbpb24_test_sample.sh` VER_TAG **Sep2026.v4**, **jediTaskID 52785715**, outDS
`user.yuhang.NTUP.Pythia_5p36TeV_pp_hQCD_DiMu_pTH125_300.FullSimHIJINGOverlayPbPb24.vtxz-71_2mm_b0_5fm.Sep2026.v4.`
No skim-code change needed (r17913 is an r17864 clone → same 2024 HI menu, same RUNMODE
`ppmcfullsim_hioverlay24`); the branch check on the downloaded NTUP confirms it.

### R11. FIRST VERDICT on r17913: the calorimeter is FIXED (100-event local test, 2026-09-29)

Ahead of the grid skim, one AOD file of the new dataset (`AOD.52733789._000008.pool.root.1`,
5.06 GB, downloaded from DESY-ZN_DATADISK to `…/test_aod_tid52733789/`) was skimmed locally with the
standard 100-event test (`run_test.sh r17913 orig <local file>`,
`TRIGRATES_RUNMODE=ppmcfullsim_hioverlay24`, athena exit 0, zero ERROR/FATAL lines). Raw per-event
skim branches only (`FCal_Et` [MeV], `L1TE` [**GeV**], `centrality`, `trk_numqual`) — same
provenance justification as `plot_overlay_rtag_event_level.cxx`.

| | FCal ΣE_T q16/50/84 [TeV] | negative | L1 ΣE_T median (mean) [TeV] | cen. undefined | 0–10 % | N_trk | HITight |
|---|---|---|---|---|---|---|---|
| **NEW r17913** (100 ev) | **3.94 / 4.44 / 5.12** | **0.0 %** | **16.22 (16.69)** | **0 %** | **100 %** | 5721 | **2362** |
| OLD r17864 tid52629267 (10 000 ev) | −0.77 / 0.41 / 2.07 | 36.2 % | 3.21 (3.47) | 37.4 % | 7.6 % | 5830 | 2183 |
| REF pbpb23 r17618 (10 000 ev) | 3.79 / **4.40** / 5.11 | **0 %** | **16.14 (16.26)** | 0 % | 99.3 % | 5435 | 2400 |

**Every symptom of R1 / R4 is gone**, and each quantity lands on the 2023-conditions reference rather
than merely moving in the right direction: FCal median 4.44 vs 4.40 TeV with ZERO negative events,
L1 ΣE_T 16.22 vs 16.14 TeV, and the centrality calibration now puts 100 % of a b < 5 fm sample in
0–10 % (reference 99.3 %) with nothing off the end. So the digitisation fix works on BOTH the offline
cell energies and the L1Calo tower path, as it must if the cause was upstream of both (R4).

**The milder second symptom recovered too.** The HITight track count — R4's "extra low-quality, fewer
tight tracks", the ID-side companion of R2's −11 % Tight muon reco efficiency — goes 2183 (−9.0 % vs
reference) → **2362 (−1.6 %)**, and the untagged count 5830 (+7.3 %) → 5721 (+5.3 %). On 100 events
these carry ~±3 % statistical uncertainty each, so this is indicative, not yet a measurement; the
10 000-event NTUP settles it (R12) and makes the R2 efficiency comparison studiable again.

Caveat on the statistics: 100 events give a binomial σ of ~5 % on a negative-event fraction, so
"0.0 %" excludes the r17864 value (36.2 %) overwhelmingly but only bounds a residual pathology at the
few-percent level; the quantiles and the medians are already precise enough to be conclusive.

### R12. FULL-STATISTICS VERDICT on r17913: the calorimeter is repaired and the R9 signature is GONE (2026-09-29)

Skim `jediTaskID 52785715` (VER_TAG Sep2026.v4) completed 2026-09-29 01:40: 5 grid NTUPs hadd-ed to
**10 000 entries, 21.66 GB** (r17864: 22.44 GB), pre- and post-merge counts agree, `merging-record.txt`
updated. All seven event-level modes of `plot_overlay_rtag_event_level.cxx` regenerated into
**`plots/r17913_rtag_sanity/`** (the r17864 set is untouched in `plots/r17864_rtag_sanity/`).

| check (10 000 events each) | **pbpb24 r17913** | pbpb24 r17864 | pbpb23 r17618 | verdict |
|---|---|---|---|---|
| FCal ΣE_T q16/50/84 [TeV] | **3.74 / 4.29 / 5.04** | −0.77 / 0.41 / 2.07 | 3.79 / **4.40** / 5.11 | **PASS** |
| FCal ΣE_T < 0 | **0.0 %** | 36.2 % | 0 % | **PASS** |
| L1 ΣE_T median (mean) [TeV] | **15.70 (15.94)** | 3.20 (3.47) | 16.14 (16.26) | **PASS** |
| centrality undefined | **0.00 %** | 37.4 % | 0.00 % | **PASS** |
| centrality 0–10 % / 10–20 % | **99.32 % / 0.68 %** | 7.6 % / 7.9 % | 99.32 % / 0.68 % | **PASS — identical to the reference** |
| (0–5 % / 5–10 % sub-split) | 78.23 / 21.09 % | 4.20 / 3.39 % | 81.54 / 17.78 % | usable bins again |
| reco PV residual [µm] | −4.33, RMS 9.58 | −5.43, RMS 11.09 | −0.22, RMS 9.10 | improved |
| N_trk (p_T > 400 MeV) | 5285 ± 1520 | 5830 ± 1740 | 5435 ± 1460 | **−2.8 %** (was +7.3 %) |
| N_trk HITight | **2257 ± 345** | 2183 ± 359 | 2400 ± 362 | **−6.0 %** (was −9.0 %) — see below |
| filled bcids (of 3564) | 1010 | 1010 | 1010 | identical pattern |

**R9 is gone — this is the decisive check.** `trainpos`, in-train bunches only, ⟨FCal ΣE_T⟩ [TeV],
positions 0–9: r17913 **4.316 / 4.444 / 4.310 / 4.368 / 4.355 / 4.362 / 4.432 / 4.391 / 4.344 / 4.341**
vs r17618 4.444 / 4.508 / 4.465 / 4.458 / 4.526 / 4.398 / 4.439 / 4.527 / 4.406 / 4.417. Folding
positions 8–55 modulo 8 (the sub-burst period that r17864 was locked to):

| phase | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|---|---|
| **r17913** | +4.385 | +4.363 | +4.346 | +4.343 | +4.362 | +4.369 | +4.370 | +4.417 |
| uncertainty | ±0.020 | ±0.019 | ±0.019 | ±0.019 | ±0.019 | ±0.019 | ±0.019 | ±0.019 |
| r17864 | +0.371 | **+2.836** | +1.388 | +0.458 | +0.011 | −0.086 | −0.160 | **−0.743** |
| r17618 | +4.443 | +4.437 | +4.419 | +4.419 | +4.435 | +4.443 | +4.405 | +4.409 |

r17864 swung by **3.58 TeV** across the phases; r17913 swings by **0.074 TeV**, i.e. flat to the same
±0.02 TeV precision as the 2023 reference (0.038 TeV). And the control that gave R9 away has closed:
the **2 isolated bunches** (⟨FCal⟩ 4.432, ⟨L1TE⟩ 16.161) now agree with the in-train bunches instead of
being the sample's only healthy crossings. Removing the out-of-time crossings removed the failure — the
R9 hypothesis (a crossing sitting on its predecessors' accumulated LAr undershoot) is therefore
*supported* by this production, though still not *proven*, since r17913 changes six flags at once.

**Two residuals, both small and neither a repeat of R1/R4** — recorded because they must not be
silently absorbed into the 2024 efficiencies:

1. **Calorimeter ~2 % low vs the 2023 reference, uniformly.** FCal median 4.29 vs 4.40 TeV (−2.5 %),
   phase-folded mean 4.37 vs 4.43 (−1.3 %), L1TE median 15.70 vs 16.14 (−2.7 %). It is flat in train
   position, so it is a normalisation, not the R9 modulation. Expected sources of a difference of just
   this size: (a) the pinned vertex at **z = −71.2 mm vs −3.3 mm** for r17618 — a 7 cm shift moves the
   FCal η acceptance by Δη ≈ 0.015, a percent-level ΣE_T change, as already estimated in R2; (b) the
   2024 conditions genuinely differ from 2023. **This is a cross-year comparison, not a closure test**,
   so a percent-level offset is not by itself evidence of a defect.
2. **HITight track multiplicity −6.0 %** (2257 vs 2400), while the untagged count is now −2.8 %
   (r17864: +7.3 % untagged, −9.0 % HITight). So r17864's *extra low-quality tracks* are gone — the
   8.5–11.5 k untagged shoulder of R4 has disappeared, consistent with the out-of-time crossings no
   longer contributing — but a HITight deficit of −6 % survives. **The 100-event test (R11) read −1.6 %;
   that was a fluctuation** (±3 % per sample at that size) and the 10 000-event number supersedes it.
   Whether this propagates to the muon reco efficiency is exactly what the R13 comparison measures:
   R2 found −11 % Tight reco ε alongside −9 % HITight under r17864, so a −6 % HITight would predict a
   partial, not complete, recovery.

**The 400 MeV track threshold cannot be lifted from these NTUPs (2026-10-06).** `m_trk_numqual`
(`TrigRates.cxx:1833`) has 8 slots: `[0-3]` = {no quality, PP_MIN_BIAS, HI_LOOSE, HI_TIGHT} **with**
`p_T > 400 MeV`, `[4-7]` = the same four **without** it (added Dec 2023 for exactly that purpose).
The two halves are equal in **10 000 of 10 000 events in both samples, for all four quality levels** —
the loop runs over `InDetTrackParticles` (`TrigRates.cxx:85`) and the HI track reconstruction of these
AODs stores nothing below 400 MeV. So the threshold is **upstream of the skim**: a re-skim cannot lift
it, only a different reconstruction can, and there is no fallback (the per-track `track_pt` branch is
not stored in the overlay NTUPs, `StoreBasic` off). The `trk_multiplicity` plot was accordingly
reduced from 2×2 to **1×2** (`[0]` and `[3]` only) — an earlier label called `[4]`/`[7]` the "stored
collection", which misread their purpose and invited the reading "the p_T cut removes nothing".

**Branch check** (`scripts/check_skim_output.C` vs the pbpb23 pTH125_300 NTUP): 265 vs 225 branches,
**24 missing / 64 extra — byte-for-byte the same bookkeeping difference as r17864**, i.e. the three
2023 `_VTE50` muon chains absent from the 2024 HI menu vs the 2024 `mu8` / `mu10` /
`mu4noL1_hi_uccTh*` / `mu4noL1_L1ZDC_HELT*` chains added. 14 always-empty branches: `ActIntPerXing`,
`AvgIntPerXing` and the 12 L1ZDC-seeded `mu4noL1` decisions (no ZDC simulation in the overlay) — the
known set. All muon / truth / vertex / FCal / trigger branches filled (analysis chains mu4, 2mu4 fill
0.9994 / 0.9952). **No skim-code change was needed for r17913**, as expected for an r17864 clone.

Files: `fcal_sum_et_r17913_vs_r17618_pTH125_300.png`, `l1te_…png`, `trk_multiplicity_…png`,
`vertex_z_reco_vs_truth_r17913.png`, `bcid_structure_…png`, `train_position_…png`,
`centrality_population_…txt`.

### R13. Efficiency chain on r17913 (the R2 redo) — the reconstruction deficit is GONE (2026-09-29)

The chain the FIXED FCal unblocks (user ruling 2026-09-17 made R2 not studiable only because the
r17864 centrality bins were not the same event populations as r17618's; R12 shows r17913 reproduces
the reference centrality population to two decimals, 99.32 % / 0.68 %).

Commands, unchanged from the r17864 round — the run scripts read the canonical, configuration-keyed
NTUP name, which now IS r17913, so nothing in the chain needed editing:
- `OVERLAY_YEAR=24 pipelines/pipeline_pythia_fullsim_overlay.sh hijing` — Condor NTP (cluster 8,
  01:46→01:52) → RDF hist filling (01:52→01:54) → reco-eff / det-response plotters for Medium + Tight.
  All five stages exit 0.
- `OVERLAY_YEAR=24 NTupleProcessingCode/run_pythia_fullsim_overlay_{single_muon,mc_trig,single_muon_mc_trig}.sh`.

Comparison against the existing pbpb23 pTH125_300 slice run
(`..._pbpb23/r17618_pTH125_300_run/`) via `plot_overlay_rtag_effcy_comparison.cxx`
(`RTAG24 = "r17913"`), modes `recoeff` / `trigeff` / `detresp` / `samb`.

**Two bookkeeping notes from the run:**
- **`w_factor` = 2.09729e-05 for pbpb24, now IDENTICAL to the pbpb23 slice.** Under r17864 it read
  2.09733e-05 because the sample's own `ami_info/` (the AOD record, σ 89.541 nb) was still being
  preferred; after the 2026-09-17 evgen ruling both overlays read the same nPDF evgen
  (802776 e8599, σ 89.529 nb × ε_filt 2.342579e-3 / 10 000), which is exactly right — same evgen,
  same beam, same slice must give the same weight. The 2e-5 relative difference was in any case
  irrelevant to R2 (single-slice constant weight, and it cancels in every efficiency ratio).
- The `single_muon` mode prints "the first bin … of the same-sign cut-acceptance histogram has ZERO
  bin content". **Pre-existing and mode-specific**, not an r17913 symptom: the identical line is in
  the archived r17864 log (`logs/ntp_single_muon_r17864.log`) and in both pbpb23 slice logs. The
  single-muon mode builds no pairs, so the SS cut-acceptance histogram is legitimately empty.

**Result (Tight WP, the nominal one; all four modes exit 0).** R2's r17864 column is repeated so the
recovery is readable in one place:

| quantity | **r17913** | r17864 (R2) | r17618 | **r17913 / r17618** |
|---|---|---|---|---|
| single-µ reco ε, 0–5 %, bins 4.5-8 / 8-14 / 14-25 / 25-100 GeV | **0.580 / 0.655 / 0.690 / 0.707** | 0.563 / 0.624 / 0.660 / 0.640 | 0.588 / 0.661 / 0.675 / 0.709 | **0.986 / 0.992 / 1.023 / 0.997** (±0.017) |
| single-µ reco ε, 0–5 %, integrated | **0.6473** (9292/14356) | — | 0.6496 (9679/14899) | **0.996** |
| single-µ reco ε, 5–10 %, integrated | 0.6679 (2554/3824) | (0.83–0.83 ratios) | 0.6994 (2283/3264) | 0.955 — see the migration note |
| single-µ reco ε, **centrality-blind** | **0.6518** (11930/18304) | 0.584 (10681/18304) | 0.6595 (12058/18284) | **0.988** (r17864: **0.885**) |
| mu4 ε, Step-1 selection, 0–5 % | **0.785** (6635/8452) | 0.741 (314/424) | 0.773 (6882/8902) | **1.016** (r17864: 0.958) |
| mu4 ε, centrality-blind | **0.7944** (8630/10864) | 0.761 (7388/9710) | 0.7808 (8654/11084) | **1.017** (r17864: 0.975) |
| (p_T^reco − p_T^truth)/p_T^truth, 0–5 %: mean / RMS | **−0.12 % / 3.03 %** | −0.25 % / 3.24 % | −0.13 % / 2.99 % | RMS **+1.3 %** (r17864: +8 %) |
| η_reco − η_truth RMS | 0.74e-3 | 0.74e-3 | 0.74e-3 | = |
| from_same_b, all OS pairs | 0.3135 ± 0.0057 (6689) | 0.3135 ± 0.0057 (6689) | 0.3200 ± 0.0057 (6593) | consistent |
| from_same_b, Tight OS pairs | 0.2913 ± 0.0083 (**2997**) | 0.2947 (**2389**) | 0.2982 ± 0.0084 (**2998**) | consistent; **pair count 2997 vs 2998** (r17864: −20 %) |

Medium WP, same trees: reco ε centrality-blind 0.6707 vs 0.6775 (0.990), 0–5 % 0.6657 vs 0.6680
(0.997), mu4 ε centrality-blind 0.7898 vs 0.7767 (1.017) — the same picture.

**What this settles.** R2's three reconstruction-level deficits were conditions artefacts, not
physics, and all three close on r17913:
- Tight reco efficiency **−11.5 % → −1.2 %** (centrality-blind, the migration-free number), with
  every coarse-p_T ratio panel consistent with unity.
- mu4 efficiency **−2.5 % → +1.7 %** (r17913 now sits marginally *above* the 2023 reference).
- p_T-resolution width **+8 % → +1.3 %**.
- And the cleanest single statement: the **Tight OS pair count is 2997 vs 2998**, against 2389 under
  r17864. R2 had noted that its −20 % pair deficit was ≈ 0.885² — the square of the single-muon
  efficiency ratio. With the single-muon ratio back at ~0.99, the square is back at ~1.

**Why `all OS reco pairs` (6689) and its `from_same_b` (0.3135) are byte-identical to r17864.** That
row is TRUTH-seeded — the overlay pair trees are built from Pythia truth pairs and then matched to
reco — so its denominator and its `from_same_b` are properties of the generator record, which the
r-tag cannot change (R2 already established truth is r-tag-independent). It is a consistency check,
not a staleness symptom: the *reco-dependent* row in the same table, Tight OS pairs, moved
2389 → 2997 on the same file path with the same macro, which is only possible if the file is the new
one.

**The one non-unity row, 5–10 % at 0.955, is centrality migration, not inefficiency.** R12 measured the
r17913 FCal ~2 % below the 2023 reference; a 2 % downward shift moves events across the percentile
edges, and the populations show exactly that — r17913 has FEWER events in 0–5 % (14356 vs 14899 truth
muons) and MORE in 5–10 % (3824 vs 3264). So r17913's "5–10 %" is contaminated from the *central*
side, and more-central events have lower efficiency, which is the sign observed. The 0–5 % bin, which
loses only its most peripheral events, is unaffected (ratio 0.996), and the centrality-blind row —
immune to migration by construction — is 0.988. **Before any pbpb24 efficiency is used for physics the
FCal centrality calibration must be re-derived for 2024 conditions anyway**; it is currently the 2023
calibration applied to a 2024 sample, which is what this residual is.

Outputs (all in `plots/r17913_rtag_sanity/`): `single_muon_reco_effcy_vs_pt_ctr{0_5,5_10}_r17913_vs_r17618.png`,
`single_muon_mu4_effcy_vs_pt{,_coarse}_ctr0_5_r17913_vs_r17618.png`,
`single_muon_pt_eta_response_ctr0_5_r17913_vs_r17618.png`, `from_same_b_fraction_r17913_vs_r17618.csv`,
plus the nominal pipeline plot sets `plots/hijing_overlay_pbpb24_{reco_effcy_plots,reco_effcy_plots_require_signal_cuts,det_resp_plots}/`.

### R14. Is the 2024-vs-2023 low-side shift also in Pb+Pb DATA? (2026-10-08) — YES, same sign

Plots: `plots/r17913_rtag_compr_PbPb_data/{fcal_l1te,trk_multiplicity}_data_vs_hijing_pbpb23_vs_pbpb24.png`
(2x2, DATA top row / HIJING overlay bottom row, 2023 vs 2024 overlaid, each curve unit-normalised).
Macro `plotting_codes/overlay_rtag_checks/plot_overlay_rtag_compr_pbpb_data.cxx`; data histograms
cached in `pbpb_data_event_level_hists.root` (one pass over the skims, merged from per-part files).

Ratios **2023 / 2024** (a ratio > 1 means 2024 sits LOWER):

| quantity | DATA mean | DATA median | MC mean | MC median |
|---|---|---|---|---|
| FCal Sigma E_T      | 1.0055 | 1.0104 | 1.0143 | 1.0242 |
| L1 Sigma E_T        | n/a (see below) | n/a | 1.0206 | 1.0278 |
| N_trk (no HITight)  | 1.0447 | 1.0494 | 1.0284 | 1.0454 |
| N_trk (HITight)     | 1.0500 | 1.0541 | 1.0634 | 1.0662 |

Absolute values — DATA: FCal mean 2.429 / 2.416 TeV, N_trk 2344 / 2244, N_trk^HITight 1519 / 1447
(2023 / 2024); L1TE 2023 only, mean 7.80 TeV. MC: FCal 4.432 / 4.370 TeV, L1TE 16.26 / 15.94 TeV,
N_trk 5435 / 5285, N_trk^HITight 2400 / 2257.

**Conclusion.** The 2024-over-2023 deficit is NOT an artefact of the cross-year overlay comparison
alone: Pb+Pb data shows the same sign in every quantity available. The track multiplicities agree
in magnitude between data and overlay (4.5-5.4 % in data vs 2.8-6.6 % in MC). The FCal shift is
present in data but roughly HALF the overlay's (0.6-1.0 % vs 1.4-2.4 %) — consistent with the R12
reading that part of the overlay's FCal offset is the cross-vertex (-71.2 vs -3.3 mm) difference
on top of a real year-to-year detector change. Nothing here contradicts r17913's acceptance (R10-R13);
it reinforces the standing requirement for a 2024-conditions FCal centrality calibration.

Procedure (USER DECISION 2026-10-06): data = `HLT_mu4_L1MU3V` + `PbPbBadRuns(yr)` + the full 5-cut
Pb+Pb event selection from each year's own `event_sel_cuts_pbpb_20yy.root` (the
`plot_pbpb_fcal_comparison.cxx::FillFCalHist()` recipe, reproduced cut for cut); **no centrality
cut** (centrality is FCal-derived with the fixed 2023 calibration, so cutting on it would partly
define away the shift under test). Overlay row: no selection, unweighted, as in R12.

**L1TE does not exist in PbPb24 data.** Measured: `L1TE` and `L1TE24` are identically zero in both
2024 parts (0 non-zero in 40.9 M selected events of part1), while 2023 is filled (mean 7.80 TeV).
The 2024 skim carries no L1 sum-E_T VALUE branch, only `L1jTE*` trigger-decision flags: TrigRates
reads the legacy `LVL1EnergySumRoI` (`SkimCode/.../TrigRates.cxx:2099-2104`), which the 2024 HI menu
stopped filling after the move to Phase-I jFEX (same for 2025/2026, see
`pbpb2026_skim_and_lgd_storage.md`). Per user decision the panel keeps the 2023 curve alone. A data
L1TE year comparison needs a re-skim that reads the jFEX energy object.

**Statistics caveat (open).** `data_pbpb23_part2.root` (49.65 GB, /pnfs LOCALGROUPDISK) could not be
read on 2026-10-08: its namespace entry stats fine over both the POSIX mount and `xrdfs`, but every
read wedges the process in uninterruptible disk sleep (0 bytes in 20 min, POSIX and xrootd alike) —
a dCache pool problem, not a missing file. The 2023 row is therefore 4 of 5 parts, 64.2 M selected
events (2024 is complete, 84.2 M). The four available 2023 parts agree to **<= 0.4 %** in every mean
(FCal 2.4277-2.4371, N_trk 2337-2353, N_trk^HITight 1517-1523), i.e. an order of magnitude below the
year difference, so the conclusion is safe; the numbers should still be regenerated once the pool
recovers (`fill:23:1`, then re-hadd and replot).


## Remaining Work

- **The r-tag question is answered: r17913 is acceptable as the recipe.** What is left is production
  scale-up, not validation — the remaining 4 vertex-z points (−28.5 / −4.8 / 18.9 / 61.6 mm) and the
  b intervals 5–9 / 9–12 / 12+ fm, then the 4 beams × 6 slices. `FullSimOverlayConfigTag` becomes a
  user-decided combination rule at that point, never a silent glob.
- **A 2024-conditions FCal centrality calibration is needed before any pbpb24 efficiency is used for
  physics.** Everything above applies the 2023 calibration to a 2024 sample; that is the origin of
  the one non-unity row in R13 (5–10 %, ratio 0.955, pure bin migration) and of R12's uniform ~2 %
  FCal offset. Not a defect of the r-tag.
- Delete the TEMPORARY `allow_missing_slices = (overlay_year == 24)` in the five overlay run scripts
  + the Condor worker once a 6-slice pbpb24 production exists.
- **Report to the production contact** — now a *confirmation* rather than a bug report: the
  r17864 → r17913 preExec change fixed it, the R9 train-position modulation is gone (3.58 → 0.07 TeV
  swing), and the efficiencies are back on the 2023 reference. Worth sending so the six-flag recipe
  is locked in for the full production.
- User decisions still open from R5: whether ε_acc (R5 #2) and the generator-level composition plots
  (R5 #4) move to the pp-suitable PDF evgen via the fullsim truth block.
- The r17864 material is kept, not deleted: `broken_calo_tid52480369/` (first production NTUP +
  all derived products), `broken_calo_r17864_rerun_tid52629267/` (the bit-identical re-run NTUP),
  `plots/r17864_rtag_sanity/`.

## Latest Stage

**2026-10-08 (session 7) DONE — R14 answered: the 2024-vs-2023 low-side shift IS present in Pb+Pb
data**, same sign in FCal and in both track multiplicities; magnitudes in the table in R14. Two
PNGs in `plots/r17913_rtag_compr_PbPb_data/`. Open: regenerate once `data_pbpb23_part2.root` is
readable again (dCache pool wedged; 2023 is 4/5 parts, part-to-part spread <= 0.4 %), and a data
L1TE year comparison remains impossible without a re-skim that reads the Phase-I jFEX energy.
