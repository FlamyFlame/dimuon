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
| N_trk, p_T > 400 MeV (`[0]` ≡ `[4]`: the stored track collection already has p_T > 400 MeV) | 5830 ± 1740 (RMS) | 5435 ± 1460 | 1.073 |
| N_trk, HITight (`[3]` ≡ `[7]`) | 2183 ± 359 | 2400 ± 362 | 0.910 |
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

## Remaining Work

- User decision on the verification plan / on reporting R1 to the production contact (message draft offered).
- The verification itself (items 1-6 of the proposal) is NOT started; item 1 (calorimeter) is already answered negatively by R1.

## Latest Stage

2026-09-17 (session 3) — three user tasks, in progress:
1. **NTUP naming.** `FullSimHIJINGOverlayPP24` → `FullSimHIJINGOverlayPbPb<yy>.<cfg>` with
   `<cfg> = vtxz<z>mm_b<lo>_<hi>fm` (z in mm, '.'→'_', sign kept; b = HITS `ip` range):
   pbpb24 `FullSimHIJINGOverlayPbPb24.vtxz-71_2mm_b0_5fm`, pbpb23 `FullSimHIJINGOverlayPbPb23.vtxz-3_3mm_b0_5fm`
   (+ `_r17662` variant). Rename the on-disk files (both dirs; `.bak_*` LGD symlinks untouched — remote names),
   repoint the diagnostic run-dir symlinks, `FullSimSampleFileTag(t, pbpb_year)` + `FullSimOverlayConfigTag`,
   `PythiaAlgCoreT.c` InitInputFullsim, `grid_monitor.sh` outDS→file parsing, the three overlay `grid_sub*.sh`
   (OUT_TAG + `_cfg_tag` helper + comment), `plot_r17864_event_level.cxx`, SkimCode/README + pythia_fullsim_overlay.md.
   Verify: recompile, rerun the pbpb24 single-muon NTP (1 slice, minutes) and the pbpb23 slice script → same entry counts.
2. **Track multiplicity plot** `trk_numqual[0,3,4,7]` (4 panels, 2024 vs 2023 pTH125_300, no centrality cut) →
   `plots/r17864_rtag_sanity/trk_multiplicity_r17864_vs_r17618_pTH125_300.png` (new mode `ntrk`).
3. **L1TE plot** (1 panel) → `plots/r17864_rtag_sanity/l1te_r17864_vs_r17618_pTH125_300.png` (mode `l1te`; L1TE is GeV,
   2023 extends to ~25 TeV → axis 0-25 TeV).
Reviews: /review-analysis-code for (1), /review-plot for (2)+(3). Then commit.
