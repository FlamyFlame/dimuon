# Hand-off: skim session → analysis-code session (pbpb26 lumi & code update)

Written by the skim/storage session (owner: `pbpb2026_skim_and_lgd_storage.md`). The analysis
session owns every code change listed here; this file only states what is on disk and why.
Updated as each grid task lands. **Read the status column — only DONE rows are on disk.**

## 1. Part counts (`file_batch_max`, `.sub` `queue N`, `QUEUE_COUNTS`, `ScrambGen::NParts`, plot part lists)

| year | old | new | status | why |
|---|---|---|---|---|
| PbPb23 | 4 | **5** | **DONE 2026-09-15 21:49 UTC** — `pbpb_2023/data_pbpb23_part5.root` (117 025 entries, on LGD + symlink, verified) | May-2026 skim abandoned 78 files of run 462969 (SiGNET "Service not available"); recovered by task 52569154, `done` 78/78, 0 duplicate events vs parts 1–4 (exact check on the shared run) |
| PbPb25 | 6 | **7** | PENDING — task 52569156 running | 685 abandoned files: run 510510 (352, RPD-only skim bug D6 — ZDC fine), 511020 (55), 511035 (237), 512013 (1), 512049 (40) |
| pp24 | 12 | **13** | PENDING — task 52569161 | 109 abandoned files, run 488534 (CSCS batch `LRMS error`) |
| PbPb26 | 0 (unset) | **7** | **DONE 2026-09-16 13:40 UTC** — `pbpb_2026/data_pbpb26_part{1..7}.root`, 270 087 106 entries; **every run 100 %** (35 runs, 111 279/111 279 files across the union of the 7 tasks) ⇒ **lift P14 (provisional-lumi warning): GRL total 2.62316 nb⁻¹ is final**; remove the `ScrambGen::NParts(26)` throw; part set 1–7 contiguous | **Correction to the 2026-09-14 note "there is no part 7": there IS a part 7** (task 52568862, 244 files: 169 of run 522546 lost to SiGNET disk failures + 75 of runs 522200/522949 lost to the absent-ZDC crashes, now skipped per D9). 2026 parts are LOCAL (not on LGD) |

Code sites (from the skim doc D10 entry; verify against the tree): `PbPbExtras.c`
`run_year_to_file_batch_max_map` `{23,4}→{23,5}`, `{25,6}→{25,7}`, `{26,0}→{26,7}`;
`run_pbpb_23*.sub` `queue 4→5`, `run_pbpb_25*.sub` `queue 6→7`, `run_pbpb_26*.sub` `queue 7`;
`PPExtras.c` pp24 `file_batch_max 12→13`, `run_pp_24*.sub` `queue 12→13`. Set each year's
sites together and run `pipelines/preflight_pbpb_year.sh <yr>`.

**Blast radius:** every downstream PbPb23 (now), PbPb25 / pp24 / PbPb26 (when their rows go
DONE) result is stale until the NTuple processing is rerun with the new part — the missing
files were part of the luminosity denominator all along (`signal_selection_change_impact.md`).

## 2. Branch difference of the recovery parts

Recovery parts carry **169 branches = the May parts' 168 + `muon_match_L1MU3V`** (added
2026-07-22, purely additive). Anything that reads a branch list from the first part of a year
is unaffected; anything that requires identical branch sets across parts must tolerate this.

## 3. Luminosity — USER RULINGS 2026-09-15 (nothing to implement)

- **No LB exclusion, luminosity unchanged.** Events with missing required ZDC data are skipped
  at skim time (D9); the luminosity stays the full GRL value because the loss is negligible at
  dataset level. LBs with a *considerable* missing-ZDC fraction are reported to the user for
  escalation to the ZDC experts, not cut. ⇒ `PbPbMu4SampledLumiNb(26)` = GRL total 2.62316 nb⁻¹
  is final once part 7 lands (lift the provisional warning then).
- **No `zdc_ZdcModuleMask == 255` requirement, in any year.** A missing module pulse is shower
  physics (neutron interaction depth ~e^(−L/λ); energy depleted before the last module).
- PbPb23 part5: all 35 run jobs exit 0, **0 ZDC-skipped events**.
- 2026 run 522200 LBs ~250–254: ZDC TTC restart (ZdcCalib rate ~0), known to the ZDC experts;
  events without ZDC are dropped by the skim; luminosity influence to be ignored (user, 2026-09-15).
  Same for 2026 run 522949 LBs ~180–186 and 2025 run 510510 LB 174 (checked, known issues).
  Standing rule: a contiguous few-LB cluster of missing-ZDC events whose fraction of the TOTAL
  dataset is negligible = expected detector-incident pattern → no report, no lumi change.

## 4. Provenance of the recovery parts

AthAnalysis **25.2.89** (same as the May skim of each year) with the ZDC-fixed `TrigRates`
(bit-identical output on 2025 data, 169/169 branches). `--inputFileList` = exactly the files
PanDA marked not-finished; file-level disjointness vs every merged part verified before
submission; event-level `(RunNumber, eventNumber)` uniqueness verified after merge.
