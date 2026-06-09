# HIJING Overlay r17618: Grid Reprocessing & NTUP Replacement

## Objective

Reprocess all 6 pT slices of the Pythia fullsim HIJING overlay sample
using the corrected r17618 AODs (beamspot fix, clone of r17044).  Replace
the old (wrong vertex-z) NTUPs with the new skimmed output.

## Context

The existing 6 NTUP files in `pythia_fullsim_hijing_overlay_test_sample/`
were produced from r17044 AODs, which had incorrect beam-spot z position.
New r17618 AODs fix this.  The r-tag comparison investigation
(`hijing_overlay_truth_barcode_duplicate_investigation.md`) established
that r17618 (with full HIJING truth) is the production configuration —
r17662 (StandardSignalOnlyTruth) eliminated truth contamination but
r17618 is the baseline since downstream code handles barcode duplication
via `emplace` first-writer-wins.

## Datasets

All 6 pT slices, r17618 reco over r15970 HIJING digitization:

| pT slice | AOD dataset | Grid task ID |
|----------|-------------|-------------|
| pTH8_14 | `mc23_5p36TeV.802781.…e8599_s4614_r17618_r15970` | 50774690 |
| pTH14_24 | `mc23_5p36TeV.802777.…e8599_s4614_r17618_r15970` | 50774664 |
| pTH24_40 | `mc23_5p36TeV.802778.…e8599_s4614_r17618_r15970` | 50774671 |
| pTH40_70 | `mc23_5p36TeV.802779.…e8599_s4614_r17618_r15970` | 50774677 |
| pTH70_125 | `mc23_5p36TeV.802780.…e8599_s4614_r17618_r15970` | 50774683 |
| pTH125_300 | `mc23_5p36TeV.802776.…e8599_s4614_r17618_r15970` | 50774658 |

Output container pattern: `user.yuhang.NTUP.Pythia_5p36TeV_pp_hQCD_DiMu_<pT>.FullSimHIJINGOverlayPP24.June2026.v1.`

## Scope

- Grid submission: `SkimCode/run_pythia_fullsim_HIJING_overlay/grid_sub.sh` (updated r17044→r17618, VER_TAG→June2026.v1)
- Grid monitoring: `SkimCode/scripts/grid_monitor.sh --mode overlay <task_ids>`
- Hadd: `/gpfs/mnt/atlasgpfs01/usatlas/data/yuhanguo/hadd_merge_datasets.sh` (configured for overlay)
- NTUP destination: `/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_hijing_overlay_test_sample/`
- Old NTUPs: backed up to `backup_wrong_vtx_z/` subdirectory

## Design Decisions

1. **grid_monitor.sh --mode overlay** (2026-06-07): Added `--mode` flag to
   `grid_monitor.sh` instead of writing a standalone overlay monitor script.
   Mode controls: DATA_BASE path, RECORD_FILE name, `map_outds()` naming
   logic, code-update skip. All shared infrastructure (BigPanDA polling,
   rucio download, hadd, validation, parallel coordination) is mode-agnostic.

## Implementation Plan

1. ✅ Update `grid_sub.sh`: r17044→r17618 AODs, VER_TAG=June2026.v1
2. ✅ Submit 6 grid jobs via pathena (all submitted 2026-06-07)
3. ✅ Back up old NTUPs to `backup_wrong_vtx_z/`
4. ✅ Add `--mode overlay` to `grid_monitor.sh` (map_outds, apply_mode_config, skip chunked fallback & code-update)
5. ✅ Configure `hadd_merge_datasets.sh` for overlay (DIR, SUBDIRS, SUBDIR_OUTPUT_NAMES)
6. ✅ Run `grid_monitor.sh --mode overlay` to poll, download, hadd, validate
7. ✅ Verify 6 new NTUP files exist with correct entry counts (all 10,000 entries)
8. ✅ Update `merging-record.txt` (auto by grid_monitor, 6 groups sorted alphabetically)

## Progress Log

### 2026-06-07: Steps 1–5

- Updated `grid_sub.sh`: changed 6 inDS from r17044 to r17618, VER_TAG April2026.v1→June2026.v1.
- Submitted 6 grid jobs: 50774658, 50774664, 50774671, 50774677, 50774683, 50774690. All confirmed `running` (0%) on BigPanDA.
- Moved 6 old `*.FullSimHIJINGOverlayPP24.NTUP.root` files to `backup_wrong_vtx_z/`.
- Added `--mode` flag to `grid_monitor.sh`:
  - `apply_mode_config()`: sets DATA_BASE, RECORD_FILE, log/state paths per mode.
  - `map_outds()`: overlay mode extracts `Pythia_5p36TeV_<sample>.<tag>.NTUP` from outDS (flat dir, no subdirs).
  - `get_code_update_info()`: returns empty for overlay (no C++ source auto-update).
  - `process_task()`: overlay skips chunked_hadd_fallback (fails immediately if recursive_hadd fails).
- Configured `hadd_merge_datasets.sh`: DIR→overlay repo, 6 SUBDIRS with June2026.v1 outDS names, SUBDIR_OUTPUT_NAMES mapping each to final NTUP filename.

## Remaining Work

None — all steps complete, plots reviewed and approved.

## Latest Stage

**Complete.** All 8 implementation steps done. Pipeline produced 964 plots across 4 directories. `/review-plot` passed at iteration 1 with zero CRITICAL/WARNING issues (2 INFO cosmetics: raw var names in det_resp titles, mild legend overlap in minv comparison).

### 2026-06-07: Pipeline execution

- First attempt: Condor job OOM at 4 GB after ~40 min (126 GB input too large for 4 GB limit).
- Fixed: `run_pythia_fullsim_overlay.sub` `request_memory` 4GB → 8GB.
- Second attempt: `pipeline_pythia_fullsim_overlay.sh hijing` completed successfully (22:02–22:54).
- Stage 1 (NTP): Condor cluster 2840, ~50 min runtime → `muon_pairs_pythia_fullsim_hijing_overlay_pp24_no_data_resonance_cuts.root` (38 MB).
- Stage 3 (RDF): `histograms_pythia_fullsim_hijing_overlay_pp24_no_data_resonance_cuts.root` (3.6 MB).
- Stage 5 (Plotting): 964 plots in 4 dirs (reco_effcy: 300, det_resp: 72, signal_cuts: 350, no_ip: 150).

### 2026-06-07: Steps 6–8

- Ran `grid_monitor.sh --mode overlay` from 16:58 to 21:18 UTC.
- All 6 tasks completed: 50774658 (pTH125_300, first done ~17:45), 50774677 (pTH40_70), 50774683 (pTH70_125), 50774664 (pTH14_24), 50774690 (pTH8_14), 50774671 (pTH24_40, last done ~21:16).
- Each task: 5 grid output files → rucio download → hadd → 10,000 entries validated → cleanup.
- 6 NTUP files in `/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_hijing_overlay_test_sample/`:
  - `Pythia_5p36TeV_pp_hQCD_DiMu_pTH{8_14,14_24,24_40,40_70,70_125,125_300}.FullSimHIJINGOverlayPP24.NTUP.root`
- `merging-record.txt` auto-updated and reorganized (6 groups, sorted alphabetically).
