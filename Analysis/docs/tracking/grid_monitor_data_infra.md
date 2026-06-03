# Grid Monitoring & Data Infrastructure

## Objective
Build robust infrastructure for May 2026 skim: grid task monitoring, parallel downloading, hadd with fallbacks, and automatic NTuple processing code updates when datasets are split into chunks.

## Technical Procedure
(Tooling/infrastructure — substitutes for Physics Procedure.)

1. **Grid monitoring**: Poll BigPanDA for task status, download via rucio, hadd, validate entry counts, update data-merging-record.txt.
2. **File naming convention**: All data files use `data_{prefix}_part{N}.root` (e.g., `data_pbpb23_part1.root`, `data_pp24_part1.root`). No subdirectories. PP and PbPb use identical conventions.
3. **hadd strategy**:
   a. Try direct hadd of all files.
   b. If fails, try recursive_hadd (split to ≤50-file chunks, merge bottom-up into single output).
   c. If recursive_hadd fails, use chunked fallback: split into separate partN files, each from ≤HADD_CHUNK_MAX_FILES inputs.
4. **Chunked fallback rules**:
   a. First chunk keeps original part number; additional chunks get numbers above current max → ensures contiguous 1..max range.
   b. data-merging-record.txt entries note "(chunk X/N of partM)".
   c. C++ source file_batch_max and .sub queue count are updated via deterministic sed (no claude -p).
   d. Code changes are git committed automatically.
5. **Code-to-dataset mapping** (for auto-update):
   - pp24 → PPExtras.c `{24, N}`, run_pp_24.sub `queue N`
   - pbpb23 → PbPbExtras.c `{23, N}`, run_pbpb_23.sub `queue N`
   - pbpb24 → PbPbExtras.c `{24, N}`, run_pbpb_24.sub `queue N`
   - pbpb25 → PbPbExtras.c `{25, N}`, run_pbpb_25.sub `queue N`
6. **NTuple processing**: Each file_batch index maps to exactly one ROOT file: `data_{prefix}_part{file_batch}.root`. The .sub `queue N` submits N condor jobs (Process 0..N-1 → file_batch 1..N).

## Context
- May 2026 skim: 15 grid tasks, 137 datasets, 4 dataset groups (PbPb 23/24/25, pp24)
- Task 50267518 (pp24 part1, 500 files) previously failed hadd — triggered recursive hadd and chunked fallback
- Old pp24 code used hardcoded subdirectory structure from old skim — replaced with uniform partN convention

## Scope
- `SkimCode/scripts/grid_monitor.sh` — monitoring, downloading, hadd, fallbacks, code auto-updates
- `Analysis/NTupleProcessingCode/PPExtras.c` — PP input file pattern and file_batch_max
- `Analysis/NTupleProcessingCode/PbPbExtras.c` — PbPb file_batch_max (pattern already correct)
- `Analysis/NTupleProcessingCode/run_{pp_24,pbpb_23,pbpb_24,pbpb_25}.sub` — queue counts

## Design Decisions
1. **Uniform partN convention (Step 5)**: PP now uses same `data_pp{YY}_part{N}.root` as PbPb. Old subdirectory switch block removed from PPExtras.c. Reason: consistency enables shared fallback logic; old subdirs were from a deprecated skim.
2. **Map-based file_batch_max (Step 5)**: PP changed from ternary to `std::map<int,int>` (matching PbPb). Reason: enables deterministic sed-based auto-update from grid_monitor.sh.
3. **Recursive hadd before chunked fallback (Step 7)**: Try recursive merge first (produces single output). Only split into separate files if recursive merge fails entirely. Reason: single file preferred for simpler downstream handling.
4. **Deterministic bash for code updates (Step 7-8)**: No claude -p calls on critical path. Renaming, record update, code update all via sed/bash. Reason: reliability over convenience.

## Implementation Plan
1. ✅ grid_monitor.sh: basic monitoring, BigPanDA polling, download, hadd, validation, record update
2. ✅ grid_monitor.sh: multi-node parallel support with flock coordination, compare-and-swap states
3. ✅ grid_monitor.sh: recursive hadd (≤50-file chunks, bottom-up merge)
4. ✅ Bug fix: `local` outside function → false "completed" marking
5. ✅ PPExtras.c: uniform partN pattern + map-based file_batch_max — per §2, §6
6. ✅ run_pp_24.sub: queue 4 → queue 3 — per §6
7. ✅ grid_monitor.sh: chunked fallback function — per §3c, §4
8. ✅ grid_monitor.sh: auto-update of C++ source and .sub — per §4c, §5
9. [ ] Clean up orphaned pp24 files (data_pp24_part1_1.root, data_pp24_part1__1.root — 366-368 bytes, corrupt)

## Progress Log

### Prior conversations (2026-05-13 to 2026-05-19)
- **Step 1**: Built grid_monitor.sh with BigPanDA polling, rucio download, hadd, entry-count validation, data-merging-record.txt updates. Initially sequential single-node.
- **Step 2**: Rewrote for multi-node parallel execution. flock-based file locking on GPFS. Compare-and-swap state transitions (pending→ready→downloading→completed|failed). Stale download recovery (3h timeout). Worker claims one task at a time.
- **Step 3**: Added half-merge fallback → later replaced with full recursive hadd (CHUNK=50).
- **Step 4 (bug fix)**: `local` used outside function → `success` variable empty → `if $success` always true → 9 tasks falsely marked completed. Fix: removed `local`, explicit string comparison.
- Task 50267518 (pp24 part1, 500 files) failed hadd 3x. Three corrupt partial outputs identified (data_pp24_part1_1.root, __1.root — 366-368 bytes).
- Old download dir (500 files, 396 GB) preserved; rucio download skips via checksum.
- Separated download from hadd in process_task so retries don't re-download.
- Removed outer retry loop (recursive hadd handles splitting internally).
- 8 tasks completed, 50267518 reset to pending, 6 tasks still on grid.

### Current session (2026-05-19)
- User identified: (1) retry restarts from download, (2) half-merge insufficient, (3) no further fallback.
- Fixed retry/fallback: download runs once, recursive_hadd handles splitting, no outer retry.
- User requested: uniform partN for PP, chunked fallback keeping separate files, auto-update of code & .sub.
- **Step 5**: PPExtras.c — replaced hardcoded switch/subdir with `data_pp{YY}_part{N}.root` (matches PbPb). Map-based `file_batch_max = {{17, 3}, {24, 3}}`. AccessPathName check added.
- **Step 6**: run_pp_24.sub — `queue 4` → `queue 3` (new skim has 3 parts).
- **Step 7**: grid_monitor.sh — added `chunked_hadd_fallback()`: splits files into ≤100-file chunks, uses `recursive_hadd` per chunk, assigns contiguous part numbers (original part keeps its number, extras go above max), validates entry counts per chunk and total, records "(chunk X/N of partM)" in record file.
- **Step 8**: grid_monitor.sh — added `get_code_update_info()` (maps target_subdir → Extras file, .sub, run_year_key) and `update_source_for_new_max()` (sed-based update of `{KEY, N}` in C++ and `queue N` in .sub, with git commit). Called automatically at end of chunked fallback.

## Results & Observations
- **State** (2026-05-19): 8 completed, 1 downloading (50267518 on attsub05), 6 pending
- **Completed**: 50267236, 50267250, 50267265, 50267371, 50267396, 50267482, 50267490, 50267539
- **Downloading**: 50267518 (pp24 part1, 500 files)
- **Pending**: 50270743 (pbpb23 part1), 50267423/435/444/472 (pbpb25 parts 1-4), 50267529 (pp24 part2)
- **Orphaned files**: data_pp24_part1_1.root (366B), data_pp24_part1__1.root (368B) — corrupt, to clean up

## Remaining Work
- Step 9: clean up orphaned files after monitor finishes current task
- Monitor remaining grid tasks to completion
- Verify all data files after all downloads complete

## Latest Stage
**Steps 5-8 complete.** All code changes implemented:
- PPExtras.c: uniform `data_pp{YY}_part{N}.root` pattern, map-based `file_batch_max` `{{17, 3}, {24, 3}}`
- run_pp_24.sub: `queue 3`
- grid_monitor.sh: `chunked_hadd_fallback()` + `update_source_for_new_max()` + `get_code_update_info()`
- Remaining: Step 9 (clean orphaned files after monitor finishes current task)
