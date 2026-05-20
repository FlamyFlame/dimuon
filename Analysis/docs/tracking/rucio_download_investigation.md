# Investigation: rucio download failure

## Objective
Determine why `rucio download` fails with `[Errno 5] Input/output error` on
Python stdlib files, while all other rucio commands (`list-dids`, `whoami`)
and ROOT/hadd work fine in the same environment. Fix it so
`grid_monitor.sh` can complete the May 2026 skim downloads.

## Context
- 15 grid tasks (May 2026 skim) completed on the grid, need downloading
- `grid_monitor.sh` calls `rucio download` which fails consistently
- Error: `[Errno 5] Input/output error: '/cvmfs/sft.cern.ch/.../python3.11/getopt.py'`
- Same failure observed on attsub05 and attsub06
- `source ~/setup.sh && lsetup rucio` works; `rucio list-dids` works; ROOT/hadd works
- Only `rucio download` fails

## Sub-steps
1. Reproduce `rucio download` failure in a controlled test
2. Compare what `rucio download` does differently from `rucio list-dids`
3. Test with verbose flag and `--ndownloader 1`
4. Identify fix and verify with a single-dataset test download

## Accumulated Findings

- 2026-05-19: `rucio download` of pbpb23 part4 (33 files) succeeded at
  ~400 MB/s from BNL-OSG2_SCRATCHDISK (dCache). All 33 files downloaded,
  verified non-empty.
- The `[Errno 5]` error was **transient** — not persistent CVMFS corruption.
  All CVMFS-dependent tools (ROOT, hadd, rucio list-dids) worked throughout.
- `rucio -v` and `rucio --verbose` are NOT valid for the download subcommand
  in rucio 35.6.0. Both produce "unrecognized arguments" errors.
- The operational problem was in `grid_monitor.sh`: any download failure
  (including transient) was marked as permanent "failed" with no retry.

## Ruled Out

- **CVMFS cache corruption**: ROOT, hadd, rucio list-dids all worked fine
  on the same nodes (attsub05, attsub06). BNL SDCC nodes are well maintained.
- **Rucio authentication/config issue**: `rucio list-dids`, `rucio whoami`
  work. Proxy valid (85+ hours remaining).
- **Storage quota issue**: ~909 GB free, first dataset (~25 GB) would fit.
  Error was I/O on Python stdlib, not disk space.

## Resolution

**Root cause**: Transient CVMFS I/O error during `rucio download`'s heavier
Python import path (multiprocessing, subprocess). Not reproducible on retry.

**Fix applied**: Added retry logic to `grid_monitor.sh` — up to 3 attempts
per task with 60s delay between retries, instead of permanent failure on
first error.

**State file** (`grid_monitor_state.txt`): Reset to all-pending for 15 tasks.

**Next step**: Rerun `grid_monitor.sh may2026_skim.txt` in tmux.
