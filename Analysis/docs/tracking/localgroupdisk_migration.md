# LOCALGROUPDISK Storage Migration

## Objective

Migrate ~1.4 TB of analysis data files (pythia_truth, powheg_fullsim, dimuon_data) from personal pnfs/dCache storage (`/pnfs/usatlas.bnl.gov/users/yuhanguo/`) to the RSE `BNL-OSG2_LOCALGROUPDISK`, then replace the current `~/dcachearea` symlink with one pointing at the LOCALGROUPDISK pnfs path so analysis code continues to work transparently.

## Context

- Current storage: `~/dcachearea` → symlink → `/pnfs/usatlas.bnl.gov/users/yuhanguo/`
- Three top-level directories: `pythia_truth_full_sample` (359 GB), `powheg_full_sample` (439 GB), `dimuon_data` (594 GB)
- Analysis code in `Analysis/NTupleProcessingCode/` references paths under `~/dcachearea` (files: PythiaAlgCoreT.c, PowhegAlgCoreT.c, check_muon_count_per_event.cxx)
- Source instructions:
  - Upload procedure: `rucio_upload_to_bnl_localgroupdisk.pdf` (Luke, Oct 2025), extracted to `rucio_upload_to_bnl_localgroupdisk_steps.md` and `_reference.md`
  - Access/PFN procedure: `bnl_localgroupdisk_pfn_root_instructions.md` (Luke, May 19 2026), with Condor templates in `mattermost_luke_may19_attachments/`

## Combined Procedure Summary (from Luke's instructions)

### Part 1: Upload to LOCALGROUPDISK

Users cannot write directly to LOCALGROUPDISK. The path is: upload to scratchdisk first → add Rucio replication rule to LOCALGROUPDISK.

**Workflow A — data still on grid scratchdisk (e.g. recent skim):**
1. Check for expired files: `rucio list-files <container_name>`
2. Add a Rucio replication rule (UI at rucio-ui.cern.ch or CLI):
   - UI: enter `user.<username>:<container_name>`, select RSE `BNL-OSG2_LOCALGROUPDISK`, group all datasets, leave lifetime empty (permanent)
   - CLI: `rucio add-rule <scope>:<did> 1 BNL-OSG2_LOCALGROUPDISK`
3. Wait for rule state → `OK`

**Workflow B — local files not on grid (expired or self-produced):**
1. Check scratchdisk quota: `rucio list-account-usage <username>`
2. Upload to scratchdisk:
   ```
   rucio upload --rse BNL_OSG2-SCRATCHDISK --scope user.yuhang --register-after-upload --recursive <local_path>
   ```
3. If `Database error` (duplicate DID): rename files with `_COPY` suffix, retry upload
4. Create a dataset: `rucio add-dataset user.yuhang:<new_dataset_name>`
5. Attach files: `rucio attach user.yuhang:<new_dataset_name> *.root`
6. Add replication rule to LOCALGROUPDISK (same as Workflow A step 2)

### Part 2: Access files on LOCALGROUPDISK

1. Set up Rucio + valid grid proxy (`voms-proxy-init -voms atlas`)
2. Generate PFN list:
   ```
   rucio list-file-replicas <scope:DID> --protocols root --pfns --rses BNL-OSG2_LOCALGROUPDISK > pfns.txt
   ```
   Each line: `root://dcgftp.usatlas.bnl.gov:1094//pnfs/<path>.root`
3. Use PFNs as if local paths in TChain (requires valid proxy for full runtime)
4. Optional: `export SITE_NAME=BNL-ATLAS` may speed up Rucio operations
5. For Condor: forward grid proxy via `transfer_input_files`, use `getenv=True`

## Decisions

- **Workflow B** confirmed: original grid job outputs have expired from scratchdisk. Must re-upload from local pnfs → scratchdisk → replicate to LOCALGROUPDISK.
- **One Rucio dataset per logical file group** (e.g. powheg_bb_truth, powheg_cc_truth, dimuon_pbpb23, etc.). Not one giant dataset.
- **Access method:** Symlink farm (preferred). For each dataset, create a directory of symlinks where each symlink points to the file's actual pnfs path on LOCALGROUPDISK (extracted from PFN by stripping `root://dcgftp.usatlas.bnl.gov:1094/`). No grid proxy needed for reading, no code changes. Prerequisite: verify `/pnfs/usatlas.bnl.gov/atlaslocalgroupdisk/` is mounted on BNL SDCC nodes.
- **Luke's EosSubmit note** is CERN lxplus-specific (special Condor schedd for EOS file transfer). Does not apply to BNL SDCC.

## Remaining Open Questions

1. **RSE name inconsistency.** Upload PDF uses `BNL_OSG2-LOCALGROUPDISK` in one place but `BNL-OSG2_LOCALGROUPDISK` in step slides. Access instructions consistently use `BNL-OSG2_LOCALGROUPDISK`. Need to verify via `rucio list-rses`.

2. **Quota.** PDF says 50 TB per user. Our total is ~1.4 TB, should be fine. Worth verifying before starting.

3. **Duplicate DID risk.** If any filenames were previously registered in Rucio (from the original grid jobs), `rucio upload` will fail with `Database error`. Need to check before uploading. The `_COPY` rename hack from Luke is a last resort — changes filenames, could break downstream code.

4. **Grid proxy lifetime.** Default 12h may not be enough for 1.4 TB upload. May need `voms-proxy-init -voms atlas -valid 96:00`.

5. **LOCALGROUPDISK pnfs mount path.** PFNs contain `/pnfs/...` after `root://dcgftp.usatlas.bnl.gov:1094/`. Need to confirm this pnfs path is directly readable on BNL nodes (for symlink farm approach).

## Pre-flight Results (Step 2–3)

- RSE names confirmed: `BNL-OSG2_LOCALGROUPDISK`, `BNL-OSG2_SCRATCHDISK`
- VOMS proxy: ~96h remaining (ample)
- Scratchdisk quota: 50 TB, 110 GB used → plenty of room
- LOCALGROUPDISK pnfs mount: **`/pnfs/usatlas.bnl.gov/LOCALGROUPDISK/`** (not `atlaslocalgroupdisk`). Rucio user files land under `/pnfs/usatlas.bnl.gov/LOCALGROUPDISK/rucio/user/yuhang/<hash>/<filename>`. No `yuhang` dir yet (first upload).
- DID conflict check: `mc_truth_bb_01.root` through `_03` not registered under `user.yuhang` scope → no duplicate DID risk for this dataset
- Test dataset: `bb_evgen_truth_full_sample/`, 25 files, 98 GB total
- **Upload does NOT delete/overwrite source**: `rucio upload` copies files to scratchdisk; original pnfs files are untouched → **no backup needed**

### Analysis code path structure (PowhegAlgCoreT.c:27-39,91-95)

Non-local truth path constructed as:
```
/pnfs/usatlas.bnl.gov/users/yuhanguo/powheg_full_sample/bb_evgen_truth_full_sample/mc_truth_bb_XX.root
```
= `powheg_dcache_dir` + `data_subdir` + filename

Symlink farm must create: a directory at the same relative path under `~/dcachearea` containing symlinks named `mc_truth_bb_XX.root` → actual LOCALGROUPDISK pnfs paths.

## Pilot Test: `bb_evgen_truth_full_sample` (25 files, 98 GB)

### Implementation Plan

**Step P1: Upload to scratchdisk**
```bash
rucio upload --rse BNL-OSG2_SCRATCHDISK --scope user.yuhang \
  /pnfs/usatlas.bnl.gov/users/yuhanguo/powheg_full_sample/bb_evgen_truth_full_sample/*.root
```
STOP if any upload error. Source files are NOT modified.

**Step P2: Create dataset & attach files**
```bash
rucio add-dataset user.yuhang:powheg_bb_evgen_truth_pilot
rucio attach user.yuhang:powheg_bb_evgen_truth_pilot user.yuhang:mc_truth_bb_*.root
```

**Step P3: Add replication rule to LOCALGROUPDISK**
```bash
rucio add-rule user.yuhang:powheg_bb_evgen_truth_pilot 1 BNL-OSG2_LOCALGROUPDISK
```
Check status: `rucio rule-info <rule_id>` — wait for state `OK`.

**Step P4: Get PFNs and extract pnfs paths**
```bash
rucio list-file-replicas user.yuhang:powheg_bb_evgen_truth_pilot \
  --protocols root --pfns --rses BNL-OSG2_LOCALGROUPDISK > pfns_bb_truth.txt
```
Each line: `root://dcgftp.usatlas.bnl.gov:1094//pnfs/usatlas.bnl.gov/LOCALGROUPDISK/rucio/user/yuhang/<hash>/mc_truth_bb_XX.root`
Extract pnfs path: strip `root://dcgftp.usatlas.bnl.gov:1094/` prefix.

**Step P5: Build symlink farm**
```bash
FARM_DIR=~/dcachearea/powheg_full_sample/bb_evgen_truth_full_sample_lgd
mkdir -p "$FARM_DIR"
while read -r pfn; do
  pnfs_path="${pfn#root://dcgftp.usatlas.bnl.gov:1094/}"
  filename=$(basename "$pnfs_path")
  ln -s "$pnfs_path" "$FARM_DIR/$filename"
done < pfns_bb_truth.txt
```
Verify: `ls -la "$FARM_DIR"` — each symlink should resolve.

**Step P6: Test analysis read**
- Temporarily point `powheg_dcache_dir` in `PowhegAlgCoreT.c` to use the `_lgd` directory
- Run a quick test: open one file via ROOT, check TTree entries match original
- Full test: run powheg truth analysis with `file_batch=1` (reads 5 files)

### Rollback / Recovery

| Failure point | What happened | Recovery |
|---|---|---|
| P1 upload fails | Files not on scratchdisk | Nothing to clean up; source files untouched |
| P2 dataset/attach fails | Files on scratchdisk but no dataset | `rucio delete-did user.yuhang:mc_truth_bb_*.root`; retry |
| P3 rule fails | Dataset exists, no LOCALGROUPDISK copy | Delete rule, retry or investigate |
| P4 PFNs empty | Replication not complete or failed | Wait longer; check `rucio rule-info` |
| P5 symlinks broken | Wrong pnfs path | `rm -rf $FARM_DIR`; re-derive paths |
| P6 analysis fails | Path mismatch or file corruption | Revert `PowhegAlgCoreT.c`; original data still at old path |

**At no point is original data at `~/dcachearea/powheg_full_sample/bb_evgen_truth_full_sample/` modified or deleted.** Full rollback = delete Rucio rule + dataset + `$FARM_DIR`, revert code changes.

### Replication Recipe (for applying to other datasets)

Once pilot succeeds, apply to any dataset with:
```bash
# 1. Upload
rucio upload --rse BNL-OSG2_SCRATCHDISK --scope user.yuhang /path/to/files/*.root

# 2. Dataset
rucio add-dataset user.yuhang:<dataset_name>
rucio attach user.yuhang:<dataset_name> user.yuhang:<file_pattern>

# 3. Replicate
RULE_ID=$(rucio add-rule user.yuhang:<dataset_name> 1 BNL-OSG2_LOCALGROUPDISK)
# Wait: rucio rule-info $RULE_ID → state OK

# 4. PFNs → symlink farm
rucio list-file-replicas user.yuhang:<dataset_name> \
  --protocols root --pfns --rses BNL-OSG2_LOCALGROUPDISK > pfns.txt
mkdir -p <farm_dir>
while read -r pfn; do
  pnfs_path="${pfn#root://dcgftp.usatlas.bnl.gov:1094/}"
  ln -s "$pnfs_path" "<farm_dir>/$(basename "$pnfs_path")"
done < pfns.txt
```

## Sub-steps (updated)

1. **Read & summarize** Luke's instructions ✓ DONE
2. **Pre-flight checks** ✓ DONE (RSE, quota, DID, proxy, pnfs mount all verified)
3. **Pilot P1**: Upload bb_evgen_truth to scratchdisk ✓ DONE
4. **Pilot P2**: Create dataset, attach files ✓ DONE
5. **Pilot P3**: Add replication rule to LOCALGROUPDISK ✓ DONE (rule `86adb3d14cfe489ca6311c0ce587a123`)
6. **Pilot P4**: Wait for replication, get PFNs ✓ DONE
7. **Pilot P5**: Build symlink farm ✓ DONE
8. **Pilot P6**: Test analysis read ✓ DONE
9. **Plugin**: Create `bnl-localgroupdisk` Claude Code plugin ✓ DONE
10. **Plugin v1.1**: Add decision-point flow (upload-only / same-path swap / full integration) ✓ DONE
11. **Plugin v1.2**: Add autonomous mode, safety features, README settings section ✓ DONE
12. **Plugin publish**: Standalone repo at `FlamyFlame/claude-bnl-localgroupdisk`, marketplace PR #61 ✓ DONE (PR pending review)
13. **Plugin marketplace compliance**: Frontmatter, section order, cursor/codex manifests ✓ DONE
14. **Plugin test**: Isolated environment at `/tmp/plugin-test/` — NOT RUN YET
15. **Apply** to remaining datasets (powheg cc truth, pythia truth, powheg fullsim) ← NEARLY DONE
    - 15a: Deleted `bb_evgen_truth_full_sample_orig` (pilot backup) ✓ DONE
    - 15b: CC truth ✓ DONE — dataset `powheg_cc_evgen_truth` (26 files, 108 GB). Old rules (ghost/stuck) cleaned up; fresh rule `9aa0d72bc282413a85bf6658e7e753e4` → OK (26/26 replicated in ~10 min on 2026-06-04). Symlink farm at `~/dcachearea/powheg_full_sample/cc_evgen_truth_full_sample/` rebuilt with all LGD targets. Original → `cc_evgen_truth_full_sample_orig`.
    - 15c: Pythia truth 5TeV ✓ DONE — uploaded to scratchdisk (24 files, ~172 GB, ~1 min/file), dataset `pythia_truth_5TeV`, rule `94b3464b4d3746128a466f5a6c26144f` → OK. Symlink farm at `~/dcachearea/pythia_truth_full_sample/pythia_5TeV/`. Original → `pythia_5TeV_orig`.
    - 15d: Pythia truth 5.36TeV ✓ DONE — uploaded (24 files, ~187 GB), dataset `pythia_truth_5p36TeV`, rule `f5d11de7f63f4bc795d484f36629a906` → OK. Symlink farm at `~/dcachearea/pythia_truth_full_sample/pythia_5p36TeV/`. Original → `pythia_5p36TeV_orig`.
    - 15e: Powheg fullsim bb ✓ DONE — uploaded with clean names (51 files, 119 GB), dataset `powheg_fullsim_bb`, rule `06110deffc8049e88ce24f9d1543689c` → OK. Symlink farm at `~/dcachearea/powheg_full_sample/user.yuhang.TrigRates.dimuon.PowhegPythia.fullsim.pp17.bb.Feb2026.v1._MYSTREAM/`. Symlink names have `user.yuhang.` prefix restored → LGD pnfs paths use stripped names. Original → `*_MYSTREAM_orig`.
    - 15f: Powheg fullsim cc ✓ DONE — uploaded with clean names (52 files, 125 GB), dataset `powheg_fullsim_cc`, rule `b07e82704a564688b18010b4af2a2eea` → OK. Symlink farm at `~/dcachearea/powheg_full_sample/user.yuhang.TrigRates.dimuon.PowhegPythia.fullsim.pp17.cc.Feb2026.v1._MYSTREAM/`. Same name mapping as bb. Original → `*_MYSTREAM_orig`.

### CC Truth Issues — Investigation & Resolution (2026-06-04)

**Original problems (from 2026-06-03):**
- `mc_truth_cc_01.root`: Ghost replica — Rucio catalog showed LGD entry at `9d/ca/` but pnfs path did not exist
- `mc_truth_cc_02.root`: FTS transfer stuck at REPLICATING for >12h; had duplicate file-level rule `6588c2674dbe4947bc540fe6426b8860` from prior session

**Investigation findings:**
- `declare-bad-file-replicas` via Python API → AccessDenied (non-admin can't declare bad replicas)
- `delete_replicas` via Python API → AccessDenied
- `add_bad_pfns` via Python API → AccessDenied
- `update-rule --stuck --boost-rule` → AccessDenied
- `rucio download` of ghost cc_01 from LGD → failed (both davs and root protocols), but did NOT auto-mark suspicious
- `declare_suspicious_file_replicas` via Python API → **succeeded** (only method available to regular users)
- Cannot write directly to LOCALGROUPDISK pnfs (`mkdir` → Permission denied)

**Resolution steps:**
1. Declared cc_01 ghost replica as suspicious via `ReplicaClient.declare_suspicious_file_replicas()`
2. Deleted duplicate file-level rule for cc_02 (`rucio delete-rule 6588c267...`)
3. Deleted old dataset rule `f228a61b` (`rucio delete-rule`) — intended to recreate fresh
4. **LESSON LEARNED**: `rucio update-rule --lifetime 0` means "expire in 0 seconds" (immediate expiration), NOT "set to permanent/no expiration". Both rules were purged by daemon, and all 24 previously-replicated LGD files were also deleted.
5. Ghost replica for cc_01 was cleaned up (suspicious flag + rule deletion caused catalog cleanup)
6. Created fresh rule `9aa0d72bc282413a85bf6658e7e753e4` — clean slate with 0/26/0 (all 26 need replication from SCRATCHDISK)
7. All 26 symlinks in farm temporarily point to `_orig` directory while replication proceeds

**Deduplication result:** Now only 1 rule for CC truth dataset. No file-level rules. Clean.

**After replication completes:** Get PFNs for all 26 files, rebuild symlink farm pointing to LGD paths.

### Fullsim Upload Workaround (for reference)

Files named `user.yuhang.<jobid>.MYSTREAM._NNNNNN.root` cause Rucio server 500 "Database exception" when uploading with `--scope user.yuhang`. Root cause: Rucio server-side fails on POST `/replicas` when DID name starts with the same string as the scope (`user.yuhang.`).

**Workaround**: create symlinks in `/tmp/fullsim_{bb,cc}_upload/` with prefix stripped (e.g., `48591366.MYSTREAM._000001.root`), upload those. Rucio DID names are the stripped versions. Symlink farm maps original names back: `user.yuhang.48591366.MYSTREAM._000001.root` → `/pnfs/.../48591366.MYSTREAM._000001.root`.

### `_orig` directories to delete

All at `~/dcachearea/`:
1. `pythia_truth_full_sample/pythia_5TeV_orig/` (24 files, ~172 GB)
2. `pythia_truth_full_sample/pythia_5p36TeV_orig/` (24 files, ~187 GB)
3. `powheg_full_sample/cc_evgen_truth_full_sample_orig/` (26 files, 101 GB)
4. `powheg_full_sample/user.yuhang.TrigRates.dimuon.PowhegPythia.fullsim.pp17.bb.Feb2026.v1._MYSTREAM_orig/` (51 files, 111 GB)
5. `powheg_full_sample/user.yuhang.TrigRates.dimuon.PowhegPythia.fullsim.pp17.cc.Feb2026.v1._MYSTREAM_orig/` (52 files, 117 GB)

**Do not delete until**: analysis code verified to work with all symlink farms. The `bb_evgen_truth_full_sample_orig` was already deleted (pilot was verified).

## Accumulated Findings

- Step 1: Luke's two-phase procedure: (1) upload to scratchdisk + replicate to LOCALGROUPDISK, (2) access via PFNs. Workflow B confirmed (local re-upload needed). One Rucio dataset per logical file group. LOCALGROUPDISK does not preserve directory structure — files stored in hash-based pnfs subdirs. Symlink farm (one symlink per file) is the most code-compatible access method.
- Steps 2–3: RSE confirmed `BNL-OSG2_LOCALGROUPDISK`. LOCALGROUPDISK pnfs mount at `/pnfs/usatlas.bnl.gov/LOCALGROUPDISK/rucio/user/<username>/<hash>/`. No DID conflicts for test filenames. Upload does not modify source files — no backup needed. Analysis code path: `powheg_dcache_dir + data_subdir + filename` (PowhegAlgCoreT.c:28,39,95).
- Pilot P1: All 25 files (98 GB) uploaded to `BNL-OSG2_SCRATCHDISK` successfully (~30s per file). Time: 19:47–20:28 on 2026-05-27.
- Pilot P2: Dataset `user.yuhang:powheg_bb_evgen_truth_pilot` created, all 25 files attached.
- Pilot P3 **BLOCKED then RESOLVED**: Initially failed with "insufficient quota" — `yuhang` not in `/atlas/usatlas` VOMS group. After VOMS membership granted (2026-05-28), `rucio list-account-limits yuhang` shows `BNL-OSG2_LOCALGROUPDISK = 50.000 TB`. Rule added successfully: `rucio add-rule user.yuhang:powheg_bb_evgen_truth_pilot 1 BNL-OSG2_LOCALGROUPDISK` → rule ID `86adb3d14cfe489ca6311c0ce587a123`. Replication took ~9.5 hours (queued ~9h, transferred in ~5 min once FTS picked it up). State: OK (25/25) at 20:42 UTC 2026-05-28.
- Pilot P4: PFNs retrieved for all 25 files. Each PFN: `root://dcgftp.usatlas.bnl.gov:1094//pnfs/usatlas.bnl.gov/LOCALGROUPDISK/rucio/user/yuhang/<hash>/<filename>.root`. Saved to `BNL_LOCALGROUPDISK_backup/pfns_bb_truth.txt`. pnfs paths confirmed directly readable on BNL SDCC nodes (no grid proxy needed).
- Pilot P5: Symlink farm built at `~/dcachearea/powheg_full_sample/bb_evgen_truth_full_sample/` (25 symlinks). Original directory renamed to `bb_evgen_truth_full_sample_orig`. Verified: file size and TTree entries (200,000 per file) match exactly between LGD symlinks and original. TChain test with 5 files (file_batch=1): 1,000,000 entries, all OK.
- Pilot P6: **SUCCESS**. Full Powheg truth analysis (`PowhegTruthAnalysis pw(1, "bb", 0)`) ran to completion on LOCALGROUPDISK symlinks. 1,000,000 events processed in 625s. Output: `muon_pairs_powheg_bb_truth_part1.root` (131 MB), `hists_powheg_ntuple_processing_powheg_bb_truth_part1.root` (948 KB). No code changes needed — symlink farm is fully transparent.
- Plugin: Created `bnl-localgroupdisk` Claude Code plugin at `BNL_LOCALGROUPDISK_backup/bnl-localgroupdisk-plugin/` with 4 skills (`preflight`, `migrate`, `check-rule`, `build-symlinks`). Encodes full procedure, troubleshooting, and key facts from pilot. Publicly shareable.
- Plugin v1.1: Restructured `migrate` skill with 3 decision points: (1) upload only vs. upload + symlink swap (default), (2) same-path swap (default, no code changes) vs. different path, (3) symlink only vs. full integration (codebase scan, path updates, test, git-branch rollback). `farm_dir` dropped from required argument to conditional. `disable-model-invocation` set to `false`, `Edit` added to `allowed-tools` for full integration path. User is asked for explicit approval before code modification begins.
- Plugin v1.2 (2026-06-03): Added autonomous mode (recognizes "no confirmation" etc.), non-.root file warning, TTree smoke test for same-path swap, `_orig` cleanup prompt, 96h VOMS validity note, BNL storage docs link. Published standalone repo at `github.com:FlamyFlame/claude-bnl-localgroupdisk` (4 commits on main). README updated with recommended settings, quick start with autonomous mode, worked example.
- Marketplace PR #61 (2026-06-03): Submitted to `usatlas/marketplace` from fork `FlamyFlame/marketplace` branch `add-bnl-localgroupdisk-plugin`. CodeRabbit review required compliance fixes: frontmatter (only `name`+`description` allowed, "Use when..." prefix), canonical section order (Overview→When to Use→Key Concepts→Canonical Patterns→Gotchas→Interop→Docs), `## Docs` section, `.cursor-plugin/plugin.json`, `marketplace.json` entry, `.codex/INSTALL.md` update. All fixes committed and pushed (`8958f34`). Standalone repo back-ported content improvements (Gotchas, Interop, When to Use, Key Concepts, Docs) while keeping standalone-specific frontmatter (`arguments`, `argument-hint`, `allowed-tools`) — commit `ecdf841`.
- Plugin test environment: set up at `/tmp/plugin-test/` with isolated `.claude/settings.json` and minimal CLAUDE.md (no project tracking docs). Plugin copied to `.claude/plugins/`. **Test was never executed** — no upload_scratch.log or test-log.md generated. The environment is ready for a future test run.
- Step 15 batch migration (2026-06-03): Uploaded 4 new datasets to scratchdisk in parallel (~35 min total). CC truth was already on scratchdisk and had a dataset from a prior session — most files (24/26) already on LGD. Pythia truth 5TeV and 5.36TeV replicated almost instantly after upload (rules went OK within minutes — much faster than pilot's 9h queue). Fullsim bb/cc also replicated within ~30 min. Upload rate: ~1 min/file for 7–15 GB pythia files, ~20s/file for 2 GB fullsim files. All 4 new symlink farms verified (pnfs paths resolve, `ls -lL` returns real file sizes). Fullsim naming workaround required (see above). CC truth has 2 problem files: cc_01 ghost replica, cc_02 still replicating.
- **PFN files saved**: `/tmp/pfns_pythia_5TeV.txt`, `/tmp/pfns_pythia_5p36TeV.txt`, `/tmp/pfns_cc_truth.txt`, `/tmp/pfns_fullsim_bb.txt`, `/tmp/pfns_fullsim_cc.txt`. Also temp upload dirs at `/tmp/fullsim_{bb,cc}_upload/` (symlinks with stripped names).
- CC truth fix (2026-06-04): Ghost replica (cc_01) and stuck transfer (cc_02) resolved by: (1) `declare_suspicious_file_replicas` on ghost, (2) deleting all old rules, (3) creating fresh rule `9aa0d72bc282413a85bf6658e7e753e4`. **Critical lesson**: `rucio update-rule --lifetime 0` = expire immediately, NOT "set permanent". Old rules were purged and all 24 LGD files were deleted. Fortunately all 26 source files intact on scratchdisk. Fresh rule triggered clean re-replication. Regular users cannot: `declare-bad-file-replicas`, `delete_replicas`, `add_bad_pfns`, `update-rule --stuck/--boost` on LOCALGROUPDISK. Only `declare_suspicious_file_replicas` is available.

## Ruled Out

- Workflow A (add-rule only): original grid scratchdisk replicas have expired
- Single symlink swap of `~/dcachearea`: LOCALGROUPDISK doesn't preserve directory structure
- Backup before upload: `rucio upload` copies, does not delete source
- `rucio update-rule --lifetime 0` for undoing deletion: causes immediate expiration + file purge. Don't use.
- `rucio declare-bad-file-replicas` / `delete_replicas` / `add_bad_pfns` / `update-rule --stuck/--boost`: all require admin, not available to regular users
- `rucio download` of a ghost replica does NOT auto-mark it suspicious on server side

## Latest Stage

**Step 15 complete (2026-06-04).** All 6 datasets migrated with same-path-swap symlink farms pointing to LOCALGROUPDISK. CC truth fixed: ghost/stuck issues resolved by deleting old rules + creating fresh rule → 26/26 replicated in ~10 min.

### Current symlink farm status (all under `~/dcachearea/`)

| Directory | Symlinks | Status |
|---|---|---|
| `powheg_full_sample/bb_evgen_truth_full_sample/` | 25/25 | ✓ Complete (LGD targets) |
| `powheg_full_sample/cc_evgen_truth_full_sample/` | 26/26 | ✓ Complete (LGD targets, rebuilt 2026-06-04) |
| `pythia_truth_full_sample/pythia_5TeV/` | 24/24 | ✓ Complete (LGD targets) |
| `pythia_truth_full_sample/pythia_5p36TeV/` | 24/24 | ✓ Complete (LGD targets) |
| `powheg_full_sample/user.yuhang.TrigRates...bb..._MYSTREAM/` | 51/51 | ✓ Complete (LGD targets) |
| `powheg_full_sample/user.yuhang.TrigRates...cc..._MYSTREAM/` | 52/52 | ✓ Complete (LGD targets) |

### TO-DO

1. **Verify analysis code** works with all symlink farms (run quick tests for each dataset type)
2. **Delete `_orig` directories** (5 dirs, ~688 GB) after verification
3. **`dimuon_data` migration** not yet started (594 GB, 3rd top-level dir under `~/dcachearea/`) — separate batch
4. **Marketplace PR #61** — still pending CodeRabbit re-review after compliance fixes
5. **Plugin test** at `/tmp/plugin-test/` — never executed
