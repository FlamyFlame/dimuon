# SUB: LOCALGROUPDISK migration #1 (scratch doc — append-only)

Subagent scratch doc for the ~293 GB GPFS→`BNL-OSG2_LOCALGROUPDISK` MOVE migration.
Orchestrator owns git + canonical docs. **This subagent NEVER runs git.**

## Autonomy Contract
- Mandate: run autonomously to DONE; do NOT pause to confirm progress.
- Done = groups A–D uploaded to LGD (via SCRATCHDISK + replication rule), every file verified
  (replica on LGD, size match, .root entry-count match), local originals deleted, quota
  before/after recorded, bb/cc-POWHEG-fullsim-on-LGD question answered, code-path grep +
  smoke test done, FINAL REPORT section written.
- Stop-and-ask = ANY physics-results-bending ambiguity (no fixed list; use judgment;
  when unsure whether an ambiguity is blocking, treat it as blocking → AskUserQuestion).

## Recipe (from CLOSED `localgroupdisk_migration.md`)
1. `source /usatlas/u/yuhanguo/setup.sh; lsetup rucio`
2. `rucio upload --rse BNL-OSG2_SCRATCHDISK --scope user.yuhang <files>`
   - **Gotcha:** DID names starting with `user.yuhang.` → server 500. Workaround: upload via
     symlinks with the prefix stripped.
3. `rucio add-dataset user.yuhang:<ds>` + `rucio attach user.yuhang:<ds> user.yuhang:<file>...`
4. `rucio add-rule user.yuhang:<ds> 1 BNL-OSG2_LOCALGROUPDISK` → wait state OK
5. `rucio list-file-replicas ... --pfns --rses BNL-OSG2_LOCALGROUPDISK`
6. Verify (size + TTree entries via the LGD pnfs path, directly readable on SDCC) → delete local.
   **NEVER** `rucio update-rule --lifetime 0` (= expire immediately, purges files).

## Step 0 — quota BEFORE (2026-07-13)
`mmlsquota atlasgpfs01`, fileset `usatlast3-data` (KB, allocated = 2× real):
- used 3085734912 KB → **1471.2 GiB real**
- soft quota 3221225472 KB → 1536 GiB ; hard 3758096384 KB → 1792 GiB
- headroom to soft: **64.8 GiB**

## Step 1 — inventory (in progress)

### COURSE CORRECTION (mid-task): follow the `bnl-localgroupdisk` plugin skill
Read `.claude/plugins/bnl-localgroupdisk/skills/{preflight,migrate}/SKILL.md`.
Procedure = upload to SCRATCHDISK → add-rule to LGD → **same-path symlink swap**
(`<dir>` → `<dir>_orig`, farm at original path) → verify → delete `_orig`.
Deleting `_orig` is what frees the quota. NO code edits needed.

## Step 1 — PREFLIGHT (2026-07-13, all PASS)
| Check | Status | Details |
|---|---|---|
| Rucio account | OK | `yuhang` → scope `user.yuhang` |
| Grid proxy | OK | 94:59 left; VOMS `/atlas/Role=NULL` **and `/atlas/usatlas/Role=NULL`** present |
| RSE names | OK | `BNL-OSG2_LOCALGROUPDISK`, `BNL-OSG2_SCRATCHDISK` |
| LGD quota | OK | limit 50 TB, used 841.469 GB, free 49.159 TB |
| SCRATCHDISK | OK | limit 50 TB, used 0 B |
| pnfs mount | OK | `/pnfs/usatlas.bnl.gov/LOCALGROUPDISK/rucio/user/yuhang/` exists |

DID-conflict probes (all **FREE**, no collision): `Pythia_..._pTH8_14...NTUP.bak_20260709.root`,
`pytree_10.root`, `muon_pairs_powheg_bbcc_fullsim_mixed_batch100.root`, `AOD.41716150._000001.pool.root.1`.

## Step 2 — C-QUESTION (read-only, ANSWERED): bb/cc POWHEG fullsim on LGD?
**YES — both are already on BNL-OSG2_LOCALGROUPDISK, rules OK:**
| DID | files | size | rule | state |
|---|---|---|---|---|
| `user.yuhang:powheg_fullsim_bb` | 51 | 119.024 GB | `06110deffc8049e88ce24f9d1543689c` | OK[51/0/0] |
| `user.yuhang:powheg_fullsim_cc` | 52 | 124.814 GB | `b07e82704a564688b18010b4af2a2eea` | OK[52/0/0] |
Also on LGD from June: `powheg_bb_evgen_truth_pilot` (25, 105.0 GB, OK), `powheg_cc_evgen_truth`
(26, 108.1 GB, OK), `pythia_truth_5TeV` (24, 189.7 GB, OK), `pythia_truth_5p36TeV` (24, 194.8 GB, OK).
These correspond to the `~/dcachearea/powheg_full_sample/*_MYSTREAM/` symlink farms — a DIFFERENT
location from the in-scope `~/usatlasdata/powheg_full_sample/*_MYSTREAM/` (OUT OF SCOPE, untouched).

## Step 3 — Group A pre-check: canonical counterparts (BLOCKING GATE) — PASS
All 7 `.bak_20260709` files have a canonical replacement on disk with **10 000 entries** in
`HeavyIonD3PD` (verified by re-opening each canonical NTUP in ROOT):
pTH125_300, pTH14_24, pTH24_40, pTH40_70, pTH70_125, pTH8_14 (r17618, 22.1 GB each) + pTH8_14_r17662.
⇒ all 7 `.bak` are safe to migrate. **No `.bak_20260710.root` exists** (checked). The r17662
`.bak_20260709` (171.9 MB) IS superseded → included.

## Step 4 — INVENTORY (real bytes, `du -sb` / `stat -c %s`)
| grp | source | files | bytes |
|---|---|---|---|
| A | `pythia_fullsim_hijing_overlay_test_sample/*.bak_20260709.root` | 7 | 132 640 132 413 (132.6 GB) |
| B | `pythia_private_sample/{0317,0318,0325,0401,0429}_all_k/` | 11 241 `.root` (+39 log/err/out) | 113 273 015 135 (113.3 GB) |
| C | `powheg_full_sample/{mixed, mixed_mass_1_3GeV}/` | 960 | 26 825 195 979 (26.8 GB) |
| D | `dimuon_data/{test_aod, test_AOD_files_for_skimming, backups_before_test_20260222_152110}/` | 65 | 18 913 567 029 (18.9 GB) |
| | **TOTAL** | **12 273** | **291.65 GB** |

### Naming hazards found (would have caused DID collisions — handled)
- **B: basenames COLLIDE.** Only **720 unique** `pytree_N.root` names across 11 241 files (repeated
  in every `<date>_all_k/k<K>/<iso>/`). Rucio DID names are flat within a scope ⇒ uploaded via
  staging symlinks renamed `pytree_<N>_<date>_<k>_<iso>.root` (11 241 unique, verified 0 dupes).
- **C: `mixed/` and `mixed_mass_1_3GeV/` have IDENTICAL basenames** (`muon_pairs_powheg_bbcc_
  fullsim_mixed_batch*.root`, 480 each) but DIFFERENT content (batch100: 28 820 369 vs 27 090 835 B)
  ⇒ suffixed `_pwmix` / `_pwmixm13` for upload.
- **D: nested** ⇒ flattened to `<subdir>__<subdir>__<file>` for upload.
- **D: one 0-byte file** `backups_.../pbpb_2023/muon_pairs_pbpb_2023_single_mu4_backup.root`
  → SKIPPED from upload (rucio cannot register an empty file; 0 bytes = no data to preserve).
- **B `log_files_example/` (39 .log/.err/.out, <1 MB)** → not uploaded; will be recreated in the
  symlink farm as real files (copied), so the dir stays complete.

### Code dependency (drives the farm layout) — `PythiaAlgCoreT.c`
`Analysis/NTupleProcessingCode/PythiaAlgCoreT.c:41` `py_dir = ".../pythia_private_sample/"`;
`:94,:104` `job_dirs = {"0317_all_k/","0318_all_k/","0318_k0/"}` / `{"0322_k0_k1/","0323_k0_k1/",
"0325_all_k/","0401_all_k/","0429_all_k/"}`; `:162` `job_path = py_dir + job_dir + kin_dir +
beam_dir`; `:169` `fpath = job_path + "pytree_" + N + ".root"`; `:170` existence test via
`std::ifstream::good()` (symlinks OK); `:177` `evChain->Add(fpath + "?#PyTree")`.
⇒ **the B symlink farm MUST reproduce the exact nested layout** `<date>_all_k/k<K>/<iso>/pytree_<N>.root`.
⇒ `0318_k0/`, `0322_k0_k1/`, `0323_k0_k1/` are ALSO in `job_dirs` but are NOT in the approved
migration list → they stay on disk untouched (small).

## Step 5 — upload benchmark
50 small B files → 114 s = **2.3 s/file** (dominated by per-file registration, not bandwidth).
11 241 files single-stream ≈ 7.1 h ⇒ **run 6 parallel workers** (`upload_worker.sh`, chunks of 200).

## Step 6 — UPLOADS LAUNCHED (2026-07-13 ~19:30)
- A: single `rucio upload` of the 7 `.bak` files (154 GB read + upload), background.
- B: 6 parallel workers over 11 191 remaining files (50 already done by the benchmark).
- C: 1 worker, 960 files. D: 1 worker, 64 files.

## Step 7 — UPLOAD ISSUE FOUND & FIXED (group D)
`rucio upload` aborted the whole D chunk (rc=1):
> `ERROR An unknown exception occurred. Details: Trying to upload ROOT files but
> pool_extractFileIdentifier tool can not be found. Setup your ATHENA environment and try again.`

**Root cause:** the rucio client calls `pool_extractFileIdentifier` (to extract the POOL GUID) for
any file whose **basename contains `pool.root`** — the two test AODs
`test_aod/data23_hi/AOD.41716150._00000{1,6}.pool.root.1`. Athena is not set up in the rucio env.
**Fix (no Athena needed):** rename the *staging* DID name `pool.root` → `poolroot`
(`...AOD.41716150._000001.poolroot.1`). Only the Rucio DID name changes; the symlink farm restores
the original filename, so nothing downstream sees it. Re-ran → 64/64 OK.
(Trade-off: no POOL GUID registered for those 2 AODs — irrelevant for archive+path-read use.)

## Step 8 — DATASETS + REPLICATION RULES CREATED
| grp | Rucio DID (`user.yuhang:`) | files | rule ID | state |
|---|---|---|---|---|
| A | `hijing_overlay_ntup_bak_20260709` | 7 | `2069511936064e7b933370a10ad6cfcc` | **OK 7/7** |
| B | `pythia_private_sample_raw` | 11 241 | `3d283f8fc80f440c98e3b4bc435ae650` | replicating |
| B2 | `pythia_private_sample_raw_0323` | 2 500 | (pending) | — |
| C | `powheg_fullsim_mixed_muon_pairs` | 960 | `dd707a82cddf4856a01c632d2dfc078e` | **OK 960/960** |
| D | `dimuon_data_test_aod_and_backups` | 64 | `cdf34bde2cf44c04b4c6cfda998751d6` | **OK 64/64** |

Gotcha hit: `rucio upload` log contains **ANSI colour escapes** → naive
`grep -oE "Successfully uploaded file [^ ]+"` yields DID names with a trailing `\x1b[0m`, and
`rucio attach` then silently attaches 0 files. Must `sed -r 's/\x1b\[[0-9;]*m//g'` first.
(Caught it: "files in dataset: 0" after the first attach.)

## Step 9 — SCOPE ADDITION (user-approved mid-task): `0323_k0_k1/`
2 500 `.root`, 744 275 866 B (0.74 GB). Same basename collision (**720 unique** names) → same
renaming scheme `pytree_<N>_0323_<k>_<iso>.root`; 2 500 staging symlinks, **0 dupes**. Uploaded via
6 parallel workers. Layout `0323_k0_k1/k{0,1}/{pp,pn,np,nn}/` must be reproduced in the farm
(`PythiaAlgCoreT.c:104` lists it in `job_dirs`).
**Still OUT OF SCOPE, left as real files:** `0322_k0_k1/` (0.75 GB), `0318_k0/` (0.11 GB).
⇒ `pythia_private_sample/` ends up a MIX of symlink-farm dirs and real dirs — expected.

## Step 10 — GROUP A VERIFIED + SWAPPED — **132 640 136 413 B (132.64 GB) FREED**
Per-file verification (all 7 PASS):
| file (`Pythia_5p36TeV_pp_hQCD_DiMu_*.FullSimHIJINGOverlayPP24*.NTUP.bak_20260709.root`) | local B | LGD B | size | local ent | LGD ent |
|---|---|---|---|---|---|
| pTH125_300 | 22 111 495 840 | 22 111 495 840 | OK | 10 000 | 10 000 |
| pTH14_24 | 22 059 597 395 | 22 059 597 395 | OK | 10 000 | 10 000 |
| pTH24_40 | 22 067 568 388 | 22 067 568 388 | OK | 10 000 | 10 000 |
| pTH40_70 | 22 081 267 337 | 22 081 267 337 | OK | 10 000 | 10 000 |
| pTH70_125 | 22 094 315 530 | 22 094 315 530 | OK | 10 000 | 10 000 |
| pTH8_14 | 22 053 987 366 | 22 053 987 366 | OK | 10 000 | 10 000 |
| pTH8_14_r17662 | 171 904 557 | 171 904 557 | OK | 10 000 | 10 000 |
(tree `HeavyIonD3PD`; LGD copies opened directly from `/pnfs/.../LOCALGROUPDISK/rucio/user/yuhang/<hh>/<hh>/`)
Swap: each local file replaced **in place** by a symlink to its LGD pnfs path (file-level same-path
swap — these are loose files in a shared dir, not a whole dir). Post-swap: `ls -lL` returns the real
sizes and ROOT reads **10 000** entries *through* the symlink. ✔

## Step 11 — GROUPS C & D VERIFIED (size)
- C: **960/960 SIZE_OK**, 0 mismatch. D: **64/64 SIZE_OK**, 0 mismatch (byte-for-byte vs LGD pnfs).
- Entry-count (local vs LGD) comparison running.
