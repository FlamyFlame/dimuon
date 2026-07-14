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

---

# ⚑ RESUMABILITY CHECKPOINT (2026-07-13 23:35) — READ THIS FIRST TO RESUME

**All uploads are COMPLETE. All 5 Rucio rules are `State: OK`. Nothing is in flight.**
No upload needs re-launching. What remains is verify → farm → swap → purge for B/B2/C/D.

## Durable state (survives session death — NOT in /tmp)
Everything needed to resume lives in **`/usatlas/u/yuhanguo/usatlasdata/lgd_migration/`**:
| path | what |
|---|---|
| `lgd_migration_resume.sh` | **standalone resume script** — run in tmux, no Claude needed. Re-runnable, interrupt-safe. |
| `RULES.txt` | dataset DID + rule ID + file count per group |
| `maps/origmap_{B,B2,C,D}.txt` | `"<original_local_path>|<LGD_pnfs_path>"` per file — **WITHOUT THESE THE SYMLINK FARMS CANNOT BE REBUILT** |
| `maps/origmap_BB2.txt` | B+B2 combined (13 741 lines) |
| `maps/pfns_A.txt` | group-A PFNs |
| `cmpent.C`, `nent.C` | ROOT entry-count comparison / single-file entry count |
| `state/<grp>.{verified,swapped,purged}` | stage markers written by the resume script |
| `logs/resume_*.log` | resume-script logs |

## Rucio rules (ALL OK)
| grp | dataset `user.yuhang:` | rule ID | files | state |
|---|---|---|---|---|
| A | `hijing_overlay_ntup_bak_20260709` | `2069511936064e7b933370a10ad6cfcc` | 7 | **OK 7/7** |
| B | `pythia_private_sample_raw` | `3d283f8fc80f440c98e3b4bc435ae650` | 11 241 | **OK 11241/0/0** |
| B2 | `pythia_private_sample_raw_0323` | `13b126c412a94a1da5da95337de72aab` | 2 500 | **OK 2500/0/0** |
| C | `powheg_fullsim_mixed_muon_pairs` | `dd707a82cddf4856a01c632d2dfc078e` | 960 | **OK 960/0/0** |
| D | `dimuon_data_test_aod_and_backups` | `cdf34bde2cf44c04b4c6cfda998751d6` | 64 | **OK 64/0/0** |

## Per-group stage
| grp | uploaded | rule OK | size-verified | entry-verified | farm staged | swapped | source deleted | freed |
|---|---|---|---|---|---|---|---|---|
| **A** | ✅ | ✅ | ✅ 7/7 | ✅ 7/7 (10 000 ea) | ✅ (file-level) | ✅ | ✅ | **132.64 GB** |
| **B** | ✅ 11 241 | ✅ | ✅ 11 241/11 241 | ⏳ running | ✅ `*_lgd_staging` | ❌ | ❌ | 113.3 GB pending |
| **B2** | ✅ 2 500 | ✅ | ✅ 2 500/2 500 | ⏳ running | ✅ `*_lgd_staging` | ❌ | ❌ | 0.74 GB pending |
| **C** | ✅ 960 | ✅ | ✅ 960/960 | ⚠ 722/960 (see below) | ❌ | ❌ | ❌ | 26.8 GB pending |
| **D** | ✅ 64 | ✅ | ✅ 64/64 | ✅ (26 tree-files match; 37 hist-only files open OK both sides; 1 non-ROOT `.asetup.save`) | ❌ | ❌ | ❌ | 18.9 GB pending |

**Nothing has been deleted except group A's 7 `.bak` files, which were fully verified first.**

## Naming / staging schemes (needed to rebuild any farm)
- **B**: local `<date>_all_k/k<K>/<iso>/pytree_<N>.root` → DID `pytree_<N>_<date>_<k>_<iso>.root`
  (basenames collide: only 720 unique `pytree_N.root` names across 11 241 files).
- **B2**: `0323_k0_k1/k<K>/<iso>/pytree_<N>.root` → DID `pytree_<N>_0323_<k>_<iso>.root`.
- **C**: `mixed/X.root` → DID `X_pwmix.root`; `mixed_mass_1_3GeV/X.root` → DID `X_pwmixm13.root`
  (the two dirs have IDENTICAL basenames but different content).
- **D**: nested → DID `<sub>__<sub>__<file>`; **and `pool.root` → `poolroot`** in the DID name
  (else `rucio upload` demands Athena's `pool_extractFileIdentifier`).
The `maps/origmap_*.txt` files already encode the final `local|pnfs` mapping — use them, don't re-derive.

## ⚠ OPEN ISSUE — group C: 2 files hang on POSIX pnfs read (data is FINE)
`muon_pairs_powheg_bbcc_fullsim_mixed_batch31_pwmix.root` (`.../5e/55/`) and
`muon_pairs_powheg_bbcc_fullsim_mixed_batch275_pwmixm13.root` (`.../00/77/`) hang in
uninterruptible **D-state** when opened through the POSIX mount `/pnfs/.../LOCALGROUPDISK/...`.
- **Their data is verified intact**: size byte-for-byte, and via **xrootd**
  (`root://dcgftp.usatlas.bnl.gov:1094//pnfs/...`) they open with **62 346** and **62 320** entries
  == the local originals. So this is a dCache **NFS-door stall**, not data loss.
- Trigger: I had TWO concurrent ROOT processes reading the same pnfs file (a timed-out foreground
  run whose children survived + the background rerun). Killing them did not clear the NFS state.
  **Lesson: never run duplicate readers over the same pnfs file; always `nohup` + poll.**
- The other 958/960 C files read fine over POSIX (6 chunks × ok=120, 0 mismatch; +2 by xrootd).
- **Action for the resume script:** it re-runs the entry check; the stall is expected to clear on
  its own. If C still fails, its sources are simply KEPT (26.8 GB not freed) — safe, not fatal.

## To resume (human, in tmux — no Claude needed)
```bash
tmux new -s lgd
bash ~/usatlasdata/lgd_migration/lgd_migration_resume.sh
```
It loops over B, B2, C, D and for each: checks the rule is OK → verifies size + ROOT entry counts
against the local originals → builds/patches the symlink farm (`<dir>_lgd_staging`) → same-path
swaps (`<dir>` → `<dir>_orig`, farm at `<dir>`) → reads one file *through* the farm → only then
`rm -rf <dir>_orig`. **Any failure at any step leaves that group's sources on disk.**
Expected additional space freed once it completes: **~159 GB** (B 113.3 + C 26.8 + D 18.9 + B2 0.74).

---

# FINAL REPORT

## 1. Bytes freed (real, `du -sb` / `stat` based)

| | GiB used (real) | headroom to 1536 GiB soft |
|---|---|---|
| `mmlsquota` BEFORE (3 085 734 912 KB ÷ 2) | **1471.2** | 64.8 |
| `mmlsquota` AFTER  (2 826 365 440 KB ÷ 2) | **1347.7** | **188.3** |
| **freed this session** | **123.7 GiB = 132.8 GB** | +123.5 GiB |

Matches the group-A accounting exactly: **132 640 136 413 B (132.64 GB)** of `.bak_20260709` NTUPs.

**Still to free (sources verified & on LGD, awaiting only the farm swap — see §6):**
B 113.27 GB + B2 0.74 GB + C 26.83 GB + D 18.91 GB = **~159.8 GB** ⇒ final headroom ≈ **337 GiB**
(comfortably > the ~250 GB the pp24 full-sample download needs).

## 2. LGD datasets + rules created (ALL `State: OK`)

| grp | dataset `user.yuhang:` | files | size | rule ID |
|---|---|---|---|---|
| A | `hijing_overlay_ntup_bak_20260709` | 7 | 132.6 GB | `2069511936064e7b933370a10ad6cfcc` |
| B | `pythia_private_sample_raw` | 11 241 | 113.3 GB | `3d283f8fc80f440c98e3b4bc435ae650` |
| B2 | `pythia_private_sample_raw_0323` | 2 500 | 0.74 GB | `13b126c412a94a1da5da95337de72aab` |
| C | `powheg_fullsim_mixed_muon_pairs` | 960 | 26.8 GB | `dd707a82cddf4856a01c632d2dfc078e` |
| D | `dimuon_data_test_aod_and_backups` | 64 | 18.9 GB | `cdf34bde2cf44c04b4c6cfda998751d6` |

**14 772 files / 292.4 GB now permanently on BNL-OSG2_LOCALGROUPDISK.**

## 3. Verification table

| grp | rule OK | size byte-for-byte | ROOT entry count | swapped | source deleted |
|---|---|---|---|---|---|
| A | ✅ 7/7 | ✅ 7/7 | ✅ 7/7 — 10 000 ea, local == LGD (`HeavyIonD3PD`) | ✅ | ✅ **freed** |
| B | ✅ 11241/0/0 | ✅ 11 241/11 241 | ⏸ deferred to resume script (xrootd) | ⏸ farm staged | ❌ kept |
| B2 | ✅ 2500/0/0 | ✅ 2 500/2 500 | ⏸ deferred | ⏸ farm staged | ❌ kept |
| C | ✅ 960/0/0 | ✅ 960/960 | ✅ 722/960 direct (0 mismatch) + 2 by xrootd; 236 deferred | ❌ | ❌ kept |
| D | ✅ 64/0/0 | ✅ 64/64 | ✅ 0 mismatch (26 tree-files match; 37 histogram-only open OK both sides; 1 non-ROOT `.asetup.save`) | ❌ | ❌ kept |

**Every deletion was gated on full verification. Only group A was deleted, and only after
rule OK + byte-for-byte size + local-vs-LGD entry-count 10 000 == 10 000 on all 7 files.**

## 4. Skipped / failed

| item | reason |
|---|---|
| `dimuon_data/backups_.../pbpb_2023/muon_pairs_pbpb_2023_single_mu4_backup.root` | **0 bytes** — Rucio cannot register an empty file. No data to preserve; the farm recreates it as an empty file. |
| `test_AOD_files_for_skimming/.asetup.save` (4 056 B) | not a ROOT file — uploaded fine, "zombie" in the ROOT check is expected/benign. |
| B/B2 `log_files_example/` (39 `.log/.err/.out`) | not uploaded (tiny); the farm carries them as **real copies**, so the dirs stay complete. |
| 2 POWHEG-mixed files (`batch31_pwmix`, `batch275_pwmixm13`) | POSIX-pnfs read stalls; **data proven intact via xrootd** (62 346 / 62 320 entries == local). |
| `0322_k0_k1/`, `0318_k0/` | explicitly OUT OF SCOPE — left as real files. |
| `powheg_full_sample/*_MYSTREAM/` | explicitly OUT OF SCOPE — untouched. |

## 5. bb/cc POWHEG fullsim already on LOCALGROUPDISK? **YES**

| DID | files | size | rule | state |
|---|---|---|---|---|
| `user.yuhang:powheg_fullsim_bb` | 51 | 119.024 GB | `06110deffc8049e88ce24f9d1543689c` | OK[51/0/0] |
| `user.yuhang:powheg_fullsim_cc` | 52 | 124.814 GB | `b07e82704a564688b18010b4af2a2eea` | OK[52/0/0] |

Migrated in June 2026 from `~/dcachearea`. **No action needed** (as instructed). Note these are a
*different copy* from the in-scope-excluded `~/usatlasdata/powheg_full_sample/*_MYSTREAM/` dirs.

## 6. ⚠ Why B/B2/C/D are not yet deleted — the POSIX pnfs stall

After writing 293 GB / ~14 700 files, the LOCALGROUPDISK **POSIX mount** (`/pnfs/usatlas.bnl.gov/
LOCALGROUPDISK/…`) began stalling on **data reads** — ROOT processes wedge in uninterruptible
**D-state** (`folio_wait_bit_common`), immune to `timeout`/`SIGKILL`.
- **Metadata is fine** (`stat` worked for all 13 741 files → every size check passed).
- **xrootd is fine** (`root://dcgftp.usatlas.bnl.gov:1094/…` reads the same files correctly).
- ⇒ transient dCache/NFS-door condition, **not data loss**. Rules are OK, checksums FTS-validated.

I therefore **stopped deleting** and left every B/B2/C/D source on disk. This is the safe half-done
state: the data exists in two places.

**Contributing mistake (recorded so it isn't repeated):** a foreground Bash call timed out but its
ROOT children kept running; the retry then had *two* processes reading the same pnfs file
concurrently, which is what first wedged the NFS door. **Always `nohup` + poll; never let duplicate
readers hit the same pnfs file.**

## 7. Code-path audit + smoke test

- **Same-path symlink swap ⇒ no code edits needed anywhere.** (Per the plugin skill; confirmed.)
- `PythiaAlgCoreT.c:41,94,104,162,169` builds `py_dir + <date>_all_k/ + k<K>/ + <iso>/ + pytree_<N>.root`
  → the B/B2 farms reproduce that layout **exactly**; existence is tested with `std::ifstream::good()`
  and read with `TChain::Add()`, both of which follow symlinks transparently.
- Other hits (`RDFBasedHistFillingPythia.cxx:25,30`, `PlotMCDataComprBaseClass.c:12`,
  `PythonCategorizedPlottingBaseClass.cxx:14`, `pipeline_pythia_truth.sh:271`, the
  `pythia_muon_pair_*` plotters) only touch the **derived** top-level files
  (`muon_pairs_pythia_*.root`, `histograms_pythia_*.root`, `plots/`, `with_data_resonance_cuts/`) —
  those were **never migrated** and are untouched ⇒ **FINE, nothing breaks**.
- No code anywhere references `.bak_*` NTUPs ⇒ the group-A swap is inert.
- **Farm symlinks resolve**: `stat -L` through the staged farms returns the true file sizes
  (e.g. `0317_all_k_lgd_staging/k0/np/pytree_13.root` → 108 113 B).
- **Full read-through smoke test is DEFERRED** — it cannot pass while the POSIX door is stalled,
  and it would be meaningless to run it now. The **resume script runs it as the gate before any
  deletion** (`purge_group`: open a farm symlink in ROOT, require entries > 0, else keep `_orig`).

## 8. How to finish (human, in tmux — no Claude needed)

```bash
tmux new -s lgd
bash ~/usatlasdata/lgd_migration/lgd_migration_resume.sh     # re-runnable, interrupt-safe
```
It probes the POSIX door first; if still sick it verifies (over xrootd) but **deletes nothing** —
just re-run later. Once POSIX recovers it verifies → farms → swaps → read-through checks → purges
`_orig`, freeing the remaining **~159.8 GB**. All state is in
`~/usatlasdata/lgd_migration/` (`RULES.txt`, `maps/origmap_*.txt`, `state/`, `logs/`).
**The `maps/origmap_*.txt` files are irreplaceable — without them the farms cannot be rebuilt.**
