# SkimCode quick README

## What this code does

This package runs the `TrigRates` Athena algorithm to produce a skim ROOT ntuple (`myfile.root`) for heavy-ion / pp dimuon studies.

Main features in `TrigRates`:

- Event cleaning (GRL, detector-error checks, vertex cut)
- Trigger decisions and prescales (including L1 info)
- Muon trigger matching (single- and dimuon paths)
- Muon MCP tools:
	- selection
	- momentum calibration
	- efficiency corrections (Medium + Tight)
	- trigger matching
- Optional storage of tracks, vertex, MET, truth, and MC event info
- ZDC information (energy, time, status, per-module PreSampleAmp; RPD centroid optional)


## Config scripts (JO vs CA)

Under `scripts/`:

- `TrigRates_JO.py`:
	- Legacy job options configuration
	- Use for **R21 / pre-R24** workflows
- `TrigRates_CA.py`:
	- ComponentAccumulator-based configuration
	- Use for **R24/R25 (recommended for modern setup)**
	- **This is the canonical master copy** — always edit here, then sync to run directories (see below)
- `TrigRates.py`:
	- Original legacy script (kept for compatibility with existing run folders)
- `TruthTrigRates.py`:
	- Truth-only skimming configuration


## !! IMPORTANT: How to update TrigRates_CA.py !!

`scripts/TrigRates_CA.py` is the **single source of truth**.
After every change, always sync it to all run directories in one command:

```bash
cp scripts/TrigRates_CA.py run_23hi/TrigRates_CA.py
cp scripts/TrigRates_CA.py run_24hi/TrigRates_CA.py
cp scripts/TrigRates_CA.py run_24pp/TrigRates_CA.py
cp scripts/TrigRates_CA.py run_25hi/TrigRates_CA.py
```

The script auto-selects the correct dataset based on the directory name it is run from (via `_run_dir = os.path.basename(os.getcwd())`), so the identical file works in every run folder without manual edits.


## How to run

### R25 (CA) — standard workflow

#### 1. Setup (once per login session)

```bash
cd /afs/cern.ch/user/y/yuhang/eos/dimuon/SkimCode
source setup_25.sh
```

`setup_25.sh` runs `asetup AthAnalysis,25.2.55` and configures the `acm` workarea in `build_25/`.

#### 2. Compile (after any C++ change in `source/`)

```bash
cd /afs/cern.ch/user/y/yuhang/eos/dimuon/SkimCode
source setup_25.sh && cd build_25 && acm compile
cd ..
```

(`setup_25.sh` must be sourced in the same shell call as `acm compile` because environment variables don't persist between separate shell invocations on this cluster.)

#### 3. Test run (100 events)

```bash
cd /afs/cern.ch/user/y/yuhang/eos/dimuon/SkimCode
source setup_25.sh && cd run_24pp && athena TrigRates_CA.py --evtMax=100
```

Or for a HI run folder:

```bash
source setup_25.sh && cd run_24hi && athena TrigRates_CA.py --evtMax=100
source setup_25.sh && cd run_25hi && athena TrigRates_CA.py --evtMax=100
```

If `athena` is not accepted by your environment wrapper, fallback:

```bash
python TrigRates_CA.py --evtMax=100
```

#### 4. Full run

```bash
cd <run_folder>
athena TrigRates_CA.py
```

Default `m_EvtMax = 1000` (set in the script). Override with `--evtMax=<N>`.

---

### R21 (legacy JO)

Important: run `setupATLAS -c centos7` first, then run setup inside the container shell.

```bash
cd /afs/cern.ch/user/y/yuhang/eos/dimuon/SkimCode
setupATLAS -c centos7
```

Inside the container:

```bash
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh
cd build_21
acmSetup --sourcedir=../source AthAnalysis,21.2.200
cd ../run_powheg_pp17fullsim
cp ../scripts/TrigRates_JO.py ./TrigRates_JO.py
athena TrigRates_JO.py --evtMax=100
```

If you still use an old tested local script in the run folder, you can run:

```bash
athena TrigRates.py --evtMax=100
```


## Output

- Default output file: `myfile.root`
- Tree name: `HeavyIonD3PD`


## Notes

- Keep R21 and R25 environments separate (do not mix setup scripts in one shell).
- If C++ in `source/HFtrigValidation/` changes, always rebuild (`acm compile`) in the matching release build directory.
- `setup_25.sh` and `acm compile` must be called in the **same shell invocation** (environment variables are lost between calls on this cluster).
- Recommended workflow:
	- R21 -> `TrigRates_JO.py`
	- R25 -> `TrigRates_CA.py`


## Dataset toggle cheat-sheet

Use exactly one dataset toggle as `True` in the chosen script, **or rely on automatic directory-based selection** (recommended).

### Common run folders

- `run_powheg_pp17fullsim` (R21 JO):
	- script: `TrigRates_JO.py` (or legacy local `TrigRates.py`)
	- toggle: `do_pp_MC_fullsim_17 = True`

- `run_23hi` (R25 CA):
	- script: `TrigRates_CA.py`
	- toggle: `do_hi2023 = True` (auto-detected from directory name)

- `run_24hi` (R25 CA):
	- script: `TrigRates_CA.py`
	- toggle: `do_hi2024 = True` (auto-detected from directory name)

- `run_24pp` (R25 CA):
	- script: `TrigRates_CA.py`
	- toggle: `do_pp2024 = True` (auto-detected from directory name)

- `run_25hi` (R25 CA):
	- script: `TrigRates_CA.py`
	- toggle: `do_hi2025 = True` (auto-detected from directory name)
	- GRL: `physics_HI2025_50ns_PbPb_IgnoreBSPOT_INVALID.xml` (2025 Pb+Pb run)


### Available toggles

- Data:
	- `do_hi2015`, `do_hi2018`, `do_hi2023`, `do_hi2024`, `do_hi2025`
	- `do_pp2015`, `do_pp2017`, `do_pp2024`
- MC fullsim:
	- `do_pp_MC_fullsim_17`, `do_pp_MC_fullsim_24`


### Quick sanity checks before running

- Exactly one `do_*` toggle is `True` (or rely on auto-detection from directory name)
- Input file path for that toggle exists
- Correct script for release:
	- R21: JO (`TrigRates_JO.py`)
	- R25: CA (`TrigRates_CA.py`)


## StoreTracks bitmap reference

`alg.StoreTracks` controls what track information is stored:

| Value | Branches stored |
|-------|----------------|
| `0`   | nothing (tracks disabled) |
| `1`   | `trk_numqual[8]` only — per-event track counts for: all/PPMinBias/HILoose/HITight (with and without pt > 400 MeV cut) |
| `1+2` | + per-track vectors: `trk_pt`, `trk_eta`, `trk_phi`, `trk_charge`, `trk_qual` |
| `1+2+4` | + hit details: `trk_d0`, `trk_z0_wrtPV`, pixel/SCT hits, chi2, etc. |
| `1+2+4+8` | + truth links: `trk_truth_index`, `trk_truth_prob`, etc. (MC only) |

Current setting: `StoreTracks = 1` for all runs (HI and pp).


## StoreZdc bitmap reference

`alg.StoreZdc` controls ZDC output:

| Value | Branches stored |
|-------|----------------|
| `0`   | ZDC disabled |
| `1`   | Basic: `zdc_ZdcAmp[2]`, `zdc_ZdcAmpErr[2]`, `zdc_ZdcEnergy[2]`, `zdc_ZdcEnergyErr[2]`, `zdc_ZdcTime[2]`, `zdc_ZdcStatus[2]`, `zdc_ZdcModuleMask`, `zdc_ZdcModulePreSampleAmp[2][4]` |
| `1+2` | + RPD/centroid: `zdc_RpdSubAmpSum[2]`, `zdc_xDetCentroid[2]`, `zdc_yDetCentroid[2]`, `zdc_xCentroid[2]`, `zdc_yCentroid[2]`, `zdc_xDetCentroidUnsub[2]`, `zdc_yDetCentroidUnsub[2]`, `zdc_xDetRowCentroidStdev[2]`, `zdc_yDetColCentroidStdev[2]`, `zdc_reactionPlaneAngle[2]`, `zdc_cosDeltaReactionPlaneAngle`, `zdc_centroidStatus[2]` |

Array index convention: `[0]` = A-side (`zdcSide < 0`), `[1]` = C-side (`zdcSide > 0`).
Current setting: `StoreZdc = 1` for HI runs; `0` for pp runs.


---

## Grid Submission Guide

### Key PanDA/pathena parameters

| Parameter | Description |
|-----------|-------------|
| `--nFilesPerJob N` | N input AOD files per grid job. Primary splitting knob. |
| `--nGBPerJob MAX` | Alternative: let PanDA auto-size to the max allowed per site (site-dependent, ~30–200 files for large PbPb AOD). |
| `--nEventsPerJob N` | Split by events rather than files (needs `--useAMIEventLevelSplit` for accurate AMI metadata). |
| `--mergeOutput` | Merge per-job ROOT outputs into one file per task at job completion. |
| `--extOutFile myfile.root` | Required with `--trf`-style CA submission to declare the output filename. |
| `--inDsTxt FILE` | File listing one `scope:dataset` per line; enables multi-dataset tasks. |

**Maximum jobs per task:** ATLAS PanDA enforces ~500 jobs per user-analysis task as a practical grid-quota guideline. Larger tasks experience reduced scheduling priority. Use `--nGBPerJob MAX` to let PanDA auto-optimize within the limit when replanning.

**Job wall-time estimate:** For our fast skimming algorithm (trigger bits + ZDC + track counts only):
- PbPb HardProbes AOD: ~11 GB/file, ~5800 events/file, read via XRootD
- Estimated processing: ~10–30 min per file over WAN XRootD
- At nFilesPerJob=60: expected wall time ~10–30 h/job (within typical 24–48 h site limit)
- Reduce to nFilesPerJob=30 if jobs frequently timeout

### Submitting grid jobs (CA-based, R25)

Setup before running `grid_sub.sh`:

```bash
cd SkimCode
source setup_25.sh          # sets up the Athena release
lsetup panda                # sets up panda client
voms-proxy-init -voms atlas # initialize grid certificate proxy
cd run_25hi
bash grid_sub.sh            # submits all 6 tasks
```

Monitor with `pbook` or [BigPanDA](https://bigpanda.cern.ch/?user=yuhang).

---

## Dataset inventory and partition plan — 2025 Pb+Pb HardProbes

**Source container:** `data25_hi.00512049.physics_HardProbes.merge.AOD.f1671_m2272`
**Query used:** `rucio list-dids "data25_hi:data25_hi.*physics_HardProbes.merge.AOD.f*" --filter type=DATASET`

### Per-dataset summary

| Run | f-tag | Files | Size (TB) | Events |
|-----|-------|------:|----------:|-------:|
| 510493 | f1655 | 177 | 1.7 | 913,064 |
| 510501 | f1655 | 653 | 4.6 | 2,380,959 |
| 510510 | f1655 | 1,452 | 16.4 | 8,208,685 |
| 510547 | f1655 | 1,281 | 13.5 | 6,826,787 |
| 510583 | f1655 | 1,077 | 11.0 | 5,541,343 |
| 510587 | f1658 | 1,899 | 20.5 | 10,772,329 |
| 510594 | f1658 | 1,045 | 11.3 | 5,951,748 |
| 510601 | f1658 | 104 | 0.4 | 244,062 |
| 510639 | f1658 | 4,373 | 47.1 | 24,671,037 |
| 510703 | f1658 | 5,850 | 65.2 | 33,617,500 |
| 510732 | f1658 | 1,951 | 22.6 | 11,672,220 |
| 510816 | f1658 | 5,863 | 65.7 | 34,138,968 |
| 510878 | f1661 | 6,401 | 72.6 | 37,713,469 |
| 510992 | f1661 | 3,080 | 34.7 | 18,041,617 |
| ~~511013~~ | f1661 | ~~1~~ | ~~0.004 GB~~ | ~~8~~ — **EXCLUDED** (empty run) |
| 511020 | f1661 | 4,211 | 46.5 | 24,251,231 |
| 511026 | f1661 | 4,338 | 48.0 | 25,066,852 |
| 511035 | f1661 | 4,032 | 43.7 | 22,827,689 |
| 511094 | f1661 | 1,933 | 21.8 | 11,362,830 |
| 511115 | f1661 | 2,635 | 28.7 | 14,989,285 |
| 511132 | f1661 | 456 | 5.0 | 2,628,288 |
| 511151 | f1668 | 3,632 | 41.1 | 21,378,815 |
| 511244 | f1668 | 3,907 | 43.2 | 22,538,346 |
| 511278 | f1668 | 3,738 | 41.4 | 21,576,139 |
| 511360 | f1668 | 13 | 0.1 | 30,842 |
| 511382 | f1668 | 391 | 4.4 | 2,292,474 |
| 511399 | f1668 | 2,184 | 24.6 | 12,808,315 |
| 511436 | f1668 | 4,782 | 53.1 | 27,695,803 |
| 511463 | f1668 | 4,579 | 50.9 | 26,524,848 |
| 511506 | f1668 | 4,794 | 53.1 | 27,687,139 |
| 511519 | f1668 | 4,044 | 44.9 | 23,378,777 |
| 511552 | f1668 | 4,425 | 48.3 | 25,117,779 |
| 511617 | f1668 | 2,182 | 24.8 | 12,913,701 |
| 511649 | f1668 | 102 | 1.1 | 567,288 |
| 511650 | f1668 | 1,440 | 16.2 | 8,135,377 |
| 511658 | f1671 | 4,798 | 51.7 | 27,424,097 |
| 511703 | f1674 | 76 | 0.7 | 411,890 |
| 511783 | f1671 | 4,342 | 48.1 | 25,288,410 |
| 511857 | f1671 | 5,026 | 55.3 | 29,109,824 |
| 511866 | f1671 | 4,915 | 53.9 | 28,416,691 |
| 511894 | f1671 | 4,555 | 50.1 | 26,424,662 |
| 511902 | f1671 | 5,319 | 58.7 | 30,938,864 |
| 511920 | f1671 | 4,733 | 52.3 | 27,549,147 |
| 511967 | f1671 | 5,286 | 58.2 | 30,689,624 |
| 511978 | f1671 | 4,335 | 47.8 | 25,158,610 |
| 512013 | f1671 | 5,008 | 58.3 | 29,376,001 |
| 512028 | f1671 | 1,001 | 11.2 | 5,911,121 |
| 512049 | f1671 | 1,110 | 12.6 | 6,646,963 |

**Totals (47 active datasets):** 143,528 files | 1587 TB | ~827.8 M events
**Average per file:** 11.3 GB | 5,768 events

### Partition plan (nFilesPerJob=60, ≤500 jobs/task)

**Strategy:** Greedy sequential grouping by run number. Each task keeps ≤ 500 × 60 = 30,000 files to stay within the ~500 jobs/task ATLAS PanDA user-analysis guideline. This maximises events per task at the given per-job limit.

| Part | Datasets (runs) | Files | Jobs | Size (TB) | Events |
|------|-----------------|------:|-----:|----------:|-------:|
| 1 | 510493–510816 (12 ds) | 25,725 | 429 | 279.9 | 144.9 M |
| 2 | 510878–511132 (8 ds)  | 27,086 | 452 | 301.1 | 156.9 M |
| 3 | 511151–511506 (9 ds)  | 28,020 | 467 | 311.8 | 162.5 M |
| 4 | 511519–511857 (9 ds)  | 26,435 | 441 | 291.2 | 152.3 M |
| 5 | 511866–511978 (6 ds)  | 29,143 | 486 | 321.0 | 169.2 M |
| 6 | 512013–512049 (3 ds)  |  7,119 | 119 |  82.2 |  41.9 M |
| **Total** | **47 datasets** | **143,528** | **2,394** | **1587** | **827.8 M** |

**Events per task (max):** ~346,000 events/job × 500 jobs = **173 M events/task**
**Total expected grid output:** ~827.8 M events across 6 tasks

Partition text files are in `run_25hi/`: `InDstxt_PbPb2025_5p36TeV_part{1..6}.txt`
Grid submission script: `run_25hi/grid_sub.sh`

### Notes on dataset scale (2025 vs 2024)

The 2025 Pb+Pb run (Oct–Nov 2025) produced dramatically more data than 2024:
- 2024 PbPb HardProbes: ~41 datasets, processed in 3 tasks with nFilesPerJob=60
- 2025 PbPb HardProbes: 47 active datasets, **requires 6 tasks** with the same per-job setting
- The largest individual runs (510878, 510703, 510816) each have 5,800–6,400 files (~65–73 TB), comparable to the entire 2024 dataset
- Run 511013 (1 file, 8 events) is excluded as an empty/corrupt run

### Re-running with different nFilesPerJob

If jobs time out (wall-time exceeded), reduce `nFilesPerJob`:
- `nFilesPerJob=30` → ~4,785 total jobs → ~10 tasks
- `nFilesPerJob=40` → ~3,588 total jobs → ~8 tasks

If you want fewer tasks and jobs do not time out, increase `nFilesPerJob`:
- `nFilesPerJob=100` → ~1,436 total jobs → ~3 tasks (Part 5 alone: 292 jobs; all 6 parts fit in 3 tasks)

Alternative: use `--nGBPerJob MAX` to let PanDA auto-determine per-job input size based on site capabilities. This is the most adaptive option but gives less predictable job counts.
