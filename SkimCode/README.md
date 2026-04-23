# SkimCode — ATLAS dimuon AOD skim

Runs the `TrigRates` Athena algorithm to produce a flat ntuple (`HeavyIonD3PD`
TTree) from heavy-ion and pp AODs, for dimuon trigger + physics studies.

Main features of `TrigRates`:

- Event cleaning (GRL, detector-error checks, vertex cut)
- Trigger decisions and prescales (L1 + HLT)
- Muon trigger matching (single-muon and dimuon chains)
- Muon MCP tools: selection, momentum calibration, efficiency corrections (Medium + Tight), trigger matching
- Optional storage of tracks, vertex, MET, truth, and MC event info
- ZDC information (energy, time, status, per-module PreSampleAmp; optional RPD centroid)
- HI event-shape quantities (FCal Et A/C, Q-vectors, centrality)


## Repository layout

```
SkimCode/
├── setup_21.sh / setup_25.sh    # base AthAnalysis release setup
├── source/                      # C++ package (shared across releases)
│   └── HFtrigValidation/
│       ├── HFtrigValidation/    # public headers (TrigRates.h, Module_*.h)
│       └── src/                 # TrigRates.cxx + Module_EventShape.cxx etc.
├── build_21/ build_25/          # per-release cmake build dirs (gitignored)
├── scripts/                     # canonical Python/bash templates
│   ├── TrigRates_CA.py          # CA config (master copy; R25)
│   ├── TrigRates_JO.py          # legacy JO config (R21)
│   ├── TruthTrigRates.py        # truth-only skim config
│   ├── grid_sub.sh              # pathena submission template
│   ├── print_ds_size.sh         # rucio helper: list datasets + sizes
│   ├── add_replication_rule_data23.sh  # batch rucio add-rule example
│   └── useful_rucio_commands.txt
├── run_<year/sample>/           # one per dataset family (see "Run modes")
├── xmls/                        # GRL XMLs (data-quality lumiblock lists)
└── datasetnames/                # txt files listing rucio datasets per year
```

Available AthAnalysis releases: see
<https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/AnalysisRelease>.


## Build

```bash
cd SkimCode
source setup_25.sh
```

`setup_25.sh` runs `asetup AthAnalysis,25.2.89` and configures the `acm` work
area in `build_25/`. After any C++ change in `source/`, in the same shell:

```bash
cd build_25 && acm compile && cd ..
```

If you hit mysterious build problems, delete `build_25/` and re-run
`setup_25.sh`.


## Run modes (`TRIGRATES_RUNMODE`)

`scripts/TrigRates_CA.py` selects per-year / per-sample configuration via a
single global run mode. Dispatch priority (first match wins):

1. Environment variable `TRIGRATES_RUNMODE=<mode>` (most robust; used in grid submission).
2. `--filesInput=...` substring match (`data23_hi` → `hi2023`, etc.).
3. `basename($PWD)` (running from `run_25hi/` implies `hi2025`, etc.).

Supported modes:

| Mode                       | Year  | Sample                                | Input     | GRL                                  | Triggers | ZDC readout |
|----------------------------|------:|---------------------------------------|-----------|--------------------------------------|----------|-------------|
| `hi2023`                   | 2023  | data23_hi Pb+Pb                       | data      | v120-pro33-03                        | yes      | yes         |
| `hi2024`                   | 2024  | data24_hi Pb+Pb                       | data      | HI2024_50ns                          | yes      | yes         |
| `hi2025`                   | 2025  | data25_hi Pb+Pb                       | data      | HI2025_50ns                          | yes      | yes         |
| `pp2024`                   | 2024  | data24_5p36TeV pp reference           | data      | 2024ppRef_25ns                       | yes      | no          |
| `ppmcfullsim2024`          | –     | Pythia8 pp fullsim MC                 | geant4    | (none)                               | off      | no          |
| `ppmcfullsim_hioverlay24`* | 2024  | Pythia8 pp fullsim + HIJING overlay   | geant4    | (none)                               | off      | **off** — overlay AOD has no ZDC |

*Lives in `run_pythia_fullsim_HIJING_overlay/TrigRates_CA.py`, not in
`scripts/TrigRates_CA.py`. Uses a specialized config that adds the overlay run
mode and an environment-driven output filename (`TRIGRATES_OUTPUT_FILE`,
see below). Keep in mind when syncing from `scripts/`.

Centrality calibration (`Module_EventShape.cxx`, as of 2026-04): Run-3 years
2023/24/25/26 all map to the PbPb2023 FCal-ET thresholds. Update the switch
when per-year calibrations are finalised.


### Sync `scripts/TrigRates_CA.py` → run dirs

`scripts/TrigRates_CA.py` is the canonical copy for the data-skim run modes
(`hi*`, `pp*`, `ppmcfullsim2024`). After every change:

```bash
for d in run_23hi run_24hi run_24pp run_25hi run_pythia_fullsim; do
  cp scripts/TrigRates_CA.py $d/TrigRates_CA.py
done
```

`run_pythia_fullsim_HIJING_overlay/` has its own specialization and should be
synced manually when shared config (triggers, GRL paths, muon-tool settings)
changes.


## Workflow

### Local test (one AOD, 20–100 events)

```bash
cd run_<dir>/
source ../setup_25.sh
athena.py TrigRates_CA.py --evtMax=20
```

Output: `myfile.root` (TTree `HeavyIonD3PD`). Default `m_EvtMax = 1000` if
`--evtMax` is omitted. If `athena.py` isn't found, source of `setup_25.sh`
failed silently — re-run without piping its output (see Pitfalls).

### Grid submission

```bash
cd run_<dir>/
source ../setup_25.sh && lsetup panda
voms-proxy-init -voms atlas        # if proxy expired
bash grid_sub.sh                   # submits one pathena task per line
```

Bump the `.v<N>.` tag in `grid_sub.sh` before every resubmission; task names
are immutable and reusing a tag fails at submission.

Monitor:
- `pbook` (interactive)
- BigPanDA web UI: <https://bigpanda.cern.ch/?user=yuhang>


## Output conventions

### Dataset naming

```
user.yuhang.<type>.<sample>.<campaign>.v<N>.part<K>.
```
Trailing dot is required by pathena (it appends a 10-char hash).

### NTUP file naming (`TRIGRATES_OUTPUT_FILE`)

The `THistSvc` output filename in `TrigRates_CA.py` is (in the HIJING-overlay
run dir) taken from `TRIGRATES_OUTPUT_FILE`, defaulting to `myfile.root` for
local tests. In `grid_sub.sh`, both the env var and `--extOutFile` are set to
the same sample-specific name so that pathena picks up the produced file:

```bash
pathena --trf "TRIGRATES_RUNMODE=... TRIGRATES_OUTPUT_FILE=${fname} athena.py TrigRates_CA.py --filesInput=%IN --evtMax=%MAXEVENTS" \
        --inDS ... --outDS ... --extOutFile "${fname}" --mergeOutput ...
```

`scripts/TrigRates_CA.py` currently hardcodes `myfile.root`. Standard
data-skim run dirs therefore produce `myfile.root` inside their output
datasets (downloaded files land as
`user.yuhang.<...>._EXT0.myfile.root`).


## Output branch reference

### `StoreTracks` bitmap (`alg.StoreTracks`)

| Value     | Branches stored |
|-----------|-----------------|
| `0`       | nothing |
| `1`       | `trk_numqual[8]` only — per-event track counts for all / PPMinBias / HILoose / HITight, with and without pt > 400 MeV cut |
| `1+2`     | + per-track vectors: `trk_pt`, `trk_eta`, `trk_phi`, `trk_charge`, `trk_qual` |
| `1+2+4`   | + hit details: `trk_d0`, `trk_z0_wrtPV`, pixel/SCT hits, chi² |
| `1+2+4+8` | + truth links (MC only): `trk_truth_index`, `trk_truth_prob` |

Current setting: `StoreTracks = 1` in all run modes.

`trk_numqual` layout: `[0–3]` with pt > 400 MeV (all / PPMinBias / HILoose /
HITight); `[4–7]` no pt cut.

### `StoreZdc` bitmap (`alg.StoreZdc`)

| Value | Branches stored |
|-------|-----------------|
| `0`   | ZDC disabled |
| `1`   | Basic: `zdc_ZdcAmp[2]`, `zdc_ZdcAmpErr[2]`, `zdc_ZdcEnergy[2]`, `zdc_ZdcEnergyErr[2]`, `zdc_ZdcTime[2]`, `zdc_ZdcStatus[2]`, `zdc_ZdcModuleMask`, `zdc_ZdcModulePreSampleAmp[2][4]` |
| `1+2` | + RPD / centroid: `zdc_RpdSubAmpSum[2]`, `zdc_xDetCentroid[2]`, `zdc_yDetCentroid[2]`, `zdc_xCentroid[2]`, `zdc_yCentroid[2]`, `zdc_xDetCentroidUnsub[2]`, `zdc_yDetCentroidUnsub[2]`, `zdc_xDetRowCentroidStdev[2]`, `zdc_yDetColCentroidStdev[2]`, `zdc_reactionPlaneAngle[2]`, `zdc_cosDeltaReactionPlaneAngle`, `zdc_centroidStatus[2]` |

Array index: `[0]` = A-side (`zdcSide < 0`), `[1]` = C-side (`zdcSide > 0`).

Current setting:
- HI data modes (`hi2023/24/25`): `StoreZdc = 1`.
- pp modes and Pythia fullsim: `StoreZdc = 0`.
- HIJING-overlay MC: `StoreZdc = 0` — the overlay AOD has no
  `ZdcSums`/`ZdcModules` containers even though `is_HION = True`.

### `StoreEventInfo` bitmap

Flag `& 1` always writes `RunNumber / lbn / bcid / eventNumber / (Act|Avg)IntPerXing`. Flag `& 4` additionally writes `NumTrackNoCuts`, `NumTrackPPMinBias`, `NumTrackHILoose`, `NumTrackHITight`, `FCalET_aux`, `FCalETP_aux`, `FCalETN_aux` — **but** those auxdata entries come from a private derivation and don't exist on production AODs. Leave `StoreEventInfo = 1` unless you know you're reading a derivation that has them.


## Non-obvious details

- **HIJING-overlay MC has no ZDC containers.** The reco chain skips ZDC
  digitization. `StoreZdc = 1` will fail with
  `No valid proxy for object ZdcSums`. The overlay-specific config forces
  `StoreZdc = 0`.
- **Run-3 centrality uses 2023 bins for all years.** See
  `Module_EventShape.cxx` — the year switch falls through for 2023–2026.
  Events in a mode without a year case (e.g. pp-only MC with
  `RunYear = 0`) get the default sentinel and should be treated as
  unreliable centrality.
- **`--nGBPerJob MAX` is mandatory for large PbPb AOD tasks.** Without it
  PanDA assumes local staging (~490 GB scratch/job for 60 × ~8 GB files) and
  no site qualifies. `MAX` enables XRootD remote-read brokerage; only
  output + workdir disk is required (~50 GB).
- **`HLT_MuonsCB_RoI` / `HLT_MuonsCB_FS` keys are data-only.** For MC modes
  the HLT muon container keys are left empty; setting them on a MC AOD
  causes retrieval failures.
- **`asetup` must run in the same shell that later calls `athena.py`.**
  Piping its output (e.g. `source foo.sh | tail -3`) puts `source` in a
  subshell — env changes never reach the caller. Use
  `source foo.sh > /tmp/log 2>&1` instead, then inspect the log.


## External resources

- **BigPanDA** (no auth for JSON API):
  - Task info: `https://bigpanda.cern.ch/task/<tid>/?json`
  - Failed jobs: `https://bigpanda.cern.ch/jobs/?jeditaskid=<tid>&jobstatus=failed&json&limit=2000`
  - User dashboard: <https://bigpanda.cern.ch/?user=yuhang>
- **Rucio scopes**:
  - user outputs: `user.yuhang`
  - data inputs: `data23_hi`, `data24_hi`, `data25_hi`, `data24_5p36TeV`
  - MC inputs: `mc23_5p36TeV`, …
- **pathena client reference**:
  <https://panda-wms.readthedocs.io/en/latest/client/pathena.html>
  (full option list also via `pathena --help`).
- **ATLAS AnalysisRelease list**:
  <https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/AnalysisRelease>


## Data-staging helpers (tape-only datasets)

If the Rucio query shows an input container has replicas only on tape
(common for older `data23_hi` reprocessings), use the helpers in `scripts/`
before submitting:

- `scripts/print_ds_size.sh` — list datasets with size (TB/GB/MB) for a
  given year/collision-type/stream. Copy into the target run dir and edit
  the filter.
- `scripts/add_replication_rule_data23.sh` — batched `rucio add-rule`
  commands to pull datasets to a SCRATCHDISK RSE. See the script header
  for the expected form.
- `scripts/useful_rucio_commands.txt` — grab-bag of common rucio
  invocations (list-rules, list-dataset-replicas, update-rule, quota
  checks).

Workflow sketch:
1. `rucio list-dataset-replicas <scope:container>` — check disk vs tape.
2. If disk-only: submit directly.
3. If partial or tape-only: add a replication rule to a SCRATCHDISK RSE
   with enough free quota, wait for OK, then submit. For stuck rules, a
   user can only edit `--source-replica-expression`; `--boost-rule` and
   `--priority` are privileged.


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

### Partition plan (`--nGBPerJob MAX`, ≤500 jobs/task)

**Strategy:** sequential greedy grouping by run number. Each task targets
~30 k input files (≤500 jobs × ~60 files/job) to stay within the
~500-jobs/user-task PanDA practical guideline.

| Part | Datasets (runs) | Files | Jobs | Size (TB) | Events |
|------|-----------------|------:|-----:|----------:|-------:|
| 1 | 510493–510816 (12 ds) | 25,725 | 429 | 279.9 | 144.9 M |
| 2 | 510878–511132 (8 ds)  | 27,086 | 452 | 301.1 | 156.9 M |
| 3 | 511151–511506 (9 ds)  | 28,020 | 467 | 311.8 | 162.5 M |
| 4 | 511519–511857 (9 ds)  | 26,435 | 441 | 291.2 | 152.3 M |
| 5 | 511866–511978 (6 ds)  | 29,143 | 486 | 321.0 | 169.2 M |
| 6 | 512013–512049 (3 ds)  |  7,119 | 119 |  82.2 |  41.9 M |
| **Total** | **47 datasets** | **143,528** | **2,394** | **1587** | **827.8 M** |

Partition text files: `run_25hi/InDstxt_PbPb2025_5p36TeV_part{1..6}.txt`.
Submission script: `run_25hi/grid_sub.sh`.

### Notes on dataset scale (2025 vs 2024)

- 2024 PbPb HardProbes: ~41 datasets, processed in 3 tasks with `nFilesPerJob=60`.
- 2025 PbPb HardProbes: 47 active datasets, **requires 6 tasks** at the same
  per-job setting. The largest individual runs (510878, 510703, 510816) each
  have 5,800–6,400 files (~65–73 TB), comparable to the entire 2024 dataset.
- Run 511013 (1 file, 8 events) is excluded as an empty/corrupt run.

### Re-running with different `nFilesPerJob`

| `nFilesPerJob` | Total jobs | Tasks needed |
|---------------:|-----------:|-------------:|
| 30 | ~4,785 | ~10 |
| 40 | ~3,588 |  ~8 |
| 60 | ~2,394 |   6 |
| 100 | ~1,436 |  ~3 |

Alternative: `--nGBPerJob MAX` lets PanDA auto-size per site; most adaptive
but per-task job count is less predictable.


## Legacy R21 JO workflow

For legacy samples that still need R21 (`run_powheg_pp17fullsim/`, etc.):

```bash
setupATLAS -c centos7      # enter CentOS7 container
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

Legacy local scripts such as `TrigRates.py` are kept in some run dirs for
compatibility; new work should use `TrigRates_JO.py` (R21) or
`TrigRates_CA.py` (R25). Do not mix R21 and R25 setups in the same shell.
