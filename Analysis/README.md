# Dimuon Continuum Analysis

Low-mass dimuon continuum measurement in pp and Pb+Pb collisions at the LHC (ATLAS).
Covers muon-pair production in the mass range ~1-3 GeV from open heavy-flavour (b/c-quark) decays.

**This is a Run 3 analysis.** Run 2 code paths (pp17, PbPb15/18) exist for cross-checks but are no longer maintained for new analysis decisions as of March 2026. All active development targets Run 3 datasets (pp24, PbPb23/24/25).

## Environment

- Cluster: BNL SDCC (`ssh spar*.sdcc.bnl.gov`)
- Setup: `source ~/setup.sh` (provides ROOT, hadd, condor tools via LCG_107a_ATLAS_2)
- C++ compilation: ROOT ACLiC (`.L MyClass.cxx+`), not CMake
- No `pip install` - use packages in the LCG release

## Repository layout

```
Analysis/
|-- NTupleProcessingCode/   Stage 1: NTuple processing (CRTP algorithm cores)
|-- RDFBasedHistFilling/    Stage 2: RDataFrame histogram filling
|-- plotting_codes/         Stage 3: Plotting
|   |-- reco_effcy/           Reco efficiency & detector response plotters (fullsim)
|   |-- trig_effcy/           Trigger efficiency plotters (data)
|   |-- pythia_plotting_codes/  Pythia truth flavor/origin plots
|   |-- mc_data_compr/        MC-data comparison plots
|   |-- pbpb_pp_compr/        PbPb vs pp comparison
|   |-- event_selection/      PbPb event selection plots
|   +-- ...
|-- pipelines/              End-to-end pipeline scripts
|-- MuonObjectsParamsAndHelpers/  Data types (Muon, MuonPair), params, base classes
|-- Utilities/              Plot helpers, general utilities
|-- EfficiencyCorrs/        Efficiency correction derivation
|-- SingleBAnalysis/        Single-b extraction analysis
|-- docs/                   Per-pipeline documentation
+-- ScrambGen/              Scrambled generator utilities
```

## Three-stage analysis pattern

Every analysis pipeline follows the same three stages:

| Stage | What | Where | Typical output |
|-------|------|-------|----------------|
| 1. NTuple processing | Read NTUP/xAOD, apply selections, write flat trees | `NTupleProcessingCode/` | `muon_pairs_*.root` |
| 2. Histogram filling | Read flat trees, fill histograms via RDataFrame | `RDFBasedHistFilling/` | `histograms_*.root` |
| 3. Plotting | Read histograms, make efficiency/response/comparison plots | `plotting_codes/` | `.png` files |

Pipelines in `pipelines/` automate all three stages with inter-stage validation. Each stage's output is checked before the next stage proceeds; corrupted or missing output causes an immediate exit with an error.

## Class architecture

### NTuple processing (Stage 1)

Uses the CRTP (Curiously Recurring Template Pattern) with variadic Extras mixins:

```
DimuonAlgCoreT<PairT, MuonT, Derived, ...Extras>     (base algorithm template)
  |-- DimuonDataAlgCoreT                               (pp / PbPb data)
  |-- PythiaAlgCoreT                                   (Pythia MC)
  +-- PowhegAlgCoreT                                   (Powheg MC)
```

**Extras mixins** (composed via multiple inheritance):
- `PPExtras` / `PbPbExtras` -- collision-system specifics
- `PythiaTruthExtras` -- truth ancestry tracing
- `PythiaFullSimExtras` -- reco-truth matching for fullsim
- `PythiaFullSimOverlayExtras` -- overlay-specific (centrality, FCal)
- `PowhegTruthExtras` / `PowhegFullSimExtras` / `PowhegFullSimOverlayExtras`

**Concrete analysis classes** (in `*AnalysisClasses.h`):

| Class | Template args | Purpose |
|-------|--------------|---------|
| `PPAnalysis` | `DimuonDataAlgCoreT` + `PPExtras` | pp data |
| `PbPbAnalysis` | `DimuonDataAlgCoreT` + `PbPbExtras` | Pb+Pb data |
| `PythiaTruthAnalysis` | `PythiaAlgCoreT` + `PythiaTruthExtras` | Pythia truth-level |
| `PythiaFullSimAnalysis` | `PythiaAlgCoreT` + `PythiaFullSimExtras` + `PythiaTruthExtras` | Pythia fullsim pp |
| `PythiaFullSimOverlayAnalysis` | `PythiaAlgCoreT` + `PythiaFullSimExtras` + `PythiaFullSimOverlayExtras` + `PythiaTruthExtras` | Pythia fullsim overlay |
| `PowhegTruthAnalysis` | `PowhegAlgCoreT` + `PowhegTruthExtras` | Powheg truth-level |
| `PowhegFullSimAnalysisWTruth` | `PowhegAlgCoreT` + `PowhegFullSimExtras` + `PowhegTruthExtras` | Powheg fullsim |

### Histogram filling (Stage 2)

```
RDFBasedHistFillingBaseClass                (base: variable registration, hist booking, I/O)
  |-- RDFBasedHistFillingPP                 (pp data)
  |-- RDFBasedHistFillingPbPb               (PbPb data)
  |-- RDFBasedHistFillingPythiaTruth        (Pythia truth)
  |-- RDFBasedHistFillingPythia             (Pythia reco, base for fullsim)
  |   +-- RDFBasedHistFillingPythiaFullsim  (Pythia fullsim pp)
  |       +-- RDFBasedHistFillingPythiaFullsimOverlay  (+ PbPbBaseClass, per-centrality)
  |-- RDFBasedHistFillingPowhegTruth        (Powheg truth)
  +-- RDFBasedHistFillingPowheg             (Powheg reco, base for fullsim)
      +-- RDFBasedHistFillingPowhegFullsim  (Powheg fullsim)
```

### Plotting (Stage 3)

**Reco efficiency / detector response** (`plotting_codes/reco_effcy/`):
- `PythiaFullsimRecoEffPlotterBase` -> `PythiaFullsimRecoEffPlotter` (pp) / `PythiaFullsimRecoEffPlotterOverlay` (overlay, per-centrality)
- `PowhegFullsimRecoEffPlotter` -- Powheg fullsim pair reco efficiency
- `PowhegFullsimDetRespPlotterSingleMuon` -- Powheg single-muon detector response

**Trigger efficiency** (`plotting_codes/trig_effcy/`):
- `TrigEffPlotterBaseClass` -> `TrigEffPlotterPP` / `TrigEffPlotterPbPb`

### Shared data types (`MuonObjectsParamsAndHelpers/`)

- `Muon.h` -- muon struct variants (`MuonBase`, `MuonPythiaExtra`, `MuonFullsimExtra`, etc.)
- `MuonPair*.h` -- pair structs (truth, reco, fullsim, overlay variants)
- `PbPbBaseClass.h` -- centrality binning (CRTP), cross-section factors
- `PPBaseClass.h` -- pp-specific params
- `FullSimSampleType.h` -- enum `{pp, hijing, zmumu, data}` with I/O path helpers
- `ParamsSet.h` -- analysis parameter sets (mass cuts, eta ranges, etc.)

## Analysis pipelines

| Pipeline | Script | Modes | Doc |
|----------|--------|-------|-----|
| Pythia truth | `pipelines/pipeline_pythia_truth.sh` | `private`, `nonprivate_5p36`, `nonprivate_5p02` | [docs/pythia_truth.md](docs/pythia_truth.md) |
| Pythia fullsim pp | (manual, see doc) | pp24 | [docs/pythia_fullsim_pp.md](docs/pythia_fullsim_pp.md) |
| Pythia fullsim overlay | `pipelines/pipeline_pythia_fullsim_overlay.sh` | `hijing`, `zmumu`, `data` | [docs/pythia_fullsim_overlay.md](docs/pythia_fullsim_overlay.md) |
| Powheg truth | (manual) | bb, cc | [docs/powheg.md](docs/powheg.md) |
| Powheg fullsim single-muon | `pipelines/pipeline_powheg_fullsim_single_muon.sh` | bb+cc | [docs/powheg.md](docs/powheg.md) |
| Powheg fullsim mixed pairs | `pipelines/pipeline_powheg_fullsim_mixed_pairs.sh` | `--mass-filter`/`--no-mass-filter` | [docs/powheg.md](docs/powheg.md) |
| pp data | (manual) | pp17, pp24 | [docs/data_analysis.md](docs/data_analysis.md) |
| PbPb data | (manual) | PbPb15, PbPb18, PbPb23, PbPb24 | [docs/data_analysis.md](docs/data_analysis.md) |

## Quick start

```bash
# Run Pythia fullsim overlay pipeline (hijing test sample)
cd Analysis
./pipelines/pipeline_pythia_fullsim_overlay.sh hijing --dry-run

# Run Pythia truth pipeline (nonprivate 5.36 TeV)
./pipelines/pipeline_pythia_truth.sh nonprivate_5p36

# Compile and run a single stage manually
source ~/setup.sh
cd RDFBasedHistFilling
root -l -b -q '.L RDFBasedHistFillingPP.cxx+'
```
