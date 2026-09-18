# Pythia Truth Analysis

Generator-level Pythia8 dimuon analysis: traces muon ancestry to tag heavy-flavour origin (b, c, light), forms opposite-sign and same-sign pairs, and fills kinematic distributions categorized by flavour and origin.

## Sample modes

| Mode | Energy | Sample | Condor submit |
|------|--------|--------|---------------|
| `private` | 5.36 TeV | Private production (2 batches) | `run_pythia_truth_kn_private.sub` |
| `nonprivate_5p36` | 5.36 TeV | Central production (6 kn ranges) | `run_pythia_truth_kn_nonprivate.sub` |
| `nonprivate_5p02` | 5.02 TeV | Central production (6 kn ranges) | `run_pythia_truth_kn_nonprivate.sub` (args overridden by pipeline) |

Each kinematic range (kn 0-5) covers a different pT-hat slice of the Pythia generation. Non-private samples split into 6 Condor jobs (one per kn), then hadd the results.

**Which evgen, which weights.** The 5.36 TeV central sample IS the **nPDF** Pythia evgen (e8599,
DSIDs 802758–802781: four isospin beams pp/pn/np/nn combined 4:6:6:9, nuclear PDF
`nNNPDF30_nlo_as_0118_A208_Z82`) — the evgen made for **Pb+Pb** conditions, which is also the
simulation input of every HIJING overlay. Its AMI weights are read from the evgen's own
`pythia_truth_full_sample/pythia_5p36TeV/ami_info_nPDF/` (a missing file is fatal;
`docs/ami_weights.md` table A). The pp-conditions counterpart — the pp-only, proton-PDF `_pdf`
evgen (803015–803020, `ami_info_PDF/`) — is the input of the pp24 fullsim FULL sample and has
**no truth-only skim**: the fullsim NTUP carries its complete truth block (all events stored), so
generator-level pp-conditions quantities can be taken from the fullsim's truth pair trees (same `PythiaTruthExtras` ancestry tracing) without a new skim.

## Key files

### NTuple processing (Stage 1)

| File | Role |
|------|------|
| `NTupleProcessingCode/PythiaAnalysisClasses.h` | Defines `PythiaTruthAnalysis` |
| `NTupleProcessingCode/PythiaAlgCoreT.{h,c}` | Template algorithm core |
| `NTupleProcessingCode/PythiaTruthExtras.{h,c}` | Truth ancestry tracing, HF origin tagging |
| `NTupleProcessingCode/run_pythia_truth_kn.sh` | Condor worker: `./run_pythia_truth_kn.sh <kn> [is_private] [e_com] [use_local]` |

### Histogram filling (Stage 2)

| File | Role |
|------|------|
| `RDFBasedHistFilling/RDFBasedHistFillingPythiaTruth.cxx` | RDF filler: flavor/origin categorized histograms |
| `RDFBasedHistFilling/run_rdf_pythia_truth.sh` | Runner: `./run_rdf_pythia_truth.sh [private|nonprivate] [5.36|5.02] [withcuts|nocuts]` |

### Plotting (Stage 3)

| File | Role |
|------|------|
| `plotting_codes/pythia_plotting_codes/plot_flavor_categorized_kinematics.cxx` | Flavor-tagged kinematic plots |
| `plotting_codes/pythia_plotting_codes/plot_origin_categorized_kinematics.cxx` | Origin-tagged kinematic plots |

## Class hierarchy

```
PythiaTruthAnalysis
  inherits PythiaAlgCoreT<MuonPairPythiaTruth, MuonPythiaTruth, Self, PythiaTruthExtras>
  inherits PythiaTruthExtras  -- truth ancestry, HF origin tagging

RDFBasedHistFillingPythiaTruth
  inherits RDFBasedHistFillingBaseClass
```

## Pipeline

```bash
./pipelines/pipeline_pythia_truth.sh nonprivate_5p36
```

### Stage 1: NTuple processing (Condor)

Submits 6 jobs (kn 0-5) for non-private, 2 jobs for private. Each job:
- Reads Pythia NTUP files for one kinematic range
- Traces truth muon ancestry to determine HF origin
- Forms muon pairs, writes `muon_pair_tree_sign1/sign2`

Output per batch: `muon_pairs_pythia_{ecom}_kn{N}_no_data_resonance_cuts.root`

### Stage 2: Validate + hadd + RDF

1. Validates all per-batch ROOT files
2. `hadd` into combined file: `muon_pairs_pythia_{ecom}_no_data_resonance_cuts.root`
3. Validates combined file (both sign trees non-empty)
4. Runs RDF histogram filling

### Stage 3: Plotting

Runs flavor- and origin-categorized kinematic plotters.

## Running manually

```bash
source ~/setup.sh

# Single kn range
cd NTupleProcessingCode
root -b -l <<'EOF'
.L PythiaAnalysisClasses.h
PythiaTruthAnalysis py(0, false, 5.36);  // kn=0, nonprivate, 5.36 TeV
py.Run();
.q
EOF

# RDF filling
cd ../RDFBasedHistFilling
./run_rdf_pythia_truth.sh nonprivate 5.36 nocuts
```
