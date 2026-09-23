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

## External truth input, raw counts (opt-in mode)

**Purpose.** Run the *identical* truth processing (same pair building, same `PassCuts` truth
fiducial cut pT > 4.5 GeV and |η| < 2.4 on both legs, same pair kinematics, same
`PythiaTruthExtras` ancestry tracing, same output trees `muon_pair_tree_sign{1,2}`,
`muon_pair_tree_kin<K>_sign<d>`, `meta_tree_out`) on truth-only NTUPs **other than** the central
nPDF sample — e.g. new test samples in the same SkimCode `TruthTrigRates.py` format that have
**no AMI record yet** — and get unweighted pair counts from which per-event rates can be formed.

**What differs from a nominal run: only the input, the output, and the weight.**

| | nominal | external mode |
|---|---|---|
| input | `pythia_5p36TeV_<beam>_hQCD_DiMu_pTH<lo>_<hi>.TRUTH0.NTUP.root`, 4 beams | `truth_input_override`: file(s) of ONE slice `kn` (glob and/or space/comma list), chained as the pp slot (beam 0) |
| weight | σ·ε_filt(AMI) · isospin ratio / N_proc (nb) | **exactly 1**; AMI read and DSID guard skipped |
| `crossx` | σ·ε_filt·ratio | N_proc (same formula, weight·N/eff) |
| output | `.../pythia_truth_full_sample/pythia_5p36TeV/muon_pairs_pythia_5p36TeV_kn<N>_...` | `<truth_output_dir_override>/muon_pairs_pythia_truth_<label>_kn<N>_no_data_resonance_cuts.root` (+ `hists_pythia_ntuple_processing_truth_<label>_kn<N>_...`, diagnostics in `<dir>/ancestor_tracing_output_diagnostic_files/*_truth_<label>_kn<N>.txt`) |
| `meta_tree_out` | not filled (0 entries) on the truth path | 1 entry: `nproc_kin<N>_beam0` = events processed, `nbeam_kin<N>_beam0` = events in the input |

The kinematic-slice definition is unchanged (kn 0..5 = pT-hat 8-14, 14-24, 24-40, 40-70, 70-125,
125-300 GeV); `kn` only chooses which `muon_pair_tree_kin<kn>_*` tree is filled and the output name.

**Guards (all throw before any output is written).** `truth_raw_counts` without
`truth_input_override` (would alter a nominal run); `truth_input_override` without
`truth_raw_counts` (would weight a new dataset with the central sample's AMI cross-sections —
forbidden by `docs/ami_weights.md`); either on the private/fullsim path; missing output dir or
label; a non-existent literal input file, a glob matching nothing, or 0 entries.

**How to run** (per slice):

```bash
source ~/setup.sh
cd NTupleProcessingCode
root -b -l <<'EOF'
.L PythiaAnalysisClasses.h
PythiaTruthAnalysis py(5, false, 5.36);        // kn = 5 (pTH125_300), non-private
py.truth_input_override      = "/path/to/new_sample/*.NTUP.root";
py.truth_output_dir_override = "/path/to/output_dir";   // must exist
py.truth_output_label        = "ss_pTH125_300_test";
py.truth_raw_counts          = true;
py.Run();
.q
EOF
```

Per-event rate of a category = (pairs in the category) / `nproc_kin<N>_beam0`. Validated
2026-09-23 on the 10-event pTH125_300 local test NTUP (`PythiaCentralProdJO/evgen_high_pT_additional/local_test/aod_test/myfile.root`):
3 SS + 11 OS pairs, identical to a hand count of `truth_mupair_*` with both legs pT > 4.5 GeV,
|η| < 2.4; all weights 1; `nproc_kin5_beam0 = 10`; ancestry tracing ran without error.
