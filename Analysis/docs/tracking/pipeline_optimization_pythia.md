# Pipeline Optimization: Pythia Group

## Objective

Review and improve the two Pythia pipeline scripts to meet all 6 quality
dimensions. `pipeline_pythia_truth.sh` produces generator-level baseline
results; `pipeline_pythia_fullsim_overlay.sh` adds detector simulation with
HIJING overlay. Truth serves as the reference for validating fullsim overlay.

## Physics Procedure

### Motivation

The Pythia pipelines produce the MC reference distributions used to validate
detector effects and extract corrections. Truth-level provides the generator
baseline; fullsim overlay includes detector response and heavy-ion environment.

### What the pipelines produce

1. **Truth pipeline**: generator-level dimuon distributions (mass, pT, eta,
   rapidity) for Pythia b-bbar samples.
2. **Fullsim overlay pipeline**: reconstructed dimuon distributions after
   detector simulation and HIJING overlay, for comparison with truth and data.

### Review criteria

The 6-dimension checklist in `.claude/commands/review-pipeline.md` is the
authoritative reference. Key cross-pipeline checks:
- Truth and overlay pipelines should use consistent sample definitions and
  binning where they produce comparable distributions.
- Overlay pipeline must handle the HIJING BC structure correctly (generator
  bc 1-47k, Geant4 bc>200k).

### Negative constraints

- Do NOT modify the C++ analysis code or RDF classes during pipeline review.
  Flag issues for `/review-analysis-code` instead.

## Scope

- `Analysis/pipelines/pipeline_pythia_truth.sh`
- `Analysis/pipelines/pipeline_pythia_fullsim_overlay.sh`
- Reference docs: `Analysis/docs/pythia_truth.md`, `Analysis/docs/pythia_fullsim_overlay.md`

## Design Decisions

(none yet)

## Implementation Plan

1. Read both pipeline scripts and produce stage maps — per review-pipeline Step 1
2. Spawn reviewer on both scripts — per review-pipeline Step 2
3. Parse verdict, amend if needed — iterate until PASS or max iterations
4. Update this tracking doc with final results

**Status**: Complete

## Progress Log

### Step 1-2: Initial review (Iteration 1)
- Read both pipeline scripts, produced stage maps
- Truth pipeline: 7 stages (condor submit, wait, validate batch, hadd, RDF, validate RDF, plot)
- Overlay pipeline: 5 stages (condor submit/dry-run, validate NTP, RDF, validate hist, plot)
- Reviewed all 29 items across 6 dimensions
- Found 2 issues: (1) truth pipeline missing SKIP_CONDOR, (2) docs/pythia_truth.md wrong .sub file for 5.02

### Step 3: Amend + re-review (Iteration 2)
- Added SKIP_CONDOR env var to pipeline_pythia_truth.sh (lines 21, 37, 344-356)
- Updated docs/pythia_truth.md sample modes table for nonprivate_5p02
- Re-reviewed: all 29 items pass

## Results & Observations

- Both pipelines are well-structured with proper `set -Eeuo pipefail`, ERR traps,
  per-stage validation, and configurable timeouts.
- Helper functions (validate_root_file_quick, wait_for_cluster_completion, etc.) are
  consistent across both scripts and match the pattern used in other pipelines.
- Truth pipeline delegates RDF to `run_rdf_pythia_truth.sh`; overlay compiles+runs
  RDF inline. Both approaches are valid for their use cases.
- Overlay's `--dry-run` mode correctly skips Condor but runs remaining stages.
- Cross-pipeline consistency: both use the same data types (MuonPairPythia variants),
  sign conventions (sign1=SS, sign2=OS), and validation patterns.
- No C++ bugs found at the orchestration level.

## Remaining Work

(none)

## Latest Stage

Complete. Review approved at iteration 2 with 0 remaining issues.
