# Pipeline Optimization: POWHEG Group

## Objective

Review and improve the two POWHEG fullsim pipeline scripts to meet all 6
quality dimensions. `pipeline_powheg_fullsim_single_muon.sh` processes
single-muon efficiency studies; `pipeline_powheg_fullsim_mixed_pairs.sh`
processes mixed dimuon pair distributions. Both use the same POWHEG generator
with full detector simulation.

## Physics Procedure

### Motivation

The POWHEG pipelines produce NLO-accurate MC predictions for b-bbar dimuon
production. Single-muon pipeline extracts per-muon efficiencies; mixed-pairs
pipeline produces dimuon distributions for comparison with data and Pythia.

### What the pipelines produce

1. **Single muon pipeline**: single-muon efficiency distributions from POWHEG
   fullsim samples.
2. **Mixed pairs pipeline**: dimuon pair distributions (mass, pT, dR, etc.)
   from POWHEG fullsim.

### Review criteria

The 6-dimension checklist in `.claude/commands/review-pipeline.md` is the
authoritative reference. Key cross-pipeline checks:
- Both pipelines should use consistent POWHEG sample paths and event selection.
- Trigger configurations should match the physics goal of each pipeline.

### Negative constraints

- Do NOT modify the C++ analysis code or RDF classes during pipeline review.
  Flag issues for `/review-analysis-code` instead.

## Scope

- `Analysis/pipelines/pipeline_powheg_fullsim_single_muon.sh`
- `Analysis/pipelines/pipeline_powheg_fullsim_mixed_pairs.sh`
- Reference doc: `Analysis/docs/powheg.md`

## Design Decisions

1. Added `--skip-condor` to single-muon pipeline for parity with mixed-pairs
   pipeline's `--skip-mixing` flag. Rationale: 22 Condor jobs is expensive;
   re-running post-Condor stages (hadd, RDF, plot) should not require
   re-submitting jobs.

## Implementation Plan

1. [DONE] Read both pipeline scripts and produce stage maps -- per review-pipeline Step 1
2. [DONE] Review both scripts against all 29 criteria across 6 dimensions
3. [DONE] Fix issues: added `--skip-condor` flag and Usage section to single-muon pipeline
4. [DONE] Re-review: all 29 criteria pass

**Status**: Complete

## Progress Log

### Step 1-2: Initial review (2026-06-02)
- Read both pipeline scripts, all Condor .sub and .sh wrapper files
- Produced stage maps for both pipelines
- Reviewed all 29 criteria across 6 dimensions
- Found 2 WARNING issues in single-muon pipeline:
  - Missing `--skip-condor` flag (Dim 3, item 14)
  - Missing Usage section in header (Dim 5, item 22/25)

### Step 3: Fixes applied (2026-06-02)
- Added `--skip-condor` flag parsing: lines 40-47 of single-muon pipeline
- Added Usage section: lines 17-24 of single-muon pipeline
- Wrapped Steps 1-2 in `if (( SKIP_CONDOR ))` conditional
- Added pipeline mode logging line
- Verified bash syntax: both scripts pass `bash -n`

### Step 4: Re-review (2026-06-02)
- All 29 criteria pass across both pipelines
- Cross-pipeline consistency verified:
  - Sample paths identical (POWHEG_BASE, BB_DIR, CC_DIR)
  - Combined tree paths match (BB_COMBINED, CC_COMBINED)
  - Helper functions identical (now, log, fail, on_error, source_env_once, etc.)

## Results & Observations

### Cross-pipeline consistency
- Both pipelines use identical POWHEG_BASE, BB_DIR, CC_DIR paths
- Combined single-muon tree paths (BB_COMBINED, CC_COMBINED) are identical
  between both pipelines, confirming the dependency chain
- Helper functions (now, log, fail, on_error, require_cmd, source_env_once,
  extract_cluster_id, validate_root_file_quick) are identical across both scripts

### Quality assessment
- Both pipelines follow the same patterns as `pipeline_pythia_fullsim_overlay.sh`
- Error handling is thorough: ERR trap, per-stage validation, held-job detection
- Mixed-pairs pipeline has additional sophistication: lone-job timeout, skip-batch
  propagation to RDF, smoke testing, mass-filter mode selection
- No C++ issues flagged for `/review-analysis-code`

## Remaining Work

None. Review complete.

## Latest Stage

Complete. Review passed at iteration 1 after 2 WARNING fixes.
