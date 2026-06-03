# Pipeline Optimization: PbPb Group

## Objective

Review and improve the two PbPb analysis pipeline scripts to meet all 6
quality dimensions (stage completeness, safety, idempotency, environment,
documentation, physics integration). These pipelines are tightly coupled:
trig eff outputs feed into the crossx pipeline.

## Physics Procedure

### Motivation

The PbPb pipelines orchestrate the full data analysis chain for heavy-ion
dimuon measurements. `pipeline_pbpb_trig_eff.sh` produces single-muon trigger
efficiency corrections (Pipelines 2+3: no-correlation efficiency via Fermi+log
pT fitting, then dR corrections via inverse weighting). `pipeline_pbpb_crossx.sh`
produces the cross-section measurement using nominal NTuple processing and
RDF hist filling.

### What the pipelines produce

1. **Trig eff pipeline**: per-year TF1 pT turn-on fits, TH2D fallback
   efficiency maps, Pipeline 2 efficiency plots, Pipeline 3 dR correction
   plots with plateau normalization.
2. **Crossx pipeline**: per-year RDF histograms for cross-section extraction,
   combined cross-section plots.

### Review criteria

The 6-dimension checklist in `.claude/commands/review-pipeline.md` is the
authoritative reference. Key cross-pipeline checks:
- Trig eff output filenames must match what crossx (or downstream analysis)
  expects as input.
- Both pipelines must use consistent QUEUE_COUNTS, year lists, trigger_mode,
  mindR settings, and sign conventions.
- Shared helper functions (validate_root_file_quick, source_env_once, etc.)
  should be identical or factored out.

### Negative constraints

- Do NOT modify the C++ analysis code or RDF classes during pipeline review.
  Flag issues for `/review-analysis-code` instead.
- Do NOT change physics parameters (trigger_mode, mindR, binning) unless a
  cross-pipeline inconsistency is found.

## Scope

- `Analysis/pipelines/pipeline_pbpb_trig_eff.sh`
- `Analysis/pipelines/pipeline_pbpb_crossx.sh`
- Reference doc: `Analysis/docs/pbpb_pipelines.md`

## Design Decisions

(none yet)

## Implementation Plan

1. Read both pipeline scripts and produce stage maps — per review-pipeline Step 1
2. Spawn reviewer on both scripts — per review-pipeline Step 2
3. Parse verdict, amend if needed — iterate until PASS or max iterations
4. Update this tracking doc with final results

**Status**: DONE

## Progress Log

### 2026-06-02: Review iteration 1 — PASS

**Issues found and fixed:**

1. **CRITICAL: `root_rc=$?` unreachable under `set -Eeuo pipefail`** (trig eff pipeline,
   4 call sites: Stages 5, 6, 7, 9). Under `set -e`, a failing `root` heredoc causes
   immediate shell exit via ERR trap before `root_rc=$?` is reached. `popd` also
   unreachable in the failure path, leaving directory stack unbalanced (cosmetic since
   script exits). Fixed by using `root_rc=0; root ... <<EOF || root_rc=$?` pattern which
   puts the command in a conditional context, suppressing `set -e`.

2. **WARNING: Crossx pipeline stage numbering gap** (header lists stages 1-7 but code
   jumps from "Stage 5" to "Stage 7"). Fixed by merging header stages 5+6 into
   "Stage 5: run RDF crossx hist filling per year + validate RDF outputs" and renumbering
   "Stage 7" to "Stage 6". Updated `docs/pbpb_pipelines.md` to match.

**Files modified:**
- `pipelines/pipeline_pbpb_trig_eff.sh` — 4 heredoc call sites fixed
- `pipelines/pipeline_pbpb_crossx.sh` — header + inline stage comment
- `docs/pbpb_pipelines.md` — stage numbering alignment

**All 29 review criteria pass across 6 dimensions.** See log file for details.

## Results & Observations

### Cross-pipeline consistency
- **QUEUE_COUNTS**: identical `[23]=4 [24]=2 [25]=6` in both scripts, matches all .sub files.
- **trigger_mode**: both use `1` (single mu4). Crossx adds `mu4_nominal_pbpb_NO_trig_calc=true`.
- **mindR_trig**: both use `0.02`.
- **Filename patterns**: trig eff uses `_res_cut_v2` suffix, crossx uses no suffix. Consistent with NTuple run scripts (`force_nominal=true` skips resonance veto).
- **DATA_BASE**: identical hardcoded path.
- **Shared helpers**: `source_env_once`, `validate_root_file_quick`, `validate_files_or_fail`, `validate_combined_muon_pair_trees_nonempty_or_fail`, `extract_cluster_id`, `wait_for_cluster_completion` are duplicated verbatim. Could be factored to a shared file but not flagged as an issue.

### `set -e` + heredoc pattern
The correct pattern for capturing `root` heredoc exit codes under `set -e` is:
```bash
root_rc=0
root -l -b <<EOF || root_rc=$?
...
EOF
if [[ $root_rc -ne 0 ]]; then fail "..."; fi
```
The `|| root_rc=$?` makes the command part of a compound list that is being tested,
so `set -e` does not trigger. The initialization `root_rc=0` is needed because
`|| root_rc=$?` only fires on failure.

## Remaining Work

(none)

## Latest Stage

**Review complete** (2026-06-02). All 29 criteria pass. 2 CRITICAL + 1 WARNING fixed.
