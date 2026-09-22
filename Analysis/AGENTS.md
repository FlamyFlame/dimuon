# Analysis-specific Codex guidance

This file extends the repository-root `AGENTS.md` for work under `Analysis/`.

## Before changing analysis artifacts

1. Read `docs/analysis_overview.md`, `README.md`, the relevant pipeline doc, and
   all active tracking documents listed by the root context hook.
2. For an implementation tracking document, re-read `Objective`, `Physics
   Procedure`, `Remaining Work`, and `Latest Stage` before the step.
3. Update the tracking plan before work and verify that the planned variables,
   cuts, weights, binning, sample roles, and normalization follow the Physics
   Procedure. Stop and ask if they conflict.

Do not copy a code pattern merely because another pipeline looks similar.
Require the physics procedure to justify reuse. Derive names from physical
quantities, not implementation convenience.

## Required review lanes

- C++/ROOT/RDataFrame changes: use `review-analysis-code`.
- Plots, fits, labels, or visual outputs: use `review-plot`.
- Discrepancies and unexpected results: use `review-investigation`.
- Multi-stage pipeline execution: use `steer-pipeline`.
- Internal-note or paper work: use `academic-writing-gates`.

These skills provide Codex-native checklists. A separate reviewer agent is used
only when the user explicitly requests subagents; otherwise perform the full
checklist in the main thread and state that it was a self-review.

## Scientific verification

- Compile changed ROOT/C++ classes with the project's ACLiC pattern when the
  environment supports ROOT.
- Verify numerical outputs from source ROOT files or logs independently of the
  code that produced them when a result will be quoted.
- Review generated plots visually and check normalization, units, axes,
  selections, centrality, year/sample labels, uncertainty treatment, and output
  provenance.
- Never promote placeholder reconstruction efficiency, TAA, luminosity,
  hadronic cross section, trigger configuration, or preliminary output to a
  final result.
- Append validation evidence and exact paths to the active tracking document
  after each step.
