---
name: steer-pipeline
description: Preflight, execute, monitor, recover, validate, and document dimuon analysis pipelines. Use when the user asks to run, steer, resume, or autonomously monitor a multi-stage NTuple, RDataFrame, fitting, plotting, Condor, grid, or end-to-end analysis pipeline.
---

# Steer an analysis pipeline

1. Confirm the requested script, arguments, environment, starting stage, and
   authority for state-changing cluster/grid operations.
2. Read `.claude/commands/steer-pipeline.md` as the latest read-only procedure,
   then read the full script, its pipeline documentation, active tracking
   document, `Analysis/README.md`, and relevant KB entries. Translate all logs to
   `.codex/logs/`; never write `.claude/logs/`.
3. Before execution, record in the tracking document and run log:
   - task, script, arguments, timestamp, Git SHA, and dirty files
   - input paths and expected outputs
   - dependency/source changes since the last successful run
   - stages to run/skip and validation for each
4. Do a dry-run or syntax/preflight pass where supported. Verify executables,
   environment setup, disk paths, input freshness, and output collision policy.
5. Run with captured logs. Monitor until terminal completion; do not infer
   success from a submitted job. Record job IDs and recovery attempts.
6. On failure, preserve evidence, identify the failed stage, apply only an
   authorized minimal recovery, and resume from the earliest invalid stage.
7. Validate output existence, timestamps, object counts, warnings, numerical
   sanity, and plots. Invoke `review-analysis-code` for source fixes and
   `review-plot` for generated visual outputs.
8. Update tracking status, run log, downstream rerun obligations, and final
   limitations. Commit or push only when explicitly requested.

Never start grid submission, Rucio rule mutation, or `pbook` mutation without
explicit user authority.
