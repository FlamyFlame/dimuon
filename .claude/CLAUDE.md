- Remote cluster (BNL SDCC); run `/usatlas/u/yuhanguo/setup.sh` for ROOT environment
- No pip install — use what's in the release
- For C++ class testing, use ROOT ACLiC: `.L MyClass.cxx+`
- Read `/usatlas/u/yuhanguo/workarea/dimuon_codes/Analysis/README.md` and `/usatlas/u/yuhanguo/workarea/dimuon_codes/Analysis/docs/` for analysis context (class hierarchy, pipelines, sample types)
- For any analysis change: always update and maintain the relevant documentation in those files

## Tracking Documents

Create a tracking doc when: (a) investigating an unknown root cause with
multiple hypotheses, or (b) the user requests documentation for a new
analysis or major rework spanning many editing cycles.

**Start:** Create `/usatlas/u/yuhanguo/workarea/dimuon_codes/Analysis/docs/tracking/<name>.md`. Register in Active Tracking
Docs below. Structure depends on mode:
- *Investigation:* Objective, Context, Sub-steps, Accumulated Findings
  (append-only), Ruled Out (append-only), Latest Stage.
- *Implementation:* Objective, Scope, Design Decisions (with rationale),
  Implementation Plan (numbered, with status), Progress Log (append-only),
  Remaining Work, Latest Stage.

**Before each step:** update Latest Stage with current step, plan, and
files involved. Write BEFORE doing any work.

**After each step:** append results to Accumulated Findings (investigation)
or Progress Log (implementation) with step number. Keep exact values:
paths, line numbers, names, numbers. For investigations, add ruled-out
approaches/hypotheses to Ruled Out with reason. For implementations, mark step done
and update Remaining Work.

**Design changes (implementation):** record old approach, new approach,
and reason in Design Decisions before proceeding.

**Continuity (CRITICAL):** At the start of every new conversation, and
before resuming after any context compression, check Active Tracking Docs.
If any exist, Read each fully before doing anything else. The doc is
ground truth — if conversation history conflicts, trust the doc.

**Completion:** Write final summary, clear Latest Stage, remove from
Active Tracking Docs.

**Never:** keep findings/progress only in conversation (write to doc
before next action); start a step without writing plan to doc first;
declare complete without re-reading Objective; delete from append-only
sections; revisit Ruled Out without a new reason.

## Documentation References

Before working on any task, check these existing docs:
- Knowledge base: `.claude/kb/index.md` (analysis overview, decisions, samples, variables, gotchas)
- Analysis overview, class hierarchy, pipelines: `Analysis/README.md`
- Per-pipeline docs: `Analysis/docs/` (pythia_truth, pythia_fullsim_pp, pythia_fullsim_overlay, powheg, data_analysis)
- Skimming code, branches, grid workflow: `SkimCode/README.md`
- Data paths and directory layout: root `README.md`
- Internal note structure: `IntNote/tex/` (section files), `IntNote/mydocument.tex` (master)

## Auto-Dispatch Rules

When the user asks to create, fix, or modify a plot → invoke `/review-plot`
When the user asks to write, modify, or fix C++/ROOT/RDF analysis code → invoke `/review-analysis-code`
When the user asks to investigate a discrepancy, debug, or understand an unexpected result → invoke `/review-investigation`
When the user asks to write or edit an internal note section → invoke `/review-note`
When the user asks to write or polish paper text for publication → invoke `/review-paper`

Each command runs an executor-reviewer loop: the main agent executes, then
spawns a reviewer subagent (via Agent tool) for independent evaluation. The
reviewer returns a structured VERDICT: PASS/FAIL. Loop continues until PASS
or MAX_ITERATIONS (5).

## Agent Reference Files

These files define review criteria. They are embedded in subagent prompts by
the commands above — you do not need to read them directly unless working
outside of a command loop.

- `.claude/agents/executor.md` — executor behavior
- `.claude/agents/plot-reviewer.md` — plot review criteria
- `.claude/agents/code-reviewer.md` — code review criteria
- `.claude/agents/physics-reviewer.md` — physics correctness criteria
- `.claude/agents/numerical-verifier.md` — number verification procedure
- `.claude/agents/note-reviewer.md` — note/paper review criteria
- `.claude/conventions/atlas-plotting.md` — plotting conventions
- `.claude/conventions/codebase-patterns.md` — codebase conventions
- `.claude/conventions/rdf-patterns.md` — RDataFrame conventions

## Agent Logs

All executor-reviewer loops write logs to `.claude/logs/`. Log files are named
`review-<type>-YYYYMMDD-HHMMSS-<task-slug>.md`. The active session's log path
is also written to `.claude/logs/active-session.txt` for compaction recovery.

## Active Tracking Docs