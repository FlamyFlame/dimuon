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
approaches to Ruled Out with reason. For implementations, mark step done
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

## Agent Dispatch

When working on plotting: read `.claude/agents/plot-reviewer.md` + `.claude/conventions/atlas-plotting.md`
When working on analysis code (C++/RDF): read `.claude/agents/code-reviewer.md` + `.claude/conventions/codebase-patterns.md` + `.claude/conventions/rdf-patterns.md`
When working on physics / event selection: read `.claude/agents/physics-reviewer.md`
When verifying numbers or yields: read `.claude/agents/numerical-verifier.md`
When writing internal notes or papers: read `.claude/agents/note-reviewer.md`

## Available Commands

- `/review-plot` — Executor-reviewer loop for creating/fixing plots
- `/review-analysis-code` — Executor-reviewer loop for analysis code changes
- `/review-investigation` — Executor-reviewer loop for physics investigations
- `/review-note` — Executor-reviewer loop for internal note writing
- `/review-paper` — Executor-reviewer loop for paper writing

## Agent Logs

All executor-reviewer loops write logs to `.claude/logs/`. Review to understand what was attempted and what issues arose.

## Active Tracking Docs