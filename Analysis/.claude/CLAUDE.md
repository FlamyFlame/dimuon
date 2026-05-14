- Remote cluster (BNL SDCC); run `/usatlas/u/yuhanguo/setup.sh` for ROOT environment
- No pip install — use what's in the release
- For C++ class testing, use ROOT ACLiC: `.L MyClass.cxx+`
- Read `README.md` and `docs/` for analysis context (class hierarchy, pipelines, sample types)
- For any analysis change: always update and maintain the relevant documentation in `README.md` / `docs/`

## Tracking Documents

Create a tracking doc when: (a) investigating an unknown root cause with
multiple hypotheses, or (b) the user requests documentation for a new
analysis or major rework spanning many editing cycles.

**Start:** Create `docs/tracking/<name>.md`. Register in Active Tracking
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

## Active Tracking Docs