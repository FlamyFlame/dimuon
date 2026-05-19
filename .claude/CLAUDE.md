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
- *Implementation:* Objective, **Physics Procedure** (REQUIRED — see below),
  Context, Scope, Design Decisions (with rationale), Implementation Plan
  (numbered, with status), Progress Log (append-only), Results &
  Observations (organized, mutable), Remaining Work, Latest Stage.

### Physics Procedure requirement (implementation docs)

Every implementation doc MUST have a **Physics Procedure** section
immediately after Objective. This is the most important section in the
document — all design decisions, implementation choices, variable names,
and filter conditions must follow from it.

Contents:
1. **Motivation:** Why this measurement/correction is needed.
2. **Top-level equation:** The final result that will be applied to the
   analysis. Define every symbol.
3. **Step-by-step method:** Each pipeline/step with (a) what physical
   quantity it measures, (b) the mathematical procedure, (c) explicit
   statements of what is and is NOT required (e.g., "no trigger
   requirement on the other muon").
4. **Negative constraints:** Explicitly state what the code must NOT do
   where confusion with similar-but-different procedures is likely (e.g.,
   "this differs from mu4_mu4noL1 where the other muon must pass mu4").

This section is the **authoritative reference**. If a code spec, design
decision, or task prompt contradicts the Physics Procedure, flag to the
user before proceeding.

### Implementation doc sections (not investigation)

- **Results & Observations** (replaces Accumulated Findings): organize by
  topic, combine related items, delete obsolete entries when issues are
  fixed. Not append-only.
- **No Ruled Out section** in implementation docs. (Keep Ruled Out for
  investigation docs only.)

### Before each step

Update Latest Stage with current step, plan, and files involved. Write
BEFORE doing any work. **Verify the planned step is consistent with the
Physics Procedure section.** Each implementation plan step should specify
which reviewer command applies (e.g., /review-analysis-code, /review-plot)
and which Physics Procedure sections to include in the task prompt.

### After each step

Append results to Progress Log (implementation) or Accumulated Findings
(investigation) with step number. Keep exact values: paths, line numbers,
names, numbers. For investigations, add ruled-out approaches/hypotheses to
Ruled Out with reason. For implementations, mark step done and update
Remaining Work.

### Design changes (implementation)

Record old approach, new approach, and reason in Design Decisions before
proceeding. When a physics motivation applies, reference the Physics
Procedure section and include a concise physics reason (not just a code
rationale). Some design decisions are technical (e.g., code refactoring,
structure changes) — these need only a technical rationale. Verify the
new approach is consistent with the Physics Procedure. If it contradicts
the procedure, update the procedure first (with user approval).

### Continuity (CRITICAL)

At the start of every new conversation, and before resuming after any
context compression, check Active Tracking Docs. If any exist, Read each
fully before doing anything else. **For implementation docs, re-read the
Physics Procedure section first — it is the authoritative reference for
what the code should do.** The doc is ground truth — if conversation
history or compaction summaries conflict, trust the doc.

### Anti-drift rules

- Implementation Plan steps must reference the Physics Procedure section
  they implement (e.g., "per §3a"). This cross-reference is mandatory.
- Naming conventions in code (variable names, histogram suffixes, function
  names) must derive from the physics terminology in the Physics Procedure,
  not from code convenience or other pipeline patterns.
- If a code pattern from one pipeline is being reused in another, explicitly
  verify in the Physics Procedure that the physics justifies reuse.
  Different physics → different code, even if the structure looks similar.
- When writing task prompts for /review-analysis-code or similar, include
  the relevant Physics Procedure sections. The reviewer must check physics
  correctness against the procedure, not just code correctness against a
  code spec.

### Mandatory reviewer loops for implementation plan steps

The Auto-Dispatch Rules (above) apply to user requests. For **agent-
initiated** work within an implementation plan, the same rules apply:
- Steps that write or modify C++/ROOT/RDF code → /review-analysis-code.
  Include the relevant Physics Procedure section(s) in the task prompt.
- Steps that create or modify plots or fitting output → /review-plot.
- This applies even after compaction — check the implementation plan
  for which reviewer each step requires.

### Completion

Write final summary, clear Latest Stage, remove from Active Tracking Docs.

### Never

Keep findings/progress only in conversation (write to doc before next
action); start a step without writing plan to doc first; declare complete
without re-reading Objective and Physics Procedure; implement code that
contradicts the Physics Procedure without user approval.

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

- `Analysis/docs/tracking/mu4_trig_effcy_implementation.md` — MU4 trigger efficiency: 3-pipeline separation, no-corr efficiency, inverse-weighted dR corrections
- `Analysis/docs/tracking/mu4_trig_effcy_investigation.md` — Investigation: Pipeline 3 physics procedure bugs, plot directory structure, isBNL toggle